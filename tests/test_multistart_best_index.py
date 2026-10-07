"""Regression tests for multi-start result reporting.

Bug (fixed): in parallel mode the winning restart was identified by the
position of its error inside ``all_errors``, which is filled in *completion*
order (``as_completed``) and interleaved with ``None`` for failed fits. The
reported ``best_start_index`` therefore did not correspond to the restart that
actually produced the best fit. The fix tracks ``start_idx`` in a list aligned
with ``all_results``.
"""

import time

import numpy as np
import pytest

import eis_analysis.fitting.multistart as ms
from eis_analysis.fitting.circuit import FitResult


class _FakeCircuit:
    """Minimal R|C circuit stub with an impedance off the (unit) data by the
    controlled error, so a start's weighted RSS - the multistart selection
    criterion - ranks the same as its error."""

    err = 0.0

    def get_param_labels(self):
        return ['R', 'C']

    def update_params(self, params):
        return 2

    def impedance(self, freq, params):
        return np.full(len(freq), 1.0 + self.err, dtype=complex)


def _install_deterministic_mocks(monkeypatch, errors, sleeps):
    """Patch perturbation + fit so each restart idx is identifiable.

    The perturbation encodes the restart index in ``params[0]`` (perturbations
    are generated sequentially in the main thread, so a simple counter is
    safe). The fit mock then reads that index back to assign a controlled
    error and an optional sleep to force a particular completion order.
    """
    base_params = np.array([100.0, 1e-6])
    idx_counter = {"i": 2}  # perturbed restarts are numbered 2, 3, ...

    def fake_perturb(params, factor=3.0, bounds=None, rng=None):
        i = idx_counter["i"]
        idx_counter["i"] += 1
        return np.array([float(i), 1e-6])

    monkeypatch.setattr(ms, "perturb_log_uniform", fake_perturb)

    def fake_fit(frequencies, Zdata, circuit, weighting=None,
                 initial_guess=None, use_analytic_jacobian=True):
        if initial_guess is None:
            idx = 1  # the initial fit is restart #1
            params = base_params
        else:
            idx = int(round(initial_guess[0]))
            params = np.asarray(initial_guess, dtype=float)
            time.sleep(sleeps.get(idx, 0.0))
        circuit.err = errors[idx]
        res = FitResult(
            circuit=circuit,
            params_opt=params,
            params_stderr=np.array([np.inf, np.inf]),  # forces log_uniform path
            fit_error_rel=errors[idx],
            cov=None,
            is_well_conditioned=False,
            n_free_params=2,
        )
        # fit_equivalent_circuit() writes its parameters into the shared
        # circuit object; mimic that so the circuit-sync test is meaningful.
        circuit.update_params(list(params))
        Z_fit = np.ones(len(frequencies), dtype=complex)
        return res, Z_fit

    monkeypatch.setattr(ms, "fit_equivalent_circuit", fake_fit)


def test_best_start_index_parallel_completion_order(monkeypatch):
    """Winning restart must be named correctly even when it finishes last."""
    freq = np.logspace(5, -1, 20)
    Z = np.ones_like(freq) + 0j

    # Restart #3 is the unique global best, but is forced to complete LAST.
    errors = {1: 0.50, 2: 0.40, 3: 0.10, 4: 0.30}
    sleeps = {2: 0.0, 3: 0.20, 4: 0.0}
    _install_deterministic_mocks(monkeypatch, errors, sleeps)

    result, _ = ms.fit_circuit_multistart(
        _FakeCircuit(), freq, Z, n_restarts=4, parallel=True, max_workers=4
    )

    assert result.diagnostics.best_start_index == 3
    assert result.diagnostics.best_error == pytest.approx(0.10)


def test_best_start_index_sequential(monkeypatch):
    """Sequential mode reports the correct winning restart."""
    freq = np.logspace(5, -1, 20)
    Z = np.ones_like(freq) + 0j

    errors = {1: 0.50, 2: 0.10, 3: 0.30, 4: 0.40}  # restart #2 wins
    _install_deterministic_mocks(monkeypatch, errors, sleeps={})

    result, _ = ms.fit_circuit_multistart(
        _FakeCircuit(), freq, Z, n_restarts=4, parallel=False
    )

    assert result.diagnostics.best_start_index == 2
    assert result.diagnostics.best_error == pytest.approx(0.10)


def test_best_start_index_initial_fit_wins(monkeypatch):
    """When the initial fit is best, restart #1 is reported."""
    freq = np.logspace(5, -1, 20)
    Z = np.ones_like(freq) + 0j

    errors = {1: 0.05, 2: 0.40, 3: 0.30, 4: 0.50}  # initial fit wins
    _install_deterministic_mocks(monkeypatch, errors, sleeps={})

    result, _ = ms.fit_circuit_multistart(
        _FakeCircuit(), freq, Z, n_restarts=4, parallel=True, max_workers=4
    )

    assert result.diagnostics.best_start_index == 1
    assert result.diagnostics.best_error == pytest.approx(0.05)


class _RecordingCircuit(_FakeCircuit):
    """Circuit stub that stores the parameters written into it, as the real
    Circuit does via update_params()."""

    def __init__(self):
        self.params = None

    def update_params(self, params):
        self.params = list(params)
        return len(self.params)


def test_best_result_circuit_holds_best_params(monkeypatch):
    """best_result.circuit must carry the BEST parameters, not the last fit's.

    Bug (fixed): every restart fits the same Circuit object, and each fit
    writes its own parameters into it, so the circuit ended up holding the
    parameters of the restart that happened to run last. Consumers reading
    parameters from the circuit tree (oxide capacitance -> thickness /
    permittivity) then used the wrong fit.
    """
    freq = np.logspace(5, -1, 20)
    Z = np.ones_like(freq) + 0j

    errors = {1: 0.50, 2: 0.10, 3: 0.30, 4: 0.40}  # restart #2 wins, #4 runs last
    _install_deterministic_mocks(monkeypatch, errors, sleeps={})

    circuit = _RecordingCircuit()
    result, _ = ms.fit_circuit_multistart(
        circuit, freq, Z, n_restarts=4, parallel=False
    )

    assert result.best_result.params_opt[0] == pytest.approx(2.0)
    assert circuit.params[0] == pytest.approx(2.0)  # not 4.0 (the last restart)
    assert result.best_result.circuit.params == circuit.params


def test_parallel_restarts_do_not_share_circuit(monkeypatch):
    """Each restart must fit its own circuit copy.

    Bug (fixed): all restarts fitted the one circuit object the caller passed
    in, and each fit writes its parameters into it - so in parallel mode
    several threads wrote (and impedance() read) the same parameters at once,
    and every result's circuit ended up holding whatever ran last.
    """
    freq = np.logspace(5, -1, 20)
    Z = np.ones_like(freq) + 0j

    errors = {1: 0.50, 2: 0.30, 3: 0.10, 4: 0.40}  # restart #3 wins
    sleeps = {2: 0.0, 3: 0.15, 4: 0.0}             # ... and completes last
    _install_deterministic_mocks(monkeypatch, errors, sleeps)

    circuit = _RecordingCircuit()
    result, _ = ms.fit_circuit_multistart(
        circuit, freq, Z, n_restarts=4, parallel=True, max_workers=4
    )

    circuits = [r.circuit for r in result.all_results]
    assert len({id(c) for c in circuits}) == len(circuits)  # no shared object

    # Each restart's circuit holds that restart's own parameters
    for r in result.all_results:
        assert r.circuit.params[0] == pytest.approx(r.params_opt[0])

    assert result.best_result.params_opt[0] == pytest.approx(3.0)
    assert circuit.params[0] == pytest.approx(3.0)  # caller's circuit = best fit


def test_rng_makes_restarts_reproducible():
    """The same seed gives the same restarts; the global np.random is untouched.

    Regression: the perturbations drew from the global np.random, so a fit was
    not reproducible and depended on whatever else had consumed that stream.
    """
    from eis_analysis.cli.utils import parse_circuit_expression

    f = np.logspace(5, -1, 40)
    truth = parse_circuit_expression("R(10) - (R(1000)|C(1e-6)) - (R(500)|Q(1e-4,0.8))")
    noise = np.random.default_rng(0).standard_normal((2, f.size))
    Z = truth.impedance(f, truth.get_all_params())
    Z = Z + 0.01 * np.abs(Z) * (noise[0] + 1j * noise[1])

    def restarts(rng):
        circuit = parse_circuit_expression("R(20) - (R(2000)|C(3e-6)) - (R(200)|Q(3e-5,0.7))")
        result, _ = ms.fit_circuit_multistart(circuit, f, Z, n_restarts=4, rng=rng)
        return np.array([r.params_opt for r in result.all_results])

    np.random.seed(123)
    state = np.random.get_state()[1].copy()
    first = restarts(7)
    assert np.array_equal(np.random.get_state()[1], state)
    assert np.array_equal(first, restarts(np.random.default_rng(7)))
    assert not np.array_equal(first, restarts(8))


@pytest.mark.parametrize('perturb', ['stderr', 'covariance', 'log_uniform'])
def test_perturbed_conductance_keeps_its_magnitude(perturb):
    """A G of 1e-17 S (bounds 0..1e4 S) stays near 1e-17, never below 0.

    Regression: every perturbation ended in max(x, 1e-15), which lifted an
    oxide-film conductance below 1e-15 S up to 1e-15 and G = 0 off its
    valid lower bound.
    """
    G = np.array([1e-17])
    bounds = (np.array([0.0]), np.array([1e4]))
    rng = np.random.default_rng(0)
    if perturb == 'stderr':
        draws = [ms.perturb_from_stderr(G, np.array([2e-18]), bounds=bounds, rng=rng)
                 for _ in range(200)]
    elif perturb == 'covariance':
        draws = [ms.perturb_from_covariance(G, np.array([[4e-36]]), bounds=bounds, rng=rng)
                 for _ in range(200)]
    else:
        draws = [ms.perturb_log_uniform(G, bounds=bounds, rng=rng) for _ in range(200)]
    draws = np.concatenate(draws)
    assert np.all(draws >= 0)
    assert np.max(draws) < 1e-16


def test_perturbation_reaches_zero_but_not_below():
    """With and without bounds, a perturbation crossing 0 stops at exactly 0."""
    rng = np.random.default_rng(0)
    G, stderr = np.array([1e-17]), np.array([1e-16])   # most draws cross 0
    for bounds in [(np.array([0.0]), np.array([1e4])), None]:
        draws = np.concatenate([ms.perturb_from_stderr(G, stderr, bounds=bounds, rng=rng)
                                for _ in range(200)])
        assert np.all(draws >= 0)
        assert np.any(draws == 0)

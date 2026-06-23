"""Regression tests for multi-start best_start_index reporting.

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
    """Minimal circuit stub: no real parameter labels, trivial impedance."""

    def get_param_labels(self):
        return None

    def impedance(self, freq, params):
        return np.ones(len(freq), dtype=complex)


def _install_deterministic_mocks(monkeypatch, errors, sleeps):
    """Patch perturbation + fit so each restart idx is identifiable.

    The perturbation encodes the restart index in ``params[0]`` (perturbations
    are generated sequentially in the main thread, so a simple counter is
    safe). The fit mock then reads that index back to assign a controlled
    error and an optional sleep to force a particular completion order.
    """
    base_params = np.array([100.0, 1e-6])
    idx_counter = {"i": 2}  # perturbed restarts are numbered 2, 3, ...

    def fake_perturb(params, factor=3.0, bounds=None):
        i = idx_counter["i"]
        idx_counter["i"] += 1
        return np.array([float(i), 1e-6])

    monkeypatch.setattr(ms, "perturb_log_uniform", fake_perturb)

    def fake_fit(frequencies, Zdata, circuit, weighting=None,
                 initial_guess=None, plot=False, use_analytic_jacobian=True):
        if initial_guess is None:
            idx = 1  # the initial fit is restart #1
            params = base_params
        else:
            idx = int(round(initial_guess[0]))
            params = np.asarray(initial_guess, dtype=float)
            time.sleep(sleeps.get(idx, 0.0))
        res = FitResult(
            circuit=circuit,
            params_opt=params,
            params_stderr=np.array([np.inf, np.inf]),  # forces log_uniform path
            fit_error_rel=errors[idx],
            cov=None,
            is_well_conditioned=False,
        )
        Z_fit = np.ones(len(frequencies), dtype=complex)
        return res, Z_fit, None

    monkeypatch.setattr(ms, "fit_equivalent_circuit", fake_fit)


def test_best_start_index_parallel_completion_order(monkeypatch):
    """Winning restart must be named correctly even when it finishes last."""
    freq = np.logspace(5, -1, 20)
    Z = np.ones_like(freq) + 0j

    # Restart #3 is the unique global best, but is forced to complete LAST.
    errors = {1: 0.50, 2: 0.40, 3: 0.10, 4: 0.30}
    sleeps = {2: 0.0, 3: 0.20, 4: 0.0}
    _install_deterministic_mocks(monkeypatch, errors, sleeps)

    result, _, _ = ms.fit_circuit_multistart(
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

    result, _, _ = ms.fit_circuit_multistart(
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

    result, _, _ = ms.fit_circuit_multistart(
        _FakeCircuit(), freq, Z, n_restarts=4, parallel=True, max_workers=4
    )

    assert result.diagnostics.best_start_index == 1
    assert result.diagnostics.best_error == pytest.approx(0.05)

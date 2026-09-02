"""Test G element (conductance, Y = G).

G exists so that a parallel resistance the data cannot determine from above
lands on a regular interior point (G = 0) instead of R's artificial upper
bound. The tests below pin the three things that has to mean:

1. G is electrically R(1/G) - same impedance, same fit quality.
2. G = 0 stays finite everywhere (impedance, parallel combination, Jacobian).
3. A fit of an unresolvably large parallel resistance reports an interpretable
   G = (0 +- s) S instead of a spurious-precision R at its bound.
"""

import numpy as np
import pytest

from eis_analysis.fitting import C, G, Q, R, fit_equivalent_circuit
from eis_analysis.fitting.bounds import (PARAMETER_BOUNDS, classify_bound_status,
                                         log_scale_ci_mask, log_search_bounds)
from eis_analysis.fitting.jacobian import circuit_jacobian, element_jacobian

FREQ = np.logspace(-2, 5, 60)


# --- 1. G is R by another name ---

def test_impedance_is_reciprocal_of_conductance():
    Z = G(2e-3).impedance(FREQ, [2e-3])
    assert np.allclose(Z, 500.0 + 0j)


def test_matches_equivalent_resistor_in_a_parallel_branch():
    circuit_G = R(10) - (G(1e-4) | Q(1e-6, 0.85))
    circuit_R = R(10) - (R(1e4) | Q(1e-6, 0.85))
    Z_G = circuit_G.impedance(FREQ, circuit_G.get_all_params())
    Z_R = circuit_R.impedance(FREQ, circuit_R.get_all_params())
    assert np.allclose(Z_G, Z_R, rtol=1e-12)


def test_string_argument_fixes_the_parameter():
    g = G("1e-9")
    assert g.fixed_params == [True]
    assert g.G == 1e-9
    assert repr(g) == 'G("1e-09")'


# --- 2. G = 0 is finite everywhere ---

def test_zero_conductance_gives_finite_impedance():
    Z = G(0.0).impedance(FREQ, [0.0])
    assert np.all(np.isfinite(Z))
    # 1/(0+0j) would be inf+nanj in numpy; the clip must avoid that entirely
    assert not np.any(np.isnan(Z))


def test_zero_conductance_vanishes_from_a_parallel_combination():
    """G = 0 is an open branch: (G(0) | C) must equal C alone."""
    with_G = G(0.0) | C(1e-6)
    Z_with_G = with_G.impedance(FREQ, with_G.get_all_params())
    Z_C_only = C(1e-6).impedance(FREQ, [1e-6])
    assert np.all(np.isfinite(Z_with_G))
    assert np.allclose(Z_with_G, Z_C_only, rtol=1e-12)


def test_zero_conductance_parallel_jacobian_is_minus_Z_squared():
    """dZ/dG = -Z_total^2 at G = 0, the point the clip has to get right."""
    circuit = G(0.0) | C(1e-6)
    params = circuit.get_all_params()
    Z, dZ = circuit_jacobian(circuit, FREQ, params)
    assert np.all(np.isfinite(dZ))
    assert np.allclose(dZ[:, 0], -Z ** 2, rtol=1e-8)


def test_analytic_jacobian_matches_numeric():
    g, h = 3e-4, 1e-9
    Z, dZ = element_jacobian(G(g), FREQ, [g])
    numeric = (G(g + h).impedance(FREQ, [g + h]) - Z) / h
    assert np.allclose(dZ[:, 0], numeric, rtol=1e-4)


# --- 3. bounds semantics: zero is a result, not a warning ---

def test_lower_bound_is_exactly_zero():
    """Anything above 0 would make G a log-scale parameter, undoing the point."""
    lo, hi = PARAMETER_BOUNDS['G']
    assert lo == 0.0
    assert hi == 1.0 / PARAMETER_BOUNDS['R'][0]


def test_small_conductance_is_not_flagged_as_at_bound():
    lo, hi = PARAMETER_BOUNDS['G']
    assert classify_bound_status(1e-10, lo, hi) == ''
    assert classify_bound_status(0.0, lo, hi) == ''
    # the upper bound is a real box edge and must still be reported
    assert classify_bound_status(hi * 0.999, lo, hi) == 'upper'


def test_nonzero_lower_bound_still_flagged():
    """The zero-lower-bound rule must not disable the diagnostic generally."""
    assert classify_bound_status(0.301, 0.3, 1.0) == 'lower'


# --- 4. DE must still search G over decades ---

def test_de_searches_G_in_log_space():
    """Regression: G's zero lower bound must not drop it out of DE's log search.

    DE draws its population uniformly over the search bounds. Linear on
    [0, 1e4], essentially every member shorts the branch, the energies
    collapse and DE stops after a handful of generations at ~100 % error.
    log_search_bounds exists to keep that from happening, so G must not
    follow log_scale_ci_mask here even though it does for the CI.
    """
    labels = ['R', 'G', 'C', 'n', 'α_CC']
    lo = [PARAMETER_BOUNDS[k][0] for k in labels]
    hi = [PARAMETER_BOUNDS[k][1] for k in labels]

    mask, search_lo, search_hi = log_search_bounds(lo, hi, labels)
    assert mask == [True, True, True, False, False]
    # G is searched from 1e-10 S = 1e10 Ohm, mirroring R's upper bound
    assert search_lo[1] == pytest.approx(-10.0)
    assert search_hi[1] == pytest.approx(4.0)
    # ...and it is the one parameter where search and CI semantics differ
    assert log_scale_ci_mask(lo, hi) == [True, False, True, False, False]


def test_de_converges_with_G_like_it_does_with_R():
    """The DE stage itself must do the work, not leave it all to refinement."""
    from eis_analysis.fitting import fit_circuit_diffevo

    rng = np.random.default_rng(1)
    truth = R(10) - (R(1e5) | C(1e-8)) - (R(1e3) | C(1e-5))
    Z0 = truth.impedance(FREQ, truth.get_all_params())
    Z = Z0 + 0.005 * np.abs(Z0) * (rng.standard_normal(len(FREQ))
                                   + 1j * rng.standard_normal(len(FREQ)))

    start = R(1) - (G(1) | C(1e-12)) - (G(1) | C(1e-12))
    result, _, _ = fit_circuit_diffevo(start, FREQ, Z, seed=42)

    assert [lab for lab in result.diagnostics.log_search_params
            if lab.startswith('G')] == ['G0', 'G1']
    # Before the fix this was ~100 % after 4-5 generations
    assert result.de_error < 5.0
    assert result.diagnostics.de_iterations > 20


# --- 5. the payoff: an unresolvable parallel R becomes an interpretable G ---

def _blocking_coating_data():
    """L-free R0 - (R_huge | Q): R is far too large for the window to resolve."""
    rng = np.random.default_rng(0)
    circuit = R(12.0) - (R(5e9) | Q(3e-7, 0.88))
    Z = circuit.impedance(FREQ, circuit.get_all_params())
    noise = 0.01 * np.abs(Z)
    return Z + noise * (rng.standard_normal(len(FREQ))
                        + 1j * rng.standard_normal(len(FREQ)))


def test_unresolvable_parallel_resistance_reports_a_usable_interval():
    Z = _blocking_coating_data()
    circuit = R(10.0) - (G(1e-6) | Q(1e-7, 0.9))
    result, _, _ = fit_equivalent_circuit(FREQ, Z, circuit, weighting='modulus', plot=False)

    idx = next(i for i, lab in enumerate(result.param_labels)
               if lab.startswith('G'))
    g, sigma = result.params_opt[idx], result.params_stderr[idx]

    assert np.isfinite(g) and np.isfinite(sigma)
    assert result.bound_status[idx] == '', "G near zero must not read as at-bound"
    # The interval has to actually cover the truth (G_true = 2e-10 S) and
    # admit zero - that is what makes it a one-sided lower bound on R.
    assert g - 2 * sigma <= 0.0
    assert g - 2 * sigma <= 2e-10 <= g + 2 * sigma
    # ...and it must be an interval, not the collapsed covariance R produces.
    assert sigma > 0


def test_fit_quality_matches_the_R_parametrization():
    """Reparametrization must not change the fit, only what can be said about it."""
    Z = _blocking_coating_data()
    fit_G, _, _ = fit_equivalent_circuit(
        FREQ, Z, R(10.0) - (G(1e-6) | Q(1e-7, 0.9)), weighting='modulus', plot=False)
    fit_R, _, _ = fit_equivalent_circuit(
        FREQ, Z, R(10.0) - (R(1e6) | Q(1e-7, 0.9)), weighting='modulus', plot=False)
    assert fit_G.fit_error_rel == pytest.approx(fit_R.fit_error_rel, rel=1e-3)


# --- 6. integration ---

def test_parses_from_circuit_expression():
    from eis_analysis.cli.utils import parse_circuit_expression
    circuit = parse_circuit_expression("R(10) - (G(1e-9) | C(1e-6))")
    assert circuit.get_param_labels() == ['R', 'G', 'C']


def test_old_gerischer_syntax_gives_a_migration_error():
    from eis_analysis.cli.utils import parse_circuit_expression
    with pytest.raises(ValueError, match="renamed to GE"):
        parse_circuit_expression("R(10) - G(100, 1e-3)")


def test_oxide_analysis_reads_G_as_the_parallel_resistance():
    from eis_analysis.analysis.oxide import _find_capacitive_elements
    with_G = _find_capacitive_elements(R(10) - (G(1e-8) | C(1e-6)), FREQ)
    with_R = _find_capacitive_elements(R(10) - (R(1e8) | C(1e-6)), FREQ)
    assert with_G[0]['R'] == pytest.approx(with_R[0]['R'], rel=1e-9)


def test_oxide_analysis_reads_a_series_G_as_series_resistance():
    from eis_analysis.analysis.oxide import _find_series_resistance
    assert _find_series_resistance(G(0.1) - (R(1e5) | C(1e-6))) == pytest.approx(10.0)
    # G = 0 is an open series branch and contributes nothing
    assert _find_series_resistance(G(0.0) - (R(1e5) | C(1e-6))) is None


def test_migration_hint_stays_off_unrelated_parse_errors():
    from eis_analysis.cli.utils import parse_circuit_expression
    with pytest.raises(ValueError, match="name 'X' is not defined") as exc:
        parse_circuit_expression("R(10) - (G(1e-9) | X(3))")
    assert "renamed to GE" not in str(exc.value)


def test_oxide_analysis_treats_zero_conductance_as_no_DC_path():
    from eis_analysis.analysis.oxide import _find_capacitive_elements
    elements = _find_capacitive_elements(R(10) - (G(0.0) | C(1e-6)), FREQ)
    assert elements[0]['R'] is None

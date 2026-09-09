#!/usr/bin/env python3
"""
Tests for --fit-on: the Z-HIT reconstruction used as a data correction.

Two things are checked here.

1. The scientific claim. A spectrum whose low-frequency modulus drifts during
   the measurement while the phase stays sound - a coating taking up water,
   the case Zahner's Original/Smoothed/Z-HIT switch is built for - is fitted
   both as measured and against the Z-HIT reconstruction of |Z| from the
   phase. The reconstruction must recover the true resistances; the raw fit
   must not.

2. The plumbing. `apply_zhit_reconstruction` attaches or substitutes the
   reconstruction per --fit-on, the frequency filter carries it under the same
   mask, and `LoadedData.Z_for_fit` hands the right array to the fit.
"""

import argparse

import numpy as np
import pytest

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from eis_analysis.cli.data_handling import filter_by_frequency
from eis_analysis.cli.handlers.validation import apply_zhit_reconstruction
from eis_analysis.cli.utils import EISAnalysisError, LoadedData
from eis_analysis.fitting import fit_equivalent_circuit, R, Q
from eis_analysis.validation import zhit_validation
from eis_analysis.validation.zhit import ZHITResult


# =============================================================================
# Reference spectrum
# =============================================================================
# Rs - (R0||Q0) - (R1||Q1) with both relaxations well inside 1 mHz - 100 kHz
# (tau0 ~ 0.1 ms, tau1 ~ 1 s), so both arcs close within the window and every
# resistance is identifiable. Without that the low-frequency arc is open at the
# edge, R1 runs to its bound on any perturbation, and the comparison below
# would measure the bound rather than the correction.
FREQUENCIES = np.logspace(-3, 5, 81)


def _truth():
    """The circuit the synthetic spectrum is generated from."""
    return R(10) - (R(100) | Q(1e-6, 0.9)) - (R(1000) | Q(1e-3, 0.85))


def _initial_guess():
    """Deliberately off the truth by ~3x, so a passing fit means convergence."""
    return R(30) - (R(300) | Q(3e-6, 0.8)) - (R(3000) | Q(3e-3, 0.8))


# Indices of Rs, R0 and R1 in the flat parameter vector
R_INDICES = (0, 1, 4)


def _clean_spectrum():
    circuit = _truth()
    params = np.array(circuit.get_all_params(), dtype=float)
    return circuit.impedance(FREQUENCIES, params), params


def _apply_lf_drift(Z):
    """
    Shrink |Z| by up to 30%, ramped over the two decades below 0.1 Hz, leaving
    the phase untouched.

    This is the drift Z-HIT is meant to undo: the magnitude is wrong, the phase
    - which the transform integrates - is not.
    """
    ramp = np.clip((np.log10(0.1) - np.log10(FREQUENCIES)) / 2.0, 0.0, 1.0)
    return np.abs(Z) * (1.0 - 0.30 * ramp) * np.exp(1j * np.angle(Z))


def _fit_max_resistance_error(Z, params_true):
    """Fit the reference circuit to Z, return the worst relative R error [%]."""
    result, _, _ = fit_equivalent_circuit(
        FREQUENCIES, Z, _initial_guess(), weighting='modulus'
    )
    plt.close('all')
    params = np.array(result.params_opt, dtype=float)
    return max(abs(params[i] - params_true[i]) / params_true[i]
               for i in R_INDICES) * 100.0


# =============================================================================
# The scientific claim
# =============================================================================

def test_reconstruction_recovers_resistances_from_drifted_modulus():
    """A 30% low-frequency modulus drift is undone by fitting the Z-HIT curve."""
    Z_clean, params_true = _clean_spectrum()
    Z_drifted = _apply_lf_drift(Z_clean)

    reconstruction = zhit_validation(FREQUENCIES, Z_drifted)
    plt.close('all')
    assert reconstruction.success

    err_original = _fit_max_resistance_error(Z_drifted, params_true)
    err_reconstructed = _fit_max_resistance_error(reconstruction.Z_fit,
                                                  params_true)

    # Measured on this spectrum: 19.6% raw, 0.12% reconstructed (factor ~160).
    # The thresholds keep an order of magnitude of headroom on each side.
    assert err_original > 10.0, (
        f"drift did not bias the raw fit ({err_original:.2f}%) - the test "
        "no longer exercises what it claims to")
    assert err_reconstructed < 1.0, (
        f"reconstruction did not recover the resistances ({err_reconstructed:.2f}%)")
    assert err_original / err_reconstructed > 20.0


def test_reconstruction_costs_accuracy_on_noisy_stationary_data():
    """The switch is not free: on noise without drift it makes the fit worse.

    The second-order term differentiates the phase, so phase noise is amplified
    rather than smoothed (open point 2 of doc/ZHIT_AUDIT_2026-04-26.md). This
    pins the trade-off the README warns about, so it cannot quietly change.

    Measured over seeds 0-4 at 1% noise, max resistance error raw vs
    reconstructed: 0.45/2.88, 0.59/7.27, 0.33/1.48, 0.58/2.91, 1.60/1.45 -
    usually several times worse, occasionally a wash. Seed 1 is used here for
    its clear margin. Should the reconstruction ever learn to smooth the phase,
    this test is meant to fail and be rewritten, not silently kept passing.
    """
    Z_clean, params_true = _clean_spectrum()
    rng = np.random.default_rng(1)
    sigma = np.abs(Z_clean) * 0.01
    Z_noisy = Z_clean + rng.normal(0, sigma) + 1j * rng.normal(0, sigma)

    reconstruction = zhit_validation(FREQUENCIES, Z_noisy)
    plt.close('all')

    err_raw = _fit_max_resistance_error(Z_noisy, params_true)
    err_reconstructed = _fit_max_resistance_error(reconstruction.Z_fit,
                                                  params_true)

    assert err_raw < 2.0, f"1% noise alone should barely move the fit ({err_raw:.2f}%)"
    assert err_reconstructed > 3 * err_raw, (
        f"reconstruction no longer costs accuracy on noisy data "
        f"({err_raw:.2f}% raw vs {err_reconstructed:.2f}% reconstructed) - "
        "if that is a deliberate improvement, rewrite this test")


def test_reconstruction_is_least_accurate_at_the_high_frequency_edge():
    """Where --fit-on all's cost lands: np.gradient degrades at the edges.

    R_inf and the high-frequency end of the DRT read only that edge, which is
    the opposite end of the spectrum from the drift the switch corrects. The
    CLI warns about this under --fit-on all.
    """
    Z_clean, _ = _clean_spectrum()
    residuals = abs(zhit_validation(FREQUENCIES, Z_clean).residuals_mag)
    plt.close('all')

    lowest_decade = residuals[FREQUENCIES <= 1e-2].mean()
    highest_decade = residuals[FREQUENCIES >= 1e4].mean()

    # Measured: 0.08% vs 1.02% on this exactly K-K compliant spectrum.
    assert lowest_decade < 0.3
    assert highest_decade > 3 * lowest_decade


def test_reconstruction_error_floor_on_undisturbed_data():
    """On a stationary spectrum the switch must cost almost nothing."""
    Z_clean, params_true = _clean_spectrum()

    reconstruction = zhit_validation(FREQUENCIES, Z_clean)
    plt.close('all')

    # Z-HIT has its own error floor (numerical integration plus the edge
    # behavior of np.gradient in the second-order term). Measured: 0.58% mean
    # magnitude residual, 0.12% resistance error.
    assert reconstruction.mean_residual_mag < 2.0
    assert _fit_max_resistance_error(reconstruction.Z_fit, params_true) < 1.0


# =============================================================================
# Plumbing
# =============================================================================

def _args(fit_on='original', f_min=None, f_max=None):
    return argparse.Namespace(fit_on=fit_on, f_min=f_min, f_max=f_max)


def _loaded():
    Z, _ = _clean_spectrum()
    return LoadedData(frequencies=FREQUENCIES, Z=Z, title="test",
                      metadata=None)


def _reconstruction_of(data):
    result = zhit_validation(data.frequencies, data.Z)
    plt.close('all')
    return result


def test_fit_on_original_leaves_data_untouched():
    data = _loaded()
    out = apply_zhit_reconstruction(data, _reconstruction_of(data),
                                    _args('original'))
    assert out is data
    assert out.Z_zhit is None
    assert out.Z_for_fit is out.Z


def test_fit_on_zhit_attaches_reconstruction_without_replacing_z():
    data = _loaded()
    reconstruction = _reconstruction_of(data)
    out = apply_zhit_reconstruction(data, reconstruction, _args('zhit'))

    np.testing.assert_array_equal(out.Z, data.Z)          # R_inf/DRT untouched
    np.testing.assert_array_equal(out.Z_zhit, reconstruction.Z_fit)
    np.testing.assert_array_equal(out.Z_for_fit, reconstruction.Z_fit)


def test_fit_on_all_replaces_z_and_marks_the_title():
    data = _loaded()
    reconstruction = _reconstruction_of(data)
    out = apply_zhit_reconstruction(data, reconstruction, _args('all'))

    np.testing.assert_array_equal(out.Z, reconstruction.Z_fit)
    assert out.Z_zhit is None
    np.testing.assert_array_equal(out.Z_for_fit, reconstruction.Z_fit)
    assert "Z-HIT" in out.title


def _failed_reconstruction():
    """What zhit_validation returns when the integration blew up."""
    empty = np.array([])
    return ZHITResult(
        Z_mag_reconstructed=empty,
        Z_fit=np.array([], dtype=np.complex128),
        residuals_mag=empty, residuals_real=empty, residuals_imag=empty,
        pseudo_chisqr=0.0, noise_estimate=0.0, quality=0.0, ref_freq=1.0,
    )


@pytest.mark.parametrize('unusable', [None, _failed_reconstruction()])
def test_unusable_reconstruction_raises_instead_of_falling_back(unusable):
    """Falling back to the original would fit something else than asked for."""
    with pytest.raises(EISAnalysisError, match="Z-HIT"):
        apply_zhit_reconstruction(_loaded(), unusable, _args('zhit'))


def test_frequency_filter_masks_the_reconstruction_alongside_z():
    """Both arrays must survive the same mask or the fit pairs wrong points."""
    data = _loaded()
    attached = apply_zhit_reconstruction(data, _reconstruction_of(data),
                                         _args('zhit'))

    filtered = filter_by_frequency(attached, _args('zhit', f_min=1.0, f_max=1e3))

    assert len(filtered.Z_zhit) == len(filtered.frequencies) < len(FREQUENCIES)
    mask = (FREQUENCIES >= 1.0) & (FREQUENCIES <= 1e3)
    np.testing.assert_array_equal(filtered.Z_zhit, attached.Z_zhit[mask])
    np.testing.assert_array_equal(filtered.Z_for_fit, filtered.Z_zhit)


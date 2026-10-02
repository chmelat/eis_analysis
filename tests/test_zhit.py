#!/usr/bin/env python3
"""
Tests for the Z-HIT integration constant (doc/ZHIT_REVIEW.md, 1.1 and 1.2).

The phase fixes ln|Z| only up to a constant. It used to be taken from a single
point at the geometric-mean frequency, which carried that point's noise and
the local error of the second-order term into the whole reconstruction.
"""

import warnings

import numpy as np
import pytest

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from eis_analysis.validation import zhit_validation, zhit_reconstruct_magnitude


FREQUENCIES = np.logspace(-2, 5, 71)
OMEGA = 2 * np.pi * FREQUENCIES


def _r_rc(tau):
    """R0 = 10 Ohm in series with R = 100 Ohm || C, exactly K-K compliant."""
    return 10 + 100 / (1 + 1j * OMEGA * tau)


@pytest.mark.parametrize("tau", [1e-2, 1.0])
def test_clean_rc_is_not_biased_by_the_anchor(tau):
    """tau = 1e-2 puts the relaxation at the geometric-mean frequency, where
    the old single-point anchor gave 1.72%. Measured now: 0.63% and 0.60%."""
    result = zhit_validation(FREQUENCIES, _r_rc(tau))
    plt.close('all')
    assert result.mean_residual_mag < 1.0


def test_outlier_at_the_middle_does_not_shift_the_rest():
    """A 5% outlier at the old anchor point moved every other point by ~3%."""
    Z = _r_rc(1.0)
    mid = int(np.argmin(np.abs(FREQUENCIES - np.sqrt(FREQUENCIES[0] * FREQUENCIES[-1]))))
    Z[mid] *= 1.05
    residuals = np.abs(zhit_validation(FREQUENCIES, Z).residuals_mag)
    plt.close('all')
    # Measured: 0.29% median over the other points (3.3% with the old anchor)
    assert np.median(np.delete(residuals, mid)) < 0.5


def test_offset_is_the_median_difference():
    Z = _r_rc(1e-2)
    ln_Z_exp = np.log(np.abs(Z))
    ln_Z = zhit_reconstruct_magnitude(FREQUENCIES, np.angle(Z), ln_Z_exp)
    assert np.median(ln_Z_exp - ln_Z) == pytest.approx(0.0, abs=1e-12)


def test_validation_leaves_no_open_figures():
    plt.close('all')
    for _ in range(3):
        zhit_validation(FREQUENCIES, _r_rc(1e-2))
    assert plt.get_fignums() == []


# A sub-band measured twice, as in an up/down sweep that starts mid-range
# (Zahner EisGenerateJob, Thales) or a file holding two sweeps.
REMEASURED = FREQUENCIES[40:50]


def _noisy(frequencies, seed=0):
    rng = np.random.default_rng(seed)
    Z = 10 + 100 / (1 + 1j * 2 * np.pi * frequencies * 1e-2)
    return Z * (1 + 1e-3 * (rng.standard_normal(Z.size) + 1j * rng.standard_normal(Z.size)))


@pytest.mark.parametrize("shift", [0.0, 1e-4], ids=["exact", "near"])
def test_repeated_frequencies_leave_the_reconstruction_unchanged(shift):
    """Used to give NaN everywhere: np.gradient divided by a zero step."""
    frequencies = np.concatenate([FREQUENCIES, REMEASURED * (1 + shift)])
    Z = 10 + 100 / (1 + 1j * 2 * np.pi * frequencies * 1e-2)
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        residuals = zhit_validation(frequencies, Z).residuals_mag
    single = zhit_validation(FREQUENCIES, _r_rc(1e-2)).residuals_mag
    plt.close('all')
    assert len(residuals) == len(frequencies)
    # Point by point, not the mean: the re-measured band counts twice there.
    # Measured: 0.008 percentage points, the shift of the median offset
    np.testing.assert_allclose(residuals[:len(FREQUENCIES)], single, atol=0.02)
    np.testing.assert_allclose(residuals[len(FREQUENCIES):], residuals[40:50], atol=0.02)


@pytest.mark.parametrize("shift", [1e-4, 0.009, 0.011, 0.049])
def test_close_frequencies_do_not_spike_on_noise(shift):
    """A 1e-4 step turns 1e-3 rad of phase noise into a derivative of 10 rad.
    Plain np.gradient gave 1.4e6 %, 10.6 %, 9.2 % and 4.2 % max residual; the
    steps around 1 % guard against a fixed merge threshold, which is a cliff."""
    frequencies = np.concatenate([FREQUENCIES, REMEASURED * (1 + shift)])
    Z = _noisy(frequencies)
    residuals = np.abs(zhit_validation(frequencies, Z).residuals_mag)
    single = np.abs(zhit_validation(FREQUENCIES, Z[:len(FREQUENCIES)]).residuals_mag)
    plt.close('all')
    # Measured: max 2.89-2.90 % against 2.85 % for the same points without
    # the re-measured band
    assert residuals.max() < single.max() + 0.1


@pytest.mark.parametrize("frequencies", [
    np.logspace(-2, 5, 2100),
    np.unique(np.concatenate([FREQUENCIES, np.logspace(0, 1, 300)])),
], ids=["dense", "dense-band"])
def test_dense_sweeps_are_reconstructed(frequencies):
    """Every step here is below MIN_DERIVATIVE_STEP; the derivative must still
    come from real neighbours (merging close points chained such a grid into
    one group: an empty result, or 138 % residuals on the dense band)."""
    Z = 10 + 100 / (1 + 1j * 2 * np.pi * frequencies * 1e-2)
    residuals = np.abs(zhit_validation(frequencies, Z).residuals_mag)
    plt.close('all')
    assert len(residuals) == len(frequencies)
    # Measured: 3.23 % and 3.27 %, against 2.84 % on the 10 points/decade grid
    assert residuals.max() < 3.5


def test_dense_noisy_sweep_is_not_dominated_by_noise():
    """The minimum step averages over neighbours that np.gradient would use
    one by one; at 300 points/decade and 0.1 % noise it gave 17.2 % max."""
    frequencies = np.logspace(-2, 5, 2100)
    residuals = np.abs(zhit_validation(frequencies, _noisy(frequencies)).residuals_mag)
    plt.close('all')
    # Measured: 5.2 %
    assert residuals.max() < 6.0


@pytest.mark.parametrize("frequencies", [np.array([]), np.array([1.0]), np.array([1.0, 1.01])])
def test_too_narrow_a_spectrum_is_a_clear_error(frequencies):
    with pytest.raises(ValueError, match="Z-HIT needs frequencies spanning"):
        zhit_reconstruct_magnitude(frequencies, np.zeros_like(frequencies), np.zeros_like(frequencies))

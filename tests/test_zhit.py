#!/usr/bin/env python3
"""
Tests for the Z-HIT integration constant (doc/ZHIT_REVIEW.md, 1.1 and 1.2).

The phase fixes ln|Z| only up to a constant. It used to be taken from a single
point at the geometric-mean frequency, which carried that point's noise and
the local error of the second-order term into the whole reconstruction.
"""

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


def test_plot_zhit_validation_marks_flagged_frequencies():
    from eis_analysis.visualization import plot_zhit_validation

    Z = _r_rc(1e-2)
    result = zhit_validation(FREQUENCIES, Z)
    plt.close('all')
    fig = plot_zhit_validation(FREQUENCIES, Z, result, flagged_frequencies=[FREQUENCIES[5]])
    try:
        assert len(fig.axes) == 2
        assert len(fig.axes[1].lines) == 2 + 3 + 1
    finally:
        plt.close(fig)


def test_validation_leaves_no_open_figures():
    plt.close('all')
    for _ in range(3):
        zhit_validation(FREQUENCIES, _r_rc(1e-2))
    assert plt.get_fignums() == []

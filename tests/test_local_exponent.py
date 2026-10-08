"""Local CPE exponent map: flat for one CPE, flagged for two processes."""

import numpy as np

from eis_analysis.analysis import local_exponent
from eis_analysis.fitting import Wa


def test_flat_for_one_cpe_and_flagged_for_two_processes():
    f = np.logspace(-3, 5, 81)
    w = 2 * np.pi * f
    R_inf, L = 1.2, 8e-8
    series = R_inf + 1j * w * L

    # One CPE parallel to C: n(f) is the CPE exponent everywhere, no warning
    Z_cpe = series + 1 / (2e-7 * (1j * w) ** 0.6 + 1j * w * 7e-8)
    res = local_exponent(f, Z_cpe, R_inf, L)
    assert res.valid.sum() > 50
    assert np.allclose(res.n[res.valid], 0.6, atol=1e-3)
    assert res.warnings == []

    # A blocked anomalous channel beside the CPE (the M136 ZrO2 model):
    # n(f) dips between the two exponents and the span is reported
    wa = Wa(2e7, 20, 0.68)
    Y = 1 / wa.impedance(f, [2e7, 20, 0.68]) + 1.5e-7 * (1j * w) ** 0.74 + 1j * w * 7e-8
    res = local_exponent(f, series + 1 / Y, R_inf, L)
    assert res.span > 0.1
    assert len(res.warnings) == 1 and 'span' in res.warnings[0]


def test_rinf_bound_leaves_no_point_wrongly_determined():
    """estimate_rinf falls back to Re Z(f_max), an upper bound, when R_inf is
    not determined; on an oxide that can be 100x R_s. +-5 % of it marked
    points determined to 0.02 that were 0.2 off (stress test, oxide/119;
    here 0.15). With R_inf_range = (0, bound) every determined point holds
    n (measured: 0.019)."""
    f = np.logspace(-2, 6, 81)
    w = 2 * np.pi * f
    Z = 3.0 + 1 / (1e-9 * (1j * w) ** 0.85 + 1j * w * 2e-11)
    bound = float(Z[-1].real)
    assert bound > 50 * 3.0     # 271 Ohm, 90x R_s

    res = local_exponent(f, Z, bound, R_inf_range=(0.0, bound))
    assert res.valid.sum() > 30
    assert np.max(np.abs(res.n[res.valid] - 0.85)) <= 2 * res.uncertainty_max

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

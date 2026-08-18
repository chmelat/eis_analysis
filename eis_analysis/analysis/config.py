"""
Configuration constants for oxide layer analysis.
"""

# Permittivity of vacuum [F/cm]
# Standard: ε₀ = 8.854×10⁻¹² F/m = 8.854×10⁻¹⁴ F/cm
EPSILON_0 = 8.854e-14

# Relative permittivity assumed when the caller does not supply one.
# 22 is the commonly quoted value for ZrO₂, the oxide this toolkit is
# primarily used on. Other oxides differ widely (Al₂O₃ ~ 9, TiO₂ ~ 80,
# SiO₂ ~ 3.9), so the default is only a starting point.
DEFAULT_EPSILON_R = 22.0

# Electrode area assumed when neither the caller nor the DTA metadata
# supplies one. Unit area means the specific capacitance equals the
# measured capacitance, so results are reported per whatever area the
# measurement actually had - correct only if that area really was 1 cm².
DEFAULT_AREA_CM2 = 1.0

# CPE exponent below which the effective-capacitance conversion
# (Hsu-Mansfeld) is not well-defined: the distribution of time constants
# is too broad for a single C_eff to be meaningful.
CPE_N_RELIABLE_MIN = 0.8

# Width of the top frequency band (one decade, f >= f_max/10) used for the
# median-based high-frequency capacitance estimate in fallback mode.
HF_ESTIMATE_DECADE_FACTOR = 10.0

# Maximum spread (max/min) of the per-point estimates C_i = -1/(ω·Z'')
# across the top frequency decade. C_i is frequency-independent when the
# capacitance dominates (ωRC ≫ 1 for a parallel R||C); a spread above this
# ratio means that assumption does not hold within the decade and the
# fallback estimate is unreliable. 1.2 corresponds to ωRC >= ~2 at the
# bottom of the decade (C_i = C·(1 + 1/(ωRC)²) for R||C).
HF_C_SPREAD_MAX_RATIO = 1.2

# Series resistance below which the Brug (2D) effective capacitance is not
# reported. The default optimizer floor for an R parameter is 0.1 mOhm, and a
# CPE with n < 1 mimics a series resistance at high frequency, so a degenerate
# fit drives R_s to that floor instead of to the true ohmic resistance. Brug
# gives C_eff ~ R_s^((1-n)/n), so a floored R_s yields an arbitrarily small
# capacitance (and permittivity) with no warning. 10 mOhm sits far below the
# ohmic resistance of any real electrochemical cell (electrolyte + leads), so
# a value under it means "not identified by the fit", not "highly conductive".
BRUG_RS_MIN_OHM = 1e-2

# Ratio C_eff(Hsu-Mansfeld) / C_eff(Brug) above which the two CPE models are
# flagged as diverging rather than bracketing C_eff. The ratio is exactly
# (1 + R_ct/R_s)^((1-n)/n), so it grows both with a blocking layer (large
# R_ct) and with a poorly determined R_s. Within one order of magnitude the
# pair still brackets a useful C_eff; beyond that the spread is dominated by
# how well R_s is known and the "comparison" carries no information.
BRUG_HM_DIVERGENCE_MAX = 10.0

__all__ = [
    'EPSILON_0',
    'DEFAULT_EPSILON_R',
    'DEFAULT_AREA_CM2',
    'CPE_N_RELIABLE_MIN',
    'HF_ESTIMATE_DECADE_FACTOR',
    'HF_C_SPREAD_MAX_RATIO',
    'BRUG_RS_MIN_OHM',
    'BRUG_HM_DIVERGENCE_MAX',
]

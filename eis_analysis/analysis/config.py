"""
Configuration constants for oxide layer analysis.
"""

# Permittivity of vacuum [F/cm]
# Standard: ε₀ = 8.854×10⁻¹² F/m = 8.854×10⁻¹⁴ F/cm
EPSILON_0 = 8.854e-14

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

__all__ = [
    'EPSILON_0',
    'CPE_N_RELIABLE_MIN',
    'HF_ESTIMATE_DECADE_FACTOR',
    'HF_C_SPREAD_MAX_RATIO',
]

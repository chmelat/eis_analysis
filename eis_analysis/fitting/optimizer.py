"""
least_squares with unit-free stopping criteria, shared by the LM fit and the
DE polish.
"""

import numpy as np
from scipy.optimize import least_squares

# Scale of a parameter whose start is exactly 0. A zero carries no magnitude,
# and the bounds give none either when the lower one is 0 (G: 0..1e4 S would
# make the first step 1e4 S). 1e-10 keeps the first steps small; it is the
# floor of the former x_scale, so zero starts behave as before. It only sizes
# the steps: the parameter itself may still take any value within its bounds.
ZERO_START_SCALE = 1e-10


def least_squares_normalized(residual, x0, jac, bounds, **kwargs):
    """least_squares over u = x / |x0|, returned in the original variables.

    least_squares stops on criteria that mix units. xtol compares norm(step)
    with norm(x) over R ~ 1e7 Ohm and Q ~ 1e-10 alike; its x_scale only
    steers the trust region, not that test. gtol is an absolute threshold on
    the gradient, so it depends on the cost's units, and on a dimensionless
    cost it ends the fit early wherever a parameter is small relative to |Z|
    (a 0.05 Ohm R_s in front of 3e3 Ohm stopped at 0.064 Ohm). A fit of the
    same spectrum in other units therefore stopped elsewhere (stress test:
    R_inf moved under Z -> Z/2, an exact binary scaling).

    Here every normalized parameter starts at +-1, which makes xtol unit-free,
    and gtol is disabled, so the fit ends on the relative ftol or xtol. The
    scale is |x0| itself, without a floor: the former floor of 1e-10 rescaled
    a 1e-12 F capacitor differently from a 1e-9 F one. Only an exact zero
    start gets ZERO_START_SCALE.

    Returns the least_squares result with `x` and `jac` in the original
    variables; `fun` and `cost` do not depend on them. `grad`, `optimality`
    and `active_mask` stay in the normalized variables u.
    """
    x0 = np.asarray(x0, dtype=float)
    scale = np.abs(x0)
    scale[scale == 0] = ZERO_START_SCALE
    jac_u = jac if isinstance(jac, str) else (lambda u: jac(u * scale) * scale)
    result = least_squares(
        lambda u: residual(u * scale), x0 / scale, jac=jac_u,
        bounds=(np.asarray(bounds[0], dtype=float) / scale,
                np.asarray(bounds[1], dtype=float) / scale),
        gtol=None, **kwargs)
    result.x = result.x * scale
    result.jac = result.jac / scale
    return result

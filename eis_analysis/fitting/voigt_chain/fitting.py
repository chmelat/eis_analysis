"""
Main Voigt chain fitting functions (Lin-KK compatible).

This module provides the core fitting functions for Voigt chain models
using linear regression.
"""

import numpy as np
import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Tuple, List, Optional
from numpy.typing import NDArray

if TYPE_CHECKING:  # circular at runtime: mu_optimization imports from here
    from .mu_optimization import MuOptimization

from .validation import validate_eis_data, validate_tau
from .tau_grid import generate_tau_grid
from .solvers import robust_nnls

from ..circuit_elements import R, K, L
from ..diagnostics import compute_weights
from ..circuit_builder import Series, Circuit

logger = logging.getLogger(__name__)

# Type alias for Voigt chain circuit
VoigtChain = Series


@dataclass
class VoigtChainDiagnostics:
    """
    How fit_voigt_chain_linear() arrived at its circuit.

    Mirrors the four steps the CLI reports: how the tau grid was chosen,
    what the regression produced, what pruning removed, and what the
    circuit ended up as. A caller that only wants the circuit can ignore
    all of it; the CLI cannot compose its section without it.
    """
    # Step 1 - tau grid: either the mu search or a fixed grid
    auto_optimize_M: bool
    mu_optimization: Optional["MuOptimization"] = None
    mu_threshold: float = 0.85
    max_M: int = 50
    n_per_decade: int = 3
    extend_decades: float = 0.0
    n_tau_generated: int = 0

    # Step 2 - linear regression
    fit_type: str = 'complex'
    weighting: str = 'modulus'
    allow_negative: bool = False
    include_Rs: bool = True
    include_L: bool = True
    R_s: float = 0.0
    R_i_min: Optional[float] = None   # from the regression, before pruning
    R_i_max: Optional[float] = None
    residual: float = 0.0
    L_value: Optional[float] = None

    # Step 3 - pruning
    pruning_enabled: bool = False
    prune_threshold: float = 0.0
    threshold_relative: float = 0.0
    threshold_absolute: float = 0.0
    threshold_effective: float = 0.0
    n_before_prune: int = 0
    n_after_prune: int = 0

    # Step 4 - the circuit that survived pruning
    n_elements: int = 0
    R_i_min_kept: Optional[float] = None
    R_i_max_kept: Optional[float] = None
    tau_min: Optional[float] = None
    tau_max: Optional[float] = None
    n_params: int = 0

    warnings: List[str] = field(default_factory=list)


@dataclass
class VoigtChainFit:
    """Result of fit_voigt_chain_linear(): the circuit and how it was built."""
    circuit: Circuit
    initial_params: List[float]
    diagnostics: VoigtChainDiagnostics


def estimate_R_linear(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    tau: NDArray[np.float64],
    include_Rs: bool = True,
    include_L: bool = True,
    include_C: bool = False,
    fit_type: str = 'complex',
    allow_negative: bool = False,
    weighting: str = 'modulus'
) -> Tuple[NDArray[np.float64], float, Optional[float], Optional[float]]:
    """
    Estimate R_i values using least squares (Lin-KK compatible).

    Implements the linear Kramers-Kronig test from Schonleber et al. (2014).
    The model includes series resistance R_s, M Voigt elements (R_k, tau_k),
    and optionally series inductance L and series capacitance C.

    Model:
        Z(omega) = R_s + sum R_k/(1 + j*omega*tau_k) + j*omega*L + 1/(j*omega*C)

    Parameters
    ----------
    frequencies : ndarray of float
        Measured frequencies [Hz] (M points)
    Z : ndarray of complex
        Measured impedance [Ohm] (M points)
    tau : ndarray of float
        Fixed time constants [s] (N points)
    include_Rs : bool, optional
        If True, estimate R_s (series resistance) (default: True)
    include_L : bool, optional
        If True, include series inductance L in the model (default: True)
        This captures measurement system artifacts (Lin-KK standard).
    include_C : bool, optional
        If True, include series capacitance C in the model (default: False).
        Captures blocking (capacitive) low-frequency behavior, e.g. in
        two-electrode cells (Schonleber Lin-KK 'add_cap'). Like L, a series
        C has zero real part, so in 'real' mode it is extracted from the
        imaginary residual.
    fit_type : str, optional
        Which components to fit:
        - 'real': Fit only real part Z', then extract L (and C) from
          imaginary residual
        - 'imag': Fit only imaginary part Z''
        - 'complex': Fit both parts simultaneously (default)
    allow_negative : bool, optional
        If False (default), use NNLS to enforce R_i >= 0 (physical constraint)
        If True, use pseudoinverse (allows negative R_i like Lin-KK test)
    weighting : str, optional
        Point weighting scheme (default: 'modulus' = Lin-KK standard):
        - 'uniform': all points equal weight (w = 1)
        - 'sqrt': compromise weighting (w = 1/sqrt|Z|)
        - 'modulus': Lin-KK standard (w = 1/|Z|) - DEFAULT
        - 'proportional': strong low-Z emphasis (w = 1/|Z|^2)

    Returns
    -------
    R : ndarray of float
        Estimated resistances [Ohm]: [R_s, R_1, R_2, ..., R_N] or [R_1, ..., R_N]
        May contain negative values if allow_negative=True
    residual : float
        Residual norm
    L_value : float or None
        Estimated inductance [H] if include_L=True, else None
    C_value : float or None
        Estimated series capacitance [F] if include_C=True, else None.
        The linear fit coefficient is D = 1/C; C_value = 1/D (None if D=0).
        Never stored inside the elements array (unlike L).

    Notes
    -----
    The default 'modulus' weighting (1/|Z|) is the Lin-KK standard.
    This ensures equal relative weighting across the frequency range,
    which is critical for wide dynamic range EIS data.

    References
    ----------
    Schonleber, M. et al. "A Method for Improving the Robustness of linear
    Kramers-Kronig Validity Tests." Electrochimica Acta 131, 20-27 (2014)
    """
    # Validate inputs
    validate_eis_data(frequencies, Z, context="estimate_R_linear")
    validate_tau(tau, context="estimate_R_linear")

    # Validate weighting parameter
    valid_weightings = ['uniform', 'sqrt', 'proportional', 'modulus']
    if weighting not in valid_weightings:
        raise ValueError(
            f"estimate_R_linear: weighting must be one of {valid_weightings}, "
            f"got '{weighting}'"
        )

    # Extract components
    Z_real = Z.real
    Z_imag = Z.imag
    omega = 2 * np.pi * frequencies
    n_freq = len(frequencies)
    n_tau = len(tau)

    # Weights (normalized to mean = 1 for numerical stability). Shared with
    # circuit fitting so both paths weight a spectrum identically.
    weights = compute_weights(Z, weighting)

    # Determine matrix dimensions
    # Columns: [R_s (optional), R_1, R_2, ..., R_N, L (optional), D=1/C (optional)]
    n_cols = n_tau
    col_offset = 0

    if include_Rs:
        n_cols += 1
        col_offset = 1

    # Number of resistive columns (R_s + R_i) - everything before L/C tail
    n_R_cols = n_cols

    L_col = None
    if include_L:
        L_col = n_cols
        n_cols += 1

    C_col = None
    if include_C:
        C_col = n_cols
        n_cols += 1

    # Build design matrices with weighting applied
    A_real = np.zeros((n_freq, n_cols))
    A_imag = np.zeros((n_freq, n_cols))

    # R_s column: contributes only to real part
    if include_Rs:
        A_real[:, 0] = 1.0 * weights
        A_imag[:, 0] = 0.0

    # Voigt element columns: K(R, tau) = R / (1 + j*omega*tau)
    for i, tau_i in enumerate(tau):
        tau_omega = tau_i * omega
        tau_omega_sq = tau_omega ** 2
        denom = 1 + tau_omega_sq

        # Real part: R / (1 + (omega*tau)^2)
        A_real[:, col_offset + i] = (1.0 / denom) * weights
        # Imaginary part: -R*omega*tau / (1 + (omega*tau)^2)
        A_imag[:, col_offset + i] = (-tau_omega / denom) * weights

    # Inductance column: Z_L = j*omega*L -> real=0, imag=omega*L
    if include_L:
        A_real[:, L_col] = 0.0
        A_imag[:, L_col] = omega * weights

    # Series capacitance column: Z_C = -j*D/omega with D = 1/C
    # -> real=0, imag=-D/omega (linear in D)
    if include_C:
        A_real[:, C_col] = 0.0
        A_imag[:, C_col] = (-1.0 / omega) * weights

    # Weighted target vectors
    b_real = Z_real * weights
    b_imag = Z_imag * weights

    # Solve based on fit_type
    if fit_type == 'real':
        # Fit only real part (Lin-KK default mode)
        # L and C don't contribute to real part, so exclude them from real fit
        A_real_fit = A_real[:, :n_R_cols]

        if allow_negative:
            elements = np.linalg.pinv(A_real_fit) @ b_real
        else:
            elements, _ = robust_nnls(A_real_fit, b_real)

        # Extract L and C from imaginary residual (Schonleber et al. approach)
        L_value = None
        C_value = None
        if include_L or include_C:
            # Compute Z_fit imaginary part from R_k only (unweighted)
            A_imag_unweighted = np.zeros((n_freq, n_R_cols))
            for i, tau_i in enumerate(tau):
                tau_omega = tau_i * omega
                tau_omega_sq = tau_omega ** 2
                denom = 1 + tau_omega_sq
                A_imag_unweighted[:, col_offset + i] = -tau_omega / denom
            if include_Rs:
                A_imag_unweighted[:, 0] = 0.0

            Z_fit_imag = A_imag_unweighted @ elements

            # Fit tail terms from residual:
            # Z_imag - Z_fit_imag = omega*L - D/omega  (D = 1/C)
            tail_cols = []
            if include_L:
                tail_cols.append(omega * weights)
            if include_C:
                tail_cols.append((-1.0 / omega) * weights)
            A_tail = np.column_stack(tail_cols)
            b_tail = (Z_imag - Z_fit_imag) * weights

            tail = np.linalg.pinv(A_tail) @ b_tail

            if include_L:
                L_value = float(tail[0])
                # Add L to elements (C is kept out of the elements array)
                elements = np.append(elements, L_value)
            if include_C:
                D = float(tail[-1])
                C_value = 1.0 / D if D != 0.0 else None

        residual = np.linalg.norm(A_real_fit @ elements[:n_R_cols] - b_real)

    elif fit_type == 'imag':
        # Fit only imaginary part
        if allow_negative:
            elements = np.linalg.pinv(A_imag) @ b_imag
        else:
            elements, _ = robust_nnls(A_imag, b_imag)

        # For imag fit, find R_s from real part residual (Boukamp approach)
        if include_Rs:
            # Unweighted real-part prediction of the fitted R_k: dividing by
            # `weights` undoes the weighting baked into A_real columns (the
            # R_s column is excluded; the L column of A_real is zero, so L
            # drops out). The previous `* Z_mag_safe` un-weighting was only
            # valid for unnormalized modulus weights (audit K1 2026-07-03).
            z_re_fit = (A_real[:, 1:] / weights[:, None]) @ elements[1:]
            ws = 1 / (Z_real**2 + Z_imag**2 + 1e-30)
            elements[0] = np.sum(ws * (Z_real - z_re_fit)) / np.sum(ws)

        residual = np.linalg.norm(A_imag @ elements - b_imag)
        L_value = float(elements[L_col]) if include_L else None

        # Pop D = 1/C off the solution vector (C is kept out of elements)
        C_value = None
        if include_C:
            D = float(elements[C_col])
            C_value = 1.0 / D if D != 0.0 else None
            elements = elements[:C_col]

    elif fit_type == 'complex':
        # Fit both parts simultaneously (analytical normal equations)
        if allow_negative:
            try:
                ATA = A_real.T @ A_real + A_imag.T @ A_imag
                ATb = A_real.T @ b_real + A_imag.T @ b_imag
                elements = np.linalg.solve(ATA, ATb)
            except np.linalg.LinAlgError:
                # Fallback to pseudoinverse if singular
                A_combined = np.vstack([A_real, A_imag])
                b_combined = np.hstack([b_real, b_imag])
                elements = np.linalg.pinv(A_combined) @ b_combined
        else:
            # NNLS requires stacked form
            A_combined = np.vstack([A_real, A_imag])
            b_combined = np.hstack([b_real, b_imag])
            elements, _ = robust_nnls(A_combined, b_combined)

        # Compute residual
        res_real = A_real @ elements - b_real
        res_imag = A_imag @ elements - b_imag
        residual = np.sqrt(np.sum(res_real**2) + np.sum(res_imag**2))
        L_value = float(elements[L_col]) if include_L else None

        # Pop D = 1/C off the solution vector (C is kept out of elements)
        C_value = None
        if include_C:
            D = float(elements[C_col])
            C_value = 1.0 / D if D != 0.0 else None
            elements = elements[:C_col]

    else:
        raise ValueError(f"Unknown fit_type: {fit_type}. Use 'real', 'imag', or 'complex'.")

    # Log results
    method_str = "pinv" if allow_negative else "NNLS"
    logger.debug(f"Linear regression ({fit_type}, {method_str}):")
    logger.debug(f"  Parameters: {len(elements)}, Residual: {residual:.3e}")
    if include_Rs:
        logger.debug(f"  R_s: {elements[0]:.3e} Ohm")
    if include_L and L_value is not None:
        logger.debug(f"  L: {L_value:.3e} H")
    if include_C and C_value is not None:
        logger.debug(f"  C (series): {C_value:.3e} F")

    # Check for negative R_i (overfit indicator)
    R_start = 1 if include_Rs else 0
    R_end = L_col if include_L else len(elements)
    R_i = elements[R_start:R_end]
    n_negative = np.sum(R_i < 0)
    if n_negative > 0:
        logger.debug(f"  Negative R_i: {n_negative}/{len(R_i)}")

    return elements, residual, L_value, C_value


def fit_voigt_chain_linear(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    n_per_decade: int = 3,
    extend_decades: float = 0.0,
    include_Rs: bool = True,
    include_L: bool = True,
    fit_type: str = 'complex',
    prune_threshold: float = 0.01,
    allow_negative: bool = False,
    auto_optimize_M: bool = False,
    mu_threshold: float = 0.85,
    max_M: int = 50,
    weighting: str = 'modulus'
) -> VoigtChainFit:
    """
    Fit Voigt chain to EIS data using linear regression (Lin-KK method).

    Implements the Lin-KK approach from Schonleber et al. (2014). This is a
    complete fitting method that returns a ready-to-use circuit with optimized
    parameters.

    Model:
        Z(omega) = R_s + sum R_k/(1 + j*omega*tau_k) + j*omega*L

    Parameters
    ----------
    frequencies : ndarray of float
        Measured frequencies [Hz] (N points)
    Z : ndarray of complex
        Measured impedance [Ohm] (N points)
    n_per_decade : int, optional
        Number of tau values per decade (default: 3)
    extend_decades : float, optional
        Extend tau range toward lower frequencies (default: 0.0 for Lin-KK)
    include_Rs : bool, optional
        Include series resistance R_s (default: True)
    include_L : bool, optional
        Include series inductance L in the model (default: True)
    fit_type : str, optional
        Which components to fit: 'real', 'imag', or 'complex' (default)
    prune_threshold : float, optional
        Relative threshold for pruning small R_i elements (default: 0.01 = 1%)
    allow_negative : bool, optional
        Allow negative R_i values (default: False)
    auto_optimize_M : bool, optional
        Automatically find optimal number of elements using mu metric (default: False)
    mu_threshold : float, optional
        Threshold for mu metric when auto_optimize_M=True (default: 0.85)
    max_M : int, optional
        Maximum number of elements to try when auto_optimize_M=True (default: 50)
    weighting : str, optional
        Point weighting scheme (default: 'modulus' = Lin-KK standard)

    Returns
    -------
    VoigtChainFit
        The circuit, its initial parameters, and the diagnostics behind the
        four steps. Nothing is logged; the CLI composes its section from
        the diagnostics.

    References
    ----------
    Schonleber, M. et al. "A Method for Improving the Robustness of linear
    Kramers-Kronig Validity Tests." Electrochimica Acta 131, 20-27 (2014)
    """
    # Validate inputs
    validate_eis_data(frequencies, Z, context="fit_voigt_chain_linear")

    diag = VoigtChainDiagnostics(
        auto_optimize_M=auto_optimize_M,
        mu_threshold=mu_threshold, max_M=max_M,
        n_per_decade=n_per_decade, extend_decades=extend_decades,
        fit_type=fit_type, weighting=weighting,
        allow_negative=allow_negative,
        include_Rs=include_Rs, include_L=include_L,
        prune_threshold=prune_threshold)

    L_value = None

    # Step 1: Generate tau grid OR optimize M using mu metric
    if auto_optimize_M:
        # Import here to avoid circular dependency
        from .mu_optimization import find_optimal_M_mu

        # Use mu optimization to find optimal M (Lin-KK style)
        mu_opt = find_optimal_M_mu(
            frequencies, Z,
            mu_threshold=mu_threshold,
            max_M=max_M,
            extend_decades=extend_decades,
            include_Rs=include_Rs,
            include_L=include_L,
            fit_type=fit_type,
            allow_negative=True,  # mu metric requires negative R detection
            weighting=weighting
        )
        tau = mu_opt.tau
        diag.mu_optimization = mu_opt

        # Step 2: Refit with NNLS to get physically meaningful R values
        elements, residual, L_value, _ = estimate_R_linear(
            frequencies, Z, tau,
            include_Rs=include_Rs,
            include_L=include_L,
            fit_type=fit_type,
            allow_negative=False,  # Force non-negative R for physical circuit
            weighting=weighting
        )

        # Extract R_s and R_i from elements
        if include_Rs:
            R_s = elements[0]
            R_i_end = -1 if include_L else len(elements)
            R_i = elements[1:R_i_end]
        else:
            R_s = 0.0
            R_i_end = -1 if include_L else len(elements)
            R_i = elements[:R_i_end]
    else:
        # Use fixed n_per_decade
        tau = generate_tau_grid(frequencies, n_per_decade, extend_decades)
        diag.n_tau_generated = len(tau)

        # Warn if too many tau (can cause NNLS convergence issues)
        if len(tau) > 20 and not allow_negative:
            diag.warnings.append(
                f"{len(tau)} tau is a lot - may cause convergence issues; "
                f"try --voigt-n-per-decade 2 or --voigt-extend-decades 0.5")

        # Step 2: Linear regression for R values
        elements, residual, L_value, _ = estimate_R_linear(
            frequencies, Z, tau,
            include_Rs=include_Rs,
            include_L=include_L,
            fit_type=fit_type,
            allow_negative=allow_negative,
            weighting=weighting
        )

        # Split R_s, R_i, and L from elements
        if include_Rs:
            R_s = elements[0]
            R_i_start = 1
        else:
            R_s = 0.0
            R_i_start = 0

        if include_L:
            R_i = elements[R_i_start:-1]
        else:
            R_i = elements[R_i_start:]

    diag.R_s = R_s
    diag.residual = residual
    diag.L_value = L_value
    if len(R_i) > 0:
        diag.R_i_min, diag.R_i_max = float(R_i.min()), float(R_i.max())

    # Step 3: Prune small R_i values
    if prune_threshold > 0 and len(R_i) > 0:
        R_max = np.max(np.abs(R_i))  # Use abs for allow_negative case

        # Relative threshold: fraction of max R_i
        threshold_relative = prune_threshold * R_max

        # Absolute minimum threshold: small fraction of total R_pol
        R_pol_total = np.sum(np.abs(R_i))
        threshold_absolute = 0.001 * R_pol_total  # 0.1% of total polarization resistance

        # Use the SMALLER of the two thresholds (more conservative pruning)
        threshold_effective = min(threshold_relative, threshold_absolute)

        diag.pruning_enabled = True
        diag.threshold_relative = threshold_relative
        diag.threshold_absolute = threshold_absolute
        diag.threshold_effective = threshold_effective

        keep_mask = np.abs(R_i) >= threshold_effective

        diag.n_before_prune = len(R_i)
        diag.n_after_prune = int(np.sum(keep_mask))

        if diag.n_after_prune < diag.n_before_prune:
            R_i = R_i[keep_mask]
            tau = tau[keep_mask]

    # Step 4: Build circuit using K(R, tau) elements

    # Start with R_s if included
    circuit: Optional[Circuit]
    if include_Rs and R_s > 0:
        circuit = R(R_s)
    else:
        circuit = None

    # Add K elements: K(R_i, tau_i)
    for i, (r, t) in enumerate(zip(R_i, tau)):
        k_element = K(r, t)

        if circuit is None:
            circuit = k_element
        else:
            circuit = circuit - k_element

    if circuit is None:
        raise ValueError(
            "Cannot build circuit: no series resistance and no K elements"
        )

    # Add inductance if fitted
    if include_L and L_value is not None:
        circuit = circuit - L(L_value)

    # Extract initial parameters
    initial_params = circuit.get_all_params()

    # The summary describes the pruned chain, so these are re-read here
    diag.n_elements = len(R_i)
    if len(R_i) > 0:
        diag.R_i_min_kept, diag.R_i_max_kept = float(R_i.min()), float(R_i.max())
        diag.tau_min, diag.tau_max = float(tau.min()), float(tau.max())
    diag.n_params = len(initial_params)

    return VoigtChainFit(circuit, initial_params, diag)


__all__ = ['estimate_R_linear', 'fit_voigt_chain_linear', 'VoigtChain']

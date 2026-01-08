"""
Main Voigt chain fitting functions (Lin-KK compatible).

This module provides the core fitting functions for Voigt chain models
using linear regression.
"""

import numpy as np
import logging
from typing import Tuple, List, Optional
from numpy.typing import NDArray

from .validation import validate_eis_data, validate_tau
from .tau_grid import generate_tau_grid
from .solvers import robust_nnls

from ..circuit_elements import R, K, L
from ..circuit_builder import Series

logger = logging.getLogger(__name__)

# Type alias for Voigt chain circuit
VoigtChain = Series


def estimate_R_linear(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    tau: NDArray[np.float64],
    include_Rs: bool = True,
    include_L: bool = True,
    fit_type: str = 'complex',
    allow_negative: bool = False,
    weighting: str = 'modulus'
) -> Tuple[NDArray[np.float64], float, Optional[float]]:
    """
    Estimate R_i values using least squares (Lin-KK compatible).

    Implements the linear Kramers-Kronig test from Schonleber et al. (2014).
    The model includes series resistance R_s, M Voigt elements (R_k, tau_k),
    and optionally series inductance L.

    Model:
        Z(omega) = R_s + sum R_k/(1 + j*omega*tau_k) + j*omega*L

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
    fit_type : str, optional
        Which components to fit:
        - 'real': Fit only real part Z', then extract L from imaginary residual
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

    # Compute weights based on weighting scheme
    Z_mag = np.abs(Z)
    Z_mag_safe = np.maximum(Z_mag, 1e-15)

    if weighting == 'uniform':
        weights = np.ones_like(Z_mag)
    elif weighting == 'sqrt':
        weights = 1.0 / np.sqrt(Z_mag_safe)
    elif weighting == 'modulus':
        weights = 1.0 / Z_mag_safe  # Lin-KK standard
    elif weighting == 'proportional':
        weights = 1.0 / (Z_mag_safe ** 2)
    else:
        weights = 1.0 / Z_mag_safe  # Fallback to modulus

    # Normalize weights so mean = 1 (for numerical stability)
    weights = weights / np.mean(weights)

    # Determine matrix dimensions
    # Columns: [R_s (optional), R_1, R_2, ..., R_N, L (optional)]
    n_cols = n_tau
    col_offset = 0

    if include_Rs:
        n_cols += 1
        col_offset = 1

    L_col = None
    if include_L:
        L_col = n_cols
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

    # Weighted target vectors
    b_real = Z_real * weights
    b_imag = Z_imag * weights

    # Solve based on fit_type
    if fit_type == 'real':
        # Fit only real part (Lin-KK default mode)
        # L doesn't contribute to real part, so exclude it from real fit
        if include_L:
            A_real_fit = A_real[:, :L_col]
        else:
            A_real_fit = A_real

        if allow_negative:
            elements = np.linalg.pinv(A_real_fit) @ b_real
        else:
            elements, _ = robust_nnls(A_real_fit, b_real)

        # Extract L from imaginary residual (Schonleber et al. approach)
        L_value = None
        if include_L:
            # Reconstruct fit without L
            R_full = np.zeros(n_cols)
            R_full[:L_col] = elements

            # Compute Z_fit imaginary part from R_k only (unweighted)
            A_imag_unweighted = np.zeros((n_freq, L_col))
            for i, tau_i in enumerate(tau):
                tau_omega = tau_i * omega
                tau_omega_sq = tau_omega ** 2
                denom = 1 + tau_omega_sq
                A_imag_unweighted[:, col_offset + i] = -tau_omega / denom
            if include_Rs:
                A_imag_unweighted[:, 0] = 0.0

            Z_fit_imag = A_imag_unweighted @ elements

            # Fit L from residual: Z_imag - Z_fit_imag = omega*L
            a_L = omega * weights
            b_L = (Z_imag - Z_fit_imag) * weights

            L_value = np.linalg.pinv(a_L.reshape(-1, 1)) @ b_L
            L_value = float(L_value)

            # Add L to elements
            elements = np.append(elements, L_value)

        residual = np.linalg.norm(A_real_fit @ elements[:L_col] - b_real)

    elif fit_type == 'imag':
        # Fit only imaginary part
        if allow_negative:
            elements = np.linalg.pinv(A_imag) @ b_imag
        else:
            elements, _ = robust_nnls(A_imag, b_imag)

        # For imag fit, find R_s from real part residual (Boukamp approach)
        if include_Rs:
            z_re_fit = A_real[:, 1:] @ elements[1:]
            ws = 1 / (Z_real**2 + Z_imag**2 + 1e-30)
            elements[0] = np.sum(ws * (Z_real - z_re_fit * Z_mag_safe)) / np.sum(ws)

        residual = np.linalg.norm(A_imag @ elements - b_imag)
        L_value = float(elements[L_col]) if include_L else None

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

    # Check for negative R_i (overfit indicator)
    R_start = 1 if include_Rs else 0
    R_end = L_col if include_L else len(elements)
    R_i = elements[R_start:R_end]
    n_negative = np.sum(R_i < 0)
    if n_negative > 0:
        logger.debug(f"  Negative R_i: {n_negative}/{len(R_i)}")

    return elements, residual, L_value


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
    refit_positive: bool = False,
    weighting: str = 'modulus'
) -> Tuple[VoigtChain, List[float]]:
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
    refit_positive : bool, optional
        Unused parameter for backward compatibility
    weighting : str, optional
        Point weighting scheme (default: 'modulus' = Lin-KK standard)

    Returns
    -------
    circuit : VoigtChain
        Circuit object: R_s - K(R1, tau1) - K(R2, tau2) - ... [- L]
    initial_params : list of float
        Initial parameter values [R_s, R_1, tau_1, R_2, tau_2, ..., [L]]

    References
    ----------
    Schonleber, M. et al. "A Method for Improving the Robustness of linear
    Kramers-Kronig Validity Tests." Electrochimica Acta 131, 20-27 (2014)
    """
    # Validate inputs
    validate_eis_data(frequencies, Z, context="fit_voigt_chain_linear")

    logger.info("="*60)
    logger.info("Voigt chain initial guess estimation (Lin-KK compatible)")
    logger.info("="*60)

    L_value = None

    # Step 1: Generate tau grid OR optimize M using mu metric
    if auto_optimize_M:
        # Import here to avoid circular dependency
        from .mu_optimization import find_optimal_M_mu

        # Use mu optimization to find optimal M (Lin-KK style)
        logger.info("Step 1: Auto-optimizing M using mu metric (Lin-KK)")
        M_opt, mu_final, tau, elements_mu, L_value_mu = find_optimal_M_mu(
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
        logger.info(f"  Optimal M: {M_opt}, mu: {mu_final:.4f}")

        # Step 2: Refit with NNLS to get physically meaningful R values
        logger.info("Step 2: Refit with NNLS (R_i >= 0) for physical circuit")
        elements, residual, L_value = estimate_R_linear(
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
        logger.info(f"Step 1: Generating tau grid ({n_per_decade} per decade, +{extend_decades} dec extension)")
        tau = generate_tau_grid(frequencies, n_per_decade, extend_decades)
        logger.info(f"  Generated {len(tau)} time constants")

        # Warn if too many tau (can cause NNLS convergence issues)
        if len(tau) > 20 and not allow_negative:
            logger.warning(f"  NOTE: {len(tau)} tau is a lot - may cause convergence issues.")
            logger.warning("  Recommendation: try --voigt-n-per-decade 2 or --voigt-extend-decades 0.5")

        # Step 2: Linear regression for R values
        method_str = "NNLS (R_i >= 0)" if not allow_negative else "pseudoinverse (allows R_i < 0)"
        weighting_labels = {
            'uniform': 'uniform (w=1)',
            'sqrt': 'sqrt (w=1/sqrt|Z|)',
            'modulus': 'modulus (w=1/|Z|, Lin-KK standard)',
            'proportional': 'proportional (w=1/|Z|^2)'
        }
        logger.info(f"Step 2: Linear regression (method: {method_str}, fit_type: {fit_type})")
        logger.info(f"  Weighting: {weighting_labels.get(weighting, weighting)}")
        elements, residual, L_value = estimate_R_linear(
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

    # Log results
    if include_Rs:
        logger.info(f"  R_s (series): {R_s:.3e} Ohm")
    else:
        logger.info("  R_s not included (set to 0)")

    if len(R_i) > 0:
        logger.info(f"  R_i range: [{R_i.min():.3e}, {R_i.max():.3e}] Ohm")
    logger.info(f"  Residual: {residual:.3e}")
    if include_L and L_value is not None:
        logger.info(f"  L (inductance): {L_value:.3e} H")

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

        logger.info("Step 3: Pruning small R_i")
        logger.info(f"  Relative threshold: {prune_threshold * 100:.1f}% of max = {threshold_relative:.3e} Ohm")
        logger.info(f"  Absolute threshold: 0.1% of R_pol = {threshold_absolute:.3e} Ohm")
        logger.info(f"  Effective threshold: {threshold_effective:.3e} Ohm")

        keep_mask = np.abs(R_i) >= threshold_effective

        n_before = len(R_i)
        n_after = np.sum(keep_mask)
        n_pruned = n_before - n_after

        if n_pruned > 0:
            logger.info(f"  Removed {n_pruned} elements (kept {n_after}/{n_before})")
            R_i = R_i[keep_mask]
            tau = tau[keep_mask]
        else:
            logger.info("  No elements pruned (all above threshold)")
    else:
        logger.info("Step 3: Pruning disabled")

    # Step 4: Build circuit using K(R, tau) elements
    logger.info(f"Step 4: Building circuit with {len(R_i)} K elements")

    # Start with R_s if included
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

    # Add inductance if fitted
    if include_L and L_value is not None:
        circuit = circuit - L(L_value)
        logger.info(f"  Added inductance: L = {L_value:.3e} H")

    # Extract initial parameters
    initial_params = circuit.get_all_params()

    logger.info("")
    logger.info("Initial guess summary:")
    if include_Rs:
        logger.info(f"  R_s = {R_s:.3e} Ohm")
    if len(R_i) > 0:
        logger.info(f"  {len(R_i)} K elements:")
        logger.info(f"    R_i in [{R_i.min():.3e}, {R_i.max():.3e}] Ohm")
        logger.info(f"    tau_i in [{tau.min():.3e}, {tau.max():.3e}] s")
    if include_L and L_value is not None:
        logger.info(f"  Inductance: L = {L_value:.3e} H")
    logger.info(f"  Total parameters: {len(initial_params)}")
    logger.info("="*60)

    return circuit, initial_params


__all__ = ['estimate_R_linear', 'fit_voigt_chain_linear', 'VoigtChain']

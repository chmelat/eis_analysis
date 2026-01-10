"""
Analytic Jacobian computation for EIS circuit fitting.

This module provides analytic derivatives for circuit elements,
which can significantly improve numerical stability of least_squares
refinement, especially for ill-conditioned problems with parameters
spanning many orders of magnitude.

Theory
------
For least_squares minimization, the Jacobian J has elements:
    J[i,j] = d(residual_i) / d(param_j)

For complex impedance with real/imag residuals:
    residual = [Re(Z - Z_fit) * w, Im(Z - Z_fit) * w]

So:
    J[i,j] = -w * d(Re(Z_fit)) / d(param_j)  for real part
    J[i,j] = -w * d(Im(Z_fit)) / d(param_j)  for imag part

Element Derivatives
-------------------
R:   Z = R           -> dZ/dR = 1
C:   Z = 1/(jwC)     -> dZ/dC = -Z/C
L:   Z = jwL         -> dZ/dL = jw
Q:   Z = 1/(Q(jw)^n) -> dZ/dQ = -Z/Q, dZ/dn = -Z*ln(jw)
W:   Z = s(1-j)/sqrt(w) -> dZ/ds = Z/s
Wo:  Z = Rw*tanh(u)/u   -> dZ/dRw = Z/Rw, dZ/dtau = complex formula
K:   Z = R/(1+jwt)   -> dZ/dR = Z/R, dZ/dtau = -jw*Z^2/R
G:   Z = s/sqrt(1+jwt)  -> dZ/ds = 1/sqrt(1+jwt), dZ/dtau = -s*jw/(2*(1+jwt)^1.5)

Circuit Composition
-------------------
Series:   Z = sum(Zi)  -> dZ/dp = dZi/dp (for param p in element i)
Parallel: Z = 1/sum(1/Zi) -> dZ/dp = (Z/Zi)^2 * dZi/dp

Author: EIS Analysis Toolkit
"""

import numpy as np
from typing import List, Tuple, Union
from numpy.typing import NDArray

from .circuit_elements import R, C, L, Q, W, Wo, K, G, CircuitElement
from .circuit_builder import Series, Parallel, CompositeCircuit


def element_jacobian(
    element: CircuitElement,
    freq: NDArray[np.float64],
    params: List[float]
) -> Tuple[NDArray[np.complex128], NDArray[np.complex128]]:
    """
    Compute impedance and Jacobian for a single circuit element.

    Parameters
    ----------
    element : CircuitElement
        Circuit element (R, C, L, Q, W, Wo, K, G)
    freq : ndarray of float
        Frequencies [Hz]
    params : list of float
        Element parameters

    Returns
    -------
    Z : ndarray of complex, shape (n_freq,)
        Impedance at each frequency
    dZ : ndarray of complex, shape (n_freq, n_params)
        Jacobian: dZ/dp for each parameter

    Raises
    ------
    NotImplementedError
        If element type is not supported
    """
    omega = 2 * np.pi * freq
    n_freq = len(freq)

    # Resistor: Z = R
    if isinstance(element, R):
        R_val = params[0]
        Z = R_val * np.ones(n_freq, dtype=complex)
        dZ = np.ones((n_freq, 1), dtype=complex)
        return Z, dZ

    # Capacitor: Z = 1/(jwC)
    if isinstance(element, C):
        C_val = params[0]
        Z = 1 / (1j * omega * C_val)
        dZ_dC = -Z / C_val
        return Z, dZ_dC.reshape(-1, 1)

    # Inductor: Z = jwL
    if isinstance(element, L):
        L_val = params[0]
        Z = 1j * omega * L_val
        dZ_dL = 1j * omega
        return Z, dZ_dL.reshape(-1, 1)

    # Q (CPE): Z = 1/(Q*(jw)^n)
    if isinstance(element, Q):
        Q_val, n_val = params[0], params[1]
        jw_n = (1j * omega) ** n_val
        Z = 1 / (Q_val * jw_n)

        dZ_dQ = -Z / Q_val
        # dZ/dn = -Z * ln(jw)
        # ln(jw) = ln(w) + j*pi/2
        ln_jw = np.log(omega) + 1j * np.pi / 2
        dZ_dn = -Z * ln_jw

        dZ = np.column_stack([dZ_dQ, dZ_dn])
        return Z, dZ

    # Warburg: Z = sigma*(1-j)/sqrt(w)
    if isinstance(element, W):
        sigma_val = params[0]
        Z = sigma_val / np.sqrt(omega) * (1 - 1j)
        dZ_dsigma = Z / sigma_val
        return Z, dZ_dsigma.reshape(-1, 1)

    # Warburg open: Z = Rw * tanh(u) / u, where u = sqrt(jw*tau)
    if isinstance(element, Wo):
        Rw_val, tau_val = params[0], params[1]
        u = np.sqrt(1j * omega * tau_val)
        tanh_u = np.tanh(u)
        Z = Rw_val * tanh_u / u

        # dZ/dRw = Z/Rw
        dZ_dRw = Z / Rw_val

        # dZ/dtau: use chain rule
        # Z = Rw * tanh(u)/u
        # dZ/dtau = Rw * d/dtau[tanh(u)/u]
        # du/dtau = jw / (2*sqrt(jw*tau)) = jw / (2u)
        # d/du[tanh(u)/u] = (sech^2(u)*u - tanh(u)) / u^2
        #                 = (1 - tanh^2(u))/u - tanh(u)/u^2
        sech2_u = 1 - tanh_u**2
        d_tanh_u_over_u = (sech2_u * u - tanh_u) / u**2
        du_dtau = 1j * omega / (2 * u)
        dZ_dtau = Rw_val * d_tanh_u_over_u * du_dtau

        dZ = np.column_stack([dZ_dRw, dZ_dtau])
        return Z, dZ

    # K element (Voigt): Z = R/(1 + jw*tau)
    if isinstance(element, K):
        R_val, tau_val = params[0], params[1]
        denom = 1 + 1j * omega * tau_val
        Z = R_val / denom

        # dZ/dR = 1/(1 + jw*tau) = Z/R
        dZ_dR = Z / R_val

        # dZ/dtau = -R * jw / (1 + jw*tau)^2 = -jw * Z^2 / R
        dZ_dtau = -1j * omega * Z**2 / R_val

        dZ = np.column_stack([dZ_dR, dZ_dtau])
        return Z, dZ

    # Gerischer element: Z = sigma / sqrt(1 + jw*tau)
    if isinstance(element, G):
        sigma_val, tau_val = params[0], params[1]
        sqrt_term = np.sqrt(1 + 1j * omega * tau_val)
        Z = sigma_val / sqrt_term

        # dZ/dsigma = 1/sqrt(1 + jw*tau)
        dZ_dsigma = 1.0 / sqrt_term

        # dZ/dtau = -sigma * jw / (2 * (1 + jw*tau)^(3/2))
        dZ_dtau = -sigma_val * 1j * omega / (2 * (1 + 1j * omega * tau_val) ** 1.5)

        dZ = np.column_stack([dZ_dsigma, dZ_dtau])
        return Z, dZ

    raise NotImplementedError(
        f"Analytic Jacobian not implemented for {type(element).__name__}. "
        f"Use numeric Jacobian (use_analytic_jacobian=False)."
    )


def circuit_jacobian(
    circuit: Union[CircuitElement, CompositeCircuit],
    freq: NDArray[np.float64],
    params: List[float]
) -> Tuple[NDArray[np.complex128], NDArray[np.complex128]]:
    """
    Compute impedance and Jacobian for a circuit (element or composite).

    Recursively computes Jacobian for Series and Parallel circuits using
    chain rule.

    Parameters
    ----------
    circuit : CircuitElement or CompositeCircuit
        Circuit to evaluate
    freq : ndarray of float
        Frequencies [Hz]
    params : list of float
        All circuit parameters (flattened)

    Returns
    -------
    Z : ndarray of complex, shape (n_freq,)
        Total impedance at each frequency
    dZ : ndarray of complex, shape (n_freq, n_params)
        Jacobian: dZ/dp for each parameter

    Notes
    -----
    For Series:   Z = sum(Zi), so dZ/dp = dZi/dp
    For Parallel: Z = 1/Y where Y = sum(Yi), Yi = 1/Zi
                  dZ/dp = -Z^2 * dY/dp = -Z^2 * (-1/Zi^2) * dZi/dp
                        = (Z/Zi)^2 * dZi/dp
    """
    n_freq = len(freq)

    # Single element
    if isinstance(circuit, CircuitElement):
        return element_jacobian(circuit, freq, params)

    # Series connection: Z = sum(Zi)
    if isinstance(circuit, Series):
        n_params_total = len(params)
        Z_total = np.zeros(n_freq, dtype=complex)
        dZ_total = np.zeros((n_freq, n_params_total), dtype=complex)

        param_idx = 0
        for elem in circuit.elements:
            n_elem_params = len(elem.get_all_params())
            elem_params = params[param_idx:param_idx + n_elem_params]

            Z_elem, dZ_elem = circuit_jacobian(elem, freq, elem_params)

            Z_total += Z_elem
            dZ_total[:, param_idx:param_idx + n_elem_params] = dZ_elem

            param_idx += n_elem_params

        return Z_total, dZ_total

    # Parallel connection: Z = 1/Y, Y = sum(1/Zi)
    if isinstance(circuit, Parallel):
        n_params_total = len(params)
        Y_total = np.zeros(n_freq, dtype=complex)
        dZ_total = np.zeros((n_freq, n_params_total), dtype=complex)

        # First pass: compute Y_total and store element impedances/jacobians
        elem_data = []
        param_idx = 0
        for elem in circuit.elements:
            n_elem_params = len(elem.get_all_params())
            elem_params = params[param_idx:param_idx + n_elem_params]

            Z_elem, dZ_elem = circuit_jacobian(elem, freq, elem_params)
            Y_total += 1 / Z_elem
            elem_data.append((param_idx, n_elem_params, Z_elem, dZ_elem))

            param_idx += n_elem_params

        Z_total = 1 / Y_total

        # Second pass: compute dZ/dp using chain rule
        # dZ/dp = (Z/Zi)^2 * dZi/dp
        for param_idx, n_elem_params, Z_elem, dZ_elem in elem_data:
            factor = (Z_total / Z_elem) ** 2
            # factor is (n_freq,), dZ_elem is (n_freq, n_elem_params)
            dZ_total[:, param_idx:param_idx + n_elem_params] = (
                factor[:, np.newaxis] * dZ_elem
            )

        return Z_total, dZ_total

    raise NotImplementedError(
        f"Analytic Jacobian not implemented for {type(circuit).__name__}"
    )


def compute_analytic_jacobian(
    circuit: Union[CircuitElement, CompositeCircuit],
    freq: NDArray[np.float64],
    params: List[float],
    weights: NDArray[np.float64],
    fixed_params: List[bool] = None
) -> NDArray[np.float64]:
    """
    Compute the full Jacobian matrix for least_squares residuals.

    This function computes the Jacobian of the residual vector:
        residual = [Re(Z - Z_fit) * w, Im(Z - Z_fit) * w]

    The Jacobian is:
        J[i,j] = d(residual_i) / d(param_j)

    For the real part (first half):
        J[i,j] = -w[i] * d(Re(Z_fit)) / d(param_j)

    For the imaginary part (second half):
        J[i,j] = -w[i] * d(Im(Z_fit)) / d(param_j)

    Parameters
    ----------
    circuit : CircuitElement or CompositeCircuit
        Circuit to evaluate
    freq : ndarray of float
        Frequencies [Hz]
    params : list of float
        All circuit parameters (only free parameters if fixed_params given)
    weights : ndarray of float
        Weights for each frequency point
    fixed_params : list of bool, optional
        Which parameters are fixed (True = fixed, False = free)
        If given, params should contain only free parameters

    Returns
    -------
    J : ndarray of float, shape (2*n_freq, n_free_params)
        Jacobian matrix suitable for scipy.optimize.least_squares
    """
    # If we have fixed params, reconstruct full param list for evaluation
    if fixed_params is not None and any(fixed_params):
        # params contains only free params, need full list for circuit
        # This case is handled by caller - we receive full params here
        pass

    # Compute impedance and complex Jacobian
    _, dZ = circuit_jacobian(circuit, freq, list(params))

    # dZ is (n_freq, n_params) complex
    # Residual = [Re(Z_data - Z_fit)*w, Im(Z_data - Z_fit)*w]
    # d(Residual)/d(param) = [-Re(dZ/dp)*w, -Im(dZ/dp)*w]

    # Build real Jacobian from complex derivatives
    dZ_real = -dZ.real * weights[:, np.newaxis]
    dZ_imag = -dZ.imag * weights[:, np.newaxis]

    J = np.vstack([dZ_real, dZ_imag])

    # Filter out fixed parameters if needed
    if fixed_params is not None and any(fixed_params):
        free_mask = [not f for f in fixed_params]
        J = J[:, free_mask]

    return J


def make_jacobian_function(
    circuit: Union[CircuitElement, CompositeCircuit],
    freq: NDArray[np.float64],
    weights: NDArray[np.float64],
    fixed_params: List[bool] = None,
    full_initial_guess: List[float] = None
):
    """
    Create a Jacobian function compatible with scipy.optimize.least_squares.

    Returns a callable that takes free parameters and returns the Jacobian
    matrix.

    Parameters
    ----------
    circuit : CircuitElement or CompositeCircuit
        Circuit to fit
    freq : ndarray of float
        Frequencies [Hz]
    weights : ndarray of float
        Weights for each frequency point
    fixed_params : list of bool, optional
        Which parameters are fixed
    full_initial_guess : list of float, optional
        Full parameter list (needed to reconstruct fixed params)

    Returns
    -------
    jac_func : callable
        Function that takes free_params array and returns Jacobian matrix
    """
    def jacobian_func(free_params):
        # Reconstruct full params
        if fixed_params is None or not any(fixed_params):
            full_params = list(free_params)
        else:
            full_params = []
            free_idx = 0
            for i, is_fixed in enumerate(fixed_params):
                if is_fixed:
                    full_params.append(full_initial_guess[i])
                else:
                    full_params.append(free_params[free_idx])
                    free_idx += 1

        return compute_analytic_jacobian(
            circuit, freq, full_params, weights, fixed_params
        )

    return jacobian_func


__all__ = [
    'element_jacobian',
    'circuit_jacobian',
    'compute_analytic_jacobian',
    'make_jacobian_function',
]

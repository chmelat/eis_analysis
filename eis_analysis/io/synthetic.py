"""
Synthetic data generation for testing and demonstration.
"""

import numpy as np
import logging
from typing import Tuple
from numpy.typing import NDArray

logger = logging.getLogger(__name__)


def generate_synthetic_data(
    Rs: float = 10,
    R0: float = 1e5,
    Q0: Tuple[float, float] = (1e-6, 0.6),
    R1: float = 8e5,
    Q1: Tuple[float, float] = (3e-5, 0.43),
    noise: float = 0.01
) -> Tuple[NDArray[np.float64], NDArray[np.complex128]]:
    """
    Generate synthetic test data for Rs-(R0||CPE0)-(R1||CPE1) circuit.
    Use for testing the script without real data.

    Circuit structure: Rs - (R0||Q0) - (R1||Q1)
    where Q is a CPE (Constant Phase Element).

    Parameters
    ----------
    Rs : float
        Series resistance (electrolyte) [Ω]
    R0 : float
        First parallel resistance [Ω]
    Q0 : tuple (Y0, n)
        First CPE parameters: Y0 [S·s^n], n [-] (0 < n <= 1)
    R1 : float
        Second parallel resistance [Ω]
    Q1 : tuple (Y0, n)
        Second CPE parameters: Y0 [S·s^n], n [-] (0 < n <= 1)
    noise : float
        Noise level (0.01 = 1%)

    Returns
    -------
    frequencies : ndarray of float
        Frequency values [Hz]
    Z : ndarray of complex
        Complex impedance values [Ω]
    """
    logger.info("="*60)
    logger.info("Generating synthetic data")
    logger.info("="*60)

    Y0_0, n0 = Q0
    Y0_1, n1 = Q1

    logger.info("Circuit: Rs - (R0||Q0) - (R1||Q1)")
    logger.info(f"Rs = {Rs} Ω")
    logger.info(f"R0 = {R0:.2e} Ω, Q0 = ({Y0_0:.2e} S·s^n, n={n0})")
    logger.info(f"R1 = {R1:.2e} Ω, Q1 = ({Y0_1:.2e} S·s^n, n={n1})")

    # Frequencies
    frequencies = np.logspace(-2, 5, 70)  # 10 mHz - 100 kHz
    omega = 2 * np.pi * frequencies

    # CPE impedance: Z_CPE = 1 / (Y0 * (jω)^n)
    # Parallel R||CPE: Z = R / (1 + R * Y0 * (jω)^n)
    Z_RQ0 = R0 / (1 + R0 * Y0_0 * (1j * omega)**n0)
    Z_RQ1 = R1 / (1 + R1 * Y0_1 * (1j * omega)**n1)
    Z = Rs + Z_RQ0 + Z_RQ1

    # Add some noise
    noise_level = noise
    Z += noise_level * np.abs(Z) * (np.random.randn(len(Z)) +
                                      1j * np.random.randn(len(Z)))

    return frequencies, Z

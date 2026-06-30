"""
Circuit elements for EIS analysis using operator overloading.

This module replaces the old string parser from impedance.py with a modern
operator overloading approach inspired by EISAnalysis.jl (Julia).

Usage:
    from eis_analysis.fitting.circuit_elements import R, C, Q, W

    # Build circuits using operators
    circuit = R(100) - (R(5000) | C(1e-6))  # Voigt element
    randles = R(10) - (R(100) - W(50)) | Q(1e-4, 0.8)

    # Values serve as initial guess for fitting
    params = circuit.get_all_params()  # [100, 5000, 1e-6]

Operators:
    -  : Series connection (Z_total = Z1 + Z2)
    |  : Parallel connection (1/Z_total = 1/Z1 + 1/Z2)
    *  : Scale parameter (2*R(100) = R(200))
    ** : Exponent for Q/CPE (Q()**0.9)

Author: EIS Analysis Toolkit
Date: 2025-12-14

The element classes are organized into submodules:
    base        - CircuitElement abstract base class
    basic       - R, C, L (lumped ideal elements)
    distributed - Q, W, Wo (CPE and Warburg diffusion)
    composite   - K, G (Voigt R-τ and Gerischer)
"""
from .base import CircuitElement
from .basic import R, C, L
from .distributed import Q, W, Wo
from .composite import K, G

__all__ = ['R', 'C', 'Q', 'L', 'W', 'Wo', 'K', 'G', 'CircuitElement']

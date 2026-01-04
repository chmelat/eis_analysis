"""
Circuit combinators (Series and Parallel) for building complex circuits.

This module provides Series and Parallel classes that allow combining
circuit elements into complex topologies using operator overloading.

Author: EIS Analysis Toolkit
Date: 2025-12-14
"""
from __future__ import annotations

import numpy as np
from abc import ABC, abstractmethod
from typing import List, Union, Iterator, TYPE_CHECKING
from numpy.typing import NDArray

if TYPE_CHECKING:
    from .circuit_elements import CircuitElement


class CompositeCircuit(ABC):
    """
    Abstract base class for composite circuits (Series and Parallel).

    Implements shared functionality to eliminate code duplication between
    Series and Parallel classes.

    Attributes
    ----------
    elements : list
        List of circuit elements or sub-circuits
    """

    def __init__(self, elements: List):
        """
        Initialize composite circuit.

        Parameters
        ----------
        elements : list
            List of circuit elements or sub-circuits
        """
        self.elements = elements

    def get_all_params(self) -> List[float]:
        """
        Extract all parameters from all elements.

        Returns flat list of all parameters (initial guess for fitting).

        Returns
        -------
        params : list of float
            All parameters in order
        """
        params = []
        for elem in self.elements:
            if hasattr(elem, 'get_all_params'):
                params.extend(elem.get_all_params())
            else:
                raise TypeError(f"Element {elem} does not have get_all_params()")
        return params

    def get_all_fixed_params(self) -> List[bool]:
        """
        Extract fixed parameter flags from all elements.

        Returns flat list of booleans indicating which parameters are fixed.

        Returns
        -------
        fixed : list of bool
            True for fixed parameters, False for free parameters
        """
        fixed = []
        for elem in self.elements:
            if hasattr(elem, 'get_all_fixed_params'):
                fixed.extend(elem.get_all_fixed_params())
            elif hasattr(elem, 'fixed_params'):
                fixed.extend(elem.fixed_params)
            else:
                # Fallback: assume all params are free
                n_params = len(elem.get_all_params())
                fixed.extend([False] * n_params)
        return fixed

    def get_param_labels(self) -> List[str]:
        """
        Get parameter labels from all elements.

        Returns
        -------
        labels : list of str
            Parameter labels in order
        """
        labels = []
        for elem in self.elements:
            if hasattr(elem, 'get_param_labels'):
                labels.extend(elem.get_param_labels())
            else:
                # Fallback for elements without labels
                n_params = len(elem.get_all_params())
                labels.extend([f'p{i}' for i in range(n_params)])
        return labels

    def update_params(self, params: List[float]) -> int:
        """
        Update all element parameters with fitted values.

        Recursively updates parameters for all elements in order.

        Parameters
        ----------
        params : list of float
            New parameter values (flattened, same order as get_all_params)

        Returns
        -------
        n_consumed : int
            Total number of parameters consumed
        """
        param_idx = 0
        for elem in self.elements:
            if hasattr(elem, 'update_params'):
                n_consumed = elem.update_params(params[param_idx:])
                param_idx += n_consumed
            else:
                # Fallback: skip elements without update_params
                n_params = len(elem.get_all_params())
                param_idx += n_params
        return param_idx

    def _iter_element_impedances(
        self,
        freq: NDArray[np.float64],
        params: List[float]
    ) -> Iterator[NDArray[np.complex128]]:
        """
        Iterate over element impedances, extracting appropriate parameters.

        This generator handles parameter slicing and impedance calculation
        for each element, eliminating duplicate iteration logic.

        Parameters
        ----------
        freq : ndarray of float
            Frequencies [Hz]
        params : list of float
            All parameters (flattened)

        Yields
        ------
        Z_elem : ndarray of complex
            Impedance of each element [Ω]
        """
        param_idx = 0

        for elem in self.elements:
            # Get number of parameters for this element
            n_params = len(elem.get_all_params())

            # Extract parameters for this element
            elem_params = params[param_idx:param_idx + n_params]

            # Calculate and yield impedance
            yield elem.impedance(freq, elem_params)

            param_idx += n_params

    @abstractmethod
    def impedance(self, freq: NDArray[np.float64],
                  params: List[float]) -> NDArray[np.complex128]:
        """
        Calculate total impedance.

        Must be implemented by subclasses (Series, Parallel).

        Parameters
        ----------
        freq : ndarray of float
            Frequencies [Hz]
        params : list of float
            All parameters (flattened)

        Returns
        -------
        Z : ndarray of complex
            Total impedance [Ω]
        """
        pass


class Series(CompositeCircuit):
    """
    Series connection of circuit elements.

    Z_total = Z1 + Z2 + ... + Zn

    Inherits from CompositeCircuit which provides:
    - __init__, get_all_params, get_param_labels
    - _iter_element_impedances for parameter slicing

    Examples
    --------
    >>> from eis_analysis.fitting.circuit_elements import R, C
    >>> # Manual creation
    >>> series = Series([R(100), C(1e-6)])
    >>> # Using operator
    >>> series = R(100) - C(1e-6)
    """

    def impedance(self, freq: NDArray[np.float64],
                  params: List[float]) -> NDArray[np.complex128]:
        """
        Calculate total impedance of series connection.

        Parameters
        ----------
        freq : ndarray of float
            Frequencies [Hz]
        params : list of float
            All parameters (flattened)

        Returns
        -------
        Z : ndarray of complex
            Total impedance [Ω]
        """
        Z_total = np.zeros_like(freq, dtype=complex)

        # Use base class iterator for parameter handling
        for Z_elem in self._iter_element_impedances(freq, params):
            Z_total += Z_elem

        return Z_total

    def __sub__(self, other: Union['CircuitElement', 'Series', 'Parallel']) -> 'Series':
        """Add another element in series"""
        return Series(self.elements + [other])

    def __or__(self, other: Union['CircuitElement', 'Series', 'Parallel']) -> 'Parallel':
        """Combine in parallel"""
        return Parallel([self, other])

    def __repr__(self) -> str:
        """String representation"""
        parts = [str(elem) for elem in self.elements]
        return ' - '.join(parts)


class Parallel(CompositeCircuit):
    """
    Parallel connection of circuit elements.

    1/Z_total = 1/Z1 + 1/Z2 + ... + 1/Zn

    Inherits from CompositeCircuit which provides:
    - __init__, get_all_params, get_param_labels
    - _iter_element_impedances for parameter slicing

    Examples
    --------
    >>> from eis_analysis.fitting.circuit_elements import R, C
    >>> # Manual creation
    >>> parallel = Parallel([R(100), C(1e-6)])
    >>> # Using operator
    >>> parallel = R(100) | C(1e-6)
    """

    def impedance(self, freq: NDArray[np.float64],
                  params: List[float]) -> NDArray[np.complex128]:
        """
        Calculate total impedance of parallel connection.

        Parameters
        ----------
        freq : ndarray of float
            Frequencies [Hz]
        params : list of float
            All parameters (flattened)

        Returns
        -------
        Z : ndarray of complex
            Total impedance [Ω]
        """
        Y_total = np.zeros_like(freq, dtype=complex)  # admittance

        # Use base class iterator for parameter handling
        for Z_elem in self._iter_element_impedances(freq, params):
            Y_total += 1 / Z_elem

        return 1 / Y_total

    def __sub__(self, other: Union['CircuitElement', 'Series', 'Parallel']) -> 'Series':
        """Combine in series"""
        return Series([self, other])

    def __or__(self, other: Union['CircuitElement', 'Series', 'Parallel']) -> 'Parallel':
        """Add another element in parallel"""
        return Parallel(self.elements + [other])

    def __repr__(self) -> str:
        """String representation"""
        parts = [str(elem) for elem in self.elements]
        return '(' + ' | '.join(parts) + ')'


# Type alias for any circuit type
Circuit = Union[Series, Parallel, 'CircuitElement']


__all__ = ['CompositeCircuit', 'Series', 'Parallel', 'Circuit']

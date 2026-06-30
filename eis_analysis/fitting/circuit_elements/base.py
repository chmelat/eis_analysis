"""
Abstract base class for circuit elements.

Defines ``CircuitElement`` with the operator overloading (``-``, ``|``, ``*``)
used to build circuits. Concrete elements live in the sibling modules
``basic``, ``distributed`` and ``composite``.
"""
from __future__ import annotations

import numpy as np
from typing import List, Union, TYPE_CHECKING
from abc import ABC, abstractmethod
from numpy.typing import NDArray

if TYPE_CHECKING:
    from ..circuit_builder import Series, Parallel
    Circuit = Union[Series, Parallel, 'CircuitElement']


class CircuitElement(ABC):
    """
    Abstract base class for all circuit elements.

    Attributes
    ----------
    params : list of float
        Element parameters (serve as initial guess for fitting)
    fixed_params : list of bool
        Whether each parameter is fixed (True) or free to fit (False)

    Notes
    -----
    Parameters passed as strings are automatically treated as fixed.
    Example: R("0.86") creates a fixed resistor, R(0.86) is free to fit.
    """

    def __init__(self, *params):
        # Convert string parameters to float and mark as fixed
        self.params = []
        self.fixed_params = []
        for param in params:
            if isinstance(param, str):
                self.params.append(float(param))
                self.fixed_params.append(True)
            else:
                self.params.append(float(param))
                self.fixed_params.append(False)

    @abstractmethod
    def impedance(self, freq: NDArray[np.float64],
                  params: List[float]) -> NDArray[np.complex128]:
        """
        Calculate impedance Z(ω) with given parameters.

        Parameters
        ----------
        freq : ndarray of float
            Frequencies [Hz]
        params : list of float
            Element parameters

        Returns
        -------
        Z : ndarray of complex
            Impedance [Ω]
        """
        pass

    @abstractmethod
    def _scale(self, scalar: float) -> 'CircuitElement':
        """Return scaled version of element (for * operator)"""
        pass

    @abstractmethod
    def __repr__(self) -> str:
        """String representation"""
        pass

    # Operator overloading for circuit building
    def __sub__(self, other: Union['CircuitElement', 'Circuit']) -> 'Circuit':
        """Series connection: Z_total = Z1 + Z2"""
        from ..circuit_builder import Series
        return Series([self, other])

    def __or__(self, other: Union['CircuitElement', 'Circuit']) -> 'Circuit':
        """Parallel connection: 1/Z_total = 1/Z1 + 1/Z2"""
        from ..circuit_builder import Parallel
        return Parallel([self, other])

    def __mul__(self, scalar: float) -> 'CircuitElement':
        """Scale element parameter"""
        return self._scale(scalar)

    def __rmul__(self, scalar: float) -> 'CircuitElement':
        """Allow 2*R syntax"""
        return self.__mul__(scalar)

    def get_all_params(self) -> List[float]:
        """Get all parameters (for compatibility with Circuit class)"""
        return self.params

    def get_all_fixed_params(self) -> List[bool]:
        """Get fixed parameter flags (for compatibility with Circuit class)"""
        return self.fixed_params

    def update_params(self, params: List[float]) -> int:
        """
        Update element parameters with fitted values.

        Parameters
        ----------
        params : list of float
            New parameter values (must match number of element parameters)

        Returns
        -------
        n_consumed : int
            Number of parameters consumed from the list
        """
        n_params = len(self.params)
        self.params = list(params[:n_params])
        return n_params

    @abstractmethod
    def get_param_labels(self) -> List[str]:
        """
        Get parameter labels for display.

        Returns
        -------
        labels : list of str
            Parameter labels (e.g., ['R'], ['Q', 'n'], etc.)
        """
        pass

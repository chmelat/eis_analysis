"""
Data structures for DRT (Distribution of Relaxation Times) analysis.

All result containers returned by the DRT pipeline. No logic — diagnostics are
returned as structured data, the CLI layer is responsible for user output.
"""

import matplotlib.pyplot as plt
import numpy as np
from dataclasses import dataclass, field
from typing import Optional, List, Dict
from numpy.typing import NDArray


@dataclass
class RinfEstimate:
    """Result of R_inf estimation."""
    R_inf: float
    method: str  # 'preset', 'median', 'rl_fit', 'voigt_fit'
    R_inf_median: Optional[float] = None  # For comparison
    figure: Optional[plt.Figure] = None

    # Fit diagnostics (if applicable)
    behavior: Optional[str] = None  # 'capacitive', 'inductive', 'mixed'
    n_points_used: Optional[int] = None
    R_squared: Optional[float] = None
    L_nH: Optional[float] = None
    R_ct: Optional[float] = None
    C_nF: Optional[float] = None
    f_characteristic: Optional[float] = None

    warnings: List[str] = field(default_factory=list)


@dataclass
class LambdaSelection:
    """Result of lambda selection."""
    lambda_value: float
    method: str  # 'user', 'default', 'gcv', 'hybrid', 'fallback'
    lambda_gcv: Optional[float] = None  # If L-curve correction was applied
    gcv_score: Optional[float] = None
    corner_at_edge: bool = False   # L-curve corner landed at edge of search range (F7)
    lambda_at_edge: bool = False   # selected lambda hit a bound of the GCV range (F3/F7)


@dataclass
class NNLSSolution:
    """Result of NNLS solver."""
    gamma: Optional[NDArray[np.float64]]
    success: bool
    n_inductive_points: int = 0
    inductive_fraction: float = 0.0
    max_inductive_imag: float = 0.0
    gamma_max: Optional[float] = None
    gamma_min_nonzero: Optional[float] = None
    condition_number: float = 0.0
    warnings: List[str] = field(default_factory=list)


@dataclass
class DRTDiagnostics:
    """Comprehensive DRT diagnostics."""
    # Frequency info
    freq_min: float
    freq_max: float
    freq_range_ratio: float
    n_points: int

    # Matrix info
    n_tau: int
    condition_number: float  # of the regularized system [A; sqrt(lambda)*L]
    d_ln_tau: float

    # R_inf estimation
    rinf: RinfEstimate

    # Lambda selection
    lambda_sel: LambdaSelection

    # NNLS solution
    nnls: NNLSSolution

    # R_pol info
    R_pol_from_data: float
    R_pol_from_gamma: float
    normalized: bool

    # Reconstruction
    reconstruction_error_rel: float

    # Peak detection
    peak_method: str
    n_peaks: int
    scipy_peaks: Optional[List[Dict]] = None  # For scipy method

    # Shape diagnostics (F3): effective number of gamma bins (participation ratio)
    n_effective_bins: Optional[float] = None


@dataclass
class DRTMatrices:
    """Container for DRT matrices."""
    A: NDArray[np.float64]
    A_re: NDArray[np.float64]
    A_im: NDArray[np.float64]
    b: NDArray[np.float64]
    L: NDArray[np.float64]
    tau: NDArray[np.float64]
    d_ln_tau: float
    condition_number: float


@dataclass
class DRTResult:
    """
    Complete DRT analysis result.

    All information previously logged is now available as structured data.
    """
    # Core results
    tau: Optional[NDArray[np.float64]] = None
    gamma: Optional[NDArray[np.float64]] = None
    gamma_original: Optional[NDArray[np.float64]] = None

    # Figures
    figure: Optional[plt.Figure] = None
    figure_rinf: Optional[plt.Figure] = None

    # GMM peaks (if GMM method used)
    peaks: Optional[List[Dict]] = None
    bic_scores: Optional[List[float]] = None

    # Key values for quick access
    R_inf: Optional[float] = None
    R_pol: Optional[float] = None
    lambda_used: Optional[float] = None
    reconstruction_error: Optional[float] = None

    # Full diagnostics
    diagnostics: Optional[DRTDiagnostics] = None

    @property
    def success(self) -> bool:
        """Check if DRT calculation was successful."""
        return self.tau is not None and self.gamma is not None

    @property
    def warnings(self) -> List[str]:
        """Collect all warnings from diagnostics."""
        if self.diagnostics is None:
            return []
        warnings = []
        warnings.extend(self.diagnostics.rinf.warnings)
        warnings.extend(self.diagnostics.nnls.warnings)
        return warnings

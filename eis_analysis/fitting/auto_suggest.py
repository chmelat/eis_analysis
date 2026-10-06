"""
Voigt element analysis from DRT spectra.

This module analyzes DRT (Distribution of Relaxation Times) spectra and
identifies Voigt elements (R||C) with parameter estimates for each peak.
Does NOT generate circuit strings - provides information for manual circuit building.
"""

import numpy as np
import logging
from dataclasses import dataclass, field
from typing import List, Optional, Dict, Sequence
from numpy.typing import NDArray

from ..drt.estimation import _estimate_peak_resistance, refine_peak_tau
from ..utils.impedance import calculate_rpol
from .config import (
    MAX_VOIGT_ELEMENTS,
    RPOL_RATIO_WARNING_THRESHOLD_LOW,
    RPOL_RATIO_WARNING_THRESHOLD_HIGH,
)

logger = logging.getLogger(__name__)


@dataclass
class VoigtElement:
    """One relaxation process read off a DRT peak."""
    id: int              # 1-based, ordered by tau
    tau: float           # time constant [s]
    freq: float          # characteristic frequency [Hz]
    R: float             # resistance from the peak area [Ohm]
    C: float             # tau / R [F]
    warnings: List[str] = field(default_factory=list)


@dataclass
class VoigtSuggestion:
    """
    Voigt elements read off a DRT spectrum, with the diagnostics behind them.

    Carries the counts the report is built from (`n_peaks_raw` before the
    edge and height filters, `n_peaks_valid` after) because the reader
    cannot tell from the element list alone that peaks were dropped.

    Attributes
    ----------
    elements : list of VoigtElement
        The suggested elements, ordered by tau
    quality : str
        'good', 'acceptable', 'uncertain' or 'poor'
    total_R : float
        Sum of the element resistances [Ohm]
    R_pol, R_inf : float
        Polarization and high-frequency resistance from the data [Ohm]
    ratio : float
        total_R / R_pol; inf when R_pol is zero
    method : str
        'gmm' or 'scipy' - which peak detection produced the elements
    n_peaks_raw, n_peaks_valid : int
        Peaks detected, and peaks the suggestion is built from
    warnings : list of str
        Caveats about the analysis. `quality` is derived from how many
        there are, which is why a dropped peak is not one of them - see
        `excluded_peaks`.
    excluded_peaks : list of str
        Why individual peaks were dropped, one entry each
    """
    elements: List[VoigtElement]
    quality: str
    total_R: float
    R_pol: float
    R_inf: float
    ratio: float
    method: str
    n_peaks_raw: int = 0
    n_peaks_valid: int = 0
    warnings: List[str] = field(default_factory=list)
    excluded_peaks: List[str] = field(default_factory=list)


def analyze_voigt_elements(
    tau: NDArray[np.float64],
    gamma: NDArray[np.float64],
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    peaks_gmm: Optional[List[Dict]] = None,
    peak_indices: Optional[Sequence[int]] = None
) -> VoigtSuggestion:
    """
    Analyze Voigt elements (R||C) from DRT spectrum.

    Identifies relaxation processes in DRT and estimates R and C parameters
    for each peak. Provides quality diagnostics and warnings.

    NOTE: This function does NOT generate circuit strings (breaking change from v3.0.0).
    It only reports estimated parameters for individual elements.

    Parameter estimates based on:
    - R_i: area of gamma over the peak's valley-to-valley segment (on the
      scipy path the DRT R_estimate; GMM's R_estimate is R_pol * weight)
    - tau_i: peak position
    - C_i = tau_i / R_i

    Parameters
    ----------
    tau : ndarray of float
        Time constants from DRT [s] (M points)
    gamma : ndarray of float
        Distribution function from DRT [Ohm] (M points)
    frequencies : ndarray of float
        Original frequencies [Hz] (N points)
    Z : ndarray of complex
        Original impedance [Ohm] (N points)
    peaks_gmm : list of dict, optional
        GMM peaks from gmm_peak_detection(). If provided (and non-empty),
        used instead of peak_indices.
    peak_indices : sequence of int, optional
        Grid indices of the DRT's significant maxima,
        ``[p['index'] for p in drt.diagnostics.scipy_peaks]``. The peaks are
        decided once, by calculate_drt's significance test, so the suggestion
        cannot differ from the peaks it reports. One of peaks_gmm and
        peak_indices is required.

    Returns
    -------
    VoigtSuggestion
        The elements, the counts behind them and any caveat. Nothing is
        logged; the CLI composes its section from this.

    Notes
    -----
    Analysis is based on empirical rules:
    1. Each peak in DRT corresponds to a relaxation process (R||C)
    2. Maximum number of elements is MAX_VOIGT_ELEMENTS (4) - more = overfit
    3. Peaks at tau range edges are filtered (artifacts)

    Examples
    --------
    >>> drt = calculate_drt(freq, Z)
    >>> suggestion = analyze_voigt_elements(
    ...     drt.tau, drt.gamma, freq, Z,
    ...     peak_indices=[p['index'] for p in drt.diagnostics.scipy_peaks])
    >>> print(f"Found {len(suggestion.elements)} Voigt elements")
    >>> print(f"Analysis quality: {suggestion.quality}")

    See Also
    --------
    config.MAX_VOIGT_ELEMENTS : Maximum number of parallel RC elements
    drt.significance : how calculate_drt decides which maxima are peaks
    """
    warnings: List[str] = []
    excluded_peaks: List[str] = []
    method = 'gmm' if peaks_gmm is not None else 'scipy'

    # Basic characteristics from data (DRY: uses utils.impedance)
    n_avg = min(5, max(1, len(frequencies) // 10))
    R_pol_data, R_inf, R_dc = calculate_rpol(frequencies, Z, n_avg)

    # Find peaks - either from GMM or scipy
    # Peak tau by grid index, used for every message and element: GMM centers
    # stay as fitted, scipy maxima are refined between nodes (tau[idx] is off
    # by up to half a grid step)
    peak_tau: Dict[int, float] = {}
    if peaks_gmm is not None and len(peaks_gmm) > 0:
        # Convert GMM peaks to format compatible with rest of function
        # Find nearest index in tau for each GMM peak
        gmm_indices = []
        for peak_gmm in peaks_gmm:
            tau_center = peak_gmm['tau_center']
            idx = int(np.argmin(np.abs(tau - tau_center)))
            gmm_indices.append(idx)
            peak_tau[idx] = float(tau_center)
        peaks = np.array(gmm_indices)
    elif peak_indices is not None:
        peaks = np.asarray(peak_indices, dtype=int)
    else:
        raise ValueError("analyze_voigt_elements needs peaks_gmm or peak_indices "
                         "(the DRT's scipy_peaks indices)")

    n_peaks_raw = len(peaks)

    if len(peaks) == 0:
        warnings.append("DRT contains no distinct peaks")

        # Fallback: single Voigt element estimated from -Z'' maximum
        idx_max_zimag = np.argmax(-Z.imag)
        f_char = frequencies[idx_max_zimag]
        tau_char = 1 / (2 * np.pi * f_char)
        C_est = tau_char / R_pol_data if R_pol_data > 0 else 1e-6

        return VoigtSuggestion(
            elements=[VoigtElement(
                id=1, tau=tau_char, freq=f_char, R=R_pol_data, C=C_est,
                warnings=["estimated from -Z'' maximum (no DRT peaks)"])],
            quality='poor',
            total_R=R_pol_data,
            R_pol=R_pol_data,
            R_inf=R_inf,
            ratio=1.0,
            method=method,
            n_peaks_raw=n_peaks_raw,
            n_peaks_valid=0,
            warnings=warnings,
            excluded_peaks=excluded_peaks)

    for peak in peaks:
        if int(peak) not in peak_tau:
            peak_tau[int(peak)] = refine_peak_tau(tau, gamma, int(peak))

    # Filter peaks at edges (may be artifacts or truncated)
    edge_margin = max(2, len(tau) // 20)  # 5% from edge
    valid_peaks = []

    for peak in peaks:
        peak_info = {
            'index': peak,
            'tau': peak_tau[peak],
            'gamma': gamma[peak],
            'valid': True,
            'warnings': []
        }

        # Edge check
        if peak < edge_margin:
            peak_info['warnings'].append('near left edge (high f)')
            peak_info['valid'] = False
        elif peak > len(tau) - edge_margin:
            peak_info['warnings'].append('near right edge (low f)')
            peak_info['valid'] = False

        if peak_info['valid']:
            valid_peaks.append(peak)
        else:
            # Surface why a detected DRT peak is dropped from the circuit
            # suggestion (otherwise the count silently shrinks, e.g. 2 -> 1).
            f_peak = 1 / (2 * np.pi * peak_tau[peak])
            reason = '; '.join(peak_info['warnings'])
            excluded_peaks.append(
                f"Peak at tau = {peak_tau[peak]:.2e} s (f = {f_peak:.2e} Hz) "
                f"excluded: {reason}"
            )

    # If all peaks were at edges, use at least the highest one
    if len(valid_peaks) == 0 and len(peaks) > 0:
        highest_peak = peaks[np.argmax(gamma[peaks])]
        valid_peaks = [highest_peak]
        warnings.append(f"All peaks at tau range edges, using the highest "
                        f"at tau = {peak_tau[highest_peak]:.2e} s")

    n_peaks_valid = len(valid_peaks)

    # Sort peaks by tau (smallest to largest)
    valid_peaks = sorted(valid_peaks, key=lambda p: peak_tau[p])

    # Limit number of Voigt elements (config.MAX_VOIGT_ELEMENTS)
    if len(valid_peaks) > MAX_VOIGT_ELEMENTS:
        warnings.append(
            f"Too many peaks ({len(valid_peaks)}), limited to "
            f"{MAX_VOIGT_ELEMENTS} most prominent"
        )
        # Select most prominent peaks
        peak_heights = [gamma[p] for p in valid_peaks]
        top_indices = np.argsort(peak_heights)[-MAX_VOIGT_ELEMENTS:]
        valid_peaks = sorted([valid_peaks[i] for i in top_indices], key=lambda p: peak_tau[p])

    # Calculate Voigt elements from peaks
    n_voigt = len(valid_peaks)

    # R_i over the valley partition of all detected peaks: the R_estimate the
    # DRT peak list prints. Walking out to a fraction of the peak height
    # instead climbed across a shallow valley into a taller neighbour.
    ordered = np.sort(peaks)
    peak_R = dict(zip(ordered.tolist(), _estimate_peak_resistance(tau, gamma, ordered)))
    total_R_from_peaks = 0
    elements = []

    for i, peak in enumerate(valid_peaks):
        tau_i = peak_tau[peak]
        f_i = 1 / (2 * np.pi * tau_i)

        R_i = peak_R[int(peak)]

        # Element warnings
        elem_warnings = []

        # Fallback only if the segment holds no mass at all. Not a floor in
        # Ohm: the significance test keeps small resolved peaks, and a fixed
        # 1 Ohm turned a 0.5 Ohm arc into R_pol / n (and mOhm cells always).
        if R_i <= 0:
            R_i = R_pol_data / n_voigt
            elem_warnings.append('heuristic R estimate (integration failed)')
            warnings.append(f"Peak {i+1}: used heuristic R estimate")

        total_R_from_peaks += R_i

        # C_i = tau_i / R_i
        C_i = tau_i / R_i

        # Clamp C to reasonable range
        if C_i < 1e-12:
            C_i = 1e-12
            elem_warnings.append('C clamped to lower limit (1e-12 F)')
        elif C_i > 1e-1:
            C_i = 1e-1
            elem_warnings.append('C clamped to upper limit (1e-1 F)')

        elements.append(VoigtElement(
            id=i + 1, tau=tau_i, freq=f_i, R=R_i, C=C_i,
            warnings=elem_warnings))

    # Consistency check: sum of R_i should be close to R_pol
    if total_R_from_peaks > 0:
        ratio = R_pol_data / total_R_from_peaks
        if ratio < RPOL_RATIO_WARNING_THRESHOLD_LOW or ratio > RPOL_RATIO_WARNING_THRESHOLD_HIGH:
            warnings.append(
                f"Inconsistency: sum(R_i) = {total_R_from_peaks:.1f} Ohm "
                f"vs R_pol = {R_pol_data:.1f} Ohm"
            )

    # Quality assessment
    if len(warnings) == 0:
        quality = 'good'
    elif len(warnings) <= 2:
        quality = 'acceptable'
    else:
        quality = 'uncertain'

    ratio = total_R_from_peaks / R_pol_data if R_pol_data > 0 else float('inf')

    return VoigtSuggestion(
        elements=elements,
        quality=quality,
        total_R=total_R_from_peaks,
        R_pol=R_pol_data,
        R_inf=R_inf,
        ratio=ratio,
        method=method,
        n_peaks_raw=n_peaks_raw,
        n_peaks_valid=n_peaks_valid,
        warnings=warnings,
        excluded_peaks=excluded_peaks)


__all__ = ['analyze_voigt_elements', 'VoigtSuggestion', 'VoigtElement']

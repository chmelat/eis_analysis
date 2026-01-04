"""
Voigt element analysis from DRT spectra.

This module analyzes DRT (Distribution of Relaxation Times) spectra and
identifies Voigt elements (R||C) with parameter estimates for each peak.
Does NOT generate circuit strings - provides information for manual circuit building.
"""

import numpy as np
import logging
from typing import List, Optional, Dict
from numpy.typing import NDArray
from scipy.signal import find_peaks

from ..utils.compat import np_trapz
from ..utils.impedance import calculate_rpol
from ..drt.term_classification import classify_all_peaks
from .config import (
    DRT_PEAK_HEIGHT_THRESHOLD,
    DRT_PEAK_PROMINENCE_THRESHOLD,
    GMM_PEAK_HEIGHT_FACTOR,
    MAX_VOIGT_ELEMENTS,
    PEAK_INTEGRATION_TOLERANCE,
    RPOL_RATIO_WARNING_THRESHOLD_LOW,
    RPOL_RATIO_WARNING_THRESHOLD_HIGH,
)

logger = logging.getLogger(__name__)


def analyze_voigt_elements(
    tau: NDArray[np.float64],
    gamma: NDArray[np.float64],
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    peaks_gmm: Optional[List[Dict]] = None,
    classify_terms: bool = False
) -> Dict:
    """
    Analyze Voigt elements (R||C) from DRT spectrum.

    Identifies relaxation processes in DRT and estimates R and C parameters
    for each peak. Provides quality diagnostics and warnings.

    NOTE: This function does NOT generate circuit strings (breaking change from v3.0.0).
    It only reports estimated parameters for individual elements.

    Parameter estimates based on:
    - R_i: peak area (integral of gamma over peak, or R_estimate from GMM)
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
        GMM peaks from gmm_peak_detection()
        If provided, used instead of scipy.find_peaks
    classify_terms : bool, optional
        If True, classifies term types (Voigt/Q/other) based on peak shape.
        Requires peaks_gmm (GMM detection). Default: False

    Returns
    -------
    voigt_info : dict
        Information about Voigt elements:
        - 'elements': list of dict, each contains:
            - 'id': int (1-based index)
            - 'tau': float [s]
            - 'freq': float [Hz]
            - 'R': float [Ohm]
            - 'C': float [F]
            - 'warnings': list of str
            - 'classification': dict or None (if classify_terms=True):
                - 'type': str ('voigt', 'cpe', 'multiple', 'uncertain')
                - 'confidence': str ('high', 'medium', 'low')
                - 'cpe_n_estimate': float or None
                - 'reasoning': str
        - 'quality': str ('good', 'acceptable', 'uncertain', 'poor')
        - 'total_R': float (sum R_i) [Ohm]
        - 'R_pol': float (from data) [Ohm]
        - 'R_inf': float (from data) [Ohm]
        - 'ratio': float (total_R / R_pol)
        - 'warnings': list of str (global warnings)
        - 'method': str ('gmm' or 'scipy')

    Notes
    -----
    Analysis is based on empirical rules:
    1. Each peak in DRT corresponds to a relaxation process (R||C)
    2. Maximum number of elements is MAX_VOIGT_ELEMENTS (4) - more = overfit
    3. Peaks at tau range edges are filtered (artifacts)

    Examples
    --------
    >>> tau, gamma, fig, peaks, _ = calculate_drt(freq, Z)
    >>> voigt_info = analyze_voigt_elements(tau, gamma, freq, Z, peaks)
    >>> print(f"Found {len(voigt_info['elements'])} Voigt elements")
    >>> print(f"Analysis quality: {voigt_info['quality']}")

    See Also
    --------
    config.MAX_VOIGT_ELEMENTS : Maximum number of parallel RC elements
    config.DRT_PEAK_HEIGHT_THRESHOLD : Minimum peak height (10% of maximum)
    format_voigt_report : Format output for terminal
    """
    logger.info("="*60)
    logger.info("Automatic circuit suggestion from DRT")
    logger.info("="*60)

    diagnostics = {
        'n_peaks_raw': 0,
        'n_peaks_valid': 0,
        'peaks_info': [],
        'warnings': [],
        'quality': 'unknown',
        'method': 'gmm' if peaks_gmm is not None else 'scipy'
    }

    # Basic characteristics from data (DRY: uses utils.impedance)
    n_avg = min(5, max(1, len(frequencies) // 10))
    R_pol_data, R_inf, R_dc = calculate_rpol(frequencies, Z, n_avg)

    logger.info(f"R_inf (from data) = {R_inf:.2f} Ohm")
    logger.info(f"R_pol (from data) = {R_pol_data:.2f} Ohm")

    # Find peaks - either from GMM or scipy
    if peaks_gmm is not None and len(peaks_gmm) > 0:
        # Use GMM peaks
        logger.info(f"Using GMM peak detection ({len(peaks_gmm)} peaks)")

        # Convert GMM peaks to format compatible with rest of function
        # Find nearest index in tau for each GMM peak
        peaks = []
        for peak_gmm in peaks_gmm:
            tau_center = peak_gmm['tau_center']
            idx = np.argmin(np.abs(tau - tau_center))
            peaks.append(idx)
        peaks = np.array(peaks)
        properties = {}  # GMM doesn't need properties from find_peaks
    else:
        # Use scipy.find_peaks (original method)
        logger.info("Using scipy.find_peaks peak detection")
        min_distance = max(3, len(tau) // 20)  # At least 5% of spectrum width
        peaks, properties = find_peaks(
            gamma,
            height=np.max(gamma) * DRT_PEAK_HEIGHT_THRESHOLD,
            distance=min_distance,
            prominence=np.max(gamma) * DRT_PEAK_PROMINENCE_THRESHOLD
        )

    diagnostics['n_peaks_raw'] = len(peaks)
    logger.info(f"Found {len(peaks)} peaks in DRT spectrum")

    if len(peaks) == 0:
        diagnostics['warnings'].append("No peaks found")
        diagnostics['quality'] = 'poor'
        logger.warning("DRT contains no distinct peaks")

        # Fallback: single Voigt element estimated from -Z'' maximum
        idx_max_zimag = np.argmax(-Z.imag)
        f_char = frequencies[idx_max_zimag]
        tau_char = 1 / (2 * np.pi * f_char)
        C_est = tau_char / R_pol_data if R_pol_data > 0 else 1e-6

        return {
            'elements': [
                {
                    'id': 1,
                    'tau': tau_char,
                    'freq': f_char,
                    'R': R_pol_data,
                    'C': C_est,
                    'warnings': ['estimated from -Z\'\' maximum (no DRT peaks)'],
                    'classification': None
                }
            ],
            'quality': 'poor',
            'total_R': R_pol_data,
            'R_pol': R_pol_data,
            'R_inf': R_inf,
            'ratio': 1.0,
            'warnings': diagnostics['warnings'],
            'method': diagnostics['method']
        }

    # Filter peaks at edges (may be artifacts or truncated)
    edge_margin = max(2, len(tau) // 20)  # 5% from edge
    valid_peaks = []

    for peak in peaks:
        peak_info = {
            'index': peak,
            'tau': tau[peak],
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

        # Height check (very small peaks may be noise)
        if gamma[peak] < np.max(gamma) * GMM_PEAK_HEIGHT_FACTOR:
            peak_info['warnings'].append(f'low height (<{GMM_PEAK_HEIGHT_FACTOR*100:.0f}% of max)')
            # Don't mark as invalid, just warn

        diagnostics['peaks_info'].append(peak_info)

        if peak_info['valid']:
            valid_peaks.append(peak)

    diagnostics['n_peaks_valid'] = len(valid_peaks)

    # If all peaks were at edges, use at least the highest one
    if len(valid_peaks) == 0 and len(peaks) > 0:
        diagnostics['warnings'].append("All peaks at edges, using highest")
        highest_peak = peaks[np.argmax(gamma[peaks])]
        valid_peaks = [highest_peak]
        logger.warning("All peaks are at tau range edges")
        logger.warning(f"Using highest peak at tau = {tau[highest_peak]:.2e} s")

    # Sort peaks by tau (smallest to largest)
    valid_peaks = sorted(valid_peaks, key=lambda p: tau[p])

    logger.info(f"Valid peaks for circuit suggestion: {len(valid_peaks)}")

    # Limit number of Voigt elements (config.MAX_VOIGT_ELEMENTS)
    if len(valid_peaks) > MAX_VOIGT_ELEMENTS:
        diagnostics['warnings'].append(
            f"Too many peaks ({len(valid_peaks)}), limited to {MAX_VOIGT_ELEMENTS}"
        )
        logger.warning(f"Found {len(valid_peaks)} peaks, limiting to {MAX_VOIGT_ELEMENTS} most prominent")
        # Select most prominent peaks
        peak_heights = [gamma[p] for p in valid_peaks]
        top_indices = np.argsort(peak_heights)[-MAX_VOIGT_ELEMENTS:]
        valid_peaks = sorted([valid_peaks[i] for i in top_indices], key=lambda p: tau[p])

    # Term type classification (if requested and we have GMM peaks)
    classifications = None
    if classify_terms:
        if peaks_gmm is not None and len(peaks_gmm) > 0:
            logger.info("Term type classification activated")
            classifications = classify_all_peaks(peaks_gmm, tau, gamma)
        else:
            logger.warning("Term classification requires GMM peak detection (--peak-method gmm)")
            logger.warning("Skipping term classification")

    # Calculate Voigt elements from peaks
    n_voigt = len(valid_peaks)
    logger.info(f"Analyzing {n_voigt} Voigt elements")

    ln_tau = np.log(tau)
    total_R_from_peaks = 0
    elements = []

    for i, peak in enumerate(valid_peaks):
        tau_i = tau[peak]
        f_i = 1 / (2 * np.pi * tau_i)

        # Estimate R_i from peak area
        # Find peak boundaries (where gamma drops to PEAK_INTEGRATION_TOLERANCE of peak height)
        peak_height = gamma[peak]
        threshold = peak_height * PEAK_INTEGRATION_TOLERANCE

        # Left boundary
        left = peak
        while left > 0 and gamma[left] > threshold:
            left -= 1

        # Right boundary
        right = peak
        while right < len(gamma) - 1 and gamma[right] > threshold:
            right += 1

        # Integral over peak (trapezoidal method)
        R_i = np_trapz(gamma[left:right+1], ln_tau[left:right+1])

        # Element warnings
        elem_warnings = []

        # Fallback if R_i is too small or negative
        if R_i < 1:
            R_i = R_pol_data / n_voigt
            elem_warnings.append('heuristic R estimate (integration failed)')
            diagnostics['warnings'].append(f"Peak {i+1}: used heuristic R estimate")

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

        # Create element dict
        element = {
            'id': i + 1,
            'tau': tau_i,
            'freq': f_i,
            'R': R_i,
            'C': C_i,
            'warnings': elem_warnings
        }

        # Add classification if available
        if classifications is not None and i < len(classifications):
            element['classification'] = classifications[i]
        else:
            element['classification'] = None

        elements.append(element)

        logger.info(f"  Element {i+1}: tau = {tau_i:.2e} s, f = {f_i:.2e} Hz, R = {R_i:.1f} Ohm, C = {C_i:.2e} F")

    # Consistency check: sum of R_i should be close to R_pol
    if total_R_from_peaks > 0:
        ratio = R_pol_data / total_R_from_peaks
        if ratio < RPOL_RATIO_WARNING_THRESHOLD_LOW or ratio > RPOL_RATIO_WARNING_THRESHOLD_HIGH:
            diagnostics['warnings'].append(
                f"Inconsistency: sum(R_i) = {total_R_from_peaks:.1f} Ohm vs R_pol = {R_pol_data:.1f} Ohm"
            )
            logger.warning(
                f"Sum of resistances from peaks ({total_R_from_peaks:.1f} Ohm) "
                f"differs from R_pol ({R_pol_data:.1f} Ohm)"
            )

    # Quality assessment
    if len(diagnostics['warnings']) == 0:
        quality = 'good'
    elif len(diagnostics['warnings']) <= 2:
        quality = 'acceptable'
    else:
        quality = 'uncertain'

    logger.info(f"Analysis quality: {quality}")
    if diagnostics['warnings']:
        logger.info("Warnings:")
        for w in diagnostics['warnings']:
            logger.info(f"  - {w}")

    # Build return dict
    ratio = total_R_from_peaks / R_pol_data if R_pol_data > 0 else float('inf')

    return {
        'elements': elements,
        'quality': quality,
        'total_R': total_R_from_peaks,
        'R_pol': R_pol_data,
        'R_inf': R_inf,
        'ratio': ratio,
        'warnings': diagnostics['warnings'],
        'method': diagnostics['method']
    }


def format_voigt_report(voigt_info: dict) -> str:
    """
    Format Voigt analysis into readable report.

    Parameters
    ----------
    voigt_info : dict
        Result from analyze_voigt_elements()

    Returns
    -------
    report : str
        Formatted report for terminal
    """
    lines = []
    lines.append("=" * 60)
    lines.append("VOIGT ELEMENT ANALYSIS (R||C) FROM DRT")
    lines.append("=" * 60)

    # Detection method
    method = voigt_info['method'].upper()
    lines.append(f"Peak detection method: {method}")
    lines.append("")

    # Elements
    elements = voigt_info['elements']
    if len(elements) == 0:
        lines.append("No Voigt elements found")
        lines.append(f"Quality: {voigt_info['quality']}")
        lines.append("=" * 60)
        return "\n".join(lines)

    lines.append(f"Found {len(elements)} Voigt elements:")
    lines.append("")

    # Element table
    lines.append("  ID | tau [s]    | f [Hz]     | R [Ohm]   | C [F]      | Warnings")
    lines.append("  " + "-" * 72)

    for elem in elements:
        warnings_str = ", ".join(elem['warnings']) if elem['warnings'] else "-"
        if len(warnings_str) > 20:
            warnings_str = warnings_str[:17] + "..."

        line = (f"  {elem['id']:2d} | "
                f"{elem['tau']:10.2e} | "
                f"{elem['freq']:10.2e} | "
                f"{elem['R']:9.1f} | "
                f"{elem['C']:10.2e} | "
                f"{warnings_str}")
        lines.append(line)

    lines.append("")

    # Term classification (if available)
    has_classification = any(elem.get('classification') is not None for elem in elements)
    if has_classification:
        lines.append("Term type classification:")
        lines.append("")
        lines.append("  ID | Type     | Confidence   | Q n    | Reason")
        lines.append("  " + "-" * 72)

        for elem in elements:
            cls = elem.get('classification')
            if cls is not None:
                type_str = cls['type'].upper()
                conf_str = cls['confidence']
                n_str = f"{cls['cpe_n_estimate']:.3f}" if cls['cpe_n_estimate'] is not None else "-"
                reason = cls['reasoning']
                if len(reason) > 30:
                    reason = reason[:27] + "..."

                line = (f"  {elem['id']:2d} | "
                        f"{type_str:8s} | "
                        f"{conf_str:12s} | "
                        f"{n_str:6s} | "
                        f"{reason}")
                lines.append(line)
            else:
                lines.append(f"  {elem['id']:2d} | N/A      | -            | -      | (classification unavailable)")

        lines.append("")

    # R_pol validation
    lines.append("Consistency validation:")
    lines.append(f"  Sum R_i (from peaks): {voigt_info['total_R']:9.1f} Ohm")
    lines.append(f"  R_pol (from data):    {voigt_info['R_pol']:9.1f} Ohm")
    ratio_val = voigt_info['ratio']
    if ratio_val == float('inf'):
        lines.append("  Ratio:                INF (R_pol = 0)")
    else:
        lines.append(f"  Ratio:                {ratio_val:9.2f}")

    if ratio_val < 0.5 or ratio_val > 2.0:
        lines.append("  WARNING: Large difference between sum(R_i) and R_pol!")

    lines.append("")

    # Quality
    quality_map = {
        'good': 'GOOD',
        'acceptable': 'ACCEPTABLE',
        'uncertain': 'UNCERTAIN',
        'poor': 'POOR'
    }
    quality_en = quality_map.get(voigt_info['quality'], voigt_info['quality'])
    lines.append(f"Analysis quality: {quality_en}")

    # Global warnings
    if voigt_info['warnings']:
        lines.append("")
        lines.append("Warnings:")
        for warning in voigt_info['warnings']:
            lines.append(f"  - {warning}")

    lines.append("")

    # Recommendations for manual circuit building
    lines.append("Recommendations for manual circuit building:")
    lines.append("  1. Start with R_inf (series resistance):")
    lines.append("     R(R_inf)")

    # Decide based on classification whether to use C or Q
    if has_classification:
        lines.append("  2. Add elements according to classification:")
    else:
        lines.append("  2. Add Voigt elements (R||C) for each peak:")

    for i, elem in enumerate(elements, 1):
        cls = elem.get('classification')
        if cls is not None and cls['type'] == 'cpe':
            # Q element - n is not estimated, suggest R||C with comment
            lines.append(f"     Element {i}: (R(R{i}) | C(C{i}))  [{cls['type'].upper()} - try Q(Q{i}, n{i})]")
        else:
            # Standard Voigt (R||C)
            type_label = f" [{cls['type'].upper()}]" if cls is not None else ""
            lines.append(f"     Element {i}: (R(R{i}) | C(C{i})){type_label}")

    lines.append("  3. Connect elements in series with '-' operator:")

    # Example circuit (symbolic)
    example_parts = ["R(R_inf)"]
    for i in range(1, min(len(elements) + 1, 4)):  # Max 3 elements
        example_parts.append(f"(R(R{i}) | C(C{i}))")
    if len(elements) > 3:
        example_parts.append("...")
    example = " - ".join(example_parts)
    lines.append(f"     Example: {example}")

    lines.append("=" * 60)

    return "\n".join(lines)


__all__ = ['analyze_voigt_elements', 'format_voigt_report']

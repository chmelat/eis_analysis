"""
Classification of circuit element types from DRT peak characteristics.

This module analyzes the shape and properties of DRT peaks to classify
whether they correspond to ideal Voigt elements (R||C), Q-based elements
(R||Q), or other types.

Theoretical background:
- Voigt (R||C): Symmetric peak with characteristic width ~1.14 decades FWHM
- CPL/Q (R||Q, n<1): Broader, potentially asymmetric peak (distributed relaxation)
- Multiple processes: Multiple peaks close together
"""

import numpy as np
import logging
from typing import Dict, List, Optional, Literal
from numpy.typing import NDArray

logger = logging.getLogger(__name__)

# Classification thresholds (empirical, based on theoretical DRT)
VOIGT_WIDTH_THRESHOLD = 1.5  # decades in log10(tau) space
Q_WIDTH_THRESHOLD = 2.0       # decades - broader than ideal Voigt
ASYMMETRY_THRESHOLD = 1.5     # ratio (upper-center)/(center-lower)

# Confidence levels
Confidence = Literal['high', 'medium', 'low']


def classify_peak_type(
    peak: Dict,
    tau_full: NDArray[np.float64],
    gamma_full: NDArray[np.float64],
    adjacent_peaks: Optional[List[Dict]] = None
) -> Dict:
    """
    Klasifikuje typ obvodu odpovídající DRT píku.

    Analyzuje tvar píku v DRT spektru a rozhoduje, zda odpovídá:
    - 'voigt': Ideální Voigt element (R||C)
    - 'cpe': Q-based element (R||Q s n < 1)
    - 'multiple': Vícenásobné procesy blízko sebe
    - 'uncertain': Nelze spolehlivě určit

    Parameters
    ----------
    peak : dict
        GMM peak z gmm_peak_detection() obsahující:
        - 'tau_center': float [s]
        - 'tau_bounds': (lower, upper) [s]
        - 'log_tau_std': float (σ v log10 prostoru)
        - 'weight': float
        - 'R_estimate': float [Ω]
    tau_full : ndarray of float
        Kompletní τ osa z DRT
    gamma_full : ndarray of float
        Kompletní γ spektrum z DRT
    adjacent_peaks : list of dict, optional
        Seznam sousedních píků pro detekci vícenásobných procesů

    Returns
    -------
    classification : dict
        - 'type': str ('voigt', 'cpe', 'multiple', 'uncertain')
        - 'confidence': str ('high', 'medium', 'low')
        - 'cpe_n_estimate': float or None (1.0 pro Voigt, None jinak)
        - 'metrics': dict s analytickými metrikami
        - 'reasoning': str (vysvětlení klasifikace)
    """
    tau_center = peak['tau_center']
    tau_lower, tau_upper = peak['tau_bounds']
    log_tau_std = peak['log_tau_std']

    # === Metrika 1: Šířka píku v dekádách ===
    peak_width_decades = np.log10(tau_upper) - np.log10(tau_lower)
    log_tau_center = np.log10(tau_center)

    # === Metrika 2: Blízkost sousedních píků ===
    has_close_neighbors = False
    if adjacent_peaks:
        for adj_peak in adjacent_peaks:
            # Vzdálenost v log prostoru
            log_distance = abs(np.log10(adj_peak['tau_center']) - log_tau_center)
            if log_distance < 0.5:  # Méně než půl dekády
                has_close_neighbors = True
                break

    # === KLASIFIKAČNÍ LOGIKA ===

    metrics = {
        'peak_width_decades': peak_width_decades,
        'log_tau_std': log_tau_std,
        'has_close_neighbors': has_close_neighbors
    }

    # Priorita: detekce vícenásobných procesů
    if has_close_neighbors:
        return {
            'type': 'multiple',
            'confidence': 'medium',
            'cpe_n_estimate': None,
            'metrics': metrics,
            'reasoning': 'Vícenásobné píky blízko sebe (< 0.5 dekády) naznačují sériové procesy'
        }

    # Voigt element (R||C) - ideální kondenzátor
    # Očekávaná šířka: ~1.14 dekád (FWHM pro single RC)
    if peak_width_decades < VOIGT_WIDTH_THRESHOLD:
        confidence = 'high' if peak_width_decades < 1.3 else 'medium'
        return {
            'type': 'voigt',
            'confidence': confidence,
            'cpe_n_estimate': 1.0,  # Ideální kondenzátor
            'metrics': metrics,
            'reasoning': f'Úzký pík šířky {peak_width_decades:.2f} dekád odpovídá R||C'
        }

    # Q element (R||Q, n < 1) - distribuovaná relaxace
    # Širší pík kvůli distribuci relaxačních časů
    elif peak_width_decades > VOIGT_WIDTH_THRESHOLD:
        # POZNÁMKA: n < 0.6 (Warburg-like) se projevuje jako divergující ocas, ne pík
        # Pro wide peaky (Q s n typicky 0.7-0.9) ponecháváme určení n na circuit fitting

        if peak_width_decades > Q_WIDTH_THRESHOLD:
            confidence = 'high'
            reasoning = f'Široký pík ({peak_width_decades:.2f} dekád) naznačuje Q element - určete n z circuit fittingu'
        else:
            confidence = 'medium'
            reasoning = f'Mírně širší pík ({peak_width_decades:.2f} dekád) možná Q nebo rozložený R||C'

        return {
            'type': 'cpe',
            'confidence': confidence,
            'cpe_n_estimate': None,
            'metrics': metrics,
            'reasoning': reasoning
        }

    # Nejednoznačný případ
    else:
        return {
            'type': 'uncertain',
            'confidence': 'low',
            'cpe_n_estimate': None,
            'metrics': metrics,
            'reasoning': f'Nelze jednoznačně klasifikovat (šířka {peak_width_decades:.2f} dekád na hranici)'
        }


def classify_all_peaks(
    peaks_gmm: List[Dict],
    tau: NDArray[np.float64],
    gamma: NDArray[np.float64]
) -> List[Dict]:
    """
    Klasifikuje všechny GMM píky.

    Parameters
    ----------
    peaks_gmm : list of dict
        Seznam GMM píků z gmm_peak_detection()
    tau : ndarray of float
        Kompletní τ osa z DRT
    gamma : ndarray of float
        Kompletní γ spektrum z DRT

    Returns
    -------
    classifications : list of dict
        Seznam klasifikací pro každý pík (stejné pořadí jako peaks_gmm)
    """
    if not peaks_gmm:
        logger.warning("Žádné GMM píky k klasifikaci")
        return []

    logger.info("="*60)
    logger.info("Klasifikace typů termů z DRT píků")
    logger.info("="*60)

    classifications = []

    for i, peak in enumerate(peaks_gmm):
        # Pro detekci sousedů: všechny ostatní píky
        adjacent_peaks = [p for j, p in enumerate(peaks_gmm) if j != i]

        classification = classify_peak_type(peak, tau, gamma, adjacent_peaks)
        classifications.append(classification)

        # Log výsledek
        logger.info(f"\nPík {i+1} (τ = {peak['tau_center']:.2e} s):")
        logger.info(f"  Typ: {classification['type'].upper()}")
        logger.info(f"  Spolehlivost: {classification['confidence']}")
        if classification['cpe_n_estimate'] is not None:
            logger.info(f"  Q n estimate: {classification['cpe_n_estimate']:.3f}")
        logger.info(f"  Důvod: {classification['reasoning']}")
        logger.info("  Metriky:")
        logger.info(f"    - Šířka: {classification['metrics']['peak_width_decades']:.2f} dekád")
        logger.info(f"    - σ(log τ): {classification['metrics']['log_tau_std']:.3f}")

    # Shrnutí
    type_counts = {}
    for cls in classifications:
        type_counts[cls['type']] = type_counts.get(cls['type'], 0) + 1

    logger.info("\n" + "="*60)
    logger.info("Shrnutí klasifikace:")
    for term_type, count in sorted(type_counts.items()):
        logger.info(f"  {term_type.upper()}: {count} píků")
    logger.info("="*60)

    return classifications

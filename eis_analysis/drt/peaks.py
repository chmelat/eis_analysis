"""
Peak detection methods for DRT spectra.
"""

import numpy as np
import logging
from typing import Tuple, List, Dict, Optional
from numpy.typing import NDArray

logger = logging.getLogger(__name__)

# GMM for robust peak detection
try:
    from sklearn.mixture import GaussianMixture
    GMM_AVAILABLE = True
except ImportError:
    GMM_AVAILABLE = False


def gmm_peak_detection(
    tau: NDArray[np.float64],
    gamma: NDArray[np.float64],
    n_components_range: Tuple[int, int] = (1, 6),
    bic_threshold: float = 10.0
) -> Tuple[List[Dict], Optional[object], List[float]]:
    """
    Detekuje píky v DRT spektru pomocí Gaussian Mixture Models.

    Poskytuje explicitní hranice píků pro přesný výpočet R_i v suggest_circuit_from_drt.
    GMM automaticky dekonvolvuje překrývající se píky a dává pravděpodobnostní rámec.

    Model selection strategie:
    Používá BIC (Bayesian Information Criterion) s early stopping. Začíná s nejmenším
    počtem komponent a přidává další, dokud BIC zlepšení překročí bic_threshold.
    Toto implementuje Occamovu břitvu - preferuje jednodušší modely, pokud komplexnější
    model nepřináší signifikantní zlepšení.

    Parametry:
    - tau: časové konstanty [s]
    - gamma: distribuční funkce [Ω]
    - n_components_range: rozsah pro hledání optimálního počtu píků (default: 1-6)
    - bic_threshold: minimální BIC zlepšení pro přidání další komponenty (default: 10.0)
                     Vyšší hodnota = konzervativnější (méně píků)
                     Typické hodnoty: 2-20 (10 je rozumný kompromis)

    Vrací:
    - peaks: list of dict obsahující informace o pících
    - gmm_model: natrénovaný GaussianMixture model (nebo None při selhání)
    - bic_scores: BIC skóre pro různé počty komponent
    """
    if not GMM_AVAILABLE:
        logger.error("sklearn není dostupný - GMM detekce není možná")
        return [], None, []

    logger.info("="*60)
    logger.info("GMM detekce píků")
    logger.info("="*60)

    # Log transformace pro lepší separaci
    log_tau = np.log10(tau)

    # Pouze kladné hodnoty
    mask_pos = gamma > 0
    if not np.any(mask_pos):
        logger.warning("Všechny hodnoty γ jsou ≤0, nemohu fitovat GMM")
        return [], None, []

    X = log_tau[mask_pos].reshape(-1, 1)
    weights = gamma[mask_pos]

    # Replikace bodů podle váhy
    weights_normalized = weights / np.mean(weights)
    weights_int = np.maximum(1, np.round(weights_normalized).astype(int))
    X_weighted = np.repeat(X, weights_int, axis=0)

    logger.info(f"Hledám optimální počet komponent v rozsahu {n_components_range}")
    logger.debug(f"Data: {len(X)} bodů → {len(X_weighted)} vážených bodů")

    # Model selection pomocí BIC
    bic_scores = []
    models = []

    for n in range(n_components_range[0], n_components_range[1] + 1):
        try:
            gmm = GaussianMixture(
                n_components=n,
                covariance_type='full',
                random_state=42,
                max_iter=200,
                tol=1e-4
            )
            gmm.fit(X_weighted)
            # BIC počítat na stejných datech jako fit (X_weighted)
            bic = gmm.bic(X_weighted)
            bic_scores.append(bic)
            models.append(gmm)
            logger.info(f"n={n} komponenty: BIC={bic:.2f}")
        except Exception as e:
            logger.warning(f"GMM fit selhal pro n={n}: {e}")
            bic_scores.append(np.inf)
            models.append(None)

    if not models or all(m is None for m in models):
        logger.error("Všechny GMM fity selhaly")
        return [], None, bic_scores

    # Vyber model pomocí BIC s early stopping (bic_threshold)
    # Strategie: Přidávej komponenty dokud BIC zlepšení > threshold
    valid_bic = [(i, bic) for i, bic in enumerate(bic_scores) if bic != np.inf]

    if not valid_bic:
        logger.error("Všechny BIC hodnoty jsou neplatné")
        return [], None, bic_scores

    # Seřaď podle indexu (počtu komponent)
    valid_bic_sorted = sorted(valid_bic, key=lambda x: x[0])

    # Early stopping: začni s nejmenším n, přidávej dokud zlepšení > threshold
    best_idx = valid_bic_sorted[0][0]  # Začni s prvním validním modelem

    logger.info(f"BIC threshold pro přidání komponenty: {bic_threshold:.1f}")

    for i in range(1, len(valid_bic_sorted)):
        prev_idx = valid_bic_sorted[i-1][0]
        curr_idx = valid_bic_sorted[i][0]
        prev_bic = bic_scores[prev_idx]
        curr_bic = bic_scores[curr_idx]

        improvement = prev_bic - curr_bic
        n = curr_idx + n_components_range[0]
        logger.info(f"n={n}: BIC zlepšení {improvement:.1f}")

        if improvement > bic_threshold:
            # Signifikantní zlepšení - přijmi tento model
            best_idx = curr_idx
        else:
            # Zlepšení pod prahem - zastav (Occamova břitva)
            logger.info(f"  → Zlepšení {improvement:.1f} < threshold {bic_threshold:.1f}, zastavuji")
            break

    best_n = best_idx + n_components_range[0]
    best_gmm = models[best_idx]
    n_peaks = best_n

    # Porovnání s absolutním minimem BIC
    abs_min_idx = min(valid_bic, key=lambda x: x[1])[0]
    abs_min_n = abs_min_idx + n_components_range[0]

    if abs_min_idx != best_idx:
        logger.info(f"Early stopping: n={n_peaks} (absolutní BIC minimum by bylo n={abs_min_n})")
        total_improvement = bic_scores[0] - bic_scores[best_idx]
        logger.info(f"BIC selection: n={n_peaks} komponent (celkové BIC zlepšení: {total_improvement:.1f})")
    elif best_idx > 0:
        total_improvement = bic_scores[0] - bic_scores[best_idx]
        logger.info(f"BIC selection: n={n_peaks} komponent (BIC zlepšení: {total_improvement:.1f})")

    # Varování pokud optimum je na okraji rozsahu
    if n_peaks == n_components_range[0]:
        logger.warning(f"⚠ Optimum na spodní hranici rozsahu (n={n_peaks})")
        logger.warning("   Zvažte rozšíření rozsahu nebo kontrolu dat")
    elif n_peaks == n_components_range[1]:
        logger.warning(f"⚠ Optimum na horní hranici rozsahu (n={n_peaks})")
        logger.warning("   Zvažte rozšíření rozsahu směrem nahoru")

    if best_gmm is None:
        logger.error(f"Vybraný model (n={n_peaks}) selhal")
        return [], None, bic_scores

    logger.info(f"✓ Optimální počet píků: {n_peaks} (BIC={bic_scores[best_idx]:.2f})")

    # Extrahuj parametry píků
    peaks = []
    ln_tau = np.log(tau)

    for i in range(best_gmm.n_components):
        mu = best_gmm.means_[i, 0]  # v log10(tau)
        sigma = np.sqrt(best_gmm.covariances_[i, 0, 0])
        weight = best_gmm.weights_[i]

        # Převod zpět na tau (lineární prostor)
        tau_center = 10**mu
        tau_lower = 10**(mu - 2*sigma)  # 95% confidence interval
        tau_upper = 10**(mu + 2*sigma)

        f_center = 1 / (2 * np.pi * tau_center)

        # Odhad R_i přímou integrací gamma v oblasti píku (tau_bounds)
        # Najdi indexy odpovídající tau_bounds
        idx_lower = np.searchsorted(tau, tau_lower)
        idx_upper = np.searchsorted(tau, tau_upper)

        # Ošetření hranic
        idx_lower = max(0, idx_lower)
        idx_upper = min(len(tau) - 1, idx_upper)

        # Integrál gamma v této oblasti
        if idx_upper > idx_lower:
            from ..utils.compat import np_trapz
            R_estimate = np_trapz(gamma[idx_lower:idx_upper+1], ln_tau[idx_lower:idx_upper+1])
        else:
            R_estimate = 0.0

        peaks.append({
            'tau_center': tau_center,
            'tau_bounds': (tau_lower, tau_upper),
            'log_tau_std': sigma,
            'weight': weight,
            'R_estimate': R_estimate,
            'f_center': f_center
        })

    # Seřaď píky podle tau_center
    peaks = sorted(peaks, key=lambda p: p['tau_center'])

    logger.info("\nDetekované píky (seřazeno podle τ):")
    for i, peak in enumerate(peaks):
        logger.info(f"  Pík {i+1}:")
        logger.info(f"    τ = {peak['tau_center']:.2e} s (f = {peak['f_center']:.2e} Hz)")
        logger.info(f"    Hranice τ: [{peak['tau_bounds'][0]:.2e}, {peak['tau_bounds'][1]:.2e}] s")
        logger.info(f"    Šířka (σ): {peak['log_tau_std']:.3f} dekád")
        logger.info(f"    Váha: {peak['weight']:.3f}")
        logger.info(f"    R ~ {peak['R_estimate']:.2f} Ω")

    return peaks, best_gmm, bic_scores

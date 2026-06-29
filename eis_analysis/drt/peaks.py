"""
Peak detection methods for DRT spectra.
"""

import numpy as np
import logging
from dataclasses import dataclass
from typing import Tuple, List, Dict, Optional
from numpy.typing import NDArray
from scipy.special import logsumexp

logger = logging.getLogger(__name__)

# Regularizace kovariance (přičteno k rozptylu) — brání singularitě, když se
# komponenta smrští na jediný bod. Stejná hodnota jako default sklearn reg_covar.
_REG_COVAR = 1e-6


@dataclass
class WeightedGMMResult:
    """Výsledek váženého EM fitu 1D Gaussovské směsi.

    Atributy zrcadlí rozhraní sklearn ``GaussianMixture`` (tvary means_,
    covariances_, weights_), aby navazující kód nemusel rozlišovat původ modelu.
    """
    n_components: int
    means_: NDArray[np.float64]        # (K, 1) — v log10(tau)
    covariances_: NDArray[np.float64]  # (K, 1, 1)
    weights_: NDArray[np.float64]      # (K,) — mixing weights, Σ = 1


def _weighted_gaussian_mixture_1d(
    x: NDArray[np.float64],
    w: NDArray[np.float64],
    n_components: int,
    max_iter: int = 200,
    tol: float = 1e-4,
    var_floor: float = 1e-12,
) -> Tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], float]:
    """Fituje 1D Gaussovskou směs na vážená data váženým EM.

    DRT γ(τ) je vážená hustota, ne vzorek. Místo celočíselné replikace bodů
    (metodicky nekorektní, viz audit F1) vstupují váhy ``w`` přímo do EM jako
    frakční počty. Odhady μ, σ², π jsou na škále vah invariantní; absolutní škála
    ``w`` ovlivní jen velikost log-likelihood (zde využito pro BIC s N = počet
    skutečných měření).

    Parametry
    ---------
    x : (M,) souřadnice bodů (log10 τ)
    w : (M,) nezáporné váhy (γ)
    n_components : počet komponent K
    var_floor : dolní mez rozptylu (kvadrát rozlišení mřížky), brání kolapsu

    Vrací (means (K,), variances (K,), mix_weights (K,), weighted_logL).
    """
    K = n_components
    w_total = w.sum()

    # Deterministický init: μ na rovnoměrně rozmístěných vážených kvantilech x.
    # Nahrazuje sklearn random_state=42 — pro 1D je kvantilový init stabilní.
    order = np.argsort(x)
    x_sorted = x[order]
    cdf = np.cumsum(w[order]) / w_total
    quantile_pts = (np.arange(K) + 0.5) / K
    means = np.interp(quantile_pts, cdf, x_sorted)

    # Vážený celkový rozptyl / K jako počáteční šířka komponent.
    x_mean = float(np.sum(w * x) / w_total)
    var_total = float(np.sum(w * (x - x_mean) ** 2) / w_total)
    variances = np.full(K, max(var_total / K, var_floor))
    mix = np.full(K, 1.0 / K)

    prev_logL = -np.inf
    weighted_logL = -np.inf
    for _ in range(max_iter):
        # E-krok: log-responsibility přes logsumexp (numericky stabilní).
        # log N(x|μ,σ²) = -0.5*ln(2πσ²) - (x-μ)²/(2σ²)
        log_norm = -0.5 * np.log(2 * np.pi * variances)            # (K,)
        sq = (x[:, None] - means[None, :]) ** 2                    # (M, K)
        log_prob = log_norm[None, :] - sq / (2 * variances[None, :])
        log_weighted = log_prob + np.log(mix)[None, :]             # (M, K)
        log_denom = logsumexp(log_weighted, axis=1)                # (M,)
        weighted_logL = float(np.sum(w * log_denom))
        resp = np.exp(log_weighted - log_denom[:, None])           # (M, K)

        # M-krok (vážený)
        wr = w[:, None] * resp                                     # (M, K)
        Nk = wr.sum(axis=0)                                        # (K,)
        Nk = np.maximum(Nk, 1e-300)  # ochrana proti dělení nulou u prázdné komponenty
        mix = Nk / w_total
        means = (wr * x[:, None]).sum(axis=0) / Nk
        diff_sq = (x[:, None] - means[None, :]) ** 2
        variances = (wr * diff_sq).sum(axis=0) / Nk + _REG_COVAR
        variances = np.maximum(variances, var_floor)

        if abs(weighted_logL - prev_logL) < tol:
            break
        prev_logL = weighted_logL

    return means, variances, mix, weighted_logL


def gmm_peak_detection(
    tau: NDArray[np.float64],
    gamma: NDArray[np.float64],
    n_components_range: Tuple[int, int] = (1, 6),
    bic_threshold: float = 10.0,
    n_data: Optional[int] = None
) -> Tuple[List[Dict], Optional[WeightedGMMResult], List[float]]:
    """
    Detekuje píky v DRT spektru pomocí vážené Gaussovské směsi (vážený EM).

    Poskytuje explicitní hranice píků pro přesný výpočet R_i v suggest_circuit_from_drt.
    Směs automaticky dekonvolvuje překrývající se píky a dává pravděpodobnostní rámec.

    DRT γ(τ) je vážená hustota — váhou bodu je přímo γ. Směs se fituje váženým EM
    (žádná replikace bodů, viz audit F1). Pro model selection se používá BIC, jehož
    N = počet skutečných měření (frekvencí), takže ``bic_threshold`` je srovnatelný
    napříč datasety.

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
    - n_data: počet skutečných měření (frekvencí) pro penalizaci BIC. Pokud None,
              použije se počet binů s γ>0 (fallback pro přímé volání).

    Vrací:
    - peaks: list of dict obsahující informace o pících
    - gmm_model: WeightedGMMResult (nebo None při selhání)
    - bic_scores: BIC skóre pro různé počty komponent
    """
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

    x = log_tau[mask_pos]
    w = gamma[mask_pos].astype(np.float64)

    if n_data is None:
        n_data = len(x)

    # Normalizace vah na Σw = n_data: log-likelihood pak má škálu n_data měření
    # a ladí s penaltou k·ln(n_data) v BIC (jádro nápravy F1).
    w = w * (n_data / w.sum())

    # Podlaha rozptylu = (polovina mediánu rozestupu mřížky)² — brání tomu, aby se
    # komponenta smrskla na jediný bin mřížky (singularita / nereálně úzký pík).
    if len(x) > 1:
        spacing = np.median(np.diff(np.sort(x)))
        var_floor = (0.5 * spacing) ** 2 if spacing > 0 else 1e-12
    else:
        var_floor = 1e-12

    logger.info(f"Hledám optimální počet komponent v rozsahu {n_components_range}")
    logger.debug(f"Data: {len(x)} binů γ>0, N (měření) = {n_data}")

    # Model selection pomocí BIC
    bic_scores = []
    models: List[Optional[WeightedGMMResult]] = []

    for n in range(n_components_range[0], n_components_range[1] + 1):
        try:
            means, variances, mix, logL = _weighted_gaussian_mixture_1d(
                x, w, n_components=n, var_floor=var_floor
            )
            # Počet volných parametrů 1D směsi: μ(K) + σ²(K) + π(K-1) = 3K-1.
            k_params = 3 * n - 1
            bic = -2.0 * logL + k_params * np.log(n_data)
            model = WeightedGMMResult(
                n_components=n,
                means_=means.reshape(-1, 1),
                covariances_=variances.reshape(-1, 1, 1),
                weights_=mix,
            )
            bic_scores.append(bic)
            models.append(model)
            logger.info(f"n={n} komponenty: BIC={bic:.2f}")
        except (ValueError, np.linalg.LinAlgError) as e:
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

    # Celkový R_pol = Σγ·Δlnτ (obdélníkové pravidlo, konzistentní s jádrem DRT,
    # viz audit F10). Dílčí odpory rozdělíme podle vah GMM (Σ weight = 1), takže
    # Σ R_i = R_pol přesně. Tím se vyhneme double-countingu, který vzniká
    # integrací CELKOVÉHO γ přes ±2σ okno každého píku.
    d_ln_tau = float(np.mean(np.diff(np.log(tau))))
    R_pol = float(np.sum(gamma) * d_ln_tau)

    for i in range(best_gmm.n_components):
        mu = best_gmm.means_[i, 0]  # v log10(tau)
        sigma = np.sqrt(best_gmm.covariances_[i, 0, 0])
        weight = best_gmm.weights_[i]

        # Převod zpět na tau (lineární prostor)
        tau_center = 10**mu
        tau_lower = 10**(mu - 2*sigma)  # 95% confidence interval
        tau_upper = 10**(mu + 2*sigma)

        f_center = 1 / (2 * np.pi * tau_center)

        # Odhad R_i rozkladem R_pol podle váhy komponenty (Σ weight = 1)
        R_estimate = float(weight * R_pol)

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

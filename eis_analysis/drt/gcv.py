"""
Regularization parameter selection for DRT analysis.

Methods:
- GCV (Generalized Cross-Validation) - fast initial estimate
- L-curve - robust correction for NNLS problems
- Hybrid (GCV + L-curve) - recommended for DRT with NNLS
"""

import numpy as np
import logging
from typing import Tuple, Optional
from numpy.typing import NDArray
from scipy.optimize import nnls

logger = logging.getLogger(__name__)


def compute_gcv_score(lambda_val: float, A: NDArray[np.float64],
                      b: NDArray[np.float64], L: NDArray[np.float64]) -> float:
    """
    Vypočítá GCV (Generalized Cross Validation) score pro danou λ.

    GCV(λ) = n * ||b - A·x(λ)||² / (trace(I - K(λ)))²

    kde K(λ) = A @ inv(A^T·A + λ·L^T·L) @ A^T

    Menší GCV → lepší volba λ

    Parametry:
    - lambda_val: testovaná hodnota regularizačního parametru
    - A: matice systému [2N × N_tau]
    - b: pravá strana [2N]
    - L: regularizační matice [(N_tau-2) × N_tau] pro druhou derivaci

    Vrací:
    - float: GCV score (menší je lepší)

    Reference:
    - Wahba (1985), Annals of Statistics
    - Saccoccio et al. (2014), Electrochimica Acta 147, 470-482
    """
    n = len(b)

    # Regularizovaná matice pro NNLS
    L_scaled = np.sqrt(lambda_val) * L
    A_reg = np.vstack([A, L_scaled])
    b_reg = np.concatenate([b, np.zeros(L.shape[0])])

    # Řeš NNLS
    try:
        x, residual_nnls = nnls(A_reg, b_reg)
    except Exception as e:
        logger.debug(f"NNLS selhalo pro λ={lambda_val:.2e}: {e}")
        return np.inf

    # Reziduum původního problému (bez regularizace)
    residual = b - A @ x
    residual_norm_sq = np.sum(residual**2)

    # Výpočet trace(I - K) pomocí přímé inverze
    # K = A @ inv(A^T @ A + λ·L^T·L) @ A^T
    # trace(I - K) = n - trace(K)
    #
    # POZNÁMKA: Pro NNLS je toto pouze aproximace, protože NNLS má nelineární
    # constraint (γ ≥ 0). GCV je striktně definováno pro lineární LSQ.
    # Pro přesnější výběr λ u NNLS zvažte L-curve metodu.

    try:
        # Přímý výpočet trace(K) - matematicky korektní pro LSQ
        M = A.T @ A + lambda_val * L.T @ L

        # K = A @ inv(M) @ A^T, ale potřebujeme pouze trace(K)
        # trace(K) = trace(A^T @ A @ inv(M))
        # Pro symetrickou pozitivně definitní M použijeme solve
        K = A @ np.linalg.solve(M, A.T)
        trace_K = np.trace(K)

        trace_I_minus_K = n - trace_K

        if trace_I_minus_K <= 1e-10:
            logger.debug(f"Neplatný trace(I-K)={trace_I_minus_K:.2e} pro λ={lambda_val:.2e}")
            return np.inf

        # GCV score (s normalizací n)
        gcv = n * residual_norm_sq / trace_I_minus_K**2

        return gcv

    except (np.linalg.LinAlgError, ValueError) as e:
        logger.debug(f"Chyba při výpočtu GCV pro λ={lambda_val:.2e}: {e}")
        return np.inf


def compute_lcurve_point(lambda_val: float, A: NDArray[np.float64],
                          b: NDArray[np.float64], L: NDArray[np.float64]
                          ) -> Tuple[float, float, Optional[NDArray]]:
    """
    Vypočítá bod na L-křivce pro danou λ.

    L-křivka je graf: log||Ax - b|| vs log||Lx||
    (reziduum vs regularizační norma)

    Parametry:
    - lambda_val: regularizační parametr
    - A: matice systému
    - b: pravá strana
    - L: regularizační matice

    Vrací:
    - tuple: (log_residual, log_reg_norm, x) nebo (inf, inf, None) při selhání
    """
    # Regularizovaná matice pro NNLS
    L_scaled = np.sqrt(lambda_val) * L
    A_reg = np.vstack([A, L_scaled])
    b_reg = np.concatenate([b, np.zeros(L.shape[0])])

    try:
        x, _ = nnls(A_reg, b_reg)
    except Exception as e:
        logger.debug(f"NNLS selhalo pro λ={lambda_val:.2e}: {e}")
        return np.inf, np.inf, None

    # Norma rezidua ||Ax - b||
    residual = b - A @ x
    residual_norm = np.linalg.norm(residual)

    # Regularizační norma ||Lx||
    reg_norm = np.linalg.norm(L @ x)

    # Log transformace (s ochranou proti log(0))
    log_res = np.log10(max(residual_norm, 1e-15))
    log_reg = np.log10(max(reg_norm, 1e-15))

    return log_res, log_reg, x


def compute_lcurve_curvature(rho: NDArray[np.float64],
                              eta: NDArray[np.float64]) -> NDArray[np.float64]:
    """
    Vypočítá křivost L-křivky v každém bodě.

    Křivost v log-log prostoru:
    κ = (ρ'η'' - ρ''η') / (ρ'² + η'²)^(3/2)

    kde ρ = log||Ax-b||, η = log||Lx||

    Používá np.gradient() pro numerickou derivaci (vektorizováno).

    Parametry:
    - rho: log||Ax - b|| pro různé λ
    - eta: log||Lx|| pro různé λ

    Vrací:
    - curvature: křivost v každém bodě
    """
    n = len(rho)
    if n < 3:
        return np.zeros(n)

    # První derivace pomocí np.gradient (centrální diference)
    d_rho = np.gradient(rho)
    d_eta = np.gradient(eta)

    # Druhá derivace
    dd_rho = np.gradient(d_rho)
    dd_eta = np.gradient(d_eta)

    # Křivost (vektorizovaně)
    numerator = d_rho * dd_eta - dd_rho * d_eta
    denominator = (d_rho**2 + d_eta**2)**1.5

    # Podmíněné přiřazení pomocí np.where
    curvature = np.where(np.abs(denominator) > 1e-15, numerator / denominator, 0.0)

    return curvature


def find_lcurve_corner(lambda_values: NDArray[np.float64],
                        rho: NDArray[np.float64],
                        eta: NDArray[np.float64]) -> Tuple[float, int, NDArray]:
    """
    Najde roh L-křivky (bod maximální křivosti).

    Parametry:
    - lambda_values: testované hodnoty λ
    - rho: log||Ax - b|| pro každé λ
    - eta: log||Lx|| pro každé λ

    Vrací:
    - tuple: (lambda_corner, corner_index, curvature_array)
    """
    curvature = compute_lcurve_curvature(rho, eta)

    # Najdi maximum křivosti (ignoruj krajní body)
    # Hledáme kladnou křivost (konvexní část L-křivky)
    valid_range = slice(1, len(curvature) - 1)
    curvature_inner = curvature[valid_range]

    if len(curvature_inner) == 0:
        # Fallback: prostřední bod
        corner_idx = len(lambda_values) // 2
    else:
        # Maximum kladné křivosti
        corner_idx_inner = np.argmax(curvature_inner)
        corner_idx = corner_idx_inner + 1  # Kompenzace za slice

    lambda_corner = lambda_values[corner_idx]

    return lambda_corner, corner_idx, curvature


def find_optimal_lambda_gcv(A: NDArray[np.float64], b: NDArray[np.float64],
                            L: NDArray[np.float64],
                            lambda_range: Tuple[float, float] = (1e-5, 1.0),
                            n_search: int = 20) -> Tuple[float, float]:
    """
    Najde optimální lambda minimalizací GCV score.

    Strategie: Nejprve hrubé prohledání, pak jemné doladění.

    Parametry:
    - A: matice systému [2N × N_tau]
    - b: pravá strana [2N]
    - L: regularizační matice
    - lambda_range: tuple (min, max) pro prohledávání λ
    - n_search: počet bodů pro hrubé prohledání

    Vrací:
    - tuple: (lambda_optimal, gcv_optimal)
    """
    logger.info("Automatický výběr λ pomocí GCV...")

    # Fáze 1: Hrubé prohledání v log-prostoru
    lambda_values = np.logspace(np.log10(lambda_range[0]),
                                 np.log10(lambda_range[1]),
                                 n_search)

    gcv_scores = []
    for lam in lambda_values:
        score = compute_gcv_score(lam, A, b, L)
        gcv_scores.append(score)

    gcv_scores = np.array(gcv_scores)

    # Najdi minimum
    min_idx = np.argmin(gcv_scores)
    lambda_coarse = lambda_values[min_idx]
    gcv_coarse = gcv_scores[min_idx]

    logger.debug(f"Hrubé prohledání: λ = {lambda_coarse:.4e}, GCV = {gcv_coarse:.4e}")

    # Fáze 2: Jemné doladění kolem minima
    if min_idx == 0:
        # Minimum na okraji - rozšiř prohledávání dolů
        fine_range = (lambda_range[0] / 10, lambda_values[min(min_idx + 2, len(lambda_values)-1)])
    elif min_idx == len(lambda_values) - 1:
        # Minimum na okraji - rozšiř prohledávání nahoru
        fine_range = (lambda_values[max(min_idx - 2, 0)], lambda_range[1] * 10)
    else:
        # Minimum uvnitř - jemné prohledání v okolí
        fine_range = (lambda_values[max(min_idx - 1, 0)],
                      lambda_values[min(min_idx + 1, len(lambda_values)-1)])

    # Jemné prohledání
    lambda_fine = np.logspace(np.log10(fine_range[0]),
                              np.log10(fine_range[1]),
                              n_search)

    gcv_fine = []
    for lam in lambda_fine:
        score = compute_gcv_score(lam, A, b, L)
        gcv_fine.append(score)

    gcv_fine = np.array(gcv_fine)
    min_fine_idx = np.argmin(gcv_fine)
    lambda_optimal = lambda_fine[min_fine_idx]
    gcv_optimal = gcv_fine[min_fine_idx]

    logger.info(f"✓ Optimální λ = {lambda_optimal:.4e} (GCV = {gcv_optimal:.4e})")

    return lambda_optimal, gcv_optimal


def find_optimal_lambda_hybrid(A: NDArray[np.float64], b: NDArray[np.float64],
                                L: NDArray[np.float64],
                                lambda_range: Tuple[float, float] = (1e-5, 1.0),
                                n_search: int = 20,
                                lcurve_decades: float = 1.5
                                ) -> Tuple[float, float, dict]:
    """
    Hybridní metoda pro výběr λ: GCV jako initial guess, L-curve jako korektor.

    Strategie:
    1. GCV pro rychlý hrubý odhad λ_gcv
    2. L-curve prohledání v okolí λ_gcv (±lcurve_decades dekád)
    3. Výběr λ z rohu L-křivky

    Toto řeší problém GCV pro NNLS: GCV předpokládá lineární řešení,
    ale NNLS má nelineární constraint (x ≥ 0). L-curve je robustnější
    pro NNLS problémy, ale pomalejší pro celý rozsah.

    Parametry:
    - A: matice systému [2N × N_tau]
    - b: pravá strana [2N]
    - L: regularizační matice
    - lambda_range: tuple (min, max) pro počáteční GCV prohledávání
    - n_search: počet bodů pro prohledávání
    - lcurve_decades: rozsah L-curve kolem GCV odhadu v dekádách (default: 1.5)

    Vrací:
    - tuple: (lambda_optimal, score, diagnostics)
      - diagnostics obsahuje: lambda_gcv, lambda_lcurve, method_used, curvature, rho, eta
    """
    logger.info("Hybridní výběr λ (GCV + L-curve korekce)...")

    diagnostics = {
        'lambda_gcv': None,
        'lambda_lcurve': None,
        'method_used': None,
        'curvature': None,
        'rho': None,
        'eta': None,
        'lambda_values': None
    }

    # === Fáze 1: GCV pro initial guess ===
    logger.info("Fáze 1: GCV initial guess...")

    lambda_values = np.logspace(np.log10(lambda_range[0]),
                                 np.log10(lambda_range[1]),
                                 n_search)

    gcv_scores = []
    for lam in lambda_values:
        score = compute_gcv_score(lam, A, b, L)
        gcv_scores.append(score)

    gcv_scores = np.array(gcv_scores)
    gcv_min_idx = np.argmin(gcv_scores)
    lambda_gcv = lambda_values[gcv_min_idx]

    diagnostics['lambda_gcv'] = lambda_gcv
    logger.info(f"  GCV minimum: λ_gcv = {lambda_gcv:.4e}")

    # === Fáze 2: L-curve v okolí GCV odhadu ===
    logger.info("Fáze 2: L-curve korekce...")

    # Rozsah pro L-curve: ±lcurve_decades dekád kolem λ_gcv
    lcurve_min = max(lambda_gcv / (10**lcurve_decades), lambda_range[0] / 10)
    lcurve_max = min(lambda_gcv * (10**lcurve_decades), lambda_range[1] * 10)

    lambda_lcurve_range = np.logspace(np.log10(lcurve_min),
                                       np.log10(lcurve_max),
                                       n_search)

    # Výpočet L-curve bodů
    rho = []  # log||Ax - b||
    eta = []  # log||Lx||

    for lam in lambda_lcurve_range:
        log_res, log_reg, _ = compute_lcurve_point(lam, A, b, L)
        rho.append(log_res)
        eta.append(log_reg)

    rho = np.array(rho)
    eta = np.array(eta)

    diagnostics['rho'] = rho
    diagnostics['eta'] = eta
    diagnostics['lambda_values'] = lambda_lcurve_range

    # Najdi roh L-křivky
    lambda_lcurve, corner_idx, curvature = find_lcurve_corner(lambda_lcurve_range, rho, eta)

    diagnostics['lambda_lcurve'] = lambda_lcurve
    diagnostics['curvature'] = curvature

    logger.info(f"  L-curve roh: λ_lcurve = {lambda_lcurve:.4e}")

    # === Fáze 3: Rozhodnutí ===
    # Porovnej GCV a L-curve výsledky
    ratio = lambda_lcurve / lambda_gcv

    if 0.1 < ratio < 10:
        # Výsledky jsou konzistentní (v rámci 1 dekády) - použij L-curve
        lambda_optimal = lambda_lcurve
        diagnostics['method_used'] = 'lcurve'
        logger.info(f"  Konzistentní výsledky (ratio={ratio:.2f}), použita L-curve")
    elif ratio >= 10:
        # L-curve chce výrazně vyšší λ - pravděpodobně NNLS efekt
        # L-curve je spolehlivější pro NNLS
        lambda_optimal = lambda_lcurve
        diagnostics['method_used'] = 'lcurve_correction'
        logger.warning(f"  L-curve korekce: λ_lcurve >> λ_gcv (ratio={ratio:.1f})")
        logger.warning("  GCV pravděpodobně underestimuje λ kvůli NNLS constraint")
    else:
        # L-curve chce výrazně nižší λ - neobvyklé, buď opatrný
        # Použij geometrický průměr jako kompromis
        lambda_optimal = np.sqrt(lambda_gcv * lambda_lcurve)
        diagnostics['method_used'] = 'geometric_mean'
        logger.warning(f"  Nekonzistentní výsledky (ratio={ratio:.2f})")
        logger.warning(f"  Použit geometrický průměr: λ = {lambda_optimal:.4e}")

    # Spočítej finální GCV score pro reporting
    final_gcv = compute_gcv_score(lambda_optimal, A, b, L)

    logger.info(f"✓ Hybridní optimum: λ = {lambda_optimal:.4e} (GCV = {final_gcv:.4e})")

    return lambda_optimal, final_gcv, diagnostics

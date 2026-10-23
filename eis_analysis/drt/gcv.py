"""
Regularization parameter selection for DRT analysis.

Methods:
- GCV (Generalized Cross-Validation) - fast initial estimate
- L-curve - corner of the residual / roughness trade-off
- Hybrid (GCV + L-curve, the larger lambda wins) - recommended for DRT with NNLS

CONTRACT

Problem. For a lambda > 0, x(lambda) is the NNLS solution (x >= 0) of
min ||A x - b||^2 + lambda ||L x||^2, solved as the stacked system
[A; sqrt(lambda) L] x = [b; 0] with DRT_NNLS_MAXITER_FACTOR x (columns)
iterations. An NNLS that raises (no convergence, NaN or inf in the system)
is a failed point: it gives no number, not a partial one.

Callers. calculate_drt through linear_system._select_lambda calls
find_optimal_lambda_hybrid(A, b, L, n_search=20, lcurve_decades=1.5) with
the default lambda_range = DRT_LAMBDA_RANGE. A is (2N, k), the real rows
of the N points, then the imaginary ones, weighted; b (2N,) is
weighted Re(Z) - R_inf, then Im(Z); L is (k - 2, k), the second difference
in ln tau scaled by sqrt(2N) / d_ln_tau^1.5, with a zero last column when
the model has a series inductance. All finite (calculate_drt rejects NaN
and inf in f and Z). k is n_tau (default 100, set by --n-tau) plus the tau
extension, so k > 2N is normal. It reads lambda_gcv, lambda_lcurve,
method_used, corner_at_edge and corner_below_gcv from the diagnostics,
and sets lambda_at_edge itself where the returned lambda or lambda_gcv is
outside the open DRT_LAMBDA_RANGE. On ValueError or LinAlgError it falls
back to find_optimal_lambda_gcv, on a second one to DRT_LAMBDA_DEFAULT;
any other exception fails the DRT. linear_system and drt/significance.py
use DRT_NNLS_MAXITER_FACTOR for their own NNLS. Public API
(eis_analysis.drt; compute_gcv_score and find_optimal_lambda_gcv also
eis_analysis): any caller, nothing guaranteed beyond the rules below.

compute_gcv_score(lambda_val, A, b, L) -> float.
  GCV = m ||b - A x||^2 / (m - trace K)^2, m = len(b), x = x(lambda), and
  K = A (A^T A + lambda L^T L)^-1 A^T from the unconstrained problem over
  all columns (deliberately, see the comment in the code). >= 0, smaller
  is better. inf for a failed NNLS, for a matrix np.linalg.solve finds
  singular, and for m - trace K <= 1e-10. The condition of the matrix is
  not checked: on DRT systems trace K agrees with an exact reference to 6+
  digits from lambda 1e-12 up (the search starts at 1e-11); far below
  (1e-20 and less) rounding dominates it and the score is finite but
  wrong (measured 3x too low at 1e-22).

compute_lcurve_point(lambda_val, A, b, L) -> (rho, eta, x).
  rho = log10 ||A x - b||, eta = log10 ||L x||, each norm floored at 1e-15;
  x = x(lambda). (inf, inf, None) for a failed NNLS.

compute_lcurve_curvature(rho, eta) -> array, same length.
  Signed curvature (rho' eta'' - rho'' eta') / (rho'^2 + eta'^2)^1.5, the
  derivatives by np.gradient over the index (the lambdas are log-spaced).
  0 where the denominator is <= 1e-15, all zeros for fewer than 3 points.
  NaN within two points of a failed point (rho or eta not finite): the
  second difference reaches that far.

find_lcurve_corner(lambda_values, rho, eta) -> (lambda_corner, index,
curvature).
  The corner is the interior point (not the first or last) of largest
  finite curvature, at a tie the lowest index; with no interior point, the
  middle one; with no finite interior curvature, the first interior point
  (its curvature NaN). Increasing lambda runs the curve from small residual
  to small norm, so the corner turns counter-clockwise and its curvature is
  positive; the caller checks the sign. A failed point (inf) never decides
  the corner, neither itself nor through a neighbour's finite differences.
  The caller passes at least one point.

find_optimal_lambda_gcv(A, b, L, lambda_range=DRT_LAMBDA_RANGE, n_search=20)
-> (lambda, score).
  Coarse search: n_search log-spaced lambdas over lambda_range. Fine
  search: n_search log-spaced lambdas between the coarse neighbours of the
  coarse minimum; at an end of the range, from that end one decade beyond
  it to two coarse steps inside. The result is the lambda of lowest score
  on the fine grid, at a tie the lowest lambda. The fine grid does not
  hold the coarse minimum (n_search even), so its best can score slightly
  above it (1 of 40 random spectra, by 3.4e-4 relative); with no finite
  score on the fine grid, the coarse minimum. lambda lies in
  [lambda_range[0] / 10, lambda_range[1] * 10]. ValueError when no lambda
  of the coarse search gives a finite score (also where b does not match
  A, or A or b holds NaN or inf). The caller passes 0 < lo < hi finite and
  n_search >= 3; other values are not checked.

find_optimal_lambda_hybrid(A, b, L, lambda_range=DRT_LAMBDA_RANGE,
n_search=20, lcurve_decades=1.5) -> (lambda, score, diagnostics).
  lambda_gcv from find_optimal_lambda_gcv. L-curve: n_search log-spaced
  lambdas from c / 10^lcurve_decades to c x 10^lcurve_decades, c =
  lambda_gcv clipped to lambda_range, the window clipped to
  [lambda_range[0] / 10, lambda_range[1] * 10]; lambda_lcurve is its
  corner. lambda = max(lambda_gcv, lambda_lcurve), method_used 'lcurve'
  if lambda_lcurve >= lambda_gcv, else 'gcv'. So lambda >= lambda_gcv
  always. Where the corner's curvature is not positive (a flat curve, as
  for b = 0, or only failed points nearby) there is no corner: lambda =
  lambda_gcv, lambda_lcurve None, method_used 'gcv', both flags False.
  score: the GCV score at the returned lambda.
  diagnostics: lambda_gcv, lambda_lcurve, method_used, curvature, rho,
  eta, lambda_values (the window), corner_at_edge (L-curve won and its
  corner is within one point of a window end), corner_below_gcv (GCV won
  and lambda_lcurve < lambda_gcv / 10). ValueError as for
  find_optimal_lambda_gcv. The caller passes lcurve_decades finite and
  > 0; it is not checked.

Outside what the caller guarantees nothing is promised: lambda_val
finite and > 0 (0 returns numbers, though x is defined for lambda > 0
only), entries of A and b small enough that A^T A does not overflow
(~1e150: score inf or NaN).

What must hold:
  - The inputs are not changed. The result is deterministic.
  - Raises only ValueError and LinAlgError (the caller's fallback catches
    these), and only where stated; no RuntimeWarning (the NNLS of a
    negative or non-finite lambda_val aside).
  - Units: b scaled by c > 0 (Z in other units) leaves every lambda,
    index, curvature, method_used and flag unchanged (to rounding), scales
    scores by c^2 and shifts rho and eta by log10 c; except where a norm
    is below the 1e-15 floor.
  - The rows of (A, b) permuted together change nothing (to rounding).
  - A shape mismatch between A, b and L is a ValueError.
"""

import numpy as np
import logging
from typing import Any, Dict, Tuple, Optional
from numpy.typing import NDArray
from scipy.optimize import nnls

from ..fitting.config import DRT_LAMBDA_RANGE

logger = logging.getLogger(__name__)

# NNLS iteration budget per unknown. scipy's default 3*n gives up on weighted
# DRT systems with |Z| over many decades at small lambda (needed ~5*n on a
# 7-decade oxide spectrum); Lawson-Hanson only runs longer when it needs to.
DRT_NNLS_MAXITER_FACTOR = 50


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
    - L: regularizační matice [(N_tau-2) × N_tau]: druhá derivace podle ln τ,
         škálovaná tak, aby λ‖Lx‖² byl λ·∫(γ'')² d ln τ vůči průměrnému
         čtverci rezidua (viz _build_drt_matrices)

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
        x, residual_nnls = nnls(A_reg, b_reg, maxiter=DRT_NNLS_MAXITER_FACTOR * A_reg.shape[1])
    except (RuntimeError, ValueError) as e:
        logger.debug(f"NNLS selhalo pro λ={lambda_val:.2e}: {e}")
        return np.inf

    # Reziduum původního problému (bez regularizace)
    residual = b - A @ x
    residual_norm_sq = float(np.sum(residual**2))

    # Výpočet trace(I - K) pomocí přímé inverze
    # K = A @ inv(A^T @ A + λ·L^T·L) @ A^T
    # trace(I - K) = n - trace(K)
    #
    # Trace se záměrně bere z neomezeného řešení přes VŠECHNY sloupce, ačkoli
    # reziduum je z NNLS (γ ≥ 0). "Správnější" K jen nad aktivní množinou
    # (γ > 0) bylo změřeno a vybírá hůř: aktivní množinu určují sama data, takže
    # stupně volnosti z ní vycházejí podhodnocené a GCV pak přeje malé λ.
    # Naměřeno (2026-09-29): real_gamry_example λ 9.4e-6 -> 2.4e-10 a 2 -> 5
    # píků; syntetika se 2 procesy a 2 % šumu 2 -> 3 píky; λ_gcv navíc
    # přestala růst se šumem (0.5 % -> 6.2e-7, 2 % -> 2.3e-7). Plný trace dává
    # konzervativní (větší) λ; hybrid ho rohem L-křivky může jen zvýšit.

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
        x, _ = nnls(A_reg, b_reg, maxiter=DRT_NNLS_MAXITER_FACTOR * A_reg.shape[1])
    except (RuntimeError, ValueError) as e:
        logger.debug(f"NNLS selhalo pro λ={lambda_val:.2e}: {e}")
        return np.inf, np.inf, None

    # Norma rezidua ||Ax - b||
    residual = b - A @ x
    residual_norm = float(np.linalg.norm(residual))

    # Regularizační norma ||Lx||
    reg_norm = float(np.linalg.norm(L @ x))

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

    # Selhaný bod (inf) dává inf/NaN i v derivacích sousedů; np.where dělí
    # všude, i kde se výsledek zahodí
    with np.errstate(invalid='ignore', divide='ignore'):
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

    # Druhá derivace sahá dva body daleko: tam křivost ze selhaného bodu
    # není číslo, i když vyjde konečná
    failed = ~(np.isfinite(rho) & np.isfinite(eta))
    if failed.any():
        curvature[np.convolve(failed, np.ones(5), mode='full')[2:-2] > 0] = np.nan

    return curvature


def find_lcurve_corner(lambda_values: NDArray[np.float64],
                        rho: NDArray[np.float64],
                        eta: NDArray[np.float64]) -> Tuple[float, int, NDArray]:
    """
    Najde roh L-křivky (bod maximální křivosti).

    Konvence znaménka: L-křivka (ρ = log||Ax-b||, η = log||Lx||) procházená
    rostoucím λ jde z malého rezidua / velké normy řešení do velkého rezidua /
    malé normy (ρ roste, η klesá). Roh je konvexní směrem k počátku, tedy
    levotočivá (CCW) zatáčka, která má v této orientaci **kladnou** znaménkovou
    křivost. Proto se roh hledá jako `argmax` (kladné) křivosti, ne `argmin`.
    Tuto konvenci hlídají testy v `tests/test_lcurve_corner.py`.

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
        # Maximum kladné křivosti; NaN (selhaný bod) by argmax vybral
        corner_idx_inner = np.argmax(np.where(np.isfinite(curvature_inner),
                                              curvature_inner, -np.inf))
        corner_idx = int(corner_idx_inner) + 1  # Kompenzace za slice

    lambda_corner = lambda_values[corner_idx]

    return lambda_corner, corner_idx, curvature


def find_optimal_lambda_gcv(A: NDArray[np.float64], b: NDArray[np.float64],
                            L: NDArray[np.float64],
                            lambda_range: Tuple[float, float] = DRT_LAMBDA_RANGE,
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
    logger.debug("Automatický výběr λ pomocí GCV...")

    # Fáze 1: Hrubé prohledání v log-prostoru
    lambda_values = np.logspace(np.log10(lambda_range[0]),
                                 np.log10(lambda_range[1]),
                                 n_search)

    gcv_scores = []
    for lam in lambda_values:
        score = compute_gcv_score(lam, A, b, L)
        gcv_scores.append(score)

    gcv_scores_arr = np.array(gcv_scores)

    # Najdi minimum
    min_idx = int(np.argmin(gcv_scores_arr))
    lambda_coarse = lambda_values[min_idx]
    gcv_coarse = gcv_scores_arr[min_idx]

    logger.debug(f"Hrubé prohledání: λ = {lambda_coarse:.4e}, GCV = {gcv_coarse:.4e}")
    # argmin nad samými inf by vrátil první λ jako platný výsledek; volající
    # na ValueError přejde na záložní λ
    if not np.isfinite(gcv_coarse):
        raise ValueError("no lambda in the GCV search gives a finite score")

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

    gcv_fine_arr = np.array(gcv_fine)
    min_fine_idx = int(np.argmin(gcv_fine_arr))
    lambda_optimal = lambda_fine[min_fine_idx]
    gcv_optimal = gcv_fine_arr[min_fine_idx]
    if not np.isfinite(gcv_optimal):
        return lambda_coarse, gcv_coarse

    logger.debug(f"✓ Optimální λ = {lambda_optimal:.4e} (GCV = {gcv_optimal:.4e})")

    return lambda_optimal, gcv_optimal


def find_optimal_lambda_hybrid(A: NDArray[np.float64], b: NDArray[np.float64],
                                L: NDArray[np.float64],
                                lambda_range: Tuple[float, float] = DRT_LAMBDA_RANGE,
                                n_search: int = 20,
                                lcurve_decades: float = 1.5
                                ) -> Tuple[float, float, dict]:
    """
    Hybridní metoda pro výběr λ: GCV a roh L-křivky, vyhrává větší λ.

    Strategie:
    1. λ_gcv z find_optimal_lambda_gcv (hrubé + jemné hledání)
    2. L-curve prohledání ±lcurve_decades dekád kolem λ_gcv oříznutého do
       lambda_range (jemné GCV může skončit dekádu za mezí)
    3. λ = max(λ_gcv, λ_lcurve) - obě kritéria na NNLS spíš podregularizují

    GCV předpokládá lineární řešení, NNLS má nelineární constraint (x ≥ 0).
    L-křivka to nepředpokládá, ale na celém rozsahu může mít falešný roh
    (na real_gamry_example roh u λ = 0.18, zjevně přehlazený), proto se
    hledá jen v okně kolem λ_gcv.

    Parametry:
    - A: matice systému [2N × N_tau]
    - b: pravá strana [2N]
    - L: regularizační matice
    - lambda_range: tuple (min, max) pro počáteční GCV prohledávání
    - n_search: počet bodů pro prohledávání
    - lcurve_decades: rozsah L-curve kolem GCV odhadu v dekádách (default: 1.5)

    Vrací:
    - tuple: (lambda_optimal, score, diagnostics)
      - diagnostics obsahuje: lambda_gcv, lambda_lcurve, method_used, curvature,
        rho, eta, corner_at_edge, corner_below_gcv
    """
    logger.debug("Hybridní výběr λ (GCV a roh L-křivky, platí větší)...")

    diagnostics: Dict[str, Any] = {
        'lambda_gcv': None,
        'lambda_lcurve': None,
        'method_used': None,
        'curvature': None,
        'rho': None,
        'eta': None,
        'lambda_values': None
    }

    # === Fáze 1: GCV ===
    # Zpřesněné, ne jen z hrubé mřížky (krok 0.42 dekády): když GCV vyhraje,
    # je λ_gcv přímo výsledkem.
    logger.debug("Fáze 1: GCV...")
    lambda_gcv, gcv_at_gcv = find_optimal_lambda_gcv(A, b, L, lambda_range, n_search)

    diagnostics['lambda_gcv'] = lambda_gcv
    logger.debug(f"  GCV minimum: λ_gcv = {lambda_gcv:.4e}")

    # === Fáze 2: L-curve v okolí GCV odhadu ===
    logger.debug("Fáze 2: roh L-křivky kolem λ_gcv...")

    # Rozsah pro L-curve: ±lcurve_decades dekád kolem λ_gcv. Jemné hledání GCV
    # může skončit až dekádu za mezí rozsahu; střed se proto ořízne do rozsahu,
    # jinak by ořezání okna níže nechalo roh jen na jedné straně GCV.
    center = float(np.clip(lambda_gcv, lambda_range[0], lambda_range[1]))
    lcurve_min = max(center / (10**lcurve_decades), lambda_range[0] / 10)
    lcurve_max = min(center * (10**lcurve_decades), lambda_range[1] * 10)

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

    rho_arr = np.array(rho)
    eta_arr = np.array(eta)

    diagnostics['rho'] = rho_arr
    diagnostics['eta'] = eta_arr
    diagnostics['lambda_values'] = lambda_lcurve_range

    # Najdi roh L-křivky
    lambda_lcurve, corner_idx, curvature = find_lcurve_corner(lambda_lcurve_range, rho_arr, eta_arr)

    diagnostics['lambda_lcurve'] = lambda_lcurve
    diagnostics['curvature'] = curvature

    # Háček: roh na okraji prohledávaného L-rozsahu znamená, že skutečný roh
    # leží nejspíš mimo okno (úzké lcurve_decades) — λ_lcurve je pak
    # nespolehlivé a stojí za to rozsah rozšířit.
    n_lc = len(lambda_lcurve_range)
    corner_at_edge = corner_idx <= 1 or corner_idx >= n_lc - 2
    if corner_at_edge:
        logger.debug(
            f"  L-curve roh na okraji rozsahu (index {corner_idx}/{n_lc - 1}) "
            f"— zvaž větší lcurve_decades"
        )

    logger.debug(f"  L-curve roh: λ_lcurve = {lambda_lcurve:.4e}")

    # Bez kladné křivosti (rovná křivka, např. b = 0) nebo jen se selhanými
    # body v okolí není roh, jen největší ze zbytku: λ nezvyšuje
    if not curvature[corner_idx] > 0:
        logger.debug("  L-křivka bez kladné křivosti: bez rohu, platí GCV")
        diagnostics['lambda_lcurve'] = None
        diagnostics['method_used'] = 'gcv'
        diagnostics['corner_at_edge'] = False
        diagnostics['corner_below_gcv'] = False
        return lambda_gcv, gcv_at_gcv, diagnostics

    # === Fáze 3: Rozhodnutí ===
    # Větší z obou λ. Obě kritéria na NNLS-DRT chybují stejným směrem - λ
    # příliš malá: změřeno (2026-09-29) na dvou-ZARC syntetice s 0.1-0.5 %
    # šumu je hybridní λ o 1-3 dekády pod hodnotou, která dá správný počet
    # píků (4-5 místo 2). Kompromis (geometrický průměr) by nebyl optimem
    # ani jednoho kritéria; skutečnou neshodu hlásí corner_below_gcv níže.
    lambda_optimal = max(lambda_gcv, lambda_lcurve)
    diagnostics['method_used'] = 'lcurve' if lambda_lcurve >= lambda_gcv else 'gcv'
    logger.debug(f"  Větší z obou: {diagnostics['method_used']} "
                 f"(ratio lcurve/gcv = {lambda_lcurve / lambda_gcv:.2f})")

    # Roh na okraji okna je nespolehlivý jen když dal λ; když vyhrálo GCV,
    # roh se zahodil.
    diagnostics['corner_at_edge'] = corner_at_edge and diagnostics['method_used'] == 'lcurve'

    # Roh víc než dekádu pod GCV: kritéria se opravdu rozcházejí a λ z GCV
    # nemá oporu v L-křivce. Roh nad GCV je očekávaný efekt NNLS, žádná
    # neshoda. (Roh na spodním okraji neoříznutého okna je vždy > 1.3 dekády
    # pod GCV; na oříznutém je GCV samo na mezi a hlásí to lambda_at_edge.)
    diagnostics['corner_below_gcv'] = (
        diagnostics['method_used'] == 'gcv' and lambda_lcurve * 10 < lambda_gcv
    )

    # Finální GCV score pro reporting; když vyhrálo GCV, už je spočítané
    final_gcv = (gcv_at_gcv if diagnostics['method_used'] == 'gcv'
                 else compute_gcv_score(lambda_optimal, A, b, L))

    logger.debug(f"✓ Hybridní optimum: λ = {lambda_optimal:.4e} (GCV = {final_gcv:.4e})")

    return lambda_optimal, final_gcv, diagnostics

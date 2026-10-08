"""
Differential Evolution global optimization for circuit fitting.

Clean design: No logging in core functions, all diagnostics returned as data.

Strategy options:
    1 = 'randtobest1bin' (default) - balanced exploration/exploitation
    2 = 'best1bin' - faster convergence, may miss global minimum
    3 = 'rand1bin' - better exploration, slower convergence
"""

import numpy as np
import logging
import multiprocessing
import warnings
from typing import Tuple, List, Optional, Any
from numpy.typing import NDArray
from dataclasses import dataclass, field
from scipy.optimize import differential_evolution, OptimizeWarning

from .circuit import FitResult, FitDiagnostics, Circuit
from .bounds import (generate_simple_bounds, build_bound_status, log_scale_ci_mask,
                     log_search_bounds,
                     validate_fixed_params)
from .covariance import compute_covariance_matrix
from .diagnostics import compute_residual_weights, compute_fit_metrics, compute_significance
from .optimizer import least_squares_normalized
from .de_archive import Refinement, choose, select_archive_candidates, selection_warnings
from .jacobian import make_jacobian_function
from .config import DE_STALLED_ERROR_PCT, DE_STALLED_IMPROVEMENT_FACTOR

logger = logging.getLogger(__name__)

# Strategy mapping: number -> scipy strategy name
DE_STRATEGIES = {
    1: 'randtobest1bin',
    2: 'best1bin',
    3: 'rand1bin',
}

# rand1bin by default: randtobest1bin and best1bin pull the whole population
# toward the current best member, and when an early leader sits in a
# degenerate basin - an arc collapsed to R -> 0, a Wo whose tau runs past the
# window so it acts as W - DE stops there after ~45 generations instead of
# ~270, and least_squares cannot leave it. Measured on 7 circuits (the four
# ZScope references, R-(R|Wo), R-(R|Ws), R-((R-Wo)|Q)), noise-free and 1 %,
# 10 seeds each (140 fits): randtobest1bin failed 7 (two time constants 4/20,
# 5-11 % error), rand1bin 0, for about twice the DE time. A smaller population
# (popsize 10) halves that cost but already missed one Randles-Wo fit.
DEFAULT_DE_STRATEGY = 3

# How far inside its bounds the DE starting point is held, as a fraction of
# the bound span. differential_evolution rescales x0 to [0, 1] as
# (x - midpoint) / span + 0.5 and rejects the result if it falls outside, so a
# value sitting exactly ON a bound can come back as -1.1e-16 and raise
# "Some entries in x0 lay outside the specified bounds". 1e-9 of the span is
# some nine orders of magnitude above that rounding error and still far below
# any parameter's physical resolution - on the CPE exponent's (0.3, 1.0) range
# it moves the start by 7e-10.
DE_X0_BOUND_MARGIN = 1e-9


class _DECostFunction:
    """Picklable cost function for differential evolution with workers > 1."""

    def __init__(self, circuit: Circuit, frequencies: NDArray[np.float64],
                 Z: NDArray[np.complexfloating], weights: NDArray[np.float64],
                 fixed_params: Optional[List[bool]] = None,
                 full_initial_guess: Optional[List[float]] = None):
        self.circuit = circuit
        self.frequencies = frequencies
        self.Z = Z
        self.weights = weights
        self.fixed_params = fixed_params
        self.full_initial_guess = full_initial_guess
        # A list to record every (cost, params) evaluation into, or None. Only
        # set for workers == 1: with more, DE calls pickled copies in worker
        # processes, and the map in fit_circuit_diffevo records in the parent.
        self.archive: Optional[List[Tuple[float, NDArray[np.float64]]]] = None

    def _reconstruct_params(self, free_params):
        if self.fixed_params is None or not any(self.fixed_params):
            return list(free_params)
        full, idx = [], 0
        for i, is_fixed in enumerate(self.fixed_params):
            if is_fixed:
                full.append(self.full_initial_guess[i])
            else:
                full.append(free_params[idx])
                idx += 1
        return full

    def __call__(self, params):
        full_params = self._reconstruct_params(params)
        Z_pred = self.circuit.impedance(self.frequencies, full_params)
        residuals_real = (self.Z.real - Z_pred.real) * self.weights
        residuals_imag = (self.Z.imag - Z_pred.imag) * self.weights
        cost = np.sum(residuals_real**2 + residuals_imag**2)
        if self.archive is not None:
            self.archive.append((float(cost), np.array(params, dtype=float)))
        return cost


def _to_linear(x, log_mask: NDArray[np.bool_]) -> NDArray[np.float64]:
    """Map a DE search vector back to physical parameters (10**x where masked)."""
    x = np.asarray(x, dtype=float)
    return np.where(log_mask, 10.0 ** x, x)


class _LogSpaceCost:
    """Picklable adapter that lets DE search log10 of the masked parameters.

    A closure would not survive pickling for workers > 1, hence a class -
    same reason _DECostFunction is one.
    """

    def __init__(self, cost_function: _DECostFunction, log_mask: NDArray[np.bool_]):
        self.cost_function = cost_function
        self.log_mask = log_mask

    def __call__(self, params):
        return self.cost_function(_to_linear(params, self.log_mask))


@dataclass
class DiffEvoDiagnostics:
    """Diagnostics from differential evolution optimization."""
    # DE settings
    strategy: str
    popsize: int
    maxiter: int
    tol: float
    workers: int
    weighting: str
    jacobian_type: str

    # DE results
    de_converged: bool
    de_iterations: int
    de_evaluations: int
    de_error: float

    # Refinement results
    refined_error: float
    refinement_improved: bool
    total_evaluations: int

    # Optimized objective values (weighted SSR, S = sum w^2 |dZ|^2).
    # These drive the DE-vs-refinement selection; the *_error fields above
    # are the human-readable weighted mean relative error (%) used for display.
    de_cost: float = 0.0
    refined_cost: float = 0.0

    # Archive check (fitting/de_archive.py): whether it ran (archive_check),
    # how many early-generation candidates refined, and whether the result is
    # one of them.
    archive_checked: bool = False
    archive_candidates: int = 0
    archive_used: bool = False

    # Fixed params info
    n_fixed_params: int = 0
    fixed_param_indices: List[int] = field(default_factory=list)

    # Free parameters DE searched as log10(value) - those whose bounds span
    # enough decades that a linear population would sample only the top one.
    log_search_params: List[str] = field(default_factory=list)

    # Initial guess passed to DE (full parameter vector, including fixed).
    # Captured before the optimizer runs so it survives circuit.update_params()
    # at the end of fit_circuit_diffevo.
    initial_guess: List[float] = field(default_factory=list)

    # Warnings
    warnings: List[str] = field(default_factory=list)


@dataclass
class DiffEvoResult:
    """
    Result from differential evolution optimization.

    Attributes
    ----------
    best_result : FitResult
        Best fitting result after least_squares refinement
    de_result : scipy.optimize.OptimizeResult
        Raw result from differential_evolution
    de_error : float
        Fit error after DE (before refinement) [%]
    final_error : float
        Fit error after least_squares refinement [%]
    n_evaluations : int
        Number of function evaluations
    strategy : str
        DE strategy used
    improvement : float
        Relative improvement from DE to refined [%]
    diagnostics : DiffEvoDiagnostics
        Detailed diagnostics
    """
    best_result: FitResult
    de_result: Any
    de_error: float
    final_error: float
    n_evaluations: int
    strategy: str
    improvement: float
    diagnostics: Optional[DiffEvoDiagnostics] = None


def fit_circuit_diffevo(
    circuit: Circuit,
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complexfloating],
    strategy: int = DEFAULT_DE_STRATEGY,
    popsize: int = 15,
    maxiter: int = 1000,
    tol: float = 0.01,
    workers: int = 1,
    weighting: str = 'modulus',
    use_analytic_jacobian: bool = True,
    seed: Optional[int] = None,
    archive_check: bool = True
) -> Tuple[DiffEvoResult, NDArray[np.complexfloating]]:
    """
    Fit circuit using Differential Evolution global optimization.

    Parameters
    ----------
    circuit : Circuit
        Circuit object with initial parameter guesses
    frequencies : ndarray
        Frequency array [Hz]
    Z : ndarray
        Complex impedance data [Ohm]
    strategy : int, optional
        DE strategy: 1='randtobest1bin', 2='best1bin', 3='rand1bin'
        (default: 3, see DEFAULT_DE_STRATEGY)
    popsize : int, optional
        Population size multiplier (default: 15)
    maxiter : int, optional
        Maximum number of generations (default: 1000)
    tol : float, optional
        Relative tolerance for convergence (default: 0.01)
    workers : int, optional
        Number of parallel workers (default: 1)
    weighting : str, optional
        Weighting scheme
    use_analytic_jacobian : bool, optional
        Use analytic Jacobian for refinement (default: True)
    seed : int, optional
        Seed for differential_evolution's random generator. Default None
        (non-deterministic). Set an int for reproducible runs (e.g. tests).
    archive_check : bool, optional
        Refine early-generation candidates from DE's evaluation archive
        against local minima and ambiguous models (default True; ~30 extra
        least_squares runs, see fitting/de_archive.py).

    Returns
    -------
    diffevo_result : DiffEvoResult
        Differential evolution result with all diagnostics
        (plot: `plot_circuit_fit(frequencies, Z, diffevo_result.best_result)`)
    Z_fit : ndarray
        Best fit impedance
    """
    strategy_name = DE_STRATEGIES.get(strategy, DE_STRATEGIES[DEFAULT_DE_STRATEGY])
    diag_warnings = []

    # Get initial guess from circuit definition
    initial_guess_full = list(circuit.get_all_params())

    param_labels = circuit.get_param_labels()
    lower_bounds_full, upper_bounds_full = generate_simple_bounds(param_labels)
    fixed_params = circuit.get_all_fixed_params()
    fixed_param_indices = [i for i, f in enumerate(fixed_params) if f]

    # Indexed labels R0, R1, Q0, ... (same as _prepare_optimization)
    param_labels_indexed = [f"{label}{param_labels[:i].count(label)}"
                            for i, label in enumerate(param_labels)]

    # Raises when every parameter is fixed (empty optimization vector)
    diag_warnings.extend(validate_fixed_params(
        initial_guess_full, lower_bounds_full, upper_bounds_full,
        fixed_params, param_labels_indexed
    ))

    # Filter to free parameters only
    if any(fixed_params):
        initial_guess = np.array([v for v, f in zip(initial_guess_full, fixed_params) if not f])
        lower_bounds = [lb for lb, f in zip(lower_bounds_full, fixed_params) if not f]
        upper_bounds = [ub for ub, f in zip(upper_bounds_full, fixed_params) if not f]
        free_labels = [lab for lab, f in zip(param_labels, fixed_params) if not f]
    else:
        initial_guess = np.array(initial_guess_full)
        lower_bounds = lower_bounds_full
        upper_bounds = upper_bounds_full
        free_labels = param_labels

    # Clip initial guess to bounds
    initial_guess = np.clip(initial_guess, lower_bounds, upper_bounds)

    # DE samples its population uniformly over the bounds, so a scale parameter
    # whose bounds span decades (R: 1e-4..1e10, C: 1e-15..1e-1) would be drawn
    # almost exclusively from its top decade: the parallel branches are then
    # shorted, every member predicts the series resistance alone, and the
    # population energies are so nearly equal that DE's convergence test
    # (std <= tol * |mean|) can fire after the first generation. Searching
    # log10 of those parameters spreads the population over the decades
    # instead. The CPE exponent n (0.3-1.0) keeps its linear scale; G too is
    # searched in log space, despite its zero lower bound - see
    # log_search_bounds.
    log_mask_list, de_lower, de_upper = log_search_bounds(
        list(lower_bounds), list(upper_bounds), free_labels
    )
    log_mask = np.array(log_mask_list, dtype=bool)
    # Labels of the log-searched parameters, for the diagnostics line. The mask
    # is in free-parameter space; map it back to full-space labels.
    free_indices = [i for i, f in enumerate(fixed_params) if not f]
    log_search_params = [
        param_labels_indexed[full_i]
        for free_i, full_i in enumerate(free_indices) if log_mask[free_i]
    ]

    # Precompute weights
    weights = compute_residual_weights(Z, weighting)

    # Cost function for DE
    cost_function = _DECostFunction(
        circuit, frequencies, Z, weights,
        fixed_params=fixed_params,
        full_initial_guess=initial_guess_full
    )

    reconstruct_params = cost_function._reconstruct_params

    # Residual function for least_squares
    def residual_function(params):
        full_params = reconstruct_params(params)
        Z_pred = circuit.impedance(frequencies, full_params)
        return np.concatenate([
            (Z.real - Z_pred.real) * weights,
            (Z.imag - Z_pred.imag) * weights
        ])

    # Step 1: Run Differential Evolution, in the search space chosen above
    de_objective = (_LogSpaceCost(cost_function, log_mask) if log_mask.any()
                    else cost_function)
    de_bounds = list(zip(de_lower, de_upper))
    # G's initial guess may be exactly 0, which has no logarithm. Floor it
    # before the transform; the clip puts it back inside the search bounds.
    #
    # Inside, not onto: initial_guess was clipped to these same bounds above,
    # so any guess outside them lands exactly on one - a CPE exponent n <= 0.3
    # written into --circuit, for instance - and scipy rejects such an x0 (see
    # DE_X0_BOUND_MARGIN).
    x0_positive = np.maximum(initial_guess, np.finfo(float).tiny)
    de_lower_arr = np.asarray(de_lower, dtype=float)
    de_upper_arr = np.asarray(de_upper, dtype=float)
    x0_margin = DE_X0_BOUND_MARGIN * (de_upper_arr - de_lower_arr)
    de_x0 = np.clip(np.where(log_mask, np.log10(x0_positive), initial_guess),
                    de_lower_arr + x0_margin, de_upper_arr - x0_margin)

    archive: Optional[list] = [] if archive_check else None
    pool = None
    de_workers: Any = workers
    if archive is not None:
        if workers == 1:
            cost_function.archive = archive
        else:
            map_func: Any = workers
            if not callable(workers):
                pool = multiprocessing.Pool(None if workers == -1 else workers)
                map_func = pool.map

            def de_workers(func, xs):
                # Worker processes keep their own archives; record in the parent
                xs = [np.asarray(x, dtype=float) for x in xs]
                costs = list(map_func(func, xs))
                archive.extend((float(c), _to_linear(x, log_mask)) for c, x in zip(costs, xs))
                return costs
    try:
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            de_result = differential_evolution(
                de_objective,
                de_bounds,
                x0=de_x0,
                strategy=strategy_name,
                popsize=popsize,
                maxiter=maxiter,
                tol=tol,
                workers=de_workers,
                polish=False,
                seed=seed,
                disp=False,
                updating='deferred' if workers != 1 else 'immediate',
            )
    except Exception as e:
        raise RuntimeError(f"DE optimization failed: {e}") from e
    finally:
        if pool is not None:
            pool.close()
            pool.join()
    # Stop recording: the cost evaluations below are not part of the search.
    cost_function.archive = None

    # Back to physical parameters right away: everything below (refinement
    # start, cost comparison, diagnostics, DiffEvoResult.de_result) works in
    # the linear parameter space and stays unaware of the DE transform.
    de_result.x = _to_linear(de_result.x, log_mask)

    # Compute DE error
    de_params_full = reconstruct_params(de_result.x)
    Z_fit_de = circuit.impedance(frequencies, de_params_full)
    de_metrics = compute_fit_metrics(Z, Z_fit_de, weighting)
    de_error_rel = de_metrics[0]

    # Step 2: Refine with least_squares
    jacobian_type = 'analytic'
    if use_analytic_jacobian:
        try:
            jac_func = make_jacobian_function(
                circuit, frequencies, weights,
                fixed_params=fixed_params,
                full_initial_guess=initial_guess_full
            )
        except NotImplementedError:
            jac_func = '2-point'
            jacobian_type = 'numeric'
    else:
        jac_func = '2-point'
        jacobian_type = 'numeric'

    # One least_squares setup for every start - DE's point and the archive
    # candidates - so their costs compare like for like.
    lb_arr, ub_arr = np.asarray(lower_bounds, float), np.asarray(upper_bounds, float)

    def refine(x0):
        x0 = np.clip(np.asarray(x0, dtype=float), lb_arr, ub_arr)
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always", OptimizeWarning)
            return least_squares_normalized(
                residual_function, x0, jac_func, (lower_bounds, upper_bounds),
                method='trf', ftol=1e-10, xtol=1e-10, max_nfev=5000,
            )

    def spectrum(x_free):
        full = list(reconstruct_params(x_free))
        return full, circuit.impedance(frequencies, full)

    refine_error = None
    try:
        r = refine(de_result.x)
        full, Z_r = spectrum(r.x)
        from_de: Optional[Refinement] = Refinement(float(np.sum(r.fun ** 2)), Z_r, full, r)
    except Exception as e:
        from_de, refine_error = None, e

    # Step 2b: second look through the evaluation archive (see de_archive.py):
    # the best distinct point of each early window of generations is refined
    # the same way. A failed refinement must not masquerade as a successful
    # one, so a start that raises simply contributes nothing.
    from_archive: List[Refinement] = []
    if archive:
        arch_costs = np.array([c for c, _ in archive])
        arch_params = np.array([x for _, x in archive])
        archive = None                       # the arrays above are all that is needed
        n_pop = max(5, popsize * len(free_labels))   # scipy's population size
        candidates = select_archive_candidates(arch_costs, arch_params, n_pop, ~log_mask)[1:]  # [0] = DE's best
        for i in candidates:
            try:
                r = refine(arch_params[i])
            except Exception:
                continue
            full, Z_r = spectrum(r.x)
            from_archive.append(Refinement(float(np.sum(r.fun ** 2)), Z_r, full, r))
    archive_nfev = sum(r.result.nfev for r in from_archive)

    # Selection and improvement use the *optimized* objective (weighted SSR,
    # S = sum w^2 |dZ|^2), not the weighted mean relative error: DE and
    # least_squares both minimize S, so choosing on a different metric could
    # discard a genuinely better refined fit. One decision among all of them,
    # before anything derived from the choice is computed.
    de_cost = float(cost_function(de_result.x))
    sel = choose(de_cost, from_de, from_archive, weights, len(free_labels))
    best = sel.best_refined
    best_metrics = compute_fit_metrics(Z, best.Z, weighting) if best is not None else de_metrics
    ls_error_rel = best_metrics[0]
    refined_cost = best.cost if best is not None else de_cost
    improvement = (de_cost - refined_cost) / de_cost * 100 if de_cost > 0 else 0

    chosen = sel.chosen
    final_ls = chosen.result if chosen is not None else None
    used_refinement = chosen is not None
    archive_used = chosen is not None and chosen is not from_de
    if chosen is not None:
        params_opt_free = np.array(chosen.result.x)
        params_opt = np.array(chosen.params)
        Z_fit = chosen.Z
        fit_metrics = best_metrics          # the chosen refinement is the best one
    else:
        params_opt_free = np.array(de_result.x)
        params_opt = np.array(de_params_full)
        Z_fit = Z_fit_de
        fit_metrics = de_metrics

    if refine_error is not None:
        diag_warnings.append(f"Refinement failed: {refine_error}, using DE result" if not used_refinement
                             else f"Refinement from the DE result failed: {refine_error}")
    elif not used_refinement:
        diag_warnings.append("Refinement worsened fit, using DE result")
    fit_error_rel, fit_error_abs, quality = fit_metrics
    diag_warnings.extend(selection_warnings(sel, fit_error_rel, param_labels_indexed))

    # The global stage is only useful if it actually explored. When it ends far
    # from the data and a local refinement then improves by an order of
    # magnitude, the reported fit came from that local run, not from DE. (An
    # archive repair into another model has its own warning above.)
    if (used_refinement and not sel.local_minimum and de_error_rel > DE_STALLED_ERROR_PCT
            and fit_error_rel * DE_STALLED_IMPROVEMENT_FACTOR < de_error_rel):
        diag_warnings.append(
            f"Global search contributed nothing: DE stopped after "
            f"{de_result.nit} iteration(s) at {de_error_rel:.1f}% error and the "
            f"local refinement reached {fit_error_rel:.1f}% on its own. The fit "
            "rests on that single local run. Check that the circuit suits the "
            "data, then raise --de-maxiter or lower --de-tol"
        )

    # Step 3: Compute covariance
    # Both the residuals and the Jacobian must be evaluated at the *chosen*
    # point (params_opt_free). When the DE result is kept, the refinement's jac is the
    # Jacobian at the LS point, not the chosen one, so it must not be reused.
    final_residuals = residual_function(params_opt_free)
    n_params_full = len(params_opt)

    if jacobian_type == 'analytic':
        jac_at_opt = jac_func(params_opt_free)
    elif final_ls is not None and getattr(final_ls, 'jac', None) is not None:
        # Numeric Jacobian; LS point coincides with the chosen point.
        jac_at_opt = final_ls.jac
    else:
        jac_at_opt = None

    cov_result = None
    if jac_at_opt is not None:
        cov_result = compute_covariance_matrix(
            jacobian=jac_at_opt,
            residuals=final_residuals,
            n_params=n_params_full,
            fixed_params=fixed_params
        )
        params_stderr = cov_result.stderr
        cov = cov_result.cov
        condition_number = cov_result.condition_number
        is_well_conditioned = cov_result.is_well_conditioned
    else:
        params_stderr = np.full(n_params_full, np.inf)
        cov = None
        condition_number = np.inf
        is_well_conditioned = False

    # Step 4: Update circuit with fitted parameters
    circuit.update_params(list(params_opt))

    # Per-parameter bound status and derived warnings — same contract as
    # fit_equivalent_circuit (full-space indices, classify_bound_status
    # criterion).
    bound_status = build_bound_status(
        params_opt, lower_bounds_full, upper_bounds_full, fixed_params
    )
    bounds_warnings = []
    params_at_bounds = []
    for i, status in enumerate(bound_status):
        if status not in ('lower', 'upper'):
            continue
        params_at_bounds.append(i)
        bound_val = lower_bounds_full[i] if status == 'lower' else upper_bounds_full[i]
        bounds_warnings.append(
            f"Parameter {param_labels_indexed[i]} = {params_opt[i]:.3e} near {status} "
            f"bound {bound_val:.1e}"
        )

    # Build FitDiagnostics
    # Optimizer metadata belongs to the least_squares run that produced the
    # result, or to the one from DE's point when DE's own point is kept; with
    # no refinement at all it would be misattributed.
    ls_meta = final_ls if final_ls is not None else (from_de.result if from_de is not None else None)
    ls_nfev = (from_de.result.nfev if from_de is not None else 0) + archive_nfev
    fit_diagnostics = FitDiagnostics(
        optimizer_status=ls_meta.status if ls_meta is not None else -1,
        optimizer_message=ls_meta.message if ls_meta is not None else 'DE only (refinement failed)',
        optimizer_success=ls_meta.success if ls_meta is not None else de_result.success,
        n_function_evals=de_result.nfev + ls_nfev,
        jacobian_type=jacobian_type,
        condition_number=condition_number,
        covariance_rank=cov_result.rank if cov_result else 0,
        covariance_warning=cov_result.warning_message if cov_result else None,
        params_at_bounds=params_at_bounds,
        bounds_warnings=bounds_warnings,
        warnings=diag_warnings
    )

    # Step 5: Create FitResult
    # When cov_result is None, stderr is inf so the CI is +/-inf regardless of dof.
    fit_result = FitResult(
        circuit=circuit,
        params_opt=params_opt,
        params_stderr=params_stderr,
        fit_error_rel=fit_error_rel,
        fit_error_abs=fit_error_abs,
        quality=quality,
        condition_number=condition_number,
        is_well_conditioned=is_well_conditioned,
        cov=cov,
        diagnostics=fit_diagnostics,
        param_labels=param_labels_indexed,
        n_free_params=len(free_indices),
        bound_status=bound_status,
        _dof=cov_result.dof if cov_result is not None else 0,
        _ci_log_scale=log_scale_ci_mask(lower_bounds_full, upper_bounds_full),
        params_significance=compute_significance(circuit, frequencies, params_opt)
    )

    # Build DiffEvoDiagnostics
    de_diagnostics = DiffEvoDiagnostics(
        strategy=strategy_name,
        popsize=popsize,
        maxiter=maxiter,
        tol=tol,
        workers=workers,
        weighting=weighting,
        jacobian_type=jacobian_type,
        de_converged=de_result.success,
        de_iterations=de_result.nit,
        de_evaluations=de_result.nfev,
        de_error=de_error_rel,
        refined_error=ls_error_rel,
        de_cost=de_cost,
        refined_cost=refined_cost,
        refinement_improved=used_refinement,
        total_evaluations=de_result.nfev + ls_nfev,
        archive_checked=archive_check,
        archive_candidates=len(from_archive),
        archive_used=archive_used,
        n_fixed_params=len(fixed_param_indices),
        log_search_params=log_search_params,
        fixed_param_indices=fixed_param_indices,
        initial_guess=list(initial_guess_full),
        warnings=diag_warnings
    )

    # Create result object
    diffevo_result = DiffEvoResult(
        best_result=fit_result,
        de_result=de_result,
        de_error=de_error_rel,
        final_error=fit_error_rel,
        n_evaluations=de_diagnostics.total_evaluations,
        strategy=strategy_name,
        improvement=improvement,
        diagnostics=de_diagnostics
    )

    return diffevo_result, Z_fit


__all__ = [
    'DiffEvoResult',
    'DiffEvoDiagnostics',
    'DE_STRATEGIES',
    'fit_circuit_diffevo',
]

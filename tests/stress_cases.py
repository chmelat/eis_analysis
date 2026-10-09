"""Random, physically faithful EIS cases for the stress test (tests/stress.py).

Every case is reproducible from (family, index) alone: each purpose draws from
its own stream np.random.default_rng([family_id, index, purpose]), so adding a
draw to the circuit generator moves neither the noise nor the fit starts.
See doc/STRESS_TEST_PLAN.md for the design.
"""

from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np

from eis_analysis.cli.utils import parse_circuit_expression
from eis_analysis.fitting import Wa, Wat
from eis_analysis.fitting.bounds import classify_bound_status, generate_simple_bounds
from eis_analysis.validation.kramers_kronig import low_frequency_slope

# Fixed ids: a new family gets a new number, existing ones are never renumbered
# (hash(str) would differ between processes under PYTHONHASHSEED).
FAMILY_IDS = {'rc': 1, 'cpe': 2, 'diffusion': 3, 'blocking': 4, 'oxide': 5,
              'anomalous': 6}

# Random stream per purpose
CIRCUIT, NOISE, START, MULTISTART, DE, ORDER, DE_START = range(7)

# How a parameter scales when Z -> k*Z (power of k), by parameter label
K_POWER = {'R': 1, 'L': 1, 'σ': 1, 'R_W': 1, 'C': -1, 'Q': -1, 'G': -1,
           'n': 0, 'τ_W': 0, 'γ_W': 0}

# Each arc's tau stays this many decades inside the measured window
# [1/(2 pi f_max), 1/(2 pi f_min)], so its whole arc is visible.
TAU_MARGIN_DEC = 0.5

# Spread of all resistances of one case (Rs and the arcs), decades
R_SPREAD_DEC = 3.0

NOISE_LEVELS = (0.0, 1e-3, 1e-2, 3e-2)
CONSTANT_NOISE_SHARE = 0.2
# Constant noise is level x max|Z|, but at most min|Z| / this: on spectra
# spanning many decades it otherwise left points with |Z| 1e4x below the
# noise (oxide/177), data no instrument returns and no method can read.
# 10 keeps the weakest point at 10 % noise, 3x the top proportional level.
CONSTANT_NOISE_MIN_SNR = 10.0

# A draw that keeps failing the bound margin points to a broken family
# definition, not to bad luck.
MAX_REJECTIONS = 10_000


@dataclass
class Case:
    """One generated spectrum with its ground truth."""
    family: str
    index: int
    expression: str            # circuit with the true values
    labels: List[str]
    truth: np.ndarray
    frequencies: np.ndarray    # descending, as an instrument writes them
    Z: np.ndarray              # with noise
    noise_level: float
    noise_kind: str            # 'proportional' or 'constant'
    n_arcs: int
    min_sep: Optional[float]   # smallest tau separation of neighbouring arcs, decades
    min_frac: Optional[float]  # smallest arc R / sum of arc R
    has_L: bool
    n_rejected: int            # draws discarded for the bound margin
    Z_clean: np.ndarray        # noise-free truth
    sigma: np.ndarray          # noise standard deviation per point, on Re and Im each

    def rng(self, purpose: int) -> np.random.Generator:
        return case_rng(self.family, self.index, purpose)

    def circuit(self, params=None):
        """Fresh circuit of the true structure, at `params` (default: truth)."""
        circuit = parse_circuit_expression(self.expression)
        if params is not None:
            circuit.update_params(list(params))
        return circuit


def case_rng(family: str, index: int, purpose: int) -> np.random.Generator:
    return np.random.default_rng([FAMILY_IDS[family], index, purpose])


def k_powers(labels: List[str]) -> np.ndarray:
    return np.array([K_POWER[label] for label in labels])


def bound_violations(labels: List[str], params) -> List[str]:
    """Parameters closer to PARAMETER_BOUNDS than classify_bound_status allows.

    Uses the library's own rule (1 decade on log scale, 1 % of the range on
    linear scale), so a truth that passes here draws no bound warning.
    """
    lower, upper = generate_simple_bounds(labels)
    return [f'{label}={value:.3g} {status}'
            for label, value, lo, hi in zip(labels, params, lower, upper)
            if (status := classify_bound_status(value, lo, hi))]


def _fmt(value: float) -> str:
    return repr(float(value))


def _log_uniform(rng, lo: float, hi: float) -> float:
    return float(10 ** rng.uniform(np.log10(lo), np.log10(hi)))


def _arc_taus(rng, tau_lo: float, tau_hi: float, n_arcs: int) -> Optional[np.ndarray]:
    """log10 tau of `n_arcs` arcs, neighbours 0.3-3 decades apart, all inside
    the window shrunk by TAU_MARGIN_DEC. None if the separations do not fit."""
    lo = np.log10(tau_lo) + TAU_MARGIN_DEC
    hi = np.log10(tau_hi) - TAU_MARGIN_DEC
    # 0.3-3 decades: tight pairs on purpose (no fixed lower limit beyond
    # what a DRT can still show), up to well separated
    seps = rng.uniform(0.3, 3.0, n_arcs - 1)
    span = float(np.sum(seps))
    if span > hi - lo:
        return None
    start = rng.uniform(lo, hi - span)
    return start + np.concatenate([[0.0], np.cumsum(seps)])


def _resistances(rng, n: int) -> np.ndarray:
    """n resistances inside one R_SPREAD_DEC band; the band lies in 1e-3..1e9
    Ohm (mOhm batteries to GOhm films, 1 decade inside RESISTANCE_RANGE)."""
    base = rng.uniform(-3.0, 9.0 - R_SPREAD_DEC)
    return 10 ** rng.uniform(base, base + R_SPREAD_DEC, n)


def _separation_stats(log_taus: np.ndarray, arc_R: np.ndarray):
    min_sep = float(np.min(np.diff(log_taus))) if len(log_taus) > 1 else None
    return min_sep, float(np.min(arc_R) / np.sum(arc_R))


# --- families ---------------------------------------------------------------
# Each returns (expression, n_arcs, min_sep, min_frac, has_L) or None to redraw.
Draw = Optional[Tuple[str, int, Optional[float], Optional[float], bool]]


def _rc(rng, tau_lo, tau_hi, f_min) -> Draw:
    n_arcs = int(rng.integers(1, 5))
    log_taus = _arc_taus(rng, tau_lo, tau_hi, n_arcs)
    if log_taus is None:
        return None
    R = _resistances(rng, n_arcs + 1)
    arcs = [f'(R({_fmt(r)})|C({_fmt(10**lt / r)}))' for r, lt in zip(R[1:], log_taus)]
    expr = ' - '.join([f'R({_fmt(R[0])})'] + arcs)
    has_L = bool(rng.random() < 0.3)
    if has_L:
        # Cable / fixture inductance, 10 nH - 1 uH
        expr += f' - L({_fmt(_log_uniform(rng, 1e-8, 1e-6))})'
    return (expr, n_arcs, *_separation_stats(log_taus, R[1:]), has_L)


def _cpe_arcs(rng, log_taus, R, n_lo, n_hi) -> List[str]:
    arcs = []
    for r, lt in zip(R, log_taus):
        n = rng.uniform(n_lo, n_hi)
        # tau = (R Q)^(1/n)  ->  Q = tau^n / R
        arcs.append(f'(R({_fmt(r)})|Q({_fmt(10**(lt * n) / r)},{_fmt(n)}))')
    return arcs


def _cpe(rng, tau_lo, tau_hi, f_min) -> Draw:
    n_arcs = int(rng.integers(1, 4))
    log_taus = _arc_taus(rng, tau_lo, tau_hi, n_arcs)
    if log_taus is None:
        return None
    R = _resistances(rng, n_arcs + 1)
    # n <= 0.98 keeps the truth off the n = 1 bound; ideal n = 1 is C in `rc`
    arcs = _cpe_arcs(rng, log_taus, R[1:], 0.6, 0.98)
    expr = ' - '.join([f'R({_fmt(R[0])})'] + arcs)
    return (expr, n_arcs, *_separation_stats(log_taus, R[1:]), False)


def _randles(rng, tau_lo, tau_hi, warburg: Callable[[float], str]) -> Draw:
    """Rs - (Q | (R_ct - Warburg)); `warburg(R_ct)` draws the element."""
    log_tau = _arc_taus(rng, tau_lo, tau_hi, 1)
    Rs, Rct = _resistances(rng, 2)
    n = rng.uniform(0.6, 0.98)
    Q = 10 ** (log_tau[0] * n) / Rct
    expr = f'R({_fmt(Rs)}) - (Q({_fmt(Q)},{_fmt(n)})|(R({_fmt(Rct)}) - {warburg(Rct)}))'
    return (expr, 1, None, None, False)


def _finite_warburg(rng, kind: str, Rct: float, tau_lo, tau_hi, *extra: float) -> str:
    """R_W 0.1-10x R_ct, tau_W inside the window like an arc's tau."""
    R_W = Rct * 10 ** rng.uniform(-1.0, 1.0)
    tau_W = 10 ** _arc_taus(rng, tau_lo, tau_hi, 1)[0]
    return f'{kind}(' + ','.join(_fmt(v) for v in (R_W, tau_W, *extra)) + ')'


def _diffusion(rng, tau_lo, tau_hi, f_min) -> Draw:
    def warburg(Rct: float) -> str:
        kind = rng.choice(['W', 'Ws', 'Wo'])
        if kind != 'W':
            return _finite_warburg(rng, kind, Rct, tau_lo, tau_hi)
        # |Z_W(f_min)| = sigma*sqrt(2/omega) from 0.3x to 10x R_ct: the tail is
        # visible but does not have to dominate
        omega_min = 2 * np.pi * f_min
        sigma = Rct * 10 ** rng.uniform(-0.5, 1.0) / np.sqrt(2 / omega_min)
        return f'W({_fmt(sigma)})'
    return _randles(rng, tau_lo, tau_hi, warburg)


def _blocking(rng, tau_lo, tau_hi, f_min) -> Draw:
    log_tau = _arc_taus(rng, tau_lo, tau_hi, 1)
    Rs, R = _resistances(rng, 2)
    arc = _cpe_arcs(rng, log_tau, [R], 0.6, 0.98)[0]
    omega_min = 2 * np.pi * f_min
    # |Z_block(f_min)| = 1x-100x R: the capacitive end rises inside the window
    Z_block = R * 10 ** rng.uniform(0.0, 2.0)
    if rng.random() < 0.5:
        block = f'C({_fmt(1 / (omega_min * Z_block))})'
    else:
        n = rng.uniform(0.85, 0.98)
        block = f'Q({_fmt(1 / (omega_min**n * Z_block))},{_fmt(n)})'
    return (f'R({_fmt(Rs)}) - {arc} - {block}', 1, None, None, False)


def _oxide(rng, tau_lo, tau_hi, f_min) -> Draw:
    n_arcs = int(rng.integers(1, 3))
    log_taus = _arc_taus(rng, tau_lo, tau_hi, n_arcs)
    if log_taus is None:
        return None
    # Zr oxide films: electrolyte 1-100 Ohm, film 1 MOhm - 1 GOhm; Q follows
    # from tau and R (the bound margin rejects Q below 1e-11)
    Rs = _log_uniform(rng, 1.0, 100.0)
    R = 10 ** rng.uniform(6.0, 9.0, n_arcs)
    arcs = _cpe_arcs(rng, log_taus, R, 0.8, 0.98)
    expr = ' - '.join([f'R({_fmt(Rs)})'] + arcs)
    return (expr, n_arcs, *_separation_stats(log_taus, R), False)


def _unit_admittance(kind: str, omega: float, tau_W: float, gamma: float) -> complex:
    """Admittance of Wa or Wat at R_W = 1, from the library's own element."""
    element = Wa if kind == 'Wa' else Wat
    params = [1.0, tau_W, gamma]
    return complex(1 / element(*params).impedance(np.array([omega / (2 * np.pi)]), params)[0])


def _anomalous(rng, tau_lo, tau_hi, f_min) -> Draw:
    """Wa or Wat: half Randles, half the ZrO2 layer model they were made for."""
    kind = str(rng.choice(['Wa', 'Wat']))
    if rng.random() < 0.5:
        # gamma <= 0.95 keeps the truth off gamma = 1 (exactly Wo/Ws, the
        # bound); 0.5 is below every gamma the ZrO2 fits gave (0.60-0.72)
        gamma = rng.uniform(0.5, 0.95)
        return _randles(rng, tau_lo, tau_hi,
                        lambda Rct: _finite_warburg(rng, kind, Rct, tau_lo, tau_hi, gamma))

    # ZrO2 layer (M136 fits, CHANGELOG `Wa`), drawn by its crossovers so that
    # every element shows in the window: 1/omega where Q meets C, where the
    # Warburg meets Q, and tau_W, ascending
    log_taus = _arc_taus(rng, tau_lo, tau_hi, 3)
    if log_taus is None:
        return None
    tau_QC, tau_WQ, tau_W = 10 ** log_taus
    # Electrolyte 1-100 Ohm, as in `oxide` (the fits: 1.2 Ohm)
    Rs = _log_uniform(rng, 1.0, 100.0)
    # Film capacitance: 0.1-10 um of ZrO2 on mm^2 - cm^2
    C = _log_uniform(rng, 1e-9, 1e-7)
    # The fits gave n 0.73-0.79, widened by ~0.15 to either side. gamma < n,
    # as in every fit: the Warburg's low-frequency CPE (exponent gamma) then
    # decays slower than Q and owns the low-frequency end, where the
    # blocking shows; with gamma > n, Q would take it back below 1/tau_W.
    n = rng.uniform(0.6, 0.9)
    gamma = rng.uniform(0.5, n)
    Q = C * tau_QC ** (n - 1)                                   # Q w^n = C w
    # |Y_W(1/tau_WQ)| = Q w^n, exactly, not by the high-frequency asymptote:
    # tau_W can be only 0.3 decade above tau_WQ
    R_W = abs(_unit_admittance(kind, 1 / tau_WQ, tau_W, gamma)) * tau_WQ ** n / Q
    element = f'{kind}({_fmt(R_W)},{_fmt(tau_W)},{_fmt(gamma)})'
    branches = [element, f'Q({_fmt(Q)},{_fmt(n)})', f'C({_fmt(C)})']
    if kind == 'Wa':
        # A DC conductance that bends the low-frequency end without hiding
        # the blocking: 0.1-1x the Warburg admittance at f_min. Wat has its
        # own DC path, the (Wat|Q|C) model.
        Y_min = abs(_unit_admittance(kind, 2 * np.pi * f_min, tau_W, gamma)) / R_W
        branches.insert(0, f'G({_fmt(Y_min * 10 ** rng.uniform(-1.0, 0.0))})')
    return (f'R({_fmt(Rs)}) - (' + '|'.join(branches) + ')', 1, None, None, False)


FAMILIES: Dict[str, Callable[..., Draw]] = {
    'rc': _rc, 'cpe': _cpe, 'diffusion': _diffusion,
    'blocking': _blocking, 'oxide': _oxide, 'anomalous': _anomalous,
}


def generate_case(family: str, index: int) -> Case:
    rng = case_rng(family, index, CIRCUIT)
    for n_rejected in range(MAX_REJECTIONS):
        f_max = 10 ** rng.uniform(3, 7)
        f_min = 10 ** rng.uniform(-3, 1)
        per_decade = int(rng.integers(5, 16))
        n_points = int(round(np.log10(f_max / f_min) * per_decade)) + 1
        frequencies = np.logspace(np.log10(f_max), np.log10(f_min), n_points)
        tau_lo, tau_hi = 1 / (2 * np.pi * f_max), 1 / (2 * np.pi * f_min)

        draw = FAMILIES[family](rng, tau_lo, tau_hi, f_min)
        if draw is None:
            continue
        expression, n_arcs, min_sep, min_frac, has_L = draw
        circuit = parse_circuit_expression(expression)
        labels = circuit.get_param_labels()
        truth = np.array(circuit.get_all_params(), dtype=float)
        if not bound_violations(labels, truth):
            break
    else:
        raise RuntimeError(f'{family}/{index}: no draw within the bound margin')

    Z_clean = circuit.impedance(frequencies, list(truth))

    noise_rng = case_rng(family, index, NOISE)
    noise_level = float(noise_rng.choice(NOISE_LEVELS))
    constant = noise_level > 0 and noise_rng.random() < CONSTANT_NOISE_SHARE
    Z_mag = np.abs(Z_clean)
    sigma = (min(noise_level * np.max(Z_mag), np.min(Z_mag) / CONSTANT_NOISE_MIN_SNR)
             if constant else noise_level * Z_mag)
    Z = (Z_clean + noise_rng.normal(0.0, 1.0, n_points) * sigma
         + 1j * noise_rng.normal(0.0, 1.0, n_points) * sigma)

    return Case(family, index, expression, labels, truth, frequencies, Z,
                noise_level, 'constant' if constant else 'proportional',
                n_arcs, min_sep, min_frac, has_L, n_rejected, Z_clean,
                np.broadcast_to(sigma, Z_clean.shape).astype(float))


def capacitive_end(case: Case) -> bool:
    """True when the truth's low-frequency end is capacitive (-Z'' still
    rising): the rule by which the CLI suggests Lin-KK's series C
    (cli/handlers/validation.py), applied to the noise-free spectrum. A DC
    path far below f_min (a small G) does not count: in the window the end
    is capacitive all the same."""
    return bool(low_frequency_slope(case.frequencies, case.Z_clean) < 0)


# Parameters on a linear scale, capped at 1: a start factor would throw them
# out of their bounds
LINEAR_LABELS = ('n', 'γ_W')

# DE start: scale parameters up to this many decades from the truth, as a
# rough guess for an unfamiliar sample can be
DE_START_DEC = 2.0


def de_start(case: Case) -> np.ndarray:
    """Start of the DE fit: truth x 10^U(-2, 2), the exponents n and gamma
    uniform within their bounds, clipped into the bounds.

    DE takes the start as one member of its population, so a start near the
    truth (fit_start) would test little more than its polish.
    """
    rng = case.rng(DE_START)
    lower, upper = (np.asarray(b, dtype=float) for b in generate_simple_bounds(case.labels))
    start = case.truth * 10 ** rng.uniform(-DE_START_DEC, DE_START_DEC, len(case.truth))
    linear = np.isin(case.labels, LINEAR_LABELS)
    start[linear] = rng.uniform(lower[linear], upper[linear])
    return np.clip(start, lower, upper)


def fit_start(case: Case) -> np.ndarray:
    """Start of the LM fit: truth x U(0.3, 3), clipped into the bounds.

    The exponents n and gamma are linear and capped at 1, where a factor
    0.3-3 would throw them out of the bounds; they get x U(0.9, 1.1) instead.
    """
    rng = case.rng(START)
    factors = 10 ** rng.uniform(np.log10(0.3), np.log10(3.0), len(case.truth))
    linear = np.isin(case.labels, LINEAR_LABELS)
    factors[linear] = rng.uniform(0.9, 1.1, int(linear.sum()))
    lower, upper = generate_simple_bounds(case.labels)
    return np.clip(case.truth * factors, lower, upper)

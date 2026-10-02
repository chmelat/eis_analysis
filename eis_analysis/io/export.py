"""
Export of a circuit fit result to machine-readable files.

The console output of a fit is for reading; this is for the next program:
a JSON record of the parameters and statistics, and a CSV of the data with
the fitted curve that opens in any spreadsheet and loads back into eis.
"""

import json
from typing import Any, Dict, Optional, Tuple

import numpy as np
from numpy.typing import NDArray

from ..fitting import FitResult
from ..fitting.diagnostics import compute_information_criteria
from ..version import __version__

# Raise on any change a reader of an older record would misread
FIT_EXPORT_FORMAT_VERSION = 1

# Headers chosen so load_csv_data reads the data columns and skips the fit:
# 'fit' is one of its words that rule a column out.
_CSV_HEADER = "freq_Hz,Z_real_Ohm,Z_imag_Ohm,Z_fit_real_Ohm,Z_fit_imag_Ohm"


def _finite(x: float) -> Optional[float]:
    """The value as float, or None (JSON null) where it is inf or NaN."""
    return float(x) if np.isfinite(x) else None


def fit_result_record(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    Z_fit: NDArray[np.complex128],
    result: FitResult,
    weighting: str,
    context: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    """
    JSON-ready record of a fit result.

    Parameters
    ----------
    frequencies : ndarray of float
        Frequencies the fit ran on [Hz]
    Z : ndarray of complex
        Impedance the fit ran on [Ohm]
    Z_fit : ndarray of complex
        Fitted impedance at `frequencies` [Ohm]
    result : FitResult
        The fit
    weighting : str
        Weighting the fit used; AIC/BIC are computed with it
    context : dict, optional
        How the fit came about (input file, optimizer, ...). Stored as given.

    Returns
    -------
    record : dict
        Serializable with ``json.dumps(..., allow_nan=False)`` as long as
        `context` is

    Notes
    -----
    Uncertainties are None where they carry no meaning: for a parameter at a
    bound or fixed (the Jacobian CI does not apply there) and where the
    standard error is not finite. Every inf and NaN is None, so the record
    is strict JSON.
    """
    params = np.asarray(result.params_opt, dtype=float)
    n = len(params)
    labels = result.param_labels
    if labels is None:
        # The linear Voigt chain sets none; index the circuit's own labels
        # the way the nonlinear fitters do (R0, R1, τ0, ...)
        raw = result.circuit.get_param_labels()
        labels = [f"{label}{raw[:i].count(label)}" for i, label in enumerate(raw)]
    status = result.bound_status or [''] * n
    significance = result.params_significance
    ci_low, ci_high = result.params_ci_95

    parameters = []
    for i in range(n):
        has_ci = status[i] == '' and np.isfinite(result.params_stderr[i])
        parameters.append({
            'label': labels[i],
            'value': _finite(params[i]),
            'stderr': _finite(result.params_stderr[i]) if has_ci else None,
            'ci95': [_finite(ci_low[i]), _finite(ci_high[i])] if has_ci else None,
            'significance': None if significance is None else _finite(significance[i]),
            'status': status[i],
        })

    # The linear Voigt chain leaves n_free_params at 0: there is then no
    # free-parameter count, so no dof and no AIC/BIC (which refuse k = 0).
    # dof is the one the CIs above were computed with, not 2N - k again:
    # the fit clamps it to >= 1, and DE leaves 0 (no CIs) when the
    # covariance fails - None then, not a dof of 0.
    k = int(result.n_free_params)
    dof = aic = bic = None
    if k > 0:
        dof = int(result._dof) if result._dof > 0 else None
        _, aic_value, bic_value = compute_information_criteria(Z, Z_fit, weighting, k)
        aic, bic = _finite(aic_value), _finite(bic_value)

    return {
        'format': 'eis_analysis.fit',
        'format_version': FIT_EXPORT_FORMAT_VERSION,
        'eis_analysis_version': __version__,
        # Parses back with parse_circuit_expression: the fitted circuit, free
        # values rounded for display. With the full-precision values below as
        # params it is the fitted model.
        'circuit': str(result.circuit),
        'weighting': weighting,
        'context': context or {},
        'n_points': len(frequencies),
        'f_range_Hz': [_finite(float(np.min(frequencies))), _finite(float(np.max(frequencies)))],
        'parameters': parameters,
        'metrics': {
            'fit_error_rel_pct': _finite(result.fit_error_rel),
            'fit_error_abs_ohm': _finite(result.fit_error_abs),
            'quality': result.quality,
            'condition_number': _finite(result.condition_number),
            'well_conditioned': bool(result.is_well_conditioned),
            'n_free_params': k if k > 0 else None,
            'dof': dof,
            'aic': aic,
            'bic': bic,
        },
        'covariance': (None if result.cov is None
                       else [[_finite(c) for c in row] for row in result.cov]),
        'warnings': result.all_warnings,
    }


def save_fit_result(
    prefix: str,
    record: Dict[str, Any],
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    Z_fit: NDArray[np.complex128]
) -> Tuple[str, str]:
    """
    Write a fit record to ``{prefix}.json`` and the curves to ``{prefix}.csv``.

    Parameters
    ----------
    prefix : str
        Path without extension
    record : dict
        From fit_result_record
    frequencies : ndarray of float
        Frequencies the fit ran on [Hz]
    Z : ndarray of complex
        Impedance the fit ran on [Ohm]
    Z_fit : ndarray of complex
        Fitted impedance [Ohm]

    Returns
    -------
    json_path, csv_path : str
        Paths of the written files

    Raises
    ------
    TypeError, ValueError
        If the record does not serialize (a context value that is not JSON,
        or a non-finite number); nothing is written then
    OSError
        If a file cannot be written
    """
    # Serialize before touching the disk: a failure halfway through
    # json.dump would leave a truncated file for the next reader to choke on.
    # allow_nan=False keeps the output strict JSON.
    text = json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + '\n'

    json_path, csv_path = f"{prefix}.json", f"{prefix}.csv"
    # The CSV first: if it cannot be written, the JSON of an earlier run is
    # not replaced either. A failure writing the JSON itself still leaves a
    # new CSV beside an old or partial JSON; the caller is told by the raise.
    # %.17g: enough digits for every float64 to read back exactly
    np.savetxt(csv_path,
               np.column_stack([frequencies, Z.real, Z.imag, Z_fit.real, Z_fit.imag]),
               delimiter=',', fmt='%.17g', header=_CSV_HEADER, comments='')
    with open(json_path, 'w', encoding='utf-8') as f:
        f.write(text)
    return json_path, csv_path

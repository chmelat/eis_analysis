"""Tests for the fit result export (eis_analysis.io.export) and its --save hook."""

import json
import pathlib
import sys

import numpy as np
import pytest

from eis_analysis.cli.handlers.fitting import run_circuit_fitting
from eis_analysis.cli.parser import parse_arguments
from eis_analysis.cli.utils import parse_circuit_expression
from eis_analysis.fitting import FitResult, fit_equivalent_circuit
from eis_analysis.io import fit_result_record, load_csv_data, save_fit_result

EXPR = 'R(10) - (R("1000") | Q(1e-5, 0.9))'


def _cli_args(monkeypatch, *argv):
    """Namespace from the real parser, so every option has its real default."""
    monkeypatch.setattr(sys, 'argv', ['eis', 'x.csv', '--no-show', *argv])
    return parse_arguments()


def _strict(record):
    """Round-trip through strict JSON: fails on NaN/inf and non-JSON types."""
    return json.loads(json.dumps(record, allow_nan=False))


def _read_json(path):
    def reject(token):
        raise ValueError(f"non-standard JSON constant {token}")
    with open(path, encoding='utf-8') as f:
        return json.load(f, parse_constant=reject)


@pytest.fixture
def spectrum():
    circuit = parse_circuit_expression(EXPR)
    f = np.logspace(-1, 5, 40)
    Z = circuit.impedance(f, list(circuit.get_all_params()))
    Z = Z * (1 + 0.001 * np.random.default_rng(0).normal(size=f.size))
    return f, Z


@pytest.fixture
def fit(spectrum):
    f, Z = spectrum
    result, Z_fit = fit_equivalent_circuit(f, Z, parse_circuit_expression(EXPR))
    return result, Z_fit


@pytest.fixture
def record(spectrum, fit):
    f, Z = spectrum
    result, Z_fit = fit
    return _strict(fit_result_record(f, Z, Z_fit, result, 'modulus',
                                     {'input': 'x.DTA'}))


def test_values_are_exact_and_context_kept(fit, record):
    result, _ = fit
    assert [p['value'] for p in record['parameters']] == list(result.params_opt)
    assert record['circuit'] == str(result.circuit)
    assert '"1000.0"' in record['circuit']  # the fixed value, exact
    assert record['context'] == {'input': 'x.DTA'}
    metrics = record['metrics']
    assert metrics['n_free_params'] == 3 and metrics['dof'] == 2 * 40 - 3
    assert metrics['dof'] == result._dof
    assert metrics['aic'] is not None and metrics['bic'] is not None
    assert np.array(record['covariance']).shape == (4, 4)


def test_fixed_parameter_has_no_uncertainty(record):
    fixed = record['parameters'][1]
    assert fixed['status'] == 'fixed'
    assert fixed['stderr'] is None and fixed['ci95'] is None
    free = record['parameters'][0]
    assert free['status'] == '' and free['stderr'] > 0
    assert free['ci95'][0] < free['value'] < free['ci95'][1]


def test_non_finite_values_become_null(spectrum):
    # As the linear Voigt chain leaves it: stderr unknown, n_free_params unset
    f, Z = spectrum
    circuit = parse_circuit_expression(EXPR)
    params = np.array(circuit.get_all_params())
    result = FitResult(circuit=circuit, params_opt=params,
                       params_stderr=np.full(len(params), np.inf),
                       fit_error_rel=0.1, condition_number=np.nan,
                       is_well_conditioned=False, _dof=76)
    Z_fit = circuit.impedance(f, list(params))
    record = _strict(fit_result_record(f, Z, Z_fit, result, 'modulus'))
    assert all(p['stderr'] is None and p['ci95'] is None for p in record['parameters'])
    metrics = record['metrics']
    assert metrics['condition_number'] is None and metrics['well_conditioned'] is False
    assert all(metrics[k] is None for k in ('n_free_params', 'dof', 'aic', 'bic'))
    # No labels on the result: the circuit's own, indexed like the fitters'
    assert [p['label'] for p in record['parameters']] == ['R0', 'R1', 'Q0', 'n0']
    assert record['covariance'] is None


def test_non_finite_parameter_value_becomes_null(spectrum, fit):
    f, Z = spectrum
    result, Z_fit = fit
    result.params_opt = result.params_opt.copy()
    result.params_opt[1] = np.inf
    record = _strict(fit_result_record(f, Z, Z_fit, result, 'modulus'))
    assert record['parameters'][1]['value'] is None


def test_dof_is_the_one_the_fit_used(spectrum, fit):
    # The fit clamps dof to >= 1; the record must not recompute 2N - k
    f, Z = spectrum
    result, Z_fit = fit
    result._dof = 1
    record = fit_result_record(f, Z, Z_fit, result, 'modulus')
    assert record['metrics']['dof'] == 1


def test_dof_without_covariance_is_null(spectrum, fit):
    # DE leaves dof 0 when the covariance fails: no CIs, and no dof either
    f, Z = spectrum
    result, Z_fit = fit
    result._dof = 0
    assert fit_result_record(f, Z, Z_fit, result, 'modulus')['metrics']['dof'] is None


def test_singular_fit_is_reported_ill_conditioned(spectrum, fit):
    # inf is a singular Jacobian - definitely not well conditioned, unlike
    # NaN, which means not computed
    f, Z = spectrum
    result, Z_fit = fit
    result.condition_number, result.is_well_conditioned = np.inf, False
    metrics = _strict(fit_result_record(f, Z, Z_fit, result, 'modulus'))['metrics']
    assert metrics['condition_number'] is None and metrics['well_conditioned'] is False


def test_files_load_back(spectrum, fit, tmp_path):
    f, Z = spectrum
    result, Z_fit = fit
    record = fit_result_record(f, Z, Z_fit, result, 'modulus')
    json_path, csv_path = save_fit_result(str(tmp_path / 'run_fit'), record, f, Z, Z_fit)

    assert _read_json(json_path) == _strict(record)

    # The data columns load back as input; the fit columns are skipped
    loaded = load_csv_data(csv_path)
    order = np.argsort(loaded.frequencies)  # the loader may reorder
    np.testing.assert_array_equal(loaded.frequencies[order], np.sort(f))
    np.testing.assert_allclose(loaded.Z[order], Z[np.argsort(f)], rtol=1e-15)

    # The fit columns are the model rebuilt from the record
    values = [p['value'] for p in record['parameters']]
    Z_rebuilt = parse_circuit_expression(record['circuit']).impedance(f, values)
    data = np.loadtxt(csv_path, delimiter=',', skiprows=1)
    np.testing.assert_array_equal(data[:, 3] + 1j * data[:, 4], Z_rebuilt)


def test_unserializable_record_writes_nothing(spectrum, fit, tmp_path):
    f, Z = spectrum
    result, Z_fit = fit
    record = fit_result_record(f, Z, Z_fit, result, 'modulus',
                               {'input': pathlib.Path('x.DTA')})
    with pytest.raises(TypeError):
        save_fit_result(str(tmp_path / 'run_fit'), record, f, Z, Z_fit)
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize('circuits, suffixes', [
    ([EXPR], ['fit']),
    ([EXPR, 'R(10) - (R(1000) | C(1e-5))'], ['fit_1', 'fit_2']),
])
def test_save_writes_one_record_per_fit(spectrum, tmp_path, monkeypatch,
                                        circuits, suffixes):
    f, Z = spectrum
    prefix = str(tmp_path / 'run')
    circuit_args = [a for c in circuits for a in ('--circuit', c)]
    run_circuit_fitting(f, Z, _cli_args(monkeypatch, *circuit_args,
                                        '--optimizer', 'single', '--save', prefix))
    for suffix in suffixes:
        record = _read_json(f"{prefix}_{suffix}.json")
        assert record['context']['optimizer'] == 'single'
        assert (tmp_path / f"run_{suffix}.csv").stat().st_size > 0


def test_voigt_chain_record_rebuilds_the_model(spectrum, tmp_path, monkeypatch):
    f, Z = spectrum
    prefix = str(tmp_path / 'run')
    run_circuit_fitting(f, Z, _cli_args(monkeypatch, '--voigt-chain', '--save', prefix))
    record = _read_json(f"{prefix}_fit.json")
    assert record['context']['optimizer'] == 'linear'
    # The linear fit computes no conditioning; it must not claim one
    assert record['metrics']['condition_number'] is None
    assert record['metrics']['well_conditioned'] is False
    assert {p['label'][:1] for p in record['parameters']} >= {'R', 'τ'}
    values = [p['value'] for p in record['parameters']]
    Z_rebuilt = parse_circuit_expression(record['circuit']).impedance(f, values)
    data = np.loadtxt(f"{prefix}_fit.csv", delimiter=',', skiprows=1)
    np.testing.assert_allclose(data[:, 3] + 1j * data[:, 4], Z_rebuilt, rtol=1e-12)


def test_export_error_keeps_the_fit(spectrum, tmp_path, monkeypatch):
    # A caller building the namespace by hand may pass a Path
    f, Z = spectrum
    args = _cli_args(monkeypatch, '--circuit', EXPR, '--optimizer', 'single',
                     '--save', str(tmp_path / 'run'))
    args.input = pathlib.Path('x.csv')
    result, _ = run_circuit_fitting(f, Z, args)
    assert result is not None
    assert not any(p.suffix in ('.json', '.csv') for p in tmp_path.iterdir())


def test_namespace_without_input_options_keeps_the_fit(spectrum, tmp_path, monkeypatch):
    # A caller building the namespace by hand need not set f_min/f_max/input
    f, Z = spectrum
    args = _cli_args(monkeypatch, '--circuit', EXPR, '--optimizer', 'single',
                     '--save', str(tmp_path / 'run'))
    del args.f_min, args.f_max
    result, _ = run_circuit_fitting(f, Z, args)
    assert result is not None
    assert _read_json(tmp_path / 'run_fit.json')['context']['f_min'] is None


def test_unexpected_export_error_keeps_the_fit(spectrum, tmp_path, monkeypatch):
    import eis_analysis.cli.handlers.fitting as handler

    def broken(*args, **kwargs):
        raise IndexError("boom")
    monkeypatch.setattr(handler, 'fit_result_record', broken)
    f, Z = spectrum
    args = _cli_args(monkeypatch, '--voigt-chain', '--save', str(tmp_path / 'run'))
    result, _ = run_circuit_fitting(f, Z, args)
    assert result is not None


def test_no_save_writes_nothing(spectrum, tmp_path, monkeypatch):
    f, Z = spectrum
    monkeypatch.chdir(tmp_path)
    run_circuit_fitting(f, Z, _cli_args(monkeypatch, '--circuit', EXPR,
                                        '--optimizer', 'single'))
    assert not any(p.suffix in ('.json', '.csv') for p in tmp_path.iterdir())

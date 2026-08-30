"""
Tests for CLI argument parsing (cli/parser.py).

Regression tests for audit finding P1 (2026-07-08): --multistart N without
--optimizer multistart used to be silently ignored (the fit ran DE), and
the help text claimed "0 = disabled" while the handler mapped 0 to a hidden
default of 16 restarts.
"""

import sys

import pytest

from eis_analysis.cli.parser import parse_arguments


def parse(monkeypatch, *argv):
    monkeypatch.setattr(sys, 'argv', ['eis', *argv])
    return parse_arguments()


def test_default_optimizer_is_de(monkeypatch):
    args = parse(monkeypatch)
    assert args.optimizer == 'de'
    assert args.multistart is None


def test_multistart_implies_multistart_optimizer(monkeypatch):
    args = parse(monkeypatch, '--multistart', '20')
    assert args.optimizer == 'multistart'
    assert args.multistart == 20


def test_explicit_multistart_optimizer_without_n(monkeypatch):
    # N defaults in the handler (16); parser leaves None
    args = parse(monkeypatch, '--optimizer', 'multistart')
    assert args.optimizer == 'multistart'
    assert args.multistart is None


def test_multistart_conflicts_with_other_optimizer(monkeypatch):
    for optimizer in ('de', 'single'):
        with pytest.raises(SystemExit):
            parse(monkeypatch, '--optimizer', optimizer, '--multistart', '20')


def test_multistart_must_be_positive(monkeypatch):
    for n in ('0', '-5'):
        with pytest.raises(SystemExit):
            parse(monkeypatch, '--multistart', n)


def test_lambda_probe_flag(monkeypatch):
    assert parse(monkeypatch).lambda_probe is False
    assert parse(monkeypatch, '--lambda-probe').lambda_probe is True


def test_oxide_defaults_are_none(monkeypatch):
    """None defaults let the handler tell "not given" from an explicit value."""
    args = parse(monkeypatch)
    assert args.thickness is None
    assert args.epsilon_r is None


def test_thickness_and_epsilon_r_parse_together(monkeypatch):
    """Conflicting flags are not a parse error - the handler warns instead."""
    args = parse(monkeypatch, '--analyze-oxide', '--thickness', '20', '--epsilon-r', '9')
    assert args.thickness == 20.0
    assert args.epsilon_r == 9.0


# --- --circuit accepts several candidates (action='append') -----------------

def test_no_circuit_stays_none(monkeypatch):
    """The auto-skip in run_circuit_fitting keys on None; keep that contract."""
    args = parse(monkeypatch)
    assert args.circuit is None


def test_single_circuit_is_a_one_element_list(monkeypatch):
    args = parse(monkeypatch, '-c', 'R(10)-(R(200)|C(2e-5))')
    assert args.circuit == ['R(10)-(R(200)|C(2e-5))']


def test_repeated_circuit_collects_all_candidates(monkeypatch):
    """Before action='append' the last --circuit silently won."""
    args = parse(monkeypatch,
                 '-c', 'R(10)-(R(200)|C(2e-5))',
                 '--circuit', 'R(10)-(R(200)|Q(2e-5,0.85))')
    assert args.circuit == ['R(10)-(R(200)|C(2e-5))',
                            'R(10)-(R(200)|Q(2e-5,0.85))']

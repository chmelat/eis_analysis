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

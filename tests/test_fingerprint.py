"""Compare logic of tests/fingerprint.py on records built in memory (no analyses run)."""

import pickle
from dataclasses import dataclass, field

import numpy as np

from tests.fingerprint import _flags, classify, compare_files

META = {'commit': 'x', 'version': '0', 'python': '3', 'numpy': '2', 'scipy': '1', 'n': 1,
        'package': 'p'}


def _entry(value, kind='num', atol=None):
    return {'value': value, 'power': 0, 'kind': kind, 'finite': True, 'atol': atol}


def _sig(shift=0.0, flag='a'):
    return {'gamma': _entry(np.array([1.0, 2.0, 0.0]) * (1 + shift)),
            'params': _entry(np.array([1e-6, 3e4]) * (1 + shift), 'rel'),
            'flags': _entry([flag], 'exact')}


def test_classes_and_exit_code(tmp_path):
    assert classify(_sig(), _sig()) == ('identical', '')
    cls, detail = classify(_sig(), _sig(1e-12))
    assert cls == 'within' and 'rel 1.0e-12' in detail
    assert classify(_sig(), _sig(1e-3))[0] == 'beyond'
    assert classify(_sig(), _sig(flag='b'))[0] == 'beyond'
    error = ('exception', 'ValueError', 'bad')
    assert classify(error, error) == ('identical', '')
    assert classify(error, ('exception', 'ValueError', 'worse'))[0] == 'exception'
    assert classify(_sig(), error)[0] == 'exception'
    assert classify(None, _sig())[0] == 'missing'

    # A stored atol is the tolerance (C moved by half of itself, within its atol),
    # and a 'rel' change is relative to the element, not to the largest one
    ref = {'params': _entry(np.array([1e-6, 3e4]), 'rel', atol=np.array([1e-6, 0.0]))}
    got = {'params': _entry(np.array([1.5e-6, 3e4]), 'rel', atol=np.array([1e-6, 0.0]))}
    assert classify(ref, got) == ('within', 'params: rel 5.0e-01')

    def write(name, cases):
        path = tmp_path / name
        path.write_bytes(pickle.dumps({'meta': META, 'cases': cases}))
        return path

    same = write('a.pkl', {('rc', 0): {'drt': _sig()}})
    assert compare_files(same, write('b.pkl', {('rc', 0): {'drt': _sig()}})) == 0
    assert compare_files(same, write('c.pkl', {('rc', 0): {'drt': _sig(1e-12)}})) == 1
    assert compare_files(same, write('e.pkl', {('rc', 1): {'drt': _sig()}})) == 1
    # Another input: reported as such, its analyses are not compared
    one = write('f.pkl', {('rc', 0): {'_input': 'a', 'drt': _sig()}})
    assert compare_files(one, write('g.pkl', {('rc', 0): {'_input': 'a', 'drt': _sig()}})) == 0
    assert compare_files(one, write('h.pkl', {('rc', 0): {'_input': 'b', 'drt': _sig()}})) == 1


@dataclass
class _Result:
    noise: float
    method: str
    gamma: np.ndarray
    warnings: list = field(default_factory=list)
    scores: list = field(default_factory=list)
    _dof: int = 3

    @property
    def is_valid(self):
        return np.float64(self.noise) < 1   # np.bool_, as numpy comparisons give


def test_flags_take_flag_and_text_fields_only():
    """Numbers stay with the signatures, whose entries carry their tolerances."""
    flags = _flags(_Result(0.5, 'fit', np.ones(3), ['low f'], [1.0, 2.0]))
    assert set(flags) == {'raw.method', 'raw.warnings', 'raw.is_valid'}
    assert flags['raw.is_valid']['value'] is True
    changed = _flags(_Result(5.0, 'fit', np.ones(3), ['low f'], [1.0, 2.0]))
    assert classify(flags, changed)[0] == 'beyond'

#!/usr/bin/env python3
"""Measured spectra: a batch of 54 nominally identical Li-ion cells.

The data (`tests/data/nmc811_zenodo_15422339`, CC BY 4.0, source and
description in its README) are 54 cylindrical 21700 NMC811 cells at 30 %
SoC, 44 points from 10 kHz to 50 mHz. The true processes are not known, so
the checks are invariants of a batch, not a ground truth: cells this alike
must get the same analysis, and a method that splits them (KAIRIS found 1,
2 or 3 DRT peaks as its L-curve lambda jumped 9 decades between cells,
doc/KAIRIS_COMPARISON.md, untracked) is unstable, not the cells.

The spectra are hard at both ends: an inductive high-frequency end where
Re(Z) rises with f (skin effect, not representable by the DRT), and a
diffusion tail that continues below 50 mHz. That is where a change in the
lambda choice, the inductance column or the R_inf estimate shows up.

About 15 s for the batch, so marked `slow`.
"""

from pathlib import Path

import numpy as np
import pytest

from eis_analysis import calculate_drt, kramers_kronig_validation, load_csv_data
from eis_analysis.validation.outliers import find_outliers
from eis_analysis.validation.zhit import zhit_validation

pytestmark = pytest.mark.slow

DATA = Path(__file__).parent / 'data' / 'nmc811_zenodo_15422339'
N_CELLS = 54
N_POINTS = 44
# Measured 0.99986-1.0074 (2026-11-02). Re(Z) at f_max lies 4-5 % above the
# minimum (skin effect), so an R_inf taken from the top point fails this.
RINF_TO_MIN_RE_TOL = 0.02
# Measured 89.8-96.0 nH; the band is the cells' lead and winding inductance
L_SERIES_RANGE = (80e-9, 120e-9)
# Measured max/min 1.069 across the batch
L_SERIES_MAX_SPREAD = 1.15


@pytest.fixture(scope='module')
def batch():
    files = sorted(DATA.glob('ID*.csv'))
    out = []
    for path in files:
        data = load_csv_data(str(path))
        out.append((path.name, data,
                    calculate_drt(data.frequencies, data.Z),
                    kramers_kronig_validation(data.frequencies, data.Z)))
    return out


@pytest.fixture(scope='module')
def zhit(batch):
    return [zhit_validation(data.frequencies, data.Z) for _, data, _, _ in batch]


def test_all_cells_load(batch):
    assert len(batch) == N_CELLS
    for name, data, _, _ in batch:
        assert len(data.frequencies) == N_POINTS, name
        assert not data.warnings, (name, data.warnings)
        # Imag_Ohm is signed: the high-frequency end is inductive
        assert data.Z.imag[np.argmax(data.frequencies)] > 0, name


def test_kk_valid_for_all_cells(batch):
    invalid = [name for name, _, _, kk in batch if not kk.is_valid]
    assert not invalid


def test_same_peak_count_across_batch(batch):
    counts = {name: drt.diagnostics.n_peaks for name, _, drt, _ in batch}
    assert len(set(counts.values())) == 1, counts


def test_inductance_found_and_consistent(batch):
    L = np.array([drt.L_series for _, _, drt, _ in batch])
    assert np.all((L >= L_SERIES_RANGE[0]) & (L <= L_SERIES_RANGE[1])), L
    assert L.max() / L.min() < L_SERIES_MAX_SPREAD


def test_rinf_tracks_minimum_of_real_part(batch):
    for name, data, drt, _ in batch:
        ratio = drt.R_inf / data.Z.real.min()
        assert abs(ratio - 1) < RINF_TO_MIN_RE_TOL, (name, ratio)


def test_parallel_rl_found_in_every_cell(batch):
    # Re(Z) rises ~5 % toward 10 kHz in every cell (skin effect): the R_inf
    # check reports an R|L, with R_L and L alike across the batch
    found = [drt.diagnostics.rinf.estimate.parallel_RL for _, _, drt, _ in batch]
    assert all(rl is not None for rl in found)
    R_L = np.array([rl.R_L for rl in found])
    assert R_L.max() / R_L.min() < 1.5


# The spectra pass KK; on clean measured data the validations must not raise
# false alarms. Measured 2026-11-02: Z-HIT mean |res_mag| 0.17-0.48 %, max
# 0.84-2.3 %; Lin-KK max point residual 1.7-2.7 %; no point flagged.

def test_zhit_valid_for_all_cells(zhit):
    assert all(z.is_valid for z in zhit)


def test_no_point_flagged_as_outlier(batch, zhit):
    for (name, data, _, kk), z in zip(batch, zhit):
        report = find_outliers(data.frequencies, kk, z)
        assert not report.points, (name, [p.frequency for p in report.points])
        assert not report.skipped, (name, report.skipped)

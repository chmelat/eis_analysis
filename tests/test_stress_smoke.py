"""Smoke run of the stress test: the first cases of every family may fail
only where tests/stress_baseline.json says they do.

Marked `stress` and excluded by default (~1.5 min on 4 processes):

    python3 -m pytest tests/ -m stress

A subprocess, because tests/stress.py pins BLAS to one thread before numpy
is imported (invariant K and the baseline need bit-identical results),
which cannot be done inside pytest once numpy is loaded.
"""

import subprocess
import sys
from pathlib import Path

import pytest

STRESS = Path(__file__).with_name('stress.py')

# Cases per family: index 0 also runs DE (invariant G, every 5th case)
SMOKE_N = 3


@pytest.mark.stress
def test_stress_smoke_matches_baseline():
    run = subprocess.run([sys.executable, str(STRESS), '--n', str(SMOKE_N), '--check'],
                         capture_output=True, text=True)
    assert run.returncode == 0, run.stdout[-3000:] + run.stderr[-3000:]

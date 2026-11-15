"""
Test that a storm with rainfall losses gives the same answer, and the same
volumes of excess and lost rain, on one process and on three.
"""

import os
import platform
import subprocess
import sys

import pytest

try:
    import mpi4py  # noqa: F401
except ImportError:
    pass

from anuga.file.sww import sww_files_are_equal

run_filename = os.path.join(os.path.dirname(__file__),
                            'run_parallel_rainfall_losses.py')


def run(cmd, cwd):
    result = subprocess.run(cmd.split(), capture_output=True, cwd=cwd, text=True)
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr)
        raise Exception(result.stderr)
    line = [ln for ln in result.stdout.splitlines() if ln.startswith('TOTALS')][-1]
    return [float(v) for v in line.split()[1:]]


def mpiexec_options():
    if platform.system() == 'Windows':
        return ''
    result = subprocess.run('mpiexec -np 3 --oversubscribe echo'.split(),
                            capture_output=True)
    return '--oversubscribe' if result.returncode == 0 else ''


@pytest.mark.skipif('mpi4py' not in sys.modules,
                    reason="requires the mpi4py module")
def test_that_sequential_and_parallel_outputs_are_identical(tmp_path):
    cwd = str(tmp_path)
    sequential = run('%s %s' % (sys.executable, run_filename), cwd)
    parallel = run('mpiexec -np 3 %s %s %s' % (mpiexec_options(), sys.executable,
                                               run_filename), cwd)

    assert sww_files_are_equal(os.path.join(cwd, 's_rainfall_losses.sww'),
                               os.path.join(cwd, 'p_rainfall_losses.sww'))

    excess, loss = sequential
    # 60 mm on 200 m x 100 m, all of it either excess or loss
    assert excess + loss == pytest.approx(0.060*200.0*100.0, rel=1e-9)
    assert 0.0 < loss < excess + loss
    assert parallel == pytest.approx(sequential, rel=1e-9)

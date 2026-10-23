"""Tracers survive distribute() and give the same answer in parallel (#278).

Runs run_parallel_tracers.py sequentially and under mpiexec, then compares.
The tracers are registered before distribute(), which is the case that used to
raise NotImplementedError because nothing carried them to the sub-domains.
"""

import os
import subprocess
import sys
import unittest

import numpy as num
import pytest

try:
    import mpi4py
except ImportError:
    pass


verbose = False

path = os.path.dirname(__file__)
run_filename = os.path.join(path, 'run_parallel_tracers.py')

paths_run_filename = os.path.join(path, 'run_parallel_tracers_paths.py')

sequential_file = 'tracers_sequential.txt'
parallel_file = 'tracers_parallel.txt'


def _mpi_prefix(np=3):
    """mpiexec plus --oversubscribe where the local MPI understands it."""
    import platform
    cmd = ['mpiexec', '-np', str(np)]
    if platform.system() == 'Windows':
        return cmd
    probe = subprocess.run(cmd + ['--oversubscribe', 'echo'],
                           capture_output=True)
    if probe.returncode == 0:
        cmd.append('--oversubscribe')
    return cmd


def _run(cmd):
    result = subprocess.run(cmd, capture_output=True)
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr)
        raise Exception(result.stderr)


def _read(filename):
    with open(filename) as fid:
        n = int(fid.readline())
        header = [float(x) for x in fid.readline().split()]
        rows = num.loadtxt(fid)
    return n, header, rows


@pytest.mark.skipif('mpi4py' not in sys.modules,
                    reason="requires the mpi4py module")
class Test_parallel_tracers(unittest.TestCase):
    def setUp(self):
        _run([sys.executable, run_filename])
        _run(_mpi_prefix(3) + [sys.executable, run_filename])

    def tearDown(self):
        for f in (sequential_file, parallel_file):
            try:
                os.remove(f)
            except OSError:
                pass

    def test_the_tracer_fields_are_the_same_in_parallel(self):
        n_seq, _, seq = _read(sequential_file)
        n_par, _, par = _read(parallel_file)

        assert n_seq == n_par, \
            'parallel run owns %d cells, sequential %d' % (n_par, n_seq)
        assert num.allclose(seq[:, :2], par[:, :2]), \
            'the two runs are not describing the same cells'

        # A structured tracer: if the partition mixed cells up, this is where
        # it shows. The run script already checked the field is not flat.
        assert num.allclose(seq[:, 2], par[:, 2], atol=1e-12), \
            'salinity differs: max %g' % num.abs(seq[:, 2] - par[:, 2]).max()
        assert num.allclose(seq[:, 3], par[:, 3], atol=1e-12), \
            'dye differs: max %g' % num.abs(seq[:, 3] - par[:, 3]).max()

    def test_the_field_actually_varies(self):
        """Otherwise the comparison above would pass on any partition at all."""
        _, _, seq = _read(sequential_file)
        assert num.ptp(seq[:, 2]) > 1e-3

    def test_the_mass_is_the_same_in_parallel(self):
        _, hs, _ = _read(sequential_file)
        _, hp, _ = _read(parallel_file)
        mass_seq, mass_par = hs[0], hp[0]
        assert abs(mass_seq - mass_par) < 1e-12 * max(abs(mass_seq), 1.0), \
            'total tracer mass differs: %g vs %g' % (mass_seq, mass_par)

    def test_the_budget_balances_across_ranks(self):
        """get_tracer_mass reduces; the flux integral must reduce to match."""
        _, hp, _ = _read(parallel_file)
        change, flux = hp[1], hp[2]
        assert abs(change) > 1e-6, 'no tracer left, so this proves nothing'
        assert abs(change - flux) < 1e-10 * abs(change), \
            'parallel budget does not balance: change %g vs flux %g' % (change, flux)


@pytest.mark.skipif('mpi4py' not in sys.modules,
                    reason="requires the mpi4py module")
class Test_parallel_tracers_other_paths(unittest.TestCase):
    """distribute_collaborative and the dump/load partition files.

    Both share partition_mesh with distribute(), so the values ride along; each
    has its own point where the reserved entries must be taken back out before
    they are mistaken for quantities. The run file asserts internally.
    """

    def test_the_other_distribution_paths_carry_tracers(self):
        _run(_mpi_prefix(3) + [sys.executable, paths_run_filename])


if __name__ == "__main__":
    runner = unittest.TextTestRunner()
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_parallel_tracers)
    runner.run(suite)

"""Sediment mass is conserved across ranks, and the parallel run agrees with
the sequential one (#423, #424).

The sediment kernels update ghost cells too, but from their state before the
step's ghost exchange, and elevation is not in that exchange, so ghost beds
drifted from their owners' and sediment mass leaked across partition
boundaries (about 0.5% of the eroded volume on 2 ranks). The sediment operator
now copies each owner's bed into its ghosts after every step. Bedload forms the
flux on a shared edge from the ghost's transport vector, so in mode 1 the ghosts
are refreshed before the step too; mode 2's device halo now carries the bed.

See run_parallel_sediment_mass.py for the set-up.
"""
import os
import subprocess
import sys
import unittest

import pytest

try:
    import mpi4py  # noqa: F401
except ImportError:
    pass

path = os.path.dirname(__file__)
run_filename = os.path.join(path, 'run_parallel_sediment_mass.py')


def _mpi_prefix(np=3):
    import platform
    cmd = ['mpiexec', '-np', str(np)]
    if platform.system() == 'Windows':
        return cmd
    probe = subprocess.run(cmd + ['--oversubscribe', 'echo'], capture_output=True)
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
    with open(filename) as f:
        n = int(f.readline())
        return n, [float(x) for x in f.readline().split()]


@pytest.mark.skipif('mpi4py' not in sys.modules, reason='requires mpi4py')
class Test_parallel_sediment_mass(unittest.TestCase):

    CASES = {'suspended': [], 'composition': ['--composition'],
             'bedload': ['--bedload'], 'bedload_unified': ['--bedload', '--unified'],
             'bedload_thin': ['--bedload', '--thin']}

    @classmethod
    def setUpClass(cls):
        cls.files = []
        cls.results = {}
        for name, flags in cls.CASES.items():
            seq, par = 'sediment_mass_%s_seq.txt' % name, 'sediment_mass_%s_par.txt' % name
            cls.files += [seq, par]
            _run([sys.executable, run_filename, seq] + flags)
            _run(_mpi_prefix(3) + [sys.executable, run_filename, par] + flags)
            cls.results[name] = (_read(seq), _read(par))

    @classmethod
    def tearDownClass(cls):
        for f in cls.files:
            if os.path.exists(f):
                os.remove(f)

    def _check(self, name, balance_tol=1e-10):
        (n_seq, seq), (n_par, par) = self.results[name]
        self.assertEqual((n_seq, n_par), (1, 3))
        eroded = abs(seq[3])
        self.assertGreater(eroded, 0.1)                 # the bed really moved
        # closed box: the solid balance stays at zero on every partitioning
        self.assertLess(abs(seq[0]), 1e-10 * eroded)
        self.assertLess(abs(par[0]), balance_tol * eroded)
        # and the parallel run is the sequential one
        for a, b in zip(par[1:], seq[1:]):
            self.assertAlmostEqual(a, b, delta=1e-9 * max(abs(b), 1.0))

    def test_suspended_sediment_and_bed(self):
        self._check('suspended')

    def test_with_bed_composition(self):
        self._check('composition')

    def test_bedload(self):
        self._check('bedload')

    def test_bedload_unified(self):
        self._check('bedload_unified')

    def test_bedload_on_an_exhausting_bed(self):
        """Where [L-5] binds, each rank sets a ghost's exhaustion flag from that
        ghost's own neighbourhood, which an outer ghost lacks part of, so the
        ranks can disagree about one edge: ~1e-10 of the eroded volume."""
        self._check('bedload_thin', balance_tol=1e-8)


if __name__ == '__main__':
    unittest.main()

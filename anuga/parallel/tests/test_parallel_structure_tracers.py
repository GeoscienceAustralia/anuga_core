"""Structures carry tracers in parallel: the tracer given up upstream arrives
downstream across sub-domains, exactly, and parallel agrees with sequential.

See run_parallel_structure_tracers.py for the set-up.
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
run_filename = os.path.join(path, 'run_parallel_structure_tracers.py')


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
class Test_parallel_structure_tracers(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.files = ['structure_tracers_seq.txt', 'structure_tracers_par.txt']
        _run([sys.executable, run_filename, cls.files[0]])
        _run(_mpi_prefix(3) + [sys.executable, run_filename, cls.files[1]])
        cls.results = {1: _read(cls.files[0]), 3: _read(cls.files[1])}

    @classmethod
    def tearDownClass(cls):
        for f in cls.files:
            if os.path.exists(f):
                os.remove(f)

    def test_the_runs_used_the_ranks_asked_for(self):
        self.assertEqual(self.results[1][0], 1)
        self.assertEqual(self.results[3][0], 3)

    def test_tracer_is_conserved_and_goes_through(self):
        for n, (m0, m1, up0, up1, dn0, dn1, flow) in self.results.values():
            self.assertGreater(flow, 1.0)
            self.assertAlmostEqual(m1, m0, delta=1e-12 * m0)
            self.assertAlmostEqual(dn1 - dn0, up0 - up1, delta=1e-9 * m0)
            self.assertGreater(dn1, 0.0)
            # arrives at the upstream concentration (to the small reverse flow)
            self.assertAlmostEqual((dn1 - dn0) / flow, 0.02, delta=5e-3 * 0.02)

    def test_parallel_matches_sequential(self):
        s, p = self.results[1][1], self.results[3][1]
        for a, b in zip(s, p):
            self.assertAlmostEqual(a, b, delta=1e-8 * max(abs(a), 1.0))


if __name__ == '__main__':
    unittest.main()

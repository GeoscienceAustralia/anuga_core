"""Inlet operators carry tracers in parallel: the budget closes on every rank
and the parallel run agrees with the sequential one.

See run_parallel_inlet_tracers.py for the set-up.
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
run_filename = os.path.join(path, 'run_parallel_inlet_tracers.py')
C0, C_IN, Q_IN, Q_OUT, T = 0.01, 0.05, 20.0, -30.0, 60.0


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
class Test_parallel_inlet_tracers(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.files = ['inlet_tracers_seq.txt', 'inlet_tracers_par.txt',
                     'inlet_tracers_seq_mixed.txt', 'inlet_tracers_par_mixed.txt']
        seq, par, seq_m, par_m = cls.files
        _run([sys.executable, run_filename, seq])
        _run(_mpi_prefix(3) + [sys.executable, run_filename, par])
        _run([sys.executable, run_filename, seq_m, '--mixed'])
        _run(_mpi_prefix(3) + [sys.executable, run_filename, par_m, '--mixed'])
        cls.results = {1: _read(seq), 3: _read(par)}
        cls.mixed = {1: _read(seq_m), 3: _read(par_m)}

    @classmethod
    def tearDownClass(cls):
        for f in cls.files:
            if os.path.exists(f):
                os.remove(f)

    def test_the_runs_used_the_ranks_asked_for(self):
        self.assertEqual(self.results[1][0], 1)
        self.assertEqual(self.results[3][0], 3)

    def test_the_budget_closes(self):
        """Tracer change = c_in x water in - c0 x water out, in both runs."""
        for n, (w0, w1, m0, m1, t_in, t_out) in self.results.values():
            water_in, water_out = Q_IN * T, -Q_OUT * T
            self.assertAlmostEqual(w1 - w0, water_in - water_out, delta=1e-6 * water_in)
            self.assertAlmostEqual(t_in, C_IN * water_in, delta=1e-9 * C_IN * water_in)
            self.assertAlmostEqual(t_out, C0 * water_out, delta=1e-9 * C0 * water_out)
            self.assertAlmostEqual(m1 - m0, t_in - t_out, delta=1e-9 * m0)

    def test_parallel_matches_sequential(self):
        for results in (self.results, self.mixed):
            s, p = results[1][1], results[3][1]
            for a, b in zip(s, p):
                self.assertAlmostEqual(a, b, delta=1e-9 * max(abs(a), 1.0))

    def test_a_pool_across_a_front_conserves_tracer(self):
        """The outlet straddles a concentration front and the sub-domains, so
        water moves between ranks inside the pool at unequal concentrations.
        The domain's tracer changes by exactly what the inlets moved."""
        for n, (w0, w1, m0, m1, t_in, t_out) in self.mixed.values():
            self.assertAlmostEqual(m1 - m0, t_in - t_out, delta=1e-9 * m0)
            self.assertGreater(t_out, 0.01 * 1800.0)   # more than c0 per unit water


if __name__ == '__main__':
    unittest.main()

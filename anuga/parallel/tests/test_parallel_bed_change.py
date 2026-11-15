"""A bed changed mid-run on owned cells only gives the sequential result in
parallel, in both compute modes.

Elevation is not in the per-step ghost exchange, so before the domain copied
the owners' bed into the ghosts after a host-side change, the ghosts kept the
old bed and the run drifted from serial (3e-4 to 1e-3 m here, in both modes).
With Set_elevation_operator in mode 2 the ghosts also kept the old stage into
the next step, and water was created. The water volume is checked too: mode 2
used to count ghost cells in it (+14% on 3 ranks).

See run_parallel_bed_change.py for the set-up.
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
run_filename = os.path.join(path, 'run_parallel_bed_change.py')


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
        volume = float(f.readline())
        stage = [float(v) for v in f.readline().split()]
    return n, volume, stage


@pytest.mark.skipif('mpi4py' not in sys.modules, reason='requires mpi4py')
class Test_parallel_bed_change(unittest.TestCase):

    CHANGES = ['none', 'set_quantity', 'direct', 'operator']

    @classmethod
    def setUpClass(cls):
        cls.files = []
        cls.results = {}
        for change in cls.CHANGES:
            seq = 'bed_change_%s_seq.txt' % change
            _run([sys.executable, run_filename, seq, change])
            cls.files.append(seq)
            cls.results[change, 'seq'] = _read(seq)
            for flags, key in (([], 'legacy'), (['--unified'], 'unified')):
                par = 'bed_change_%s_%s_par.txt' % (change, key)
                _run(_mpi_prefix(3) + [sys.executable, run_filename, par, change] + flags)
                cls.files.append(par)
                cls.results[change, key] = _read(par)

    @classmethod
    def tearDownClass(cls):
        for f in cls.files:
            if os.path.exists(f):
                os.remove(f)

    def _check(self, change, mode):
        n_seq, v_seq, seq = self.results[change, 'seq']
        n_par, v_par, par = self.results[change, mode]
        self.assertEqual((n_seq, n_par), (1, 3))
        self.assertEqual(len(seq), len(par))
        if change != 'none':
            ref = self.results['none', 'seq'][2]
            self.assertGreater(max(abs(a - b) for a, b in zip(seq, ref)), 1e-4,
                               'the bed change should alter the flow')
        self.assertLess(max(abs(a - b) for a, b in zip(par, seq)), 1e-12)
        self.assertAlmostEqual(v_par, v_seq, delta=1e-9 * v_seq)

    def test_no_change_legacy(self):
        self._check('none', 'legacy')

    def test_no_change_unified(self):
        self._check('none', 'unified')

    def test_set_quantity_legacy(self):
        self._check('set_quantity', 'legacy')

    def test_set_quantity_unified(self):
        self._check('set_quantity', 'unified')

    def test_direct_write_legacy(self):
        self._check('direct', 'legacy')

    def test_direct_write_unified(self):
        self._check('direct', 'unified')

    def test_set_elevation_operator_legacy(self):
        self._check('operator', 'legacy')

    def test_set_elevation_operator_unified(self):
        self._check('operator', 'unified')


if __name__ == '__main__':
    unittest.main()

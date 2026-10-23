"""
Tests for Parallel_domain.load_balance_statistics() and
print_load_balance_statistics().

Serial path (numprocs==1) is tested inline; the MPI path is smoke-tested by
spawning run_parallel_load_balance.py via mpiexec.
"""

import os
import sys
import subprocess
import unittest

import numpy as num
import pytest

try:
    import mpi4py
    MPI4PY_AVAILABLE = True
except ImportError:
    MPI4PY_AVAILABLE = False

path = os.path.dirname(os.path.abspath(__file__))
run_script = os.path.join(path, 'run_parallel_load_balance.py')


def _make_domain():
    """Small 10x10 wet domain, run one yieldstep so timing is populated."""
    import anuga
    from anuga import rectangular_cross_domain, Reflective_boundary

    domain = rectangular_cross_domain(10, 10)
    domain.set_flow_algorithm('DE0')
    domain.store = False
    domain.set_quantity('elevation', 0.0)
    domain.set_quantity('stage', 1.0)
    domain.set_boundary({'left':   Reflective_boundary(domain),
                         'right':  Reflective_boundary(domain),
                         'top':    Reflective_boundary(domain),
                         'bottom': Reflective_boundary(domain)})
    for _ in domain.evolve(yieldstep=0.1, finaltime=0.1):
        pass
    return domain


def _make_parallel_domain():
    """Wrap _make_domain() as a Parallel_domain with numprocs=1."""
    from anuga.parallel.parallel_shallow_water import Parallel_domain
    import anuga
    from anuga import rectangular_cross_domain, Reflective_boundary

    # Build a minimal Parallel_domain in serial (trivial decomposition)
    base = rectangular_cross_domain(10, 10)
    base.set_flow_algorithm('DE0')
    base.store = False
    base.set_quantity('elevation', 0.0)
    base.set_quantity('stage', 1.0)

    from anuga import distribute
    pdomain = distribute(base)

    pdomain.set_boundary({'left':   Reflective_boundary(pdomain),
                          'right':  Reflective_boundary(pdomain),
                          'top':    Reflective_boundary(pdomain),
                          'bottom': Reflective_boundary(pdomain)})
    for _ in pdomain.evolve(yieldstep=0.1, finaltime=0.1):
        pass
    return pdomain


class TestLoadBalanceStatisticsSerial(unittest.TestCase):
    """Tests that run without MPI (numprocs == 1)."""

    def setUp(self):
        try:
            import mpi4py  # noqa: F401
        except ImportError:
            self.skipTest('mpi4py not available')
        self.domain = _make_parallel_domain()

    def test_returns_dict_on_rank0(self):
        stats = self.domain.load_balance_statistics()
        # In serial, processor==0 always
        self.assertIsInstance(stats, dict)

    def test_expected_keys(self):
        stats = self.domain.load_balance_statistics()
        expected_keys = {
            'n_full', 'n_ghost', 'n_wet_full', 'wet_fraction',
            'ghost_fraction', 'wall_time', 'comm_time',
            'reduce_wait_time', 'compute_time',
            'imbalance_ratio', 'wet_compute_corr',
        }
        self.assertEqual(set(stats.keys()), expected_keys)

    def test_array_lengths_equal_numproc(self):
        stats = self.domain.load_balance_statistics()
        np_ = self.domain.numproc
        for key in ('n_full', 'n_ghost', 'n_wet_full', 'wet_fraction',
                    'ghost_fraction', 'wall_time', 'comm_time',
                    'reduce_wait_time', 'compute_time'):
            self.assertEqual(len(stats[key]), np_,
                             f'{key} has wrong length')

    def test_n_full_positive(self):
        stats = self.domain.load_balance_statistics()
        self.assertTrue(all(stats['n_full'] > 0))

    def test_wet_fraction_in_unit_interval(self):
        stats = self.domain.load_balance_statistics()
        for wf in stats['wet_fraction']:
            self.assertGreaterEqual(wf, 0.0)
            self.assertLessEqual(wf, 1.0)

    def test_all_wet_domain_wet_fraction_near_one(self):
        """Domain is fully wet — every full triangle should be wet."""
        stats = self.domain.load_balance_statistics()
        for wf in stats['wet_fraction']:
            self.assertGreater(wf, 0.9,
                msg=f'Expected wet fraction near 1 for fully-wet domain, got {wf}')

    def test_compute_time_nonnegative(self):
        stats = self.domain.load_balance_statistics()
        self.assertTrue(all(stats['compute_time'] >= 0.0))

    def test_imbalance_ratio_at_least_one(self):
        stats = self.domain.load_balance_statistics()
        self.assertGreaterEqual(stats['imbalance_ratio'], 1.0)

    def test_print_load_balance_statistics_runs(self):
        """print_load_balance_statistics should not raise."""
        self.domain.print_load_balance_statistics()

    def test_custom_minimum_height(self):
        """Custom minimum_height changes wet count."""
        stats_low  = self.domain.load_balance_statistics(minimum_height=0.0)
        stats_high = self.domain.load_balance_statistics(minimum_height=2.0)
        # With stage=1, elev=0 → depth=1. threshold 2.0 → all dry
        self.assertGreater(stats_low['n_wet_full'][0],
                           stats_high['n_wet_full'][0])

    def test_partially_dry_domain(self):
        """Half-dry domain should have wet fraction < 1."""
        import anuga
        from anuga import rectangular_cross_domain, Reflective_boundary
        from anuga import distribute

        base = rectangular_cross_domain(20, 10, len1=2.0, len2=1.0)
        base.set_flow_algorithm('DE0')
        base.store = False
        # Left half wet (stage=1, elev=0), right half dry (stage=elev=0.5)
        base.set_quantity('elevation',
                          lambda x, y: num.where(x < 1.0, 0.0, 0.5))
        base.set_quantity('stage',
                          lambda x, y: num.where(x < 1.0, 1.0, 0.5))

        pdomain = distribute(base)
        pdomain.set_boundary({'left':   Reflective_boundary(pdomain),
                              'right':  Reflective_boundary(pdomain),
                              'top':    Reflective_boundary(pdomain),
                              'bottom': Reflective_boundary(pdomain)})
        for _ in pdomain.evolve(yieldstep=0.05, finaltime=0.05):
            pass

        stats = pdomain.load_balance_statistics()
        total_wet = sum(stats['n_wet_full'])
        total_full = sum(stats['n_full'])
        overall_wet_frac = total_wet / total_full
        # About half the triangles should be wet (±25% tolerance)
        self.assertGreater(overall_wet_frac, 0.25)
        self.assertLess(overall_wet_frac, 0.75)


def _extra_mpiexec_options():
    import platform
    if platform.system() == 'Windows':
        return ''
    result = subprocess.run(
        ['mpiexec', '-np', '2', '--oversubscribe', 'echo', 'ok'],
        capture_output=True)
    return '--oversubscribe' if result.returncode == 0 else ''


@pytest.mark.slow
@pytest.mark.skipif(not MPI4PY_AVAILABLE, reason='requires mpi4py')
class TestLoadBalanceStatisticsMPI(unittest.TestCase):
    """Smoke test: run_parallel_load_balance.py exits 0 under mpiexec."""

    def _run(self, nprocs):
        extra = _extra_mpiexec_options()
        cmd = (['mpiexec', '-np', str(nprocs)]
               + (extra.split() if extra else [])
               + [sys.executable, run_script])
        result = subprocess.run(cmd, capture_output=True)
        if result.returncode != 0:
            raise AssertionError(
                f'mpiexec -np {nprocs} failed:\n'
                + result.stdout.decode()
                + result.stderr.decode())

    def test_mpi_2ranks(self):
        self._run(2)

    def test_mpi_4ranks(self):
        self._run(4)


if __name__ == '__main__':
    unittest.main()

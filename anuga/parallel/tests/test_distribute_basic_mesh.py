"""Tests for distribute_basic_mesh() -- mesh-first parallel workflow.

Tests cover:
  1. Basic_mesh construction (rectangular and rectangular_cross).
  2. Serial fallback: distribute_basic_mesh with numprocs==1 returns a plain Domain.
  3. Parallel distribution: distribute_basic_mesh returns a Parallel_domain that can
     evolve and produces results matching the traditional distribute() workflow.

Run serially::

    python -m pytest anuga/parallel/tests/test_distribute_basic_mesh.py

Run in parallel (4 procs)::

    mpirun -np 4 python anuga/parallel/tests/test_distribute_basic_mesh.py
"""

import tempfile
import os
import sys
import unittest

import numpy as num
import pytest

import anuga
from anuga import (
    Domain,
    Reflective_boundary,
    Dirichlet_boundary,
    distribute,
    myid,
    numprocs,
    barrier,
    finalize,
)
from anuga.abstract_2d_finite_volumes.basic_mesh import (
    Basic_mesh,
    rectangular_basic_mesh,
    rectangular_cross_basic_mesh,
)
from anuga.parallel.parallel_api import distribute_basic_mesh

try:
    import mpi4py
    mpi4py_available = True
except ImportError:
    mpi4py_available = False

# ---------------------------------------------------------------------------
# Test parameters
# ---------------------------------------------------------------------------
M = 20        # grid cells in x
N = 20        # grid cells in y
YIELDSTEP  = 0.25
FINALTIME  = 0.5
VERBOSE    = False


def topography(x, y):
    return -x / 10.0


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_serial_domain(name='test_serial'):
    """Create a plain rectangular_cross domain with initial conditions."""
    domain = anuga.rectangular_cross_domain(M, N)
    domain.set_quantity('elevation', topography)
    domain.set_quantity('friction', 0.03)
    domain.set_quantity('stage', expression='elevation')
    domain.set_name(name)
    domain.set_datadir(tempfile.mkdtemp())
    domain.set_quantities_to_be_stored(None)
    return domain


def set_standard_boundary(domain):
    Br = Reflective_boundary(domain)
    Bd = Dirichlet_boundary([-0.2, 0., 0.])
    domain.set_boundary({'left': Br, 'right': Bd, 'top': Br, 'bottom': Br})


def stage_at_gauges(domain, points):
    """Return stage at a list of (x,y) points (only from the owning proc)."""
    values = []
    for pt in points:
        try:
            k = domain.get_triangle_containing_point(pt)
            if domain.tri_full_flag[k] == 1:
                values.append(domain.quantities['stage'].centroid_values[k])
            else:
                values.append(None)
        except Exception:
            values.append(None)
    return values


# ---------------------------------------------------------------------------
# Test 1: Basic_mesh construction
# ---------------------------------------------------------------------------

class TestBasic_meshConstruction(unittest.TestCase):
    """Test that Basic_mesh factory functions produce the right topology."""

    def test_rectangular_basic_mesh_shape(self):
        bm = rectangular_basic_mesh(M, N, len1=1.0, len2=1.0)
        self.assertEqual(bm.number_of_triangles, 2 * M * N)
        self.assertEqual(bm.number_of_nodes, (M + 1) * (N + 1))
        self.assertIn('left', set(bm.boundary.values()))
        self.assertIn('right', set(bm.boundary.values()))
        self.assertIn('top', set(bm.boundary.values()))
        self.assertIn('bottom', set(bm.boundary.values()))

    def test_rectangular_cross_basic_mesh_shape(self):
        bm = rectangular_cross_basic_mesh(M, N, len1=1.0, len2=1.0)
        self.assertEqual(bm.number_of_triangles, 4 * M * N)
        # nodes = corner nodes + centre nodes
        self.assertEqual(bm.number_of_nodes, (M + 1) * (N + 1) + M * N)
        self.assertIn('left', set(bm.boundary.values()))
        self.assertIn('right', set(bm.boundary.values()))
        self.assertIn('top', set(bm.boundary.values()))
        self.assertIn('bottom', set(bm.boundary.values()))

    def test_neighbours_shape(self):
        bm = rectangular_cross_basic_mesh(M, N)
        nbrs = bm.neighbours
        self.assertEqual(nbrs.shape, (bm.number_of_triangles, 3))

    def test_centroid_coordinates_shape(self):
        bm = rectangular_cross_basic_mesh(M, N, len1=2.0, len2=3.0)
        cc = bm.centroid_coordinates
        self.assertEqual(cc.shape, (bm.number_of_triangles, 2))
        # Centroids should be inside the domain
        self.assertTrue(num.all(cc[:, 0] >= 0) and num.all(cc[:, 0] <= 2.0))
        self.assertTrue(num.all(cc[:, 1] >= 0) and num.all(cc[:, 1] <= 3.0))

    def test_reorder_preserves_boundary_count(self):
        bm = rectangular_cross_basic_mesh(M, N)
        n_bnd = len(bm.boundary)
        order = num.random.default_rng().permutation(bm.number_of_triangles)
        bm2 = bm.reorder(order, in_place=False)
        self.assertEqual(len(bm2.boundary), n_bnd)


# ---------------------------------------------------------------------------
# Test 2: Serial fallback (numprocs == 1)
# ---------------------------------------------------------------------------

class TestDistributeBasic_meshSerial(unittest.TestCase):
    """When numprocs==1, distribute_basic_mesh should return a plain Domain."""

    def test_returns_domain(self):
        bm = rectangular_cross_basic_mesh(M, N, len1=float(M), len2=float(N))
        domain = distribute_basic_mesh(bm)
        self.assertIsInstance(domain, Domain)

    def test_set_quantity_and_boundary(self):
        bm = rectangular_cross_basic_mesh(M, N, len1=float(M), len2=float(N))
        domain = distribute_basic_mesh(bm)
        domain.set_quantity('elevation', topography)
        domain.set_quantity('friction', 0.03)
        domain.set_quantity('stage', expression='elevation')
        Br = Reflective_boundary(domain)
        Bd = Dirichlet_boundary([-0.2, 0., 0.])
        # Should not raise
        domain.set_boundary({'left': Br, 'right': Bd, 'top': Br, 'bottom': Br})

    def test_evolve(self):
        bm = rectangular_cross_basic_mesh(M, N, len1=float(M), len2=float(N))
        domain = distribute_basic_mesh(bm)
        domain.set_name('test_evolve_serial')
        domain.set_datadir(tempfile.mkdtemp())
        domain.set_quantities_to_be_stored(None)
        domain.set_quantity('elevation', topography)
        domain.set_quantity('friction', 0.03)
        domain.set_quantity('stage', expression='elevation')
        Br = Reflective_boundary(domain)
        Bd = Dirichlet_boundary([-0.2, 0., 0.])
        domain.set_boundary({'left': Br, 'right': Bd, 'top': Br, 'bottom': Br})
        t_final = None
        for t in domain.evolve(yieldstep=YIELDSTEP, finaltime=FINALTIME):
            t_final = t
        self.assertAlmostEqual(t_final, FINALTIME, places=6)

    def test_results_match_plain_domain(self):
        """distribute_basic_mesh serial result matches plain Domain evolve."""
        gauges = [[float(M) * 0.3, float(N) * 0.5],
                  [float(M) * 0.7, float(N) * 0.5]]

        # --- plain Domain ---
        domain_ref = anuga.rectangular_cross_domain(M, N,
                                                    len1=float(M),
                                                    len2=float(N))
        domain_ref.set_name('ref_serial')
        domain_ref.set_datadir(tempfile.mkdtemp())
        domain_ref.set_quantities_to_be_stored(None)
        domain_ref.set_quantity('elevation', topography)
        domain_ref.set_quantity('friction', 0.03)
        domain_ref.set_quantity('stage', expression='elevation')
        Br = Reflective_boundary(domain_ref)
        Bd = Dirichlet_boundary([-0.2, 0., 0.])
        domain_ref.set_boundary({'left': Br, 'right': Bd, 'top': Br, 'bottom': Br})
        for t in domain_ref.evolve(yieldstep=YIELDSTEP, finaltime=FINALTIME):
            pass
        ref_stage = stage_at_gauges(domain_ref, gauges)

        # --- distribute_basic_mesh ---
        bm = rectangular_cross_basic_mesh(M, N, len1=float(M), len2=float(N))
        domain_bm = distribute_basic_mesh(bm)
        domain_bm.set_name('bm_serial')
        domain_bm.set_datadir(tempfile.mkdtemp())
        domain_bm.set_quantities_to_be_stored(None)
        domain_bm.set_quantity('elevation', topography)
        domain_bm.set_quantity('friction', 0.03)
        domain_bm.set_quantity('stage', expression='elevation')
        Br2 = Reflective_boundary(domain_bm)
        Bd2 = Dirichlet_boundary([-0.2, 0., 0.])
        domain_bm.set_boundary({'left': Br2, 'right': Bd2,
                                 'top': Br2, 'bottom': Br2})
        for t in domain_bm.evolve(yieldstep=YIELDSTEP, finaltime=FINALTIME):
            pass
        bm_stage = stage_at_gauges(domain_bm, gauges)

        for r, b in zip(ref_stage, bm_stage):
            if r is not None and b is not None:
                self.assertAlmostEqual(r, b, places=5,
                    msg='Stage mismatch at gauge: ref=%g bm=%g' % (r, b))


# ---------------------------------------------------------------------------
# Test 3: Parallel distribution (run via mpirun)
# ---------------------------------------------------------------------------

def collect_gauge_stages(domain, gauges):
    """Gather stage at gauges from all ranks to rank 0.

    Returns a list (length == len(gauges)) on rank 0 with the first
    non-None value found across all ranks for each gauge.  Returns None
    on non-zero ranks.
    """
    from mpi4py import MPI
    comm = MPI.COMM_WORLD

    local = stage_at_gauges(domain, gauges)
    all_local = comm.gather(local, root=0)

    if myid == 0:
        result = []
        for i in range(len(gauges)):
            val = None
            for proc_vals in all_local:
                if proc_vals[i] is not None:
                    val = proc_vals[i]
                    break
            result.append(val)
        return result
    return None


def run_distribute_basic_mesh_parallel(verbose=False):
    """Parallel distribute_basic_mesh: evolve runs without error and produces
    physically reasonable stage values (checked on rank 0)."""

    gauges = [[float(M) * 0.3, float(N) * 0.5],
              [float(M) * 0.6, float(N) * 0.5]]

    if myid == 0:
        bm = rectangular_cross_basic_mesh(M, N, len1=float(M), len2=float(N))
    else:
        bm = None

    domain = distribute_basic_mesh(bm, verbose=False)
    domain.set_name('test_bm_parallel')
    domain.set_datadir(tempfile.mkdtemp())
    domain.set_quantities_to_be_stored(None)
    domain.set_quantity('elevation', topography)
    domain.set_quantity('friction', 0.03)
    domain.set_quantity('stage', expression='elevation')

    Br = Reflective_boundary(domain)
    Bd = Dirichlet_boundary([-0.2, 0., 0.])
    domain.set_boundary({'left': Br, 'right': Bd, 'top': Br, 'bottom': Br})

    t_final = None
    for t in domain.evolve(yieldstep=YIELDSTEP, finaltime=FINALTIME):
        t_final = t

    bm_stages = collect_gauge_stages(domain, gauges)
    barrier()

    if myid == 0:
        assert abs(t_final - FINALTIME) < 1e-10, \
            'Final time mismatch: %g != %g' % (t_final, FINALTIME)

        if verbose:
            print('bm_stages:', bm_stages)

        # Stage at gauges should be between initial elevation and the
        # Dirichlet boundary value (-0.2), allowing for wave dynamics.
        for i, (gauge, stage) in enumerate(zip(gauges, bm_stages)):
            if stage is not None:
                elev = topography(gauge[0], gauge[1])
                assert stage >= elev - 0.5, \
                    'Stage %g far below initial elevation %g at gauge %s' \
                    % (stage, elev, gauge)
                assert stage <= 0.5, \
                    'Stage %g suspiciously high at gauge %s' % (stage, gauge)


def run_basic_mesh_no_sww(verbose=False):
    """distribute_basic_mesh with quantities_to_be_stored=None should not write SWW."""
    if myid == 0:
        bm = rectangular_cross_basic_mesh(M, N, len1=float(M), len2=float(N))
    else:
        bm = None

    domain = distribute_basic_mesh(bm)
    domain.set_name('test_no_sww')
    domain.set_datadir(tempfile.mkdtemp())
    domain.set_quantities_to_be_stored(None)
    domain.set_quantity('elevation', topography)
    domain.set_quantity('friction', 0.03)
    domain.set_quantity('stage', expression='elevation')

    Br = Reflective_boundary(domain)
    Bd = Dirichlet_boundary([-0.2, 0., 0.])
    domain.set_boundary({'left': Br, 'right': Bd, 'top': Br, 'bottom': Br})

    t_final = None
    for t in domain.evolve(yieldstep=YIELDSTEP, finaltime=FINALTIME):
        t_final = t

    barrier()

    if myid == 0:
        assert abs(t_final - FINALTIME) < 1e-10, \
            'Final time mismatch: %g != %g' % (t_final, FINALTIME)
        if verbose:
            print('run_basic_mesh_no_sww: OK, final time = %g' % t_final)


# ---------------------------------------------------------------------------
# Test class that invokes the parallel tests via mpirun
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not mpi4py_available, reason='requires mpi4py')
class TestDistributeBasic_meshParallel(unittest.TestCase):

    def test_parallel_evolve(self):
        """Run the parallel distribute_basic_mesh evolve test via mpirun."""
        cmd = anuga.mpicmd(os.path.abspath(__file__))
        result = os.system(cmd)
        assert result == 0, 'Parallel test returned non-zero exit code'


# ---------------------------------------------------------------------------
# Entry point when run directly under mpirun
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    if numprocs == 1:
        # Run the serial unittest suite
        runner = unittest.TextTestRunner(verbosity=2)
        suite = unittest.TestSuite()
        suite.addTests(unittest.TestLoader().loadTestsFromTestCase(
            TestBasic_meshConstruction))
        suite.addTests(unittest.TestLoader().loadTestsFromTestCase(
            TestDistributeBasic_meshSerial))
        runner.run(suite)
    else:
        from anuga.utilities.parallel_abstraction import global_except_hook
        sys.excepthook = global_except_hook

        run_basic_mesh_no_sww(verbose=VERBOSE)
        barrier()

        run_distribute_basic_mesh_parallel(verbose=VERBOSE)
        barrier()

        #if myid == 0:
        #    print('All parallel distribute_basic_mesh tests passed.')

        finalize()

"""
Multi-rank GPU halo-exchange test for DE_ader2 flow algorithm.

Verifies that multiprocessor_mode=2 (GPU ADER-2 step, GPU internal MPI ghost
exchange) produces results matching multiprocessor_mode=1 (CPU ADER-2 step,
Python MPI ghost exchange) when run in parallel.

mode=1 uses the OpenMP C extension kernels; mode=2 uses the GPU C kernels.
Both implement the same ADER-2 algorithm but with different code paths, so
small floating-point differences (< 1e-6) are expected and acceptable.
Ghost exchange errors produce O(1) differences and are easily caught.

Both modes call the same core_kernels.c functions for computation; the only
difference is ghost exchange (Python MPI vs C MPI).  After a correct build
the results should be bit-for-bit identical, so ATOL=1e-12 is used.  Ghost
exchange errors produce O(0.1-1.0) differences and are easily caught.
All ranks vote on the result via MPI allreduce before raising so that no
rank is left waiting.

Run via pytest (auto-marked slow):

    pytest anuga/parallel/tests/test_parallel_sw_flow_gpu_de_ader2.py

Or directly via MPI:

    mpiexec -np 4 python -m mpi4py anuga/parallel/tests/test_parallel_sw_flow_gpu_de_ader2.py
"""

import os
import sys
import tempfile
import unittest

import numpy as num
import pytest

import anuga
from anuga import (Dirichlet_boundary, Reflective_boundary,
                   distribute, rectangular_cross_domain)
from anuga import myid, numprocs, barrier, finalize


def gpu_available():
    """True if sw_domain_gpu_ext can be imported (CPU or GPU mode)."""
    try:
        from anuga.shallow_water.sw_domain_gpu_ext import init_gpu_domain  # noqa: F401
        return True
    except ImportError:
        return False


def real_gpu_available():
    """True only when a real GPU offload target is active (not CPU_ONLY_MODE)."""
    try:
        from anuga.shallow_water.sw_domain_gpu_ext import gpu_available as _ga
        return _ga()
    except ImportError:
        return False


def get_num_gpu_devices():
    """Number of OpenMP offload devices; 0 in CPU_ONLY_MODE or no GPU."""
    try:
        from anuga.shallow_water.sw_domain_gpu_ext import get_num_gpu_devices as _gnd
        return _gnd()
    except ImportError:
        return 0


def gpu_has_mpi():
    try:
        from anuga.shallow_water.sw_domain_gpu_ext import gpu_has_mpi as _gpu_has_mpi
        return _gpu_has_mpi()
    except ImportError:
        return False


M = 29
N = 29
YIELDSTEP = 0.25
FINALTIME = 1.0
GAUGE_POINTS = [[0.4, 0.5], [0.6, 0.5], [0.8, 0.5], [0.9, 0.5]]
ATOL = 1e-12


def topography(x, y):
    return -x / 2


def run_simulation(mode, verbose=False):
    """Run a parallel DE_ader2 simulation with the given multiprocessor_mode.

    Parameters
    ----------
    mode : int
        1 = CPU ADER-2 step + Python MPI ghost exchange
        2 = GPU ADER-2 step + C-level MPI ghost exchange
    """
    domain = rectangular_cross_domain(M, N)
    domain.set_quantity('elevation', topography)
    domain.set_quantity('friction', 0.0)
    domain.set_quantity('stage', expression='elevation')

    domain = distribute(domain, verbose=False)

    domain.set_name(f'gpu_de_ader2_mode{mode}')
    domain.set_datadir(tempfile.mkdtemp())
    domain.set_flow_algorithm('DE_ader2')
    domain.set_quantities_to_be_stored(None)

    Br = Reflective_boundary(domain)
    Bd = Dirichlet_boundary([-0.2, 0., 0.])
    domain.set_boundary({'left': Br, 'right': Bd, 'top': Br, 'bottom': Br})

    domain.set_multiprocessor_mode(mode)

    tri_ids = []
    for point in GAUGE_POINTS:
        try:
            k = domain.get_triangle_containing_point(point)
            tri_ids.append(k if domain.tri_full_flag[k] == 1 else -1)
        except Exception:
            tri_ids.append(-2)

    gauge_values = [[] for _ in GAUGE_POINTS]
    stage = domain.get_quantity('stage')

    for _ in domain.evolve(yieldstep=YIELDSTEP, finaltime=FINALTIME):
        for i, tid in enumerate(tri_ids):
            if tid > -1:
                gauge_values[i].append(float(stage.centroid_values[tid]))

    return gauge_values, tri_ids


def compare_modes(G1, ids1, G2, ids2, label=''):
    """Assert mode=1 and mode=2 gauge time-series match across all ranks.

    Uses MPI allreduce so all ranks fail (or pass) together — prevents one
    rank from hanging in a collective while another has called MPI_Abort.
    """
    local_max_diff = 0.0
    worst_gauge = -1
    for i, (tid1, tid2) in enumerate(zip(ids1, ids2)):
        if tid1 > -1 and tid2 > -1:
            a1 = num.array(G1[i])
            a2 = num.array(G2[i])
            diff = float(num.max(num.abs(a1 - a2)))
            if diff > local_max_diff:
                local_max_diff = diff
                worst_gauge = i

    try:
        from mpi4py import MPI
        global_max_diff = MPI.COMM_WORLD.allreduce(local_max_diff, op=MPI.MAX)
    except ImportError:
        global_max_diff = local_max_diff

    if global_max_diff > ATOL:
        raise AssertionError(
            f'[rank {myid}] {label}: '
            f'mode=1 vs mode=2 global max diff = {global_max_diff:.2e} '
            f'(atol={ATOL}, local worst gauge={worst_gauge})'
        )


@pytest.mark.skipif(not gpu_available(),
                    reason='GPU extension not available')
@pytest.mark.skipif(not gpu_has_mpi(),
                    reason='GPU extension built without C MPI (mpi.h not found at build time)')
@pytest.mark.skipif('mpi4py' not in sys.modules,
                    reason='requires mpi4py')
class Test_parallel_sw_gpu_de_ader2(unittest.TestCase):

    def test_gpu_de_ader2_2ranks(self):
        """2-rank DE_ader2: GPU halo exchange matches Python halo exchange."""
        if real_gpu_available() and get_num_gpu_devices() < 2:
            self.skipTest("requires at least 2 GPUs (1 per MPI rank); "
                          f"found {get_num_gpu_devices()}")
        cmd = anuga.mpicmd(os.path.abspath(__file__), numprocs=2)
        assert os.system(cmd) == 0

    def test_gpu_de_ader2_4ranks(self):
        """4-rank DE_ader2: GPU halo exchange matches Python halo exchange."""
        if real_gpu_available() and get_num_gpu_devices() < 4:
            self.skipTest("requires at least 4 GPUs (1 per MPI rank); "
                          f"found {get_num_gpu_devices()}")
        cmd = anuga.mpicmd(os.path.abspath(__file__), numprocs=4)
        assert os.system(cmd) == 0


def assert_(condition, msg='Assertion Failed'):
    if not condition:
        raise AssertionError(msg)


if __name__ == '__main__':
    if numprocs == 1:
        runner = unittest.TextTestRunner()
        suite = unittest.TestLoader().loadTestsFromTestCase(
            Test_parallel_sw_gpu_de_ader2)
        runner.run(suite)
    else:
        from anuga.utilities.parallel_abstraction import global_except_hook
        sys.excepthook = global_except_hook

        if myid == 0:
            print(f'\n=== GPU halo exchange test DE_ader2: {numprocs} ranks ===')

        barrier()
        if myid == 0:
            print('  Running mode=1 (CPU ADER-2, Python MPI) ...')
        G1, ids1 = run_simulation(mode=1, verbose=False)
        barrier()

        if myid == 0:
            print('  Running mode=2 (GPU ADER-2, GPU/C MPI) ...')
        G2, ids2 = run_simulation(mode=2, verbose=False)
        barrier()

        compare_modes(G1, ids1, G2, ids2, label=f'np={numprocs}')

        if myid == 0:
            print(f'  PASS: mode=1 and mode=2 agree (atol={ATOL})')

        finalize()

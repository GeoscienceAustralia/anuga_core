"""
Simple water flow example using ANUGA  (DE2 flow algorithm)

Water driven up a linear slope and time varying boundary,
similar to a beach environment.

Tests that the parallel DE2 evolve produces the same gauge-station
time series as the equivalent sequential run.
"""

import tempfile
import unittest
import os
import sys
import numpy as num

import anuga

from anuga import Reflective_boundary, Dirichlet_boundary
from anuga import rectangular_cross_domain
from anuga import distribute, myid, numprocs, send, receive, barrier, finalize

try:
    import mpi4py
except ImportError:
    pass

import pytest

#--------------------------------------------------------------------------
# Setup parameters
#--------------------------------------------------------------------------
yieldstep = 0.25
finaltime = 1.0
nprocs = 4
N = 29
M = 29
verbose = False

#---------------------------------
# Setup Functions
#---------------------------------
def topography(x, y):
    return -x / 2

###########################################################################
# Setup Test
##########################################################################
def run_simulation(parallel=False, G=None, seq_interpolation_points=None, verbose=False):

    #--------------------------------------------------------------------------
    # Setup computational domain and quantities
    #--------------------------------------------------------------------------
    domain = rectangular_cross_domain(M, N)
    domain.set_quantity('elevation', topography)
    domain.set_quantity('friction', 0.0)
    domain.set_quantity('stage', expression='elevation')

    #--------------------------------------------------------------------------
    # Create the parallel domain
    #--------------------------------------------------------------------------
    if parallel:
        if myid == 0 and verbose:
            print('DISTRIBUTING PARALLEL DOMAIN')
        domain = distribute(domain, verbose=False)

    #--------------------------------------------------------------------------
    # Setup domain parameters — flow algorithm set after distribute
    #--------------------------------------------------------------------------
    domain.set_name('runup')
    domain.set_datadir(tempfile.mkdtemp())
    domain.set_flow_algorithm('DE2')
    domain.set_quantities_to_be_stored(None)

    #--------------------------------------------------------------------------
    # Setup boundary conditions  (must happen after distribute)
    #--------------------------------------------------------------------------
    Br = Reflective_boundary(domain)
    Bd = Dirichlet_boundary([-0.2, 0., 0.])
    domain.set_boundary({'left': Br, 'right': Bd, 'top': Br, 'bottom': Br})

    #--------------------------------------------------------------------------
    # Locate gauge triangles
    # Interpolation points may straddle partition boundaries in the parallel
    # run, so we snap them to the centroid found in the sequential run.
    #--------------------------------------------------------------------------
    interpolation_points = [[0.4, 0.5], [0.6, 0.5], [0.8, 0.5], [0.9, 0.5]]

    gauge_values = []
    tri_ids = []
    for i, point in enumerate(interpolation_points):
        gauge_values.append([])
        try:
            k = domain.get_triangle_containing_point(point)
            if domain.tri_full_flag[k] == 1:
                tri_ids.append(k)
            else:
                tri_ids.append(-1)
        except Exception:
            tri_ids.append(-2)

    if verbose:
        print('P%d has points = %s' % (myid, tri_ids))

    c_coord = domain.get_centroid_coordinates()
    interpolation_points = []
    for tid in tri_ids:
        if tid < 1:
            if verbose:
                print('WARNING: Interpolation point not within the domain!')
        interpolation_points.append(c_coord[tid, :])

    #--------------------------------------------------------------------------
    # Evolve
    #--------------------------------------------------------------------------
    if parallel:
        if myid == 0 and verbose:
            print('PARALLEL EVOLVE')
    else:
        if myid == 0 and verbose:
            print('SEQUENTIAL EVOLVE')

    time = []
    for t in domain.evolve(yieldstep=yieldstep, finaltime=finaltime):
        if myid == 0 and verbose:
            domain.write_time()

        time.append(domain.get_time())

        stage = domain.get_quantity('stage')
        for i in range(4):
            if tri_ids[i] > -1:
                gauge_values[i].append(stage.centroid_values[tri_ids[i]])

    #--------------------------------------------------------------------------
    # Compare against sequential reference on parallel runs
    #--------------------------------------------------------------------------
    if not parallel:
        G = []
        for i in range(4):
            G.append(gauge_values[i])

    success = True
    for i in range(4):
        if tri_ids[i] > -1:
            success = success and num.allclose(gauge_values[i], G[i])

    assert_(success)

    return G, interpolation_points


@pytest.mark.skipif('mpi4py' not in sys.modules,
                    reason="requires the mpi4py module")
class Test_parallel_sw_flow_de2(unittest.TestCase):
    def test_parallel_sw_flow_de2(self):
        cmd = anuga.mpicmd(os.path.abspath(__file__))
        result = os.system(cmd)
        assert_(result == 0)


def assert_(condition, msg="Assertion Failed"):
    if condition == False:
        raise AssertionError(msg)


if __name__ == "__main__":
    if numprocs == 1:
        runner = unittest.TextTestRunner()
        suite = unittest.TestLoader().loadTestsFromTestCase(
            Test_parallel_sw_flow_de2)
        runner.run(suite)
    else:
        barrier()
        if myid == 0 and verbose:
            print('SEQUENTIAL START')

        G, interpolation_points = run_simulation(parallel=False, verbose=verbose)
        G = num.array(G, float)

        barrier()

        if myid == 0 and verbose:
            print('PARALLEL START')

        from anuga.utilities.parallel_abstraction import global_except_hook
        sys.excepthook = global_except_hook

        run_simulation(parallel=True, G=G,
                       seq_interpolation_points=interpolation_points,
                       verbose=verbose)

        finalize()

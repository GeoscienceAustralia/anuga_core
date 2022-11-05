
import os.path
import sys

from anuga.utilities.system_tools import get_pathname_from_package
from anuga.geometry.polygon_function import Polygon_function

from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
from anuga.abstract_2d_finite_volumes.quantity import Quantity
from anuga.abstract_2d_finite_volumes.util import file_function

import anuga

import warnings
warnings.simplefilter("ignore")

#from anuga.structures.boyd_box_operator import Boyd_box_operator
#from anuga.structures.inlet_operator import Inlet_operator

#from anuga.culvert_flows.culvert_routines import boyd_generalised_culvert_model

from anuga.utilities import parallel_abstraction as pypar

from math import pi, pow, sqrt

import numpy as num

import anuga
from anuga.parallel.parallel_inlet_operator import Parallel_Inlet_operator
from anuga.parallel import distribute, myid, numprocs, finalize
from anuga.geometry.polygon import inside_polygon, is_inside_polygon, line_intersect

from anuga.parallel.parallel_operator_factory import Inlet_operator, Boyd_box_operator

import random
import unittest
import pickle
import pathlib

# Setup to skip test if mpi4py not available
import sys
try:
    import mpi4py
except ImportError:
    pass

import pytest

run_path = pathlib.Path(__file__).parent.resolve()
serial_output_path = os.path.join(run_path, 'serial_boyd.p')

"""

This test exercises the parallel culvert and checks values
"""
verbose = False
nprocs = 3


length = 40.
width = 15.

dx = dy = 0.5           # Resolution: Length of subdivisions on both axes

#----------------------------------------------------------------------
# Setup initial conditions
#----------------------------------------------------------------------

def topography(x, y):
    """Set up a weir

    A culvert will connect either side
    """
    # General Slope of Topography
    z=-x/1000

    N = len(x)
    for i in range(N):

        # Sloping Embankment Across Channel
        if 5.0 < x[i] < 10.1:
            # Cut Out Segment for Culvert face
            if  1.0+(x[i]-5.0)/5.0 <  y[i]  < 4.0 - (x[i]-5.0)/5.0:
                z[i]=z[i]
            else:
                z[i] +=  0.5*(x[i] -5.0)    # Sloping Segment  U/S Face
        if 10.0 < x[i] < 12.1:
            z[i] +=  2.5                    # Flat Crest of Embankment
        if 12.0 < x[i] < 14.5:
            # Cut Out Segment for Culvert face
            if  2.0-(x[i]-12.0)/2.5 <  y[i]  < 3.0 + (x[i]-12.0)/2.5:
                z[i]=z[i]
            else:
                z[i] +=  2.5-1.0*(x[i] -12.0) # Sloping D/S Face


    return z

def stage(x,y):

    z = topography(x,y)

    N = len(x)
    for i in range(N):
        if x[i] < 10:
            z[i] += 2
    return z
#filename=os.path.join(path, 'example_rating_curve.csv')

#mod_path = get_pathname_from_package('anuga.parallel')

def run_simulation(parallel, verbose=False):

##-----------------------------------------------------------------------
## Setup domain
##-----------------------------------------------------------------------
    if parallel:
        if myid == 0:
            points, vertices, boundary = rectangular_cross(int(length/dx),
                                                           int(width/dy),
                                                           len1=length,
                                                           len2=width)

            domain = anuga.Domain(points, vertices, boundary)
            domain.set_store(False)
            domain.set_flow_algorithm('DE0')
            domain.set_store_vertices_uniquely(True)
            domain.set_quantity('elevation', topography)
            domain.set_quantity('friction', 0.01)         # Constant friction
            # domain.set_quantity('stage',
            #                     expression='elevation')   # Dry initial condition
            domain.set_quantity('stage',stage)
        else:
            domain = None
    else:
        points, vertices, boundary = rectangular_cross(int(length/dx),
                                                       int(width/dy),
                                                       len1=length,
                                                       len2=width)

        domain = anuga.Domain(points, vertices, boundary)
        domain.set_store(False)
        domain.set_flow_algorithm('DE0')
        domain.set_store_vertices_uniquely(True)
        domain.set_quantity('elevation', topography)
        domain.set_quantity('friction', 0.01)         # Constant friction
        # domain.set_quantity('stage',
        #                     expression='elevation')   # Dry initial condition
        domain.set_quantity('stage',stage)

##-----------------------------------------------------------------------
## Distribute domain
##-----------------------------------------------------------------------

    if parallel:
        domain = distribute(domain)
        #domain.dump_triangulation("frac_op_domain.png")


##-----------------------------------------------------------------------
## Setup boundary conditions
##-----------------------------------------------------------------------

    Br = anuga.Reflective_boundary(domain)              # Solid reflective wall
    domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})



    gate = Boyd_box_operator(domain,
                            end_points=[[9.0, 2.5],[13.0, 2.5]],
                            losses=1.5,
                            width=1.5,
                            height = 0.0001,
                            apron=5.0,
                            use_momentum_jet=True,
                            use_velocity_head=False,
                            manning=0.013,
                            verbose=False)

    if gate is not None:
        gate.set_culvert_height(10.0)

    ##-----------------------------------------------------------------------
    ## Evolve system through time
    ##-----------------------------------------------------------------------

    if verbose:
        if parallel:
            if myid == 0:
                print("========== Parallel ==========")
        else:
            print("========== Serial ==========")

    if parallel:
        pypar.barrier()

    return_output = 0
    for t in domain.evolve(yieldstep = 0.01, finaltime = 0.02):
        if myid == 0 and verbose:
            domain.write_time()

        if gate is not None:
            output = gate.discharge_routine()
            if myid == gate.get_master_proc():
                return_output = output
                #print('myid ',myid,output)

        if gate is not None and verbose:
            if myid == gate.get_master_proc():
                print('master_proc_id', myid, 'discharge', return_output)

    if parallel:
        pypar.barrier()

    if parallel:
        if gate is not None:
            return gate.get_master_proc(), return_output
        else:
            return -1, return_output
    else:
        pickle.dump(return_output, open(serial_output_path,"wb"))



@pytest.mark.skipif('mpi4py' not in sys.modules,
                    reason="requires the mpi4py module")
class Test_parallel_boyd_box_operator_consistency(unittest.TestCase):
    def test_parallel_operator(self):
        #print "Expect this test to fail if not run from the parallel/test directory."

        os.system(anuga.mpicmd(os.path.abspath(__file__),1))
        cmd = anuga.mpicmd(os.path.abspath(__file__))
        exitstatus = os.system(cmd)
        os.remove(serial_output_path)
        assert_(exitstatus == 0)


# Because we are doing assertions outside of the TestCase class
# the PyUnit defined assert_ function can't be used.
def assert_(condition, msg="Assertion Failed"):
    if condition == False:
        #pypar.finalize()
        raise (AssertionError, msg)

if __name__=="__main__":
    if numprocs == 1:
        run_simulation(parallel=False, verbose=verbose)

    else:
        master_proc, parallel_output = run_simulation(parallel=True, verbose=verbose)
        success = True
        if myid == master_proc:
            serial_output = pickle.load(open(serial_output_path,"rb"))
            success = num.allclose(parallel_output, serial_output)

        finalize()
        import sys
        if success:
            sys.exit(0)
        else:
            sys.exit(1)

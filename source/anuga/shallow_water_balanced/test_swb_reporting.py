#!/usr/bin/env python

import unittest, os
import os.path
from math import pi, sqrt
import tempfile

from anuga.config import g, epsilon
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a
from anuga.utilities.numerical_tools import mean
from anuga.utilities.polygon import is_inside_polygon
from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.abstract_2d_finite_volumes.quantity import Quantity
from anuga.geospatial_data.geospatial_data import Geospatial_data
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross

from anuga.utilities.system_tools import get_pathname_from_package
from swb_domain import *

import numpy as num

# Get gateway to C implementation of flux function for direct testing
from shallow_water_ext import flux_function_central as flux_function




class Test_swb_clean(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass


    def test_extrema(self):
        """Test that extrema of quantities are computed correctly
        Extrema are updated at every *internal* timestep
        """

        from anuga.abstract_2d_finite_volumes.mesh_factory \
                import rectangular_cross

        initial_runup_height = -0.4
        final_runup_height = -0.3

        #--------------------------------------------------------------
        # Setup computational domain
        #--------------------------------------------------------------
        N = 5
        points, vertices, boundary = rectangular_cross(N, N)
        domain = Domain(points, vertices, boundary)
        domain.set_name('extrema_test')

        #--------------------------------------------------------------
        # Setup initial conditions
        #--------------------------------------------------------------
        def topography(x,y):
            return -x/2                             # linear bed slope

        domain.set_quantity('elevation', topography)    # function for elevation
        domain.set_quantity('friction', 0.)             # Zero friction
        # Constant negative initial stage
        domain.set_quantity('stage', initial_runup_height)
        domain.set_quantities_to_be_monitored(['stage', 'stage-elevation'],
                                              time_interval=[0.5, 2.7],
                                              polygon=[[0,0], [0,1],
                                                       [1,1], [1,0]])

        assert len(domain.quantities_to_be_monitored) == 2
        assert domain.quantities_to_be_monitored.has_key('stage')
        assert domain.quantities_to_be_monitored.has_key('stage-elevation')
        for key in domain.quantities_to_be_monitored['stage'].keys():
            assert domain.quantities_to_be_monitored['stage'][key] is None

        #--------------------------------------------------------------
        # Setup boundary conditions
        #--------------------------------------------------------------
        Br = Reflective_boundary(domain)              # Reflective wall
        # Constant inflow
        Bd = Dirichlet_boundary([final_runup_height, 0, 0])

        # All reflective to begin with (still water)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        #--------------------------------------------------------------
        # Let triangles adjust and check extrema
        #--------------------------------------------------------------
        for t in domain.evolve(yieldstep=0.1, finaltime=1.0):
            domain.quantity_statistics() # Run it silently

        #--------------------------------------------------------------
        # Test extrema
        #--------------------------------------------------------------
        stage = domain.quantities_to_be_monitored['stage']
        assert stage['min'] <= stage['max']

        assert num.allclose(stage['min'], initial_runup_height,
                            rtol=1.0/N)    # First order accuracy

        depth = domain.quantities_to_be_monitored['stage-elevation']
        assert depth['min'] <= depth['max']
        assert depth['min'] >= 0.0
        assert depth['max'] >= 0.0

        #--------------------------------------------------------------
        # Update boundary to allow inflow
        #--------------------------------------------------------------
        domain.set_boundary({'right': Bd})

        #--------------------------------------------------------------
        # Evolve system through time
        #--------------------------------------------------------------
        for t in domain.evolve(yieldstep=0.1, finaltime=3.0):
            domain.quantity_statistics()    # Run it silently

        #--------------------------------------------------------------
        # Test extrema again
        #--------------------------------------------------------------
        stage = domain.quantities_to_be_monitored['stage']
        assert stage['min'] <= stage['max']

        assert num.allclose(stage['min'], initial_runup_height,
                            rtol = 1.0/N) # First order accuracy

        depth = domain.quantities_to_be_monitored['stage-elevation']
        assert depth['min'] <= depth['max']
        assert depth['min'] >= 0.0
        assert depth['max'] >= 0.0

        # Cleanup
        os.remove(domain.get_name() + '.sww')





#################################################################################

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_swb_clean, 'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)

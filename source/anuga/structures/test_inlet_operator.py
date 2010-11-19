#!/usr/bin/env python


import unittest
import os.path
import sys

import numpy
import anuga

from anuga.utilities.system_tools import get_pathname_from_package
from anuga.structures.boyd_box_operator import Boyd_box_operator
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
from anuga.shallow_water.shallow_water_domain import Domain
from anuga.shallow_water.forcing import Rainfall, Inflow

from anuga.structures.inlet_operator import Inlet_operator

class Test_inlet_operator(unittest.TestCase):
    """
	Test the boyd box operator, in particular the discharge_routine!
    """

    def setUp(self):
        pass

    def tearDown(self):
        pass
    
    
    def _create_domain(self,d_length,
                            d_width,
                            dx,
                            dy,
                            elevation_0,
                            elevation_1,
                            stage_0,
                            stage_1):
        
        points, vertices, boundary = rectangular_cross(int(d_length/dx), int(d_width/dy),
                                                        len1=d_length, len2=d_width)
        domain = Domain(points, vertices, boundary)   
        domain.set_name('Test_Outlet_Inlet')                 # Output name
        domain.set_store()
        domain.set_default_order(2)
        domain.H0 = 0.01
        domain.tight_slope_limiters = 1

        #print 'Size', len(domain)

        #------------------------------------------------------------------------------
        # Setup initial conditions
        #------------------------------------------------------------------------------

        def elevation(x, y):
            """Set up a elevation
            """
            
            z = numpy.zeros(x.shape,dtype='d')
            z[:] = elevation_0
            
            numpy.putmask(z, x > d_length/2, elevation_1)
    
            return z
            
        def stage(x,y):
            """Set up stage
            """
            z = numpy.zeros(x.shape,dtype='d')
            z[:] = stage_0
            
            numpy.putmask(z, x > d_length/2, stage_1)

            return z
            
        #print 'Setting Quantities....'
        domain.set_quantity('elevation', elevation)  # Use function for elevation
        domain.set_quantity('stage',  stage)   # Use function for elevation

        Br = anuga.Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})
        
        return domain

    def test_inlet_Q(self):
        """test_inlet_Q
        
        This tests that the inlet operator adds the correct amount of water
        """

        stage_0 = 11.0
        stage_1 = 10.0
        elevation_0 = 10.0
        elevation_1 = 10.0

        domain_length = 200.0
        domain_width = 200.0
        

        domain = self._create_domain(d_length=domain_length,
                                     d_width=domain_width,
                                     dx = 10.0,
                                     dy = 10.0,
                                     elevation_0 = elevation_0,
                                     elevation_1 = elevation_1,
                                     stage_0 = stage_0,
                                     stage_1 = stage_1)

        vol0 = domain.compute_total_volume()

        finaltime = 3.0
        line1 = [[95.0, 10.0], [105.0, 10.0]]
        Q1 = 5.00
        
        line2 = [[10.0, 90.0], [20.0, 90.0]]
        Q2 = 10.0
        
        Inlet_operator(domain, line1, Q1)
        Inlet_operator(domain, line2, Q2)

        for t in domain.evolve(yieldstep = 1.0, finaltime = finaltime):
            #domain.write_time()
            #print domain.volumetric_balance_statistics()
            pass
 

        vol1 = domain.compute_total_volume()

        assert numpy.allclose((Q1+Q2)*finaltime, vol1-vol0, rtol=1.0e-8) 
        
        


# =========================================================================
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_inlet_operator, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)

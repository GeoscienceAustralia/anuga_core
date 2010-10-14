#!/usr/bin/env python


import unittest
import os.path
import sys

from anuga.utilities.system_tools import get_pathname_from_package
from anuga.structures.boyd_box_operator import Boyd_box_operator
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
from anuga.shallow_water.shallow_water_domain import Domain
from anuga.shallow_water.forcing import Rainfall, Inflow
import numpy


class Test_boyd_box_operator(unittest.TestCase):
    """
	Test the boyd box operator, in particular the discharge_routine!
    """

    def setUp(self):
        pass

    def tearDown(self):
        pass


    def test_boyd_1(self):
        """test_boyd_1
        
        This tests the Boyd routine with data obtained from culvertw application 1.1 by IceMindserer  BD Parkinson, 
        calculation code by MJ Boyd 
        """

        stage_0 = 11.0
        stage_1 = 10.0
        elevation_0 = 10.0
        elevation_1 = 10.0

        culvert_length = 20.0
        culvert_width = 3.66
        culvert_height = 3.66
        culvert_losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
        culvert_mannings = 0.013
        
        culvert_apron = 0.001
        enquiry_gap = 1.0

        domain_length = 200.  #x-Dir
        domain_width  = 200.  #y-dir
        dx = dy = 10.0          # Resolution: Length of subdivisions on both axes


        points, vertices, boundary = rectangular_cross(int(domain_length/dx), int(domain_width/dy),
                                                        len1=domain_length, len2=domain_width)
        domain = Domain(points, vertices, boundary)   
        domain.set_name('Test_Outlet_Inlet')                 # Output name
        domain.set_default_order(2)
        domain.H0 = 0.01
        domain.tight_slope_limiters = 1

        print 'Size', len(domain)

        #------------------------------------------------------------------------------
        # Setup initial conditions
        #------------------------------------------------------------------------------

        def elevation(x, y):
            """Set up a elevation
            """
            
            z = numpy.zeros(x.shape,dtype='d')
            z[:] = elevation_0
            
            numpy.putmask(z, x > domain_length/2, elevation_1)
    
            return z
            
        def stage(x,y):
            """Set up stage
            """
            z = numpy.zeros(x.shape,dtype='d')
            z[:] = stage_0
            
            numpy.putmask(z, x > domain_length/2, stage_1)

            return z
            
        print 'Setting Quantities....'
        domain.set_quantity('elevation', elevation)  # Use function for elevation
        domain.set_quantity('stage',  stage)   # Use function for elevation


        print 'Defining Structures'
        
        culvert = Boyd_box_operator(domain,
            label='3.6x3.6RCBC',
            end_point0=[domain_length/2-culvert_length/2, 100.0],
            end_point1=[domain_length/2+culvert_length/2, 100.0],
            losses=culvert_losses,
            width=culvert_width,
            height=culvert_height,
            apron=culvert_apron,
            enquiry_gap=enquiry_gap,
            use_momentum_jet=False,
            use_velocity_head=False,
            manning=culvert_mannings,
            verbose=False)


        culvert.determine_inflow_outflow()
        
        ( Q, v, d ) = culvert.discharge_routine()
        

        print Q
        print v
        print d
        

        assert numpy.allclose(Q, 4.55, rtol=1.0e-1) #inflow
        assert numpy.allclose(v, 2.3,  rtol=1.0e-1) #outflow velocity
        assert numpy.allclose(d, 0.54, rtol=1.0e-1) #depth at outlet used to calc v 
        

# =========================================================================
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_boyd_box_operator, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)

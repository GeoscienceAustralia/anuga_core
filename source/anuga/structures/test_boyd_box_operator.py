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
        
        return domain

    def test_boyd_non_skew(self):
        """test_boyd_non_skew
        
        This tests the Boyd routine with data obtained from culvertw application 1.1 by IceMindserer  BD Parkinson, 
        calculation code by MJ Boyd 
        """

        stage_0 = 11.0
        stage_1 = 10.0
        elevation_0 = 10.0
        elevation_1 = 10.0

        domain_length = 200.0
        domain_width = 200.0

        culvert_length = 20.0
        culvert_width = 3.66
        culvert_height = 3.66
        culvert_losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
        culvert_mannings = 0.013
        
        culvert_apron = 0.0
        enquiry_gap = 10.0

        
        expected_Q = 6.23
        expected_v = 2.55
        expected_d = 0.66
        

        domain = self._create_domain(d_length=domain_length,
                                     d_width=domain_width,
                                     dx = 5.0,
                                     dy = 5.0,
                                     elevation_0 = elevation_0,
                                     elevation_1 = elevation_1,
                                     stage_0 = stage_0,
                                     stage_1 = stage_1)
 

        #print 'Defining Structures'
        
        ep0 = numpy.array([domain_length/2-culvert_length/2, 100.0])
        ep1 = numpy.array([domain_length/2+culvert_length/2, 100.0])
        
        
        culvert = Boyd_box_operator(domain,
                                    losses=culvert_losses,
                                    width=culvert_width,
                                    end_points=[ep0, ep1],
                                    height=culvert_height,
                                    apron=culvert_apron,
                                    enquiry_gap=enquiry_gap,
                                    use_momentum_jet=False,
                                    use_velocity_head=False,
                                    manning=culvert_mannings,
                                    label='3.6x3.6RCBC',
                                    verbose=False)

        #culvert.determine_inflow_outflow()
        
        ( Q, v, d ) = culvert.discharge_routine()
        
        #print 'test_boyd_non_skew'
        #print 'Q: ', Q, 'expected_Q: ', expected_Q
        

        assert numpy.allclose(Q, expected_Q, rtol=1.0e-2) #inflow
        assert numpy.allclose(v, expected_v, rtol=1.0e-2) #outflow velocity
        assert numpy.allclose(d, expected_d, rtol=1.0e-2) #depth at outlet used to calc v 
        
        
    def test_boyd_skew(self):
        """test_boyd_skew
        
        This tests the Boyd routine with data obtained from culvertw application 1.1 by IceMindserer  BD Parkinson, 
        calculation code by MJ Boyd 
        """

        stage_0 = 11.0
        stage_1 = 10.0
        elevation_0 = 10.0
        elevation_1 = 10.0

        domain_length = 200.0
        domain_width = 200.0

        culvert_length = 20.0
        culvert_width = 3.66
        culvert_height = 3.66
        culvert_losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
        culvert_mannings = 0.013
        
        culvert_apron = 0.0
        enquiry_gap = 10.0
        
        expected_Q = 6.23
        expected_v = 2.55
        expected_d = 0.66
        

        domain = self._create_domain(d_length=domain_length,
                                     d_width=domain_width,
                                     dx = 5.0,
                                     dy = 5.0,
                                     elevation_0 = elevation_0,
                                     elevation_1 = elevation_1,
                                     stage_0 = stage_0,
                                     stage_1 = stage_1)

        #print 'Defining Structures'
        
        a = domain_length/2 - culvert_length/2
        b = domain_length/2 + culvert_length/2
        
        el0 = numpy.array([[a, 100.0 - culvert_width/2], [a, 100.0 + culvert_width/2]])
        el1 = numpy.array([[b, 100.0 - culvert_width/2], [b, 100.0 + culvert_width/2]])
        
        culvert = Boyd_box_operator(domain,
                                    losses=culvert_losses,
                                    width=culvert_width,
                                    exchange_lines=[el0, el1],
                                    height=culvert_height,
                                    apron=culvert_apron,
                                    enquiry_gap=enquiry_gap,
                                    use_momentum_jet=False,
                                    use_velocity_head=False,
                                    manning=culvert_mannings,
                                    label='3.6x3.6RCBC',
                                    verbose=False)

        #culvert.determine_inflow_outflow()
        
        ( Q, v, d ) = culvert.discharge_routine()
        
        #print 'test_boyd_skew'
        #print 'Q: ', Q, 'expected_Q: ', expected_Q

        assert numpy.allclose(Q, expected_Q, rtol=1.0e-2) #inflow
        assert numpy.allclose(v, expected_v, rtol=1.0e-2) #outflow velocity
        assert numpy.allclose(d, expected_d, rtol=1.0e-2) #depth at outlet used to calc v         


    def test_boyd_non_skew_enquiry_points(self):
        """test_boyd_skew_enquiry_points
        
        This tests the Boyd routine with data obtained from culvertw application 1.1 by IceMindserer  BD Parkinson, 
        calculation code by MJ Boyd 
        """

        stage_0 = 11.0
        stage_1 = 10.0
        elevation_0 = 10.0
        elevation_1 = 10.0

        domain_length = 200.0
        domain_width = 200.0

        culvert_length = 20.0
        culvert_width = 3.66
        culvert_height = 3.66
        culvert_losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
        culvert_mannings = 0.013
        
        culvert_apron = 0.0
        enquiry_gap = 10.0
        
        expected_Q = 6.23
        expected_v = 2.55
        expected_d = 0.66
        

        # Probably no need to change below here
        
        domain = self._create_domain(d_length=domain_length,
                                     d_width=domain_width,
                                     dx = 5.0,
                                     dy = 5.0,
                                     elevation_0 = elevation_0,
                                     elevation_1 = elevation_1,
                                     stage_0 = stage_0,
                                     stage_1 = stage_1)


        #print 'Defining Structures'
        
        a = domain_length/2 - culvert_length/2
        b = domain_length/2 + culvert_length/2
        
        el0 = numpy.array([[a, 100.0 - culvert_width/2], [a, 100.0 + culvert_width/2]])
        el1 = numpy.array([[b, 100.0 - culvert_width/2], [b, 100.0 + culvert_width/2]])
        
        enquiry_points = (numpy.array([85, 100]), numpy.array([115, 100]))
        
        culvert = Boyd_box_operator(domain,
                                    losses=culvert_losses,
                                    width=culvert_width,
                                    exchange_lines=[el0, el1],
                                    enquiry_points=enquiry_points,
                                    height=culvert_height,
                                    apron=culvert_apron,
                                    enquiry_gap=enquiry_gap,
                                    use_momentum_jet=False,
                                    use_velocity_head=False,
                                    manning=culvert_mannings,
                                    label='3.6x3.6RCBC',
                                    verbose=False)

        #culvert.determine_inflow_outflow()
        
        ( Q, v, d ) = culvert.discharge_routine()
        
        #print 'test_boyd_non_skew_enquiry_points'
        #print 'Q: ', Q, 'expected_Q: ', expected_Q
        #print 'v: ', v, 'expected_v: ', expected_v
        #print 'd: ', d, 'expected_d: ', expected_d

        assert numpy.allclose(Q, expected_Q, rtol=1.0e-2) #inflow
        assert numpy.allclose(v, expected_v, rtol=1.0e-2) #outflow velocity
        assert numpy.allclose(d, expected_d, rtol=1.0e-2) #depth at outlet used to calc v         


# =========================================================================
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_boyd_box_operator, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)

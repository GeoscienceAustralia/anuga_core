#!/usr/bin/env python


import unittest
import os.path
import sys

import numpy
import anuga

from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
from anuga.shallow_water.shallow_water_domain import Domain
from anuga.abstract_2d_finite_volumes.util import file_function
from anuga.utilities.system_tools import get_pathname_from_package

from anuga.structures.inlet_operator import Inlet_operator

import warnings
warnings.simplefilter("ignore")

class Test_inlet_operator(unittest.TestCase):
    """
	Test the boyd box operator, in particular the discharge_routine!
    """

    def setUp(self):
        pass

    def tearDown(self):
        try:
            os.remove('Test_Outlet_Inlet.sww')
        except:
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

    def test_inlet_constant_Q(self):
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
        
        Inlet_operator(domain, line1, Q1, logging=False)
        Inlet_operator(domain, line2, Q2)

        for t in domain.evolve(yieldstep = 1.0, finaltime = finaltime):
            #domain.write_time()
            #print domain.volumetric_balance_statistics()
            pass
 

        vol1 = domain.compute_total_volume()

        assert numpy.allclose((Q1+Q2)*finaltime, vol1-vol0, rtol=1.0e-8) 
        assert numpy.allclose((Q1+Q2)*finaltime, domain.fractional_step_volume_integral, rtol=1.0e-8) 



    def test_inlet_constant_Q_polygon(self):
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
        poly1 = [[95.0, 10.0], [105.0, 10.0], [105, 20.0], [95.0, 20.0]]
        Q1 = 5.00


        Inlet_operator(domain, poly1, Q1, logging=False)


        for t in domain.evolve(yieldstep = 1.0, finaltime = finaltime):
            #domain.write_time()
            #print domain.volumetric_balance_statistics()
            pass


        vol1 = domain.compute_total_volume()

        assert numpy.allclose((Q1)*finaltime, vol1-vol0, rtol=1.0e-8)
        assert numpy.allclose((Q1)*finaltime, domain.fractional_step_volume_integral, rtol=1.0e-8) 



    def test_inlet_variable_Q(self):
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

        #Make sure we are inthe right directory to find the
        #time series data for the inlets
        import os
        
        path = get_pathname_from_package('anuga.structures')
        filename1 = os.path.join(path, 'tests', 'data', 'inlet_operator_test1.tms')
        filename2 = os.path.join(path, 'tests', 'data', 'inlet_operator_test2.tms')

        line1 = [[95.0, 10.0], [105.0, 10.0]]
        Q1 = file_function(filename=filename1, quantities=['hydrograph'])
        
        line2 = [[10.0, 90.0], [20.0, 90.0]]
        Q2 = file_function(filename=filename2, quantities=['hydrograph'])

        
        Inlet_operator(domain, line1, Q1)
        Inlet_operator(domain, line2, Q2)

        for t in domain.evolve(yieldstep = 1.0, finaltime = finaltime):
            #domain.write_time()
            #print domain.volumetric_balance_statistics()
            pass
 

        vol1 = domain.compute_total_volume()

        #print vol1-vol0
        
        assert numpy.allclose(13.5, vol1-vol0, rtol=1.0e-8) 
        assert numpy.allclose(vol1-vol0, domain.fractional_step_volume_integral, rtol=1.0e-8) 
                
    def test_inlet_variable_Q_default(self):
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

        finaltime = 5.0

        #Make sure we are inthe right directory to find the
        #time series data for the inlets
        import os
        baseDir = os.getcwd()

        path = get_pathname_from_package('anuga.structures')
        filename1 = os.path.join(path, 'tests', 'data', 'inlet_operator_test1.tms')
        filename2 = os.path.join(path, 'tests', 'data', 'inlet_operator_test2.tms')

        line1 = [[95.0, 10.0], [105.0, 10.0]]
        Q1 = file_function(filename=filename1, quantities=['hydrograph'])

        line2 = [[10.0, 90.0], [20.0, 90.0]]
        Q2 = file_function(filename=filename2, quantities=['hydrograph'])

        os.chdir(baseDir)

        Inlet_operator(domain, line1, Q1, default=6)
        Inlet_operator(domain, line2, Q2, default=3)

        for t in domain.evolve(yieldstep = 1.0, finaltime = finaltime):
            #domain.write_time()
            #print domain.volumetric_balance_statistics()
            pass


        vol1 = domain.compute_total_volume()

        #print vol1-vol0

        assert numpy.allclose(31.5, vol1-vol0, rtol=1.0e-8)
        assert numpy.allclose(vol1-vol0, domain.fractional_step_volume_integral, rtol=1.0e-8) 

# =========================================================================
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_inlet_operator, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)

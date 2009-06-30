#!/usr/bin/env python


import unittest
import os.path
import sys

from anuga.utilities.system_tools import get_pathname_from_package
from anuga.culvert_flows.culvert_routines import boyd_generalised_culvert_model
import numpy as num


class Test_culvert_routines_pipe_10pct(unittest.TestCase):
    """
	This unit test sets up 6 tests for various culvert conditions for a Circular Pipe Culvert on a 10% Slope
    """

    def setUp(self):
        pass

    def tearDown(self):
        pass


    def test_boyd_1(self):
        """test_boyd_1
        
        This tests the Boyd routine with data obtained from ??? by Petar Milevski    
        """
        # FIXME(Ole): This test fails (20 Feb 2009)

        g=9.81


        inlet_depth=0.150
        outlet_depth=0.15
        inlet_velocity=1.00
        outlet_velocity=0.5
        
        culvert_length=10.0
        culvert_width=0.0
        culvert_height=1.20
        
        culvert_type='circle'
        manning=0.013
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 
        culvert_slope=10.0  # % Downward
        z_in = 10.0
        z_out = -culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out
        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 

        Q, v, d = boyd_generalised_culvert_model(inlet_depth,
                                                 outlet_depth,
                                                 inlet_velocity,
                                                 outlet_velocity,
                                                 inlet_specific_energy, 
                                                 delta_total_energy, 
                                                 g,
                                                 culvert_length,
                                                 culvert_width,
                                                 culvert_height,
                                                 culvert_type,
                                                 manning,
                                                 sum_loss)
        
        #print ('%s,%.2f,%.2f,%.2f' %('ANUGAcalcsTEST01 Q-v-d',Q,v,d))
        #print('%s,%.2f,%.2f,%.2f' %('Spreadsheet_Boydcalcs', 0.113, 0.297, 0.443))
        assert num.allclose(Q, 0.113, rtol=1.0e-1) #inflow
        assert num.allclose(v, 0.297, rtol=1.0e-1) #outflow velocity
        assert num.allclose(d, 0.443, rtol=1.0e-1) #depth at outlet used to calc v 
        
    def test_boyd_2(self):
        """test_boyd_2
        
        This tests the Boyd routine with data obtained from ??? by Petar Milevski    
        """
        # FIXME(Ole): This test fails (20 Feb 2009)

        g=9.81
        culvert_slope=10  # Downward

        inlet_depth=0.500
        outlet_depth=0.700
        inlet_velocity=1.0
        outlet_velocity=0.50
        
        culvert_length=10.0
        culvert_width=0.0
        culvert_height=1.20
        culvert_width=0.0
        culvert_type='circle'
        manning=0.013
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 
        z_in = 0.0
        z_out = -culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out

        Q, v, d = boyd_generalised_culvert_model(inlet_depth,
                                                 outlet_depth,
                                                 inlet_velocity,
                                                 outlet_velocity,
                                                 inlet_specific_energy, 
                                                 delta_total_energy, 
                                                 g,
                                                 culvert_length,
                                                 culvert_width,
                                                 culvert_height,
                                                 culvert_type,
                                                 manning,
                                                 sum_loss)
        
        #print ('%s,%.2f,%.2f,%.2f' %('ANUGAcalcsTEST02 Q-v-d',Q,v,d))
        #print ('%s,%.2f,%.2f,%.2f' %('Spreadsheet_Boydcalcs', 0.585, 0.873, 0.687))
        assert num.allclose(Q, 0.585, rtol=1.0e-1) #inflow
        assert num.allclose(v, 0.873, rtol=1.0e-1) #outflow velocity
        assert num.allclose(d, 0.687, rtol=1.0e-1) #depth at outlet used to calc v  

    def test_boyd_3(self):
        """test_boyd_3
        
        This tests the Boyd routine with data obtained from ??? by Petar Milevski    
        """
        # FIXME(Ole): This test fails (20 Feb 2009)

        g=9.81
        culvert_slope=10  # Downward

        inlet_depth=1.800
        outlet_depth=0.80
        inlet_velocity=1.0
        outlet_velocity=0.5
        
        culvert_length=10.0
        
        culvert_height=1.20
        culvert_width=0.0
        culvert_type='circle'
        manning=0.013
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 
        z_in = 0.0
        z_out = -culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out

        Q, v, d = boyd_generalised_culvert_model(inlet_depth,
                                                 outlet_depth,
                                                 inlet_velocity,
                                                 outlet_velocity,
                                                 inlet_specific_energy, 
                                                 delta_total_energy, 
                                                 g,
                                                 culvert_length,
                                                 culvert_width,
                                                 culvert_height,
                                                 culvert_type,
                                                 manning,
                                                 sum_loss)
        #print ('%s,%.2f'%('SPEC_E = ',inlet_specific_energy))
        #print ('%s,%.2f'%('Delta E = ',delta_total_energy))
        #print ('%s,%.2f,%.2f,%.2f' %('ANUGAcalcsTEST03 Q-v-d',Q,v,d))
        #print ('%s,%.2f,%.2f,%.2f' %('Spreadsheet_Boydcalcs', 3.441, 3.042, 1.2))
        assert num.allclose(Q, 3.441, rtol=1.0e-2) #inflow
        assert num.allclose(v, 3.042, rtol=1.0e-2) #outflow velocity
        assert num.allclose(d, 1.2, rtol=1.0e-2) #depth at outlet used to calc v

#NOTE FROM HERE DOWN THE UNITS TEST HAVE NOT BEEN AMENDED TO ALLOW VELOCITY COMPONENT TO BE USED. ONLY ABOVE 3 TESTS WORK. PM WILL FIX THE ONES BELOW WHEN THE ABOVE 3 ARE WORKING
    def test_boyd_4(self):
        """test_boyd_4
        
        This tests the Boyd routine with data obtained from ??? by Petar Milevski    
        """
        # FIXME(Ole): This test fails (20 Feb 2009)

        g=9.81
        culvert_slope=10  # Downward

        inlet_depth=1.00
        outlet_depth=0.8
        inlet_velocity=1.0
        outlet_velocity=0.5 
        culvert_length=10.0
        culvert_width=0.0
        culvert_height=1.20
       
        culvert_type='circle'
        manning=0.013
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 
        z_in = 10.0
        z_out = 10.0-culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out



        Q, v, d = boyd_generalised_culvert_model(inlet_depth,
                                                 outlet_depth,
                                                 inlet_velocity,
                                                 outlet_velocity,
                                                 inlet_specific_energy, 
                                                 delta_total_energy, 
                                                 g,
                                                 culvert_length,
                                                 culvert_width,
                                                 culvert_height,
                                                 culvert_type,
                                                 manning,
                                                 sum_loss)        
        #print ('%s,%.2f'%('SPEC_E = ',inlet_specific_energy))
        #print ('%s,%.2f'%('Delta E = ',delta_total_energy))
        #print ('%s,%.2f,%.2f,%.2f' %('ANUGAcalcsTEST04 Q-v-d',Q,v,d))
        #print ('%s,%.2f,%.2f,%.2f' %('Spreadsheet_Boydcalcs', 1.676, 1.821, 0.91))
        assert num.allclose(Q, 1.676, rtol=1.0e-2) #inflow
        assert num.allclose(v, 1.821, rtol=1.0e-2) #outflow velocity
        assert num.allclose(d, 0.91, rtol=1.0e-2) #depth at outlet used to calc v

    def test_boyd_5(self):
        """test_boyd_5
        
        This tests the Boyd routine with data obtained from ??? by Petar Milevski    
        """
        # FIXME(Ole): This test fails (20 Feb 2009)

        g=9.81
        culvert_slope=10  # Downward

        inlet_depth=1.50
        inlet_velocity= 1.0
        outlet_depth=2.5
        outlet_velocity=0.5
        culvert_length=10.0
        culvert_width=0.0
        culvert_height=1.20
        
        culvert_type='circle'
        manning=0.013
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 
        z_in = 10.0
        z_out = 10.0-culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out



        Q, v, d = boyd_generalised_culvert_model(inlet_depth,
                                                 outlet_depth,
                                                 inlet_velocity,
                                                 outlet_velocity,
                                                 inlet_specific_energy, 
                                                 delta_total_energy, 
                                                 g,
                                                 culvert_length,
                                                 culvert_width,
                                                 culvert_height,
                                                 culvert_type,
                                                 manning,
                                                 sum_loss)        
        #print ('%s,%.3f'%('SPEC_E = ',inlet_specific_energy))
        #print ('%s,%.3f'%('Delta E = ',delta_total_energy))
        
        #print ('%s,%.3f,%.3f,%.3f' %('ANUGAcalcsTEST05 Q-v-d',Q,v,d))
        #print ('%s,%.3f,%.3f,%.3f' %('Spreadsheet_Boydcalcs',0.759, 0.671, 1.20))
        assert num.allclose(Q, 0.759, rtol=1.0e-2) #inflow
        assert num.allclose(v, 0.671, rtol=1.0e-2) #outflow velocity
        assert num.allclose(d, 1.20, rtol=1.0e-2) #depth at outlet used to calc v

    def test_boyd_6(self):
        """test_boyd_6
        
        This tests the Boyd routine with data obtained from ??? by Petar Milevski    
        """
        # FIXME(Ole): This test fails (20 Feb 2009)

        g=9.81
        culvert_slope=10  # Downward

        inlet_depth=1.50
        inlet_velocity= 4.0
        outlet_depth=0.80
        outlet_velocity=4.0
        culvert_length=10.0
        culvert_height=1.20
        culvert_width=0.0
        culvert_type='circle'
        manning=0.013
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 
        z_in = 10.0
        z_out = 10.0-culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out



        Q, v, d = boyd_generalised_culvert_model(inlet_depth,
                                                 outlet_depth,
                                                 inlet_velocity,
                                                 outlet_velocity,
                                                 inlet_specific_energy, 
                                                 delta_total_energy, 
                                                 g,
                                                 culvert_length,
                                                 culvert_width,
                                                 culvert_height,
                                                 culvert_type,
                                                 manning,
                                                 sum_loss)        
        #print ('%s,%.3f,%.3f,%.3f' %('ANUGAcalcsTEST06 Q-v-d',Q,v,d))
        #print ('%s,%.3f,%.3f,%.3f' %('Spreadsheet_Boydcalcs',3.962, 3.503, 1.20))
        assert num.allclose(Q, 3.962, rtol=1.0e-2) #inflow
        assert num.allclose(v, 3.503, rtol=1.0e-2) #outflow velocity
        assert num.allclose(d, 1.20, rtol=1.0e-2) #depth at outlet used to calc v
# =========================================================================
# =========================================================================

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_culvert_routines_pipe_10pct, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)

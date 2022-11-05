#!/usr/bin/env python

import unittest
import os.path
import sys

from anuga.utilities.system_tools import get_pathname_from_package
from anuga.culvert_flows.culvert_routines import boyd_generalised_culvert_model
import numpy as num


class Test_culvert_routines(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass


    def test_boyd_0(self):
        """test_boyd_0
        
        This tests the Boyd routine with data obtained from ??? by Petar Milevski
        This test is the only one that passed in late February 2009
        """
      
        g=9.81
        culvert_slope=0.1  # Downward

        inlet_depth=2.0
        outlet_depth=0.0
        
        inlet_velocity=0.0,
        outlet_velocity=0.0,        

        culvert_length=4.0
        culvert_width=1.2
        culvert_height=0.75

        culvert_type='box'
        manning=0.013
        sum_loss=0.0

        inlet_specific_energy=inlet_depth #+0.5*v**2/g 
        z_in = 0.0
        z_out = -culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth # + 
        E_out = z_out+outlet_depth # +
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
        
        #print Q, v, d
        assert num.allclose(Q, 3.118, rtol=1.0e-3)
        

        #assert num.allclose(v, 0.93)
        #assert num.allclose(d, 0.0)
        

    def Xtest_boyd_00(self):
        """test_boyd_00
        
        This tests the Boyd routine with data obtained from ??? by Petar Milevski    
        """
        # FIXME(Ole): This test fails (20 Feb 2009)

        g=9.81
        culvert_slope=0.1  # Downward

        inlet_depth=0.2
        outlet_depth=0.0

        inlet_velocity=0.0,
        outlet_velocity=0.0,                
        
        culvert_length=4.0
        culvert_width=1.2
        culvert_height=0.75

        culvert_type='box'
        manning=0.013
        sum_loss=0.0

        inlet_specific_energy=inlet_depth #+0.5*v**2/g 
        z_in = 0.0
        z_out = -culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth # + 
        E_out = z_out+outlet_depth # +
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
        
        #print Q, v, d
        assert num.allclose(Q, 0.185, rtol=1.0e-3)
        #assert num.allclose(v, 0.93)
        #assert num.allclose(d, 0.0)
        
    def Xtest_boyd_1(self):
        """test_boyd_1
        
        This tests the Boyd routine with data obtained from ??? by Petar Milevski    
        """
        # FIXME(Ole): This test fails (20 Feb 2009)

        g=9.81
        culvert_slope=0.01  # Downward

        inlet_depth=0.263
        outlet_depth=0.0

        culvert_length=4.0
        culvert_width=0.75
        culvert_height=0.75
        
        culvert_type='pipe'
        manning=0.013
        sum_loss=1.5

        inlet_specific_energy=inlet_depth #+0.5*v**2/g 
        z_in = 0.0
        z_out = -culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth  #+ 0.5*v**2/g
        E_out = z_out+outlet_depth  #+ 0.5*v**2/g
        delta_total_energy = E_in-E_out

        Q, v, d = boyd_generalised_culvert_model(inlet_depth, 
                                                 outlet_depth,
                                                 inlet_specific_energy, 
                                                 delta_total_energy, 
                                                 g,
                                                 culvert_length,
                                                 culvert_width,
                                                 culvert_height,
                                                 culvert_type,
                                                 manning,
                                                 sum_loss)
        
        print(Q, v, d)
        assert num.allclose(Q, 0.10, rtol=1.0e-2) #inflow
        assert num.allclose(v, 1.13, rtol=1.0e-2) #outflow velocity
        assert num.allclose(d, 0.15, rtol=1.0e-2) #depth at outlet used to calc v 
        
    def Xtest_boyd_2(self):
        """test_boyd_2
        
        This tests the Boyd routine with data obtained from ??? by Petar Milevski    
        """
        # FIXME(Ole): This test fails (20 Feb 2009)

        g=9.81
        culvert_slope=0.01  # Downward

        inlet_depth=1.135
        outlet_depth=0.0

        culvert_length=4.0
        culvert_width=0.75
        culvert_height=0.75
        
        culvert_type='pipe'
        manning=0.013
        sum_loss=1.5

        inlet_specific_energy=inlet_depth #+0.5*v**2/g 
        z_in = 0.0
        z_out = -culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth  #+ 0.5*v**2/g
        E_out = z_out+outlet_depth  #+ 0.5*v**2/g
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
        
        print(Q, v, d)
        assert num.allclose(Q, 1.00, rtol=1.0e-2) #inflow
        assert num.allclose(v, 2.59, rtol=1.0e-2) #outflow velocity
        assert num.allclose(d, 0.563, rtol=1.0e-2) #depth at outlet used to calc v  

    def Xtest_boyd_3(self):
        """test_boyd_3
        
        This tests the Boyd routine with data obtained from ??? by Petar Milevski    
        """
        # FIXME(Ole): This test fails (20 Feb 2009)

        g=9.81
        culvert_slope=0.01  # Downward

        inlet_depth=12.747
        outlet_depth=0.0

        culvert_length=4.0
        culvert_width=0.75
        culvert_height=0.75
        
        culvert_type='pipe'
        manning=0.013
        sum_loss=1.5

        inlet_specific_energy=inlet_depth #+0.5*v**2/g 
        z_in = 0.0
        z_out = -culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth  #+ 0.5*v**2/g
        E_out = z_out+outlet_depth  #+ 0.5*v**2/g
        delta_total_energy = E_in-E_out

        Q, v, d = boyd_generalised_culvert_model(inlet_depth, 
                                                 outlet_depth,
                                                 inlet_specific_energy, 
                                                 delta_total_energy, 
                                                 g,
                                                 culvert_length,
                                                 culvert_width,
                                                 culvert_height,
                                                 culvert_type,
                                                 manning,
                                                 sum_loss)
        
        print(Q, v, d)
        assert num.allclose(Q, 5.00, rtol=1.0e-2) #inflow
        assert num.allclose(v, 11.022, rtol=1.0e-2) #outflow velocity
        assert num.allclose(d, 0.72, rtol=1.0e-2) #depth at outlet used to calc v

    def Xtest_boyd_4(self):
        """test_boyd_4
        
        This tests the Boyd routine with data obtained from ??? by Petar Milevski    
        """
        # FIXME(Ole): This test fails (20 Feb 2009)

        g=9.81
        culvert_slope=0.01  # Downward

        inlet_depth=1.004
        outlet_depth=1.00

        culvert_length=4.0
        culvert_width=0.75
        culvert_height=0.75
        
        culvert_type='pipe'
        manning=0.013
        sum_loss=1.5

        inlet_specific_energy=inlet_depth #+0.5*v**2/g 
        z_in = 0.0
        z_out = -culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth  #+ 0.5*v**2/g
        E_out = z_out+outlet_depth  #+ 0.5*v**2/g
        delta_total_energy = E_in-E_out

        Q, v, d = boyd_generalised_culvert_model(inlet_depth, 
                                                 outlet_depth,
                                                 inlet_specific_energy, 
                                                 delta_total_energy, 
                                                 g,
                                                 culvert_length,
                                                 culvert_width,
                                                 culvert_height,
                                                 culvert_type,
                                                 manning,
                                                 sum_loss)
        
        print(Q, v, d)
        assert num.allclose(Q, 0.10, rtol=1.0e-2) #inflow
        assert num.allclose(v, 0.22, rtol=1.0e-2) #outflow velocity
        assert num.allclose(d, 0.76, rtol=1.0e-2) #depth at outlet used to calc v

    def Xtest_boyd_5(self):
        """test_boyd_5
        
        This tests the Boyd routine with data obtained from ??? by Petar Milevski    
        """
        # FIXME(Ole): This test fails (20 Feb 2009)

        g=9.81
        culvert_slope=0.01  # Downward

        inlet_depth=1.401
        outlet_depth=1.00

        culvert_length=4.0
        culvert_width=0.75
        culvert_height=0.75
        
        culvert_type='pipe'
        manning=0.013
        sum_loss=1.5

        inlet_specific_energy=inlet_depth #+0.5*v**2/g 
        z_in = 0.0
        z_out = -culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth  #+ 0.5*v**2/g
        E_out = z_out+outlet_depth  #+ 0.5*v**2/g
        delta_total_energy = E_in-E_out

        Q, v, d = boyd_generalised_culvert_model(inlet_depth, 
                                                 outlet_depth,
                                                 inlet_specific_energy, 
                                                 delta_total_energy, 
                                                 g,
                                                 culvert_length,
                                                 culvert_width,
                                                 culvert_height,
                                                 culvert_type,
                                                 manning,
                                                 sum_loss)
        
        print(Q, v, d)
        assert num.allclose(Q, 1.00, rtol=1.0e-2) #inflow
        assert num.allclose(v, 2.204, rtol=1.0e-2) #outflow velocity
        assert num.allclose(d, 0.76, rtol=1.0e-2) #depth at outlet used to calc v


    def Xtest_boyd_6(self):
        """test_boyd_5
        
        This tests the Boyd routine with data obtained from ??? by Petar Milevski    
        """
        # FIXME(Ole): This test fails (20 Feb 2009)

        g=9.81
        culvert_slope=0.01  # Downward

        inlet_depth=12.747
        outlet_depth=1.00

        culvert_length=4.0
        culvert_width=0.75
        culvert_height=0.75
        
        culvert_type='pipe'
        manning=0.013
        sum_loss=1.5

        inlet_specific_energy=inlet_depth #+0.5*v**2/g 
        z_in = 0.0
        z_out = -culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth  #+ 0.5*v**2/g
        E_out = z_out+outlet_depth  #+ 0.5*v**2/g
        delta_total_energy = E_in-E_out

        Q, v, d = boyd_generalised_culvert_model(inlet_depth, 
                                                 outlet_depth,
                                                 inlet_specific_energy, 
                                                 delta_total_energy, 
                                                 g,
                                                 culvert_length,
                                                 culvert_width,
                                                 culvert_height,
                                                 culvert_type,
                                                 manning,
                                                 sum_loss)
        
        print(Q, v, d)
        assert num.allclose(Q, 5.00, rtol=1.0e-2) #inflow
        assert num.allclose(v, 11.022, rtol=1.0e-2) #outflow velocity
        assert num.allclose(d, 0.76, rtol=1.0e-2) #depth at outlet used to calc v


    def Xtest_boyd_7(self):
        """test_boyd_7
        
        This tests the Boyd routine with data obtained from ??? by Petar Milevski    
        """
        # FIXME(Ole): This test fails (20 Feb 2009)

        g=9.81
        culvert_slope=0.1  # Downward

        inlet_depth=0.303
        outlet_depth=0.00

        culvert_length=4.0
        culvert_width=0.75
        culvert_height=0.75
        
        culvert_type='pipe'
        manning=0.013
        sum_loss=1.5

        inlet_specific_energy=inlet_depth #+0.5*v**2/g 
        z_in = 0.0
        z_out = -culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth  #+ 0.5*v**2/g
        E_out = z_out+outlet_depth  #+ 0.5*v**2/g
        delta_total_energy = E_in-E_out

        Q, v, d = boyd_generalised_culvert_model(inlet_depth, 
                                                 outlet_depth,
                                                 inlet_specific_energy, 
                                                 delta_total_energy, 
                                                 g,
                                                 culvert_length,
                                                 culvert_width,
                                                 culvert_height,
                                                 culvert_type,
                                                 manning,
                                                 sum_loss)
        
        print(Q, v, d)
        assert num.allclose(Q, 0.10, rtol=1.0e-2) #inflow
        assert num.allclose(v, 1.13, rtol=1.0e-2) #outflow velocity
        assert num.allclose(d, 0.19, rtol=1.0e-2) #depth at outlet used to calc v


    def Xtest_boyd_8(self):
        """test_boyd_8
        
        This tests the Boyd routine with data obtained from ??? by Petar Milevski    
        """
        # FIXME(Ole): This test fails (20 Feb 2009)

        g=9.81
        culvert_slope=0.1  # Downward

        inlet_depth=1.135
        outlet_depth=0.00

        culvert_length=4.0
        culvert_width=0.75
        culvert_height=0.75
        
        culvert_type='pipe'
        manning=0.013
        sum_loss=1.5

        inlet_specific_energy=inlet_depth #+0.5*v**2/g 
        z_in = 0.0
        z_out = -culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth  #+ 0.5*v**2/g
        E_out = z_out+outlet_depth  #+ 0.5*v**2/g
        delta_total_energy = E_in-E_out

        Q, v, d = boyd_generalised_culvert_model(inlet_depth, 
                                                 outlet_depth,
                                                 inlet_specific_energy, 
                                                 delta_total_energy, 
                                                 g,
                                                 culvert_length,
                                                 culvert_width,
                                                 culvert_height,
                                                 culvert_type,
                                                 manning,
                                                 sum_loss)
        
        print(Q, v, d)
        assert num.allclose(Q, 1.00, rtol=1.0e-2) #inflow
        assert num.allclose(v, 2.204, rtol=1.0e-2) #outflow velocity
        assert num.allclose(d, 0.76, rtol=1.0e-2) #depth at outlet used to calc v

    def Xtest_boyd_9(self):
        """test_boyd_9
        
        This tests the Boyd routine with data obtained from ??? by Petar Milevski    
        """
        # FIXME(Ole): This test fails (20 Feb 2009)

        g=9.81
        culvert_slope=0.1  # Downward

        inlet_depth=1.1504
        outlet_depth=1.5

        culvert_length=4.0
        culvert_width=0.75
        culvert_height=0.75
        
        culvert_type='pipe'
        manning=0.013
        sum_loss=1.5

        inlet_specific_energy=inlet_depth #+0.5*v**2/g 
        z_in = 0.0
        z_out = -culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth  #+ 0.5*v**2/g
        E_out = z_out+outlet_depth  #+ 0.5*v**2/g
        delta_total_energy = E_in-E_out

        Q, v, d = boyd_generalised_culvert_model(inlet_depth, 
                                                 outlet_depth,
                                                 inlet_specific_energy, 
                                                 delta_total_energy, 
                                                 g,
                                                 culvert_length,
                                                 culvert_width,
                                                 culvert_height,
                                                 culvert_type,
                                                 manning,
                                                 sum_loss)
        
        print(Q, v, d)
        assert num.allclose(Q, 0.10, rtol=1.0e-2) #inflow
        assert num.allclose(v, 0.22, rtol=1.0e-2) #outflow velocity
        assert num.allclose(d, 0.76, rtol=1.0e-2) #depth at outlet used to calc v


    def Xtest_boyd_10(self):
        """test_boyd_9
        
        This tests the Boyd routine with data obtained from ??? by Petar Milevski    
        """
        # FIXME(Ole): This test fails (20 Feb 2009)

        g=9.81
        culvert_slope=0.1  # Downward

        inlet_depth=1.901
        outlet_depth=1.5

        culvert_length=4.0
        culvert_width=0.75
        culvert_height=0.75
        
        culvert_type='pipe'
        manning=0.013
        sum_loss=1.5

        inlet_specific_energy=inlet_depth #+0.5*v**2/g 
        z_in = 0.0
        z_out = -culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth  #+ 0.5*v**2/g
        E_out = z_out+outlet_depth  #+ 0.5*v**2/g
        delta_total_energy = E_in-E_out

        Q, v, d = boyd_generalised_culvert_model(inlet_depth, 
                                                 outlet_depth,
                                                 inlet_specific_energy, 
                                                 delta_total_energy, 
                                                 g,
                                                 culvert_length,
                                                 culvert_width,
                                                 culvert_height,
                                                 culvert_type,
                                                 manning,
                                                 sum_loss)
        
        print(Q, v, d)
        assert num.allclose(Q, 1.00, rtol=1.0e-2) #inflow
        assert num.allclose(v, 2.204, rtol=1.0e-2) #outflow velocity
        assert num.allclose(d, 0.76, rtol=1.0e-2) #depth at outlet used to calc v
    
               
#-------------------------------------------------------------
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_culvert_routines, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)


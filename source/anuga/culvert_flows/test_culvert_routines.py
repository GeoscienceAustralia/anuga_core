#!/usr/bin/env python


import unittest
import os.path
import sys

from anuga.utilities.system_tools import get_pathname_from_package
from anuga.culvert_flows.culvert_routines import boyd_generalised_culvert_model
import Numeric as num


class Test_culvert_routines(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass


    def NOtest_boyd_1(self):
        """test_boyd_1
        
        This tests the Boyd routine with data obtained from ??? by Petar Milevski    
        """
        # FIXME(Ole): This test fails (20 Feb 2009)

        g=9.81
        culvert_slope=0.1  # Downward

        inlet_depth=0.1
        outlet_depth=0.09

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
                                                 inlet_specific_energy, 
                                                 delta_total_energy, 
                                                 g,
                                                 culvert_length,
                                                 culvert_width,
                                                 culvert_height,
                                                 culvert_type,
                                                 manning,
                                                 sum_loss)
        
        print Q, v, d
        assert num.allclose(Q, 0.1)
        assert num.allclose(v, 0.93)
        assert num.allclose(d, 0.09)
        

    
               
#-------------------------------------------------------------
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_culvert_routines, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)


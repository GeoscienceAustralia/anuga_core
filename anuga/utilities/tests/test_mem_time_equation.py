#!/usr/bin/env python


import unittest
import numpy as num
import random


# Import standard shallow water domain and standard boundaries.
import anuga

from anuga.utilities.mem_time_equation import *

class Test_mem_time_equation(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_estimate_time_mem(self):
        points, vertices, boundary = anuga.rectangular_cross(10, 5,
                                               len1=10.0, len2=5.0) # Mesh

        domain = anuga.Domain(points, vertices, boundary)  # Create domain
        tri_num =  len(domain)
        yieldstep = 360
        finaltime = 3600
        
        time, memory = estimate_time_mem(domain, yieldstep, finaltime, 
                          use_test_constants=True, log_results=False)

        
        actual = system_constants[TEST_CON]['tri_a_T'] * tri_num ** 2 + \
             system_constants[TEST_CON]['tri_b_T'] * tri_num + \
              system_constants[TEST_CON]['tim_a_T'] * finaltime + \
              system_constants[TEST_CON]['fil_a_T'] * finaltime/yieldstep + \
               system_constants[TEST_CON]['cons_T']
        self.assertEqual(time, actual)

################################################################################

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_mem_time_equation, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)


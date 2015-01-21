#!/usr/bin/env python

import unittest
from math import sqrt

from anuga.abstract_2d_finite_volumes.generic_domain import Generic_Domain
from anuga.config import epsilon

import numpy as num


class Test_Domain(unittest.TestCase):
    def setUp(self):
        pass


    def tearDown(self):
        pass


    def test_simple(self):
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe, daf, dae
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        conserved_quantities = ['stage', 'xmomentum', 'ymomentum']
        other_quantities = ['elevation', 'friction']

        domain = Generic_Domain(points, vertices, None,
                        conserved_quantities, None, other_quantities)
        domain.check_integrity()

        for name in conserved_quantities + other_quantities:
            assert domain.quantities.has_key(name)


        assert num.alltrue(domain.get_conserved_quantities(0, edge=1) == 0.)


#-------------------------------------------------------------

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_Domain,'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)

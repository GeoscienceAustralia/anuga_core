#!/usr/bin/env python

import unittest
from math import sqrt, pi
import numpy

from anuga_1d.base.generic_domain import *



class Test_Generic_Domain(unittest.TestCase):
    def setUp(self):
        self.points = [0.0, 1.0, 2.0, 3.0]
        self.vertex_values = [[1.0,2.0],[4.0,5.0],[-1.0,2.0]]
        self.points2 = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
        
    def tearDown(self):
        pass
        #print "  Tearing down"


    def test_creation(self):
        domain = Generic_domain(self.points)

        assert numpy.allclose(domain.neighbours, \
            [[-1 , 1],[ 0 , 2],[ 1 ,-1]])
 
        assert numpy.allclose(domain.neighbour_vertices, \
            [[-1 , 0],[ 1 , 0],[ 1 ,-1]])

        assert numpy.allclose(domain.number_of_boundaries, \
            [1, 0, 1])

        assert numpy.allclose(domain.surrogate_neighbours, \
            [[0, 1],[0, 2],[1, 2]])

        assert numpy.allclose(domain.vertices, \
            [[ 0.,  1.],[ 1.,  2.],[ 2.,  3.]])


        assert numpy.allclose(domain.centroids, [0.5, 1.5, 2.5])

        assert numpy.allclose(domain.areas, [ 1.,  1.,  1.])

        assert numpy.allclose(domain.normals, \
            [[-1.,  1.],[-1.,  1.],[-1.,  1.]])

    def test_set_timestepping_method(self):

        domain = Generic_domain(self.points)

        domain.set_timestepping_method('euler')
        assert domain.timestepping_method == 'euler'

        domain.set_timestepping_method('rk2')
        assert domain.timestepping_method == 'rk2'

        domain.set_timestepping_method('rk3')
        assert domain.timestepping_method == 'rk3'

        domain.set_timestepping_method(1)
        assert domain.timestepping_method == 'euler'

        domain.set_timestepping_method(2)
        assert domain.timestepping_method == 'rk2'

        domain.set_timestepping_method(3)
        assert domain.timestepping_method == 'rk3'

        try:
            domain.set_timestepping_method(4)
        except:
            pass
        else:
            raise Exception,  'Should have raised "wrong method" error'

    def test_set_spatial_order(self):

        domain = Generic_domain(self.points)

        domain.set_spatial_order(1)
        assert domain.order == 1

        domain.set_spatial_order(2)
        assert  domain.order == 2


        try:
            domain.set_spatial_order(3)
        except:
            pass
        else:
            raise Exception,  'Should have raised "wrong order" error'


#-------------------------------------------------------------
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_Generic_Domain, 'test')
    #suite = unittest.makeSuite(Test_Shallow_Water, 'test_evolve_first_order')


    runner = unittest.TextTestRunner()
    runner.run(suite)

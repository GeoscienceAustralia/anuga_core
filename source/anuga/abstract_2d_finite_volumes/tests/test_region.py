#!/usr/bin/env python

import unittest
from math import sqrt

from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular
from anuga.shallow_water.shallow_water_domain import Domain
from anuga.abstract_2d_finite_volumes.generic_domain import Generic_Domain
from anuga.abstract_2d_finite_volumes.region import *
#from anuga.config import epsilon

import numpy as num


"""
This is what the mesh in these tests look like;

3---7
|5 /|
| /6|
2---6
|3 /|
| /2|
1---5
|1 /|
| /0|
0---4
"""

def add_x_y(x, y):
    return x+y

def give_me_23(x, y):
    return 23.0

class Test_region(unittest.TestCase):
    def setUp(self):
        pass


    def tearDown(self):
        pass


    def test_region_indices(self):
        """create region based on triangle lists."""

        #Create basic mesh
        points, vertices, boundary = rectangular(1, 3)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        
        region = Region(domain, indices=[0,2,3])

        expected_indices = [0,2,3]  
        assert num.allclose(region.indices, expected_indices)
        

    def test_region_polygon(self):
        """create region based on triangle lists."""

        #Create basic mesh
        points, vertices, boundary = rectangular(1, 3)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        
        poly = [[0.0,0.0], [0.5,0.0], [0.5,0.5]]
        
        #print poly
        region = Region(domain, polygon=poly)
        
        expected_indices = [1]
        assert num.allclose(region.indices, expected_indices)


    def test_region_polygon_expanded(self):
        """create region based on triangle lists."""

        #Create basic mesh
        points, vertices, boundary = rectangular(1, 3)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        
        region = Region(domain, polygon=[[0.0,0.0], [0.5,0.0], [0.5,0.5]], expand_polygon=True)
        
        expected_indices = [0,1,2,3]
        assert num.allclose(region.indices, expected_indices)

        
                  
#-------------------------------------------------------------

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_region, 'test')    
    runner = unittest.TextTestRunner()
    runner.run(suite)

#!/usr/bin/env python

#FIXME (Ole): Maxe this test independent of anything that inherits from General_mesh (namely shallow_water)


import unittest
from math import sqrt, pi


from quantity import *
from config import epsilon
from Numeric import allclose, array, ones, Float


class Test_General_Mesh(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_get_vertex_values(self):
        """
        get connectivity based on triangle lists.
        """
        from mesh_factory import rectangular
        from shallow_water import Domain
        from Numeric import zeros, Float

        #Create basic mesh
        points, vertices, boundary = rectangular(1, 3)
        domain = Domain(points, vertices, boundary)

        value = [7]
        #indexes = [1]  #FIXME (Ole): Should this be used
        assert  domain.get_vertices() == domain.triangles
        assert domain.get_vertices([0,4]) == [domain.triangles[0],
                                              domain.triangles[4]]
    def test_areas(self):
        from mesh_factory import rectangular
        from shallow_water import Domain
        from Numeric import zeros, Float

        #Create basic mesh
        points, vertices, boundary = rectangular(1, 3)
        domain = Domain(points, vertices, boundary)

	assert domain.get_area() == 1.0


    def test_get_unique_vertex_values(self):
        """
        get unique_vertex based on triangle lists.
        """
        from mesh_factory import rectangular
        from shallow_water import Domain
        from Numeric import zeros, Float

        #Create basic mesh
        points, vertices, boundary = rectangular(1, 3)
        domain = Domain(points, vertices, boundary)

        assert  domain.get_unique_vertices() == [0,1,2,3,4,5,6,7]
        unique_vertices = domain.get_unique_vertices([0,1,4])
        unique_vertices.sort()
        assert unique_vertices == [0,1,2,4,5,6,7]

        unique_vertices = domain.get_unique_vertices([0,4])
        unique_vertices.sort()
        assert unique_vertices == [0,2,4,5,6,7]

#-------------------------------------------------------------
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_General_Mesh,'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)


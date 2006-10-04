#!/usr/bin/env python


import unittest
from math import sqrt, pi


from anuga.config import epsilon
from Numeric import allclose, array, ones, Float
from general_mesh import General_mesh



class Test_General_Mesh(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass


    def test_get_vertex_coordinates(self):
        from mesh_factory import rectangular
        from Numeric import zeros, Float

        #Create basic mesh
        points, vertices, boundary = rectangular(1, 3)
        domain = General_mesh(points, vertices, boundary)

        assert allclose(domain.get_vertex_coordinates(unique=True), domain.coordinates)

        #assert allclose(domain.get_vertex_coordinates(), ...TODO
        #assert allclose(domain.get_vertex_coordinates(absolute=True), ...TODO
        
        

    def test_get_vertex_values(self):
        """Get connectivity based on triangle lists.
        """
        from mesh_factory import rectangular
        from Numeric import zeros, Float

        #Create basic mesh
        points, vertices, boundary = rectangular(1, 3)
        domain = General_mesh(points, vertices, boundary)

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
        domain = General_mesh(points, vertices, boundary)        

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
        domain = General_mesh(points, vertices, boundary)                

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


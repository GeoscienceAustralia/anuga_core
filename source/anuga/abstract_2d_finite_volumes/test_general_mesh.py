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
        nodes, triangles, _ = rectangular(1, 3)
        domain = General_mesh(nodes, triangles)


        assert allclose(domain.get_nodes(), nodes)


        M = domain.number_of_triangles        

        V = domain.get_vertex_coordinates()
        assert V.shape[0] == 3*M

        for i in range(M):
            for j in range(3):
                k = triangles[i,j]  #Index of vertex j in triangle i
                assert allclose(V[3*i+j,:], nodes[k])


    def test_get_vertex_coordinates_triangle_id(self):
        """test_get_vertex_coordinates_triangle_id
        Test that vertices for one triangle can be returned.
        """
        from mesh_factory import rectangular
        from Numeric import zeros, Float

        #Create basic mesh
        nodes, triangles, _ = rectangular(1, 3)
        domain = General_mesh(nodes, triangles)


        assert allclose(domain.get_nodes(), nodes)


        M = domain.number_of_triangles        

        for i in range(M):
            V = domain.get_vertex_coordinates(triangle_id=i)
            assert V.shape[0] == 3

            for j in range(3):
                k = triangles[i,j]  #Index of vertex j in triangle i
                assert allclose(V[j,:], nodes[k])


        
        

    def test_get_vertex_values(self):
        """Get connectivity based on triangle lists.
        """
        from mesh_factory import rectangular
        from Numeric import zeros, Float

        #Create basic mesh
        nodes, triangles, _ = rectangular(1, 3)
        domain = General_mesh(nodes, triangles)

        value = [7]
        assert allclose(domain.get_triangles(), triangles)
        assert allclose(domain.get_triangles([0,4]),
                        [triangles[0], triangles[4]])
        

    def test_vertex_value_indices(self):
        """Check that structures are correct.
        """
        from mesh_factory import rectangular
        from Numeric import zeros, Float, array

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        nodes = array([a, b, c, d, e, f])
        #bac, bce, ecf, dbe, daf, dae
        triangles = array([[1,0,2], [1,2,4], [4,2,5], [3,1,4]])

        domain1 = General_mesh(nodes, triangles)
        
        #Create larger mesh
        nodes, triangles, _ = rectangular(3, 6)
        domain2 = General_mesh(nodes, triangles)

        # Test both meshes
        for domain in [domain1, domain2]:
            assert sum(domain.number_of_triangles_per_node) ==\
                   len(domain.vertex_value_indices)

            # Check number of triangles per node
            count = [0]*domain.number_of_nodes
            for triangle in domain.triangles:
                for i in triangle:
                    count[i] += 1

            #print count
            #
            assert allclose(count, domain.number_of_triangles_per_node)
            
            # Check indices
            current_node = 0
            k = 0 # Track triangles touching on node
            for index in domain.vertex_value_indices:
                k += 1
                
                triangle = index / 3
                vertex = index % 3

                assert domain.triangles[triangle, vertex] == current_node

                if domain.number_of_triangles_per_node[current_node] == k:
                    # Move on to next node
                    k = 0
                    current_node += 1
                

    def test_get_triangles_and_vertices_per_node(self):
        """test_get_triangles_and_vertices_per_node -

        Test that tuples of triangle, vertex can be extracted
        from inverted triangles structure

        """

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        nodes = array([a, b, c, d, e, f])
        #bac, bce, ecf, dbe, daf, dae
        triangles = array([[1,0,2], [1,2,4], [4,2,5], [3,1,4]])

        domain = General_mesh(nodes, triangles)

        # One node
        L = domain.get_triangles_and_vertices_per_node(node=2)

        assert allclose(L[0], [0, 2])
        assert allclose(L[1], [1, 1])
        assert allclose(L[2], [2, 1])

        # All nodes
        ALL = domain.get_triangles_and_vertices_per_node()

        assert len(ALL) == 6
        for i, Lref in enumerate(ALL):
            L = domain.get_triangles_and_vertices_per_node(node=i)
            assert allclose(L, Lref)
            

        

        


    def test_areas(self):
        from mesh_factory import rectangular
        from shallow_water import Domain
        from Numeric import zeros, Float

        #Create basic mesh
        points, vertices, boundary = rectangular(1, 3)
        domain = General_mesh(points, vertices)        

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
        domain = General_mesh(points, vertices)                

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
    #suite = unittest.makeSuite(Test_General_Mesh,'test_vertex_value_indices')
    runner = unittest.TextTestRunner()
    runner.run(suite)


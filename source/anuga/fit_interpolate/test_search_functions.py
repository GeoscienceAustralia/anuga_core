#!/usr/bin/env python


import unittest
from search_functions import search_tree_of_vertices, set_last_triangle
from search_functions import _search_triangles_of_vertices

from anuga.abstract_2d_finite_volumes.neighbour_mesh import Mesh
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular
from anuga.utilities.polygon import is_inside_polygon
from anuga.utilities.quad import build_quadtree, Cell


class Test_search_functions(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def NOtest_that_C_extension_compiles(self):
        FN = 'search_functions_ext.c'
        try:
            import search_functions_ext
        except:
            from compile import compile

            try:
                compile(FN)
            except:
                raise 'Could not compile %s' %FN
            else:
                import search_functions_ext



    def test_small(self):
        """test_small: Two triangles
        """

        points, vertices, boundary = rectangular(1, 1, 1, 1)
        mesh = Mesh(points, vertices, boundary)

        #Test that points are arranged in a counter clock wise order
        mesh.check_integrity()

        root = build_quadtree(mesh, max_points_per_cell = 1)
        set_last_triangle()

        x = [0.2, 0.7]
        found, s0, s1, s2, k = search_tree_of_vertices(root, mesh, x)
        assert k == 1 # Triangle one
        assert found is True

    def test_bigger(self):
        """test_larger mesh
        """

        points, vertices, boundary = rectangular(4, 4, 1, 1)
        mesh = Mesh(points, vertices, boundary)

        #Test that points are arranged in a counter clock wise order
        mesh.check_integrity()

        root = build_quadtree(mesh, max_points_per_cell = 4)
        set_last_triangle()

        for x in [[0.6, 0.3], [0.1, 0.2], [0.7,0.7],
                  [0.1,0.9], [0.4,0.6], [0.9,0.1],
                  [10, 3]]:
                
            found, s0, s1, s2, k = search_tree_of_vertices(root, mesh, x)

            if k >= 0:
                V = mesh.get_vertex_coordinates(k) # nodes for triangle k
                assert is_inside_polygon(x, V)
                assert found is True
                #print k, x
            else:
                assert found is False                

        

    def test_large(self):
        """test_larger mesh and different quad trees
        """

        points, vertices, boundary = rectangular(10, 12, 1, 1)
        mesh = Mesh(points, vertices, boundary)

        #Test that points are arranged in a counter clock wise order
        mesh.check_integrity()

        
        for m in range(8):
            root = build_quadtree(mesh, max_points_per_cell = m)
            set_last_triangle()
            #print m, root.show()

            for x in [[0.6, 0.3], [0.1, 0.2], [0.7,0.7],
                      [0.1,0.9], [0.4,0.6], [0.9,0.1],
                      [10, 3]]:
                
                found, s0, s1, s2, k = search_tree_of_vertices(root, mesh, x)

                if k >= 0:
                    V = mesh.get_vertex_coordinates(k) # nodes for triangle k
                    assert is_inside_polygon(x, V)
                    assert found is True
                else:
                    assert found is False                

                

    def test_underlying_function(self):
        """test_larger mesh and different quad trees
        """

        points, vertices, boundary = rectangular(2, 2, 1, 1)
        mesh = Mesh(points, vertices, boundary)

        root = build_quadtree(mesh, max_points_per_cell = 4)
        set_last_triangle()

        # One point
        x = [0.5, 0.5]
        candidate_vertices = root.search(x[0], x[1])

        #print x, candidate_vertices
        found, sigma0, sigma1, sigma2, k = \
               _search_triangles_of_vertices(mesh,
                                             candidate_vertices,
                                             x)

        if k >= 0:
            V = mesh.get_vertex_coordinates(k) # nodes for triangle k
            assert is_inside_polygon(x, V)
            assert found is True
        else:
            assert found is False                

        

        # More points    
        for x in [[0.6, 0.3], [0.1, 0.2], [0.7,0.7],
                  [0.1,0.9], [0.4,0.6], [0.9,0.1],
                  [10, 3]]:
                
            candidate_vertices = root.search(x[0], x[1])

            #print x, candidate_vertices
            found, sigma0, sigma1, sigma2, k = \
                   _search_triangles_of_vertices(mesh,
                                                 candidate_vertices,
                                                 x)
            if k >= 0:
                V = mesh.get_vertex_coordinates(k) # nodes for triangle k
                assert is_inside_polygon(x, V)
                assert found is True
            else:
                assert found is False

                

    def expanding_search(self):
        """test_larger mesh and different quad trees
        """
        
        p0 = [2,1]
        p1 = [4,1]
        p2 = [4.,4]
        p3 = [2,4]
        p4 = [5,4]

        p5 = [-1,-1]
        p6 = [1,-1]
        p7 = [1,1]
        p8 = [-1,1]

        points = [p0,p1,p2, p3,p4,p5,p6,p7,p8]
        #
        vertices = [[0,1,2],[0,2,3],[1,4,2],[5,6,7], [5,7,8]]
        mesh = Mesh(points, vertices)

        # Don't do this, want to control the max and mins
        #root = build_quadtree(mesh, max_points_per_cell=4)
    

        root = Cell(-3, 9, -3, 9, mesh,
                    max_points_per_cell = 4)
        #Insert indices of all vertices
        root.insert( range(mesh.number_of_nodes) )

        #Build quad tree and return
        root.split()
        
        # One point
        #x = [3.5, 1.5]
        x = [2.5, 1.5]
        element_found, sigma0, sigma1, sigma2, k = \
                       search_tree_of_vertices(root, mesh, x)
        # One point
        x = [3.00005, 2.999994]
        element_found, sigma0, sigma1, sigma2, k = \
                       search_tree_of_vertices(root, mesh, x)
        assert element_found is True
        assert k == 1
        

#-------------------------------------------------------------
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_search_functions,'test')
    #suite = unittest.makeSuite(Test_search_functions,'expanding_search')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
    

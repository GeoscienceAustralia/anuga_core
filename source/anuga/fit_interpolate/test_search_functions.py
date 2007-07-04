#!/usr/bin/env python


import unittest
from search_functions import *

from Numeric import zeros, array, allclose

from anuga.abstract_2d_finite_volumes.neighbour_mesh import Mesh
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

from anuga.utilities.quad import build_quadtree





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

        #print mesh.nodes
        #print mesh.triangles
        root = build_quadtree(mesh, max_points_per_cell = 1)
        #print 'root', root.show()

        x = [0.7, 0.7]
        found, s0, s1, s2, k = search_tree_of_vertices(root, mesh, x)
        #print k
        # What is k??
        assert found is True

        
        

        
#-------------------------------------------------------------
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_search_functions,'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
    

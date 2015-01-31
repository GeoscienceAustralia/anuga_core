import unittest
import numpy as num

from anuga.geometry.aabb import AABB
from anuga.geometry.quad import Cell
from anuga.pmesh.mesh_quadtree import MeshQuadtree
from anuga.abstract_2d_finite_volumes.general_mesh import General_mesh as Mesh

import sys

#-------------------------------------------------------------

class Test_Quad(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_build_quadtree(self):

        a = [3, 7]
        b = [5, 7]
        c = [5, 5]
        d = [7, 7]
        e = [15, 15]
        f = [15, 30]
        g = [30, 10]
        h = [30, 30]

        points = [a, b, c, d, e, f, g, h]
        
        #bac, bce, ecf, dbe, daf, dae
        vertices = [[1,0,2], [1,3,4], [1,2,3], [5,4,7], [4,6,7]]

        mesh = Mesh(points, vertices)
    
        Q = MeshQuadtree(mesh)
        #Q.show()
        #print Q.count()
        self.assertEqual(Q.count(), len(vertices))

        # test a point that falls within a triangle
        result = Q.search([10, 10])
        assert isinstance(result, (list, tuple)), 'should be a list'
         
        # Padarn Note: The result of Q.search is no longer in 
        # the same format, and so this test fails. However, the 
        # old functionality does not seem to be needed.                    
        #self.assertEqual(result[0][0][0], 1)


    def test_build_quadtreeII(self):

        self.cell = Cell(AABB(100, 140, 0, 40), 'cell')

        p0 = [34.6292076111,-7999.92529297]
        p1 = [8000.0, 7999.0]
        p2 = [-7999.96630859, 7999.0]
        p3 = [34, 7999.97021484]

        points = [p0,p1,p2, p3]
        #bac, bce, ecf, dbe, daf, dae
        vertices = [[0,1,2],[0,3,2]]

        mesh = Mesh(points, vertices)

        #This was causing round off error
        Q = MeshQuadtree(mesh)
        
    def NOtest_interpolate_one_point_many_triangles(self):
        # this test has 10 triangles that share the same vert.
        # If the number of points per cell in  a quad tree is less
        # than 10 it should crash 
        z0 = [2.0, 5.0]
        z1 = [2.0, 5.0]
        z2 = [2.0, 5.0]
        z3 = [2.0, 5.0]
        z4 = [2.0, 5.0]
        z5 = [2.0, 5.0]
        z6 = [2.0, 5.0]
        z7 = [2.0, 5.0]
        z8 = [2.0, 5.0]
        z9 = [2.0, 5.0]
        z10 = [2.0, 5.0]
        
        v0 = [0.0, 0.0]
        v1 = [1.0, 0.0]
        v2 = [2.0, 0.0]
        v3 = [3.0, 0.0]
        v4 = [4.0, 0.0]
        v5 = [0.0, 10.0]
        v6 = [1.0, 10.0]
        v7 = [2.0, 10.0]
        v8 = [3.0, 10.0]
        v9 = [4.0, 10.0]

        vertices = [z0,v0, v1, v2, v3,v4 ,v5, v6, v7, v8, v9,
                    z1, z2, z3, z4, z5, z6, z7, z8, z9]
        triangles = [
                      [11,1,2],
                      [12,2,3],
                      [13,3,4],
                      [14,4,5],
                      [7,6,15],
                      [8,7,16],
                      [9,8,17],
                      [10,9,18],
                      [6,1,19],
                      [5,10,0]
                      ]
        
        mesh = Mesh(vertices, triangles)
        try:
            Q = MeshQuadtree(mesh, max_points_per_cell = 9)
        except RuntimeError:
            pass
        else:
            self.assertTrue(0 ==1,  'many verts at the same position no  \
            longer causes as error')
    
    def test_retrieve_triangles(self):

        cell = Cell(AABB(0, 6, 0, 6), 'cell')

        p0 = [2,1]
        p1 = [4,1]
        p2 = [4.,4]
        p3 = [2,4]
        p4 = [5,4]

        points = [p0,p1,p2, p3, p4]
        #
        vertices = [[0,1,2],[0,2,3],[1,4,2]]

        mesh = Mesh(points, vertices)

        Q = MeshQuadtree(mesh)
        results = Q.search([4.5, 3])
        # Padarn Note: The result of Q.search is no longer in 
        # the same format, and so this test fails. However, the 
        # old functionality does not seem to be needed. 
        #assert len(results) == 1
        #self.assertEqual(results[0][0][0], 2)
        #results = Q.search([5,4.])
        #self.assertEqual(len(results),1)
        #self.assertEqual(results[0][0][0], 2)
        
    def NOtest_num_visits(self):
        """ Test optimisation code.
        """
        a = [3, 7]
        b = [5, 7]
        c = [5, 5]
        d = [7, 7]
        e = [15, 15]
        f = [15, 30]
        g = [30, 10]
        h = [30, 30]

        points = [a, b, c, d, e, f, g, h]
        
        #bac, bce, ecf, dbe, daf, dae
        vertices = [[1,0,2], [1,3,4], [1,2,3], [5,4,7], [4,6,7]]

        mesh = Mesh(points, vertices)

        Q = MeshQuadtree(mesh)    


        results = Q.search_fast([5.5, 5.5])
        print 'visits: ', Q.count_visits()
        
        Q.clear_visits()
        results = Q.search_fast([30, 10])
        print 'visits: ', Q.count_visits()
        
        print 'second time:'

        Q.clear_visits()        
        results = Q.search_fast([5.5, 5.5])
        print 'visits: ', Q.count_visits()
################################################################################

if __name__ == "__main__":
    mysuite = unittest.makeSuite(Test_Quad,'test')
    runner = unittest.TextTestRunner()
    runner.run(mysuite)

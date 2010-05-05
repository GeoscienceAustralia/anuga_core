import unittest
import numpy as num

from quad import AABB, Cell, build_quadtree
from anuga.abstract_2d_finite_volumes.general_mesh import General_mesh as Mesh

import types, sys

#-------------------------------------------------------------

class Test_Quad(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_AABB_contains(self):
        box = AABB(1, 21, 1, 11)
        assert box.contains(10, 5)
        assert box.contains(1, 1)
        assert box.contains(20, 6)
        assert not box.contains(-1, -1)
        assert not box.contains(5, 70)
        assert not box.contains(6, -70)
        assert not box.contains(-1, 6)
        assert not box.contains(50, 6)        
        
    def test_AABB_split_vert(self):
        parent = AABB(1, 21, 1, 11)
        
        child1, child2 = parent.split(0.6)

        self.assertEqual(child1.xmin, 1)
        self.assertEqual(child1.xmax, 13)
        self.assertEqual(child1.ymin, 1)
        self.assertEqual(child1.ymax, 11)
        
        self.assertEqual(child2.xmin, 9)
        self.assertEqual(child2.xmax, 21)
        self.assertEqual(child2.ymin, 1)
        self.assertEqual(child2.ymax, 11)    

    def test_AABB_split_horiz(self):
        parent = AABB(1, 11, 1, 41)
        
        child1, child2 = parent.split(0.6)

        self.assertEqual(child1.xmin, 1)
        self.assertEqual(child1.xmax, 11)
        self.assertEqual(child1.ymin, 1)
        self.assertEqual(child1.ymax, 25)
        
        self.assertEqual(child2.xmin, 1)
        self.assertEqual(child2.xmax, 11)
        self.assertEqual(child2.ymin, 17)
        self.assertEqual(child2.ymax, 41)          
        
    def test_add_data(self):
        cell = Cell(AABB(0,10, 0,5))
        cell.insert([(AABB(1,3, 1, 3), 111), (AABB(8,9, 1, 2), 222),  \
                     (AABB(7, 8, 3, 4), 333), (AABB(1, 10, 0, 1), 444)])

        result = cell.retrieve()
        assert type(result) in [types.ListType,types.TupleType],\
                            'should be a list'

        self.assertEqual(len(result),4)
        
    def test_search(self):
        test_region = (AABB(8,9, 1, 2), 222)
        cell = Cell(AABB(0,10, 0,5))
        cell.insert([(AABB(1,3, 1, 3), 111), test_region,  \
                     (AABB(7, 8, 3, 4), 333), (AABB(1, 10, 0, 1), 444)])

        result =  cell.search(x = 8.5, y = 1.5, get_vertices=True)
        assert type(result) in [types.ListType,types.TupleType],\
                            'should be a list'
        self.assertEqual(result, [test_region], 'only 1 point should intersect')


    def test_clear_1(self):
        cell = Cell(AABB(0,10, 0,5))    
        cell.insert([(AABB(1,3, 1, 3), 111), (AABB(8,9, 1, 2), 222),  \
                     (AABB(7, 8, 3, 4), 333), (AABB(1, 10, 0, 1), 444)])
                     
        assert len(cell.retrieve()) == 4
        cell.clear()

        assert len(cell.retrieve()) == 0

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
    
        Q = build_quadtree(mesh)
        #Q.show()
        #print Q.count()
        self.assertEqual(Q.count(), len(vertices))

        # test a point that falls within a triangle
        result = Q.search(10, 10, get_vertices=True)
        assert type(result) in [types.ListType,types.TupleType],\
                            'should be a list'
        pos, index = result[0]
        self.assertEqual(index, 1)


    def test_build_quadtreeII(self):

        self.cell = Cell(AABB(100, 140, 0, 40), 'cell')

        p0 = [34.6292076111,-7999.92529297]
        p1 = [8000.0, 7999.0]
        p2 = [-7999.96630859, 7999.0]
        p3 = [34, 7999.97021484]

        points = [p0,p1,p2, p3]
        #bac, bce, ecf, dbe, daf, dae
        vertices = [[0,1,2],[0,2,3]]

        mesh = Mesh(points, vertices)

        #This was causing round off error
        Q = build_quadtree(mesh)
        
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
            Q = build_quadtree(mesh, max_points_per_cell = 9)
        except RuntimeError:
            pass
        else:
            self.failUnless(0 ==1,  'many verts at the same position no  \
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

        Q = build_quadtree(mesh)
        results = Q.search(4.5, 3)
        assert len(results) == 1
        self.assertEqual(results[0], 2)
        results = Q.search(5,4.)
        self.assertEqual(len(results),1)
        self.assertEqual(results[0], 2)
################################################################################

if __name__ == "__main__":
    mysuite = unittest.makeSuite(Test_Quad,'test')
    runner = unittest.TextTestRunner()
    runner.run(mysuite)

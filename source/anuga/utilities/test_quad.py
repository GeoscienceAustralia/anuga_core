import unittest
##import numpy

from quad import Cell, build_quadtree
from anuga.abstract_2d_finite_volumes.general_mesh import General_mesh as Mesh

import types, sys

#-------------------------------------------------------------

class Test_Quad(unittest.TestCase):

    def setUp(self):

        a = [3, 107]
        b = [5, 107]
        c = [5, 105]
        d = [7, 107]
        e = [15, 115]
        f = [15, 130]
        g = [30, 110]
        h = [30, 130]

        points = [a, b, c, d, e, f, g, h]
        
        #bac, bce, ecf, dbe, daf, dae
        vertices = [[1,0,2], [1,3,4], [1,2,3], [5,4,7], [4,6,7]]

	mesh = Mesh(points, vertices)
        self.mesh = mesh
        self.cell = Cell(100, 140, 0, 40, mesh, 'cell')

    def tearDown(self):
        pass

    def test_add_points_2_cell(self):
        self.cell.insert(0)
        self.cell.insert(1)

        result = self.cell.retrieve()
        assert type(result) in [types.ListType,types.TupleType],\
                                'should be a list'
        self.assertEqual(len(result),2)

    def test_add_points_2_cellII(self):
        self.cell.insert([0,1,2,3,4,5,6,7])

        result = self.cell.retrieve()
        assert type(result) in [types.ListType,types.TupleType],\
	                        'should be a list'
        self.assertEqual(len(result),8)


    def test_search(self):
        self.cell.insert([0,1,2,3,4,5,6,7])
	self.cell.split(4)

        result =  self.cell.search(x = 1, y = 101, get_vertices=True)
        assert type(result) in [types.ListType,types.TupleType],\
	                        'should be a list'
        self.assertEqual(result, [0,1,2,3])


    def test_clear_1(self):
        self.cell.insert([0,1,2,3,4,5,6,7])
	assert self.cell.count() == 8
	self.cell.clear()

        #This one actually revealed a bug :-)
	assert self.cell.count() == 0

    def test_clear_2(self):
        self.cell.insert([0,1,2,3,4,5,6,7])
	assert self.cell.count() == 8
	self.cell.split(2)
	assert self.cell.count() == 8

	self.cell.clear()
	assert self.cell.count() == 0



    def test_split(self):
        self.cell.insert([0,1,2,3,4,5,6,7], split = False)

	#No children yet
	assert self.cell.children is None
	assert self.cell.count() == 8

        #Split
	self.cell.split(4)
	#self.cell.show()
	#self.cell.show_all()


	#Now there are children
	assert self.cell.children is not None
	assert self.cell.count() == 8



    def test_collapse(self):
        self.cell.insert([0,1,2,3,4,5,6,7], split = False)

        #Split maximally
	self.cell.split(1)

	#Now there are children
	assert self.cell.children is not None
	assert self.cell.count() == 8

	#Collapse
	self.cell.collapse(8)

	#No children
	assert self.cell.children is None
	assert self.cell.count() == 8

    def test_build_quadtree(self):

        Q = build_quadtree(self.mesh)
        #Q.show()
        #print Q.count()
	assert Q.count() == 8



        result = Q.search(3, 105, get_vertices=True)
        assert type(result) in [types.ListType,types.TupleType],\
	                        'should be a list'
        #print "result",result
        self.assertEqual(result, [0,1,2,3])


    def test_build_quadtreeII(self):

        self.cell = Cell(100, 140, 0, 40, 'cell')

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
        
    def test_interpolate_one_point_many_triangles(self):
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

        cell = Cell(0, 6, 0, 6, 'cell', max_points_per_cell=4)

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
        results = Q.search(5,1)
        assert len(results),2
        #print "results", results
        #print "results[0][0]", results[0][0]
        assert results[0],0
        assert results[1],2
        assert results[0][1],[[ 2.,  1.],
                     [ 4.,  1.],
                     [ 4.,  4.]]
        assert results[1][1],[[ 4.,  1.],
                     [ 5.,  4.],
                     [ 4.,  4.]]
        # this is the normals
        assert results[0][1][1],[[1.,  0.],
                     [-0.83205029,  0.5547002],
                     [ 0.,  -1.]]
                     
        # assert numpy.allclose(numpy.array(results),[[[ 2.,  1.],
        #[ 4.,  1.], [ 4.,  4.]], [[ 4.,  1.],[ 5.,  4.],[ 4.,  4.]]] )
        results = Q.search(5,4.)
        ### print "results",results 
        # results_dic={}
        # results_dic.update(results)
        assert len(results),3
        #print "results_dic[0]", results_dic[0]
        assert results[0][1],[[ 2.,  1.],
                     [ 4.,  1.],
                     [ 4.,  4.]]
        assert results[1][1],[[ 2.,  1.],
                     [ 4.,  4.],
                     [ 2.,  4.]]
        assert results[2][1],[[ 4.,  1.],
                     [ 5.,  4.],
                     [ 4.,  4.]]
        #assert numpy.allclose(numpy.array(results),[[[ 2.,  1.],[ 4.,  1.], [ 4.,  4.]]
         #                               ,[[ 2.,  1.],[ 4.,  4.], [ 2.,  4.]],
        #[[ 4.,  1.],  [ 5.,  4.], [ 4.,  4.]],
         #                               [[ 4.,  1.], [ 5.,  4.], [ 4.,  4.]]])
        
        
#-------------------------------------------------------------
if __name__ == "__main__":

    mysuite = unittest.makeSuite(Test_Quad,'test')
    #mysuite = unittest.makeSuite(Test_Quad,'test_interpolate_one_point_many_triangles')
    runner = unittest.TextTestRunner()
    runner.run(mysuite)

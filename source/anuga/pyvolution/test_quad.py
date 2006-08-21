import unittest
from quad import Cell, build_quadtree

#from domain import *
from general_mesh import General_mesh as Mesh

import types, sys

#-------------------------------------------------------------

class Test_Quad(unittest.TestCase):

    def setUp(self):
        self.cell = Cell(100, 140, 0, 40, 'cell')

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
        Cell.initialise(mesh)

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

        result =  self.cell.search(x = 1, y = 101)
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
	assert Q.count() == 8

        #Q.show()

        result = Q.search(3, 105)
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

#-------------------------------------------------------------
if __name__ == "__main__":

    mysuite = unittest.makeSuite(Test_Quad,'test')
    runner = unittest.TextTestRunner()
    runner.run(mysuite)

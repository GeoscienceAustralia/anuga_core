""" Test for the geometry classes.

    Pylint quality rating as of June 2010: 8.51/10.
"""

import unittest

from anuga.geometry.aabb import AABB
from anuga.geometry.quad import Cell

#-------------------------------------------------------------

class Test_Geometry(unittest.TestCase):
    """ Test geometry classes. """
    def setUp(self):
        """ Generic set up for geometry tests. """
        pass

    def tearDown(self):
        """ Generic shut down for geometry tests. """
        pass

    def test_aabb_contains(self):
        """ Test if point is correctly classified as falling inside or
            outside of bounding box. """
        box = AABB(1, 21, 1, 11)
        assert box.contains([10, 5])
        assert box.contains([1, 1])
        assert box.contains([20, 6])
        assert not box.contains([-1, -1])
        assert not box.contains([5, 70])
        assert not box.contains([6, -70])
        assert not box.contains([-1, 6])
        assert not box.contains([50, 6])        
        
    def test_aabb_split_vert(self):
        """ Test that a bounding box can be split correctly along an axis.
        """
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

    def test_aabb_split_horiz(self):
        """ Test that a bounding box will be split along the horizontal axis
        correctly. """
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
        """ Test add and retrieve arbitrary data from tree structure. """
        cell = Cell(AABB(0, 10, 0, 5), None)
        cell.insert([(AABB(1, 3, 1, 3), 111), (AABB(8, 9, 1, 2), 222),  \
                     (AABB(7, 8, 3, 4), 333), (AABB(1, 10, 0, 1), 444)])

        result = cell.retrieve()
        assert isinstance(result, (list, tuple)), 'should be a list'

        self.assertEqual(len(result), 4)
        
    def test_search(self):
        """ Test search tree for an intersection. """
        test_tag = 222
        cell = Cell(AABB(0, 10, 0,5), None)
        cell.insert([(AABB(1, 3, 1, 3), 111), (AABB(8, 9, 1, 2), test_tag),  \
                     (AABB(7, 8, 3, 4), 333), (AABB(1, 10, 0, 1), 444)])

        result = cell.search([8.5, 1.5])
        assert isinstance(result, (list, tuple)), 'should be a list'
        assert(len(result) == 1)
        data, _ = result[0]
        self.assertEqual(data, test_tag, 'only 1 point should intersect')

    def test_get_siblings(self):
        """ Make sure children know their parent. """
        cell = Cell(AABB(0, 10, 0, 5), None)
        cell.insert([(AABB(1, 3, 1, 3), 111), (AABB(8, 9, 1, 2), 222)])

        assert len(cell.children) == 2
        assert cell.parent == None
        assert cell.children[0].parent == cell
        assert cell.children[1].parent == cell


################################################################################

if __name__ == "__main__":
    mysuite = unittest.makeSuite(Test_Geometry, 'test')
    runner = unittest.TextTestRunner()
    runner.run(mysuite)

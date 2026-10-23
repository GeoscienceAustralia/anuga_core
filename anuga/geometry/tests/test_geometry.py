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
        assert cell.parent is None
        assert cell.children[0].parent == cell
        assert cell.children[1].parent == cell


class Test_Geometry_extra(unittest.TestCase):
    """Additional coverage tests for quad.py and aabb.py."""

    # ---- AABB tests ----

    def test_aabb_repr(self):
        """AABB.__repr__ covers line 51."""
        a = AABB(1.0, 3.0, 2.0, 4.0)
        s = repr(a)
        self.assertIn('AABB', s)
        self.assertIn('xmin', s)

    def test_aabb_grow(self):
        """AABB.grow covers lines 60-63."""
        a = AABB(0.0, 10.0, 0.0, 5.0)
        a.grow(0.1)
        # After growing, extents should expand
        self.assertGreater(a.xmax, 10.0)
        self.assertLess(a.xmin, 0.0)

    def test_aabb_from_bad_single_arg_raises(self):
        """AABB from non-list single arg raises (lines 38-39)."""
        with self.assertRaises(Exception) as cm:
            AABB(42)  # single non-list arg
        self.assertIn('Single parameter', str(cm.exception))

    def test_aabb_include_expands_bounds(self):
        """AABB.include covers line 127 (new xmin smaller)."""
        a = AABB(2.0, 8.0, 2.0, 8.0)
        a.include([[-1.0, 5.0], [5.0, 12.0]])
        self.assertAlmostEqual(a.xmin, -1.0)
        self.assertAlmostEqual(a.ymax, 12.0)

    # ---- Cell/quad tests ----

    def test_cell_repr_with_children(self):
        """Cell.__repr__ with children covers lines 40-44.
        Note: Cell.__repr__ has a pre-existing bug (self.name not set in __init__).
        """
        cell = Cell(AABB(0, 10, 0, 5), None)
        cell.name = 'root'  # work around bug: name not stored in __init__
        cell.insert([(AABB(1, 3, 1, 3), 111), (AABB(8, 9, 1, 2), 222)])
        # Give children a name too
        for i, child in enumerate(cell.children):
            child.name = f'child{i}'
        s = repr(cell)
        self.assertIn('children', s)

    def test_cell_insert_single_item(self):
        """Cell.insert with single non-list item covers line 56."""
        cell = Cell(AABB(0, 10, 0, 5), None)
        cell.insert((AABB(1, 3, 1, 3), 42))  # single tuple, not a list
        self.assertEqual(len(cell.retrieve()), 1)

    def test_cell_count(self):
        """Cell.count covers lines 123-131."""
        cell = Cell(AABB(0, 10, 0, 5), None)
        cell.insert([(AABB(1, 3, 1, 3), 111), (AABB(8, 9, 1, 2), 222)])
        count = cell.count()
        self.assertEqual(count, 2)

    def _set_name_recursive(self, cell, prefix='cell'):
        """Recursively set .name on all cells to work around the bug in __init__."""
        cell.name = prefix
        if cell.children:
            for i, child in enumerate(cell.children):
                self._set_name_recursive(child, f'{prefix}.{i}')

    def test_cell_show(self):
        """Cell.show covers lines 136-143 (just ensure it runs)."""
        cell = Cell(AABB(0, 10, 0, 5), None)
        cell.insert([(AABB(1, 3, 1, 3), 111), (AABB(8, 9, 1, 2), 222)])
        self._set_name_recursive(cell)
        cell.show(depth=0)  # output goes to stdout; just ensure it runs

    def test_cell_get_siblings(self):
        """Cell.get_siblings covers lines 183-188."""
        cell = Cell(AABB(0, 10, 0, 5), None)
        cell.insert([(AABB(1, 3, 1, 3), 111), (AABB(8, 9, 1, 2), 222)])
        # cell has two children; each child should have one sibling
        siblings = cell.children[0].get_siblings()
        self.assertEqual(len(siblings), 1)
        self.assertIs(siblings[0], cell.children[1])

    def test_cell_get_siblings_no_parent(self):
        """Cell.get_siblings with no parent returns [] (line 184)."""
        cell = Cell(AABB(0, 10, 0, 5), None)
        self.assertEqual(cell.get_siblings(), [])


################################################################################

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_Geometry)
    runner = unittest.TextTestRunner()
    runner.run(suite)

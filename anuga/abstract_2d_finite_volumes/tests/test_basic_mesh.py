"""Tests for basic_mesh.py — Basic_mesh and factory functions."""
import unittest
import numpy as num

from anuga.abstract_2d_finite_volumes.basic_mesh import (
    Basic_mesh,
    rectangular_basic_mesh,
    rectangular_cross_basic_mesh,
)
from anuga.coordinate_transforms.geo_reference import Geo_reference


def _simple_mesh():
    """Return a minimal 2-triangle mesh for testing."""
    nodes = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]
    triangles = [[0, 1, 2], [1, 3, 2]]
    boundary = {(0, 0): 'left', (0, 1): 'bottom', (1, 0): 'right', (1, 2): 'top'}
    return nodes, triangles, boundary


class Test_Basic_mesh_init(unittest.TestCase):

    def test_basic_construction(self):
        nodes, triangles, boundary = _simple_mesh()
        bm = Basic_mesh(nodes, triangles, boundary=boundary)
        self.assertEqual(bm.number_of_triangles, 2)
        self.assertEqual(bm.number_of_nodes, 4)

    def test_default_geo_reference(self):
        nodes, triangles, _ = _simple_mesh()
        bm = Basic_mesh(nodes, triangles)
        self.assertIsInstance(bm.geo_reference, Geo_reference)

    def test_custom_geo_reference(self):
        nodes, triangles, _ = _simple_mesh()
        gr = Geo_reference(zone=55, xllcorner=100.0, yllcorner=200.0)
        bm = Basic_mesh(nodes, triangles, geo_reference=gr)
        self.assertIs(bm.geo_reference, gr)

    def test_default_boundary_is_empty_dict(self):
        nodes, triangles, _ = _simple_mesh()
        bm = Basic_mesh(nodes, triangles)
        self.assertEqual(bm.boundary, {})

    def test_nodes_stored_as_float_array(self):
        nodes = [[0, 0], [1, 0], [0, 1]]
        triangles = [[0, 1, 2]]
        bm = Basic_mesh(nodes, triangles)
        self.assertEqual(bm.nodes.dtype, float)

    def test_triangles_stored_as_int_array(self):
        nodes = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]
        triangles = [[0, 1, 2]]
        bm = Basic_mesh(nodes, triangles)
        self.assertTrue(num.issubdtype(bm.triangles.dtype, num.integer))

    def test_precomputed_neighbours(self):
        nodes, triangles, _ = _simple_mesh()
        # Pass explicit triangle_neighbours
        pre_nbrs = num.array([[-1, -1, 1], [-1, -1, 0]])
        bm = Basic_mesh(nodes, triangles, triangle_neighbours=pre_nbrs)
        num.testing.assert_array_equal(bm.neighbours, pre_nbrs)


class Test_Basic_mesh_protocol(unittest.TestCase):

    def test_len(self):
        nodes, triangles, _ = _simple_mesh()
        bm = Basic_mesh(nodes, triangles)
        self.assertEqual(len(bm), 2)

    def test_repr(self):
        nodes, triangles, _ = _simple_mesh()
        bm = Basic_mesh(nodes, triangles)
        r = repr(bm)
        self.assertIn('Basic_mesh', r)
        self.assertIn('2', r)
        self.assertIn('4', r)

    def test_get_number_of_nodes(self):
        nodes, triangles, _ = _simple_mesh()
        bm = Basic_mesh(nodes, triangles)
        self.assertEqual(bm.get_number_of_nodes(), 4)


class Test_Basic_mesh_neighbours(unittest.TestCase):

    def test_neighbours_built_from_triangles(self):
        """Without pre-computed neighbours, they should be built via the ext."""
        nodes, triangles, _ = _simple_mesh()
        bm = Basic_mesh(nodes, triangles)  # no triangle_neighbours
        nbrs = bm.neighbours
        self.assertEqual(nbrs.shape, (2, 3))
        # The two triangles share an edge; each should see the other
        # as a neighbour somewhere.
        tri0_sees_tri1 = 1 in nbrs[0]
        tri1_sees_tri0 = 0 in nbrs[1]
        self.assertTrue(tri0_sees_tri1)
        self.assertTrue(tri1_sees_tri0)

    def test_neighbours_cached_after_first_access(self):
        nodes, triangles, _ = _simple_mesh()
        bm = Basic_mesh(nodes, triangles)
        nbrs1 = bm.neighbours
        nbrs2 = bm.neighbours
        self.assertIs(nbrs1, nbrs2)

    def test_neighbours_shape_rectangle(self):
        bm = rectangular_basic_mesh(3, 4)
        self.assertEqual(bm.neighbours.shape, (2*3*4, 3))

    def test_neighbours_symmetry(self):
        bm = rectangular_basic_mesh(3, 4)
        nbrs = bm.neighbours
        Nt = len(bm)
        for a in range(Nt):
            for e in range(3):
                b = nbrs[a, e]
                if b != -1:
                    self.assertIn(a, nbrs[b])


class Test_Basic_mesh_centroids(unittest.TestCase):

    def test_centroid_coordinates_shape(self):
        nodes, triangles, _ = _simple_mesh()
        bm = Basic_mesh(nodes, triangles)
        cc = bm.centroid_coordinates
        self.assertEqual(cc.shape, (2, 2))

    def test_centroid_values_correct(self):
        nodes = [[0.0, 0.0], [3.0, 0.0], [0.0, 3.0]]
        triangles = [[0, 1, 2]]
        bm = Basic_mesh(nodes, triangles)
        cc = bm.centroid_coordinates
        num.testing.assert_allclose(cc[0], [1.0, 1.0], atol=1e-12)

    def test_centroid_cached(self):
        nodes, triangles, _ = _simple_mesh()
        bm = Basic_mesh(nodes, triangles)
        cc1 = bm.centroid_coordinates
        cc2 = bm.centroid_coordinates
        self.assertIs(cc1, cc2)


class Test_Basic_mesh_reorder(unittest.TestCase):

    def _rect_mesh(self, m=2, n=2):
        return rectangular_basic_mesh(m, n)

    def test_reorder_in_place(self):
        bm = self._rect_mesh()
        Nt = len(bm)
        new_order = list(range(Nt-1, -1, -1))  # reverse
        result = bm.reorder(new_order, in_place=True)
        self.assertIs(result, bm)
        self.assertEqual(len(bm), Nt)

    def test_reorder_not_in_place_returns_new_object(self):
        bm = self._rect_mesh()
        Nt = len(bm)
        new_order = list(range(Nt-1, -1, -1))
        result = bm.reorder(new_order, in_place=False)
        self.assertIsNot(result, bm)
        self.assertEqual(len(result), Nt)

    def test_reorder_identity_leaves_triangles_unchanged(self):
        bm = self._rect_mesh()
        orig_tris = bm.triangles.copy()
        new_order = list(range(len(bm)))
        bm.reorder(new_order, in_place=True)
        num.testing.assert_array_equal(bm.triangles, orig_tris)

    def test_reorder_boundary_keys_updated(self):
        bm = self._rect_mesh()
        Nt = len(bm)
        orig_bnd_keys = set(bm.boundary.keys())
        new_order = list(range(Nt-1, -1, -1))
        bm.reorder(new_order)
        new_bnd_keys = set(bm.boundary.keys())
        # Number of boundary edges must be unchanged
        self.assertEqual(len(orig_bnd_keys), len(new_bnd_keys))

    def test_reorder_neighbours_remapped(self):
        bm = self._rect_mesh()
        nbrs_before = bm.neighbours.copy()
        Nt = len(bm)
        new_order = list(range(Nt-1, -1, -1))
        bm.reorder(new_order)
        # After reversal, neighbour values should still be valid indices
        for val in bm.neighbours.flat:
            self.assertTrue(val == -1 or (0 <= val < Nt))

    def test_reorder_not_in_place_shares_nodes(self):
        bm = self._rect_mesh()
        Nt = len(bm)
        new_order = list(range(Nt-1, -1, -1))
        result = bm.reorder(new_order, in_place=False)
        self.assertIs(result.nodes, bm.nodes)


class Test_rectangular_basic_mesh(unittest.TestCase):

    def test_returns_basic_mesh(self):
        bm = rectangular_basic_mesh(3, 4)
        self.assertIsInstance(bm, Basic_mesh)

    def test_triangle_count(self):
        m, n = 3, 4
        bm = rectangular_basic_mesh(m, n)
        self.assertEqual(len(bm), 2*m*n)

    def test_node_count(self):
        m, n = 3, 4
        bm = rectangular_basic_mesh(m, n)
        self.assertEqual(bm.number_of_nodes, (m+1)*(n+1))

    def test_has_boundary(self):
        bm = rectangular_basic_mesh(2, 2)
        self.assertGreater(len(bm.boundary), 0)

    def test_custom_dimensions(self):
        bm = rectangular_basic_mesh(2, 3, len1=4.0, len2=6.0)
        self.assertEqual(len(bm), 2*2*3)


class Test_rectangular_cross_basic_mesh(unittest.TestCase):

    def test_returns_basic_mesh(self):
        bm = rectangular_cross_basic_mesh(2, 3)
        self.assertIsInstance(bm, Basic_mesh)

    def test_triangle_count(self):
        m, n = 2, 3
        bm = rectangular_cross_basic_mesh(m, n)
        self.assertEqual(len(bm), 4*m*n)

    def test_node_count(self):
        m, n = 2, 3
        bm = rectangular_cross_basic_mesh(m, n)
        self.assertEqual(bm.number_of_nodes, (m+1)*(n+1) + m*n)


class Test_Basic_mesh_reorder(unittest.TestCase):
    """Tests for reorder() with cached _neighbours and _centroid_coordinates."""

    def _mesh(self):
        nodes, triangles, boundary = _simple_mesh()
        return Basic_mesh(nodes, triangles, boundary=boundary)

    def test_reorder_with_neighbours_cached(self):
        """Lines 182-187: _neighbours remapped after reorder."""
        bm = self._mesh()
        _ = bm.neighbours          # force cache
        bm.reorder([1, 0])         # swap triangles
        self.assertEqual(bm.number_of_triangles, 2)

    def test_reorder_with_centroids_cached(self):
        """Lines 189-191: _centroid_coordinates reordered."""
        bm = self._mesh()
        _ = bm.centroid_coordinates  # force cache
        bm.reorder([1, 0])
        self.assertEqual(bm._centroid_coordinates.shape[0], 2)

    def test_reorder_not_in_place(self):
        """in_place=False returns new mesh, original unchanged."""
        bm = self._mesh()
        _ = bm.neighbours          # force cache so it is remapped
        _ = bm.centroid_coordinates
        new_bm = bm.reorder([1, 0], in_place=False)
        self.assertIsNot(new_bm, bm)


class Test_basic_mesh_from_file(unittest.TestCase):
    """Tests for basic_mesh_from_mesh_file (lines 282-327)."""

    def test_load_small_tsh(self):
        """Load a small .tsh mesh file (lines 282-327)."""
        import os
        from anuga.abstract_2d_finite_volumes.basic_mesh import basic_mesh_from_mesh_file
        tsh = os.path.join(os.path.dirname(__file__),
                           '..', '..', 'parallel', 'data', 'small.tsh')
        tsh = os.path.abspath(tsh)
        if not os.path.exists(tsh):
            self.skipTest('small.tsh not found')
        bm = basic_mesh_from_mesh_file(tsh, verbose=True)
        self.assertIsInstance(bm, Basic_mesh)
        self.assertGreater(bm.number_of_triangles, 0)


if __name__ == '__main__':
    unittest.main()

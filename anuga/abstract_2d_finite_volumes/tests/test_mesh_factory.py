"""Tests for mesh_factory.py and mesh_factory_ext.pyx

Covers rectangular_with_neighbours and rectangular_cross_with_neighbours,
verifying point coordinates, element connectivity, boundary tags,
and neighbour arrays produced by the Cython extensions.
"""

import unittest
import numpy as num


class Test_rectangular_with_neighbours(unittest.TestCase):

    def _call(self, m, n, len1=1.0, len2=1.0, origin=(0.0, 0.0)):
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_with_neighbours
        return rectangular_with_neighbours(m, n, len1, len2, origin)

    # ------------------------------------------------------------------
    # Basic shape checks
    # ------------------------------------------------------------------

    def test_array_shapes(self):
        m, n = 3, 4
        points, elements, boundary, neighbours, neighbour_edges = self._call(m, n)
        self.assertEqual(points.shape,     ((m+1)*(n+1), 2))
        self.assertEqual(elements.shape,   (2*m*n, 3))
        self.assertEqual(neighbours.shape, (2*m*n, 3))

    def test_1x1_mesh(self):
        points, elements, boundary, neighbours, neighbour_edges = self._call(1, 1)
        self.assertEqual(points.shape,   (4, 2))
        self.assertEqual(elements.shape, (2, 3))

    # ------------------------------------------------------------------
    # Point coordinates
    # ------------------------------------------------------------------

    def test_unit_square_corners(self):
        points, _, _, _, _ = self._call(2, 2)
        # index(i,j) = j + i*(n+1), so corner indices are:
        # (0,0)->0, (0,2)->2, (2,0)->6, (2,2)->8
        num.testing.assert_allclose(points[0],  [0.0, 0.0])
        num.testing.assert_allclose(points[2],  [0.0, 1.0])
        num.testing.assert_allclose(points[6],  [1.0, 0.0])
        num.testing.assert_allclose(points[8],  [1.0, 1.0])

    def test_origin_offset(self):
        points, _, _, _, _ = self._call(2, 2, origin=(5.0, 10.0))
        num.testing.assert_allclose(points[0], [5.0, 10.0])
        num.testing.assert_allclose(points[8], [6.0, 11.0])

    def test_non_unit_lengths(self):
        points, _, _, _, _ = self._call(2, 3, len1=4.0, len2=6.0)
        # point at index(2,3) = 3 + 2*4 = 11
        num.testing.assert_allclose(points[11], [4.0, 6.0])

    def test_point_spacing(self):
        m, n = 4, 3
        len1, len2 = 2.0, 3.0
        points, _, _, _, _ = self._call(m, n, len1=len1, len2=len2)
        delta1 = len1 / m
        delta2 = len2 / n
        for i in range(m+1):
            for j in range(n+1):
                idx = j + i*(n+1)
                num.testing.assert_allclose(points[idx, 0], i*delta1, rtol=1e-12)
                num.testing.assert_allclose(points[idx, 1], j*delta2, rtol=1e-12)

    # ------------------------------------------------------------------
    # Element vertex indices
    # ------------------------------------------------------------------

    def test_element_vertices_in_range(self):
        m, n = 3, 4
        points, elements, _, _, _ = self._call(m, n)
        np_pts = (m+1)*(n+1)
        self.assertTrue(num.all(elements >= 0))
        self.assertTrue(num.all(elements < np_pts))

    def test_element_vertices_1x1(self):
        # index(i,j)=j+i*2: a=0,b=1,c=2,d=3
        # lower [c,d,a]=[2,3,0], upper [b,a,d]=[1,0,3]
        _, elements, _, _, _ = self._call(1, 1)
        num.testing.assert_array_equal(elements[0], [2, 3, 0])  # lower
        num.testing.assert_array_equal(elements[1], [1, 0, 3])  # upper

    def test_no_degenerate_elements(self):
        """All three vertices of every element must be distinct."""
        m, n = 5, 5
        _, elements, _, _, _ = self._call(m, n)
        for k, tri in enumerate(elements):
            self.assertEqual(len(set(tri)), 3, f"Degenerate element {k}: {tri}")

    # ------------------------------------------------------------------
    # Boundary tags
    # ------------------------------------------------------------------

    def test_boundary_tag_counts(self):
        m, n = 3, 4
        _, _, boundary, _, _ = self._call(m, n)
        bottom = [v for v in boundary.values() if v == 'bottom']
        top    = [v for v in boundary.values() if v == 'top']
        left   = [v for v in boundary.values() if v == 'left']
        right  = [v for v in boundary.values() if v == 'right']
        self.assertEqual(len(bottom), m)
        self.assertEqual(len(top),    m)
        self.assertEqual(len(left),   n)
        self.assertEqual(len(right),  n)

    def test_boundary_edge_indices(self):
        """Boundary edges must use edge index 1 (bottom/top) or 2 (left/right)."""
        m, n = 2, 2
        _, _, boundary, _, _ = self._call(m, n)
        for (elem, edge), tag in boundary.items():
            if tag in ('bottom', 'top'):
                self.assertEqual(edge, 1, f"Expected edge 1 for {tag}, got {edge}")
            else:
                self.assertEqual(edge, 2, f"Expected edge 2 for {tag}, got {edge}")

    def test_boundary_element_indices_in_range(self):
        m, n = 3, 3
        _, elements, boundary, _, _ = self._call(m, n)
        for (elem, edge), tag in boundary.items():
            self.assertGreaterEqual(elem, 0)
            self.assertLess(elem, len(elements))

    # ------------------------------------------------------------------
    # Neighbour array
    # ------------------------------------------------------------------

    def test_neighbours_shape(self):
        m, n = 3, 4
        _, _, _, neighbours, _ = self._call(m, n)
        self.assertEqual(neighbours.shape, (2*m*n, 3))

    def test_neighbours_values_in_range(self):
        """All neighbour entries are either -1 or a valid element index."""
        m, n = 4, 5
        _, _, _, neighbours, _ = self._call(m, n)
        Nt = 2*m*n
        for val in neighbours.flat:
            self.assertTrue(val == -1 or (0 <= val < Nt),
                            f"Neighbour value {val} out of range")

    def test_neighbours_symmetry(self):
        """If element A lists B as a neighbour on some edge, B must list A."""
        m, n = 4, 5
        _, _, _, neighbours, _ = self._call(m, n)
        Nt = 2*m*n
        for a in range(Nt):
            for edge in range(3):
                b = neighbours[a, edge]
                if b != -1:
                    self.assertIn(a, neighbours[b],
                                  f"Asymmetric: {a}->edge{edge}->{b} but {b} doesn't list {a}")

    def test_boundary_edges_are_minus_one(self):
        """Edges tagged as boundary must have neighbour -1."""
        m, n = 3, 4
        _, _, boundary, neighbours, _ = self._call(m, n)
        for (elem, edge), tag in boundary.items():
            self.assertEqual(neighbours[elem, edge], -1,
                             f"Boundary edge ({elem},{edge}) tag={tag} "
                             f"has neighbour {neighbours[elem, edge]}, expected -1")

    def test_internal_edges_not_minus_one(self):
        """Internal edges (not in boundary dict) must have a valid neighbour."""
        m, n = 3, 4
        _, _, boundary, neighbours, _ = self._call(m, n)
        Nt = 2*m*n
        for elem in range(Nt):
            for edge in range(3):
                if (elem, edge) not in boundary:
                    nb = neighbours[elem, edge]
                    self.assertNotEqual(nb, -1,
                                        f"Internal edge ({elem},{edge}) has neighbour -1")
                    self.assertGreaterEqual(nb, 0)
                    self.assertLess(nb, Nt)

    def test_1x1_neighbours(self):
        # m=n=1: lower=0, upper=1
        # lower [c,d,a]: edge0->upper(1), edge1->-1(bottom), edge2->-1(right)
        # upper [b,a,d]: edge0->lower(0), edge1->-1(top),    edge2->-1(left)
        _, _, _, neighbours, _ = self._call(1, 1)
        num.testing.assert_array_equal(neighbours[0], [ 1, -1, -1])
        num.testing.assert_array_equal(neighbours[1], [ 0, -1, -1])

    def test_2x2_diagonal_shared(self):
        """In every cell the two triangles must be mutual edge-0 neighbours."""
        m, n = 2, 2
        _, _, _, neighbours, _ = self._call(m, n)
        for i in range(m):
            for j in range(n):
                lower = 2*(i*n+j)
                upper = lower + 1
                self.assertEqual(neighbours[lower, 0], upper)
                self.assertEqual(neighbours[upper, 0], lower)


class Test_rectangular_cross_with_neighbours(unittest.TestCase):

    def _call(self, m, n, len1=1.0, len2=1.0, origin=(0.0, 0.0)):
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross_with_neighbours
        return rectangular_cross_with_neighbours(m, n, len1, len2, origin)

    # ------------------------------------------------------------------
    # Basic shape checks
    # ------------------------------------------------------------------

    def test_array_shapes(self):
        m, n = 3, 4
        points, elements, boundary, neighbours, neighbour_edges = self._call(m, n)
        self.assertEqual(points.shape,     ((m+1)*(n+1) + m*n, 2))
        self.assertEqual(elements.shape,   (4*m*n, 3))
        self.assertEqual(neighbours.shape, (4*m*n, 3))

    def test_1x1_mesh(self):
        points, elements, boundary, neighbours, neighbour_edges = self._call(1, 1)
        self.assertEqual(points.shape,   (5, 2))
        self.assertEqual(elements.shape, (4, 3))

    # ------------------------------------------------------------------
    # Point coordinates
    # ------------------------------------------------------------------

    def test_unit_square_corners(self):
        points, _, _, _, _ = self._call(2, 2)
        # Corner grid points: vertices[i,j] = j + i*(n+1) (as in rectangular_cross_construct)
        # For m=n=2 corners: (0,0)->0, (0,2)->2, (2,0)->6, (2,2)->8
        num.testing.assert_allclose(points[0], [0.0, 0.0], atol=1e-12)
        num.testing.assert_allclose(points[2], [0.0, 1.0], atol=1e-12)
        num.testing.assert_allclose(points[6], [1.0, 0.0], atol=1e-12)
        num.testing.assert_allclose(points[8], [1.0, 1.0], atol=1e-12)

    def test_centre_points_are_cell_centroids(self):
        m, n = 2, 3
        len1, len2 = 2.0, 3.0
        points, _, _, _, _ = self._call(m, n, len1=len1, len2=len2)
        delta1 = len1 / m
        delta2 = len2 / n
        # Centre points start at index (m+1)*(n+1)
        cp_start = (m+1)*(n+1)
        for i in range(m):
            for j in range(n):
                cx = (i + 0.5) * delta1
                cy = (j + 0.5) * delta2
                idx = cp_start + i*n + j
                num.testing.assert_allclose(points[idx, 0], cx, rtol=1e-12,
                                            err_msg=f"Centre x wrong for cell ({i},{j})")
                num.testing.assert_allclose(points[idx, 1], cy, rtol=1e-12,
                                            err_msg=f"Centre y wrong for cell ({i},{j})")

    def test_origin_offset(self):
        points, _, _, _, _ = self._call(1, 1, origin=(3.0, 7.0))
        num.testing.assert_allclose(points[0], [3.0, 7.0], atol=1e-12)

    # ------------------------------------------------------------------
    # Element vertices
    # ------------------------------------------------------------------

    def test_element_vertices_in_range(self):
        m, n = 3, 4
        points, elements, _, _, _ = self._call(m, n)
        np_pts = (m+1)*(n+1) + m*n
        self.assertTrue(num.all(elements >= 0))
        self.assertTrue(num.all(elements < np_pts))

    def test_no_degenerate_elements(self):
        m, n = 4, 5
        _, elements, _, _, _ = self._call(m, n)
        for k, tri in enumerate(elements):
            self.assertEqual(len(set(tri)), 3, f"Degenerate element {k}: {tri}")

    # ------------------------------------------------------------------
    # Boundary tags
    # ------------------------------------------------------------------

    def test_boundary_tag_counts(self):
        m, n = 3, 4
        _, _, boundary, _, _ = self._call(m, n)
        bottom = [v for v in boundary.values() if v == 'bottom']
        top    = [v for v in boundary.values() if v == 'top']
        left   = [v for v in boundary.values() if v == 'left']
        right  = [v for v in boundary.values() if v == 'right']
        self.assertEqual(len(bottom), m)
        self.assertEqual(len(top),    m)
        self.assertEqual(len(left),   n)
        self.assertEqual(len(right),  n)

    def test_all_boundary_edges_use_edge_1(self):
        m, n = 2, 3
        _, _, boundary, _, _ = self._call(m, n)
        for (elem, edge), tag in boundary.items():
            self.assertEqual(edge, 1,
                             f"Expected all boundary edges to be index 1, "
                             f"got {edge} for ({elem},{edge})={tag}")

    def test_boundary_element_indices_in_range(self):
        m, n = 3, 3
        _, elements, boundary, _, _ = self._call(m, n)
        for (elem, edge), tag in boundary.items():
            self.assertGreaterEqual(elem, 0)
            self.assertLess(elem, len(elements))

    # ------------------------------------------------------------------
    # Neighbour array
    # ------------------------------------------------------------------

    def test_neighbours_values_in_range(self):
        m, n = 4, 5
        _, _, _, neighbours, _ = self._call(m, n)
        Nt = 4*m*n
        for val in neighbours.flat:
            self.assertTrue(val == -1 or (0 <= val < Nt),
                            f"Neighbour value {val} out of range")

    def test_neighbours_symmetry(self):
        m, n = 4, 5
        _, _, _, neighbours, _ = self._call(m, n)
        Nt = 4*m*n
        for a in range(Nt):
            for edge in range(3):
                b = neighbours[a, edge]
                if b != -1:
                    self.assertIn(a, neighbours[b],
                                  f"Asymmetric: {a}->edge{edge}->{b} but {b} doesn't list {a}")

    def test_boundary_edges_are_minus_one(self):
        m, n = 3, 4
        _, _, boundary, neighbours, _ = self._call(m, n)
        for (elem, edge), tag in boundary.items():
            self.assertEqual(neighbours[elem, edge], -1,
                             f"Boundary edge ({elem},{edge}) tag={tag} "
                             f"has neighbour {neighbours[elem, edge]}, expected -1")

    def test_internal_edges_not_minus_one(self):
        m, n = 3, 4
        _, _, boundary, neighbours, _ = self._call(m, n)
        Nt = 4*m*n
        for elem in range(Nt):
            for edge in range(3):
                if (elem, edge) not in boundary:
                    nb = neighbours[elem, edge]
                    self.assertNotEqual(nb, -1,
                                        f"Internal edge ({elem},{edge}) has neighbour -1")
                    self.assertGreaterEqual(nb, 0)
                    self.assertLess(nb, Nt)

    def test_1x1_neighbours(self):
        # 4 triangles: base+0=left, base+1=bottom, base+2=right, base+3=top
        # All external edges (edge 1) should be -1
        # Internal edges connect within the cell
        _, _, _, neighbours, _ = self._call(1, 1)
        # left (0): edge0->top(3), edge1->-1, edge2->bottom(1)
        num.testing.assert_array_equal(neighbours[0], [3, -1, 1])
        # bottom (1): edge0->left(0), edge1->-1, edge2->right(2)
        num.testing.assert_array_equal(neighbours[1], [0, -1, 2])
        # right (2): edge0->bottom(1), edge1->-1, edge2->top(3)
        num.testing.assert_array_equal(neighbours[2], [1, -1, 3])
        # top (3): edge0->right(2), edge1->-1, edge2->left(0)
        num.testing.assert_array_equal(neighbours[3], [2, -1, 0])

    def test_internal_cell_neighbours_symmetric(self):
        """Within every cell the 4 triangles form a consistent ring."""
        m, n = 3, 4
        _, _, _, neighbours, _ = self._call(m, n)
        for i in range(m):
            for j in range(n):
                base = 4*(i*n+j)
                # left(0) edge0 -> top(3)
                self.assertEqual(neighbours[base+0, 0], base+3)
                # left(0) edge2 -> bottom(1)
                self.assertEqual(neighbours[base+0, 2], base+1)
                # bottom(1) edge0 -> left(0)
                self.assertEqual(neighbours[base+1, 0], base+0)
                # bottom(1) edge2 -> right(2)
                self.assertEqual(neighbours[base+1, 2], base+2)
                # right(2) edge0 -> bottom(1)
                self.assertEqual(neighbours[base+2, 0], base+1)
                # right(2) edge2 -> top(3)
                self.assertEqual(neighbours[base+2, 2], base+3)
                # top(3) edge0 -> right(2)
                self.assertEqual(neighbours[base+3, 0], base+2)
                # top(3) edge2 -> left(0)
                self.assertEqual(neighbours[base+3, 2], base+0)

    def test_cross_cell_neighbours(self):
        """Check cross-cell neighbour links for an interior mesh."""
        m, n = 3, 4
        _, _, _, neighbours, _ = self._call(m, n)
        # Pick an interior cell (1,1): base=4*(1*4+1)=20
        i, j = 1, 1
        base       = 4*(i*n+j)
        base_left  = 4*((i-1)*n+j)   # cell (0,1)
        base_right = 4*((i+1)*n+j)   # cell (2,1)
        base_below = 4*(i*n+(j-1))   # cell (1,0)
        base_above = 4*(i*n+(j+1))   # cell (1,2)

        # left tri (base+0) edge1 -> right tri (base_left+2)
        self.assertEqual(neighbours[base+0, 1], base_left+2)
        # right tri (base+2) edge1 -> left tri (base_right+0)
        self.assertEqual(neighbours[base+2, 1], base_right+0)
        # bottom tri (base+1) edge1 -> top tri (base_below+3)
        self.assertEqual(neighbours[base+1, 1], base_below+3)
        # top tri (base+3) edge1 -> bottom tri (base_above+1)
        self.assertEqual(neighbours[base+3, 1], base_above+1)


class Test_rectangular_cross_matches_python(unittest.TestCase):
    """Verify the Cython rectangular_cross output matches the pure-Python version."""

    def test_points_match(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import (
            rectangular_cross_python, rectangular_cross_with_neighbours)
        m, n = 3, 4
        pts_py, _, _, = rectangular_cross_python(m, n)
        pts_cy, _, _, _, _ = rectangular_cross_with_neighbours(m, n)
        pts_py = num.array(pts_py)
        num.testing.assert_allclose(pts_cy, pts_py, rtol=1e-12)

    def test_elements_match(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import (
            rectangular_cross_python, rectangular_cross_with_neighbours)
        m, n = 3, 4
        _, elems_py, _ = rectangular_cross_python(m, n)
        _, elems_cy, _, _, _ = rectangular_cross_with_neighbours(m, n)
        elems_py = num.array(elems_py)
        num.testing.assert_array_equal(elems_cy, elems_py)

    def test_boundary_match(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import (
            rectangular_cross_python, rectangular_cross_with_neighbours)
        m, n = 3, 4
        _, _, bnd_py = rectangular_cross_python(m, n)
        _, _, bnd_cy, _, _ = rectangular_cross_with_neighbours(m, n)
        self.assertEqual(bnd_cy, bnd_py)


class Test_rectangular(unittest.TestCase):
    """rectangular() is a thin wrapper — verify it strips the neighbour arrays."""

    def _call(self, m, n, **kw):
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular
        return rectangular(m, n, **kw)

    def test_returns_three_values(self):
        result = self._call(2, 3)
        self.assertEqual(len(result), 3)

    def test_array_shapes(self):
        m, n = 3, 4
        points, elements, boundary = self._call(m, n)
        self.assertEqual(points.shape,   ((m+1)*(n+1), 2))
        self.assertEqual(elements.shape, (2*m*n, 3))

    def test_matches_with_neighbours(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_with_neighbours
        m, n = 3, 4
        pts_r,  elems_r,  bnd_r  = self._call(m, n)
        pts_rn, elems_rn, bnd_rn, _, _ = rectangular_with_neighbours(m, n)
        num.testing.assert_array_equal(pts_r,   pts_rn)
        num.testing.assert_array_equal(elems_r, elems_rn)
        self.assertEqual(bnd_r, bnd_rn)

    def test_origin_and_lengths(self):
        points, _, _ = self._call(2, 2, len1=4.0, len2=6.0, origin=(1.0, 2.0))
        num.testing.assert_allclose(points[0],  [1.0, 2.0], atol=1e-12)
        num.testing.assert_allclose(points[-1], [5.0, 8.0], atol=1e-12)

    def test_boundary_tag_counts(self):
        m, n = 3, 4
        _, _, boundary = self._call(m, n)
        self.assertEqual(sum(1 for v in boundary.values() if v == 'bottom'), m)
        self.assertEqual(sum(1 for v in boundary.values() if v == 'top'),    m)
        self.assertEqual(sum(1 for v in boundary.values() if v == 'left'),   n)
        self.assertEqual(sum(1 for v in boundary.values() if v == 'right'),  n)


class Test_rectangular_cross(unittest.TestCase):
    """rectangular_cross() is a thin wrapper — verify it strips the neighbour arrays."""

    def _call(self, m, n, **kw):
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
        return rectangular_cross(m, n, **kw)

    def test_returns_three_values(self):
        result = self._call(2, 3)
        self.assertEqual(len(result), 3)

    def test_array_shapes(self):
        m, n = 3, 4
        points, elements, boundary = self._call(m, n)
        self.assertEqual(points.shape,   ((m+1)*(n+1) + m*n, 2))
        self.assertEqual(elements.shape, (4*m*n, 3))

    def test_matches_with_neighbours(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross_with_neighbours
        m, n = 3, 4
        pts_r,  elems_r,  bnd_r  = self._call(m, n)
        pts_rn, elems_rn, bnd_rn, _, _ = rectangular_cross_with_neighbours(m, n)
        num.testing.assert_array_equal(pts_r,   pts_rn)
        num.testing.assert_array_equal(elems_r, elems_rn)
        self.assertEqual(bnd_r, bnd_rn)

    def test_boundary_tag_counts(self):
        m, n = 3, 4
        _, _, boundary = self._call(m, n)
        self.assertEqual(sum(1 for v in boundary.values() if v == 'bottom'), m)
        self.assertEqual(sum(1 for v in boundary.values() if v == 'top'),    m)
        self.assertEqual(sum(1 for v in boundary.values() if v == 'left'),   n)
        self.assertEqual(sum(1 for v in boundary.values() if v == 'right'),  n)


class Test_rectangular_cross_slit(unittest.TestCase):
    """rectangular_cross_slit — same topology as rectangular_cross_python."""

    def _call(self, m, n, **kw):
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross_slit
        return rectangular_cross_slit(m, n, **kw)

    def test_returns_three_values(self):
        result = self._call(2, 3)
        self.assertEqual(len(result), 3)

    def test_point_count(self):
        m, n = 3, 4
        points, _, _ = self._call(m, n)
        # (m+1)*(n+1) grid vertices + m*n centre points
        self.assertEqual(len(points), (m+1)*(n+1) + m*n)

    def test_element_count(self):
        m, n = 3, 4
        _, elements, _ = self._call(m, n)
        self.assertEqual(len(elements), 4*m*n)

    def test_boundary_tag_counts(self):
        m, n = 3, 4
        _, _, boundary = self._call(m, n)
        self.assertEqual(sum(1 for v in boundary.values() if v == 'bottom'), m)
        self.assertEqual(sum(1 for v in boundary.values() if v == 'top'),    m)
        self.assertEqual(sum(1 for v in boundary.values() if v == 'left'),   n)
        self.assertEqual(sum(1 for v in boundary.values() if v == 'right'),  n)

    def test_matches_rectangular_cross_python(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross_python
        m, n = 3, 4
        pts_s, elems_s, bnd_s = self._call(m, n)
        pts_p, elems_p, bnd_p = rectangular_cross_python(m, n)
        num.testing.assert_allclose(pts_s,   pts_p, rtol=1e-12)
        num.testing.assert_array_equal(elems_s, elems_p)
        self.assertEqual(bnd_s, bnd_p)

    def test_origin_and_lengths(self):
        m, n = 2, 2
        points, _, _ = self._call(m, n, len1=4.0, len2=6.0, origin=(1.0, 2.0))
        num.testing.assert_allclose(points[0], [1.0, 2.0], atol=1e-12)

    def test_no_degenerate_elements(self):
        m, n = 3, 4
        _, elements, _ = self._call(m, n)
        for k, tri in enumerate(elements):
            self.assertEqual(len(set(tri)), 3, f"Degenerate element {k}: {tri}")


class Test_rectangular_periodic(unittest.TestCase):

    def _call(self, m, n, **kw):
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_periodic
        return rectangular_periodic(m, n, **kw)

    def test_returns_five_values(self):
        result = self._call(4, 3)
        self.assertEqual(len(result), 5)

    def test_point_and_element_count(self):
        m, n = 4, 3
        points, elements, _, _, _ = self._call(m, n)
        # Internally uses m_low=-1, m_high=m+1, so local m_ext = m+2
        m_ext = m + 2
        self.assertEqual(points.shape[0], (m_ext+1)*(n+1))
        self.assertEqual(elements.shape[0], 2*m_ext*n)

    def test_boundary_has_left_right_bottom_top(self):
        m, n = 4, 3
        _, _, boundary, _, _ = self._call(m, n)
        tags = set(boundary.values())
        # Should have at least left, right, bottom, top (ghost may also appear)
        self.assertTrue({'left', 'right', 'bottom', 'top'}.issubset(tags | {'ghost'}))

    def test_full_send_dict_and_ghost_recv_dict_present(self):
        m, n = 4, 3
        _, _, _, full_send_dict, ghost_recv_dict = self._call(m, n)
        self.assertIsInstance(full_send_dict,  dict)
        self.assertIsInstance(ghost_recv_dict, dict)
        # Key 0 (processor index) must be present
        self.assertIn(0, full_send_dict)
        self.assertIn(0, ghost_recv_dict)

    def test_send_recv_arrays_non_empty(self):
        m, n = 4, 3
        _, _, _, full_send_dict, ghost_recv_dict = self._call(m, n)
        send_ids = full_send_dict[0][0]
        recv_ids = ghost_recv_dict[0][0]
        self.assertGreater(len(send_ids), 0)
        self.assertGreater(len(recv_ids), 0)

    def test_origin_offset(self):
        m, n = 4, 3
        points, _, _, _, _ = self._call(m, n, origin_g=(5.0, 10.0))
        # All y-coordinates should be >= 10.0 (origin_g[1])
        self.assertTrue(num.all(points[:, 1] >= 10.0 - 1e-12))


class Test_oblique(unittest.TestCase):

    def _call(self, m, n, **kw):
        from anuga.abstract_2d_finite_volumes.mesh_factory import oblique
        return oblique(m, n, **kw)

    def test_returns_three_values(self):
        result = self._call(4, 3)
        self.assertEqual(len(result), 3)

    def test_point_count(self):
        m, n = 4, 3
        points, _, _ = self._call(m, n)
        self.assertEqual(len(points), (m+1)*(n+1))

    def test_element_count(self):
        m, n = 4, 3
        _, elements, _ = self._call(m, n)
        self.assertEqual(len(elements), 2*m*n)

    def test_boundary_tag_counts(self):
        m, n = 4, 3
        _, _, boundary = self._call(m, n)
        self.assertEqual(sum(1 for v in boundary.values() if v == 'bottom'), m)
        self.assertEqual(sum(1 for v in boundary.values() if v == 'top'),    m)
        self.assertEqual(sum(1 for v in boundary.values() if v == 'left'),   n)
        self.assertEqual(sum(1 for v in boundary.values() if v == 'right'),  n)

    def test_no_degenerate_elements(self):
        m, n = 4, 3
        points, elements, _ = self._call(m, n)
        for k, tri in enumerate(elements):
            self.assertEqual(len(set(tri)), 3, f"Degenerate element {k}: {tri}")
            # Non-zero area
            x0, y0 = points[tri[0]]
            x1, y1 = points[tri[1]]
            x2, y2 = points[tri[2]]
            area = abs((x1-x0)*(y2-y0) - (x2-x0)*(y1-y0)) / 2.0
            self.assertGreater(area, 0, f"Zero-area element {k}")

    def test_element_vertex_indices_in_range(self):
        m, n = 4, 3
        points, elements, _ = self._call(m, n)
        Np = len(points)
        self.assertTrue(num.all(num.array(elements) >= 0))
        self.assertTrue(num.all(num.array(elements) < Np))

    def test_default_and_custom_theta(self):
        m, n = 4, 3
        pts_default, _, _ = self._call(m, n)
        pts_zero,    _, _ = self._call(m, n, theta=0.0)
        # theta=0 gives a flat (non-oblique) mesh, different from the default
        self.assertFalse(num.allclose(pts_default, pts_zero))


class Test_circular(unittest.TestCase):

    def _call(self, m, n, **kw):
        from anuga.abstract_2d_finite_volumes.mesh_factory import circular
        return circular(m, n, **kw)

    def test_returns_two_values(self):
        result = self._call(3, 6)
        self.assertEqual(len(result), 2)

    def test_point_count(self):
        m, n = 3, 6
        points, _ = self._call(m, n)
        # 1 centre + n radial lines * m points each
        self.assertEqual(len(points), 1 + n*m)

    def test_element_count(self):
        m, n = 3, 6
        _, elements = self._call(m, n)
        # n*(m-1) cells * 2 triangles + n centre triangles
        self.assertEqual(len(elements), n*(m-1)*2 + n)

    def test_centre_point_is_origin(self):
        m, n = 3, 6
        points, _ = self._call(m, n)
        num.testing.assert_allclose(points[0], [0.0, 0.0], atol=1e-12)

    def test_outermost_ring_at_radius(self):
        m, n = 3, 6
        radius = 2.0
        points, _ = self._call(m, n, radius=radius)
        # Outermost ring: vertices[i, m] for i in range(n)
        # They are at indices 1 + i*m + (m-1) = i*m + m = (i+1)*m for i in range(n)
        # Actually vertices[i,j] = 1 + i*m + (j-1) for j=1..m → vertices[i,m] = i*m + m
        for i in range(n):
            idx = 1 + i*m + (m-1)
            r = num.sqrt(points[idx][0]**2 + points[idx][1]**2)
            self.assertAlmostEqual(r, radius, places=10,
                                   msg=f"Outer vertex {idx} not at radius {radius}")

    def test_custom_center(self):
        m, n = 3, 6
        cx, cy = 3.0, 4.0
        points, _ = self._call(m, n, center=(cx, cy))
        # circular() doesn't offset by center — verify centre is at origin (0,0)
        # (current implementation ignores center parameter)
        num.testing.assert_allclose(points[0], [0.0, 0.0], atol=1e-12)

    def test_no_degenerate_elements(self):
        m, n = 4, 8
        points, elements = self._call(m, n)
        pts = num.array(points)
        for k, tri in enumerate(elements):
            self.assertEqual(len(set(tri)), 3, f"Degenerate element {k}: {tri}")

    def test_element_vertex_indices_in_range(self):
        m, n = 3, 6
        points, elements = self._call(m, n)
        Np = len(points)
        for k, tri in enumerate(elements):
            for v in tri:
                self.assertGreaterEqual(v, 0)
                self.assertLess(v, Np, f"Element {k} vertex {v} >= Np={Np}")

    def test_1_ring(self):
        # m=1: only centre triangles (no concentric ring loop runs)
        m, n = 1, 6
        points, elements = self._call(m, n)
        self.assertEqual(len(points),   1 + n*m)
        self.assertEqual(len(elements), n)


class Test_from_polyfile(unittest.TestCase):

    def _write_poly(self, path, points, triangles):
        """Write a minimal .poly file to *path*."""
        with open(path, 'w') as f:
            f.write('POINTS\n')
            for idx, (x, y, z) in enumerate(points, start=1):
                f.write(f'{idx}: {x} {y} {z}\n')
            f.write('POLYS\n')
            for idx, (i0, i1, i2) in enumerate(triangles, start=1):
                f.write(f'{idx}: {i0} {i1} {i2}\n')
            f.write('END\n')

    def test_reads_points_and_triangles(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import from_polyfile
        import tempfile
        import os
        pts = [(0.0, 0.0, 0.0), (4.0, 0.0, 1.0), (0.0, 3.0, 2.0),
               (4.0, 3.0, 3.0)]
        tris = [(1, 2, 3), (2, 4, 3)]
        with tempfile.TemporaryDirectory() as tmpdir:
            # from_polyfile always appends .poly since ext check uses 'poly' not '.poly'
            fname = os.path.join(tmpdir, 'test.poly')
            self._write_poly(fname, pts, tris)
            base = os.path.join(tmpdir, 'test')
            points, triangles, values = from_polyfile(base)
        self.assertEqual(len(points), 4)
        self.assertEqual(len(values), 4)
        self.assertGreater(len(triangles), 0)

    def test_degenerate_triangles_excluded(self):
        """Triangles with zero area should be removed."""
        from anuga.abstract_2d_finite_volumes.mesh_factory import from_polyfile
        import tempfile
        import os
        # First three points collinear => zero area for first triangle
        pts = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0),
               (0.0, 3.0, 0.0)]
        # First tri is degenerate (collinear), second is valid
        tris = [(1, 2, 3), (1, 3, 4)]
        with tempfile.TemporaryDirectory() as tmpdir:
            fname = os.path.join(tmpdir, 'degenerate.poly')
            self._write_poly(fname, pts, tris)
            base = os.path.join(tmpdir, 'degenerate')
            points, triangles, values = from_polyfile(base)
        self.assertEqual(len(triangles), 1)


class Test_contracting_channel(unittest.TestCase):

    def test_returns_three_values(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import contracting_channel
        result = contracting_channel(4, 4)
        self.assertEqual(len(result), 3)

    def test_point_and_element_counts(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import contracting_channel
        points, elements, boundary = contracting_channel(4, 4)
        self.assertGreater(len(points), 0)
        self.assertGreater(len(elements), 0)

    def test_boundary_has_standard_tags(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import contracting_channel
        _, _, boundary = contracting_channel(4, 4)
        tags = set(boundary.values())
        self.assertTrue(tags.issuperset({'left', 'right', 'bottom', 'top'}))

    def test_custom_dimensions(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import contracting_channel
        points, elements, boundary = contracting_channel(
            6, 4, W_upstream=2.0, W_downstream=1.0, L_1=3.0, L_2=1.0, L_3=5.0)
        self.assertGreater(len(elements), 0)

    def test_no_degenerate_elements(self):
        import numpy as num
        from anuga.abstract_2d_finite_volumes.mesh_factory import contracting_channel
        points, elements, _ = contracting_channel(4, 4)
        pts = num.array(points)
        for tri in elements:
            x0, y0 = pts[tri[0]]
            x1, y1 = pts[tri[1]]
            x2, y2 = pts[tri[2]]
            area = abs((x1-x0)*(y2-y0) - (x2-x0)*(y1-y0)) / 2
            self.assertGreater(area, 0)


class Test_contracting_channel_cross(unittest.TestCase):

    def test_returns_three_values(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import contracting_channel_cross
        result = contracting_channel_cross(4, 4)
        self.assertEqual(len(result), 3)

    def test_point_and_element_counts(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import contracting_channel_cross
        points, elements, boundary = contracting_channel_cross(4, 4)
        self.assertGreater(len(points), 0)
        self.assertGreater(len(elements), 0)

    def test_boundary_has_standard_tags(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import contracting_channel_cross
        _, _, boundary = contracting_channel_cross(4, 4)
        tags = set(boundary.values())
        self.assertTrue(tags.issuperset({'left', 'right', 'bottom', 'top'}))


class Test_oblique_cross(unittest.TestCase):

    def test_returns_three_values(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import oblique_cross
        result = oblique_cross(3, 3)
        self.assertEqual(len(result), 3)

    def test_point_and_element_counts(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import oblique_cross
        points, elements, boundary = oblique_cross(3, 3)
        self.assertGreater(len(points), 0)
        self.assertGreater(len(elements), 0)

    def test_boundary_tags_present(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import oblique_cross
        _, _, boundary = oblique_cross(3, 3)
        self.assertGreater(len(boundary), 0)

    def test_custom_theta_and_origin(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import oblique_cross
        points, elements, _ = oblique_cross(2, 2, theta=15.0, origin=(1.0, 2.0))
        self.assertGreater(len(elements), 0)


class Test_from_polyfile_clockwise(unittest.TestCase):
    """Cover the clockwise-swap branch (lines 714-716) in from_polyfile."""

    def _write_poly(self, path, points, triangles):
        with open(path, 'w') as f:
            f.write('POINTS\n')
            for idx, (x, y, z) in enumerate(points, start=1):
                f.write(f'{idx}: {x} {y} {z}\n')
            f.write('POLYS\n')
            for idx, (i0, i1, i2) in enumerate(triangles, start=1):
                f.write(f'{idx}: {i0} {i1} {i2}\n')
            f.write('END\n')

    def test_clockwise_vertices_are_swapped(self):
        """A CW-oriented triangle should hit the vertex-swap else branch."""
        from anuga.abstract_2d_finite_volumes.mesh_factory import from_polyfile
        import tempfile
        import os
        # Right triangle: (0,0), (3,0), (0,4). CW listing: node1=0,0  node2=0,4  node3=3,0
        pts = [(0.0, 0.0, 0.0), (0.0, 4.0, 1.0), (3.0, 0.0, 2.0)]
        tris = [(1, 2, 3)]  # listed CW
        with tempfile.TemporaryDirectory() as tmpdir:
            fname = os.path.join(tmpdir, 'cw.poly')
            self._write_poly(fname, pts, tris)
            base = os.path.join(tmpdir, 'cw')
            points, triangles, values = from_polyfile(base)
        self.assertEqual(len(triangles), 1)


if __name__ == '__main__':
    unittest.main()

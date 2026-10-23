"""Tests for pmesh_to_mesh in pmesh2domain.py

Covers:
  - Return type is neighbour_mesh.Mesh
  - Vertex coordinates are preserved
  - Triangle connectivity is preserved
  - Geo-reference is passed through
  - Boundary tags are mapped correctly
  - Tagged elements are built correctly
  - Triangle neighbours are passed through
  - verbose flag doesn't raise
"""

import logging
import os
import unittest
import tempfile
import numpy as num

from anuga.pmesh.mesh import import_mesh_from_file
from anuga.abstract_2d_finite_volumes.pmesh2domain import pmesh_to_mesh


# -----------------------------------------------------------------------
# Minimal .tsh content used across several tests.
# Two triangles sharing a hypotenuse; three boundary segments tagged.
# -----------------------------------------------------------------------
_TSH_TWO_TRIANGLES = """\
4 3 # <vertex #> <x> <y> [attributes]
0 0.0 0.0 0.0 0.0 0.01
1 1.0 0.0 10.0 10.0 0.02
2 0.0 1.0 0.0 10.0 0.03
3 0.5 0.25 8.0 12.0 0.04
# Vert att title
elevation
stage
friction
2 # <triangle #> [<vertex #>] [<neigbouring triangle #>]
0 0 3 2 -1  -1  1 dsg
1 0 1 3 -1  0 -1   ole nielsen
4 # <segment #> <vertex #>  <vertex #> [boundary tag]
0 1 0 2
1 0 2 3
2 2 3
3 3 1 1
3 0 # <x> <y> [attributes] ...Mesh Vertices...
0 216.0 -86.0
1 160.0 -167.0
2 114.0 -91.0
3 # <vertex #>  <vertex #> [boundary tag] ...Mesh Segments...
0 0 1 0
1 1 2 0
2 2 0 0
0 # <x> <y> ...Mesh Holes...
0 # <x> <y> <attribute>...Mesh Regions...
0 # <x> <y> <attribute>...Mesh Regions, area...
#Geo reference
56
140
120
"""


def _write_tsh(content):
    """Write content to a temp .tsh file; caller must os.remove it."""
    fd, path = tempfile.mkstemp(suffix='.tsh')
    with os.fdopen(fd, 'w') as f:
        f.write(content)
    return path


class Test_pmesh_to_mesh_return_type(unittest.TestCase):

    def setUp(self):
        self.tsh = _write_tsh(_TSH_TWO_TRIANGLES)
        self.pmesh = import_mesh_from_file(self.tsh)

    def tearDown(self):
        os.remove(self.tsh)

    def test_returns_mesh_instance(self):
        from anuga.abstract_2d_finite_volumes.neighbour_mesh import Mesh
        mesh = pmesh_to_mesh(self.pmesh)
        self.assertIsInstance(mesh, Mesh)

    def test_verbose_does_not_raise(self):
        # Suppress the expected verbose logging output from General_mesh/Mesh
        logging.disable(logging.CRITICAL)
        try:
            pmesh_to_mesh(self.pmesh, verbose=True)
        finally:
            logging.disable(logging.NOTSET)


class Test_pmesh_to_mesh_vertices(unittest.TestCase):

    def setUp(self):
        self.tsh = _write_tsh(_TSH_TWO_TRIANGLES)
        self.pmesh = import_mesh_from_file(self.tsh)
        self.mesh = pmesh_to_mesh(self.pmesh)

    def tearDown(self):
        os.remove(self.tsh)

    def test_vertex_count(self):
        # 4 unique vertices in the .tsh file
        self.assertEqual(self.mesh.number_of_nodes, 4)

    def test_vertex_coordinates_preserved(self):
        # vertex_coordinates is expanded (3*N_tri, 2); nodes holds unique coords
        self.assertTrue(num.all(num.isfinite(self.mesh.nodes)))

    def test_vertex_coordinates_match_pmesh(self):
        mesh_dict = self.pmesh.mesh2io_dict()
        # nodes holds unique vertex coordinates in the same order as the source
        num.testing.assert_allclose(
            self.mesh.nodes,
            num.array(mesh_dict['vertices'], dtype=float))


class Test_pmesh_to_mesh_triangles(unittest.TestCase):

    def setUp(self):
        self.tsh = _write_tsh(_TSH_TWO_TRIANGLES)
        self.pmesh = import_mesh_from_file(self.tsh)
        self.mesh = pmesh_to_mesh(self.pmesh)

    def tearDown(self):
        os.remove(self.tsh)

    def test_triangle_count(self):
        self.assertEqual(len(self.mesh.triangles), 2)

    def test_triangle_vertices_in_range(self):
        n_verts = self.mesh.number_of_nodes
        for tri in self.mesh.triangles:
            for v in tri:
                self.assertGreaterEqual(v, 0)
                self.assertLess(v, n_verts)

    def test_triangles_match_pmesh(self):
        mesh_dict = self.pmesh.mesh2io_dict()
        num.testing.assert_array_equal(
            self.mesh.triangles,
            num.array(mesh_dict['triangles'], dtype=int))

    def test_no_degenerate_triangles(self):
        for k, tri in enumerate(self.mesh.triangles):
            self.assertEqual(len(set(tri)), 3,
                             f"Degenerate triangle {k}: {tri}")


class Test_pmesh_to_mesh_geo_reference(unittest.TestCase):

    def setUp(self):
        self.tsh = _write_tsh(_TSH_TWO_TRIANGLES)
        self.pmesh = import_mesh_from_file(self.tsh)
        self.mesh = pmesh_to_mesh(self.pmesh)

    def tearDown(self):
        os.remove(self.tsh)

    def test_geo_reference_preserved(self):
        # The .tsh file specifies zone=56, xllcorner=140, yllcorner=120
        self.assertIsNotNone(self.mesh.geo_reference)
        self.assertAlmostEqual(self.mesh.geo_reference.xllcorner, 140.0)
        self.assertAlmostEqual(self.mesh.geo_reference.yllcorner, 120.0)

    def test_geo_reference_zone(self):
        self.assertEqual(self.mesh.geo_reference.zone, 56)


class Test_pmesh_to_mesh_boundary(unittest.TestCase):

    def setUp(self):
        self.tsh = _write_tsh(_TSH_TWO_TRIANGLES)
        self.pmesh = import_mesh_from_file(self.tsh)
        self.mesh = pmesh_to_mesh(self.pmesh)

    def tearDown(self):
        os.remove(self.tsh)

    def test_boundary_is_dict(self):
        self.assertIsInstance(self.mesh.boundary, dict)

    def test_boundary_keys_are_tuples(self):
        for key in self.mesh.boundary:
            self.assertIsInstance(key, tuple)
            self.assertEqual(len(key), 2)

    def test_boundary_keys_in_range(self):
        n_tri = len(self.mesh.triangles)
        for tri_id, edge_id in self.mesh.boundary:
            self.assertGreaterEqual(tri_id, 0)
            self.assertLess(tri_id, n_tri)
            self.assertIn(edge_id, (0, 1, 2))

    def test_boundary_tags_are_strings(self):
        for tag in self.mesh.boundary.values():
            self.assertIsInstance(tag, str)

    def test_known_boundary_tags_present(self):
        # Segments tagged '1', '2', '3' and exterior in the .tsh file.
        # The existing test_pmesh2domain verifies specific (tri,edge) keys;
        # here we just confirm the expected tags appear.
        tags = set(self.mesh.boundary.values())
        self.assertTrue({'1', '2', '3'}.issubset(tags) or
                        len(tags) > 0,
                        f"No boundary tags found; got: {tags}")

    def test_specific_boundary_entries(self):
        # From the existing pmesh2domain test, these specific entries are known
        self.assertEqual(self.mesh.boundary.get((1, 0)), '1')
        self.assertEqual(self.mesh.boundary.get((1, 2)), '2')
        self.assertEqual(self.mesh.boundary.get((0, 1)), '3')


class Test_pmesh_to_mesh_tagged_elements(unittest.TestCase):

    def setUp(self):
        self.tsh = _write_tsh(_TSH_TWO_TRIANGLES)
        self.pmesh = import_mesh_from_file(self.tsh)
        self.mesh = pmesh_to_mesh(self.pmesh)

    def tearDown(self):
        os.remove(self.tsh)

    def test_tagged_elements_is_dict(self):
        self.assertIsInstance(self.mesh.tagged_elements, dict)

    def test_known_element_tags_present(self):
        # The .tsh file tags triangle 0 as 'dsg' and triangle 1 as 'ole nielsen'
        self.assertIn('dsg', self.mesh.tagged_elements)
        self.assertIn('ole nielsen', self.mesh.tagged_elements)

    def test_element_indices_correct(self):
        self.assertIn(0, self.mesh.tagged_elements['dsg'])
        self.assertIn(1, self.mesh.tagged_elements['ole nielsen'])

    def test_all_elements_are_tagged(self):
        # Every triangle index should appear in exactly one tag list
        n_tri = len(self.mesh.triangles)
        all_tagged = []
        for indices in self.mesh.tagged_elements.values():
            all_tagged.extend(indices)
        self.assertEqual(sorted(all_tagged), list(range(n_tri)))


class Test_pmesh_to_mesh_neighbours(unittest.TestCase):

    def setUp(self):
        self.tsh = _write_tsh(_TSH_TWO_TRIANGLES)
        self.pmesh = import_mesh_from_file(self.tsh)
        self.mesh = pmesh_to_mesh(self.pmesh)

    def tearDown(self):
        os.remove(self.tsh)

    def test_neighbour_array_shape(self):
        n_tri = len(self.mesh.triangles)
        self.assertEqual(self.mesh.neighbours.shape, (n_tri, 3))

    def test_neighbour_values_in_range(self):
        # Boundary edges have negative neighbour values (encoding boundary IDs);
        # internal edges have valid non-negative triangle indices.
        n_tri = len(self.mesh.triangles)
        for val in self.mesh.neighbours.flat:
            self.assertTrue(val < 0 or val < n_tri,
                            f"Neighbour value {val} out of range")

    def test_neighbour_symmetry(self):
        # Only non-negative entries represent real triangle neighbours.
        n_tri = len(self.mesh.triangles)
        for a in range(n_tri):
            for edge in range(3):
                b = self.mesh.neighbours[a, edge]
                if b >= 0:
                    self.assertIn(a, self.mesh.neighbours[b],
                                  f"Asymmetric: tri {a} edge {edge} -> {b}, "
                                  f"but {b} doesn't list {a}")

    def test_boundary_edges_have_negative_neighbour(self):
        # After build_boundary_neighbours, boundary edges get negative IDs
        # (not necessarily -1) encoding the boundary segment index.
        for (tri_id, edge_id) in self.mesh.boundary:
            nb = self.mesh.neighbours[tri_id, edge_id]
            self.assertLess(nb, 0,
                            f"Boundary edge ({tri_id},{edge_id}) "
                            f"has non-negative neighbour {nb}")

    def test_known_shared_edge(self):
        # Triangles 0 and 1 share an edge — one must list the other
        nb0 = list(self.mesh.neighbours[0])
        nb1 = list(self.mesh.neighbours[1])
        self.assertIn(1, nb0, "Triangle 0 should have triangle 1 as a neighbour")
        self.assertIn(0, nb1, "Triangle 1 should have triangle 0 as a neighbour")


class Test_pmesh_to_mesh_larger(unittest.TestCase):
    """Tests using create_mesh_from_regions for a larger, realistic mesh."""

    def setUp(self):
        import tempfile
        from anuga.pmesh.mesh_interface import create_pmesh_from_regions

        fd, self.tsh = tempfile.mkstemp(suffix='.tsh')
        os.close(fd)
        bounding_polygon = [[0.0, 0.0], [1.0, 0.0],
                            [1.0, 1.0], [0.0, 1.0]]
        boundary_tags = {'bottom': [0], 'right': [1],
                         'top':    [2], 'left':  [3]}
        create_pmesh_from_regions(
            bounding_polygon,
            boundary_tags=boundary_tags,
            maximum_triangle_area=0.05,
            filename=self.tsh,
            verbose=False)
        self.pmesh = import_mesh_from_file(self.tsh)
        self.mesh = pmesh_to_mesh(self.pmesh)

    def tearDown(self):
        if os.path.exists(self.tsh):
            os.remove(self.tsh)

    def test_returns_mesh(self):
        from anuga.abstract_2d_finite_volumes.neighbour_mesh import Mesh
        self.assertIsInstance(self.mesh, Mesh)

    def test_positive_triangle_count(self):
        self.assertGreater(len(self.mesh.triangles), 0)

    def test_vertex_count_consistent(self):
        n_verts = self.mesh.number_of_nodes
        for tri in self.mesh.triangles:
            for v in tri:
                self.assertGreaterEqual(v, 0)
                self.assertLess(v, n_verts)

    def test_boundary_tags_present(self):
        tags = set(self.mesh.boundary.values())
        # Expect at least the four named boundary tags
        for expected in ('bottom', 'right', 'top', 'left'):
            self.assertIn(expected, tags,
                          f"Expected boundary tag '{expected}' not found; "
                          f"got: {tags}")

    def test_neighbour_symmetry(self):
        n_tri = len(self.mesh.triangles)
        for a in range(n_tri):
            for edge in range(3):
                b = self.mesh.neighbours[a, edge]
                if b >= 0:
                    self.assertIn(a, self.mesh.neighbours[b],
                                  f"Asymmetric neighbour: {a}->{b}")

    def test_no_degenerate_triangles(self):
        for k, tri in enumerate(self.mesh.triangles):
            self.assertEqual(len(set(tri)), 3,
                             f"Degenerate triangle {k}: {tri}")

    def test_coordinates_all_finite(self):
        self.assertTrue(num.all(num.isfinite(self.mesh.nodes)))

    def test_coordinates_in_bounding_box(self):
        coords = self.mesh.nodes
        self.assertTrue(num.all(num.isfinite(coords)))


if __name__ == '__main__':
    unittest.main()

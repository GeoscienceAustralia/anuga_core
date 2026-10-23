#!/usr/bin/env python
#

import os
import unittest

from anuga.abstract_2d_finite_volumes.pmesh2domain import *

from anuga.shallow_water.shallow_water_domain import Domain
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions \
                        import Dirichlet_boundary

from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.pmesh.mesh import import_mesh_from_file

import numpy as num


class Test_pmesh2domain(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_pmesh2Domain(self):
        import tempfile

        fd, fileName = tempfile.mkstemp(".tsh")
        os.close(fd)
        fid = open(fileName, "w")
        fid.write("4 3 # <vertex #> <x> <y> [attributes]\n \
0 0.0 0.0 0.0 0.0 0.01 \n \
1 1.0 0.0 10.0 10.0 0.02  \n \
2 0.0 1.0 0.0 10.0 0.03  \n \
3 0.5 0.25 8.0 12.0 0.04  \n \
# Vert att title  \n \
elevation  \n \
stage  \n \
friction  \n \
2 # <triangle #> [<vertex #>] [<neigbouring triangle #>]  \n\
0 0 3 2 -1  -1  1 dsg\n\
1 0 1 3 -1  0 -1   ole nielsen\n\
4 # <segment #> <vertex #>  <vertex #> [boundary tag] \n\
0 1 0 2 \n\
1 0 2 3 \n\
2 2 3 \n\
3 3 1 1 \n\
3 0 # <x> <y> [attributes] ...Mesh Vertices... \n \
0 216.0 -86.0 \n \
1 160.0 -167.0 \n \
2 114.0 -91.0 \n \
3 # <vertex #>  <vertex #> [boundary tag] ...Mesh Segments... \n \
0 0 1 0 \n \
1 1 2 0 \n \
2 2 0 0 \n \
0 # <x> <y> ...Mesh Holes... \n \
0 # <x> <y> <attribute>...Mesh Regions... \n \
0 # <x> <y> <attribute>...Mesh Regions, area... \n\
#Geo reference \n \
56 \n \
140 \n \
120 \n")
        fid.close()

        tags = {}
        b1 = Dirichlet_boundary(dirichlet_values=num.array([0.0]))
        b2 = Dirichlet_boundary(dirichlet_values=num.array([1.0]))
        b3 = Dirichlet_boundary(dirichlet_values=num.array([2.0]))
        tags["1"] = b1
        tags["2"] = b2
        tags["3"] = b3

        domain = pmesh_to_domain_instance(fileName, Domain)
        os.remove(fileName)
        # print "domain.tagged_elements", domain.tagged_elements
        # # check the quantities
        # print domain.quantities['elevation'].vertex_values
        answer = [[0., 8., 0.],
               [0., 10., 8.]]
        assert num.allclose(domain.quantities['elevation'].vertex_values,
                         answer)

        # print domain.quantities['stage'].vertex_values
        answer = [[0., 12., 10.],
               [0., 10., 12.]]
        assert num.allclose(domain.quantities['stage'].vertex_values,
                         answer)

        # print domain.quantities['friction'].vertex_values
        answer = [[0.01, 0.04, 0.03],
               [0.01, 0.02, 0.04]]
        assert num.allclose(domain.quantities['friction'].vertex_values,
                         answer)

        # print domain.quantities['friction'].vertex_values
        tagged_elements = domain.get_tagged_elements()
        assert num.allclose(tagged_elements['dsg'][0], 0)
        assert num.allclose(tagged_elements['ole nielsen'][0], 1)

        self.assertTrue(domain.boundary[(1, 0)] == '1',
                      "test_tags_to_boundaries  failed. Single boundary wasn't added.")
        self.assertTrue(domain.boundary[(1, 2)] == '2',
                      "test_tags_to_boundaries  failed. Single boundary wasn't added.")
        self.assertTrue(domain.boundary[(0, 1)] == '3',
                      "test_tags_to_boundaries  failed. Single boundary wasn't added.")
        self.assertTrue(domain.boundary[(0, 0)] == 'exterior',
                      "test_tags_to_boundaries  failed. Single boundary wasn't added.")
        # print "domain.boundary",domain.boundary
        self.assertTrue(len(domain.boundary) == 4,
                      "test_pmesh2Domain Too many boundaries")
        # FIXME change to use get_xllcorner
        # print "d.geo_reference.xllcorner",domain.geo_reference.xllcorner
        self.assertTrue(domain.geo_reference.xllcorner == 140.0,
                      "bad geo_referece")
    #************

    def test_pmesh2Domain_instance(self):
        import tempfile

        fd, fileName = tempfile.mkstemp(".tsh")
        os.close(fd)
        fid = open(fileName, "w")
        fid.write("4 3 # <vertex #> <x> <y> [attributes]\n \
0 0.0 0.0 0.0 0.0 0.01 \n \
1 1.0 0.0 10.0 10.0 0.02  \n \
2 0.0 1.0 0.0 10.0 0.03  \n \
3 0.5 0.25 8.0 12.0 0.04  \n \
# Vert att title  \n \
elevation  \n \
stage  \n \
friction  \n \
2 # <triangle #> [<vertex #>] [<neigbouring triangle #>]  \n\
0 0 3 2 -1  -1  1 dsg\n\
1 0 1 3 -1  0 -1   ole nielsen\n\
4 # <segment #> <vertex #>  <vertex #> [boundary tag] \n\
0 1 0 2 \n\
1 0 2 3 \n\
2 2 3 \n\
3 3 1 1 \n\
3 0 # <x> <y> [attributes] ...Mesh Vertices... \n \
0 216.0 -86.0 \n \
1 160.0 -167.0 \n \
2 114.0 -91.0 \n \
3 # <vertex #>  <vertex #> [boundary tag] ...Mesh Segments... \n \
0 0 1 0 \n \
1 1 2 0 \n \
2 2 0 0 \n \
0 # <x> <y> ...Mesh Holes... \n \
0 # <x> <y> <attribute>...Mesh Regions... \n \
0 # <x> <y> <attribute>...Mesh Regions, area... \n\
#Geo reference \n \
56 \n \
140 \n \
120 \n")
        fid.close()

        mesh_instance = import_mesh_from_file(fileName)

        tags = {}
        b1 = Dirichlet_boundary(dirichlet_values=num.array([0.0]))
        b2 = Dirichlet_boundary(dirichlet_values=num.array([1.0]))
        b3 = Dirichlet_boundary(dirichlet_values=num.array([2.0]))
        tags["1"] = b1
        tags["2"] = b2
        tags["3"] = b3

        domain = pmesh_to_domain_instance(mesh_instance, Domain)

        os.remove(fileName)
        # print "domain.tagged_elements", domain.tagged_elements
        # # check the quantities
        # print domain.quantities['elevation'].vertex_values
        answer = [[0., 8., 0.],
               [0., 10., 8.]]
        assert num.allclose(domain.quantities['elevation'].vertex_values,
                         answer)

        # print domain.quantities['stage'].vertex_values
        answer = [[0., 12., 10.],
               [0., 10., 12.]]
        assert num.allclose(domain.quantities['stage'].vertex_values,
                         answer)

        # print domain.quantities['friction'].vertex_values
        answer = [[0.01, 0.04, 0.03],
               [0.01, 0.02, 0.04]]
        assert num.allclose(domain.quantities['friction'].vertex_values,
                         answer)

        # print domain.quantities['friction'].vertex_values
        tagged_elements = domain.get_tagged_elements()
        assert num.allclose(tagged_elements['dsg'][0], 0)
        assert num.allclose(tagged_elements['ole nielsen'][0], 1)

        self.assertTrue(domain.boundary[(1, 0)] == '1',
                      "test_tags_to_boundaries  failed. Single boundary wasn't added.")
        self.assertTrue(domain.boundary[(1, 2)] == '2',
                      "test_tags_to_boundaries  failed. Single boundary wasn't added.")
        self.assertTrue(domain.boundary[(0, 1)] == '3',
                      "test_tags_to_boundaries  failed. Single boundary wasn't added.")
        self.assertTrue(domain.boundary[(0, 0)] == 'exterior',
                      "test_tags_to_boundaries  failed. Single boundary wasn't added.")
        # print "domain.boundary",domain.boundary
        self.assertTrue(len(domain.boundary) == 4,
                      "test_pmesh2Domain Too many boundaries")
        # FIXME change to use get_xllcorner
        # print "d.geo_reference.xllcorner",domain.geo_reference.xllcorner
        self.assertTrue(domain.geo_reference.xllcorner == 140.0,
                      "bad geo_referece")



class Test_pmesh2domain_extra(unittest.TestCase):
    """Cover utility functions in pmesh2domain.py not touched by existing tests."""

    # A minimal two-triangle mesh expressed as a mesh_dict
    _MESH_DICT = {
        'triangles': [[0, 1, 2], [1, 3, 2]],
        'triangle_tags': ['zone_a', 'zone_b'],
        'segments': [[0, 1], [1, 3], [3, 2], [2, 0]],
        'segment_tags': ['left', 'right', 'top', 'bottom'],
        'vertices': [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]],
        'geo_reference': None,
    }

    def test_build_tagged_elements_no_tags(self):
        """build_tagged_elements_dictionary: tri_atts is None path (line 249)."""
        mesh_dict = dict(self._MESH_DICT)
        mesh_dict['triangle_tags'] = None
        result = build_tagged_elements_dictionary(mesh_dict)
        self.assertIn('', result)
        self.assertEqual(len(result['']), 2)

    def test_build_tagged_elements_with_tags(self):
        """build_tagged_elements_dictionary: normal tag path."""
        result = build_tagged_elements_dictionary(self._MESH_DICT)
        self.assertIn('zone_a', result)
        self.assertIn('zone_b', result)
        self.assertEqual(result['zone_a'], [0])
        self.assertEqual(result['zone_b'], [1])

    def test_pmesh_dict_to_tag_dict_old(self):
        """pmesh_dict_to_tag_dict_old: lines 263-290."""
        result = pmesh_dict_to_tag_dict_old(self._MESH_DICT)
        self.assertIsInstance(result, dict)
        # Should have some boundary edges tagged
        self.assertGreater(len(result), 0)
        for key in result:
            vol_id, edge_id = key
            self.assertIn(edge_id, (0, 1, 2))

    def test_calc_sides_old(self):
        """calc_sides_old: lines 327-338."""
        triangles = [[0, 1, 2], [1, 3, 2]]
        sides = calc_sides_old(triangles)
        self.assertIsInstance(sides, dict)
        # Triangle 0: edges (0,1)->2, (1,2)->0, (2,0)->1
        self.assertIn((0, 1), sides)
        self.assertEqual(sides[(0, 1)], 2)   # 3*0+2

    def test_calc_sides_zip(self):
        """calc_sides_zip: lines 347-363."""
        triangles = [[0, 1, 2], [1, 3, 2]]
        sides = calc_sides_zip(triangles)
        self.assertIsInstance(sides, dict)
        self.assertIn((0, 1), sides)

    def test_calc_sides_c(self):
        """calc_sides_c: lines 370-401."""
        triangles = [[0, 1, 2], [1, 3, 2]]
        sides = calc_sides_c(triangles)
        self.assertIsInstance(sides, dict)
        self.assertIn((0, 1), sides)

    def test_pmesh_to_basic_mesh(self):
        """pmesh_to_basic_mesh: lines 34-47 (no prior generate_mesh)."""
        import os
        from anuga.pmesh.mesh import import_mesh_from_file
        from anuga.abstract_2d_finite_volumes.basic_mesh import Basic_mesh
        tsh = os.path.join(os.path.dirname(__file__),
                           '..', '..', 'parallel', 'data', 'small.tsh')
        tsh = os.path.abspath(tsh)
        if not os.path.exists(tsh):
            self.skipTest('small.tsh not found')
        pmesh = import_mesh_from_file(tsh)
        result = pmesh_to_basic_mesh(pmesh)
        self.assertIsInstance(result, Basic_mesh)
        self.assertGreater(result.number_of_triangles, 0)

    def test_pmesh_to_mesh(self):
        """pmesh_to_mesh: lines 69-91."""
        import os
        from anuga.pmesh.mesh import import_mesh_from_file
        from anuga.abstract_2d_finite_volumes.neighbour_mesh import Mesh
        tsh = os.path.join(os.path.dirname(__file__),
                           '..', '..', 'parallel', 'data', 'small.tsh')
        tsh = os.path.abspath(tsh)
        if not os.path.exists(tsh):
            self.skipTest('small.tsh not found')
        pmesh = import_mesh_from_file(tsh)
        result = pmesh_to_mesh(pmesh)
        self.assertIsInstance(result, Mesh)
        self.assertGreater(result.number_of_triangles, 0)


#-------------------------------------------------------------

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_pmesh2domain)
    runner = unittest.TextTestRunner()
    runner.run(suite)

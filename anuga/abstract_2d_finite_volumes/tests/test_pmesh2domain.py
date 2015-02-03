#!/usr/bin/env python
#

import unittest

from anuga.abstract_2d_finite_volumes.pmesh2domain import *

from anuga.shallow_water.shallow_water_domain import Domain
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions \
                        import Dirichlet_boundary

from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.pmesh.mesh import importMeshFromFile

import numpy as num


class Test_pmesh2domain(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_pmesh2Domain(self):
        import os
        import tempfile

        fileName = tempfile.mktemp(".tsh")
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
        import os
        import tempfile

        fileName = tempfile.mktemp(".tsh")
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

        mesh_instance = importMeshFromFile(fileName)

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



#-------------------------------------------------------------

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_pmesh2domain, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)

#!/usr/bin/env python
#

import unittest
import numpy

#from anuga.pyvolution.pmesh2domain import *
from pmesh2domain import *

from anuga.shallow_water import Domain,\
     Reflective_boundary, Dirichlet_boundary,\
     Transmissive_boundary

from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.pmesh.mesh import importMeshFromFile

class Test_pmesh2domain(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_pmesh2Domain(self):
         import os
         import tempfile

         fileName = tempfile.mktemp(".tsh")
         file = open(fileName,"w")
         file.write("4 3 # <vertex #> <x> <y> [attributes]\n \
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
         file.close()

         tags = {}
         b1 =  Dirichlet_boundary(conserved_quantities = numpy.array([0.0]))
         b2 =  Dirichlet_boundary(conserved_quantities = numpy.array([1.0]))
         b3 =  Dirichlet_boundary(conserved_quantities = numpy.array([2.0]))
         tags["1"] = b1
         tags["2"] = b2
         tags["3"] = b3

         domain = pmesh_to_domain_instance(fileName, Domain)
         os.remove(fileName)
         #print "domain.tagged_elements", domain.tagged_elements
         ## check the quantities
         #print domain.quantities['elevation'].vertex_values
         answer = [[0., 8., 0.],
                   [0., 10., 8.]]
         assert numpy.allclose(domain.quantities['elevation'].vertex_values,
                               answer)

         #print domain.quantities['stage'].vertex_values
         answer = [[0., 12., 10.],
                   [0., 10., 12.]]
         assert numpy.allclose(domain.quantities['stage'].vertex_values,
                               answer)

         #print domain.quantities['friction'].vertex_values
         answer = [[0.01, 0.04, 0.03],
                   [0.01, 0.02, 0.04]]
         assert numpy.allclose(domain.quantities['friction'].vertex_values,
                               answer)

         #print domain.quantities['friction'].vertex_values
         assert numpy.allclose(domain.tagged_elements['dsg'][0],0)
         assert numpy.allclose(domain.tagged_elements['ole nielsen'][0],1)

         self.failUnless( domain.boundary[(1, 0)]  == '1',
                          "test_tags_to_boundaries  failed. Single boundary wasn't added.")
         self.failUnless( domain.boundary[(1, 2)]  == '2',
                          "test_tags_to_boundaries  failed. Single boundary wasn't added.")
         self.failUnless( domain.boundary[(0, 1)]  == '3',
                          "test_tags_to_boundaries  failed. Single boundary wasn't added.")
         self.failUnless( domain.boundary[(0, 0)]  == 'exterior',
                          "test_tags_to_boundaries  failed. Single boundary wasn't added.")
         #print "domain.boundary",domain.boundary
         self.failUnless( len(domain.boundary)  == 4,
                          "test_pmesh2Domain Too many boundaries")
         #FIXME change to use get_xllcorner
         #print "d.geo_reference.xllcorner",domain.geo_reference.xllcorner 
         self.failUnless(domain.geo_reference.xllcorner  == 140.0,
                          "bad geo_referece")
    #************
    
    def test_pmesh2Domain_instance(self):
         import os
         import tempfile

         fileName = tempfile.mktemp(".tsh")
         file = open(fileName,"w")
         file.write("4 3 # <vertex #> <x> <y> [attributes]\n \
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
         file.close()

         mesh_instance = importMeshFromFile(fileName)
        
         tags = {}
         b1 =  Dirichlet_boundary(conserved_quantities = numpy.array([0.0]))
         b2 =  Dirichlet_boundary(conserved_quantities = numpy.array([1.0]))
         b3 =  Dirichlet_boundary(conserved_quantities = numpy.array([2.0]))
         tags["1"] = b1
         tags["2"] = b2
         tags["3"] = b3

         domain = pmesh_instance_to_domain_instance(mesh_instance, Domain)

         os.remove(fileName)
         #print "domain.tagged_elements", domain.tagged_elements
         ## check the quantities
         #print domain.quantities['elevation'].vertex_values
         answer = [[0., 8., 0.],
                   [0., 10., 8.]]
         assert numpy.allclose(domain.quantities['elevation'].vertex_values,
                               answer)

         #print domain.quantities['stage'].vertex_values
         answer = [[0., 12., 10.],
                   [0., 10., 12.]]
         assert numpy.allclose(domain.quantities['stage'].vertex_values,
                               answer)

         #print domain.quantities['friction'].vertex_values
         answer = [[0.01, 0.04, 0.03],
                   [0.01, 0.02, 0.04]]
         assert numpy.allclose(domain.quantities['friction'].vertex_values,
                               answer)

         #print domain.quantities['friction'].vertex_values
         assert numpy.allclose(domain.tagged_elements['dsg'][0],0)
         assert numpy.allclose(domain.tagged_elements['ole nielsen'][0],1)

         self.failUnless( domain.boundary[(1, 0)]  == '1',
                          "test_tags_to_boundaries  failed. Single boundary wasn't added.")
         self.failUnless( domain.boundary[(1, 2)]  == '2',
                          "test_tags_to_boundaries  failed. Single boundary wasn't added.")
         self.failUnless( domain.boundary[(0, 1)]  == '3',
                          "test_tags_to_boundaries  failed. Single boundary wasn't added.")
         self.failUnless( domain.boundary[(0, 0)]  == 'exterior',
                          "test_tags_to_boundaries  failed. Single boundary wasn't added.")
         #print "domain.boundary",domain.boundary
         self.failUnless( len(domain.boundary)  == 4,
                          "test_pmesh2Domain Too many boundaries")
         #FIXME change to use get_xllcorner
         #print "d.geo_reference.xllcorner",domain.geo_reference.xllcorner 
         self.failUnless(domain.geo_reference.xllcorner  == 140.0,
                          "bad geo_referece")
         
    #***********
    def old_test_tags_to_boundaries (self):
         meshDict = {}
         p0 = [0.0, 0.0]
         p1 = [1.0, 0.0]
         p2 = [0.0, 1.0]
         p3 = [0.646446609407, 0.353553390593]
         meshDict['vertices'] = [p0,p1,p2,p3]
         meshDict['vertex_attributes'] = [[0.0, 0.0,7.0],[10.0, 0.0,7.0],[0.0, 10.0,7.0],[6.46446609407, 3.53553390593,7.0]]
         meshDict['triangles'] = [[0,3,2],[0,1,3]]
         meshDict['triangle_tags'] = [6.6,6.6]
         meshDict['triangle_neighbors'] = [[-1,-1,1],[-1,0,-1]]
         meshDict['segments'] = [[1,0],[0,2],[2,3],[3,1]]
         meshDict['segment_tags'] = [2,3,1,1]

         domain = Domain.pmesh_dictionary_to_domain(meshDict)

         #domain.set_tag_dict(tag_dict)
         #Boundary tests
         b1 =  Dirichlet_boundary(conserved_quantities = numpy.array([0.0]))
         b2 =  Dirichlet_boundary(conserved_quantities = numpy.array([1.0]))
         b3 =  Dirichlet_boundary(conserved_quantities = numpy.array([1.0]))
         #test adding a boundary
         tags = {}
         tags[1] = b1
         tags[2] = b2
         tags[3] = b3
         domain.set_boundary(tags)
         inverted_id = Volume.instances[0].neighbours[0]
         id = -(inverted_id+1)
         boundary_value = Boundary_value.instances[id]
         boundary_obj = boundary_value.boundary_object

         self.failUnless( boundary_obj  == b1,
                          "test_tags_to_boundaries  failed. Single boundary wasn't added.")

         inverted_id = Volume.instances[0].neighbours[1]
         id = -(inverted_id+1)
         boundary_value = Boundary_value.instances[id]
         boundary_obj = boundary_value.boundary_object
         self.failUnless( boundary_obj  == b3,
                          "test_tags_to_boundaries  failed. Single boundary wasn't added.")

#-------------------------------------------------------------
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_pmesh2domain, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)





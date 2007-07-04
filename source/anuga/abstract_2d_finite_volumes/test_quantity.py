#!/usr/bin/env python

import unittest
from math import sqrt, pi
import tempfile

from quantity import *
from anuga.config import epsilon
from Numeric import allclose, array, ones, Float

from anuga.fit_interpolate.fit import fit_to_mesh
#from anuga.pyvolution.least_squares import fit_to_mesh         
from domain import Domain
from anuga.geospatial_data.geospatial_data import Geospatial_data
from anuga.coordinate_transforms.geo_reference import Geo_reference

#Aux for fit_interpolate.fit example
def linear_function(point):
    point = array(point)
    return point[:,0]+point[:,1]


class Test_Quantity(unittest.TestCase):
    def setUp(self):
        from domain import Domain

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]

        #bac, bce, ecf, dbe
        elements = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4] ]

        self.mesh1 = Domain(points[:3], [elements[0]])
        self.mesh1.check_integrity()

        self.mesh4 = Domain(points, elements)
        self.mesh4.check_integrity()

        # UTM round Onslow
        a = [240000, 7620000]
        b = [240000, 7680000]
        c = [300000, 7620000]

        points = [a, b, c]
        elements = [[0,2,1]]
        
        self.mesh_onslow = Domain(points, elements)
        self.mesh_onslow.check_integrity()
        
    def tearDown(self):
        pass
        #print "  Tearing down"


    def test_creation(self):

        quantity = Quantity(self.mesh1, [[1,2,3]])
        assert allclose(quantity.vertex_values, [[1.,2.,3.]])

        try:
            quantity = Quantity()
        except:
            pass
        else:
            raise 'Should have raised empty quantity exception'


        try:
            quantity = Quantity([1,2,3])
        except AssertionError:
            pass
        except:
            raise 'Should have raised "mising mesh object" error'


    def test_creation_zeros(self):

        quantity = Quantity(self.mesh1)
        assert allclose(quantity.vertex_values, [[0.,0.,0.]])


        quantity = Quantity(self.mesh4)
        assert allclose(quantity.vertex_values, [[0.,0.,0.], [0.,0.,0.],
                                                 [0.,0.,0.], [0.,0.,0.]])


    def test_interpolation(self):
        quantity = Quantity(self.mesh1, [[1,2,3]])
        assert allclose(quantity.centroid_values, [2.0]) #Centroid

        assert allclose(quantity.edge_values, [[2.5, 2.0, 1.5]])


    def test_interpolation2(self):
        quantity = Conserved_quantity(self.mesh4,
                            [[1,2,3], [5,5,5], [0,0,9], [-6, 3, 3]])
        assert allclose(quantity.centroid_values, [2., 5., 3., 0.]) #Centroid


        quantity.extrapolate_second_order()

        #print quantity.vertex_values
        #assert allclose(quantity.vertex_values, [[2., 2., 2.],
        #                                         [3.+2./3, 6.+2./3, 4.+2./3],
        #                                         [7.5, 0.5, 1.],
        #                                         [-5, -2.5, 7.5]])

        assert allclose(quantity.vertex_values[1,:],[3.+2./3, 6.+2./3, 4.+2./3])
        #FIXME: Work out the others


        #print quantity.edge_values
        assert allclose(quantity.edge_values, [[2.5, 2.0, 1.5],
                                               [5., 5., 5.],
                                               [4.5, 4.5, 0.],
                                               [3.0, -1.5, -1.5]])

    def test_get_maximum_1(self):
        quantity = Conserved_quantity(self.mesh4,
                                      [[1,2,3], [5,5,5], [0,0,9], [-6, 3, 3]])
        assert allclose(quantity.centroid_values, [2., 5., 3., 0.]) #Centroids

        v = quantity.get_maximum_value()
        assert v == 5

        i = quantity.get_maximum_index()
        assert i == 1
        
        x,y = quantity.get_maximum_location()
        xref, yref = 4.0/3, 4.0/3
        assert x == xref
        assert y == yref

        v = quantity.get_values(interpolation_points = [[x,y]])
        assert allclose(v, 5)

    def test_get_maximum_2(self):

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [[1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices)

        quantity = Quantity(domain)
        quantity.set_values(lambda x, y: x+2*y) #2 4 4 6
        
        v = quantity.get_maximum_value()
        assert v == 6

        i = quantity.get_maximum_index()
        assert i == 3
        
        x,y = quantity.get_maximum_location()
        xref, yref = 2.0/3, 8.0/3
        assert x == xref
        assert y == yref

        v = quantity.get_values(interpolation_points = [[x,y]])
        assert allclose(v, 6)


        #Multiple locations for maximum -
        #Test that the algorithm picks the first occurrence        
        v = quantity.get_maximum_value(indices=[0,1,2])
        assert allclose(v, 4)

        i = quantity.get_maximum_index(indices=[0,1,2])
        assert i == 1
        
        x,y = quantity.get_maximum_location(indices=[0,1,2])
        xref, yref = 4.0/3, 4.0/3
        assert x == xref
        assert y == yref

        v = quantity.get_values(interpolation_points = [[x,y]])
        assert allclose(v, 4)        

        # More test of indices......
        v = quantity.get_maximum_value(indices=[2,3])
        assert allclose(v, 6)

        i = quantity.get_maximum_index(indices=[2,3])
        assert i == 3
        
        x,y = quantity.get_maximum_location(indices=[2,3])
        xref, yref = 2.0/3, 8.0/3
        assert x == xref
        assert y == yref

        v = quantity.get_values(interpolation_points = [[x,y]])
        assert allclose(v, 6)        

        

    def test_boundary_allocation(self):
        quantity = Conserved_quantity(self.mesh4,
                            [[1,2,3], [5,5,5], [0,0,9], [-6, 3, 3]])

        assert quantity.boundary_values.shape[0] == len(self.mesh4.boundary)


    def test_set_values(self):
        quantity = Quantity(self.mesh4)


        quantity.set_values([[1,2,3], [5,5,5], [0,0,9], [-6, 3, 3]],
                            location = 'vertices')
        assert allclose(quantity.vertex_values,
                        [[1,2,3], [5,5,5], [0,0,9], [-6, 3, 3]])
        assert allclose(quantity.centroid_values, [2., 5., 3., 0.]) #Centroid
        assert allclose(quantity.edge_values, [[2.5, 2.0, 1.5],
                                               [5., 5., 5.],
                                               [4.5, 4.5, 0.],
                                               [3.0, -1.5, -1.5]])


        #Test default
        quantity.set_values([[1,2,3], [5,5,5], [0,0,9], [-6, 3, 3]])
        assert allclose(quantity.vertex_values,
                        [[1,2,3], [5,5,5], [0,0,9], [-6, 3, 3]])
        assert allclose(quantity.centroid_values, [2., 5., 3., 0.]) #Centroid
        assert allclose(quantity.edge_values, [[2.5, 2.0, 1.5],
                                               [5., 5., 5.],
                                               [4.5, 4.5, 0.],
                                               [3.0, -1.5, -1.5]])

        #Test centroids
        quantity.set_values([1,2,3,4], location = 'centroids')
        assert allclose(quantity.centroid_values, [1., 2., 3., 4.]) #Centroid

        #Test edges
        quantity.set_values([[1,2,3], [5,5,5], [0,0,9], [-6, 3, 3]],
                            location = 'edges')
        assert allclose(quantity.edge_values,
                        [[1,2,3], [5,5,5], [0,0,9], [-6, 3, 3]])

        #Test exceptions
        try:
            quantity.set_values([[1,2,3], [5,5,5], [0,0,9], [-6, 3, 3]],
                                location = 'bas kamel tuba')
        except:
            pass


        try:
            quantity.set_values([[1,2,3], [0,0,9]])
        except AssertionError:
            pass
        except:
            raise 'should have raised Assertionerror'



    def test_set_values_const(self):
        quantity = Quantity(self.mesh4)

        quantity.set_values(1.0, location = 'vertices')
        assert allclose(quantity.vertex_values,
                        [[1,1,1], [1,1,1], [1,1,1], [1, 1, 1]])

        assert allclose(quantity.centroid_values, [1, 1, 1, 1]) #Centroid
        assert allclose(quantity.edge_values, [[1, 1, 1],
                                               [1, 1, 1],
                                               [1, 1, 1],
                                               [1, 1, 1]])


        quantity.set_values(2.0, location = 'centroids')
        assert allclose(quantity.centroid_values, [2, 2, 2, 2])

        quantity.set_values(3.0, location = 'edges')
        assert allclose(quantity.edge_values, [[3, 3, 3],
                                               [3, 3, 3],
                                               [3, 3, 3],
                                               [3, 3, 3]])


    def test_set_values_func(self):
        quantity = Quantity(self.mesh4)

        def f(x, y):
            return x+y

        quantity.set_values(f, location = 'vertices')
        #print "quantity.vertex_values",quantity.vertex_values
        assert allclose(quantity.vertex_values,
                        [[2,0,2], [2,2,4], [4,2,4], [4,2,4]])
        assert allclose(quantity.centroid_values,
                        [4.0/3, 8.0/3, 10.0/3, 10.0/3])
        assert allclose(quantity.edge_values,
                        [[1,2,1], [3,3,2], [3,4,3], [3,4,3]])


        quantity.set_values(f, location = 'centroids')
        assert allclose(quantity.centroid_values,
                        [4.0/3, 8.0/3, 10.0/3, 10.0/3])


    def test_integral(self):
        quantity = Quantity(self.mesh4)

        #Try constants first
        const = 5
        quantity.set_values(const, location = 'vertices')
        #print 'Q', quantity.get_integral()

        assert allclose(quantity.get_integral(), self.mesh4.get_area() * const)

        #Try with a linear function
        def f(x, y):
            return x+y

        quantity.set_values(f, location = 'vertices')


        ref_integral = (4.0/3 + 8.0/3 + 10.0/3 + 10.0/3) * 2

        assert allclose (quantity.get_integral(), ref_integral)



    def test_set_vertex_values(self):
        quantity = Quantity(self.mesh4)
        quantity.set_vertex_values([0,1,2,3,4,5])

        assert allclose(quantity.vertex_values,
                        [[1,0,2], [1,2,4], [4,2,5], [3,1,4]])
        assert allclose(quantity.centroid_values,
                        [1., 7./3, 11./3, 8./3]) #Centroid
        assert allclose(quantity.edge_values, [[1., 1.5, 0.5],
                                               [3., 2.5, 1.5],
                                               [3.5, 4.5, 3.],
                                               [2.5, 3.5, 2]])


    def test_set_vertex_values_subset(self):
        quantity = Quantity(self.mesh4)
        quantity.set_vertex_values([0,1,2,3,4,5])
        quantity.set_vertex_values([0,20,30,50], indices = [0,2,3,5])

        assert allclose(quantity.vertex_values,
                        [[1,0,20], [1,20,4], [4,20,50], [30,1,4]])


    def test_set_vertex_values_using_general_interface(self):
        quantity = Quantity(self.mesh4)


        quantity.set_values([0,1,2,3,4,5])


        assert allclose(quantity.vertex_values,
                        [[1,0,2], [1,2,4], [4,2,5], [3,1,4]])

        #Centroid
        assert allclose(quantity.centroid_values, [1., 7./3, 11./3, 8./3])

        assert allclose(quantity.edge_values, [[1., 1.5, 0.5],
                                               [3., 2.5, 1.5],
                                               [3.5, 4.5, 3.],
                                               [2.5, 3.5, 2]])



    def test_set_vertex_values_using_general_interface_with_subset(self):
        """test_set_vertex_values_using_general_interface_with_subset(self):
        Test that indices and polygon works
        """
        
        quantity = Quantity(self.mesh4)


        quantity.set_values([0,2,3,5], indices=[0,2,3,5])
        assert allclose(quantity.vertex_values,
                        [[0,0,2], [0,2,0], [0,2,5], [3,0,0]])


        # Constant
        quantity.set_values(0.0)
        quantity.set_values(3.14, indices=[0,2], location='vertices')

        # Indices refer to triangle numbers here - not vertices (why?)
        assert allclose(quantity.vertex_values,
                        [[3.14,3.14,3.14], [0,0,0],
                         [3.14,3.14,3.14], [0,0,0]])        
        


        # Now try with polygon (pick points where y>2)
        polygon = [[0,2.1], [4,2.1], [4,7], [0,7]]
        quantity.set_values(0.0)
        quantity.set_values(3.14, polygon=polygon)
        
        assert allclose(quantity.vertex_values,
                        [[0,0,0], [0,0,0], [0,0,0],
                         [3.14,3.14,3.14]])                


        # Another polygon (pick triangle 1 and 2 (rightmost triangles)
        polygon = [[2.1, 0.0], [3.5,0.1], [2,2.2], [0.2,2]]
        quantity.set_values(0.0)
        quantity.set_values(3.14, polygon=polygon)

        assert allclose(quantity.vertex_values,
                        [[0,0,0],
                         [3.14,3.14,3.14],
                         [3.14,3.14,3.14],                         
                         [0,0,0]])                



        # Test input checking
        try:
            quantity.set_values(3.14, polygon=polygon, indices = [0,2])
        except:
            pass
        else:
            msg = 'Should have caught this'
            raise msg

            

    def test_set_values_using_fit(self):


        quantity = Quantity(self.mesh4)

        #Get (enough) datapoints
        data_points = [[ 0.66666667, 0.66666667],
                       [ 1.33333333, 1.33333333],
                       [ 2.66666667, 0.66666667],
                       [ 0.66666667, 2.66666667],
                       [ 0.0, 1.0],
                       [ 0.0, 3.0],
                       [ 1.0, 0.0],
                       [ 1.0, 1.0],
                       [ 1.0, 2.0],
                       [ 1.0, 3.0],
                       [ 2.0, 1.0],
                       [ 3.0, 0.0],
                       [ 3.0, 1.0]]

        z = linear_function(data_points)

        #Use built-in fit_interpolate.fit
        quantity.set_values( Geospatial_data(data_points, z), alpha = 0 )
        #quantity.set_values(points = data_points, values = z, alpha = 0)


        answer = linear_function(quantity.domain.get_vertex_coordinates())
        #print quantity.vertex_values, answer
        assert allclose(quantity.vertex_values.flat, answer)


        #Now try by setting the same values directly
        vertex_attributes = fit_to_mesh(quantity.domain.get_nodes(),
                                        quantity.domain.triangles, #FIXME
                                        data_points,
                                        z,
                                        alpha = 0,
                                        verbose=False)

        #print vertex_attributes
        quantity.set_values(vertex_attributes)
        assert allclose(quantity.vertex_values.flat, answer)





    def test_test_set_values_using_fit_w_geo(self):


        #Mesh
        vertex_coordinates = [[0.76, 0.76],
                              [0.76, 5.76],
                              [5.76, 0.76]]
        triangles = [[0,2,1]]

        mesh_georef = Geo_reference(56,-0.76,-0.76)
        mesh1 = Domain(vertex_coordinates, triangles,
                       geo_reference = mesh_georef)
        mesh1.check_integrity()

        #Quantity
        quantity = Quantity(mesh1)

        #Data
        data_points = [[ 201.0, 401.0],
                       [ 201.0, 403.0],
                       [ 203.0, 401.0]]

        z = [2, 4, 4]

        data_georef = Geo_reference(56,-200,-400)


        #Reference
        ref = fit_to_mesh(vertex_coordinates, triangles, data_points, z,
                          data_origin = data_georef.get_origin(),
                          mesh_origin = mesh_georef.get_origin(),
                          alpha = 0)

        assert allclose( ref, [0,5,5] )


        #Test set_values

        quantity.set_values( Geospatial_data(data_points, z, data_georef), alpha = 0 )

        #quantity.set_values(points = data_points,
        #                    values = z,
        #                    data_georef = data_georef,
        #                    alpha = 0)


        #quantity.set_values(points = data_points,
        #                    values = z,
        #                    data_georef = data_georef,
        #                    alpha = 0)
        assert allclose(quantity.vertex_values.flat, ref)



        #Test set_values using geospatial data object
        quantity.vertex_values[:] = 0.0

        geo = Geospatial_data(data_points, z, data_georef)


        quantity.set_values(geospatial_data = geo, alpha = 0)
        assert allclose(quantity.vertex_values.flat, ref)



    def test_set_values_from_file1(self):
        quantity = Quantity(self.mesh4)

        #Get (enough) datapoints
        data_points = [[ 0.66666667, 0.66666667],
                       [ 1.33333333, 1.33333333],
                       [ 2.66666667, 0.66666667],
                       [ 0.66666667, 2.66666667],
                       [ 0.0, 1.0],
                       [ 0.0, 3.0],
                       [ 1.0, 0.0],
                       [ 1.0, 1.0],
                       [ 1.0, 2.0],
                       [ 1.0, 3.0],
                       [ 2.0, 1.0],
                       [ 3.0, 0.0],
                       [ 3.0, 1.0]]

        data_geo_spatial = Geospatial_data(data_points,
                         geo_reference = Geo_reference(56, 0, 0))
        data_points_absolute = data_geo_spatial.get_data_points(absolute=True)
        attributes = linear_function(data_points_absolute)
        att = 'spam_and_eggs'
        
        #Create .txt file
        ptsfile = tempfile.mktemp(".txt")
        file = open(ptsfile,"w")
        file.write(" x,y," + att + " \n")
        for data_point, attribute in map(None, data_points_absolute
                                         ,attributes):
            row = str(data_point[0]) + ',' + str(data_point[1]) \
                  + ',' + str(attribute)
            file.write(row + "\n")
        file.close()


        #Check that values can be set from file
        quantity.set_values(filename = ptsfile,
                            attribute_name = att, alpha = 0)
        answer = linear_function(quantity.domain.get_vertex_coordinates())

        #print quantity.vertex_values.flat
        #print answer


        assert allclose(quantity.vertex_values.flat, answer)


        #Check that values can be set from file using default attribute
        quantity.set_values(filename = ptsfile, alpha = 0)
        assert allclose(quantity.vertex_values.flat, answer)

        #Cleanup
        import os
        os.remove(ptsfile)

    def Cache_cache_test_set_values_from_file(self):
        quantity = Quantity(self.mesh4)

        #Get (enough) datapoints
        data_points = [[ 0.66666667, 0.66666667],
                       [ 1.33333333, 1.33333333],
                       [ 2.66666667, 0.66666667],
                       [ 0.66666667, 2.66666667],
                       [ 0.0, 1.0],
                       [ 0.0, 3.0],
                       [ 1.0, 0.0],
                       [ 1.0, 1.0],
                       [ 1.0, 2.0],
                       [ 1.0, 3.0],
                       [ 2.0, 1.0],
                       [ 3.0, 0.0],
                       [ 3.0, 1.0]]

        data_geo_spatial = Geospatial_data(data_points,
                         geo_reference = Geo_reference(56, 0, 0))
        data_points_absolute = data_geo_spatial.get_data_points(absolute=True)
        attributes = linear_function(data_points_absolute)
        att = 'spam_and_eggs'
        
        #Create .txt file
        ptsfile = tempfile.mktemp(".txt")
        file = open(ptsfile,"w")
        file.write(" x,y," + att + " \n")
        for data_point, attribute in map(None, data_points_absolute
                                         ,attributes):
            row = str(data_point[0]) + ',' + str(data_point[1]) \
                  + ',' + str(attribute)
            file.write(row + "\n")
        file.close()


        #Check that values can be set from file
        quantity.set_values(filename = ptsfile,
                            attribute_name = att, alpha = 0, use_cache=True,
                            verbose=True)
        answer = linear_function(quantity.domain.get_vertex_coordinates())

        #print quantity.vertex_values.flat
        #print answer


        assert allclose(quantity.vertex_values.flat, answer)


        #Check that values can be set from file using default attribute
        quantity.set_values(filename = ptsfile, alpha = 0)
        assert allclose(quantity.vertex_values.flat, answer)

        #checking cache
        #quantity.set_values(filename = ptsfile,
         #                   attribute_name = att, alpha = 0, use_cache=True,
          #                  verbose=True)
        #Cleanup
        import os
        os.remove(ptsfile)

    def test_set_values_from_lat_long(self):
        quantity = Quantity(self.mesh_onslow)

        #Get (enough) datapoints
        data_points = [[-21.5, 114.5],[-21.4, 114.6],[-21.45,114.65],
                       [-21.35, 114.65],[-21.45, 114.55],[-21.45,114.6]]

        data_geo_spatial = Geospatial_data(data_points,
                                           points_are_lats_longs=True)
        points_UTM = data_geo_spatial.get_data_points(absolute=True)
        attributes = linear_function(points_UTM)
        att = 'elevation'
        
        #Create .txt file
        txt_file = tempfile.mktemp(".txt")
        file = open(txt_file,"w")
        file.write(" lat,long," + att + " \n")
        for data_point, attribute in map(None, data_points, attributes):
            row = str(data_point[0]) + ',' + str(data_point[1]) \
                  + ',' + str(attribute)
            #print "row", row 
            file.write(row + "\n")
        file.close()


        #Check that values can be set from file
        quantity.set_values(filename = txt_file,
                            attribute_name = att, alpha = 0)
        answer = linear_function(quantity.domain.get_vertex_coordinates())

        #print "quantity.vertex_values.flat", quantity.vertex_values.flat
        #print "answer",answer

        assert allclose(quantity.vertex_values.flat, answer)


        #Check that values can be set from file using default attribute
        quantity.set_values(filename = txt_file, alpha = 0)
        assert allclose(quantity.vertex_values.flat, answer)

        #Cleanup
        import os
        os.remove(txt_file)
         
    def test_set_values_from_lat_long(self):
        quantity = Quantity(self.mesh_onslow)

        #Get (enough) datapoints
        data_points = [[-21.5, 114.5],[-21.4, 114.6],[-21.45,114.65],
                       [-21.35, 114.65],[-21.45, 114.55],[-21.45,114.6]]

        data_geo_spatial = Geospatial_data(data_points,
                                           points_are_lats_longs=True)
        points_UTM = data_geo_spatial.get_data_points(absolute=True)
        attributes = linear_function(points_UTM)
        att = 'elevation'
        
        #Create .txt file
        txt_file = tempfile.mktemp(".txt")
        file = open(txt_file,"w")
        file.write(" lat,long," + att + " \n")
        for data_point, attribute in map(None, data_points, attributes):
            row = str(data_point[0]) + ',' + str(data_point[1]) \
                  + ',' + str(attribute)
            #print "row", row 
            file.write(row + "\n")
        file.close()


        #Check that values can be set from file
        quantity.set_values(filename = txt_file,
                            attribute_name = att, alpha = 0)
        answer = linear_function(quantity.domain.get_vertex_coordinates())

        #print "quantity.vertex_values.flat", quantity.vertex_values.flat
        #print "answer",answer

        assert allclose(quantity.vertex_values.flat, answer)


        #Check that values can be set from file using default attribute
        quantity.set_values(filename = txt_file, alpha = 0)
        assert allclose(quantity.vertex_values.flat, answer)

        #Cleanup
        import os
        os.remove(txt_file)
        
    def test_set_values_from_UTM_pts(self):
        quantity = Quantity(self.mesh_onslow)

        #Get (enough) datapoints
        data_points = [[-21.5, 114.5],[-21.4, 114.6],[-21.45,114.65],
                       [-21.35, 114.65],[-21.45, 114.55],[-21.45,114.6]]

        data_geo_spatial = Geospatial_data(data_points,
                                           points_are_lats_longs=True)
        points_UTM = data_geo_spatial.get_data_points(absolute=True)
        attributes = linear_function(points_UTM)
        att = 'elevation'
        
        #Create .txt file
        txt_file = tempfile.mktemp(".txt")
        file = open(txt_file,"w")
        file.write(" x,y," + att + " \n")
        for data_point, attribute in map(None, points_UTM, attributes):
            row = str(data_point[0]) + ',' + str(data_point[1]) \
                  + ',' + str(attribute)
            #print "row", row 
            file.write(row + "\n")
        file.close()


        pts_file = tempfile.mktemp(".pts")        
        convert = Geospatial_data(txt_file)
        convert.export_points_file(pts_file)

        #Check that values can be set from file
        quantity.set_values_from_file(pts_file, att, 0,
                                      'vertices', None)
        answer = linear_function(quantity.domain.get_vertex_coordinates())
        #print "quantity.vertex_values.flat", quantity.vertex_values.flat
        #print "answer",answer
        assert allclose(quantity.vertex_values.flat, answer)

        #Check that values can be set from file
        quantity.set_values(filename = pts_file,
                            attribute_name = att, alpha = 0)
        answer = linear_function(quantity.domain.get_vertex_coordinates())
        #print "quantity.vertex_values.flat", quantity.vertex_values.flat
        #print "answer",answer
        assert allclose(quantity.vertex_values.flat, answer)


        #Check that values can be set from file using default attribute
        quantity.set_values(filename = txt_file, alpha = 0)
        assert allclose(quantity.vertex_values.flat, answer)

        #Cleanup
        import os
        os.remove(txt_file)
        os.remove(pts_file)
        
    def verbose_test_set_values_from_UTM_pts(self):
        quantity = Quantity(self.mesh_onslow)

        #Get (enough) datapoints
        data_points = [[-21.5, 114.5],[-21.4, 114.6],[-21.45,114.65],
                       [-21.35, 114.65],[-21.45, 114.55],[-21.45,114.6],
                       [-21.5, 114.5],[-21.4, 114.6],[-21.45,114.65],
                       [-21.35, 114.65],[-21.45, 114.55],[-21.45,114.6],
                       [-21.5, 114.5],[-21.4, 114.6],[-21.45,114.65],
                       [-21.35, 114.65],[-21.45, 114.55],[-21.45,114.6],
                       [-21.5, 114.5],[-21.4, 114.6],[-21.45,114.65],
                       [-21.35, 114.65],[-21.45, 114.55],[-21.45,114.6],
                       [-21.5, 114.5],[-21.4, 114.6],[-21.45,114.65],
                       [-21.35, 114.65],[-21.45, 114.55],[-21.45,114.6],
                       [-21.5, 114.5],[-21.4, 114.6],[-21.45,114.65],
                       [-21.35, 114.65],[-21.45, 114.55],[-21.45,114.6],
                       [-21.5, 114.5],[-21.4, 114.6],[-21.45,114.65],
                       [-21.35, 114.65],[-21.45, 114.55],[-21.45,114.6],
                       [-21.35, 114.65],[-21.45, 114.55],[-21.45,114.6],
                       [-21.5, 114.5],[-21.4, 114.6],[-21.45,114.65],
                       [-21.35, 114.65],[-21.45, 114.55],[-21.45,114.6],
                       [-21.5, 114.5],[-21.4, 114.6],[-21.45,114.65],
                       [-21.35, 114.65],[-21.45, 114.55],[-21.45,114.6],
                       [-21.5, 114.5],[-21.4, 114.6],[-21.45,114.65],
                       [-21.35, 114.65],[-21.45, 114.55],[-21.45,114.6],
                       [-21.5, 114.5],[-21.4, 114.6],[-21.45,114.65],
                       [-21.35, 114.65],[-21.45, 114.55],[-21.45,114.6],
                       [-21.5, 114.5],[-21.4, 114.6],[-21.45,114.65],
                       [-21.35, 114.65],[-21.45, 114.55],[-21.45,114.6],
                       [-21.5, 114.5],[-21.4, 114.6],[-21.45,114.65],
                       [-21.35, 114.65],[-21.45, 114.55],[-21.45,114.6],
                       ]

        data_geo_spatial = Geospatial_data(data_points,
                                           points_are_lats_longs=True)
        points_UTM = data_geo_spatial.get_data_points(absolute=True)
        attributes = linear_function(points_UTM)
        att = 'elevation'
        
        #Create .txt file
        txt_file = tempfile.mktemp(".txt")
        file = open(txt_file,"w")
        file.write(" x,y," + att + " \n")
        for data_point, attribute in map(None, points_UTM, attributes):
            row = str(data_point[0]) + ',' + str(data_point[1]) \
                  + ',' + str(attribute)
            #print "row", row 
            file.write(row + "\n")
        file.close()


        pts_file = tempfile.mktemp(".pts")        
        convert = Geospatial_data(txt_file)
        convert.export_points_file(pts_file)

        #Check that values can be set from file
        quantity.set_values_from_file(pts_file, att, 0,
                                      'vertices', None, verbose = True,
                                      max_read_lines=2)
        answer = linear_function(quantity.domain.get_vertex_coordinates())
        #print "quantity.vertex_values.flat", quantity.vertex_values.flat
        #print "answer",answer
        assert allclose(quantity.vertex_values.flat, answer)

        #Check that values can be set from file
        quantity.set_values(filename = pts_file,
                            attribute_name = att, alpha = 0)
        answer = linear_function(quantity.domain.get_vertex_coordinates())
        #print "quantity.vertex_values.flat", quantity.vertex_values.flat
        #print "answer",answer
        assert allclose(quantity.vertex_values.flat, answer)


        #Check that values can be set from file using default attribute
        quantity.set_values(filename = txt_file, alpha = 0)
        assert allclose(quantity.vertex_values.flat, answer)

        #Cleanup
        import os
        os.remove(txt_file)
        os.remove(pts_file)
        
    def test_set_values_from_file_with_georef1(self):

        #Mesh in zone 56 (absolute coords)

        x0 = 314036.58727982
        y0 = 6224951.2960092

        a = [x0+0.0, y0+0.0]
        b = [x0+0.0, y0+2.0]
        c = [x0+2.0, y0+0.0]
        d = [x0+0.0, y0+4.0]
        e = [x0+2.0, y0+2.0]
        f = [x0+4.0, y0+0.0]

        points = [a, b, c, d, e, f]

        #bac, bce, ecf, dbe
        elements = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4] ]

        #absolute going in ..
        mesh4 = Domain(points, elements,
                       geo_reference = Geo_reference(56, 0, 0))
        mesh4.check_integrity()
        quantity = Quantity(mesh4)

        #Get (enough) datapoints (relative to georef)
        data_points_rel = [[ 0.66666667, 0.66666667],
                       [ 1.33333333, 1.33333333],
                       [ 2.66666667, 0.66666667],
                       [ 0.66666667, 2.66666667],
                       [ 0.0, 1.0],
                       [ 0.0, 3.0],
                       [ 1.0, 0.0],
                       [ 1.0, 1.0],
                       [ 1.0, 2.0],
                       [ 1.0, 3.0],
                       [ 2.0, 1.0],
                       [ 3.0, 0.0],
                       [ 3.0, 1.0]]

        data_geo_spatial = Geospatial_data(data_points_rel,
                         geo_reference = Geo_reference(56, x0, y0))
        data_points_absolute = data_geo_spatial.get_data_points(absolute=True)
        attributes = linear_function(data_points_absolute)
        att = 'spam_and_eggs'
        
        #Create .txt file
        ptsfile = tempfile.mktemp(".txt")
        file = open(ptsfile,"w")
        file.write(" x,y," + att + " \n")
        for data_point, attribute in map(None, data_points_absolute
                                         ,attributes):
            row = str(data_point[0]) + ',' + str(data_point[1]) \
                  + ',' + str(attribute)
            file.write(row + "\n")
        file.close()

        #file = open(ptsfile, 'r')
        #lines = file.readlines()
        #file.close()
     

        #Check that values can be set from file
        quantity.set_values(filename = ptsfile,
                            attribute_name = att, alpha = 0)
        answer = linear_function(quantity.domain.get_vertex_coordinates())

        assert allclose(quantity.vertex_values.flat, answer)


        #Check that values can be set from file using default attribute
        quantity.set_values(filename = ptsfile, alpha = 0)
        assert allclose(quantity.vertex_values.flat, answer)

        #Cleanup
        import os
        os.remove(ptsfile)


    def test_set_values_from_file_with_georef2(self):

        #Mesh in zone 56 (relative coords)

        x0 = 314036.58727982
        y0 = 6224951.2960092

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]

        #bac, bce, ecf, dbe
        elements = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4] ]

        mesh4 = Domain(points, elements,
                       geo_reference = Geo_reference(56, x0, y0))
        mesh4.check_integrity()
        quantity = Quantity(mesh4)

        #Get (enough) datapoints
        data_points = [[ x0+0.66666667, y0+0.66666667],
                       [ x0+1.33333333, y0+1.33333333],
                       [ x0+2.66666667, y0+0.66666667],
                       [ x0+0.66666667, y0+2.66666667],
                       [ x0+0.0, y0+1.0],
                       [ x0+0.0, y0+3.0],
                       [ x0+1.0, y0+0.0],
                       [ x0+1.0, y0+1.0],
                       [ x0+1.0, y0+2.0],
                       [ x0+1.0, y0+3.0],
                       [ x0+2.0, y0+1.0],
                       [ x0+3.0, y0+0.0],
                       [ x0+3.0, y0+1.0]]


        data_geo_spatial = Geospatial_data(data_points,
                         geo_reference = Geo_reference(56, 0, 0))
        data_points_absolute = data_geo_spatial.get_data_points(absolute=True)
        attributes = linear_function(data_points_absolute)
        att = 'spam_and_eggs'
        
        #Create .txt file
        ptsfile = tempfile.mktemp(".txt")
        file = open(ptsfile,"w")
        file.write(" x,y," + att + " \n")
        for data_point, attribute in map(None, data_points_absolute
                                         ,attributes):
            row = str(data_point[0]) + ',' + str(data_point[1]) \
                  + ',' + str(attribute)
            file.write(row + "\n")
        file.close()


        #Check that values can be set from file
        quantity.set_values(filename = ptsfile,
                            attribute_name = att, alpha = 0)
        answer = linear_function(quantity.domain. \
                                 get_vertex_coordinates(absolute=True))


        assert allclose(quantity.vertex_values.flat, answer)


        #Check that values can be set from file using default attribute
        quantity.set_values(filename = ptsfile, alpha = 0)
        assert allclose(quantity.vertex_values.flat, answer)

        #Cleanup
        import os
        os.remove(ptsfile)




    def test_set_values_from_quantity(self):

        quantity1 = Quantity(self.mesh4)
        quantity1.set_vertex_values([0,1,2,3,4,5])

        assert allclose(quantity1.vertex_values,
                        [[1,0,2], [1,2,4], [4,2,5], [3,1,4]])


        quantity2 = Quantity(self.mesh4)
        quantity2.set_values(quantity = quantity1)
        assert allclose(quantity2.vertex_values,
                        [[1,0,2], [1,2,4], [4,2,5], [3,1,4]])

        quantity2.set_values(quantity = 2*quantity1)
        assert allclose(quantity2.vertex_values,
                        [[2,0,4], [2,4,8], [8,4,10], [6,2,8]])

        quantity2.set_values(quantity = 2*quantity1 + 3)
        assert allclose(quantity2.vertex_values,
                        [[5,3,7], [5,7,11], [11,7,13], [9,5,11]])


        #Check detection of quantity as first orgument
        quantity2.set_values(2*quantity1 + 3)
        assert allclose(quantity2.vertex_values,
                        [[5,3,7], [5,7,11], [11,7,13], [9,5,11]])





    def test_overloading(self):

        quantity1 = Quantity(self.mesh4)
        quantity1.set_vertex_values([0,1,2,3,4,5])

        assert allclose(quantity1.vertex_values,
                        [[1,0,2], [1,2,4], [4,2,5], [3,1,4]])


        quantity2 = Quantity(self.mesh4)
        quantity2.set_values([[1,2,3], [5,5,5], [0,0,9], [-6, 3, 3]],
                             location = 'vertices')



        quantity3 = Quantity(self.mesh4)
        quantity3.set_values([[2,2,2], [7,8,9], [7,6,3], [3, 8, -8]],
                             location = 'vertices')


        #Negation
        Q = -quantity1
        assert allclose(Q.vertex_values, -quantity1.vertex_values)
        assert allclose(Q.centroid_values, -quantity1.centroid_values)
        assert allclose(Q.edge_values, -quantity1.edge_values)

        #Addition
        Q = quantity1 + 7
        assert allclose(Q.vertex_values, quantity1.vertex_values + 7)
        assert allclose(Q.centroid_values, quantity1.centroid_values + 7)
        assert allclose(Q.edge_values, quantity1.edge_values + 7)

        Q = 7 + quantity1
        assert allclose(Q.vertex_values, quantity1.vertex_values + 7)
        assert allclose(Q.centroid_values, quantity1.centroid_values + 7)
        assert allclose(Q.edge_values, quantity1.edge_values + 7)

        Q = quantity1 + quantity2
        assert allclose(Q.vertex_values,
                        quantity1.vertex_values + quantity2.vertex_values)
        assert allclose(Q.centroid_values,
                        quantity1.centroid_values + quantity2.centroid_values)
        assert allclose(Q.edge_values,
                        quantity1.edge_values + quantity2.edge_values)


        Q = quantity1 + quantity2 - 3
        assert allclose(Q.vertex_values,
                        quantity1.vertex_values + quantity2.vertex_values - 3)

        Q = quantity1 - quantity2
        assert allclose(Q.vertex_values,
                        quantity1.vertex_values - quantity2.vertex_values)

        #Scaling
        Q = quantity1*3
        assert allclose(Q.vertex_values, quantity1.vertex_values*3)
        assert allclose(Q.centroid_values, quantity1.centroid_values*3)
        assert allclose(Q.edge_values, quantity1.edge_values*3)
        Q = 3*quantity1
        assert allclose(Q.vertex_values, quantity1.vertex_values*3)

        #Multiplication
        Q = quantity1 * quantity2
        #print Q.vertex_values
        #print Q.centroid_values
        #print quantity1.centroid_values
        #print quantity2.centroid_values

        assert allclose(Q.vertex_values,
                        quantity1.vertex_values * quantity2.vertex_values)

        #Linear combinations
        Q = 4*quantity1 + 2
        assert allclose(Q.vertex_values,
                        4*quantity1.vertex_values + 2)

        Q = quantity1*quantity2 + 2
        assert allclose(Q.vertex_values,
                        quantity1.vertex_values * quantity2.vertex_values + 2)

        Q = quantity1*quantity2 + quantity3
        assert allclose(Q.vertex_values,
                        quantity1.vertex_values * quantity2.vertex_values +
                        quantity3.vertex_values)
        Q = quantity1*quantity2 + 3*quantity3
        assert allclose(Q.vertex_values,
                        quantity1.vertex_values * quantity2.vertex_values +
                        3*quantity3.vertex_values)
        Q = quantity1*quantity2 + 3*quantity3 + 5.0
        assert allclose(Q.vertex_values,
                        quantity1.vertex_values * quantity2.vertex_values +
                        3*quantity3.vertex_values + 5)

        Q = quantity1*quantity2 - quantity3
        assert allclose(Q.vertex_values,
                        quantity1.vertex_values * quantity2.vertex_values -
                        quantity3.vertex_values)
        Q = 1.5*quantity1*quantity2 - 3*quantity3 + 5.0
        assert allclose(Q.vertex_values,
                        1.5*quantity1.vertex_values * quantity2.vertex_values -
                        3*quantity3.vertex_values + 5)

        #Try combining quantities and arrays and scalars
        Q = 1.5*quantity1*quantity2.vertex_values -\
            3*quantity3.vertex_values + 5.0
        assert allclose(Q.vertex_values,
                        1.5*quantity1.vertex_values * quantity2.vertex_values -
                        3*quantity3.vertex_values + 5)


        #Powers
        Q = quantity1**2
        assert allclose(Q.vertex_values, quantity1.vertex_values**2)

        Q = quantity1**2 +quantity2**2
        assert allclose(Q.vertex_values,
                        quantity1.vertex_values**2 + \
                        quantity2.vertex_values**2)

        Q = (quantity1**2 +quantity2**2)**0.5
        assert allclose(Q.vertex_values,
                        (quantity1.vertex_values**2 + \
                        quantity2.vertex_values**2)**0.5)







    def test_gradient(self):
        quantity = Conserved_quantity(self.mesh4)

        #Set up for a gradient of (2,0) at mid triangle
        quantity.set_values([2.0, 4.0, 6.0, 2.0],
                            location = 'centroids')


        #Gradients
        a, b = quantity.compute_gradients()

        #print self.mesh4.centroid_coordinates
        #print a, b

        #The central triangle (1)
        #(using standard gradient based on neigbours controid values)
        assert allclose(a[1], 2.0)
        assert allclose(b[1], 0.0)


        #Left triangle (0) using two point gradient
        #q0 = q1 + a*(x0-x1) + b*(y0-y1)  <=>
        #2  = 4  + a*(-2/3)  + b*(-2/3)
        assert allclose(a[0] + b[0], 3)
        #From orthogonality (a*(y0-y1) + b*(x0-x1) == 0)
        assert allclose(a[0] - b[0], 0)


        #Right triangle (2) using two point gradient
        #q2 = q1 + a*(x2-x1) + b*(y2-y1)  <=>
        #6  = 4  + a*(4/3)  + b*(-2/3)
        assert allclose(2*a[2] - b[2], 3)
        #From orthogonality (a*(y1-y2) + b*(x2-x1) == 0)
        assert allclose(a[2] + 2*b[2], 0)


        #Top triangle (3) using two point gradient
        #q3 = q1 + a*(x3-x1) + b*(y3-y1)  <=>
        #2  = 4  + a*(-2/3)  + b*(4/3)
        assert allclose(a[3] - 2*b[3], 3)
        #From orthogonality (a*(y1-y3) + b*(x3-x1) == 0)
        assert allclose(2*a[3] + b[3], 0)



        #print a, b
        quantity.extrapolate_second_order()

        #Apply q(x,y) = qc + a*(x-xc) + b*(y-yc)
        assert allclose(quantity.vertex_values[0,:], [3., 0.,  3.])
        assert allclose(quantity.vertex_values[1,:], [4./3, 16./3,  16./3])


        #a = 1.2, b=-0.6
        #q(4,0) = 6 + a*(4 - 8/3) + b*(-2/3)
        assert allclose(quantity.vertex_values[2,2], 8)


    def test_second_order_extrapolation2(self):
        quantity = Conserved_quantity(self.mesh4)

        #Set up for a gradient of (3,1), f(x) = 3x+y
        quantity.set_values([2.0+2.0/3, 4.0+4.0/3, 8.0+2.0/3, 2.0+8.0/3],
                            location = 'centroids')

        #Gradients
        a, b = quantity.compute_gradients()

        #print a, b

        assert allclose(a[1], 3.0)
        assert allclose(b[1], 1.0)

        #Work out the others

        quantity.extrapolate_second_order()

        #print quantity.vertex_values
        assert allclose(quantity.vertex_values[1,0], 2.0)
        assert allclose(quantity.vertex_values[1,1], 6.0)
        assert allclose(quantity.vertex_values[1,2], 8.0)



    def test_first_order_extrapolator(self):
        quantity = Conserved_quantity(self.mesh4)

        #Test centroids
        quantity.set_values([1.,2.,3.,4.], location = 'centroids')
        assert allclose(quantity.centroid_values, [1, 2, 3, 4]) #Centroid

        #Extrapolate
        quantity.extrapolate_first_order()

        #Check vertices but not edge values
        assert allclose(quantity.vertex_values,
                        [[1,1,1], [2,2,2], [3,3,3], [4, 4, 4]])


    def test_second_order_extrapolator(self):
        quantity = Conserved_quantity(self.mesh4)

        #Set up for a gradient of (3,0) at mid triangle
        quantity.set_values([2.0, 4.0, 8.0, 2.0],
                            location = 'centroids')



        quantity.extrapolate_second_order()
        quantity.limit()


        #Assert that central triangle is limited by neighbours
        assert quantity.vertex_values[1,0] >= quantity.vertex_values[0,0]
        assert quantity.vertex_values[1,0] >= quantity.vertex_values[3,1]

        assert quantity.vertex_values[1,1] <= quantity.vertex_values[2,1]
        assert quantity.vertex_values[1,1] >= quantity.vertex_values[0,2]

        assert quantity.vertex_values[1,2] <= quantity.vertex_values[2,0]
        assert quantity.vertex_values[1,2] >= quantity.vertex_values[3,1]


        #Assert that quantities are conserved
        from Numeric import sum
        for k in range(quantity.centroid_values.shape[0]):
            assert allclose (quantity.centroid_values[k],
                             sum(quantity.vertex_values[k,:])/3)





    def test_limiter(self):
        quantity = Conserved_quantity(self.mesh4)

        #Create a deliberate overshoot (e.g. from gradient computation)
        quantity.set_values([[3,0,3], [2,2,6], [5,3,8], [8,3,5]])


        #Limit
        quantity.limit()

        #Assert that central triangle is limited by neighbours
        assert quantity.vertex_values[1,0] >= quantity.vertex_values[0,0]
        assert quantity.vertex_values[1,0] <= quantity.vertex_values[3,1]

        assert quantity.vertex_values[1,1] <= quantity.vertex_values[2,1]
        assert quantity.vertex_values[1,1] >= quantity.vertex_values[0,2]

        assert quantity.vertex_values[1,2] <= quantity.vertex_values[2,0]
        assert quantity.vertex_values[1,2] <= quantity.vertex_values[3,1]



        #Assert that quantities are conserved
        from Numeric import sum
        for k in range(quantity.centroid_values.shape[0]):
            assert allclose (quantity.centroid_values[k],
                             sum(quantity.vertex_values[k,:])/3)


    def test_limiter2(self):
        """Taken from test_shallow_water
        """
        quantity = Conserved_quantity(self.mesh4)
        quantity.domain.beta_w = 0.9
        
        #Test centroids
        quantity.set_values([2.,4.,8.,2.], location = 'centroids')
        assert allclose(quantity.centroid_values, [2, 4, 8, 2]) #Centroid


        #Extrapolate
        quantity.extrapolate_second_order()

        assert allclose(quantity.vertex_values[1,:], [0.0, 6, 6])

        #Limit
        quantity.limit()

        # limited value for beta_w = 0.9
        
        assert allclose(quantity.vertex_values[1,:], [2.2, 4.9, 4.9])
        # limited values for beta_w = 0.5
        #assert allclose(quantity.vertex_values[1,:], [3.0, 4.5, 4.5])


        #Assert that quantities are conserved
        from Numeric import sum
        for k in range(quantity.centroid_values.shape[0]):
            assert allclose (quantity.centroid_values[k],
                             sum(quantity.vertex_values[k,:])/3)





    def test_distribute_first_order(self):
        quantity = Conserved_quantity(self.mesh4)

        #Test centroids
        quantity.set_values([1.,2.,3.,4.], location = 'centroids')
        assert allclose(quantity.centroid_values, [1, 2, 3, 4]) #Centroid


        #Extrapolate
        quantity.extrapolate_first_order()

        #Interpolate
        quantity.interpolate_from_vertices_to_edges()

        assert allclose(quantity.vertex_values,
                        [[1,1,1], [2,2,2], [3,3,3], [4, 4, 4]])
        assert allclose(quantity.edge_values, [[1,1,1], [2,2,2],
                                               [3,3,3], [4, 4, 4]])


    def test_distribute_second_order(self):
        quantity = Conserved_quantity(self.mesh4)

        #Test centroids
        quantity.set_values([2.,4.,8.,2.], location = 'centroids')
        assert allclose(quantity.centroid_values, [2, 4, 8, 2]) #Centroid


        #Extrapolate
        quantity.extrapolate_second_order()

        assert allclose(quantity.vertex_values[1,:], [0.0, 6, 6])


    def test_update_explicit(self):
        quantity = Conserved_quantity(self.mesh4)

        #Test centroids
        quantity.set_values([1.,2.,3.,4.], location = 'centroids')
        assert allclose(quantity.centroid_values, [1, 2, 3, 4]) #Centroid

        #Set explicit_update
        quantity.explicit_update = array( [1.,1.,1.,1.] )

        #Update with given timestep
        quantity.update(0.1)

        x = array([1, 2, 3, 4]) + array( [.1,.1,.1,.1] )
        assert allclose( quantity.centroid_values, x)

    def test_update_semi_implicit(self):
        quantity = Conserved_quantity(self.mesh4)

        #Test centroids
        quantity.set_values([1.,2.,3.,4.], location = 'centroids')
        assert allclose(quantity.centroid_values, [1, 2, 3, 4]) #Centroid

        #Set semi implicit update
        quantity.semi_implicit_update = array([1.,1.,1.,1.])

        #Update with given timestep
        timestep = 0.1
        quantity.update(timestep)

        sem = array([1.,1.,1.,1.])/array([1, 2, 3, 4])
        denom = ones(4, Float)-timestep*sem

        x = array([1, 2, 3, 4])/denom
        assert allclose( quantity.centroid_values, x)


    def test_both_updates(self):
        quantity = Conserved_quantity(self.mesh4)

        #Test centroids
        quantity.set_values([1.,2.,3.,4.], location = 'centroids')
        assert allclose(quantity.centroid_values, [1, 2, 3, 4]) #Centroid

        #Set explicit_update
        quantity.explicit_update = array( [4.,3.,2.,1.] )

        #Set semi implicit update
        quantity.semi_implicit_update = array( [1.,1.,1.,1.] )

        #Update with given timestep
        timestep = 0.1
        quantity.update(0.1)

        sem = array([1.,1.,1.,1.])/array([1, 2, 3, 4])
        denom = ones(4, Float)-timestep*sem

        x = array([1., 2., 3., 4.])
        x /= denom
        x += timestep*array( [4.0, 3.0, 2.0, 1.0] )

        assert allclose( quantity.centroid_values, x)




    #Test smoothing
    def test_smoothing(self):

        from mesh_factory import rectangular
        from shallow_water import Domain, Transmissive_boundary
        from Numeric import zeros, Float
        from anuga.utilities.numerical_tools import mean

        #Create basic mesh
        points, vertices, boundary = rectangular(2, 2)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.default_order=2
        domain.reduction = mean


        #Set some field values
        domain.set_quantity('elevation', lambda x,y: x)
        domain.set_quantity('friction', 0.03)


        ######################
        # Boundary conditions
        B = Transmissive_boundary(domain)
        domain.set_boundary( {'left': B, 'right': B, 'top': B, 'bottom': B})


        ######################
        #Initial condition - with jumps

        bed = domain.quantities['elevation'].vertex_values
        stage = zeros(bed.shape, Float)

        h = 0.03
        for i in range(stage.shape[0]):
            if i % 2 == 0:
                stage[i,:] = bed[i,:] + h
            else:
                stage[i,:] = bed[i,:]

        domain.set_quantity('stage', stage)

        stage = domain.quantities['stage']

        #Get smoothed stage
        A, V = stage.get_vertex_values(xy=False, smooth=True)
        Q = stage.vertex_values


        assert A.shape[0] == 9
        assert V.shape[0] == 8
        assert V.shape[1] == 3

        #First four points
        assert allclose(A[0], (Q[0,2] + Q[1,1])/2)
        assert allclose(A[1], (Q[1,0] + Q[3,1] + Q[2,2])/3)
        assert allclose(A[2], Q[3,0])
        assert allclose(A[3], (Q[0,0] + Q[5,1] + Q[4,2])/3)

        #Center point
        assert allclose(A[4], (Q[0,1] + Q[1,2] + Q[2,0] +\
                               Q[5,0] + Q[6,2] + Q[7,1])/6)


        #Check V
        assert allclose(V[0,:], [3,4,0])
        assert allclose(V[1,:], [1,0,4])
        assert allclose(V[2,:], [4,5,1])
        assert allclose(V[3,:], [2,1,5])
        assert allclose(V[4,:], [6,7,3])
        assert allclose(V[5,:], [4,3,7])
        assert allclose(V[6,:], [7,8,4])
        assert allclose(V[7,:], [5,4,8])

        #Get smoothed stage with XY
        X, Y, A1, V1 = stage.get_vertex_values(xy=True, smooth=True)

        assert allclose(A, A1)
        assert allclose(V, V1)

        #Check XY
        assert allclose(X[4], 0.5)
        assert allclose(Y[4], 0.5)

        assert allclose(X[7], 1.0)
        assert allclose(Y[7], 0.5)




    def test_vertex_values_no_smoothing(self):

        from mesh_factory import rectangular
        from shallow_water import Domain, Transmissive_boundary
        from Numeric import zeros, Float
        from anuga.utilities.numerical_tools import mean


        #Create basic mesh
        points, vertices, boundary = rectangular(2, 2)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.default_order=2
        domain.reduction = mean


        #Set some field values
        domain.set_quantity('elevation', lambda x,y: x)
        domain.set_quantity('friction', 0.03)


        ######################
        #Initial condition - with jumps

        bed = domain.quantities['elevation'].vertex_values
        stage = zeros(bed.shape, Float)

        h = 0.03
        for i in range(stage.shape[0]):
            if i % 2 == 0:
                stage[i,:] = bed[i,:] + h
            else:
                stage[i,:] = bed[i,:]

        domain.set_quantity('stage', stage)

        #Get stage
        stage = domain.quantities['stage']
        A, V = stage.get_vertex_values(xy=False, smooth=False)
        Q = stage.vertex_values.flat

        for k in range(8):
            assert allclose(A[k], Q[k])


        for k in range(8):
            assert V[k, 0] == 3*k
            assert V[k, 1] == 3*k+1
            assert V[k, 2] == 3*k+2



        X, Y, A1, V1 = stage.get_vertex_values(xy=True, smooth=False)


        assert allclose(A, A1)
        assert allclose(V, V1)

        #Check XY
        assert allclose(X[1], 0.5)
        assert allclose(Y[1], 0.5)
        assert allclose(X[4], 0.0)
        assert allclose(Y[4], 0.0)
        assert allclose(X[12], 1.0)
        assert allclose(Y[12], 0.0)



    def set_array_values_by_index(self):

        from mesh_factory import rectangular
        from shallow_water import Domain
        from Numeric import zeros, Float

        #Create basic mesh
        points, vertices, boundary = rectangular(1, 1)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        #print "domain.number_of_elements ",domain.number_of_elements
        quantity = Quantity(domain,[[1,1,1],[2,2,2]])
        value = [7]
        indices = [1]
        quantity.set_array_values_by_index(value,
                                           location = 'centroids',
                                           indices = indices)
        #print "quantity.centroid_values",quantity.centroid_values

        assert allclose(quantity.centroid_values, [1,7])

        quantity.set_array_values([15,20,25], indices = indices)
        assert allclose(quantity.centroid_values, [1,20])

        quantity.set_array_values([15,20,25], indices = indices)
        assert allclose(quantity.centroid_values, [1,20])

    def test_setting_some_vertex_values(self):
        """
        set values based on triangle lists.
        """
        from mesh_factory import rectangular
        from shallow_water import Domain
        from Numeric import zeros, Float

        #Create basic mesh
        points, vertices, boundary = rectangular(1, 3)
        #print "vertices",vertices
        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        #print "domain.number_of_elements ",domain.number_of_elements
        quantity = Quantity(domain,[[1,1,1],[2,2,2],[3,3,3],
                                    [4,4,4],[5,5,5],[6,6,6]])


        # Check that constants work
        value = 7
        indices = [1]
        quantity.set_values(value,
                            location = 'centroids',
                            indices = indices)
        #print "quantity.centroid_values",quantity.centroid_values
        assert allclose(quantity.centroid_values, [1,7,3,4,5,6])
        
        value = [7]
        indices = [1]
        quantity.set_values(value,
                            location = 'centroids',
                            indices = indices)
        #print "quantity.centroid_values",quantity.centroid_values
        assert allclose(quantity.centroid_values, [1,7,3,4,5,6])

        value = [[15,20,25]]
        quantity.set_values(value, indices = indices)
        #print "1 quantity.vertex_values",quantity.vertex_values
        assert allclose(quantity.vertex_values[1], value[0])


        #print "quantity",quantity.vertex_values
        values = [10,100,50]
        quantity.set_values(values, indices = [0,1,5], location = 'centroids')
        #print "2 quantity.vertex_values",quantity.vertex_values
        assert allclose(quantity.vertex_values[0], [10,10,10])
        assert allclose(quantity.vertex_values[5], [50,50,50])
        #quantity.interpolate()
        #print "quantity.centroid_values",quantity.centroid_values
        assert allclose(quantity.centroid_values, [10,100,3,4,5,50])


        quantity = Quantity(domain,[[1,1,1],[2,2,2],[3,3,3],
                                    [4,4,4],[5,5,5],[6,6,6]])
        values = [10,100,50]
        #this will be per unique vertex, indexing the vertices
        #print "quantity.vertex_values",quantity.vertex_values
        quantity.set_values(values, indices = [0,1,5])
        #print "quantity.vertex_values",quantity.vertex_values
        assert allclose(quantity.vertex_values[0], [1,50,10])
        assert allclose(quantity.vertex_values[5], [6,6,6])
        assert allclose(quantity.vertex_values[1], [100,10,50])

        quantity = Quantity(domain,[[1,1,1],[2,2,2],[3,3,3],
                                    [4,4,4],[5,5,5],[6,6,6]])
        values = [[31,30,29],[400,400,400],[1000,999,998]]
        quantity.set_values(values, indices = [3,3,5])
        quantity.interpolate()
        assert allclose(quantity.centroid_values, [1,2,3,400,5,999])

        values = [[1,1,1],[2,2,2],[3,3,3],
                                    [4,4,4],[5,5,5],[6,6,6]]
        quantity.set_values(values)

        # testing the standard set values by vertex
        # indexed by vertex_id in general_mesh.coordinates
        values = [0,1,2,3,4,5,6,7]

        quantity.set_values(values)
        #print "1 quantity.vertex_values",quantity.vertex_values
        assert allclose(quantity.vertex_values,[[ 4.,  5.,  0.],
                                                [ 1.,  0.,  5.],
                                                [ 5.,  6.,  1.],
                                                [ 2.,  1.,  6.],
                                                [ 6.,  7.,  2.],
                                                [ 3.,  2.,  7.]])

    def test_setting_unique_vertex_values(self):
        """
        set values based on unique_vertex lists.
        """
        from mesh_factory import rectangular
        from shallow_water import Domain
        from Numeric import zeros, Float

        #Create basic mesh
        points, vertices, boundary = rectangular(1, 3)
        #print "vertices",vertices
        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        #print "domain.number_of_elements ",domain.number_of_elements
        quantity = Quantity(domain,[[0,0,0],[1,1,1],[2,2,2],[3,3,3],
                                    [4,4,4],[5,5,5]])
        value = 7
        indices = [1,5]
        quantity.set_values(value,
                            location = 'unique vertices',
                            indices = indices)
        #print "quantity.centroid_values",quantity.centroid_values
        assert allclose(quantity.vertex_values[0], [0,7,0])
        assert allclose(quantity.vertex_values[1], [7,1,7])
        assert allclose(quantity.vertex_values[2], [7,2,7])


    def test_get_values(self):
        """
        get values based on triangle lists.
        """
        from mesh_factory import rectangular
        from shallow_water import Domain
        from Numeric import zeros, Float

        #Create basic mesh
        points, vertices, boundary = rectangular(1, 3)

        #print "points",points
        #print "vertices",vertices
        #print "boundary",boundary

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        #print "domain.number_of_elements ",domain.number_of_elements
        quantity = Quantity(domain,[[0,0,0],[1,1,1],[2,2,2],[3,3,3],
                                    [4,4,4],[5,5,5]])

        #print "quantity.get_values(location = 'unique vertices')", \
        #      quantity.get_values(location = 'unique vertices')

        #print "quantity.get_values(location = 'unique vertices')", \
        #      quantity.get_values(indices=[0,1,2,3,4,5,6,7], \
        #                          location = 'unique vertices')

        answer = [0.5,2,4,5,0,1,3,4.5]
        assert allclose(answer,
                        quantity.get_values(location = 'unique vertices'))

        indices = [0,5,3]
        answer = [0.5,1,5]
        assert allclose(answer,
                        quantity.get_values(indices=indices, \
                                            location = 'unique vertices'))
        #print "quantity.centroid_values",quantity.centroid_values
        #print "quantity.get_values(location = 'centroids') ",\
        #      quantity.get_values(location = 'centroids')




    def test_get_values_2(self):
        """Different mesh (working with domain object) - also check centroids.
        """

        
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices)

        quantity = Quantity(domain)
        quantity.set_values(lambda x, y: x+2*y) #2 4 4 6
        
        assert allclose(quantity.get_values(location='centroids'), [2,4,4,6])
        assert allclose(quantity.get_values(location='centroids', indices=[1,3]), [4,6])


        assert allclose(quantity.get_values(location='vertices'), [[4,0,2],
                                                                   [4,2,6],
                                                                   [6,2,4],
                                                                   [8,4,6]])
        
        assert allclose(quantity.get_values(location='vertices', indices=[1,3]), [[4,2,6],
                                                                                  [8,4,6]])


        assert allclose(quantity.get_values(location='edges'), [[1,3,2],
                                                                [4,5,3],
                                                                [3,5,4],
                                                                [5,7,6]])
        assert allclose(quantity.get_values(location='edges', indices=[1,3]),
                        [[4,5,3],
                         [5,7,6]])        

        # Check averaging over vertices
        #a: 0
        #b: (4+4+4)/3
        #c: (2+2+2)/3
        #d: 8
        #e: (6+6+6)/3        
        #f: 4
        assert(quantity.get_values(location='unique vertices'), [0, 4, 2, 8, 6, 4])        
                                                                                  
        




    def test_get_interpolated_values(self):

        from mesh_factory import rectangular
        from shallow_water import Domain
        from Numeric import zeros, Float

        #Create basic mesh
        points, vertices, boundary = rectangular(1, 3)
        domain = Domain(points, vertices, boundary)

        #Constant values
        quantity = Quantity(domain,[[0,0,0],[1,1,1],[2,2,2],[3,3,3],
                                    [4,4,4],[5,5,5]])

        

        # Get interpolated values at centroids
        interpolation_points = domain.get_centroid_coordinates()
        answer = quantity.get_values(location='centroids')

        
        #print quantity.get_values(points=interpolation_points)
        assert allclose(answer, quantity.get_values(interpolation_points=interpolation_points))


        #Arbitrary values
        quantity = Quantity(domain,[[0,1,2],[3,1,7],[2,1,2],[3,3,7],
                                    [1,4,-9],[2,5,0]])


        # Get interpolated values at centroids
        interpolation_points = domain.get_centroid_coordinates()
        answer = quantity.get_values(location='centroids')
        #print answer
        #print quantity.get_values(points=interpolation_points)
        assert allclose(answer, quantity.get_values(interpolation_points=interpolation_points))        
                        

        #FIXME TODO
        #indices = [0,5,3]
        #answer = [0.5,1,5]
        #assert allclose(answer,
        #                quantity.get_values(indices=indices, \
        #                                    location = 'unique vertices'))




    def test_get_interpolated_values_2(self):
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [[1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices)

        quantity = Quantity(domain)
        quantity.set_values(lambda x, y: x+2*y) #2 4 4 6

        #First pick one point
        x, y = 2.0/3, 8.0/3
        v = quantity.get_values(interpolation_points = [[x,y]])
        assert allclose(v, 6)        

        # Then another to test that algorithm won't blindly
        # reuse interpolation matrix
        x, y = 4.0/3, 4.0/3
        v = quantity.get_values(interpolation_points = [[x,y]])
        assert allclose(v, 4)        




    def test_getting_some_vertex_values(self):
        """
        get values based on triangle lists.
        """
        from mesh_factory import rectangular
        from shallow_water import Domain
        from Numeric import zeros, Float

        #Create basic mesh
        points, vertices, boundary = rectangular(1, 3)

        #print "points",points
        #print "vertices",vertices
        #print "boundary",boundary

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        #print "domain.number_of_elements ",domain.number_of_elements
        quantity = Quantity(domain,[[1,1,1],[2,2,2],[3,3,3],
                                    [4,4,4],[5,5,5],[6,6,6]])
        value = [7]
        indices = [1]
        quantity.set_values(value,
                            location = 'centroids',
                            indices = indices)
        #print "quantity.centroid_values",quantity.centroid_values
        #print "quantity.get_values(location = 'centroids') ",\
        #      quantity.get_values(location = 'centroids')
        assert allclose(quantity.centroid_values,
                        quantity.get_values(location = 'centroids'))


        value = [[15,20,25]]
        quantity.set_values(value, indices = indices)
        #print "1 quantity.vertex_values",quantity.vertex_values
        assert allclose(quantity.vertex_values, quantity.get_values())

        assert allclose(quantity.edge_values,
                        quantity.get_values(location = 'edges'))

        # get a subset of elements
        subset = quantity.get_values(location='centroids', indices=[0,5])
        answer = [quantity.centroid_values[0],quantity.centroid_values[5]]
        assert allclose(subset, answer)


        subset = quantity.get_values(location='edges', indices=[0,5])
        answer = [quantity.edge_values[0],quantity.edge_values[5]]
        #print "subset",subset
        #print "answer",answer
        assert allclose(subset, answer)

        subset = quantity.get_values( indices=[1,5])
        answer = [quantity.vertex_values[1],quantity.vertex_values[5]]
        #print "subset",subset
        #print "answer",answer
        assert allclose(subset, answer)




#-------------------------------------------------------------
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_Quantity, 'test')

    #suite = unittest.makeSuite(Test_Quantity, 'test_set_values_from_UTM')
    #print "restricted test"
    #suite = unittest.makeSuite(Test_Quantity,'verbose_test_set_values_from_UTM_pts')
    runner = unittest.TextTestRunner()
    runner.run(suite)

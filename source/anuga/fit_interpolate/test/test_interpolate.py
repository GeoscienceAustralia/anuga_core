#!/usr/bin/env python

#TEST

#import time, os


import sys
import os
import unittest
import tempfile
import csv
import time

import numpy as num
from math import sqrt

# ANUGA code imports
import anuga
from anuga.fit_interpolate.interpolate import Interpolate
from anuga.fit_interpolate.interpolate import Interpolation_function
from anuga.fit_interpolate.interpolate import interpolate
from anuga.fit_interpolate.interpolate import interpolate_sww2csv

from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.utilities.numerical_tools import mean, NAN
from anuga.file.sww import SWW_file
from anuga.shallow_water.shallow_water_domain import Domain
from anuga.abstract_2d_finite_volumes.quantity import Quantity
from anuga.geospatial_data.geospatial_data import Geospatial_data
from anuga.pmesh.mesh import Mesh
from anuga.file.netcdf import NetCDFFile

def distance(x, y):
    return sqrt(num.sum((num.array(x)-num.array(y))**2))

def linear_function(point):
    point = num.array(point)
    return point[:,0]+point[:,1]


class Test_Interpolate(unittest.TestCase):

    def setUp(self):

        import time
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular


        #Create basic mesh
        points, vertices, boundary = rectangular(2, 2)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.default_order=2


        #Set some field values
        domain.set_quantity('elevation', lambda x,y: -x)
        domain.set_quantity('friction', 0.03)


        ######################
        # Boundary conditions
        B = anuga.Transmissive_boundary(domain)
        domain.set_boundary( {'left': B, 'right': B, 'top': B, 'bottom': B})


        ######################
        #Initial condition - with jumps

        bed = domain.quantities['elevation'].vertex_values
        stage = num.zeros(bed.shape, num.float)

        h = 0.3
        for i in range(stage.shape[0]):
            if i % 2 == 0:
                stage[i,:] = bed[i,:] + h
            else:
                stage[i,:] = bed[i,:]

        domain.set_quantity('stage', stage)

        domain.distribute_to_vertices_and_edges()


        self.domain = domain

        C = domain.get_vertex_coordinates()
        self.X = C[:,0:6:2].copy()
        self.Y = C[:,1:6:2].copy()

        self.F = bed



    def tearDown(self):
        pass

    def test_datapoint_at_centroid(self):
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        points = [a, b, c]
        vertices = [ [1,0,2] ]   #bac

        data = [ [2.0/3, 2.0/3] ] #Use centroid as one data point

        interp = Interpolate(points, vertices)
        A, _, _, _ = interp._build_interpolation_matrix_A(data)
        assert num.allclose(A.todense(), [[1./3, 1./3, 1./3]])

    def test_datapoint_in_hole(self):
        # create 3 right-angled triangles arranged in a bigger triangle       
        a = [0.0, 0.0] #0
        b = [0.0, 2.0] #1
        c = [2.0,0.0]  #2
        d = [0.0,4.0]  #3
        e = [2.0,2.0]  #4
        f = [4.0,0.0]  #5
        points = [a, b, c, d, e, f]
        vertices = [ [1,0,2], [3,1,4], [4,2,5] ]   #bac dbe ecf

        point_in_hole = [1.5, 1.5] # a point in the hole
        data = [ [20, 20], [0.3, 0.3], point_in_hole, [2.5, 0.3], [30, 30] ] #some points inside and outside the hole

        # any function for the vertices, we don't care about the result
        f = num.array([linear_function(points), 2*linear_function(points)])
        f = num.transpose(f)

        interp = Interpolate(points, vertices)
        interp.interpolate(f, data)
        
        assert not set(interp.inside_poly_indices).intersection(set(interp.outside_poly_indices)), \
               'Some points are in both lists!'    
        assert len(interp.inside_poly_indices) == 2
        assert len(interp.outside_poly_indices) == 3
        
        interp.outside_poly_indices.sort()
        assert interp.outside_poly_indices[1] == 2, \
               'third outside point should be inside the hole!'

    def test_simple_interpolation_example(self):
        
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(1, 3)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)

        #----------------
        #Constant values
        #----------------        
        quantity = Quantity(domain,[[0,0,0],[1,1,1],[2,2,2],[3,3,3],
                                    [4,4,4],[5,5,5]])


        x, y, vertex_values, triangles = quantity.get_vertex_values(xy=True, smooth=False)
        vertex_coordinates = num.concatenate( (x[:, num.newaxis], y[:, num.newaxis]), axis=1 )
        # FIXME: This concat should roll into get_vertex_values


        # Get interpolated values at centroids
        interpolation_points = domain.get_centroid_coordinates()
        answer = quantity.get_values(location='centroids')

        I = Interpolate(vertex_coordinates, triangles)
        result = I.interpolate(vertex_values, interpolation_points)
        assert num.allclose(result, answer)


        #----------------
        # Variable values
        #----------------
        quantity = Quantity(domain,[[0,1,2],[3,1,7],[2,1,2],[3,3,7],
                                    [1,4,-9],[2,5,0]])
        
        x, y, vertex_values, triangles = quantity.get_vertex_values(xy=True, smooth=False)
        vertex_coordinates = num.concatenate( (x[:, num.newaxis], y[:, num.newaxis]), axis=1 )
        # FIXME: This concat should roll into get_vertex_values


        # Get interpolated values at centroids
        interpolation_points = domain.get_centroid_coordinates()
        answer = quantity.get_values(location='centroids')

        I = Interpolate(vertex_coordinates, triangles)
        result = I.interpolate(vertex_values, interpolation_points)
        assert num.allclose(result, answer)        
        

    def test_simple_interpolation_example_using_direct_interface(self):
        
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(1, 3)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)

        #----------------
        # Constant values
        #----------------        
        quantity = Quantity(domain,[[0,0,0],[1,1,1],[2,2,2],[3,3,3],
                                    [4,4,4],[5,5,5]])


        x, y, vertex_values, triangles = quantity.get_vertex_values(xy=True, smooth=False)
        vertex_coordinates = num.concatenate( (x[:, num.newaxis], y[:, num.newaxis]), axis=1 )
        # FIXME: This concat should roll into get_vertex_values


        # Get interpolated values at centroids
        interpolation_points = domain.get_centroid_coordinates()
        answer = quantity.get_values(location='centroids')

        result = interpolate(vertex_coordinates, triangles, vertex_values, interpolation_points)
        assert num.allclose(result, answer)


        #----------------
        # Variable values
        #----------------
        quantity = Quantity(domain,[[0,1,2],[3,1,7],[2,1,2],[3,3,7],
                                    [1,4,-9],[2,5,0]])
        
        x, y, vertex_values, triangles = quantity.get_vertex_values(xy=True, smooth=False)
        vertex_coordinates = num.concatenate( (x[:, num.newaxis], y[:, num.newaxis]), axis=1 )
        # FIXME: This concat should roll into get_vertex_values


        # Get interpolated values at centroids
        interpolation_points = domain.get_centroid_coordinates()
        answer = quantity.get_values(location='centroids')

        result = interpolate(vertex_coordinates, triangles,
                             vertex_values, interpolation_points)
        assert num.allclose(result, answer)        
        
        
    def test_simple_interpolation_example_using_direct_interface_and_caching(self):
        
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(1, 3)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)

        #----------------
        # First call
        #----------------
        quantity = Quantity(domain,[[0,1,2],[3,1,7],[2,1,2],[3,3,7],
                                    [1,4,-9],[2,5,0]])
        
        x, y, vertex_values, triangles = quantity.get_vertex_values(xy=True, smooth=False)
        vertex_coordinates = num.concatenate( (x[:, num.newaxis], y[:, num.newaxis]), axis=1 )
        # FIXME: This concat should roll into get_vertex_values

        # Get interpolated values at centroids
        interpolation_points = domain.get_centroid_coordinates()
        answer = quantity.get_values(location='centroids')
        result = interpolate(vertex_coordinates, triangles,
                             vertex_values, interpolation_points,
                             use_cache=True,
                             verbose=False)
        assert num.allclose(result, answer)                

        # Second call using the cache
        result = interpolate(vertex_coordinates, triangles,
                             vertex_values, interpolation_points,
                             use_cache=True,
                             verbose=False)
        assert num.allclose(result, answer)                        
        
        
    def test_quad_tree(self):
        p0 = [-10.0, -10.0]
        p1 = [20.0, -10.0]
        p2 = [-10.0, 20.0]
        p3 = [10.0, 50.0]
        p4 = [30.0, 30.0]
        p5 = [50.0, 10.0]
        p6 = [40.0, 60.0]
        p7 = [60.0, 40.0]
        p8 = [-66.0, 20.0]
        p9 = [10.0, -66.0]

        points = [p0, p1, p2, p3, p4, p5, p6, p7, p8, p9]
        triangles = [ [0, 1, 2],
                      [3, 2, 4],
                      [4, 2, 1],
                      [4, 1, 5],
                      [3, 4, 6],
                      [6, 4, 7],
                      [7, 4, 5],
                      [8, 0, 2],
                      [0, 9, 1]]

        data = [ [4,4] ]
        interp = Interpolate(points, triangles)
        #print "PDSG - interp.get_A()", interp.get_A()
        answer =  [ [ 0.06666667,  0.46666667,  0.46666667,  0.,
                      0., 0. , 0., 0., 0., 0.]]

        A,_,_,_ = interp._build_interpolation_matrix_A(data)
        assert num.allclose(A.todense(), answer)
        
        #interp.set_point_coordinates([[-30, -30]]) #point outside of mesh
        #print "PDSG - interp.get_A()", interp.get_A()
        data = [[-30, -30]]
        answer =  [ [ 0.0,  0.0,  0.0,  0.,
                      0., 0. , 0., 0., 0., 0.]]
        
        A,_,_,_ = interp._build_interpolation_matrix_A(data)        
        assert num.allclose(A.todense(), answer)


        #point outside of quad tree root cell
        #interp.set_point_coordinates([[-70, -70]])
        #print "PDSG - interp.get_A()", interp.get_A()
        data = [[-70, -70]]
        answer =  [ [ 0.0,  0.0,  0.0,  0.,
                      0., 0. , 0., 0., 0., 0.]]
                      
        A,_,_,_ = interp._build_interpolation_matrix_A(data)        
        assert num.allclose(A.todense(), answer)


    def test_datapoints_at_vertices(self):
        #Test that data points coinciding with vertices yield a diagonal matrix
        

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        points = [a, b, c]
        vertices = [ [1,0,2] ]   #bac

        data = points #Use data at vertices

        interp = Interpolate(points, vertices)
        answer = [[1., 0., 0.],
                   [0., 1., 0.],
                   [0., 0., 1.]]
                   
        A,_,_,_ = interp._build_interpolation_matrix_A(data)
        assert num.allclose(A.todense(), answer)


    def test_datapoints_on_edge_midpoints(self):
        #Try datapoints midway on edges -
        #each point should affect two matrix entries equally
        

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        points = [a, b, c]
        vertices = [ [1,0,2] ]   #bac

        data = [ [0., 1.], [1., 0.], [1., 1.] ]
        answer =  [[0.5, 0.5, 0.0],  #Affects vertex 1 and 0
                    [0.5, 0.0, 0.5],  #Affects vertex 0 and 2
                    [0.0, 0.5, 0.5]]
        interp = Interpolate(points, vertices)

        A,_,_,_ = interp._build_interpolation_matrix_A(data)
        assert num.allclose(A.todense(), answer)

    def test_datapoints_on_edges(self):
        #Try datapoints on edges -
        #each point should affect two matrix entries in proportion
        

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        points = [a, b, c]
        vertices = [ [1,0,2] ]   #bac

        data = [ [0., 1.5], [1.5, 0.], [1.5, 0.5] ]
        answer =  [[0.25, 0.75, 0.0],  #Affects vertex 1 and 0
                   [0.25, 0.0, 0.75],  #Affects vertex 0 and 2
                   [0.0, 0.25, 0.75]]

        interp = Interpolate(points, vertices)

        A,_,_,_ = interp._build_interpolation_matrix_A(data)
        assert num.allclose(A.todense(), answer)


    def test_arbitrary_datapoints(self):
        #Try arbitrary datapoints
        

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        points = [a, b, c]
        vertices = [ [1,0,2] ]   #bac

        data = [ [0.2, 1.5], [0.123, 1.768], [1.43, 0.44] ]

        interp = Interpolate(points, vertices)
        #print "interp.get_A()", interp.get_A()
        
        A,_,_,_ = interp._build_interpolation_matrix_A(data)
        results = A.todense()
        assert num.allclose(num.sum(results, axis=1), 1.0)

        
    def test_arbitrary_datapoints_return_centroids(self):
        #Try arbitrary datapoints, making sure they return
        #an interpolation matrix for the intersected triangle's
        #centroid.
        
        a = [1.0, 0.0]
        b = [0.0, 3.0]
        c = [4.0,0.0]
        points = [a, b, c]
        vertices = [ [1,0,2] ]

        data = [ [1.2, 1.5], [1.123, 1.768], [2.43, 0.44] ]

        interp = Interpolate(points, vertices)
        
        third = [1.0/3.0, 1.0/3.0, 1.0/3.0]
        answer = [third, third, third]
        
        A,_,_,_ = interp._build_interpolation_matrix_A(data, output_centroids=True)
        results = A.todense()
        assert num.allclose(results, answer)        
        
        
    def test_arbitrary_datapoints_some_outside(self):
        #Try arbitrary datapoints one outside the triangle.
        #That one should be ignored
        

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        points = [a, b, c]
        vertices = [ [1,0,2] ]   #bac

        data = [ [0.2, 1.5], [0.123, 1.768], [1.43, 0.44], [5.0, 7.0]]

        interp = Interpolate(points, vertices)
        
        A,_,_,_ = interp._build_interpolation_matrix_A(data)
        results = A.todense()
        assert num.allclose(num.sum(results, axis=1), [1,1,1,0])



    # this causes a memory error in scipy.sparse
    def test_more_triangles(self):

        a = [-1.0, 0.0]
        b = [3.0, 4.0]
        c = [4.0,1.0]
        d = [-3.0, 2.0] #3
        e = [-1.0,-2.0]
        f = [1.0, -2.0] #5

        points = [a, b, c, d,e,f]
        triangles = [[0,1,3],[1,0,2],[0,4,5], [0,5,2]] #abd bac aef afc

        #Data points
        data = [ [-3., 2.0], [-2, 1], [0.0, 1], [0, 3], [2, 3], [-1.0/3,-4./3] ]
        interp = Interpolate(points, triangles)

        answer = [[0.0, 0.0, 0.0, 1.0, 0.0, 0.0],    #Affects point d
                  [0.5, 0.0, 0.0, 0.5, 0.0, 0.0],    #Affects points a and d
                  [0.75, 0.25, 0.0, 0.0, 0.0, 0.0],  #Affects points a and b
                  [0.0, 0.5, 0.0, 0.5, 0.0, 0.0],    #Affects points a and d
                  [0.25, 0.75, 0.0, 0.0, 0.0, 0.0],  #Affects points a and b
                  [1./3, 0.0, 0.0, 0.0, 1./3, 1./3]] #Affects points a, e and f


        A,_,_,_ = interp._build_interpolation_matrix_A(data)
        A = A.todense()
        for i in range(A.shape[0]):
            for j in range(A.shape[1]):
                if not num.allclose(A[i,j], answer[i][j]):
                    print i,j,':',A[i,j], answer[i][j]


        #results = interp._build_interpolation_matrix_A(data).todense()

        assert num.allclose(A, answer)
    
    def test_geo_ref(self):
        v0 = [0.0, 0.0]
        v1 = [0.0, 5.0]
        v2 = [5.0, 0.0]

        vertices_absolute = [v0, v1, v2]
        triangles = [ [1,0,2] ]   #bac

        geo = Geo_reference(57,100, 500)

        vertices = geo.change_points_geo_ref(vertices_absolute)
        #print "vertices",vertices 
        
        d0 = [1.0, 1.0]
        d1 = [1.0, 2.0]
        d2 = [3.0, 1.0]
        point_coords = [ d0, d1, d2]

        interp = Interpolate(vertices, triangles, mesh_origin=geo)
        f = linear_function(vertices_absolute)
        z = interp.interpolate(f, point_coords)
        answer = linear_function(point_coords)

        #print "z",z 
        #print "answer",answer 
        assert num.allclose(z, answer)

        
        z = interp.interpolate(f, point_coords, start_blocking_len = 2)
        answer = linear_function(point_coords)

        #print "z",z 
        #print "answer",answer 
        assert num.allclose(z, answer)
        
     
    def test_sigma_epsilon(self):
        """
        def test_sigma_epsilon(self):
            Testing ticket 168. I could not reduce the bug to this small
            test though.
        
        """
        v0 = [22031.25, 59687.5]
        v1 = [22500., 60000.]
        v2 = [22350.31640625, 59716.71484375]

        vertices = [v0, v1, v2]
        triangles = [ [1,0,2] ]   #bac

        
        point_coords = [[22050., 59700.]]

        interp = Interpolate(vertices, triangles)
        f = linear_function(vertices)
        z = interp.interpolate(f, point_coords)
        answer = linear_function(point_coords)

        #print "z",z 
        #print "answer",answer 
        assert num.allclose(z, answer)

        
        z = interp.interpolate(f, point_coords, start_blocking_len = 2)
        answer = linear_function(point_coords)

        #print "z",z 
        #print "answer",answer 
        assert num.allclose(z, answer)

        
    def test_Geospatial_verts(self):
        v0 = [0.0, 0.0]
        v1 = [0.0, 5.0]
        v2 = [5.0, 0.0]

        vertices_absolute = [v0, v1, v2]
        triangles = [ [1,0,2] ]   #bac

        geo = Geo_reference(57,100, 500)
        vertices = geo.change_points_geo_ref(vertices_absolute)
        geopoints = Geospatial_data(vertices,geo_reference = geo)
        #print "vertices",vertices 
        
        d0 = [1.0, 1.0]
        d1 = [1.0, 2.0]
        d2 = [3.0, 1.0]
        point_coords = [ d0, d1, d2]

        interp = Interpolate(geopoints, triangles)
        f = linear_function(vertices_absolute)
        z = interp.interpolate(f, point_coords)
        answer = linear_function(point_coords)

        #print "z",z 
        #print "answer",answer 
        assert num.allclose(z, answer)
        
        z = interp.interpolate(f, point_coords, start_blocking_len = 2)
        answer = linear_function(point_coords)

        #print "z",z 
        #print "answer",answer 
        assert num.allclose(z, answer)
        
    def test_interpolate_attributes_to_points(self):
        v0 = [0.0, 0.0]
        v1 = [0.0, 5.0]
        v2 = [5.0, 0.0]

        vertices = [v0, v1, v2]
        triangles = [ [1,0,2] ]   #bac

        d0 = [1.0, 1.0]
        d1 = [1.0, 2.0]
        d2 = [3.0, 1.0]
        point_coords = [ d0, d1, d2]

        interp = Interpolate(vertices, triangles)
        f = linear_function(vertices)
        z = interp.interpolate(f, point_coords)
        answer = linear_function(point_coords)

        #print "z",z 
        #print "answer",answer 
        assert num.allclose(z, answer)


        z = interp.interpolate(f, point_coords, start_blocking_len = 2)
        answer = linear_function(point_coords)

        #print "z",z 
        #print "answer",answer 
        assert num.allclose(z, answer)

    def test_interpolate_attributes_to_pointsII(self):
        a = [-1.0, 0.0]
        b = [3.0, 4.0]
        c = [4.0, 1.0]
        d = [-3.0, 2.0] #3
        e = [-1.0, -2.0]
        f = [1.0, -2.0] #5

        vertices = [a, b, c, d,e,f]
        triangles = [[0,1,3], [1,0,2], [0,4,5], [0,5,2]] #abd bac aef afc


        point_coords = [[-2.0, 2.0],
                        [-1.0, 1.0],
                        [0.0, 2.0],
                        [1.0, 1.0],
                        [2.0, 1.0],
                        [0.0, 0.0],
                        [1.0, 0.0],
                        [0.0, -1.0],
                        [-0.2, -0.5],
                        [-0.9, -1.5],
                        [0.5, -1.9],
                        [3.0, 1.0]]

        interp = Interpolate(vertices, triangles)
        f = linear_function(vertices)
        z = interp.interpolate(f, point_coords)
        answer = linear_function(point_coords)
        #print "z",z
        #print "answer",answer
        assert num.allclose(z, answer)

        z = interp.interpolate(f, point_coords, start_blocking_len = 2)
        answer = linear_function(point_coords)

        #print "z",z 
        #print "answer",answer 
        assert num.allclose(z, answer)
        
    def test_interpolate_attributes_to_pointsIII(self):
        #Test linear interpolation of known values at vertices to
        #new points inside a triangle
        
        a = [0.0, 0.0]
        b = [0.0, 5.0]
        c = [5.0, 0.0]
        d = [5.0, 5.0]

        vertices = [a, b, c, d]
        triangles = [ [1,0,2], [2,3,1] ]   #bac, cdb

        #Points within triangle 1
        d0 = [1.0, 1.0]
        d1 = [1.0, 2.0]
        d2 = [3.0, 1.0]

        #Point within triangle 2
        d3 = [4.0, 3.0]

        #Points on common edge
        d4 = [2.5, 2.5]
        d5 = [4.0, 1.0]

        #Point on common vertex
        d6 = [0., 5.]
        
        point_coords = [d0, d1, d2, d3, d4, d5, d6]

        interp = Interpolate(vertices, triangles)

        #Known values at vertices
        #Functions are x+y, x+2y, 2x+y, x-y-5
        f = [ [0., 0., 0., -5.],        # (0,0)
              [5., 10., 5., -10.],      # (0,5)
              [5., 5., 10.0, 0.],       # (5,0)
              [10., 15., 15., -5.]]     # (5,5)

        z = interp.interpolate(f, point_coords)
        answer = [ [2., 3., 3., -5.],   # (1,1)
                   [3., 5., 4., -6.],   # (1,2)
                   [4., 5., 7., -3.],   # (3,1)
                   [7., 10., 11., -4.], # (4,3)
                   [5., 7.5, 7.5, -5.], # (2.5, 2.5)
                   [5., 6., 9., -2.],   # (4,1)
                   [5., 10., 5., -10.]]  # (0,5)

        #print "***********"
        #print "z",z
        #print "answer",answer
        #print "***********"

        assert num.allclose(z, answer)


        z = interp.interpolate(f, point_coords, start_blocking_len = 2)

        #print "z",z 
        #print "answer",answer 
        assert num.allclose(z, answer)
        
    def test_interpolate_point_outside_of_mesh(self):
        #Test linear interpolation of known values at vertices to
        #new points inside a triangle
        
        a = [0.0, 0.0]
        b = [0.0, 5.0]
        c = [5.0, 0.0]
        d = [5.0, 5.0]

        vertices = [a, b, c, d]
        triangles = [ [1,0,2], [2,3,1] ]   #bac, cdb

        #Far away point
        d7 = [-1., -1.]
        
        point_coords = [ d7]
        interp = Interpolate(vertices, triangles)

        #Known values at vertices
        #Functions are x+y, x+2y, 2x+y, x-y-5
        f = [ [0., 0., 0., -5.],        # (0,0)
              [5., 10., 5., -10.],      # (0,5)
              [5., 5., 10.0, 0.],       # (5,0)
              [10., 15., 15., -5.]]     # (5,5)

        z = interp.interpolate(f, point_coords) #, verbose=True)
        answer = num.array([ [NAN, NAN, NAN, NAN]]) # (-1,-1)

        #print "***********"
        #print "z",z
        #print "answer",answer
        #print "***********"

        #Should an error message be returned if points are outside
        # of the mesh?
        # A warning message is printed, if verbose is on.

        for i in range(4):
            self.assertTrue( z[0,i] == answer[0,i], 'Fail!')
        
        z = interp.interpolate(f, point_coords, start_blocking_len = 2)

        #print "z",z 
        #print "answer",answer
        
        for i in range(4):
            self.assertTrue( z[0,i] == answer[0,i], 'Fail!')
        
        
    def test_interpolate_attributes_to_pointsIV(self):
        a = [-1.0, 0.0]
        b = [3.0, 4.0]
        c = [4.0, 1.0]
        d = [-3.0, 2.0] #3
        e = [-1.0, -2.0]
        f = [1.0, -2.0] #5

        vertices = [a, b, c, d,e,f]
        triangles = [[0,1,3], [1,0,2], [0,4,5], [0,5,2]] #abd bac aef afc


        point_coords = [[-2.0, 2.0],
                        [-1.0, 1.0],
                        [0.0, 2.0],
                        [1.0, 1.0],
                        [2.0, 1.0],
                        [0.0, 0.0],
                        [1.0, 0.0],
                        [0.0, -1.0],
                        [-0.2, -0.5],
                        [-0.9, -1.5],
                        [0.5, -1.9],
                        [3.0, 1.0]]

        interp = Interpolate(vertices, triangles)
        f = num.array([linear_function(vertices),2*linear_function(vertices)])
        f = num.transpose(f)
        #print "f",f
        z = interp.interpolate(f, point_coords)
        answer = [linear_function(point_coords),
                  2*linear_function(point_coords) ]
        answer = num.transpose(answer)
        #print "z",z
        #print "answer",answer
        assert num.allclose(z, answer)

        z = interp.interpolate(f, point_coords, start_blocking_len = 2)

        #print "z",z 
        #print "answer",answer 
        assert num.allclose(z, answer)

    def test_interpolate_blocking(self):
        a = [-1.0, 0.0]
        b = [3.0, 4.0]
        c = [4.0, 1.0]
        d = [-3.0, 2.0] #3
        e = [-1.0, -2.0]
        f = [1.0, -2.0] #5

        vertices = [a, b, c, d,e,f]
        triangles = [[0,1,3], [1,0,2], [0,4,5], [0,5,2]] #abd bac aef afc


        point_coords = [[-2.0, 2.0],
                        [-1.0, 1.0],
                        [0.0, 2.0],
                        [1.0, 1.0],
                        [2.0, 1.0],
                        [0.0, 0.0],
                        [1.0, 0.0],
                        [0.0, -1.0],
                        [-0.2, -0.5],
                        [-0.9, -1.5],
                        [0.5, -1.9],
                        [3.0, 1.0]]

        interp = Interpolate(vertices, triangles)
        f = num.array([linear_function(vertices),2*linear_function(vertices)])
        f = num.transpose(f)
        #print "f",f
        for blocking_max in range(len(point_coords)+2):
        #if True:
         #   blocking_max = 5
            z = interp.interpolate(f, point_coords,
                                   start_blocking_len=blocking_max)
            answer = [linear_function(point_coords),
                      2*linear_function(point_coords) ]
            answer = num.transpose(answer)
            #print "z",z
            #print "answer",answer
            assert num.allclose(z, answer)
            
        f = num.array([linear_function(vertices),2*linear_function(vertices),
                       2*linear_function(vertices) - 100])
        f = num.transpose(f)
        #print "f",f
        for blocking_max in range(len(point_coords)+2):
        #if True:
         #   blocking_max = 5
            z = interp.interpolate(f, point_coords,
                                   start_blocking_len=blocking_max)
            answer = num.array([linear_function(point_coords),
                                2*linear_function(point_coords) ,
                                2*linear_function(point_coords)-100])
            z = num.transpose(z)
            #print "z",z
            #print "answer",answer
            assert num.allclose(z, answer)

    def test_interpolate_geo_spatial(self):
        a = [-1.0, 0.0]
        b = [3.0, 4.0]
        c = [4.0, 1.0]
        d = [-3.0, 2.0] #3
        e = [-1.0, -2.0]
        f = [1.0, -2.0] #5

        vertices = [a, b, c, d,e,f]
        triangles = [[0,1,3], [1,0,2], [0,4,5], [0,5,2]] #abd bac aef afc


        point_coords_absolute = [[-2.0, 2.0],
                        [-1.0, 1.0],
                        [0.0, 2.0],
                        [1.0, 1.0],
                        [2.0, 1.0],
                        [0.0, 0.0],
                        [1.0, 0.0],
                        [0.0, -1.0],
                        [-0.2, -0.5],
                        [-0.9, -1.5],
                        [0.5, -1.9],
                        [3.0, 1.0]]

        geo = Geo_reference(57,100, 500)
        point_coords = geo.change_points_geo_ref(point_coords_absolute)
        point_coords = Geospatial_data(point_coords,geo_reference = geo)
        
        interp = Interpolate(vertices, triangles)
        f = num.array([linear_function(vertices),2*linear_function(vertices)])
        f = num.transpose(f)
        #print "f",f
        for blocking_max in range(14):
        #if True:
        #   blocking_max = 5
            z = interp.interpolate(f, point_coords,
                                   start_blocking_len=blocking_max)
            answer = [linear_function(point_coords.get_data_points( \
                      absolute = True)),
                      2*linear_function(point_coords.get_data_points( \
                      absolute = True)) ]
            answer = num.transpose(answer)
            #print "z",z
            #print "answer",answer
            assert num.allclose(z, answer)
            
        f = num.array([linear_function(vertices),2*linear_function(vertices),
                       2*linear_function(vertices) - 100])
        f = num.transpose(f)
        #print "f",f
        for blocking_max in range(14):
        #if True:
        #   blocking_max = 5
            z = interp.interpolate(f, point_coords,
                                   start_blocking_len=blocking_max)
            answer = num.array([linear_function(point_coords.get_data_points( \
                                                              absolute = True)),
                                                              2*linear_function(point_coords.get_data_points( \
                                                              absolute = True)) ,
                                                              2*linear_function(point_coords.get_data_points( \
                                                              absolute = True))-100])
            z = num.transpose(z)
            #print "z",z
            #print "answer",answer
            assert num.allclose(z, answer)

        #z = interp.interpolate(f, point_coords, start_blocking_len = 2)

        #print "z",z 
        #print "answer",answer 
        #assert num.allclose(z, answer)
        
    def test_interpolate_geo_spatial_2(self):
        a = [-1.0, 0.0]
        b = [3.0, 4.0]
        c = [4.0, 1.0]
        d = [-3.0, 2.0] #3
        e = [-1.0, -2.0]
        f = [1.0, -2.0] #5

        vertices = [a, b, c, d,e,f]
        triangles = [[0,1,3], [1,0,2], [0,4,5], [0,5,2]] #abd bac aef afc

        point_coords_absolute = [[-2.0,  2.0],
                                 [-1.0,  1.0],
                                 [ 0.0,  2.0],
                                 [ 1.0,  1.0],
                                 [ 2.0,  1.0],
                                 [ 0.0,  0.0],
                                 [ 1.0,  0.0],
                                 [ 0.0, -1.0],
                                 [-0.2, -0.5],
                                 [-0.9, -1.5],
                                 [ 0.5, -1.9],
                                 [ 3.0,  1.0]]

        geo = Geo_reference(57, 100, 500)
        point_coords = geo.change_points_geo_ref(point_coords_absolute)
        point_coords = Geospatial_data(point_coords, geo_reference=geo)
        
        interp = Interpolate(vertices, triangles)
        f = num.array([linear_function(vertices), 2*linear_function(vertices)])
        f = num.transpose(f)
        z = interp.interpolate_block(f, point_coords)
        answer = [linear_function(point_coords.get_data_points(absolute=True)),
                  2*linear_function(point_coords.get_data_points(absolute=True))
                 ]
        answer = num.transpose(answer)
        msg = ('Expected z\n%s\nto be close to answer\n%s'
               % (str(z), str(answer)))
        assert num.allclose(z, answer), msg
            
        z = interp.interpolate(f, point_coords, start_blocking_len = 2)

        msg = ('Expected z\n%s\nto be close to answer\n%s'
               % (str(z), str(answer)))
        assert num.allclose(z, answer)

        
    def test_interpolate_reuse_if_None(self):
        a = [-1.0, 0.0]
        b = [3.0, 4.0]
        c = [4.0, 1.0]
        d = [-3.0, 2.0] #3
        e = [-1.0, -2.0]
        f = [1.0, -2.0] #5

        vertices = [a, b, c, d,e,f]
        triangles = [[0,1,3], [1,0,2], [0,4,5], [0,5,2]] #abd bac aef afc


        point_coords = [[-2.0,  2.0],
                        [-1.0,  1.0],
                        [ 0.0,  2.0],
                        [ 1.0,  1.0],
                        [ 2.0,  1.0],
                        [ 0.0,  0.0],
                        [ 1.0,  0.0],
                        [ 0.0, -1.0],
                        [-0.2, -0.5],
                        [-0.9, -1.5],
                        [ 0.5, -1.9],
                        [ 3.0,  1.0]]

        interp = Interpolate(vertices, triangles)
        f = num.array([linear_function(vertices),2*linear_function(vertices)])
        f = num.transpose(f)
        z = interp.interpolate(f, point_coords,
                               start_blocking_len=20)
        answer = [linear_function(point_coords),
                  2*linear_function(point_coords) ]
        answer = num.transpose(answer)
        #print "z",z
        #print "answer",answer
        assert num.allclose(z, answer)
        assert num.allclose(interp._A_can_be_reused, True)

        z = interp.interpolate(f)
        assert num.allclose(z, answer)
        
        # This causes blocking to occur. 
        z = interp.interpolate(f, start_blocking_len=10)
        assert num.allclose(z, answer)
        assert num.allclose(interp._A_can_be_reused, False)

        #A is recalculated
        z = interp.interpolate(f)
        assert num.allclose(z, answer)
        assert num.allclose(interp._A_can_be_reused, True)
        
        interp = Interpolate(vertices, triangles)
        #Must raise an exception, no points specified
        try:
            z = interp.interpolate(f)
        except:
            pass
        
    def xxtest_interpolate_reuse_if_same(self):

        # This on tests that repeated identical interpolation
        # points makes use of precomputed matrix (Ole)
        # This is not really a test and is disabled for now
        
        a = [-1.0, 0.0]
        b = [3.0, 4.0]
        c = [4.0, 1.0]
        d = [-3.0, 2.0] #3
        e = [-1.0, -2.0]
        f = [1.0, -2.0] #5

        vertices = [a, b, c, d,e,f]
        triangles = [[0,1,3], [1,0,2], [0,4,5], [0,5,2]] #abd bac aef afc


        point_coords = [[-2.0,  2.0],
                        [-1.0,  1.0],
                        [ 0.0,  2.0],
                        [ 1.0,  1.0],
                        [ 2.0,  1.0],
                        [ 0.0,  0.0],
                        [ 1.0,  0.0],
                        [ 0.0, -1.0],
                        [-0.2, -0.5],
                        [-0.9, -1.5],
                        [ 0.5, -1.9],
                        [ 3.0,  1.0]]

        interp = Interpolate(vertices, triangles)
        f = num.array([linear_function(vertices), 2*linear_function(vertices)])
        f = num.transpose(f)
        z = interp.interpolate(f, point_coords)
        answer = [linear_function(point_coords),
                  2*linear_function(point_coords) ]
        answer = num.transpose(answer)

        assert num.allclose(z, answer)
        assert num.allclose(interp._A_can_be_reused, True)


        z = interp.interpolate(f)    # None
        assert num.allclose(z, answer)        
        z = interp.interpolate(f, point_coords) # Repeated (not really a test)        
        assert num.allclose(z, answer)
        


    def test_interpolation_interface_time_only(self):

        # Test spatio-temporal interpolation
        # Test that spatio temporal function performs the correct
        # interpolations in both time and space
        


        #Three timesteps
        time = [1.0, 5.0, 6.0]
        

        #One quantity
        Q = num.zeros( (3,6), num.float )

        #Linear in time and space
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        
        for i, t in enumerate(time):
            Q[i, :] = t*linear_function(points)

            
        #Check basic interpolation of one quantity using averaging
        #(no interpolation points or spatial info)
        I = Interpolation_function(time, [mean(Q[0,:]),
                                          mean(Q[1,:]),
                                          mean(Q[2,:])])



        #Check temporal interpolation
        for i in [0,1,2]:
            assert num.allclose(I(time[i]), mean(Q[i,:]))

        #Midway    
        assert num.allclose(I( (time[0] + time[1])/2 ),
                            (I(time[0]) + I(time[1]))/2 )

        assert num.allclose(I( (time[1] + time[2])/2 ),
                            (I(time[1]) + I(time[2]))/2 )

        assert num.allclose(I( (time[0] + time[2])/2 ),
                            (I(time[0]) + I(time[2]))/2 )                 

        #1/3
        assert num.allclose(I( (time[0] + time[2])/3 ),
                            (I(time[0]) + I(time[2]))/3 )                         


        #Out of bounds checks
        try:
            I(time[0]-1) 
        except:
            pass
        else:
            raise Exception('Should raise exception')

        try:
            I(time[-1]+1) 
        except:
            pass
        else:
            raise Exception('Should raise exception')


        

    def test_interpolation_interface_spatial_only(self):
        # Test spatio-temporal interpolation with constant time
        
        #Three timesteps
        time = [1.0, 5.0, 6.0]
               
        #Setup mesh used to represent fitted function
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        triangles = [[1,0,2], [1,2,4], [4,2,5], [3,1,4]]


        #New datapoints where interpolated values are sought
        interpolation_points = [[ 0.0, 0.0],
                                [ 0.5, 0.5],
                                [ 0.7, 0.7],
                                [ 1.0, 0.5],
                                [ 2.0, 0.4],
                                [ 2.8, 1.2]]


        #One quantity linear in space
        Q = linear_function(points)


        #Check interpolation of one quantity using interpolaton points
        I = Interpolation_function(time, Q,
                                   vertex_coordinates = points,
                                   triangles = triangles, 
                                   interpolation_points = interpolation_points,
                                   verbose = False)


        answer = linear_function(interpolation_points)

        t = time[0]
        for j in range(50): #t in [1, 6]
            for id in range(len(interpolation_points)):
                assert num.allclose(I(t, id), answer[id])
            t += 0.1    

        try:    
            I(1)
        except:
            pass
        else:
            raise Exception('Should raise exception')

            
    def test_interpolation_interface(self):
        # Test spatio-temporal interpolation
        # Test that spatio temporal function performs the correct
        # interpolations in both time and space
    
        #Three timesteps
        time = [1.0, 5.0, 6.0]    

        #Setup mesh used to represent fitted function
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        triangles = [[1,0,2], [1,2,4], [4,2,5], [3,1,4]]


        #New datapoints where interpolated values are sought
        interpolation_points = [[ 0.0, 0.0],
                                [ 0.5, 0.5],
                                [ 0.7, 0.7],
                                [ 1.0, 0.5],
                                [ 2.0, 0.4],
                                [ 2.8, 1.2]]

        #One quantity
        Q = num.zeros( (3,6), num.float )

        #Linear in time and space
        for i, t in enumerate(time):
            Q[i, :] = t*linear_function(points)

        #Check interpolation of one quantity using interpolaton points)
        I = Interpolation_function(time, Q,
                                   vertex_coordinates = points,
                                   triangles = triangles, 
                                   interpolation_points = interpolation_points,
                                   verbose = False)

        answer = linear_function(interpolation_points)
        
        t = time[0]
        for j in range(50): #t in [1, 6]
            for id in range(len(interpolation_points)):
                assert num.allclose(I(t, id), t*answer[id])
            t += 0.1    

        try:    
            I(1)
        except:
            pass
        else:
            raise Exception('Should raise exception')



    def test_interpolation_interface_with_time_thinning(self):
        # Test spatio-temporal interpolation
        # Test that spatio temporal function performs the correct
        # interpolations in both time and space
    
        # Eight timesteps
        time = [1.0, 2.0, 4.0, 5.0, 7.0, 8.0, 9.0, 10.0]    

        # Setup mesh used to represent fitted function
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        # bac, bce, ecf, dbe
        triangles = [[1,0,2], [1,2,4], [4,2,5], [3,1,4]]


        # New datapoints where interpolated values are sought
        interpolation_points = [[ 0.0, 0.0],
                                [ 0.5, 0.5],
                                [ 0.7, 0.7],
                                [ 1.0, 0.5],
                                [ 2.0, 0.4],
                                [ 2.8, 1.2]]

        # One quantity
        Q = num.zeros((8,6), num.float)

        # Linear in time and space
        for i, t in enumerate(time):
            Q[i, :] = t*linear_function(points)

        # Check interpolation of one quantity using interpolaton points) using default
        # time_thinning of 1
        I = Interpolation_function(time, Q,
                                   vertex_coordinates=points,
                                   triangles=triangles, 
                                   interpolation_points=interpolation_points,
                                   verbose=False)

        answer = linear_function(interpolation_points)

        
        t = time[0]
        for j in range(90): #t in [1, 10]
            for id in range(len(interpolation_points)):
                assert num.allclose(I(t, id), t*answer[id])
            t += 0.1    


        # Now check time_thinning
        I = Interpolation_function(time, Q,
                                   vertex_coordinates=points,
                                   triangles=triangles, 
                                   interpolation_points=interpolation_points,
                                   time_thinning=2,
                                   verbose=False)


        assert len(I.time) == 4
        assert( num.allclose(I.time, [1.0, 4.0, 7.0, 9.0] ))    

        answer = linear_function(interpolation_points)

        t = time[0]
        for j in range(80): #t in [1, 9]
            for id in range(len(interpolation_points)):
                assert num.allclose(I(t, id), t*answer[id])
            t += 0.1    



    def test_interpolation_precompute_points(self):
        # looking at a discrete mesh
        #
    
        #Three timesteps
        time = [0.0, 60.0]    

        #Setup mesh used to represent fitted function
        points = [[ 15., -20.],
                  [ 15.,  10.],
                  [  0., -20.],
                  [  0.,  10.],
                  [  0., -20.],
                  [ 15.,  10.]]
        
        triangles = [[0, 1, 2],
                     [3, 4, 5]]

        #New datapoints where interpolated values are sought
        interpolation_points = [[ 1.,  0.], [0.,1.]]

        #One quantity
        Q = num.zeros( (2,6), num.float )

        #Linear in time and space
        for i, t in enumerate(time):
            Q[i, :] = t*linear_function(points)
        #print "Q", Q


        
        interp = Interpolate(points, triangles)
        f = num.array([linear_function(points),2*linear_function(points)])
        f = num.transpose(f)
        #print "f",f
        z = interp.interpolate(f, interpolation_points)
        answer = [linear_function(interpolation_points),
                  2*linear_function(interpolation_points) ]
        answer = num.transpose(answer)
        #print "z",z
        #print "answer",answer
        assert num.allclose(z, answer)


        #Check interpolation of one quantity using interpolaton points)
        I = Interpolation_function(time, Q,
                                   vertex_coordinates = points,
                                   triangles = triangles, 
                                   interpolation_points = interpolation_points,
                                   verbose = False)
        
        #print "I.precomputed_values", I.precomputed_values

        msg = 'Interpolation failed'
        assert num.allclose(I.precomputed_values['Attribute'][1], [60, 60]), msg
        #self.assertTrue( I.precomputed_values['Attribute'][1] == 60.0,
        #                ' failed')
        
    def test_interpolation_function_outside_point(self):
        # Test spatio-temporal interpolation
        # Test that spatio temporal function performs the correct
        # interpolations in both time and space
    
        # Three timesteps
        time = [1.0, 5.0, 6.0]    

        # Setup mesh used to represent fitted function
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        triangles = [[1,0,2], [1,2,4], [4,2,5], [3,1,4]]


        # New datapoints where interpolated values are sought
        interpolation_points = [[ 0.0, 0.0],
                                [ 0.5, 0.5],
                                [ 0.7, 0.7],
                                [ 1.0, 0.5],
                                [ 2.0, 0.4],
                                [ 545354534, 4354354353]] # outside the mesh

        # One quantity
        Q = num.zeros( (3,6), num.float )

        # Linear in time and space
        for i, t in enumerate(time):
            Q[i, :] = t*linear_function(points)

        # Check interpolation of one quantity using interpolaton points)

        I = Interpolation_function(time, Q,
                                   vertex_coordinates = points,
                                   triangles = triangles, 
                                   interpolation_points = interpolation_points,
                                   verbose = False)
        
        
        assert num.alltrue(I.precomputed_values['Attribute'][:,4] != NAN)
        assert num.sometrue(I.precomputed_values['Attribute'][:,5] == NAN)

        #X = I.precomputed_values['Attribute'][1,:]
        #print X
        #print take(X, X == NAN)
        #print where(X == NAN, range(len(X)), 0)        
        
        answer = linear_function(interpolation_points)
          
        t = time[0]
        for j in range(50): #t in [1, 6]
            for id in range(len(interpolation_points)-1):
                assert num.allclose(I(t, id), t*answer[id])
            t += 0.1
            
        # Now test the point outside the mesh
        t = time[0]
        for j in range(50): #t in [1, 6]
            self.assertTrue(I(t, 5) == NAN, 'Fail!')
            t += 0.1  
            
        try:    
            I(1)
        except:
            pass
        else:
            raise Exception('Should raise exception')


    def test_interpolation_function_time(self):
        #Test a long time series with an error in it (this did cause an
        #error once)
        

        time = num.array(\
        [0.00000000e+00, 5.00000000e-02, 1.00000000e-01,   1.50000000e-01,
        2.00000000e-01,   2.50000000e-01,   3.00000000e-01,   3.50000000e-01,
        4.00000000e-01,   4.50000000e-01,   5.00000000e-01,   5.50000000e-01,
        6.00000000e-01,   6.50000000e-01,   7.00000000e-01,   7.50000000e-01,
        8.00000000e-01,   8.50000000e-01,   9.00000000e-01,   9.50000000e-01,
        1.00000000e-00,   1.05000000e+00,   1.10000000e+00,   1.15000000e+00,
        1.20000000e+00,   1.25000000e+00,   1.30000000e+00,   1.35000000e+00,
        1.40000000e+00,   1.45000000e+00,   1.50000000e+00,   1.55000000e+00,
        1.60000000e+00,   1.65000000e+00,   1.70000000e+00,   1.75000000e+00,
        1.80000000e+00,   1.85000000e+00,   1.90000000e+00,   1.95000000e+00,
        2.00000000e+00,   2.05000000e+00,   2.10000000e+00,   2.15000000e+00,
        2.20000000e+00,   2.25000000e+00,   2.30000000e+00,   2.35000000e+00,
        2.40000000e+00,   2.45000000e+00,   2.50000000e+00,   2.55000000e+00,
        2.60000000e+00,   2.65000000e+00,   2.70000000e+00,   2.75000000e+00,
        2.80000000e+00,   2.85000000e+00,   2.90000000e+00,   2.95000000e+00,
        3.00000000e+00,   3.05000000e+00,   9.96920997e+36,   3.15000000e+00,
        3.20000000e+00,   3.25000000e+00,   3.30000000e+00,   3.35000000e+00,
        3.40000000e+00,   3.45000000e+00,   3.50000000e+00,   3.55000000e+00,
        3.60000000e+00,   3.65000000e+00,   3.70000000e+00,   3.75000000e+00,
        3.80000000e+00,   3.85000000e+00,   3.90000000e+00,   3.95000000e+00,
        4.00000000e+00,   4.05000000e+00,   4.10000000e+00,   4.15000000e+00,
        4.20000000e+00,   4.25000000e+00,   4.30000000e+00,   4.35000000e+00,
        4.40000000e+00,   4.45000000e+00,   4.50000000e+00,   4.55000000e+00,
        4.60000000e+00,   4.65000000e+00,   4.70000000e+00,   4.75000000e+00,
        4.80000000e+00,   4.85000000e+00,   4.90000000e+00,   4.95000000e+00,
        5.00000000e+00,   5.05000000e+00,   5.10000000e+00,   5.15000000e+00,
        5.20000000e+00,   5.25000000e+00,   5.30000000e+00,   5.35000000e+00,
        5.40000000e+00,   5.45000000e+00,   5.50000000e+00,   5.55000000e+00,
        5.60000000e+00,   5.65000000e+00,   5.70000000e+00,   5.75000000e+00,
        5.80000000e+00,   5.85000000e+00,   5.90000000e+00,   5.95000000e+00,
        6.00000000e+00,   6.05000000e+00,   6.10000000e+00,   6.15000000e+00,
        6.20000000e+00,   6.25000000e+00,   6.30000000e+00,   6.35000000e+00,
        6.40000000e+00,   6.45000000e+00,   6.50000000e+00,   6.55000000e+00,
        6.60000000e+00,   6.65000000e+00,   6.70000000e+00,   6.75000000e+00,
        6.80000000e+00,   6.85000000e+00,   6.90000000e+00,   6.95000000e+00,
        7.00000000e+00,   7.05000000e+00,   7.10000000e+00,   7.15000000e+00,
        7.20000000e+00,   7.25000000e+00,   7.30000000e+00,   7.35000000e+00,
        7.40000000e+00,   7.45000000e+00,   7.50000000e+00,   7.55000000e+00,
        7.60000000e+00,   7.65000000e+00,   7.70000000e+00,   7.75000000e+00,
        7.80000000e+00,   7.85000000e+00,   7.90000000e+00,   7.95000000e+00,
        8.00000000e+00,   8.05000000e+00,   8.10000000e+00,   8.15000000e+00,
        8.20000000e+00,   8.25000000e+00,   8.30000000e+00,   8.35000000e+00,
        8.40000000e+00,   8.45000000e+00,   8.50000000e+00,   8.55000000e+00,
        8.60000000e+00,   8.65000000e+00,   8.70000000e+00,   8.75000000e+00,
        8.80000000e+00,   8.85000000e+00,   8.90000000e+00,   8.95000000e+00,
        9.00000000e+00,   9.05000000e+00,   9.10000000e+00,   9.15000000e+00,
        9.20000000e+00,   9.25000000e+00,   9.30000000e+00,   9.35000000e+00,
        9.40000000e+00,   9.45000000e+00,   9.50000000e+00,   9.55000000e+00,
        9.60000000e+00,   9.65000000e+00,   9.70000000e+00,   9.75000000e+00,
        9.80000000e+00,   9.85000000e+00,   9.90000000e+00,   9.95000000e+00,
        1.00000000e+01,   1.00500000e+01,   1.01000000e+01,   1.01500000e+01,
        1.02000000e+01,   1.02500000e+01,   1.03000000e+01,   1.03500000e+01,
        1.04000000e+01,   1.04500000e+01,   1.05000000e+01,   1.05500000e+01,
        1.06000000e+01,   1.06500000e+01,   1.07000000e+01,   1.07500000e+01,
        1.08000000e+01,   1.08500000e+01,   1.09000000e+01,   1.09500000e+01,
        1.10000000e+01,   1.10500000e+01,   1.11000000e+01,   1.11500000e+01,
        1.12000000e+01,   1.12500000e+01,   1.13000000e+01,   1.13500000e+01,
        1.14000000e+01,   1.14500000e+01,   1.15000000e+01,   1.15500000e+01,
        1.16000000e+01,   1.16500000e+01,   1.17000000e+01,   1.17500000e+01,
        1.18000000e+01,   1.18500000e+01,   1.19000000e+01,   1.19500000e+01,
        1.20000000e+01,   1.20500000e+01,   1.21000000e+01,   1.21500000e+01,
        1.22000000e+01,   1.22500000e+01,   1.23000000e+01,   1.23500000e+01,
        1.24000000e+01,   1.24500000e+01,   1.25000000e+01,   1.25500000e+01,
        1.26000000e+01,   1.26500000e+01,   1.27000000e+01,   1.27500000e+01,
        1.28000000e+01,   1.28500000e+01,   1.29000000e+01,   1.29500000e+01,
        1.30000000e+01,   1.30500000e+01,   1.31000000e+01,   1.31500000e+01,
        1.32000000e+01,   1.32500000e+01,   1.33000000e+01,   1.33500000e+01,
        1.34000000e+01,   1.34500000e+01,   1.35000000e+01,   1.35500000e+01,
        1.36000000e+01,   1.36500000e+01,   1.37000000e+01,   1.37500000e+01,
        1.38000000e+01,   1.38500000e+01,   1.39000000e+01,   1.39500000e+01,
        1.40000000e+01,   1.40500000e+01,   1.41000000e+01,   1.41500000e+01,
        1.42000000e+01,   1.42500000e+01,   1.43000000e+01,   1.43500000e+01,
        1.44000000e+01,   1.44500000e+01,   1.45000000e+01,   1.45500000e+01,
        1.46000000e+01,   1.46500000e+01,   1.47000000e+01,   1.47500000e+01,
        1.48000000e+01,   1.48500000e+01,   1.49000000e+01,   1.49500000e+01,
        1.50000000e+01,   1.50500000e+01,   1.51000000e+01,   1.51500000e+01,
        1.52000000e+01,   1.52500000e+01,   1.53000000e+01,   1.53500000e+01,
        1.54000000e+01,   1.54500000e+01,   1.55000000e+01,   1.55500000e+01,
        1.56000000e+01,   1.56500000e+01,   1.57000000e+01,   1.57500000e+01,
        1.58000000e+01,   1.58500000e+01,   1.59000000e+01,   1.59500000e+01,
        1.60000000e+01,   1.60500000e+01,   1.61000000e+01,   1.61500000e+01,
        1.62000000e+01,   1.62500000e+01,   1.63000000e+01,   1.63500000e+01,
        1.64000000e+01,   1.64500000e+01,   1.65000000e+01,   1.65500000e+01,
        1.66000000e+01,   1.66500000e+01,   1.67000000e+01,   1.67500000e+01,
        1.68000000e+01,   1.68500000e+01,   1.69000000e+01,   1.69500000e+01,
        1.70000000e+01,   1.70500000e+01,   1.71000000e+01,   1.71500000e+01,
        1.72000000e+01,   1.72500000e+01,   1.73000000e+01,   1.73500000e+01,
        1.74000000e+01,   1.74500000e+01,   1.75000000e+01,   1.75500000e+01,
        1.76000000e+01,   1.76500000e+01,   1.77000000e+01,   1.77500000e+01,
        1.78000000e+01,   1.78500000e+01,   1.79000000e+01,   1.79500000e+01,
        1.80000000e+01,   1.80500000e+01,   1.81000000e+01,   1.81500000e+01,
        1.82000000e+01,   1.82500000e+01,   1.83000000e+01,   1.83500000e+01,
        1.84000000e+01,   1.84500000e+01,   1.85000000e+01,   1.85500000e+01,
        1.86000000e+01,   1.86500000e+01,   1.87000000e+01,   1.87500000e+01,
        1.88000000e+01,   1.88500000e+01,   1.89000000e+01,   1.89500000e+01,
        1.90000000e+01,   1.90500000e+01,   1.91000000e+01,   1.91500000e+01,
        1.92000000e+01,   1.92500000e+01,   1.93000000e+01,   1.93500000e+01,
        1.94000000e+01,   1.94500000e+01,   1.95000000e+01,   1.95500000e+01,
        1.96000000e+01,   1.96500000e+01,   1.97000000e+01,   1.97500000e+01,
        1.98000000e+01,   1.98500000e+01,   1.99000000e+01,   1.99500000e+01,
        2.00000000e+01,   2.00500000e+01,   2.01000000e+01,   2.01500000e+01,
        2.02000000e+01,   2.02500000e+01,   2.03000000e+01,   2.03500000e+01,
        2.04000000e+01,   2.04500000e+01,   2.05000000e+01,   2.05500000e+01,
        2.06000000e+01,   2.06500000e+01,   2.07000000e+01,   2.07500000e+01,
        2.08000000e+01,   2.08500000e+01,   2.09000000e+01,   2.09500000e+01,
        2.10000000e+01,   2.10500000e+01,   2.11000000e+01,   2.11500000e+01,
        2.12000000e+01,   2.12500000e+01,   2.13000000e+01,   2.13500000e+01,
        2.14000000e+01,   2.14500000e+01,   2.15000000e+01,   2.15500000e+01,
        2.16000000e+01,   2.16500000e+01,   2.17000000e+01,   2.17500000e+01,
        2.18000000e+01,   2.18500000e+01,   2.19000000e+01,   2.19500000e+01,
        2.20000000e+01,   2.20500000e+01,   2.21000000e+01,   2.21500000e+01,
        2.22000000e+01,   2.22500000e+01,   2.23000000e+01,   2.23500000e+01,
        2.24000000e+01,   2.24500000e+01,   2.25000000e+01])

        #print 'Diff', time[1:] - time[:-1] 

        #Setup mesh used to represent fitted function
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        triangles = [[1,0,2], [1,2,4], [4,2,5], [3,1,4]]


        #New datapoints where interpolated values are sought
        interpolation_points = [[ 0.0, 0.0],
                                [ 0.5, 0.5],
                                [ 0.7, 0.7],
                                [ 1.0, 0.5],
                                [ 2.0, 0.4],
                                [ 545354534, 4354354353]] # outside the mesh

        #One quantity
        Q = num.zeros( (len(time),6), num.float )

        #Linear in time and space
        for i, t in enumerate(time):
            Q[i, :] = t*linear_function(points)

        #Check interpolation of one quantity using interpolaton points)
        try:
            I = Interpolation_function(time, Q,
                                       vertex_coordinates = points,
                                       triangles = triangles, 
                                       interpolation_points = interpolation_points,
                                       verbose = False)
        except:
            pass
        else:
            raise Exception('Should raise exception due to time being non-monotoneous')
      

    def test_points_outside_the_polygon(self):
        a = [-1.0, 0.0]
        b = [3.0, 4.0]
        c = [4.0, 1.0]
        d = [-3.0, 2.0] #3
        e = [-1.0, -2.0]
        f = [1.0, -2.0] #5

        vertices = [a, b, c, d,e,f]
        triangles = [[0,1,3], [1,0,2], [0,4,5], [0,5,2]] #abd bac aef afc

        point_coords = [[-2.0, 2.0],
                        [-1.0, 1.0],
                        [9999.0, 9999.0], # point Outside poly
                        [-9999.0, 1.0], # point Outside poly
                        [2.0, 1.0],
                        [0.0, 0.0],
                        [1.0, 0.0],
                        [0.0, -1.0],
                        [-0.2, -0.5],
                        [-0.9, -1.5],
                        [0.5, -1.9],
                        [999999, 9999999]] # point Outside poly
        geo_data = Geospatial_data(data_points = point_coords)

        interp = Interpolate(vertices, triangles)
        f = num.array([linear_function(vertices),2*linear_function(vertices)])
        f = num.transpose(f)
        #print "f",f
        z = interp.interpolate(f, geo_data)
        #z = interp.interpolate(f, point_coords)
        answer = [linear_function(point_coords),
                  2*linear_function(point_coords) ]
        answer = num.transpose(answer)
        answer[2,:] = [NAN, NAN]
        answer[3,:] = [NAN, NAN]
        answer[11,:] = [NAN, NAN]
        #print "z",z
        #print "answer _ fixed",answer
        assert num.allclose(z[0:1], answer[0:1])
        assert num.allclose(z[4:10], answer[4:10])
        for i in [2,3,11]:
            self.assertTrue( z[i,1] == answer[11,1], 'Fail!')
            self.assertTrue( z[i,0] == answer[11,0], 'Fail!')




    def test_points_in_hole(self):

        v0  = [0.0,     0.0]
        v1  = [1.0/3.0, 0.0]
        v2  = [2.0/3.0, 0.0]
        v3  = [1.0,     0.0]
        v4  = [0.0,     1.0/3.0]
        v5  = [1.0/3.0, 1.0/3.0]
        v6  = [2.0/3.0, 1.0/3.0]
        v7  = [1.0,     1.0/3.0]
        v8  = [0.0,     2.0/3.0]
        v9  = [1.0/3.0, 2.0/3.0]
        v10 = [2.0/3.0, 2.0/3.0]
        v11 = [1.0,     2.0/3.0]
        v12 = [0.0,     1.0]
        v13 = [1.0/3.0, 1.0]
        v14 = [2.0/3.0, 1.0]
        v15 = [1.0,     1.0]


        vertices = [v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15]

        triangles = [[1,4,0], [1,5,4],    [2,5,1], [2,6,5],       [3,6,2],  [3,7,6],
                     [5,8,4], [5,9,8],                            [7,10,6], [7,11,10],
                     [9,12,8], [9,13,12], [10,13,9], [10,14,13], [11,14,10], [11,15,14]]




        point_coords = [[0.25, 0.25],
                        [0.25, 0.2],
                        [10.0, 12.0], # point Outside poly
                        [0.5, 0.5], # point in hole
                        [0.5, 0.4], # point in hole
                        [0.0, 0.0],
                        [1.0, 0.0]]

        geo_data = Geospatial_data(data_points = point_coords)

        interp = Interpolate(vertices, triangles)
        f = num.array([linear_function(vertices),2*linear_function(vertices)])
        f = num.transpose(f)
        #print "f",f
        z = interp.interpolate(f, geo_data)
        #z = interp.interpolate(f, point_coords)
        answer = [linear_function(point_coords),
                  2*linear_function(point_coords) ]
        answer = num.transpose(answer)
        answer[2,:] = [NAN, NAN]
        answer[3,:] = [NAN, NAN]
        answer[4,:] = [NAN, NAN]
        #print "z",z
        #print "answer _ fixed",answer
        assert num.allclose(z[0:1], answer[0:1])
        assert num.allclose(z[5:6], answer[5:6])
        for i in [2,3,4]:
            self.assertTrue( z[i,1] == answer[2,1], 'Fail!')
            self.assertTrue( z[i,0] == answer[2,0], 'Fail!')


    def test_interpolate_sww2csv(self):

        def elevation_function(x, y):
            return -x
        
        # Create mesh
        mesh_file = tempfile.mktemp(".tsh")    
        points = [[0.0,0.0],[6.0,0.0],[6.0,6.0],[0.0,6.0]]
        m = Mesh()
        m.add_vertices(points)
        m.auto_segment()
        m.generate_mesh(verbose=False)
        m.export_mesh_file(mesh_file)
        
        # Create shallow water domain
        domain = Domain(mesh_file)
        os.remove(mesh_file)
        
        domain.default_order = 2

        # This test was made before tight_slope_limiters were introduced
        # Since were are testing interpolation values this is OK
        domain.tight_slope_limiters = 0 

        # Set some field values
        domain.set_quantity('elevation', elevation_function)
        domain.set_quantity('friction', 0.03)
        domain.set_quantity('xmomentum', 3.0)
        domain.set_quantity('ymomentum', 4.0)

        ######################
        # Boundary conditions
        B = anuga.Transmissive_boundary(domain)
        domain.set_boundary( {'exterior': B})

        # This call mangles the stage values.
        domain.distribute_to_vertices_and_edges()
        domain.set_quantity('stage', 1.0)


        domain.set_name('datatest' + str(time.time()))
        domain.smooth = True
        domain.reduction = mean

        sww = SWW_file(domain)
        sww.store_connectivity()
        sww.store_timestep()
        domain.set_quantity('stage', 10.0) # This is automatically limited
        # So it will not be less than the elevation
        domain.time = 2.
        sww.store_timestep()

        # Test the function
        points = [[5.0,1.],[0.5,2.]]
        depth_file = tempfile.mktemp(".csv") 
        velocity_x_file = tempfile.mktemp(".csv") 
        velocity_y_file = tempfile.mktemp(".csv") 
        interpolate_sww2csv(sww.filename, points, depth_file,
                            velocity_x_file,
                            velocity_y_file,
                            verbose=False)

        depth_answers_array = [[0.0, 6.0, 1.5], [2.0, 15., 10.5]] 
        velocity_x_answers_array = [[0.0, 3./6.0, 3./1.5],
                                    [2.0, 3./15., 3/10.5]]
        velocity_y_answers_array = [[0.0, 4./6.0, 4./1.5],
                                    [2.0, 4./15., 4./10.5]]
        depth_file_handle = file(depth_file)
        depth_reader = csv.reader(depth_file_handle)
        depth_reader.next()
        velocity_x_file_handle = file(velocity_x_file)
        velocity_x_reader = csv.reader(velocity_x_file_handle)
        velocity_x_reader.next()
        for depths, velocitys, depth_answers, velocity_answers in map(None,
                                              depth_reader,
                                              velocity_x_reader,
                                              depth_answers_array,
                                              velocity_x_answers_array):
            for i in range(len(depths)):
                #print "depths",depths[i] 
                #print "depth_answers",depth_answers[i]
                #print "velocitys",velocitys[i] 
                #print "velocity_answers_array", velocity_answers[i]
                msg = 'Interpolation failed'
                assert num.allclose(float(depths[i]), depth_answers[i]), msg
                assert num.allclose(float(velocitys[i]), velocity_answers[i]), msg

        velocity_y_file_handle = file(velocity_y_file)
        velocity_y_reader = csv.reader(velocity_y_file_handle)
        velocity_y_reader.next()
        for velocitys, velocity_answers in map(None,
                                              velocity_y_reader,
                                              velocity_y_answers_array):
            for i in range(len(depths)):
                #print "depths",depths[i] 
                #print "depth_answers",depth_answers[i]
                #print "velocitys",velocitys[i] 
                #print "velocity_answers_array", velocity_answers[i]
                msg = 'Interpolation failed'
                assert num.allclose(float(depths[i]), depth_answers[i]), msg
                assert num.allclose(float(velocitys[i]), velocity_answers[i]), msg
                
        # clean up
        depth_file_handle.close()
        velocity_y_file_handle.close()
        velocity_x_file_handle.close()
        #print "sww.filename",sww.filename 
        os.remove(sww.filename)
        os.remove(depth_file)
        os.remove(velocity_x_file)
        os.remove(velocity_y_file)

        
    def test_interpolate_one_point_many_triangles(self):
        z0 = [2.0, 5.0]

        v0 = [0.0, 0.0]
        v1 = [1.0, 0.0]
        v2 = [2.0, 0.0]
        v3 = [3.0, 0.0]
        v4 = [4.0, 0.0]
        v5 = [5.0, 0.0]
        v6 = [6.0, 0.0]
        v7 = [0.0, 10.0]
        v8 = [1.0, 10.0]
        v9 = [2.0, 10.0]
        v10= [3.0, 10.0]
        v11= [4.0, 10.0]
        v12= [5.0, 10.0]
        v13= [6.0, 10.0]
        
        vertices = [z0,v0, v1, v2, v3,v4 ,v5, v6, v7, v8, v9, v10, v11,
                    v12, v13]
        triangles = [
                      [0,1,2],
                      [0,2,3],
                      [0,3,4],
                      [0,4,5],
                      [0,5,6],
                      [0,6,7],
                      [0,9,8],
                      [0,10,9],
                      [0,11,10],
                      [0,12,11],
                      [0,13,12],
                      [0,14,13]
                      ]

        d0 = [1.0, 1.0]
        d1 = [1.0, 2.0]
        d2 = [3.0, 1.0]
        point_coords = [ d0, d1, d2]
        try:
            interp = Interpolate(vertices, triangles)
        except RuntimeError:
            self.assertTrue(0 ==1,  'quad fails with 14 verts at the same \
            position. Should be able to handle any number.')
        f = linear_function(vertices)
        z = interp.interpolate(f, point_coords)
        answer = linear_function(point_coords)

        #print "z",z 
        #print "answer",answer 
        assert num.allclose(z, answer)

################################################################################

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_Interpolate,'test')
    runner = unittest.TextTestRunner() #verbosity=1)
    runner.run(suite)


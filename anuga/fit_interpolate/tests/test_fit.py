#!/usr/bin/env python

"""
Note, fitting/blocking is also tested in
test_quantity.test_set_values_from_UTM_pts.
"""

#TEST

import unittest
from math import sqrt
import tempfile
import os

from anuga.fit_interpolate.fit import *
from anuga.abstract_2d_finite_volumes.neighbour_mesh import Mesh
from anuga.utilities.sparse import Sparse, Sparse_CSR
from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.utilities.numerical_tools import ensure_numeric
from anuga.geospatial_data.geospatial_data import Geospatial_data
from anuga.shallow_water.shallow_water_domain import Domain

import numpy as num


def distance(x, y):
    return sqrt(num.sum((num.array(x)-num.array(y))**2))

def linear_function(point):
    point = num.array(point)
    return point[:,0]+point[:,1]


class Test_Fit(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass


    def test_smooth_attributes_to_mesh(self):
        a = [0.0, 0.0]
        b = [0.0, 5.0]
        c = [5.0, 0.0]
        points = [a, b, c]
        triangles = [ [1,0,2] ]   #bac

        d1 = [1.0, 1.0]
        d2 = [1.0, 3.0]
        d3 = [3.0,1.0]
        z1 = 2
        z2 = 4
        z3 = 4
        data_coords = [d1, d2, d3]

        z = [z1, z2, z3]
        fit = Fit(points, triangles, alpha=0)
        #print "interp.get_A()", interp.get_A()
        fit._build_matrix_AtA_Atz(ensure_numeric(data_coords),
                                  ensure_numeric(z))
        #print "Atz - from fit", fit.Atz
        #print "AtA - from fit", fit.AtA.todense()
        #print "z",z

        assert num.allclose(fit.Atz, [2.8, 3.6, 3.6], atol=1e-7)

        f = fit.fit()

        answer = [0, 5., 5.]

        #print "f\n",f
        #print "answer\n",answer

        assert num.allclose(f, answer, atol=1e-7)

    def test_smooth_att_to_meshII(self):

        a = [0.0, 0.0]
        b = [0.0, 5.0]
        c = [5.0, 0.0]
        points = [a, b, c]
        triangles = [ [1,0,2] ]   #bac

        d1 = [1.0, 1.0]
        d2 = [1.0, 2.0]
        d3 = [3.0,1.0]
        data_coords = [d1, d2, d3]
        z = linear_function(data_coords)
        #print "z",z

        interp = Fit(points, triangles, alpha=0.0)
        f = interp.fit(data_coords, z)
        answer = linear_function(points)
        #print "f\n",f
        #print "answer\n",answer

        assert num.allclose(f, answer)

    def test_smooth_attributes_to_meshIII(self):

        a = [-1.0, 0.0]
        b = [3.0, 4.0]
        c = [4.0,1.0]
        d = [-3.0, 2.0] #3
        e = [-1.0,-2.0]
        f = [1.0, -2.0] #5

        vertices = [a, b, c, d,e,f]
        triangles = [[0,1,3], [1,0,2], [0,4,5], [0,5,2]] #abd bac aef afc

        point_coords = [[-2.0, 2.0],
                        [-1.0, 1.0],
                        [0.0,2.0],
                        [1.0, 1.0],
                        [2.0, 1.0],
                        [0.0,0.0],
                        [1.0, 0.0],
                        [0.0, -1.0],
                        [-0.2,-0.5],
                        [-0.9, -1.5],
                        [0.5, -1.9],
                        [3.0,1.0]]

        z = linear_function(point_coords)
        interp = Fit(vertices, triangles,
                                alpha=0.0)

        #print 'z',z
        f = interp.fit(point_coords,z)
        answer = linear_function(vertices)
        #print "f\n",f
        #print "answer\n",answer
        assert num.allclose(f, answer)


    def test_smooth_attributes_to_meshIV(self):
        # Testing 2 attributes smoothed to the mesh


        a = [0.0, 0.0]
        b = [0.0, 5.0]
        c = [5.0, 0.0]
        points = [a, b, c]
        triangles = [ [1,0,2] ]   #bac

        d1 = [1.0, 1.0]
        d2 = [1.0, 3.0]
        d3 = [3.0, 1.0]
        z1 = [2, 4]
        z2 = [4, 8]
        z3 = [4, 8]
        data_coords = [d1, d2, d3]

        z = [z1, z2, z3]
        fit = Fit(points, triangles, alpha=0)

        f =  fit.fit(data_coords,z)
        answer = [[0,0], [5., 10.], [5., 10.]]
        assert num.allclose(f, answer)

    def test_smooth_attributes_to_mesh_build_fit_subset(self):

        a = [-1.0, 0.0]
        b = [3.0, 4.0]
        c = [4.0,1.0]
        d = [-3.0, 2.0] #3
        e = [-1.0,-2.0]
        f = [1.0, -2.0] #5

        vertices = [a, b, c, d,e,f]
        triangles = [[0,1,3], [1,0,2], [0,4,5], [0,5,2]] #abd bac aef afc

        interp = Fit(vertices, triangles,
                                alpha=0.0)

        point_coords = [[-2.0, 2.0],
                        [-1.0, 1.0],
                        [0.0,2.0],
                        [1.0, 1.0],
                       ]

        z = linear_function(point_coords)

        f = interp.build_fit_subset(point_coords,z)

        point_coords = [
                        [2.0, 1.0],
                        [0.0,0.0],
                        [1.0, 0.0],
                        [0.0, -1.0],
                        [-0.2,-0.5],
                        [-0.9, -1.5],
                        [0.5, -1.9],
                        [3.0,1.0]]

        z = linear_function(point_coords)

        f = interp.build_fit_subset(point_coords,z)

        #print 'z',z
        f = interp.fit()
        answer = linear_function(vertices)
        #print "f\n",f
        #print "answer\n",answer
        assert num.allclose(f, answer)

        # test fit 2 mesh as well.
    def test_fit_file_blocking(self):

        a = [-1.0, 0.0]
        b = [3.0, 4.0]
        c = [4.0,1.0]
        d = [-3.0, 2.0] #3
        e = [-1.0,-2.0]
        f = [1.0, -2.0] #5

        vertices = [a, b, c, d,e,f]
        triangles = [[0,1,3], [1,0,2], [0,4,5], [0,5,2]] #abd bac aef afc

        interp = Fit(vertices, triangles,
                     alpha=0.0)


        fileName = tempfile.mktemp(".ddd")
        file = open(fileName,"w")
        file.write(" x,y, elevation \n\
-2.0, 2.0, 0.\n\
-1.0, 1.0, 0.\n\
0.0, 2.0 , 2.\n\
1.0, 1.0 , 2.\n\
 2.0,  1.0 ,3. \n\
 0.0,  0.0 , 0.\n\
 1.0,  0.0 , 1.\n\
 0.0,  -1.0, -1.\n\
 -0.2, -0.5, -0.7\n\
 -0.9, -1.5, -2.4\n\
 0.5,  -1.9, -1.4\n\
 3.0,  1.0 , 4.\n")
        file.close()

        f = interp.fit(fileName, max_read_lines=2)
        answer = linear_function(vertices)
        #print "f\n",f
        #print "answer\n",answer
        assert num.allclose(f, answer)
        os.remove(fileName)

    def test_fit_to_mesh_UTM_file(self):
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
        for data_point, attribute in zip(points_UTM, attributes):
            row = str(data_point[0]) + ',' + str(data_point[1]) \
                  + ',' + str(attribute)
            #print "row", row
            file.write(row + "\n")
        file.close()

        # setting up the mesh
        a = [240000, 7620000]
        b = [240000, 7680000]
        c = [300000, 7620000]
        points = [a, b, c]
        elements = [[0,2,1]]
        f = fit_to_mesh(txt_file, points, elements,
                        alpha=0.0, max_read_lines=2)
        answer = linear_function(points)
        #print "f",f
        #print "answer",answer
        assert num.allclose(f, answer)

        # Delete file!
        os.remove(txt_file)


    def test_fit_to_mesh_pts(self):
        a = [-1.0, 0.0]
        b = [3.0, 4.0]
        c = [4.0,1.0]
        d = [-3.0, 2.0] #3
        e = [-1.0,-2.0]
        f = [1.0, -2.0] #5

        vertices = [a, b, c, d,e,f]
        triangles = [[0,1,3], [1,0,2], [0,4,5], [0,5,2]] #abd bac aef afc


        fileName = tempfile.mktemp(".txt")
        file = open(fileName,"w")
        file.write(" x, y, elevation \n\
-2.0, 2.0, 0.\n\
-1.0, 1.0, 0.\n\
0.0, 2.0 , 2.\n\
1.0, 1.0 , 2.\n\
 2.0,  1.0 ,3. \n\
 0.0,  0.0 , 0.\n\
 1.0,  0.0 , 1.\n\
 0.0,  -1.0, -1.\n\
 -0.2, -0.5, -0.7\n\
 -0.9, -1.5, -2.4\n\
 0.5,  -1.9, -1.4\n\
 3.0,  1.0 , 4.\n")
        file.close()
        geo = Geospatial_data(fileName)
        fileName_pts = tempfile.mktemp(".pts")
        geo.export_points_file(fileName_pts)
        f = fit_to_mesh(fileName_pts, vertices, triangles,
                                alpha=0.0, max_read_lines=2)
        answer = linear_function(vertices)
        #print "f\n",f
        #print "answer\n",answer
        assert num.allclose(f, answer)
        os.remove(fileName)
        os.remove(fileName_pts)

    def test_fit_to_mesh_pts_passing_mesh_in(self):
        a = [-1.0, 0.0]
        b = [3.0, 4.0]
        c = [4.0,1.0]
        d = [-3.0, 2.0] #3
        e = [-1.0,-2.0]
        f = [1.0, -2.0] #5

        vertices = [a, b, c, d,e,f]
        triangles = [[0,1,3], [1,0,2], [0,4,5], [0,5,2]] #abd bac aef afc


        fileName = tempfile.mktemp(".txt")
        file = open(fileName,"w")
        file.write(" x, y, elevation \n\
-2.0, 2.0, 0.\n\
-1.0, 1.0, 0.\n\
0.0, 2.0 , 2.\n\
1.0, 1.0 , 2.\n\
 2.0,  1.0 ,3. \n\
 0.0,  0.0 , 0.\n\
 1.0,  0.0 , 1.\n\
 0.0,  -1.0, -1.\n\
 -0.2, -0.5, -0.7\n\
 -0.9, -1.5, -2.4\n\
 0.5,  -1.9, -1.4\n\
 3.0,  1.0 , 4.\n")
        file.close()
        geo = Geospatial_data(fileName)
        fileName_pts = tempfile.mktemp(".pts")
        geo.export_points_file(fileName_pts)
        f = fit_to_mesh(fileName_pts, vertices, triangles,
                                alpha=0.0, max_read_lines=2)
        answer = linear_function(vertices)
        #print "f\n",f
        #print "answer\n",answer
        assert num.allclose(f, answer)
        os.remove(fileName)
        os.remove(fileName_pts)

    def test_fit_to_mesh(self):

        a = [-1.0, 0.0]
        b = [3.0, 4.0]
        c = [4.0,1.0]
        d = [-3.0, 2.0] #3
        e = [-1.0,-2.0]
        f = [1.0, -2.0] #5

        vertices = [a, b, c, d,e,f]
        triangles = [[0,1,3], [1,0,2], [0,4,5], [0,5,2]] #abd bac aef afc


        fileName = tempfile.mktemp(".ddd")
        file = open(fileName,"w")
        file.write(" x,y, elevation \n\
-2.0, 2.0, 0.\n\
-1.0, 1.0, 0.\n\
0.0, 2.0 , 2.\n\
1.0, 1.0 , 2.\n\
 2.0,  1.0 ,3. \n\
 0.0,  0.0 , 0.\n\
 1.0,  0.0 , 1.\n\
 0.0,  -1.0, -1.\n\
 -0.2, -0.5, -0.7\n\
 -0.9, -1.5, -2.4\n\
 0.5,  -1.9, -1.4\n\
 3.0,  1.0 , 4.\n")
        file.close()

        f = fit_to_mesh(fileName, vertices, triangles,
                                alpha=0.0, max_read_lines=2)
                        #use_cache=True, verbose=True)
        answer = linear_function(vertices)
        #print "f\n",f
        #print "answer\n",answer
        assert num.allclose(f, answer)

        os.remove(fileName)

    def test_fit_to_mesh_using_jacobi_precon(self):

        a = [-1.0, 0.0]
        b = [3.0, 4.0]
        c = [4.0,1.0]
        d = [-3.0, 2.0] #3
        e = [-1.0,-2.0]
        f = [1.0, -2.0] #5

        vertices = [a, b, c, d,e,f]
        triangles = [[0,1,3], [1,0,2], [0,4,5], [0,5,2]] #abd bac aef afc


        fileName = tempfile.mktemp(".ddd")
        file = open(fileName,"w")
        file.write(" x,y, elevation \n\
-2.0, 2.0, 0.\n\
-1.0, 1.0, 0.\n\
0.0, 2.0 , 2.\n\
1.0, 1.0 , 2.\n\
 2.0,  1.0 ,3. \n\
 0.0,  0.0 , 0.\n\
 1.0,  0.0 , 1.\n\
 0.0,  -1.0, -1.\n\
 -0.2, -0.5, -0.7\n\
 -0.9, -1.5, -2.4\n\
 0.5,  -1.9, -1.4\n\
 3.0,  1.0 , 4.\n")
        file.close()
        f = fit_to_mesh(fileName, vertices, triangles,
                                alpha=0.0, max_read_lines=2, use_c_cg=True, cg_precon='Jacobi')
                        #use_cache=True, verbose=True)
        answer = linear_function(vertices)
        #print "f\n",f
        #print "answer\n",answer

        assert num.allclose(f, answer)

        os.remove(fileName)

    def test_fit_to_mesh_using_jacobi_precon_no_c_cg(self):

        a = [-1.0, 0.0]
        b = [3.0, 4.0]
        c = [4.0,1.0]
        d = [-3.0, 2.0] #3
        e = [-1.0,-2.0]
        f = [1.0, -2.0] #5

        vertices = [a, b, c, d,e,f]
        triangles = [[0,1,3], [1,0,2], [0,4,5], [0,5,2]] #abd bac aef afc


        fileName = tempfile.mktemp(".ddd")
        file = open(fileName,"w")
        file.write(" x,y, elevation \n\
-2.0, 2.0, 0.\n\
-1.0, 1.0, 0.\n\
0.0, 2.0 , 2.\n\
1.0, 1.0 , 2.\n\
 2.0,  1.0 ,3. \n\
 0.0,  0.0 , 0.\n\
 1.0,  0.0 , 1.\n\
 0.0,  -1.0, -1.\n\
 -0.2, -0.5, -0.7\n\
 -0.9, -1.5, -2.4\n\
 0.5,  -1.9, -1.4\n\
 3.0,  1.0 , 4.\n")
        file.close()
        f = fit_to_mesh(fileName, vertices, triangles,
                                alpha=0.0, max_read_lines=2, use_c_cg=False, cg_precon='Jacobi')
                        #use_cache=True, verbose=True)
        answer = linear_function(vertices)
        #print "f\n",f
        #print "answer\n",answer

        assert num.allclose(f, answer)

        os.remove(fileName)

    def test_fit_to_mesh_2_atts(self):

        a = [-1.0, 0.0]
        b = [3.0, 4.0]
        c = [4.0,1.0]
        d = [-3.0, 2.0] #3
        e = [-1.0,-2.0]
        f = [1.0, -2.0] #5

        vertices = [a, b, c, d,e,f]
        triangles = [[0,1,3], [1,0,2], [0,4,5], [0,5,2]] #abd bac aef afc


        fileName = tempfile.mktemp(".ddd")
        file = open(fileName,"w")
        # the 2nd att name is wacky so it's the first key off a hash table
        file.write(" x,y, elevation, afriqction \n\
-2.0, 2.0, 0., 0.\n\
-1.0, 1.0, 0., 0.\n\
0.0, 2.0 , 2., 20.\n\
1.0, 1.0 , 2., 20.\n\
 2.0,  1.0 ,3., 30. \n\
 0.0,  0.0 , 0., 0.\n\
 1.0,  0.0 , 1., 10.\n\
 0.0,  -1.0, -1., -10.\n\
 -0.2, -0.5, -0.7, -7.\n\
 -0.9, -1.5, -2.4, -24. \n\
 0.5,  -1.9, -1.4, -14. \n\
 3.0,  1.0 , 4., 40. \n")
        file.close()

        f = fit_to_mesh(fileName, vertices, triangles,
                        alpha=0.0,
                        attribute_name='elevation', max_read_lines=2)
        answer = linear_function(vertices)
        #print "f\n",f
        #print "answer\n",answer
        assert num.allclose(f, answer)
        os.remove(fileName)





    def test_fit_and_interpolation(self):

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe, daf, dae
        triangles = [[1,0,2], [1,2,4], [4,2,5], [3,1,4]]

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
        interp = Fit(points, triangles,
                                alpha=0.0)

        answer = linear_function(points)

        f = interp.fit(data_points, z)

        #print "f",f
        #print "answer",answer
        assert num.allclose(f, answer)

    def test_fit_and_interpolation_using_jacobi_precon(self):

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe, daf, dae
        triangles = [[1,0,2], [1,2,4], [4,2,5], [3,1,4]]

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
        interp = Fit(points, triangles,
                                alpha=0.0,cg_precon='Jacobi')

        answer = linear_function(points)

        f = interp.fit(data_points, z)

        #print "f",f
        #print "answer",answer
        assert num.allclose(f, answer,atol=5)


    def test_smoothing_and_interpolation(self):

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe, daf, dae
        triangles = [[1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        #Get (too few!) datapoints
        data_points = [[ 0.66666667, 0.66666667],
                       [ 1.33333333, 1.33333333],
                       [ 2.66666667, 0.66666667],
                       [ 0.66666667, 2.66666667]]

        z = linear_function(data_points)
        answer = linear_function(points)

        #Make interpolator with too few data points and no smoothing
        interp = Fit(points, triangles, alpha=0.0)
        #Must raise an exception
        try:
            f = interp.fit(data_points,z)
        except TooFewPointsError:
            pass

        #Now try with smoothing parameter
        interp = Fit(points, triangles, alpha=1.0e-13)

        f = interp.fit(data_points,z)
        #f will be different from answer due to smoothing
        assert num.allclose(f, answer,atol=5)


    #Tests of smoothing matrix
    def test_smoothing_matrix_one_triangle(self):
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        points = [a, b, c]

        vertices = [ [1,0,2] ]   #bac

        interp = Fit(points, vertices)

        assert num.allclose(interp.get_D(), [[1, -0.5, -0.5],
                                             [-0.5, 0.5, 0],
                                             [-0.5, 0, 0.5]])

        #Define f(x,y) = x
        f = num.array([0,0,2]) #Value at global vertex 2

        #Check that int (df/dx)**2 + (df/dy)**2 dx dy =
        #           int 1 dx dy = area = 2
        assert num.dot(num.dot(f, interp.get_D()), f) == 2

        #Define f(x,y) = y
        f = num.array([0,2,0])  #Value at global vertex 1

        #Check that int (df/dx)**2 + (df/dy)**2 dx dy =
        #           int 1 dx dy = area = 2
        assert num.dot(num.dot(f, interp.get_D()), f) == 2

        #Define f(x,y) = x+y
        f = num.array([0,2,2])  #Values at global vertex 1 and 2

        #Check that int (df/dx)**2 + (df/dy)**2 dx dy =
        #           int 2 dx dy = 2*area = 4
        assert num.dot(num.dot(f, interp.get_D()), f) == 4


    def test_smoothing_matrix_more_triangles(self):
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe, daf, dae
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        interp = Fit(points, vertices)


        #assert allclose(interp.get_D(), [[1, -0.5, -0.5],
        #                           [-0.5, 0.5, 0],
        #                           [-0.5, 0, 0.5]])

        #Define f(x,y) = x
        f = num.array([0,0,2,0,2,4]) #f evaluated at points a-f

        #Check that int (df/dx)**2 + (df/dy)**2 dx dy =
        #           int 1 dx dy = total area = 8
        assert num.dot(num.dot(f, interp.get_D()), f) == 8

        #Define f(x,y) = y
        f = num.array([0,2,0,4,2,0]) #f evaluated at points a-f

        #Check that int (df/dx)**2 + (df/dy)**2 dx dy =
        #           int 1 dx dy = area = 8
        assert num.dot(num.dot(f, interp.get_D()), f) == 8

        #Define f(x,y) = x+y
        f = num.array([0,2,2,4,4,4])  #f evaluated at points a-f

        #Check that int (df/dx)**2 + (df/dy)**2 dx dy =
        #           int 2 dx dy = 2*area = 16
        assert num.dot(num.dot(f, interp.get_D()), f) == 16



    def test_fit_and_interpolation_with_different_origins(self):
        """Fit a surface to one set of points. Then interpolate that surface
        using another set of points.
        This test tests situtaion where points and mesh belong to a different
        coordinate system as defined by origin.
        """

        #Setup mesh used to represent fitted function
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe, daf, dae
        triangles = [[1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        #Datapoints to fit from
        data_points1 = [[ 0.66666667, 0.66666667],
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


        #First check that things are OK when using same origin
        mesh_origin = (56, 290000, 618000) #zone, easting, northing
        data_origin = (56, 290000, 618000) #zone, easting, northing


        #Fit surface to mesh
        interp = Fit(points, triangles,
                               alpha=0.0,
                               mesh_origin = mesh_origin)

        data_geo_spatial = Geospatial_data(data_points1,
                         geo_reference = Geo_reference(56, 290000, 618000))
        z = linear_function(data_points1) #Example z-values
        f = interp.fit(data_geo_spatial, z)   #Fitted values at vertices

        #Shift datapoints according to new origins
        for k in range(len(data_points1)):
            data_points1[k][0] += mesh_origin[1] - data_origin[1]
            data_points1[k][1] += mesh_origin[2] - data_origin[2]



        #Fit surface to mesh
        interp = Fit(points, triangles, alpha=0.0)
        #Fitted values at vertices (using same z as before)
        f1 = interp.fit(data_points1,z)

        assert num.allclose(f,f1), 'Fit should have been unaltered'


    def test_smooth_attributes_to_mesh_function(self):
        #Testing 2 attributes smoothed to the mesh


        a = [0.0, 0.0]
        b = [0.0, 5.0]
        c = [5.0, 0.0]
        points = [a, b, c]
        triangles = [ [1,0,2] ]   #bac

        d1 = [1.0, 1.0]
        d2 = [1.0, 3.0]
        d3 = [3.0, 1.0]
        z1 = [2, 4]
        z2 = [4, 8]
        z3 = [4, 8]
        data_coords = [d1, d2, d3]
        z = [z1, z2, z3]

        f = fit_to_mesh(data_coords, vertex_coordinates=points,
                        triangles=triangles,
                        point_attributes=z, alpha=0.0)
        answer = [[0, 0], [5., 10.], [5., 10.]]

        assert num.allclose(f, answer)

    def test_fit_to_mesh_w_georef(self):
        """Simple check that georef works at the fit_to_mesh level
        """

        from anuga.coordinate_transforms.geo_reference import Geo_reference

        #Mesh
        vertex_coordinates = [[0.76, 0.76],
                              [0.76, 5.76],
                              [5.76, 0.76]]
        triangles = [[0,2,1]]

        mesh_geo = Geo_reference(56,-0.76,-0.76)
        #print "mesh_geo.get_absolute(vertex_coordinates)", \
         #     mesh_geo.get_absolute(vertex_coordinates)

        #Data
        data_points = [[ 201.0, 401.0],
                       [ 201.0, 403.0],
                       [ 203.0, 401.0]]

        z = [2, 4, 4]

        data_geo = Geo_reference(56,-200,-400)

        #print "data_geo.get_absolute(data_points)", \
        #      data_geo.get_absolute(data_points)

        #Fit
        zz = fit_to_mesh(data_points, vertex_coordinates=vertex_coordinates,
                         triangles=triangles,
                         point_attributes=z,
                         data_origin = data_geo.get_origin(),
                         mesh_origin = mesh_geo.get_origin(),
                         alpha = 0)
        assert num.allclose( zz, [0,5,5] )


    def test_fit_to_mesh_file2domain(self):
        from anuga.load_mesh.loadASCII import import_mesh_file, \
             export_mesh_file
        import tempfile
        import os

        # create a .tsh file, no user outline
        mesh_dic = {}
        mesh_dic['vertices'] = [[0.0, 0.0],
                                          [0.0, 5.0],
                                          [5.0, 0.0]]
        mesh_dic['triangles'] =  [[0, 2, 1]]
        mesh_dic['segments'] = [[0, 1], [2, 0], [1, 2]]
        mesh_dic['triangle_tags'] = ['']
        mesh_dic['vertex_attributes'] = [[], [], []]
        mesh_dic['vertiex_attribute_titles'] = []
        mesh_dic['triangle_neighbors'] = [[-1, -1, -1]]
        mesh_dic['segment_tags'] = ['external',
                                                  'external',
                                                  'external']
        mesh_file = tempfile.mktemp(".tsh")
        export_mesh_file(mesh_file,mesh_dic)

        # create a points .csv file
        point_file = tempfile.mktemp(".csv")
        fd = open(point_file,'w')
        fd.write("x,y, elevation, stage \n\
        1.0, 1.0,2.,4 \n\
        1.0, 3.0,4,8 \n\
        3.0,1.0,4.,8 \n")
        fd.close()

        mesh_output_file = tempfile.mktemp(".tsh")
        fit_to_mesh_file(mesh_file,
                         point_file,
                         mesh_output_file,
                         alpha = 0.0)
        # load in the .tsh file we just wrote
        mesh_dic = import_mesh_file(mesh_output_file)
        #print "mesh_dic",mesh_dic
        ans =[[0.0, 0.0],
              [5.0, 10.0],
              [5.0,10.0]]
        assert num.allclose(mesh_dic['vertex_attributes'],ans)

        self.assertTrue(mesh_dic['vertex_attribute_titles']  ==
                        ['elevation','stage'],
                        'test_fit_to_mesh_file failed')
        domain = Domain(mesh_output_file, use_cache=True, verbose=False)

        answer = [0., 5., 5.]
        assert num.allclose(domain.quantities['elevation'].vertex_values,
                            answer)
        #clean up
        os.remove(mesh_file)
        os.remove(point_file)
        os.remove(mesh_output_file)

    def test_fit_to_mesh_file3(self):
        from anuga.load_mesh.loadASCII import import_mesh_file, \
             export_mesh_file
        import tempfile
        import os

        # create a .tsh file, no user outline
        mesh_dic = {}
        mesh_dic['vertices'] = [[0.76, 0.76],
                                          [0.76, 5.76],
                                          [5.76, 0.76]]
        mesh_dic['triangles'] =  [[0, 2, 1]]
        mesh_dic['segments'] = [[0, 1], [2, 0], [1, 2]]
        mesh_dic['triangle_tags'] = ['']
        mesh_dic['vertex_attributes'] = [[], [], []]
        mesh_dic['vertiex_attribute_titles'] = []
        mesh_dic['triangle_neighbors'] = [[-1, -1, -1]]
        mesh_dic['segment_tags'] = ['external',
                                                  'external',
                                                  'external']
        mesh_dic['geo_reference'] = Geo_reference(56,-0.76,-0.76)
        mesh_file = tempfile.mktemp(".tsh")
        export_mesh_file(mesh_file,mesh_dic)

        # create a points .csv file
        point_file = tempfile.mktemp(".csv")
        fd = open(point_file,'w')
        fd.write("x,y, elevation, stage \n\
        1.0, 1.0,2.,4 \n\
        1.0, 3.0,4,8 \n\
        3.0,1.0,4.,8 \n")
        fd.close()

        mesh_output_file = tempfile.mktemp(".tsh")
        fit_to_mesh_file(mesh_file,
                         point_file,
                         mesh_output_file,
                         alpha = 0.0)
        # load in the .tsh file we just wrote
        mesh_dic = import_mesh_file(mesh_output_file)
        #print "mesh_dic",mesh_dic
        ans =[[0.0, 0.0],
              [5.0, 10.0],
              [5.0,10.0]]
        assert num.allclose(mesh_dic['vertex_attributes'],ans)

        self.assertTrue(mesh_dic['vertex_attribute_titles']  ==
                        ['elevation','stage'],
                        'test_fit_to_mesh_file failed')

        #clean up
        os.remove(mesh_file)
        os.remove(point_file)
        os.remove(mesh_output_file)

    def test_fit_to_mesh_fileII(self):
        from anuga.load_mesh.loadASCII import import_mesh_file, \
             export_mesh_file
        import tempfile
        import os

        # create a .tsh file, no user outline
        mesh_dic = {}
        mesh_dic['vertices'] = [[0.0, 0.0],
                                [0.0, 5.0],
                                [5.0, 0.0]]
        mesh_dic['triangles'] =  [[0, 2, 1]]
        mesh_dic['segments'] = [[0, 1], [2, 0], [1, 2]]
        mesh_dic['triangle_tags'] = ['']
        mesh_dic['vertex_attributes'] = [[1,2], [1,2], [1,2]]
        mesh_dic['vertex_attribute_titles'] = ['density', 'temp']
        mesh_dic['triangle_neighbors'] = [[-1, -1, -1]]
        mesh_dic['segment_tags'] = ['external',
                                                  'external',
                                                  'external']
        mesh_file = tempfile.mktemp(".tsh")
        export_mesh_file(mesh_file,mesh_dic)

        # create a points .csv file
        point_file = tempfile.mktemp(".csv")
        fd = open(point_file,'w')
        fd.write("x,y,elevation, stage \n\
        1.0, 1.0,2.,4 \n\
        1.0, 3.0,4,8 \n\
        3.0,1.0,4.,8 \n")
        fd.close()

        mesh_output_file = "new_triangle.tsh"
        fit_to_mesh_file(mesh_file,
                         point_file,
                         mesh_output_file,
                         alpha = 0.0)
        # load in the .tsh file we just wrote
        mesh_dic = import_mesh_file(mesh_output_file)

        assert num.allclose(mesh_dic['vertex_attributes'],
                            [[1.0, 2.0,0.0, 0.0],
                             [1.0, 2.0,5.0, 10.0],
                             [1.0, 2.0,5.0,10.0]])

        self.assertTrue(mesh_dic['vertex_attribute_titles']  ==
                        ['density', 'temp','elevation','stage'],
                        'test_fit_to_mesh_file failed')

        #clean up
        os.remove(mesh_file)
        os.remove(mesh_output_file)
        os.remove(point_file)

    def test_fit_to_mesh_file_errors(self):
        from anuga.load_mesh.loadASCII import import_mesh_file, export_mesh_file
        import tempfile
        import os

        # create a .tsh file, no user outline
        mesh_dic = {}
        mesh_dic['vertices'] = [[0.0, 0.0],[0.0, 5.0],[5.0, 0.0]]
        mesh_dic['triangles'] =  [[0, 2, 1]]
        mesh_dic['segments'] = [[0, 1], [2, 0], [1, 2]]
        mesh_dic['triangle_tags'] = ['']
        mesh_dic['vertex_attributes'] = [[1,2], [1,2], [1,2]]
        mesh_dic['vertex_attribute_titles'] = ['density', 'temp']
        mesh_dic['triangle_neighbors'] = [[-1, -1, -1]]
        mesh_dic['segment_tags'] = ['external', 'external','external']
        mesh_file = tempfile.mktemp(".tsh")
        export_mesh_file(mesh_file,mesh_dic)

        # create a bad points .csv file
        point_file = tempfile.mktemp(".csv")
        fd = open(point_file,'w')
        fd.write("x,y,elevation stage \n\
        1.0, 1.0,2.,4 \n\
        1.0, 3.0,4,8 \n\
        3.0,1.0,4.,8 \n")
        fd.close()

        mesh_output_file = "new_triangle.tsh"
        try:
            fit_to_mesh_file(mesh_file, point_file,
                             mesh_output_file, display_errors = False)
        except SyntaxError:
            pass
        else:
            #self.assertTrue(0 ==1,  'Bad file did not raise error!')
            raise Exception('Bad file did not raise error!')

        #clean up
        os.remove(mesh_file)
        os.remove(point_file)

    def test_fit_to_mesh_file_errorsII(self):
        from anuga.load_mesh.loadASCII import import_mesh_file, export_mesh_file
        import tempfile
        import os

        # create a .tsh file, no user outline
        mesh_file = tempfile.mktemp(".tsh")
        fd = open(mesh_file,'w')
        fd.write("unit testing a bad .tsh file \n")
        fd.close()

        # create a points .csv file
        point_file = tempfile.mktemp(".csv")
        fd = open(point_file,'w')
        fd.write("x,y,elevation, stage \n\
        1.0, 1.0,2.,4 \n\
        1.0, 3.0,4,8 \n\
        3.0,1.0,4.,8 \n")
        fd.close()

        mesh_output_file = "new_triangle.tsh"
        try:
            fit_to_mesh_file(mesh_file, point_file,
                             mesh_output_file, display_errors = False)
        except IOError:
            pass
        else:
            raise Exception('Bad file did not raise error!')

        #clean up
        try:
            os.remove(mesh_file)
            os.remove(point_file)
        except PermissionError:
            pass

    def test_fit_to_mesh_file_errorsIII(self):
        from anuga.load_mesh.loadASCII import import_mesh_file, export_mesh_file
        import tempfile
        import os

        # create a .tsh file, no user outline
        mesh_dic = {}
        mesh_dic['vertices'] = [[0.0, 0.0],[0.0, 5.0],[5.0, 0.0]]
        mesh_dic['triangles'] =  [[0, 2, 1]]
        mesh_dic['segments'] = [[0, 1], [2, 0], [1, 2]]
        mesh_dic['triangle_tags'] = ['']
        mesh_dic['vertex_attributes'] = [[1,2], [1,2], [1,2]]
        mesh_dic['vertex_attribute_titles'] = ['density', 'temp']
        mesh_dic['triangle_neighbors'] = [[-1, -1, -1]]
        mesh_dic['segment_tags'] = ['external', 'external','external']
        mesh_file = tempfile.mktemp(".tsh")
        export_mesh_file(mesh_file,mesh_dic)


        # create a points .csv file
        point_file = tempfile.mktemp(".csv")
        fd = open(point_file,'w')
        fd.write("x,y,elevation, stage \n\
        1.0, 1.0,2.,4 \n\
        1.0, 3.0,4,8 \n\
        3.0,1.0,4.,8 \n")
        fd.close()

        #This a deliberately illegal filename to invoke the error.
        mesh_output_file = r".../\z\z:ya.tsh"

        try:
            fit_to_mesh_file(mesh_file, point_file,
                             mesh_output_file, display_errors = False)
        except IOError:
            pass
        else:
            raise Exception('Bad file did not raise error!')

        #clean up
        os.remove(mesh_file)
        os.remove(point_file)


#-------------------------------------------------------------
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_Fit,'test')
    runner = unittest.TextTestRunner() #verbosity=1)
    runner.run(suite)

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


        fd, fileName = tempfile.mkstemp(".ddd")
        os.close(fd)
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
        fd, txt_file = tempfile.mkstemp(".txt")
        os.close(fd)
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


        fd, fileName = tempfile.mkstemp(".txt")
        os.close(fd)
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
        fd, fileName_pts = tempfile.mkstemp(".pts")
        os.close(fd)
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


        fd, fileName = tempfile.mkstemp(".txt")
        os.close(fd)
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
        fd, fileName_pts = tempfile.mkstemp(".pts")
        os.close(fd)
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


        fd, fileName = tempfile.mkstemp(".ddd")
        os.close(fd)
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


        fd, fileName = tempfile.mkstemp(".ddd")
        os.close(fd)
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


        fd, fileName = tempfile.mkstemp(".ddd")
        os.close(fd)
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


        fd, fileName = tempfile.mkstemp(".ddd")
        os.close(fd)
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
        fd, mesh_file = tempfile.mkstemp(".tsh")
        os.close(fd)
        export_mesh_file(mesh_file,mesh_dic)

        # create a points .csv file
        fd, point_file = tempfile.mkstemp(".csv")
        os.close(fd)
        fd = open(point_file,'w')
        fd.write("x,y, elevation, stage \n\
        1.0, 1.0,2.,4 \n\
        1.0, 3.0,4,8 \n\
        3.0,1.0,4.,8 \n")
        fd.close()

        fd, mesh_output_file = tempfile.mkstemp(".tsh")
        os.close(fd)
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
        fd, mesh_file = tempfile.mkstemp(".tsh")
        os.close(fd)
        export_mesh_file(mesh_file,mesh_dic)

        # create a points .csv file
        fd, point_file = tempfile.mkstemp(".csv")
        os.close(fd)
        fd = open(point_file,'w')
        fd.write("x,y, elevation, stage \n\
        1.0, 1.0,2.,4 \n\
        1.0, 3.0,4,8 \n\
        3.0,1.0,4.,8 \n")
        fd.close()

        fd, mesh_output_file = tempfile.mkstemp(".tsh")
        os.close(fd)
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
        fd, mesh_file = tempfile.mkstemp(".tsh")
        os.close(fd)
        export_mesh_file(mesh_file,mesh_dic)

        # create a points .csv file
        fd, point_file = tempfile.mkstemp(".csv")
        os.close(fd)
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
        fd, mesh_file = tempfile.mkstemp(".tsh")
        os.close(fd)
        export_mesh_file(mesh_file,mesh_dic)

        # create a bad points .csv file
        fd, point_file = tempfile.mkstemp(".csv")
        os.close(fd)
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

        # create a .tsh file, no user outline
        fd, mesh_file = tempfile.mkstemp(".tsh")
        os.close(fd)
        fd = open(mesh_file,'w')
        fd.write("unit testing a bad .tsh file \n")
        fd.close()

        # create a points .csv file
        fd, point_file = tempfile.mkstemp(".csv")
        os.close(fd)
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
        except OSError:
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
        fd, mesh_file = tempfile.mkstemp(".tsh")
        os.close(fd)
        export_mesh_file(mesh_file,mesh_dic)


        # create a points .csv file
        fd, point_file = tempfile.mkstemp(".csv")
        os.close(fd)
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
        except OSError:
            pass
        else:
            raise Exception('Bad file did not raise error!')

        #clean up
        os.remove(mesh_file)
        os.remove(point_file)


#-------------------------------------------------------------

class Test_Fit_coverage(unittest.TestCase):
    """Targeted tests for previously-uncovered branches in fit.py and
    general_fit_interpolate.py."""

    # ------------------------------------------------------------------
    # Shared mesh: 3 vertices, 1 triangle
    # ------------------------------------------------------------------
    _pts  = num.array([[0.,0.],[0.,5.],[5.,0.]], float)
    _tris = num.array([[1,0,2]], int)
    _data = num.array([[1.,1.],[1.,3.],[3.,1.]], float)
    _z    = num.array([2., 4., 4.], float)

    def _fit(self, **kw):
        f = Fit(self._pts, self._tris, **kw)
        f._build_matrix_AtA_Atz(self._data, self._z)
        return f

    # ------------------------------------------------------------------
    # Fit.__init__ verbose path
    # ------------------------------------------------------------------
    def test_fit_init_verbose(self):
        """Fit.__init__ with verbose=True covers the log.info smoothing path."""
        f = Fit(self._pts, self._tris, verbose=True)
        assert f.D is not None

    # ------------------------------------------------------------------
    # Fit.fit verbose path and no-data warning
    # ------------------------------------------------------------------
    def test_fit_fit_verbose(self):
        """Fit.fit with verbose=True covers the print('Initializing') branch."""
        f = self._fit()
        result = f.fit(self._data, self._z, verbose=True)
        assert result is not None

    def test_fit_fit_no_data_verbose(self):
        """fit() with no data and verbose=True logs the no-data warning."""
        f = self._fit()
        # Force a second call with None so the no-data branch fires
        result = f.fit(None, verbose=True)
        assert result is not None

    # ------------------------------------------------------------------
    # TooFewPointsError when alpha=0 and not enough points
    # ------------------------------------------------------------------
    def test_too_few_points_raises(self):
        """Fewer data points than mesh vertices with alpha=0 raises TooFewPointsError."""
        f = Fit(self._pts, self._tris, alpha=0.0)
        # 1 data point < 3 vertices
        f._build_matrix_AtA_Atz(num.array([[1., 1.]]), num.array([2.]))
        with self.assertRaises(TooFewPointsError):
            f.fit()

    # ------------------------------------------------------------------
    # fit_to_mesh use_cache path
    # ------------------------------------------------------------------
    def test_fit_to_mesh_use_cache(self):
        """fit_to_mesh with use_cache=True exercises the cache() branch."""
        result = fit_to_mesh(self._data, self._pts, self._tris,
                             point_attributes=self._z,
                             alpha=0.1, use_cache=True, verbose=False)
        assert result is not None
        assert len(result) == len(self._pts)

    # ------------------------------------------------------------------
    # _build_matrix_AtA_Atz verbose dot path
    # ------------------------------------------------------------------
    def test_build_matrix_verbose_dot(self):
        """_build_matrix_AtA_Atz with verbose=True, output='dot' hits print path."""
        f = Fit(self._pts, self._tris)
        f._build_matrix_AtA_Atz(self._data, self._z, verbose=True, output='dot')
        assert f.AtA is not None

    # ------------------------------------------------------------------
    # general_fit_interpolate verbose and build_quad_tree
    # ------------------------------------------------------------------
    def test_fitinterpolate_verbose(self):
        """FitInterpolate.__init__ with verbose=True covers the log.info paths."""
        from anuga.fit_interpolate.general_fit_interpolate import FitInterpolate
        fi = FitInterpolate(vertex_coordinates=self._pts,
                            triangles=self._tris, verbose=True)
        assert fi.mesh is not None

    def test_fitinterpolate_build_quad_tree(self):
        """build_quad_tree() rebuilds self.root."""
        from anuga.fit_interpolate.general_fit_interpolate import FitInterpolate
        fi = FitInterpolate(vertex_coordinates=self._pts, triangles=self._tris)
        old_root = fi.root
        fi.build_quad_tree()
        assert fi.root is not None

    def test_fitinterpolate_mesh_none(self):
        """FitInterpolate with no vertices/triangles leaves mesh as None."""
        from anuga.fit_interpolate.general_fit_interpolate import FitInterpolate
        fi = FitInterpolate()
        assert fi.mesh is None

    def test_get_build_quadtree_time(self):
        """get_build_quadtree_time() returns a non-negative float."""
        from anuga.fit_interpolate.general_fit_interpolate import get_build_quadtree_time, FitInterpolate
        FitInterpolate(vertex_coordinates=self._pts, triangles=self._tris)
        assert get_build_quadtree_time() >= 0.0

    def test_fitinterpolate_repr(self):
        """__repr__ returns a non-empty string."""
        from anuga.fit_interpolate.general_fit_interpolate import FitInterpolate
        fi = FitInterpolate(vertex_coordinates=self._pts, triangles=self._tris)
        assert 'Interpolation' in repr(fi)

    # ------------------------------------------------------------------
    # interpolate.py — output_centroids=True path
    # ------------------------------------------------------------------
    def test_interpolate_output_centroids(self):
        """interpolate_block with output_centroids=True returns centroid values."""
        from anuga.fit_interpolate.interpolate import Interpolate
        vertices = num.array([[-1.,0.],[3.,4.],[4.,1.],
                               [-3.,2.],[-1.,-2.],[1.,-2.]], float)
        triangles = num.array([[0,1,3],[1,0,2],[0,4,5],[0,5,2]])
        interp = Interpolate(vertices, triangles)
        f = linear_function(vertices)
        pts = num.array([[0.,0.],[1.,1.]])
        z = interp.interpolate(f, pts, output_centroids=True)
        assert z is not None
        assert len(z) == len(pts)
        # centroids list is populated
        assert len(interp.centroids) > 0

    # ------------------------------------------------------------------
    # interpolate.py — Interpolate verbose path in _build_interpolation_matrix_A
    # ------------------------------------------------------------------
    def test_interpolate_block_verbose(self):
        """interpolate_block with verbose=True runs log.info paths."""
        from anuga.fit_interpolate.interpolate import Interpolate
        vertices = num.array([[-1.,0.],[3.,4.],[4.,1.],
                               [-3.,2.],[-1.,-2.],[1.,-2.]], float)
        triangles = num.array([[0,1,3],[1,0,2],[0,4,5],[0,5,2]])
        interp = Interpolate(vertices, triangles)
        f = linear_function(vertices)
        pts = num.array([[0.,0.],[1.,1.]])
        z = interp.interpolate_block(f, pts, verbose=True)
        assert z is not None

    # ------------------------------------------------------------------
    # fit_to_mesh_file — verbose paths and display_errors=True error paths
    # ------------------------------------------------------------------
    def _make_mesh_file(self):
        """Helper: write a minimal .tsh mesh file, return path."""
        from anuga.load_mesh.loadASCII import export_mesh_file
        mesh_dic = {
            'vertices': [[0.0, 0.0], [0.0, 5.0], [5.0, 0.0]],
            'triangles': [[0, 2, 1]],
            'segments': [[0, 1], [2, 0], [1, 2]],
            'triangle_tags': [''],
            'vertex_attributes': [[], [], []],
            'vertex_attribute_titles': [],
            'triangle_neighbors': [[-1, -1, -1]],
            'segment_tags': ['external', 'external', 'external'],
        }
        import tempfile
        import os
        fd, path = tempfile.mkstemp('.tsh')
        os.close(fd)
        export_mesh_file(path, mesh_dic)
        return path

    def _make_points_file(self):
        """Helper: write a minimal points .csv file, return path."""
        import tempfile
        import os
        fd, path = tempfile.mkstemp('.csv')
        os.close(fd)
        with open(path, 'w') as f:
            f.write('x,y,elevation\n1.0,1.0,2.0\n1.0,3.0,4.0\n3.0,1.0,4.0\n')
        return path

    def test_fit_to_mesh_file_verbose(self):
        """fit_to_mesh_file with verbose=True covers all log.info branches."""
        import tempfile
        import os
        mesh_file = self._make_mesh_file()
        point_file = self._make_points_file()
        fd, out_file = tempfile.mkstemp('.tsh')
        os.close(fd)
        try:
            fit_to_mesh_file(mesh_file, point_file, out_file,
                             alpha=0.0, verbose=True)
            from anuga.load_mesh.loadASCII import import_mesh_file
            result = import_mesh_file(out_file)
            assert result['vertex_attribute_titles'] == ['elevation']
        finally:
            for p in (mesh_file, point_file, out_file):
                try:
                    os.remove(p)
                except OSError:
                    pass

    def test_fit_to_mesh_file_bad_mesh_display_errors_true(self):
        """fit_to_mesh_file logs warning when mesh file is bad and display_errors=True."""
        import tempfile
        import os
        fd, bad_mesh = tempfile.mkstemp('.tsh')
        os.close(fd)
        with open(bad_mesh, 'w') as f:
            f.write('not a valid tsh file\n')
        fd, point_file = tempfile.mkstemp('.csv')
        os.close(fd)
        with open(point_file, 'w') as f:
            f.write('x,y,elevation\n1.0,1.0,2.0\n')
        try:
            with self.assertRaises(OSError):
                fit_to_mesh_file(bad_mesh, point_file, '/tmp/out.tsh',
                                 display_errors=True)
        finally:
            for p in (bad_mesh, point_file):
                try:
                    os.remove(p)
                except OSError:
                    pass

    def test_fit_to_mesh_file_bad_points_display_errors_true(self):
        """fit_to_mesh_file logs warning when points file is bad and display_errors=True."""
        import tempfile
        import os
        mesh_file = self._make_mesh_file()
        # Bad points file: header missing comma between last two columns → SyntaxError
        fd, bad_pts = tempfile.mkstemp('.csv')
        os.close(fd)
        with open(bad_pts, 'w') as f:
            f.write('x,y,elevation stage\n1.0,1.0,2.0,4.0\n')
        fd, out_file = tempfile.mkstemp('.tsh')
        os.close(fd)
        try:
            with self.assertRaises(Exception):
                fit_to_mesh_file(mesh_file, bad_pts, out_file,
                                 display_errors=True)
        finally:
            for p in (mesh_file, bad_pts, out_file):
                try:
                    os.remove(p)
                except OSError:
                    pass

    def test_fit_to_mesh_file_bad_output_display_errors_true(self):
        """fit_to_mesh_file logs warning when output file path is bad and display_errors=True."""
        import os
        mesh_file = self._make_mesh_file()
        point_file = self._make_points_file()
        bad_out = r'.../\z\z:ya.tsh'
        try:
            with self.assertRaises(OSError):
                fit_to_mesh_file(mesh_file, point_file, bad_out,
                                 display_errors=True)
        finally:
            for p in (mesh_file, point_file):
                try:
                    os.remove(p)
                except OSError:
                    pass


#-------------------------------------------------------------

class Test_Alpha_Selection(unittest.TestCase):
    """Tests for L-curve automatic alpha selection."""

    # Medium mesh: 6 vertices, 4 triangles — enough triangles to see smoothing effect
    _pts = num.array([[0., 0.], [0., 2.], [2., 0.],
                      [0., 4.], [2., 2.], [4., 0.]], float)
    _tris = num.array([[1, 0, 2], [1, 2, 4], [4, 2, 5], [3, 1, 4]], int)

    # Dense data that exactly matches a linear function (z = x + y)
    _data = num.array([[0.5, 0.5], [0.5, 1.5], [1.5, 0.5],
                       [0.5, 2.5], [1.5, 1.5], [2.5, 0.5],
                       [0.5, 3.5], [1.5, 2.5], [2.5, 1.5], [3.5, 0.5]], float)
    _z = _data[:, 0] + _data[:, 1]  # exact linear: z = x + y

    def _make_fit(self):
        f = Fit(self._pts, self._tris, alpha=DEFAULT_ALPHA)
        f._build_matrix_AtA_Atz(self._data, self._z)
        return f

    # ------------------------------------------------------------------
    # select_alpha returns a float in the candidate set
    # ------------------------------------------------------------------
    def test_select_alpha_returns_valid(self):
        """select_alpha returns a positive float."""
        f = self._make_fit()
        alpha = f.select_alpha()
        assert isinstance(alpha, float)
        assert alpha > 0, 'alpha must be positive, got %g' % alpha

    # ------------------------------------------------------------------
    # return_curve gives three arrays of the right length
    # ------------------------------------------------------------------
    def test_select_alpha_return_curve(self):
        """select_alpha(return_curve=True) returns (alpha, (log_rss, log_s, kappa))."""
        f = self._make_fit()
        alpha, (log_rss, log_s, kappa) = f.select_alpha(return_curve=True)
        assert len(log_rss) == 20
        assert len(log_s) == 20
        assert len(kappa) == 20
        # log_rss should decrease as alpha increases (more smoothing → worse fit)
        assert log_rss[-1] > log_rss[0], 'rss should grow with alpha'
        # log_s should decrease as alpha increases (more smoothing → smaller f^T D f)
        assert log_s[0] > log_s[-1], 'smoothing norm should shrink with alpha'

    # ------------------------------------------------------------------
    # z_sq is accumulated correctly across blocking calls
    # ------------------------------------------------------------------
    def test_z_sq_accumulated(self):
        """z_sq equals the total sum of z² over all data."""
        f = self._make_fit()
        expected = float(self._z @ self._z)
        assert num.isclose(f.z_sq, expected), \
            'z_sq=%g expected %g' % (f.z_sq, expected)

    def test_z_sq_accumulated_blocking(self):
        """z_sq accumulates correctly when _build_matrix_AtA_Atz is called in blocks."""
        f = Fit(self._pts, self._tris)
        half = len(self._data) // 2
        f._build_matrix_AtA_Atz(self._data[:half], self._z[:half])
        f._build_matrix_AtA_Atz(self._data[half:], self._z[half:])
        expected = float(self._z @ self._z)
        assert num.isclose(f.z_sq, expected)

    # ------------------------------------------------------------------
    # alpha='auto' via fit_to_mesh produces finite results
    # ------------------------------------------------------------------
    def test_fit_to_mesh_alpha_auto(self):
        """fit_to_mesh with alpha='auto' runs end-to-end and produces finite results."""
        result = fit_to_mesh(self._data, self._pts, self._tris,
                             point_attributes=self._z,
                             alpha='auto', verbose=False)
        assert result is not None
        assert result.shape == (len(self._pts),)
        assert num.all(num.isfinite(result)), 'auto-alpha result contains non-finite values'

    # ------------------------------------------------------------------
    # L-curve finds a non-trivial corner when data is noisy
    # ------------------------------------------------------------------
    def test_select_alpha_noisy_data_picks_interior(self):
        """With noisy data the L-curve selects an interior alpha with positive kappa.

        Uses a 6-node mesh with 60 noisy points; with the 20-candidate log-scale
        (1e-6 … 100) the L-corner falls at alpha ~ 1–10, giving positive kappa
        and hitting the opt_alpha = float(alphas[best]) branch.
        """
        rng = num.random.default_rng(7)
        n_pts = 60
        xy = rng.uniform(0.5, 3.5, (n_pts, 2))
        z_noisy = xy[:, 0] + xy[:, 1] + rng.standard_normal(n_pts) * 1.5

        pts  = num.array([[0.,0.],[0.,2.],[2.,0.],[0.,4.],[2.,2.],[4.,0.]], float)
        tris = num.array([[1,0,2],[1,2,4],[4,2,5],[3,1,4]], int)

        f = Fit(pts, tris, alpha='auto')
        f._build_matrix_AtA_Atz(xy, z_noisy)
        alpha_opt, (_, _, kappa) = f.select_alpha(return_curve=True)

        best = int(num.argmax(kappa))
        assert kappa[best] > 0, \
            'Expected L-curve corner (positive kappa) for noisy data'
        assert best < 19, \
            'Expected interior corner, not last candidate'
        assert alpha_opt > 0, \
            'alpha_opt should be positive, got %g' % alpha_opt

    def test_select_alpha_multi_attribute(self):
        """select_alpha handles 2D Atz (multi-attribute) by using only first column."""
        pts  = num.array([[0.,0.],[0.,2.],[2.,0.],[0.,4.],[2.,2.],[4.,0.]], float)
        tris = num.array([[1,0,2],[1,2,4],[4,2,5],[3,1,4]], int)
        rng  = num.random.default_rng(3)
        n_pts = 30
        xy   = rng.uniform(0.5, 3.5, (n_pts, 2))
        # Two attributes: z1 = x+y, z2 = x-y
        z_2d = num.column_stack([xy[:,0]+xy[:,1], xy[:,0]-xy[:,1]])

        f = Fit(pts, tris, alpha='auto')
        f._build_matrix_AtA_Atz(xy, z_2d)
        assert f.Atz.ndim == 2, 'Expected 2D Atz for multi-attribute input'
        alpha = f.select_alpha()
        assert isinstance(alpha, float) and alpha > 0

    # ------------------------------------------------------------------
    # alpha='auto' verbose path
    # ------------------------------------------------------------------
    def test_fit_alpha_auto_verbose(self):
        """Fit.fit with alpha='auto' and verbose=True logs the selected alpha."""
        f = Fit(self._pts, self._tris, alpha='auto')
        result = f.fit(self._data, self._z, verbose=True)
        assert result is not None
        # alpha should now be a concrete float (resolved during fit)
        assert isinstance(f.alpha, float)

    def test_select_alpha_sparse_data_row_ptr_extended(self):
        """select_alpha works when data doesn't cover the highest-indexed nodes.

        When data only lands in triangles that don't reference node 5 (index 5
        out of 6), sparse_dok.num_rows < m-1, so row_ptr is shorter than m+1.
        The _to_scipy_csr helper must extend it before constructing the CSR matrix.
        """
        pts  = num.array([[0.,0.],[0.,2.],[2.,0.],[0.,4.],[2.,2.],[4.,0.]], float)
        tris = num.array([[1,0,2],[1,2,4],[4,2,5],[3,1,4]], int)
        # All points are away from triangle [4,2,5], so node 5 is never addressed
        xy = num.array([[0.5, 0.5], [0.5, 1.0], [0.8, 0.8], [0.3, 0.3], [0.7, 1.2]], float)
        z  = xy[:, 0] + xy[:, 1]

        f = Fit(pts, tris, alpha='auto')
        f._build_matrix_AtA_Atz(xy, z)
        alpha = f.select_alpha()
        assert isinstance(alpha, float) and alpha > 0

    def test_select_alpha_degenerate_falls_back_to_default(self):
        """Degenerate fallback: when kappa has no valid interior corner, DEFAULT_ALPHA is returned."""
        # 3 data points for 6 nodes — underdetermined.  On some numpy versions the
        # L-curve computation places the best-kappa index at the last candidate
        # (triggering the fallback); on others it finds an apparent interior corner.
        # We use return_curve=True to see what kappa says and then verify the
        # correct outcome for *that* platform, exercising the fallback logic itself.
        pts  = num.array([[0.,0.],[0.,2.],[2.,0.],[0.,4.],[2.,2.],[4.,0.]], float)
        tris = num.array([[1,0,2],[1,2,4],[4,2,5],[3,1,4]], int)
        xy   = num.array([[0.5, 0.5], [1.5, 1.5], [2.5, 0.5]], float)
        z    = num.array([1.0, 3.0, 3.0], float)

        f = Fit(pts, tris, alpha='auto')
        f._build_matrix_AtA_Atz(xy, z)
        alpha, (_, _, kappa) = f.select_alpha(return_curve=True)

        best = int(num.argmax(kappa))
        fallback_expected = (kappa[best] <= 0 or best == len(kappa) - 1)

        if fallback_expected:
            assert alpha == DEFAULT_ALPHA, \
                'Degenerate fallback should return DEFAULT_ALPHA, got %g' % alpha
        else:
            # L-curve found a genuine interior corner on this platform;
            # verify the returned alpha is within the candidate range.
            assert 1e-6 <= alpha <= 100, \
                'select_alpha returned out-of-range alpha %g' % alpha

    # ------------------------------------------------------------------
    # select_alpha fails gracefully before data is loaded
    # ------------------------------------------------------------------
    def test_select_alpha_requires_data(self):
        """select_alpha raises AssertionError if called before _build_matrix_AtA_Atz."""
        f = Fit(self._pts, self._tris)
        with self.assertRaises(AssertionError):
            f.select_alpha()


#-------------------------------------------------------------
if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_Fit)
    runner = unittest.TextTestRunner() #verbosity=1)
    runner.run(suite)

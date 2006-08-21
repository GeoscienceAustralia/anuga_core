#!/usr/bin/env python

#TEST
import sys
import unittest
from math import sqrt


from anuga.pyvolution.least_squares import *
from anuga.pyvolution.neighbour_mesh import Mesh

from Numeric import allclose, array, transpose

from anuga.coordinate_transforms.geo_reference import Geo_reference

def distance(x, y):
    return sqrt( sum( (array(x)-array(y))**2 ))

def linear_function(point):
    point = array(point)
    return point[:,0]+point[:,1]


class Test_Least_Squares(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_datapoint_at_centroid(self):
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        points = [a, b, c]
        vertices = [ [1,0,2] ]   #bac

        data = [ [2.0/3, 2.0/3] ] #Use centroid as one data point

        
        interp = Interpolation(points, vertices, data)
        assert allclose(interp.get_A(), [[1./3, 1./3, 1./3]])


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
        interp = Interpolation(points, triangles, data, alpha = 0.0,
                               max_points_per_cell = 4)
        #print "PDSG - interp.get_A()", interp.get_A()
        answer =  [ [ 0.06666667,  0.46666667,  0.46666667,  0.,
                      0., 0. , 0., 0., 0., 0.]]
        assert allclose(interp.get_A(), answer)
        interp.set_point_coordinates([[-30, -30]]) #point outside of mesh
        #print "PDSG - interp.get_A()", interp.get_A()
        answer =  [ [ 0.0,  0.0,  0.0,  0.,
                      0., 0. , 0., 0., 0., 0.]]
        assert allclose(interp.get_A(), answer)


        #point outside of quad tree root cell
        interp.set_point_coordinates([[-70, -70]])
        #print "PDSG - interp.get_A()", interp.get_A()
        answer =  [ [ 0.0,  0.0,  0.0,  0.,
                      0., 0. , 0., 0., 0., 0.]]
        assert allclose(interp.get_A(), answer)

    def test_expand_search(self):
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

        data = [ [4,4],
                 [-30,10],
                 [-20,0],
                 [-20,10],
                 [0,30],
                 [10,-40],
                 [10,-30],
                 [10,-20],
                 [10,10],
                 [10,20],
                 [10,30],
                 [10,40],
                 [20,10],
                 [25,45],
                 [30,0],
                 [30,10],
                 [30,30],
                 [30,40],
                 [30,50],
                 [40,10],
                 [40,30],
                 [40,40],
                 [40,50],
                 [50,20],
                 [50,30],
                 [50,40],
                 [50,50],
                 [30,0],
                 [-20,-20]]
        point_attributes = [ -400000,
                     10,
                     10,
                     10,
                    10,
                     10,
                     10,
                     10,
                     10,
                     10,
                     10,
                     10,
                     10,
                     10,
                    10,
                     10,
                     10,
                     10,
                     10,
                     10,
                     10,
                     10,
                     10,
                     10,
                     10,
                     10,
                     10,
                   10,
                   99]

        interp = Interpolation(points, triangles, data,
                               alpha=0.0, expand_search=False, #verbose = True, #False,
                               max_points_per_cell = 4)
        calc = interp.fit_points(point_attributes, )
        #print "calc",calc

        # the point at 4,4 is ignored.  An expanded search has to be done
        # to fine which triangel it's in.
        # An expanded search isn't done to find that the last point
        # isn't in the mesh.  But this isn't tested.
        answer= [ 10,
                     10,
                     10,
                     10,
                    10,
                     10,
                     10,
                     10,
                     10,
                     10]
        assert allclose(calc, answer)

    def test_quad_treeII(self):
        p0 = [-66.0, 14.0]
        p1 = [14.0, -66.0]
        p2 = [14.0, 14.0]
        p3 = [60.0, 20.0]
        p4 = [10.0, 60.0]
        p5 = [60.0, 60.0]

        points = [p0, p1, p2, p3, p4, p5]
        triangles = [ [0, 1, 2],
                      [3, 2, 1],
                      [0, 2, 4],
                      [4, 2, 5],
                      [5, 2, 3]]

        data = [ [-26.0,-26.0] ]
        interp = Interpolation(points, triangles, data, alpha = 0.0,
                 max_points_per_cell = 4)
        #print "PDSG - interp.get_A()", interp.get_A()
        answer =  [ [ 0.5,  0.5,  0.0,  0.,
                      0., 0.]]
        assert allclose(interp.get_A(), answer)
        interp.set_point_coordinates([[-30, -30]]) #point outside of mesh
        #print "PDSG -30,-30 - interp.get_A()", interp.get_A()
        answer =  [ [ 0.0,  0.0,  0.0,  0.,
                      0., 0.]]
        assert allclose(interp.get_A(), answer)


        #point outside of quad tree root cell
        interp.set_point_coordinates([[-70, -70]])
        #print "PDSG -70,-70 interp.get_A()", interp.get_A()
        answer =  [ [ 0.0,  0.0,  0.0,  0.,
                      0., 0. ]]
        assert allclose(interp.get_A(), answer)


    def test_datapoints_at_vertices(self):
        """Test that data points coinciding with vertices yield a diagonal matrix
        """

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        points = [a, b, c]
        vertices = [ [1,0,2] ]   #bac

        data = points #Use data at vertices

        interp = Interpolation(points, vertices, data)
        assert allclose(interp.get_A(), [[1., 0., 0.],
                                   [0., 1., 0.],
                                   [0., 0., 1.]])



    def test_datapoints_on_edge_midpoints(self):
        """Try datapoints midway on edges -
        each point should affect two matrix entries equally
        """

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        points = [a, b, c]
        vertices = [ [1,0,2] ]   #bac

        data = [ [0., 1.], [1., 0.], [1., 1.] ]

        interp = Interpolation(points, vertices, data)

        assert allclose(interp.get_A(), [[0.5, 0.5, 0.0],  #Affects vertex 1 and 0
                                   [0.5, 0.0, 0.5],  #Affects vertex 0 and 2
                                   [0.0, 0.5, 0.5]]) #Affects vertex 1 and 2


    def test_datapoints_on_edges(self):
        """Try datapoints on edges -
        each point should affect two matrix entries in proportion
        """

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        points = [a, b, c]
        vertices = [ [1,0,2] ]   #bac

        data = [ [0., 1.5], [1.5, 0.], [1.5, 0.5] ]

        interp = Interpolation(points, vertices, data)

        assert allclose(interp.get_A(), [[0.25, 0.75, 0.0],  #Affects vertex 1 and 0
                                   [0.25, 0.0, 0.75],  #Affects vertex 0 and 2
                                   [0.0, 0.25, 0.75]]) #Affects vertex 1 and 2

    def test_arbitrary_datapoints(self):
        """Try arbitrary datapoints
        """

        from Numeric import sum

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        points = [a, b, c]
        vertices = [ [1,0,2] ]   #bac

        data = [ [0.2, 1.5], [0.123, 1.768], [1.43, 0.44] ]

        interp = Interpolation(points, vertices, data)
        #print "interp.get_A()", interp.get_A()
        assert allclose(sum(interp.get_A(), axis=1), 1.0)

    def test_arbitrary_datapoints_some_outside(self):
        """Try arbitrary datapoints one outside the triangle.
        That one should be ignored
        """

        from Numeric import sum

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        points = [a, b, c]
        vertices = [ [1,0,2] ]   #bac

        data = [ [0.2, 1.5], [0.123, 1.768], [1.43, 0.44], [5.0, 7.0]]


        interp = Interpolation(points, vertices, data, precrop = True)
        assert allclose(sum(interp.get_A(), axis=1), 1.0)

        interp = Interpolation(points, vertices, data, precrop = False)
        assert allclose(sum(interp.get_A(), axis=1), [1,1,1,0])



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
        data_points = [ [-3., 2.0], [-2, 1], [0.0, 1], [0, 3], [2, 3], [-1.0/3,-4./3] ]
        interp = Interpolation(points, triangles, data_points)

        answer = [[0.0, 0.0, 0.0, 1.0, 0.0, 0.0],    #Affects point d
                  [0.5, 0.0, 0.0, 0.5, 0.0, 0.0],    #Affects points a and d
                  [0.75, 0.25, 0.0, 0.0, 0.0, 0.0],  #Affects points a and b
                  [0.0, 0.5, 0.0, 0.5, 0.0, 0.0],    #Affects points a and d
                  [0.25, 0.75, 0.0, 0.0, 0.0, 0.0],  #Affects points a and b
                  [1./3, 0.0, 0.0, 0.0, 1./3, 1./3]] #Affects points a, e and f


        A = interp.get_A()
        for i in range(A.shape[0]):
            for j in range(A.shape[1]):
                if not allclose(A[i,j], answer[i][j]):
                    print i,j,':',A[i,j], answer[i][j]


        assert allclose(interp.get_A(), answer)




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

        interp = Interpolation(points, triangles, data_coords, alpha=5.0e-20)
        z = [z1, z2, z3]
        f = interp.fit(z)
        answer = [0, 5., 5.]

        #print "f\n",f
        #print "answer\n",answer

        assert allclose(f, answer, atol=1e-7)


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
        interp = Interpolation(points, triangles, data_coords, alpha=0.0)
        f = interp.fit(z)
        answer = linear_function(points)

        assert allclose(f, answer)

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
        interp = Interpolation(vertices, triangles, point_coords, alpha=0.0)

        #print 'z',z
        f = interp.fit(z)
        answer = linear_function(vertices)
        #print "f\n",f
        #print "answer\n",answer
        assert allclose(f, answer)


    def test_smooth_attributes_to_meshIV(self):
        """ Testing 2 attributes smoothed to the mesh
        """

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

        interp = Interpolation(points, triangles, data_coords, alpha=0.0)
        z = [z1, z2, z3]
        f =  interp.fit_points(z)
        answer = [[0,0], [5., 10.], [5., 10.]]
        assert allclose(f, answer)

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

        interp = Interpolation(vertices, triangles, point_coords)
        f = linear_function(vertices)
        #z = interp.interpolate(f)
        answer = linear_function(point_coords)


        #assert allclose(z, answer)


    def test_interpolate_attributes_to_points_interp_only(self):
        v0 = [0.0, 0.0]
        v1 = [0.0, 5.0]
        v2 = [5.0, 0.0]

        vertices = [v0, v1, v2]
        triangles = [ [1,0,2] ]   #bac

        d0 = [1.0, 1.0]
        d1 = [1.0, 2.0]
        d2 = [3.0, 1.0]
        point_coords = [ d0, d1, d2]

        interp = Interpolation(vertices, triangles, point_coords,
                               interp_only = True)
        
        f = linear_function(vertices)
        #z = interp.interpolate(f)
        answer = linear_function(point_coords)
        #print "answer", answer
        #print "z", z 

        #assert allclose(z, answer)

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

        interp = Interpolation(vertices, triangles, point_coords)
        f = linear_function(vertices)
        #z = interp.interpolate(f)
        answer = linear_function(point_coords)
        #print "z",z
        #print "answer",answer
        #assert allclose(z, answer)

    def test_interpolate_attributes_to_pointsIII(self):
        """Test linear interpolation of known values at vertices to
        new points inside a triangle
        """
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

        interp = Interpolation(vertices, triangles, point_coords)

        #Known values at vertices
        #Functions are x+y, x+2y, 2x+y, x-y-5
        f = [ [0., 0., 0., -5.],        # (0,0)
              [5., 10., 5., -10.],      # (0,5)
              [5., 5., 10.0, 0.],       # (5,0)
              [10., 15., 15., -5.]]     # (5,5)

        #z = interp.interpolate(f)
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

        #Should an error message be returned if points are outside
        # of the mesh? Not currently.  

        #assert allclose(z, answer)


    def test_interpolate_point_outside_of_mesh(self):
        """Test linear interpolation of known values at vertices to
        new points inside a triangle
        """
        a = [0.0, 0.0]
        b = [0.0, 5.0]
        c = [5.0, 0.0]
        d = [5.0, 5.0]

        vertices = [a, b, c, d]
        triangles = [ [1,0,2], [2,3,1] ]   #bac, cdb

        #Far away point
        d7 = [-1., -1.]
        
        point_coords = [ d7]

        interp = Interpolation(vertices, triangles, point_coords)

        #Known values at vertices
        #Functions are x+y, x+2y, 2x+y, x-y-5
        f = [ [0., 0., 0., -5.],        # (0,0)
              [5., 10., 5., -10.],      # (0,5)
              [5., 5., 10.0, 0.],       # (5,0)
              [10., 15., 15., -5.]]     # (5,5)


        #z = interp.interpolate(f)
        #answer = [ [0., 0., 0., 0.]] # (-1,-1)

        #print "***********"
        #print "z",z
        #print "answer",answer
        #print "***********"

        #Should an error message be returned if points are outside
        # of the mesh? Not currently.  

        #assert allclose(z, answer)

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

        interp = Interpolation(vertices, triangles, point_coords)
        f = array([linear_function(vertices),2*linear_function(vertices) ])
        f = transpose(f)
        #print "f",f
        #z = interp.interpolate(f)
        answer = [linear_function(point_coords),
                  2*linear_function(point_coords) ]
        answer = transpose(answer)
        #print "z",z
        #print "answer",answer
        #assert allclose(z, answer)

    def test_smooth_attributes_to_mesh_function(self):
        """ Testing 2 attributes smoothed to the mesh
        """

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

        f = fit_to_mesh(points, triangles, data_coords, z, alpha=0.0)
        answer = [[0, 0], [5., 10.], [5., 10.]]

        assert allclose(f, answer)



    def test_pts2rectangular(self):

    	import time, os
        FN = 'xyatest' + str(time.time()) + '.xya'
    	fid = open(FN, 'w')
	fid.write(' %s \n' %('elevation'))
	fid.write('%f %f %f\n' %(1,1,2) )
	fid.write('%f %f %f\n' %(1,3,4) )
	fid.write('%f %f %f\n' %(3,1,4) )
	fid.close()

	points, triangles, boundary, attributes =\
                pts2rectangular(FN, 4, 4)


        data_coords = [ [1,1], [1,3], [3,1] ]
        z = [2, 4, 4]

        ref = fit_to_mesh(points, triangles, data_coords, z, verbose=False)

        #print attributes
        #print ref
	assert allclose(attributes, ref)

        os.remove(FN)


    #Tests of smoothing matrix
    def test_smoothing_matrix_one_triangle(self):
        from Numeric import dot
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        points = [a, b, c]

        vertices = [ [1,0,2] ]   #bac

        interp = Interpolation(points, vertices)

        assert allclose(interp.get_D(), [[1, -0.5, -0.5],
                                   [-0.5, 0.5, 0],
                                   [-0.5, 0, 0.5]])

        #Define f(x,y) = x
        f = array([0,0,2]) #Value at global vertex 2

        #Check that int (df/dx)**2 + (df/dy)**2 dx dy =
        #           int 1 dx dy = area = 2
        assert dot(dot(f, interp.get_D()), f) == 2

        #Define f(x,y) = y
        f = array([0,2,0])  #Value at global vertex 1

        #Check that int (df/dx)**2 + (df/dy)**2 dx dy =
        #           int 1 dx dy = area = 2
        assert dot(dot(f, interp.get_D()), f) == 2

        #Define f(x,y) = x+y
        f = array([0,2,2])  #Values at global vertex 1 and 2

        #Check that int (df/dx)**2 + (df/dy)**2 dx dy =
        #           int 2 dx dy = 2*area = 4
        assert dot(dot(f, interp.get_D()), f) == 4



    def test_smoothing_matrix_more_triangles(self):
        from Numeric import dot

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe, daf, dae
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        interp = Interpolation(points, vertices)


        #assert allclose(interp.get_D(), [[1, -0.5, -0.5],
        #                           [-0.5, 0.5, 0],
        #                           [-0.5, 0, 0.5]])

        #Define f(x,y) = x
        f = array([0,0,2,0,2,4]) #f evaluated at points a-f

        #Check that int (df/dx)**2 + (df/dy)**2 dx dy =
        #           int 1 dx dy = total area = 8
        assert dot(dot(f, interp.get_D()), f) == 8

        #Define f(x,y) = y
        f = array([0,2,0,4,2,0]) #f evaluated at points a-f

        #Check that int (df/dx)**2 + (df/dy)**2 dx dy =
        #           int 1 dx dy = area = 8
        assert dot(dot(f, interp.get_D()), f) == 8

        #Define f(x,y) = x+y
        f = array([0,2,2,4,4,4])  #f evaluated at points a-f

        #Check that int (df/dx)**2 + (df/dy)**2 dx dy =
        #           int 2 dx dy = 2*area = 16
        assert dot(dot(f, interp.get_D()), f) == 16


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

        interp = Interpolation(points, triangles, data_points, alpha=0.0)

        z = linear_function(data_points)
        answer = linear_function(points)

        f = interp.fit(z)

        #print "f",f
        #print "answer",answer
        assert allclose(f, answer)

        #Map back
        #z1 = interp.interpolate(f)
        #print "z1\n", z1
        #print "z\n",z
        #assert allclose(z, z1)


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
        interp = Interpolation(points, triangles, data_points, alpha=0.0)
        #Must raise an exception
        try:
            f = interp.fit(z)
        except:
            pass

        #Now try with smoothing parameter
        interp = Interpolation(points, triangles, data_points, alpha=1.0e-13)

        f = interp.fit(z)
        #f will be different from answer due to smoothing
        assert allclose(f, answer,atol=5)

        #Map back
        #z1 = interp.interpolate(f)
        #assert allclose(z, z1)



    def test_fit_and_interpolation_with_new_points(self):
        """Fit a surface to one set of points. Then interpolate that surface
        using another set of points.
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
                        [ 15, -17],   #Outside mesh
                        [ 1.0, 2.0],
                        [ 1.0, 3.0],
                        [ 2.0, 1.0],
                        [ 3.0, 0.0],
                        [ 3.0, 1.0]]

        #Fit surface to mesh
        interp = Interpolation(points, triangles, data_points1, alpha=0.0,
                               precrop = True, verbose=False)
        z = linear_function(data_points1) #Example z-values
        f = interp.fit(z)                 #Fitted values at vertices



        #New datapoints where interpolated values are sought
        data_points2 = [[ 0.0, 0.0],
                        [ 0.5, 0.5],
                        [ 0.7, 0.7],
                        [-13, 65],  #Outside
                        [ 1.0, 0.5],
                        [ 2.0, 0.4],
                        [ 2.8, 1.2]]



        #Build new A matrix based on new points (without precrop)
        interp.build_interpolation_matrix_A(data_points2, precrop = False)

        #Interpolate using fitted surface
        #z1 = interp.interpolate(f)

        #import Numeric
        #data_points2 = Numeric.take(data_points2, interp.point_indices)

        #Desired result (OK for points inside)

        answer = linear_function(data_points2)
        import Numeric
        #z1 = Numeric.take(z1, [0,1,2,4,5,6])
        answer = Numeric.take(answer, [0,1,2,4,5,6])
        #assert allclose(z1, answer)

        #Build new A matrix based on new points (with precrop)
        interp.build_interpolation_matrix_A(data_points2, precrop = True)

        #Interpolate using fitted surface
        #z1 = interp.interpolate(f)

        import Numeric
        data_points2 = Numeric.take(data_points2, interp.point_indices)

        #Desired result
        answer = linear_function(data_points2)
        #assert allclose(z1, answer)






    def test_interpolation_from_discontinuous_vertex_values(self):
        """test_interpolation_from_discontinuous_vertex_values.
        This will test the format used internally in pyvolution and also
        interpolation from sww files
        """
        
        #Setup mesh used to represent discontinuous function
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [b, a, c,
                  b, c, e,
                  e, c, f,
                  d, b, e]
        
        #bac, bce, ecf, dbe
        triangles = [[0,1,2], [3,4,5], [6,7,8], [9,10,11]]

        
        vertex_values = [0.,0.,0.,1.,1.,1.,2.,2.,2.,7.,3.,3.]
                          
        

        #New datapoints where interpolated values are sought
        data_points = [[0.0, 0.0],  #T0
                       [0.5, 0.5],  #T0
                       [1.5, 1.5],  #T1
                       [2.5, 0.5],  #T2
                       [0.0, 3.0],  #T3
                       [1.0, 2.0],  #In between T1 and T3 (T1 is used) FIXME?
                       [2.0, 1.0],  #In between T1 and T2 (T1 is used) FIXME?
                       [1.0, 1.0]]  #In between T1 and T0 (T0 is used) FIXME?
        



        #Build interpolation matrix
        interp = Interpolation(points, triangles, data_points)
                               #, alpha=0.0, precrop = True)

        #print interp.A.todense()
        #print vertex_values

        #Interpolate using fitted surface
        #z = interp.interpolate(vertex_values)

        #print z

        #assert allclose(z, [0,0,1,2,5,1,1,0])




    def test_interpolation_function_time_only(self):
        """Test spatio-temporal interpolation
        Test that spatio temporal function performs the correct
        interpolations in both time and space
        """


        #Three timesteps
        time = [1.0, 5.0, 6.0]
        

        #One quantity
        Q = zeros( (3,6), Float )

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
        from anuga.utilities.numerical_tools import mean        
        I = Interpolation_function(time, [mean(Q[0,:]),
                                          mean(Q[1,:]),
                                          mean(Q[2,:])])



        #Check temporal interpolation
        for i in [0,1,2]:
            assert allclose(I(time[i]), mean(Q[i,:]))

        #Midway    
        assert allclose(I( (time[0] + time[1])/2 ),
                        (I(time[0]) + I(time[1]))/2 )

        assert allclose(I( (time[1] + time[2])/2 ),
                        (I(time[1]) + I(time[2]))/2 )

        assert allclose(I( (time[0] + time[2])/2 ),
                        (I(time[0]) + I(time[2]))/2 )                 

        #1/3
        assert allclose(I( (time[0] + time[2])/3 ),
                        (I(time[0]) + I(time[2]))/3 )                         


        #Out of bounds checks
        try:
            I(time[0]-1) 
        except:
            pass
        else:
            raise 'Should raise exception'

        try:
            I(time[-1]+1) 
        except:
            pass
        else:
            raise 'Should raise exception'        


        

    def interpolation_test_interpolation_function_spatial_only(self):
        """Test spatio-temporal interpolation with constant time
        """

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
                assert allclose(I(t, id), answer[id])

            t += 0.1    


        try:    
            I(1)
        except:
            pass
        else:
            raise 'Should raise exception'
            


    def interpolation_test_interpolation_function(self):
        """Test spatio-temporal interpolation
        Test that spatio temporal function performs the correct
        interpolations in both time and space
        """


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
        Q = zeros( (3,6), Float )

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
                assert allclose(I(t, id), t*answer[id])

            t += 0.1    


        try:    
            I(1)
        except:
            pass
        else:
            raise 'Should raise exception'
            
        #
        #interpolation_points = [[ 0.0, 0.0],
        #                        [ 0.5, 0.5],
        #                        [ 0.7, 0.7],
        #                        [-13, 65],  #Outside
        #                        [ 1.0, 0.5],
        #                        [ 2.0, 0.4],
        #                        [ 2.8, 1.2]]
        #
        #try:
        #    I = Interpolation_function(time, Q,
        #                               vertex_coordinates = points,
        #                               triangles = triangles, 
        #                               interpolation_points = interpolation_points,
        #                               verbose = False)
        #except:
        #    pass
        #else:
        #    raise 'Should raise exception'





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
        interp = Interpolation(points, triangles, data_points1,
                               alpha=0.0,
                               data_origin = data_origin,
                               mesh_origin = mesh_origin)

        z = linear_function(data_points1) #Example z-values
        f = interp.fit(z)                 #Fitted values at vertices


        #New datapoints where interpolated values are sought
        data_points2 = [[ 0.0, 0.0],
                        [ 0.5, 0.5],
                        [ 0.7, 0.7],
                        [ 1.0, 0.5],
                        [ 2.0, 0.4],
                        [ 2.8, 1.2]]


        #Build new A matrix based on new points
        interp.build_interpolation_matrix_A(data_points2)

        #Interpolate using fitted surface
        #z1 = interp.interpolate(f)

        #Desired result
        #answer = linear_function(data_points2)
        #assert allclose(z1, answer)


        ##############################################

        #Then check situation where points are relative to a different
        #origin (same zone, though, until we figure that out (FIXME))

        mesh_origin = (56, 290000, 618000) #zone, easting, northing
        data_origin = (56, 10000, 10000) #zone, easting, northing

        #Shift datapoints according to new origin

        for k in range(len(data_points1)):
            data_points1[k][0] += mesh_origin[1] - data_origin[1]
            data_points1[k][1] += mesh_origin[2] - data_origin[2]

        for k in range(len(data_points2)):
            data_points2[k][0] += mesh_origin[1] - data_origin[1]
            data_points2[k][1] += mesh_origin[2] - data_origin[2]



        #Fit surface to mesh
        interp = Interpolation(points, triangles, data_points1,
                               alpha=0.0,
                               data_origin = data_origin,
                               mesh_origin = mesh_origin)

        f1 = interp.fit(z) #Fitted values at vertices (using same z as before)

        assert allclose(f,f1), 'Fit should have been unaltered'


        #Build new A matrix based on new points
        interp.build_interpolation_matrix_A(data_points2)

        #Interpolate using fitted surface
        #z1 = interp.interpolate(f)
        #assert allclose(z1, answer)


        #########################################################
        #Finally try to relate data_points2 to new origin without
        #rebuilding matrix

        data_origin = (56, 2000, 2000) #zone, easting, northing
        for k in range(len(data_points2)):
            data_points2[k][0] += 8000
            data_points2[k][1] += 8000

        #Build new A matrix based on new points
        interp.build_interpolation_matrix_A(data_points2,
                                            data_origin = data_origin)

        #Interpolate using fitted surface
        #z1 = interp.interpolate(f)
        #assert allclose(z1, answer)



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


        #Data                       
        data_points = [[ 201.0, 401.0],
                       [ 201.0, 403.0],
                       [ 203.0, 401.0]]

        z = [2, 4, 4]
        
        data_geo = Geo_reference(56,-200,-400)
        
        #Fit
        zz = fit_to_mesh(vertex_coordinates, triangles, data_points, z,
                         data_origin = data_geo.get_origin(),
                         mesh_origin = mesh_geo.get_origin(),
                         alpha = 0)
        assert allclose( zz, [0,5,5] )


    def test_fit_to_mesh_file(self):
        from load_mesh.loadASCII import import_mesh_file, \
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

        # create an .xya file
        point_file = tempfile.mktemp(".xya")
        fd = open(point_file,'w')
        fd.write("elevation, stage \n 1.0, 1.0,2.,4 \n 1.0, 3.0,4,8 \n 3.0,1.0,4.,8 \n")
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
        assert allclose(mesh_dic['vertex_attributes'],ans)

        self.failUnless(mesh_dic['vertex_attribute_titles']  ==
                        ['elevation','stage'],
                        'test_fit_to_mesh_file failed')

        #clean up
        os.remove(mesh_file)
        os.remove(point_file)
        os.remove(mesh_output_file)

    def test_fit_to_mesh_file3(self):
        from load_mesh.loadASCII import import_mesh_file, \
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

        # create an .xya file
        point_file = tempfile.mktemp(".xya")
        fd = open(point_file,'w')
        fd.write("elevation, stage \n 1.0, 1.0,2.,4 \n 1.0, 3.0,4,8 \n 3.0,1.0,4.,8 \n")
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
        assert allclose(mesh_dic['vertex_attributes'],ans)

        self.failUnless(mesh_dic['vertex_attribute_titles']  ==
                        ['elevation','stage'],
                        'test_fit_to_mesh_file failed')

        #clean up
        os.remove(mesh_file)
        os.remove(point_file)
        os.remove(mesh_output_file)

    def test_fit_to_mesh_file4(self):
        from load_mesh.loadASCII import import_mesh_file, \
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

        geo_ref = Geo_reference(56,-200,-400)
        # create an .xya file
        point_file = tempfile.mktemp(".xya")
        fd = open(point_file,'w')
        fd.write("elevation, stage \n 201.0, 401.0,2.,4 \n 201.0, 403.0,4,8 \n 203.0, 401.0,4.,8 \n")
        geo_ref.write_ASCII(fd)
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
              [5.0, 10.0]]
        assert allclose(mesh_dic['vertex_attributes'],ans)

        self.failUnless(mesh_dic['vertex_attribute_titles']  ==
                        ['elevation','stage'],
                        'test_fit_to_mesh_file failed')

        #clean up
        os.remove(mesh_file)
        os.remove(point_file)
        os.remove(mesh_output_file)

    def test_fit_to_mesh_fileII(self):
        from load_mesh.loadASCII import import_mesh_file, \
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

        # create an .xya file
        point_file = tempfile.mktemp(".xya")
        fd = open(point_file,'w')
        fd.write("elevation, stage \n 1.0, 1.0,2.,4 \n 1.0, 3.0,4,8 \n 3.0,1.0,4.,8 \n")
        fd.close()

        mesh_output_file = "new_triangle.tsh"
        fit_to_mesh_file(mesh_file,
                         point_file,
                         mesh_output_file,
                         alpha = 0.0)
        # load in the .tsh file we just wrote
        mesh_dic = import_mesh_file(mesh_output_file)

        assert allclose(mesh_dic['vertex_attributes'],
                        [[1.0, 2.0,0.0, 0.0],
                         [1.0, 2.0,5.0, 10.0],
                         [1.0, 2.0,5.0,10.0]])

        self.failUnless(mesh_dic['vertex_attribute_titles']  ==
                        ['density', 'temp','elevation','stage'],
                        'test_fit_to_mesh_file failed')

        #clean up
        os.remove(mesh_file)
        os.remove(mesh_output_file)
        os.remove(point_file)

    def test_fit_to_mesh_file_errors(self):
        from load_mesh.loadASCII import import_mesh_file, export_mesh_file
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

        # create an .xya file
        point_file = tempfile.mktemp(".xya")
        fd = open(point_file,'w')
        fd.write("elevation stage \n 1.0, 1.0,2.,4 \n 1.0, 3.0,4,8 \n 3.0,1.0,4.,8 \n")
        fd.close()

        mesh_output_file = "new_triangle.tsh"
        try:
            fit_to_mesh_file(mesh_file, point_file,
                             mesh_output_file, display_errors = False)
        except IOError:
            pass
        else:
            #self.failUnless(0 ==1,  'Bad file did not raise error!')
            raise 'Bad file did not raise error!'
            
        #clean up
        os.remove(mesh_file)
        os.remove(point_file)

    def test_fit_to_mesh_file_errorsII(self):
        from load_mesh.loadASCII import import_mesh_file, export_mesh_file
        import tempfile
        import os

        # create a .tsh file, no user outline
        mesh_file = tempfile.mktemp(".tsh")
        fd = open(mesh_file,'w')
        fd.write("unit testing a bad .tsh file \n")
        fd.close()

        # create an .xya file
        point_file = tempfile.mktemp(".xya")
        fd = open(point_file,'w')
        fd.write("elevation, stage \n 1.0, 1.0,2.,4 \n 1.0, 3.0,4,8 \n 3.0,1.0,4.,8 \n")
        fd.close()

        mesh_output_file = "new_triangle.tsh"
        try:
            fit_to_mesh_file(mesh_file, point_file,
                             mesh_output_file, display_errors = False)
        except IOError:
            pass
        else:
            raise 'Bad file did not raise error!'
            
        #clean up
        os.remove(mesh_file)
        os.remove(point_file)

    def test_fit_to_mesh_file_errorsIII(self):
        from load_mesh.loadASCII import import_mesh_file, export_mesh_file
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

        # create an .xya file
        point_file = tempfile.mktemp(".xya")
        fd = open(point_file,'w')
        fd.write("elevation, stage \n 1.0, 1.0,2.,4 \n 1.0, 3.0,4,8 \n 3.0,1.0,4.,8 \n")
        fd.close()

        #This a deliberately illegal filename to invoke the error.
        mesh_output_file = ".../\z\z:ya.tsh"        

        try:
            fit_to_mesh_file(mesh_file, point_file,
                             mesh_output_file, display_errors = False)
        except IOError:
            pass
        else:
            raise 'Bad file did not raise error!'
        
        #clean up
        os.remove(mesh_file)
        os.remove(point_file)
  
## FIXME?  Running from the Comand line isn't in vogue these days
#  The test was breaking when test_all at the inundation level was running
# was running it.issue - not running the test in this directory
    def Bad_test_fit_to_mesh_file_errorsIV(self):
        import os
        command = '%s least_squares.py q q q e n 0.9 n'  %(sys.executable)
        status = os.system(command) 
        self.failUnless(status%255  == 1,
                        'command prompt least_squares.py failed.  Incorect exit status.')
        
    def test_fit_to_msh_netcdf_fileII(self):
        from load_mesh.loadASCII import import_mesh_file, export_mesh_file
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
        mesh_file = tempfile.mktemp(".msh")
        export_mesh_file(mesh_file,mesh_dic)

        # create an .xya file
        point_file = tempfile.mktemp(".xya")
        fd = open(point_file,'w')
        fd.write("elevation, stage \n 1.0, 1.0,2.,4 \n 1.0, 3.0,4,8 \n 3.0,1.0,4.,8 \n")
        fd.close()

        mesh_output_file = "new_triangle.msh"
        fit_to_mesh_file(mesh_file,
                         point_file,
                         mesh_output_file,
                         alpha = 0.0)
        # load in the .tsh file we just wrote
        mesh_dic = import_mesh_file(mesh_output_file)

        assert allclose(mesh_dic['vertex_attributes'],
                        [[1.0, 2.0,0.0, 0.0],
                         [1.0, 2.0,5.0, 10.0],
                         [1.0, 2.0,5.0,10.0]])

        self.failUnless(mesh_dic['vertex_attribute_titles']  ==
                        ['density', 'temp','elevation','stage'],
                        'test_fit_to_mesh_file failed')

        #clean up
        os.remove(mesh_file)
        os.remove(mesh_output_file)
        os.remove(point_file)

	
	
    def test_fit_using_fit_to_mesh(self):
        """Fit a surface to one set of points. Then interpolate that surface
        using another set of points.
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
                        [ 15, -17],   #Outside mesh
                        [ 1.0, 2.0],
                        [ 1.0, 3.0],
                        [ 2.0, 1.0],
                        [ 3.0, 0.0],
                        [ 3.0, 1.0]]

        #Fit surface to mesh
        z = linear_function(data_points1) #Example z-values
	v = fit_to_mesh(points, triangles, data_points1, z, alpha=0.0,
                        precrop=True, verbose=False)

	assert allclose(linear_function(points), v)


	
    def test_acceptable_overshoot(self):
        """Fit a surface to one set of points. Then interpolate that surface
        using another set of points.
        Check that exceedance in fitted values are caught.
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
                        [ 15, -17],   #Outside mesh
                        [ 1.0, 2.0],
                        [ 1.0, 3.0],
                        [ 2.0, 1.0],
                        [ 3.0, 0.0],
                        [ 3.0, 1.0]]

        #Fit surface to mesh
        z = linear_function(data_points1) #Example z-values

        try:
            v = fit_to_mesh(points, triangles, data_points1, z, alpha=0.0,
                            acceptable_overshoot = 0.2,
                            precrop=True, verbose=False)
        except FittingError, e:
            pass
        else:
            raise 'Should have raised exception'
            

	#assert allclose(linear_function(points), v)
        

	
#-------------------------------------------------------------
if __name__ == "__main__":
    #suite = unittest.makeSuite(Test_Least_Squares,'test_smooth_attributes_to_mesh_function')
    #suite = unittest.makeSuite(Test_Least_Squares,'test_datapoint_at_centroid')
    suite = unittest.makeSuite(Test_Least_Squares,'test')

    #suite = unittest.makeSuite(Test_Least_Squares,'test_fit_to_msh_netcdf_fileII')
    #suite = unittest.makeSuite(Test_Least_Squares,'test_fit_to_mesh_fileII')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)






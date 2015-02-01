#!/usr/bin/env python



#FIXME: Seperate the tests for mesh and general_mesh

#FIXME (Ole): Maxe this test independent of anything that inherits from General_mesh (namely shallow_water)

import unittest
from math import sqrt


from anuga.abstract_2d_finite_volumes.neighbour_mesh import *
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_periodic
from anuga.config import epsilon

from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.geometry.polygon import is_inside_polygon
from anuga.utilities.numerical_tools import ensure_numeric

import numpy as num


def distance(x, y):
    return sqrt(num.sum((num.array(x)-num.array(y))**2))

class Test_Mesh(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_triangle_inputs(self):
        points = [[0.0, 0.0], [4.0, 0.0], [0.0, 3.0]]
        vertices = [0,1,2] #Wrong

        try:
            mesh = Mesh(points, vertices)
        except:
            pass
        else:
            msg = 'Should have raised exception'
            raise Exception(msg)


    def test_basic_triangle(self):

        a = [0.0, 0.0]
        b = [4.0, 0.0]
        c = [0.0, 3.0]

        points = [a, b, c]
        vertices = [[0,1,2]]
        mesh = Mesh(points, vertices)

        #Centroid
        centroid = mesh.centroid_coordinates[0]
        assert centroid[0] == 4.0/3
        assert centroid[1] == 1.0

        #Area
        assert mesh.areas[0] == 6.0,\
               'Area was %f, should have been 6.0' %mesh.areas[0]

        #Normals
        normals = mesh.get_normals()
        assert num.allclose(normals[0, 0:2], [3.0/5, 4.0/5])
        assert num.allclose(normals[0, 2:4], [-1.0, 0.0])
        assert num.allclose(normals[0, 4:6], [0.0, -1.0])

        assert num.allclose(mesh.get_normal(0,0), [3.0/5, 4.0/5])
        assert num.allclose(mesh.get_normal(0,1), [-1.0, 0.0])
        assert num.allclose(mesh.get_normal(0,2), [0.0, -1.0])

        #Edge lengths
        assert num.allclose(mesh.edgelengths[0], [5.0, 3.0, 4.0])


        #Vertex coordinates
        #V = mesh.get_vertex_coordinates()
        #assert allclose(V[0], [0.0, 0.0, 4.0, 0.0, 0.0, 3.0])
        

        V = mesh.get_vertex_coordinates()
        assert num.allclose(V, [ [0.0, 0.0],
                             [4.0, 0.0],
                             [0.0, 3.0] ])

        V0 = mesh.get_vertex_coordinate(0, 0)
        assert num.allclose(V0, [0.0, 0.0])

        V1 = mesh.get_vertex_coordinate(0, 1)
        assert num.allclose(V1, [4.0, 0.0])

        V2 = mesh.get_vertex_coordinate(0, 2)
        assert num.allclose(V2, [0.0, 3.0])


        #General tests:

        #Test that points are arranged in a counter clock wise order etc
        mesh.check_integrity()


        #Test that the centroid is located 2/3 of the way
        #from each vertex to the midpoint of the opposite side

        V = mesh.get_vertex_coordinates()
        x0 = V[0, 0]; y0 = V[0, 1]
        x1 = V[1, 0]; y1 = V[1, 1]
        x2 = V[2, 0]; y2 = V[2, 1]
        #x0 = V[0,0]
        #y0 = V[0,1]
        #x1 = V[0,2]
        #y1 = V[0,3]
        #x2 = V[0,4]
        #y2 = V[0,5]

        m0 = [(x1 + x2)/2, (y1 + y2)/2]
        m1 = [(x0 + x2)/2, (y0 + y2)/2]
        m2 = [(x1 + x0)/2, (y1 + y0)/2]

        d0 = distance(centroid, [x0, y0])
        d1 = distance(m0, [x0, y0])
        assert d0 == 2*d1/3
        #
        d0 = distance(centroid, [x1, y1])
        d1 = distance(m1, [x1, y1])
        assert abs(d0 - 2*d1/3) < epsilon, '%e, %e' %(d0, 2*d1/3)

        d0 = distance(centroid, [x2, y2])
        d1 = distance(m2, [x2, y2])
        assert abs(d0 - 2*d1/3) < epsilon, '%e, %e' %(d0, 2*d1/3)

        #Radius
        d0 = distance(centroid, m0)
        assert d0 == 5.0/6

        d1 = distance(centroid, m1)
        assert d1 == sqrt(73.0/36)

        d2 = distance(centroid, m2)
        assert d2 == sqrt(13.0/9)

        assert mesh.radii[0] == min(d0, d1, d2)
        assert mesh.radii[0] == 5.0/6


        #Let x be the centroid of triangle abc.
        #Test that areas of the three triangles axc, cxb, and bxa are equal.
        points = [a, b, c, centroid]
        vertices = [[0,3,2], [2,3,1], [1,3,0]]
        new_mesh = Mesh(points, vertices)

        assert new_mesh.areas[0] == new_mesh.areas[1]
        assert new_mesh.areas[1] == new_mesh.areas[2]
        assert new_mesh.areas[1] == new_mesh.areas[2]

        assert new_mesh.areas[1] == mesh.areas[0]/3



    def test_general_triangle(self):
        a = [2.0, 1.0]
        b = [6.0, 2.0]
        c = [1.0, 3.0]

        points = [a, b, c]
        vertices = [[0,1,2]]

        mesh = Mesh(points, vertices)
        centroid = mesh.centroid_coordinates[0]


        #Test that the centroid is located 2/3 of the way
        #from each vertex to the midpoint of the opposite side

        V = mesh.get_vertex_coordinates()
        x0 = V[0, 0]; y0 = V[0, 1]
        x1 = V[1, 0]; y1 = V[1, 1]
        x2 = V[2, 0]; y2 = V[2, 1]        

        #x0 = V[0,0]
        #y0 = V[0,1]
        #x1 = V[0,2]
        #y1 = V[0,3]
        #x2 = V[0,4]
        #y2 = V[0,5]

        m0 = [(x1 + x2)/2, (y1 + y2)/2]
        m1 = [(x0 + x2)/2, (y0 + y2)/2]
        m2 = [(x1 + x0)/2, (y1 + y0)/2]

        d0 = distance(centroid, [x0, y0])
        d1 = distance(m0, [x0, y0])
        assert abs(d0 - 2*d1/3) < epsilon, '%e, %e' %(d0, 2*d1/3)
        #
        d0 = distance(centroid, [x1, y1])
        d1 = distance(m1, [x1, y1])
        assert abs(d0 - 2*d1/3) < epsilon, '%e, %e' %(d0, 2*d1/3)

        d0 = distance(centroid, [x2, y2])
        d1 = distance(m2, [x2, y2])
        assert abs(d0 - 2*d1/3) < epsilon, '%e, %e' %(d0, 2*d1/3)

        #Radius
        d0 = distance(centroid, m0)
        d1 = distance(centroid, m1)
        d2 = distance(centroid, m2)
        assert mesh.radii[0] == min(d0, d1, d2)



        #Let x be the centroid of triangle abc.
        #Test that areas of the three triangles axc, cxb, and bxa are equal.

        points = [a, b, c, centroid]
        vertices = [[0,3,2], [2,3,1], [1,3,0]]
        new_mesh = Mesh(points, vertices)

        assert new_mesh.areas[0] == new_mesh.areas[1]
        assert new_mesh.areas[1] == new_mesh.areas[2]
        assert new_mesh.areas[1] == new_mesh.areas[2]

        assert new_mesh.areas[1] == mesh.areas[0]/3


        #Test that points are arranged in a counter clock wise order
        mesh.check_integrity()

    def test_inscribed_circle_equilateral(self):
        """test that the radius is calculated correctly by mesh in the case of an equilateral triangle"""
        a = [0.0, 0.0]
        b = [2.0, 0.0]
        c = [1.0, sqrt(3.0)]

        points = [a, b, c]
        vertices = [[0,1,2]]

        mesh = Mesh(points, vertices,use_inscribed_circle=False)
        assert num.allclose(mesh.radii[0],sqrt(3.0)/3),'Steve''s doesn''t work'

        mesh = Mesh(points, vertices,use_inscribed_circle=True)
        assert num.allclose(mesh.radii[0],sqrt(3.0)/3),'inscribed circle doesn''t work'

    def test_inscribed_circle_rightangle_triangle(self):
        """test that the radius is calculated correctly by mesh in the case of a right-angled triangle"""
        a = [0.0, 0.0]
        b = [4.0, 0.0]
        c = [0.0, 3.0]

        points = [a, b, c]
        vertices = [[0,1,2]]

        mesh = Mesh(points, vertices,use_inscribed_circle=False)
        assert num.allclose(mesh.radii[0],5.0/6),'Steve''s doesn''t work'

        mesh = Mesh(points, vertices,use_inscribed_circle=True)
        assert num.allclose(mesh.radii[0],1.0),'inscribed circle doesn''t work'


    def test_two_triangles(self):
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        e = [2.0, 2.0]
        points = [a, b, c, e]
        vertices = [ [1,0,2], [1,2,3] ]   #bac, bce
        mesh = Mesh(points, vertices)

        assert mesh.areas[0] == 2.0

        assert num.allclose(mesh.centroid_coordinates[0], [2.0/3, 2.0/3])


        #Test that points are arranged in a counter clock wise order
        mesh.check_integrity()



    def test_more_triangles(self):

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe, daf, dae
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]
        mesh = Mesh(points, vertices)

        #Test that points are arranged in a counter clock wise order
        mesh.check_integrity()

        assert mesh.areas[0] == 2.0
        assert mesh.areas[1] == 2.0
        assert mesh.areas[2] == 2.0
        assert mesh.areas[3] == 2.0

        assert mesh.edgelengths[1,0] == 2.0
        assert mesh.edgelengths[1,1] == 2.0
        assert mesh.edgelengths[1,2] == sqrt(8.0)

        assert num.allclose(mesh.centroid_coordinates[0], [2.0/3, 2.0/3])
        assert num.allclose(mesh.centroid_coordinates[1], [4.0/3, 4.0/3])
        assert num.allclose(mesh.centroid_coordinates[2], [8.0/3, 2.0/3])
        assert num.allclose(mesh.centroid_coordinates[3], [2.0/3, 8.0/3])

    def test_mesh_and_neighbours(self):
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]


        points = [a, b, c, d, e, f]

        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4] ]
        mesh = Mesh(points, vertices)

        mesh.check_integrity()


        T = mesh
        tid = 0
        assert T.number_of_boundaries[tid] == 2
        assert T.neighbours[tid, 0] < 0  #Opposite point b (0,2)
        assert T.neighbours[tid, 1] == 1 #Opposite point a (0,0)
        assert T.neighbours[tid, 2] < 0  #Opposite point c (2,0)

        tid = 1
        assert T.number_of_boundaries[tid] == 0
        assert T.neighbours[tid, 0] == 2  #Opposite point b (0,2)
        assert T.neighbours[tid, 1] == 3  #Opposite point c (2,0)
        assert T.neighbours[tid, 2] == 0  #Opposite point e (2,2)

        tid = 2
        assert T.number_of_boundaries[tid] == 2
        assert T.neighbours[tid, 0] < 0   #Opposite point e (2,2)
        assert T.neighbours[tid, 1] < 0   #Opposite point c (2,0)
        assert T.neighbours[tid, 2] == 1  #Opposite point f (4,0)

        tid = 3
        assert T.number_of_boundaries[tid] == 2
        assert T.neighbours[tid, 0] == 1  #Opposite point d (0,4)
        assert T.neighbours[tid, 1] < 0   #Opposite point b (0,3)
        assert T.neighbours[tid, 2] < 0   #Opposite point e (2,2)

        #Neighbouring edges
        tid = 0
        assert T.neighbour_edges[tid, 0] < 0  #Opposite point b (0,2)
        assert T.neighbour_edges[tid, 1] == 2 #Opposite point a (0,0)
        assert T.neighbour_edges[tid, 2] < 0  #Opposite point c (2,0)

        tid = 1
        assert T.neighbour_edges[tid, 0] == 2 #Opposite point b (0,2)
        assert T.neighbour_edges[tid, 1] == 0 #Opposite point c (2,0)
        assert T.neighbour_edges[tid, 2] == 1 #Opposite point e (2,2)

        tid = 2
        assert T.neighbour_edges[tid, 0] < 0  #Opposite point e (2,2)
        assert T.neighbour_edges[tid, 1] < 0  #Opposite point c (2,0)
        assert T.neighbour_edges[tid, 2] == 0 #Opposite point f (4,0)

        tid = 3
        assert T.neighbour_edges[tid, 0] == 1 #Opposite point d (0,4)
        assert T.neighbour_edges[tid, 1] < 0  #Opposite point b (0,3)
        assert T.neighbour_edges[tid, 2] < 0  #Opposite point e (2,2)


    def test_build_neighbour_structure_duplicates(self):
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
                      [0, 2, 4],
                      [4, 2, 5],
                      [5, 2, 3]]
        try:
            mesh = Mesh(points, triangles)
        except:
            pass
        else:
            raise Exception("triangle edge duplicates not caught")


    def test_rectangular_mesh_basic(self):
        M=1
        N=1

        points, vertices, boundary = rectangular(M, N)
        mesh = Mesh(points, vertices, boundary)

        #Test that points are arranged in a counter clock wise order
        mesh.check_integrity()

        M=2
        N=2
        points, vertices, boundary = rectangular(M, N)
        mesh = Mesh(points, vertices, boundary)

        #Test that points are arranged in a counter clock wise order
        mesh.check_integrity()

        #assert mesh.boundary[(7,1)] == 2 # top
        assert mesh.boundary[(7,1)] == 'top' # top
        assert mesh.boundary[(3,1)] == 'top' # top






    def test_boundary_tags(self):


        points, vertices, boundary = rectangular(4, 4)
        mesh = Mesh(points, vertices, boundary)


        #Test that points are arranged in a counter clock wise order
        mesh.check_integrity()

        #print mesh.get_boundary_tags()
        #print mesh.boundary

        for k in [1, 3, 5, 7]:
            assert mesh.boundary[(k,2)] == 'left'

        for k in [24, 26, 28, 30]:
            assert mesh.boundary[(k,2)] == 'right'

        for k in [7, 15, 23, 31]:
            assert mesh.boundary[(k,1)] == 'top'
        for k in [0, 8, 16, 24]:
            assert mesh.boundary[(k,1)] == 'bottom'



    def test_rectangular_mesh(self):
        M=4
        N=16
        len1 = 100.0
        len2 = 17.0

        points, vertices, boundary = rectangular(M, N, len1, len2)
        mesh = Mesh(points, vertices, boundary)

        assert len(mesh) == 2*M*N

        for i in range(len(mesh)):
            assert mesh.areas[i] == len1*len2/(2*M*N)

            hypo = sqrt((len1/M)**2 + (len2/N)**2) #hypothenuse
            assert mesh.edgelengths[i, 0] == hypo
            assert mesh.edgelengths[i, 1] == len1/M #x direction
            assert mesh.edgelengths[i, 2] == len2/N #y direction

        #Test that points are arranged in a counter clock wise order
        mesh.check_integrity()


    def test_rectangular_mesh2(self):
        #Check that integers don't cause trouble
        N = 16

        points, vertices, boundary = rectangular(2*N, N, len1=10, len2=10)
        mesh = Mesh(points, vertices, boundary)


    def test_surrogate_neighbours(self):
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]

        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4] ]
        mesh = Mesh(points, vertices)
        mesh.check_integrity()


        T = mesh
        tid = 0
        assert T.number_of_boundaries[tid] == 2
        assert T.surrogate_neighbours[tid, 0] == tid
        assert T.surrogate_neighbours[tid, 1] == 1
        assert T.surrogate_neighbours[tid, 2] == tid

        tid = 1
        assert T.number_of_boundaries[tid] == 0
        assert T.surrogate_neighbours[tid, 0] == 2
        assert T.surrogate_neighbours[tid, 1] == 3
        assert T.surrogate_neighbours[tid, 2] == 0

        tid = 2
        assert T.number_of_boundaries[tid] == 2
        assert T.surrogate_neighbours[tid, 0] == tid
        assert T.surrogate_neighbours[tid, 1] == tid
        assert T.surrogate_neighbours[tid, 2] == 1

        tid = 3
        assert T.number_of_boundaries[tid] == 2
        assert T.surrogate_neighbours[tid, 0] == 1
        assert T.surrogate_neighbours[tid, 1] == tid
        assert T.surrogate_neighbours[tid, 2] == tid


    def test_boundary_inputs(self):
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]

        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4] ]

        boundary = { (0, 0): 'First',
                     (0, 2): 'Second',
                     (2, 0): 'Third',
                     (2, 1): 'Fourth',
                     (3, 1): 'Fifth',
                     (3, 2): 'Sixth'}


        mesh = Mesh(points, vertices, boundary)
        mesh.check_integrity()


        #Check enumeration
        #for k, (vol_id, edge_id) in enumerate(mesh.boundary_segments):
        #    b = -k-1
        #    assert mesh.neighbours[vol_id, edge_id] == b



    def test_boundary_inputs_using_one_default(self):
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]

        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4] ]

        boundary = { (0, 0): 'First',
                     (0, 2): 'Second',
                     (2, 0): 'Third',
                     (2, 1): 'Fourth',
                     #(3, 1): 'Fifth',  #Skip this
                     (3, 2): 'Sixth'}


        mesh = Mesh(points, vertices, boundary)
        mesh.check_integrity()

        from anuga.config import default_boundary_tag
        assert mesh.boundary[ (3, 1) ] == default_boundary_tag


        #Check enumeration
        #for k, (vol_id, edge_id) in enumerate(mesh.boundary_segments):
        #    b = -k-1
        #    assert mesh.neighbours[vol_id, edge_id] == b

    def test_boundary_inputs_using_all_defaults(self):
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]

        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4] ]

        boundary = { (0, 0): 'First',
                     (0, 2): 'Second',
                     (2, 0): 'Third',
                     (2, 1): 'Fourth',
                     #(3, 1): 'Fifth',  #Skip this
                     (3, 2): 'Sixth'}


        mesh = Mesh(points, vertices) #, boundary)
        mesh.check_integrity()

        from anuga.config import default_boundary_tag
        assert mesh.boundary[ (0, 0) ] == default_boundary_tag
        assert mesh.boundary[ (0, 2) ] == default_boundary_tag
        assert mesh.boundary[ (2, 0) ] == default_boundary_tag
        assert mesh.boundary[ (2, 1) ] == default_boundary_tag
        assert mesh.boundary[ (3, 1) ] == default_boundary_tag
        assert mesh.boundary[ (3, 2) ] == default_boundary_tag


        #Check enumeration
        #for k, (vol_id, edge_id) in enumerate(mesh.boundary_segments):
        #    b = -k-1
        #    assert mesh.neighbours[vol_id, edge_id] == b






    def test_inputs(self):
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]

        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4] ]

        #Too few points
        try:
            mesh = Mesh([points[0]], vertices)
        except AssertionError:
            pass
        else:
            raise Exception('Should have raised an exception')

        #Too few points - 1 element
        try:
            mesh = Mesh([points[0]], [vertices[0]])
        except AssertionError:
            pass
        else:
            raise Exception('Should have raised an exception')

        #Wrong dimension of vertices
        try:
            mesh = Mesh(points, vertices[0])
        except AssertionError:
            pass
        else:
            raise Exception('Should have raised an exception')

        #Unsubscriptable coordinates object raises exception
        try:
            mesh = Mesh(points[0], [vertices[0]])
        except AssertionError:
            pass
        else:
            raise Exception('Should have raised an exception')

        #FIXME: This has been commented out pending a decision
        #whether to allow partial boundary tags or not
        #
        #Not specifying all boundary tags
        #try:
        #    mesh = Mesh(points, vertices, {(3,0): 'x'})
        #except AssertionError:
        #    pass
        #else:
        #    raise Exception('Should have raised an exception')

        #Specifying wrong non existing segment       
        try:
            mesh = Mesh(points, vertices, {(5,0): 'x'})
        except AssertionError:
            pass
        except RuntimeError:
            pass
        else:
            raise Exception('Should have raised an exception')




    def test_internal_boundaries(self):
        """
        get values based on triangle lists.
        """
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        #Create basic mesh
        points, vertices, boundary = rectangular(1, 3)

        # Add an internal boundary
        boundary[(2,0)] = 'internal'
        boundary[(1,0)] = 'internal'

        #Create shallow water domain
        domain = Mesh(points, vertices, boundary)
        domain.build_tagged_elements_dictionary({'bottom':[0,1],
                                                 'top':[4,5],
                                                 'all':[0,1,2,3,4,5]})


    def test_boundary_polygon(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular
        #from mesh import Mesh

        #Create basic mesh
        points, vertices, boundary = rectangular(2, 2)
        mesh = Mesh(points, vertices, boundary)


        P = mesh.get_boundary_polygon()

        assert len(P) == 8
        assert num.allclose(P, [[0.0, 0.0], [0.5, 0.0], [1.0, 0.0],
                                [1.0, 0.5], [1.0, 1.0], [0.5, 1.0],
                                [0.0, 1.0], [0.0, 0.5]])
        for p in points:
            #print p, P
            assert is_inside_polygon(p, P)


    def test_boundary_polygon_II(self):

        #Points
        a = [0.0, 0.0] #0
        b = [0.0, 0.5] #1
        c = [0.0, 1.0] #2
        d = [0.5, 0.0] #3
        e = [0.5, 0.5] #4
        f = [1.0, 0.0] #5
        g = [1.0, 0.5] #6
        h = [1.0, 1.0] #7
        i = [1.5, 0.5] #8

        points = [a, b, c, d, e, f, g, h, i]

        #dea, bae, bec, fgd,
        #edg, ghe, gfi, gih
        vertices = [ [3,4,0], [1,0,4], [1,4,2], [5,6,3],
                     [4,3,6], [6,7,4], [6,5,8], [6,8,7]]

        mesh = Mesh(points, vertices)

        mesh.check_integrity()

        P = mesh.get_boundary_polygon()

        assert len(P) == 8
        assert num.allclose(P, [a, d, f, i, h, e, c, b])

        for p in points:
            #print p, P
            assert is_inside_polygon(p, P)


    def test_boundary_polygon_III(self):
        """Same as II but vertices ordered differently
        """


        #Points
        a = [0.0, 0.0] #0
        b = [0.0, 0.5] #1
        c = [0.0, 1.0] #2
        d = [0.5, 0.0] #3
        e = [0.5, 0.5] #4
        f = [1.0, 0.0] #5
        g = [1.0, 0.5] #6
        h = [1.0, 1.0] #7
        i = [1.5, 0.5] #8

        points = [a, b, c, d, e, f, g, h, i]

        #edg, ghe, gfi, gih
        #dea, bae, bec, fgd,
        vertices = [[4,3,6], [6,7,4], [6,5,8], [6,8,7],
                    [3,4,0], [1,0,4], [1,4,2], [5,6,3]]


        mesh = Mesh(points, vertices)
        mesh.check_integrity()


        P = mesh.get_boundary_polygon()

        assert len(P) == 8
        assert num.allclose(P, [a, d, f, i, h, e, c, b])

        for p in points:
            assert is_inside_polygon(p, P)

	    
    def test_boundary_polygon_IIIa(self):
        """test_boundary_polygon_IIIa - Check pathological situation where
	one triangle has no neighbours. This may be the case if a mesh
	is partitioned using pymetis.
        """


        #Points
        a = [0.0, 0.0] #0
        b = [0.0, 0.5] #1
        c = [0.0, 1.0] #2
        d = [0.5, 0.0] #3
        e = [0.5, 0.5] #4
        f = [1.0, 0.0] #5
        g = [1.0, 0.5] #6
        h = [1.0, 1.0] #7
       
	# Add pathological triangle with no neighbours to an otherwise
	# trivial mesh
	
	points = [a, b, c, d, e, f, g, h]

        #cbe, aeb, dea, fed, ghe (pathological triangle)
        vertices = [[2,1,4], [0,4,1], [3,4,0], [5,4,3],
     	            [6,7,4]]	
		    
        mesh = Mesh(points, vertices)
        mesh.check_integrity()

        P = mesh.get_boundary_polygon(verbose=False)

        
        assert len(P) == 9
	
	# Note that point e appears twice!
        assert num.allclose(P, [a, d, f, e, g, h, e, c, b])

        for p in points:
	    msg = 'Point %s is not inside polygon %s'\
	    %(p, P)	
            assert is_inside_polygon(p, P), msg		    


				 	
	    
	    

    def test_boundary_polygon_IV(self):
        """Reproduce test test_spatio_temporal_file_function_time
        from test_util.py that looked as if it produced the wrong boundary
        """

        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular        

        #Create a domain to hold test grid
        #(0:15, -20:10)
        points, vertices, boundary =\
                rectangular(4, 4, 15, 30, origin = (0, -20))        

        #####
        mesh = Mesh(points, vertices)
        mesh.check_integrity()

        P = mesh.get_boundary_polygon()

        #print P
        assert len(P) == 16
        for p in points:
            assert is_inside_polygon(p, P)



        #####
        mesh = Mesh(points, vertices, boundary)
        mesh.check_integrity()

        P = mesh.get_boundary_polygon()

        
        #print P, len(P)
        assert len(P) == 16

        for p in points:
            assert is_inside_polygon(p, P)

        #print mesh.statistics()    



    def test_boundary_polygon_V(self):
        """Create a discontinuous mesh (duplicate vertices)
        and check that boundary is as expected
        
        """

        #Points
        a = [0.0, 0.0] #0
        b = [0.0, 0.5] #1
        c = [0.0, 1.0] #2
        d = [0.5, 0.0] #3
        e = [0.5, 0.5] #4
        f = [1.0, 0.0] #5
        g = [1.0, 0.5] #6
        h = [1.0, 1.0] #7
        i = [1.5, 0.5] #8

        #Duplicate points for triangles edg [4,3,6] (central) and
        #gid [6,8,7] (top right boundary) to them disconnected
        #from the others

        e0 = [0.5, 0.5] #9 
        d0 = [0.5, 0.0] #10
        g0 = [1.0, 0.5] #11
        i0 = [1.5, 0.5] #12        
        

        points = [a, b, c, d, e, f, g, h, i, e0, d0, g0, i0]



        #dea, bae, bec, fgd,
        #edg, ghe, gfi, gih
        #vertices = [ [3,4,0], [1,0,4], [1,4,2], [5,6,3],
        #             [4,3,6], [6,7,4], [6,5,8], [6,8,7]]


        #dea, bae, bec, fgd,
        #e0d0g0, ghe, gfi, g0i0h
        vertices = [[3,4,0], [1,0,4], [1,4,2], [5,6,3],
                    [9,10,11], [6,7,4], [6,5,8], [11,12,7]]        

        mesh = Mesh(points, vertices)

        mesh.check_integrity()

        P = mesh.get_boundary_polygon()

        #print P
        
        assert len(P) == 8
        assert num.allclose(P, [a, d, f, i, h, e, c, b])
        assert num.allclose(P, [(0.0, 0.0), (0.5, 0.0), (1.0, 0.0), (1.5, 0.5), (1.0, 1.0), (0.5, 0.5), (0.0, 1.0), (0.0, 0.5)])
        

        for p in points:
            #print p, P
            assert is_inside_polygon(p, P)



    def test_boundary_polygon_VI(self):
        """test_boundary_polygon_VI(self)

        Create a discontinuous mesh (duplicate vertices) from a real situation that failed
        and check that boundary is as expected
        """

        # First do the continuous version of mesh

        points = [[   6626.85400391,      0.        ],
                  [      0.        ,  38246.4140625 ],
                  [   9656.2734375 ,  68351.265625  ],
                  [  20827.25585938,  77818.203125  ],
                  [  32755.59375   ,  58126.9765625 ],
                  [  35406.3359375 ,  79332.9140625 ],
                  [  31998.23828125,  88799.84375   ],
                  [  23288.65820313, 104704.296875  ],
                  [  32187.57617188, 109816.4375    ],
                  [  50364.08984375, 110763.1328125 ],
                  [  80468.9453125 ,  96184.0546875 ],
                  [  86149.1015625 , 129886.34375   ],
                  [ 118715.359375  , 129886.34375   ],
                  [ 117768.6640625 ,  85770.4296875 ],
                  [ 101485.5390625 ,  45251.9453125 ],
                  [  49985.4140625 ,   2272.06396484],
                  [  51737.94140625,  90559.2109375 ],
                  [  56659.0703125 ,  65907.6796875 ],
                  [  75735.4765625 ,  23762.00585938],
                  [  52341.70703125,  38563.39453125]]

        triangles = [[19, 0,15],
                     [ 2, 4, 3],
                     [ 4, 2, 1],
                     [ 1,19, 4],
                     [15,18,19],
                     [18,14,17],
                     [19, 1, 0],
                     [ 6, 8, 7],
                     [ 8, 6,16],
                     [10, 9,16],
                     [17, 5, 4],
                     [16,17,10],
                     [17,19,18],
                     [ 5,17,16],
                     [10,14,13],
                     [10,17,14],
                     [ 8,16, 9],
                     [12,11,10],
                     [10,13,12],
                     [19,17, 4],
                     [16, 6, 5]]


        triangles = num.array(triangles,num.int)
        points = num.array(points,num.float)

        mesh = Mesh(points, triangles)
        mesh.check_integrity()
        Pref = mesh.get_boundary_polygon()

        #plot_polygons([ensure_numeric(Pref)], 'goodP')

        for p in points:
            assert is_inside_polygon(p, Pref)
        
        
        # Then do the discontinuous version
        import warnings
        warnings.filterwarnings('ignore')

        
        points = [[  52341.70703125,  38563.39453125],
                  [   6626.85400391,      0.        ],
                  [  49985.4140625 ,   2272.06396484],
                  [   9656.2734375 ,  68351.265625  ],
                  [  32755.59375   ,  58126.9765625 ],
                  [  20827.25585938,  77818.203125  ],
                  [  32755.59375   ,  58126.9765625 ],
                  [   9656.2734375 ,  68351.265625  ],
                  [      0.        ,  38246.4140625 ],
                  [      0.        ,  38246.4140625 ],
                  [  52341.70703125,  38563.39453125],
                  [  32755.59375   ,  58126.9765625 ],
                  [  49985.4140625 ,   2272.06396484],
                  [  75735.4765625 ,  23762.00585938],
                  [  52341.70703125,  38563.39453125],
                  [  75735.4765625 ,  23762.00585938],
                  [ 101485.5390625 ,  45251.9453125 ],
                  [  56659.0703125 ,  65907.6796875 ],
                  [  52341.70703125,  38563.39453125],
                  [      0.        ,  38246.4140625 ],
                  [   6626.85400391,      0.        ],
                  [  31998.23828125,  88799.84375   ],
                  [  32187.57617188, 109816.4375    ],
                  [  23288.65820313, 104704.296875  ],
                  [  32187.57617188, 109816.4375    ],
                  [  31998.23828125,  88799.84375   ],
                  [  51737.94140625,  90559.2109375 ],
                  [  80468.9453125 ,  96184.0546875 ],
                  [  50364.08984375, 110763.1328125 ],
                  [  51737.94140625,  90559.2109375 ],
                  [  56659.0703125 ,  65907.6796875 ],
                  [  35406.3359375 ,  79332.9140625 ],
                  [  32755.59375   ,  58126.9765625 ],
                  [  51737.94140625,  90559.2109375 ],
                  [  56659.0703125 ,  65907.6796875 ],
                  [  80468.9453125 ,  96184.0546875 ],
                  [  56659.0703125 ,  65907.6796875 ],
                  [  52341.70703125,  38563.39453125],
                  [  75735.4765625 ,  23762.00585938],
                  [  35406.3359375 ,  79332.9140625 ],
                  [  56659.0703125 ,  65907.6796875 ],
                  [  51737.94140625,  90559.2109375 ],
                  [  80468.9453125 ,  96184.0546875 ],
                  [ 101485.5390625 ,  45251.9453125 ],
                  [ 117768.6640625 ,  85770.4296875 ],
                  [  80468.9453125 ,  96184.0546875 ],
                  [  56659.0703125 ,  65907.6796875 ],
                  [ 101485.5390625 ,  45251.9453125 ],
                  [  32187.57617188, 109816.4375    ],
                  [  51737.94140625,  90559.2109375 ],
                  [  50364.08984375, 110763.1328125 ],
                  [ 118715.359375  , 129886.34375   ],
                  [  86149.1015625 , 129886.34375   ],
                  [  80468.9453125 ,  96184.0546875 ],
                  [  80468.9453125 ,  96184.0546875 ],
                  [ 117768.6640625 ,  85770.4296875 ],
                  [ 118715.359375  , 129886.34375   ],
                  [  52341.70703125,  38563.39453125],
                  [  56659.0703125 ,  65907.6796875 ],
                  [  32755.59375   ,  58126.9765625 ],
                  [  51737.94140625,  90559.2109375 ],
                  [  31998.23828125,  88799.84375   ],
                  [  35406.3359375 ,  79332.9140625 ]]

        scaled_points = ensure_numeric(points, num.int)/1000  # Simplify for ease of interpretation

        triangles = [[ 0, 1, 2],
                     [ 3, 4, 5],
                     [ 6, 7, 8],
                     [ 9,10,11],
                     [12,13,14],
                     [15,16,17],
                     [18,19,20],
                     [21,22,23],
                     [24,25,26],
                     [27,28,29],
                     [30,31,32],
                     [33,34,35],
                     [36,37,38],
                     [39,40,41],
                     [42,43,44],
                     [45,46,47],
                     [48,49,50],
                     [51,52,53],
                     [54,55,56],
                     [57,58,59],
                     [60,61,62]]


        # First use scaled points for ease of debugging
        mesh = Mesh(scaled_points, triangles)
        mesh.check_integrity()
        P = mesh.get_boundary_polygon()

        for p in scaled_points:
            assert is_inside_polygon(p, P)            

        # Then use original points and test        
        mesh = Mesh(points, triangles)
        mesh.check_integrity()
        P = mesh.get_boundary_polygon()

        for p in points:
            assert is_inside_polygon(p, P)            

        assert num.allclose(P, Pref)    

    def test_lone_vertices(self):
        a = [2.0, 1.0]
        b = [6.0, 2.0]
        c = [1.0, 3.0]
        d = [2.0, 4.0]
        e = [4.0, 3.0]

        points = [a, b, d, c, e]
        vertices = [[0,1,3]]

        mesh = Mesh(points, vertices)
        mesh.check_integrity()
        loners = mesh.get_lone_vertices()

        #print loners
        self.assertTrue(loners==[2,4],
                        'FAILED!')

        
        a = [2.0, 1.0]
        b = [6.0, 2.0]
        c = [1.0, 3.0]
        d = [2.0, 4.0]

        points = [d, a, b, c]
        vertices = [[3,1,2]]

        mesh = Mesh(points, vertices)
        mesh.check_integrity()
        loners = mesh.get_lone_vertices()
        self.assertTrue(loners==[0],
                        'FAILED!') 

    def test_mesh_get_boundary_polygon_with_georeferencing(self):
        """test_mesh_get_boundary_polygon_with_georeferencing
        
        Test that get_boundary_polygon returns absolute coordinates
        """
        
        # test
        a = [0.0, 0.0]
        b = [4.0, 0.0]
        c = [0.0, 4.0]

        absolute_points = [a, b, c]
        vertices = [[0, 1, 2]]
        
        geo = Geo_reference(56, 67, -56)

        relative_points = geo.change_points_geo_ref(absolute_points)

        mesh = Mesh(relative_points, vertices, geo_reference=geo)
        boundary_polygon = mesh.get_boundary_polygon()

        assert num.allclose(absolute_points, boundary_polygon)

    def test_get_triangle_containing_point(self):

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]
        mesh = Mesh(points, vertices)
        
        mesh.check_integrity()


        try:
            id = mesh.get_triangle_containing_point([3.0, 5.0])
        except:
            pass
        else:
            msg = 'Should have caught point outside polygon (Non)'            
            raise Exception(msg)
            
        id = mesh.get_triangle_containing_point([0.5, 1.0])
        assert id == 0

        id = mesh.get_triangle_containing_point([1.0, 3.0])
        assert id == 3        

        for i, point in enumerate(mesh.get_centroid_coordinates()):
            id = mesh.get_triangle_containing_point(point)
            assert id == i        

    def test_get_triangle_neighbours(self):
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        e = [2.0, 2.0]
        points = [a, b, c, e]
        vertices = [ [1,0,2], [1,2,3] ]   #bac, bce
        mesh = Mesh(points, vertices)

        neighbours = mesh.get_triangle_neighbours(0)
        assert num.allclose(neighbours, [-1,1,-2])
        neighbours = mesh.get_triangle_neighbours(-10)
        assert neighbours == []
        neighbours = mesh.get_triangle_neighbours(2)
        assert neighbours == []


    def test_get_intersecting_segments1(self):
        """test_get_intersecting_segments(self):

        Very simple test (horizontal lines)
        """

        # Build test mesh
        
        # Create basic mesh
        # 9 points at (0,0), (0, 1), ..., (2,2)
        # 8 triangles enumerated from left bottom to right top.
        points, vertices, boundary = rectangular(2, 2, 2, 2)
        mesh = Mesh(points, vertices, boundary)

        # Very simple horizontal line intersecting
        #


        for y_line in [0.1, 0.2, 0.314159, 0.41, 0.6, 0.99, 1.01, 1.5, 1.77, 1.9]:
            if y_line < 1:
                ceiling = 1
                floor = 0
                intersected_triangles = [0,1,4,5]
            elif y_line > 1:
                ceiling = 2
                floor = 1
                intersected_triangles = [2,3,6,7]
            else:
                raise Exception('this test is not for parallel lines')


            line = [[-1,y_line], [3,y_line]]

            L = mesh.get_intersecting_segments(line)
            assert len(L) == 4

            

            # Check all normals point straight down etc
            total_length = 0
            for x in L:
                if x.triangle_id % 2 == 0:
                    assert num.allclose(x.length, ceiling-y_line)
                else:
                    assert num.allclose(x.length, y_line-floor)                

                
                assert num.allclose(x.normal, [0,-1])

                assert num.allclose(x.segment[1][0], x.segment[0][0] + x.length)
                assert num.allclose(x.segment[0][1], y_line)
                assert num.allclose(x.segment[1][1], y_line)                

                assert x.triangle_id in intersected_triangles

                total_length += x.length

            msg = 'Segments do not add up'
            assert num.allclose(total_length, 2), msg
            

    def test_get_intersecting_segments_coinciding(self):
        """test_get_intersecting_segments_coinciding(self):

        Test that lines coinciding with triangle edges work.
        """

        # Build test mesh
        
        # Create basic mesh
        # 9 points at (0,0), (0, 1), ..., (2,2)
        # 8 triangles enumerated from left bottom to right top.
        points, vertices, boundary = rectangular(2, 2, 2, 2)
        mesh = Mesh(points, vertices, boundary)
        intersected_triangles = [1,5]

        # Very simple horizontal line intersecting
        #

        y_line = 1.0
        
        line = [[-1,y_line], [3,y_line]]

        L = mesh.get_intersecting_segments(line)


        msg = 'Only two triangles should be returned'    
        assert len(L) == 2, msg    
            

        # Check all 
        total_length = 0
        for x in L:
            assert num.allclose(x.length, 1.0)
            assert num.allclose(x.normal, [0,-1])

            assert num.allclose(x.segment[1][0], x.segment[0][0] + x.length)
            assert num.allclose(x.segment[0][1], y_line)
            assert num.allclose(x.segment[1][1], y_line)                            



            assert x.triangle_id in intersected_triangles

            total_length += x.length

        msg = 'Segments do not add up'
        assert num.allclose(total_length, 2), msg
        

    def test_get_intersecting_segments_partially_coinciding(self):
        """test_get_intersecting_segments_partially_coinciding(self):

        Test that line coinciding with triangle edges work.
        But this ones only coincide with parts of the edge. 
        """

        # Build test mesh
        
        # Create basic mesh
        # 9 points at (0,0), (0, 1), ..., (2,2)
        # 8 triangles enumerated from left bottom to right top.
        points, vertices, boundary = rectangular(2, 2, 2, 2)
        mesh = Mesh(points, vertices, boundary)
        intersected_triangles = [1,5]

        # Horizontal line intersecting along center but stopping
        # parway through second triangle's edge
        #

        y_line = 1.0
        
        #line = [[0, y_line], [2, y_line]]
        line = [[0, y_line], [1.5, y_line]]

        L = mesh.get_intersecting_segments(line)
        #for x in L:
        #    print x

        msg = 'Two triangles should be returned'    
        assert len(L) == 2, msg    
            

        # Check all 
        total_length = 0
        for x in L:
            if x.triangle_id == 1:
                assert num.allclose(x.length, 1)        
                assert num.allclose(x.normal, [0, -1])
                
            if x.triangle_id == 5:
                assert num.allclose(x.length, 0.5)
                assert num.allclose(x.normal, [0, -1])


            assert x.triangle_id in intersected_triangles

            total_length += x.length

        msg = 'Segments do not add up'
        assert num.allclose(total_length, 1.5), msg            



    def test_get_intersecting_segments2(self):
        """test_get_intersecting_segments(self):

        Lines with a slope
        """

        s2 = sqrt(2.0)/2
        

        # Build test mesh
        
        # Create basic mesh
        # 9 points at (0,0), (0, 1), ..., (2,2)
        # 8 triangles enumerated from left bottom to right top.
        points, vertices, boundary = rectangular(2, 2, 2, 2)
        mesh = Mesh(points, vertices, boundary)
        

        # Diagonal cutting through a vertex and hypothenuses
        line = [[0, 2], [2, 0]]
        intersected_triangles = [3,2,5,4]                

        L = mesh.get_intersecting_segments(line)
        assert len(L) == 4

        #print L
        
        # Check all segments
        total_length = 0
        for i, x in enumerate(L): 
            assert num.allclose(x.length, s2)
            assert num.allclose(x.normal, [-s2, -s2])
            assert num.allclose(sum(x.normal**2), 1)
            
            assert x.triangle_id in intersected_triangles

            total_length += x.length

        msg = 'Segments do not add up'
        assert num.allclose(total_length, 4*s2), msg


        # Diagonal cutting through a vertex and hypothenuses (reversed)
        line = [[2, 0], [0, 2]]
        intersected_triangles = [3,2,5,4]                

        L = mesh.get_intersecting_segments(line)
        assert len(L) == 4

        #print L
        
        # Check all segments
        total_length = 0
        for i, x in enumerate(L): 
            assert num.allclose(x.length, s2)
            assert num.allclose(x.normal, [s2, s2])
            assert num.allclose(sum(x.normal**2), 1)
            
            assert x.triangle_id in intersected_triangles

            total_length += x.length

        msg = 'Segments do not add up'
        assert num.allclose(total_length, 4*s2), msg                    



        # Diagonal coinciding with hypothenuses
        line = [[2, 2], [0, 0]]
        intersected_triangles = [6,0]                

        L = mesh.get_intersecting_segments(line)
        assert len(L) == 2

        #print L
        
        # Check all segments
        total_length = 0
        for i, x in enumerate(L): 
            assert num.allclose(x.length, 2*s2)
            assert num.allclose(x.normal, [-s2, s2])
            assert num.allclose(sum(x.normal**2), 1)
            
            assert x.triangle_id in intersected_triangles

            total_length += x.length

        msg = 'Segments do not add up'
        assert num.allclose(total_length, 4*s2), msg                        


        # Diagonal coinciding with hypothenuses (reversed)
        line = [[0, 0], [2, 2]]
        intersected_triangles = [6,0]                

        L = mesh.get_intersecting_segments(line)
        assert len(L) == 2

        #print L
        
        # Check all segments
        total_length = 0
        for i, x in enumerate(L): 
            assert num.allclose(x.length, 2*s2)
            assert num.allclose(x.normal, [s2, -s2])
            assert num.allclose(sum(x.normal**2), 1)
            
            assert x.triangle_id in intersected_triangles

            total_length += x.length

        msg = 'Segments do not add up'
        assert num.allclose(total_length, 4*s2), msg                        



        # line with slope [1, -1] cutting through vertices of tri 7 and 6
        line = [[1, 2], [2, 1]]
        intersected_triangles = [7,6]                

        L = mesh.get_intersecting_segments(line)
        assert len(L) == 2

        #print L
        
        # Check all segments
        total_length = 0
        for i, x in enumerate(L): 
            assert num.allclose(x.length, s2)
            assert num.allclose(x.normal, [-s2, -s2])
            assert num.allclose(sum(x.normal**2), 1)
            
            assert x.triangle_id in intersected_triangles

            total_length += x.length

        msg = 'Segments do not add up'
        assert num.allclose(total_length, 2*s2), msg


        # Arbitrary line with slope [1, -1] cutting through tri 7 and 6
        line = [[1.1, 2], [2.1, 1]]
        intersected_triangles = [7,6]                

        L = mesh.get_intersecting_segments(line)
        assert len(L) == 2
        
        # Check all segments
        total_length = 0
        for i, x in enumerate(L): 
            assert num.allclose(x.normal, [-s2, -s2])
            assert num.allclose(sum(x.normal**2), 1)

            msg = 'Triangle %d' %x.triangle_id + ' is not in %s' %(intersected_triangles)
            assert x.triangle_id in intersected_triangles, msg
            


    def test_get_intersecting_segments3(self):
        """test_get_intersecting_segments(self):

        Check that line can stop inside a triangle
        
        """



        s2 = sqrt(2.0)/2
        

        # Build test mesh
        
        # Create basic mesh
        # 9 points at (0,0), (0, 1), ..., (2,2)
        # 8 triangles enumerated from left bottom to right top.
        points, vertices, boundary = rectangular(2, 2, 2, 2)
        mesh = Mesh(points, vertices, boundary)
        

        # Line cutting through one triangle and ending on its edge
        line = [[0.5, 3], [0.5, 1.5]]
        intersected_triangles = [3]                

        L = mesh.get_intersecting_segments(line)
        assert len(L) == 1
        assert L[0].triangle_id == 3
        assert num.allclose(L[0].length, 0.5)        
        assert num.allclose(L[0].normal, [-1,0])                



        # Now try to shorten it so that its endpoint falls short of the far edge
        line = [[0.5, 3], [0.5, 1.6]]
        intersected_triangles = [3]                

        L = mesh.get_intersecting_segments(line)
        assert len(L) == 1
        assert L[0].triangle_id == 3
        assert num.allclose(L[0].length, 0.4)
        assert num.allclose(L[0].normal, [-1,0])

        intersected_triangles = [3]

        # Now the same, but with direction changed
        line = [[0.5, 3], [0.5, 1.6]]
        line = [[0.5, 1.6], [0.5, 3]]        
        intersected_triangles = [3]                

        L = mesh.get_intersecting_segments(line)
        assert len(L) == 1
        assert L[0].triangle_id == 3
        assert num.allclose(L[0].length, 0.4)
        assert num.allclose(L[0].normal, [1,0])                
        

            

    def test_get_intersecting_segments4(self):
        """test_get_intersecting_segments(self):

        Simple poly line
        
        """



        s2 = sqrt(2.0)/2
        

        # Build test mesh
        
        # Create basic mesh
        # 9 points at (0,0), (0, 1), ..., (2,2)
        # 8 triangles enumerated from left bottom to right top.
        points, vertices, boundary = rectangular(2, 2, 2, 2)
        mesh = Mesh(points, vertices, boundary)
        

        # Polyline with three segments cutting through mesh
        line = [[0.5, 3], [0.5, 1.5], [1,1]]

        L = mesh.get_intersecting_segments(line)
        assert len(L) == 2

        for x in L:
            if x.triangle_id == 3:
                assert num.allclose(x.length, 0.5)        
                assert num.allclose(x.normal, [-1,0])
                
            if x.triangle_id == 2:
                assert num.allclose(x.length, s2)
                assert num.allclose(x.normal, [-s2,-s2])



    def test_get_intersecting_segments5(self):
        """test_get_intersecting_segments(self):

        More complex poly line
        
        """



        s2 = sqrt(2.0)/2
        

        # Build test mesh
        
        # Create basic mesh
        # 9 points at (0,0), (0, 1), ..., (2,2)
        # 8 triangles enumerated from left bottom to right top.
        points, vertices, boundary = rectangular(2, 2, 2, 2)
        mesh = Mesh(points, vertices, boundary)
        

        # Polyline with three segments cutting through mesh
        line = [[0.5, 3], [0.5, 1.5], [1.25, 0.75]] 

        L = mesh.get_intersecting_segments(line)
        assert len(L) == 3

        for x in L:
            if x.triangle_id == 3:
                assert num.allclose(x.length, 0.5)        
                assert num.allclose(x.normal, [-1,0])
                
            if x.triangle_id == 2:
                msg = str(x.length)
                assert num.allclose(x.length, s2), msg
                assert num.allclose(x.normal, [-s2,-s2])

            if x.triangle_id == 5:
                segvec = num.array([line[2][0]-1,
                                    line[2][1]-1])
                msg = str(x.length)
                assert num.allclose(x.length, sqrt(sum(segvec**2))), msg
                assert num.allclose(x.normal, [-s2,-s2])                                                


    def test_get_intersecting_segments6(self):
        """test_get_intersecting_segments(self):

        Even more complex poly line, where line breaks within triangle 5

        5 segments are returned even though only four triangles [3,2,5,6] are touched.
        Triangle 5 therefor has two segments in it.
        
        """



        s2 = sqrt(2.0)/2
        

        # Build test mesh
        
        # Create basic mesh
        # 9 points at (0,0), (0, 1), ..., (2,2)
        # 8 triangles enumerated from left bottom to right top.
        points, vertices, boundary = rectangular(2, 2, 2, 2)
        mesh = Mesh(points, vertices, boundary)
        

        # Polyline with three segments cutting through mesh
        line = [[0.5, 3], [0.5, 1.5], [1.25, 0.75], [2.25, 1.75]]

        L = mesh.get_intersecting_segments(line)
        #for x in L:
        #    print x

        assert len(L) == 5

        for x in L:
            if x.triangle_id == 3:
                assert num.allclose(x.length, 0.5)        
                assert num.allclose(x.normal, [-1,0])
                
            if x.triangle_id == 2:
                msg = str(x.length)
                assert num.allclose(x.length, s2), msg
                assert num.allclose(x.normal, [-s2,-s2])

            if x.triangle_id == 5:
                if x.segment == ((1.0, 1.0), (1.25, 0.75)):                    
                    segvec = num.array([line[2][0]-1,
                                        line[2][1]-1])
                    msg = str(x.length)
                    assert num.allclose(x.length, sqrt(sum(segvec**2))), msg
                    assert num.allclose(x.normal, [-s2,-s2])
                elif x.segment == ((1.25, 0.75), (1.5, 1.0)):
                    segvec = num.array([1.5-line[2][0],
                                        1.0-line[2][1]])
                    
                    assert num.allclose(x.length, sqrt(sum(segvec**2))), msg
                    assert num.allclose(x.normal, [s2,-s2])
                else:
                    msg = 'Unknown segment: %s' %x.segment
                    raise Exception(msg)
                

                    
            if x.triangle_id == 6:
                assert num.allclose(x.normal, [s2,-s2])
                assert num.allclose(x.segment, ((1.5, 1.0), (2, 1.5)))


      # Internal test that sum of line segments add up
      # to length of input line
      #
      # Could be useful perhaps
      #
      #xi1 = line[1][0]
      #eta1 = line[1][1]
      #linevector = num.array([xi1-xi0, eta1-eta0])
      #linelength = sqrt(sum(linevector**2))
      #
      #segmentlength = 0
      #for segment in triangle_intersections:
      #    vector = array([segment[1][0] - segment[0][0],
      #                    segment[1][1] - segment[0][1]])
      #    length = sqrt(sum(vector**2))      
      #    segmentlength += length
      #
      #msg = 'Sum of intersecting segments do not add up'    
      #assert allclose(segmentlength, linelength), msg    




    def test_get_intersecting_segments7(self):
        """test_get_intersecting_segments(self):

        Check that line can stop inside a triangle - this is from
        flow throug a cross sections example in test_datamanager.
        
        """

        # Build test mesh
        width = 5
        length = 100
        t_end = 1
        points, vertices, boundary = rectangular(length, width,
                                                 length, width)

        mesh = Mesh(points, vertices, boundary)
        

        # A range of partial lines
        x = length/2.
        for i in range(10):
            start_point = [length/2., i*width/10.]
            #print 
            #print start_point
                            
            line = [start_point, [length/2., width]]
 
            L = mesh.get_intersecting_segments(line)

            if start_point[1] < 1:
                assert len(L) == 5
                
            
            total_length = 0    
            for x in L:
                total_length += x.length
                

            ref_length = line[1][1] - line[0][1]
            #print ref_length, total_length
            assert num.allclose(total_length, ref_length)


#-------------------------------------------------------------

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_Mesh, 'test_mesh_and_neighbours')
    runner = unittest.TextTestRunner()#verbosity=2)
    runner.run(suite)

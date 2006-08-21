#!/usr/bin/env python



#FIXME: Seperate the tests for mesh and general_mesh

#FIXME (Ole): Maxe this test independent of anything that inherits from General_mesh (namely shallow_water)

import unittest
from math import sqrt

from neighbour_mesh import *
from mesh_factory import rectangular
from anuga.config import epsilon
from Numeric import allclose, array

from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.utilities.polygon import is_inside_polygon

def distance(x, y):
    return sqrt( sum( (array(x)-array(y))**2 ))

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
            raise msg


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
        assert allclose(normals[0, 0:2], [3.0/5, 4.0/5])
        assert allclose(normals[0, 2:4], [-1.0, 0.0])
        assert allclose(normals[0, 4:6], [0.0, -1.0])

        assert allclose(mesh.get_normal(0,0), [3.0/5, 4.0/5])
        assert allclose(mesh.get_normal(0,1), [-1.0, 0.0])
        assert allclose(mesh.get_normal(0,2), [0.0, -1.0])

        #Edge lengths
        assert allclose(mesh.edgelengths[0], [5.0, 3.0, 4.0])


        #Vertex coordinates
        V = mesh.get_vertex_coordinates()
        assert allclose(V[0], [0.0, 0.0, 4.0, 0.0, 0.0, 3.0])

        V = mesh.get_vertex_coordinates(obj=True)
        assert allclose(V, [ [0.0, 0.0],
                             [4.0, 0.0],
                             [0.0, 3.0] ])

        V0 = mesh.get_vertex_coordinate(0, 0)
        assert allclose(V0, [0.0, 0.0])

        V1 = mesh.get_vertex_coordinate(0, 1)
        assert allclose(V1, [4.0, 0.0])

        V2 = mesh.get_vertex_coordinate(0, 2)
        assert allclose(V2, [0.0, 3.0])


        #General tests:

        #Test that points are arranged in a counter clock wise order etc
        mesh.check_integrity()


        #Test that the centroid is located 2/3 of the way
        #from each vertex to the midpoint of the opposite side

        V = mesh.get_vertex_coordinates()

        x0 = V[0,0]
        y0 = V[0,1]
        x1 = V[0,2]
        y1 = V[0,3]
        x2 = V[0,4]
        y2 = V[0,5]

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

        x0 = V[0,0]
        y0 = V[0,1]
        x1 = V[0,2]
        y1 = V[0,3]
        x2 = V[0,4]
        y2 = V[0,5]

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
        assert allclose(mesh.radii[0],sqrt(3.0)/3),'Steve''s doesn''t work'

        mesh = Mesh(points, vertices,use_inscribed_circle=True)
        assert allclose(mesh.radii[0],sqrt(3.0)/3),'inscribed circle doesn''t work'

    def test_inscribed_circle_rightangle_triangle(self):
        """test that the radius is calculated correctly by mesh in the case of a right-angled triangle"""
        a = [0.0, 0.0]
        b = [4.0, 0.0]
        c = [0.0, 3.0]

        points = [a, b, c]
        vertices = [[0,1,2]]

        mesh = Mesh(points, vertices,use_inscribed_circle=False)
        assert allclose(mesh.radii[0],5.0/6),'Steve''s doesn''t work'

        mesh = Mesh(points, vertices,use_inscribed_circle=True)
        assert allclose(mesh.radii[0],1.0),'inscribed circle doesn''t work'


    def test_two_triangles(self):
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        e = [2.0, 2.0]
        points = [a, b, c, e]
        vertices = [ [1,0,2], [1,2,3] ]   #bac, bce
        mesh = Mesh(points, vertices)

        assert mesh.areas[0] == 2.0

        assert allclose(mesh.centroid_coordinates[0], [2.0/3, 2.0/3])


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

        assert allclose(mesh.centroid_coordinates[0], [2.0/3, 2.0/3])
        assert allclose(mesh.centroid_coordinates[1], [4.0/3, 4.0/3])
        assert allclose(mesh.centroid_coordinates[2], [8.0/3, 2.0/3])
        assert allclose(mesh.centroid_coordinates[3], [2.0/3, 8.0/3])

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
            raise "triangle edge duplicates not caught"

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
            raise 'Should have raised an exception'

        #Too few points - 1 element
        try:
            mesh = Mesh([points[0]], [vertices[0]])
        except AssertionError:
            pass
        else:
            raise 'Should have raised an exception'

        #Wrong dimension of vertices
        try:
            mesh = Mesh(points, vertices[0])
        except AssertionError:
            pass
        else:
            raise 'Should have raised an exception'

        #Unsubscriptable coordinates object raises exception
        try:
            mesh = Mesh(points[0], [vertices[0]])
        except AssertionError:
            pass
        else:
            raise 'Should have raised an exception'

        #FIXME: This has been commented out pending a decision
        #whether to allow partial boundary tags or not
        #
        #Not specifying all boundary tags
        #try:
        #    mesh = Mesh(points, vertices, {(3,0): 'x'})
        #except AssertionError:
        #    pass
        #else:
        #    raise 'Should have raised an exception'

        #Specifying wrong non existing segment
        try:
            mesh = Mesh(points, vertices, {(5,0): 'x'})
        except AssertionError:
            pass
        else:
            raise 'Should have raised an exception'




    def test_internal_boundaries(self):
        """
        get values based on triangle lists.
        """
        from mesh_factory import rectangular
        from shallow_water import Domain
        from Numeric import zeros, Float

        #Create basic mesh
        points, vertices, boundary = rectangular(1, 3)

        # Add an internal boundary
        boundary[(2,0)] = 'internal'
        boundary[(1,0)] = 'internal'

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.build_tagged_elements_dictionary({'bottom':[0,1],
                                                 'top':[4,5],
                                                 'all':[0,1,2,3,4,5]})


    def test_boundary_polygon(self):
        from mesh_factory import rectangular
        #from mesh import Mesh
        from Numeric import zeros, Float

        #Create basic mesh
        points, vertices, boundary = rectangular(2, 2)
        mesh = Mesh(points, vertices, boundary)


        P = mesh.get_boundary_polygon()

        assert len(P) == 8
        assert allclose(P, [[0.0, 0.0], [0.5, 0.0], [1.0, 0.0],
                            [1.0, 0.5], [1.0, 1.0], [0.5, 1.0],
                            [0.0, 1.0], [0.0, 0.5]])
        for p in points:
            #print p, P
            assert is_inside_polygon(p, P)


    def test_boundary_polygon_II(self):
        from Numeric import zeros, Float
        

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
        assert allclose(P, [a, d, f, i, h, e, c, b])

        for p in points:
            #print p, P
            assert is_inside_polygon(p, P)


    def test_boundary_polygon_III(self):
        """Same as II but vertices ordered differently
        """

        from Numeric import zeros, Float


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
        assert allclose(P, [a, d, f, i, h, e, c, b])

        for p in points:
            assert is_inside_polygon(p, P)


    def test_boundary_polygon_IV(self):
        """Reproduce test test_spatio_temporal_file_function_time
        from test_util.py that looked as if it produced the wrong boundary
        """

        from Numeric import zeros, Float
        from mesh_factory import rectangular        

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
        from Numeric import zeros, Float
        

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
        assert allclose(P, [a, d, f, i, h, e, c, b])
        assert allclose(P, [(0.0, 0.0), (0.5, 0.0), (1.0, 0.0), (1.5, 0.5), (1.0, 1.0), (0.5, 0.5), (0.0, 1.0), (0.0, 0.5)])
        

        for p in points:
            #print p, P
            assert is_inside_polygon(p, P)



    def test_lone_vertices(self):
        a = [2.0, 1.0]
        b = [6.0, 2.0]
        c = [1.0, 3.0]
        d = [2.0, 4.0]

        points = [a, b, c, d]
        vertices = [[0,1,2]]

        mesh = Mesh(points, vertices)
        mesh.check_integrity()
        loners = mesh.get_lone_vertices()
        self.failUnless(loners==[3],
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
        self.failUnless(loners==[0],
                        'FAILED!') 

    def test_mesh_get_boundary_polygon_with_georeferencing(self):
     
        # test
        a = [0.0, 0.0]
        b = [4.0, 0.0]
        c = [0.0, 4.0]

        absolute_points = [a, b, c]
        vertices = [[0,1,2]]
        
        geo = Geo_reference(56,67,-56)

        relative_points = geo.change_points_geo_ref(absolute_points)

        #print 'Relative', relative_points
        #print 'Absolute', absolute_points       

        mesh = Mesh(relative_points, vertices, geo_reference=geo)
        boundary_polygon = mesh.get_boundary_polygon()

        assert allclose(absolute_points, boundary_polygon)
        
#-------------------------------------------------------------
if __name__ == "__main__":
    #suite = unittest.makeSuite(Test_Mesh,'test_mesh_get_boundary_polygon_with_georeferencing')
    suite = unittest.makeSuite(Test_Mesh,'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)





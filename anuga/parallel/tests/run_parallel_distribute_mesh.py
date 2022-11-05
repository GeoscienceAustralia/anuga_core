#!/usr/bin/env python


import unittest
import sys
import os
from math import sqrt



import anuga

from anuga import Domain
from anuga import rectangular_cross

from anuga.parallel.distribute_mesh import pmesh_divide_metis
from anuga.parallel.distribute_mesh import build_submesh
from anuga.parallel.distribute_mesh import submesh_full, submesh_ghost, submesh_quantities
from anuga.parallel.distribute_mesh import extract_submesh, rec_submesh, send_submesh

from anuga import myid, numprocs, barrier, finalize

from pprint import pprint, pformat

import numpy
from numpy import array

#----------------------------------
# set up MPI to abort on error
#----------------------------------
from anuga.utilities.parallel_abstraction import global_except_hook
import sys
sys.excepthook = global_except_hook


verbose = False
verbose0 = False
verbose1 = False
verbose2 = False


def topography(x, y):
    return -x/2.0


def xcoord(x, y):
    return x


def ycoord(x, y):
    return y


def distibute_three_processors():
    """
    Do a parallel test of distributing a rectangle onto 3 processors

    """


    # FIXME: Need to update expected values on macos
    if sys.platform == 'darwin':
        return

    # FIXME: Need to update expected values on macos
    #if sys.platform == 'win32':
    #	return


    from anuga.utilities import parallel_abstraction as pypar

    myid = anuga.myid
    numprocs = anuga.numprocs

    if not numprocs == 3:
        if verbose: print('Should be run with 3 processes')
        return

    if sys.platform == 'win32':
        from pymetis import part_graph
        metis_version = "5_part_graph"
    else:
        from pymetis import part_graph
        metis_version = "5_part_graph"
    

    if myid == 0 and verbose0: print('metis version = ', metis_version)

    #barrier()


    if myid == 0:

        nodes_0, triangles_0, boundary_0 = rectangular_cross(2, 2)

        domain = Domain(nodes_0, triangles_0, boundary_0)

        domain.set_quantity('elevation', topography)  # Use function for elevation
        domain.set_quantity('friction', 0.0)         # Constant friction
        domain.set_quantity('stage', expression='elevation')  # Dry initial stage
        domain.set_quantity('xmomentum', expression='friction + 2.0')
        domain.set_quantity('ymomentum', ycoord)

        #----------------------------------------------------------------------------------
        # Test pmesh_divide_metis
        #----------------------------------------------------------------------------------
        vertices, triangles, boundary, triangles_per_proc, quantities = pmesh_divide_metis(
            domain, numprocs)

        if verbose0:
            print_seq_values(vertices, triangles, triangles_per_proc)

        true_seq_values = get_true_seq_values(metis_version=metis_version)
        
        if verbose0:
            print("True Seq Values = \\")
            pprint(true_seq_values)

        numpy.allclose(vertices, true_seq_values['vertices'] )
        numpy.allclose(triangles, true_seq_values['triangles'] )
        numpy.allclose(triangles_per_proc, true_seq_values['triangles_per_proc'])


        #----------------------------------------------------------------------------------
        # Test build_submesh
        #----------------------------------------------------------------------------------
        submesh = build_submesh(vertices, triangles, boundary,
                                quantities, triangles_per_proc)


        if verbose0: 
            print('submesh_values = \\')
            print_submesh_values(submesh)

        true_values = get_true_submesh_values(metis_version)



        numpy.allclose(submesh['full_nodes'][0], true_values['full_nodes_0'])
        numpy.allclose(submesh['full_nodes'][1], true_values['full_nodes_1'])
        numpy.allclose(submesh['full_nodes'][2], true_values['full_nodes_2'])

        numpy.allclose(submesh['ghost_nodes'][0], true_values['ghost_nodes_0'])
        numpy.allclose(submesh['ghost_nodes'][1], true_values['ghost_nodes_1'])
        numpy.allclose(submesh['ghost_nodes'][2], true_values['ghost_nodes_2'])

        numpy.allclose(submesh['full_triangles'][0], true_values['full_triangles_0'])
        numpy.allclose(submesh['full_triangles'][1], true_values['full_triangles_1'])
        numpy.allclose(submesh['full_triangles'][2], true_values['full_triangles_2'])

        numpy.allclose(submesh['ghost_triangles'][0], true_values['ghost_triangles_0'])
        numpy.allclose(submesh['ghost_triangles'][1], true_values['ghost_triangles_1'])
        numpy.allclose(submesh['ghost_triangles'][2], true_values['ghost_triangles_2'])
        
        numpy.allclose(submesh['ghost_commun'][0], true_values['ghost_commun_0'])
        numpy.allclose(submesh['ghost_commun'][1], true_values['ghost_commun_1'])
        numpy.allclose(submesh['ghost_commun'][2], true_values['ghost_commun_2'])

        assert submesh['full_commun'] == true_values['full_commun']

    barrier()
    #--------------------------------
    # Now do the comunnication part
    #--------------------------------

    if myid == 0:

        points, vertices, boundary, quantities, \
                    ghost_recv_dict, full_send_dict, tri_map, node_map, tri_l2g, node_l2g, \
                    ghost_layer_width =\
                    extract_submesh(submesh, triangles_per_proc)

        #----------------------------------------------------------------------------------
        # Test send_submesh
        #----------------------------------------------------------------------------------
        for p in range(1, numprocs):
            send_submesh(submesh, triangles_per_proc, p, verbose=False)
    else:
        #----------------------------------------------------------------------------------
        # Test rec_submesh
        #----------------------------------------------------------------------------------
        points, triangles, boundary, quantities, \
                ghost_recv_dict, full_send_dict, \
                no_full_nodes, no_full_trigs, tri_map, node_map, tri_l2g, node_l2g, \
                ghost_layer_width = \
                rec_submesh(0, verbose=False)

    barrier()

    #--------------------------------
    # Now do the test
    #--------------------------------
    if myid == 0:

        if verbose0:
            print('extract_values = \\')
            print_extract_submesh(points, triangles, ghost_recv_dict, \
                                  full_send_dict, tri_map, node_map, ghost_layer_width)

        true_values  = get_true_extract_submesh(metis_version)



        numpy.allclose(points,   true_values['points'])
        numpy.allclose(triangles, true_values['triangles'])
        numpy.allclose(ghost_recv_dict[1], true_values['ghost_recv_dict_1'])
        numpy.allclose(ghost_recv_dict[2], true_values['ghost_recv_dict_2'])
        numpy.allclose(full_send_dict[1], true_values['full_send_dict_1'])
        numpy.allclose(full_send_dict[2], true_values['full_send_dict_2'])
        numpy.allclose(tri_map, true_values['tri_map'])
        numpy.allclose(node_map, true_values['node_map'])
        numpy.allclose(ghost_layer_width,  true_values['ghost_layer_width'])


    if myid == 1:

        if verbose1:
            print("rec_submesh_1 = \\")
            print_rec_submesh_1(points, triangles, ghost_recv_dict, full_send_dict, \
                         tri_map, node_map, ghost_layer_width)


        true_values = get_true_rec_submesh_1(metis_version)

        if verbose1:
            print('true_rec_values_1 = \\')
            pprint(true_values)

        numpy.allclose(points,   true_values['points'])
        numpy.allclose(triangles, true_values['triangles'])
        numpy.allclose(ghost_recv_dict[0], true_values['ghost_recv_dict_0'])
        numpy.allclose(ghost_recv_dict[2], true_values['ghost_recv_dict_2'])
        numpy.allclose(full_send_dict[0], true_values['full_send_dict_0'])
        numpy.allclose(full_send_dict[2], true_values['full_send_dict_2'])
        numpy.allclose(tri_map, true_values['tri_map'])
        numpy.allclose(node_map, true_values['node_map'])
        numpy.allclose(ghost_layer_width,  true_values['ghost_layer_width'])

    if myid == 2:

        if verbose2:
            print("rec_submesh_2 = \\")
            print_rec_submesh_2(points, triangles, ghost_recv_dict, full_send_dict, \
                         tri_map, node_map, ghost_layer_width)

        true_values = get_true_rec_submesh_2(metis_version)

        if verbose2:
            print('true_rec_values_2 = \\')
            pprint(true_values)

        numpy.allclose(points,   true_values['points'])
        numpy.allclose(triangles, true_values['triangles'])
        numpy.allclose(ghost_recv_dict[0], true_values['ghost_recv_dict_0'])
        numpy.allclose(ghost_recv_dict[1], true_values['ghost_recv_dict_1'])
        numpy.allclose(full_send_dict[0], true_values['full_send_dict_0'])
        numpy.allclose(full_send_dict[1], true_values['full_send_dict_1'])
        numpy.allclose(tri_map, true_values['tri_map'])
        numpy.allclose(node_map, true_values['node_map'])
        numpy.allclose(ghost_layer_width,  true_values['ghost_layer_width'])


    finalize()


#==============================================================================================


def get_true_seq_values(metis_version=4):

    from numpy import array
    import sys

    if sys.platform == 'win32' and metis_version == 4:
        true_seq_values = dict(
            triangles = array([[ 0,  9,  1],
                [ 1, 10,  2],
                [ 4, 10,  1],
                [ 2, 10,  5],
                [ 3,  9,  0],
                [ 6, 11,  3],
                [ 7, 11,  6],
                [ 4, 11,  7],
                [ 7, 12,  4],
                [ 8, 12,  7],
                [ 4,  9,  3],
                [ 1,  9,  4],
                [ 5, 10,  4],
                [ 3, 11,  4],
                [ 4, 12,  5],
                [ 5, 12,  8]]),
            triangles_per_proc = array([4, 6, 6]),
            vertices = array([[0.  , 0.  ],
                [0.  , 0.5 ],
                [0.  , 1.  ],
                [0.5 , 0.  ],
                [0.5 , 0.5 ],
                [0.5 , 1.  ],
                [1.  , 0.  ],
                [1.  , 0.5 ],
                [1.  , 1.  ],
                [0.25, 0.25],
                [0.25, 0.75],
                [0.75, 0.25],
                [0.75, 0.75]]))

        return true_seq_values

    
    if sys.platform == 'win32' and metis_version == '5_part_graph':
        true_seq_values = dict(
            triangles = array([[ 0,  9,  1],
            [ 3,  9,  0],
            [ 4,  9,  3],
            [ 1,  9,  4],
            [ 4, 10,  1],
            [ 3, 11,  4],
            [ 6, 11,  3],
            [ 7, 11,  6],
            [ 4, 11,  7],
            [ 7, 12,  4],
            [ 1, 10,  2],
            [ 5, 10,  4],
            [ 2, 10,  5],
            [ 4, 12,  5],
            [ 8, 12,  7],
            [ 5, 12,  8]]),
            triangles_per_proc = array([5, 5, 6]),
            vertices = array([[0.  , 0.  ],
            [0.  , 0.5 ],
            [0.  , 1.  ],
            [0.5 , 0.  ],
            [0.5 , 0.5 ],
            [0.5 , 1.  ],
            [1.  , 0.  ],
            [1.  , 0.5 ],
            [1.  , 1.  ],
            [0.25, 0.25],
            [0.25, 0.75],
            [0.75, 0.25],
            [0.75, 0.75]]))

        return true_seq_values


    if metis_version == "5_part_mesh":
        true_seq_values =\
            {'triangles': array([[7, 11,  6],
                                 [7, 12,  4],
                                 [8, 12,  7],
                                 [5, 12,  8],
                                 [0,  9,  1],
                                 [3,  9,  0],
                                 [4,  9,  3],
                                 [1,  9,  4],
                                 [4, 10,  1],
                                 [4, 11,  7],
                                 [4, 12,  5],
                                 [1, 10,  2],
                                 [5, 10,  4],
                                 [2, 10,  5],
                                 [3, 11,  4],
                                 [6, 11,  3]]),
                'triangles_per_proc': array([4, 7, 5]),
                'vertices': array([[0., 0.],
                                   [0., 0.5],
                                   [0., 1.],
                                   [0.5, 0.],
                                   [0.5, 0.5],
                                   [0.5, 1.],
                                   [1., 0.],
                                   [1., 0.5],
                                   [1., 1.],
                                   [0.25, 0.25],
                                   [0.25, 0.75],
                                   [0.75, 0.25],
                                   [0.75, 0.75]])}

        return true_seq_values



    vertices = array([[0.  , 0.  ],
                [0.  , 0.5 ],
                [0.  , 1.  ],
                [0.5 , 0.  ],
                [0.5 , 0.5 ],
                [0.5 , 1.  ],
                [1.  , 0.  ],
                [1.  , 0.5 ],
                [1.  , 1.  ],
                [0.25, 0.25],
                [0.25, 0.75],
                [0.75, 0.25],
                [0.75, 0.75]])


    if False:
        from pprint import pformat
        true_values = dict(
        true_ghost_layer_width = ghost_layer_width,
        true_points = points,
        true_vertices = vertices,
        true_ghost_recv_dict_1 = ghost_recv_dict[1],
        true_ghost_recv_dict_2 = ghost_recv_dict[2],
        true_full_send_dict_1 = full_send_dict[1],
        true_full_send_dict_2 = full_send_dict[2])
        for key,item in true_values.items():
            msg = key + '=' + pformat(item)
            print (msg)

    if metis_version == 4:
        triangles = [[4, 9, 3], [4, 12, 5], [7, 12, 4], [8, 12, 7], [5, 12, 8],
                    [0, 9, 1], [1, 9, 4], [1, 10, 2], [4, 10, 1], [5, 10, 4],
                    [2, 10, 5], [3, 9, 0], [3, 11, 4], [6, 11, 3], [7, 11, 6], [4, 11, 7]]
        part = [5,6,5]

        true_ghost_layer_width=2
        true_ghost_recv_dict_1=[array([5, 6, 7, 8, 9]), array([ 5,  6,  8,  9, 10])]
        true_ghost_recv_dict_2=[array([10, 11, 12, 13, 14]), array([11, 12, 13, 14, 15])]
        true_vertices=array([[ 1,  5,  0],
            [ 1,  6,  2],
            [ 3,  6,  1],
            [ 4,  6,  3],
            [ 2,  6,  4],
            [ 7,  5,  8],
            [ 8,  5,  1],
            [ 1, 11,  8],
            [ 2, 11,  1],
            [ 9, 11,  2],
            [ 0,  5,  7],
            [ 0, 12,  1],
            [10, 12,  0],
            [ 3, 12, 10],
            [ 1, 12,  3]])
        true_points=array([[ 0.5 ,  0.  ],
            [ 0.5 ,  0.5 ],
            [ 0.5 ,  1.  ],
            [ 1.  ,  0.5 ],
            [ 1.  ,  1.  ],
            [ 0.25,  0.25],
            [ 0.75,  0.75],
            [ 0.  ,  0.  ],
            [ 0.  ,  0.5 ],
            [ 0.  ,  1.  ],
            [ 1.  ,  0.  ],
            [ 0.25,  0.75],
            [ 0.75,  0.25]])
        true_full_send_dict_1=[array([0, 1, 2, 4]), array([0, 1, 2, 4])]
        true_full_send_dict_2=[array([0, 1, 2, 3]), array([0, 1, 2, 3])]

    if metis_version == '5_part_graph':
        triangles = [[0,  9,  1], [3,  9,  0], [4,  9,  3], [1,  9,  4],
                      [3, 11,  4], [1, 10,  2], [4, 10,  1], [5, 10,  4],
                      [2, 10,  5], [4, 12,  5], [6, 11,  3], [7, 11,  6],
                      [4, 11,  7], [7, 12,  4], [8, 12,  7], [5, 12,  8]]
        part = [5,5,6]

        true_ghost_recv_dict_1=[array([5, 6, 7]), array([5, 6, 7])]
        true_ghost_layer_width=2
        true_ghost_recv_dict_2=[array([ 8,  9, 10, 11]), array([10, 11, 12, 13])]
        true_vertices=array([[ 0,  4,  1],
            [ 2,  4,  0],
            [ 3,  4,  2],
            [ 1,  4,  3],
            [ 2,  5,  3],
            [ 1, 10,  6],
            [ 3, 10,  1],
            [ 7, 10,  3],
            [ 8,  5,  2],
            [ 9,  5,  8],
            [ 3,  5,  9],
            [ 9, 11,  3]])
        true_points=array([[ 0.  ,  0.  ],
            [ 0.  ,  0.5 ],
            [ 0.5 ,  0.  ],
            [ 0.5 ,  0.5 ],
            [ 0.25,  0.25],
            [ 0.75,  0.25],
            [ 0.  ,  1.  ],
            [ 0.5 ,  1.  ],
            [ 1.  ,  0.  ],
            [ 1.  ,  0.5 ],
            [ 0.25,  0.75],
            [ 0.75,  0.75]])
        true_full_send_dict_1=[array([0, 2, 3]), array([0, 2, 3])]
        true_full_send_dict_2=[array([2, 4]), array([2, 4])]


        # numpy.allclose(ghost_layer_width,  true_ghost_layer_width)
        # numpy.allclose(points,   true_points)
        # numpy.allclose(vertices, true_vertices)
        # numpy.allclose(ghost_recv_dict[1], true_ghost_recv_dict_1)
        # numpy.allclose(ghost_recv_dict[2], true_ghost_recv_dict_2)
        # numpy.allclose(full_send_dict[1], true_full_send_dict_1)
        # numpy.allclose(full_send_dict[2], true_full_send_dict_2)

        #print triangles_per_proc


    if myid == 1:

        from numpy import array
        
        if metis_version == 4:
            true_vertices=array([[ 0,  5,  1],
                [ 1,  5,  3],
                [ 1,  6,  2],
                [ 3,  6,  1],
                [ 4,  6,  3],
                [ 2,  6,  4],
                [ 3,  5,  7],
                [ 3, 11,  4],
                [ 8, 11,  3],
                [ 4, 11,  9],
                [ 7,  5,  0],
                [ 7, 10,  3]])
            true_points=array([[ 0.  ,  0.  ],
                [ 0.  ,  0.5 ],
                [ 0.  ,  1.  ],
                [ 0.5 ,  0.5 ],
                [ 0.5 ,  1.  ],
                [ 0.25,  0.25],
                [ 0.25,  0.75],
                [ 0.5 ,  0.  ],
                [ 1.  ,  0.5 ],
                [ 1.  ,  1.  ],
                [ 0.75,  0.25],
                [ 0.75,  0.75]])
            true_full_send_dict_0=[array([0, 1, 3, 4, 5]), array([ 5,  6,  8,  9, 10])]
            true_node_map=array([ 0,  1,  2,  7,  3,  4, -1,  8,  9,  5,  6, 10, 11])
            true_full_send_dict_2=[array([0, 1]), array([5, 6])]
            true_ghost_recv_dict_0=[array([6, 7, 8, 9]), array([0, 1, 2, 4])]
            true_ghost_recv_dict_2=[array([10, 11]), array([11, 12])]
            true_ghost_layer_width=2
            true_tri_map=array([ 6,  7,  8, -1,  9,  0,  1,  2,  3,  4,  5, 10, 11])

        if metis_version == '5_part_graph':
            true_vertices=array([[ 0,  4,  1],
                [ 2,  4,  0],
                [ 3,  4,  2],
                [ 1,  4,  3],
                [ 2,  5,  3],
                [ 6, 10,  0],
                [ 2, 10,  7],
                [ 0, 10,  2],
                [ 2, 11,  8],
                [ 8,  5,  2],
                [ 9,  5,  8],
                [ 3,  5,  9]])
            true_points=array([[ 0.  ,  0.5 ],
                [ 0.  ,  1.  ],
                [ 0.5 ,  0.5 ],
                [ 0.5 ,  1.  ],
                [ 0.25,  0.75],
                [ 0.75,  0.75],
                [ 0.  ,  0.  ],
                [ 0.5 ,  0.  ],
                [ 1.  ,  0.5 ],
                [ 1.  ,  1.  ],
                [ 0.25,  0.25],
                [ 0.75,  0.25]])
            true_full_send_dict_0=[array([0, 1, 2]), array([5, 6, 7])]
            true_node_map=array([ 6,  0,  1,  7,  2,  3, -1,  8,  9, 10,  4, 11,  5])
            true_full_send_dict_2=[array([2, 4]), array([7, 9])]
            true_ghost_recv_dict_0=[array([5, 6, 7]), array([0, 2, 3])]
            true_ghost_recv_dict_2=[array([ 8,  9, 10, 11]), array([12, 13, 14, 15])]
            true_ghost_layer_width=2
            true_tri_map=array([ 5, -1,  6,  7, -1,  0,  1,  2,  3,  4, -1, -1,  8,  9, 10, 11])

    if myid == 2:

        from numpy import array

        if metis_version == 4:
            true_vertices=array([[ 1,  5,  0],
                [ 1,  6,  2],
                [ 3,  6,  1],
                [ 4,  6,  3],
                [ 2,  6,  4],
                [ 2,  5,  1],
                [ 2, 10,  8],
                [ 4, 10,  2],
                [ 9, 10,  4],
                [ 0,  5,  7],
                [ 7,  5,  2]])
            true_points=array([[ 0.  ,  0.  ],
                [ 0.5 ,  0.  ],
                [ 0.5 ,  0.5 ],
                [ 1.  ,  0.  ],
                [ 1.  ,  0.5 ],
                [ 0.25,  0.25],
                [ 0.75,  0.25],
                [ 0.  ,  0.5 ],
                [ 0.5 ,  1.  ],
                [ 1.  ,  1.  ],
                [ 0.75,  0.75]])
            true_full_send_dict_0=[array([0, 1, 2, 3, 4]), array([11, 12, 13, 14, 15])]
            true_full_send_dict_1=[array([0, 1]), array([11, 12])]
            true_node_map=array([ 0,  7, -1,  1,  2,  8,  3,  4,  9,  5, -1,  6, 10])
            true_ghost_recv_dict_1=[array([ 9, 10]), array([5, 6])]
            true_ghost_recv_dict_0=[array([5, 6, 7, 8]), array([0, 1, 2, 3])]
            true_ghost_layer_width=2
            true_tri_map=array([ 5,  6,  7,  8, -1,  9, 10, -1, -1, -1, -1,  0,  1,  2,  3,  4, -1])

        if metis_version == '5_part_graph':
            true_vertices=array([[3, 6, 0],
                [4, 6, 3],
                [1, 6, 4],
                [4, 7, 1],
                [5, 7, 4],
                [2, 7, 5],
                [1, 8, 0],
                [0, 6, 1],
                [2, 9, 1],
                [1, 7, 2]])
            true_points=array([[ 0.5 ,  0.  ],
                [ 0.5 ,  0.5 ],
                [ 0.5 ,  1.  ],
                [ 1.  ,  0.  ],
                [ 1.  ,  0.5 ],
                [ 1.  ,  1.  ],
                [ 0.75,  0.25],
                [ 0.75,  0.75],
                [ 0.25,  0.25],
                [ 0.25,  0.75]])
            true_full_send_dict_0=[array([0, 1, 2, 3]), array([10, 11, 12, 13])]
            true_full_send_dict_1=[array([2, 3, 4, 5]), array([12, 13, 14, 15])]
            true_node_map=array([-1, -1, -1,  0,  1,  2,  3,  4,  5,  8,  9,  6,  7])
            true_ghost_recv_dict_1=[array([8, 9]), array([7, 9])]
            true_ghost_recv_dict_0=[array([6, 7]), array([2, 4])]
            true_ghost_layer_width=2
            true_tri_map=array([-1, -1,  6, -1,  7, -1, -1,  8, -1,  9,  0,  1,  2,  3,  4,  5, -1])

    true_seq_values = dict(vertices = vertices, triangles = triangles, triangles_per_proc = part)
    return true_seq_values


def print_seq_values(vertices, triangles, triangles_per_proc):

    values = dict(vertices = vertices, triangles = triangles, triangles_per_proc = triangles_per_proc)
    print("seq_values")
    pprint(values)

    print('end of seq_values')

    return


def print_submesh_values(submesh):
    from pprint import pformat
    for i in [0,1,2]:
        parms = [ 'full_nodes',
            'ghost_nodes',
            'full_triangles',
            'ghost_triangles',
            'ghost_commun'
        ]
        for parm in parms:
            name = "submesh['"+parm+"']["+str(i)+"]"
            value = eval(name)
            msg = parm + '_'+str(i)+'='+ pformat(value) + ','
            msg = str(msg)
            print(msg)
    value = submesh['full_commun']
    msg = 'full_commun='+ pformat(value)
    print(msg)

def get_true_submesh_values(metis_version = 4):
    
    if sys.platform == 'win32'  and metis_version == 4:

        true_values = dict( \
            full_nodes_0=array([[ 0.  ,  0.  ,  0.  ],
                [ 1.  ,  0.  ,  0.5 ],
                [ 2.  ,  0.  ,  1.  ],
                [ 4.  ,  0.5 ,  0.5 ],
                [ 5.  ,  0.5 ,  1.  ],
                [ 9.  ,  0.25,  0.25],
                [10.  ,  0.25,  0.75]]),
            ghost_nodes_0=array([[ 3.  ,  0.5 ,  0.  ],
                [12.  ,  0.75,  0.75]]),
            full_triangles_0=array([[ 0,  9,  1],
                [ 1, 10,  2],
                [ 4, 10,  1],
                [ 2, 10,  5]]),
            ghost_triangles_0=array([[ 4,  3,  9,  0],
                [10,  4,  9,  3],
                [11,  1,  9,  4],
                [12,  5, 10,  4],
                [14,  4, 12,  5]]),
            ghost_commun_0=array([[ 4,  1],
                [10,  2],
                [11,  2],
                [12,  2],
                [14,  2]]),
            full_nodes_1=array([[ 0.  ,  0.  ,  0.  ],
                [ 3.  ,  0.5 ,  0.  ],
                [ 4.  ,  0.5 ,  0.5 ],
                [ 6.  ,  1.  ,  0.  ],
                [ 7.  ,  1.  ,  0.5 ],
                [ 8.  ,  1.  ,  1.  ],
                [ 9.  ,  0.25,  0.25],
                [11.  ,  0.75,  0.25],
                [12.  ,  0.75,  0.75]]),
            ghost_nodes_1=array([[ 1.  ,  0.  ,  0.5 ],
                [ 5.  ,  0.5 ,  1.  ],
                [10.  ,  0.25,  0.75]]),
            full_triangles_1=array([[ 3,  9,  0],
                [ 6, 11,  3],
                [ 7, 11,  6],
                [ 4, 11,  7],
                [ 7, 12,  4],
                [ 8, 12,  7]]),
            ghost_triangles_1=array([[ 0,  0,  9,  1],
                [10,  4,  9,  3],
                [11,  1,  9,  4],
                [12,  5, 10,  4],
                [13,  3, 11,  4],
                [14,  4, 12,  5],
                [15,  5, 12,  8]]),
            ghost_commun_1=array([[ 0,  0],
                [10,  2],
                [11,  2],
                [12,  2],
                [13,  2],
                [14,  2],
                [15,  2]]),
            full_nodes_2=array([[ 1.  ,  0.  ,  0.5 ],
                [ 3.  ,  0.5 ,  0.  ],
                [ 4.  ,  0.5 ,  0.5 ],
                [ 5.  ,  0.5 ,  1.  ],
                [ 8.  ,  1.  ,  1.  ],
                [ 9.  ,  0.25,  0.25],
                [10.  ,  0.25,  0.75],
                [11.  ,  0.75,  0.25],
                [12.  ,  0.75,  0.75]]),
            ghost_nodes_2=array([[0. , 0. , 0. ],
                [2. , 0. , 1. ],
                [6. , 1. , 0. ],
                [7. , 1. , 0.5]]),
            full_triangles_2=array([[ 4,  9,  3],
                [ 1,  9,  4],
                [ 5, 10,  4],
                [ 3, 11,  4],
                [ 4, 12,  5],
                [ 5, 12,  8]]),
            ghost_triangles_2=array([[ 0,  0,  9,  1],
                [ 1,  1, 10,  2],
                [ 2,  4, 10,  1],
                [ 3,  2, 10,  5],
                [ 4,  3,  9,  0],
                [ 5,  6, 11,  3],
                [ 6,  7, 11,  6],
                [ 7,  4, 11,  7],
                [ 8,  7, 12,  4],
                [ 9,  8, 12,  7]]),
            ghost_commun_2=array([[0, 0],
                [1, 0],
                [2, 0],
                [3, 0],
                [4, 1],
                [5, 1],
                [6, 1],
                [7, 1],
                [8, 1],
                [9, 1]]),
            full_commun=[{0: [1, 2], 1: [2], 2: [2], 3: [2]},
            {4: [0, 2], 5: [2], 6: [2], 7: [2], 8: [2], 9: [2]},
            {10: [0, 1], 11: [0, 1], 12: [0, 1], 13: [1], 14: [0, 1], 15: [1]}] )

        return true_values


    if sys.platform == 'win32'  and metis_version == '5_part_graph':

        true_values = dict( \
            full_nodes_0=array([[ 0.  ,  0.  ,  0.  ],
                [ 1.  ,  0.  ,  0.5 ],
                [ 3.  ,  0.5 ,  0.  ],
                [ 4.  ,  0.5 ,  0.5 ],
                [ 9.  ,  0.25,  0.25],
                [10.  ,  0.25,  0.75]]),
            ghost_nodes_0=array([[ 2.  ,  0.  ,  1.  ],
                [ 5.  ,  0.5 ,  1.  ],
                [ 6.  ,  1.  ,  0.  ],
                [ 7.  ,  1.  ,  0.5 ],
                [11.  ,  0.75,  0.25],
                [12.  ,  0.75,  0.75]]),
            full_triangles_0=array([[ 0,  9,  1],
                [ 3,  9,  0],
                [ 4,  9,  3],
                [ 1,  9,  4],
                [ 4, 10,  1]]),
            ghost_triangles_0=array([[ 5,  3, 11,  4],
                [ 6,  6, 11,  3],
                [ 8,  4, 11,  7],
                [10,  1, 10,  2],
                [11,  5, 10,  4],
                [12,  2, 10,  5],
                [13,  4, 12,  5]]),
            ghost_commun_0=array([[ 5,  1],
                [ 6,  1],
                [ 8,  1],
                [10,  2],
                [11,  2],
                [12,  2],
                [13,  2]]),
            full_nodes_1=array([[ 3.  ,  0.5 ,  0.  ],
                [ 4.  ,  0.5 ,  0.5 ],
                [ 6.  ,  1.  ,  0.  ],
                [ 7.  ,  1.  ,  0.5 ],
                [11.  ,  0.75,  0.25],
                [12.  ,  0.75,  0.75]]),
            ghost_nodes_1=array([[ 0.  ,  0.  ,  0.  ],
                [ 1.  ,  0.  ,  0.5 ],
                [ 5.  ,  0.5 ,  1.  ],
                [ 8.  ,  1.  ,  1.  ],
                [ 9.  ,  0.25,  0.25],
                [10.  ,  0.25,  0.75]]),
            full_triangles_1=array([[ 3, 11,  4],
                [ 6, 11,  3],
                [ 7, 11,  6],
                [ 4, 11,  7],
                [ 7, 12,  4]]),
            ghost_triangles_1=array([[ 1,  3,  9,  0],
                [ 2,  4,  9,  3],
                [ 3,  1,  9,  4],
                [11,  5, 10,  4],
                [13,  4, 12,  5],
                [14,  8, 12,  7],
                [15,  5, 12,  8]]),
            ghost_commun_1=array([[ 1,  0],
                [ 2,  0],
                [ 3,  0],
                [11,  2],
                [13,  2],
                [14,  2],
                [15,  2]]),
            full_nodes_2=array([[ 1.  ,  0.  ,  0.5 ],
                [ 2.  ,  0.  ,  1.  ],
                [ 4.  ,  0.5 ,  0.5 ],
                [ 5.  ,  0.5 ,  1.  ],
                [ 7.  ,  1.  ,  0.5 ],
                [ 8.  ,  1.  ,  1.  ],
                [10.  ,  0.25,  0.75],
                [12.  ,  0.75,  0.75]]),
            ghost_nodes_2=array([[ 9.  ,  0.25,  0.25],
                [11.  ,  0.75,  0.25]]),
            full_triangles_2=array([[ 1, 10,  2],
                [ 5, 10,  4],
                [ 2, 10,  5],
                [ 4, 12,  5],
                [ 8, 12,  7],
                [ 5, 12,  8]]),
            ghost_triangles_2=array([[ 3,  1,  9,  4],
                [ 4,  4, 10,  1],
                [ 8,  4, 11,  7],
                [ 9,  7, 12,  4]]),
            ghost_commun_2=array([[3, 0],
                [4, 0],
                [8, 1],
                [9, 1]]),
            full_commun=[{0: [], 1: [1], 2: [1], 3: [1, 2], 4: [2]},
            {5: [0], 6: [0], 7: [], 8: [0, 2], 9: [2]},
            {10: [0], 11: [0, 1], 12: [0], 13: [0, 1], 14: [1], 15: [1]}] )

        return true_values

    if metis_version == 4:

        true_values = dict(
        full_nodes_0=array([[  3.  ,   0.5 ,   0.  ],
            [  4.  ,   0.5 ,   0.5 ],
            [  5.  ,   0.5 ,   1.  ],
            [  7.  ,   1.  ,   0.5 ],
            [  8.  ,   1.  ,   1.  ],
            [  9.  ,   0.25,   0.25],
            [ 12.  ,   0.75,   0.75]]),
        ghost_nodes_0=array([[  0.  ,   0.  ,   0.  ],
            [  1.  ,   0.  ,   0.5 ],
            [  2.  ,   0.  ,   1.  ],
            [  6.  ,   1.  ,   0.  ],
            [ 10.  ,   0.25,   0.75],
            [ 11.  ,   0.75,   0.25]]),
        full_triangles_0=array([[ 4,  9,  3],
            [ 4, 12,  5],
            [ 7, 12,  4],
            [ 8, 12,  7],
            [ 5, 12,  8]]),
        ghost_triangles_0=array([[ 5,  0,  9,  1],
            [ 6,  1,  9,  4],
            [ 8,  4, 10,  1],
            [ 9,  5, 10,  4],
            [10,  2, 10,  5],
            [11,  3,  9,  0],
            [12,  3, 11,  4],
            [13,  6, 11,  3],
            [14,  7, 11,  6],
            [15,  4, 11,  7]]),
        ghost_commun_0=array([[ 5,  1],
            [ 6,  1],
            [ 8,  1],
            [ 9,  1],
            [10,  1],
            [11,  2],
            [12,  2],
            [13,  2],
            [14,  2],
            [15,  2]]),
        full_nodes_1=array([[  0.  ,   0.  ,   0.  ],
            [  1.  ,   0.  ,   0.5 ],
            [  2.  ,   0.  ,   1.  ],
            [  4.  ,   0.5 ,   0.5 ],
            [  5.  ,   0.5 ,   1.  ],
            [  9.  ,   0.25,   0.25],
            [ 10.  ,   0.25,   0.75]]),
        ghost_nodes_1=array([[  3.  ,   0.5 ,   0.  ],
            [  7.  ,   1.  ,   0.5 ],
            [  8.  ,   1.  ,   1.  ],
            [ 11.  ,   0.75,   0.25],
            [ 12.  ,   0.75,   0.75]]),
        full_triangles_1=array([[ 0,  9,  1],
            [ 1,  9,  4],
            [ 1, 10,  2],
            [ 4, 10,  1],
            [ 5, 10,  4],
            [ 2, 10,  5]]),
        ghost_triangles_1=array([[ 0,  4,  9,  3],
            [ 1,  4, 12,  5],
            [ 2,  7, 12,  4],
            [ 4,  5, 12,  8],
            [11,  3,  9,  0],
            [12,  3, 11,  4]]),
        ghost_commun_1=array([[ 0,  0],
            [ 1,  0],
            [ 2,  0],
            [ 4,  0],
            [11,  2],
            [12,  2]]),
        full_nodes_2=array([[  0.  ,   0.  ,   0.  ],
            [  3.  ,   0.5 ,   0.  ],
            [  4.  ,   0.5 ,   0.5 ],
            [  6.  ,   1.  ,   0.  ],
            [  7.  ,   1.  ,   0.5 ],
            [  9.  ,   0.25,   0.25],
            [ 11.  ,   0.75,   0.25]]),
        ghost_nodes_2=array([[  1.  ,   0.  ,   0.5 ],
            [  5.  ,   0.5 ,   1.  ],
            [  8.  ,   1.  ,   1.  ],
            [ 12.  ,   0.75,   0.75]]),
        full_triangles_2=array([[ 3,  9,  0],
            [ 3, 11,  4],
            [ 6, 11,  3],
            [ 7, 11,  6],
            [ 4, 11,  7]]),
        ghost_triangles_2=array([[ 0,  4,  9,  3],
            [ 1,  4, 12,  5],
            [ 2,  7, 12,  4],
            [ 3,  8, 12,  7],
            [ 5,  0,  9,  1],
            [ 6,  1,  9,  4]]),
        ghost_commun_2=array([[0, 0],
            [1, 0],
            [2, 0],
            [3, 0],
            [5, 1],
            [6, 1]]),
        full_commun = [{0: [1, 2], 1: [1, 2], 2: [1, 2], 3: [2], 4: [1]}, {5: [0, 2], 6: [
        0, 2], 7: [], 8: [0], 9: [0], 10: [0]}, {11: [0, 1], 12: [0, 1], 13: [0], 14: [0], 15: [0]}]
        )
        return true_values

    #===============================================
    if metis_version == '5_part_graph': # get_true_submesh_values
        true_values = dict(
        full_nodes_0=array([[  0.  ,   0.  ,   0.  ],
            [  1.  ,   0.  ,   0.5 ],
            [  3.  ,   0.5 ,   0.  ],
            [  4.  ,   0.5 ,   0.5 ],
            [  9.  ,   0.25,   0.25],
            [ 11.  ,   0.75,   0.25]]),
        ghost_nodes_0=array([[  2.  ,   0.  ,   1.  ],
            [  5.  ,   0.5 ,   1.  ],
            [  6.  ,   1.  ,   0.  ],
            [  7.  ,   1.  ,   0.5 ],
            [ 10.  ,   0.25,   0.75],
            [ 12.  ,   0.75,   0.75]]),
        full_triangles_0=array([[ 0,  9,  1],
            [ 3,  9,  0],
            [ 4,  9,  3],
            [ 1,  9,  4],
            [ 3, 11,  4]]),
        ghost_triangles_0=array([[ 5,  1, 10,  2],
            [ 6,  4, 10,  1],
            [ 7,  5, 10,  4],
            [10,  6, 11,  3],
            [11,  7, 11,  6],
            [12,  4, 11,  7],
            [13,  7, 12,  4]]),
        ghost_commun_0=array([[ 5,  1],
            [ 6,  1],
            [ 7,  1],
            [10,  2],
            [11,  2],
            [12,  2],
            [13,  2]]),
        full_nodes_1=array([[  1.  ,   0.  ,   0.5 ],
            [  2.  ,   0.  ,   1.  ],
            [  4.  ,   0.5 ,   0.5 ],
            [  5.  ,   0.5 ,   1.  ],
            [ 10.  ,   0.25,   0.75],
            [ 12.  ,   0.75,   0.75]]),
        ghost_nodes_1=array([[  0.  ,   0.  ,   0.  ],
            [  3.  ,   0.5 ,   0.  ],
            [  7.  ,   1.  ,   0.5 ],
            [  8.  ,   1.  ,   1.  ],
            [  9.  ,   0.25,   0.25],
            [ 11.  ,   0.75,   0.25]]),
        full_triangles_1=array([[ 1, 10,  2],
            [ 4, 10,  1],
            [ 5, 10,  4],
            [ 2, 10,  5],
            [ 4, 12,  5]]),
        ghost_triangles_1=array([[ 0,  0,  9,  1],
            [ 2,  4,  9,  3],
            [ 3,  1,  9,  4],
            [12,  4, 11,  7],
            [13,  7, 12,  4],
            [14,  8, 12,  7],
            [15,  5, 12,  8]]),
        ghost_commun_1=array([[ 0,  0],
            [ 2,  0],
            [ 3,  0],
            [12,  2],
            [13,  2],
            [14,  2],
            [15,  2]]),
        full_nodes_2=array([[  3.  ,   0.5 ,   0.  ],
            [  4.  ,   0.5 ,   0.5 ],
            [  5.  ,   0.5 ,   1.  ],
            [  6.  ,   1.  ,   0.  ],
            [  7.  ,   1.  ,   0.5 ],
            [  8.  ,   1.  ,   1.  ],
            [ 11.  ,   0.75,   0.25],
            [ 12.  ,   0.75,   0.75]]),
        ghost_nodes_2=array([[  9.  ,   0.25,   0.25],
            [ 10.  ,   0.25,   0.75]]),
        full_triangles_2=array([[ 6, 11,  3],
            [ 7, 11,  6],
            [ 4, 11,  7],
            [ 7, 12,  4],
            [ 8, 12,  7],
            [ 5, 12,  8]]),
        ghost_triangles_2=array([[ 2,  4,  9,  3],
            [ 4,  3, 11,  4],
            [ 7,  5, 10,  4],
            [ 9,  4, 12,  5]]),
        ghost_commun_2=array([[2, 0],
            [4, 0],
            [7, 1],
            [9, 1]]),
        full_commun=[{0: [1], 1: [], 2: [1, 2], 3: [1], 4: [2]}, {5: [0], 
            6: [0], 7: [0, 2], 8: [], 9: [2]}, {10: [0], 11: [0], 12: [0, 1], 13: [0, 1], 14: [1], 15: [1]}]    
        )
        return true_values


    if metis_version == '5_part_mesh':
        true_values = \
            dict(
                full_nodes_0=array([[4.,  0.5,  0.5],
                                    [5.,  0.5,  1.],
                                    [6.,  1.,  0.],
                                    [7.,  1.,  0.5],
                                    [8.,  1.,  1.],
                                    [11.,  0.75,  0.25],
                                    [12.,  0.75,  0.75]]),
                ghost_nodes_0=array([[3.,  0.5,  0.],
                                     [10.,  0.25,  0.75]]),
                full_triangles_0=array([[7, 11,  6],
                                        [7, 12,  4],
                                        [8, 12,  7],
                                        [5, 12,  8]]),
                ghost_triangles_0=array([[9,  4, 11,  7],
                                         [10,  4, 12,  5],
                                         [12,  5, 10,  4],
                                         [14,  3, 11,  4],
                                         [15,  6, 11,  3]]),
                ghost_commun_0=array([[9,  1],
                                      [10,  1],
                                      [12,  2],
                                      [14,  2],
                                      [15,  2]]),
                full_nodes_1=array([[0.,  0.,  0.],
                                    [1.,  0.,  0.5],
                                    [3.,  0.5,  0.],
                                    [4.,  0.5,  0.5],
                                    [5.,  0.5,  1.],
                                    [7.,  1.,  0.5],
                                    [9.,  0.25,  0.25],
                                    [10.,  0.25,  0.75],
                                    [11.,  0.75,  0.25],
                                    [12.,  0.75,  0.75]]),
                ghost_nodes_1=array([[2., 0., 1.],
                                     [6., 1., 0.],
                                     [8., 1., 1.]]),
                full_triangles_1=array([[0,  9,  1],
                                        [3,  9,  0],
                                        [4,  9,  3],
                                        [1,  9,  4],
                                        [4, 10,  1],
                                        [4, 11,  7],
                                        [4, 12,  5]]),
                ghost_triangles_1=array([[0,  7, 11,  6],
                                         [1,  7, 12,  4],
                                         [2,  8, 12,  7],
                                         [3,  5, 12,  8],
                                         [11,  1, 10,  2],
                                         [12,  5, 10,  4],
                                         [13,  2, 10,  5],
                                         [14,  3, 11,  4],
                                         [15,  6, 11,  3]]),
                ghost_commun_1=array([[0,  0],
                                      [1,  0],
                                      [2,  0],
                                      [3,  0],
                                      [11,  2],
                                      [12,  2],
                                      [13,  2],
                                      [14,  2],
                                      [15,  2]]),
                full_nodes_2=array([[1.,  0.,  0.5],
                                    [2.,  0.,  1.],
                                    [3.,  0.5,  0.],
                                    [4.,  0.5,  0.5],
                                    [5.,  0.5,  1.],
                                    [6.,  1.,  0.],
                                    [10.,  0.25,  0.75],
                                    [11.,  0.75,  0.25]]),
                ghost_nodes_2=array([[0.,  0.,  0.],
                                     [7.,  1.,  0.5],
                                     [8.,  1.,  1.],
                                     [9.,  0.25,  0.25],
                                     [12.,  0.75,  0.75]]),
                full_triangles_2=array([[1, 10,  2],
                                        [5, 10,  4],
                                        [2, 10,  5],
                                        [3, 11,  4],
                                        [6, 11,  3]]),
                ghost_triangles_2=array([[0,  7, 11,  6],
                                         [1,  7, 12,  4],
                                         [3,  5, 12,  8],
                                         [5,  3,  9,  0],
                                         [6,  4,  9,  3],
                                         [7,  1,  9,  4],
                                         [8,  4, 10,  1],
                                         [9,  4, 11,  7],
                                         [10,  4, 12,  5]]),
                ghost_commun_2=array([[0,  0],
                                      [1,  0],
                                      [3,  0],
                                      [5,  1],
                                      [6,  1],
                                      [7,  1],
                                      [8,  1],
                                      [9,  1],
                                      [10,  1]]),
                full_commun=[{0: [1, 2], 1: [1, 2], 2: [1], 3: [1, 2]},
                             {4: [], 5: [2], 6: [2], 7: [2],
                             8: [2], 9: [0, 2], 10: [0, 2]},
                             {11: [1], 12: [0, 1], 13: [1], 14: [0, 1], 15: [0, 1]}])
    
        return true_values
                



def print_extract_submesh(points, triangles, ghost_recv_dict, full_send_dict, \
                         tri_map, node_map, ghost_layer_width):

    values = dict(
        points = points,
        triangles = triangles,
        ghost_layer_width=ghost_layer_width,
        ghost_recv_dict_1= ghost_recv_dict[1],
        ghost_recv_dict_2= ghost_recv_dict[2],
        full_send_dict_1=full_send_dict[1],
        full_send_dict_2=full_send_dict[2],
        tri_map=tri_map,
        node_map=node_map)

    pprint(values)

def print_extract_submesh_1(points, triangles, ghost_recv_dict, full_send_dict, \
                         tri_map, node_map, ghost_layer_width):

    values = dict(
        points = points,
        triangles = triangles,
        ghost_layer_width=ghost_layer_width,
        ghost_recv_dict_1= ghost_recv_dict[1],
        ghost_recv_dict_2= ghost_recv_dict[2],
        full_send_dict_1=full_send_dict[1],
        full_send_dict_2=full_send_dict[2],
        tri_map=tri_map,
        node_map=node_map)

    pprint(values)

def get_true_extract_submesh(metis_version=4):

    if sys.platform == 'win32' and metis_version == 4:
        true_values = \
            {'full_send_dict_1': [array([0]), array([0])],
            'full_send_dict_2': [array([0, 1, 2, 3]),
                                array([0, 1, 2, 3])],
            'ghost_layer_width': 2,
            'ghost_recv_dict_1': [array([4]), array([4])],
            'ghost_recv_dict_2': [array([5, 6, 7, 8]),
                                array([10, 11, 12, 14])],
            'node_map': array([ 0,  1,  2,  7,  3,  4, -1, -1, -1,  5,  6, -1,  8]),
            'points': array([[0.  , 0.  ],
                [0.  , 0.5 ],
                [0.  , 1.  ],
                [0.5 , 0.5 ],
                [0.5 , 1.  ],
                [0.25, 0.25],
                [0.25, 0.75],
                [0.5 , 0.  ],
                [0.75, 0.75]]),
            'tri_map': array([ 0,  1,  2,  3,  4, -1, -1, -1, -1, -1,  5,  6,  7, -1,  8],
                dtype=int),
            'triangles': array([[ 0,  9,  1],
                [ 1, 10,  2],
                [ 4, 10,  1],
                [ 2, 10,  5],
                [ 3,  9,  0],
                [ 6, 11,  3],
                [ 7, 11,  6],
                [ 4, 11,  7],
                [ 7, 12,  4],
                [ 8, 12,  7],
                [ 4,  9,  3],
                [ 1,  9,  4],
                [ 5, 10,  4],
                [ 3, 11,  4],
                [ 4, 12,  5],
                [ 5, 12,  8]])}

        return true_values

    if sys.platform == 'win32' and metis_version == '5_part_graph':
        true_values = \
            {'full_send_dict_1': [array([1, 2, 3]), array([1, 2, 3])],
             'full_send_dict_2': [array([3, 4]), array([3, 4])],
             'ghost_layer_width': 2,
             'ghost_recv_dict_1': [array([5, 6, 7]), array([5, 6, 8])],
             'ghost_recv_dict_2': [array([8,  9, 10, 11]),
                                   array([10, 11, 12, 13])],
             'node_map': array([0,  1,  6,  2,  3,  7,  8,  9, -1,  4,  5, 10, 11]),
             'points': array([[0., 0.],
                              [0., 0.5],
                              [0.5, 0.],
                              [0.5, 0.5],
                              [0.25, 0.25],
                              [0.25, 0.75],
                              [0., 1.],
                              [0.5, 1.],
                              [1., 0.],
                              [1., 0.5],
                              [0.75, 0.25],
                              [0.75, 0.75]]),
                'tri_map': array([0,  1,  2,  3,  4,  5,  6, -1,  7, -1,  8,  9, 10, 11]),
                'triangles': array([[0,  9,  1],
                                    [3,  9,  0],
                                    [4,  9,  3],
                                    [1,  9,  4],
                                    [4, 10,  1],
                                    [3, 11,  4],
                                    [6, 11,  3],
                                    [7, 11,  6],
                                    [4, 11,  7],
                                    [7, 12,  4],
                                    [1, 10,  2],
                                    [5, 10,  4],
                                    [2, 10,  5],
                                    [4, 12,  5],
                                    [8, 12,  7],
                                    [5, 12,  8]])}

        return true_values


    if metis_version == 4:
        true_values = \
            {'full_send_dict_1': [array([0, 1, 2, 4]), array([0, 1, 2, 4])],
            'full_send_dict_2': [array([0, 1, 2, 3]), array([0, 1, 2, 3])],
            'ghost_layer_width': 2,
            'ghost_recv_dict_1': [array([5, 6, 7, 8, 9]), array([ 5,  6,  8,  9, 10])],
            'ghost_recv_dict_2': [array([10, 11, 12, 13, 14]), array([11, 12, 13, 14, 15])],
            'node_map': array([ 7,  8,  9,  0,  1,  2, 10,  3,  4,  5, 11, 12,  6]),
            'points': array([[0.5 , 0.  ],
                [0.5 , 0.5 ],
                [0.5 , 1.  ],
                [1.  , 0.5 ],
                [1.  , 1.  ],
                [0.25, 0.25],
                [0.75, 0.75],
                [0.  , 0.  ],
                [0.  , 0.5 ],
                [0.  , 1.  ],
                [1.  , 0.  ],
                [0.25, 0.75],
                [0.75, 0.25]]),
            'tri_map': array([ 0,  1,  2,  3,  4,  5,  6, -1,  7,  8,  9, 10, 11, 12, 13, 14]),
            'triangles': array([[ 4,  9,  3],
                [ 4, 12,  5],
                [ 7, 12,  4],
                [ 8, 12,  7],
                [ 5, 12,  8],
                [ 0,  9,  1],
                [ 1,  9,  4],
                [ 1, 10,  2],
                [ 4, 10,  1],
                [ 5, 10,  4],
                [ 2, 10,  5],
                [ 3,  9,  0],
                [ 3, 11,  4],
                [ 6, 11,  3],
                [ 7, 11,  6],
                [ 4, 11,  7]])}

    if metis_version == '5_part_graph': # get_true_extract_submesh
        true_values = \
            {'full_send_dict_1': [array([0, 2, 3]), array([0, 2, 3])],
            'full_send_dict_2': [array([2, 4]), array([2, 4])],
            'ghost_layer_width': 2,
            'ghost_recv_dict_1': [array([5, 6, 7]), array([5, 6, 7])],
            'ghost_recv_dict_2': [array([ 8,  9, 10, 11]), array([10, 11, 12, 13])],
            'node_map': array([ 0,  1,  6,  2,  3,  7,  8,  9, -1,  4, 10,  5, 11]),
            'points': array([[0.  , 0.  ],
                [0.  , 0.5 ],
                [0.5 , 0.  ],
                [0.5 , 0.5 ],
                [0.25, 0.25],
                [0.75, 0.25],
                [0.  , 1.  ],
                [0.5 , 1.  ],
                [1.  , 0.  ],
                [1.  , 0.5 ],
                [0.25, 0.75],
                [0.75, 0.75]]),
            'tri_map': array([ 0,  1,  2,  3,  4,  5,  6,  7, -1, -1,  8,  9, 10, 11]),
            'triangles': array([[ 0,  9,  1],
                [ 3,  9,  0],
                [ 4,  9,  3],
                [ 1,  9,  4],
                [ 3, 11,  4],
                [ 1, 10,  2],
                [ 4, 10,  1],
                [ 5, 10,  4],
                [ 2, 10,  5],
                [ 4, 12,  5],
                [ 6, 11,  3],
                [ 7, 11,  6],
                [ 4, 11,  7],
                [ 7, 12,  4],
                [ 8, 12,  7],
                [ 5, 12,  8]])}


        return true_values

    if metis_version == '5_part_mesh':

        true_values = \
            {'full_send_dict_1': [array([0, 1, 2, 3]), array([0, 1, 2, 3])],
             'full_send_dict_2': [array([0, 1, 3]), array([0, 1, 3])],
             'ghost_layer_width': 2,
             'ghost_recv_dict_1': [array([4, 5]), array([9, 10])],
             'ghost_recv_dict_2': [array([6, 7, 8]), array([12, 14, 15])],
             'node_map': array([-1, -1, -1,  7,  0,  1,  2,  3,  4, -1,  8,  5,  6]),
             'points': array([[0.5, 0.5],
                              [0.5, 1.],
                              [1., 0.],
                              [1., 0.5],
                              [1., 1.],
                              [0.75, 0.25],
                              [0.75, 0.75],
                              [0.5, 0.],
                              [0.25, 0.75]]),
                'tri_map': array([0,  1,  2,  3, -1, -1, -1, -1, -1,  4,  5, -1,  6, -1,  7,  8]),
                'triangles': array([[7, 11,  6],
                                    [7, 12,  4],
                                    [8, 12,  7],
                                    [5, 12,  8],
                                    [0,  9,  1],
                                    [3,  9,  0],
                                    [4,  9,  3],
                                    [1,  9,  4],
                                    [4, 10,  1],
                                    [4, 11,  7],
                                    [4, 12,  5],
                                    [1, 10,  2],
                                    [5, 10,  4],
                                    [2, 10,  5],
                                    [3, 11,  4],
                                    [6, 11,  3]])}

        return true_values

def print_rec_submesh_1(points, triangles, ghost_recv_dict, full_send_dict, \
                         tri_map, node_map, ghost_layer_width):

    values = dict(
        points = points,
        triangles = triangles,
        ghost_layer_width=ghost_layer_width,
        ghost_recv_dict_0= ghost_recv_dict[0],
        ghost_recv_dict_2= ghost_recv_dict[2],
        full_send_dict_0=full_send_dict[0],
        full_send_dict_2=full_send_dict[2],
        tri_map=tri_map,
        node_map=node_map)

    pprint(values)


def get_true_rec_submesh_1(metis_version=4):

    if sys.platform == 'win32' and metis_version == 4:
        true_values = \
        {'full_send_dict_0': [array([0]), array([4])],
        'full_send_dict_2': [array([0, 1, 2, 3, 4, 5]),
                            array([4, 5, 6, 7, 8, 9])],
        'ghost_layer_width': 2,
        'ghost_recv_dict_0': [array([6]), array([0])],
        'ghost_recv_dict_2': [array([ 7,  8,  9, 10, 11, 12]),
                            array([10, 11, 12, 13, 14, 15])],
        'node_map': array([ 0,  9, -1,  1,  2, 10,  3,  4,  5,  6, 11,  7,  8]),
        'points': array([[0.  , 0.  ],
            [0.5 , 0.  ],
            [0.5 , 0.5 ],
            [1.  , 0.  ],
            [1.  , 0.5 ],
            [1.  , 1.  ],
            [0.25, 0.25],
            [0.75, 0.25],
            [0.75, 0.75],
            [0.  , 0.5 ],
            [0.5 , 1.  ],
            [0.25, 0.75]]),
        'tri_map': array([ 6, -1, -1, -1,  0,  1,  2,  3,  4,  5,  7,  8,  9, 10, 11, 12],
            dtype=int),
        'triangles': array([[ 1,  6,  0],
            [ 3,  7,  1],
            [ 4,  7,  3],
            [ 2,  7,  4],
            [ 4,  8,  2],
            [ 5,  8,  4],
            [ 0,  6,  9],
            [ 2,  6,  1],
            [ 9,  6,  2],
            [10, 11,  2],
            [ 1,  7,  2],
            [ 2,  8, 10],
            [10,  8,  5]])}

        return true_values

    if sys.platform == 'win32' and metis_version == '5_part_graph':
        rec_submesh_1 = \
            {'full_send_dict_0': [array([0, 1, 3]), array([5, 6, 8])],
             'full_send_dict_2': [array([3, 4]), array([8, 9])],
             'ghost_layer_width': 2,
             'ghost_recv_dict_0': [array([5, 6, 7]), array([1, 2, 3])],
             'ghost_recv_dict_2': [array([8,  9, 10, 11]), array([11, 13, 14, 15])],
             'node_map': array([6,  7, -1,  0,  1,  8,  2,  3,  9, 10, 11,  4,  5]),
             'points': array([[0.5, 0.],
                              [0.5, 0.5],
                              [1., 0.],
                              [1., 0.5],
                              [0.75, 0.25],
                              [0.75, 0.75],
                              [0., 0.],
                              [0., 0.5],
                              [0.5, 1.],
                              [1., 1.],
                              [0.25, 0.25],
                              [0.25, 0.75]]),
                'tri_map': array([-1,  5,  6,  7, -1,  0,  1,  2,  3,  4, -1,  8, -1,  9, 10, 11]),
                'triangles': array([[0,  4,  1],
                                    [2,  4,  0],
                                    [3,  4,  2],
                                    [1,  4,  3],
                                    [3,  5,  1],
                                    [0, 10,  6],
                                    [1, 10,  0],
                                    [7, 10,  1],
                                    [8, 11,  1],
                                    [1,  5,  8],
                                    [9,  5,  3],
                                    [8,  5,  9]])}

    if metis_version == 4:
        true_values = \
            {'full_send_dict_0': [array([0, 1, 3, 4, 5]), array([ 5,  6,  8,  9, 10])],
            'full_send_dict_2': [array([0, 1]), array([5, 6])],
            'ghost_layer_width': 2,
            'ghost_recv_dict_0': [array([6, 7, 8, 9]), array([0, 1, 2, 4])],
            'ghost_recv_dict_2': [array([10, 11]), array([11, 12])],
            'node_map': array([ 0,  1,  2,  7,  3,  4, -1,  8,  9,  5,  6, 10, 11]),
            'points': array([[0.  , 0.  ],
                [0.  , 0.5 ],
                [0.  , 1.  ],
                [0.5 , 0.5 ],
                [0.5 , 1.  ],
                [0.25, 0.25],
                [0.25, 0.75],
                [0.5 , 0.  ],
                [1.  , 0.5 ],
                [1.  , 1.  ],
                [0.75, 0.25],
                [0.75, 0.75]]),
            'tri_map': array([ 6,  7,  8, -1,  9,  0,  1,  2,  3,  4,  5, 10, 11]),
            'triangles': array([[ 0,  5,  1],
                [ 1,  5,  3],
                [ 1,  6,  2],
                [ 3,  6,  1],
                [ 4,  6,  3],
                [ 2,  6,  4],
                [ 3,  5,  7],
                [ 3, 11,  4],
                [ 8, 11,  3],
                [ 4, 11,  9],
                [ 7,  5,  0],
                [ 7, 10,  3]])}

    if metis_version == '5_part_graph': # get_true_rec_submesh_1
        true_values = \
            {'full_send_dict_0': [array([0, 1, 2]), array([5, 6, 7])],
            'full_send_dict_2': [array([2, 4]), array([7, 9])],
            'ghost_layer_width': 2,
            'ghost_recv_dict_0': [array([5, 6, 7]), array([0, 2, 3])],
            'ghost_recv_dict_2': [array([ 8,  9, 10, 11]), array([12, 13, 14, 15])],
            'node_map': array([ 6,  0,  1,  7,  2,  3, -1,  8,  9, 10,  4, 11,  5]),
            'points': array([[0.  , 0.5 ],
                [0.  , 1.  ],
                [0.5 , 0.5 ],
                [0.5 , 1.  ],
                [0.25, 0.75],
                [0.75, 0.75],
                [0.  , 0.  ],
                [0.5 , 0.  ],
                [1.  , 0.5 ],
                [1.  , 1.  ],
                [0.25, 0.25],
                [0.75, 0.25]]),
            'tri_map': array([ 5, -1,  6,  7, -1,  0,  1,  2,  3,  4, -1, -1,  8,  9, 10, 11]),
            'triangles': array([[ 0,  4,  1],
                [ 2,  4,  0],
                [ 3,  4,  2],
                [ 1,  4,  3],
                [ 2,  5,  3],
                [ 6, 10,  0],
                [ 2, 10,  7],
                [ 0, 10,  2],
                [ 2, 11,  8],
                [ 8,  5,  2],
                [ 9,  5,  8],
                [ 3,  5,  9]])}

        return true_values

    if metis_version == '5_part_mesh':
        true_values = \
        {'full_send_dict_0': [array([5, 6]), array([ 9, 10])],
        'full_send_dict_2': [array([1, 2, 3, 4, 5, 6]),
                            array([ 5,  6,  7,  8,  9, 10])],
        'ghost_layer_width': 2,
        'ghost_recv_dict_0': [array([ 7,  8,  9, 10]), array([0, 1, 2, 3])],
        'ghost_recv_dict_2': [array([11, 12, 13, 14, 15]),
                            array([11, 12, 13, 14, 15])],
        'node_map': array([ 0,  1, 10,  2,  3,  4, 11,  5, 12,  6,  7,  8,  9]),
        'points': array([[0.  , 0.  ],
            [0.  , 0.5 ],
            [0.5 , 0.  ],
            [0.5 , 0.5 ],
            [0.5 , 1.  ],
            [1.  , 0.5 ],
            [0.25, 0.25],
            [0.25, 0.75],
            [0.75, 0.25],
            [0.75, 0.75],
            [0.  , 1.  ],
            [1.  , 0.  ],
            [1.  , 1.  ]]),
        'tri_map': array([ 7,  8,  9, 10,  0,  1,  2,  3,  4,  5,  6, 11, 12, 13, 14, 15]),
        'triangles': array([[ 0,  6,  1],
            [ 2,  6,  0],
            [ 3,  6,  2],
            [ 1,  6,  3],
            [ 3,  7,  1],
            [ 3,  8,  5],
            [ 3,  9,  4],
            [ 5,  8, 11],
            [ 5,  9,  3],
            [12,  9,  5],
            [ 4,  9, 12],
            [ 1,  7, 10],
            [ 4,  7,  3],
            [10,  7,  4],
            [ 2,  8,  3],
            [11,  8,  2]])}
    return true_values



def print_rec_submesh_2(points, triangles, ghost_recv_dict, full_send_dict, \
                         tri_map, node_map, ghost_layer_width):

    values = dict(
        points = points,
        triangles = triangles,
        ghost_layer_width=ghost_layer_width,
        ghost_recv_dict_0= ghost_recv_dict[0],
        ghost_recv_dict_1= ghost_recv_dict[1],
        full_send_dict_0=full_send_dict[0],
        full_send_dict_1=full_send_dict[1],
        tri_map=tri_map,
        node_map=node_map)

    pprint(values)

def get_true_rec_submesh_2(metis_version=4):

    if sys.platform == 'win32' and metis_version == 4:
        true_values = \
        {'full_send_dict_0': [array([0, 1, 2, 4]),
                            array([10, 11, 12, 14])],
        'full_send_dict_1': [array([0, 1, 2, 3, 4, 5]),
                            array([10, 11, 12, 13, 14, 15])],
        'ghost_layer_width': 2,
        'ghost_recv_dict_0': [array([6, 7, 8, 9]),
                            array([0, 1, 2, 3])],
        'ghost_recv_dict_1': [array([10, 11, 12, 13, 14, 15]),
                            array([4, 5, 6, 7, 8, 9])],
        'node_map': array([ 9,  0, 10,  1,  2,  3, 11, 12,  4,  5,  6,  7,  8]),
        'points': array([[0.  , 0.5 ],
            [0.5 , 0.  ],
            [0.5 , 0.5 ],
            [0.5 , 1.  ],
            [1.  , 1.  ],
            [0.25, 0.25],
            [0.25, 0.75],
            [0.75, 0.25],
            [0.75, 0.75],
            [0.  , 0.  ],
            [0.  , 1.  ],
            [1.  , 0.  ],
            [1.  , 0.5 ]]),
        'tri_map': array([ 6,  7,  8,  9, 10, 11, 12, 13, 14, 15,  0,  1,  2,  3,  4,  5, -1],
            dtype=int),
        'triangles': array([[ 2,  5,  1],
            [ 0,  5,  2],
            [ 3,  6,  2],
            [ 1,  7,  2],
            [ 2,  8,  3],
            [ 3,  8,  4],
            [ 9,  5,  0],
            [ 0,  6, 10],
            [ 2,  6,  0],
            [10,  6,  3],
            [ 1,  5,  9],
            [11,  7,  1],
            [12,  7, 11],
            [ 2,  7, 12],
            [12,  8,  2],
            [ 4,  8, 12]])}

        return true_values

    if sys.platform == 'win32' and metis_version == 4:
        rec_submesh_2 = \
            {'full_send_dict_0': [array([0, 1, 2, 3]), array([10, 11, 12, 13])],
             'full_send_dict_1': [array([1, 3, 4, 5]), array([11, 13, 14, 15])],
             'ghost_layer_width': 2,
             'ghost_recv_dict_0': [array([6, 7]), array([3, 4])],
             'ghost_recv_dict_1': [array([8, 9]), array([8, 9])],
             'node_map': array([-1,  0,  1, -1,  2,  3, -1,  4,  5,  8,  6,  9,  7]),
             'points': array([[0., 0.5],
                              [0., 1.],
                              [0.5, 0.5],
                              [0.5, 1.],
                              [1., 0.5],
                              [1., 1.],
                              [0.25, 0.75],
                              [0.75, 0.75],
                              [0.25, 0.25],
                              [0.75, 0.25]]),
                'tri_map': array([-1, -1, -1,  6,  7, -1, -1, -1,  8,  9,  0,  1,  2,  3,  4,  5, -1]),
                'triangles': array([[0, 6, 1],
                                    [3, 6, 2],
                                    [1, 6, 3],
                                    [2, 7, 3],
                                    [5, 7, 4],
                                    [3, 7, 5],
                                    [0, 8, 2],
                                    [2, 6, 0],
                                    [2, 9, 4],
                                    [4, 7, 2]])}

    if metis_version == 4:  # get_true_rec_submesh_2
        true_values = dict(
        triangles=array([[ 1,  5,  0],
            [ 1,  6,  2],
            [ 3,  6,  1],
            [ 4,  6,  3],
            [ 2,  6,  4],
            [ 2,  5,  1],
            [ 2, 10,  8],
            [ 4, 10,  2],
            [ 9, 10,  4],
            [ 0,  5,  7],
            [ 7,  5,  2]]),
        points=array([[0.  , 0.  ],
            [0.5 , 0.  ],
            [0.5 , 0.5 ],
            [1.  , 0.  ],
            [1.  , 0.5 ],
            [0.25, 0.25],
            [0.75, 0.25],
            [0.  , 0.5 ],
            [0.5 , 1.  ],
            [1.  , 1.  ],
            [0.75, 0.75]]),
        full_send_dict_0=[array([0, 1, 2, 3, 4]), array([11, 12, 13, 14, 15])],
        full_send_dict_1=[array([0, 1]), array([11, 12])],
        node_map=array([ 0,  7, -1,  1,  2,  8,  3,  4,  9,  5, -1,  6, 10]),
        ghost_recv_dict_1=[array([ 9, 10]), array([5, 6])],
        ghost_recv_dict_0=[array([5, 6, 7, 8]), array([0, 1, 2, 3])],
        ghost_layer_width=2,
        tri_map=array([ 5,  6,  7,  8, -1,  9, 10, -1, -1, -1, -1,  0,  1,  2,  3,  4, -1]))

    if metis_version == '5_part_graph': # get_true_rec_submesh_2
        true_values = dict(
        triangles=array([[3, 6, 0],
            [4, 6, 3],
            [1, 6, 4],
            [4, 7, 1],
            [5, 7, 4],
            [2, 7, 5],
            [1, 8, 0],
            [0, 6, 1],
            [2, 9, 1],
            [1, 7, 2]]),
        points=array([[0.5 , 0.  ],
            [0.5 , 0.5 ],
            [0.5 , 1.  ],
            [1.  , 0.  ],
            [1.  , 0.5 ],
            [1.  , 1.  ],
            [0.75, 0.25],
            [0.75, 0.75],
            [0.25, 0.25],
            [0.25, 0.75]]),
        full_send_dict_0=[array([0, 1, 2, 3]), array([10, 11, 12, 13])],
        full_send_dict_1=[array([2, 3, 4, 5]), array([12, 13, 14, 15])],
        node_map=array([-1, -1, -1,  0,  1,  2,  3,  4,  5,  8,  9,  6,  7]),
        ghost_recv_dict_1=[array([8, 9]), array([7, 9])],
        ghost_recv_dict_0=[array([6, 7]), array([2, 4])],
        ghost_layer_width=2,
        tri_map=array([-1, -1,  6, -1,  7, -1, -1,  8, -1,  9,  0,  1,  2,  3,  4,  5, -1]))    

        return true_values

    if metis_version == '5_part_mesh':

        true_values = \
            {'full_send_dict_0': [array([1, 3, 4]), array([12, 14, 15])],
             'full_send_dict_1': [array([0, 1, 2, 3, 4]), array([11, 12, 13, 14, 15])],
             'ghost_layer_width': 2,
             'ghost_recv_dict_0': [array([5, 6, 7]), array([0, 1, 3])],
             'ghost_recv_dict_1': [array([8,  9, 10, 11, 12, 13]),
                                   array([5,  6,  7,  8,  9, 10])],
             'node_map': array([8,  0,  1,  2,  3,  4,  5,  9, 10, 11,  6,  7, 12]),
             'points': array([[0., 0.5],
                              [0., 1.],
                              [0.5, 0.],
                              [0.5, 0.5],
                              [0.5, 1.],
                              [1., 0.],
                              [0.25, 0.75],
                              [0.75, 0.25],
                              [0., 0.],
                              [1., 0.5],
                              [1., 1.],
                              [0.25, 0.25],
                              [0.75, 0.75]]),
                'tri_map': array([5,  6, -1,  7, -1,  8,  9, 10, 11, 12, 13,  0,  1,  2,  3,  4, -1]),
                'triangles': array([[0,  6,  1],
                                    [4,  6,  3],
                                    [1,  6,  4],
                                    [2,  7,  3],
                                    [5,  7,  2],
                                    [9,  7,  5],
                                    [9, 12,  3],
                                    [4, 12, 10],
                                    [2, 11,  8],
                                    [3, 11,  2],
                                    [0, 11,  3],
                                    [3,  6,  0],
                                    [3,  7,  9],
                                    [3, 12,  4]])}

    return true_values

#-------------------------------------------------------------
if __name__ == "__main__":
    
    distibute_three_processors()

#!/usr/bin/env python
#########################################################
#
#
# Calculate and print the norms of the domain
#
#
#  The routines defined here are intended for debugging
# use. They print the norms of the quantities in the
# domain. As opposed to the definitions given
# in utiltites.norms these calculations take a
# parallel domain into account.
#
#  Authors: Linda Stals and April 2006
#
#
#########################################################


from builtins import range
import sys

from anuga.utilities import parallel_abstraction as pypar

from numpy import array, zeros, ones, take, nonzero, float
from anuga.utilities.norms import l1_norm, l2_norm, linf_norm

#########################################################
#
# Find out which triangles are full triangles (only these
# triangles should be included in the norm calculations)
#
# *)
#
# -------------------------------------------------------
#
# *) A 1-D array, tri_full_flag, is returned.
#
# *) The size of tri_full_flag is the same as the number
# of vertices in the domain
#
# *) If tri_full_flag[i] = 1, then triangle number i is
# a full triangle, if tri_full_flag[i] = 0  the triangle
# is a ghost triangle
#
#########################################################

def build_full_flag(domain, ghost_recv_dict):

    tri_full_flag = ones(len(domain.get_triangles()), Int8)
    for i in list(ghost_recv_dict.keys()):
        for id in ghost_recv_dict[i][0]:
            tri_full_flag[id] = 0
        

    return tri_full_flag

#########################################################
#
# Print the l1 norm of the given quantity
#
# *) The quantity is an array containing three columns
#
# -------------------------------------------------------
#
# *) The l1 norm is calculated along each axis
#
# *) The l1 norm is printed
#
# *) Processor 0 prints the results
#
#
#########################################################
def print_l1_stats(full_edge):

    numprocs = pypar.size()
    myid = pypar.rank()
    
    tri_norm = zeros(3, Float)
    recv_norm = zeros(3, Float)
    tri_norm[0] = l1_norm(full_edge[:, 0])
    tri_norm[1] = l1_norm(full_edge[:, 1])
    tri_norm[2] = l1_norm(full_edge[:, 2])
    if myid == 0:
        for p in range(numprocs-1):
            pypar.receive(p+1, recv_norm)
            tri_norm[0] = tri_norm[0]+recv_norm[0]
            tri_norm[1] = tri_norm[1]+recv_norm[1]
            tri_norm[2] = tri_norm[2]+recv_norm[2]
        print('l1_norm along each axis : [', tri_norm[0],', ', tri_norm[1], ', ', tri_norm[2], ']')

    else:
        pypar.send(tri_norm, 0)

#########################################################
#
# Print the l2 norm of the given quantity
#
# *) The quantity is an array containing three columns
#
# -------------------------------------------------------
#
# *) The l2 norm is calculated along each axis
#
# *) The l2 norm is printed
#
# *) Processor 0 prints the results
#
#
#########################################################
def print_l2_stats(full_edge):

    numprocs = pypar.size()
    myid = pypar.rank()
    
    tri_norm = zeros(3, Float)
    recv_norm = zeros(3, Float)
    tri_norm[0] = pow(l2_norm(full_edge[:, 0]), 2)
    tri_norm[1] = pow(l2_norm(full_edge[:, 1]), 2)
    tri_norm[2] = pow(l2_norm(full_edge[:, 2]), 2)
    if myid == 0:
        for p in range(numprocs-1):
            pypar.receive(p+1, recv_norm)
            tri_norm[0] = tri_norm[0]+recv_norm[0]
            tri_norm[1] = tri_norm[1]+recv_norm[1]
            tri_norm[2] = tri_norm[2]+recv_norm[2]
        print('l2_norm along each axis : [', pow(tri_norm[0], 0.5),', ', pow(tri_norm[1], 0.5), \
              ', ', pow(tri_norm[2], 0.5), ']')
    else:
        pypar.send(tri_norm, 0)


#########################################################
#
# Print the linf norm of the given quantity
#
# *) The quantity is an array containing three columns
#
# -------------------------------------------------------
#
# *) The linf norm is calculated along each axis
#
# *) The linf norm is printed
#
# *) Processor 0 prints the results
#
#
#########################################################
def print_linf_stats(full_edge):

    numprocs = pypar.size()
    myid = pypar.rank()
    
    tri_norm = zeros(3, Float)
    recv_norm = zeros(3, Float)
        
    tri_norm[0] = linf_norm(full_edge[:, 0])
    tri_norm[1] = linf_norm(full_edge[:, 1])
    tri_norm[2] = linf_norm(full_edge[:, 2])
    if myid == 0:
        for p in range(numprocs-1):
            pypar.receive(p+1, recv_norm)
            tri_norm[0] = max(tri_norm[0], recv_norm[0])
            tri_norm[1] = max(tri_norm[1], recv_norm[1])
            tri_norm[2] = max(tri_norm[2], recv_norm[2])
        print('linf_norm along each axis : [', tri_norm[0],', ', tri_norm[1], ', ', tri_norm[2], ']')
    else:
        pypar.send(tri_norm, 0)

       
#########################################################
#
# Print the norms of the quantites assigned to the domain
# (this is useful for checking the numerical results
# in the parallel computation)
#
# *) tri_full_flag states which of the triangles are
# full triangles
#
#
# -------------------------------------------------------
#
# *) For each quantity, the l1, l2 and linf norms are
# printed
#
# *) The size of tri_full_flag is the same as the number
# of vertices in the domain
#
# *) Only the full triangles are used in the norm
# calculations
#
# *) Processor 0 prints the results
#
#########################################################
def print_test_stats(domain, tri_full_flag):

    myid = pypar.rank()

    for k in list(domain.quantities.keys()):
        TestStage = domain.quantities[k]
        if myid == 0:
            print(" ===== ", k, " ===== ")
        full_edge = take(TestStage.edge_values, nonzero(tri_full_flag))
        print_l1_stats(full_edge)
        print_l2_stats(full_edge)
        print_linf_stats(full_edge)


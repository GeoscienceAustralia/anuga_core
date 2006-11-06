#########################################################
#
# Handle the communication between the host machine
# (processor 0) and the processors. The host machine is
# responsible for the doing the initial grid partitioning.
#
# The routines given below should be moved to the
# build_submesh.py and build_local.py file to allow
# overlapping of  communication and computation.
# This should be done after more debugging.
#
#
#  Author: Linda Stals, June 2005
#  Modified: Linda Stals, Nov 2005 (optimise python code)
#
#
#########################################################

from Numeric import array, Int, Float, zeros
import logging, logging.config
logger = logging.getLogger('parallel')
logger.setLevel(logging.WARNING)

try:
    logging.config.fileConfig('log.ini')
except:
    pass

import sys

import pypar

from build_local import build_local_mesh

#########################################################
#
# Send the submesh to processor p.
#
# *) The order and form is strongly coupled with
# rec_submesh.
#
# -------------------------------------------------------
#
# *) All of the information has been sent to processor p.
#
#########################################################

def send_submesh(submesh, triangles_per_proc, p):

    myid = pypar.rank()
    print 'process %d sending submesh to process %d' %(myid, p)
    
    # build and send the tagmap for the boundary conditions
    
    tagmap = {}
    counter = 1
    for b in submesh["full_boundary"][p]:
         bkey = submesh["full_boundary"][p][b]
         if not tagmap.has_key(bkey):
             tagmap[bkey] = counter
             counter = counter+1
    for b in submesh["ghost_boundary"][p]:
         bkey = submesh["ghost_boundary"][p][b]
         if not tagmap.has_key(bkey):
             tagmap[bkey] = counter
             counter = counter+1
    pypar.send(tagmap, p)

    # send the quantities key information
    
    pypar.send(submesh["full_quan"].keys(), p)
    
    # send the number of triangles per processor

    pypar.send(triangles_per_proc, p, use_buffer=True)

    # compress full_commun

    flat_full_commun = []

    for c in submesh["full_commun"][p]:
        for i in range(len(submesh["full_commun"][p][c])):
            flat_full_commun.append([c,submesh["full_commun"][p][c][i]])

    # send the array sizes so memory can be allocated

    setup_array = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    setup_array[0] = len(submesh["full_nodes"][p])
    setup_array[1] = len(submesh["ghost_nodes"][p])
    setup_array[2] = len(submesh["full_triangles"][p])
    setup_array[3] = len(submesh["ghost_triangles"][p])
    setup_array[4] = len(submesh["full_boundary"][p])
    setup_array[5] = len(submesh["ghost_boundary"][p])
    setup_array[6] = len(submesh["ghost_commun"][p])
    setup_array[7] = len(flat_full_commun)

    pypar.send(setup_array, p)
    
    # send the nodes
    
    pypar.send(submesh["full_nodes"][p], p, use_buffer=True)
    pypar.send(submesh["ghost_nodes"][p], p, use_buffer=True)

    # send the triangles

    pypar.send(array(submesh["full_triangles"][p], Int), p, use_buffer=True)
    pypar.send(submesh["ghost_triangles"][p], p, use_buffer=True)

    # send the boundary

    bc = []
    for b in submesh["full_boundary"][p]:
        bc.append([b[0], b[1], tagmap[submesh["full_boundary"][p][b]]])
    pypar.send(bc, p, use_buffer=True)
    bc = []
    for b in submesh["ghost_boundary"][p]:
        bc.append([b[0], b[1], tagmap[submesh["ghost_boundary"][p][b]]])
    pypar.send(bc, p, use_buffer=True)

    # send the communication pattern

    pypar.send(submesh["ghost_commun"][p], p, use_buffer=True)
    pypar.send(flat_full_commun, p, use_buffer=True)

    # send the quantities
    
    for k in submesh["full_quan"]:
        pypar.send(submesh["full_quan"][k][p], p, use_buffer=True)
        
    for k in submesh["ghost_quan"]:
        pypar.send(submesh["ghost_quan"][k][p], p, use_buffer=True)
        

#########################################################
#
# Receive the submesh from processor p.
#
# *) The order and form is strongly coupled with
# send_submesh.
#
# -------------------------------------------------------
#
# *) All of the information has been received by the
# processor p and passed into build_local.
#
# *) The information is returned in a form needed by the
# GA datastructure.
#
#########################################################

def rec_submesh_flat(p):

    numproc = pypar.size()
    myid = pypar.rank()

    submesh_cell = {}
    
    print 'process %d receiving submesh from process %d' %(myid, p)

    # receive the tagmap for the boundary conditions
    
    tagmap = pypar.receive(p)
    itagmap = {}
    for t in tagmap:
        itagmap[tagmap[t]]=t

    # receive the quantities key information
    
    qkeys = pypar.receive(p)

    # receive the number of triangles per processor

    triangles_per_proc = []
    for i in range(numproc):
        triangles_per_proc.append([0])

    triangles_per_proc = pypar.receive(p, triangles_per_proc)

    # recieve information about the array sizes

    setup_array = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    setup_array = pypar.receive(p, setup_array)

    # receive the full nodes

    no_full_nodes = setup_array[0]
    full_nodes = zeros((no_full_nodes, 3), Float)
    submesh_cell["full_nodes"] = pypar.receive(p, full_nodes)
    
    # receive the ghost nodes

    no_ghost_nodes = setup_array[1]
    ghost_nodes = zeros((no_ghost_nodes, 3), Float)
    submesh_cell["ghost_nodes"] = pypar.receive(p, ghost_nodes)

    
    # receive the full triangles

    no_full_triangles = setup_array[2]
    full_triangles = zeros((no_full_triangles, 3), Int)
    submesh_cell["full_triangles"] = pypar.receive(p, full_triangles)
    
    # receive the ghost triangles

    no_ghost_triangles = setup_array[3]
    ghost_triangles = zeros((no_ghost_triangles, 4), Int)
    submesh_cell["ghost_triangles"] = pypar.receive(p, ghost_triangles)
    
    # receive the full boundary

    no_full_boundary = setup_array[4]
    bc = []
    for i in range(no_full_boundary):
        bc.append([0.0, 0.0, 0.0])
    bnd_c = pypar.receive(p, bc)

    submesh_cell["full_boundary"] = {}
    for b in bnd_c:
        submesh_cell["full_boundary"][b[0],b[1]]=itagmap[b[2]]

    # receive the ghost boundary

    no_ghost_boundary = setup_array[5]
    bc = []
    for i in range(no_ghost_boundary):
        bc.append([0.0, 0.0, 0.0])
    bnd_c = pypar.receive(p, bc)

    submesh_cell["ghost_boundary"] = {}
    for b in bnd_c:
        submesh_cell["ghost_boundary"][b[0],b[1]]=itagmap[b[2]]

    # receive the ghost communication pattern

    no_ghost_commun = setup_array[6]
    ghost_commun = zeros((no_ghost_commun, 2), Int)
    submesh_cell["ghost_commun"] = pypar.receive(p, ghost_commun)
    
    # receive the full communication pattern

    no_full_commun = setup_array[7]
    full_commun = []
    for i in range(no_full_commun):
        full_commun.append([0.0, 0.0])

    full_commun = pypar.receive(p, full_commun)

    submesh_cell["full_commun"] = {}
    for c in full_commun:
        submesh_cell["full_commun"][c[0]] = []
    for c in full_commun:
        submesh_cell["full_commun"][c[0]].append(c[1])

    # receive the quantities

    no_quantities = len(qkeys)
    new_quan = zeros((no_full_triangles, 3), Float)
    submesh_cell["full_quan"]={}
    
    for i in range(no_quantities):
        tmp = pypar.receive(p, new_quan)
        submesh_cell["full_quan"][qkeys[i]]=zeros((no_full_triangles,3), Float)
        submesh_cell["full_quan"][qkeys[i]][:] = tmp[:]

    new_quan = zeros((no_ghost_triangles, 3), Float)
    submesh_cell["ghost_quan"]={}
    for i in range(no_quantities):
        tmp = pypar.receive(p, new_quan)
        submesh_cell["ghost_quan"][qkeys[i]]= zeros((no_ghost_triangles,3), Float)
        submesh_cell["ghost_quan"][qkeys[i]][:] = tmp[:]
    
    return submesh_cell, triangles_per_proc,\
           no_full_nodes, no_full_triangles



#########################################################
#
# Receive the submesh from processor p.
#
# *) The order and form is strongly coupled with
# send_submesh.
#
# -------------------------------------------------------
#
# *) All of the information has been received by the
# processor p and passed into build_local.
#
# *) The information is returned in a form needed by the
# GA datastructure.
#
#########################################################

def rec_submesh(p):

    numproc = pypar.size()
    myid = pypar.rank()

    [submesh_cell, triangles_per_proc,\
     number_of_full_nodes, number_of_full_triangles] = rec_submesh_flat(p)
    
    # find the full triangles assigned to this processor

    lower_t = 0
    for i in range(myid):
        lower_t = lower_t+triangles_per_proc[i]
    upper_t = lower_t+triangles_per_proc[myid]

    # convert the information into a form needed by the GA
    # datastructure

    [GAnodes, GAtriangles, boundary, quantities, ghost_rec, full_send] = \
              build_local_mesh(submesh_cell, lower_t, upper_t, \
                               numproc)
    
    return GAnodes, GAtriangles, boundary, quantities,\
           ghost_rec, full_send,\
           number_of_full_nodes, number_of_full_triangles


#########################################################
#
# Extract the submesh that will belong to the
# "host processor" (i.e. processor zero)
#
#  *) See the documentation for build_submesh
#
# -------------------------------------------------------
#
#  *) A dictionary containing the full_triangles,
# full_nodes, full_boundary, ghost_triangles, ghost_nodes,
# ghost_boundary, ghost_commun and full_commun belonging
# to processor zero are returned.
#
#########################################################
def extract_hostmesh(submesh, triangles_per_proc):

    submesh_cell = {}
    submesh_cell["full_nodes"] = submesh["full_nodes"][0]
    submesh_cell["ghost_nodes"] = submesh["ghost_nodes"][0]
    submesh_cell["full_triangles"] = submesh["full_triangles"][0]
    submesh_cell["ghost_triangles"] = submesh["ghost_triangles"][0]
    submesh_cell["full_boundary"] = submesh["full_boundary"][0]
    submesh_cell["ghost_boundary"] = submesh["ghost_boundary"][0]
    submesh_cell["ghost_commun"] = submesh["ghost_commun"][0]
    submesh_cell["full_commun"] = submesh["full_commun"][0]
    submesh_cell["full_quan"] ={}
    submesh_cell["ghost_quan"]={}
    for k in submesh["full_quan"]:
        submesh_cell["full_quan"][k] = submesh["full_quan"][k][0]
        submesh_cell["ghost_quan"][k] = submesh["ghost_quan"][k][0]

    numprocs = pypar.size()
    points, vertices, boundary, quantities, ghost_recv_dict, full_send_dict = \
            build_local_mesh(submesh_cell, 0, triangles_per_proc[0], numprocs)
    return  points, vertices, boundary, quantities, ghost_recv_dict, \
           full_send_dict




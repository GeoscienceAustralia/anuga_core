#########################################################
#
#
#  Read in a data file and subdivide the triangle list
#
#
#  The final routine, pmesh_divide_metis, does automatic
# grid partitioning. Once testing has finished on this
# routine the others should be removed.
#
#  Authors: Linda Stals and Matthew Hardy, June 2005
#  Modified: Linda Stals, Nov 2005
#            Jack Kelly, Nov 2005
#            Steve Roberts, Aug 2009 (updating to numpy)
#
#
#########################################################


import sys
from os import sep
from sys import path
from math import floor

import numpy as num

from anuga.abstract_2d_finite_volumes.neighbour_mesh import Mesh
from anuga import indent


#########################################################
#
# If the triangles list is reordered, the quantities
# assigned to the triangles must also be reorded.
#
# *) quantities contain the quantites in the old ordering
# *) proc_sum[i] contains the number of triangles in
# processor i
# *) tri_index is a map from the old triangle ordering to
# the new ordering, where the new number for triangle
# i is proc_sum[tri_index[i][0]]+tri_index[i][1]
#
# -------------------------------------------------------
#
# *) The quantaties are returned in the new ordering
#
#########################################################

def reorder(quantities, tri_index, proc_sum):

    # Find the number triangles

    N = len(tri_index)

    # Temporary storage area

    index = num.zeros(N, num.int)
    q_reord = {}

    # Find the new ordering of the triangles

    for i in range(N):
        bin = tri_index[i][0]
        bin_off_set = tri_index[i][1]
        index[i] = proc_sum[bin]+bin_off_set

    # Reorder each quantity according to the new ordering

    for k in quantities:
        q_reord[k] = num.zeros((N, 3), num.float)
        for i in range(N):
            q_reord[k][index[i]]=quantities[k].vertex_values[i]
    del index

    return q_reord


#########################################################
#
# Divide the mesh using a call to metis, through pymetis.
#
# -------------------------------------------------------
#
# *)  The nodes, triangles, boundary, and quantities are
# returned. triangles_per_proc defines the subdivision.
# The first triangles_per_proc[0] triangles are assigned
# to processor 0, the next triangles_per_proc[1] are
# assigned to processor 1 etc. The boundary and quantites
# are ordered the same way as the triangles
#
#########################################################

#path.append('..' + sep + 'pymetis')

try:
    from pymetis.metis import partMeshNodal
except ImportError:
    print "***************************************************"
    print "         Metis is probably not compiled."
    print "         Read \anuga_core\source\pymetis\README"
    print "***************************************************"
    raise ImportError

def pmesh_divide_metis(domain, n_procs):
    # Wrapper for old pmesh_divide_metis which does not return tri_index or r_tri_index
    nodes, ttriangles, boundary, triangles_per_proc, quantities, tri_index, r_tri_index = pmesh_divide_metis_helper(domain, n_procs)
    return nodes, ttriangles, boundary, triangles_per_proc, quantities

def pmesh_divide_metis_with_map(domain, n_procs):
    return pmesh_divide_metis_helper(domain, n_procs)

def pmesh_divide_metis_helper(domain, n_procs):
    
    # Initialise the lists
    # List, indexed by processor of # triangles.
    
    triangles_per_proc = []
    
    # List of lists, indexed by processor of vertex numbers
    
    tri_list = []

    # Serial to Parallel and Parallel to Serial Triangle index maps
    tri_index = {}
    r_tri_index = {} # reverse tri index, parallel to serial triangle index mapping
    
    # List indexed by processor of cumulative total of triangles allocated.
    
    proc_sum = []
    for i in range(n_procs):
        tri_list.append([])
        triangles_per_proc.append(0)
        proc_sum.append([])

    # Prepare variables for the metis call
    
    n_tri = len(domain.triangles)
    if n_procs != 1: #Because metis chokes on it...
        n_vert = domain.get_number_of_nodes()
        t_list = domain.triangles.copy()
        t_list = num.reshape(t_list, (-1,))
    
        # The 1 here is for triangular mesh elements.
        edgecut, epart, npart = partMeshNodal(n_tri, n_vert, t_list, 1, n_procs)
        # print edgecut
        # print npart
        # print epart
        del edgecut
        del npart

        # Sometimes (usu. on x86_64), partMeshNodal returns an array of zero
        # dimensional arrays. Correct this.
        if type(epart[0]) == num.ndarray:
            epart_new = num.zeros(len(epart), num.int)
            for i in range(len(epart)):
                epart_new[i] = epart[i][0]
            epart = epart_new
            del epart_new
        # Assign triangles to processes, according to what metis told us.
        
        # tri_index maps triangle number -> processor, new triangle number
        # (local to the processor)
        
        triangles = []        
        for i in range(n_tri):
            triangles_per_proc[epart[i]] = triangles_per_proc[epart[i]] + 1
            tri_list[epart[i]].append(domain.triangles[i])
            tri_index[i] = ([epart[i], len(tri_list[epart[i]]) - 1])
            r_tri_index[epart[i], len(tri_list[epart[i]]) - 1] = i
        
        # Order the triangle list so that all of the triangles belonging
        # to processor i are listed before those belonging to processor
        # i+1

        for i in range(n_procs):
            for t in tri_list[i]:
                triangles.append(t)
            
        # The boundary labels have to changed in accoradance with the
        # new triangle ordering, proc_sum and tri_index help with this

        proc_sum[0] = 0
        for i in range(n_procs - 1):
            proc_sum[i+1]=proc_sum[i]+triangles_per_proc[i]

        # Relabel the boundary elements to fit in with the new triangle
        # ordering

        boundary = {}
        for b in domain.boundary:
            t =  tri_index[b[0]]
            boundary[proc_sum[t[0]]+t[1], b[1]] = domain.boundary[b]

        quantities = reorder(domain.quantities, tri_index, proc_sum)
    else:
        boundary = domain.boundary.copy()
        triangles_per_proc[0] = n_tri
        triangles = domain.triangles.copy()
        
        # This is essentially the same as a chunk of code from reorder.
        
        quantities = {}
        for k in domain.quantities:
            quantities[k] = num.zeros((n_tri, 3), num.float)
            for i in range(n_tri):
                quantities[k][i] = domain.quantities[k].vertex_values[i]
        
    # Extract the node list
    
    nodes = domain.get_nodes().copy()
    
    # Convert the triangle datastructure to be an array type,
    # this helps with the communication

    ttriangles = num.zeros((len(triangles), 3), num.int)
    for i in range(len(triangles)):
        ttriangles[i] = triangles[i]

    #return nodes, ttriangles, boundary, triangles_per_proc, quantities
    
    return nodes, ttriangles, boundary, triangles_per_proc, quantities, tri_index, r_tri_index


#########################################################
#
# Subdivide the domain. This module is primarily
# responsible for building the ghost layer and
# communication pattern
#
#
#  Author: Linda Stals, June 2005
#  Modified: Linda Stals, Nov 2005 (optimise python code)
#            Steve Roberts, Aug 2009 (convert to numpy)
#
#
#########################################################



#########################################################
#
# Subdivide the triangles into non-overlapping domains.
#
#  *)  The subdivision is controlled by triangles_per_proc.
# The first triangles_per_proc[0] triangles are assigned
# to the first processor, the second triangles_per_proc[1]
# are assigned to the second processor etc.
#
#  *) nodes, triangles and boundary contains all of the
# nodes, triangles and boundary tag information for the
# whole domain. The triangles should be orientated in the
# correct way and the nodes number consecutively from 0.
#
# -------------------------------------------------------
#
#  *) A dictionary containing the full_nodes, full_triangles
# and full_boundary information for each processor is
# returned. The node information consists of
# [global_id, x_coord, y_coord].
#
#########################################################

def submesh_full(nodes, triangles, boundary, triangles_per_proc):

    # Initialise

    tlower = 0
    nproc = len(triangles_per_proc)
    nnodes = len(nodes)
    node_list = []
    triangle_list = []
    boundary_list = []
    submesh = {}
    node_range = num.reshape(num.arange(nnodes),(nnodes,1))

    #print node_range
    tsubnodes = num.concatenate((node_range, nodes), 1)


    # Loop over processors

    for p in range(nproc):

        # Find triangles on processor p

        tupper = triangles_per_proc[p]+tlower
        subtriangles = triangles[tlower:tupper]
        triangle_list.append(subtriangles)

        # Find the boundary edges on processor p

        subboundary = {}
        for k in boundary:
            if (k[0] >=tlower and k[0] < tupper):
                subboundary[k]=boundary[k]
        boundary_list.append(subboundary)

        # Find nodes in processor p

        nodemap = num.zeros(nnodes, 'i')
        for t in subtriangles:
            nodemap[t[0]]=1
            nodemap[t[1]]=1
            nodemap[t[2]]=1

        
        node_list.append(tsubnodes.take(num.flatnonzero(nodemap),axis=0))

        # Move to the next processor

        tlower = tupper

    # Put the results in a dictionary

    submesh["full_nodes"] = node_list
    submesh["full_triangles"] = triangle_list
    submesh["full_boundary"] = boundary_list

    # Clean up before exiting

    del (nodemap)

    return submesh


#########################################################
#
# Build the ghost layer of triangles
#
#  *) Given the triangle subpartion for the processor
# build a ghost layer of triangles. The ghost layer
# consists of two layers of neighbouring triangles.
#
#  *) The vertices in the ghost triangles must also
# be added to the node list for the current processor
#
#
# -------------------------------------------------------
#
#  *) The extra triangles and nodes are returned.
#
#  *)  The node information consists of
# [global_id, x_coord, y_coord].
#
#  *) The triangle information consists of
# [triangle number, t], where t = [v1, v2, v3].
#
#########################################################

def ghost_layer(submesh, mesh, p, tupper, tlower):

    ncoord = mesh.number_of_nodes
    ntriangles = mesh.number_of_triangles

    # Find the first layer of boundary triangles

    trianglemap = num.zeros(ntriangles, 'i')
    for t in range(tlower, tupper):
        
        n = mesh.neighbours[t, 0]

        if n >= 0:
            if n < tlower or n >= tupper:
                trianglemap[n] = 1
        n = mesh.neighbours[t, 1]
        if n >= 0:
            if n < tlower or n >= tupper:
                trianglemap[n] = 1
        n = mesh.neighbours[t, 2]
        if n >= 0:
            if n < tlower or n >= tupper:
                trianglemap[n] = 1

    # Find the second layer of boundary triangles

    for t in range(len(trianglemap)):
        if trianglemap[t]==1:
            n = mesh.neighbours[t, 0]
            if n >= 0:
                if (n < tlower or n >= tupper) and trianglemap[n] == 0:
                    trianglemap[n] = 2
            n = mesh.neighbours[t, 1]
            if n >= 0:
                if (n < tlower or n >= tupper) and trianglemap[n] == 0:
                    trianglemap[n] = 2
            n = mesh.neighbours[t, 2]
            if n >= 0:
                if (n < tlower or n >= tupper) and trianglemap[n] == 0:
                    trianglemap[n] = 2

    # Build the triangle list and make note of the vertices

    nodemap = num.zeros(ncoord, 'i')
    fullnodes = submesh["full_nodes"][p]

    subtriangles = []
    for i in range(len(trianglemap)):
        if trianglemap[i] != 0:
            t = list(mesh.triangles[i])
            nodemap[t[0]] = 1
            nodemap[t[1]] = 1
            nodemap[t[2]] = 1

    trilist = num.reshape(num.arange(ntriangles),(ntriangles,1))
    tsubtriangles = num.concatenate((trilist, mesh.triangles), 1)
    subtriangles = tsubtriangles.take(num.flatnonzero(trianglemap),axis=0)

    
    # Keep a record of the triangle vertices, if they are not already there

    subnodes = []
    for n in fullnodes:
        nodemap[int(n[0])] = 0

    nodelist = num.reshape(num.arange(ncoord),(ncoord,1))
    tsubnodes = num.concatenate((nodelist, mesh.get_nodes()), 1)
    subnodes = tsubnodes.take(num.flatnonzero(nodemap),axis=0)

    # Clean up before exiting

    del (nodelist)
    del (trilist)
    del (tsubnodes)
    del (nodemap)
    del (trianglemap)

    # Return the triangles and vertices sitting on the boundary layer

    return subnodes, subtriangles

#########################################################
#
# Find the edges of the ghost trianlges that do not
# have a neighbour in the current cell. These are
# treated as a special type of boundary edge.
#
#  *) Given the ghost triangles in a particular
# triangle, use the mesh to find its neigbours. If
# the neighbour is not in the processor set it to
# be a boundary edge
#
#  *) The vertices in the ghost triangles must also
# be added to the node list for the current processor
#
#  *) The boundary edges for the ghost triangles are
# ignored.
#
# -------------------------------------------------------
#
#  *) The type assigned to the ghost boundary edges is 'ghost'
#
#  *)  The boundary information is returned as a directorier
# with the key = (triangle id, edge no) and the values
# assigned to the key is 'ghost'
#
#
#########################################################
def is_in_processor(ghost_list, tlower, tupper, n):

    return num.equal(ghost_list,n).any() or (tlower <= n and tupper > n)


def ghost_bnd_layer(ghosttri, tlower, tupper, mesh, p):

    ghost_list = []
    subboundary = {}


    for t in ghosttri:
        ghost_list.append(t[0])
    
    for t in ghosttri:

        n = mesh.neighbours[t[0], 0]
        if not is_in_processor(ghost_list, tlower, tupper, n):
            subboundary[t[0], 0] = 'ghost'

        n = mesh.neighbours[t[0], 1]
        if not is_in_processor(ghost_list, tlower, tupper, n):
            subboundary[t[0], 1] = 'ghost'

        n = mesh.neighbours[t[0], 2]
        if not is_in_processor(ghost_list, tlower, tupper, n):
            subboundary[t[0], 2] = 'ghost'
            
    return subboundary

#########################################################
#
# The ghost triangles on the current processor will need
# to get updated information from the neighbouring
# processor containing the corresponding full triangles.
#
#  *) The tri_per_proc is used to determine which
# processor contains the full node copy.
#
# -------------------------------------------------------
#
#  *) The ghost communication pattern consists of
# [global node number, neighbour processor number].
#
#########################################################

def ghost_commun_pattern(subtri, p, tri_per_proc):

    # Loop over the ghost triangles

    ghost_commun = num.zeros((len(subtri), 2), num.int)

    for i in range(len(subtri)):
        global_no = subtri[i][0]

        # Find which processor contains the full triangle

        nproc = len(tri_per_proc)
        neigh = nproc-1
        sum = 0
        for q in range(nproc-1):
            if (global_no < sum+tri_per_proc[q]):
                neigh = q
                break
            sum = sum+tri_per_proc[q]

        # Keep a copy of the neighbour processor number

        ghost_commun[i] = [global_no, neigh]

    return ghost_commun

#########################################################
#
# The full triangles in this processor must communicate
# updated information to neighbouring processor that
# contain ghost triangles
#
#  *) The ghost communication pattern for all of the
# processor must be built before calling this processor.
#
#  *) The full communication pattern is found by looping
# through the ghost communication pattern for all of the
# processors. Recall that this information is stored in
# the form [global node number, neighbour processor number].
# The full communication for the neighbour processor is
# then updated.
#
# -------------------------------------------------------
#
#  *) The full communication pattern consists of
# [global id, [p1, p2, ...]], where p1, p2 etc contain
# a ghost node copy of the triangle global id.
#
#########################################################

def full_commun_pattern(submesh, tri_per_proc):
    tlower = 0
    nproc = len(tri_per_proc)
    full_commun = []

    # Loop over the processor

    for p in range(nproc):

        # Loop over the full triangles in the current processor
        # and build an empty dictionary

        fcommun = {}
        tupper = tri_per_proc[p]+tlower
        for i in range(tlower, tupper):
            fcommun[i] = []
        full_commun.append(fcommun)
        tlower = tupper

    # Loop over the processor again

    for p in range(nproc):

        # Loop over the ghost triangles in the current processor,
        # find which processor contains the corresponding full copy
        # and note that the processor must send updates to this
        # processor

        for g in submesh["ghost_commun"][p]:
            neigh = g[1]
            full_commun[neigh][g[0]].append(p)

    return full_commun


#########################################################
#
# Given the non-overlapping grid partition, an extra layer
# of triangles are included to help with the computations.
# The triangles in this extra layer are not updated by
# the processor, their updated values must be sent by the
# processor containing the original, full, copy of the
# triangle. The communication pattern that controls these
# updates must also be built.
#
#  *) Assumes that full triangles, nodes etc have already
# been found and stored in submesh
#
#  *) See the documentation for ghost_layer,
# ghost_commun_pattern and full_commun_pattern
#
# -------------------------------------------------------
#
#  *) The additional information is added to the submesh
# dictionary. See the documentation for ghost_layer,
# ghost_commun_pattern and full_commun_pattern
#
#  *) The ghost_triangles, ghost_nodes, ghost_boundary,
# ghost_commun and full_commun is added to submesh
#########################################################

def submesh_ghost(submesh, mesh, triangles_per_proc):

    nproc = len(triangles_per_proc)
    tlower = 0
    ghost_triangles = []
    ghost_nodes = []
    ghost_commun = []
    ghost_bnd = []

    # Loop over the processors

    for p in range(nproc):

        # Find the full triangles in this processor

        tupper = triangles_per_proc[p]+tlower

        # Build the ghost boundary layer

        [subnodes, subtri] = \
                   ghost_layer(submesh, mesh, p, tupper, tlower)
        ghost_triangles.append(subtri)
        ghost_nodes.append(subnodes)


        # Find the boundary layer formed by the ghost triangles
        
        subbnd = ghost_bnd_layer(subtri, tlower, tupper, mesh, p)
        ghost_bnd.append(subbnd)
        
        # Build the communication pattern for the ghost nodes

        gcommun = \
                ghost_commun_pattern(subtri, p, triangles_per_proc)
        ghost_commun.append(gcommun)

        # Move to the next processor

        tlower = tupper


    # Record the ghost layer and communication pattern

    submesh["ghost_nodes"] = ghost_nodes
    submesh["ghost_triangles"] = ghost_triangles
    submesh["ghost_commun"] = ghost_commun
    submesh["ghost_boundary"] = ghost_bnd
    
    # Build the communication pattern for the full triangles

    full_commun = full_commun_pattern(submesh, triangles_per_proc)
    submesh["full_commun"] = full_commun

    # Return the submesh

    return submesh


#########################################################
#
# Certain quantities may be assigned to the triangles,
# these quantities must be subdivided in the same way
# as the triangles
#
#  *) The quantities are ordered in the same way as the
# triangles
#
# -------------------------------------------------------
#
#  *) The quantites attached to the full triangles are
# stored in full_quan
#
#  *) The quantities attached to the ghost triangles are
# stored in ghost_quan
#########################################################

def submesh_quantities(submesh, quantities, triangles_per_proc):

    nproc = len(triangles_per_proc)

    lower = 0

    # Build an empty dictionary to hold the quantites

    submesh["full_quan"] = {}
    submesh["ghost_quan"] = {}
    for k in quantities:
        submesh["full_quan"][k] = []
        submesh["ghost_quan"][k] = []

    # Loop trough the subdomains

    for p in range(nproc):
        upper =   lower+triangles_per_proc[p]

        # Find the global ID of the ghost triangles

        global_id = []
        M = len(submesh["ghost_triangles"][p])
        for j in range(M):
            global_id.append(submesh["ghost_triangles"][p][j][0])

        # Use the global ID to extract the quantites information from
        # the full domain

        for k in quantities:
            submesh["full_quan"][k].append(quantities[k][lower:upper])
            submesh["ghost_quan"][k].append(num.zeros( (M,3) , num.float))
            for j in range(M):
                submesh["ghost_quan"][k][p][j] = \
                                               quantities[k][global_id[j]]

        lower = upper

    return submesh

#########################################################
#
# Build the grid partition on the host.
#
#  *) See the documentation for submesh_ghost and
# submesh_full
#
# -------------------------------------------------------
#
#  *) A dictionary containing the full_triangles,
# full_nodes, full_boundary, ghost_triangles, ghost_nodes,
# ghost_boundary, ghost_commun and full_commun and true boundary polygon is returned.
#
#########################################################

def build_submesh(nodes, triangles, edges, quantities,
                  triangles_per_proc):

    # Temporarily build the mesh to find the neighbouring
    # triangles and true boundary polygon

    mesh = Mesh(nodes, triangles)
    boundary_polygon = mesh.get_boundary_polygon()
    

    # Subdivide into non-overlapping partitions

    submeshf = submesh_full(nodes, triangles, edges, \
                            triangles_per_proc)
    
    # Add any extra ghost boundary layer information

    submeshg = submesh_ghost(submeshf, mesh, triangles_per_proc)

    # Order the quantities information to be the same as the triangle
    # information

    submesh = submesh_quantities(submeshg, quantities, \
                                 triangles_per_proc)

    submesh["boundary_polygon"] = boundary_polygon
    return submesh

#########################################################
#
#  Given the subdivision of the grid assigned to the
# current processor convert it into a form that is
# appropriate for the GA datastructure.
#
#  The main function of these modules is to change the
# node numbering. The GA datastructure assumes they
# are numbered consecutively from 0.
#
#  The module also changes the communication pattern
# datastructure into a form needed by parallel_advection
#
#  Authors: Linda Stals and Matthew Hardy, June 2005
#  Modified: Linda Stals, Nov 2005 (optimise python code)
#            Steve Roberts, Aug 2009 (updating to numpy)
#
#
#########################################################

 
#########################################################
# Convert the format of the data to that used by ANUGA
#
#
# *) Change the nodes global ID's to an integer value,
#starting from 0.
#
# *) The triangles and boundary edges must also be
# updated accordingly.
#
# -------------------------------------------------------
#
# *) The nodes, triangles and boundary edges defined by
# the new numbering scheme are returned
#
#########################################################

def build_local_GA(nodes, triangles, boundaries, tri_map):

    Nnodes =len(nodes)
    Ntriangles = len(triangles)
    
    # Extract the nodes (using the local ID)
    
    GAnodes = num.take(nodes, (1, 2), 1)

    # Build a global ID to local ID mapping

    NGlobal = 0
    for i in range(Nnodes):
        if nodes[i][0] > NGlobal:
            NGlobal = nodes[i][0]

    node_map = -1*num.ones(int(NGlobal)+1, num.int)

    num.put(node_map, num.take(nodes, (0,), 1).astype(num.int), \
        num.arange(Nnodes))
        
    # Change the global IDs in the triangles to the local IDs

    GAtriangles = num.zeros((Ntriangles, 3), num.int)
    GAtriangles[:,0] = num.take(node_map, triangles[:,0])
    GAtriangles[:,1] = num.take(node_map, triangles[:,1])
    GAtriangles[:,2] = num.take(node_map, triangles[:,2])

    # Change the triangle numbering in the boundaries

    GAboundaries = {}
    for b in boundaries:
        GAboundaries[tri_map[b[0]], b[1]] = boundaries[b]
        
    
    return GAnodes, GAtriangles, GAboundaries, node_map


#########################################################
# Change the communication format to that needed by the
# parallel advection file.
#
# *) The index contains [global triangle no,
# local triangle no.]
#
# -------------------------------------------------------
#
# *) The ghost_recv and full_send dictionaries are
# returned.
#
# *) ghost_recv dictionary is local id, global id, value
#
# *) full_recv dictionary is local id, global id, value
#
# *) The information is ordered by the global id. This
# means that the communication order is predetermined and
# local and global id do not need to be
# compared when the information is sent/received.
#
#########################################################

def build_local_commun(tri_map, ghostc, fullc, nproc):

    # Initialise

    full_send = {}
    ghost_recv = {}

    # Build the ghost_recv dictionary (sort the
    # information by the global numbering)
    
    ghostc = num.sort(ghostc, 0)
    
    for c in range(nproc):
        s = ghostc[:,0]
        d = num.compress(num.equal(ghostc[:,1],c), s)
        if len(d) > 0:
            ghost_recv[c] = [0, 0]
            ghost_recv[c][1] = d
            ghost_recv[c][0] = num.take(tri_map, d)
            
    # Build a temporary copy of the full_send dictionary
    # (this version allows the information to be stored
    # by the global numbering)

    tmp_send = {}
    for global_id in fullc:
        for i in range(len(fullc[global_id])):
            neigh = fullc[global_id][i]
            if not tmp_send.has_key(neigh):
                tmp_send[neigh] = []
            tmp_send[neigh].append([global_id, \
                                    tri_map[global_id]])

    # Extract the full send information and put it in the form
    # required for the full_send dictionary

    for neigh in tmp_send:
        neigh_commun = num.sort(tmp_send[neigh], 0)
        full_send[neigh] = [0, 0]
        full_send[neigh][0] = neigh_commun[:,1]
        full_send[neigh][1] = neigh_commun[:,0]

    return ghost_recv, full_send


#########################################################
# Convert the format of the data to that used by ANUGA
#
#
# *) Change the nodes global ID's to an integer value,
# starting from 0. The node numbering in the triangles
# must also be updated to take this into account.
#
# *) The triangle number will also change, which affects
# the boundary tag information and the communication
# pattern.
#
# -------------------------------------------------------
#
# *) The nodes, triangles, boundary edges and communication
# pattern defined by the new numbering scheme are returned
#
#########################################################

def build_local_mesh(submesh, lower_t, upper_t, nproc):

    # Combine the full nodes and ghost nodes

    nodes = num.concatenate((submesh["full_nodes"], \
                         submesh["ghost_nodes"]))
    
    # Combine the full triangles and ghost triangles

    gtri =  num.take(submesh["ghost_triangles"],(1, 2, 3),1)
    triangles = num.concatenate((submesh["full_triangles"], gtri))

    # Combine the full boundaries and ghost boundaries

    boundaries = submesh["full_boundary"]
    for b in submesh["ghost_boundary"]:
        boundaries[b]=submesh["ghost_boundary"][b]

    # Make note of the new triangle numbers, including the ghost
    # triangles

    NGlobal = upper_t
    for i in range(len(submesh["ghost_triangles"])):
        id = submesh["ghost_triangles"][i][0]
        if id > NGlobal:
            NGlobal = id
    #index = num.zeros(int(NGlobal)+1, num.int)
    tri_map = -1*num.ones(int(NGlobal)+1, num.int)
    tri_map[lower_t:upper_t]=num.arange(upper_t-lower_t)
    for i in range(len(submesh["ghost_triangles"])):
        tri_map[submesh["ghost_triangles"][i][0]] = i+upper_t-lower_t
    
    # Change the node numbering (and update the numbering in the
    # triangles)

    [GAnodes, GAtriangles, GAboundary, node_map] = \
    build_local_GA(nodes, triangles, boundaries, tri_map)

    # Extract the local quantities
    
    quantities ={}
    for k in submesh["full_quan"]:
        Nf = len(submesh["full_quan"][k])
        Ng = len(submesh["ghost_quan"][k])
        quantities[k] = num.zeros((Nf+Ng, 3), num.float)
        quantities[k][0:Nf] = submesh["full_quan"][k] 
        quantities[k][Nf:Nf+Ng] = submesh["ghost_quan"][k]
                             
    # Change the communication pattern into a form needed by
    # the parallel_adv

    gcommun = submesh["ghost_commun"]
    fcommun = submesh["full_commun"]
    [ghost_rec, full_send] = \
                build_local_commun(tri_map, gcommun, fcommun, nproc)


    return GAnodes, GAtriangles, GAboundary, quantities, ghost_rec, \
           full_send, tri_map, node_map


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
#            Steve Roberts, Aug 2009 (update to numpy)
#
#
#########################################################


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

def send_submesh(submesh, triangles_per_proc, p, verbose=True):

    import pypar
    
    myid = pypar.rank()
    
    if verbose: print 'process %d sending submesh to process %d' %(myid, p)
    
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

    pypar.send(triangles_per_proc, p)

    # compress full_commun

    flat_full_commun = []

    for c in submesh["full_commun"][p]:
        for i in range(len(submesh["full_commun"][p][c])):
            flat_full_commun.append([c,submesh["full_commun"][p][c][i]])

    # send the array sizes so memory can be allocated

    setup_array = num.zeros((9,),num.int)
    setup_array[0] = len(submesh["full_nodes"][p])
    setup_array[1] = len(submesh["ghost_nodes"][p])
    setup_array[2] = len(submesh["full_triangles"][p])
    setup_array[3] = len(submesh["ghost_triangles"][p])
    setup_array[4] = len(submesh["full_boundary"][p])
    setup_array[5] = len(submesh["ghost_boundary"][p])
    setup_array[6] = len(submesh["ghost_commun"][p])
    setup_array[7] = len(flat_full_commun)
    setup_array[8] = len(submesh["full_quan"])
    
    pypar.send(num.array(setup_array, num.int), p)
    
    # send the nodes

    pypar.send(num.array(submesh["full_nodes"][p], num.float), p)
    pypar.send(num.array(submesh["ghost_nodes"][p], num.float),p)

    # send the triangles

    pypar.send(num.array(submesh["full_triangles"][p],  num.int), p)
    pypar.send(num.array(submesh["ghost_triangles"][p], num.int), p)

    # send the boundary

    bc = []
    for b in submesh["full_boundary"][p]:
        bc.append([b[0], b[1], tagmap[submesh["full_boundary"][p][b]]])


    pypar.send(num.array(bc, num.int), p)

    bc = []
    for b in submesh["ghost_boundary"][p]:
        bc.append([b[0], b[1], tagmap[submesh["ghost_boundary"][p][b]]])

    pypar.send(num.array(bc, num.int), p)

    # send the communication pattern

    pypar.send(submesh["ghost_commun"][p], p)

    pypar.send(num.array(flat_full_commun, num.int), p)

    # send the quantities
    
    for k in submesh["full_quan"]:
        pypar.send(num.array(submesh["full_quan"][k][p], num.float), p)
        
    for k in submesh["ghost_quan"]:
        pypar.send(num.array(submesh["ghost_quan"][k][p], num.float),p)
        

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

def rec_submesh_flat(p, verbose=True):

    import pypar
    
    numproc = pypar.size()
    myid = pypar.rank()

    submesh_cell = {}
    
    if verbose: print indent+'P%d: Receiving submesh from P%d' %(myid, p)

    # receive the tagmap for the boundary conditions
    
    tagmap = pypar.receive(p)

    itagmap = {}
    for t in tagmap:
        itagmap[tagmap[t]]=t

    # receive the quantities key information

    qkeys = pypar.receive(p)

    # receive the number of triangles per processor

    triangles_per_proc = pypar.receive(p)

    # recieve information about the array sizes

    setup_array = pypar.receive(p)

    no_full_nodes      = setup_array[0]
    no_ghost_nodes     = setup_array[1]
    no_full_triangles  = setup_array[2]
    no_ghost_triangles = setup_array[3]
    no_full_boundary   = setup_array[4]
    no_ghost_boundary  = setup_array[5]
    no_ghost_commun    = setup_array[6]
    no_full_commun     = setup_array[7]
    no_quantities      = setup_array[8]
    
    # receive the full nodes

    submesh_cell["full_nodes"] = pypar.receive(p)

    # receive the ghost nodes

    submesh_cell["ghost_nodes"] = pypar.receive(p)
    
    # receive the full triangles

    submesh_cell["full_triangles"] = pypar.receive(p)
    
    # receive the ghost triangles

    submesh_cell["ghost_triangles"] = pypar.receive(p)

    # receive the full boundary

    bnd_c = pypar.receive(p)

    submesh_cell["full_boundary"] = {}
    for b in bnd_c:
        submesh_cell["full_boundary"][b[0],b[1]]=itagmap[b[2]]

    # receive the ghost boundary

    bnd_c = pypar.receive(p)

    submesh_cell["ghost_boundary"] = {}
    for b in bnd_c:
        submesh_cell["ghost_boundary"][b[0],b[1]]=itagmap[b[2]]

    # receive the ghost communication pattern

    submesh_cell["ghost_commun"] = pypar.receive(p)
    
    # receive the full communication pattern

    full_commun = pypar.receive(p)

    submesh_cell["full_commun"] = {}
    for c in full_commun:
        submesh_cell["full_commun"][c[0]] = []
    for c in full_commun:
        submesh_cell["full_commun"][c[0]].append(c[1])

    # receive the quantities

    submesh_cell["full_quan"]={}
    
    for i in range(no_quantities):
        tmp = pypar.receive(p)
        submesh_cell["full_quan"][qkeys[i]]=num.zeros((no_full_triangles,3), num.float)
        submesh_cell["full_quan"][qkeys[i]][:] = tmp[:]

    submesh_cell["ghost_quan"]={}
    for i in range(no_quantities):
        tmp = pypar.receive(p)
        submesh_cell["ghost_quan"][qkeys[i]]= num.zeros((no_ghost_triangles,3), num.float)
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

def rec_submesh(p, verbose=True):

    import pypar
    
    numproc = pypar.size()
    myid = pypar.rank()

    [submesh_cell, triangles_per_proc,\
     number_of_full_nodes, number_of_full_triangles] = rec_submesh_flat(p,verbose)
    
    # find the full triangles assigned to this processor

    lower_t = 0
    for i in range(myid):
        lower_t = lower_t+triangles_per_proc[i]
    upper_t = lower_t+triangles_per_proc[myid]

    # convert the information into a form needed by the GA
    # datastructure

    [GAnodes, GAtriangles, boundary, quantities, \
     ghost_rec, full_send, tri_map, node_map] = \
              build_local_mesh(submesh_cell, lower_t, upper_t, \
                               numproc)
    
    return GAnodes, GAtriangles, boundary, quantities,\
           ghost_rec, full_send,\
           number_of_full_nodes, number_of_full_triangles, tri_map, node_map
          


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

    numprocs = len(triangles_per_proc)
    points, vertices, boundary, quantities, ghost_recv_dict, \
            full_send_dict, tri_map, node_map = \
            build_local_mesh(submesh_cell, 0, triangles_per_proc[0], numprocs)


    return  points, vertices, boundary, quantities, ghost_recv_dict, \
           full_send_dict, tri_map, node_map
           





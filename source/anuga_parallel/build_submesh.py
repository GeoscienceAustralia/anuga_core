#########################################################
#
# Subdivide the GA domain. This module is primarily
# responsible for building the ghost layer and
# communication pattern
#
#
#  Author: Linda Stals, June 2005
#  Modified: Linda Stals, Nov 2005 (optimise python code)
#
#
#########################################################

import sys

from Numeric import zeros, Float, Int, concatenate, \
     reshape, arrayrange, take, nonzero

from anuga.abstract_2d_finite_volumes.neighbour_mesh import Mesh



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
    node_range = reshape(arrayrange(nnodes),(nnodes,1))
    tsubnodes = concatenate((node_range, nodes), 1)

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

        nodemap = zeros(nnodes, 'i')
        for t in subtriangles:
            nodemap[t[0]]=1
            nodemap[t[1]]=1
            nodemap[t[2]]=1

        node_list.append(take(tsubnodes,nonzero(nodemap)))

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

    trianglemap = zeros(ntriangles, 'i')
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

    nodemap = zeros(ncoord, 'i')
    fullnodes = submesh["full_nodes"][p]

    subtriangles = []
    for i in range(len(trianglemap)):
        if trianglemap[i] != 0:
            t = list(mesh.triangles[i])
            nodemap[t[0]] = 1
            nodemap[t[1]] = 1
            nodemap[t[2]] = 1

    trilist = reshape(arrayrange(ntriangles),(ntriangles,1))
    tsubtriangles = concatenate((trilist, mesh.triangles), 1)
    subtriangles = take(tsubtriangles, nonzero(trianglemap))

    # Keep a record of the triangle vertices, if they are not already there

    subnodes = []
    for n in fullnodes:
        nodemap[int(n[0])] = 0

    nodelist = reshape(arrayrange(ncoord),(ncoord,1))
    tsubnodes = concatenate((nodelist, mesh.get_nodes()), 1)
    subnodes = take(tsubnodes, nonzero(nodemap))

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
    return (n in ghost_list) or (tlower <= n and tupper > n)

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

    ghost_commun = zeros((len(subtri), 2), Int)

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
            submesh["ghost_quan"][k].append(zeros( (M,3) , Float))
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


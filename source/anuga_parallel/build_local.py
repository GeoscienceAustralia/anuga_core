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
#
#
#########################################################

from Numeric import  zeros, Float, Int, concatenate, \
     take, arrayrange, put, sort, compress, equal


#########################################################
# Convert the format of the data to that used by
# pyvolution
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

def build_local_GA(nodes, triangles, boundaries, tri_index):

    Nnodes =len(nodes)
    Ntriangles = len(triangles)
    
    # Extract the nodes (using the local ID)
    
    GAnodes = take(nodes, (1, 2), 1)

    # Build a global ID to local ID mapping

    NGlobal = 0
    for i in range(Nnodes):
        if nodes[i][0] > NGlobal:
            NGlobal = nodes[i][0]
    index = zeros(int(NGlobal)+1, Int)
    put(index, take(nodes, (0,), 1).astype(Int), \
        arrayrange(Nnodes))
        
    # Change the global IDs in the triangles to the local IDs

    GAtriangles = zeros((Ntriangles, 3), Int)
    GAtriangles[:,0] = take(index, triangles[:,0])
    GAtriangles[:,1] = take(index, triangles[:,1])
    GAtriangles[:,2] = take(index, triangles[:,2])

    # Change the triangle numbering in the boundaries

    GAboundaries = {}
    for b in boundaries:
        GAboundaries[tri_index[b[0]], b[1]] = boundaries[b]
        
    del (index)
    
    return GAnodes, GAtriangles, GAboundaries


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

def build_local_commun(index, ghostc, fullc, nproc):

    # Initialise

    full_send = {}
    ghost_recv = {}

    # Build the ghost_recv dictionary (sort the
    # information by the global numbering)
    
    ghostc = sort(ghostc, 0)
    
    for c in range(nproc):
        s = ghostc[:,0]
        d = compress(equal(ghostc[:,1],c), s)
        if len(d) > 0:
            ghost_recv[c] = [0, 0]
            ghost_recv[c][1] = d
            ghost_recv[c][0] = take(index, d)
            
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
                                    index[global_id]])

    # Extract the full send information and put it in the form
    # required for the full_send dictionary

    for neigh in tmp_send:
        neigh_commun = sort(tmp_send[neigh], 0)
        full_send[neigh] = [0, 0]
        full_send[neigh][0] = neigh_commun[:,1]
        full_send[neigh][1] = neigh_commun[:,0]

    return ghost_recv, full_send


#########################################################
# Convert the format of the data to that used by
# pyvolution
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

    nodes = concatenate((submesh["full_nodes"], \
                         submesh["ghost_nodes"]))
    
    # Combine the full triangles and ghost triangles

    gtri =  take(submesh["ghost_triangles"],(1, 2, 3),1)
    triangles = concatenate((submesh["full_triangles"], gtri))

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
    index = zeros(int(NGlobal)+1, Int)
    index[lower_t:upper_t]=arrayrange(upper_t-lower_t)
    for i in range(len(submesh["ghost_triangles"])):
        index[submesh["ghost_triangles"][i][0]] = i+upper_t-lower_t
    
    # Change the node numbering (and update the numbering in the
    # triangles)

    [GAnodes, GAtriangles, GAboundary] = build_local_GA(nodes, triangles, boundaries, index)

    # Extract the local quantities
    
    quantities ={}
    for k in submesh["full_quan"]:
        Nf = len(submesh["full_quan"][k])
        Ng = len(submesh["ghost_quan"][k])
        quantities[k] = zeros((Nf+Ng, 3), Float)
        quantities[k][0:Nf] = submesh["full_quan"][k] 
        quantities[k][Nf:Nf+Ng] = submesh["ghost_quan"][k]
                             
    # Change the communication pattern into a form needed by
    # the parallel_adv

    gcommun = submesh["ghost_commun"]
    fcommun = submesh["full_commun"]
    [ghost_rec, full_send] = \
                build_local_commun(index, gcommun, fcommun, nproc)

    # Clean up before exiting

    del(index)

    return GAnodes, GAtriangles, GAboundary, quantities, ghost_rec, \
           full_send

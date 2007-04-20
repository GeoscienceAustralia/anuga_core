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
#
#
#########################################################

from os import sep
from sys import path
from math import floor

from Numeric import zeros, Float, Int, reshape, argsort, ArrayType


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

    index = zeros(N, Int)
    q_reord = {}

    # Find the new ordering of the triangles

    for i in range(N):
        bin = tri_index[i][0]
        bin_off_set = tri_index[i][1]
        index[i] = proc_sum[bin]+bin_off_set

    # Reorder each quantity according to the new ordering

    for k in quantities:
        q_reord[k] = zeros((N, 3), Float)
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
    
    # Initialise the lists
    # List, indexed by processor of # triangles.
    
    triangles_per_proc = []
    
    # List of lists, indexed by processor of vertex numbers
    
    tri_list = []
    
    # List indexed by processor of cumulative total of triangles allocated.
    
    proc_sum = []
    for i in range(n_procs):
        tri_list.append([])
        triangles_per_proc.append(0)
        proc_sum.append([])

    # Prepare variables for the metis call
    
    n_tri = len(domain.triangles)
    if n_procs != 1: #Because metis chokes on it...
        n_vert = domain.number_of_nodes
        t_list = domain.triangles.copy()
        t_list = reshape(t_list, (-1,))
    
        # The 1 here is for triangular mesh elements.
        edgecut, epart, npart = partMeshNodal(n_tri, n_vert, t_list, 1, n_procs)
        # print edgecut
        # print npart
        # print epart
        del edgecut
        del npart

        # Sometimes (usu. on x86_64), partMeshNodal returnes an array of zero
        # dimensional arrays. Correct this.
        if type(epart[0]) == ArrayType:
            epart_new = zeros(len(epart), Int)
            for i in range(len(epart)):
                epart_new[i] = epart[i][0]
            epart = epart_new
            del epart_new
        # Assign triangles to processes, according to what metis told us.
        
        # tri_index maps triangle number -> processor, new triangle number
        # (local to the processor)
        
        tri_index = {}
        triangles = []        
        for i in range(n_tri):
            triangles_per_proc[epart[i]] = triangles_per_proc[epart[i]] + 1
            tri_list[epart[i]].append(domain.triangles[i])
            tri_index[i] = ([epart[i], len(tri_list[epart[i]]) - 1])
        
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
            quantities[k] = zeros((n_tri, 3), Float)
            for i in range(n_tri):
                quantities[k][i] = domain.quantities[k].vertex_values[i]
        
    # Extract the node list
    
    nodes = domain.get_nodes().copy()
    
    # Convert the triangle datastructure to be an array type,
    # this helps with the communication

    ttriangles = zeros((len(triangles), 3), Int)
    for i in range(len(triangles)):
        ttriangles[i] = triangles[i]
    
    return nodes, ttriangles, boundary, triangles_per_proc, quantities

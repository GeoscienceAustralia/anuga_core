"""Trying to lump parallel stuff into simpler interface


"""

# Parallelism
import pypar   # The Python-MPI interface
from anuga_parallel.pmesh_divide  import pmesh_divide_metis
from anuga_parallel.build_submesh import build_submesh
from anuga_parallel.build_local   import build_local_mesh
from anuga_parallel.build_commun  import send_submesh, rec_submesh, extract_hostmesh
from anuga_parallel.parallel_shallow_water import Parallel_Domain


#------------------------------------------------------------------------------
# Read in processor information
#------------------------------------------------------------------------------

numprocs = pypar.size()
myid = pypar.rank()
processor_name = pypar.Get_processor_name()
print 'I am processor %d of %d on node %s' %(myid, numprocs, processor_name)




def distribute_mesh(domain):

    numprocs = pypar.size()

    
    # Subdivide the mesh
    print 'Subdivide mesh'
    nodes, triangles, boundary, triangles_per_proc, quantities = \
           pmesh_divide_metis(domain, numprocs)

    # Build the mesh that should be assigned to each processor,
    # this includes ghost nodes and the communicaiton pattern
    print 'Build submeshes'    
    submesh = build_submesh(nodes, triangles, boundary,\
                            quantities, triangles_per_proc)

    # Send the mesh partition to the appropriate processor
    print 'Distribute submeshes'        
    for p in range(1, numprocs):
      send_submesh(submesh, triangles_per_proc, p)

    # Build the local mesh for processor 0
    points, vertices, boundary, quantities, ghost_recv_dict, full_send_dict = \
              extract_hostmesh(submesh, triangles_per_proc)

    # Return stuff
    return points, vertices, boundary, quantities, \
           ghost_recv_dict, full_send_dict
    




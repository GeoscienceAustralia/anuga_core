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




def distribute(domain, verbose=False):
    """ Distribute the domain to all processes
    """

    # For some obscure reason this communication must happen prior to
    # the more complex mesh distribution - Oh Well!
    if myid == 0:
        domain_name = domain.get_name()
        for p in range(1, numprocs):
            print 'p', p            
            pypar.send(domain_name, p)
    else:
        if verbose: print 'Receiving'

        domain_name = pypar.receive(0)


    if myid == 0:
        # Partition and distribute mesh.
        # Structures returned is in the
        # correct form for the ANUGA data structure


        points, vertices, boundary, quantities,\
                ghost_recv_dict, full_send_dict,\
                = distribute_mesh(domain)

        if verbose: print 'Communication done'
        
    else:
        # Read in the mesh partition that belongs to this
        # processor
        points, vertices, boundary, quantities,\
                ghost_recv_dict, full_send_dict,\
                = rec_submesh(0)



    #------------------------------------------------------------------------
    # Build the domain for this processor using partion structures
    #------------------------------------------------------------------------
    domain = Parallel_Domain(points, vertices, boundary,
                             full_send_dict  = full_send_dict,
                             ghost_recv_dict = ghost_recv_dict)

    #------------------------------------------------------------------------
    # Transfer initial conditions to each subdomain
    #------------------------------------------------------------------------
    for q in quantities:
        domain.set_quantity(q, quantities[q]) 


    #------------------------------------------------------------------------
    # Transfer other attributes to each subdomain
    #------------------------------------------------------------------------

    # FIXME Do them all
    domain.set_name(domain_name)    

    #------------------------------------------------------------------------
    # Return parallel domain to all nodes
    #------------------------------------------------------------------------
    return domain    







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

    # Return structures necessary for building the parallel domain
    return points, vertices, boundary, quantities, \
           ghost_recv_dict, full_send_dict
    




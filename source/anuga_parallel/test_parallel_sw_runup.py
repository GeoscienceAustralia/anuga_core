"""Simple water flow example using ANUGA

Water driven up a linear slope and time varying boundary,
similar to a beach environment

This is a very simple test of the parallel algorithm
"""


#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------

from anuga.pmesh.mesh_interface import create_mesh_from_regions
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
from anuga.shallow_water import Domain
from anuga.shallow_water import Reflective_boundary
from anuga.shallow_water import Dirichlet_boundary
from anuga.shallow_water import Time_boundary
from anuga.shallow_water import Transmissive_boundary

from parallel_api import *

#------------------------------------------------------------------------------
# Read in processor information
#------------------------------------------------------------------------------

numprocs = pypar.size()
myid = pypar.rank()
processor_name = pypar.Get_processor_name()
print 'I am processor %d of %d on node %s' %(myid, numprocs, processor_name)


#------------------------------------------------------------------------------
# Initialise 
#------------------------------------------------------------------------------

if myid == 0:
    #--------------------------------------------------------------------------
    # Setup computational domain
    #--------------------------------------------------------------------------


    points, vertices, boundary = rectangular_cross(10, 10) # Basic mesh
    domain = Domain(points, vertices, boundary) # Create domain


    # Unstructured mesh
    #polygon = [[1,1],[0,1],[0,0],[1,0]]
    #meshname = 'runup.msh'
    #create_mesh_from_regions(polygon,
    #                         boundary_tags={'top': [0],
    #                                        'left': [1],
    #                                        'bottom': [2],
    #                                        'right': [3]},
    #                         maximum_triangle_area=0.01, 
    #                         filename=meshname)
    #domain = Domain(meshname, use_cache=True, verbose = True)


    
    domain.set_name('runup')                    # Set sww filename
    

    #------------ -------------------------------------------------------------
    # Setup initial conditions
    #--------------------------------------------------------------------------

    def topography(x,y): 
        return -x/2                              # linear bed slope

    domain.set_quantity('elevation', topography) # Use function for elevation
    domain.set_quantity('friction', 0.1)         # Constant friction 
    domain.set_quantity('stage', -.4)            # Constant initial stage


    #------------ -------------------------------------------------------------
    # Distribute the domain
    #--------------------------------------------------------------------------
    
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

    print 'Communication done'        
    
else:
    # Read in the mesh partition that belongs to this
    # processor (note that the information is in the
    # correct form for the GA data structure)

    points, vertices, boundary, quantities, ghost_recv_dict, full_send_dict, \
            = rec_submesh(0)



#------------------------------------------------------------------------------
# Start the computations on each subpartion
#------------------------------------------------------------------------------

# Build the domain for this processor
domain = Parallel_Domain(points, vertices, boundary,
                         full_send_dict  = full_send_dict,
                         ghost_recv_dict = ghost_recv_dict)


# Name and dir, etc currently has to be set here as they are not
# transferred from the original domain
domain.set_name('runup')                    # Set sww filename


#------------------------------------------------------------------------------
# Setup initial conditions
#------------------------------------------------------------------------------
for q in quantities:
    domain.set_quantity(q, quantities[q]) # Distribute all quantities    



#------------------------------------------------------------------------------
# Setup parallel boundary conditions
#------------------------------------------------------------------------------

Br = Reflective_boundary(domain)      # Solid reflective wall
Bd = Dirichlet_boundary([-0.2,0.,0.]) # Constant boundary values

# Associate boundary tags with boundary objects
domain.set_boundary({'left': Br, 'right': Bd, 'top': Br, 'bottom': Br,
                     'ghost': None, 'exterior': Bd})



#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------

for t in domain.evolve(yieldstep = 0.1, finaltime = 10.0):
    domain.write_time()
    



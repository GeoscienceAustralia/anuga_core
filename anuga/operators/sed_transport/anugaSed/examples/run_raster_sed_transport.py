"""

Example of use of sediment transport and vegetation drag operators over
a raster-derived topography

M. Perignon
perignon@colorado.edu
July 2014

"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------

import anuga
from anuga import rectangular_cross
from anuga import Domain
from anuga import Dirichlet_boundary, Reflective_boundary


#===============================================================================
# Setup Functions
#===============================================================================

filename_root = 'topo'

# Convert an elevation raster into a point file
anuga.asc2dem(filename_root + '.asc', use_cache = False, verbose = True)
anuga.dem2pts(filename_root + '.dem', use_cache = False, verbose = True)


"""
Include the process-specific quantities when creating the domain
"""        


# import bounding polygon text file, set boundary tags
bounding_polygon = anuga.read_polygon('outline.csv')
boundary_tags = {'bottom':[0],
                  'side1':[1],
				  'side2':[2],
				    'top':[3],
			      'side3':[4],
			      'side4':[5]}

"""
Create the domain
"""

evolved_quantities =  ['stage', 'xmomentum', 'ymomentum', 'concentration']
								
from anuga.pmesh.mesh_interface import create_mesh_from_regions

create_mesh_from_regions(bounding_polygon = bounding_polygon,
                         boundary_tags = boundary_tags,
                         maximum_triangle_area = 200,
                         filename = filename_root + '.msh')

domain = Domain(filename_root + '.msh', evolved_quantities = evolved_quantities)

#------------------------------------------------------------------------------
# Setup parameters of computational domain
#------------------------------------------------------------------------------

# Print some stats about mesh and domain
print 'Number of triangles = ', len(domain)
print 'The extent is ', domain.get_extent()
print domain.statistics()


domain.set_quantity('elevation', 
					filename = filename_root + '.pts',
					use_cache = False,
					verbose = True,
					alpha = 0.1)
                    


domain.set_flow_algorithm('DE0')
domain.set_name('run_raster_sed_transport') # Output name
domain.set_store_vertices_uniquely(True)
domain.set_quantity('stage', expression='elevation')   # Dry initial condition


"""
Store process-specific quantities with same functions
""" 
domain.set_quantities_to_be_stored({'elevation': 2,
                                    'stage': 2,
                                    'xmomentum': 2,
                                    'ymomentum': 2,
                                    'concentration': 2})
                                    
#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
min_elev = domain.quantities['elevation'].vertex_values.min()

Bd = Dirichlet_boundary([1528, 0., 0.])
Bi = Dirichlet_boundary([min_elev - 1, 0., 0.])
Br = Reflective_boundary(domain)

domain.set_boundary({'bottom':Bi,
                      'side1':Br,
                      'side2':Br,
                        'top':Bd,
                      'side3':Br,
                      'side4':Br,
                   'exterior':Br})

#------------------------------------------------------------------------------
# Setup operators
#------------------------------------------------------------------------------

domain.set_quantity('concentration', 0.01)

from anuga.operators.sed_transport_operator import Sed_transport_operator

sed_op = Sed_transport_operator(domain)

#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------
for t in domain.evolve(yieldstep = 1, finaltime = 30.0):
    domain.print_timestepping_statistics()











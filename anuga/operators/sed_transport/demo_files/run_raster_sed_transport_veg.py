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
from anuga import Dirichlet_boundary

"""
Import operators
"""
from anuga.operators.sed_transport.sed_transport_operator import Sed_transport_operator, Vegetation_operator

"""
Import operator-specific boundaries
"""
from anuga.operators.sed_transport.sed_transport_utils import Reflective_boundary_Sed, Dirichlet_boundary_Sed

"""
Import operator-specific version of domain function
"""
from anuga.operators.sed_transport.sed_transport_utils import create_domain_from_regions_sed

"""
Import file conversion and quantity setting functions for vegetation file
"""
from anuga.operators.sed_transport.file_conversion.generic_asc2dem import generic_asc2dem 
from anuga.operators.sed_transport.file_conversion.generic_dem2pts import generic_dem2pts
from anuga.operators.sed_transport.sed_transport_utils import set_quantity_NNeigh



#===============================================================================
# Setup Functions
#===============================================================================

# Convert an elevation raster into a point file
anuga.asc2dem('topo.asc', use_cache = False, verbose = True)
anuga.dem2pts('topo.dem', use_cache = False, verbose = True)


"""
Include the process-specific quantities when creating the domain
"""        
evolved_quantities =  ['stage', 'xmomentum', 'ymomentum', 'concentration']

other_quantities=['elevation', 'friction', 'height', 'xvelocity', \
                  'yvelocity', 'x', 'y', 'vegetation', 'diffusivity']


# import bounding polygon text file, set boundary tags
bounding_polygon = anuga.read_polygon('outline.csv')
boundary_tags = {'bottom':[0],
                  'side1':[1],
				  'side2':[2],
				    'top':[3],
			      'side3':[4],
			      'side4':[5]}

"""
Create the domain with operator-specific function (accepts quantities)
"""
domain = create_domain_from_regions_sed(bounding_polygon,
								boundary_tags = boundary_tags,
								maximum_triangle_area = 200,
								mesh_filename = 'topo.msh',
								interior_regions = {},
								evolved_quantities = evolved_quantities,
								other_quantities = other_quantities,
								use_cache = False,
								verbose = True)


#------------------------------------------------------------------------------
# Setup parameters of computational domain
#------------------------------------------------------------------------------

domain.set_name('run_raster_sed_transport_veg') # Name of sww file

# Print some stats about mesh and domain
print 'Number of triangles = ', len(domain)
print 'The extent is ', domain.get_extent()
print domain.statistics()


domain.set_quantity('elevation', 
					filename = 'topo.pts',
					use_cache = False,
					verbose = True,
					alpha = 0.1)
                    
domain.set_quantity('stage', expression='elevation')

			
#------------------------------------------------------------------------------
# Sediment transport and vegetation operators
#------------------------------------------------------------------------------

"""
Convert a raster of vegetation types into a point file

Set the values of quantity 'vegetation' to values of point file
with Nearest Neighbour algorithm
""" 
generic_asc2dem('veg.asc',
                quantity_name = 'vegetation',
                use_cache = False,
                verbose = True)
generic_dem2pts('veg.dem',
                quantity_name = 'vegetation',
                use_cache = False,
                verbose = True)

set_quantity_NNeigh(domain, 'vegetation', filename = 'veg.pts')




op1 = Sed_transport_operator(domain,
                             erosion = True,
                             deposition = True,
                             turbulence = True,
                             momentum_sinks = True,
                             verbose = True)
                             
op2 = Vegetation_operator(domain,
                          vegfile = 'vegcodes.txt',
                          verbose = True)


domain.set_flow_algorithm('1_75')
domain.set_quantities_to_be_stored({'elevation': 2,'stage': 2,'xmomentum': 2, 'concentration': 2, 'vegetation': 1})



#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
max_elev = domain.quantities['elevation'].vertex_values.max()
min_elev = domain.quantities['elevation'].vertex_values.min()

Bd = Dirichlet_boundary_Sed([1528, 0., 0., 0.2])
Bi = anuga.Dirichlet_boundary([min_elev - 1, 0., 0.])
Br = Reflective_boundary_Sed(domain)

domain.set_boundary({'bottom':Bi,
                      'side1':Br,
                      'side2':Br,
                        'top':Bd,
                      'side3':Br,
                      'side4':Br,
                   'exterior':Br})

#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------

for t in domain.evolve(yieldstep = 5, finaltime = 100):
	print domain.timestepping_statistics()
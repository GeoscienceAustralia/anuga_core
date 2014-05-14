"""Script for running a flood scenario of the WRL Merewether model in Newcastle

Simulation is of 2007 Pasha Bulka flood.

DSR 17/02/2012
WRL2012003.01

Water Research Laboratory, UNSW
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------

# Standard modules
import os
import time
import sys
import numpy as num


# Related major packages
import anuga
from anuga_parallel import myid, distribute, finalize

from anuga_parallel.parallel_operator_factory import Inlet_operator

# Application specific imports
import project                 # Definition of file names and polygons

verbose = project.verbose
use_cache = project.use_cache

#------------------------------------------------------------------------------
# Preparation of topographic data
# Convert ASC 2 DEM 2 PTS using source data and store result in source data
#------------------------------------------------------------------------------

# Filenames
asc_name = 'topography1.asc' 
meshname = 'merewether.msh'
dem_name = 'topography1.dem'

if myid == 0:
    # Create DEM from asc data
    anuga.asc2dem(asc_name, verbose=verbose, use_cache=use_cache)

    # Create pts file for onshore DEM
    anuga.dem2pts(dem_name, verbose=verbose, use_cache=use_cache)

#------------------------------------------------------------------------------
# Create the triangular mesh based on overall clipping polygon with a tagged
# boundary and interior regions defined in project.py along with
# resolutions (maximal area of per triangle) for each polygon
#------------------------------------------------------------------------------

remainder_res = 2.0
houses_res = 1.0
merewether_res = 1.0
holes_res = 1.0
interior_regions = [[project.poly_merewether, merewether_res]]
holes = project.holes

if myid == 0:
    domain = anuga.create_domain_from_regions(project.bounding_polygon,
   		boundary_tags={'bottom': [0],
   			       'right': [1],
                 	       'top': [2],
                 	       'left': [3]},
  				maximum_triangle_area=remainder_res,
				interior_holes=holes,
				mesh_filename=meshname,
   				interior_regions=interior_regions,
   				use_cache=use_cache,
   				verbose=verbose)

    #------------------------------------------------------------------------------
    # Setup initial conditions
    #------------------------------------------------------------------------------
    domain.set_quantity('stage', 0.0)
    domain.set_quantity('friction', 0.02)
    domain.set_quantity('elevation', filename='topography1.pts',
			   	          use_cache=use_cache,
	              			  verbose=verbose)

else:
    domain = None

#------------------------------------------------------------------------------
# Now the sequential domain on processor 0 is distribued to parellel domains
# on each of the processors
#------------------------------------------------------------------------------
domain = distribute(domain)

#------------------------------------------------------------------------------
# Setup computational domain
#------------------------------------------------------------------------------
domain.set_name('merewether_1m') # Name of sww file
domain.set_datadir('.') # Store sww output here
domain.set_minimum_storable_height(0.001) # Store only depth > 1cm



#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------

if myid == 0: print 'Available boundary tags', domain.get_boundary_tags()

Br = anuga.Reflective_boundary(domain)
Bt = anuga.Transmissive_boundary(domain)

domain.set_boundary({'interior': Br,
		     'bottom':   Br,
                     'right':    Bt, # outflow
                     'top':      Bt, # outflow
                     'left':     Br})

line0 = [[382300.0,6354280.], [382300.0,6354300.]]
fixed_inflow = Inlet_operator(domain, line0, 19.7, verbose = False)

#fixed_inflow = anuga.Inflow(domain,
#			center=(382300.0,6354290.0),
#			radius=15.00,
#			rate=19.7)
#domain.forcing_terms.append(fixed_inflow)
#hydrograph = anuga.Inflow(center=(382300.0,6354290.0),radius=30.0,rate=anuga.file_function('test_hydrograph2.tms', quantities=['hydrograph']) 
#domain.forcing_terms.append(hydrograph)

#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------

for t in domain.evolve(yieldstep=10,finaltime=2000):
	if myid == 0: print domain.timestepping_statistics()


domain.sww_merge()


finalize()
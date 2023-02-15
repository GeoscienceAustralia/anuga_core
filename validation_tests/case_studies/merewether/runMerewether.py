"""Script for running a flood scenario of the WRL Merewether model in Newcastle

Simulation is of 2007 Pasha Bulka flood.

DSR 17/02/2012
WRL2012003.01

Water Research Laboratory, UNSW

Edits by Steve Roberts & Gareth Davies, 2014
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
from anuga import myid, distribute, finalize, barrier

from anuga import Inlet_operator
from anuga.utilities import quantity_setting_functions as qs
from anuga import plot_utils as util

# Application specific imports
import project                 # Definition of file names and polygons

#verbose = project.verbose
use_cache = project.use_cache

#------------------------------------------------------------------------------
# Preparation of topographic data
# Convert ASC 2 DEM 2 PTS using source data and store result in source data
#------------------------------------------------------------------------------

# Filenames
zip_name = 'topography1.zip' 
asc_name = 'topography1.asc' 
meshname = 'merewether.msh'
dem_name = 'topography1.dem'

args = anuga.get_args()
alg = args.alg
verbose = args.verbose

if myid == 0:
    # Unzip asc from zip file
    import zipfile as zf
    if verbose: print('Reading ASC from ' + zip_name)
    zf.ZipFile(zip_name).extract(asc_name)
    
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

# Either use houses as holes, or as part of the mesh (with breaklines)
houses_as_holes = False
if houses_as_holes:
    holes = project.holes
    breaklines = []
else:
    # Houses as elevation
    house_height = 3.0
    holes = []
    breaklines = project.holes
    house_addition_poly_fun_pairs = []
    for i in range(len(breaklines)):
        breaklines[i] = breaklines[i] + [breaklines[i][0]]
        house_addition_poly_fun_pairs.append( 
            [ breaklines[i], house_height])
    house_addition_poly_fun_pairs.append(['All', 0.])

#------------------------------------------------------------------------------
# Make domain
#------------------------------------------------------------------------------

if myid == 0:
    domain = anuga.create_domain_from_regions(
            project.bounding_polygon,
            boundary_tags={'bottom': [0],
                           'right': [1],
                           'top': [2],
                           'left': [3]},
            maximum_triangle_area=remainder_res,
            interior_holes=holes,
            breaklines=breaklines,
            mesh_filename=meshname,
            interior_regions=interior_regions,
            use_cache=use_cache,
            verbose=verbose)

    domain.set_zone(project.zone)

    domain.set_flow_algorithm(alg)
    if verbose: print(domain.get_extent(absolute=True))

    #------------------------------------------------------------------------------
    # Setup initial conditions
    #------------------------------------------------------------------------------
    domain.set_quantity('stage', 0.0)

    # Friction -- 2 options
    variable_friction = True
    if not variable_friction:
        # Constant friction
        domain.set_quantity('friction', 0.02)

    else:
        # Set friction to 0.02 on roads, 0.04 elsewhere
        road_polygon = anuga.read_polygon('Road/RoadPolygon.csv')
        friction_function = qs.composite_quantity_setting_function(
            [ [road_polygon, 0.02], ['All', 0.04] ], 
            domain)
        domain.set_quantity('friction', friction_function)

    # Elevation
    if houses_as_holes:
        domain.set_quantity('elevation', filename='topography1.pts',
                              use_cache=use_cache,
                                  verbose=verbose)

    else:
        domain.set_quantity('elevation', filename='topography1.pts',
                              use_cache=use_cache,
                              verbose=verbose, location='vertices')
        # Add house_height inside houses    
        house_addition_function = qs.composite_quantity_setting_function(
            house_addition_poly_fun_pairs, domain)
        domain.add_quantity('elevation', house_addition_function, 
                            location='centroids')
else:
    domain = None

#------------------------------------------------------------------------------
# Now the sequential domain on processor 0 is distribued to parellel domains
# on each of the processors
#------------------------------------------------------------------------------
domain = distribute(domain)

#domain.set_store_vertices_uniquely()

#------------------------------------------------------------------------------
# Setup computational domain
#------------------------------------------------------------------------------
domain.set_name('merewether_1m') # Name of sww file
domain.set_datadir('.') # Store sww output here
#domain.set_minimum_storable_height(0.001) # Store only depth > 1cm

#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------

if myid == 0: print('Available boundary tags', domain.get_boundary_tags())

Br = anuga.Reflective_boundary(domain)
Bt = anuga.Transmissive_boundary(domain)

boundary_dict = {'bottom':   Br,
                     'right':    Bt, # outflow
                     'top':      Bt, # outflow
                     'left':     Br}

if houses_as_holes:
    boundary_dict['interior'] = Br

domain.set_boundary(boundary_dict)

#-------------------------------------------------------------------------------
# Setup Inlet
#-------------------------------------------------------------------------------
#line0 = [[382300.0,6354280.], [382300.0,6354300.]]
line0 = [[382275.0,6354270.], [382255.0,6354290.]]
import math


#velocity = [3.5/math.sqrt(2.0), 3.5/math.sqrt(2.0)]
#fixed_inflow = Inlet_operator(domain, line0, 19.7, verbose = False)

center = (382265.0, 6354280.0)
radius = 10.0
region0 = anuga.Region(domain, center=center, radius=radius)
fixed_inflow = Inlet_operator(domain, region0 , 19.7, verbose = True)
#fixed_inflow = anuga.Inflow(domain,
#           center=(382300.0,6354290.0),
#           radius=15.00,
#           rate=19.7)
#domain.forcing_terms.append(fixed_inflow)
#hydrograph = anuga.Inflow(center=(382300.0,6354290.0),radius=30.0,rate=anuga.file_function('test_hydrograph2.tms', quantities=['hydrograph']) 
#domain.forcing_terms.append(hydrograph)

#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------

for t in domain.evolve(yieldstep=10,finaltime=1000):#1000
    if myid == 0: print(domain.timestepping_statistics())


domain.sww_merge(delete_old=True)

#------------------------------------------------------------------------------
# Make geotif output
#------------------------------------------------------------------------------

barrier()
if myid==0:
    util.Make_Geotif('merewether_1m.sww', 
                     output_quantities=['depth', 'velocity', 
                                        'friction', 'elevation'],
                     myTimeStep='last',
                     CellSize=1.0,
                     EPSG_CODE=32756,
                     bounding_polygon=project.bounding_polygon,
                     k_nearest_neighbours=1)

finalize()

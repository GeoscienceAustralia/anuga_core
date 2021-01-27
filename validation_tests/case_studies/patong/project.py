"""
This file contains all your file and directory definitions 
for elevation, meshes and outputs.
"""

import os
from anuga.utilities.system_tools import get_user_name, get_host_name
from time import localtime, strftime, gmtime
from os.path import join, exists
from anuga.geometry.polygon import read_polygon


scenario_name = 'patong'

#------------------------------------------------------------------------------
# Initial Conditions
#------------------------------------------------------------------------------

# Model specific parameters.
# One or all can be changed each time the run_model script is executed
tide = 0.8              # difference between MSL and HAT
zone = 47               # specify zone of model
model = 'small'         # different size model 'small' or 'large'
event_number = 'mux_' + model    # the event number or the mux file name
alpha = 0.1             # smoothing parameter for mesh
friction = 0.025         # manning's friction coefficient
starttime = 0           # start time for simulation
finaltime = 10000       # final time for simulation

nameFlag='xxNameFlagxx'

use_buildings = False   # Add buildings to the DEM? This is slow initially 
setup = 'basic'         # This can be one of three values
                        #    trial - coarsest mesh, fast
                        #    basic - coarse mesh
                        #    final - fine mesh, slowest
#------------------------------------------------------------------------------
# Output filename
#
# Your output filename should be unique between different runs on different
# data.
# The list of items below will be used to create a file in your output
# directory.
# Your user name and time+date will be automatically added.  For example,
#     [setup, tide, event_number]
# will result in a filename like
#     20090212_091046_run_final_0_27283_rwilson
#------------------------------------------------------------------------------

output_comment = [setup, tide, event_number, nameFlag]

#------------------------------------------------------------------------------
# Input Data
#------------------------------------------------------------------------------

# ELEVATION DATA
ascii_grid_filenames = []

# Format for point is x,y,elevation (with header)
point_filenames = ['patong_10m_grid_msl_Project.txt',    # 10m main grid
                   'patong_10m_small_grid_for_anuga_sub_Project.txt', # 10m saddle
                   'patong_bay_1s_grid_Project.txt',   # 1s grid
                   'andaman_3s_grid_Clip2_Project.txt'] #3s grid

extent_polygon_filenames = ['patong_10m.txt',
                            'saddle_10m.txt',
                            'patong_1s.txt']


# BOUNDING POLYGON - for data clipping and estimate of triangles in mesh
# Format for points easting,northing (no header)
bounding_polygon_filename = 'bounding_polygon_'+ model + '.csv'
bounding_polygon_maxarea = 100000

# INTERIOR REGIONS -  for designing the mesh
# Used in run_model.py
# Format for points easting,northing (no header)                    
interior_regions_data = [['inundation_area.csv', 75],
                         ['aoi.csv', 200],
                         ['aos.csv', 900],
                         ['sw_coast_original.csv', 7000],
                         ['building_main.csv', 20],
                         ['building_main_small.csv', 20],
                         ['building_main_south.csv', 20],
                         ['building_saddle.csv', 20]]


# LAND - used to set the initial stage/water to be offcoast only
# Used in run_model.py.  Format for points easting,northing (no header)
land_initial_conditions_filename = [['initial_conditions_' + model +'.csv', 0]]

# GAUGES - for creating timeseries at a specific point
# Used in get_timeseries.py.  
# Format easting,northing,name,elevation (with header)
gauges_filename = 'gauges/validation_gauges.csv'

# BUILDINGS POLYOGN - for elevation of buildings
# Used in run_model.py
# Format easting,northing,id,floor (with header)
building_polygon_filename = 'buildings.csv' 


# BUILDINGS EXPOSURE - for identifying inundated houses
# Used in run_building_inundation.py
# Format latitude,longitude etc (geographic)
building_exposure_filename = 'busselton_res_clip.csv' # from NEXIS

# BOUNDING POLYGON - used in build_boundary.py and run_model.py respectively
# NOTE: when files are put together the points must be in sequence
# For ease go clockwise!
# Check the run_model.py for boundary_tags

# Thinned ordering file from Hazard Map (geographic)
# Format is index,latitude,longitude (with header)
urs_order_filename = 'urs_order_'+ model +'.csv'

# Landward bounding points
# Format easting,northing (no header)
landward_boundary_filename = 'landward_boundary_'+ model +'.csv'

# MUX input filename.
# If a meta-file from EventSelection is used, set 'multi-mux' to True.
# If a single MUX stem filename (*.grd) is used, set 'multi-mux' to False.
mux_input_filename = event_number # to be found in event_folder
                                    # (ie boundaries/event_number/)
multi_mux = False
##mux_input_filename = 'event.list'
##multi_mux = True

#------------------------------------------------------------------------------
# Clipping regions for export to asc and regions for clipping data
# Final inundation maps should only be created in regions of the finest mesh
#------------------------------------------------------------------------------

#CBD extract ascii grid - coordinates from patong_1s extent
xminCBD = 417445.1119
xmaxCBD = 425601.7881
yminCBD = 870663.4547
ymaxCBD = 876965.3856



################################################################################
################################################################################
####         NOTE: NOTHING WOULD NORMALLY CHANGE BELOW THIS POINT.          ####
################################################################################
################################################################################

# Get system user and host names.
# These values can be used to distinguish between two similar runs by two
# different users or runs by the same user on two different machines.
user = get_user_name()
host = get_host_name()

#-------------------------------------------------------------------------------
# Output Elevation Data
#-------------------------------------------------------------------------------

# Output filename for elevation
# this is a combination of all the data generated in build_elevation.py
combined_elevation_basename = scenario_name + '_combined_elevation_' + model

#-------------------------------------------------------------------------------
# Directory Structure
#-------------------------------------------------------------------------------

# determines time for setting up output directories
time = strftime('%Y%m%d_%H%M%S', localtime()) 
gtime = strftime('%Y%m%d_%H%M%S', gmtime()) 
build_time = time + '_build_' + model
run_time = time + '_run_' + model

# create paths for data files.
output_dirname = os.path.dirname(__file__)
#home = os.path.join(output_dirname, 'local_data', 'data')
home = output_dirname
muxhome = home # FIXME (Ole): Get rid off

    
# check various directories/files that must exist
#anuga_folder = join(home, state, scenario_folder, 'anuga')
anuga_folder = home
topographies_folder = join(anuga_folder, 'topographies')
polygons_folder = join(anuga_folder, 'polygons')
boundaries_folder = join(anuga_folder, 'boundaries')
output_folder = join(anuga_folder, 'outputs')
gauges_folder = join(anuga_folder, 'gauges')
meshes_folder = join(anuga_folder, 'meshes')
event_folder = join(boundaries_folder, str(event_number))

# MUX data files
# Directory containing the MUX data files to be used with EventSelection.
mux_data_folder = join(muxhome, 'mux')

#------------------------------------------------------------------------------
# Location of input and output data
#------------------------------------------------------------------------------

# Convert the user output_comment to a string for run_model.py
output_comment = ('_'.join([str(x) for x in output_comment if x != user])
                  + '_' + user)

# The absolute pathname of the all elevation, generated in build_elevation.py
combined_elevation = join(topographies_folder, combined_elevation_basename)

# The absolute pathname of the mesh, generated in run_model.py
meshes = join(meshes_folder, scenario_name) + model +'.msh'

# The pathname for the urs order points, used within build_urs_boundary.py
urs_order = join(boundaries_folder, urs_order_filename)

# The absolute pathname for the landward points of the bounding polygon,
# Used within run_model.py)
landward_boundary = join(boundaries_folder, landward_boundary_filename)

# The absolute pathname for the .sts file, generated in build_boundary.py
event_sts = join(event_folder, scenario_name)

# The absolute pathname for the output folder names
# Used for build_elevation.py

output_build = join(output_folder, build_time) + '_' + user
#output_build = output_folder 

# Used for run_model.py
output_run = join(output_folder, run_time)
#output_run = output_folder

# Used by post processing
output_run_time = join(output_run, scenario_name)
#output_run_time = output_run

# The absolute pathname for the gauges file
# Used for get_timeseries.py
gauges = join(gauges_folder, gauges_filename)

# The absolute pathname for the gauges file
# Used for run_model.py
building_polygon = join(polygons_folder, building_polygon_filename)

# The absolute pathname for the building file
# Used for run_building_inundation.py
building_exposure = join(gauges_folder, building_exposure_filename)

# full path to where MUX files (or meta-files) live
mux_input = join(event_folder, mux_input_filename)

###########################################################################
#
# Code below here used to be in setup_model.py
#
###########################################################################

sanity_error=False

if setup == 'trial':
    scale_factor = 100
    yieldstep = 30.
elif setup == 'basic': 
    scale_factor = 4
    yieldstep = 120.
elif setup == 'final': 
    scale_factor = 1
    yieldstep = 30.
else:
    print ("Sorry, you must set the 'setup' variable to one of:"
           '   trial - coarsest mesh, fast\n'
           '   basic - coarse mesh\n'
           '   final - fine mesh, slowest\n'
           '\n'
           "'setup' was set to '%s'" % setup)
    sanity_error = True

#-------------------------------------------------------------------------------
# Check for errors detected above.
#-------------------------------------------------------------------------------

if sanity_error:
    msg = 'You must fix the above errors before continuing.'
    raise( Exception, msg)

#-------------------------------------------------------------------------------
# Reading polygons and creating interior regions
#-------------------------------------------------------------------------------

# Create list of land polygons with initial conditions
land_initial_conditions = []
for filename, MSL in land_initial_conditions_filename:
    polygon = read_polygon(join(polygons_folder, filename))
    land_initial_conditions.append([polygon, MSL])

# Create list of interior polygons with scaling factor
interior_regions = []
for filename, maxarea in interior_regions_data:
    polygon = read_polygon(join(polygons_folder, filename))
    interior_regions.append([polygon,
                                     maxarea*scale_factor])

# Initial bounding polygon for data clipping 
bounding_polygon = read_polygon(join(polygons_folder,
                                             bounding_polygon_filename))
bounding_maxarea = bounding_polygon_maxarea*scale_factor

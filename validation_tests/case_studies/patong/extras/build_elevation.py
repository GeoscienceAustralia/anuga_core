"""Build the elevation data to run a tsunami inundation scenario 
for Busselton, WA, Australia.

Input: elevation data from project.py
Output: pts file stored in project.topographies_dir 
The run_model.py is reliant on the output of this script.

"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------

# Standard modules
import os
from os.path import join

# Related major packages
from anuga.file_conversion.asc2dem import _convert_dem_from_ascii2netcdf as convert_dem_from_ascii2netcdf
from anuga.file_conversion.dem2pts import dem2pts
#from anuga. import start_screen_catcher
from anuga.utilities.file_utils import copy_code_files
from anuga.geospatial_data.geospatial_data import Geospatial_data
from anuga.geometry.polygon import read_polygon

# Application specific imports
import project


#------------------------------------------------------------------------------
# Copy scripts to time stamped output directory and capture screen
# output to file
#------------------------------------------------------------------------------
copy_code_files(project.output_build,__file__, 
                os.path.dirname(project.__file__)+os.sep+\
                project.__name__+'.py' )
#start_screen_catcher(project.output_build)


#------------------------------------------------------------------------------
# Preparation of topographic data
#
# Convert ASC 2 DEM 2 PTS using source data and store result in source data
# Do for coarse and fine data
# Fine pts file to be clipped to area of interest
#------------------------------------------------------------------------------

print 'project.bounding_polygon', project.bounding_polygon
print 'project.combined_elevation_basename', project.combined_elevation_basename

# Create Geospatial data from ASCII files
geospatial_data = {}

for filename in project.ascii_grid_filenames:
    absolute_filename = join(project.topographies_folder, filename)
    convert_dem_from_ascii2netcdf(absolute_filename,
                                  basename_out=absolute_filename,
                                  use_cache=True,
                                  verbose=True)
    dem2pts(absolute_filename, use_cache=True, verbose=True)

    G_grid = Geospatial_data(file_name=absolute_filename+'.pts',
                                                verbose=True)
    print 'Clip geospatial object'
    geospatial_data[filename] = G_grid.clip(project.bounding_polygon)

# Create Geospatial data from TXT files

for filename in project.point_filenames:
    absolute_filename = join(project.topographies_folder, filename)
    G_points = Geospatial_data(file_name=absolute_filename,
                                                verbose=True)
    print 'Clip geospatial object'
    geospatial_data[filename] = G_points.clip(project.bounding_polygon)

#-------------------------------------------------------------------------------
# Combine, clip and export dataset 
#-------------------------------------------------------------------------------
extent_polygons = []
for extent_polygon_filename in project.extent_polygon_filenames:
    p = read_polygon(join(project.polygons_folder, extent_polygon_filename))
    extent_polygons.append(p)
    
print 'Add geospatial objects' 
G = None
for key in geospatial_data:
    if key == project.point_filenames[0] or key == project.point_filenames[1]:
        G += geospatial_data[key]
    elif key == project.point_filenames[2]:
        D = geospatial_data[key]
        D = D.clip_outside(extent_polygons[0])
        D = D.clip_outside(extent_polygons[1])
        G += D
    elif key == project.point_filenames[3]:
        D = geospatial_data[key]
        D = D.clip_outside(extent_polygons[2])
        G += D

print 'Export combined DEM file'
G.export_points_file(project.combined_elevation + '.pts')
print 'Do txt version too'
# Use for comparision in ARC
G.export_points_file(project.combined_elevation + '.txt')


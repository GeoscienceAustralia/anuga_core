"""Create mesh and time boundary for the Okushiri island validation
"""


import numpy as num

from anuga.pmesh.mesh import *
from anuga.pmesh.mesh_interface import create_mesh_from_regions
from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.geospatial_data import Geospatial_data
from anuga.config import netcdf_float

import project

# configure my logging
import anuga.utilities.log as log
#log.console_logging_level = log.INFO
#log.log_logging_level = log.DEBUG
log.log_filename = './create_okushiri.log'


#--------------------------------------------------------------------------
# Create the triangular mesh based on overall clipping polygon with a
# tagged
# boundary and interior regions defined in project.py along with
# resolutions (maximal area of per triangle) for each polygon
#--------------------------------------------------------------------------

#base_resolution = 0.01 #  triangles
#base_resolution = 0.02 # 1,972,289 triangles
#base_resolution = 0.04 # 989,669 triangles
#base_resolution = 0.1 # 397,456 triangles
#base_resolution = 0.6 # 70703 triangles
base_resolution = 2.0 # 21884 triangles
#base_resolution = 4.0 # 11388 triangles

#print base_resolution

# Basic geometry and bounding polygon
xleft   = 0
xright  = 5.448
ybottom = 0
ytop    = 3.402


point_sw = [xleft, ybottom]
point_se = [xright, ybottom]
point_nw = [xleft, ytop]    
point_ne = [xright, ytop]

bounding_polygon = [point_se,
                    point_ne,
                    point_nw,
                    point_sw]


# Localised refined area for gulleys
xl = 4.8
xr = 5.3
yb = 1.6
yt = 2.3
p0 = [xl, yb]
p1 = [xr, yb]
p2 = [xr, yt]
p3 = [xl, yt]

gulleys = [p0, p1, p2, p3]



# Island area and drawdown region
island_0 = [xleft + 2*(xright-xleft)/3+1.2, ytop-0.5]
island_1 = [xleft + 2*(xright-xleft)/3+0.5, ybottom + 2*(ytop-ybottom)/3]
island_2 = [xleft + (xright-xleft)/2+0.4, ybottom + 2*(ytop-ybottom)/3-0.3]
island_3 = [xleft + (xright-xleft)/2+0.4, ybottom + (ytop-ybottom)/3+0.3]
island_4 = [xleft + 2*(xright-xleft)/3+0.4, ybottom + (ytop-ybottom)/3-0.3]
island_5 = [xleft + 2*(xright-xleft)/3+1.2, ybottom+0.8]
island_6 = [xl-.01, yb]  # Keep right edge just off the gulleys
island_7 = [xl-.01, yt]

island = [island_0, island_1, island_2,
          island_3, island_4, island_5,
          island_6, island_7]



# Region spanning half right hand side of domain just inside boundary
rhs_nw = [xleft + (xright-xleft)/3+1, ytop-1.4]
rhs_sw = [xleft + (xright-xleft)/3+1, ybottom+0.5]
rhs_se = [xright-0.1, ybottom+0.2]
rhs_ne = [xright-0.1, ytop-0.2]        

rhs_region = [rhs_nw, rhs_ne, rhs_se, rhs_sw]


# Interior regions and creation of mesh
interior_regions = [[rhs_region, 0.0005*base_resolution],
                    [island, 0.0003*base_resolution],
                    [gulleys, 0.00003*base_resolution]]    




def prepare_timeboundary(filename, verbose = False):
    """Convert benchmark 2 time series to NetCDF tms file.
    This is a 'throw-away' code taylor made for files like
    'Benchmark_2_input.txt' from the LWRU2 benchmark
    """

    from anuga.file.netcdf import NetCDFFile

    if verbose: print('Creating', filename)

    # Read the ascii (.txt) version of this file
    fid = open(filename[:-4] + '.txt')

    # Skip first line
    line = fid.readline()

    # Read remaining lines
    lines = fid.readlines()
    fid.close()


    N = len(lines)
    T = num.zeros(N, float)  #Time
    Q = num.zeros(N, float)  #Values

    for i, line in enumerate(lines):
        fields = line.split()

        T[i] = float(fields[0])
        Q[i] = float(fields[1])


    # Create tms NetCDF file

    fid = NetCDFFile(filename, 'w')
    fid.institution = 'Geoscience Australia'
    fid.description = 'Input wave for Benchmark 2'
    fid.starttime = 0.0
    fid.createDimension('number_of_timesteps', len(T))
    fid.createVariable('time', netcdf_float, ('number_of_timesteps',))
    fid.variables['time'][:] = T

    fid.createVariable('stage', netcdf_float, ('number_of_timesteps',))
    fid.variables['stage'][:] = Q[:]

    fid.createVariable('xmomentum', netcdf_float, ('number_of_timesteps',))
    fid.variables['xmomentum'][:] = 0.0

    fid.createVariable('ymomentum', netcdf_float, ('number_of_timesteps',))
    fid.variables['ymomentum'][:] = 0.0

    fid.close()


def create_mesh(elevation_in_mesh=False, verbose=False):
    # Prepare time boundary
    prepare_timeboundary(project.boundary_filename, verbose)

    

    

    meshname = project.mesh_filename + '.msh'
    m = create_mesh_from_regions(bounding_polygon,
                                 boundary_tags={'wall': [0, 1, 3],
                                                'wave': [2]},     
                                 maximum_triangle_area=0.01*base_resolution,
                                 interior_regions=interior_regions,
                                 filename=project.mesh_filename,
                                 use_cache=False,
                                 verbose=verbose)


    
    if elevation_in_mesh is True:
        from anuga.fit_interpolate.fit import fit_to_mesh_file


        if verbose: print('Reading xya from zip')
        import zipfile as zf
        zf.ZipFile(project.bathymetry_filename_stem+'.zip').\
          extract(project.bathymetry_filename_stem+'.xya')

        if verbose: print('Reading pts from xya')
        anuga.xya2pts(project.bathymetry_filename_stem+'.xya',\
                      verbose = verbose)
        
        fit_to_mesh_file(project.mesh_filename, project.bathymetry_filename,
                         project.mesh_filename, verbose = verbose)


#-------------------------------------------------------------
if __name__ == "__main__":

    print('Create mesh')
    create_mesh(elevation_in_mesh=True, verbose=True)



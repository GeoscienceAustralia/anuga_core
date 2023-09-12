

import numpy as num
from anuga.file.netcdf import NetCDFFile


import anuga
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular
from anuga.shallow_water.shallow_water_domain import Domain
from anuga.shallow_water.boundaries import Reflective_boundary
from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.shallow_water.forcing import *
from anuga.utilities.numerical_tools import ensure_numeric

from anuga.file.sww import Write_sww  
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a, \
                            netcdf_float

def sts2sww_mesh(basename_in, basename_out=None, 
                 spatial_thinning=1, verbose=False):
    
    from anuga.mesh_engine.mesh_engine import NoTrianglesError
    from anuga.pmesh.mesh import Mesh
    if verbose:
        print("Starting sts2sww_mesh")
    
    mean_stage=0.
    zscale=1.

    if (basename_in[:-4]=='.sts'):
        stsname = basename_in
    else: 
        stsname = basename_in + '.sts'

    if verbose: print("Reading sts NetCDF file: %s" %stsname)
    infile = NetCDFFile(stsname, netcdf_mode_r)
    cellsize = infile.cellsize
    ncols = infile.ncols
    nrows = infile.nrows
    no_data = infile.no_data
    refzone = infile.zone
    x_origin = infile.xllcorner
    y_origin = infile.yllcorner
    origin = num.array([x_origin, y_origin])
    x = infile.variables['x'][:]
    y = infile.variables['y'][:]
    times = infile.variables['time'][:]
    wind_speed_full = infile.variables['wind_speed'][:]
    wind_angle_full = infile.variables['wind_angle'][:]
    pressure_full   =   infile.variables['barometric_pressure'][:]
    infile.close()

    number_of_points = nrows*ncols
    points_utm = num.zeros((number_of_points,2),float)
    points_utm[:,0]=x+x_origin
    points_utm[:,1]=y+y_origin

    thinned_indices=[]
    for i in range(number_of_points):
        if (i//ncols==0 or i//ncols==ncols-1 or (i//ncols)%(spatial_thinning)==0):
            if ( i%(spatial_thinning)==0 or i%nrows==0 or i%nrows==nrows-1 ):  
                thinned_indices.append(i)

    #Spatial thinning
    points_utm=points_utm[thinned_indices]
    number_of_points = points_utm.shape[0]
    number_of_timesteps = wind_speed_full.shape[0]
    wind_speed = num.empty((number_of_timesteps,number_of_points),dtype=float)
    wind_angle = num.empty((number_of_timesteps,number_of_points),dtype=float)
    barometric_pressure   = num.empty((number_of_timesteps,number_of_points),dtype=float)
    if verbose:
        print("Total number of points: ", nrows*ncols)
        print("Number of thinned points: ", number_of_points)
    for i in range(number_of_timesteps):
        wind_speed[i] = wind_speed_full[i,thinned_indices]
        wind_angle[i] = wind_angle_full[i,thinned_indices]
        barometric_pressure[i]   = pressure_full[i,thinned_indices]

    #P.plot(points_utm[:,0],points_utm[:,1],'ro')
    #P.show()

    if verbose:
        print("Generating sww triangulation of gems data")

    mesh = Mesh()
    mesh.add_vertices(points_utm)
    mesh.auto_segment(smooth_indents=True, expand_pinch=True)
    mesh.auto_segment(mesh.shape.get_alpha() * 1.1)
    try:
        mesh.generate_mesh(minimum_triangle_angle=0.0, verbose=False)
    except NoTrianglesError:
        # This is a bit of a hack, going in and changing the data structure.
        mesh.holes = []
        mesh.generate_mesh(minimum_triangle_angle=0.0, verbose=False)

    mesh_dic = mesh.Mesh2MeshList()

    points_utm=ensure_numeric(points_utm)
    assert num.all(ensure_numeric(mesh_dic['generatedpointlist'])
                       == ensure_numeric(points_utm))

    volumes = mesh_dic['generatedtrianglelist']

    # Write sww intro and grid stuff.
    if (basename_out is not None and basename_out[:-4]=='.sww'): 
        swwname = basename_out
    else: 
        swwname = basename_in + '.sww'

    if verbose: 'Output to %s' % swwname

    if verbose:
        print("Writing sww wind and pressure field file")
    outfile = NetCDFFile(swwname, netcdf_mode_w)
    sww = Write_sww([], ['wind_speed','wind_angle','barometric_pressure'])
    sww.store_header(outfile, times, len(volumes), len(points_utm),
                     verbose=verbose, sww_precision='d')
    outfile.mean_stage = mean_stage
    outfile.zscale = zscale
    sww.store_triangulation(outfile, points_utm, volumes,
                            refzone,  
                            new_origin=origin, #check effect of this line
                            verbose=verbose)

    if verbose: 
        print('Converting quantities')
    
    # Read in a time slice from the sts file and write it to the SWW file

    #print wind_angle[0,:10]
    for i in range(len(times)):
        sww.store_quantities(outfile,
                             slice_index=i,
                             verbose=verbose,
                             wind_speed=wind_speed[i,:],
                             wind_angle=wind_angle[i,:],
                             barometric_pressure=barometric_pressure[i,:],
                             sww_precision=float)

    if verbose: 
        sww.verbose_quantities(outfile)
    outfile.close()

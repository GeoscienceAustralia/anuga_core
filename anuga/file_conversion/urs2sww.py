from builtins import map
from builtins import range
import os
import numpy as num
from anuga.file.netcdf import NetCDFFile

from anuga.file.urs import Read_urs

import anuga.utilities.log as log

from anuga.file_conversion.urs2nc import urs2nc
from anuga.file_conversion.ferret2sww import ferret2sww

from anuga.coordinate_transforms.redfearn import redfearn, \
     convert_from_latlon_to_utm

from anuga.geospatial_data.geospatial_data import ensure_absolute, \
                                                    Geospatial_data

from anuga.file.mux import WAVEHEIGHT_MUX_LABEL, EAST_VELOCITY_LABEL, \
                            NORTH_VELOCITY_LABEL
                            
from anuga.utilities.numerical_tools import ensure_numeric                            

from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a, \
                            netcdf_float

from anuga.file.sww import Write_sww  

###############################################################

def urs_ungridded2sww(basename_in='o', basename_out=None, verbose=False,
                      mint=None, maxt=None,
                      mean_stage=0,
                      origin=None,
                      hole_points_UTM=None,
                      zscale=1):
    r"""
    Convert URS C binary format for wave propagation to
    sww format native to abstract_2d_finite_volumes.

    Specify only basename_in and read files of the form
    basefilename-z-mux, basefilename-e-mux and
    basefilename-n-mux containing relative height,
    x-velocity and y-velocity, respectively.

    Also convert latitude and longitude to UTM. All coordinates are
    assumed to be given in the GDA94 datum. The latitude and longitude
    information is assumed ungridded grid.

    min's and max's: If omitted - full extend is used.
    To include a value min ans max may equal it.
    Lat and lon are assumed to be in decimal degrees.

    origin is a 3-tuple with geo referenced
    UTM coordinates (zone, easting, northing)
    It will be the origin of the sww file. This shouldn't be used,
    since all of anuga should be able to handle an arbitary origin.
    The mux point info is NOT relative to this origin.

    URS C binary format has data organised as TIME, LONGITUDE, LATITUDE
    which means that latitude is the fastest
    varying dimension (row major order, so to speak)

    In URS C binary the latitudes and longitudes are in assending order.

    Note, interpolations of the resulting sww file will be different
    from results of urs2sww.  This is due to the interpolation
    function used, and the different grid structure between urs2sww
    and this function.

    Interpolating data that has an underlying gridded source can
    easily end up with different values, depending on the underlying
    mesh.

    consider these 4 points
    50  -50

    0     0

    The grid can be
     -
    |\|   A
     -
     or;
      -
     |/|  B
      -

    If a point is just below the center of the midpoint, it will have a
    +ve value in grid A and a -ve value in grid B.
    """

    from anuga.mesh_engine.mesh_engine import NoTrianglesError
    from anuga.pmesh.mesh import Mesh

    files_in = [basename_in + WAVEHEIGHT_MUX_LABEL,
                basename_in + EAST_VELOCITY_LABEL,
                basename_in + NORTH_VELOCITY_LABEL]
    quantities = ['HA','UA','VA']

    # instantiate urs_points of the three mux files.
    mux = {}
    for quantity, file in zip(quantities, files_in):
        mux[quantity] = Read_urs(file)

    # Could check that the depth is the same. (hashing)

    # handle to a mux file to do depth stuff
    a_mux = mux[quantities[0]]

    # Convert to utm
    lat = a_mux.lonlatdep[:,1]
    long = a_mux.lonlatdep[:,0]
    points_utm, zone = convert_from_latlon_to_utm(latitudes=lat,
                                                  longitudes=long)

    elevation = a_mux.lonlatdep[:,2] * -1

    # grid (create a mesh from the selected points)
    # This mesh has a problem.  Triangles are streched over ungridded areas.
    # If these areas could be described as holes in pmesh, that would be great.

    # I can't just get the user to selection a point in the middle.
    # A boundary is needed around these points.
    # But if the zone of points is obvious enough auto-segment should do
    # a good boundary.
    mesh = Mesh()
    mesh.add_vertices(points_utm)
    mesh.auto_segment(smooth_indents=True, expand_pinch=True)

    # To try and avoid alpha shape 'hugging' too much
    mesh.auto_segment(mesh.shape.get_alpha() * 1.1)
    if hole_points_UTM is not None:
        point = ensure_absolute(hole_points_UTM)
        mesh.add_hole(point[0], point[1])

    try:
        mesh.generate_mesh(minimum_triangle_angle=0.0, verbose=False)
    except NoTrianglesError:
        # This is a bit of a hack, going in and changing the data structure.
        mesh.holes = []
        mesh.generate_mesh(minimum_triangle_angle=0.0, verbose=False)

    mesh_dic = mesh.Mesh2MeshList()

    #mesh.export_mesh_file(basename_in + '_168.tsh')
    #import sys; sys.exit()
    # These are the times of the mux file
    mux_times = []
    for i in range(a_mux.time_step_count):
        mux_times.append(a_mux.time_step * i)
    (mux_times_start_i, mux_times_fin_i) = read_time_from_mux(mux_times, mint, maxt)
    times = mux_times[mux_times_start_i:mux_times_fin_i]

    if mux_times_start_i == mux_times_fin_i:
        # Close the mux files
        for quantity, file in zip(quantities, files_in):
            mux[quantity].close()
        msg = "Due to mint and maxt there's no time info in the boundary SWW."
        raise Exception(msg)

    # If this raise is removed there is currently no downstream errors

    points_utm=ensure_numeric(points_utm)
    assert num.all(ensure_numeric(mesh_dic['generatedpointlist'])
                       == ensure_numeric(points_utm))

    volumes = mesh_dic['generatedtrianglelist']

    # Write sww intro and grid stuff.
    if basename_out is None:
        swwname = basename_in + '.sww'
    else:
        swwname = basename_out + '.sww'

    if verbose: log.critical('Output to %s' % swwname)

    outfile = NetCDFFile(swwname, netcdf_mode_w)

    # For a different way of doing this, check out tsh2sww
    # work out sww_times and the index range this covers
    sww = Write_sww(['elevation'], ['stage', 'xmomentum', 'ymomentum'])
    sww.store_header(outfile, times, len(volumes), len(points_utm),
                     verbose=verbose, sww_precision=netcdf_float)
    outfile.mean_stage = mean_stage
    outfile.zscale = zscale

    sww.store_triangulation(outfile, points_utm, volumes,
                            zone,  
                            new_origin=origin,
                            verbose=verbose)
    sww.store_static_quantities(outfile, elevation=elevation)

    if verbose: log.critical('Converting quantities')

    # Read in a time slice from each mux file and write it to the SWW file
    j = 0
    for ha, ua, va in zip(mux['HA'], mux['UA'], mux['VA']):
        if j >= mux_times_start_i and j < mux_times_fin_i:
            stage = zscale*ha + mean_stage
            h = stage - elevation
            xmomentum = ua*h
            ymomentum = -1 * va * h # -1 since in mux files south is positive.
            sww.store_quantities(outfile,
                                 slice_index=j-mux_times_start_i,
                                 verbose=verbose,
                                 stage=stage,
                                 xmomentum=xmomentum,
                                 ymomentum=ymomentum,
                                 sww_precision=float)
        j += 1

    if verbose: sww.verbose_quantities(outfile)

    outfile.close()

    # Do some conversions while writing the sww file


def read_time_from_mux(mux_times, mint, maxt):
    """
        Read a list of mux times.
        Return start and finish times which lie within the passed time period.
    """

    if mint is None:
        mux_times_start_i = 0
    else:
        mux_times_start_i = num.searchsorted(mux_times, mint)

    if maxt is None:
        mux_times_fin_i = len(mux_times)
    else:
        maxt += 0.5 # so if you specify a time where there is
                    # data that time will be included
        mux_times_fin_i = num.searchsorted(mux_times, maxt)

    return mux_times_start_i, mux_times_fin_i





def urs2sww(basename_in='o', basename_out=None, verbose=False,
            remove_nc_files=True,
            minlat=None, maxlat=None,
            minlon=None, maxlon=None,
            mint=None, maxt=None,
            mean_stage=0,
            origin=None,
            zscale=1,
            fail_on_NaN=True,
            NaN_filler=0):
    """Convert a URS file to an SWW file.
    Convert URS C binary format for wave propagation to
    sww format native to abstract_2d_finite_volumes.

    Specify only basename_in and read files of the form
    basefilename-z-mux2, basefilename-e-mux2 and
    basefilename-n-mux2 containing relative height,
    x-velocity and y-velocity, respectively.

    Also convert latitude and longitude to UTM. All coordinates are
    assumed to be given in the GDA94 datum. The latitude and longitude
    information is for  a grid.

    min's and max's: If omitted - full extend is used.
    To include a value min may equal it, while max must exceed it.
    Lat and lon are assumed to be in decimal degrees.
    NOTE: minlon is the most east boundary.

    origin is a 3-tuple with geo referenced
    UTM coordinates (zone, easting, northing)
    It will be the origin of the sww file. This shouldn't be used,
    since all of anuga should be able to handle an arbitary origin.

    URS C binary format has data orgainised as TIME, LONGITUDE, LATITUDE
    which means that latitude is the fastest
    varying dimension (row major order, so to speak)

    In URS C binary the latitudes and longitudes are in assending order.
    """

    if basename_out is None:
        basename_out = basename_in

    files_out = urs2nc(basename_in, basename_out)

    ferret2sww(basename_out,
               minlat=minlat,
               maxlat=maxlat,
               minlon=minlon,
               maxlon=maxlon,
               mint=mint,
               maxt=maxt,
               mean_stage=mean_stage,
               origin=origin,
               zscale=zscale,
               fail_on_NaN=fail_on_NaN,
               NaN_filler=NaN_filler,
               inverted_bathymetry=True,
               verbose=verbose)
    
    if remove_nc_files:
        for file_out in files_out:
            os.remove(file_out)


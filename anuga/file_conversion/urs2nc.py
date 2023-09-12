
import os
from struct import pack, unpack
import array as p_array
import numpy as num

from anuga.anuga_exceptions import ANUGAError
from anuga.utilities.numerical_tools import ensure_numeric   
from anuga.caching.caching import myhash

from anuga.file.netcdf import Write_nc, write_elevation_nc

import anuga.utilities.log as log
from anuga.file.mux import WAVEHEIGHT_MUX_LABEL, EAST_VELOCITY_LABEL, \
                            NORTH_VELOCITY_LABEL


def urs2nc(basename_in='o', basename_out='urs'):
    """Convert the 3 urs files to 4 nc files.

    The name of the urs file names must be;
    [basename_in]-z-mux
    [basename_in]-e-mux
    [basename_in]-n-mux
    """

    files_in = [basename_in + WAVEHEIGHT_MUX_LABEL,
                basename_in + EAST_VELOCITY_LABEL,
                basename_in + NORTH_VELOCITY_LABEL]
    files_out = [basename_out + '_ha.nc',
                 basename_out + '_ua.nc',
                 basename_out + '_va.nc']
    quantities = ['HA', 'UA', 'VA']

    #if os.access(files_in[0]+'.mux', os.F_OK) == 0 :
    for i, file_name in enumerate(files_in):
        if os.access(file_name, os.F_OK) == 0:
            if os.access(file_name + '.mux', os.F_OK) == 0 :
                msg = 'File %s does not exist or is not accessible' % file_name
                raise IOError(msg)
            else:
               files_in[i] += '.mux'
               log.critical("file_name %s" % file_name)

    hashed_elevation = None
    for file_in, file_out, quantity in zip(files_in,
                                           files_out,
                                           quantities):
        lonlatdep, lon, lat, depth = _binary_c2nc(file_in,
                                                  file_out,
                                                  quantity)
        if hashed_elevation is None:
            elevation_file = basename_out + '_e.nc'
            write_elevation_nc(elevation_file,
                               lon,
                               lat,
                               depth)
            hashed_elevation = myhash(lonlatdep)
        else:
            msg = "The elevation information in the mux files is inconsistent"
            assert hashed_elevation == myhash(lonlatdep), msg

    files_out.append(elevation_file)

    return files_out


def _binary_c2nc(file_in, file_out, quantity):
    """Reads in a quantity urs file and writes a quantity nc file.
    Additionally, returns the depth and lat, long info,
    so it can be written to a file.
    """

    columns = 3                           # long, lat , depth
    mux_file = open(file_in, 'rb')

    # Number of points/stations
    (points_num,) = unpack('i', mux_file.read(4))

    # nt, int - Number of time steps
    (time_step_count,) = unpack('i', mux_file.read(4))

    #dt, float - time step, seconds
    (time_step,) = unpack('f', mux_file.read(4))

    msg = "Bad data in the mux file."
    if points_num < 0:
        mux_file.close()
        raise ANUGAError(msg)
    if time_step_count < 0:
        mux_file.close()
        raise ANUGAError(msg)
    if time_step < 0:
        mux_file.close()
        raise ANUGAError(msg)

    lonlatdep = p_array.array('f')
    lonlatdep.read(mux_file, columns * points_num)
    lonlatdep = num.array(lonlatdep, dtype=float)
    lonlatdep = num.reshape(lonlatdep, (points_num, columns))

    lon, lat, depth = lon_lat2grid(lonlatdep)
    lon_sorted = list(lon)
    lon_sorted.sort()

    if not num.all(lon == lon_sorted):
        msg = "Longitudes in mux file are not in ascending order"
        raise IOError(msg)

    lat_sorted = list(lat)
    lat_sorted.sort()

    nc_file = Write_nc(quantity,
                       file_out,
                       time_step_count,
                       time_step,
                       lon,
                       lat)

    for i in range(time_step_count):
        #Read in a time slice from mux file
        hz_p_array = p_array.array('f')
        hz_p_array.read(mux_file, points_num)
        hz_p = num.array(hz_p_array, dtype=float)
        hz_p = num.reshape(hz_p, (len(lon), len(lat)))
        hz_p = num.transpose(hz_p)  # mux has lat varying fastest, nc has long v.f.

        #write time slice to nc file
        nc_file.store_timestep(hz_p)

    mux_file.close()
    nc_file.close()

    return lonlatdep, lon, lat, depth



def lon_lat2grid(long_lat_dep):
    """
    given a list of points that are assumed to be an a grid,
    return the long's and lat's of the grid.
    long_lat_dep is an array where each row is a position.
    The first column is longitudes.
    The second column is latitudes.

    The latitude is the fastest varying dimension - in mux files
    """

    LONG = 0
    LAT = 1
    QUANTITY = 2

    long_lat_dep = ensure_numeric(long_lat_dep, float)

    num_points = long_lat_dep.shape[0]
    this_rows_long = long_lat_dep[0,LONG]

    # Count the length of unique latitudes
    i = 0
    while long_lat_dep[i,LONG] == this_rows_long and i < num_points:
        i += 1

    # determine the lats and longsfrom the grid
    lat = long_lat_dep[:i, LAT]
    long = long_lat_dep[::i, LONG]

    lenlong = len(long)
    lenlat = len(lat)

    msg = 'Input data is not gridded'
    assert num_points % lenlat == 0, msg
    assert num_points % lenlong == 0, msg

    # Test that data is gridded
    for i in range(lenlong):
        msg = 'Data is not gridded.  It must be for this operation'
        first = i * lenlat
        last = first + lenlat

        assert num.allclose(long_lat_dep[first:last,LAT], lat), msg
        assert num.allclose(long_lat_dep[first:last,LONG], long[i]), msg

    msg = 'Out of range latitudes/longitudes'
    for l in lat:assert -90 < l < 90 , msg
    for l in long:assert -180 < l < 180 , msg

    # Changing quantity from lat being the fastest varying dimension to
    # long being the fastest varying dimension
    # FIXME - make this faster/do this a better way
    # use numeric transpose, after reshaping the quantity vector
    quantity = num.zeros(num_points, float)

    for lat_i, _ in enumerate(lat):
        for long_i, _ in enumerate(long):
            q_index = lat_i*lenlong + long_i
            lld_index = long_i*lenlat + lat_i
            temp = long_lat_dep[lld_index, QUANTITY]
            quantity[q_index] = temp

    return long, lat, quantity





##
# @brief Convert 3 URS files back to 4 NC files.
# @param basename_in Stem of the input filenames.
# @param basename_outStem of the output filenames.
# @note The name of the urs file names must be:
#          [basename_in]-z-mux
#          [basename_in]-e-mux
#          [basename_in]-n-mux
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
                raise IOError, msg
            else:
               files_in[i] += '.mux'
               log.critical("file_name %s" % file_name)

    hashed_elevation = None
    for file_in, file_out, quantity in map(None, files_in,
                                           files_out,
                                           quantities):
        lonlatdep, lon, lat, depth = _binary_c2nc(file_in,
                                                  file_out,
                                                  quantity)
        if hashed_elevation == None:
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


##
# @brief Convert a quantity URS file to a NetCDF file.
# @param file_in Path to input URS file.
# @param file_out Path to the output file.
# @param quantity Name of the quantity to be written to the output file.
# @return A tuple (lonlatdep, lon, lat, depth).
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
        raise ANUGAError, msg
    if time_step_count < 0:
        mux_file.close()
        raise ANUGAError, msg
    if time_step < 0:
        mux_file.close()
        raise ANUGAError, msg

    lonlatdep = p_array.array('f')
    lonlatdep.read(mux_file, columns * points_num)
    lonlatdep = num.array(lonlatdep, dtype=num.float)
    lonlatdep = num.reshape(lonlatdep, (points_num, columns))

    lon, lat, depth = lon_lat2grid(lonlatdep)
    lon_sorted = list(lon)
    lon_sorted.sort()

    if not num.alltrue(lon == lon_sorted):
        msg = "Longitudes in mux file are not in ascending order"
        raise IOError, msg

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
        hz_p = num.array(hz_p_array, dtype=num.float)
        hz_p = num.reshape(hz_p, (len(lon), len(lat)))
        hz_p = num.transpose(hz_p)  # mux has lat varying fastest, nc has long v.f.

        #write time slice to nc file
        nc_file.store_timestep(hz_p)

    mux_file.close()
    nc_file.close()

    return lonlatdep, lon, lat, depth



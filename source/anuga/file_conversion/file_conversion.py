""" Conversion routines.
    ANUGA needs to deal with many different file formats, and this
    module provides routines for easily converting between them.
    
    These routines are necessarily high level, sitting above the various
    ANUGA modules. They take a file as input, and output a file.
"""

#non ANUGA imports
from Scientific.IO.NetCDF import NetCDFFile
import numpy as num
import os.path
import os

#ANUGA imports
from anuga.coordinate_transforms.geo_reference import Geo_reference, \
     write_NetCDF_georeference, ensure_geo_reference
from anuga.abstract_2d_finite_volumes.pmesh2domain import \
     pmesh_to_domain_instance
from anuga.utilities.numerical_tools import ensure_numeric
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a, \
                            netcdf_float


#shallow water imports
from anuga.file.sww import Read_sww, Write_sww
from anuga.shallow_water.shallow_water_domain import Domain
from anuga.shallow_water.shallow_water_domain import Domain


def sww2obj(basefilename, size):
    """ Convert netcdf based data output to obj

        Convert SWW data to OBJ data.
        basefilename Stem of filename, needs size and extension added.
        size The number of lines to write.
    """

    from Scientific.IO.NetCDF import NetCDFFile

    # Get NetCDF
    FN = create_filename('.', basefilename, 'sww', size)
    log.critical('Reading from %s' % FN)
    fid = NetCDFFile(FN, netcdf_mode_r)  #Open existing file for read

    # Get the variables
    x = fid.variables['x']
    y = fid.variables['y']
    z = fid.variables['elevation']
    time = fid.variables['time']
    stage = fid.variables['stage']

    M = size  #Number of lines
    xx = num.zeros((M,3), num.float)
    yy = num.zeros((M,3), num.float)
    zz = num.zeros((M,3), num.float)

    for i in range(M):
        for j in range(3):
            xx[i,j] = x[i+j*M]
            yy[i,j] = y[i+j*M]
            zz[i,j] = z[i+j*M]

    # Write obj for bathymetry
    FN = create_filename('.', basefilename, 'obj', size)
    write_obj(FN,xx,yy,zz)

    # Now read all the data with variable information, combine with
    # x,y info and store as obj
    for k in range(len(time)):
        t = time[k]
        log.critical('Processing timestep %f' % t)

        for i in range(M):
            for j in range(3):
                zz[i,j] = stage[k,i+j*M]

        #Write obj for variable data
        #FN = create_filename(basefilename, 'obj', size, time=t)
        FN = create_filename('.', basefilename[:5], 'obj', size, time=t)
        write_obj(FN, xx, yy, zz)


##
# @brief 
# @param basefilename Stem of filename, needs size and extension added.
def dat2obj(basefilename):
    """Convert line based data output to obj
    FIXME: Obsolete?
    """

    import glob, os
    from anuga.config import data_dir

    # Get bathymetry and x,y's
    lines = open(data_dir+os.sep+basefilename+'_geometry.dat', 'r').readlines()

    M = len(lines)  #Number of lines
    x = num.zeros((M,3), num.float)
    y = num.zeros((M,3), num.float)
    z = num.zeros((M,3), num.float)

    for i, line in enumerate(lines):
        tokens = line.split()
        values = map(float, tokens)

        for j in range(3):
            x[i,j] = values[j*3]
            y[i,j] = values[j*3+1]
            z[i,j] = values[j*3+2]

    # Write obj for bathymetry
    write_obj(data_dir + os.sep + basefilename + '_geometry', x, y, z)

    # Now read all the data files with variable information, combine with
    # x,y info and store as obj.

    files = glob.glob(data_dir + os.sep + basefilename + '*.dat')
    for filename in files:
        log.critical('Processing %s' % filename)

        lines = open(data_dir + os.sep + filename, 'r').readlines()
        assert len(lines) == M
        root, ext = os.path.splitext(filename)

        # Get time from filename
        i0 = filename.find('_time=')
        if i0 == -1:
            #Skip bathymetry file
            continue

        i0 += 6  #Position where time starts
        i1 = filename.find('.dat')

        if i1 > i0:
            t = float(filename[i0:i1])
        else:
            raise DataTimeError, 'Hmmmm'

        for i, line in enumerate(lines):
            tokens = line.split()
            values = map(float,tokens)

            for j in range(3):
                z[i,j] = values[j]

        # Write obj for variable data
        write_obj(data_dir + os.sep + basefilename + '_time=%.4f' % t, x, y, z)

##
# @brief Convert time-series text file to TMS file.
# @param filename 
# @param quantity_names 
# @param time_as_seconds 
def timefile2netcdf(filename, quantity_names=None, time_as_seconds=False):
    """Template for converting typical text files with time series to
    NetCDF tms file.

    The file format is assumed to be either two fields separated by a comma:

        time [DD/MM/YY hh:mm:ss], value0 value1 value2 ...

    E.g

      31/08/04 00:00:00, 1.328223 0 0
      31/08/04 00:15:00, 1.292912 0 0

    or time (seconds), value0 value1 value2 ...

      0.0, 1.328223 0 0
      0.1, 1.292912 0 0

    will provide a time dependent function f(t) with three attributes

    filename is assumed to be the rootname with extenisons .txt and .sww
    """

    import time, calendar
    from anuga.config import time_format
    from anuga.utilities.numerical_tools import ensure_numeric

    file_text = filename + '.txt'
    fid = open(file_text)
    line = fid.readline()
    fid.close()

    fields = line.split(',')
    msg = "File %s must have the format 'datetime, value0 value1 value2 ...'" \
          % file_text
    assert len(fields) == 2, msg

    if not time_as_seconds:
        try:
            starttime = calendar.timegm(time.strptime(fields[0], time_format))
        except ValueError:
            msg = 'First field in file %s must be' % file_text
            msg += ' date-time with format %s.\n' % time_format
            msg += 'I got %s instead.' % fields[0]
            raise DataTimeError, msg
    else:
        try:
            starttime = float(fields[0])
        except Error:
            msg = "Bad time format"
            raise DataTimeError, msg

    # Split values
    values = []
    for value in fields[1].split():
        values.append(float(value))

    q = ensure_numeric(values)

    msg = 'ERROR: File must contain at least one independent value'
    assert len(q.shape) == 1, msg

    # Read times proper
    from anuga.config import time_format
    import time, calendar

    fid = open(file_text)
    lines = fid.readlines()
    fid.close()

    N = len(lines)
    d = len(q)

    T = num.zeros(N, num.float)       # Time
    Q = num.zeros((N, d), num.float)  # Values

    for i, line in enumerate(lines):
        fields = line.split(',')
        if not time_as_seconds:
            realtime = calendar.timegm(time.strptime(fields[0], time_format))
        else:
            realtime = float(fields[0])
        T[i] = realtime - starttime

        for j, value in enumerate(fields[1].split()):
            Q[i, j] = float(value)

    msg = 'File %s must list time as a monotonuosly ' % filename
    msg += 'increasing sequence'
    assert num.alltrue(T[1:] - T[:-1] > 0), msg

    #Create NetCDF file
    from Scientific.IO.NetCDF import NetCDFFile

    fid = NetCDFFile(filename + '.tms', netcdf_mode_w)

    fid.institution = 'Geoscience Australia'
    fid.description = 'Time series'

    #Reference point
    #Start time in seconds since the epoch (midnight 1/1/1970)
    #FIXME: Use Georef
    fid.starttime = starttime

    # dimension definitions
    #fid.createDimension('number_of_volumes', self.number_of_volumes)
    #fid.createDimension('number_of_vertices', 3)

    fid.createDimension('number_of_timesteps', len(T))

    fid.createVariable('time', netcdf_float, ('number_of_timesteps',))

    fid.variables['time'][:] = T

    for i in range(Q.shape[1]):
        try:
            name = quantity_names[i]
        except:
            name = 'Attribute%d' % i

        fid.createVariable(name, netcdf_float, ('number_of_timesteps',))
        fid.variables[name][:] = Q[:,i]

    fid.close()
    
    

##
# @brief 
# @param filename
# @param verbose
def tsh2sww(filename, verbose=False):
    """
    to check if a tsh/msh file 'looks' good.
    """

    if verbose == True: log.critical('Creating domain from %s' % filename)

    domain = pmesh_to_domain_instance(filename, Domain)

    if verbose == True: log.critical("Number of triangles = %s" % len(domain))

    domain.smooth = True
    domain.format = 'sww'   #Native netcdf visualisation format
    file_path, filename = path.split(filename)
    filename, ext = path.splitext(filename)
    domain.set_name(filename)
    domain.reduction = mean

    if verbose == True: log.critical("file_path = %s" % file_path)

    if file_path == "":
        file_path = "."
    domain.set_datadir(file_path)

    if verbose == True:
        log.critical("Output written to %s%s%s.%s"
                     % (domain.get_datadir(), sep, domain.get_name(),
                        domain.format))

    sww = SWW_file(domain)
    sww.store_connectivity()
    sww.store_timestep()


##
# @brief Convert URS file to SWW file.
# @param basename_in Stem of the input filename.
# @param basename_out Stem of the output filename.
# @param verbose True if this function is to be verbose.
# @param remove_nc_files 
# @param minlat Sets extent of area to be used.  If not supplied, full extent.
# @param maxlat Sets extent of area to be used.  If not supplied, full extent.
# @param minlon Sets extent of area to be used.  If not supplied, full extent.
# @param maxlon Sets extent of area to be used.  If not supplied, full extent.
# @param mint 
# @param maxt 
# @param mean_stage 
# @param origin A 3-tuple with geo referenced UTM coordinates
# @param zscale 
# @param fail_on_NaN 
# @param NaN_filler 
# @param elevation 
# @note Also convert latitude and longitude to UTM. All coordinates are
#       assumed to be given in the GDA94 datum.
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

    if basename_out == None:
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



##
# @brief Return max&min indexes (for slicing) of area specified.
# @param latitudes_ref ??
# @param longitudes_ref ??
# @param minlat Minimum latitude of specified area.
# @param maxlat Maximum latitude of specified area.
# @param minlon Minimum longitude of specified area.
# @param maxlon Maximum longitude of specified area.
# @return Tuple (lat_min_index, lat_max_index, lon_min_index, lon_max_index)
def get_min_max_indices(latitudes_ref, longitudes_ref,
                         minlat=None, maxlat=None,
                         minlon=None, maxlon=None):
    """
    Return max, min indexes (for slicing) of the lat and long arrays to cover
    the area specified with min/max lat/long.

    Think of the latitudes and longitudes describing a 2d surface.
    The area returned is, if possible, just big enough to cover the
    inputed max/min area. (This will not be possible if the max/min area
    has a section outside of the latitudes/longitudes area.)

    asset  longitudes are sorted,
    long - from low to high (west to east, eg 148 - 151)
    assert latitudes are sorted, ascending or decending
    """

    latitudes = latitudes_ref[:]
    longitudes = longitudes_ref[:]

    latitudes = ensure_numeric(latitudes)
    longitudes = ensure_numeric(longitudes)

    assert num.allclose(num.sort(longitudes), longitudes)

    #print latitudes[0],longitudes[0]
    #print len(latitudes),len(longitudes)
    #print latitudes[len(latitudes)-1],longitudes[len(longitudes)-1]

    lat_ascending = True
    if not num.allclose(num.sort(latitudes), latitudes):
        lat_ascending = False
        # reverse order of lat, so it's in ascending order
        latitudes = latitudes[::-1]
        assert num.allclose(num.sort(latitudes), latitudes)

    largest_lat_index = len(latitudes)-1

    #Cut out a smaller extent.
    if minlat == None:
        lat_min_index = 0
    else:
        lat_min_index = num.searchsorted(latitudes, minlat)-1
        if lat_min_index < 0:
            lat_min_index = 0

    if maxlat == None:
        lat_max_index = largest_lat_index #len(latitudes)
    else:
        lat_max_index = num.searchsorted(latitudes, maxlat)
        if lat_max_index > largest_lat_index:
            lat_max_index = largest_lat_index

    if minlon == None:
        lon_min_index = 0
    else:
        lon_min_index = num.searchsorted(longitudes, minlon)-1
        if lon_min_index < 0:
            lon_min_index = 0

    if maxlon == None:
        lon_max_index = len(longitudes)
    else:
        lon_max_index = num.searchsorted(longitudes, maxlon)

    # Reversing the indexes, if the lat array is decending
    if lat_ascending is False:
        lat_min_index, lat_max_index = largest_lat_index - lat_max_index, \
                                       largest_lat_index - lat_min_index
    lat_max_index = lat_max_index + 1 # taking into account how slicing works
    lon_max_index = lon_max_index + 1 # taking into account how slicing works

    return lat_min_index, lat_max_index, lon_min_index, lon_max_index


##
# @brief 
# @param filename 
# @param x 
# @param y 
# @param z 
def write_obj(filename, x, y, z):
    """Store x,y,z vectors into filename (obj format).

       Vectors are assumed to have dimension (M,3) where
       M corresponds to the number elements.
       triangles are assumed to be disconnected

       The three numbers in each vector correspond to three vertices,

       e.g. the x coordinate of vertex 1 of element i is in x[i,1]
    """

    import os.path

    root, ext = os.path.splitext(filename)
    if ext == '.obj':
        FN = filename
    else:
        FN = filename + '.obj'

    outfile = open(FN, 'wb')
    outfile.write("# Triangulation as an obj file\n")

    M, N = x.shape
    assert N == 3  #Assuming three vertices per element

    for i in range(M):
        for j in range(N):
            outfile.write("v %f %f %f\n" % (x[i,j], y[i,j], z[i,j]))

    for i in range(M):
        base = i * N
        outfile.write("f %d %d %d\n" % (base+1, base+2, base+3))

    outfile.close()


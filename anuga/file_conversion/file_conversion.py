""" Conversion routines.
    ANUGA needs to deal with many different file formats, and this
    module provides routines for easily converting between them.

    These routines are necessarily high level, sitting above the various
    ANUGA modules. They take a file as input, and output a file.
"""

#non ANUGA imports

from anuga.file.netcdf import NetCDFFile
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

from anuga.anuga_exceptions import *


#shallow water imports
from anuga.file.sww import Read_sww, Write_sww
from anuga.shallow_water.shallow_water_domain import Domain
from anuga.shallow_water.shallow_water_domain import Domain


def sww2obj(filename, size):
    """ Convert netcdf based data output to obj

        Convert SWW data to OBJ data.
        basefilename Stem of filename, needs size and extension added.
        size The number of lines to write.
    """

    if filename[-4:] != '.sww':
        raise IOError('Output file %s should be of type .sww.' % sww_file)

    basefilename = filename[:-4]

    # Get NetCDF
    nc_fname = create_filename('.', basefilename, 'sww', size)
    log.critical('Reading from %s' % nc_fname)
    fid = NetCDFFile(nc_fname, netcdf_mode_r)  #Open existing file for read

    # Get the variables
    x = fid.variables['x']
    y = fid.variables['y']
    z = fid.variables['elevation']
    time = fid.variables['time']
    stage = fid.variables['stage']

    M = size  #Number of lines
    xx = num.zeros((M,3), float)
    yy = num.zeros((M,3), float)
    zz = num.zeros((M,3), float)

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

def timefile2netcdf(file_text, file_out = None, quantity_names=None, \
                                time_as_seconds=False):
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

    filename is assumed to be the rootname with extensions .txt/.tms and .sww
    """

    import time, calendar
    from anuga.config import time_format
    from anuga.utilities.numerical_tools import ensure_numeric

    if file_text[-4:] != '.txt':
        raise IOError('Input file %s should be of type .txt.' % file_text)

    if file_out is None:
        file_out = file_text[:-4] + '.tms'

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
            raise DataTimeError(msg)
    else:
        try:
            starttime = float(fields[0])
        except Exception:
            msg = "Bad time format"
            raise DataTimeError(msg)

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

    T = num.zeros(N, float)       # Time
    Q = num.zeros((N, d), float)  # Values

    for i, line in enumerate(lines):
        fields = line.split(',')
        if not time_as_seconds:
            realtime = calendar.timegm(time.strptime(fields[0], time_format))
        else:
            realtime = float(fields[0])
        T[i] = realtime - starttime

        for j, value in enumerate(fields[1].split()):
            Q[i, j] = float(value)

    msg = 'File %s must list time as a monotonuosly ' % file_text
    msg += 'increasing sequence'
    assert num.alltrue(T[1:] - T[:-1] > 0), msg

    #Create NetCDF file
    fid = NetCDFFile(file_out, netcdf_mode_w)

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



def tsh2sww(filename, verbose=False):
    """
    to check if a tsh/msh file 'looks' good.
    """

    if filename[-4:] != '.tsh' and filename[-4:] != '.msh':
        raise IOError('Input file %s should be .tsh or .msh.' % name_out)

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
    if minlat is None:
        lat_min_index = 0
    else:
        lat_min_index = num.searchsorted(latitudes, minlat)-1
        if lat_min_index < 0:
            lat_min_index = 0

    if maxlat is None:
        lat_max_index = largest_lat_index #len(latitudes)
    else:
        lat_max_index = num.searchsorted(latitudes, maxlat)
        if lat_max_index > largest_lat_index:
            lat_max_index = largest_lat_index

    if minlon is None:
        lon_min_index = 0
    else:
        lon_min_index = num.searchsorted(longitudes, minlon)-1
        if lon_min_index < 0:
            lon_min_index = 0

    if maxlon is None:
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

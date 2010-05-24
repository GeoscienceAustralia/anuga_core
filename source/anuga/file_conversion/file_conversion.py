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
from anuga.shallow_water.sww_file import Read_sww, Write_sww
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
# @brief Convert CSIRO ESRI file to an SWW boundary file.
# @param bath_dir 
# @param elevation_dir 
# @param ucur_dir 
# @param vcur_dir 
# @param sww_file 
# @param minlat 
# @param maxlat 
# @param minlon 
# @param maxlon 
# @param zscale 
# @param mean_stage 
# @param fail_on_NaN 
# @param elevation_NaN_filler 
# @param bath_prefix 
# @param elevation_prefix 
# @param verbose 
# @note Also convert latitude and longitude to UTM. All coordinates are
#       assumed to be given in the GDA94 datum.
def asc_csiro2sww(bath_dir,
                  elevation_dir,
                  ucur_dir,
                  vcur_dir,
                  sww_file,
                  minlat=None, maxlat=None,
                  minlon=None, maxlon=None,
                  zscale=1,
                  mean_stage=0,
                  fail_on_NaN=True,
                  elevation_NaN_filler=0,
                  bath_prefix='ba',
                  elevation_prefix='el',
                  verbose=False):
    """
    Produce an sww boundary file, from esri ascii data from CSIRO.

    Also convert latitude and longitude to UTM. All coordinates are
    assumed to be given in the GDA94 datum.

    assume:
    All files are in esri ascii format

    4 types of information
    bathymetry
    elevation
    u velocity
    v velocity

    Assumptions
    The metadata of all the files is the same
    Each type is in a seperate directory
    One bath file with extention .000
    The time period is less than 24hrs and uniform.
    """

    from Scientific.IO.NetCDF import NetCDFFile

    from anuga.coordinate_transforms.redfearn import redfearn

    # So if we want to change the precision it's done here
    precision = netcdf_float 

    # go in to the bath dir and load the only file,
    bath_files = os.listdir(bath_dir)
    bath_file = bath_files[0]
    bath_dir_file =  bath_dir + os.sep + bath_file
    bath_metadata, bath_grid =  _read_asc(bath_dir_file)

    #Use the date.time of the bath file as a basis for
    #the start time for other files
    base_start = bath_file[-12:]

    #go into the elevation dir and load the 000 file
    elevation_dir_file = elevation_dir  + os.sep + elevation_prefix \
                         + base_start

    elevation_files = os.listdir(elevation_dir)
    ucur_files = os.listdir(ucur_dir)
    vcur_files = os.listdir(vcur_dir)
    elevation_files.sort()

    # the first elevation file should be the
    # file with the same base name as the bath data
    assert elevation_files[0] == 'el' + base_start

    number_of_latitudes = bath_grid.shape[0]
    number_of_longitudes = bath_grid.shape[1]
    number_of_volumes = (number_of_latitudes-1) * (number_of_longitudes-1) * 2

    longitudes = [bath_metadata['xllcorner'] + x*bath_metadata['cellsize'] \
                  for x in range(number_of_longitudes)]
    latitudes = [bath_metadata['yllcorner'] + y*bath_metadata['cellsize'] \
                 for y in range(number_of_latitudes)]

     # reverse order of lat, so the first lat represents the first grid row
    latitudes.reverse()

    kmin, kmax, lmin, lmax = get_min_max_indices(latitudes[:],longitudes[:],
                                                  minlat=minlat, maxlat=maxlat,
                                                  minlon=minlon, maxlon=maxlon)

    bath_grid = bath_grid[kmin:kmax,lmin:lmax]
    latitudes = latitudes[kmin:kmax]
    longitudes = longitudes[lmin:lmax]
    number_of_latitudes = len(latitudes)
    number_of_longitudes = len(longitudes)
    number_of_times = len(os.listdir(elevation_dir))
    number_of_points = number_of_latitudes * number_of_longitudes
    number_of_volumes = (number_of_latitudes-1) * (number_of_longitudes-1) * 2

    #Work out the times
    if len(elevation_files) > 1:
        # Assume: The time period is less than 24hrs.
        time_period = (int(elevation_files[1][-3:]) \
                       - int(elevation_files[0][-3:])) * 60*60
        times = [x*time_period for x in range(len(elevation_files))]
    else:
        times = [0.0]

    if verbose:
        log.critical('------------------------------------------------')
        log.critical('Statistics:')
        log.critical('  Extent (lat/lon):')
        log.critical('    lat in [%f, %f], len(lat) == %d'
                     % (min(latitudes), max(latitudes), len(latitudes)))
        log.critical('    lon in [%f, %f], len(lon) == %d'
                     % (min(longitudes), max(longitudes), len(longitudes)))
        log.critical('    t in [%f, %f], len(t) == %d'
                     % (min(times), max(times), len(times)))

    ######### WRITE THE SWW FILE #############

    # NetCDF file definition
    outfile = NetCDFFile(sww_file, netcdf_mode_w)

    #Create new file
    outfile.institution = 'Geoscience Australia'
    outfile.description = 'Converted from XXX'

    #For sww compatibility
    outfile.smoothing = 'Yes'
    outfile.order = 1

    #Start time in seconds since the epoch (midnight 1/1/1970)
    outfile.starttime = starttime = times[0]

    # dimension definitions
    outfile.createDimension('number_of_volumes', number_of_volumes)
    outfile.createDimension('number_of_vertices', 3)
    outfile.createDimension('number_of_points', number_of_points)
    outfile.createDimension('number_of_timesteps', number_of_times)

    # variable definitions
    outfile.createVariable('x', precision, ('number_of_points',))
    outfile.createVariable('y', precision, ('number_of_points',))
    outfile.createVariable('elevation', precision, ('number_of_points',))

    #FIXME: Backwards compatibility
    #outfile.createVariable('z', precision, ('number_of_points',))
    #################################

    outfile.createVariable('volumes', netcdf_int, ('number_of_volumes',
                                                   'number_of_vertices'))

    outfile.createVariable('time', precision, ('number_of_timesteps',))

    outfile.createVariable('stage', precision, ('number_of_timesteps',
                                                'number_of_points'))

    outfile.createVariable('xmomentum', precision, ('number_of_timesteps',
                                                    'number_of_points'))

    outfile.createVariable('ymomentum', precision, ('number_of_timesteps',
                                                    'number_of_points'))

    #Store
    from anuga.coordinate_transforms.redfearn import redfearn

    x = num.zeros(number_of_points, num.float)  #Easting
    y = num.zeros(number_of_points, num.float)  #Northing

    if verbose: log.critical('Making triangular grid')

    #Get zone of 1st point.
    refzone, _, _ = redfearn(latitudes[0], longitudes[0])

    vertices = {}
    i = 0
    for k, lat in enumerate(latitudes):
        for l, lon in enumerate(longitudes):
            vertices[l,k] = i

            zone, easting, northing = redfearn(lat, lon)

            #msg = 'Zone boundary crossed at longitude =', lon
            #assert zone == refzone, msg
            #print '%7.2f %7.2f %8.2f %8.2f' %(lon, lat, easting, northing)
            x[i] = easting
            y[i] = northing
            i += 1

    #Construct 2 triangles per 'rectangular' element
    volumes = []
    for l in range(number_of_longitudes-1):    #X direction
        for k in range(number_of_latitudes-1): #Y direction
            v1 = vertices[l,k+1]
            v2 = vertices[l,k]
            v3 = vertices[l+1,k+1]
            v4 = vertices[l+1,k]

            #Note, this is different to the ferrit2sww code
            #since the order of the lats is reversed.
            volumes.append([v1,v3,v2]) #Upper element
            volumes.append([v4,v2,v3]) #Lower element

    volumes = num.array(volumes, num.int)      #array default#

    geo_ref = Geo_reference(refzone, min(x), min(y))
    geo_ref.write_NetCDF(outfile)

    # This will put the geo ref in the middle
    #geo_ref = Geo_reference(refzone, (max(x)+min(x))/2., (max(x)+min(y))/2.)

    if verbose:
        log.critical('------------------------------------------------')
        log.critical('More Statistics:')
        log.critical('  Extent (/lon):')
        log.critical('    x in [%f, %f], len(lat) == %d'
                     % (min(x), max(x), len(x)))
        log.critical('    y in [%f, %f], len(lon) == %d'
                     % (min(y), max(y), len(y)))
        log.critical('geo_ref: ', geo_ref)

    z = num.resize(bath_grid,outfile.variables['elevation'][:].shape)
    outfile.variables['x'][:] = x - geo_ref.get_xllcorner()
    outfile.variables['y'][:] = y - geo_ref.get_yllcorner()
    # FIXME (Ole): Remove once viewer has been recompiled and changed
    #              to use elevation instead of z
    #outfile.variables['z'][:] = z
    outfile.variables['elevation'][:] = z
    outfile.variables['volumes'][:] = volumes.astype(num.int32) # On Opteron 64

    stage = outfile.variables['stage']
    xmomentum = outfile.variables['xmomentum']
    ymomentum = outfile.variables['ymomentum']

    outfile.variables['time'][:] = times   #Store time relative

    if verbose: log.critical('Converting quantities')

    n = number_of_times
    for j in range(number_of_times):
        # load in files
        elevation_meta, elevation_grid = \
            _read_asc(elevation_dir + os.sep + elevation_files[j])

        _, u_momentum_grid = _read_asc(ucur_dir + os.sep + ucur_files[j])
        _, v_momentum_grid = _read_asc(vcur_dir + os.sep + vcur_files[j])

        #cut matrix to desired size
        elevation_grid = elevation_grid[kmin:kmax,lmin:lmax]
        u_momentum_grid = u_momentum_grid[kmin:kmax,lmin:lmax]
        v_momentum_grid = v_momentum_grid[kmin:kmax,lmin:lmax]

        # handle missing values
        missing = (elevation_grid == elevation_meta['NODATA_value'])
        if num.sometrue (missing):
            if fail_on_NaN:
                msg = 'File %s contains missing values' \
                      % (elevation_files[j])
                raise DataMissingValuesError, msg
            else:
                elevation_grid = elevation_grid*(missing==0) \
                                 + missing*elevation_NaN_filler

        if verbose and j % ((n+10)/10) == 0: log.critical('  Doing %d of %d'
                                                          % (j, n))

        i = 0
        for k in range(number_of_latitudes):      #Y direction
            for l in range(number_of_longitudes): #X direction
                w = zscale*elevation_grid[k,l] + mean_stage
                stage[j,i] = w
                h = w - z[i]
                xmomentum[j,i] = u_momentum_grid[k,l]*h
                ymomentum[j,i] = v_momentum_grid[k,l]*h
                i += 1

    outfile.close()



def _read_asc(filename, verbose=False):
    """Read esri file from the following ASCII format (.asc)

    Example:

    ncols         3121
    nrows         1800
    xllcorner     722000
    yllcorner     5893000
    cellsize      25
    NODATA_value  -9999
    138.3698 137.4194 136.5062 135.5558 ..........

    filename Path to the file to read.
    verbose True if this function is to be verbose.
    """

    datafile = open(filename)

    if verbose: log.critical('Reading DEM from %s' % filename)

    lines = datafile.readlines()
    datafile.close()

    if verbose: log.critical('Got %d lines' % len(lines))

    ncols = int(lines.pop(0).split()[1].strip())
    nrows = int(lines.pop(0).split()[1].strip())
    xllcorner = float(lines.pop(0).split()[1].strip())
    yllcorner = float(lines.pop(0).split()[1].strip())
    cellsize = float(lines.pop(0).split()[1].strip())
    NODATA_value = float(lines.pop(0).split()[1].strip())

    assert len(lines) == nrows

    #Store data
    grid = []

    n = len(lines)
    for i, line in enumerate(lines):
        cells = line.split()
        assert len(cells) == ncols
        grid.append(num.array([float(x) for x in cells]))
    grid = num.array(grid)

    return {'xllcorner':xllcorner,
            'yllcorner':yllcorner,
            'cellsize':cellsize,
            'NODATA_value':NODATA_value}, grid


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



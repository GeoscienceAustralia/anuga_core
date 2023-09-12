"""
    Convert a ferret file to an SWW file.
"""


import numpy as num


# ANUGA modules
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_float
from anuga.file.sww import Write_sww  
from anuga.coordinate_transforms.geo_reference import Geo_reference, \
     write_NetCDF_georeference                          
import anuga.utilities.log as log

#local modules
from anuga.file_conversion.file_conversion import get_min_max_indices                            

class DataMissingValuesError(Exception):
    pass


def ferret2sww(basename_in, name_out=None,
               verbose=False,
               minlat=None, maxlat=None,
               minlon=None, maxlon=None,
               mint=None, maxt=None, mean_stage=0,
               origin=None, zscale=1,
               fail_on_NaN=True,
               NaN_filler=0,
               elevation=None,
               inverted_bathymetry=True
               ): #FIXME: Bathymetry should be obtained
                                  #from MOST somehow.
                                  #Alternatively from elsewhere
                                  #or, as a last resort,
                                  #specified here.
                                  #The value of -100 will work
                                  #for the Wollongong tsunami
                                  #scenario but is very hacky
    """Convert MOST and 'Ferret' NetCDF format for wave propagation to
    sww format native to abstract_2d_finite_volumes.

    Specify only basename_in and read files of the form
    basefilename_ha.nc, basefilename_ua.nc, basefilename_va.nc containing
    relative height, x-velocity and y-velocity, respectively.

    Also convert latitude and longitude to UTM. All coordinates are
    assumed to be given in the GDA94 datum.

    min's and max's: If omitted - full extend is used.
    To include a value min may equal it, while max must exceed it.
    Lat and lon are assuemd to be in decimal degrees

    origin is a 3-tuple with geo referenced
    UTM coordinates (zone, easting, northing)

    nc format has values organised as HA[TIME, LATITUDE, LONGITUDE]
    which means that longitude is the fastest
    varying dimension (row major order, so to speak)

    ferret2sww uses grid points as vertices in a triangular grid
    counting vertices from lower left corner upwards, then right
    """

    from anuga.file.netcdf import NetCDFFile

    _assert_lat_long(minlat, maxlat, minlon, maxlon)

    if name_out != None and name_out[-4:] != '.sww':
        raise IOError('Output file %s should be of type .sww.' % name_out)

    # Get NetCDF data
    if verbose: log.critical('Reading files %s_*.nc' % basename_in)

    # Wave amplitude (cm)
    file_h = NetCDFFile(basename_in + '_ha.nc', netcdf_mode_r) 
    
    # Velocity (x) (cm/s)
    file_u = NetCDFFile(basename_in + '_ua.nc', netcdf_mode_r)
     
    # Velocity (y) (cm/s)
    file_v = NetCDFFile(basename_in + '_va.nc', netcdf_mode_r)
    
    # Elevation (z) (m)
    file_e = NetCDFFile(basename_in + '_e.nc', netcdf_mode_r)  

    if name_out is None:
        swwname = basename_in + '.sww'
    else:
        swwname = name_out

    # Get dimensions of file_h
    for dimension in list(file_h.dimensions.keys()):
        if dimension[:3] == 'LON':
            dim_h_longitude = dimension
        if dimension[:3] == 'LAT':
            dim_h_latitude = dimension
        if dimension[:4] == 'TIME':
            dim_h_time = dimension

    times = file_h.variables[dim_h_time]
    latitudes = file_h.variables[dim_h_latitude]
    longitudes = file_h.variables[dim_h_longitude]

    kmin, kmax, lmin, lmax = get_min_max_indices(latitudes[:],
                                                  longitudes[:],
                                                  minlat, maxlat,
                                                  minlon, maxlon)
    # get dimensions for file_e
    for dimension in list(file_e.dimensions.keys()):
        if dimension[:3] == 'LON':
            dim_e_longitude = dimension
        if dimension[:3] == 'LAT':
            dim_e_latitude = dimension

    # get dimensions for file_u
    for dimension in list(file_u.dimensions.keys()):
        if dimension[:3] == 'LON':
            dim_u_longitude = dimension
        if dimension[:3] == 'LAT':
            dim_u_latitude = dimension

    # get dimensions for file_v
    for dimension in list(file_v.dimensions.keys()):
        if dimension[:3] == 'LON':
            dim_v_longitude = dimension
        if dimension[:3] == 'LAT':
            dim_v_latitude = dimension

    # Precision used by most for lat/lon is 4 or 5 decimals
    e_lat = num.around(file_e.variables[dim_e_latitude][:], 5)
    e_lon = num.around(file_e.variables[dim_e_longitude][:], 5)

    # Check that files are compatible
    assert num.allclose(latitudes, file_u.variables[dim_u_latitude])
    assert num.allclose(latitudes, file_v.variables[dim_v_latitude])
    assert num.allclose(latitudes, e_lat)

    assert num.allclose(longitudes, file_u.variables[dim_u_longitude])
    assert num.allclose(longitudes, file_v.variables[dim_v_longitude])
    assert num.allclose(longitudes, e_lon)

    if mint is None:
        jmin = 0
        mint = times[0]
    else:
        jmin = num.searchsorted(times, mint)
        
        # numpy.int32 didn't work in slicing of amplitude below
        jmin = int(jmin)

    if maxt is None:
        jmax = len(times)
        maxt = times[-1]
    else:
        jmax = num.searchsorted(times, maxt)
        
        # numpy.int32 didn't work in slicing of amplitude below
        jmax = int(jmax)        

    kmin, kmax, lmin, lmax = get_min_max_indices(latitudes[:],
                                                  longitudes[:],
                                                  minlat, maxlat,
                                                  minlon, maxlon)


    times = times[jmin:jmax]
    latitudes = latitudes[kmin:kmax]
    longitudes = longitudes[lmin:lmax]

    if verbose: log.critical('cropping')

    zname = 'ELEVATION'

    amplitudes = file_h.variables['HA'][jmin:jmax, kmin:kmax, lmin:lmax]
    uspeed = file_u.variables['UA'][jmin:jmax, kmin:kmax, lmin:lmax] #Lon
    vspeed = file_v.variables['VA'][jmin:jmax, kmin:kmax, lmin:lmax] #Lat
    elevations = file_e.variables[zname][kmin:kmax, lmin:lmax]

    # Get missing values
    nan_ha = file_h.variables['HA'].missing_value
    nan_ua = file_u.variables['UA'].missing_value
    nan_va = file_v.variables['VA'].missing_value
    if hasattr(file_e.variables[zname],'missing_value'):
        nan_e  = file_e.variables[zname].missing_value
    else:
        nan_e = None

    # Cleanup
    missing = (amplitudes == nan_ha)
    if num.any (missing):
        if fail_on_NaN:
            msg = 'NetCDFFile %s contains missing values' \
                  % basename_in + '_ha.nc'
            raise DataMissingValuesError(msg)
        else:
            amplitudes = amplitudes*(missing==0) + missing*NaN_filler

    missing = (uspeed == nan_ua)
    if num.any (missing):
        if fail_on_NaN:
            msg = 'NetCDFFile %s contains missing values' \
                  % basename_in + '_ua.nc'
            raise DataMissingValuesError(msg)
        else:
            uspeed = uspeed*(missing==0) + missing*NaN_filler

    missing = (vspeed == nan_va)
    if num.any (missing):
        if fail_on_NaN:
            msg = 'NetCDFFile %s contains missing values' \
                  % basename_in + '_va.nc'
            raise DataMissingValuesError(msg)
        else:
            vspeed = vspeed*(missing==0) + missing*NaN_filler

    missing = (elevations == nan_e)
    if num.any (missing):
        if fail_on_NaN:
            msg = 'NetCDFFile %s contains missing values' \
                  % basename_in + '_e.nc'
            raise DataMissingValuesError(msg)
        else:
            elevations = elevations*(missing==0) + missing*NaN_filler

    number_of_times = times.shape[0]
    number_of_latitudes = latitudes.shape[0]
    number_of_longitudes = longitudes.shape[0]

    assert amplitudes.shape[0] == number_of_times
    assert amplitudes.shape[1] == number_of_latitudes
    assert amplitudes.shape[2] == number_of_longitudes

    if verbose:
        _show_stats((latitudes, longitudes), times, amplitudes, \
                    (uspeed, vspeed), elevations)

    # print number_of_latitudes, number_of_longitudes
    number_of_points = number_of_latitudes * number_of_longitudes
    number_of_volumes = (number_of_latitudes-1) * (number_of_longitudes-1) * 2

    file_h.close()
    file_u.close()
    file_v.close()
    file_e.close()

    # NetCDF file definition
    outfile = NetCDFFile(swwname, netcdf_mode_w)

    description = 'Converted from Ferret files: %s, %s, %s, %s' \
                  % (basename_in + '_ha.nc',
                     basename_in + '_ua.nc',
                     basename_in + '_va.nc',
                     basename_in + '_e.nc')

    # Create new file
    starttime = times[0]

    sww = Write_sww(['elevation'], ['stage', 'xmomentum', 'ymomentum'])
    sww.store_header(outfile, times, number_of_volumes,
                     number_of_points, description=description,
                     verbose=verbose, sww_precision=netcdf_float)

    # Store
    from anuga.coordinate_transforms.redfearn import redfearn
    x = num.zeros(number_of_points, float)  #Easting
    y = num.zeros(number_of_points, float)  #Northing

    if verbose:
        log.critical('Making triangular grid')

    # Check zone boundaries
    refzone, _, _ = redfearn(latitudes[0], longitudes[0])

    vertices = {}
    i = 0
    for k, lat in enumerate(latitudes):       # Y direction
        for l, lon in enumerate(longitudes):  # X direction
            vertices[l, k] = i

            _, easting, northing = redfearn(lat, lon)

            #msg = 'Zone boundary crossed at longitude =', lon
            #assert zone == refzone, msg
            #print '%7.2f %7.2f %8.2f %8.2f' %(lon, lat, easting, northing)
            x[i] = easting
            y[i] = northing
            i += 1

    #Construct 2 triangles per 'rectangular' element
    volumes = []
    for l in range(number_of_longitudes-1):    # X direction
        for k in range(number_of_latitudes-1): # Y direction
            v1 = vertices[l, k+1]
            v2 = vertices[l, k]
            v3 = vertices[l+1, k+1]
            v4 = vertices[l+1, k]

            volumes.append([v1, v2, v3]) #Upper element
            volumes.append([v4, v3, v2]) #Lower element

    volumes = num.array(volumes, int)      #array default#

    if origin is None:
        origin = Geo_reference(refzone, min(x), min(y))
    geo_ref = write_NetCDF_georeference(origin, outfile)

    if elevation is not None:
        z = elevation
    else:
        if inverted_bathymetry:
            z = -1 * elevations
        else:
            z = elevations
    #FIXME: z should be obtained from MOST and passed in here

    #FIXME use the Write_sww instance(sww) to write this info
    z = num.resize(z, outfile.variables['elevation'][:].shape)
    outfile.variables['x'][:] = x - geo_ref.get_xllcorner()
    outfile.variables['y'][:] = y - geo_ref.get_yllcorner()
    #outfile.variables['z'][:] = z             #FIXME HACK for bacwards compat.
    outfile.variables['elevation'][:] = z
    outfile.variables['volumes'][:] = volumes.astype(num.int32) #For Opteron 64

    #Time stepping
    stage = outfile.variables['stage']
    xmomentum = outfile.variables['xmomentum']
    ymomentum = outfile.variables['ymomentum']

    if verbose:
        log.critical('Converting quantities')

    n = len(times)
    for j in range(n):
        if verbose and j % ((n+10)//10) == 0:
            log.critical('  Doing %d of %d' % (j, n))

        i = 0
        for k in range(number_of_latitudes):      # Y direction
            for l in range(number_of_longitudes): # X direction
                w = zscale * amplitudes[j, k, l]/100.0 + mean_stage
                stage[j, i] = w
                h = w - z[i]
                xmomentum[j, i] = uspeed[j, k, l]/100.0*h
                ymomentum[j, i] = vspeed[j, k, l]/100.0*h
                i += 1

    #outfile.close()

    #FIXME: Refactor using code from file_function.statistics
    #Something like print swwstats(swwname)
    if verbose:
        time_info = times, starttime, mint, maxt
        _show_sww_stats(outfile, swwname, geo_ref, time_info)

    outfile.close()


def _show_stats(latlong, times, amplitudes, speeds, elevations):
    """ Print the statistics nicely to the log file """
    
    latitudes, longitudes = latlong
    uspeed, vspeed = speeds
    
    log.critical('------------------------------------------------')
    log.critical('Statistics:')
    log.critical('  Extent (lat/lon):')
    log.critical('    lat in [%f, %f], len(lat) == %d'
                 % (num.min(latitudes), num.max(latitudes),
                    len(latitudes.flat)))
    log.critical('    lon in [%f, %f], len(lon) == %d'
                 % (num.min(longitudes), num.max(longitudes),
                    len(longitudes.flat)))
    log.critical('    t in [%f, %f], len(t) == %d'
                 % (num.min(times), num.max(times), len(times.flat)))

    name = 'Amplitudes (ha) [cm]'
    log.critical('  %s in [%f, %f]'
                 % (name, num.min(amplitudes), num.max(amplitudes)))

    name = 'Speeds (ua) [cm/s]'
    log.critical('  %s in [%f, %f]'
                 % (name, num.min(uspeed), num.max(uspeed)))

    name = 'Speeds (va) [cm/s]'
    log.critical('  %s in [%f, %f]'
                 % (name, num.min(vspeed), num.max(vspeed)))

    name = 'Elevations (e) [m]'
    log.critical('  %s in [%f, %f]'
                 % (name, num.min(elevations), num.max(elevations)))
                 
             

def _show_sww_stats(outfile, swwname, geo_ref, time_info):         
    """ Log SWW output stats. """
    times, starttime, mint, maxt = time_info
    x = outfile.variables['x'][:]
    y = outfile.variables['y'][:]
    log.critical('------------------------------------------------')
    log.critical('Statistics of output file:')
    log.critical('  Name: %s' %swwname)
    log.critical('  Reference:')
    log.critical('    Lower left corner: [%f, %f]'
                 % (geo_ref.get_xllcorner(), geo_ref.get_yllcorner()))
    log.critical(' Start time: %f' %starttime)
    log.critical('    Min time: %f' %mint)
    log.critical('    Max time: %f' %maxt)
    log.critical('  Extent:')
    log.critical('    x [m] in [%f, %f], len(x) == %d'
                 % (num.min(x), num.max(x), len(x.flat)))
    log.critical('    y [m] in [%f, %f], len(y) == %d'
                 % (num.min(y), num.max(y), len(y.flat)))
    log.critical('    t [s] in [%f, %f], len(t) == %d'
                 % (min(times), max(times), len(times)))
    log.critical('  Quantities [SI units]:')
    for name in ['stage', 'xmomentum', 'ymomentum', 'elevation']:
        q = outfile.variables[name][:]    # .flatten()
        log.critical('    %s in [%f, %f]' % \
                        (name, num.min(q), num.max(q)))        
                 

def _assert_lat_long(minlat, maxlat, minlon, maxlon):
    """Check latitudes and longitudes for validity."""
    
    msg = 'Must use latitudes and longitudes for minlat, maxlon etc'

    if minlat != None:
        assert -90 < minlat < 90 , msg
    if maxlat != None:
        assert -90 < maxlat < 90 , msg
        if minlat != None:
            assert maxlat > minlat
    if minlon != None:
        assert -180 < minlon < 180 , msg
    if maxlon != None:
        assert -180 < maxlon < 180 , msg
        if minlon != None:
            assert maxlon > minlon

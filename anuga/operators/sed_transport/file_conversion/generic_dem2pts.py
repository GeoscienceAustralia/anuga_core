# external modules
import numpy as num

# ANUGA modules
import anuga.utilities.log as log
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a, \
                            netcdf_float

from generic_asc2dem import generic_asc2dem
                            

def generic_dem2pts(name_in, name_out=None, quantity_name=None,
            easting_min=None, easting_max=None,
            northing_min=None, northing_max=None,
            use_cache=False, verbose=False,):
    """Read raster file from the following NetCDF format (.dem)
    Generic function, created from dem2pts

    Example:

    ncols         3121
    nrows         1800
    xllcorner     722000
    yllcorner     5893000
    cellsize      25
    NODATA_value  -9999
    138.3698 137.4194 136.5062 135.5558 ..........

    name_in may be a .asc or .dem file to be converted.

    Convert to NetCDF pts format which is

    points:  (Nx2) float array
    elevation: N float array
    """

    kwargs = {'name_out': name_out,
              'quantity_name': quantity_name,
              'easting_min': easting_min,
              'easting_max': easting_max,
              'northing_min': northing_min,
              'northing_max': northing_max,
              'verbose': verbose}

    if use_cache is True:
        from caching import cache
        result = cache(_generic_dem2pts, name_in, kwargs,
                       dependencies = [name_in],
                       verbose = verbose)

    else:
        result = apply(_generic_dem2pts, [name_in], kwargs)

    return result


def _generic_dem2pts(name_in, name_out=None, quantity_name=None, verbose=False,
            easting_min=None, easting_max=None,
            northing_min=None, northing_max=None):
    """Read raster from the following NetCDF format (.dem)

    Internal function. See public function generic_dem2pts for details.
    """

    # FIXME: Can this be written feasibly using write_pts?

    import os
    from anuga.file.netcdf import NetCDFFile

    root = name_in[:-4]

    if name_in[-4:] == '.asc':
        intermediate = root + '.dem'
        if verbose:
            log.critical('Preconvert %s from asc to %s' % \
                                    (name_in, intermediate))
        asc2dem(name_in)
        name_in = intermediate
    elif name_in[-4:] != '.dem':
        raise IOError('Input file %s should be of type .asc or .dem.' % name_in)

    if name_out != None and basename_out[-4:] != '.pts':
        raise IOError('Input file %s should be of type .pts.' % name_out)

    # Get NetCDF
    infile = NetCDFFile(name_in, netcdf_mode_r) 

    if verbose: log.critical('Reading raster from %s' % (name_in))

    ncols = int(infile.ncols)
    nrows = int(infile.nrows)
    xllcorner = float(infile.xllcorner)  # Easting of lower left corner
    yllcorner = float(infile.yllcorner)  # Northing of lower left corner
    cellsize = float(infile.cellsize)
    NODATA_value = float(infile.NODATA_value)

    dem_elevation = infile.variables[quantity_name]

    zone = int(infile.zone)
    false_easting = float(infile.false_easting)
    false_northing = float(infile.false_northing)

    #print ncols, nrows, xllcorner,yllcorner, cellsize, NODATA_value, zone


    # Text strings
    projection = infile.projection
    datum = infile.datum
    units = infile.units

    #print projection, datum, units

    # Get output file
    if name_out == None:
        ptsname = root + '.pts'
    else:
        ptsname = name_out

    if verbose: log.critical('Store to NetCDF file %s' % ptsname)

    # NetCDF file definition
    outfile = NetCDFFile(ptsname, netcdf_mode_w)

    # Create new file
    outfile.institution = 'Geoscience Australia'
    outfile.description = 'NetCDF pts format for compact and portable ' \
                          'storage of spatial point data'

    # Assign default values
    if easting_min is None: easting_min = xllcorner
    if easting_max is None: easting_max = xllcorner + ncols*cellsize
    if northing_min is None: northing_min = yllcorner
    if northing_max is None: northing_max = yllcorner + nrows*cellsize


    #print easting_min, easting_max, northing_min, northing_max

    # Compute offsets to update georeferencing
    easting_offset = xllcorner - easting_min
    northing_offset = yllcorner - northing_min

    # Georeferencing
    outfile.zone = zone
    outfile.xllcorner = easting_min # Easting of lower left corner
    outfile.yllcorner = northing_min # Northing of lower left corner
    outfile.false_easting = false_easting
    outfile.false_northing = false_northing

    outfile.projection = projection
    outfile.datum = datum
    outfile.units = units

    # Grid info (FIXME: probably not going to be used, but heck)
    outfile.ncols = ncols
    outfile.nrows = nrows

    dem_elevation_r = num.reshape(dem_elevation, (nrows, ncols))
    totalnopoints = nrows*ncols



    #========================================
    # Do the preceeding with numpy
    #========================================
    y = num.arange(nrows,dtype=num.float)
    y = yllcorner + (nrows-1)*cellsize - y*cellsize

    x = num.arange(ncols,dtype=num.float)
    x = xllcorner + x*cellsize

    xx,yy = num.meshgrid(x,y)

    xx = xx.flatten()
    yy = yy.flatten()

    
    flag = num.logical_and(num.logical_and((xx <= easting_max),(xx >= easting_min)),
                           num.logical_and((yy <= northing_max),(yy >= northing_min)))

    
    dem = dem_elevation[:].flatten()


    id = num.where(flag)[0]

    xx = xx[id]
    yy = yy[id]
    dem = dem[id]


    clippednopoints = len(dem)
    #print clippedpoints
    
    #print xx
    #print yy
    #print dem

    data_flag = dem != NODATA_value

    data_id = num.where(data_flag)

    xx = xx[data_id]
    yy = yy[data_id]
    dem = dem[data_id]

    nn = clippednopoints - len(dem)

    nopoints = len(dem)


    if verbose:
        log.critical('There are %d values in the raster' % totalnopoints)
        log.critical('There are %d values in the clipped raster'
                     % clippednopoints)
        log.critical('There are %d NODATA_values in the clipped raster' % nn)

    outfile.createDimension('number_of_points', nopoints)
    outfile.createDimension('number_of_dimensions', 2) #This is 2d data

    # Variable definitions
    outfile.createVariable('points', netcdf_float, ('number_of_points',
                                                    'number_of_dimensions'))
    outfile.createVariable(quantity_name, netcdf_float, ('number_of_points',))

    # Get handles to the variables
    points = outfile.variables['points']
    elevation = outfile.variables[quantity_name]

    points[:,0] = xx - easting_min
    points[:,1] = yy - northing_min
    elevation[:] = dem


    infile.close()
    outfile.close()


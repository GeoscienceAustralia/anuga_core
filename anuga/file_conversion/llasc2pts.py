# external modules

import numpy as num
import anuga.utilities.log as log
from anuga.coordinate_transforms.redfearn import convert_from_latlon_to_utm

# ANUGA modules
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a, \
                            netcdf_float

def llasc2pts(name_in, name_out=None,
            use_cache=False, show_progress=False, verbose=False,):
    """Read Digital Elevation model from the following ASCII format (.asc)

    Example:
    ncols        968
    nrows        698
    xllcorner    144.966666666667
    yllcorner    -18.500000000000
    cellsize     0.004166666667
    NODATA_value -32767
    138.3698 137.4194 136.5062 135.5558 ..........

    Convert name_in (.asc) to NetCDF format (.pts)

    """

    kwargs = {'name_out': name_out, 'verbose': verbose, 'show_progress': show_progress}

    if use_cache is True:
        from anuga.caching import cache
        result = cache(_convert_dem_from_llasc2pts, name_in, kwargs,
                       dependencies=[name_in],
                       verbose=verbose)

    else:
        result = _convert_dem_from_llasc2pts(*[name_in], **kwargs)

    return result


def _convert_dem_from_llasc2pts(name_in, name_out = None,
                                show_progress=False,
                                verbose=False,):
    """Read Digital Elevation model from the following LL ASCII format (.asc)

    Internal function. See public function convert_dem_from_ascii2netcdf
    for details.
    """

    import os
    from anuga.file.netcdf import NetCDFFile

    #Read DEM data
    datafile = open(name_in)

    if verbose: log.critical('Reading DEM from %s' % (name_in))

    lines = datafile.readlines()
    datafile.close()

    if verbose: log.critical('Got %d lines' % len(lines))

    ncols = int(lines[0].split()[1].strip())
    nrows = int(lines[1].split()[1].strip())

    # Do cellsize (line 4) before line 2 and 3
    cellsize = float(lines[4].split()[1].strip())

    xref = lines[2].split()
    if xref[0].strip() == 'xllcorner':
        xllcorner = float(xref[1].strip()) 
    elif xref[0].strip() == 'xllcenter':
        xllcorner = float(xref[1].strip()) # - 0.5*cellsize # Correct offset
    else:
        msg = 'Unknown keyword: %s' % xref[0].strip()
        raise Exception(msg)

    yref = lines[3].split()
    if yref[0].strip() == 'yllcorner':
        yllcorner = float(yref[1].strip()) 
    elif yref[0].strip() == 'yllcenter':
        yllcorner = float(yref[1].strip()) # - 0.5*cellsize # Correct offset
    else:
        msg = 'Unknown keyword: %s' % yref[0].strip()
        raise Exception(msg)

    NODATA_value = float(lines[5].split()[1].strip())

    assert len(lines) == nrows + 6

    dem_elevation = num.loadtxt(lines, skiprows=6, dtype=float)


    totalnopoints = nrows*ncols

    y = num.arange(nrows,dtype=float)
    y = yllcorner + (nrows-1)*cellsize - y*cellsize

    x = num.arange(ncols,dtype=float)
    x = xllcorner + x*cellsize

    #print(xllcorner)
    #print(x)

    #print(yllcorner)
    #print(y)

    xx,yy = num.meshgrid(x,y)

    xx = xx.flatten()
    yy = yy.flatten()
    dem = dem_elevation[:].flatten()
    
    # ====================
    # remove NODATA points
    # ====================
    data_flag = dem != NODATA_value

    data_id = num.where(data_flag)

    xx = xx[data_id]
    yy = yy[data_id]
    dem = dem[data_id]

    nn = totalnopoints - len(dem)

    nopoints = len(dem)

    # =====================================
    # Convert xx and yy to UTM
    # =====================================
    points_UTM, zone = convert_from_latlon_to_utm(latitudes=yy, longitudes=xx, show_progress=show_progress)

    points_UTM = num.asarray(points_UTM, dtype=float)

    corners, zone_c = convert_from_latlon_to_utm(latitudes=yllcorner, longitudes=xllcorner)

    xllcorner = corners[0][0]
    yllcorner = corners[0][1]

    assert zone == zone_c

    points_UTM = points_UTM - corners

    # ===============================
    # Step up for writing to pts file
    # ===============================

    if name_out is None:
        netcdfname = name_in[:-4]+'.pts'
    else:
        netcdfname = name_out + '.pts'

    if verbose: log.critical('Store to NetCDF file %s' % netcdfname)

    # NetCDF file definition
    outfile = NetCDFFile(netcdfname, netcdf_mode_w)

    # Create new file
    outfile.institution = 'Geoscience Australia'
    outfile.description = 'NetCDF pts format for compact and portable ' \
                          'storage of spatial point data'


    # Georeferencing
    outfile.zone = zone
    outfile.xllcorner = xllcorner # Easting of lower left corner
    outfile.yllcorner = yllcorner # Northing of lower left corner
    
    # Default settings
    outfile.false_easting = 500000.0
    outfile.false_northing = 10000000.0
    outfile.projection = 'UTM'
    outfile.datum = 'WGS84'
    outfile.units = 'METERS'

    # Grid info (FIXME: probably not going to be used, but heck)
    outfile.ncols = ncols
    outfile.nrows = nrows




    if verbose:
        log.critical('There are %d values in the elevation' % totalnopoints)
        log.critical('There are %d NODATA_values in the clipped elevation' % nn)

    outfile.createDimension('number_of_points', nopoints)
    outfile.createDimension('number_of_dimensions', 2) #This is 2d data

    # Variable definitions
    outfile.createVariable('points', netcdf_float, ('number_of_points',
                                                    'number_of_dimensions'))
    outfile.createVariable('elevation', netcdf_float, ('number_of_points',))

    # Get handles to the variables
    points = outfile.variables['points']
    elevation = outfile.variables['elevation']

    points[:,:]= points_UTM
    elevation[:] = dem

    outfile.close()

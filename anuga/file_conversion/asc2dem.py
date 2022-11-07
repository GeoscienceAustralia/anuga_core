# external modules

import numpy as num
import anuga.utilities.log as log

# ANUGA modules
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a, \
                            netcdf_float

def asc2dem(name_in, name_out=None,
                                  use_cache=False,
                                  verbose=False):
    """Read Digital Elevation model from the following ASCII format (.asc)

    Example:
    ncols         3121
    nrows         1800
    xllcorner     722000
    yllcorner     5893000
    cellsize      25
    NODATA_value  -9999
    138.3698 137.4194 136.5062 135.5558 ..........

    Convert name_in (.asc) to NetCDF format (.dem)
    mimicking the ASCII format closely.

    An accompanying file with same basename but extension .prj must exist
    and is used to fix the UTM zone, datum, false northings and eastings.

    The prj format is assumed to be as

    Projection    UTM
    Zone          56
    Datum         WGS84
    Zunits        NO
    Units         METERS
    Spheroid      WGS84
    Xshift        0.0000000000
    Yshift        10000000.0000000000
    Parameters
    """

    kwargs = {'name_out': name_out, 'verbose': verbose}

    if use_cache is True:
        from caching import cache
        result = cache(_convert_dem_from_ascii2netcdf, name_in, kwargs,
                       dependencies=[name_in,
                                     name_in[:-4] + '.prj'],
                       verbose=verbose)

    else:
        result = _convert_dem_from_ascii2netcdf(*[name_in], **kwargs)

    return result


def _convert_dem_from_ascii2netcdf(name_in, name_out = None,
                                   verbose = False):
    """Read Digital Elevation model from the following ASCII format (.asc)

    Internal function. See public function convert_dem_from_ascii2netcdf
    for details.
    """

    import os
    from anuga.file.netcdf import NetCDFFile

    root = name_in[:-4]

    # Read Meta data
    if verbose: log.critical('Reading METADATA from %s' % (root + '.prj'))

    metadatafile = open(root + '.prj')
    metalines = metadatafile.readlines()
    metadatafile.close()

    L = metalines[0].strip().split()
    assert L[0].strip().lower() == 'projection'
    projection = L[1].strip()                   #TEXT

    L = metalines[1].strip().split()
    assert L[0].strip().lower() == 'zone'
    zone = int(L[1].strip())

    L = metalines[2].strip().split()
    assert L[0].strip().lower() == 'datum'
    datum = L[1].strip()                        #TEXT

    L = metalines[3].strip().split()
    assert L[0].strip().lower() == 'zunits'     #IGNORE
    zunits = L[1].strip()                       #TEXT

    L = metalines[4].strip().split()
    assert L[0].strip().lower() == 'units'
    units = L[1].strip()                        #TEXT

    L = metalines[5].strip().split()
    assert L[0].strip().lower() == 'spheroid'   #IGNORE
    spheroid = L[1].strip()                     #TEXT

    L = metalines[6].strip().split()
    assert L[0].strip().lower() == 'xshift'
    false_easting = float(L[1].strip())

    L = metalines[7].strip().split()
    assert L[0].strip().lower() == 'yshift'
    false_northing = float(L[1].strip())

    if name_in[-4:] != '.asc':
        raise IOError('Input file %s should be of type .asc.' % name_in)

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

    # Checks suggested by Joaquim Luis
    # Our internal representation of xllcorner
    # and yllcorner is non-standard.
    xref = lines[2].split()
    if xref[0].strip() == 'xllcorner':
        xllcorner = float(xref[1].strip()) # + 0.5*cellsize # Correct offset
    elif xref[0].strip() == 'xllcenter':
        xllcorner = float(xref[1].strip())
    else:
        msg = 'Unknown keyword: %s' % xref[0].strip()
        raise Exception(msg)

    yref = lines[3].split()
    if yref[0].strip() == 'yllcorner':
        yllcorner = float(yref[1].strip()) # + 0.5*cellsize # Correct offset
    elif yref[0].strip() == 'yllcenter':
        yllcorner = float(yref[1].strip())
    else:
        msg = 'Unknown keyword: %s' % yref[0].strip()
        raise Exception(msg)

    NODATA_value = float(lines[5].split()[1].strip())

    assert len(lines) == nrows + 6

    if name_out is None:
        netcdfname = name_in[:-4]+'.dem'
    else:
        netcdfname = name_out + '.dem'

    if verbose: log.critical('Store to NetCDF file %s' % netcdfname)

    # NetCDF file definition
    fid = NetCDFFile(netcdfname, netcdf_mode_w)

    #Create new file
    fid.institution = 'Geoscience Australia'
    fid.description = 'NetCDF DEM format for compact and portable storage ' \
                      'of spatial point data'

    fid.ncols = ncols
    fid.nrows = nrows
    fid.xllcorner = xllcorner
    fid.yllcorner = yllcorner
    fid.cellsize = cellsize
    fid.NODATA_value = NODATA_value

    fid.zone = zone
    fid.false_easting = false_easting
    fid.false_northing = false_northing
    fid.projection = projection
    fid.datum = datum
    fid.units = units

    # dimension definitions
    fid.createDimension('number_of_rows', nrows)
    fid.createDimension('number_of_columns', ncols)

    # variable definitions
    fid.createVariable('elevation', netcdf_float, ('number_of_rows',
                                                   'number_of_columns'))

    # Get handles to the variables
    elevation = fid.variables['elevation']

    #Store data
    import numpy

    datafile = open(name_in)
    elevation[:,:] = numpy.loadtxt(datafile, skiprows=6)
    datafile.close()

#    n = len(lines[6:])
#    for i, line in enumerate(lines[6:]):
#        fields = line.split()
#        if verbose and i % ((n+10)/10) == 0:
#            log.critical('Processing row %d of %d' % (i, nrows))
#
#        if len(fields) != ncols:
#            msg = 'Wrong number of columns in file "%s" line %d\n' % (name_in, i)
#            msg += 'I got %d elements, but there should have been %d\n' % (len(fields), ncols)
#            raise Exception(msg)
#
#        elevation[i, :] = num.array([float(x) for x in fields])

    fid.close()

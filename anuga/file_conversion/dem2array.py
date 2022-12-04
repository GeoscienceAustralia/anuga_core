
# external modules

import numpy as num

# ANUGA modules
import anuga.utilities.log as log
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a, \
                            netcdf_float

from .asc2dem import asc2dem
                            

def dem2array(filename, variable_name='elevation',
            easting_min=None, easting_max=None,
            northing_min=None, northing_max=None,
            use_cache=False, verbose=False,):
    """Read Digitial Elevation model from the following NetCDF format (.dem)

    Example:

    ncols         3121
    nrows         1800
    xllcorner     722000
    yllcorner     5893000
    cellsize      25
    NODATA_value  -9999
    138.3698 137.4194 136.5062 135.5558 ..........

    name_in should be .dem file to be read.

    """




    import os
    from anuga.file.netcdf import NetCDFFile




    msg = 'Filename must be a text string'
    assert isinstance(filename, str), msg
    

        
    msg = 'Extension should be .dem'
    assert os.path.splitext(filename)[1] in ['.dem'], msg
    
    msg = 'Variable name must be a text string'
    assert isinstance(variable_name, str), msg
     


    # Get NetCDF
    infile = NetCDFFile(filename, netcdf_mode_r) 

    if verbose: log.critical('Reading DEM from %s' % (filename))

    ncols = int(infile.ncols)
    nrows = int(infile.nrows)
    xllcorner = float(infile.xllcorner)  # Easting of lower left corner
    yllcorner = float(infile.yllcorner)  # Northing of lower left corner
    cellsize = float(infile.cellsize)
    NODATA_value = float(infile.NODATA_value)


    zone = int(infile.zone)
    false_easting = float(infile.false_easting)
    false_northing = float(infile.false_northing)
    
    # Text strings
    projection = infile.projection
    datum = infile.datum
    units = infile.units
    
    Z = infile.variables[variable_name][:]
    Z = Z.reshape(nrows,ncols)
    Z = num.where(Z == NODATA_value , num.nan, Z)
    #changed the orientation of Z array to make it consistent with grd2array result
    Z = num.fliplr(Z.T)

    #print ncols, nrows, xllcorner,yllcorner, cellsize, NODATA_value, zone

    x = num.linspace(xllcorner, xllcorner+(ncols-1)*cellsize, ncols)
    y = num.linspace(yllcorner, yllcorner+(nrows-1)*cellsize, nrows)

    return x,y, Z





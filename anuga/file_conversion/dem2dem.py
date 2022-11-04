""" Load a DEM file, decimate it, and resave it.
"""

import numpy as num

from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_float

def dem2dem(name_in, stencil, cellsize_new, name_out=None,
                 verbose=False):
    """Read Digitial Elevation model from the following NetCDF format (.dem)

    Example:

    ncols         3121
    nrows         1800
    xllcorner     722000
    yllcorner     5893000
    cellsize      25
    NODATA_value  -9999
    138.3698 137.4194 136.5062 135.5558 ..........

    Decimate data to cellsize_new using stencil and write to NetCDF dem format.
    """

    import os
    from anuga.file.netcdf import NetCDFFile

    if name_in[-4:] != '.dem':
        raise IOError('Input file %s should be of type .dem.' % name_in)

    if name_out != None and basename_out[-4:] != '.dem':
        raise IOError('Input file %s should be of type .dem.' % name_out)

    #Open existing netcdf file to read
    infile = NetCDFFile(name_in, netcdf_mode_r)

    if verbose: log.critical('Reading DEM from %s' % inname)

    # Read metadata (convert from numpy.int32 to int where appropriate)
    ncols = int(infile.ncols)
    nrows = int(infile.nrows)
    xllcorner = infile.xllcorner
    yllcorner = infile.yllcorner
    cellsize = int(infile.cellsize)
    NODATA_value = int(infile.NODATA_value)
    zone = int(infile.zone)
    false_easting = infile.false_easting
    false_northing = infile.false_northing
    projection = infile.projection
    datum = infile.datum
    units = infile.units

    dem_elevation = infile.variables['elevation']

    #Get output file name
    if name_out is None:
        outname = name_in[:-4] + '_' + repr(cellsize_new) + '.dem'
    else:
        outname = name_out

    if verbose: log.critical('Write decimated NetCDF file to %s' % outname)

    #Determine some dimensions for decimated grid
    (nrows_stencil, ncols_stencil) = stencil.shape
    x_offset = ncols_stencil/ 2
    y_offset = nrows_stencil/ 2
    cellsize_ratio = int(cellsize_new/ cellsize)
    ncols_new = 1 + (ncols - ncols_stencil)// cellsize_ratio
    nrows_new = 1 + (nrows - nrows_stencil)// cellsize_ratio

    #print type(ncols_new), ncols_new
    
    #Open netcdf file for output
    outfile = NetCDFFile(outname, netcdf_mode_w)

    #Create new file
    outfile.institution = 'Geoscience Australia'
    outfile.description = 'NetCDF DEM format for compact and portable ' \
                          'storage of spatial point data'

    #Georeferencing
    outfile.zone = zone
    outfile.projection = projection
    outfile.datum = datum
    outfile.units = units

    outfile.cellsize = cellsize_new
    outfile.NODATA_value = NODATA_value
    outfile.false_easting = false_easting
    outfile.false_northing = false_northing

    outfile.xllcorner = xllcorner + (x_offset * cellsize)
    outfile.yllcorner = yllcorner + (y_offset * cellsize)
    outfile.ncols = ncols_new
    outfile.nrows = nrows_new

    # dimension definition
    #print nrows_new, ncols_new, nrows_new*ncols_new
    #print type(nrows_new), type(ncols_new), type(nrows_new*ncols_new)
    outfile.createDimension('number_of_points', nrows_new*ncols_new)

    # variable definition
    outfile.createVariable('elevation', netcdf_float, ('number_of_points',))

    # Get handle to the variable
    elevation = outfile.variables['elevation']

    dem_elevation_r = num.reshape(dem_elevation, (nrows, ncols))

    #Store data
    global_index = 0
    for i in range(nrows_new):
        if verbose: log.critical('Processing row %d of %d' % (i, nrows_new))

        lower_index = global_index
        telev = num.zeros(ncols_new, float)
        local_index = 0
        trow = i * cellsize_ratio

        for j in range(ncols_new):
            tcol = j * cellsize_ratio
            tmp = dem_elevation_r[trow:trow+nrows_stencil,
                                  tcol:tcol+ncols_stencil]

            #if dem contains 1 or more NODATA_values set value in
            #decimated dem to NODATA_value, else compute decimated
            #value using stencil
            if num.sum(num.sum(num.equal(tmp, NODATA_value))) > 0:
                telev[local_index] = NODATA_value
            else:
                telev[local_index] = num.sum(num.sum(tmp * stencil))

            global_index += 1
            local_index += 1

        upper_index = global_index

        elevation[lower_index:upper_index] = telev

    assert global_index == nrows_new*ncols_new, \
           'index not equal to number of points'

    infile.close()
    outfile.close()


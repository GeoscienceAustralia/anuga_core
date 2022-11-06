"""
    Module to convert SWW to DEM files.
"""


# external modules

import os
import numpy as num

# ANUGA modules
from anuga.abstract_2d_finite_volumes.util import remove_lone_verts     
from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.utilities.system_tools import get_vars_in_expression
import anuga.utilities.log as log
from anuga.utilities.file_utils import get_all_swwfiles


######
# formula mappings
######

quantity_formula = {'momentum':'(xmomentum**2 + ymomentum**2)**0.5',
                    'depth':'stage-elevation',
                    'speed': \
 '(xmomentum**2 + ymomentum**2)**0.5/(stage-elevation+1.e-6/(stage-elevation))'}



# Default block size for sww2dem()
DEFAULT_BLOCK_SIZE = 10000

def sww2dem(name_in, name_out,
            quantity=None, # defaults to elevation
            reduction=None,
            cellsize=10,
            number_of_decimal_places=None,
            NODATA_value=-9999.0,
            easting_min=None,
            easting_max=None,
            northing_min=None,
            northing_max=None,
            verbose=False,
            origin=None,
            datum='WGS84',
            block_size=None):
    """Read SWW file and convert to Digitial Elevation model format
    (.asc or .ers)

    Example (ASC):
    ncols         3121
    nrows         1800
    xllcorner     722000
    yllcorner     5893000
    cellsize      25
    NODATA_value  -9999
    138.3698 137.4194 136.5062 135.5558 ..........

    The number of decimal places can be specified by the user to save
    on disk space requirements by specifying in the call to sww2dem.

    Also write accompanying file with same basename_in but extension .prj
    used to fix the UTM zone, datum, false northings and eastings.

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

    The parameter quantity must be the name of an existing quantity or
    an expression involving existing quantities. The default is
    'elevation'. Quantity is not a list of quantities.

    If reduction is given and it's an index, sww2dem will output the quantity at that time-step. 
    If reduction is given and it's a built in function (eg max, min, mean), then that 
    function is used to reduce the quantity over all time-steps. If reduction is not given, 
    reduction is set to "max" by default.

    datum

    format can be either 'asc' or 'ers'
    block_size - sets the number of slices along the non-time axis to
                 process in one block.
    """

    import sys
    import types

    from anuga.geometry.polygon import inside_polygon, outside_polygon
    from anuga.abstract_2d_finite_volumes.util import \
         apply_expression_to_dictionary

    basename_in, in_ext = os.path.splitext(name_in)
    basename_out, out_ext = os.path.splitext(name_out)
    out_ext = out_ext.lower()

    if in_ext != '.sww':
        raise IOError('Input format for %s must be .sww' % name_in)

    if out_ext not in ['.asc', '.ers']:
        raise IOError('Format for %s must be either asc or ers.' % name_out)

    false_easting = 500000
    false_northing = 10000000

    if quantity is None:
        quantity = 'elevation'
    
    if reduction is None:
        reduction = max

    if quantity in quantity_formula:
        quantity = quantity_formula[quantity]

    if number_of_decimal_places is None:
        number_of_decimal_places = 3

    if block_size is None:
        block_size = DEFAULT_BLOCK_SIZE

    assert(isinstance(block_size, (int, int, float)))

    # Read sww file
    if verbose:
        log.critical('Reading from %s' % name_in)
        log.critical('Output directory is %s' % name_out)

    from anuga.file.netcdf import NetCDFFile
    fid = NetCDFFile(name_in)

    #Get extent and reference
    x = num.array(fid.variables['x'], float)
    y = num.array(fid.variables['y'], float)
    volumes = num.array(fid.variables['volumes'], int)
    if type(reduction) is not types.BuiltinFunctionType:
        times = fid.variables['time'][reduction]
    else:
        times = fid.variables['time'][:]

    number_of_timesteps = fid.dimensions['number_of_timesteps']
    number_of_points = fid.dimensions['number_of_points']

    if origin is None:
        # Get geo_reference
        # sww files don't have to have a geo_ref
        try:
            geo_reference = Geo_reference(NetCDFObject=fid)
        except AttributeError as e:
            geo_reference = Geo_reference() # Default georef object

        xllcorner = geo_reference.get_xllcorner()
        yllcorner = geo_reference.get_yllcorner()
        zone = geo_reference.get_zone()
    else:
        zone = origin[0]
        xllcorner = origin[1]
        yllcorner = origin[2]

    # FIXME: Refactor using code from Interpolation_function.statistics
    # (in interpolate.py)
    # Something like print swwstats(swwname)
    if verbose:
        log.critical('------------------------------------------------')
        log.critical('Statistics of SWW file:')
        log.critical('  Name: %s' % name_in)
        log.critical('  Reference:')
        log.critical('    Lower left corner: [%f, %f]' % (xllcorner, yllcorner))
        if type(reduction) is not types.BuiltinFunctionType:
            log.critical('    Time: %f' % times)
        else:
            log.critical('    Start time: %f' % fid.starttime[0])
        log.critical('  Extent:')
        log.critical('    x [m] in [%f, %f], len(x) == %d'
                     %(num.min(x), num.max(x), len(x.flat)))
        log.critical('    y [m] in [%f, %f], len(y) == %d'
                     % (num.min(y), num.max(y), len(y.flat)))
        if type(reduction) is not types.BuiltinFunctionType:
            log.critical('    t [s] = %f, len(t) == %d' % (times, 1))
        else:
            log.critical('    t [s] in [%f, %f], len(t) == %d'
                         % (min(times), max(times), len(times)))
        log.critical('  Quantities [SI units]:')
        
        # Comment out for reduced memory consumption
        for name in ['stage', 'xmomentum', 'ymomentum']:
            q = fid.variables[name][:].flatten()
            if type(reduction) is not types.BuiltinFunctionType:
                q = q[reduction*len(x):(reduction+1)*len(x)]
            if verbose: log.critical('    %s in [%f, %f]'
                                     % (name, min(q), max(q)))
        for name in ['elevation']:
            q = fid.variables[name][:].flatten()
            if verbose: log.critical('    %s in [%f, %f]'
                                     % (name, min(q), max(q)))

    # Get the variables in the supplied expression.
    # This may throw a SyntaxError exception.
    var_list = get_vars_in_expression(quantity)

    # Check that we have the required variables in the SWW file.
    missing_vars = []
    for name in var_list:
        try:
            _ = fid.variables[name]
        except KeyError:
            missing_vars.append(name)
    if missing_vars:
        msg = ("In expression '%s', variables %s are not in the SWW file '%s'"
               % (quantity, str(missing_vars), name_in))
        raise(Exception, msg)

    # Create result array and start filling, block by block.
    result = num.zeros(number_of_points, float)

    if verbose:
        msg = 'Slicing sww file, num points: ' + str(number_of_points)
        msg += ', block size: ' + str(block_size)
        log.critical(msg)

    for start_slice in range(0, number_of_points, block_size):
        # Limit slice size to array end if at last block
        end_slice = min(start_slice + block_size, number_of_points)
        
        # Get slices of all required variables
        q_dict = {}
        for name in var_list:
            # check if variable has time axis
            if len(fid.variables[name].shape) == 2:
                q_dict[name] = fid.variables[name][:,start_slice:end_slice]
            else:       # no time axis
                q_dict[name] = fid.variables[name][start_slice:end_slice]

        # Evaluate expression with quantities found in SWW file
        res = apply_expression_to_dictionary(quantity, q_dict)

        if len(res.shape) == 2:
            new_res = num.zeros(res.shape[1], float)
            for k in range(res.shape[1]):
                if type(reduction) is not types.BuiltinFunctionType:
                    new_res[k] = res[reduction,k]
                else:
                    new_res[k] = reduction(res[:,k])
            res = new_res

        result[start_slice:end_slice] = res
                                    
    # Post condition: Now q has dimension: number_of_points
    assert len(result.shape) == 1
    assert result.shape[0] == number_of_points

    if verbose:
        log.critical('Processed values for %s are in [%f, %f]'
                     % (quantity, min(result), max(result)))

    # Create grid and update xll/yll corner and x,y
    # Relative extent
    if easting_min is None:
        xmin = min(x)
    else:
        xmin = easting_min - xllcorner

    if easting_max is None:
        xmax = max(x)
    else:
        xmax = easting_max - xllcorner

    if northing_min is None:
        ymin = min(y)
    else:
        ymin = northing_min - yllcorner

    if northing_max is None:
        ymax = max(y)
    else:
        ymax = northing_max - yllcorner

    msg = 'xmax must be greater than or equal to xmin.\n'
    msg += 'I got xmin = %f, xmax = %f' %(xmin, xmax)
    assert xmax >= xmin, msg

    msg = 'ymax must be greater than or equal to xmin.\n'
    msg += 'I got ymin = %f, ymax = %f' %(ymin, ymax)
    assert ymax >= ymin, msg

    if verbose: log.critical('Creating grid')
    ncols = int((xmax-xmin)/cellsize) + 1
    nrows = int((ymax-ymin)/cellsize) + 1

    # New absolute reference and coordinates
    newxllcorner = xmin + xllcorner
    newyllcorner = ymin + yllcorner

    x = x + xllcorner - newxllcorner
    y = y + yllcorner - newyllcorner

    vertex_points = num.concatenate ((x[:,num.newaxis], y[:,num.newaxis]), axis=1)
    assert len(vertex_points.shape) == 2


    def calc_grid_values_old(vertex_points, volumes, result):

        grid_points = num.zeros ((ncols*nrows, 2), float)

        for i in range(nrows):
            if out_ext == '.asc':
                yg = i * cellsize
            else:
                # this will flip the order of the y values for ers
                yg = (nrows-i) * cellsize

            for j in range(ncols):
                xg = j * cellsize
                k = i*ncols + j

                grid_points[k, 0] = xg
                grid_points[k, 1] = yg

        # Interpolate
        from anuga.fit_interpolate.interpolate import Interpolate

        # Remove loners from vertex_points, volumes here
        vertex_points, volumes = remove_lone_verts(vertex_points, volumes)
        # export_mesh_file('monkey.tsh',{'vertices':vertex_points, 'triangles':volumes})


        interp = Interpolate(vertex_points, volumes, verbose = verbose)

        bprint = 0

        # Interpolate using quantity values
        if verbose: log.critical('Interpolating')
        grid_values = interp.interpolate(bprint, result, grid_points).flatten()
        outside_indices = interp.get_outside_poly_indices()

        for i in outside_indices:
            #print 'change grid_value',NODATA_value
            grid_values[i] = NODATA_value
	
        return grid_values

    def calc_grid_values(vertex_points, volumes, result):

        grid_points = num.zeros ((ncols*nrows, 2), float)

        for i in range(nrows):
            if out_ext == '.asc':
                yg = i * cellsize
            else:
                #this will flip the order of the y values for ers
                yg = (nrows-i) * cellsize
   
            for j in range(ncols):
                xg = j * cellsize
                k = i*ncols + j

                grid_points[k, 0] = xg
                grid_points[k, 1] = yg
    
        grid_values = num.zeros(ncols*nrows, float)

        eval_grid(nrows, ncols, NODATA_value, grid_points, vertex_points.flatten(), volumes, result, grid_values);
        return grid_values.flatten()

    grid_values = calc_grid_values(vertex_points, volumes, result)

    if verbose:
        log.critical('Interpolated values are in [%f, %f]'
                     % (num.min(grid_values), num.max(grid_values)))

    # Assign NODATA_value to all points outside bounding polygon (from interpolation mesh)
    
#    P = interp.mesh.get_boundary_polygon()
#    outside_indices = outside_polygon(grid_points, P, closed=True)

    if out_ext == '.ers':
        # setup ERS header information
        grid_values = num.reshape(grid_values, (nrows, ncols))
        header = {}
        header['datum'] = '"' + datum + '"'
        # FIXME The use of hardwired UTM and zone number needs to be made optional
        # FIXME Also need an automatic test for coordinate type (i.e. EN or LL)
        header['projection'] = '"UTM-' + str(zone) + '"'
        header['coordinatetype'] = 'EN'
        if header['coordinatetype'] == 'LL':
            header['longitude'] = str(newxllcorner)
            header['latitude'] = str(newyllcorner)
        elif header['coordinatetype'] == 'EN':
            header['eastings'] = str(newxllcorner)
            header['northings'] = str(newyllcorner)
        header['nullcellvalue'] = str(NODATA_value)
        header['xdimension'] = str(cellsize)
        header['ydimension'] = str(cellsize)
        header['value'] = '"' + quantity + '"'
        #header['celltype'] = 'IEEE8ByteReal'  #FIXME: Breaks unit test

        #Write
        if verbose:
            log.critical('Writing %s' % name_out)

        import ermapper_grids

        ermapper_grids.write_ermapper_grid(name_out, grid_values, header)

        fid.close()
   
    else:
        #Write to Ascii format
        #Write prj file
        prjfile = basename_out + '.prj'

        if verbose: log.critical('Writing %s' % prjfile)
        prjid = open(prjfile, 'w')
        prjid.write('Projection    %s\n' %'UTM')
        prjid.write('Zone          %d\n' %zone)
        prjid.write('Datum         %s\n' %datum)
        prjid.write('Zunits        NO\n')
        prjid.write('Units         METERS\n')
        prjid.write('Spheroid      %s\n' %datum)
        prjid.write('Xshift        %d\n' %false_easting)
        prjid.write('Yshift        %d\n' %false_northing)
        prjid.write('Parameters\n')
        prjid.close()

        if verbose: log.critical('Writing %s' % name_out)

        ascid = open(name_out, 'w')

        ascid.write('ncols         %d\n' %ncols)
        ascid.write('nrows         %d\n' %nrows)
        ascid.write('xllcorner     %d\n' %newxllcorner)
        ascid.write('yllcorner     %d\n' %newyllcorner)
        ascid.write('cellsize      %f\n' %cellsize)
        ascid.write('NODATA_value  %d\n' %NODATA_value)

        #Get bounding polygon from mesh
        #P = interp.mesh.get_boundary_polygon()
        #inside_indices = inside_polygon(grid_points, P)

        # change printoptions so that a long string of zeros in not
        # summarized as [0.0, 0.0, 0.0, ... 0.0, 0.0, 0.0]
        #printoptions = num.get_printoptions()
        #num.set_printoptions(threshold=sys.maxint)

        format = '%.'+'%g' % number_of_decimal_places +'e'
        for i in range(nrows):
            if verbose and i % ((nrows+10)//10) == 0:
                log.critical('Doing row %d of %d' % (i, nrows))

            base_index = (nrows-i-1)*ncols

            slice = grid_values[base_index:base_index+ncols]

            num.savetxt(ascid, slice.reshape(1,ncols), format, ' ' )
            
        
        #Close
        ascid.close()
        fid.close()
     
        return basename_out



def sww2dem_batch(basename_in, extra_name_out=None,
                quantities=None, # defaults to elevation
                reduction=None,
                cellsize=10,
                number_of_decimal_places=None,
                NODATA_value=-9999,
                easting_min=None,
                easting_max=None,
                northing_min=None,
                northing_max=None,
                verbose=False,
                origin=None,
                datum='WGS84',
                format='ers'):
    """Wrapper for sww2dem.
    See sww2dem to find out what most of the parameters do. Note that since this
    is a batch command, the normal filename naming conventions do not apply.

    basename_in is a path to sww file/s, without the .sww extension.
    extra_name_out is a postfix to add to the output filename.

    Quantities is a list of quantities.  Each quantity will be
    calculated for each sww file.

    This returns the basenames of the files returned, which is made up
    of the dir and all of the file name, except the extension.

    This function returns the names of the files produced.

    It will also produce as many output files as there are input sww files.
    """

    if quantities is None:
        quantities = ['elevation']

    if type(quantities) is str:
            quantities = [quantities]

    # How many sww files are there?
    dir, base = os.path.split(basename_in)

    iterate_over = get_all_swwfiles(dir, base, verbose)

    if dir == "":
        dir = "." # Unix compatibility

    files_out = []
    for sww_file in iterate_over:
        for quantity in quantities:
            if extra_name_out is None:
                basename_out = sww_file + '_' + quantity
            else:
                basename_out = sww_file + '_' + quantity + '_' + extra_name_out

            swwin = dir+os.sep+sww_file+'.sww'
            demout = dir+os.sep+basename_out+'.'+format

            if verbose:
                log.critical('sww2dem: %s => %s' % (swwin, demout))

            file_out = sww2dem(swwin,
                               demout,
                               quantity,
                               reduction,
                               cellsize,
                               number_of_decimal_places,
                               NODATA_value,
                               easting_min,
                               easting_max,
                               northing_min,
                               northing_max,
                               verbose,
                               origin,
                               datum)
                               
            files_out.append(file_out)
    return files_out


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
DEFAULT_BLOCK_SIZE = 100000

def sww2array(name_in,
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
    """Read SWW file and convert to a numpy array (can be stored to a png file later)


    The parameter quantity must be the name of an existing quantity or
    an expression involving existing quantities. The default is
    'elevation'. Quantity is not a list of quantities.

    If reduction is given and it's an index, sww2array will output the quantity at that time-step.
    If reduction is given and it's a built in function (eg max, min, mean), then that 
    function is used to reduce the quantity over all time-steps. If reduction is not given, 
    reduction is set to "max" by default.

    datum


    block_size - sets the number of slices along the non-time axis to
                 process in one block.
    """

    import sys
    import types

    from anuga.geometry.polygon import inside_polygon, outside_polygon
    from anuga.abstract_2d_finite_volumes.util import \
         apply_expression_to_dictionary

    basename_in, in_ext = os.path.splitext(name_in)

    if in_ext != '.sww':
        raise IOError('Input format for %s must be .sww' % name_in)


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
        if type(reduction) is not types.BuiltinFunctionType:
            q_dict = {}
            for name in var_list:
                # check if variable has time axis
                if len(fid.variables[name].shape) == 2:
                    print('avoiding large array')
                    q_dict[name] = fid.variables[name][reduction,start_slice:end_slice]
                else:       # no time axis
                    q_dict[name] = fid.variables[name][start_slice:end_slice]

            # Evaluate expression with quantities found in SWW file
            res = apply_expression_to_dictionary(quantity, q_dict)

#            if len(res.shape) == 2:
#                new_res = num.zeros(res.shape[1], float)
#                for k in xrange(res.shape[1]):
#                    if type(reduction) is not types.BuiltinFunctionType:
#                        new_res[k] = res[k]
#                    else:
#                        new_res[k] = reduction(res[:,k])
#                res = new_res
        else:
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
    ncols = int((xmax-xmin)//cellsize) + 1
    nrows = int((ymax-ymin)//cellsize) + 1

    # New absolute reference and coordinates
    newxllcorner = xmin + xllcorner
    newyllcorner = ymin + yllcorner

    x = x + xllcorner - newxllcorner
    y = y + yllcorner - newyllcorner


    grid_values = num.zeros( (nrows*ncols, ), float)
    #print '---',grid_values.shape

    num_tri =  len(volumes)
    norms = num.zeros(6*num_tri, float)

    #Use fasr method to calc grid values
    from .calc_grid_values_ext import calc_grid_values

    calc_grid_values(nrows, ncols, cellsize, NODATA_value,
                     x,y, norms, volumes, result, grid_values)


    fid.close()
    
    #print outside_indices

    if verbose:
        log.critical('Interpolated values are in [%f, %f]'
                     % (num.min(grid_values), num.max(grid_values)))



    return x,y, grid_values.reshape(nrows,ncols)[::-1,:]






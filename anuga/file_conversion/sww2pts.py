
from builtins import range
import numpy as num
import os

from anuga.coordinate_transforms.geo_reference import Geo_reference

def sww2pts(name_in, name_out=None,
            data_points=None,
            quantity=None,
            timestep=None,
            reduction=None,
            NODATA_value=-9999,
            verbose=False,
            origin=None):
    """Read SWW file and convert to interpolated values at selected points

    The parameter 'quantity' must be the name of an existing quantity or
    an expression involving existing quantities. The default is 'elevation'.

    if timestep (an index) is given, output quantity at that timestep.

    if reduction is given use that to reduce quantity over all timesteps.

    data_points (Nx2 array) give locations of points where quantity is to 
    be computed.
    """

    import sys
    from anuga.geometry.polygon import inside_polygon, outside_polygon
    from anuga.abstract_2d_finite_volumes.util import \
             apply_expression_to_dictionary
    from anuga.geospatial_data.geospatial_data import Geospatial_data

    if quantity is None:
        quantity = 'elevation'

    if reduction is None:
        reduction = max

    basename_in, in_ext = os.path.splitext(name_in)
    
    if name_out != None:
        basename_out, out_ext = os.path.splitext(name_out)
    else:
        basename_out = basename_in + '_%s' % quantity
        out_ext = '.pts'
        name_out = basename_out + out_ext

    if in_ext != '.sww':
        raise IOError('Input format for %s must be .sww' % name_in)

    if out_ext != '.pts':
        raise IOError('Output format for %s must be .pts' % name_out)


    # Read sww file
    if verbose: log.critical('Reading from %s' % name_in)
    from anuga.file.netcdf import NetCDFFile
    fid = NetCDFFile(name_in)

    # Get extent and reference
    x = fid.variables['x'][:]
    y = fid.variables['y'][:]
    volumes = fid.variables['volumes'][:]


    try: # works with netcdf4
        number_of_timesteps = len(fid.dimensions['number_of_timesteps'])
        number_of_points = len(fid.dimensions['number_of_points'])
    except: #works with scientific.io.netcdf
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

    # FIXME: Refactor using code from file_function.statistics
    # Something like print swwstats(swwname)
    if verbose:
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        times = fid.variables['time'][:]
        log.critical('------------------------------------------------')
        log.critical('Statistics of SWW file:')
        log.critical('  Name: %s' % swwfile)
        log.critical('  Reference:')
        log.critical('    Lower left corner: [%f, %f]' % (xllcorner, yllcorner))
        log.critical('    Start time: %f' % fid.starttime[0])
        log.critical('  Extent:')
        log.critical('    x [m] in [%f, %f], len(x) == %d'
                     % (num.min(x), num.max(x), len(x.flat)))
        log.critical('    y [m] in [%f, %f], len(y) == %d'
                     % (num.min(y), num.max(y), len(y.flat)))
        log.critical('    t [s] in [%f, %f], len(t) == %d'
                     % (min(times), max(times), len(times)))
        log.critical('  Quantities [SI units]:')
        for name in ['stage', 'xmomentum', 'ymomentum', 'elevation']:
            q = fid.variables[name][:].flat
            log.critical('    %s in [%f, %f]' % (name, min(q), max(q)))

    # Get quantity and reduce if applicable
    if verbose: log.critical('Processing quantity %s' % quantity)

    # Turn NetCDF objects into numeric arrays
    quantity_dict = {}
    for name in list(fid.variables.keys()):
        quantity_dict[name] = fid.variables[name][:]

    # Convert quantity expression to quantities found in sww file
    q = apply_expression_to_dictionary(quantity, quantity_dict)

    if len(q.shape) == 2:
        # q has a time component and needs to be reduced along
        # the temporal dimension
        if verbose: log.critical('Reducing quantity %s' % quantity)

        q_reduced = num.zeros(number_of_points, float)
        for k in range(number_of_points):
            q_reduced[k] = reduction(q[:,k])
        q = q_reduced

    # Post condition: Now q has dimension: number_of_points
    assert len(q.shape) == 1
    assert q.shape[0] == number_of_points

    if verbose:
        log.critical('Processed values for %s are in [%f, %f]'
                     % (quantity, min(q), max(q)))

    # Create grid and update xll/yll corner and x,y
    vertex_points = num.concatenate((x[:, num.newaxis], y[:, num.newaxis]), axis=1)
    assert len(vertex_points.shape) == 2

    # Interpolate
    from anuga.fit_interpolate.interpolate import Interpolate
    interp = Interpolate(vertex_points, volumes, verbose=verbose)

    # Interpolate using quantity values
    if verbose: log.critical('Interpolating')
    interpolated_values = interp.interpolate(q, data_points).flatten()

    if verbose:
        log.critical('Interpolated values are in [%f, %f]'
                     % (num.min(interpolated_values),
                        num.max(interpolated_values)))

    # Assign NODATA_value to all points outside bounding polygon
    # (from interpolation mesh)
    P = interp.mesh.get_boundary_polygon()
    outside_indices = outside_polygon(data_points, P, closed=True)

    for i in outside_indices:
        interpolated_values[i] = NODATA_value

    # Store results
    G = Geospatial_data(data_points=data_points, attributes=interpolated_values)

    G.export_points_file(name_out, absolute = True)

    fid.close()


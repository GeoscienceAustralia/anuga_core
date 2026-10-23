#!/usr/bin/env python
"""File function
Takes a file as input, and returns it as a mathematical function.
For example, you can load an arbitrary 2D heightfield mesh, and treat it as a
function like so:

F = file_function('my_mesh.sww', ...)
evaluated_point = F(x, y)

Values will be interpolated across the surface of the mesh. Holes in the mesh
have an undefined value.

"""
import numpy as num

from anuga.geospatial_data.geospatial_data import ensure_absolute
from anuga.file.netcdf import NetCDFFile
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a
from anuga.utilities.numerical_tools import ensure_numeric

import anuga.utilities.log as log


def file_function(filename,
                  domain=None,
                  quantities=None,
                  interpolation_points=None,
                  use_relative_time=True,
                  time_thinning=1,
                  time_limit=None,
                  verbose=False,
                  use_cache=False,
                  boundary_polygon=None,
                  output_centroids=False):
    """Read time history of spatial and/or temporal data from NetCDF file and return
    a callable object.

    Input variables:

    filename - Name of sww, tms or sts file

       If the file has extension 'sww' then it is assumed to be spatio-temporal
       or temporal and the callable object will have the form f(t,x,y) or f(t)
       depending on whether the file contains spatial data

       If the file has extension 'tms' then it is assumed to be temporal only
       and the callable object will have the form f(t)

       Either form will return interpolated values based on the input file
       using the underlying interpolation_function.

    domain - Associated domain object
       If domain is specified, model time (domain.starttime)
       will be checked and possibly modified.

       All times are assumed to be in UTC

       All spatial information is assumed to be in absolute UTM coordinates.

    quantities - the name of the quantity to be interpolated or a
                 list of quantity names. The resulting function will return
                 a tuple of values - one for each quantity
                 If quantities are None, the default quantities are
                 ['stage', 'xmomentum', 'ymomentum']


    interpolation_points - list of absolute UTM coordinates for points (N x 2)
    or geospatial object or points file name at which values are sought

    use_relative_time -

    time_thinning -

    verbose -

    use_cache: True means that caching of intermediate result of
               Interpolation_function is attempted

    boundary_polygon -


    See Interpolation function in anuga.fit_interpolate.interpolation for
    further documentation
    """

    if quantities is None:
        if verbose:
            msg = 'Quantities specified in file_function are None,'
            msg += ' so using stage, xmomentum, and ymomentum in that order'
            log.info(msg)
        quantities = ['stage', 'xmomentum', 'ymomentum']

    # Use domain's startime if available
    if domain is not None:
        domain_starttime = domain.get_starttime()
    else:
        domain_starttime = None

    # Build arguments and keyword arguments for use with caching or apply.
    args = (filename,)

    # domain is not passed to the cache key directly because instance hash
    # codes change; only the scalar attributes needed (domain_starttime) are
    # extracted above so the cache key stays stable.
    kwargs = {'quantities': quantities,
              'interpolation_points': interpolation_points,
              'domain_starttime': domain_starttime,
              'use_relative_time': use_relative_time,
              'time_thinning': time_thinning,
              'time_limit': time_limit,
              'verbose': verbose,
              'boundary_polygon': boundary_polygon,
              'output_centroids': output_centroids}

    # Call underlying engine with or without caching
    if use_cache is True:
        try:
            from anuga.caching import cache
        except ImportError:
            msg = 'Caching was requested, but caching module'+\
                  'could not be imported'
            raise Exception(msg)

        f, starttime = cache(_file_function,
                             args, kwargs,
                             dependencies=[filename],
                             compression=False,
                             verbose=verbose)
    else:
        f, starttime = _file_function(*args, **kwargs)

    f.starttime = starttime
    f.filename = filename

    if domain is not None:
        #Update domain.startime if it is *earlier* than starttime from file
        if starttime > domain.starttime:
            msg = 'WARNING: Start time as specified in domain (%f)' \
                  % domain.starttime
            msg += ' is earlier than the starttime of file %s (%f).' \
                     % (filename, starttime)
            msg += ' Modifying domain starttime accordingly.'

            if verbose: log.info(msg)

            domain.set_starttime(starttime) #Modifying model time

            if verbose: log.info('Domain starttime is now set to %f'
                                     % domain.starttime)
    return f


def _file_function(filename,
                   quantities=None,
                   interpolation_points=None,
                   domain_starttime=None,
                   use_relative_time=True,
                   time_thinning=1,
                   time_limit=None,
                   verbose=False,
                   boundary_polygon=None,
                   output_centroids=False):
    """Internal function

    See file_function for documentatiton
    """

    assert isinstance(filename, str),\
               'First argument to File_function must be a string'

    import os
    ext = os.path.splitext(filename)[1]
    msg = 'Extension should be csv  sww, tms or sts '
    assert ext in [".csv",  ".sww", ".tms", ".sts"], msg


    if ext in [".sww", ".tms", ".sts"]:
        return get_netcdf_file_function(filename,
                                        quantities,
                                        interpolation_points,
                                        domain_starttime,
                                        use_relative_time=use_relative_time,
                                        time_thinning=time_thinning,
                                        time_limit=time_limit,
                                        verbose=verbose,
                                        boundary_polygon=boundary_polygon,
                                        output_centroids=output_centroids)
    else:
        raise Exception('file_function requires a .sww, .tms or .sts file, got %s' % ext)


def get_netcdf_file_function(filename,
                             quantity_names=None,
                             interpolation_points=None,
                             domain_starttime=None,
                             use_relative_time=False,
                             time_thinning=1,
                             time_limit=None,
                             verbose=False,
                             boundary_polygon=None,
                             output_centroids=False):
    """Read time history of spatial data from NetCDF sww file and
    return a callable object f(t,x,y)
    which will return interpolated values based on the input file.

    Model time (domain_starttime)
    will be checked, possibly modified and returned

    All times are assumed to be in UTC

    See Interpolation function for further documentation
    """

    # Origin reconciliation: interpolation_points are supplied as absolute UTM
    # coordinates and are shifted by (xllcorner, yllcorner) below to match the
    # file's local coordinate frame.  A full zone-consistency check between
    # domain and file would require passing domain.geo_reference through the
    # call chain — left for a future enhancement.

    import time
    import calendar
    from anuga.config import time_format

    # Open NetCDF file
    if verbose: log.info('Reading %s' % filename)

    fid = NetCDFFile(filename, netcdf_mode_r)

    if isinstance(quantity_names, str):
        quantity_names = [quantity_names]

    if quantity_names is None or len(quantity_names) < 1:
        msg = 'No quantities are specified in file_function'
        raise Exception(msg)

    if interpolation_points is not None:

        #interpolation_points = num.array(interpolation_points, float)
        interpolation_points = ensure_absolute(interpolation_points)
        msg = 'Points must by N x 2. I got %d' % interpolation_points.shape[1]
        assert interpolation_points.shape[1] == 2, msg

    # Now assert that requested quantitites (and the independent ones)
    # are present in file
    missing = []
    for quantity in ['time'] + quantity_names:
        if quantity not in fid.variables:
            missing.append(quantity)

    if len(missing) > 0:
        msg = 'Quantities %s could not be found in file %s'\
              % (str(missing), filename)
        fid.close()
        raise Exception(msg)

    # Decide whether this data has a spatial dimension
    spatial = True
    for quantity in ['x', 'y']:
        if quantity not in fid.variables:
            spatial = False

    if filename[-3:] == 'tms' and spatial is True:
        msg = 'Files of type TMS must not contain spatial information'
        raise Exception(msg)

    if filename[-3:] == 'sww' and spatial is False:
        msg = 'Files of type SWW must contain spatial information'
        raise Exception(msg)

    if filename[-3:] == 'sts' and spatial is False:
        #What if mux file only contains one point
        msg = 'Files of type STS must contain spatial information'
        raise Exception(msg)

    # Get first timestep
    try:
        starttime = float(fid.starttime)
    except ValueError:
        msg = 'Could not read starttime from file %s' % filename
        raise Exception(msg)


    # Get variables
    # if verbose: log.info('Get variables'    )
    time = fid.variables['time'][:]

    if not use_relative_time:
        time = time + starttime

    msg = 'Time vector in %s is not monotonically increasing' % filename
    assert num.all(num.diff(time) > 0), msg

    # Apply time limit if requested
    upper_time_index = len(time)
    msg = 'Time vector obtained from file %s has length 0' % filename
    assert upper_time_index > 0, msg

    if time_limit is not None:
        # Adjust given time limit to given start time
        time_limit = time_limit - starttime


        # Find limit point
        for i, t in enumerate(time):
            if t > time_limit:
                upper_time_index = i
                break

        msg = 'Time vector is zero. Requested time limit is %f' % time_limit
        assert upper_time_index > 0, msg

        if time_limit < time[-1] and verbose is True:
            log.info('Limited time vector from %.2fs to %.2fs'
                         % (time[-1], time_limit))

    time = time[:upper_time_index]




    # Get time independent stuff
    if spatial:
        # Get origin — use Geo_reference so zone, hemisphere and EPSG are
        # read consistently (same as sww_merge, SWW_file, etc.)
        from anuga.coordinate_transforms.geo_reference import Geo_reference
        geo_ref = Geo_reference(NetCDFObject=fid)
        xllcorner = geo_ref.get_xllcorner()
        yllcorner = geo_ref.get_yllcorner()
        zone = geo_ref.zone

        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        if filename.endswith('sww'):
            triangles = fid.variables['volumes'][:]

        x = num.reshape(x, (len(x), 1))
        y = num.reshape(y, (len(y), 1))
        vertex_coordinates = num.concatenate((x, y), axis=1) #m x 2 array

        if boundary_polygon is not None:
            # Restrict to STS points that lie on the boundary polygon.
            # STS points are produced using an ordering file so all should
            # be present, but we filter defensively.
            boundary_polygon=ensure_numeric(boundary_polygon)
            boundary_polygon[:, 0] -= xllcorner
            boundary_polygon[:, 1] -= yllcorner
            temp=[]
            boundary_id=[]
            gauge_id=[]
            for i in range(len(boundary_polygon)):
                for j in range(len(x)):
                    if num.allclose(vertex_coordinates[j],
                                    boundary_polygon[i], rtol=1e-4, atol=1e-4):
                        # Loose tolerance because gauge lat/lon stored as
                        # float32 in STS files, causing slight coordinate drift
                        # when cast to float64.
                        temp.append(boundary_polygon[i])
                        gauge_id.append(j)
                        boundary_id.append(i)
                        break
            gauge_neighbour_id=[]
            for i in range(len(boundary_id)-1):
                if boundary_id[i]+1==boundary_id[i+1]:
                    gauge_neighbour_id.append(i+1)
                else:
                    gauge_neighbour_id.append(-1)
            if boundary_id[len(boundary_id)-1]==len(boundary_polygon)-1 \
               and boundary_id[0]==0:
                gauge_neighbour_id.append(0)
            else:
                gauge_neighbour_id.append(-1)
            gauge_neighbour_id=ensure_numeric(gauge_neighbour_id)


            if len(num.compress(gauge_neighbour_id>=0, gauge_neighbour_id)) \
               != len(temp)-1:
                msg='incorrect number of segments'
                raise Exception(msg)
            vertex_coordinates=ensure_numeric(temp)
            if len(vertex_coordinates)==0:
                msg = 'None of the sts gauges fall on the boundary'
                raise Exception(msg)
        else:
            gauge_neighbour_id=None

        if interpolation_points is not None:
            # Adjust for georef
            interpolation_points[:, 0] -= xllcorner
            interpolation_points[:, 1] -= yllcorner
    else:
        gauge_neighbour_id=None

    if domain_starttime is not None:
        # If domain_startime is *later* than starttime,
        # move time back - relative to domain's time
        if domain_starttime > starttime:
            time = time - domain_starttime + starttime

    if verbose:
        log.info('File_function data obtained from: %s' % filename)
        log.info('  References:')
        if spatial:
            log.info('    Lower left corner: [%f, %f]'
                         % (xllcorner, yllcorner))
        log.info('    Start time:   %f' % starttime)


    # Produce values for desired data points at
    # each timestep for each quantity
    quantities = {}
    for i, name in enumerate(quantity_names):
        quantities[name] = fid.variables[name][:]
        if boundary_polygon is not None:
            #removes sts points that do not lie on boundary
            quantities[name] = num.take(quantities[name], gauge_id, axis=1)

    # Close sww, tms or sts netcdf file
    fid.close()

    from anuga.fit_interpolate.interpolate import Interpolation_function

    if not spatial:
        vertex_coordinates = triangles = interpolation_points = None
    if filename[-3:] == 'sts':#added
        triangles = None
        #vertex coordinates is position of urs gauges

    if verbose:
        log.info('Calling interpolation function')

    # Return Interpolation_function instance as well as
    # starttime for use to possible modify that of domain
    return (Interpolation_function(time,
                                   quantities,
                                   quantity_names,
                                   vertex_coordinates,
                                   triangles,
                                   interpolation_points,
                                   time_thinning=time_thinning,
                                   verbose=verbose,
                                   gauge_neighbour_id=gauge_neighbour_id,
                                   output_centroids=output_centroids),
            starttime)

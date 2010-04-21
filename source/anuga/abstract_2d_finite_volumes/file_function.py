#!/usr/bin/env python
"""File function
Takes a file as input, and returns it as a mathematical function.
For example, you can load an arbitrary 2D heightfield mesh, and treat it as a function as so:

F = file_function('my_mesh.sww', ...)
evaluated_point = F(x, y)

Values will be interpolated across the surface of the mesh. Holes in the mesh have an undefined value.

"""

import numpy as num

from anuga.geospatial_data.geospatial_data import ensure_absolute
from Scientific.IO.NetCDF import NetCDFFile
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a
from anuga.utilities.numerical_tools import ensure_numeric

import anuga.utilities.log as log



##
# @brief Read time history of data from NetCDF file, return callable object.
# @param filename  Name of .sww or .tms file.
# @param domain Associated domain object.
# @param quantities Name of quantity to be interpolated or a list of names.
# @param interpolation_points List of absolute UTM coordinates for points
#                             (N x 2) or geospatial object or
#                             points file name at which values are sought.
# @param time_thinning 
# @param verbose True if this function is to be verbose.
# @param use_cache True means that caching of intermediate result is attempted.
# @param boundary_polygon 
# @param output_centroids if True, data for the centroid of the triangle will be output
# @return A callable object.
def file_function(filename,
                  domain=None,
                  quantities=None,
                  interpolation_points=None,
                  time_thinning=1,
                  time_limit=None,
                  verbose=False,
                  use_cache=False,
                  boundary_polygon=None,
                  output_centroids=False):
    """Read time history of spatial data from NetCDF file and return
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

    time_thinning - 

    verbose - 

    use_cache: True means that caching of intermediate result of
               Interpolation_function is attempted

    boundary_polygon - 

    
    See Interpolation function in anuga.fit_interpolate.interpolation for
    further documentation
    """

    # FIXME (OLE): Should check origin of domain against that of file
    # In fact, this is where origin should be converted to that of domain
    # Also, check that file covers domain fully.

    # Take into account:
    # - domain's georef
    # - sww file's georef
    # - interpolation points as absolute UTM coordinates

    if quantities is None:
        if verbose:
            msg = 'Quantities specified in file_function are None,'
            msg += ' so I will use stage, xmomentum, and ymomentum in that order'
            log.critical(msg)
        quantities = ['stage', 'xmomentum', 'ymomentum']

    # Use domain's startime if available
    if domain is not None:    
        domain_starttime = domain.get_starttime()
    else:
        domain_starttime = None

    # Build arguments and keyword arguments for use with caching or apply.
    args = (filename,)

    # FIXME (Ole): Caching this function will not work well
    # if domain is passed in as instances change hash code.
    # Instead we pass in those attributes that are needed (and return them
    # if modified)
    kwargs = {'quantities': quantities,
              'interpolation_points': interpolation_points,
              'domain_starttime': domain_starttime,
              'time_thinning': time_thinning,      
              'time_limit': time_limit,                                 
              'verbose': verbose,
              'boundary_polygon': boundary_polygon,
              'output_centroids': output_centroids}

    # Call underlying engine with or without caching
    if use_cache is True:
        try:
            from caching import cache
        except:
            msg = 'Caching was requested, but caching module'+\
                  'could not be imported'
            raise msg

        f, starttime = cache(_file_function,
                             args, kwargs,
                             dependencies=[filename],
                             compression=False,                  
                             verbose=verbose)
    else:
        f, starttime = apply(_file_function,
                             args, kwargs)

    #FIXME (Ole): Pass cache arguments, such as compression, in some sort of
    #structure

    f.starttime = starttime
    f.filename = filename
    
    if domain is not None:
        #Update domain.startime if it is *earlier* than starttime from file
        if starttime > domain.starttime:
            msg = 'WARNING: Start time as specified in domain (%f)'\
                  %domain.starttime
            msg += ' is earlier than the starttime of file %s (%f).'\
                     %(filename, starttime)
            msg += ' Modifying domain starttime accordingly.'
            
            if verbose: log.critical(msg)

            domain.set_starttime(starttime) #Modifying model time

            if verbose: log.critical('Domain starttime is now set to %f'
                                     % domain.starttime)
    return f


##
# @brief ??
# @param filename  Name of .sww or .tms file.
# @param domain Associated domain object.
# @param quantities Name of quantity to be interpolated or a list of names.
# @param interpolation_points List of absolute UTM coordinates for points
#                             (N x 2) or geospatial object or
#                             points file name at which values are sought.
# @param time_thinning 
# @param verbose True if this function is to be verbose.
# @param use_cache True means that caching of intermediate result is attempted.
# @param boundary_polygon 
def _file_function(filename,
                   quantities=None,
                   interpolation_points=None,
                   domain_starttime=None,
                   time_thinning=1,
                   time_limit=None,
                   verbose=False,
                   boundary_polygon=None,
                   output_centroids=False):
    """Internal function
    
    See file_function for documentatiton
    """

    assert type(filename) == type(''),\
               'First argument to File_function must be a string'

    try:
        fid = open(filename)
    except Exception, e:
        msg = 'File "%s" could not be opened: Error="%s"' % (filename, e)
        raise msg

    # read first line of file, guess file type
    line = fid.readline()
    fid.close()

    if line[:3] == 'CDF':
        return get_netcdf_file_function(filename,
                                        quantities,
                                        interpolation_points,
                                        domain_starttime,
                                        time_thinning=time_thinning,
                                        time_limit=time_limit,
                                        verbose=verbose,
                                        boundary_polygon=boundary_polygon,
                                        output_centroids=output_centroids)
    else:
        # FIXME (Ole): Could add csv file here to address Ted Rigby's
        # suggestion about reading hydrographs.
        # This may also deal with the gist of ticket:289 
        raise 'Must be a NetCDF File'


##
# @brief ??
# @param filename  Name of .sww or .tms file.
# @param quantity_names Name of quantity to be interpolated or a list of names.
# @param interpolation_points List of absolute UTM coordinates for points
#                             (N x 2) or geospatial object or
#                             points file name at which values are sought.
# @param domain_starttime Start time from domain object.
# @param time_thinning ??
# @param verbose True if this function is to be verbose.
# @param boundary_polygon ??
# @return A callable object.
def get_netcdf_file_function(filename,
                             quantity_names=None,
                             interpolation_points=None,
                             domain_starttime=None,                            
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

    See Interpolation function for further documetation
    """

    # FIXME: Check that model origin is the same as file's origin
    # (both in UTM coordinates)
    # If not - modify those from file to match domain
    # (origin should be passed in)
    # Take this code from e.g. dem2pts in data_manager.py
    # FIXME: Use geo_reference to read and write xllcorner...

    import time, calendar, types
    from anuga.config import time_format

    # Open NetCDF file
    if verbose: log.critical('Reading %s' % filename)

    fid = NetCDFFile(filename, netcdf_mode_r)

    if type(quantity_names) == types.StringType:
        quantity_names = [quantity_names]        

    if quantity_names is None or len(quantity_names) < 1:
        msg = 'No quantities are specified in file_function'
        raise Exception, msg
 
    if interpolation_points is not None:
        interpolation_points = ensure_absolute(interpolation_points)
        msg = 'Points must by N x 2. I got %d' % interpolation_points.shape[1]
        assert interpolation_points.shape[1] == 2, msg

    # Now assert that requested quantitites (and the independent ones)
    # are present in file 
    missing = []
    for quantity in ['time'] + quantity_names:
        if not fid.variables.has_key(quantity):
            missing.append(quantity)

    if len(missing) > 0:
        msg = 'Quantities %s could not be found in file %s'\
              % (str(missing), filename)
        fid.close()
        raise Exception, msg

    # Decide whether this data has a spatial dimension
    spatial = True
    for quantity in ['x', 'y']:
        if not fid.variables.has_key(quantity):
            spatial = False

    if filename[-3:] == 'tms' and spatial is True:
        msg = 'Files of type tms must not contain spatial  information'
        raise msg

    if filename[-3:] == 'sww' and spatial is False:
        msg = 'Files of type sww must contain spatial information'        
        raise msg

    if filename[-3:] == 'sts' and spatial is False:
        #What if mux file only contains one point
        msg = 'Files of type sts must contain spatial information'        
        raise msg

    if filename[-3:] == 'sts' and boundary_polygon is None:
        #What if mux file only contains one point
        msg = 'Files of type sts require boundary polygon'        
        raise msg

    # Get first timestep
    try:
        starttime = fid.starttime[0]
    except ValueError:
        msg = 'Could not read starttime from file %s' % filename
        raise msg

    # Get variables
    # if verbose: log.critical('Get variables'    )
    time = fid.variables['time'][:]    
    # FIXME(Ole): Is time monotoneous?

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
            log.critical('Limited time vector from %.2fs to %.2fs'
                         % (time[-1], time_limit))

    time = time[:upper_time_index]


    
    
    # Get time independent stuff
    if spatial:
        # Get origin
        xllcorner = fid.xllcorner[0]
        yllcorner = fid.yllcorner[0]
        zone = fid.zone[0]        

        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        if filename.endswith('sww'):
            triangles = fid.variables['volumes'][:]

        x = num.reshape(x, (len(x),1))
        y = num.reshape(y, (len(y),1))
        vertex_coordinates = num.concatenate((x,y), axis=1) #m x 2 array

        if boundary_polygon is not None:
            # Remove sts points that do not lie on boundary
            # FIXME(Ole): Why don't we just remove such points from the list of points and associated data?
            # I am actually convinced we can get rid of neighbour_gauge_id altogether as the sts file is produced using the ordering file.
            # All sts points are therefore always present in the boundary. In fact, they *define* parts of the boundary.
            boundary_polygon=ensure_numeric(boundary_polygon)
            boundary_polygon[:,0] -= xllcorner
            boundary_polygon[:,1] -= yllcorner
            temp=[]
            boundary_id=[]
            gauge_id=[]
            for i in range(len(boundary_polygon)):
                for j in range(len(x)):
                    if num.allclose(vertex_coordinates[j],
                                    boundary_polygon[i], 1e-4):
                        #FIXME:
                        #currently gauges lat and long is stored as float and
                        #then cast to double. This cuases slight repositioning
                        #of vertex_coordinates.
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

            if len(num.compress(gauge_neighbour_id>=0,gauge_neighbour_id)) \
               != len(temp)-1:
                msg='incorrect number of segments'
                raise msg
            vertex_coordinates=ensure_numeric(temp)
            if len(vertex_coordinates)==0:
                msg = 'None of the sts gauges fall on the boundary'
                raise msg
        else:
            gauge_neighbour_id=None

        if interpolation_points is not None:
            # Adjust for georef
            interpolation_points[:,0] -= xllcorner
            interpolation_points[:,1] -= yllcorner        
    else:
        gauge_neighbour_id=None
        
    if domain_starttime is not None:
        # If domain_startime is *later* than starttime,
        # move time back - relative to domain's time
        if domain_starttime > starttime:
            time = time - domain_starttime + starttime

        # FIXME Use method in geo to reconcile
        # if spatial:
        # assert domain.geo_reference.xllcorner == xllcorner
        # assert domain.geo_reference.yllcorner == yllcorner
        # assert domain.geo_reference.zone == zone        
        
    if verbose:
        log.critical('File_function data obtained from: %s' % filename)
        log.critical('  References:')
        if spatial:
            log.critical('    Lower left corner: [%f, %f]'
                         % (xllcorner, yllcorner))
        log.critical('    Start time:   %f' % starttime)
        
    
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
        log.critical('Calling interpolation function')
        
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

    # NOTE (Ole): Caching Interpolation function is too slow as
    # the very long parameters need to be hashed.

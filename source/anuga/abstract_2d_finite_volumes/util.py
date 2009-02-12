"""This module contains various auxiliary function used by pyvolution.

It is also a clearing house for functions that may later earn a module
of their own.
"""

import anuga.utilities.polygon
import sys
import os

from os import remove, mkdir, access, F_OK, R_OK, W_OK, sep,getcwd
from os.path import exists, basename, split,join
from warnings import warn
from shutil import copy

from anuga.utilities.numerical_tools import ensure_numeric
from Scientific.IO.NetCDF import NetCDFFile
    
from anuga.geospatial_data.geospatial_data import ensure_absolute
from math import sqrt, atan, degrees

# FIXME (Ole): Temporary short cuts -
# FIXME (Ole): remove and update scripts where they are used
from anuga.utilities.system_tools import get_revision_number
from anuga.utilities.system_tools import store_version_info

from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a

import Numeric as num


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
# @return A callable object.
def file_function(filename,
                  domain=None,
                  quantities=None,
                  interpolation_points=None,
                  time_thinning=1,
                  time_limit=None,
                  verbose=False,
                  use_cache=False,
                  boundary_polygon=None):
    """Read time history of spatial data from NetCDF file and return
    a callable object.

    Input variables:
    
    filename - Name of sww or tms file
       
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
            print msg
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
              'boundary_polygon': boundary_polygon}

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
            
            if verbose: print msg

            domain.set_starttime(starttime) #Modifying model time

            if verbose: print 'Domain starttime is now set to %f'\
               %domain.starttime
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
                   boundary_polygon=None):
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
                                        boundary_polygon=boundary_polygon)
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
                             boundary_polygon=None):
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
    if verbose: print 'Reading', filename

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
    # if verbose: print 'Get variables'    
    time = fid.variables['time'][:]    
    # FIXME(Ole): Is time monotoneous?

    # Apply time limit if requested
    upper_time_index = len(time)    
    msg = 'Time vector obtained from file %s has length 0' % filename
    assert upper_time_index > 0, msg
    
    if time_limit is not None:
        #if verbose is True:
        #    print '****** Time limit', time_limit
        #    print '****** Start time', starttime
        #    print '****** Time in ', time[0], time[-1]
        
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
            print 'Limited time vector from %.2fs to %.2fs'\
                  % (time[-1], time_limit)

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

        #print gauge_neighbour_id
        
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
        print 'File_function data obtained from: %s' %filename
        print '  References:'
        #print '    Datum ....' #FIXME
        if spatial:
            print '    Lower left corner: [%f, %f]'\
                  %(xllcorner, yllcorner)
        print '    Start time:   %f' %starttime                
        
    
    # Produce values for desired data points at
    # each timestep for each quantity
    quantities = {}
    for i, name in enumerate(quantity_names):
        quantities[name] = fid.variables[name][:]
        if boundary_polygon is not None:
            #removes sts points that do not lie on boundary
            quantities[name] = num.take(quantities[name], gauge_id, 1)
            
    # Close sww, tms or sts netcdf file         
    fid.close()

    from anuga.fit_interpolate.interpolate import Interpolation_function

    if not spatial:
        vertex_coordinates = triangles = interpolation_points = None
    if filename[-3:] == 'sts':#added
        triangles = None
        #vertex coordinates is position of urs gauges

    if verbose:
        print 'Calling interpolation function'
        
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
                                   gauge_neighbour_id=gauge_neighbour_id),
            starttime)

    # NOTE (Ole): Caching Interpolation function is too slow as
    # the very long parameters need to be hashed.


##
# @brief Replace multiple substrings in a string.
# @param text The string to operate on.
# @param dictionary A dict containing replacements, key->value.
# @return The new string.
def multiple_replace(text, dictionary):
    """Multiple replace of words in text

    text:       String to be modified
    dictionary: Mapping of words that are to be substituted

    Python Cookbook 3.14 page 88 and page 90
    http://code.activestate.com/recipes/81330/
    """

    import re
    
    #Create a regular expression from all of the dictionary keys
    #matching only entire words
    regex = re.compile(r'\b'+ \
                       r'\b|\b'.join(map(re.escape, dictionary.keys()))+ \
                       r'\b' )

    #For each match, lookup the corresponding value in the dictionary
    return regex.sub(lambda match: dictionary[match.group(0)], text)


##
# @brief Apply arbitrary expressions to the values of a dict.
# @param expression A string expression to apply.
# @param dictionary The dictionary to apply the expression to.
def apply_expression_to_dictionary(expression, dictionary):
    """Apply arbitrary expression to values of dictionary

    Given an expression in terms of the keys, replace key by the
    corresponding values and evaluate.

    expression: Arbitrary, e.g. arithmetric, expression relating keys
                from dictionary. 

    dictionary: Mapping between symbols in expression and objects that
                will be evaluated by expression.
                Values in dictionary must support operators given in
                expression e.g. by overloading

    Due to a limitation with Numeric, this can not evaluate 0/0
    In general, the user can fix by adding 1e-30 to the numerator.
    SciPy core can handle this situation.
    """

    import types
    import re

    assert isinstance(expression, basestring)
    assert type(dictionary) == types.DictType

    #Convert dictionary values to textual representations suitable for eval
    D = {}
    for key in dictionary:
        D[key] = 'dictionary["%s"]' % key

    #Perform substitution of variables    
    expression = multiple_replace(expression, D)

    #Evaluate and return
    try:
        return eval(expression)
    except NameError, e:
        msg = 'Expression "%s" could not be evaluated: %s' % (expression, e)
        raise NameError, msg
    except ValueError, e:
        msg = 'Expression "%s" could not be evaluated: %s' % (expression, e)
        raise ValueError, msg


##
# @brief Format a float into a string.
# @param value Float value to format.
# @param format The format to use (%.2f is default).
# @return The formatted float as a string.
def get_textual_float(value, format = '%.2f'):
    """Get textual representation of floating point numbers
    and accept None as valid entry

    format is a string - default = '%.2f'
    """

    if value is None:
        return 'None'
    else:
        try:
            float(value)
        except:
            # May this is a vector
            if len(value) > 1:
                s = '('
                for v in value:
                    s += get_textual_float(v, format) + ', '
                    
                s = s[:-2] + ')' # Strip trailing comma and close
                return s
            else:
                raise 'Illegal input to get_textual_float:', value
        else:
            return format % float(value)


#################################################################################
# OBSOLETE STUFF
#################################################################################

# @note TEMP
def angle(v1, v2):
    """Temporary Interface to new location"""

    import anuga.utilities.numerical_tools as NT	
    
    msg = 'angle has moved from util.py.  '
    msg += 'Please use "from anuga.utilities.numerical_tools import angle"'
    warn(msg, DeprecationWarning) 

    return NT.angle(v1, v2)

    
# @note TEMP
def anglediff(v0, v1):
    """Temporary Interface to new location"""

    import anuga.utilities.numerical_tools as NT
    
    msg = 'anglediff has moved from util.py.  '
    msg += 'Please use "from anuga.utilities.numerical_tools import anglediff"'
    warn(msg, DeprecationWarning) 

    return NT.anglediff(v0, v1)    

    
# @note TEMP
def mean(x):
    """Temporary Interface to new location"""

    import anuga.utilities.numerical_tools as NT    
    
    msg = 'mean has moved from util.py.  '
    msg += 'Please use "from anuga.utilities.numerical_tools import mean"'
    warn(msg, DeprecationWarning) 

    return NT.mean(x)    


# @note TEMP
def point_on_line(*args, **kwargs):
    """Temporary Interface to new location"""

    msg = 'point_on_line has moved from util.py.  '
    msg += 'Please use "from anuga.utilities.polygon import point_on_line"'
    warn(msg, DeprecationWarning) 

    return utilities.polygon.point_on_line(*args, **kwargs)	

    
# @note TEMP
def inside_polygon(*args, **kwargs):
    """Temporary Interface to new location"""

    print 'inside_polygon has moved from util.py.  ',
    print 'Please use "from anuga.utilities.polygon import inside_polygon"'

    return utilities.polygon.inside_polygon(*args, **kwargs)    

    
# @note TEMP
def outside_polygon(*args, **kwargs):
    """Temporary Interface to new location"""

    print 'outside_polygon has moved from util.py.  ',
    print 'Please use "from anuga.utilities.polygon import outside_polygon"'

    return utilities.polygon.outside_polygon(*args, **kwargs)    


# @note TEMP
def separate_points_by_polygon(*args, **kwargs):
    """Temporary Interface to new location"""

    print 'separate_points_by_polygon has moved from util.py.  ',
    print 'Please use "from anuga.utilities.polygon import ' \
          'separate_points_by_polygon"'

    return utilities.polygon.separate_points_by_polygon(*args, **kwargs)    


# @note TEMP
def read_polygon(*args, **kwargs):
    """Temporary Interface to new location"""

    print 'read_polygon has moved from util.py.  ',
    print 'Please use "from anuga.utilities.polygon import read_polygon"'

    return utilities.polygon.read_polygon(*args, **kwargs)    


# @note TEMP
def populate_polygon(*args, **kwargs):
    """Temporary Interface to new location"""

    print 'populate_polygon has moved from util.py.  ',
    print 'Please use "from anuga.utilities.polygon import populate_polygon"'

    return utilities.polygon.populate_polygon(*args, **kwargs)    


#################################################################################
# End of obsolete stuff ?
#################################################################################

# @note TEMP
def start_screen_catcher(dir_name, myid='', numprocs='', extra_info='',
                         verbose=False):
    """Temporary Interface to new location"""
    from anuga.shallow_water.data_manager import start_screen_catcher \
         as dm_start_screen_catcher

    print 'start_screen_catcher has moved from util.py.  ',
    print 'Please use "from anuga.shallow_water.data_manager import ' \
          'start_screen_catcher"'
    
    return dm_start_screen_catcher(dir_name, myid='', numprocs='',
                                   extra_info='', verbose=False)


##
# @brief Read a .sww file and plot the time series.
# @param swwfiles Dictionary of .sww files.
# @param gauge_filename Name of gauge file.
# @param production_dirs ??
# @param report If True, write figures to report directory.
# @param reportname Name of generated report (if any).
# @param plot_quantity List containing quantities to plot.
# @param generate_fig If True, generate figures as well as CSV files.
# @param surface If True, then generate solution surface with 3d plot.
# @param time_min Beginning of user defined time range for plotting purposes.
# @param time_max End of user defined time range for plotting purposes.
# @param time_thinning ??
# @param time_unit ??
# @param title_on If True, export standard graphics with title.
# @param use_cache If True, use caching.
# @param verbose If True, this function is verbose.
def sww2timeseries(swwfiles,
                   gauge_filename,
                   production_dirs,
                   report=None,
                   reportname=None,
                   plot_quantity=None,
                   generate_fig=False,
                   surface=None,
                   time_min=None,
                   time_max=None,
                   time_thinning=1,                   
                   time_unit=None,
                   title_on=None,
                   use_cache=False,
                   verbose=False):
    """ Read sww file and plot the time series for the
    prescribed quantities at defined gauge locations and
    prescribed time range.

    Input variables:

    swwfiles        - dictionary of sww files with label_ids (used in
                      generating latex output. It will be part of
                      the directory name of file_loc (typically the timestamp).
                      Helps to differentiate latex files for different
                      simulations for a particular scenario.  
                    - assume that all conserved quantities have been stored
                    - assume each sww file has been simulated with same timestep
    
    gauge_filename  - name of file containing gauge data
                        - easting, northing, name , elevation?
                    - OR (this is not yet done)
                        - structure which can be converted to a Numeric array,
                          such as a geospatial data object
                      
    production_dirs -  A list of list, example {20061101_121212: '1 in 10000', 
                                                'boundaries': 'urs boundary'}
                      this will use the second part as the label and the
                      first part as the ?
                      #FIXME: Is it a list or a dictionary
                      # This is probably obsolete by now
                     
    report          - if True, then write figures to report_figures directory in
                      relevant production directory
                    - if False, figures are already stored with sww file
                    - default to False

    reportname      - name for report if wishing to generate report
    
    plot_quantity   - list containing quantities to plot, they must
                      be the name of an existing quantity or one of
                      the following possibilities
                    - possibilities:
                        - stage; 'stage'
                        - depth; 'depth'
                        - speed; calculated as absolute momentum
                         (pointwise) divided by depth; 'speed'
                        - bearing; calculated as the angle of the momentum
                          vector (xmomentum, ymomentum) from the North; 'bearing'
                        - absolute momentum; calculated as
                          sqrt(xmomentum^2 + ymomentum^2); 'momentum'
                        - x momentum; 'xmomentum'
                        - y momentum; 'ymomentum'
                    - default will be ['stage', 'speed', 'bearing']

    generate_fig     - if True, generate figures as well as csv file
                     - if False, csv files created only
                     
    surface          - if True, then generate solution surface with 3d plot
                       and save to current working directory
                     - default = False
    
    time_min         - beginning of user defined time range for plotting purposes
                        - default will be first available time found in swwfile
                        
    time_max         - end of user defined time range for plotting purposes
                        - default will be last available time found in swwfile
                        
    title_on        - if True, export standard graphics with title
                    - if False, export standard graphics without title


    Output:
    
    - time series data stored in .csv for later use if required.
      Name = gauges_timeseries followed by gauge name 
    - latex file will be generated in same directory as where script is
      run (usually production scenario directory.
      Name = latexoutputlabel_id.tex

    Other important information:
    
    It is assumed that the used has stored all the conserved quantities
    and elevation during the scenario run, i.e.
    ['stage', 'elevation', 'xmomentum', 'ymomentum']
    If this has not occurred then sww2timeseries will not work.


    Usage example
    texname = sww2timeseries({project.boundary_name + '.sww': ''},
                             project.polygons_dir + sep + 'boundary_extent.csv',
                             project.anuga_dir, 
                             report = False,
                             plot_quantity = ['stage', 'speed', 'bearing'],
                             time_min = None,
                             time_max = None,
                             title_on = True,   
                             verbose = True)
    
    """

    msg = 'NOTE: A new function is available to create csv files from sww '
    msg += 'files called sww2csv_gauges in anuga.abstract_2d_finite_volumes.util'
    msg += ' PLUS another new function to create graphs from csv files called '
    msg += 'csv2timeseries_graphs in anuga.abstract_2d_finite_volumes.util'
    print msg
    
    k = _sww2timeseries(swwfiles,
                        gauge_filename,
                        production_dirs,
                        report,
                        reportname,
                        plot_quantity,
                        generate_fig,
                        surface,
                        time_min,
                        time_max,
                        time_thinning,                        
                        time_unit,
                        title_on,
                        use_cache,
                        verbose)
    return k


##
# @brief Read a .sww file and plot the time series.
# @param swwfiles Dictionary of .sww files.
# @param gauge_filename Name of gauge file.
# @param production_dirs ??
# @param report If True, write figures to report directory.
# @param reportname Name of generated report (if any).
# @param plot_quantity List containing quantities to plot.
# @param generate_fig If True, generate figures as well as CSV files.
# @param surface If True, then generate solution surface with 3d plot.
# @param time_min Beginning of user defined time range for plotting purposes.
# @param time_max End of user defined time range for plotting purposes.
# @param time_thinning ??
# @param time_unit ??
# @param title_on If True, export standard graphics with title.
# @param use_cache If True, use caching.
# @param verbose If True, this function is verbose.
def _sww2timeseries(swwfiles,
                    gauge_filename,
                    production_dirs,
                    report = None,
                    reportname = None,
                    plot_quantity = None,
                    generate_fig = False,
                    surface = None,
                    time_min = None,
                    time_max = None,
                    time_thinning = 1,                    
                    time_unit = None,
                    title_on = None,
                    use_cache = False,
                    verbose = False):   
        
    # FIXME(Ole): Shouldn't print statements here be governed by verbose?
    assert type(gauge_filename) == type(''), 'Gauge filename must be a string'
    
    try:
        fid = open(gauge_filename)
    except Exception, e:
        msg = 'File "%s" could not be opened: Error="%s"' % (gauge_filename, e)
        raise msg

    if report is None:
        report = False
        
    if plot_quantity is None:
        plot_quantity = ['depth', 'speed']
    else:
        assert type(plot_quantity) == list, 'plot_quantity must be a list'
        check_list(plot_quantity)

    if surface is None:
        surface = False

    if time_unit is None:
        time_unit = 'hours'
    
    if title_on is None:
        title_on = True
    
    if verbose: print '\n Gauges obtained from: %s \n' %gauge_filename

    gauges, locations, elev = get_gauges_from_file(gauge_filename)

    sww_quantity = ['stage', 'elevation', 'xmomentum', 'ymomentum']

    file_loc = []
    f_list = []
    label_id = []
    leg_label = []
    themaxT = 0.0
    theminT = 0.0

    for swwfile in swwfiles.keys():
        try:
            fid = open(swwfile)
        except Exception, e:
            msg = 'File "%s" could not be opened: Error="%s"' % (swwfile, e)
            raise msg

        if verbose:
            print 'swwfile', swwfile

        # Extract parent dir name and use as label
        path, _ = os.path.split(swwfile)
        _, label = os.path.split(path)        
        
        #print 'label', label
        leg_label.append(label)

        f = file_function(swwfile,
                          quantities = sww_quantity,
                          interpolation_points = gauges,
                          time_thinning = time_thinning,
                          verbose = verbose,
                          use_cache = use_cache)

        # determine which gauges are contained in sww file
        count = 0
        gauge_index = []
        for k, g in enumerate(gauges):
            if f(0.0, point_id = k)[2] > 1.0e6:
                count += 1
                if count == 1: print 'Gauges not contained here:'
                print locations[k]
            else:
                gauge_index.append(k)

        if len(gauge_index) > 0:
            print 'Gauges contained here: \n',
        else:
            print 'No gauges contained here. \n'
        for i in range(len(gauge_index)):
             print locations[gauge_index[i]]
             
        index = swwfile.rfind(sep)
        file_loc.append(swwfile[:index+1])
        label_id.append(swwfiles[swwfile])
        
        f_list.append(f)
        maxT = max(f.get_time())
        minT = min(f.get_time())
        if maxT > themaxT: themaxT = maxT
        if minT > theminT: theminT = minT

    if time_min is None:
        time_min = theminT # min(T)
    else:
        if time_min < theminT: # min(T):
            msg = 'Minimum time entered not correct - please try again'
            raise Exception, msg

    if time_max is None:
        time_max = themaxT # max(T)
    else:
        if time_max > themaxT: # max(T):
            msg = 'Maximum time entered not correct - please try again'
            raise Exception, msg

    if verbose and len(gauge_index) > 0:
         print 'Inputs OK - going to generate figures'

    if len(gauge_index) <> 0:
        texfile, elev_output = \
            generate_figures(plot_quantity, file_loc, report, reportname,
                             surface, leg_label, f_list, gauges, locations,
                             elev, gauge_index, production_dirs, time_min,
                             time_max, time_unit, title_on, label_id,
                             generate_fig, verbose)
    else:
        texfile = ''
        elev_output = []

    return texfile, elev_output


##
# @brief Read gauge info from a file.
# @param filename The name of the file to read.
# @return A (gauges, gaugelocation, elev) tuple.
def get_gauges_from_file(filename):
    """ Read in gauge information from file
    """

    from os import sep, getcwd, access, F_OK, mkdir

    # Get data from the gauge file
    fid = open(filename)
    lines = fid.readlines()
    fid.close()
    
    gauges = []
    gaugelocation = []
    elev = []

    # Check header information    
    line1 = lines[0]
    line11 = line1.split(',')

    if isinstance(line11[0], str) is True:
        # We have found text in the first line
        east_index = None
        north_index = None
        name_index = None
        elev_index = None

        for i in range(len(line11)):
            if line11[i].strip().lower() == 'easting':   east_index = i
            if line11[i].strip().lower() == 'northing':  north_index = i
            if line11[i].strip().lower() == 'name':      name_index = i
            if line11[i].strip().lower() == 'elevation': elev_index = i

        if east_index < len(line11) and north_index < len(line11):
            pass
        else:
            msg = 'WARNING: %s does not contain correct header information' \
                  % filename
            msg += 'The header must be: easting, northing, name, elevation'
            raise Exception, msg

        if elev_index is None: 
            raise Exception
    
        if name_index is None: 
            raise Exception

        lines = lines[1:] # Remove header from data
    else:
        # No header, assume that this is a simple easting, northing file

        msg = 'There was no header in file %s and the number of columns is %d' \
              % (filename, len(line11))
        msg += '- was assuming two columns corresponding to Easting and Northing'
        assert len(line11) == 2, msg

        east_index = 0
        north_index = 1

        N = len(lines)
        elev = [-9999]*N
        gaugelocation = range(N)
        
    # Read in gauge data
    for line in lines:
        fields = line.split(',')

        gauges.append([float(fields[east_index]), float(fields[north_index])])

        if len(fields) > 2:
            elev.append(float(fields[elev_index]))
            loc = fields[name_index]
            gaugelocation.append(loc.strip('\n'))

    return gauges, gaugelocation, elev


##
# @brief Check that input quantities in quantity list are legal.
# @param quantity Quantity list to check.
# @note Raises an exception of list is not legal.
def check_list(quantity):
    """ Check that input quantities in quantity list are possible
    """
    import sys
    from sets import Set as set

    all_quantity = ['stage', 'depth', 'momentum', 'xmomentum',
                    'ymomentum', 'speed', 'bearing', 'elevation']

    # convert all quanitiy names to lowercase
    for i,j in enumerate(quantity):
        quantity[i] = quantity[i].lower()

    # check that all names in 'quantity' appear in 'all_quantity'
    p = list(set(quantity).difference(set(all_quantity)))
    if len(p) != 0:
        msg = 'Quantities %s do not exist - please try again' %p
        raise Exception, msg


##
# @brief Calculate velocity bearing from North.
# @param uh ??
# @param vh ??
# @return The calculated bearing.
def calc_bearing(uh, vh):
    """ Calculate velocity bearing from North
    """
    #FIXME (Ole): I reckon we should refactor this one to use
    #             the function angle() in utilities/numerical_tools
    #
    #             It will be a simple matter of
    # * converting from radians to degrees
    # * moving the reference direction from [1,0] to North
    # * changing from counter clockwise to clocwise.
        
    angle = degrees(atan(vh/(uh+1.e-15)))

    if (0 < angle < 90.0):
        if vh > 0:
            bearing = 90.0 - abs(angle)
        if vh < 0:
            bearing = 270.0 - abs(angle)
    
    if (-90 < angle < 0):
        if vh < 0:
            bearing = 90.0 - (angle)
        if vh > 0:
            bearing = 270.0 - (angle)
    if angle == 0: bearing = 0.0

    return bearing


##
# @brief Generate figures from quantities and gauges for each sww file.
# @param plot_quantity  ??
# @param file_loc ??
# @param report ??
# @param reportname ??
# @param surface ??
# @param leg_label ??
# @param f_list ??
# @param gauges ??
# @param locations ??
# @param elev ??
# @param gauge_index ??
# @param production_dirs ??
# @param time_min ??
# @param time_max ??
# @param time_unit ??
# @param title_on ??
# @param label_id ??
# @param generate_fig ??
# @param verbose??
# @return (texfile2, elev_output)
def generate_figures(plot_quantity, file_loc, report, reportname, surface,
                     leg_label, f_list, gauges, locations, elev, gauge_index,
                     production_dirs, time_min, time_max, time_unit,
                     title_on, label_id, generate_fig, verbose):
    """ Generate figures based on required quantities and gauges for
    each sww file
    """
    from os import sep, altsep, getcwd, mkdir, access, F_OK, environ

    if generate_fig is True:
        from pylab import ion, hold, plot, axis, figure, legend, savefig, \
             xlabel, ylabel, title, close, subplot
    
        if surface is True:
            import pylab as p1
            import mpl3d.mplot3d as p3
        
    if report == True:    
        texdir = getcwd()+sep+'report'+sep
        if access(texdir,F_OK) == 0:
            mkdir (texdir)
        if len(label_id) == 1:
            label_id1 = label_id[0].replace(sep,'')
            label_id2 = label_id1.replace('_','')
            texfile = texdir + reportname + '%s' % label_id2
            texfile2 = reportname + '%s' % label_id2
            texfilename = texfile + '.tex'
            fid = open(texfilename, 'w')

            if verbose: print '\n Latex output printed to %s \n' %texfilename
        else:
            texfile = texdir+reportname 
            texfile2 = reportname
            texfilename = texfile + '.tex' 
            fid = open(texfilename, 'w')

            if verbose: print '\n Latex output printed to %s \n' %texfilename
    else:
        texfile = ''
        texfile2 = ''

    p = len(f_list)
    n = []
    n0 = 0
    for i in range(len(f_list)):
        n.append(len(f_list[i].get_time()))
        if n[i] > n0: n0 = n[i]  
    n0 = int(n0)
    m = len(locations)
    model_time = num.zeros((n0, m, p), num.Float) 
    stages = num.zeros((n0, m, p), num.Float)
    elevations = num.zeros((n0, m, p), num.Float) 
    momenta = num.zeros((n0, m, p), num.Float)
    xmom = num.zeros((n0, m, p), num.Float)
    ymom = num.zeros((n0, m, p), num.Float)
    speed = num.zeros((n0, m, p), num.Float)
    bearings = num.zeros((n0, m, p), num.Float)
    due_east = 90.0*num.ones((n0, 1), num.Float)
    due_west = 270.0*num.ones((n0, 1), num.Float)
    depths = num.zeros((n0, m, p), num.Float)
    eastings = num.zeros((n0, m, p), num.Float)
    min_stages = []
    max_stages = []
    min_momentums = []    
    max_momentums = []
    max_xmomentums = []
    max_ymomentums = []
    min_xmomentums = []
    min_ymomentums = []
    max_speeds = []
    min_speeds = []    
    max_depths = []
    model_time_plot3d = num.zeros((n0, m), num.Float)
    stages_plot3d = num.zeros((n0, m), num.Float)
    eastings_plot3d = num.zeros((n0, m),num.Float)
    if time_unit is 'mins': scale = 60.0
    if time_unit is 'hours': scale = 3600.0

    ##### loop over each swwfile #####
    for j, f in enumerate(f_list):
        if verbose: print 'swwfile %d of %d' % (j, len(f_list))

        starttime = f.starttime
        comparefile = file_loc[j] + sep + 'gauges_maxmins' + '.csv'
        fid_compare = open(comparefile, 'w')
        file0 = file_loc[j] + 'gauges_t0.csv'
        fid_0 = open(file0, 'w')

        ##### loop over each gauge #####
        for k in gauge_index:
            if verbose: print 'Gauge %d of %d' % (k, len(gauges))

            g = gauges[k]
            min_stage = 10
            max_stage = 0
            max_momentum = max_xmomentum = max_ymomentum = 0
            min_momentum = min_xmomentum = min_ymomentum = 100
            max_speed = 0
            min_speed = 0            
            max_depth = 0            
            gaugeloc = str(locations[k])
            thisfile = file_loc[j] + sep + 'gauges_time_series' + '_' \
                       + gaugeloc + '.csv'
            if j == 0:
                fid_out = open(thisfile, 'w')
                s = 'Time, Stage, Momentum, Speed, Elevation, xmom, ymom, Bearing \n'
                fid_out.write(s)            

            #### generate quantities #######
            for i, t in enumerate(f.get_time()):
                if time_min <= t <= time_max:
                    w = f(t, point_id = k)[0]
                    z = f(t, point_id = k)[1]
                    uh = f(t, point_id = k)[2]
                    vh = f(t, point_id = k)[3]
                    depth = w-z      
                    m = sqrt(uh*uh + vh*vh)
                    if depth < 0.001:
                        vel = 0.0
                    else:
                        vel = m / (depth + 1.e-6/depth) 
                    bearing = calc_bearing(uh, vh)                    
                    model_time[i,k,j] = (t + starttime)/scale #t/60.0
                    stages[i,k,j] = w
                    elevations[i,k,j] = z 
                    xmom[i,k,j] = uh 
                    ymom[i,k,j] = vh 
                    momenta[i,k,j] = m 
                    speed[i,k,j] = vel 
                    bearings[i,k,j] = bearing 
                    depths[i,k,j] = depth
                    thisgauge = gauges[k]
                    eastings[i,k,j] = thisgauge[0]
                    s = '%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f,\n' \
                            % (t, w, m, vel, z, uh, vh, bearing)
                    fid_out.write(s)
                    if t == 0:
                        s = '%.2f, %.2f, %.2f\n' % (g[0], g[1], w)
                        fid_0.write(s)
                    if t/60.0 <= 13920: tindex = i
                    if w > max_stage: max_stage = w
                    if w < min_stage: min_stage = w
                    if m > max_momentum: max_momentum = m
                    if m < min_momentum: min_momentum = m                    
                    if uh > max_xmomentum: max_xmomentum = uh
                    if vh > max_ymomentum: max_ymomentum = vh
                    if uh < min_xmomentum: min_xmomentum = uh
                    if vh < min_ymomentum: min_ymomentum = vh
                    if vel > max_speed: max_speed = vel
                    if vel < min_speed: min_speed = vel                    
                    if z > 0 and depth > max_depth: max_depth = depth
                    
                    
            s = '%.2f, %.2f, %.2f, %.2f, %s\n' \
                    % (max_stage, min_stage, z, thisgauge[0], leg_label[j])
            fid_compare.write(s)
            max_stages.append(max_stage)
            min_stages.append(min_stage)
            max_momentums.append(max_momentum)
            max_xmomentums.append(max_xmomentum)
            max_ymomentums.append(max_ymomentum)
            min_xmomentums.append(min_xmomentum)
            min_ymomentums.append(min_ymomentum)
            min_momentums.append(min_momentum)            
            max_depths.append(max_depth)
            max_speeds.append(max_speed)
            min_speeds.append(min_speed)            
            #### finished generating quantities for each swwfile #####
        
        model_time_plot3d[:,:] = model_time[:,:,j]
        stages_plot3d[:,:] = stages[:,:,j]
        eastings_plot3d[:,] = eastings[:,:,j]
            
        if surface is True:
            print 'Printing surface figure'
            for i in range(2):
                fig = p1.figure(10)
                ax = p3.Axes3D(fig)
                if len(gauges) > 80:
                    ax.plot_surface(model_time[:,:,j],
                                    eastings[:,:,j],
                                    stages[:,:,j])
                else:
                    ax.plot3D(num.ravel(eastings[:,:,j]),
                              num.ravel(model_time[:,:,j]),
                              num.ravel(stages[:,:,j]))
                ax.set_xlabel('time')
                ax.set_ylabel('x')
                ax.set_zlabel('stage')
                fig.add_axes(ax)
                p1.show()
                surfacefig = 'solution_surface%s' % leg_label[j]
                p1.savefig(surfacefig)
                p1.close()
            
    #### finished generating quantities for all swwfiles #####

    # x profile for given time
    if surface is True:
        figure(11)
        plot(eastings[tindex,:,j], stages[tindex,:,j])
        xlabel('x')
        ylabel('stage')
        profilefig = 'solution_xprofile' 
        savefig('profilefig')

    elev_output = []
    if generate_fig is True:
        depth_axis = axis([starttime/scale, time_max/scale, -0.1,
                           max(max_depths)*1.1])
        stage_axis = axis([starttime/scale, time_max/scale,
                           min(min_stages), max(max_stages)*1.1])
        vel_axis = axis([starttime/scale, time_max/scale,
                         min(min_speeds), max(max_speeds)*1.1])
        mom_axis = axis([starttime/scale, time_max/scale,
                         min(min_momentums), max(max_momentums)*1.1])
        xmom_axis = axis([starttime/scale, time_max/scale,
                          min(min_xmomentums), max(max_xmomentums)*1.1])
        ymom_axis = axis([starttime/scale, time_max/scale,
                          min(min_ymomentums), max(max_ymomentums)*1.1])
        cstr = ['g', 'r', 'b', 'c', 'm', 'y', 'k']
        nn = len(plot_quantity)
        no_cols = 2
        
        if len(label_id) > 1: graphname_report = []
        pp = 1
        div = 11.
        cc = 0
        for k in gauge_index:
            g = gauges[k]
            count1 = 0
            if report == True and len(label_id) > 1:
                s = '\\begin{figure}[ht] \n' \
                    '\\centering \n' \
                    '\\begin{tabular}{cc} \n'
                fid.write(s)
            if len(label_id) > 1: graphname_report = []

            #### generate figures for each gauge ####
            for j, f in enumerate(f_list):
                ion()
                hold(True)
                count = 0
                where1 = 0
                where2 = 0
                word_quantity = ''
                if report == True and len(label_id) == 1:
                    s = '\\begin{figure}[hbt] \n' \
                        '\\centering \n' \
                        '\\begin{tabular}{cc} \n'
                    fid.write(s)
                    
                for which_quantity in plot_quantity:
                    count += 1
                    where1 += 1
                    figure(count, frameon = False)
                    if which_quantity == 'depth':
                        plot(model_time[0:n[j]-1,k,j],
                             depths[0:n[j]-1,k,j], '-', c = cstr[j])
                        units = 'm'
                        axis(depth_axis)
                    if which_quantity == 'stage':
                        if elevations[0,k,j] <= 0:
                            plot(model_time[0:n[j]-1,k,j],
                                 stages[0:n[j]-1,k,j], '-', c = cstr[j])
                            axis(stage_axis)
                        else:
                            plot(model_time[0:n[j]-1,k,j],
                                 depths[0:n[j]-1,k,j], '-', c = cstr[j])
                            #axis(depth_axis)                 
                        units = 'm'
                    if which_quantity == 'momentum':
                        plot(model_time[0:n[j]-1,k,j],
                             momenta[0:n[j]-1,k,j], '-', c = cstr[j])
                        axis(mom_axis)
                        units = 'm^2 / sec'
                    if which_quantity == 'xmomentum':
                        plot(model_time[0:n[j]-1,k,j],
                             xmom[0:n[j]-1,k,j], '-', c = cstr[j])
                        axis(xmom_axis)
                        units = 'm^2 / sec'
                    if which_quantity == 'ymomentum':
                        plot(model_time[0:n[j]-1,k,j],
                             ymom[0:n[j]-1,k,j], '-', c = cstr[j])
                        axis(ymom_axis)
                        units = 'm^2 / sec'
                    if which_quantity == 'speed':
                        plot(model_time[0:n[j]-1,k,j],
                             speed[0:n[j]-1,k,j], '-', c = cstr[j])
                        axis(vel_axis)
                        units = 'm / sec'
                    if which_quantity == 'bearing':
                        plot(model_time[0:n[j]-1,k,j],bearings[0:n[j]-1,k,j],'-',
                             model_time[0:n[j]-1,k,j], due_west[0:n[j]-1], '-.', 
                             model_time[0:n[j]-1,k,j], due_east[0:n[j]-1], '-.')
                        units = 'degrees from North'
                        #ax = axis([time_min, time_max, 0.0, 360.0])
                        legend(('Bearing','West','East'))

                    if time_unit is 'mins': xlabel('time (mins)')
                    if time_unit is 'hours': xlabel('time (hours)')
                    #if which_quantity == 'stage' \
                    #   and elevations[0:n[j]-1,k,j] > 0:
                    #    ylabel('%s (%s)' %('depth', units))
                    #else:
                    #    ylabel('%s (%s)' %(which_quantity, units))
                        #ylabel('%s (%s)' %('wave height', units))
                    ylabel('%s (%s)' %(which_quantity, units))
                    if len(label_id) > 1: legend((leg_label),loc='upper right')

                    #gaugeloc1 = gaugeloc.replace(' ','')
                    #gaugeloc2 = gaugeloc1.replace('_','')
                    gaugeloc2 = str(locations[k]).replace(' ','')
                    graphname = '%sgauge%s_%s' %(file_loc[j],
                                                 gaugeloc2,
                                                 which_quantity)

                    if report == True and len(label_id) > 1:
                        figdir = getcwd()+sep+'report_figures'+sep
                        if access(figdir,F_OK) == 0 :
                            mkdir (figdir)
                        latex_file_loc = figdir.replace(sep,altsep) 
                        # storing files in production directory    
                        graphname_latex = '%sgauge%s%s' \
                                          % (latex_file_loc, gaugeloc2,
                                             which_quantity)
                        # giving location in latex output file
                        graphname_report_input = '%sgauge%s%s' % \
                                                 ('..' + altsep + 
                                                      'report_figures' + altsep,
                                                  gaugeloc2, which_quantity)
                        graphname_report.append(graphname_report_input)
                        
                        # save figures in production directory for report
                        savefig(graphname_latex)

                    if report == True:
                        figdir = getcwd() + sep + 'report_figures' + sep
                        if access(figdir,F_OK) == 0:
                            mkdir(figdir)
                        latex_file_loc = figdir.replace(sep,altsep)    

                        if len(label_id) == 1: 
                            # storing files in production directory  
                            graphname_latex = '%sgauge%s%s%s' % \
                                              (latex_file_loc, gaugeloc2,
                                               which_quantity, label_id2)
                            # giving location in latex output file
                            graphname_report = '%sgauge%s%s%s' % \
                                               ('..' + altsep +
                                                    'report_figures' + altsep,
                                                gaugeloc2, which_quantity,
                                                label_id2)
                            s = '\includegraphics' \
                                '[width=0.49\linewidth, height=50mm]{%s%s}' % \
                                (graphname_report, '.png')
                            fid.write(s)
                            if where1 % 2 == 0:
                                s = '\\\\ \n'
                                where1 = 0
                            else:
                                s = '& \n'
                            fid.write(s)
                            savefig(graphname_latex)
                    
                    if title_on == True:
                        title('%s scenario: %s at %s gauge' % \
                              (label_id, which_quantity, gaugeloc2))
                        #title('Gauge %s (MOST elevation %.2f, ' \
                        #      'ANUGA elevation %.2f)' % \
                        #      (gaugeloc2, elevations[10,k,0],
                        #       elevations[10,k,1]))

                    savefig(graphname) # save figures with sww file

                if report == True and len(label_id) == 1:
                    for i in range(nn-1):
                        if nn > 2:
                            if plot_quantity[i] == 'stage' \
                               and elevations[0,k,j] > 0:
                                word_quantity += 'depth' + ', '
                            else:
                                word_quantity += plot_quantity[i] + ', '
                        else:
                            if plot_quantity[i] == 'stage' \
                               and elevations[0,k,j] > 0:
                                word_quantity += 'depth' + ', '
                            else:
                                word_quantity += plot_quantity[i]
                        
                    if plot_quantity[nn-1] == 'stage' and elevations[0,k,j] > 0:
                        word_quantity += ' and ' + 'depth'
                    else:
                        word_quantity += ' and ' + plot_quantity[nn-1]
                    caption = 'Time series for %s at %s location ' \
                              '(elevation %.2fm)' % \
                              (word_quantity, locations[k], elev[k])
                    if elev[k] == 0.0:
                        caption = 'Time series for %s at %s location ' \
                                  '(elevation %.2fm)' % \
                                  (word_quantity, locations[k],
                                   elevations[0,k,j])
                        east = gauges[0]
                        north = gauges[1]
                        elev_output.append([locations[k], east, north,
                                            elevations[0,k,j]])
                    label = '%sgauge%s' % (label_id2, gaugeloc2)
                    s = '\end{tabular} \n' \
                        '\\caption{%s} \n' \
                        '\label{fig:%s} \n' \
                        '\end{figure} \n \n' % (caption, label)
                    fid.write(s)
                    cc += 1
                    if cc % 6 == 0: fid.write('\\clearpage \n')
                    savefig(graphname_latex)               
                    
            if report == True and len(label_id) > 1:
                for i in range(nn-1):
                    if nn > 2:
                        if plot_quantity[i] == 'stage' and elevations[0,k,j] > 0:
                            word_quantity += 'depth' + ','
                        else:
                            word_quantity += plot_quantity[i] + ', '
                    else:
                        if plot_quantity[i] == 'stage' and elevations[0,k,j] > 0:
                            word_quantity += 'depth'
                        else:
                            word_quantity += plot_quantity[i]
                    where1 = 0
                    count1 += 1
                    index = j*len(plot_quantity)
                    for which_quantity in plot_quantity:
                        where1 += 1
                        s = '\includegraphics' \
                            '[width=0.49\linewidth, height=50mm]{%s%s}' % \
                            (graphname_report[index], '.png')
                        index += 1
                        fid.write(s)
                        if where1 % 2 == 0:
                            s = '\\\\ \n'
                            where1 = 0
                        else:
                            s = '& \n'
                        fid.write(s)
                word_quantity += ' and ' + plot_quantity[nn-1]            
                label = 'gauge%s' %(gaugeloc2) 
                caption = 'Time series for %s at %s location ' \
                          '(elevation %.2fm)' % \
                          (word_quantity, locations[k], elev[k])
                if elev[k] == 0.0:
                        caption = 'Time series for %s at %s location ' \
                                  '(elevation %.2fm)' % \
                                  (word_quantity, locations[k],
                                   elevations[0,k,j])
                        thisgauge = gauges[k]
                        east = thisgauge[0]
                        north = thisgauge[1]
                        elev_output.append([locations[k], east, north,
                                            elevations[0,k,j]])
                        
                s = '\end{tabular} \n' \
                    '\\caption{%s} \n' \
                    '\label{fig:%s} \n' \
                    '\end{figure} \n \n' % (caption, label)
                fid.write(s)
                if float((k+1)/div - pp) == 0.:
                    fid.write('\\clearpage \n')
                    pp += 1
                #### finished generating figures ###

            close('all')
        
    return texfile2, elev_output


# FIXME (DSG): Add unit test, make general, not just 2 files,
# but any number of files.
##
# @brief ??
# @param dir_name ??
# @param filename1 ??
# @param filename2 ??
# @return ??
# @note TEMP
def copy_code_files(dir_name, filename1, filename2):
    """Temporary Interface to new location"""

    from anuga.shallow_water.data_manager import \
                    copy_code_files as dm_copy_code_files
    print 'copy_code_files has moved from util.py.  ',
    print 'Please use "from anuga.shallow_water.data_manager \
                                        import copy_code_files"'
    
    return dm_copy_code_files(dir_name, filename1, filename2)


##
# @brief Create a nested sub-directory path.
# @param root_directory The base diretory path.
# @param directories An iterable of sub-directory names.
# @return The final joined directory path.
# @note If each sub-directory doesn't exist, it will be created.
def add_directories(root_directory, directories):
    """
    Add the first sub-directory in 'directories' to root_directory.
    Then add the second sub-directory to the accumulating path and so on.

    Return the path of the final directory.

    This is handy for specifying and creating a directory where data will go.
    """
    dir = root_directory
    for new_dir in directories:
        dir = os.path.join(dir, new_dir)
        if not access(dir,F_OK):
            mkdir(dir)
    return dir


##
# @brief 
# @param filename 
# @param separator_value 
# @return 
# @note TEMP
def get_data_from_file(filename, separator_value=','):
    """Temporary Interface to new location"""
    from anuga.shallow_water.data_manager import \
                        get_data_from_file as dm_get_data_from_file
    print 'get_data_from_file has moved from util.py'
    print 'Please use "from anuga.shallow_water.data_manager \
                                     import get_data_from_file"'
    
    return dm_get_data_from_file(filename,separator_value = ',')


##
# @brief 
# @param verbose 
# @param kwargs 
# @return 
# @note TEMP
def store_parameters(verbose=False,**kwargs):
    """Temporary Interface to new location"""
    
    from anuga.shallow_water.data_manager \
                    import store_parameters as dm_store_parameters
    print 'store_parameters has moved from util.py.'
    print 'Please use "from anuga.shallow_water.data_manager \
                                     import store_parameters"'
    
    return dm_store_parameters(verbose=False,**kwargs)


##
# @brief Remove vertices that are not associated with any triangle.
# @param verts An iterable (or array) of points.
# @param triangles An iterable of 3 element tuples.
# @param number_of_full_nodes ??
# @return (verts, triangles) where 'verts' has been updated.
def remove_lone_verts(verts, triangles, number_of_full_nodes=None):
    """Removes vertices that are not associated with any triangles.

    verts is a list/array of points.
    triangles is a list of 3 element tuples.  Each tuple represents a triangle.
    number_of_full_nodes relate to parallelism when a mesh has an
        extra layer of ghost points.
    """

    verts = ensure_numeric(verts)
    triangles = ensure_numeric(triangles)
    
    N = len(verts)
    
    # initialise the array to easily find the index of the first loner
    # ie, if N=3 -> [6,5,4]
    loners=num.arange(2*N, N, -1)
    for t in triangles:
        for vert in t:
            loners[vert]= vert # all non-loners will have loners[i]=i 

    lone_start = 2*N - max(loners) # The index of the first loner

    if lone_start-1 == N:
        # no loners
        pass
    elif min(loners[lone_start:N]) > N:
        # All the loners are at the end of the vert array
        verts = verts[0:lone_start]
    else:
        # change the loners list so it can be used to modify triangles
        # Remove the loners from verts
        # Could've used X=compress(less(loners,N),loners)
        # verts=num.take(verts,X)  to Remove the loners from verts
        # but I think it would use more memory
        new_i = lone_start	# point at first loner - 'shuffle down' target
        for i in range(lone_start, N):
            if loners[i] >= N:	# [i] is a loner, leave alone
                pass
            else:		# a non-loner, move down
                loners[i] = new_i
                verts[new_i] = verts[i]
                new_i += 1
        verts = verts[0:new_i]

        # Modify the triangles
        #print "loners", loners
        #print "triangles before", triangles
        triangles = num.choose(triangles,loners)
        #print "triangles after", triangles
    return verts, triangles


##
# @brief Compute centroid values from vertex values
# @param x Values at vertices of triangular mesh.
# @param triangles Nx3 integer array pointing to vertex information.
# @return [N] array of centroid values.
def get_centroid_values(x, triangles):
    """Compute centroid values from vertex values
    
    x: Values at vertices of triangular mesh
    triangles: Nx3 integer array pointing to vertex information
    for each of the N triangels. Elements of triangles are
    indices into x
    """
        
    xc = num.zeros(triangles.shape[0], num.Float) # Space for centroid info
    
    for k in range(triangles.shape[0]):
        # Indices of vertices
        i0 = triangles[k][0]
        i1 = triangles[k][1]
        i2 = triangles[k][2]        
        
        xc[k] = (x[i0] + x[i1] + x[i2])/3

    return xc


# @note TEMP
def make_plots_from_csv_file(directories_dic={dir:['gauge', 0, 0]},
                                output_dir='',
                                base_name='',
                                plot_numbers=['3-5'],
                                quantities=['speed','stage','momentum'],
                                assess_all_csv_files=True,
                                extra_plot_name='test'):

    msg = 'make_plots_from_csv_file has been replaced by csv2timeseries_graphs '
    msg += 'Please use "from anuga.abstract_2d_finite_volumes.util import ' \
           'csv2timeseries_graphs"'
    raise Exception, msg

    return csv2timeseries_graphs(directories_dic,
                                 output_dir,
                                 base_name,
                                 plot_numbers,
                                 quantities,
                                 extra_plot_name,
                                 assess_all_csv_files)


##
# @brief Plot time series from CSV files.
# @param directories_dic 
# @param output_dir 
# @param base_name 
# @param plot_numbers 
# @param quantities 
# @param extra_plot_name 
# @param assess_all_csv_files 
# @param create_latex 
# @param verbose 
# @note Assumes that 'elevation' is in the CSV file(s).
def csv2timeseries_graphs(directories_dic={},
                          output_dir='',
                          base_name=None,
                          plot_numbers='',
                          quantities=['stage'],
                          extra_plot_name='',
                          assess_all_csv_files=True,
                          create_latex=False,
                          verbose=False):
                                
    """
    Read in csv files that have the right header information and
    plot time series such as Stage, Speed, etc. Will also plot several
    time series on one plot. Filenames must follow this convention,
    <base_name><plot_number>.csv eg gauge_timeseries3.csv
    
    NOTE: relies that 'elevation' is in the csv file!

    Each file represents a location and within each file there are
    time, quantity columns.
    
    For example:    
    if "directories_dic" defines 4 directories and in each directories
    there is a csv files corresponding to the right "plot_numbers", 
    this will create a plot with 4 lines one for each directory AND 
    one plot for each "quantities".  ??? FIXME: unclear.
    
    Usage:
        csv2timeseries_graphs(directories_dic={'slide'+sep:['Slide',0, 0],
                                       'fixed_wave'+sep:['Fixed Wave',0,0]},
                            output_dir='fixed_wave'+sep,
                            base_name='gauge_timeseries_',
                            plot_numbers='',
                            quantities=['stage','speed'],
                            extra_plot_name='',
                            assess_all_csv_files=True,                            
                            create_latex=False,
                            verbose=True)
            this will create one plot for stage with both 'slide' and 
            'fixed_wave' lines on it for stage and speed for each csv
            file with 'gauge_timeseries_' as the prefix. The graghs 
            will be in the output directory 'fixed_wave' and the graph
            axis will be determined by assessing all the 
    
    ANOTHER EXAMPLE
        new_csv2timeseries_graphs(directories_dic={'slide'+sep:['Slide',0, 0],
                                       'fixed_wave'+sep:['Fixed Wave',0,0]},
                            output_dir='fixed_wave'+sep,
                            base_name='gauge_timeseries_',
                            plot_numbers=['1-3'],
                            quantities=['stage','speed'],
                            extra_plot_name='',
                            assess_all_csv_files=False,                            
                            create_latex=False,
                            verbose=True)
        This will plot csv files called gauge_timeseries_1.csv and 
        gauge_timeseries3.csv from both 'slide' and 'fixed_wave' directories
        to 'fixed_wave'. There will be 4 plots created two speed and two stage
        one for each csv file. There will be two lines on each of these plots.
        And the axis will have been determined from only these files, had 
        assess_all_csv_files = True all csv file with 'gauges_timeseries_' prefix
        would of been assessed.
    
    ANOTHER EXAMPLE    
         csv2timeseries_graphs({'J:'+sep+'anuga_validation'+sep:['new',20,-.1],
                                   'J:'+sep+'conical_island'+sep:['test',0,0]},
                                   output_dir='',
                                   plot_numbers=['1','3'],
                                   quantities=['stage','depth','bearing'],
                                   base_name='gauge_b',
                                   assess_all_csv_files=True,
                                  verbose=True)    
        
            This will produce one plot for each quantity (therefore 3) in the
            current directory, each plot will have 2 lines on them. The first
            plot named 'new' will have the time offseted by 20secs and the stage
            height adjusted by -0.1m
        
    Inputs:
        directories_dic: dictionary of directory with values (plot 
                         legend name for directory), (start time of 
                         the time series) and the (value to add to 
                         stage if needed). For example
                         {dir1:['Anuga_ons',5000, 0],
                          dir2:['b_emoth',5000,1.5],
                          dir3:['b_ons',5000,1.5]}
                         Having multiple directories defined will plot them on 
                         one plot, therefore there will be 3 lines on each of
                         these plot. If you only want one line per plot call
                         csv2timeseries_graph separately for each directory,
                         eg only have one directory in the 'directories_dic' in
                         each call. 
                         
        output_dir: directory for the plot outputs. Only important to define when
                    you have more than one directory in your directories_dic, if
                    you have not defined it and you have multiple directories in
                    'directories_dic' there will be plots in each directory,
                    however only one directory will contain the complete
                    plot/graphs.
        
        base_name: Is used a couple of times.
                   1) to find the csv files to be plotted if there is no
                      'plot_numbers' then csv files with 'base_name' are plotted
                   2) in the title of the plots, the length of base_name is 
                      removed from the front of the filename to be used in the
                      title. 
                   This could be changed if needed. 
                   Note is ignored if assess_all_csv_files=True
        
        plot_numbers: a String list of numbers to plot. For example 
                      [0-4,10,15-17] will read and attempt to plot
                      the follow 0,1,2,3,4,10,15,16,17
                      NOTE: if no plot numbers this will create one plot per
                            quantity, per gauge

        quantities: Will get available quantities from the header in the csv
                    file.  Quantities must be one of these.
                    NOTE: ALL QUANTITY NAMES MUST BE lower case!
                    
        extra_plot_name: A string that is appended to the end of the 
                         output filename.
                    
        assess_all_csv_files: if true it will read ALL csv file with
                             "base_name", regardless of 'plot_numbers'
                              and determine a uniform set of axes for 
                              Stage, Speed and Momentum. IF FALSE it 
                              will only read the csv file within the
                             'plot_numbers'
                             
        create_latex: NOT IMPLEMENTED YET!! sorry Jane....
        
    OUTPUTS: saves the plots to 
              <output_dir><base_name><plot_number><extra_plot_name>.png
    """

    try: 
        import pylab
    except ImportError:
        msg='csv2timeseries_graphs needs pylab to be installed correctly'
        raise msg
            #ANUGA don't need pylab to work so the system doesn't 
            #rely on pylab being installed 
        return

    from os import sep
    from anuga.shallow_water.data_manager import \
                               get_all_files_with_extension, csv2dict

    seconds_in_hour = 3600
    seconds_in_minutes = 60
    
    quantities_label={}
#    quantities_label['time'] = 'time (hours)'
    quantities_label['time'] = 'time (minutes)'
    quantities_label['stage'] = 'wave height (m)'
    quantities_label['speed'] = 'speed (m/s)'
    quantities_label['momentum'] = 'momentum (m^2/sec)'
    quantities_label['depth'] = 'water depth (m)'
    quantities_label['xmomentum'] = 'momentum (m^2/sec)'
    quantities_label['ymomentum'] = 'momentum (m^2/sec)'
    quantities_label['bearing'] = 'degrees (o)'
    quantities_label['elevation'] = 'elevation (m)'
    
    if extra_plot_name != '':
        extra_plot_name = '_' + extra_plot_name

    new_plot_numbers=[]
    #change plot_numbers to list, eg ['0-4','10'] 
    #to ['0','1','2','3','4','10']
    for i, num_string in enumerate(plot_numbers):
        if '-' in num_string: 
            start = int(num_string[:num_string.rfind('-')])
            end = int(num_string[num_string.rfind('-') + 1:]) + 1
            for x in range(start, end):
                new_plot_numbers.append(str(x))
        else:
            new_plot_numbers.append(num_string)

    #finds all the files that fit the specs provided and return a list of them
    #so to help find a uniform max and min for the plots... 
    list_filenames=[]
    all_csv_filenames=[]
    if verbose: print 'Determining files to access for axes ranges.'
    
    for i,directory in enumerate(directories_dic.keys()):
        all_csv_filenames.append(get_all_files_with_extension(directory,
                                                              base_name, '.csv'))

        filenames=[]
        if plot_numbers == '': 
            list_filenames.append(get_all_files_with_extension(directory,
                                                               base_name,'.csv'))
        else:
            for number in new_plot_numbers:
                filenames.append(base_name + number)
            list_filenames.append(filenames)

    #use all the files to get the values for the plot axis
    max_start_time= -1000.
    min_start_time = 100000 
    
    if verbose: print 'Determining uniform axes' 

    #this entire loop is to determine the min and max range for the 
    #axes of the plots

#    quantities.insert(0,'elevation')
    quantities.insert(0,'time')

    directory_quantity_value={}
#    quantity_value={}
    min_quantity_value={}
    max_quantity_value={}

    for i, directory in enumerate(directories_dic.keys()):
        filename_quantity_value = {}
        if assess_all_csv_files == False:
            which_csv_to_assess = list_filenames[i]
        else:
            #gets list of filenames for directory "i"
            which_csv_to_assess = all_csv_filenames[i]
        
        for j, filename in enumerate(which_csv_to_assess):
            quantity_value = {}

            dir_filename = join(directory,filename)
            attribute_dic, title_index_dic = csv2dict(dir_filename + '.csv')
            directory_start_time = directories_dic[directory][1]
            directory_add_tide = directories_dic[directory][2]

            if verbose: print 'reading: %s.csv' %dir_filename

            #add time to get values
            for k, quantity in enumerate(quantities):
                quantity_value[quantity] = [float(x) for
                                                x in attribute_dic[quantity]]

                #add tide to stage if provided
                if quantity == 'stage':
                     quantity_value[quantity] = num.array(quantity_value[quantity], num.Float) \
                                                          + directory_add_tide

                #condition to find max and mins for all the plots
                # populate the list with something when i=0 and j=0 and
                # then compare to the other values to determine abs max and min
                if i==0 and j==0:
                    min_quantity_value[quantity], \
                        max_quantity_value[quantity] = \
                            get_min_max_values(quantity_value[quantity])

                    if quantity != 'time':
                        min_quantity_value[quantity] = \
                            min_quantity_value[quantity] *1.1
                        max_quantity_value[quantity] = \
                            max_quantity_value[quantity] *1.1
                else:
                    min, max = get_min_max_values(quantity_value[quantity])
                
                    # min and max are multipled by "1+increase_axis" to get axes
                    # that are slighty bigger than the max and mins
                    # so the plots look good.

                    increase_axis = (max-min)*0.05
                    if min <= min_quantity_value[quantity]:
                        if quantity == 'time': 
                            min_quantity_value[quantity] = min
                        else:
                            if round(min,2) == 0.00:
                                min_quantity_value[quantity] = -increase_axis
#                                min_quantity_value[quantity] = -2.
                                #min_quantity_value[quantity] = \
                                #    -max_quantity_value[quantity]*increase_axis
                            else:
#                                min_quantity_value[quantity] = \
#                                    min*(1+increase_axis)
                                min_quantity_value[quantity]=min-increase_axis
                    
                    if max > max_quantity_value[quantity]: 
                        if quantity == 'time': 
                            max_quantity_value[quantity] = max
                        else:
                            max_quantity_value[quantity] = max + increase_axis
#                            max_quantity_value[quantity]=max*(1+increase_axis)

            #set the time... ???
            if min_start_time > directory_start_time: 
                min_start_time = directory_start_time
            if max_start_time < directory_start_time: 
                max_start_time = directory_start_time
            
            filename_quantity_value[filename]=quantity_value
            
        directory_quantity_value[directory]=filename_quantity_value
    
    #final step to unifrom axis for the graphs
    quantities_axis={}
    
    for i, quantity in enumerate(quantities):
        quantities_axis[quantity] = (float(min_start_time) \
                                         / float(seconds_in_minutes),
                                     (float(max_quantity_value['time']) \
                                          + float(max_start_time)) \
                                              / float(seconds_in_minutes),
                                     min_quantity_value[quantity],
                                     max_quantity_value[quantity])

        if verbose and (quantity != 'time' and quantity != 'elevation'): 
            print 'axis for quantity %s are x:(%s to %s)%s and y:(%s to %s)%s' \
                  % (quantity, 
                     quantities_axis[quantity][0],
                     quantities_axis[quantity][1],
                     quantities_label['time'],
                     quantities_axis[quantity][2],
                     quantities_axis[quantity][3],
                     quantities_label[quantity])

    cstr = ['b', 'r', 'g', 'c', 'm', 'y', 'k']

    if verbose: print 'Now start to plot \n'
    
    i_max = len(directories_dic.keys())
    legend_list_dic = {}
    legend_list = []
    for i, directory in enumerate(directories_dic.keys()):
        if verbose: print 'Plotting in %s %s' % (directory, new_plot_numbers)

        # FIXME THIS SORT IS VERY IMPORTANT
        # Without it the assigned plot numbers may not work correctly
        # there must be a better way
        list_filenames[i].sort()
        for j, filename in enumerate(list_filenames[i]):
            if verbose: print 'Starting %s' % filename  

            directory_name = directories_dic[directory][0]
            directory_start_time = directories_dic[directory][1]
            directory_add_tide = directories_dic[directory][2]
            
            # create an if about the start time and tide height if don't exist
            attribute_dic, title_index_dic = csv2dict(directory + sep
                                                      + filename + '.csv')
            #get data from dict in to list
            #do maths to list by changing to array
            t = (num.array(directory_quantity_value[directory][filename]['time'])
                     + directory_start_time) / seconds_in_minutes

            #finds the maximum elevation, used only as a test
            # and as info in the graphs
            max_ele=-100000
            min_ele=100000
            elevation = [float(x) for x in attribute_dic["elevation"]]
            
            min_ele, max_ele = get_min_max_values(elevation)
            
            if min_ele != max_ele:
                print "Note! Elevation changes in %s" %dir_filename

            # creates a dictionary with keys that is the filename and attributes
            # are a list of lists containing 'directory_name' and 'elevation'.
            # This is used to make the contents for the legends in the graphs,
            # this is the name of the model and the elevation.  All in this
            # great one liner from DG. If the key 'filename' doesn't exist it
            # creates the entry if the entry exist it appends to the key.

            legend_list_dic.setdefault(filename,[]) \
                .append([directory_name, round(max_ele, 3)])

            # creates a LIST for the legend on the last iteration of the
            # directories which is when "legend_list_dic" has been fully
            # populated. Creates a list of strings which is used in the legend
            # only runs on the last iteration for all the gauges(csv) files
            # empties the list before creating it 

            if i == i_max - 1:
                legend_list = []
    
                for name_and_elevation in legend_list_dic[filename]:
                    legend_list.append('%s (elevation = %sm)'\
                                       % (name_and_elevation[0],
                                          name_and_elevation[1]))
            
            #skip time and elevation so it is not plotted!
            for k, quantity in enumerate(quantities):
                if quantity != 'time' and quantity != 'elevation':
                    pylab.figure(int(k*100+j))
                    pylab.ylabel(quantities_label[quantity])
                    pylab.plot(t,
                               directory_quantity_value[directory]\
                                                       [filename][quantity],
                               c = cstr[i], linewidth=1)
                    pylab.xlabel(quantities_label['time'])
                    pylab.axis(quantities_axis[quantity])
                    pylab.legend(legend_list,loc='upper right')
                    
                    pylab.title('%s at %s gauge'
                                % (quantity, filename[len(base_name):]))

                    if output_dir == '':
                        figname = '%s%s%s_%s%s.png' \
                                  % (directory, sep, filename, quantity,
                                     extra_plot_name)
                    else:
                        figname = '%s%s%s_%s%s.png' \
                                  % (output_dir, sep, filename, quantity,
                                     extra_plot_name)

                    if verbose: print 'saving figure here %s' %figname

                    pylab.savefig(figname)
           
    if verbose: print 'Closing all plots'

    pylab.close('all')
    del pylab

    if verbose: print 'Finished closing plots'

##
# @brief Return min and max of an iterable.
# @param list The iterable to return min & max of.
# @return (min, max) of 'list'.
def get_min_max_values(list=None):
    """ 
    Returns the min and max of the list it was provided.
    """

    if list == None: print 'List must be provided'
        
    return min(list), max(list)


##
# @brief Get runup around a point in a CSV file.
# @param gauge_filename gauge file name.
# @param sww_filename SWW file name.
# @param runup_filename Name of file to report into.
# @param size ??
# @param verbose ??
def get_runup_data_for_locations_from_file(gauge_filename,
                                           sww_filename,
                                           runup_filename,
                                           size=10,
                                           verbose=False):
    """this will read a csv file with the header x,y. Then look in a square
    'size'x2 around this position for the 'max_inundaiton_height' in the
    'sww_filename' and report the findings in the 'runup_filename'.
    
    WARNING: NO TESTS! 
    """

    from anuga.shallow_water.data_manager import get_all_directories_with_name,\
                                                 get_maximum_inundation_data,\
                                                 csv2dict
                                                 
    file = open(runup_filename, "w")
    file.write("easting,northing,runup \n ")
    file.close()
    
    #read gauge csv file to dictionary
    attribute_dic, title_index_dic = csv2dict(gauge_filename)
    northing = [float(x) for x in attribute_dic["y"]]
    easting = [float(x) for x in attribute_dic["x"]]

    print 'Reading %s' %sww_filename

    runup_locations=[]
    for i, x in enumerate(northing):
        poly = [[int(easting[i]+size),int(northing[i]+size)],
                [int(easting[i]+size),int(northing[i]-size)],
                [int(easting[i]-size),int(northing[i]-size)],
                [int(easting[i]-size),int(northing[i]+size)]]
        
        run_up, x_y = get_maximum_inundation_data(filename=sww_filename,
                                                  polygon=poly,
                                                  verbose=False) 

        #if no runup will return 0 instead of NONE
        if run_up==None: run_up=0
        if x_y==None: x_y=[0,0]
        
        if verbose:
            print 'maximum inundation runup near %s is %s meters' %(x_y,run_up)
        
        #writes to file
        file = open(runup_filename, "a")
        temp = '%s,%s,%s \n' % (x_y[0], x_y[1], run_up)
        file.write(temp)
        file.close()


##
# @brief ??
# @param  ??
# @param gauge_file ??
# @param out_name ??
# @param quantities ??
# @param verbose ??
# @param use_cache ??
def sww2csv_gauges(sww_file,
                   gauge_file,
                   out_name='gauge_',
                   quantities=['stage', 'depth', 'elevation',
                               'xmomentum', 'ymomentum'],
                   verbose=False,
                   use_cache=True):
    """
    
    Inputs: 
        
        NOTE: if using csv2timeseries_graphs after creating csv file,
        it is essential to export quantities 'depth' and 'elevation'.
        'depth' is good to analyse gauges on land and elevation is used
        automatically by csv2timeseries_graphs in the legend.
        
        sww_file: path to any sww file
        
        gauge_file: Assumes that it follows this format
            name, easting, northing, elevation
            point1, 100.3, 50.2, 10.0
            point2, 10.3, 70.3, 78.0
        
        NOTE: order of column can change but names eg 'easting', 'elevation' 
        must be the same! ALL lowercaps!

        out_name: prefix for output file name (default is 'gauge_')
        
    Outputs: 
        one file for each gauge/point location in the points file. They
        will be named with this format in the same directory as the 'sww_file'
            <out_name><name>.csv
        eg gauge_point1.csv if <out_name> not supplied
           myfile_2_point1.csv if <out_name> ='myfile_2_'
            
            
        They will all have a header
    
    Usage: sww2csv_gauges(sww_file='test1.sww',
                          quantities = ['stage', 'elevation','depth','bearing'],
                          gauge_file='gauge.txt')    
    
    Interpolate the quantities at a given set of locations, given
    an sww file.
    The results are written to a csv file.

    In the future let points be a points file.
    And the user choose the quantities.

    This is currently quite specific.
    If it needs to be more general, change things.

    This is really returning speed, not velocity.
    """
    
    from csv import reader,writer
    from anuga.utilities.numerical_tools import ensure_numeric, mean, NAN
    import string
    from anuga.shallow_water.data_manager import get_all_swwfiles

#    quantities =  ['stage', 'elevation', 'xmomentum', 'ymomentum']
    #print "points",points 

    assert type(gauge_file) == type(''), 'Gauge filename must be a string'
    assert type(out_name) == type(''), 'Output filename prefix must be a string'
    
    try:
        point_reader = reader(file(gauge_file))
    except Exception, e:
        msg = 'File "%s" could not be opened: Error="%s"' % (gauge_file, e)
        raise msg

    if verbose: print '\n Gauges obtained from: %s \n' %gauge_file
    
    point_reader = reader(file(gauge_file))
    points = []
    point_name = []
    
    #read point info from file
    for i,row in enumerate(point_reader):
        #read header and determine the column numbers to read correcty.
        if i==0:
            for j,value in enumerate(row):
                if value.strip()=='easting':easting=j
                if value.strip()=='northing':northing=j
                if value.strip()=='name':name=j
                if value.strip()=='elevation':elevation=j
        else:
            points.append([float(row[easting]),float(row[northing])])
            point_name.append(row[name])
        
    #convert to array for file_function
    points_array = num.array(points,num.Float)
        
    points_array = ensure_absolute(points_array)

    dir_name, base = os.path.split(sww_file)    

    #need to get current directory so when path and file
    #are "joined" below the directory is correct
    if dir_name == '':
        dir_name =getcwd()
        
    if access(sww_file,R_OK):
        if verbose: print 'File %s exists' %(sww_file)
    else:
        msg = 'File "%s" could not be opened: no read permission' % sww_file
        raise msg

    sww_files = get_all_swwfiles(look_in_dir=dir_name,
                                 base_name=base,
                                 verbose=verbose)
    print 'sww files just after get_all_swwfiles()', sww_files
    # fudge to get SWW files in 'correct' order, oldest on the left
    sww_files.sort()

    if verbose:
        print 'sww files', sww_files
    
    #to make all the quantities lower case for file_function
    quantities = [quantity.lower() for quantity in quantities]

    # what is quantities are needed from sww file to calculate output quantities
    # also 

    core_quantities = ['stage', 'elevation', 'xmomentum', 'ymomentum']

    gauge_file = out_name

    heading = [quantity for quantity in quantities]
    heading.insert(0,'time')
    heading.insert(1,'hours')

    #create a list of csv writers for all the points and write header
    points_writer = []
    for point_i,point in enumerate(points):
        points_writer.append(writer(file(dir_name + sep + gauge_file
                                         + point_name[point_i] + '.csv', "wb")))
        points_writer[point_i].writerow(heading)
    
    if verbose: print 'Writing csv files'

    quake_offset_time = None

    for sww_file in sww_files:
        sww_file = join(dir_name, sww_file+'.sww')
        #print 'sww file = ',sww_file
        callable_sww = file_function(sww_file,
                                     quantities=core_quantities,
                                     interpolation_points=points_array,
                                     verbose=verbose,
                                     use_cache=use_cache)

        if quake_offset_time is None:
            quake_offset_time = callable_sww.starttime

        for time in callable_sww.get_time():
            print 'time = ', str(time)
            for point_i, point in enumerate(points_array):
               # print 'gauge_file = ', str(point_name[point_i])
                #print 'point_i = ', str(point_i), ' point is = ', str(point) 
                #add domain starttime to relative time.
                quake_time = time + quake_offset_time
                points_list = [quake_time, quake_time/3600.]# fudge around SWW time bug
                #print 'point list = ', str(points_list)
                point_quantities = callable_sww(time,point_i)
                #print 'point quantities = ', str(point_quantities)
                
                for quantity in quantities:
                    if quantity == NAN:
                        print 'quantity does not exist in' % callable_sww.get_name
                    else:
                        if quantity == 'stage':
                            points_list.append(point_quantities[0])
                            
                        if quantity == 'elevation':
                            points_list.append(point_quantities[1])
                            
                        if quantity == 'xmomentum':
                            points_list.append(point_quantities[2])
                            
                        if quantity == 'ymomentum':
                            points_list.append(point_quantities[3])
                            
                        if quantity == 'depth':
                            points_list.append(point_quantities[0] 
                                               - point_quantities[1])

                        if quantity == 'momentum':
                            momentum = sqrt(point_quantities[2]**2 
                                            + point_quantities[3]**2)
                            points_list.append(momentum)
                            
                        if quantity == 'speed':
                            #if depth is less than 0.001 then speed = 0.0
                            if point_quantities[0] - point_quantities[1] < 0.001:
                                vel = 0.0
                            else:
                                if point_quantities[2] < 1.0e6:
                                    momentum = sqrt(point_quantities[2]**2
                                                    + point_quantities[3]**2)
        #                            vel = momentum/depth              
                                    vel = momentum / (point_quantities[0] 
                                                      - point_quantities[1])
        #                            vel = momentum/(depth + 1.e-6/depth)
                                else:
                                    momentum = 0
                                    vel = 0
                                
                            points_list.append(vel)
                            
                        if quantity == 'bearing':
                            points_list.append(calc_bearing(point_quantities[2],
                                                            point_quantities[3]))

                print 'point list before write (writer %s) = %s' % (str(point_name[point_i]), str(points_list))
                points_writer[point_i].writerow(points_list)
            

##
# @brief Get a wave height at a certain depth given wave height at another depth.
# @param d1 The first depth.
# @param d2 The second depth.
# @param h1 Wave ampitude at d1
# @param verbose True if this function is to be verbose.
# @return The wave height at d2.
def greens_law(d1, d2, h1, verbose=False):
    """Green's Law

    Green's Law allows an approximation of wave amplitude at
    a given depth based on the fourh root of the ratio of two depths
    and the amplitude at another given depth.

    Note, wave amplitude is equal to stage.
    
    Inputs:

    d1, d2 - the two depths
    h1     - the wave amplitude at d1
    h2     - the derived amplitude at d2

    h2 = h1 (d1/d2)^(1/4), where d2 cannot equal 0.
    
    """

    d1 = ensure_numeric(d1)
    d2 = ensure_numeric(d2)
    h1 = ensure_numeric(h1)

    if d1 <= 0.0:
        msg = 'the first depth, d1 (%f), must be strictly positive' % (d1)
        raise Exception(msg)

    if d2 <= 0.0:
        msg = 'the second depth, d2 (%f), must be strictly positive' % (d2)
        raise Exception(msg)
    
    if h1 <= 0.0:
        msg = 'the wave amplitude, h1 (%f), must be strictly positive' % (h1)
        raise Exception(msg)
    
    h2 = h1*(d1/d2)**0.25

    assert h2 > 0
    
    return h2
        

##
# @brief Get the square-root of a value.
# @param s The value to get the square-root of.
# @return The square-root of 's'.
def square_root(s):
    return sqrt(s)



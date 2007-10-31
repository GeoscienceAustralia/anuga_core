"""This module contains various auxiliary function used by pyvolution.

It is also a clearing house for functions that may later earn a module
of their own.
"""

import anuga.utilities.polygon
import sys
import os

from os import remove, mkdir, access, F_OK, W_OK, sep,mkdir
from os.path import exists, basename, split,join
from warnings import warn
from shutil import copy

from anuga.utilities.numerical_tools import ensure_numeric
from Numeric import arange, choose, zeros, Float, array
    
from anuga.geospatial_data.geospatial_data import ensure_absolute

def file_function(filename,
                  domain=None,
                  quantities=None,
                  interpolation_points=None,
                  time_thinning=1,
                  verbose=False,
                  use_cache=False):
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
                 If quantities are None, domain's conserved quantities
                 are used.

    interpolation_points - list of absolute UTM coordinates for points (N x 2)
    or geospatial object or points file name at which values are sought
    
    use_cache: True means that caching of intermediate result of
               Interpolation_function is attempted

    
    See Interpolation function for further documentation
    """


    #FIXME (OLE): Should check origin of domain against that of file
    #In fact, this is where origin should be converted to that of domain
    #Also, check that file covers domain fully.

    #Take into account:
    #- domain's georef
    #- sww file's georef
    #- interpolation points as absolute UTM coordinates




    # Use domain's conserved_quantity names as defaults
    if domain is not None:    
        if quantities is None: 
            quantities = domain.conserved_quantities
            
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
              'verbose': verbose}


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



def _file_function(filename,
                   quantities=None,
                   interpolation_points=None,
                   domain_starttime=None,
                   time_thinning=1,
                   verbose=False):
    """Internal function
    
    See file_function for documentatiton
    """
    

    assert type(filename) == type(''),\
               'First argument to File_function must be a string'

    try:
        fid = open(filename)
    except Exception, e:
        msg = 'File "%s" could not be opened: Error="%s"'\
                  %(filename, e)
        raise msg

    line = fid.readline()
    fid.close()


    if line[:3] == 'CDF':
        return get_netcdf_file_function(filename,
                                        quantities,
                                        interpolation_points,
                                        domain_starttime,
                                        time_thinning=time_thinning,
                                        verbose=verbose)
    else:
        raise 'Must be a NetCDF File'



def get_netcdf_file_function(filename,
                             quantity_names=None,
                             interpolation_points=None,
                             domain_starttime=None,                            
                             time_thinning=1,                             
                             verbose=False):
    """Read time history of spatial data from NetCDF sww file and
    return a callable object f(t,x,y)
    which will return interpolated values based on the input file.

    Model time (domain_starttime)
    will be checked, possibly modified and returned
    
    All times are assumed to be in UTC

    See Interpolation function for further documetation
    
    """
    
    
    #FIXME: Check that model origin is the same as file's origin
    #(both in UTM coordinates)
    #If not - modify those from file to match domain
    #(origin should be passed in)
    #Take this code from e.g. dem2pts in data_manager.py
    #FIXME: Use geo_reference to read and write xllcorner...
        

    import time, calendar, types
    from anuga.config import time_format
    from Scientific.IO.NetCDF import NetCDFFile
    from Numeric import array, zeros, Float, alltrue, concatenate, reshape

    #Open NetCDF file
    if verbose: print 'Reading', filename
    fid = NetCDFFile(filename, 'r')

    if type(quantity_names) == types.StringType:
        quantity_names = [quantity_names]        
    
    if quantity_names is None or len(quantity_names) < 1:
        #If no quantities are specified get quantities from file
        #x, y, time are assumed as the independent variables so
        #they are excluded as potentiol quantities
        quantity_names = []
        for name in fid.variables.keys():
            if name not in ['x', 'y', 'time']:
                quantity_names.append(name)

    if len(quantity_names) < 1:                
        msg = 'ERROR: At least one independent value must be specified'
        raise msg


    if interpolation_points is not None:
        interpolation_points = ensure_absolute(interpolation_points)
        msg = 'Points must by N x 2. I got %d' %interpolation_points.shape[1]
        assert interpolation_points.shape[1] == 2, msg


    #Now assert that requested quantitites (and the independent ones)
    #are present in file 
    missing = []
    for quantity in ['time'] + quantity_names:
        if not fid.variables.has_key(quantity):
            missing.append(quantity)

    if len(missing) > 0:
        msg = 'Quantities %s could not be found in file %s'\
              %(str(missing), filename)
        fid.close()
        raise Exception, msg

    #Decide whether this data has a spatial dimension
    spatial = True
    for quantity in ['x', 'y']:
        if not fid.variables.has_key(quantity):
            spatial = False

    if filename[-3:] == 'tms' and spatial is True:
        msg = 'Files of type tms must not contain spatial information'
        raise msg

    if filename[-3:] == 'sww' and spatial is False:
        msg = 'Files of type sww must contain spatial information'        
        raise msg        

    #Get first timestep
    try:
        starttime = fid.starttime[0]
    except ValueError:
        msg = 'Could not read starttime from file %s' %filename
        raise msg

    #Get variables
    #if verbose: print 'Get variables'    
    time = fid.variables['time'][:]    

    # Get time independent stuff
    if spatial:
        #Get origin
        xllcorner = fid.xllcorner[0]
        yllcorner = fid.yllcorner[0]
        zone = fid.zone[0]        

        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        triangles = fid.variables['volumes'][:]

        x = reshape(x, (len(x),1))
        y = reshape(y, (len(y),1))
        vertex_coordinates = concatenate((x,y), axis=1) #m x 2 array

        if interpolation_points is not None:
            #Adjust for georef
            interpolation_points[:,0] -= xllcorner
            interpolation_points[:,1] -= yllcorner        
        



    if domain_starttime is not None:

        #If domain_startime is *later* than starttime,
        #move time back - relative to domain's time
        if domain_starttime > starttime:
            time = time - domain_starttime + starttime

        #FIXME Use method in geo to reconcile
        #if spatial:
        #assert domain.geo_reference.xllcorner == xllcorner
        #assert domain.geo_reference.yllcorner == yllcorner
        #assert domain.geo_reference.zone == zone        
        
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
    fid.close()
    

    #from least_squares import Interpolation_function
    from anuga.fit_interpolate.interpolate import Interpolation_function

    if not spatial:
        vertex_coordinates = triangles = interpolation_points = None         


    # Return Interpolation_function instance as well as
    # starttime for use to possible modify that of domain
    return Interpolation_function(time,
                                  quantities,
                                  quantity_names,
                                  vertex_coordinates,
                                  triangles,
                                  interpolation_points,
                                  time_thinning=time_thinning,
                                  verbose=verbose), starttime

    # NOTE (Ole): Caching Interpolation function is too slow as
    # the very long parameters need to be hashed.






def multiple_replace(text, dictionary):
    """Multiple replace of words in text

    text:       String to be modified
    dictionary: Mapping of words that are to be substituted

    Python Cookbook 3.14 page 88 and page 90
    """

    import re
    
    #Create a regular expression from all of the dictionary keys
    #matching only entire words
    regex = re.compile(r'\b'+ \
                       r'\b|\b'.join(map(re.escape, dictionary.keys()))+ \
                       r'\b' )

    #For each match, lookup the corresponding value in the dictionary
    return regex.sub(lambda match: dictionary[match.group(0)], text)




def apply_expression_to_dictionary(expression, dictionary):#dictionary):
    """Apply arbitrary expression to values of dictionary

    Given an expression in terms of the keys, replace key by the
    corresponding values and evaluate.
   

    expression: Arbitrary, e.g. arithmetric, expression relating keys
                from dictionary. 

    dictionary: Mapping between symbols in expression and objects that
                will be evaluated by expression.
                Values in dictionary must support operators given in
                expression e.g. by overloading

    due to a limitation with Numeric, this can not evaluate 0/0
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
        D[key] = 'dictionary["%s"]' %key

    #Perform substitution of variables    
    expression = multiple_replace(expression, D)

    #Evaluate and return
    try:
        return eval(expression)
    except NameError, e:
        msg = 'Expression "%s" could not be evaluated: %s' %(expression, e)
        raise NameError, msg
    except ValueError, e:
        msg = 'Expression "%s" could not be evaluated: %s' %(expression, e)
        raise ValueError, msg
    

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
            return format %float(value)



####################################
####OBSOLETE STUFF


def angle(v1, v2):
    """Temporary Interface to new location"""

    import anuga.utilities.numerical_tools as NT	
    
    msg = 'angle has moved from util.py.  '
    msg += 'Please use "from anuga.utilities.numerical_tools import angle"'
    warn(msg, DeprecationWarning) 

    return NT.angle(v1, v2)
    
def anglediff(v0, v1):
    """Temporary Interface to new location"""

    import anuga.utilities.numerical_tools as NT
    
    msg = 'anglediff has moved from util.py.  '
    msg += 'Please use "from anuga.utilities.numerical_tools import anglediff"'
    warn(msg, DeprecationWarning) 

    return NT.anglediff(v0, v1)    

    
def mean(x):
    """Temporary Interface to new location"""

    import anuga.utilities.numerical_tools as NT    
    
    msg = 'mean has moved from util.py.  '
    msg += 'Please use "from anuga.utilities.numerical_tools import mean"'
    warn(msg, DeprecationWarning) 

    return NT.mean(x)    

def point_on_line(*args, **kwargs):
    """Temporary Interface to new location"""

    msg = 'point_on_line has moved from util.py.  '
    msg += 'Please use "from anuga.utilities.polygon import point_on_line"'
    warn(msg, DeprecationWarning) 

    return utilities.polygon.point_on_line(*args, **kwargs)	
    
def inside_polygon(*args, **kwargs):
    """Temporary Interface to new location"""

    print 'inside_polygon has moved from util.py.  ',
    print 'Please use "from anuga.utilities.polygon import inside_polygon"'

    return utilities.polygon.inside_polygon(*args, **kwargs)    
    
def outside_polygon(*args, **kwargs):
    """Temporary Interface to new location"""

    print 'outside_polygon has moved from util.py.  ',
    print 'Please use "from anuga.utilities.polygon import outside_polygon"'

    return utilities.polygon.outside_polygon(*args, **kwargs)    


def separate_points_by_polygon(*args, **kwargs):
    """Temporary Interface to new location"""

    print 'separate_points_by_polygon has moved from util.py.  ',
    print 'Please use "from anuga.utilities.polygon import separate_points_by_polygon"'

    return utilities.polygon.separate_points_by_polygon(*args, **kwargs)    



def read_polygon(*args, **kwargs):
    """Temporary Interface to new location"""

    print 'read_polygon has moved from util.py.  ',
    print 'Please use "from anuga.utilities.polygon import read_polygon"'

    return utilities.polygon.read_polygon(*args, **kwargs)    


def populate_polygon(*args, **kwargs):
    """Temporary Interface to new location"""

    print 'populate_polygon has moved from util.py.  ',
    print 'Please use "from anuga.utilities.polygon import populate_polygon"'

    return utilities.polygon.populate_polygon(*args, **kwargs)    

##################### end of obsolete stuff ? ############

def start_screen_catcher(dir_name, myid='', numprocs='', extra_info='',
                         print_to_screen=False, verbose=False):
    """Temporary Interface to new location"""
    from anuga.shallow_water.data_manager import start_screen_catcher as dm_start_screen_catcher

    print 'start_screen_catcher has moved from util.py.  ',
    print 'Please use "from anuga.shallow_water.data_manager import start_screen_catcher"'
    
    return dm_start_screen_catcher(dir_name, myid='', numprocs='', extra_info='',
                         print_to_screen=False, verbose=False)

def get_revision_number():
    """Get the version number of the SVN
    NOTE: This requires that the command svn is on the system PATH
    (simply aliasing svn to the binary will not work)
    """

    # Create dummy info 
    #info = 'Revision: Version info could not be obtained.'
    #info += 'A command line version of svn must be availbable '
    #info += 'on the system PATH, access to the subversion '
    #info += 'repository is necessary and the output must '
    #info += 'contain a line starting with "Revision:"'
    

    #FIXME (Ole): Change this so that svn info is attempted first.
    # If that fails, try to read a stored file with that same info (this would be created by e.g. the release script). Failing that, throw an exception.

    #FIXME (Ole): Move this and store_version_info to utilities


    try:
        from anuga.stored_version_info import version_info
    except:
	msg = 'No version info stored and command "svn" is not '
	msg += 'recognised on the system PATH.\n\n'
	msg += 'If ANUGA has been installed from a distribution e.g. as '
	msg += 'obtained from SourceForge,\n'
	msg += 'the version info should be '
	msg += 'available in the automatically generated file '
	msg += 'stored_version_info.py\n'
	msg += 'in the anuga root directory.\n'
	msg += 'If run from a Subversion sandpit, '
	msg += 'ANUGA will try to obtain the version info '
	msg += 'by using the command: "svn info".\n'
	msg += 'In this case, make sure svn is accessible on the system path. '
	msg += 'Simply aliasing svn to the binary will not work. '
	msg += 'Good luck!'

        # No file available - try using Subversion
        try:
            # The null stuff is so this section fails quitly.
            # This could cause the svn info command to fail due to
            # the redirection being bad on some platforms.
            # If that occurs then change this code.
            if sys.platform[0:3] == 'win':
                fid = os.popen('svn info 2> null')
            else:
                fid = os.popen('svn info 2>/dev/null')
	
        except:
            raise Exception(msg)
        else:
            #print 'Got version from svn'            
            version_info = fid.read()
	    
	    if version_info == '':
	        raise Exception(msg)    
    else:
        pass
        #print 'Got version from file'

            
    for line in version_info.split('\n'):
        if line.startswith('Revision:'):
            break

    fields = line.split(':')
    msg = 'Keyword "Revision" was not found anywhere in text: %s' %version_info
    assert fields[0].startswith('Revision'), msg            

    try:
        revision_number = int(fields[1])
    except:
        msg = 'Revision number must be an integer. I got %s' %fields[1]
        msg += 'Check that the command svn is on the system path' 
        raise Exception(msg)                
        
    return revision_number


def store_version_info(destination_path='.', verbose=False):
    """Obtain current version from Subversion and store it.
    
    Title: store_version_info()

    Author: Ole Nielsen (Ole.Nielsen@ga.gov.au)

    CreationDate: January 2006

    Description:
        This function obtains current version from Subversion and stores it
        is a Python file named 'stored_version_info.py' for use with
        get_version_info()

        If svn is not available on the system PATH, an Exception is thrown
    """

    # Note (Ole): This function should not be unit tested as it will only
    # work when running out of the sandpit. End users downloading the
    # ANUGA distribution would see a failure.
    #
    # FIXME: This function should really only be used by developers (
    # (e.g. for creating new ANUGA releases), so maybe it should move
    # to somewhere else.
    
    import config

    try:
        fid = os.popen('svn info')
    except:
        msg = 'Command "svn" is not recognised on the system PATH'
        raise Exception(msg)
    else:    
        txt = fid.read()
        fid.close()


        # Determine absolute filename
        if destination_path[-1] != os.sep:
            destination_path += os.sep
            
        filename = destination_path + config.version_filename

        fid = open(filename, 'w')

        docstring = 'Stored version info.\n\n'
        docstring += 'This file provides the version for distributions '
        docstring += 'that are not accessing Subversion directly.\n'
        docstring += 'The file is automatically generated and should not '
        docstring += 'be modified manually.\n'
        fid.write('"""%s"""\n\n' %docstring)
        
        fid.write('version_info = """\n%s"""' %txt)
        fid.close()


        if verbose is True:
            print 'Version info stored to %s' %filename
            
    
def sww2timeseries(swwfiles,
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
    
    """ Read sww file and plot the time series for the
    prescribed quantities at defined gauge locations and
    prescribed time range.

    Input variables:

    swwfiles        - dictionary of sww files with label_ids (used in
                      generating latex output. It will be part of
                      the directory name of file_loc (typically the timestamp).
                      Helps to differentiate latex files for different simulations
                      for a particular scenario.  
                    - assume that all conserved quantities have been stored
                    - assume each sww file has been simulated with same timestep
    
    gauge_filename  - name of file containing gauge data
                        - easting, northing, name , elevation?
                    - OR (this is not yet done)
                        - structure which can be converted to a Numeric array,
                          such as a geospatial data object
                      
    production_dirs -  A list of list, example {20061101_121212: '1 in 10000', 
                                                'boundaries': 'urs boundary'}
                      this will use the second part as the label and the first part 
                      as the ?
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

    generate_fig     - if True, generate figures as well as csv files of quantities
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
        
    assert type(gauge_filename) == type(''),\
           'Gauge filename must be a string'
    
    try:
        fid = open(gauge_filename)
    except Exception, e:
        msg = 'File "%s" could not be opened: Error="%s"'\
                  %(gauge_filename, e)
        raise msg

    if report is None:
        report = False
        
    if plot_quantity is None:
        plot_quantity = ['depth', 'speed']
    else:
        assert type(plot_quantity) == list,\
               'plot_quantity must be a list'
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
            msg = 'File "%s" could not be opened: Error="%s"'\
                  %(swwfile, e)
            raise msg

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

    if verbose and len(gauge_index) > 0: print 'Inputs OK - going to generate figures'

    if len(gauge_index) <> 0:
        texfile, elev_output = generate_figures(plot_quantity, file_loc, report, reportname, surface,
                                                leg_label, f_list, gauges, locations, elev, gauge_index,
                                                production_dirs, time_min, time_max, time_unit,
                                                title_on, label_id, generate_fig, verbose)
    else:
        texfile = ''
        elev_output = []

    return texfile, elev_output
               
def get_gauges_from_file(filename):
    """ Read in gauge information from file
    """
    from os import sep, getcwd, access, F_OK, mkdir
    fid = open(filename)
    lines = fid.readlines()
    fid.close()
    
    gauges = []
    gaugelocation = []
    elev = []

    # Check header information    
    line1 = lines[0]
    line11 = line1.split(',')

    if isinstance(line11[0],str) is True:
        # We have found text in the first line
        east_index = None
        north_index = None
        name_index = None
        elev_index = None
        for i in range(len(line11)):
            if line11[i].strip('\n').strip('\r').strip(' ').lower() == 'easting': east_index = i
            if line11[i].strip('\n').strip('\r').strip(' ').lower() == 'northing': north_index = i
            if line11[i].strip('\n').strip('\r').strip(' ').lower() == 'name': name_index = i
            if line11[i].strip('\n').strip('\r').strip(' ').lower() == 'elevation': elev_index = i

        if east_index < len(line11) and north_index < len(line11):
            pass
        else:
            msg = 'WARNING: %s does not contain correct header information' %(filename)
            msg += 'The header must be: easting, northing, name, elevation'
            raise Exception, msg

        if elev_index is None: 
            raise Exception
    
        if name_index is None: 
            raise Exception

        lines = lines[1:] # Remove header from data
    else:
        # No header, assume that this is a simple easting, northing file

        msg = 'There was no header in file %s and the number of columns is %d' %(filename, len(line11))
        msg += '- I was assuming two columns corresponding to Easting and Northing'
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



def check_list(quantity):
    """ Check that input quantities in quantity list are possible
    """
    all_quantity = ['stage', 'depth', 'momentum', 'xmomentum',
                    'ymomentum', 'speed', 'bearing', 'elevation']

	
		    
    import sys
    #if not sys.version.startswith('2.4'):
        # Backwards compatibility
    #   from sets import Set as set

    from sets import Set as set
            
    for i,j in enumerate(quantity):
        quantity[i] = quantity[i].lower()
    p = list(set(quantity).difference(set(all_quantity)))
    if len(p) <> 0:
        msg = 'Quantities %s do not exist - please try again' %p
        raise Exception, msg
        
    return 

def calc_bearing(uh, vh):
    """ Calculate velocity bearing from North
    """
    from math import atan, degrees
    
    angle = degrees(atan(vh/(uh+1.e-15)))
    if (0 < angle < 90.0):
        if vh > 0:
            bearing = 90.0 - abs(angle)
        if vh < 0:
            bearing = 270.0 - abs(angle)
    if (-90 < angle < 0):
        if vh < 0:
            bearing = 90.0 - abs(angle)
        if vh > 0:
            bearing = 270.0 - abs(angle)
    if angle == 0: bearing = 0.0
            
    return bearing

def generate_figures(plot_quantity, file_loc, report, reportname, surface,
                     leg_label, f_list, gauges, locations, elev, gauge_index,
                     production_dirs, time_min, time_max, time_unit,
                     title_on, label_id, generate_fig, verbose):
    """ Generate figures based on required quantities and gauges for
    each sww file
    """
    from math import sqrt, atan, degrees
    from Numeric import ones, allclose, zeros, Float, ravel
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
            texfile = texdir+reportname+'%s' %(label_id2)
            texfile2 = reportname+'%s' %(label_id2)
            texfilename = texfile + '.tex'
            if verbose: print '\n Latex output printed to %s \n' %texfilename
            fid = open(texfilename, 'w')
        else:
            texfile = texdir+reportname 
            texfile2 = reportname
            texfilename = texfile + '.tex' 
            if verbose: print '\n Latex output printed to %s \n' %texfilename
            fid = open(texfilename, 'w')
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
    model_time = zeros((n0,m,p), Float) 
    stages = zeros((n0,m,p), Float)
    elevations = zeros((n0,m,p), Float) 
    momenta = zeros((n0,m,p), Float)
    xmom = zeros((n0,m,p), Float)
    ymom = zeros((n0,m,p), Float)
    speed = zeros((n0,m,p), Float)
    bearings = zeros((n0,m,p), Float)
    depths = zeros((n0,m,p), Float)
    eastings = zeros((n0,m,p), Float)
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
    model_time_plot3d = zeros((n0,m), Float)
    stages_plot3d = zeros((n0,m), Float)
    eastings_plot3d = zeros((n0,m),Float)
    if time_unit is 'mins': scale = 60.0
    if time_unit is 'hours': scale = 3600.0
    ##### loop over each swwfile #####
    for j, f in enumerate(f_list):
        starttime = f.starttime
        if verbose: print 'swwfile %d of %d' %(j, len(f_list))
        comparefile = file_loc[j]+sep+'gauges_maxmins'+'.csv'
        fid_compare = open(comparefile, 'w')
        file0 = file_loc[j]+'gauges_t0.csv'
        fid_0 = open(file0, 'w')
        ##### loop over each gauge #####
        for k in gauge_index:
            g = gauges[k]
            if verbose: print 'Gauge %d of %d' %(k, len(gauges))
            min_stage = 10
            max_stage = 0
            max_momentum = max_xmomentum = max_ymomentum = 0
            min_momentum = min_xmomentum = min_ymomentum = 100
            max_speed = 0
            min_speed = 0            
            max_depth = 0            
            gaugeloc = str(locations[k])
            thisfile = file_loc[j]+sep+'gauges_time_series'+'_'\
                       +gaugeloc+'.csv'
            fid_out = open(thisfile, 'w')
            s = 'Time, Stage, Momentum, Speed, Elevation, xmom, ymom \n'
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
                    #bearing = calc_bearing(uh, vh)                    
                    model_time[i,k,j] = (t + starttime)/scale #t/60.0
                    stages[i,k,j] = w
                    elevations[i,k,j] = z 
                    xmom[i,k,j] = uh 
                    ymom[i,k,j] = vh 
                    momenta[i,k,j] = m 
                    speed[i,k,j] = vel 
                    #bearings[i,k,j] = bearing 
                    depths[i,k,j] = depth
                    thisgauge = gauges[k]
                    eastings[i,k,j] = thisgauge[0]
                    s = '%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f\n' %(t, w, m, vel, z, uh, vh)
                    fid_out.write(s)
                    if t == 0:
                        s = '%.2f, %.2f, %.2f\n' %(g[0], g[1], w)
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
                    
                    
            s = '%.2f, %.2f, %.2f, %.2f, %s\n' %(max_stage, min_stage, z, thisgauge[0], leg_label[j])
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
            
        if surface is  True:
            print 'Printing surface figure'
            for i in range(2):
                fig = p1.figure(10)
                ax = p3.Axes3D(fig)
                if len(gauges) > 80:
                    ax.plot_surface(model_time[:,:,j],eastings[:,:,j],stages[:,:,j])
                else:
                    #ax.plot_wireframe(model_time[:,:,j],eastings[:,:,j],stages[:,:,j])
                    ax.plot3D(ravel(eastings[:,:,j]),ravel(model_time[:,:,j]),ravel(stages[:,:,j]))
                ax.set_xlabel('time')
                ax.set_ylabel('x')
                ax.set_zlabel('stage')
                fig.add_axes(ax)
                p1.show()
                surfacefig = 'solution_surface%s' %leg_label[j]
                p1.savefig(surfacefig)
                p1.close()
            
    #### finished generating quantities for all swwfiles #####

    # x profile for given time
    if surface is True:
        figure(11)
        plot(eastings[tindex,:,j],stages[tindex,:,j])
        xlabel('x')
        ylabel('stage')
        profilefig = 'solution_xprofile' 
        savefig('profilefig')

    elev_output = []
    if generate_fig is True:
        depth_axis = axis([starttime/scale, time_max/scale, -0.1, max(max_depths)*1.1])
        stage_axis = axis([starttime/scale, time_max/scale, min(min_stages), max(max_stages)*1.1])
        vel_axis = axis([starttime/scale, time_max/scale, min(min_speeds), max(max_speeds)*1.1])
        mom_axis = axis([starttime/scale, time_max/scale, min(min_momentums), max(max_momentums)*1.1])
        xmom_axis = axis([starttime/scale, time_max/scale, min(min_xmomentums), max(max_xmomentums)*1.1])
        ymom_axis = axis([starttime/scale, time_max/scale, min(min_ymomentums), max(max_ymomentums)*1.1])
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
                s = '\\begin{figure}[ht] \n \\centering \n \\begin{tabular}{cc} \n'
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
                    s = '\\begin{figure}[hbt] \n \\centering \n \\begin{tabular}{cc} \n'
                    fid.write(s)
                    
                for which_quantity in plot_quantity:
                    count += 1
                    where1 += 1
                    figure(count, frameon = False)
                    if which_quantity == 'depth':
                        plot(model_time[0:n[j]-1,k,j], depths[0:n[j]-1,k,j], '-', c = cstr[j])
                        units = 'm'
                        axis(depth_axis)
                    if which_quantity == 'stage':
                        if elevations[0,k,j] <= 0:
                            plot(model_time[0:n[j]-1,k,j], stages[0:n[j]-1,k,j], '-', c = cstr[j])
                            axis(stage_axis)
                        else:
                            plot(model_time[0:n[j]-1,k,j], depths[0:n[j]-1,k,j], '-', c = cstr[j])
                            #axis(depth_axis)                 
                        units = 'm'
                    if which_quantity == 'momentum':
                        plot(model_time[0:n[j]-1,k,j], momenta[0:n[j]-1,k,j], '-', c = cstr[j])
                        axis(mom_axis)
                        units = 'm^2 / sec'
                    if which_quantity == 'xmomentum':
                        plot(model_time[0:n[j]-1,k,j], xmom[0:n[j]-1,k,j], '-', c = cstr[j])
                        axis(xmom_axis)
                        units = 'm^2 / sec'
                    if which_quantity == 'ymomentum':
                        plot(model_time[0:n[j]-1,k,j], ymom[0:n[j]-1,k,j], '-', c = cstr[j])
                        axis(ymom_axis)
                        units = 'm^2 / sec'
                    if which_quantity == 'speed':
                        plot(model_time[0:n[j]-1,k,j], speed[0:n[j]-1,k,j], '-', c = cstr[j])
                        axis(vel_axis)
                        units = 'm / sec'
                    if which_quantity == 'bearing':
                        due_east = 90.0*ones(shape(model_time[0:n[j]-1,k,j],Float))
                        due_west = 270.0*ones(shape(model_time[0:n[j]-1,k,j],Float))
                        plot(model_time[0:n[j]-1,k,j], bearings, '-', 
                             model_time[0:n[j]-1,k,j], due_west, '-.', 
                             model_time[0:n[j]-1,k,j], due_east, '-.')
                        units = 'degrees from North'
                        ax = axis([time_min, time_max, 0.0, 360.0])
                        legend(('Bearing','West','East'))

                    if time_unit is 'mins': xlabel('time (mins)')
                    if time_unit is 'hours': xlabel('time (hours)')
                    #if which_quantity == 'stage' and elevations[0:n[j]-1,k,j] > 0:
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
                        graphname_latex = '%sgauge%s%s' %(latex_file_loc, gaugeloc2, which_quantity) # storing files in production directory    
                        graphname_report_input = '%sgauge%s%s' %('..'+altsep+'report_figures'+altsep, gaugeloc2, which_quantity) # giving location in latex output file
                        graphname_report.append(graphname_report_input)
                        
                        savefig(graphname_latex) # save figures in production directory for report generation

                    if report == True:

                        figdir = getcwd()+sep+'report_figures'+sep
                        if access(figdir,F_OK) == 0 :
                            mkdir (figdir)
                        latex_file_loc = figdir.replace(sep,altsep)    

                        if len(label_id) == 1: 
                            graphname_latex = '%sgauge%s%s%s' %(latex_file_loc, gaugeloc2, which_quantity, label_id2) # storing files in production directory  
                            graphname_report = '%sgauge%s%s%s' %('..'+altsep+'report_figures'+altsep, gaugeloc2, which_quantity, label_id2) # giving location in latex output file
                            s = '\includegraphics[width=0.49\linewidth, height=50mm]{%s%s}' %(graphname_report, '.png')
                            fid.write(s)
                            if where1 % 2 == 0:
                                s = '\\\\ \n'
                                where1 = 0
                            else:
                                s = '& \n'
                            fid.write(s)
                            savefig(graphname_latex)
                    
                    if title_on == True:
                        title('%s scenario: %s at %s gauge' %(label_id, which_quantity, gaugeloc2))
                        #title('Gauge %s (MOST elevation %.2f, ANUGA elevation %.2f)' %(gaugeloc2, elevations[10,k,0], elevations[10,k,1] ))

                    savefig(graphname) # save figures with sww file

                if report == True and len(label_id) == 1:
                    for i in range(nn-1):
                        if nn > 2:
                            if plot_quantity[i] == 'stage' and elevations[0,k,j] > 0:
                                word_quantity += 'depth' + ', '
                            else:
                                word_quantity += plot_quantity[i] + ', '
                        else:
                            if plot_quantity[i] == 'stage' and elevations[0,k,j] > 0:
                                word_quantity += 'depth' + ', '
                            else:
                                word_quantity += plot_quantity[i]
                        
                    if plot_quantity[nn-1] == 'stage' and elevations[0,k,j] > 0:
                        word_quantity += ' and ' + 'depth'
                    else:
                        word_quantity += ' and ' + plot_quantity[nn-1]
                    caption = 'Time series for %s at %s location (elevation %.2fm)' %(word_quantity, locations[k], elev[k]) #gaugeloc.replace('_',' '))
                    if elev[k] == 0.0:
                        caption = 'Time series for %s at %s location (elevation %.2fm)' %(word_quantity, locations[k], elevations[0,k,j])
                        east = gauges[0]
                        north = gauges[1]
                        elev_output.append([locations[k],east,north,elevations[0,k,j]])
                    label = '%sgauge%s' %(label_id2, gaugeloc2)
                    s = '\end{tabular} \n \\caption{%s} \n \label{fig:%s} \n \end{figure} \n \n' %(caption, label)
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
                        s = '\includegraphics[width=0.49\linewidth, height=50mm]{%s%s}' %(graphname_report[index], '.png')
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
                caption = 'Time series for %s at %s location (elevation %.2fm)' %(word_quantity, locations[k], elev[k])
                if elev[k] == 0.0:
                        caption = 'Time series for %s at %s location (elevation %.2fm)' %(word_quantity, locations[k], elevations[0,k,j])
                        thisgauge = gauges[k]
                        east = thisgauge[0]
                        north = thisgauge[1]
                        elev_output.append([locations[k],east,north,elevations[0,k,j]])
                        
                s = '\end{tabular} \n \\caption{%s} \n \label{fig:%s} \n \end{figure} \n \n' %(caption, label)
                fid.write(s)
                if float((k+1)/div - pp) == 0.:
                    fid.write('\\clearpage \n')
                    pp += 1
                
                #### finished generating figures ###

            close('all')
        
    return texfile2, elev_output

# FIXME (DSG): Add unit test, make general, not just 2 files,
# but any number of files.
def copy_code_files(dir_name, filename1, filename2):
    """Temporary Interface to new location"""

    from anuga.shallow_water.data_manager import \
                    copy_code_files as dm_copy_code_files
    print 'copy_code_files has moved from util.py.  ',
    print 'Please use "from anuga.shallow_water.data_manager \
                                        import copy_code_files"'
    
    return dm_copy_code_files(dir_name, filename1, filename2)


def add_directories(root_directory, directories):
    """
    Add the first directory in directories to root_directory.
    Then add the second
    directory to the first directory and so on.

    Return the path of the final directory.

    This is handy for specifying and creating a directory 
    where data will go.
    """
    dir = root_directory
    for new_dir in directories:
        dir = os.path.join(dir, new_dir)
        if not access(dir,F_OK):
            mkdir (dir)
    return dir

def get_data_from_file(filename,separator_value = ','):
    """Temporary Interface to new location"""
    from anuga.shallow_water.data_manager import \
                        get_data_from_file as dm_get_data_from_file
    print 'get_data_from_file has moved from util.py'
    print 'Please use "from anuga.shallow_water.data_manager \
                                     import get_data_from_file"'
    
    return dm_get_data_from_file(filename,separator_value = ',')

def store_parameters(verbose=False,**kwargs):
    """Temporary Interface to new location"""
    
    from anuga.shallow_water.data_manager \
                    import store_parameters as dm_store_parameters
    print 'store_parameters has moved from util.py.'
    print 'Please use "from anuga.shallow_water.data_manager \
                                     import store_parameters"'
    
    return dm_store_parameters(verbose=False,**kwargs)

def remove_lone_verts(verts, triangles, number_of_full_nodes=None):
    """
    Removes vertices that are not associated with any triangles.

    verts is a list/array of points
    triangles is a list of 3 element tuples.  
    Each tuple represents a triangle.

    number_of_full_nodes relate to parallelism when a mesh has an
        extra layer of ghost points.
        
    """
    verts = ensure_numeric(verts)
    triangles = ensure_numeric(triangles)
    
    N = len(verts)
    
    # initialise the array to easily find the index of the first loner
    loners=arange(2*N, N, -1) # if N=3 [6,5,4]
    for t in triangles:
        for vert in t:
            loners[vert]= vert # all non-loners will have loners[i]=i 
    #print loners

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
        # verts=take(verts,X)  to Remove the loners from verts
        # but I think it would use more memory
        new_i = 0
        for i in range(lone_start, N):
            if loners[i] >= N:
                # Loner!
                pass
            else:
                loners[i] = new_i
                verts[new_i] = verts[i]
                new_i += 1
        verts = verts[0:new_i]

        # Modify the triangles
        #print "loners", loners
        #print "triangles before", triangles
        triangles = choose(triangles,loners)
        #print "triangles after", triangles
    return verts, triangles

 
def get_centroid_values(x, triangles):
    """Compute centroid values from vertex values
    
    x: Values at vertices of triangular mesh
    triangles: Nx3 integer array pointing to vertex information
    for each of the N triangels. Elements of triangles are
    indices into x
    """

        
    xc = zeros(triangles.shape[0], Float) # Space for centroid info
    
    for k in range(triangles.shape[0]):
        # Indices of vertices
        i0 = triangles[k][0]
        i1 = triangles[k][1]
        i2 = triangles[k][2]        
        
        xc[k] = (x[i0] + x[i1] + x[i2])/3


    return xc

def make_plots_from_csv_file(directories_dic={},
                            output_dir='',
                            base_name=None,
                            plot_numbers='0',
                            quantities=['Stage'],
                            extra_plot_name='',
                            assess_all_csv_files=True,                            
                            create_latex=False,
                            verbose=False):
                                
    """Read in csv files that have the right header information and
    plot time series such as Stage, Speed, etc. Will also plot several
    time series on one plot. Filenames must follow this convention,
    <base_name><plot_number>.csv eg gauge_timeseries3.csv

    Each file represents a location and within each file there are
    time, quantity columns.
    
    For example:    
    if "directories_dic" defines 4 directories and in each directories
    there is a csv files corresponding to the right "plot_numbers", 
    this will create a plot with 4 lines one for each directory AND 
    one plot for each "quantities".  
    
    Inputs:
        directories_dic: dictionary of directory with values (plot 
                         legend name for directory), (start time of 
                         the time series) and the (value to add to 
                         stage if needed). For example
                        {dir1:['Anuga_ons',5000, 0],
                         dir2:['b_emoth',5000,1.5],
                         dir3:['b_ons',5000,1.5]}
                         
        output_dir: directory for the plot outputs
        
        base_name: common name to all the csv files to be read
        
        plot_numbers: a String list of numbers to plot. For example 
                      [0-4,10,15-17] will read and attempt to plot
                       the follow 0,1,2,3,4,10,15,16,17
        quantities: Currently limited to "Stage", "Speed", and 
                     "Momentum", should be changed to incorporate 
                    any quantity read from CSV file....
                    
        extra_plot_name: A string that is appended to the end of the 
                         output filename.
                    
        assess_all_csv_files: if true it will read ALL csv file with
                             "base_name", regardless of 'plot_numbers'
                              and determine a uniform set of axes for 
                              Stage, Speed and Momentum. IF FALSE it 
                              will only read the csv file within the
                             'plot_numbers'
                             
        create_latex: NOT IMPLEMENTED YET!! sorry Jane....
        
        OUTPUTS: None, it saves the plots to 
              <output_dir><base_name><plot_number><extra_plot_name>.png
    
    """
    import pylab# import plot, xlabel, ylabel, savefig, \
                #ion, hold, axis, close, figure, legend
    from os import sep
    from anuga.shallow_water.data_manager import \
                               get_all_files_with_extension, csv2dict
    #find all the files that meet the specs
    #FIXME remove all references to word that such as Stage
    #(with a Capital letter) could use the header info from
    #the dictionary that is returned from the csv2dict this could 
    #be very good and very flexable.... but this works for now!

    #FIXME plot all available data from the csv file, currently will
    #only plot Stage,Speed and Momentum and only if the header is 
    #named like wise

#    quantities_dic={'time':'time (hour)'}
#    if 'time' in quantities:
#        quantities_dic['time']='time (hour)'
#    if 'Time' in quantities:
#        quantities_dic['Time']='time (hour)'
#    if 'stage' in quantities:
#        quantities_dic['stage']='wave height (m)'
#    if 'Stage' in quantities:
#        quantities_dic['Stage']='wave height (m)'
#    if 'speed' in quantities:
#        quantities_dic['speed']='speed (m/s)'
#    if 'Speed' in quantities:
#        quantities_dic['Speed']='speed (m/s)'
#    if 'stage' or 'Stage' in quantities:
#        quantities_dic['stage']='speed (m/s)'
#    if 'stage' or 'Stage' in quantities:
#        quantities_dic['stage']='speed (m/s)'

    seconds_in_hour = 3600
    time_label = 'time (hour)'
    stage_label = 'wave height (m)'
    speed_label = 'speed (m/s)'
    momentum_label = 'momentum (m^2/sec)'
    
    if extra_plot_name != '':
        extra_plot_name='_'+extra_plot_name

    #finds all the files that fit the specs and return a list of them
    #so to help find a uniform max and min for the plots... 
    list_filenames=[]
    if verbose: print 'Determining files to access for axes ranges \n'
    for i,directory in enumerate(directories_dic.keys()):
        list_filenames.append(get_all_files_with_extension(directory,
                              base_name,'.csv'))
#    print 'list_filenames',list_filenames

    #use all the files to get the values for the plot axis
    max_st=max_sp=max_mom=min_st=min_sp=min_mom=max_t=min_t=0.
    max_start_time= 0.
    min_start_time = 100000 

    
    new_plot_numbers=[]
    #change plot_numbers to list, eg ['0-4','10'] 
    #to ['0','1','2','3','4','10']
    for i, num_string in enumerate(plot_numbers):
        if '-' in num_string: 
            start = int(num_string[:num_string.rfind('-')])
            end = int(num_string[num_string.rfind('-')+1:])+1
            for x in range(start, end):
                new_plot_numbers.append(str(x))
        else:
            new_plot_numbers.append(num_string)
#    print 'new_plot_numbers',new_plot_numbers
    
    if verbose: print 'Determining uniform axes \n' 
    #this entire loop is to determine the min and max range for the 
    #axes of the plot
    for i, directory in enumerate(directories_dic.keys()):
        
        if assess_all_csv_files==False:
            which_csv_to_assess = new_plot_numbers
        else:
            which_csv_to_assess = list_filenames[i]

        for j, filename in enumerate(which_csv_to_assess):
            if assess_all_csv_files==False:
#                dir_filename=directory+sep+base_name+filename
                dir_filename=join(directory,base_name+filename)
            else:
                dir_filename=join(directory,filename)
#            print'dir_filename',dir_filename
            attribute_dic, title_index_dic = csv2dict(dir_filename+
                                                       '.csv')

            directory_start_time = directories_dic[directory][1]
            directory_add_tide = directories_dic[directory][2]

            time = [float(x) for x in attribute_dic["Time"]]
            min_t, max_t = get_min_max_values(time,min_t,max_t)
            
            stage = [float(x) for x in attribute_dic["Stage"]]
            stage =array(stage)+directory_add_tide
            min_st, max_st = get_min_max_values(stage,min_st,max_st)
            
            speed = [float(x) for x in attribute_dic["Speed"]]
            min_sp, max_sp = get_min_max_values(speed,min_sp,max_sp)
            
            momentum = [float(x) for x in attribute_dic["Momentum"]]
            min_mom, max_mom = get_min_max_values(momentum,
                                                  min_mom,
                                                  max_mom)
                                                                      
#            print 'min_sp, max_sp',min_sp, max_sp, 
            # print directory_start_time
            if min_start_time > directory_start_time: 
                min_start_time = directory_start_time
            if max_start_time < directory_start_time: 
                max_start_time = directory_start_time
#            print 'start_time' , max_start_time, min_start_time
    
    stage_axis = (min_start_time/seconds_in_hour,
                 (max_t+max_start_time)/seconds_in_hour,
                  min_st, max_st)
    speed_axis = (min_start_time/seconds_in_hour,
                 (max_t+max_start_time)/seconds_in_hour,
                 min_sp, max_sp)
    momentum_axis = (min_start_time/seconds_in_hour,
                    (max_t+max_start_time)/seconds_in_hour,
                     min_mom, max_mom)
#    ion()
#    pylab.hold()
    
    
    cstr = ['b', 'r', 'g', 'c', 'm', 'y', 'k']
    
#    if verbose: print 'Now start to plot \n'
    
    i_max = len(directories_dic.keys())
    legend_list_dic =[]
    legend_list =[]
    for i, directory in enumerate(directories_dic.keys()): 
        if verbose: print'Plotting in %s', directory 
        for j, number in enumerate(new_plot_numbers):
            if verbose: print'Starting %s',base_name+number  
            directory_name = directories_dic[directory][0]
            directory_start_time = directories_dic[directory][1]
            directory_add_tide = directories_dic[directory][2]
            
            #create an if about the start time and tide hieght 
            #if they don't exist
            file=directory+sep+base_name+number
            #print 'i %s,j %s, number %s, file %s' %(i,j,number,file)
            attribute_dic, title_index_dic = csv2dict(file+'.csv')
            #get data from dict in to list
            t = [float(x) for x in attribute_dic["Time"]]

            #do maths to list by changing to array
            t=(array(t)+directory_start_time)/seconds_in_hour
            stage = [float(x) for x in attribute_dic["Stage"]]
            speed = [float(x) for x in attribute_dic["Speed"]]
            momentum = [float(x) for x in attribute_dic["Momentum"]]
                        
            #Can add tide so plots are all on the same tide height
            stage =array(stage)+directory_add_tide
           
            #finds the maximum elevation
            max_ele=-100000
            min_ele=100000
            elevation = [float(x) for x in attribute_dic["Elevation"]]
            
            min_ele, max_ele = get_min_max_values(elevation,
                                                  min_ele,
                                                  max_ele)
            if min_ele != max_ele:
                print "Note! Elevation changes in %s" %dir_filename
#            print 'min_ele, max_ele',min_ele, max_ele
            

            #populates the legend_list_dic with dir_name and the elevation
            if i==0:
                legend_list_dic.append({directory_name:max_ele})
            else:
                legend_list_dic[j][directory_name]=max_ele
           
            # creates a list for the legend after at "legend_dic" has been fully populated
            # only runs on the last iteration for all the gauges(csv) files
            # empties the list before creating it 
            if i==i_max-1:
                legend_list=[]
                for k, l in legend_list_dic[j].iteritems():
                    legend_list.append('%s (elevation = %sm)'%(k,l))
                    #print k,l, legend_list_dic[j]
         
            if "Stage" in quantities:
                pylab.figure(100+j)
                pylab.plot(t, stage, c = cstr[i], linewidth=1)
                pylab.xlabel(time_label)
                pylab.ylabel(stage_label)
                pylab.axis(stage_axis)
                pylab.legend(legend_list,loc='upper right')
                figname = '%sstage_%s%s.png' %(output_dir+sep,
                                               base_name+number,
                                               extra_plot_name)
                pylab.savefig(figname)
            if "Speed" in quantities:
                pylab.figure(200+j)
                pylab.plot(t, speed, c = cstr[i], linewidth=1)
                pylab.xlabel(time_label)
                pylab.ylabel(speed_label)
                pylab.axis(speed_axis)
                pylab.legend(legend_list,loc='upper right')
                figname = '%sspeed_%s%s.png' %(output_dir+sep,
                                               base_name+number,
                                               extra_plot_name)
                pylab.savefig(figname)
            if "Momentum" in quantities:
                pylab.figure(300+j)
                pylab.plot(t, momentum, c = cstr[i], linewidth=1)
                pylab.xlabel(time_label)
                pylab.ylabel(momentum_label)
                pylab.axis(momentum_axis)
                pylab.legend(legend_list,loc='upper right')
                figname = '%smomentum_%s%s.png' %(output_dir+sep,
                                                  base_name+number,
                                                  extra_plot_name)
                pylab.savefig(figname)
    if verbose: print 'Closing all plots'
    pylab.close('all')
    del pylab
    if verbose: print 'Finished closing plots'

def get_min_max_values(list=None,min1=100,max1=-100):
    """ Returns the min and max of the list it was provided.
    NOTE: default min and max may need to change depeending on
    your list
    """
    if list == None: print 'List must be provided'
#    min = max_list = 0
    if max(list) > max1: 
        max1 = max(list)
    if min(list) < min1: 
        min1 = min(list)
        
    return min1, max1


def get_runup_data_for_locations_from_file(gauge_filename,
                                           sww_filename,
                                           runup_filename,
                                           size=10,
                                           verbose=False):

    '''this will read a csv file with the header x,y. Then look in a square 'size'x2
    around this position for the 'max_inundaiton_height' in the 'sww_filename' and 
    report the findings in the 'runup_filename
    
    WARNING: NO TESTS!
    '''

    from anuga.shallow_water.data_manager import get_all_directories_with_name,\
                                                 get_maximum_inundation_data,\
                                                 csv2dict
                                                 
    file = open(runup_filename,"w")
    file.write("easting,northing,runup \n ")
    file.close()
    
    #read gauge csv file to dictionary
    attribute_dic, title_index_dic = csv2dict(gauge_filename)
    northing = [float(x) for x in attribute_dic["y"]]
    easting = [float(x) for x in attribute_dic["x"]]

    print 'Reading %s' %sww_filename

    runup_locations=[]
    for i, x in enumerate(northing):
#        print 'easting,northing',i,easting[i],northing[i]
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
        
        if verbose:print 'maximum inundation runup near %s is %s meters' %(x_y,run_up)
        
        #writes to file
        file = open(runup_filename,"a")
        temp = '%s,%s,%s \n' %(x_y[0], x_y[1], run_up)
        file.write(temp)
        file.close()

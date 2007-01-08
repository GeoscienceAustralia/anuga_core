"""This module contains various auxiliary function used by pyvolution.

It is also a clearing house for functions that may later earn a module
of their own.
"""

import anuga.utilities.polygon
import sys
import os

from os import remove, mkdir, access, F_OK, sep
from os.path import exists, basename
from warnings import warn
from shutil import copy

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

    interpolation_points - list of absolute UTM coordinates for points (N x 2)
    or geospatial object or points file name at which values are sought
    
    use_cache: True means that caching of intermediate result of
               Interpolation_function is attempted

    
    See Interpolation function for further documentation
    """


    # Build arguments and keyword arguments for use with caching or apply.
    args = (filename,)
    
    kwargs = {'domain': domain,
              'quantities': quantities,
              'interpolation_points': interpolation_points,
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

        f = cache(_file_function,
                  args, kwargs,
                  dependencies=[filename],
                  compression=False,                  
                  verbose=verbose)

    else:
        f = apply(_file_function,
                  args, kwargs)


    #FIXME (Ole): Pass cache arguments, such as compression, in some sort of
    #structure
        

    return f



def _file_function(filename,
                   domain=None,
                   quantities=None,
                   interpolation_points=None,
                   time_thinning=1,                                                
                   verbose=False):
    """Internal function
    
    See file_function for documentatiton
    """
    

    #FIXME (OLE): Should check origin of domain against that of file
    #In fact, this is where origin should be converted to that of domain
    #Also, check that file covers domain fully.

    #Take into account:
    #- domain's georef
    #- sww file's georef
    #- interpolation points as absolute UTM coordinates


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

    if quantities is None: 
        if domain is not None:
            quantities = domain.conserved_quantities



    if line[:3] == 'CDF':
        return get_netcdf_file_function(filename, domain, quantities,
                                        interpolation_points,
                                        time_thinning=time_thinning,
                                        verbose=verbose)
    else:
        raise 'Must be a NetCDF File'



def get_netcdf_file_function(filename,
                             domain=None,
                             quantity_names=None,
                             interpolation_points=None,
                             time_thinning=1,                             
                             verbose=False):
    """Read time history of spatial data from NetCDF sww file and
    return a callable object f(t,x,y)
    which will return interpolated values based on the input file.

    If domain is specified, model time (domain.starttime)
    will be checked and possibly modified
    
    All times are assumed to be in UTC

    See Interpolation function for further documetation
    
    """
    
    
    #FIXME: Check that model origin is the same as file's origin
    #(both in UTM coordinates)
    #If not - modify those from file to match domain
    #Take this code from e.g. dem2pts in data_manager.py
    #FIXME: Use geo_reference to read and write xllcorner...
        

    #FIXME: Maybe move caching out to this level rather than at the
    #Interpolation_function level (below)

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
        



    if domain is not None:
        #Update domain.startime if it is *earlier* than starttime
        if starttime > domain.starttime:
            msg = 'WARNING: Start time as specified in domain (%f)'\
                  %domain.starttime
            msg += ' is earlier than the starttime of file %s (%f).'\
                     %(filename, starttime)
            msg += ' Modifying domain starttime accordingly.'
            
            if verbose: print msg
            domain.starttime = starttime #Modifying model time
            if verbose: print 'Domain starttime is now set to %f'\
               %domain.starttime


        #If domain.startime is *later* than starttime,
        #move time back - relative to domain's time
        if domain.starttime > starttime:
            time = time - domain.starttime + starttime

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


    return Interpolation_function(time,
                                  quantities,
                                  quantity_names,
                                  vertex_coordinates,
                                  triangles,
                                  interpolation_points,
                                  time_thinning=time_thinning,
                                  verbose=verbose)





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

def start_screen_catcher(dir_name, myid=0, numprocs=1):
    """Used to store screen output and errors to file, if run on multiple 
    processes eachprocessor will have its own output and error file.
    """

    dir_name = dir_name
    if access(dir_name,F_OK) == 0:
        print 'Make directory %s' %dir_name
        mkdir (dir_name,0777)
    screen_output_name = dir_name + "screen_output_%d_%d.txt" %(myid,numprocs)
    screen_error_name = dir_name + "screen_error_%d_%d.txt" %(myid,numprocs)
    
    #used to catch screen output to file
    sys.stdout = Screen_Catcher(screen_output_name)
    sys.stderr = Screen_Catcher(screen_error_name)

class Screen_Catcher:
    """this simply catches the screen output and stores it to file defined by
    start_screen_catcher (above)
    """
    
    def __init__(self, filename):
        self.filename = filename
 
        if exists(self.filename)is True:
            remove(self.filename)
            print'Old existing file "%s" has been deleted' %(self.filename)

    def write(self, stuff):
        fid = open(self.filename, 'a')
        fid.write(stuff)

def get_version_info():
    """gets the version number of the SVN
    NOTE: doesn't work on all systems eg GA's cyclone (64 bit linux cluster)
    """
    
    import os, sys

    # Create dummy info 
    info = 'Revision: Version info could not be obtained.'
    info += 'A command line version of svn and access to the '
    info += 'repository is necessary and the output must '
    info += 'contain a line starting with "Revision:"'

    try:
        fid = os.popen('svn info')
    except:
        msg = 'svn is not recognised'
        warn(msg, UserWarning)
    else:    
        lines = fid.readlines()
        fid.close()
        for line in lines:
            if line.startswith('Revision:'):
                info = line
                break
       
    return info
    
def sww2timeseries(swwfiles,
                   gauge_filename,
                   production_dirs,
                   report = None,
                   reportname = None,
                   plot_quantity = None,
                   surface = None,
                   time_min = None,
                   time_max = None,
                   title_on = None,
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
                      
    """

    
    k = _sww2timeseries(swwfiles,
                        gauge_filename,
                        production_dirs,
                        report,
                        reportname,
                        plot_quantity,
                        surface,
                        time_min,
                        time_max,
                        title_on,
                        verbose)

    return k

def _sww2timeseries(swwfiles,
                    gauge_filename,
                    production_dirs,
                    report = None,
                    reportname = None,
                    plot_quantity = None,
                    surface = None,
                    time_min = None,
                    time_max = None,
                    title_on = None,
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
        plot_quantity = ['depth', 'speed', 'bearing']
    else:
        assert type(plot_quantity) == list,\
               'plot_quantity must be a list'
        check_list(plot_quantity)

    if surface is None:
        surface = False
        
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

        f = file_function(swwfile,
                          quantities = sww_quantity,
                          interpolation_points = gauges,
                          verbose = True,
                          use_cache = True)

        # determine which gauges are contained in sww file
        count = 0
        gauge_index = []
        print 'swwfile', swwfile
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
        leg_label.append(production_dirs[swwfiles[swwfile]])
        
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
                                                production_dirs, time_min, time_max, title_on, label_id, verbose)
    else:
        texfile = ''
        elev_output = []

    return texfile, elev_output
                         
#Fixme - Use geospatial to read this file - it's an xya file
#Need to include other information into this filename, so xya + Name - required for report
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
    line1 = lines[0]
    line11 = line1.split(',')
    east_index = len(line11)+1
    north_index = len(line11)+1
    name_index = len(line11)+1
    elev_index = len(line11)+1
    for i in range(len(line11)):
        if line11[i].strip('\n').strip(' ').lower() == 'easting': east_index = i
        if line11[i].strip('\n').strip(' ').lower() == 'northing': north_index = i
        if line11[i].strip('\n').strip(' ').lower() == 'name': name_index = i
        if line11[i].strip('\n').strip(' ').lower() == 'elevation': elev_index = i

    for line in lines[1:]:
        fields = line.split(',')
        if east_index < len(line11) and north_index < len(line11):
            gauges.append([float(fields[east_index]), float(fields[north_index])])
        else:
            msg = 'WARNING: %s does not contain location information' %(filename)
            raise Exception, msg
        if elev_index < len(line11): elev.append(float(fields[elev_index]))
        if name_index < len(line11):
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
                     production_dirs, time_min, time_max, title_on, label_id,
                     verbose):
    """ Generate figures based on required quantities and gauges for each sww file
    """
    from math import sqrt, atan, degrees
    from Numeric import ones, allclose, zeros, Float, ravel
    from os import sep, altsep, getcwd, mkdir, access, F_OK, environ
    from pylab import ion, hold, plot, axis, figure, legend, savefig, \
         xlabel, ylabel, title, close, subplot

    import pylab as p1
    if surface is True:
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
    max_momentums = []
    max_speeds = []
    max_depths = []
    model_time_plot3d = zeros((n0,m), Float)
    stages_plot3d = zeros((n0,m), Float)
    eastings_plot3d = zeros((n0,m),Float)
    ##### loop over each swwfile #####
    for j, f in enumerate(f_list):
        if verbose: print 'swwfile %d of %d' %(j, len(f_list))
        comparefile = file_loc[j]+sep+'gauges_maxmins'+'.csv'
        fid_compare = open(comparefile, 'w')
        ##### loop over each gauge #####
        for k in gauge_index:
            g = gauges[k]
            if verbose: print 'Gauge %d of %d' %(k, len(gauges))
            min_stage = 10
            max_stage = 0
            max_momentum = 0
            max_speed = 0
            max_depth = 0
            gaugeloc = locations[k]
            thisfile = file_loc[j]+sep+'gauges_time_series'+'_'+gaugeloc+'.csv'
            fid_out = open(thisfile, 'w')
            s = 'Time, Stage, Momentum, Speed, Elevation \n'
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
                        vel = m / (depth + 1.e-30) 
                    bearing = calc_bearing(uh, vh)
                    model_time[i,k,j] = t/60.0
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
                    s = '%.2f, %.2f, %.2f, %.2f, %.2f\n' %(t, w, m, vel, z)
                    fid_out.write(s)
                    if t/60.0 <= 13920: tindex = i
                    if w > max_stage: max_stage = w
                    if w < min_stage: min_stage = w
                    if m > max_momentum: max_momentum = m
                    if vel > max_speed: max_speed = vel
                    if z > 0 and depth > max_depth: max_depth = depth
                    
                    
            s = '%.2f, %.2f, %.2f, %.2f, %s\n' %(max_stage, min_stage, z, thisgauge[0], leg_label[j])
            fid_compare.write(s)
            max_stages.append(max_stage)
            min_stages.append(min_stage)
            max_momentums.append(max_momentum)
            max_speeds.append(max_speed)
            max_depths.append(max_depth)
            #### finished generating quantities for each swwfile #####
        
        model_time_plot3d[:,:] = model_time[:,:,j]
        stages_plot3d[:,:] = stages[:,:,j]
        eastings_plot3d[:,] = eastings[:,:,j]
            
        if surface ==  True:
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
    if surface == True:
        figure(11)
        plot(eastings[tindex,:,j],stages[tindex,:,j])
        xlabel('x')
        ylabel('stage')
        profilefig = 'solution_xprofile' 
        savefig('profilefig')

    depth_axis = axis([time_min/60.0, time_max/60.0, -0.1, max(max_depths)*1.1])
    stage_axis = axis([time_min/60.0, time_max/60.0, min(min_stages), max(max_stages)*1.1])
    vel_axis = axis([time_min/60.0, time_max/60.0, min(max_speeds), max(max_speeds)*1.1])
    mom_axis = axis([time_min/60.0, time_max/60.0, min(max_momentums), max(max_momentums)*1.1])  
    
    cstr = ['g', 'r', 'b', 'c', 'm', 'y', 'k']
    nn = len(plot_quantity)
    no_cols = 2
    elev_output = []
    if len(label_id) > 1: graphname_report = []
    pp = 1
    div = 11.
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
                    if elevations[0,k,j] < 0:
                        plot(model_time[0:n[j]-1,k,j], stages[0:n[j]-1,k,j], '-', c = cstr[j])
                        axis(stage_axis)
                    else:
                        plot(model_time[0:n[j]-1,k,j], depths[0:n[j]-1,k,j], '-', c = cstr[j])
                        axis(depth_axis)                 
                    units = 'm'
                if which_quantity == 'momentum':
                    plot(model_time[0:n[j]-1,k,j], momenta[0:n[j]-1,k,j], '-', c = cstr[j])
                    axis(mom_axis)
                    units = 'm^2 / sec'
                if which_quantity == 'xmomentum':
                    plot(model_time[0:n[j]-1,k,j], xmom[0:n[j]-1,k,j], '-', c = cstr[j])
                    axis(mom_axis)
                    units = 'm^2 / sec'
                if which_quantity == 'ymomentum':
                    plot(model_time[0:n[j]-1,k,j], ymom[0:n[j]-1,k,j], '-', c = cstr[j])
                    axis(mom_axis)
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

                xlabel('time (mins)')
                if which_quantity == 'stage' and elevations[0:n[j]-1,k,j] > 0:
                    ylabel('%s (%s)' %('depth', units))
                else:
                    ylabel('%s (%s)' %(which_quantity, units))
                if len(label_id) > 1: legend((leg_label),loc='upper right')

                gaugeloc1 = gaugeloc.replace(' ','')
                #gaugeloc2 = gaugeloc1.replace('_','')
                gaugeloc2 = locations[k].replace(' ','')
                graphname = '%sgauge%s_%s' %(file_loc[j], gaugeloc2, which_quantity)

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
    """Copies "filename1" and "filename2" to "dir_name". Very useful for 
    information management """

    if access(dir_name,F_OK) == 0:
        print 'Make directory %s' %dir_name
        mkdir (dir_name,0777)
    copy(filename1, dir_name + sep + basename(filename1))
    copy(filename2, dir_name + sep + basename(filename2))
#    copy (__file__, project.output_run_time_dir + basename(__file__))
    print 'Files %s and %s copied' %(filename1, filename2)


def add_directories(root_directory, directories):
    """
    Add the first directory in directories to root_directory.
    Then add the second
    direcotory to the first directory and so on.

    Return the path of the final directory.

    This is handy for specifying and creating a directory where data will go.
    """
    dir = root_directory
    for new_dir in directories:
        dir = os.path.join(dir, new_dir)
        if not access(dir,F_OK):
            mkdir (dir)
    return dir




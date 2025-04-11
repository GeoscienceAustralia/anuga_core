"""This module contains various auxiliary function used by pyvolution.

It is also a clearing house for functions that may later earn a module
of their own.
"""
import os

from os import remove, mkdir, access, F_OK, R_OK, W_OK, sep,getcwd
from os.path import exists, basename, split,join
from warnings import warn
from shutil import copy

from anuga.utilities.numerical_tools import ensure_numeric, angle, NAN
from anuga.file.csv_file import load_csv_as_dict

from math import sqrt, atan, degrees

import anuga.utilities.log as log

import numpy as num


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
    from .file_function import file_function as file_function_new
    return file_function_new(filename, domain, quantities, interpolation_points,
                      use_relative_time, time_thinning, time_limit, verbose, use_cache,
                      boundary_polygon, output_centroids)



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
                       r'\b|\b'.join(map(re.escape, list(dictionary.keys())))+ \
                       r'\b' )

    #For each match, lookup the corresponding value in the dictionary
    return regex.sub(lambda match: dictionary[match.group(0)], text)


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

    Due to a limitation with numeric, this can not evaluate 0/0
    In general, the user can fix by adding 1e-30 to the numerator.
    SciPy core can handle this situation.
    """

    import types
    import re

    assert isinstance(expression, str)
    assert type(dictionary) == dict

    #Convert dictionary values to textual representations suitable for eval
    D = {}
    for key in dictionary:
        D[key] = 'dictionary["%s"]' % key

    #Perform substitution of variables    
    expression = multiple_replace(expression, D)

    #Evaluate and return
    try:
        return eval(expression)
    except NameError as e:
        msg = 'Expression "%s" could not be evaluated: %s' % (expression, e)
        raise NameError(msg)
    except ValueError as e:
        msg = 'Expression "%s" could not be evaluated: %s' % (expression, e)
        raise ValueError(msg)


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
                raise Exception('Illegal input to get_textual_float: %s' % str(value))
        else:
            return format % float(value)

def get_gauges_from_file(filename):
    return gauge_get_from_file(filename)


def check_list(quantity):
    """ Check that input quantities in quantity list are possible
    """
    import sys

    all_quantity = ['stage', 'depth', 'momentum', 'xmomentum',
                    'ymomentum', 'speed', 'bearing', 'elevation']

    # convert all quanitiy names to lowercase
    for i,j in enumerate(quantity):
        quantity[i] = quantity[i].lower()

    # check that all names in 'quantity' appear in 'all_quantity'
    p = list(set(quantity).difference(set(all_quantity)))
    if len(p) != 0:
        msg = 'Quantities %s do not exist - please try again' % p
        raise Exception(msg)


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

    # if indeterminate, just return 
    if uh==0 and vh==0:
        return NAN
    
    return degrees(angle([uh, vh], [0, -1]))   


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


def store_parameters(verbose=False,**kwargs):
    """Temporary Interface to new location"""
    
    from anuga.shallow_water.data_manager \
                    import store_parameters as dm_store_parameters
    log.critical('store_parameters has moved from util.py.')
    log.critical('Please use "from anuga.shallow_water.data_manager '
                 'import store_parameters"')
    
    return dm_store_parameters(verbose=False,**kwargs)


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
            try:
                loners[vert]= vert # all non-loners will have loners[i]=i
            except IndexError:
                msg = 'IndexError: t = '+str(t)+' vert = '+str(vert)+' N = '+str(N)
                raise Exception(msg)

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
        # verts=num.take(verts,X,axis=0)  to Remove the loners from verts
        # but I think it would use more memory
        new_i = lone_start    # point at first loner - 'shuffle down' target
        for i in range(lone_start, N):
            if loners[i] >= N:    # [i] is a loner, leave alone
                pass
            else:        # a non-loner, move down
                loners[i] = new_i
                verts[new_i] = verts[i]
                new_i += 1
        verts = verts[0:new_i]

        # Modify the triangles
        triangles = num.choose(triangles,loners)
    return verts, triangles


def get_centroid_values(x, triangles):
    """Compute centroid values from vertex values
    
    x: Values at vertices of triangular mesh
    triangles: Nx3 integer array pointing to vertex information
    for each of the N triangels. Elements of triangles are
    indices into x
    """
        
    xc = num.zeros(triangles.shape[0], float) # Space for centroid info
    
    for k in range(triangles.shape[0]):
        # Indices of vertices
        i0 = triangles[k][0]
        i1 = triangles[k][1]
        i2 = triangles[k][2]        
        
        xc[k] = (x[i0] + x[i1] + x[i2])/3

    return xc


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
    raise Exception(msg)

    return csv2timeseries_graphs(directories_dic,
                                 output_dir,
                                 base_name,
                                 plot_numbers,
                                 quantities,
                                 extra_plot_name,
                                 assess_all_csv_files)


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

#     try: 
#         import pylab
#     except ImportError:
#         msg='csv2timeseries_graphs needs pylab to be installed correctly'
#         raise Exception(msg)
#             #ANUGA don't need pylab to work so the system doesn't 
#             #rely on pylab being installed 
#         return
    
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as pylab
    except:
        #print "Couldn't import module from matplotlib, probably you need to update matplotlib"
        return

    from os import sep
    from anuga.utilities.file_utils import get_all_files_with_extension

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
    if verbose: log.critical('Determining files to access for axes ranges.')
    
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
    
    if verbose: log.critical('Determining uniform axes')

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
            attribute_dic, title_index_dic = load_csv_as_dict(dir_filename + '.csv')
            directory_start_time = directories_dic[directory][1]
            directory_add_tide = directories_dic[directory][2]

            if verbose: log.critical('reading: %s.csv' % dir_filename)

            #add time to get values
            for k, quantity in enumerate(quantities):
                quantity_value[quantity] = [float(x) for
                                                x in attribute_dic[quantity]]

                #add tide to stage if provided
                if quantity == 'stage':
                    quantity_value[quantity] = num.array(quantity_value[quantity],
                                                          float) + directory_add_tide

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
            log.critical('axis for quantity %s are x:(%s to %s)%s '
                         'and y:(%s to %s)%s' 
                         % (quantity, quantities_axis[quantity][0],
                            quantities_axis[quantity][1],
                            quantities_label['time'],
                            quantities_axis[quantity][2],
                            quantities_axis[quantity][3],
                            quantities_label[quantity]))

    cstr = ['b', 'r', 'g', 'c', 'm', 'y', 'k']

    if verbose: log.critical('Now start to plot')
    
    i_max = len(list(directories_dic.keys()))
    legend_list_dic = {}
    legend_list = []
    for i, directory in enumerate(directories_dic.keys()):
        if verbose: log.critical('Plotting in %s %s'
                                 % (directory, new_plot_numbers))

        # FIXME THIS SORT IS VERY IMPORTANT
        # Without it the assigned plot numbers may not work correctly
        # there must be a better way
        list_filenames[i].sort()
        for j, filename in enumerate(list_filenames[i]):
            if verbose: log.critical('Starting %s' % filename)

            directory_name = directories_dic[directory][0]
            directory_start_time = directories_dic[directory][1]
            directory_add_tide = directories_dic[directory][2]
            
            # create an if about the start time and tide height if don't exist
            attribute_dic, title_index_dic = load_csv_as_dict(directory + sep
                                                      + filename + '.csv')
            #get data from dict in to list
            #do maths to list by changing to array
            t = (num.array(directory_quantity_value[directory][filename]['time'])
                     + directory_start_time)/seconds_in_minutes

            #finds the maximum elevation, used only as a test
            # and as info in the graphs
            max_ele=-100000
            min_ele=100000
            elevation = [float(x) for x in attribute_dic["elevation"]]
            
            min_ele, max_ele = get_min_max_values(elevation)
            
            if min_ele != max_ele:
                log.critical("Note! Elevation changes in %s" % dir_filename)

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

                    if verbose: log.critical('saving figure here %s' % figname)

                    pylab.savefig(figname)
           
    if verbose: log.critical('Closing all plots')

    pylab.close('all')
    del pylab

    if verbose: log.critical('Finished closing plots')

def get_min_max_values(list=None):
    """ 
    Returns the min and max of the list it was provided.
    """

    if list is None: log.critical('List must be provided')
        
    return min(list), max(list)


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

    from anuga.shallow_water.data_manager import \
        get_maximum_inundation_data
                                                 
    file = open(runup_filename, "w")
    file.write("easting,northing,runup \n ")
    file.close()
    
    #read gauge csv file to dictionary
    attribute_dic, title_index_dic = load_csv_as_dict(gauge_filename)
    northing = [float(x) for x in attribute_dic["y"]]
    easting = [float(x) for x in attribute_dic["x"]]

    log.critical('Reading %s' % sww_filename)

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
            log.critical('maximum inundation runup near %s is %s meters'
                         % (x_y, run_up))
        
        #writes to file
        file = open(runup_filename, "a")
        temp = '%s,%s,%s \n' % (x_y[0], x_y[1], run_up)
        file.write(temp)
        file.close()

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
    from .gauge import sww2csv_gauges as sww2csv
    
    return sww2csv(sww_file, gauge_file, out_name, quantities, verbose, use_cache)
    
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
                   verbose=False,
                   output_centroids=False):
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
                        - structure which can be converted to a numeric array,
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


    from .gauge import sww2timeseries as sww2timeseries_new
    return sww2timeseries_new(swwfiles,
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
                   verbose,
                   output_centroids)                   
    
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
        

def square_root(s):
    return sqrt(s)



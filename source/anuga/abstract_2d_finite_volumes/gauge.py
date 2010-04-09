"""Gauge functions
   
   High-level functions for converting gauge and sww files into timeseries plots.


   Copyright 2010
   Ole Nielsen, Stephen Roberts, Duncan Gray, Christopher Zoppou, James Hudson
   Geoscience Australia
"""

import numpy as num

from anuga.geospatial_data.geospatial_data import ensure_absolute
from util import check_list, generate_figures
from file_function import file_function

import os

from os import remove, mkdir, access, F_OK, R_OK, W_OK, sep, getcwd
from os.path import exists, split, join
import anuga.utilities.log as log

from math import sqrt

##
# @brief Take gauge readings, given a gauge file and a sww mesh
#
#        Use this function to take a timeseries sample, given a list of gauge points
# @param  sww_file sww file to use as input
# @param gauge_file gauge file as input, containing list of gauge points to sample
# @param out_name output file prefix
# @param quantities which quantities in the sww file we want to export
# @param verbose show extra logging information for debug purposes
# @param use_cache cache requests if possible, for speed
# @param output_centroids Set to true to output the values at the centroid of the mesh triangle
def sww2csv_gauges(sww_file,
                   gauge_file,
                   out_name='gauge_',
                   quantities=['stage', 'depth', 'elevation',
                               'xmomentum', 'ymomentum'],
                   verbose=False,
                   use_cache=True,
                   output_centroids=False):
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
    from anuga.abstract_2d_finite_volumes.util import file_function	

    assert type(gauge_file) == type(''), 'Gauge filename must be a string'
    assert type(out_name) == type(''), 'Output filename prefix must be a string'
    
    try:
        point_reader = reader(file(gauge_file))
    except Exception, e:
        msg = 'File "%s" could not be opened: Error="%s"' % (gauge_file, e)
        raise msg

    if verbose: log.critical('Gauges obtained from: %s' % gauge_file)
    
    point_reader = reader(file(gauge_file))
    points = []
    point_name = []
    
    # read point info from file
    for i,row in enumerate(point_reader):
        # read header and determine the column numbers to read correctly.
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
    points_array = num.array(points,num.float)
        
    points_array = ensure_absolute(points_array)

    dir_name, base = os.path.split(sww_file)    

    #need to get current directory so when path and file
    #are "joined" below the directory is correct
    if dir_name == '':
        dir_name =getcwd()
        
    if access(sww_file,R_OK):
        if verbose: log.critical('File %s exists' % sww_file)
    else:
        msg = 'File "%s" could not be opened: no read permission' % sww_file
        raise msg

    sww_files = get_all_swwfiles(look_in_dir=dir_name,
                                 base_name=base,
                                 verbose=verbose)

    # fudge to get SWW files in 'correct' order, oldest on the left
    sww_files.sort()

    if verbose:
        log.critical('sww files=%s' % sww_files)
    
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
    
    if verbose: log.critical('Writing csv files')

    quake_offset_time = None

    for sww_file in sww_files:
        sww_file = join(dir_name, sww_file+'.sww')
        callable_sww = file_function(sww_file,
                                     quantities=core_quantities,
                                     interpolation_points=points_array,
                                     verbose=verbose,
                                     use_cache=use_cache,
                                     output_centroids = output_centroids)

        if quake_offset_time is None:
            quake_offset_time = callable_sww.starttime

        for time in callable_sww.get_time():
            for point_i, point in enumerate(points_array):
                #add domain starttime to relative time.
                quake_time = time + quake_offset_time
                points_list = [quake_time, quake_time/3600.]# fudge around SWW time bug
                point_quantities = callable_sww(time, point_i) # __call__ is overridden
                            
                for quantity in quantities:
                    if quantity == NAN:
                        log.critical('quantity does not exist in %s'
                                     % callable_sww.get_name)
                    else:
                        #core quantities that are exported from the interpolator     
                        if quantity == 'stage':
                            points_list.append(point_quantities[0])
                            
                        if quantity == 'elevation':
                            points_list.append(point_quantities[1])
                            
                        if quantity == 'xmomentum':
                            points_list.append(point_quantities[2])
                            
                        if quantity == 'ymomentum':
                            points_list.append(point_quantities[3])

                        #derived quantities that are calculated from the core ones
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
                                    vel = momentum / (point_quantities[0] 
                                                      - point_quantities[1])
                                else:
                                    momentum = 0
                                    vel = 0
                                
                            points_list.append(vel)
                            
                        if quantity == 'bearing':
                            points_list.append(calc_bearing(point_quantities[2],
                                                            point_quantities[3]))
                        if quantity == 'xcentroid':
                            points_list.append(callable_sww.centroids[point_i][0])

                        if quantity == 'ycentroid':
                            points_list.append(callable_sww.centroids[point_i][1])
							
                points_writer[point_i].writerow(points_list)

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

    msg = 'NOTE: A new function is available to create csv files from sww '
    msg += 'files called sww2csv_gauges in anuga.abstract_2d_finite_volumes.util'
    msg += ' PLUS another new function to create graphs from csv files called '
    msg += 'csv2timeseries_graphs in anuga.abstract_2d_finite_volumes.util'
    log.critical(msg)
    
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
                        verbose,
						output_centroids = output_centroids)
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
                    verbose = False,
					output_centroids = False):   
        
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
    
    if verbose: log.critical('Gauges obtained from: %s' % gauge_filename)

    gauges, locations, elev = gauge_get_from_file(gauge_filename)

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
            log.critical('swwfile = %s' % swwfile)

        # Extract parent dir name and use as label
        path, _ = os.path.split(swwfile)
        _, label = os.path.split(path)        
        
        leg_label.append(label)

        f = file_function(swwfile,
                          quantities = sww_quantity,
                          interpolation_points = gauges,
                          time_thinning = time_thinning,
                          verbose = verbose,
                          use_cache = use_cache,
						  output_centroids = output_centroids)

        # determine which gauges are contained in sww file
        count = 0
        gauge_index = []
        for k, g in enumerate(gauges):
            if f(0.0, point_id = k)[2] > 1.0e6:
                count += 1
                if count == 1: log.critical('Gauges not contained here:')
                log.critical(locations[k])
            else:
                gauge_index.append(k)

        if len(gauge_index) > 0:
            log.critical('Gauges contained here:')
        else:
            log.critical('No gauges contained here.')
        for i in range(len(gauge_index)):
             log.critical(locations[gauge_index[i]])
             
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
         log.critical('Inputs OK - going to generate figures')

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
def gauge_get_from_file(filename):
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

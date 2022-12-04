"""Gauge functions

   High-level functions for converting gauge and sww files into timeseries plots.


   Copyright 2010
   Ole Nielsen, Stephen Roberts, Duncan Gray, Christopher Zoppou, James Hudson
   Geoscience Australia
"""



from builtins import str
from six import string_types
from builtins import range
import numpy as num

from anuga.geospatial_data.geospatial_data import ensure_absolute
from .util import check_list, calc_bearing
from .file_function import file_function

import os

from os import remove, mkdir, access, F_OK, R_OK, W_OK, sep, getcwd
from os.path import exists, split, join
import anuga.utilities.log as log

from math import sqrt

def _quantities2csv(quantities, point_quantities, centroids, point_i):
    points_list = []

    for quantity in quantities:
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
                    vel = momentum/(point_quantities[0] - point_quantities[1])
                else:
                    momentum = 0
                    vel = 0

            points_list.append(vel)

        if quantity == 'bearing':
            points_list.append(calc_bearing(point_quantities[2],
                                            point_quantities[3]))
        if quantity == 'xcentroid':
            points_list.append(centroids[point_i][0])

        if quantity == 'ycentroid':
            points_list.append(centroids[point_i][1])

    return points_list


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
    from anuga.utilities.file_utils import get_all_swwfiles
    from anuga.abstract_2d_finite_volumes.util import file_function

    assert isinstance(gauge_file,string_types) or isinstance(gauge_file, str), 'Gauge filename must be a string or unicode'
    assert isinstance(out_name,string_types) or isinstance(out_name, str), 'Output filename prefix must be a string'

    try:
        gid = open(gauge_file)
        point_reader = reader(gid)
        gid.close()
    except Exception as e:
        msg = 'File "%s" could not be opened: Error="%s"' % (gauge_file, e)
        raise Exception(msg)

    if verbose: log.critical('Gauges obtained from: %s' % gauge_file)

    gid = open(gauge_file)
    point_reader = reader(gid)

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
            #points.append([float(row[easting]),float(row[northing])])
            points.append([float(row[easting]),float(row[northing])])
            point_name.append(row[name])

    gid.close()

    #convert to array for file_function
    points_array = num.array(points,float)

    points_array = ensure_absolute(points_array)

    #print 'points_array', points_array

    dir_name, base = os.path.split(sww_file)

    #need to get current directory so when path and file
    #are "joined" below the directory is correct
    if dir_name == '':
        dir_name =getcwd()

    if access(sww_file,R_OK):
        if verbose: log.critical('File %s exists' % sww_file)
    else:
        msg = 'File "%s" could not be opened: no read permission' % sww_file
        raise Exception(msg)

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

    if verbose: log.critical('Writing csv files')

    quake_offset_time = None

    is_opened = [False]*len(points_array)
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

        for point_i, point in enumerate(points_array):
            for time in callable_sww.get_time():
                # add domain starttime to relative time.
                quake_time = time + quake_offset_time
                point_quantities = callable_sww(time, point_i) # __call__ is overridden


                if point_quantities[0] != NAN:
                    if is_opened[point_i] == False:
                        points_handle = open(dir_name + sep + gauge_file
                                             + point_name[point_i] + '.csv', 'w', newline='')
                        points_writer = writer(points_handle)
                        points_writer.writerow(heading)
                        is_opened[point_i] = True
                    else:
                        points_handle = open(dir_name + sep + gauge_file
                                             + point_name[point_i] + '.csv', 'a', newline='')
                        points_writer = writer(points_handle)


                    points_list = [quake_time, quake_time/3600.] +  _quantities2csv(quantities, point_quantities, callable_sww.centroids, point_i)
                    points_writer.writerow(points_list)
                    points_handle.close()
                else:
                    if verbose:
                        msg = 'gauge' + point_name[point_i] + 'falls off the mesh in file ' + sww_file + '.'
                        log.warning(msg)

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
    except Exception as e:
        msg = 'File "%s" could not be opened: Error="%s"' % (gauge_filename, e)
        raise Exception(msg)

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

    for swwfile in list(swwfiles.keys()):
        try:
            fid = open(swwfile)
        except Exception as e:
            msg = 'File "%s" could not be opened: Error="%s"' % (swwfile, e)
            raise Exception(msg)

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
            raise Exception(msg)

    if time_max is None:
        time_max = themaxT # max(T)
    else:
        if time_max > themaxT: # max(T):
            msg = 'Maximum time entered not correct - please try again'
            raise Exception(msg)

    if verbose and len(gauge_index) > 0:
         log.critical('Inputs OK - going to generate figures')

    if len(gauge_index) != 0:
        texfile, elev_output = \
            _generate_figures(plot_quantity, file_loc, report, reportname,
                             surface, leg_label, f_list, gauges, locations,
                             elev, gauge_index, production_dirs, time_min,
                             time_max, time_unit, title_on, label_id,
                             generate_fig, verbose)
    else:
        texfile = ''
        elev_output = []

    return texfile, elev_output


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

    if isinstance(line11[0], string_types) is True:
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
            raise Exception(msg)

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
        gaugelocation = list(range(N))

    # Read in gauge data
    for line in lines:
        fields = line.split(',')

        gauges.append([float(fields[east_index]), float(fields[north_index])])

        if len(fields) > 2:
            elev.append(float(fields[elev_index]))
            loc = fields[name_index]
            gaugelocation.append(loc.strip(r'\n'))

    return gauges, gaugelocation, elev


def _generate_figures(plot_quantity, file_loc, report, reportname, surface,
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

            if verbose: log.critical('Latex output printed to %s' % texfilename)
        else:
            texfile = texdir+reportname
            texfile2 = reportname
            texfilename = texfile + '.tex'
            fid = open(texfilename, 'w')

            if verbose: log.critical('Latex output printed to %s' % texfilename)
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
    model_time = num.zeros((n0, m, p), float)
    stages = num.zeros((n0, m, p), float)
    elevations = num.zeros((n0, m, p), float)
    momenta = num.zeros((n0, m, p), float)
    xmom = num.zeros((n0, m, p), float)
    ymom = num.zeros((n0, m, p), float)
    speed = num.zeros((n0, m, p), float)
    bearings = num.zeros((n0, m, p), float)
    due_east = 90.0*num.ones((n0, 1), float)
    due_west = 270.0*num.ones((n0, 1), float)
    depths = num.zeros((n0, m, p), float)
    eastings = num.zeros((n0, m, p), float)
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
    model_time_plot3d = num.zeros((n0, m), float)
    stages_plot3d = num.zeros((n0, m), float)
    eastings_plot3d = num.zeros((n0, m),float)
    if time_unit == 'mins': scale = 60.0
    if time_unit == 'hours': scale = 3600.0

    ##### loop over each swwfile #####
    for j, f in enumerate(f_list):
        if verbose: log.critical('swwfile %d of %d' % (j, len(f_list)))

        starttime = f.starttime
        comparefile = file_loc[j] + sep + 'gauges_maxmins' + '.csv'
        fid_compare = open(comparefile, 'w')
        file0 = file_loc[j] + 'gauges_t0.csv'
        fid_0 = open(file0, 'w')

        ##### loop over each gauge #####
        for k in gauge_index:
            if verbose: log.critical('Gauge %d of %d' % (k, len(gauges)))

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
                        vel = m/ (depth + 1.e-6/depth)
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
            log.critical('Printing surface figure')
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
                s = r'\\begin{figure}[ht] \n' \
                    r'\\centering \n' \
                    r'\\begin{tabular}{cc} \n'
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
                    s = r'\\begin{figure}[hbt] \n' \
                        r'\\centering \n' \
                        r'\\begin{tabular}{cc} \n'
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

                    if time_unit == 'mins': xlabel('time (mins)')
                    if time_unit == 'hours': xlabel('time (hours)')
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
                            s = r'\includegraphics' \
                                r'[width=0.49\linewidth, height=50mm]{%s%s}' % \
                                (graphname_report, '.png')
                            fid.write(s)
                            if where1 % 2 == 0:
                                s = r'\\\\ \n'
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
                    s = r'\end{tabular} \n' \
                        r'\\caption{%s} \n' \
                        r'\label{fig:%s} \n' \
                        r'\end{figure} \n \n' % (caption, label)
                    fid.write(s)
                    cc += 1
                    if cc % 6 == 0: fid.write(r'\\clearpage \n')
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
                        s = r'\includegraphics' \
                            r'[width=0.49\linewidth, height=50mm]{%s%s}' % \
                            (graphname_report[index], '.png')
                        index += 1
                        fid.write(s)
                        if where1 % 2 == 0:
                            s = r'\\\\ \n'
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

                s = r'\end{tabular} \n' \
                    r'\\caption{%s} \n' \
                    r'\label{fig:%s} \n' \
                    r'\end{figure} \n \n' % (caption, label)
                fid.write(s)
                if float((k+1)//div - pp) == 0.:
                    fid.write(r'\\clearpage \n')
                    pp += 1
                #### finished generating figures ###

            close('all')

    return texfile2, elev_output

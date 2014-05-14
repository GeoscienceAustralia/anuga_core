"""
Replacement function that wraps the legacy function urs2sts()
and passes in the data from the EventSelection <event>.list file.
"""

import os
import os.path
from time import localtime, strftime, gmtime

import project

from anuga.file.netcdf import NetCDFFile
import numpy as num

from anuga.file_conversion.urs2sts import urs2sts


#-------------------------------------------------------------------------------
# Get gauges (timeseries of index points)
#-------------------------------------------------------------------------------
def get_sts_gauge_data(filename, verbose=False):
    print 'get_sts_gauge_data: filename=%s' % filename
    fid = NetCDFFile(filename+'.sts', 'r')      #Open existing file for read
    permutation = fid.variables['permutation'][:]
    x = fid.variables['x'][:] + fid.xllcorner   #x-coordinates of vertices
    y = fid.variables['y'][:] + fid.yllcorner   #y-coordinates of vertices
    points = num.transpose(num.asarray([x.tolist(), y.tolist()]))
    time = fid.variables['time'][:] + fid.starttime
    elevation = fid.variables['elevation'][:]
       
    basename = 'sts_gauge'
    quantity_names = ['stage', 'xmomentum', 'ymomentum']
    quantities = {}
    for i, name in enumerate(quantity_names):
        quantities[name] = fid.variables[name][:]

    #---------------------------------------------------------------------------
    # Get maximum wave height throughout timeseries at each index point
    #---------------------------------------------------------------------------
     
    maxname = 'max_sts_stage.csv'
    print 'get_sts_gauge_data: maxname=%s' % maxname
    fid_max = open(os.path.join(project.event_folder, maxname), 'w')
    fid_max.write('index, x, y, max_stage \n')    
    for j in range(len(x)):
        index = permutation[j]
        stage = quantities['stage'][:,j]
        xmomentum = quantities['xmomentum'][:,j]
        ymomentum = quantities['ymomentum'][:,j]

        fid_max.write('%d, %.6f, %.6f, %.6f\n' % (index, x[j], y[j], max(stage)))
      
    #---------------------------------------------------------------------------
    # Get minimum wave height throughout timeseries at each index point
    #---------------------------------------------------------------------------

    minname = 'min_sts_stage.csv'
    fid_min = open(os.path.join(project.event_folder, minname), 'w')
    fid_min.write('index, x, y, max_stage \n')    
    for j in range(len(x)):
        index = permutation[j]
        stage = quantities['stage'][:,j]
        xmomentum = quantities['xmomentum'][:,j]
        ymomentum = quantities['ymomentum'][:,j]

        fid_min.write('%d, %.6f, %.6f, %.6f\n' %(index, x[j], y[j], min(stage)))

        out_file = os.path.join(project.event_folder,
                                basename+'_'+str(index)+'.csv')
        fid_sts = open(out_file, 'w')
        fid_sts.write('time, stage, xmomentum, ymomentum \n')

        #-----------------------------------------------------------------------
        # End of the get gauges
        #-----------------------------------------------------------------------
        for k in range(len(time)-1):
            fid_sts.write('%.6f, %.6f, %.6f, %.6f\n'
                          % (time[k], stage[k], xmomentum[k], ymomentum[k]))

        fid_sts.close()      

    fid.close()

    return quantities,elevation,time


##
# @brief Build boundary STS files from one or more MUX files.
# @param event_file Name of mux meta-file or single mux stem.
# @param output_dir Directory to write STS data to.
# @note 'event_file' is produced by EventSelection.
def build_urs_boundary(event_file, output_dir):
    '''Build a boundary STS file from a set of MUX files.'''

    # if we are using an EventSelection multi-mux file
    if project.multi_mux:
        # get the mux+weight data from the meta-file (in <boundaries>)
        mux_event_file = os.path.join(project.event_folder, event_file)
        print 'using multi-mux file', mux_event_file
        try:
            fd = open(mux_event_file, 'r')
            mux_data = fd.readlines()
            fd.close()
        except IOError, e:
            msg = 'File %s cannot be read: %s' % (mux_event_file, str(e))
            raise Exception, msg
        except:
            raise

        # first line of file is # filenames+weight in rest of file
        num_lines = int(mux_data[0].strip())
        mux_data = mux_data[1:]
        print 'number of sources %d' % num_lines

        # quick sanity check on input mux meta-file
        if num_lines != len(mux_data):
            msg = ('Bad file %s: %d data lines, but line 1 count is %d'
                   % (event_file, len(mux_data), num_lines))
            raise Exception, msg

        # Create filename and weights lists.
        # Must chop GRD filename just after '*.grd'.
        mux_filenames = []
        for line in mux_data:
            muxname = line.strip().split()[0]
            split_index = muxname.index('.grd')
            muxname = muxname[:split_index+len('.grd')]
            muxname = os.path.join(project.mux_data_folder, muxname)
            mux_filenames.append(muxname)
        
        mux_weights = [float(line.strip().split()[1]) for line in mux_data]

        # Call legacy function to create STS file.
        print 'creating sts file'
        urs2sts(mux_filenames,
                basename_out=output_dir,
                ordering_filename=project.urs_order,
                weights=mux_weights,
                zone=project.zone,
                mean_stage=project.tide,
                verbose=True)
    else:                           # a single mux stem file, assume 1.0 weight
        mux_file = os.path.join(project.event_folder, event_file)
        mux_filenames = [mux_file]
        print 'using single-mux file', mux_file

        weight_factor = 1.0
        mux_weights = weight_factor*num.ones(len(mux_filenames), num.float)
            
        order_filename = project.urs_order

        print 'reading', order_filename
        # Create ordered sts file
        urs2sts(mux_filenames,
                basename_out=output_dir,
                ordering_filename=order_filename,
                weights=mux_weights,
                mean_stage=project.tide,
                verbose=True)

    # report on progress so far
    sts_file = os.path.join(project.event_folder, project.scenario_name)
    quantities, elevation, time = get_sts_gauge_data(sts_file, verbose=False)
    print len(elevation), len(quantities['stage'][0,:])

    


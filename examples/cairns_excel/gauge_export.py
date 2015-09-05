'''

Given an sww file and a set of gauge coordinates, export flow data from the
centroid coordinates nearest the gauge

Gareth Davies, Geoscience Australia 2014+

'''

###############################################################
#
# INPUT PARAMETERS
#
###############################################################

sww_file = 'MODEL_OUTPUTS/RUN_20150625_110925_cairns_excel/cairns_excel.sww'
gauge_sites = {'Offshore1': [666779.0, 8253357.],
               'Offshore2': [684597., 7983715.],
               'Offshore3': [395059., 8095373.]}

output_directory = 'gauge_data'

################################################################
#
# END INPUT PARAMETERS
#
################################################################

import sys

import anuga
from anuga.file.netcdf import NetCDFFile
from anuga import plot_utils as util
import numpy
import os

# Optionally take commandline arguments
l = len(sys.argv)
if l > 1:
    sww_files = [sys.argv[i] for i in range(1, l)]
else:
    sww_files = [sww_file]


for sww_file in sww_files:
    # Output directory inside sww directory
    outdir = os.path.split(sww_file)[0] + '/' + output_directory

    try:
        os.mkdir(outdir)
    except:
        pass

    # Separate labels and coordinates
    gauge_labels = gauge_sites.keys()
    gauge_coordinates = [gauge_sites[site] for site in gauge_labels]
    assert len(gauge_coordinates) == len(
        gauge_labels), 'Must have one label for each gauge coordinate'

    # Get first-timestep sww file data to help later extraction
    p = util.get_output(sww_file, timeSlices=[0])
    pc = util.get_centroids(p, timeSlices=[0])

    # Get centroid information
    xc = pc.x + pc.xllcorner
    yc = pc.y + pc.yllcorner
    elevation_c = pc.elev

    # Make file connection
    fid = NetCDFFile(sww_file)
    #gauge_indices = [((xc - p[0]) ** 2 + (yc - p[1]) ** 2).argmin()
    #                 for p in gauge_coordinates]
    gauge_indices = []
    for gc in gauge_coordinates:
        new_gauge_point = [gc[0] - p.xllcorner, gc[1] - p.yllcorner]
        gi = util.get_triangle_containing_point(p, new_gauge_point)
        gauge_indices.append(gi) 

    time = fid.variables['time'][:]
    nts = len(time)

    # Export flow data for all gauges
    for i, gi in enumerate(gauge_indices):
        #stage = util._getCentVar(fid, 'stage_c', time_indices=range(nts), space_indices=gi)
        try:
            stage = fid.variables['stage_c'][:, gi]
        except:
            vols = fid.variables['volumes'][gi, :]
            stage = (fid.variables['stage'][:, vols[0]] +
                     fid.variables['stage'][:, vols[1]] +
                     fid.variables['stage'][:, vols[2]]) / 3.

        if len(elevation_c.shape) > 1:
            elev_c = elevation_c[0, gi]
        else:
            elev_c = elevation_c[gi]
        depth = stage - elev_c
        #xmom = util._getCentVar(fid, 'xmomentum_c', time_indices=range(nts), space_indices=gi)
        try:
            xmom = fid.variables['xmomentum_c'][:, gi]
            ymom = fid.variables['ymomentum_c'][:, gi]
        except:
            vols = fid.variables['volumes'][gi, :]
            xmom = (fid.variables['xmomentum'][:, vols[0]] +
                    fid.variables['xmomentum'][:, vols[1]] +
                    fid.variables['xmomentum'][:, vols[2]]) / 3.
            ymom = (fid.variables['ymomentum'][:, vols[0]] +
                    fid.variables['ymomentum'][:, vols[1]] +
                    fid.variables['ymomentum'][:, vols[2]]) / 3.
        #ymom = util._getCentVar(fid, 'ymomentum_c', time_indices=range(nts), space_indices=gi)
        xvel = xmom / (depth + 1.0e-12)
        yvel = ymom / (depth + 1.0e-12)

        # Export the data
        export_array = numpy.vstack([time, stage.flatten(), depth.flatten(),
                                     xvel.flatten(), yvel.flatten(), xmom.flatten(), ymom.flatten()])
        export_array = export_array.transpose()

        filename = outdir + '/' + 'Gauge_' + gauge_labels[i] + '_' + \
            str(gauge_coordinates[i][0]) + '_' + \
            str(gauge_coordinates[i][1]) + '.csv'
        f = open(filename, 'w')

        header1 = '# The requested coordinate was ' + \
            str(gauge_coordinates[i][0]) + ' ' + \
            str(gauge_coordinates[i][1]) + '\n'
        header2 = '# The nearest centroid coordinate (used to get the flow values) was ' + str(
            xc[gi]) + ' ' + str(yc[gi]) + '\n'
        distance_offset = ((xc[gi] - gauge_coordinates[i][0]) ** 2 +
                           (yc[gi] - gauge_coordinates[i][1]) ** 2) ** 0.5
        header3 = '# The distance between the requested and nearest coordinate is ' + \
            str(distance_offset) + ' m' + '\n'
        f.write(header1)
        f.write(header2)
        f.write(header3)
        f.write('time, stage, depth, xvel, yvel, xvel_d, yvel_d\n')
        numpy.savetxt(f, export_array, delimiter=",")

"""
Export csv file with flow values at centroids for a chosen set of time-slices
in the sww file

Run this from ipython with 

>   %run -i points_export.py


Gareth Davies, Geoscience Australia 2014+
"""
#################################################################
# INPUT DATA

sww_file = 'MODEL_OUTPUTS/RUN_20150625_110925_cairns_excel/cairns_excel.sww'
timesteps = range(0, 5, 2) # List of indices of timesteps that you want to export
outdir = 'point_data'

################################################################


import os
outdir = os.path.split(sww_file)[0] + '/' + outdir

try:
    os.mkdir(outdir)
except:
    pass

import matplotlib
matplotlib.use('Agg')
from anuga import plot_utils as util
import numpy


for timestep in timesteps:

    p = util.get_centroids(sww_file, timeSlices=timestep)

    export_array = numpy.vstack(
        [p.x + p.xllcorner, p.y + p.yllcorner,
         p.stage[0, :], p.height[0, :],
         p.xvel[0, :], p.yvel[0, :],
         p.xmom[0, :], p.ymom[0, :],
         p.elev]).transpose()

    file_name = outdir + '/' + 'point_outputs_' + \
        str(timestep) + '_time_' + str(round(p.time[0])) + '.csv'
    f = open(file_name, 'w')
    f.write('x, y, stage, depth, xvel, yvel, xvel_d, yvel_d, elev\n')
    numpy.savetxt(f, export_array, delimiter=",")

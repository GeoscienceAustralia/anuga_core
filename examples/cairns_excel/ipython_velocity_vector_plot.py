"""
Make a velocity vector plot.

This can be viewed interactively with zooming/panning if run from ipython,
or it can be used to make a 'batch' of figures


Gareth Davies, Geoscience Australia, 2014+
"""


import matplotlib
from matplotlib import pyplot
import anuga
from anuga import plot_utils as util
import numpy

#############################################################
#
# User Inputs here
#
#############################################################

sww_file = 'MODEL_OUTPUTS/RUN_20150625_110925_cairns_excel/cairns_excel.sww'

# timesteps = range(500,1101,10) # 500, 510, 520, ... 1100
timesteps = 4

# Colors
colorbar = 'Blues'  # 'spectral' #'afmhot_r'#'summer' # 'Blues' # 'winter' #
edgecolors = 'k'
max_depth_on_colorbar = 1500.0
arrow_color = 'red'  # 'green'
arrow_scale = (0.8) ** (-1)  # 1/arrow_enlargement
variable = 'flux' #'velocity' # flux
plot_threshold = 0.1  # Smaller values of 'variable' are not plotted

# Plot bounding box
x_range = (580000, 620000)
y_range = (8160000, 8200000)

# Dimensions of figure:
figure_dim = (10, 8)  # (horizontal, vertical)
dpi = 300  # dots per inch

# Directory for figures
figdir = 'FIG'

##############################################################
#
# End user inputs
#
##############################################################

import os
figdir = os.path.split(sww_file)[0] + '/' + figdir

# Make directory
try:
    os.mkdir(figdir)
except:
    pass

# Ensure timesteps is a list
try:
    tmp = list(timesteps)
except:
    tmp = [timesteps]
timesteps = tmp

for timestep in timesteps:

    # Get vertex information
    p = util.get_output(sww_file, timeSlices=[timestep])
    # Get centroid information
    pc = util.get_centroids(p)

    # Compute maximum height at centroids
    hmax = numpy.minimum(pc.height.max(axis=0), max_depth_on_colorbar)

    # Choose colormap
    cm = matplotlib.cm.get_cmap(colorbar)

    # Make initial scatter plot
    if len(timesteps) == 1:
        pyplot.ion()

    pyplot.figure(figsize=figure_dim)
    pyplot.scatter(pc.x + p.xllcorner, pc.y + p.yllcorner,
                   c=hmax, cmap=cm, s=0.1, edgecolors='none')
    util.plot_triangles(
        p, adjustLowerLeft=True, values=hmax, values_cmap=cm, edgecolors=edgecolors)
    pyplot.gca().set_aspect('equal')
    pyplot.colorbar()

    if variable == 'velocity':
        nonzero_vel = (
            (pc.xvel[0, :] ** 2 + pc.yvel[0, :] ** 2) > plot_threshold ** 2).nonzero()[0]

        # Add arrows
        pyplot.quiver(pc.x[nonzero_vel] + p.xllcorner, pc.y[nonzero_vel] + p.yllcorner,
                      pc.xvel[0, nonzero_vel], pc.yvel[0, nonzero_vel],
                      scale=arrow_scale, scale_units='xy', color=arrow_color)
    elif variable == 'flux':
        nonzero_mom = (
            (pc.xmom[0, :] ** 2 + pc.ymom[0, :] ** 2) > plot_threshold ** 2).nonzero()[0]

        # Add arrows
        pyplot.quiver(pc.x[nonzero_mom] + p.xllcorner, pc.y[nonzero_mom] + p.yllcorner,
                      pc.xmom[0, nonzero_mom], pc.ymom[0, nonzero_mom],
                      scale=arrow_scale, scale_units='xy', color=arrow_color)

    pyplot.xlim(x_range)
    pyplot.ylim(y_range)

    figname = figdir + '/' + 'Fig_' + \
        str(timestep) + '_time_' + str(round(p.time[0])) + '.jpg'
    pyplot.savefig(figname, dpi=dpi)

    if (len(timesteps) > 1):
        pyplot.close()

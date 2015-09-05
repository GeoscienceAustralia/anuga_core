"""
Run this from ipython with 

> ipython
> %run -i raster_export.py

or in pure python with 

>   python raster_export.py


Gareth Davies, Geoscience Australia 2014+
"""
####################################################################
#
# INPUT DATA
#

# sww filename relative to the model_from_excel directory
sww_file = 'MODEL_OUTPUTS/RUN_20150625_111757_cairns_excel/cairns_excel.sww'

# xls filename in the model_from_excel directory
xls_config_file = 'cairns_excel.xls'

# Which timesteps to export?
# timesteps = [0, 13, 108, 950] # Multiple time steps at once
#
# More example choices below
# timesteps = range(500, 901, 10) # This gives 500, 510, 520, ...890, 900
timesteps = 'collected_max'  # Flow maxima
# timesteps = 12 # The 12th timestep

# Set the number of neighbours for raster interpolation
#   -- Values > 1 produce some smoothing, which can improve appearance
#   -- Larger k requires more memory
#   -- Memory can be reduced using a larger cellsize
k_nearest_neighbours = 3

# Optionally set a raster pixel size different to the one in ANUGA_setup_basic.xls
#   -- if cellsize=None then the value from xls_config_file is used
#cellsize = None
cellsize = 1000

#
# END INPUT DATA
#
##############################################################################

import os
from setup import raster_outputs as ro
from setup import prepare_data
project = prepare_data.PrepareData(xls_config_file,
                                   make_directories=False)

# Ensure timesteps is a list
if timesteps is not 'collected_max':
    try:
        tmp = list(timesteps)
    except:
        tmp = [timesteps]
    timesteps = tmp

if cellsize is None:
    cellsize = project.output_tif_cellsize

# Make the tifs
ro.make_me_some_tifs(
    sww_file,
    project.bounding_polygon,
    project.proj4string,
    my_time_step=timesteps,
    tif_output_subdir='/TIFS/',
    cell_size=cellsize,
    k_nearest_neighbours=k_nearest_neighbours)

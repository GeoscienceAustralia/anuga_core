"""Benchmark for sww2dem

Creates and exports a dem from an sww file.
"""

from anuga.shallow_water.data_manager import sww2dem
from create_test_sww import create_test_sww
import os.path

sww_name = 'test.sww'

# use existing file
if not os.path.isfile(sww_name):
	create_test_sww(sww_name)

# do export to DEM

sww2dem('test',
        basename_out='sww2dem_out',
        quantity='stage',
        cellsize=1,      
        easting_min=0,
        easting_max=100,
        northing_min=0,
        northing_max=100,        
        reduction=max, 
        verbose=True,
        format='asc')


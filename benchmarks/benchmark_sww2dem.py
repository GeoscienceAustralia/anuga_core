"""Benchmark for sww2dem

Creates and exports a dem from an sww file.
"""

from anuga.file_conversion.sww2dem import sww2dem
from create_test_sww import create_test_sww
import os.path
import cProfile
import time


sww_name = 'test.sww'

def sww2dem_test():
	# do export to DEM

	sww2dem('test',
			basename_out='sww2dem_out',
			quantity='stage',
			cellsize=0.25,      
			easting_min=0,
			easting_max=100,
			northing_min=0,
			northing_max=100,        
			reduction=max, 
			verbose=True,
			format='asc')
  
		
# use existing file
if not os.path.isfile(sww_name):
	create_test_sww(sww_name)


start_time = time.time()	
#cProfile.run('sww2dem_test()')
sww2dem_test()
stop_time = time.time()
print ('sww2dem took %.1fs\n\n\n' % (stop_time - start_time))

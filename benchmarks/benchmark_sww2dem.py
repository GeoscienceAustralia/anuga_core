"""Benchmark for sww2dem

Creates and exports a dem from an sww file.
"""

import anuga
import os.path
import cProfile
import time

from create_test_sww import create_test_sww

sww_name = 'test.sww'

def sww2dem_test():
	# do export to DEM

	anuga.sww2dem(sww_name,
			name_out='sww2dem_out.asc',
			quantity='xmomentum',
			cellsize=0.25,      
			easting_min=0,
			easting_max=100,
			northing_min=0,
			northing_max=100,        
			reduction=max, 
			verbose=True)
  
		
# use existing file
if not os.path.isfile(sww_name):
	create_test_sww(sww_name)


start_time = time.time()	
#cProfile.run('sww2dem_test()')
sww2dem_test()
stop_time = time.time()
print ('sww2dem took %.1fs\n\n\n' % (stop_time - start_time))

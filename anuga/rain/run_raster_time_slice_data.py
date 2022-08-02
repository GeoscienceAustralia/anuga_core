
import numpy as np


from anuga.rain.raster_time_slice_data import Raster_time_slice_data
import pylab as pl
from os.path import join


try:
    import ipdb as pdb
except:
    import pdb


start_time = 4700
final_time = 5500
rain = Raster_time_slice_data(start_time=start_time, final_time=final_time, verbose=True, debug=True)
rain.x = np.arange(41)/40.0*4000.0 +3000.0
rain.y = np.arange(61)/60.0*4000.0 + 100000.0
rain.extent = (rain.x.min(), rain.x.max(), rain.y.min(), rain.y.max())
rain.times = np.arange(11)*600+4000
rain.time_step = 600
nx = len(rain.x)
ny = len(rain.y)

rain.data_slices = [ i*np.ones((ny,nx),dtype="float")+4 for i, time in enumerate(rain.times)]

for i,data_slice in enumerate(rain.data_slices):
    # Rows vertical, Columns horizontal. from top
    data_slice[2*i:2*i+4,:] = 2
    data_slice[:,i:i+4] = 3

rain.data_accumulated = np.arange(nx*ny,dtype="float").reshape(ny,nx)
rain.data_accumulated[50:54,:] = 2
rain.data_accumulated[:,20:24] = 3
rain.data_max_in_period = np.max(rain.data_accumulated)
rain.radar_dir = None
l0 = [5100.0, 100000.0]
l1 = [5100.0, 103900.0]
l2 = [3400., 103400.]
l3 = [6000., 103400.]
locations = [l0,l1,l2,l3] 



print(rain.radar_dir)


pl.figure(1)

print(rain.get_extent())
rain.accumulate_data_stats()

print('time_step ',rain.time_step)     
    
import time

plot_vmax = np.max(rain.data_slices)
print('plot_vmax', plot_vmax)

for tid in range(len(rain.times)):
    rain.plot_data(tid, plot_vmax=plot_vmax, save=False, show=True)
    time.sleep(0.1)
    #ipdb.set_trace() 
    
pdb.set_trace()     

rain.plot_accumulated_data()

#pl.ioff()
pl.show()
    

rain.plot_time_hist_locations(locations) 

p_indices = rain.grid_indices_inside_polygon()



pdb.set_trace()
    
    
    
     


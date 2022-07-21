import anuga
import numpy as np 
import pylab as pl

from anuga.rain.calibrated_radar_rain import Calibrated_radar_rain

try:
    import ipdb as pdb
except:
    import pdb


scenario = 'act'
#scenario = 'grantham' # The files seem to be corrupted

HOME_DIR = "/mnt/c/Users/stoiv/Dropbox/"

if scenario == 'act':
    start_time = '20120301_0000'
    final_time = '20120302_1200'

    BASE_DIR =  anuga.join(HOME_DIR, "RAINFALL/RADAR/AUS_RADAR/Calibrated_Radar_Data/ACT_netcdf" )
    RADAR_DIR = anuga.join(BASE_DIR, 'RADAR_Rainfall/140/2012' )    
    Catchment_file     = anuga.join(BASE_DIR,'ACT_Bdy/ACT_Entire_Catchment_Poly_Simple_UTM_55.csv')
    State_boundary_file = anuga.join(BASE_DIR,'ACT_Bdy/ACT_State_Bdy_UTM_55.csv')
    Daily_plot_Vmax = 4
    rain = Calibrated_radar_rain(RADAR_DIR, start_time=start_time, final_time=final_time, verbose=True, debug=True)
    l0 = [670501.0, 6117750.0]
    l1 = [702251.0, 6157920.0]
    l2 = [748256., 5966120.]
    l3 = [600000., 5966120.] 
    locations = [l0,l1,l2,l3] 

    
        
elif scenario == 'grantham':
    start_time = '20110109_2300'
    final_time = '20110110_2300'
    
    #start_time = '20110108_2300'
    #final_time = '20110109_2300'
    
    #start_time = '20110105_2300'
    #final_time = '20110106_2300'
    
    BASE_DIR  = anuga.join(HOME_DIR, 'RAINFALL/RADAR/AUS_RADAR/Gauge_Blend_Grantham/20081211-20110223.gz_y_loc_Inverted_30min/Merged/66/2011/')
    RADAR_DIR = anuga.join(BASE_DIR, '01' ) 
    Daily_plot_Vmax = 50
    rain = Calibrated_radar_rain(RADAR_DIR, start_time=start_time, final_time=final_time, verbose=True, debug=True)
    l0 = [476160., 6889300.0]
    l1 = [476160., 6978420.0]
    l2 = [421517., 6974520.0]
    l3 = [479412.0, 6980370.0]
    locations = [l0,l1,l2,l3]

else:
    pass


print(rain.radar_dir)

try: 
    p2 = anuga.read_polygon(Catchment_file)
    print('Catchment File')
    print(Catchment_file)
    #print p2
except:
    p2 = None
    pass

try:
    p3 = anuga.read_polygon(State_boundary_file)
    print('State_boundary_file')
    print(State_boundary_file)
    #print p3
except:
    p3 = None
    pass
    
pl.figure(1)
#p2 = [[689801.0, 6091501.0], [700200.0, 6091501.0], [700200.0, 6104600.0], [689801.0, 6104600.0], [689801.0, 6091501.0]]

print(rain.get_extent())
rain.accumulate_data_stats()

print('time_step ',rain.time_step)     
    
import time
#pl.ion()


plot_vmax = np.max(rain.data_slices)
print('plot_vmax', plot_vmax)

for tid in range(len(rain.times)):
    rain.plot_data(tid, plot_vmax=plot_vmax, save=False, show=True, polygons=[p2,p3])
    time.sleep(0.1)
    #ipdb.set_trace() 
    
pdb.set_trace()     

rain.plot_accumulated_data(polygons=[p2,p3])

#pl.ioff()
pl.show()
    

rain.plot_time_hist_locations(locations) 

p_indices = rain.grid_indices_inside_polygon(polygon=p2)



pdb.set_trace()
    
    
    
     


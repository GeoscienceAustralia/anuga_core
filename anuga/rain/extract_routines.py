
"""
Basic helper routines


"""

import anuga
from anuga.fit_interpolate.interpolate2d import interpolate2d
import pylab as pl
import matplotlib.pyplot as plt
import numpy as np
import os
from os.path import join
import fnmatch
import gzip
from anuga.file.netcdf import NetCDFFile



Daily_plot_Vmax = 4

# ----------------------------------------------------------------------------------------



class Calibrated_radar_rain(object):
    
    def __init__(self, radar_dir = None, verbose=False, debug=False):
        
        self.radar_dir = radar_dir
        
        self.x = None
        self.y = None
        self.precips = []
        self.times = []
        self.precip_total = 0.0
        self.time_step = 0.0
        
        self.verbose = verbose
        self.debug = debug
        

        self.pattern = '*.nc'
        #self.pattern = '*.gz'


        # process the radar files
        if not radar_dir is None:
            self.read_data_files(radar_dir)
            

    def get_extent(self):
        
        return self.extent
 
    def read_data_files(self, radar_dir):
        """
        Given a radar_dir walk through all sub directories to find 
        radar files
        """
        
        pattern = self.pattern
        
        file_counter = 0
        first = True
        rain_max_in_period = 0.0
        
        precips = []
        times = []
        
        self.radar_dir = radar_dir
        
        for root, dirs, files in os.walk(radar_dir): 
            
            if self.debug or self.verbose:
                print 'Directory: ',dirs
                print 'Number of Files = ',len(fnmatch.filter(files, pattern))
            
            for filename in fnmatch.filter(files, pattern): 
                
                if self.debug :
                    print filename
                    
                key = filename[-16:-3]
                

                file_counter +=1
                if file_counter == 1:
                    self.radar_data_title = "RADAR_Data_"+filename[-20:-3]

                if pattern == '*.gz':
                    filename = gzip.open(filename, 'rb')
                    data = NetCDFFile(os.path.join(root,filename[:-3]), 'r') 
                else:
                    data = NetCDFFile(os.path.join(root, filename), 'r') 
                
                # RADAR NetCDF files have Dimensions, Attributes, Variables
                if self.debug:
                    print 'VARIABLES:'
                    print 'Reference LAT, LONG = ',data.reference_longitude, data.reference_latitude

                    
                # Check Time for key 
                year = int(key[0:4])
                month = int(key[4:6])
                day = int(key[6:8])
                hour = int(key[9:11])
                minute = int(key[11:13])
                
                if self.debug:
                    print year, month, day, hour, minute
                    print 'Convert to epoch'
                    
                import datetime
                valid_time = int((datetime.datetime(year,month,day,hour,minute) - datetime.datetime(1970,1,1)).total_seconds())
                
                
                file_start_time = data.variables['start_time'][0]
                file_valid_time = data.variables['valid_time'][0]
                
                if self.debug:
                    print 'VARIABLES:'
                    print 'Reference times ', file_start_time, file_valid_time, valid_time
                
                times.append(valid_time)                
        
                    
                # This handles format changes in the files from BOM !!!!
                possible_precip_names = ['precipitation',  'precip', 'rain_amount'] 
                
                # Go through each of the possible names
                for name in possible_precip_names:  # Check if name is a key in the variables dictionary
                    if name in data.variables:
                        precip_name = name
                        if self.debug:
                            print 'BOM Reference name tag in this file:'
                            print precip_name

                
                if first:
                    first = False
                    
                    self.base_filename = filename[:-20]
                    
                    precip = data.variables[precip_name][:]
                    precip_total = precip.copy() # Put into new Accumulating ARRRAY
                    
                    if self.debug: print ' Accumulate rainfall here....'
                    self.x = data.variables['x_loc'][:]
                    self.y = data.variables['y_loc'][:]
                    if self.y[0] < 0:
                        pass  # Check if y[0] = -ve if not reverse...  arr[::-1]
                    else:
                        self.y = self.y[::-1]  # Check if y[0] = -ve if not reverse...  arr[::-1]
                    
                    
                    self.reference_longitude = data.reference_longitude
                    self.reference_latitude =  data.reference_latitude
                    # Convert to UTM offsets
                    
                    self.zone, self.offset_x, self.offset_y = anuga.LLtoUTM(self.reference_latitude, self.reference_longitude)
                    
                    # Convert to UTM
                    self.x = self.x*1000 + self.offset_x
                    self.y = self.y*1000 + self.offset_y
                    
                    rain_max_in_period  = max(np.max(precip),rain_max_in_period)  
                        
    
                else:  # ---If NOT FIRST !!!
                    precip = data.variables[precip_name][:]

                    precip_total += precip
                    if self.debug: print ' Keep accumulating rainfall....'

                    rain_max_in_period  = max(np.max(precip),rain_max_in_period)

                    
                precips.append(precip)
        
        
        self.rain_max_in_period = rain_max_in_period
                  
        times = np.array(times)
        ids = np.argsort(times)
        

        self.times = times[ids]   
        self.precips = np.array([ precips[id] for id in ids ])
        self.time_step = self.times[1]-self.times[0]
        
        self.precip_total = precip_total 
        
        
        # Test sorting
        for i, id in enumerate(ids):
            np.allclose(times[id], self.times[i])
            np.allclose(precips[id], self.precips[i])
                
        self.extent = (self.x.min(), self.x.max(), self.y.min(), self.y.max())

    def extract_data_at_location(self, locations):
        """
        
        EXTRACT RAINFALL FROM GRID AT locations  
        
        """
        if self.verbose or self.debug: print 'Extract Rain Check data First'
        
        locations = np.array(locations)
        all_values = []
        x = self.x
        y = self.y
        precips = self.precips

        if self.debug: print locations
        
        for precip in precips:
            # and then do the interpolation
            values = interpolate2d(x,y,precip,locations)       
            all_values.append(values) 


        return np.array(all_values)




    def plot_grid(self, id, save=False, show=True, polygons=None):
        """
        Plot radar data at timestep id
        """

        # Aliases
        precips = self.precips
        times = self.times
        x = self.x
        y = self.y

        
        # get time from filenane
        from datetime import datetime
        date_time = datetime.utcfromtimestamp(times[id]).strftime("%Y/%m/%d %H:%M")
        
        if self.debug: print '--- Date/Time ', date_time
        
        plt.figure(1)
        plt.clf()

        
        if not polygons is None:
            for polygon in polygons:
                polygon = np.array(polygon)
                plt.plot( polygon[:,0], polygon[:,1],'--w')
            
        plot_title = 'RADAR RAINFALL '+date_time
        plt.suptitle(plot_title, size=20)
        s_title = 'Max Rainfall = %.1f mm / period' % (np.max(precips[id]))
        plt.title(s_title , size=12)
        plt.imshow(precips[id], origin='lower', interpolation='bicubic',
                   extent=self.extent,vmin=0, vmax=Daily_plot_Vmax)

        plt.colorbar()            
        
        if show: plt.draw()
        if save: plt.savefig(plot_title+'.jpg',format='jpg')
            
        return




    def plot_grid_accumulated(self, polygons=None):
        
        precip_total = self.precip_total

        time_step = self.time_step
        radar_dir = self.radar_dir
        x = self.x
        y = self.y
        
        
        if self.debug: 
            print precip_total
            print 'maximum rain =',np.max(precip_total)
            print 'mean rain =',np.mean(precip_total)
            
        Total_Rain_Vol = np.mean(precip_total)/1000.0*256.0*256.0 # Volume in Million m3 over 128km x 128km area
        Rain_Max_in_period = self.rain_max_in_period
        Peak_Intensity = Rain_Max_in_period/time_step
        extent = self.extent
        
        if self.verbose:
            print 'Total rainfall volume in Mill m3 =',Total_Rain_Vol
            print 'Peak rainfall in 1 time step = ', Rain_Max_in_period
            print 'Peak Intensity in 1 timestep =',Peak_Intensity
            print 'extent', extent
            print 'size', extent[1]-extent[0], extent[3]-extent[2]
            print time_step


        plt.figure(1)
        plt.clf()
        plt.suptitle('Accumulated Radar Rainfall', size=20)
        
        
        s_title = 'Max Int. = %.1f mm/hr, Ave. rain = %.1f mm, Tot rain Vol. = %.3f Mill. m3' \
                  % (Peak_Intensity,np.mean(precip_total),Total_Rain_Vol)
        plt.title(s_title , size=12)
        
        
        if not polygons is None:
            for polygon in polygons:
                polygon = np.array(polygon)
                plt.plot( polygon[:,0], polygon[:,1],'--w')
                
        plt.imshow(precip_total, origin='lower', interpolation='bicubic',
                   extent=(x.min(), x.max(), y.min(), y.max()))
        
        plt.colorbar()
        plt.draw()
        
        return
    
    def plot_time_hist_locations(self, locations):
        """
        Plot time histograms at selected locations i.e.
        gauge locations
        """
        

        t = (self.times-self.time_step- self.times[0])/60.0
        
        locations = np.array(locations)

        time_step = self.time_step
        
        all_values = self.extract_data_at_location(locations)

        for lid, _ in enumerate(locations):

            bar_values = [values[lid] for values in all_values]
            total_rain = sum(bar_values)
            Ave_rain = total_rain/(self.times[-1]-self.times[0])*3600
            max_Intensity = max(bar_values)/time_step*3600
            b_title = 'Tot rain = %.1f mm, Ave. rain = %.3f mm/hr, Max Int.= %.3f mm/hr' % (total_rain,Ave_rain,max_Intensity)
            #   Using list comprehension  b=[x[0] for x in a]
            #print len(t)
            #print len(bar_values)
            # Think about using...  zip(*lst)
            plt.bar(t,bar_values,width=self.time_step/60)
            plt.suptitle(' Data for Location %s:' % lid, fontsize=14, fontweight='bold')    
            plt.title(b_title)
    
            plt.xlabel('time (mins)')
            plt.ylabel('rainfall (mm)')
            plt.show()
            
        return

    
    
if __name__ == "__main__":
    

    
    BASE_DIR = "/home/steve/RAINFALL/RADAR/AUS_RADAR/Calibrated_Radar_Data/ACT_netcdf"
    RADAR_DIR          = join(BASE_DIR, 'RADAR_Rainfall/140/2012/03/01' )    
    Catchment_file     = join(BASE_DIR,'ACT_Bdy/ACT_Entire_Catchment_Poly_Simple_UTM_55.csv')
    State_boundary_file = join(BASE_DIR,'ACT_Bdy/ACT_State_Bdy_UTM_55.csv')

    
    #BASE_DIR = '/home/steve/RAINFALL/RADAR/AUS_RADAR/Gauge_Blend_Grantham/20081211-20110223.gz_y_loc_Inverted_30min/Merged/66/2011/01/'
    #RADAR_DIR          = join(BASE_DIR, '15' )  
    #RainGauge_LOC_file = join(BASE_DIR, 'ACT_Bdy/Rain_Gauge_Station_Location_subset.csv' )
    
    print RADAR_DIR
    
#     
    p2 = anuga.read_polygon(Catchment_file)
    print 'Catchment File'
    print Catchment_file
    #print p2
    
    p3 = anuga.read_polygon(State_boundary_file)
    print 'State_boundary_file'
    print State_boundary_file
    #print p3
     
    # Create object to read in and store radar data
    act_rain = Calibrated_radar_rain(RADAR_DIR, verbose=True)

    print 'time_step ',act_rain.time_step
    
    import time
    pl.ion()
    for rid in xrange(len(act_rain.times)):
        act_rain.plot_grid(rid, save=False, show=True, polygons=[p2])
        time.sleep(0.05)
    
    act_rain.plot_grid_accumulated(polygons=[p2])
    
    pl.ioff()
    pl.show()
    
    print act_rain.get_extent()

        
    l0 = [600000.0, 6000000.0]
    l1 = [800000.0, 6100000.0]
    locations = [l0,l1]
    act_rain.plot_time_hist_locations(locations) 
    
    import pdb
    #pdb.set_trace()    
    
    
    
     


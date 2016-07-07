
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



class Grid_rain_radar(object):
    
    def __init__(self,x=None,y=None,UTM_x=None,UTM_y=None,precip=None):
        
        self.precip = precip
        self.UTM_x = UTM_x
        self.UTM_y = UTM_y
        self.x = x
        self.y = y
        self.data = None
        self.offset_x = 0.0
        self.offset_y = 0.0
        self.filename = ""
        self.gauge_loc_points =[[0.0, 0.0]]
        self.gauge_loc_labels =  ["gauge0"]
        self.all_values = []
        self.precip_total = 0.0
        self.verbose = True
        self.debug = False
        
        self.all_precip = []
        self.all_times =[]
        
        self.pattern = '*.nc'
        
        # Plotting parameters
        self.plot_line = False
        
 
    def process_calibrated_radar_data(self, radar_dir):
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
        
        for root, dirs, files in os.walk(radar_dir): 
            
            if self.debug :
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
                    data = NetCDFFile(os.path.join(root,filename), 'r') 
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
        
        self.precip_total = precip_total 
        
        
        # Test sorting
        for i, id in enumerate(ids):
            np.allclose(times[id], self.times[i])
            np.allclose(precips[id], self.precips[i])
            

      


    def read_radar_UTM_offsets(self, Radar_UTM_LOC_file):
        """
        Read location of radar (in UTM coordinates)
        
        FIXME: Probably should have a geo reference!
        """
        
        fid = open(Radar_UTM_LOC_file)
        lines = fid.readlines()  # Read Entire Input File
        fid.close()  
               
        line=lines[0].strip('\n')
        fields = line.split(',')
        offset_x = float(fields[0])
        offset_y = float(fields[1])


        self.offset_x = offset_x
        self.offset_y = offset_y
        
        return 



    def convert_LLPrecip2UTM(self, x, y, precip):
        """
        OPTION TO CONVERT TO UTM Implimentation FIRST  
        # Only saves to a file not used for plots YET !!:
        """
        if self.verbose: print 'Converting to UTM...'
        outfilename = self.filename[0:-4]+'.xyz' # Output file for RADAR in UTM
        outfid = open(outfilename, 'w')
        for i in x:
            for j in y:
                if self.verbose: print x[i],y[j]
                UTM_x = x[i]*1000.0 + self.offset_x
                UTM_y = y[j]*1000.0 + self.offset_y  
                
                # FIXME: Precip not correct
                s = ' %.3f,%.3f,%.3f \n' %(UTM_x,UTM_y,precip[0][i])
                outfid.write(s)
                if self.verbose: print UTM_x,UTM_y,precip[0][i]  

        outfid.close()
        
        return                  


    def extract_radar_data_at_gauge_loc(self):
        """
        
        "OPTION TO EXTRACT RAINFALL FROM RADAR GRID AT GAUGE POINTS FIRST  
        
        """
        if self.verbose: print 'Extract Rain Check data First'
        
        data = self.data
        precip_name = self.precip_name

        x = data.variables['x_loc'][:]
        y = data.variables['y_loc'][:]
        if y[0] < 0:
            pass  # Check if y[0] = -ve if not reverse...  arr[::-1]
        else:
            y = data.variables['y_loc'][::-1]  # Check if y[0] = -ve if not reverse...  arr[::-1]
        Z = data.variables[precip_name][:]
        #print x
        #print y
        #print Z[0]
        if self.verbose: print self.Gauge_LOC_points[0]
        # and then do the interpolation
        values = interpolate2d(x,y,Z,self.Gauge_LOC_points) # This is a numpy array of Rain for this time slice at each of the gauge locations
        values.tolist()      # Convert array to list
        if self.verbose: print 'First values...'
        #print values
        self.all_values.append(values.tolist() ) # This is a time history of the List above
        if self.verbose: print 'First ALL values...'
        #print ALL_values
        # Save it to a file...??
        #raw_input('Hold at Gauge Extract....line 278')
        return


    def radar_plot_save_and_show(self, plot_line):
        # SHOW and SAVE the Plot            
    
        precip = self.precip
    
        filename = self.filename
        x = self.x
        y = self.y
        xl = self.xl
        yl = self.yl
        
        plt.figure(1)
        plt.clf()
        
        if self.plot_line : plt.plot( xl, yl,'--w')
        
        plt.suptitle('RADAR RAINFALL'+filename[-20:-3], size=20)
        s_title = 'Max Rainfall = %.1f mm / period' % (np.max(precip))
        plt.title(s_title , size=12)
        plt.imshow(precip, origin='lower', interpolation='bicubic',
                   extent=(x.min(), x.max(), y.min(), y.max()),vmin=0, vmax=Daily_plot_Vmax)
        plt.colorbar()            
        plt.show()
        plt.savefig('RADAR RAINFALL'+filename[-20:-3]+'.jpg',format='jpg')
        return



    def radar_plot_save_noshow(self):
        # DONT SHOW only SAVE the Plot
        
        precip = self.precip
        filename = self.filename
        x = self.x
        y = self.y
        xl = self.xl
        yl = self.yl
        
        if self.verbose: print 'Saving Image...'   
        plt.figure(1)
        plt.clf()
        if self.plot_line : plt.plot( xl, yl,'--w')
        plt.suptitle('RADAR RAINFALL'+ filename[-20:-3], size=20)
        s_title = 'Max Rainfall = %.1f mm / period' % (np.max(self.precip))
        plt.title(s_title , size=12)
        plt.imshow(self.precip, origin='lower', interpolation='bicubic',
                   extent=(x.min(), x.max(), y.min(), y.max()),vmin=0, vmax=Daily_plot_Vmax)
        plt.colorbar()            
        plt.savefig('RADAR RAINFALL'+self.filename[-20:-3]+'.jpg',format='jpg')
        return

    def radar_plot_show_only(self, id):
        """
        Plot radar data at timestep id
        """

    
        precips = self.precips
        base_filename = self.base_filename
        times = self.times
        x = self.x
        y = self.y
        xl = self.xl
        yl = self.yl
        plot_line = True
            
        plt.figure(1)
        plt.clf()
        
        from datetime import datetime
        date = datetime.utcfromtimestamp(times[id])
        
        if plot_line : plt.plot( xl, yl,'--w')
        plt.suptitle('RADAR RAINFALL '+date.strftime("%Y/%m/%d %H:%M"), size=20)
        s_title = 'Max Rainfall = %.1f mm / period' % (np.max(precips[id]))
        plt.title(s_title , size=12)
        plt.imshow(precips[id], origin='lower', interpolation='bicubic',
                   extent=(x.min(), x.max(), y.min(), y.max()),vmin=0, vmax=Daily_plot_Vmax)

        plt.colorbar()            
        plt.draw()
            
        return

    def plot_time_hist_radar_at_gauge_loc(self):
        """
        NOW PLOT THE DATA EXTRACTED AT GAUGE LOCATIONS 
        """
        
        precip = self.precip
        filename = self.filename
        x = self.x
        y = self.y
        xl = self.xl
        yl = self.yl
        file_counter = self.file_counter
        gauge_loc_labels = self.gauge_loc_labels
        time_step = self.time_step
        
        gauge_count = 0
        t = np.arange(file_counter)
        for item in sorted(gauge_loc_labels):
            # Access List of Lists [rows][cols]
            # need to plot gauge in columns
            #      loc1 = [[0.0 for y in range(ncols)] for x in range(nrows)]
            #bar_value = zip(*ALL_values)
            bar_values = [x[gauge_count] for x in self.all_values]
            total_rain = sum(bar_values)
            Ave_rain = total_rain/len(bar_values)
            max_Intensity = max(bar_values)/time_step*60.0
            b_title = 'Tot rain = %.1f mm/hr, Ave. rain = %.1f mm, Max Int.= %.1f mm' % (total_rain,Ave_rain,max_Intensity)
            #   Using list comprehension  b=[x[0] for x in a]
            print len(t)
            print len(bar_values)
            # Think about using...  zip(*lst)
            plt.bar(t,bar_values)
            plt.suptitle(' Rain Gauge data for Station %s:' % item, fontsize=14, fontweight='bold')    
            plt.title(b_title)
    
            plt.xlabel('time steps')
            plt.ylabel('rainfall (mm)')
            plt.show()
            gauge_count+=1
        return


    def plot_radar_accumulated(self):
        
        precip = self.precip
        precip_total = self.precip_total
        file_counter = self.file_counter
        time_step = self.time_step
        fromdir = self.fromdir
        xl = self.xl
        yl = self.yl
        x = self.x
        y = self.y
        
        plot_line = self.plot_line
        
        
        if self.verbose: 
            print precip_total
            print 'maximum rain =',np.max(precip_total)
            print 'mean rain =',np.mean(precip_total)
            
        Total_Rain_Vol = np.mean(precip_total)/1000.0*128.0*128.0 # Volume in Million m3 over 128km x 128km area
        Rain_Max_in_period = np.max(precip)
        Peak_Intensity = Rain_Max_in_period/time_step*60
        
        if self.verbose:
            print 'Total rainfall volume in Mill m3 =',Total_Rain_Vol
            print 'Peak rainfall in 1 time step = ', Rain_Max_in_period
            print 'Peak Intensity in 1 timestep =',Peak_Intensity
            
        dir_part = os.path.basename(os.path.normpath(fromdir))
        plot_sup_title = dir_part
        plt.figure(1)
        plt.clf()
        plt.suptitle('Accumulated '+str(file_counter)+' files '+plot_sup_title, size=20)
        
        
        s_title = 'Max Int. = %.1f mm/hr, Ave. rain = %.1f mm, Tot rain Vol. = %.3f Mill. m3' % (Peak_Intensity,np.mean(precip_total),Total_Rain_Vol)
        plt.title(s_title , size=12)
        
        plt.imshow(precip_total, origin='lower', interpolation='bicubic',extent=(x.min(), x.max(), y.min(), y.max()))
        if plot_line : plt.plot( xl, yl,'--w')
        plt.colorbar()
        plt.show()
        
        return
    
    
    
    
if __name__ == "__main__":
    
    act_rain = Grid_rain_radar()
    
    BASE_DIR = "/home/steve/RAINFALL/RADAR/AUS_RADAR/Calibrated_Radar_Data/ACT_netcdf"
    
    RADAR_DIR          = join(BASE_DIR, 'RADAR_Rainfall/140/2012/03/' )
    Radar_UTM_LOC_file = join(BASE_DIR, 'ACT_Bdy/Captains_Flat_RADAR_UTM_55.csv' )
    Radar_LL_LOC_file  = join(BASE_DIR, 'ACT_Bdy/Captains_Flat_RADAR_lat_long.csv' )
    Catchment_file     = join(BASE_DIR,'ACT_Bdy/ACT_Entire_Catchment_Poly_Simple_UTM_55.csv')
    
    #RainGauge_LOC_file = join(BASE_DIR, 'ACT_Bdy/Rain_Gauge_Station_Location_subset.csv' )
    
    print RADAR_DIR
    
    p0 = anuga.read_polygon(Radar_UTM_LOC_file)
    print 'radar utm loc'
    print Radar_UTM_LOC_file
    print p0
    
    p1 = anuga.read_polygon(Radar_LL_LOC_file)
    print 'radar ll loc'
    print Radar_LL_LOC_file
    print p1
    
    p2 = anuga.read_polygon(Catchment_file)
    print 'radar ll loc'
    print Catchment_file
    print p2
    
    import csv 

#     print 'raingauge loc'
#     print RainGauge_LOC_file
#     p2id = open(RainGauge_LOC_file)
#     rows = csv.reader(p2id)
#     for row in  rows:
#         print row

    act_rain.read_radar_UTM_offsets(Radar_UTM_LOC_file)
    
    print 'offsets'
    print act_rain.offset_x
    print act_rain.offset_y
    
    
    act_rain.process_calibrated_radar_data(RADAR_DIR)
    
    print 'offsets'
    print act_rain.offset_x
    print act_rain.offset_y
    
    p2 = np.array(p2)
    
    
        
    act_rain.xl = p2[:,0]
    act_rain.yl = p2[:,1]
    
    import pdb
    pdb.set_trace()    
    
    #print act_rain.all_times
    import time
    
    pl.ion()
    for id in xrange(len(act_rain.times)):
        act_rain.radar_plot_show_only(id)
        time.sleep(0.05)
    
    pl.ioff()
    pl.draw()
    
    
    import pdb
    pdb.set_trace()
    
    
    
    
     


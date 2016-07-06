
"""
Basic helper routines


"""

from anuga.fit_interpolate.interpolate2d import interpolate2d
import pylab as pl
import matplotlib.pyplot as plt
import numpy as np
import os


Daily_plot_Vmax = 15

# ----------------------------------------------------------------------------------------


class Rain_radar(object):
    
    def __init__(self,x,y,UTM_x,UTM_y,precip):
        
        self.precip = precip
        self.UTM.x = UTM_x
        self.UTM.y = UTM_y
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
        
        # Plotting parameters
        self.plot_line = False
        
        
    
    def read_radar_UTM_offsets(self, Radar_UTM_LOC_file):
        """
        CONVERT RADAR FROM LAT LONG TO UTM 
        
        PROCESS THE Lat Long to Produce UTM Grid of Radar ??
        """
        
        fid = open(Radar_UTM_LOC_file)
        lines = fid.readlines()  # Read Entire Input File
        fid.close()  
               
        for line in lines:
            if self.verbose: print line
            line=line.strip('\n')
            fields = line.split(',')
            offset_x = float(fields[0])
            offset_y = float(fields[1])
        if self.verbose: print offset_x,offset_y
        #raw_input('Hold here... line 128')
        
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

    def radar_plot_show_only(self):
        # ONLY SHOW the PLOT DONT SAVE
    
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
            
        return()

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
    return()


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
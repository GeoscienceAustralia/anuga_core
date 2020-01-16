
"""
Basic helper routines


"""
from __future__ import print_function
from __future__ import division

from builtins import range
from builtins import object
from past.utils import old_div
import anuga
#from anuga.fit_interpolate.interpolate2d import interpolate2d
from anuga.fit_interpolate.interpolate2d import interpolate_raster
import pylab as pl
import matplotlib.pyplot as plt
import numpy as np
import os
from os.path import join
import fnmatch
import gzip
from anuga.file.netcdf import NetCDFFile


#-----------------------------------------------------------------------------------------

class Raster_time_slice_data(object):
    
    
    def __init__(self, 
                 start_time = None,
                 final_time = None,
                 verbose=False, 
                 debug=False):
        """
        start_time: seconds since epoch  or string of form 20120229_1210
        final_time: seconds since epoch  or string of form 20120229_1210
        
        The data is assumed to be stored as a raster, ie. columns  
        in the x direction (eastings) and rows in the vertical y direction
        (northings) from north to south. 
        
        It is assumed that the data is in SI units. 
        
        """
        
        self.x = None
        self.y = None
        self.times = []
        self.time_step = 0.0
        
        self.verbose = verbose
        self.debug = debug
        
        self.start_time = anuga.parse_time(start_time)
        self.final_time = anuga.parse_time(final_time)
        
        
        
            

    def get_extent(self):
        
        return self.extent
    

    def read_data_files(self):
        """
        Implement method to read in specific data formats. 
        
        Convert to SI units, ie convert mm to metres when importing rain data.
        """
        
        pass


    def ungzip_data_files(self, data_dir):
        """
        Given a data_dir walk through all sub directories to find 
        gzipped files and unzip
        """
        
        pattern = '*.gz'
        
        self.radar_dir = data_dir
        
        for root, dirs, files in os.walk(data_dir): 
            
            if self.debug or self.verbose:
                print('Directory: ',dirs)
                print('Root: ',root)
                print('Number of Files = ',len(fnmatch.filter(files, pattern)))
            
            if len(fnmatch.filter(files, pattern)) > 0:
                os.chdir(root)
                os.system('gzip -d *.gz')
                
             



    def extract_data_at_locations(self, locations):
        """
        
        Extract data from Rasters at locations  
        
        """
        if self.verbose or self.debug: print('Extract data at locations', locations)
        
        locations = np.array(locations)
        all_values = []
        x = self.x
        y = self.y
        data_slices = self.data_slices

        if self.debug: print(locations)
        
        for data in data_slices:
            # and then do the interpolation
            #values = interpolate2d(x,y,np.fliplr(precip.T),locations) 
            values = interpolate_raster(x,y,data,locations)      
            all_values.append(values) 


        return np.array(all_values)



    def grid_indices_inside_polygon(self, polygon=None):
        """
        Flattened indices of grid point inside a polygon
        
        To get corresponding grid values of time slice tid use
        
        data_inside = self.data_slices[tid].flat[indices]
        """
        
        
        dx = self.extent[1]-self.extent[0]
        dy = self.extent[3]-self.extent[2]
        x = self.x
        y = self.y
        nx = len(x)
        ny = len(y)
        ldx = old_div(dx,nx)
        ldy = old_div(dy,ny)
        
        if not polygon is None:
            X,Y = np.meshgrid(x,y)
            Y = np.flipud(Y)
            
            points = np.empty((nx*ny,2), dtype="float")
            points[:,0] = X.flatten()
            points[:,1] = Y.flatten()
            from anuga import inside_polygon
            indices = inside_polygon(points,polygon)
        else:
            indices = None
            
        return indices
    

    def accumulate_data_stats(self, tid=None, polygon=None, print_stats=False):
        """
        Accumulate stats of data over slices, either for a specified tid timeslice
        or over all time slices.
        
        Can be restricted to a polygon 
        """
        
        dx = self.extent[1]-self.extent[0]
        dy = self.extent[3]-self.extent[2]
        x = self.x
        y = self.y
        nx = len(x)
        ny = len(y)
        ldx = old_div(dx,nx)
        ldy = old_div(dy,ny)
        
        time_step = self.time_step
        
        
        if tid is None:
            data = self.data_accumulated
            time_period = self.time_step*len(self.times)
        else:
            data = self.data_slices[tid]
            time = self.times[tid]
            time_period = self.time_step
            
        
        indices = self.grid_indices_inside_polygon(polygon)

        
        if indices is None:  
            data_max_in_period = np.max(data) 
            total_data_volume = np.sum(data)*ldx*ldy 
            catchment_area = len(data.flat)*ldx*ldy
        else:
            pmask = data.flat[indices]
            data_max_in_period = np.max(pmask)
            total_data_volume = np.sum(pmask)*ldx*ldy 
            catchment_area = len(pmask)*ldx*ldy

        #print indices
        peak_intensity = old_div(data_max_in_period,time_step)      

        
        if print_stats or self.verbose:
            print('Time period = ',time_period)
            print('Time step = ',time_step)
            print('Catchment Area = ', catchment_area)
            print('Total rainfall volume in cubic metres  =', total_data_volume)
            print('Peak data/time_step in time period  = ', data_max_in_period)
            print('Peak Intensity in time period (m/sec) =', peak_intensity)      
        
        return total_data_volume, data_max_in_period, peak_intensity, catchment_area, time_period


    def plot_data(self, tid, plot_vmax=None, save=False, show=True, polygons=None):
        """
        Plot data at timestep tid
        """

        data_slices = self.data_slices
        times = self.times
        x = self.x
        y = self.y
        

        # Aliases

        
        # get UTC time from epoch time
        from datetime import datetime
        date_time = datetime.utcfromtimestamp(times[tid]).strftime("%Y/%m/%d %H:%M")
        
        if self.debug: print('--- Date/Time ', date_time)
        
        plt.figure(1)
        plt.clf()

        
        if not polygons is None:
            for polygon in polygons:
                if not polygon is None:
                    polygon = np.array(polygon)
                    plt.plot( polygon[:,0], polygon[:,1],'--w')
            
        plot_title = 'RASTER DATA '+date_time
        plt.suptitle(plot_title, size=20)
        s_title = 'Max Data = %.2e / period' % (np.max(data_slices[tid]))
        plt.title(s_title , size=12)
        
        plt.imshow(data_slices[tid], origin='upper', interpolation='bicubic',
                   extent=self.extent,vmin=0, vmax=plot_vmax)
        
        
        plt.colorbar()            
        
        if show: plt.draw()
        if save: plt.savefig(plot_title+'.jpg',format='jpg')
            
        return




    def plot_accumulated_data(self, polygons=None):
        
        data_accumulated = self.data_accumulated

        time_step = self.time_step

        x = self.x
        y = self.y
        
        
        if self.debug: 
            print(data_accumulated)
            print('maximum rain =',np.max(data_accumulated))
            print('mean rain =',np.mean(data_accumulated))
            
        dx = self.extent[1]-self.extent[0]
        dy = self.extent[3]-self.extent[2]
        total_data_vol = old_div(np.mean(data_accumulated)*dx*dy,1e6) # Volume in Million m3 over 128km x 128km area
        data_max_in_period = self.data_max_in_period
        peak_intensity = old_div(data_max_in_period,time_step)
        extent = self.extent
        
        if self.verbose:
            print('Total data volume =',total_data_vol)
            print('Peak data in 1 time step = ', data_max_in_period)
            print('Peak Intensity in 1 timestep =', peak_intensity)
            print('extent', extent)
            print('size', extent[1]-extent[0], extent[3]-extent[2])
            print(time_step)


        plt.figure(1)
        plt.clf()
        plt.suptitle('Accumulated Data', size=20)
        
        
        s_title = 'Max Int. = %.2e, Average = %.2e, Tot Vol. = %.2e Mill. m3' \
                  % (peak_intensity,np.mean(data_accumulated),total_data_vol)
        plt.title(s_title , size=12)
        
        
        if not polygons is None:
            for polygon in polygons:
                if not polygon is None:
                    polygon = np.array(polygon)
                    plt.plot( polygon[:,0], polygon[:,1],'--w')
                
        plt.imshow(data_accumulated, origin='upper', interpolation='bicubic',
                   extent=(x.min(), x.max(), y.min(), y.max()))
        
        
        plt.colorbar()
        plt.draw()
        
        return
    
    def plot_time_hist_locations(self, locations):
        """
        Plot time histograms at selected locations i.e.
        gauge locations
        """
        

        t = (self.times-self.times[0])
        
        locations = np.array(locations)

        time_step = self.time_step
        
        all_values = self.extract_data_at_locations(locations)

        for lid, _ in enumerate(locations):

            bar_values = [values[lid] for values in all_values]
            total_values = sum(bar_values)
            average_values = old_div(total_values,(self.times[-1]-self.times[0]))
            max_intensity = old_div(max(bar_values),time_step)
            
            
            b_title = 'Total = %.2e, Average = %.2e, Max Int.= %.2e' % (total_values,average_values,max_intensity)

            plt.bar(t,bar_values,width=time_step,)
            plt.suptitle(' Data for Location %s:' % lid, fontsize=14, fontweight='bold')    
            plt.title(b_title)
    
            plt.xlabel('time (s)')
            plt.ylabel('Data')
            plt.show()
            
        return

    

# ----------------------------------------------------------------------------------------



class Calibrated_radar_rain(Raster_time_slice_data):
    
    def __init__(self, 
                 radar_dir = None,
                 start_time = None,
                 final_time = None,
                 verbose=False, 
                 debug=False):
        """
        start_time: seconds since epoch  or string of form 20120229_1210
        final_time: seconds since epoch  or string of form 20120229_1210
        
        The BoM data is assumed to be stored as a raster, ie. columns  
        in the x direction (eastings) and rows in the vertical y direction
        (northings) from north to south. 
        
        The data is stored in mm so we have to convert to metres. 
        
        """


        Raster_time_slice_data.__init__(self,
                                        start_time = start_time,
                                        final_time = final_time,
                                        verbose = verbose,
                                        debug = debug)

        self.radar_dir = radar_dir
        
        # process the radar files
        if not radar_dir is None:
            self.read_data_files(radar_dir)
            

    
    def read_data_files(self, radar_dir, pattern = '*.nc'):
        """
        Given a radar_dir walk through all sub directories to find 
        radar raster files
        """
        
        
        file_counter = 0
        first = True
        data_max_in_period = 0.0
        
        data_slices = []
        times = []
        
        self.radar_dir = radar_dir
        
        if self.verbose:
            import datetime
            print("READING BoM calibrated rain grid data")
        
        for root, dirs, files in os.walk(radar_dir): 
            
            if self.debug:
                print('Directory: ',dirs)
                print('Number of Files = ',len(fnmatch.filter(files, pattern)))
            
            for filename in fnmatch.filter(files, pattern): 
                    
                key = filename[-16:-3]
                
                valid_time = anuga.parse_time(key)
                
                if not self.final_time is None:
                    if valid_time > self.final_time: continue
                if not self.start_time is None:                    
                    if valid_time < self.start_time: continue

                if self.debug :
                    print(filename)
                    
                file_counter +=1
                if file_counter == 1:
                    self.radar_data_title = "RADAR_Data_"+filename[-20:-3]


                data = NetCDFFile(os.path.join(root, filename), 'r') 
                
                # RADAR NetCDF files have Dimensions, Attributes, Variables
                if self.debug:
                    print('VARIABLES:')
                    print('Reference LAT, LONG = ',data.reference_longitude, data.reference_latitude)

                    
                # Check Time for key 
                valid_time = anuga.parse_time(key)
             
                
                file_start_time = data.variables['start_time'][0]
                file_valid_time = data.variables['valid_time'][0]
                
                new_time_step = file_valid_time - file_start_time
                
                if self.debug:
                    print('VARIABLES:')
                    print('Reference times ', file_start_time, file_valid_time, valid_time)
                
                times.append(valid_time)                
        
                    
                # This handles format changes in the files from BOM !!!!
                possible_precip_names = ['precipitation',  'precip', 'rain_amount'] 
                
                # Go through each of the possible names
                for name in possible_precip_names:  # Check if name is a key in the variables dictionary
                    if name in data.variables:
                        precip_name = name
                        if self.debug:
                            print('BOM Reference name tag in this file:')
                            print(precip_name)

                
                if first:
                    first = False
                    
                    self.base_filename = filename[:-20]
                    
                    self.time_step = file_valid_time - file_start_time
                    
                    if self.debug: print(' Accumulate rainfall here....')
                    self.x = data.variables['x_loc'][:]
                    self.y = data.variables['y_loc'][:]
                    if self.y[0] < 0:
                        pass # Check if y[0] = -ve if not reverse...  arr[::-1]
                    else:
                        self.y = self.y[::-1]  # Check if y[0] = -ve if not reverse...  arr[::-1]
                    
                    self.reference_longitude = data.reference_longitude
                    self.reference_latitude =  data.reference_latitude
                    # Convert to UTM offsets
                    
                    self.zone, self.offset_x, self.offset_y = anuga.LLtoUTM(self.reference_latitude, self.reference_longitude)
                    
                    # Convert to UTM
                    self.x = self.x*1000 + self.offset_x
                    self.y = self.y*1000 + self.offset_y
                    
                    data_slice = old_div(data.variables[precip_name][:],1000)  # convert from mm to m
                    
                    data_accumulated = data_slice.copy() # Put into new Accumulating ARRRAY
                    
                    data_max_in_period  = max(np.max(data_slice),data_max_in_period)  
                        
    
                else:  # ---If NOT FIRST !!!
                    data_slice = old_div(data.variables[precip_name][:],1000) # convert from mm to m

                    data_accumulated += data_slice
                    
                    if self.debug: print(' Keep accumulating rainfall....')

                    data_max_in_period  = max(np.max(data_slice),data_max_in_period)
                    
                    assert np.allclose(self.time_step,new_time_step), "Timesteps not equal"

                    
                data_slices.append(data_slice)
        
        
        self.data_max_in_period = data_max_in_period
                  
        times = np.array(times)
        ids = np.argsort(times)
        
        #pdb.set_trace()
        if len(times) > 0:
            self.times = times[ids]   
            self.data_slices = np.array([ data_slices[tid] for tid in ids ])
            self.start_time = self.times[0]-self.time_step
            self.data_accumulated = data_accumulated 
            print("+++++", np.sum(data_accumulated))
            # Test sorting
            for i, tid in enumerate(ids):
                np.allclose(times[tid], self.times[i])
                np.allclose(data_slices[tid], self.data_slices[i])
        else:
            self.times = []  
            self.data_slices = []
            self.start_time = []
            self.data_accumulated = []
        
                
        self.extent = (self.x.min(), self.x.max(), self.y.min(), self.y.max())

        if self.verbose:
            print("    From UTC time: %s"% datetime.datetime.utcfromtimestamp(self.start_time).strftime('%c'))
            print("    To UTC time:   %s"% datetime.datetime.utcfromtimestamp(self.final_time).strftime('%c'))
            print("    Read in %g time slices" % len(self.times))

    
    
if __name__ == "__main__":
    
    try:
        import ipdb as pdb
    except:
        import pdb
    

    scenario = 'act'
    #scenario = 'grantham'
    #scenario = 'artificial'
    

    
    if scenario == 'act':
        start_time = '20120301_0000'
        final_time = '20120302_1200'
        BASE_DIR = "/home/steve/RAINFALL/RADAR/AUS_RADAR/Calibrated_Radar_Data/ACT_netcdf"
        RADAR_DIR          = join(BASE_DIR, 'RADAR_Rainfall/140/2012' )    
        Catchment_file     = join(BASE_DIR,'ACT_Bdy/ACT_Entire_Catchment_Poly_Simple_UTM_55.csv')
        State_boundary_file = join(BASE_DIR,'ACT_Bdy/ACT_State_Bdy_UTM_55.csv')
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
        
        BASE_DIR = '/home/steve/RAINFALL/RADAR/AUS_RADAR/Gauge_Blend_Grantham/20081211-20110223.gz_y_loc_Inverted_30min/Merged/66/2011/'
        RADAR_DIR          = join(BASE_DIR, '01' ) 
        Daily_plot_Vmax = 50
        rain = Calibrated_radar_rain(RADAR_DIR, start_time=start_time, final_time=final_time, verbose=True, debug=True)
        l0 = [476160., 6889300.0]
        l1 = [476160., 6978420.0]
        l2 = [421517., 6974520.0]
        l3 = [479412.0, 6980370.0]
        locations = [l0,l1,l2,l3]

        
    elif scenario == 'artificial':
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
     
    p2 = [[689801.0, 6091501.0], [700200.0, 6091501.0], [700200.0, 6104600.0], [689801.0, 6104600.0], [689801.0, 6091501.0]]

    print(rain.get_extent())
    rain.accumulate_data_stats()
    
    print('time_step ',rain.time_step)     
        
    import time
    pl.ion()
    #pdb.set_trace() 

    plot_vmax = np.max(rain.data_slices)
    print('plot_vmax', plot_vmax)
    for tid in range(len(rain.times)):
        rain.plot_data(tid, plot_vmax=plot_vmax, save=False, show=True, polygons=[p2,p3])
        time.sleep(0.05)
        #ipdb.set_trace() 
        
        
    
    rain.plot_accumulated_data(polygons=[p2])
    
    pl.ioff()
    pl.show()
     

    rain.plot_time_hist_locations(locations) 
    
    p_indices = rain.grid_indices_inside_polygon(polygon=p2)
    
    
    
    pdb.set_trace()    
    
    
    
     


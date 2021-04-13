
"""
Basic helper routines


"""
from __future__ import print_function
from __future__ import division

from builtins import range
from past.builtins import basestring
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


Daily_plot_Vmax = 4

# ----------------------------------------------------------------------------------------


class Calibrated_radar_rain(object):

    def __init__(self,
                 radar_dir=None,
                 start_time=None,
                 final_time=None,
                 verbose=False,
                 debug=False):
        """
        start_time: seconds since epoch  or string of form 20120229_1210
        final_time: seconds since epoch  or string of form 20120229_1210

        The BoM data is assumed to be stored as a raster, ie. columns
        in the x direction (eastings) and rows in the vertical y direction
        (northings) from north to south.
        """
        self.radar_dir = radar_dir

        self.x = None
        self.y = None
        self.precips = []
        self.times = []
        self.precip_total = 0.0
        self.time_step = 0.0

        self.verbose = verbose
        self.debug = debug

        self.start_time = self.parse_time(start_time)
        self.final_time = self.parse_time(final_time)

        self.pattern = '*.nc'

        # process the radar files
        if not radar_dir is None:
            self.read_data_files(radar_dir)

    def get_extent(self):

        return self.extent

    def parse_time(self, time=None):
        """
        Time: seconds since epoch  or
        string of form '20120229'  '20120229_1210' '20120229 1210' '201202291210'
        """

        if time is None:
            return None

        if not isinstance(time, basestring):

            try:
                time = float(time)
                return time
            except ValueError:
                pass

        year, month, day, hour, minute = 1970, 1, 1, 0, 0

        try:
            year = int(time[0:4])
        except:
            year, month, day = 1970, 1, 1

        #month = int(time[4:6])

        try:
            month = int(time[4:6])
        except:
            month = 1

        try:
            day = int(time[6:8])
        except:
            day = 1

        try:
            dash = time[8:9]
            assert dash == '_' or dash == ':' or dash == '/' or dash == ' '
        except:
            dash = None

        try:
            if dash is None:
                hour = int(time[8:10])
            else:
                hour = int(time[9:11])
        except:
            hour = 0

        try:
            if dash is None:
                minute = int(time[10:12])
            else:
                minute = int(time[11:13])
        except:
            minute = 0

        if self.debug:
            print(year, month, day, hour, minute)
            print('Convert to epoch')
 
                
                    
        import datetime
        time = int((datetime.datetime(year,month,day,hour,minute) - datetime.datetime(1970,1,1)).total_seconds())
 
        if self.debug: print(time)
        
        return float(time)

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
        reverse = False

        self.radar_dir = radar_dir

        for root, dirs, files in os.walk(radar_dir):

            if self.debug or self.verbose:
                print('Directory: ',dirs)
                print('Number of Files = ',len(fnmatch.filter(files, pattern)))
            
            for filename in fnmatch.filter(files, pattern): 
                    
                key = filename[-16:-3]

                valid_time = self.parse_time(key)

                if not self.final_time is None:
                    if valid_time > self.final_time:
                        continue
                if not self.start_time is None:
                    if valid_time < self.start_time:
                        continue

                if self.debug :
                    print(filename)
                    
                file_counter +=1
                if file_counter == 1:
                    self.radar_data_title = "RADAR_Data_"+filename[-20:-3]

                if pattern == '*.gz':
                    filename = gzip.open(filename, 'rb')
                    data = NetCDFFile(os.path.join(root, filename[:-3]), 'r')
                else:
                    data = NetCDFFile(os.path.join(root, filename), 'r')

                # RADAR NetCDF files have Dimensions, Attributes, Variables
                if self.debug:
                    print('VARIABLES:')
                    print('Reference LAT, LONG = ',data.reference_longitude, data.reference_latitude)

                # Check Time for key
                valid_time = self.parse_time(key)

                file_start_time = data.variables['start_time'][0]
                file_valid_time = data.variables['valid_time'][0]

                if self.debug:
                    print('VARIABLES:')
                    print('Reference times ', file_start_time, file_valid_time, valid_time)
                
                times.append(valid_time)                
        
                    
                # This handles format changes in the files from BOM !!!!
                possible_precip_names = [
                    'precipitation',  'precip', 'rain_amount']

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
                    
                    
                    if self.debug: print(' Accumulate rainfall here....')
                    self.x = data.variables['x_loc'][:]
                    self.y = data.variables['y_loc'][:]
                    if self.y[0] < 0:
                        # Check if y[0] = -ve if not reverse...  arr[::-1]
                        pass
                    else:
                        # Check if y[0] = -ve if not reverse...  arr[::-1]
                        self.y = self.y[::-1]

                    self.reference_longitude = data.reference_longitude
                    self.reference_latitude = data.reference_latitude
                    # Convert to UTM offsets

                    self.zone, self.offset_x, self.offset_y = anuga.LLtoUTM(
                        self.reference_latitude, self.reference_longitude)

                    # Convert to UTM
                    self.x = self.x*1000 + self.offset_x
                    self.y = self.y*1000 + self.offset_y

                    precip = data.variables[precip_name][:]

                    precip_total = precip.copy()  # Put into new Accumulating ARRRAY

                    rain_max_in_period = max(
                        np.max(precip), rain_max_in_period)

                else:  # ---If NOT FIRST !!!
                    precip = data.variables[precip_name][:]

                    precip_total += precip
                    
                    if self.debug: print(' Keep accumulating rainfall....')

                    rain_max_in_period = max(
                        np.max(precip), rain_max_in_period)

                precips.append(precip)

        self.rain_max_in_period = rain_max_in_period

        times = np.array(times)
        ids = np.argsort(times)

        pdb.set_trace()

        self.times = times[ids]
        self.precips = np.array([precips[tid] for tid in ids])
        self.time_step = self.times[1]-self.times[0]
        self.start_time = self.times[0]-self.time_step

        self.precip_total = precip_total

        # Test sorting
        for i, tid in enumerate(ids):
            np.allclose(times[tid], self.times[i])
            np.allclose(precips[tid], self.precips[i])

        self.extent = (self.x.min(), self.x.max(), self.y.min(), self.y.max())

    def ungzip_data_files(self, radar_dir):
        """
        Given a radar_dir walk through all sub directories to find
        radar files
        """

        pattern = '*.gz'

        self.radar_dir = radar_dir

        for root, dirs, files in os.walk(radar_dir):

            if self.debug or self.verbose:
                print('Directory: ',dirs)
                print('Root: ',root)
                print('Number of Files = ',len(fnmatch.filter(files, pattern)))
            
            if len(fnmatch.filter(files, pattern)) > 0:
                os.chdir(root)
                os.system('gzip -d *.gz')

    def extract_data_at_locations(self, locations):
        """

        EXTRACT RAINFALL FROM GRID AT locations

        """
        if self.verbose or self.debug: print('Extract Rain Check data First')
        
        locations = np.array(locations)
        all_values = []
        x = self.x
        y = self.y
        precips = self.precips

        if self.debug: print(locations)
        
        for precip in precips:
            # and then do the interpolation
            #values = interpolate2d(x,y,np.fliplr(precip.T),locations)
            values = interpolate_raster(x, y, precip, locations)
            all_values.append(values)

        return np.array(all_values)

    def grid_indices_inside_polygon(self, polygon=None):
        """
        flattened indices of grid point inside a polygon

        To get corresponding grid values of time slice tid use

        precip_inside = self.precips[tid].flat[indices]
        """

        dx = self.extent[1]-self.extent[0]
        dy = self.extent[3]-self.extent[2]
        x = self.x
        y = self.y
        nx = len(x)
        ny = len(y)
        ldx = dx/nx
        ldy = dy/ny
        
        if not polygon is None:
            X, Y = np.meshgrid(x, y)
            points = np.empty((nx*ny, 2), dtype="float")
            points[:, 0] = X.flatten()
            points[:, 1] = Y.flatten()
            from anuga import inside_polygon
            indices = inside_polygon(points, polygon)
        else:
            indices = None

        return indices

    def accumulate_grid(self, tid=None, polygon=None):
        """
        Accumulate precip over grid, either for a specified tid timeslice
        or over all time slices.

        Can be restricted to a polygon
        """

        dx = self.extent[1]-self.extent[0]
        dy = self.extent[3]-self.extent[2]
        x = self.x
        y = self.y
        nx = len(x)
        ny = len(y)
        ldx = dx/nx
        ldy = dy/ny

        time_step = self.time_step

        if tid is None:
            precip = self.precip_total
            time_period = self.times[-1]-self.times[0]
        else:
            precip = self.precips[tid]
            time = self.times[tid]
            time_period = self.time_step

        indices = self.grid_indices_inside_polygon(polygon)

        if indices is None:
            Rain_Max_in_period = np.max(precip)
            Total_Rain_Vol = np.sum(precip)/1000.0*ldx*ldy  # cubic metres
            import pdb
            pdb.set_trace()
            Catchment_Area = precip.size*ldx*ldy
        else:
            pmask = precip.flat[indices]
            Rain_Max_in_period = np.max(pmask)
            Total_Rain_Vol = np.sum(pmask)/1000.0*ldx*ldy  # cubic metres
            Catchment_Area = pmask.size*ldx*ldy

        print (indices)
        Peak_Intensity = Rain_Max_in_period/time_step  # mm/sec

        if self.verbose:
            print('Time period = ',time_period)
            print('Time step = ',time_step)
            print('Catchment Area = ', Catchment_Area)
            print('Total rainfall volume in cubic metres (m^3) =', Total_Rain_Vol)
            print('Peak rainfall/time_step in time period (mm) = ', Rain_Max_in_period)
            print('Peak Intensity in time period (mm/sec) =', Peak_Intensity)      
        
        return Total_Rain_Vol, Rain_Max_in_period, Peak_Intensity, Catchment_Area, time_period

        return Total_Rain_Vol, Rain_Max_in_period, Peak_Intensity, Catchment_Area, time_period

    def plot_grid(self, tid, save=False, show=True, polygons=None):
        """
        Plot radar data at timestep tid
        """

        # Aliases
        precips = self.precips
        times = self.times
        x = self.x
        y = self.y

        # get time from filenane
        from datetime import datetime
        date_time = datetime.utcfromtimestamp(times[tid]).strftime("%Y/%m/%d %H:%M")
        
        if self.debug: print('--- Date/Time ', date_time)
        
        plt.figure(1)
        plt.clf()

        if not polygons is None:
            for polygon in polygons:
                if not polygon is None:
                    polygon = np.array(polygon)
                    plt.plot(polygon[:, 0], polygon[:, 1], '--w')

        plot_title = 'RADAR RAINFALL '+date_time
        plt.suptitle(plot_title, size=20)
        s_title = 'Max Rainfall = %.1f mm / period' % (np.max(precips[tid]))
        plt.title(s_title, size=12)

        plt.imshow(precips[tid], origin='upper', interpolation='bicubic',
                   extent=self.extent, vmin=0, vmax=Daily_plot_Vmax)

        #plt.contourf(precips[tid], origin='lower', interpolation='bicubic',
        #           extent=self.extent, vmin=0, vmax=Daily_plot_Vmax)

        plt.colorbar()

        if show:
            plt.pause(0.001)
        if save:
            plt.savefig(plot_title+'.jpg', format='jpg')

        return

    def plot_grid_accumulated(self, polygons=None):

        precip_total = self.precip_total

        time_step = self.time_step
        radar_dir = self.radar_dir
        x = self.x
        y = self.y

        if self.debug:
            print (precip_total)
            print ('maximum rain =', np.max(precip_total))
            print ('mean rain =', np.mean(precip_total))

        dx = self.extent[1]-self.extent[0]
        dy = self.extent[3]-self.extent[2]
        
        Total_Rain_Vol = np.mean(precip_total)/1000.0*dx*dy/1e6 # Volume in Million m3 over 128km x 128km area
        Rain_Max_in_period = self.rain_max_in_period
        Peak_Intensity = old_div(Rain_Max_in_period,time_step)
        extent = self.extent

        if self.verbose:
            print('Total rainfall volume in Mill m3 =',Total_Rain_Vol)
            print('Peak rainfall in 1 time step = ', Rain_Max_in_period)
            print('Peak Intensity in 1 timestep =',Peak_Intensity)
            print('extent', extent)
            print('size', extent[1]-extent[0], extent[3]-extent[2])
            print(time_step)

        plt.figure(1)
        plt.clf()
        plt.suptitle('Accumulated Radar Rainfall', size=20)

        s_title = 'Max Int. = %.1f mm/hr, Ave. rain = %.1f mm, Tot rain Vol. = %.3f Mill. m3' \
                  % (Peak_Intensity, np.mean(precip_total), Total_Rain_Vol)
        plt.title(s_title, size=12)

        if not polygons is None:
            for polygon in polygons:
                if not polygon is None:
                    polygon = np.array(polygon)
                    plt.plot(polygon[:, 0], polygon[:, 1], '--w')

        plt.imshow(precip_total, origin='upper', interpolation='bicubic',
                   extent=(x.min(), x.max(), y.min(), y.max()))

        #plt.contourf(precip_total, origin='lower', interpolation='bicubic',
        #           extent=(x.min(), x.max(), y.min(), y.max()))

        plt.colorbar()
        #plt.show()
        plt.pause(0.01)

        return

    def plot_time_hist_locations(self, locations):
        """
        Plot time histograms at selected locations i.e.
        gauge locations
        """

        t = (self.times-self.times[0])/60.0  # minutes

        locations = np.array(locations)

        time_step = self.time_step

        all_values = self.extract_data_at_locations(locations)



        for lid, _ in enumerate(locations):

            bar_values = [values[lid] for values in all_values]
            total_rain = sum(bar_values)
            Ave_rain = total_rain/(self.times[-1]-self.times[0])*3600
            max_Intensity = max(bar_values)/time_step*3600

            b_title = 'Tot rain = %.1f mm, Ave. rain = %.3f mm/hr, Max Int.= %.3f mm/hr' % (
                total_rain, Ave_rain, max_Intensity)
            #   Using list comprehension  b=[x[0] for x in a]
            #print len(t)
            #print len(bar_values)
            # Think about using...  zip(*lst)
            #plt.figure(1)
            plt.clf()
            plt.bar(t, bar_values, width=time_step/60,)
            plt.suptitle(' Data for Location %s:' %
                         lid, fontsize=14, fontweight='bold')
            plt.title(b_title)

            plt.xlabel('time (mins)')
            plt.ylabel('rainfall (mm)')
            plt.pause(0.01)
            pdb.set_trace()

        return


if __name__ == "__main__":

    try:
        import ipdb as pdb
    except:
        import pdb

    scenario = 'act'
    #scenario = 'grantham'
    #scenario = 'artificial'
    
    HOME_DIR = "/home/anuga/RAINFALL"
    
    if scenario == 'act':
        start_time = '20120301_0000'
        final_time = '20120302_1200'
        BASE_DIR           = join(HOME_DIR, "RADAR/AUS_RADAR/Calibrated_Radar_Data/ACT_netcdf")
        RADAR_DIR          = join(BASE_DIR, 'RADAR_Rainfall/140/2012' )    
        Catchment_file     = join(BASE_DIR,'ACT_Bdy/ACT_Entire_Catchment_Poly_Simple_UTM_55.csv')
        State_boundary_file = join(BASE_DIR,'ACT_Bdy/ACT_State_Bdy_UTM_55.csv')
        Daily_plot_Vmax = 4
        rain = Calibrated_radar_rain(
            RADAR_DIR, start_time=start_time, final_time=final_time, verbose=True, debug=True)
        l0 = [670501.0, 6117750.0]
        l1 = [702251.0, 6157920.0]
        l2 = [748256., 5966120.]
        l3 = [600000., 5966120.]
        locations = [l0, l1, l2, l3]

    elif scenario == 'grantham':
        start_time = '20110109_2300'
        final_time = '20110110_2300'

        #start_time = '20110108_2300'
        #final_time = '20110109_2300'

        #start_time = '20110105_2300'
        #final_time = '20110106_2300'
        
        BASE_DIR           = join(HOME_DIR, 'RADAR/AUS_RADAR/Gauge_Blend_Grantham/20081211-20110223.gz_y_loc_Inverted_30min/Merged/66/2011/')
        RADAR_DIR          = join(BASE_DIR, '01' ) 
        Daily_plot_Vmax = 50
        rain = Calibrated_radar_rain(
            RADAR_DIR, start_time=start_time, final_time=final_time, verbose=True, debug=True)
        l0 = [476160., 6889300.0]
        l1 = [476160., 6978420.0]
        l2 = [421517., 6974520.0]
        l3 = [479412.0, 6980370.0]
        locations = [l0, l1, l2, l3]

    elif scenario == 'artificial':
        start_time = 4700
        final_time = 5500
        rain = Calibrated_radar_rain(
            start_time=start_time, final_time=final_time, verbose=True, debug=True)
        rain.x = np.arange(41)/40.0*4000.0 + 3000.0
        rain.y = np.arange(61)/60.0*4000.0 + 100000.0
        rain.extent = (rain.x.min(), rain.x.max(), rain.y.min(), rain.y.max())
        rain.times = np.arange(11)*600+4000
        rain.time_step = 600
        nx = len(rain.x)
        ny = len(rain.y)
        rain.precip_total = np.arange(nx*ny, dtype="float").reshape(ny, nx)
        rain.precips = [
            i*np.ones((ny, nx), dtype="float")+4 for i, time in enumerate(rain.times)]
        for i, precip in enumerate(rain.precips):
            # Rows vertical, Columns horizontal. from top
            precip[2*i:2*i+4, :] = 2
            precip[:, i:i+4] = 3

        rain.precip_total = np.arange(nx*ny, dtype="float").reshape(ny, nx)
        rain.precip_total[50:54, :] = 2
        rain.precip_total[:, 20:24] = 3
        rain.rain_max_in_period = np.max(rain.precip_total)
        rain.radar_dir = None
        l0 = [5100.0, 100000.0]
        l1 = [5100.0, 103900.0]
        l2 = [3400., 103400.]
        l3 = [6000., 103400.]
        locations = [l0, l1, l2, l3]

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

    print(rain.get_extent())
    rain.accumulate_grid()
    
    print('time_step ',rain.time_step)     
        
    import time
    #pl.ion()
    pdb.set_trace() 
    for tid in range(len(rain.times)):
        rain.plot_grid(tid, save=False, show=True, polygons=[p2])

        #time.sleep(0.05)
        #pdb.set_trace()


    pdb.set_trace()

    rain.plot_grid_accumulated(polygons=[p2])

    pdb.set_trace()

    #pl.ioff()
    #plt.pause(0.01)

    rain.plot_time_hist_locations(locations)

    pdb.set_trace()

    p_indices = rain.grid_indices_inside_polygon(polygon=p2)

    pdb.set_trace()

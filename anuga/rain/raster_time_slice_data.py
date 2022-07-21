
import anuga
from anuga.fit_interpolate.interpolate2d import interpolate_raster
import matplotlib.pyplot as plt
import numpy as np
import os
from os.path import join
import fnmatch


# -----------------------------------------------------------------------------------------

class Raster_time_slice_data(object):

    def __init__(self,
                 start_time=None,
                 final_time=None,
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
                print('Directory: ', dirs)
                print('Root: ', root)
                print('Number of Files = ', len(
                    fnmatch.filter(files, pattern)))

            if len(fnmatch.filter(files, pattern)) > 0:
                os.chdir(root)
                os.system('gzip -d *.gz')

    def extract_data_at_locations(self, locations):
        """
        Extract data from Rasters at locations  
        """
        
        if self.verbose or self.debug:
            print('Extract data at locations', locations)

        locations = np.array(locations)
        all_values = []
        x = self.x
        y = self.y
        data_slices = self.data_slices

        if self.debug:
            print(locations)

        for data in data_slices:
            # and then do the interpolation
            #values = interpolate2d(x,y,np.fliplr(precip.T),locations)
            values = interpolate_raster(x, y, data, locations)
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
        ldx = dx/nx
        ldy = dy/ny

        if not polygon is None:
            X, Y = np.meshgrid(x, y)
            Y = np.flipud(Y)

            points = np.empty((nx*ny, 2), dtype="float")
            points[:, 0] = X.flatten()
            points[:, 1] = Y.flatten()
            from anuga import inside_polygon
            indices = inside_polygon(points, polygon)
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
        ldx = dx/nx
        ldy = dy/ny

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
            catchment_area = data.size*ldx*ldy
        else:
            pmask = data.flat[indices]
            data_max_in_period = np.max(pmask)
            total_data_volume = np.sum(pmask)*ldx*ldy
            catchment_area = len(pmask)*ldx*ldy

        # print indices
        peak_intensity = data_max_in_period/time_step

        if print_stats or self.verbose:
            print('Time period = ', time_period)
            print('Time step = ', time_step)
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
        date_time = datetime.utcfromtimestamp(
            times[tid]).strftime("%Y/%m/%d %H:%M")

        if self.debug:
            print('--- Date/Time ', date_time)

        plt.figure(1)
        plt.clf()

        if not polygons is None:
            for polygon in polygons:
                if not polygon is None:
                    polygon = np.array(polygon)
                    plt.plot(polygon[:, 0], polygon[:, 1], '--w')

        plot_title = 'RASTER DATA '+date_time
        plt.suptitle(plot_title)
        s_title = 'Max Data = %.2e / period' % (np.max(data_slices[tid]))
        plt.title(s_title)

        plt.imshow(data_slices[tid], origin='upper', interpolation='bicubic',
                   extent=self.extent, vmin=0, vmax=plot_vmax)

        plt.colorbar()

        if show:
            plt.pause(0.001)
        if save:
            plt.savefig(plot_title+'.jpg', format='jpg')

        return

    def plot_accumulated_data(self, polygons=None):

        data_accumulated = self.data_accumulated

        time_step = self.time_step

        x = self.x
        y = self.y

        if self.debug:
            print(data_accumulated)
            print('maximum rain =', np.max(data_accumulated))
            print('mean rain =', np.mean(data_accumulated))

        dx = self.extent[1]-self.extent[0]
        dy = self.extent[3]-self.extent[2]
        # Volume in Million m3 over 128km x 128km area
        total_data_vol = np.mean(data_accumulated)*dx*dy/1e6
        data_max_in_period = self.data_max_in_period
        peak_intensity = data_max_in_period/time_step
        extent = self.extent

        if self.verbose:
            print('Total data volume =', total_data_vol)
            print('Peak data in 1 time step = ', data_max_in_period)
            print('Peak Intensity in 1 timestep =', peak_intensity)
            print('extent', extent)
            print('size', extent[1]-extent[0], extent[3]-extent[2])
            print(time_step)

        plt.figure(1)
        plt.clf()
        plt.suptitle('Accumulated Data')

        s_title = 'Max Int. = %.2e, Average = %.2e, Tot Vol. = %.2e Mill. m3' \
                  % (peak_intensity, np.mean(data_accumulated), total_data_vol)
        plt.title(s_title)

        if not polygons is None:
            for polygon in polygons:
                if not polygon is None:
                    polygon = np.array(polygon)
                    plt.plot(polygon[:, 0], polygon[:, 1], '--w')

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
            average_values = total_values/(self.times[-1]-self.times[0])
            max_intensity = max(bar_values)/time_step

            b_title = 'Total = %.2e, Average = %.2e, Max Int.= %.2e' % (
                total_values, average_values, max_intensity)

            plt.bar(t, bar_values, width=time_step,)
            #plt.suptitle(' Data for Location %s:' % lid, fontsize=12, fontweight='bold')
            plt.suptitle(' Data for Location %s:' % lid)
            plt.title(b_title)

            plt.xlabel('time (s)')
            plt.ylabel('Data')
            plt.show()

        return

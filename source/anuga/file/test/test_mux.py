import unittest
import tempfile
import numpy as num
import os
from struct import pack, unpack
from anuga.file.netcdf import NetCDFFile

from anuga.utilities.numerical_tools import ensure_numeric
from anuga.coordinate_transforms.redfearn import redfearn
from anuga.coordinate_transforms.geo_reference import Geo_reference

from anuga.file.mux import WAVEHEIGHT_MUX_LABEL, EAST_VELOCITY_LABEL, \
                            NORTH_VELOCITY_LABEL

from anuga.file.mux import WAVEHEIGHT_MUX2_LABEL, EAST_VELOCITY_MUX2_LABEL, \
                NORTH_VELOCITY_MUX2_LABEL
                
from anuga.file.mux import read_mux2_py
from anuga.file_conversion.urs2sts import urs2sts
from anuga.file.urs import Read_urs

class Test_Mux(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def write_mux(self, lat_long_points, time_step_count, time_step,
                  depth=None, ha=None, ua=None, va=None):
        """
        This will write 3 non-gridded mux files, for testing.
        If no quantities are passed in,
        na and va quantities will be the Easting values.
        Depth and ua will be the Northing value.
        
        The mux file format has south as positive so
        this function will swap the sign for va.  
        """

        #print "lat_long_points", lat_long_points
        #print "time_step_count",time_step_count
        #print "time_step",

        
        points_num = len(lat_long_points)
        lonlatdeps = []
        quantities = ['HA','UA','VA']
        
        mux_names = [WAVEHEIGHT_MUX_LABEL,
                     EAST_VELOCITY_LABEL,
                     NORTH_VELOCITY_LABEL]
        quantities_init = [[],[],[]]
        # urs binary is latitude fastest
        for point in lat_long_points:
            lat = point[0]
            lon = point[1]
            _ , e, n = redfearn(lat, lon)
            if depth is None:
                this_depth = n
            else:
                this_depth = depth
            if ha is None:
                this_ha = e
            else:
                this_ha = ha
            if ua is None:
                this_ua = n
            else:
                this_ua = ua
            if va is None:
                this_va = e   
            else:
                this_va = va         
            lonlatdeps.append([lon, lat, this_depth])
            quantities_init[0].append(this_ha) # HA
            quantities_init[1].append(this_ua) # UA
            quantities_init[2].append(this_va) # VA 
                
        file_handle, base_name = tempfile.mkstemp("")
        os.close(file_handle)
        os.remove(base_name)

        files = []        
        for i, q in enumerate(quantities): 
            quantities_init[i] = ensure_numeric(quantities_init[i])
            #print "HA_init", HA_init
            q_time = num.zeros((time_step_count, points_num), num.float)
            for time in range(time_step_count):
                q_time[time,:] = quantities_init[i] #* time * 4
            
            #Write C files
            columns = 3 # long, lat , depth
            file = base_name + mux_names[i]
            #print "base_name file",file 
            f = open(file, 'wb')
            files.append(file)
            f.write(pack('i',points_num))
            f.write(pack('i',time_step_count))
            f.write(pack('f',time_step))

            #write lat/long info
            for lonlatdep in lonlatdeps:
                for float in lonlatdep:
                    f.write(pack('f',float))
                    
            # Write quantity info
            for time in  range(time_step_count):
                for point_i in range(points_num):
                    f.write(pack('f',q_time[time,point_i]))
                    #print " mux_names[i]", mux_names[i] 
                    #print "f.write(pack('f',q_time[time,i]))", q_time[time,point_i]
            f.close()
        return base_name, files
        
    
    def delete_mux(self, files):
        for file in files:
            try:
                os.remove(file)
            except:
                pass
        
    def write_mux2(self, lat_long_points, time_step_count, time_step,
                   first_tstep, last_tstep,
                   depth=None, ha=None, ua=None, va=None):
        """
        This will write 3 non-gridded mux files, for testing.
        If no quantities are passed in,
        na and va quantities will be the Easting values.
        Depth and ua will be the Northing value.
        """
        #print "lat_long_points", lat_long_points
        #print "time_step_count",time_step_count
        #print "time_step",

        #irrelevant header information
        ig=ilon=ilat=0
        mcolat=mcolon=centerlat=centerlon=offset=az=baz=id=0.0

        points_num = len(lat_long_points)
        latlondeps = []
        quantities = ['HA','UA','VA']

        mux_names = [WAVEHEIGHT_MUX2_LABEL,
                     EAST_VELOCITY_MUX2_LABEL,
                     NORTH_VELOCITY_MUX2_LABEL]

        msg='first_tstep and last_step arrays must have same length as number of points'
        assert len(first_tstep)==points_num,msg
        assert len(last_tstep)==points_num,msg

        if depth is not None:
            depth=ensure_numeric(depth)
            assert len(depth)==points_num
        if ha is not None:
            ha=ensure_numeric(ha)
            assert ha.shape==(points_num,time_step_count)
        if ua is not None:
            ua=ensure_numeric(ua)
            assert ua.shape==(points_num,time_step_count)
        if va is not None:
            va=ensure_numeric(va)
            assert va.shape==(points_num,time_step_count)

        quantities_init = [[],[],[]]
        # urs binary is latitude fastest
        for i,point in enumerate(lat_long_points):
            lat = point[0]
            lon = point[1]
            _ , e, n = redfearn(lat, lon)
            if depth is None:
                this_depth = n
            else:
                this_depth = depth[i]
            latlondeps.append([lat, lon, this_depth])

            if ha is None:
                this_ha = e
                quantities_init[0].append(num.ones(time_step_count,num.float)*this_ha) # HA
            else:
                quantities_init[0].append(ha[i])
            if ua is None:
                this_ua = n
                quantities_init[1].append(num.ones(time_step_count,num.float)*this_ua) # UA
            else:
                quantities_init[1].append(ua[i])
            if va is None:
                this_va = e
                quantities_init[2].append(num.ones(time_step_count,num.float)*this_va) #
            else:
                quantities_init[2].append(-va[i]) # South is negative in MUX

        file_handle, base_name = tempfile.mkstemp("write_mux2")
        os.close(file_handle)
        os.remove(base_name)

        files = []        
        for i, q in enumerate(quantities):
            q_time = num.zeros((time_step_count, points_num), num.float)
            quantities_init[i] = ensure_numeric(quantities_init[i])
            for time in range(time_step_count):
                #print i, q, time, quantities_init[i][:,time]
                q_time[time,:] = quantities_init[i][:,time]
                #print i, q, time, q_time[time, :]                

            #Write C files
            columns = 3 # long, lat , depth
            file = base_name + mux_names[i]
            
            #print 'base_name file', file 
            f = open(file, 'wb')
            files.append(file)

            f.write(pack('i',points_num))
            #write mux 2 header
            for latlondep in latlondeps:
                f.write(pack('f',latlondep[0]))
                f.write(pack('f',latlondep[1]))
                f.write(pack('f',mcolat))
                f.write(pack('f',mcolon))
                f.write(pack('i',ig))
                f.write(pack('i',ilon))
                f.write(pack('i',ilat))
                f.write(pack('f',latlondep[2]))
                f.write(pack('f',centerlat))
                f.write(pack('f',centerlon))
                f.write(pack('f',offset))
                f.write(pack('f',az))
                f.write(pack('f',baz))
                f.write(pack('f',time_step))
                f.write(pack('i',time_step_count))
                for j in range(4): # identifier
                    f.write(pack('f',id))    

            #first_tstep=1
            #last_tstep=time_step_count
            for i,latlondep in enumerate(latlondeps):
                f.write(pack('i',first_tstep[i]))
            for i,latlondep in enumerate(latlondeps):
                f.write(pack('i',last_tstep[i]))

            # Find when first station starts recording
            min_tstep = min(first_tstep)
            # Find when all stations have stopped recording
            max_tstep = max(last_tstep)

            #for time in  range(time_step_count):
            for time in range(min_tstep-1,max_tstep):
                f.write(pack('f',time*time_step))                
                for point_i in range(points_num):
                    if time+1>=first_tstep[point_i] and time+1<=last_tstep[point_i]:
                        #print 'writing', time, point_i, q_time[time, point_i]
                        f.write(pack('f', q_time[time, point_i]))
            f.close()

        return base_name, files


    def test_urs2sts_read_mux2_pyI(self):
        """test_urs2sts_read_mux2_pyI(self):
        Constant stage,momentum at each gauge
        """
        tide = 1
        time_step_count = 3
        time_step = 2
        lat_long_points =[(-21.5,114.5),(-21,114.5),(-21.5,115), (-21.,115.)]
        n=len(lat_long_points)
        first_tstep=num.ones(n,num.int)
        last_tstep=time_step_count*num.ones(n,num.int)
        depth=20*num.ones(n,num.float)
        ha=2*num.ones((n,time_step_count),num.float)
        ua=5*num.ones((n,time_step_count),num.float)
        va=-10*num.ones((n,time_step_count),num.float)
        #-ve added to take into account mux file format where south is positive.
        base_name, files = self.write_mux2(lat_long_points,
                                      time_step_count, time_step,
                                      first_tstep, last_tstep,
                                      depth=depth,
                                      ha=ha,
                                      ua=ua,
                                      va=va)

        weights=num.ones(1, num.float)
        #ensure that files are indeed mux2 files
        times, latitudes, longitudes, elevation, stage, starttime = read_mux2_py([files[0]],weights)
        ua_times, ua_latitudes, ua_longitudes, ua_elevation, xvelocity,starttime_ua=read_mux2_py([files[1]],weights)
        msg='ha and ua have different gauge meta data'
        assert num.allclose(times,ua_times) and num.allclose(latitudes,ua_latitudes) and num.allclose(longitudes,ua_longitudes) and num.allclose(elevation,ua_elevation) and num.allclose(starttime,starttime_ua), msg
        va_times, va_latitudes, va_longitudes, va_elevation, yvelocity, starttime_va=read_mux2_py([files[2]],weights)
        msg='ha and va have different gauge meta data'
        assert num.allclose(times,va_times) and num.allclose(latitudes,va_latitudes) and num.allclose(longitudes,va_longitudes) and num.allclose(elevation,va_elevation) and num.allclose(starttime,starttime_va), msg

        self.delete_mux(files)

        msg='time array has incorrect length'
        assert times.shape[0]==time_step_count,msg
        
        msg = 'time array is incorrect'
        #assert allclose(times,time_step*num.arange(1,time_step_count+1)),msg
        assert num.allclose(times,time_step*num.arange(time_step_count)), msg
        
        msg='Incorrect gauge positions returned'
        for i,point in enumerate(lat_long_points):
            assert num.allclose(latitudes[i],point[0]) and num.allclose(longitudes[i],point[1]),msg

        msg='Incorrect gauge depths returned'
        assert num.allclose(elevation,-depth),msg
        msg='incorrect gauge height time series returned'
        assert num.allclose(stage,ha)
        msg='incorrect gauge ua time series returned'
        assert num.allclose(xvelocity,ua)
        msg='incorrect gauge va time series returned'
        assert num.allclose(yvelocity, -va)

    def test_urs2sts_read_mux2_pyII(self):
        """Spatially varing stage
        """
        tide = 1
        time_step_count = 3
        time_step = 2
        lat_long_points =[(-21.5,114.5),(-21,114.5),(-21.5,115), (-21.,115.)]
        n=len(lat_long_points)
        first_tstep=num.ones(n,num.int)
        last_tstep=(time_step_count)*num.ones(n,num.int)
        depth=20*num.ones(n,num.float)
        ha=2*num.ones((n,time_step_count),num.float)
        ha[0]=num.arange(0,time_step_count)+1
        ha[1]=time_step_count-num.arange(1,time_step_count+1)
        ha[1]=num.arange(time_step_count,2*time_step_count)
        ha[2]=num.arange(2*time_step_count,3*time_step_count)
        ha[3]=num.arange(3*time_step_count,4*time_step_count)
        ua=5*num.ones((n,time_step_count),num.float)
        va=-10*num.ones((n,time_step_count),num.float)
        #-ve added to take into account mux file format where south is positive.
        base_name, files = self.write_mux2(lat_long_points,
                                      time_step_count, time_step,
                                      first_tstep, last_tstep,
                                      depth=depth,
                                      ha=ha,
                                      ua=ua,
                                      va=va)

        weights=num.ones(1, num.float)
        #ensure that files are indeed mux2 files
        times, latitudes, longitudes, elevation, stage,starttime=read_mux2_py([files[0]],weights)
        ua_times, ua_latitudes, ua_longitudes, ua_elevation, xvelocity,starttime_ua=read_mux2_py([files[1]],weights)
        msg='ha and ua have different gauge meta data'
        assert num.allclose(times,ua_times) and num.allclose(latitudes,ua_latitudes) and num.allclose(longitudes,ua_longitudes) and num.allclose(elevation,ua_elevation) and num.allclose(starttime,starttime_ua), msg
        va_times, va_latitudes, va_longitudes, va_elevation, yvelocity,starttime_va=read_mux2_py([files[2]],weights)
        msg='ha and va have different gauge meta data'
        assert num.allclose(times,va_times) and num.allclose(latitudes,va_latitudes) and num.allclose(longitudes,va_longitudes) and num.allclose(elevation,va_elevation) and num.allclose(starttime,starttime_va), msg


        self.delete_mux(files)

        msg='time array has incorrect length'
        #assert times.shape[0]==time_step_count,msg
        msg = 'time array is incorrect'
        #assert allclose(times,time_step*num.arange(1,time_step_count+1)),msg
        msg='Incorrect gauge positions returned'
        for i,point in enumerate(lat_long_points):
            assert num.allclose(latitudes[i],point[0]) and num.allclose(longitudes[i],point[1]),msg

        msg='Incorrect gauge depths returned'
        assert num.allclose(elevation, -depth),msg
        msg='incorrect gauge height time series returned'
        assert num.allclose(stage, ha)
        msg='incorrect gauge ua time series returned'
        assert num.allclose(xvelocity, ua)
        msg='incorrect gauge va time series returned'
        assert num.allclose(yvelocity, -va) # South is positive in MUX

    def test_urs2sts_read_mux2_pyIII(self):
        """Varying start and finish times
        """
        tide = 1
        time_step_count = 3
        time_step = 2
        lat_long_points =[(-21.5,114.5),(-21,114.5),(-21.5,115), (-21.,115.)]
        n=len(lat_long_points)
        first_tstep=num.ones(n,num.int)
        first_tstep[0]+=1
        first_tstep[2]+=1
        last_tstep=(time_step_count)*num.ones(n,num.int)
        last_tstep[0]-=1

        depth=20*num.ones(n,num.float)
        ha=2*num.ones((n,time_step_count),num.float)
        ha[0]=num.arange(0,time_step_count)
        ha[1]=num.arange(time_step_count,2*time_step_count)
        ha[2]=num.arange(2*time_step_count,3*time_step_count)
        ha[3]=num.arange(3*time_step_count,4*time_step_count)
        ua=5*num.ones((n,time_step_count),num.float)
        va=-10*num.ones((n,time_step_count),num.float)
        #-ve added to take into account mux file format where south is positive.
        base_name, files = self.write_mux2(lat_long_points,
                                      time_step_count, time_step,
                                      first_tstep, last_tstep,
                                      depth=depth,
                                      ha=ha,
                                      ua=ua,
                                      va=va)

        weights=num.ones(1, num.float)
        #ensure that files are indeed mux2 files
        times, latitudes, longitudes, elevation, stage, starttime=read_mux2_py([files[0]],weights)
        ua_times, ua_latitudes, ua_longitudes, ua_elevation, xvelocity, starttime_ua=read_mux2_py([files[1]],weights)
        msg='ha and ua have different gauge meta data'
        assert num.allclose(times,ua_times) and num.allclose(latitudes,ua_latitudes) and num.allclose(longitudes,ua_longitudes) and num.allclose(elevation,ua_elevation) and num.allclose(starttime,starttime_ua), msg
        va_times, va_latitudes, va_longitudes, va_elevation, yvelocity,starttime_va=read_mux2_py([files[2]],weights)
        msg='ha and va have different gauge meta data'
        assert num.allclose(times,va_times) and num.allclose(latitudes,va_latitudes) and num.allclose(longitudes,va_longitudes) and num.allclose(elevation,va_elevation) and num.allclose(starttime,starttime_va), msg

        self.delete_mux(files)

        msg='time array has incorrect length'
        #assert times.shape[0]==time_step_count,msg
        msg = 'time array is incorrect'
        #assert allclose(times,time_step*num.arange(1,time_step_count+1)),msg
        msg='Incorrect gauge positions returned'
        for i,point in enumerate(lat_long_points):
            assert num.allclose(latitudes[i],point[0]) and num.allclose(longitudes[i],point[1]),msg


        # Set original data used to write mux file to be zero when gauges are 
        #not recdoring
        ha[0][0]=0.0
        ha[0][time_step_count-1]=0.0;
        ha[2][0]=0.0;
        ua[0][0]=0.0
        ua[0][time_step_count-1]=0.0;
        ua[2][0]=0.0;
        va[0][0]=0.0
        va[0][time_step_count-1]=0.0;
        va[2][0]=0.0;
        msg='Incorrect gauge depths returned'
        assert num.allclose(elevation,-depth),msg
        msg='incorrect gauge height time series returned'
        assert num.allclose(stage,ha)
        msg='incorrect gauge ua time series returned'
        assert num.allclose(xvelocity,ua)
        msg='incorrect gauge va time series returned'
        assert num.allclose(yvelocity, -va) # South is positive in mux
        

        
    def test_read_mux_platform_problem1(self):
        """test_read_mux_platform_problem1
        
        This is to test a situation where read_mux returned 
        wrong values Win32

        This test passes on Windows but test_read_mux_platform_problem2
        does not
        """
        
        from anuga.file.urs_ext import read_mux2 
        
        verbose = False
                
        tide = 1.5
        time_step_count = 10
        time_step = 0.2
        times_ref = num.arange(0, time_step_count*time_step, time_step)

        lat_long_points = [(-21.5,114.5), (-21,114.5), (-21.5,115), (-21.,115.), (-22., 117.)]
        n = len(lat_long_points)
        
        # Create different timeseries starting and ending at different times 
        first_tstep=num.ones(n, num.int)
        first_tstep[0]+=2   # Point 0 starts at 2
        first_tstep[1]+=4   # Point 1 starts at 4        
        first_tstep[2]+=3   # Point 2 starts at 3
        
        last_tstep=(time_step_count)*num.ones(n,num.int)
        last_tstep[0]-=1    # Point 0 ends 1 step early
        last_tstep[1]-=2    # Point 1 ends 2 steps early                
        last_tstep[4]-=3    # Point 4 ends 3 steps early        
        
        # Create varying elevation data (positive values for seafloor)
        gauge_depth=20*num.ones(n,num.float)
        for i in range(n):
            gauge_depth[i] += i**2
            
        # Create data to be written to first mux file        
        ha0=2*num.ones((n,time_step_count),num.float)
        ha0[0]=num.arange(0,time_step_count)
        ha0[1]=num.arange(time_step_count,2*time_step_count)
        ha0[2]=num.arange(2*time_step_count,3*time_step_count)
        ha0[3]=num.arange(3*time_step_count,4*time_step_count)
        ua0=5*num.ones((n,time_step_count),num.float)
        va0=-10*num.ones((n,time_step_count),num.float)

        # Ensure data used to write mux file to be zero when gauges are
        # not recording
        for i in range(n):
             # For each point
             for j in range(0, first_tstep[i]-1) + range(last_tstep[i], time_step_count):
                 # For timesteps before and after recording range
                 ha0[i][j] = ua0[i][j] = va0[i][j] = 0.0                                  
        
        # Write first mux file to be combined by urs2sts
        base_nameI, filesI = self.write_mux2(lat_long_points,
                                             time_step_count, time_step,
                                             first_tstep, last_tstep,
                                             depth=gauge_depth,
                                             ha=ha0,
                                             ua=ua0,
                                             va=va0)

        # Create ordering file
        permutation = ensure_numeric([4,0,2])

        _, ordering_filename = tempfile.mkstemp('')
        order_fid = open(ordering_filename, 'w')  
        order_fid.write('index, longitude, latitude\n')
        for index in permutation:
            order_fid.write('%d, %f, %f\n' %(index, 
                                             lat_long_points[index][1], 
                                             lat_long_points[index][0]))
        order_fid.close()
        
        

        # -------------------------------------
        # Now read files back and check values
        weights = ensure_numeric([1.0])

        # For each quantity read the associated list of source mux2 file with 
        # extention associated with that quantity
        file_params=-1*num.ones(3,num.float) #[nsta,dt,nt]
        OFFSET = 5

        for j, file in enumerate(filesI):
            data = read_mux2(1, [file], weights, file_params, permutation, verbose)

            number_of_selected_stations = data.shape[0]

            # Index where data ends and parameters begin
            parameters_index = data.shape[1]-OFFSET          
          
            for i in range(number_of_selected_stations):
                if j == 0: assert num.allclose(data[i][:parameters_index], ha0[permutation[i], :])
                if j == 1: assert num.allclose(data[i][:parameters_index], ua0[permutation[i], :])
                if j == 2: assert num.allclose(data[i][:parameters_index], -va0[permutation[i], :])
        
        self.delete_mux(filesI)


        
        
    def test_read_mux_platform_problem2(self):
        """test_read_mux_platform_problem2
        
        This is to test a situation where read_mux returned 
        wrong values Win32

        This test does not pass on Windows but test_read_mux_platform_problem1
        does
        """
        
        from anuga.file.urs_ext import read_mux2 
        
        from anuga.config import single_precision as epsilon        
        
        verbose = False
                
        tide = 1.5
        time_step_count = 10
        time_step = 0.2
        
        times_ref = num.arange(0, time_step_count*time_step, time_step)
        
        lat_long_points = [(-21.5,114.5), (-21,114.5), (-21.5,115),
                           (-21.,115.), (-22., 117.)]
        n = len(lat_long_points)
        
        # Create different timeseries starting and ending at different times 
        first_tstep=num.ones(n,num.int)
        first_tstep[0]+=2   # Point 0 starts at 2
        first_tstep[1]+=4   # Point 1 starts at 4        
        first_tstep[2]+=3   # Point 2 starts at 3
        
        last_tstep=(time_step_count)*num.ones(n,num.int)
        last_tstep[0]-=1    # Point 0 ends 1 step early
        last_tstep[1]-=2    # Point 1 ends 2 steps early                
        last_tstep[4]-=3    # Point 4 ends 3 steps early        
        
        # Create varying elevation data (positive values for seafloor)
        gauge_depth=20*num.ones(n,num.float)
        for i in range(n):
            gauge_depth[i] += i**2
            
        # Create data to be written to second mux file        
        ha1=num.ones((n,time_step_count),num.float)
        ha1[0]=num.sin(times_ref)
        ha1[1]=2*num.sin(times_ref - 3)
        ha1[2]=5*num.sin(4*times_ref)
        ha1[3]=num.sin(times_ref)
        ha1[4]=num.sin(2*times_ref-0.7)
                
        ua1=num.zeros((n,time_step_count),num.float)
        ua1[0]=3*num.cos(times_ref)        
        ua1[1]=2*num.sin(times_ref-0.7)   
        ua1[2]=num.arange(3*time_step_count,4*time_step_count)
        ua1[4]=2*num.ones(time_step_count)
        
        va1=num.zeros((n,time_step_count),num.float)
        va1[0]=2*num.cos(times_ref-0.87)        
        va1[1]=3*num.ones(time_step_count)
        va1[3]=2*num.sin(times_ref-0.71)        
        
        # Ensure data used to write mux file to be zero when gauges are
        # not recording
        for i in range(n):
             # For each point
             for j in range(0, first_tstep[i]-1) + range(last_tstep[i], time_step_count):
                 # For timesteps before and after recording range
                 ha1[i][j] = ua1[i][j] = va1[i][j] = 0.0 


        #print 'Second station to be written to MUX'
        #print 'ha', ha1[0,:]
        #print 'ua', ua1[0,:]
        #print 'va', va1[0,:]
        
        # Write second mux file to be combined by urs2sts 
        base_nameII, filesII = self.write_mux2(lat_long_points,
                                               time_step_count, time_step,
                                               first_tstep, last_tstep,
                                               depth=gauge_depth,
                                               ha=ha1,
                                               ua=ua1,
                                               va=va1)




        # Read mux file back and verify it's correcness

        ####################################################
        # FIXME (Ole): This is where the test should
        # verify that the MUX files are correct.

        #JJ: It appears as though
        #that certain quantities are not being stored with enough precision
        #inn muxfile or more likely that they are being cast into a
        #lower precision when read in using read_mux2 Time step and q_time
        # are equal but only to approx 1e-7
        ####################################################

        #define information as it should be stored in mus2 files
        points_num=len(lat_long_points)
        depth=gauge_depth
        ha=ha1
        ua=ua1
        va=va1
        
        quantities = ['HA','UA','VA']
        mux_names = [WAVEHEIGHT_MUX2_LABEL,
                     EAST_VELOCITY_MUX2_LABEL,
                     NORTH_VELOCITY_MUX2_LABEL]
        quantities_init = [[],[],[]]
        latlondeps = []
        #irrelevant header information
        ig=ilon=ilat=0
        mcolat=mcolon=centerlat=centerlon=offset=az=baz=id=0.0
        # urs binary is latitude fastest
        for i,point in enumerate(lat_long_points):
            lat = point[0]
            lon = point[1]
            _ , e, n = redfearn(lat, lon)
            if depth is None:
                this_depth = n
            else:
                this_depth = depth[i]
            latlondeps.append([lat, lon, this_depth])

            if ha is None:
                this_ha = e
                quantities_init[0].append(num.ones(time_step_count,num.float)*this_ha) # HA
            else:
                quantities_init[0].append(ha[i])
            if ua is None:
                this_ua = n
                quantities_init[1].append(num.ones(time_step_count,num.float)*this_ua) # UA
            else:
                quantities_init[1].append(ua[i])
            if va is None:
                this_va = e
                quantities_init[2].append(num.ones(time_step_count,num.float)*this_va) #
            else:
                quantities_init[2].append(va[i])

        for i, q in enumerate(quantities):
            #print
            #print i, q
            
            q_time = num.zeros((time_step_count, points_num), num.float)
            quantities_init[i] = ensure_numeric(quantities_init[i])
            for time in range(time_step_count):
                #print i, q, time, quantities_init[i][:,time]
                q_time[time,:] = quantities_init[i][:,time]
                #print i, q, time, q_time[time, :]

            
            filename = base_nameII + mux_names[i]
            f = open(filename, 'rb')
            assert abs(points_num-unpack('i',f.read(4))[0])<epsilon
            #write mux 2 header
            for latlondep in latlondeps:
                assert abs(latlondep[0]-unpack('f',f.read(4))[0])<epsilon
                assert abs(latlondep[1]-unpack('f',f.read(4))[0])<epsilon
                assert abs(mcolat-unpack('f',f.read(4))[0])<epsilon
                assert abs(mcolon-unpack('f',f.read(4))[0])<epsilon
                assert abs(ig-unpack('i',f.read(4))[0])<epsilon
                assert abs(ilon-unpack('i',f.read(4))[0])<epsilon
                assert abs(ilat-unpack('i',f.read(4))[0])<epsilon
                assert abs(latlondep[2]-unpack('f',f.read(4))[0])<epsilon
                assert abs(centerlat-unpack('f',f.read(4))[0])<epsilon
                assert abs(centerlon-unpack('f',f.read(4))[0])<epsilon
                assert abs(offset-unpack('f',f.read(4))[0])<epsilon
                assert abs(az-unpack('f',f.read(4))[0])<epsilon
                assert abs(baz-unpack('f',f.read(4))[0])<epsilon
                
                x = unpack('f', f.read(4))[0]
                #print time_step
                #print x
                assert abs(time_step-x)<epsilon
                assert abs(time_step_count-unpack('i',f.read(4))[0])<epsilon
                for j in range(4): # identifier
                    assert abs(id-unpack('i',f.read(4))[0])<epsilon 

            #first_tstep=1
            #last_tstep=time_step_count
            for i,latlondep in enumerate(latlondeps):
                assert abs(first_tstep[i]-unpack('i',f.read(4))[0])<epsilon
            for i,latlondep in enumerate(latlondeps):
                assert abs(last_tstep[i]-unpack('i',f.read(4))[0])<epsilon

            # Find when first station starts recording
            min_tstep = min(first_tstep)
            # Find when all stations have stopped recording
            max_tstep = max(last_tstep)

            #for time in  range(time_step_count):
            for time in range(min_tstep-1,max_tstep):
                assert abs(time*time_step-unpack('f',f.read(4))[0])<epsilon
                for point_i in range(points_num):
                    if time+1>=first_tstep[point_i] and time+1<=last_tstep[point_i]:
                        x = unpack('f',f.read(4))[0]
                        #print time, x, q_time[time, point_i]
                        if q == 'VA': x = -x # South is positive in MUX
                        assert abs(q_time[time, point_i]-x)<epsilon

            f.close()
                                               
        # Create ordering file
        permutation = ensure_numeric([4,0,2])

       #  _, ordering_filename = tempfile.mkstemp('')
#         order_fid = open(ordering_filename, 'w')  
#         order_fid.write('index, longitude, latitude\n')
#         for index in permutation:
#             order_fid.write('%d, %f, %f\n' %(index, 
#                                              lat_long_points[index][1], 
#                                              lat_long_points[index][0]))
#         order_fid.close()
        
        # -------------------------------------
        # Now read files back and check values
        weights = ensure_numeric([1.0])

        # For each quantity read the associated list of source mux2 file with 
        # extention associated with that quantity
        file_params=-1*num.ones(3,num.float) # [nsta,dt,nt]
        OFFSET = 5

        for j, file in enumerate(filesII):
            # Read stage, u, v enumerated as j
            #print 'Reading', j, file
            data = read_mux2(1, [file], weights, file_params, permutation, verbose)

            #print 'Data received by Python'
            #print data[1][8]
            number_of_selected_stations = data.shape[0]

            # Index where data ends and parameters begin
            parameters_index = data.shape[1]-OFFSET          
                 
            quantity=num.zeros((number_of_selected_stations, parameters_index), num.float)
            
            
            for i in range(number_of_selected_stations):
        
                #print i, parameters_index
                #print quantity[i][:]
                if j == 0: assert num.allclose(data[i][:parameters_index], ha1[permutation[i], :])
                if j == 1: assert num.allclose(data[i][:parameters_index], ua1[permutation[i], :])
                if j == 2:
                    # FIXME (Ole): This is where the output is wrong on Win32
                    
                    #print
                    #print j, i
                    #print 'Input'
                    #print 'u', ua1[permutation[i], 8]       
                    #print 'v', va1[permutation[i], 8]
                
                    #print 'Output'
                    #print 'v ', data[i][:parameters_index][8]  

                    # South is positive in MUX
                    #print "data[i][:parameters_index]", data[i][:parameters_index]
                    #print "-va1[permutation[i], :]", -va1[permutation[i], :]
                    assert num.allclose(data[i][:parameters_index], -va1[permutation[i], :])
        
        self.delete_mux(filesII)
           
    def test_read_mux_platform_problem3(self):
        
        # This is to test a situation where read_mux returned 
        # wrong values Win32

        
        from anuga.file.urs_ext import read_mux2 
        
        from anuga.config import single_precision as epsilon        
        
        verbose = False
                
        tide = 1.5
        time_step_count = 10
        time_step = 0.02

        '''
        Win results
        time_step = 0.2000001
        This is OK        
        '''
        
        '''
        Win results
        time_step = 0.20000001

        ======================================================================
ERROR: test_read_mux_platform_problem3 (__main__.Test_Data_Manager)
----------------------------------------------------------------------
Traceback (most recent call last):
  File "test_data_manager.py", line 6718, in test_read_mux_platform_problem3
    ha1[0]=num.sin(times_ref)
ValueError: matrices are not aligned for copy

        '''

        '''
        Win results
        time_step = 0.200000001
        FAIL
         assert num.allclose(data[i][:parameters_index],
         -va1[permutation[i], :])
        '''
        times_ref = num.arange(0, time_step_count*time_step, time_step)
        #print "times_ref", times_ref
        
        lat_long_points = [(-21.5,114.5), (-21,114.5), (-21.5,115),
                           (-21.,115.), (-22., 117.)]
        stations = len(lat_long_points)
        
        # Create different timeseries starting and ending at different times 
        first_tstep=num.ones(stations, num.int)
        first_tstep[0]+=2   # Point 0 starts at 2
        first_tstep[1]+=4   # Point 1 starts at 4        
        first_tstep[2]+=3   # Point 2 starts at 3
        
        last_tstep=(time_step_count)*num.ones(stations, num.int)
        last_tstep[0]-=1    # Point 0 ends 1 step early
        last_tstep[1]-=2    # Point 1 ends 2 steps early                
        last_tstep[4]-=3    # Point 4 ends 3 steps early        
        
        # Create varying elevation data (positive values for seafloor)
        gauge_depth=20*num.ones(stations, num.float)
        for i in range(stations):
            gauge_depth[i] += i**2
            
        # Create data to be written to second mux file        
        ha1=num.ones((stations,time_step_count), num.float)
        ha1[0]=num.sin(times_ref)
        ha1[1]=2*num.sin(times_ref - 3)
        ha1[2]=5*num.sin(4*times_ref)
        ha1[3]=num.sin(times_ref)
        ha1[4]=num.sin(2*times_ref-0.7)
                
        ua1=num.zeros((stations,time_step_count),num.float)
        ua1[0]=3*num.cos(times_ref)        
        ua1[1]=2*num.sin(times_ref-0.7)   
        ua1[2]=num.arange(3*time_step_count,4*time_step_count)
        ua1[4]=2*num.ones(time_step_count)
        
        va1=num.zeros((stations,time_step_count),num.float)
        va1[0]=2*num.cos(times_ref-0.87)        
        va1[1]=3*num.ones(time_step_count)
        va1[3]=2*num.sin(times_ref-0.71)        
        #print "va1[0]", va1[0]  # The 8th element is what will go bad.
        # Ensure data used to write mux file to be zero when gauges are
        # not recording
        for i in range(stations):
             # For each point
             for j in range(0, first_tstep[i]-1) + range(last_tstep[i],
                                                         time_step_count):
                 # For timesteps before and after recording range
                 ha1[i][j] = ua1[i][j] = va1[i][j] = 0.0 


        #print 'Second station to be written to MUX'
        #print 'ha', ha1[0,:]
        #print 'ua', ua1[0,:]
        #print 'va', va1[0,:]
        
        # Write second mux file to be combined by urs2sts 
        base_nameII, filesII = self.write_mux2(lat_long_points,
                                               time_step_count, time_step,
                                               first_tstep, last_tstep,
                                               depth=gauge_depth,
                                               ha=ha1,
                                               ua=ua1,
                                               va=va1)
        #print "filesII", filesII




        # Read mux file back and verify it's correcness

        ####################################################
        # FIXME (Ole): This is where the test should
        # verify that the MUX files are correct.

        #JJ: It appears as though
        #that certain quantities are not being stored with enough precision
        #inn muxfile or more likely that they are being cast into a
        #lower precision when read in using read_mux2 Time step and q_time
        # are equal but only to approx 1e-7
        ####################################################

        #define information as it should be stored in mus2 files
        points_num=len(lat_long_points)
        depth=gauge_depth
        ha=ha1
        ua=ua1
        va=va1
        
        quantities = ['HA','UA','VA']
        mux_names = [WAVEHEIGHT_MUX2_LABEL,
                     EAST_VELOCITY_MUX2_LABEL,
                     NORTH_VELOCITY_MUX2_LABEL]
        quantities_init = [[],[],[]]
        latlondeps = []
        #irrelevant header information
        ig=ilon=ilat=0
        mcolat=mcolon=centerlat=centerlon=offset=az=baz=id=0.0
        # urs binary is latitude fastest
        for i,point in enumerate(lat_long_points):
            lat = point[0]
            lon = point[1]
            _ , e, n = redfearn(lat, lon)
            if depth is None:
                this_depth = n
            else:
                this_depth = depth[i]
            latlondeps.append([lat, lon, this_depth])

            if ha is None:
                this_ha = e
                quantities_init[0].append(num.ones(time_step_count,
                                                   num.float)*this_ha) # HA
            else:
                quantities_init[0].append(ha[i])
            if ua is None:
                this_ua = n
                quantities_init[1].append(num.ones(time_step_count,
                                                   num.float)*this_ua) # UA
            else:
                quantities_init[1].append(ua[i])
            if va is None:
                this_va = e
                quantities_init[2].append(num.ones(time_step_count,
                                                   num.float)*this_va) #
            else:
                quantities_init[2].append(va[i])

        for i, q in enumerate(quantities):
            #print
            #print i, q
            
            q_time = num.zeros((time_step_count, points_num), num.float)
            quantities_init[i] = ensure_numeric(quantities_init[i])
            for time in range(time_step_count):
                #print i, q, time, quantities_init[i][:,time]
                q_time[time,:] = quantities_init[i][:,time]
                #print i, q, time, q_time[time, :]

            
            filename = base_nameII + mux_names[i]
            f = open(filename, 'rb')
            assert abs(points_num-unpack('i',f.read(4))[0])<epsilon
            #write mux 2 header
            for latlondep in latlondeps:
                assert abs(latlondep[0]-unpack('f',f.read(4))[0])<epsilon
                assert abs(latlondep[1]-unpack('f',f.read(4))[0])<epsilon
                assert abs(mcolat-unpack('f',f.read(4))[0])<epsilon
                assert abs(mcolon-unpack('f',f.read(4))[0])<epsilon
                assert abs(ig-unpack('i',f.read(4))[0])<epsilon
                assert abs(ilon-unpack('i',f.read(4))[0])<epsilon
                assert abs(ilat-unpack('i',f.read(4))[0])<epsilon
                assert abs(latlondep[2]-unpack('f',f.read(4))[0])<epsilon
                assert abs(centerlat-unpack('f',f.read(4))[0])<epsilon
                assert abs(centerlon-unpack('f',f.read(4))[0])<epsilon
                assert abs(offset-unpack('f',f.read(4))[0])<epsilon
                assert abs(az-unpack('f',f.read(4))[0])<epsilon
                assert abs(baz-unpack('f',f.read(4))[0])<epsilon
                
                x = unpack('f', f.read(4))[0]
                #print time_step
                #print x
                assert abs(time_step-x)<epsilon
                assert abs(time_step_count-unpack('i',f.read(4))[0])<epsilon
                for j in range(4): # identifier
                    assert abs(id-unpack('i',f.read(4))[0])<epsilon 

            #first_tstep=1
            #last_tstep=time_step_count
            for i,latlondep in enumerate(latlondeps):
                assert abs(first_tstep[i]-unpack('i',f.read(4))[0])<epsilon
            for i,latlondep in enumerate(latlondeps):
                assert abs(last_tstep[i]-unpack('i',f.read(4))[0])<epsilon

            # Find when first station starts recording
            min_tstep = min(first_tstep)
            # Find when all stations have stopped recording
            max_tstep = max(last_tstep)

            #for time in  range(time_step_count):
            for time in range(min_tstep-1,max_tstep):
                assert abs(time*time_step-unpack('f',f.read(4))[0])<epsilon
                for point_i in range(points_num):
                    if time+1>=first_tstep[point_i] and time+1<=last_tstep[point_i]:
                        x = unpack('f',f.read(4))[0]
                        #print time, x, q_time[time, point_i]
                        if q == 'VA': x = -x # South is positive in MUX
                        #print q+" q_time[%d, %d] = %f" %(time, point_i, 
                                                      #q_time[time, point_i])
                        assert abs(q_time[time, point_i]-x)<epsilon

            f.close()
                            
        permutation = ensure_numeric([4,0,2])
                   
        # Create ordering file
#         _, ordering_filename = tempfile.mkstemp('')
#         order_fid = open(ordering_filename, 'w')  
#         order_fid.write('index, longitude, latitude\n')
#         for index in permutation:
#             order_fid.write('%d, %f, %f\n' %(index, 
#                                              lat_long_points[index][1], 
#                                              lat_long_points[index][0]))
#         order_fid.close()
        
        # -------------------------------------
        # Now read files back and check values
        weights = ensure_numeric([1.0])

        # For each quantity read the associated list of source mux2 file with 
        # extention associated with that quantity
        file_params=-1*num.ones(3,num.float) # [nsta,dt,nt]
        OFFSET = 5

        for j, file in enumerate(filesII):
            # Read stage, u, v enumerated as j
            #print 'Reading', j, file
            #print "file", file
            data = read_mux2(1, [file], weights, file_params,
                             permutation, verbose)
            #print str(j) + "data", data

            #print 'Data received by Python'
            #print data[1][8]
            number_of_selected_stations = data.shape[0]
            #print "number_of_selected_stations", number_of_selected_stations
            #print "stations", stations

            # Index where data ends and parameters begin
            parameters_index = data.shape[1]-OFFSET          
                 
            for i in range(number_of_selected_stations):
        
                #print i, parameters_index
                if j == 0:
                    assert num.allclose(data[i][:parameters_index],
                                        ha1[permutation[i], :])
                    
                if j == 1: assert num.allclose(data[i][:parameters_index], ua1[permutation[i], :])
                if j == 2:
                    assert num.allclose(data[i][:parameters_index], -va1[permutation[i], :])
        
        self.delete_mux(filesII)      


          
    def test_urs2sts_nonstandard_projection_reverse(self):
        """
        Test that a point not in the specified zone can occur first
        """
        tide=0
        time_step_count = 3
        time_step = 2
        lat_long_points =[(-21.,113.5),(-21.,114.5),(-21.,114.), (-21.,115.)]
        n=len(lat_long_points)
        first_tstep=num.ones(n,num.int)
        first_tstep[0]+=1
        first_tstep[2]+=1
        last_tstep=(time_step_count)*num.ones(n,num.int)
        last_tstep[0]-=1

        gauge_depth=20*num.ones(n,num.float)
        ha=2*num.ones((n,time_step_count),num.float)
        ha[0]=num.arange(0,time_step_count)
        ha[1]=num.arange(time_step_count,2*time_step_count)
        ha[2]=num.arange(2*time_step_count,3*time_step_count)
        ha[3]=num.arange(3*time_step_count,4*time_step_count)
        ua=5*num.ones((n,time_step_count),num.float)
        va=-10*num.ones((n,time_step_count),num.float)

        base_name, files = self.write_mux2(lat_long_points,
                                      time_step_count, time_step,
                                      first_tstep, last_tstep,
                                      depth=gauge_depth,
                                      ha=ha,
                                      ua=ua,
                                      va=va)

        urs2sts(base_name,
                basename_out=base_name, 
                zone=50,
                mean_stage=tide,verbose=False)

        # now I want to check the sts file ...
        sts_file = base_name + '.sts'

        #Let's interigate the sww file
        # Note, the sww info is not gridded.  It is point data.
        fid = NetCDFFile(sts_file)

        # Make x and y absolute
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]

        geo_reference = Geo_reference(NetCDFObject=fid)
        points = geo_reference.get_absolute(map(None, x, y))
        points = ensure_numeric(points)

        x = points[:,0]
        y = points[:,1]

        # Check that all coordinate are correctly represented       
        # Using the non standard projection (50) 
        for i in range(4):
            zone, e, n = redfearn(lat_long_points[i][0], lat_long_points[i][1],
                                  zone=50) 
            assert num.allclose([x[i],y[i]], [e,n])
            assert zone==geo_reference.zone
        
        self.delete_mux(files)

            
    def test_urs2stsII(self):
        """
        Test multiple sources
        """
        tide=0
        time_step_count = 3
        time_step = 2
        lat_long_points =[(-21.5,114.5),(-21,114.5),(-21.5,115), (-21.,115.)]
        n=len(lat_long_points)
        first_tstep=num.ones(n,num.int)
        first_tstep[0]+=1
        first_tstep[2]+=1
        last_tstep=(time_step_count)*num.ones(n,num.int)
        last_tstep[0]-=1

        gauge_depth=20*num.ones(n,num.float)
        ha=2*num.ones((n,time_step_count),num.float)
        ha[0]=num.arange(0,time_step_count)
        ha[1]=num.arange(time_step_count,2*time_step_count)
        ha[2]=num.arange(2*time_step_count,3*time_step_count)
        ha[3]=num.arange(3*time_step_count,4*time_step_count)
        ua=5*num.ones((n,time_step_count),num.float)
        va=-10*num.ones((n,time_step_count),num.float)

        # Create two identical mux files to be combined by urs2sts
        base_nameI, filesI = self.write_mux2(lat_long_points,
                                             time_step_count, time_step,
                                             first_tstep, last_tstep,
                                             depth=gauge_depth,
                                             ha=ha,
                                             ua=ua,
                                             va=va)

        base_nameII, filesII = self.write_mux2(lat_long_points,
                                               time_step_count, time_step,
                                               first_tstep, last_tstep,
                                               depth=gauge_depth,
                                               ha=ha,
                                               ua=ua,
                                               va=va)

        # Call urs2sts with multiple mux files
        urs2sts([base_nameI, base_nameII], 
                basename_out=base_nameI, 
                weights=[1.0, 1.0],
                mean_stage=tide,
                verbose=False)

        # now I want to check the sts file ...
        sts_file = base_nameI + '.sts'

        #Let's interrogate the sts file
        # Note, the sts info is not gridded.  It is point data.
        fid = NetCDFFile(sts_file)

        # Make x and y absolute
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]

        geo_reference = Geo_reference(NetCDFObject=fid)
        points = geo_reference.get_absolute(map(None, x, y))
        points = ensure_numeric(points)

        x = points[:,0]
        y = points[:,1]

        #Check that first coordinate is correctly represented       
        #Work out the UTM coordinates for first point
        zone, e, n = redfearn(lat_long_points[0][0], lat_long_points[0][1]) 
        assert num.allclose([x[0],y[0]], [e,n])

        #Check the time vector
        times = fid.variables['time'][:]

        times_actual = []
        for i in range(time_step_count):
            times_actual.append(time_step * i)

        assert num.allclose(ensure_numeric(times),
                            ensure_numeric(times_actual))

        #Check first value
        stage = fid.variables['stage'][:]
        xmomentum = fid.variables['xmomentum'][:]
        ymomentum = fid.variables['ymomentum'][:]
        elevation = fid.variables['elevation'][:]

        # Set original data used to write mux file to be zero when gauges are
        # not recdoring
        
        ha[0][0]=0.0
        ha[0][time_step_count-1]=0.0
        ha[2][0]=0.0
        ua[0][0]=0.0
        ua[0][time_step_count-1]=0.0
        ua[2][0]=0.0
        va[0][0]=0.0
        va[0][time_step_count-1]=0.0
        va[2][0]=0.0;

        # The stage stored in the .sts file should be the sum of the stage
        # in the two mux2 files because both have weights = 1. In this case
        # the mux2 files are the same so stage == 2.0 * ha
        #print 2.0*num.transpose(ha) - stage 
        assert num.allclose(2.0*num.transpose(ha), stage)  #Meters

        #Check the momentums - ua
        #momentum = velocity*(stage-elevation)
        # elevation = - depth
        #momentum = velocity_ua *(stage+depth)

        depth=num.zeros((len(lat_long_points),time_step_count),num.float)
        for i in range(len(lat_long_points)):
            depth[i]=gauge_depth[i]+tide+2.0*ha[i]
            #2.0*ha necessary because using two files with weights=1 are used

        # The xmomentum stored in the .sts file should be the sum of the ua
        # in the two mux2 files multiplied by the depth.
        assert num.allclose(2.0*num.transpose(ua*depth), xmomentum) 

        #Check the momentums - va
        #momentum = velocity*(stage-elevation)
        # elevation = - depth
        #momentum = velocity_va *(stage+depth)

        # The ymomentum stored in the .sts file should be the sum of the va
        # in the two mux2 files multiplied by the depth.
        assert num.allclose(2.0*num.transpose(va*depth), ymomentum)

        # check the elevation values.
        # -ve since urs measures depth, sww meshers height,
        assert num.allclose(-elevation, gauge_depth)  #Meters

        fid.close()
        self.delete_mux(filesI)
        self.delete_mux(filesII)
        os.remove(sts_file)        


    def test_urs2sts0(self):
        """
        Test single source
        """
        tide=0
        time_step_count = 3
        time_step = 2
        lat_long_points =[(-21.5,114.5),(-21,114.5),(-21.5,115), (-21.,115.)]
        n=len(lat_long_points)
        first_tstep=num.ones(n,num.int)
        first_tstep[0]+=1
        first_tstep[2]+=1
        last_tstep=(time_step_count)*num.ones(n,num.int)
        last_tstep[0]-=1

        gauge_depth=20*num.ones(n,num.float)
        ha=2*num.ones((n,time_step_count),num.float)
        ha[0]=num.arange(0,time_step_count)
        ha[1]=num.arange(time_step_count,2*time_step_count)
        ha[2]=num.arange(2*time_step_count,3*time_step_count)
        ha[3]=num.arange(3*time_step_count,4*time_step_count)
        ua=5*num.ones((n,time_step_count),num.float)
        va=-10*num.ones((n,time_step_count),num.float)

        base_name, files = self.write_mux2(lat_long_points,
                                      time_step_count, time_step,
                                      first_tstep, last_tstep,
                                      depth=gauge_depth,
                                      ha=ha,
                                      ua=ua,
                                      va=va)

        urs2sts(base_name,
                basename_out=base_name, 
                mean_stage=tide,verbose=False)

        # now I want to check the sts file ...
        sts_file = base_name + '.sts'

        #Let's interigate the sww file
        # Note, the sww info is not gridded.  It is point data.
        fid = NetCDFFile(sts_file)

        # Make x and y absolute
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]

        geo_reference = Geo_reference(NetCDFObject=fid)
        points = geo_reference.get_absolute(map(None, x, y))
        points = ensure_numeric(points)

        x = points[:,0]
        y = points[:,1]

        #Check that first coordinate is correctly represented       
        #Work out the UTM coordinates for first point
        for i in range(4):
            zone, e, n = redfearn(lat_long_points[i][0], lat_long_points[i][1]) 
            assert num.allclose([x[i],y[i]], [e,n])

        #Check the time vector
        times = fid.variables['time'][:]

        times_actual = []
        for i in range(time_step_count):
            times_actual.append(time_step * i)

        assert num.allclose(ensure_numeric(times),
                            ensure_numeric(times_actual))

        #Check first value
        stage = fid.variables['stage'][:]
        xmomentum = fid.variables['xmomentum'][:]
        ymomentum = fid.variables['ymomentum'][:]
        elevation = fid.variables['elevation'][:]

        # Set original data used to write mux file to be zero when gauges are
        #not recdoring
        ha[0][0]=0.0
        ha[0][time_step_count-1]=0.0;
        ha[2][0]=0.0;
        ua[0][0]=0.0
        ua[0][time_step_count-1]=0.0;
        ua[2][0]=0.0;
        va[0][0]=0.0
        va[0][time_step_count-1]=0.0;
        va[2][0]=0.0;

        assert num.allclose(num.transpose(ha),stage)  #Meters

        #Check the momentums - ua
        #momentum = velocity*(stage-elevation)
        # elevation = - depth
        #momentum = velocity_ua *(stage+depth)

        depth=num.zeros((len(lat_long_points),time_step_count),num.float)
        for i in range(len(lat_long_points)):
            depth[i]=gauge_depth[i]+tide+ha[i]
        assert num.allclose(num.transpose(ua*depth),xmomentum) 

        #Check the momentums - va
        #momentum = velocity*(stage-elevation)
        # elevation = - depth
        #momentum = velocity_va *(stage+depth)

        assert num.allclose(num.transpose(va*depth),ymomentum)

        # check the elevation values.
        # -ve since urs measures depth, sww meshers height,
        assert num.allclose(-elevation, gauge_depth)  #Meters

        fid.close()
        self.delete_mux(files)
        os.remove(sts_file)

    def test_urs2sts_nonstandard_meridian(self):
        """
        Test single source using the meridian from zone 50 as a nonstandard meridian
        """
        tide=0
        time_step_count = 3
        time_step = 2
        lat_long_points =[(-21.,114.5),(-21.,113.5),(-21.,114.), (-21.,115.)]
        n=len(lat_long_points)
        first_tstep=num.ones(n,num.int)
        first_tstep[0]+=1
        first_tstep[2]+=1
        last_tstep=(time_step_count)*num.ones(n,num.int)
        last_tstep[0]-=1

        gauge_depth=20*num.ones(n,num.float)
        ha=2*num.ones((n,time_step_count),num.float)
        ha[0]=num.arange(0,time_step_count)
        ha[1]=num.arange(time_step_count,2*time_step_count)
        ha[2]=num.arange(2*time_step_count,3*time_step_count)
        ha[3]=num.arange(3*time_step_count,4*time_step_count)
        ua=5*num.ones((n,time_step_count),num.float)
        va=-10*num.ones((n,time_step_count),num.float)

        base_name, files = self.write_mux2(lat_long_points,
                                           time_step_count, time_step,
                                           first_tstep, last_tstep,
                                           depth=gauge_depth,
                                           ha=ha,
                                           ua=ua,
                                           va=va)

        urs2sts(base_name,
                basename_out=base_name, 
                central_meridian=123,
                mean_stage=tide,
                verbose=False)

        # now I want to check the sts file ...
        sts_file = base_name + '.sts'

        #Let's interigate the sww file
        # Note, the sww info is not gridded.  It is point data.
        fid = NetCDFFile(sts_file)

        # Make x and y absolute
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]

        geo_reference = Geo_reference(NetCDFObject=fid)
        points = geo_reference.get_absolute(map(None, x, y))
        points = ensure_numeric(points)

        x = points[:,0]
        y = points[:,1]

        # Check that all coordinate are correctly represented       
        # Using the non standard projection (50) 
        for i in range(4):
            zone, e, n = redfearn(lat_long_points[i][0],
                                  lat_long_points[i][1],
                                  central_meridian=123)
            assert num.allclose([x[i],y[i]], [e,n])
            assert zone==-1
        
        self.delete_mux(files)
 
    def test_Urs_points(self):
        time_step_count = 3
        time_step = 2
        lat_long_points =[(-21.5,114.5),(-21.5,115),(-21.,115)]
        base_name, files = self.write_mux(lat_long_points,
                                          time_step_count, time_step)
        for file in files:
            urs = Read_urs(file)
            assert time_step_count == urs.time_step_count
            assert time_step == urs.time_step

            for lat_lon, dep in map(None, lat_long_points, urs.lonlatdep):
                    _ , e, n = redfearn(lat_lon[0], lat_lon[1])
                    assert num.allclose(n, dep[2])
                        
            count = 0
            for slice in urs:
                count += 1
                #print slice
                for lat_lon, quantity in map(None, lat_long_points, slice):
                    _ , e, n = redfearn(lat_lon[0], lat_lon[1])
                    #print "quantity", quantity
                    #print "e", e
                    #print "n", n
                    if file[-5:] == WAVEHEIGHT_MUX_LABEL[-5:] or \
                           file[-5:] == NORTH_VELOCITY_LABEL[-5:] :
                        assert num.allclose(e, quantity)
                    if file[-5:] == EAST_VELOCITY_LABEL[-5:]:
                        assert num.allclose(n, quantity)
            assert count == time_step_count
                     
        self.delete_mux(files)        



        
################################################################################

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_Mux,'test')
    runner = unittest.TextTestRunner() #verbosity=2)
    runner.run(suite)
        

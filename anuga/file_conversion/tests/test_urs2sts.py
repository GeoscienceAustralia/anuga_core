

import numpy as num
import unittest
import tempfile
import os
import sys
from anuga.file.netcdf import NetCDFFile

from anuga.utilities.system_tools import get_pathname_from_package
from anuga.coordinate_transforms.geo_reference import Geo_reference     
from anuga.coordinate_transforms.redfearn import redfearn
from anuga.utilities.numerical_tools import ensure_numeric
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a
from anuga.file.sts import create_sts_boundary
from anuga.file.csv_file import load_csv_as_dict, load_csv_as_array

from anuga.shallow_water.shallow_water_domain import Domain

# boundary functions
from anuga.shallow_water.boundaries import Reflective_boundary, \
            Field_boundary, Transmissive_momentum_set_stage_boundary, \
            Transmissive_stage_zero_momentum_boundary
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions\
     import Transmissive_boundary, Dirichlet_boundary, \
            Time_boundary, File_boundary, AWI_boundary

from anuga.pmesh.mesh_interface import create_mesh_from_regions

from anuga.file_conversion.urs2sts import urs2sts

# Allow us to use helper methods from this test.
from anuga.file.tests.test_mux import Test_Mux

class Test_Urs2Sts(Test_Mux):
    """ A suite of tests to test urs2sts file conversion functions.
        These tests are quite coarse-grained: converting a file
        and checking that its headers and some of its contents
        are correct.
    """ 


    def tearDown(self):
        for file in ['domain.sww', 'urs_test_mesh.tsh' ]:
            try:
                os.remove(file)
            except:
                pass
                
    def test_urs2sts0(self):
        """
        Test single source
        """
        tide=0
        time_step_count = 3
        time_step = 2
        lat_long_points =[(-21.5,114.5),(-21,114.5),(-21.5,115), (-21.,115.)]
        n=len(lat_long_points)
        first_tstep=num.ones(n,int)
        first_tstep[0]+=1
        first_tstep[2]+=1
        last_tstep=(time_step_count)*num.ones(n,int)
        last_tstep[0]-=1

        gauge_depth=20*num.ones(n,float)
        ha=2*num.ones((n,time_step_count),float)
        ha[0]=num.arange(0,time_step_count)
        ha[1]=num.arange(time_step_count,2*time_step_count)
        ha[2]=num.arange(2*time_step_count,3*time_step_count)
        ha[3]=num.arange(3*time_step_count,4*time_step_count)
        ua=5*num.ones((n,time_step_count),float)
        va=-10*num.ones((n,time_step_count),float)

        base_name, files = self.write_mux2(lat_long_points,
                                      time_step_count, time_step,
                                      first_tstep, last_tstep,
                                      depth=gauge_depth,
                                      ha=ha,
                                      ua=ua,
                                      va=va)

        sts_file = base_name + '.sts'

        urs2sts(base_name,
                basename_out=sts_file, 
                mean_stage=tide,verbose=False)


        #Let's interigate the sww file
        # Note, the sww info is not gridded.  It is point data.
        fid = NetCDFFile(sts_file)

        # Make x and y absolute
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]

        geo_reference = Geo_reference(NetCDFObject=fid)
        points = geo_reference.get_absolute(list(zip(x, y)))
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

        depth=num.zeros((len(lat_long_points),time_step_count),float)
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
        first_tstep=num.ones(n,int)
        first_tstep[0]+=1
        first_tstep[2]+=1
        last_tstep=(time_step_count)*num.ones(n,int)
        last_tstep[0]-=1

        gauge_depth=20*num.ones(n,float)
        ha=2*num.ones((n,time_step_count),float)
        ha[0]=num.arange(0,time_step_count)
        ha[1]=num.arange(time_step_count,2*time_step_count)
        ha[2]=num.arange(2*time_step_count,3*time_step_count)
        ha[3]=num.arange(3*time_step_count,4*time_step_count)
        ua=5*num.ones((n,time_step_count),float)
        va=-10*num.ones((n,time_step_count),float)

        base_name, files = self.write_mux2(lat_long_points,
                                           time_step_count, time_step,
                                           first_tstep, last_tstep,
                                           depth=gauge_depth,
                                           ha=ha,
                                           ua=ua,
                                           va=va)

        sts_file = base_name + '.sts'

        urs2sts(base_name,
                basename_out=sts_file, 
                central_meridian=123,
                mean_stage=tide,
                verbose=False)

        #Let's interigate the sww file
        # Note, the sww info is not gridded.  It is point data.
        fid = NetCDFFile(sts_file)

        # Make x and y absolute
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]

        geo_reference = Geo_reference(NetCDFObject=fid)
        points = geo_reference.get_absolute(list(zip(x, y)))
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
        


    def test_urs2sts_individual_sources(self):   
        """Test that individual sources compare to actual urs output
           Test that the first recording time is the smallest
           over waveheight, easting and northing velocity
        """
        
        # Get path where this test is run
        path = get_pathname_from_package('anuga.shallow_water')        
        
        testdir = os.path.join(path, 'tests',  'urs_test_data')
        ordering_filename=os.path.join(testdir, 'thinned_bound_order_test.txt')
        
        sources = ['1-z.grd','2-z.grd','3-z.grd']
        
        # Start times by source and station taken manually from urs header files
        time_start_z = num.array([[10.0,11.5,13,14.5,17.7],
                                  [9.8,11.2,12.7,14.2,17.4],
                                  [9.5,10.9,12.4,13.9,17.1]])

        time_start_e = time_start_n = time_start_z

        # time step in urs output
        delta_t = 0.1
        
        # Make sts file for each source
        for k, source_filename in enumerate(sources):
            source_number = k + 1 # Source numbering starts at 1
            
            urs_filenames = os.path.join(testdir, source_filename)
            weights = [1.]
            sts_name_out = 'test'
            
            urs2sts(urs_filenames,
                    basename_out=sts_name_out,
                    ordering_filename=ordering_filename,
                    weights=weights,
                    mean_stage=0.0,
                    verbose=False)

            # Read in sts file for this source file
            fid = NetCDFFile(sts_name_out+'.sts', netcdf_mode_r) # Open existing file for read
            x = fid.variables['x'][:]+fid.xllcorner    # x-coordinates of vertices
            y = fid.variables['y'][:]+fid.yllcorner    # y-coordinates of vertices
            elevation = fid.variables['elevation'][:]
            time=fid.variables['time'][:]+fid.starttime

            # Get quantity data from sts file
            quantity_names=['stage','xmomentum','ymomentum']
            quantities = {}
            for i, name in enumerate(quantity_names):
                quantities[name] = fid.variables[name][:]
            
            # Make sure start time from sts file is the minimum starttime 
            # across all stations (for this source)
            #print k, time_start_z[k,:]
            starttime = min(time_start_z[k, :])
            sts_starttime = fid.starttime
            msg = 'sts starttime for source %d was %f. Should have been %f'\
                %(source_number, sts_starttime, starttime)
            assert num.allclose(sts_starttime, starttime), msg             

            # For each station, compare urs2sts output to known urs output
            for j in range(len(x)):
                index_start_urs = 0
                index_end_urs = 0
                index_start = 0
                index_end = 0
                count = 0

                # read in urs test data for stage, e and n velocity
                urs_file_name_z = 'z_'+str(source_number)+'_'+str(j)+'.csv'
                dict = load_csv_as_array(os.path.join(testdir, urs_file_name_z))
                urs_stage = dict['urs_stage']
                urs_file_name_e = 'e_'+str(source_number)+'_'+str(j)+'.csv'
                dict = load_csv_as_array(os.path.join(testdir, urs_file_name_e))
                urs_e = dict['urs_e']
                urs_file_name_n = 'n_'+str(source_number)+'_'+str(j)+'.csv'
                dict = load_csv_as_array(os.path.join(testdir, urs_file_name_n))
                urs_n = dict['urs_n']

                # find start and end time for stage             
                for i in range(len(urs_stage)):
                    if urs_stage[i] == 0.0:
                        index_start_urs_z = i+1
                    if int(urs_stage[i]) == 99 and count != 1:
                        count +=1
                        index_end_urs_z = i

                if count == 0: index_end_urs_z = len(urs_stage)

                start_times_z = time_start_z[source_number-1]

                # start times for easting velocities should be the same as stage
                start_times_e = time_start_e[source_number-1]
                index_start_urs_e = index_start_urs_z
                index_end_urs_e = index_end_urs_z

                # start times for northing velocities should be the same as stage
                start_times_n = time_start_n[source_number-1]
                index_start_urs_n = index_start_urs_z
                index_end_urs_n = index_end_urs_z
                
                # Check that actual start time matches header information for stage
                msg = 'stage start time from urs file is not the same as the '
                msg += 'header file for source %i and station %i' %(source_number,j)
                assert num.allclose(index_start_urs_z,start_times_z[j]/delta_t), msg

                msg = 'e velocity start time from urs file is not the same as the '
                msg += 'header file for source %i and station %i' %(source_number,j)
                assert num.allclose(index_start_urs_e,start_times_e[j]/delta_t), msg

                msg = 'n velocity start time from urs file is not the same as the '
                msg += 'header file for source %i and station %i' %(source_number,j)
                assert num.allclose(index_start_urs_n,start_times_n[j]/delta_t), msg
                
                # get index for start and end time for sts quantities
                index_start_stage = 0
                index_end_stage = 0
                count = 0
                sts_stage = quantities['stage'][:,j]
                for i in range(len(sts_stage)):
                    if sts_stage[i] != 0.0 and count != 1:
                        count += 1
                        index_start_stage = i
                    if int(sts_stage[i]) == 99 and count != 1:
                        count += 1
                        index_end_stage = i

                index_end_stage = index_start_stage + len(urs_stage[index_start_urs_z:index_end_urs_z])

                sts_xmom = quantities['xmomentum'][:,j]
                index_start_x = index_start_stage
                index_end_x = index_start_x + len(urs_e[index_start_urs_e:index_end_urs_e])

                sts_ymom = quantities['ymomentum'][:,j]
                index_start_y = index_start_stage
                index_end_y = index_start_y + len(urs_n[index_start_urs_n:index_end_urs_n])

                # check that urs stage and sts stage are the same
                msg = 'urs stage is not equal to sts stage for for source %i and station %i' %(source_number,j)
                assert num.allclose(urs_stage[index_start_urs_z:index_end_urs_z],
                                    sts_stage[index_start_stage:index_end_stage], 
                                    rtol=1.0e-6, atol=1.0e-5 ), msg                                
                                
                # check that urs e velocity and sts xmomentum are the same
                msg = 'urs e velocity is not equivalent to sts x momentum for for source %i and station %i' %(source_number,j)
                assert num.allclose(urs_e[index_start_urs_e:index_end_urs_e]*(urs_stage[index_start_urs_e:index_end_urs_e]-elevation[j]),
                                sts_xmom[index_start_x:index_end_x], 
                                rtol=1.0e-5, atol=1.0e-4 ), msg
                
                # check that urs n velocity and sts ymomentum are the same
                #print 'urs n velocity', urs_n[index_start_urs_n:index_end_urs_n]*(urs_stage[index_start_urs_n:index_end_urs_n]-elevation[j])
                #print 'sts momentum', sts_ymom[index_start_y:index_end_y]                                                             
                msg = 'urs n velocity is not equivalent to sts y momentum for source %i and station %i' %(source_number,j)
                assert num.allclose(urs_n[index_start_urs_n:index_end_urs_n]*(urs_stage[index_start_urs_n:index_end_urs_n]-elevation[j]),
                                -sts_ymom[index_start_y:index_end_y], 
                                rtol=1.0e-5, atol=1.0e-4 ), msg
                                                
                                
            fid.close()
            
        os.remove(sts_name_out+'.sts')

    def test_urs2sts_combined_sources(self):   
        """Test that combined sources compare to actual urs output
           Test that the first recording time is the smallest
           over waveheight, easting and northing velocity
        """

        # combined
        time_start_z = num.array([9.5,10.9,12.4,13.9,17.1])
        time_start_e = time_start_n = time_start_z
         
        # make sts file for combined sources
        weights = [1., 2., 3.]
        
        path = get_pathname_from_package('anuga.shallow_water')        
                
        testdir = os.path.join(path, 'tests', 'urs_test_data')        
        ordering_filename=os.path.join(testdir, 'thinned_bound_order_test.txt')

        urs_filenames = [os.path.join(testdir,'1-z.grd'),
                         os.path.join(testdir,'2-z.grd'),
                         os.path.join(testdir,'3-z.grd')]
        sts_name_out = 'test'
        
        urs2sts(urs_filenames,
                basename_out=sts_name_out,
                ordering_filename=ordering_filename,
                weights=weights,
                mean_stage=0.0,
                verbose=False)
        
        # read in sts file for combined source
        fid = NetCDFFile(sts_name_out+'.sts', netcdf_mode_r)    # Open existing file for read
        x = fid.variables['x'][:]+fid.xllcorner   # x-coordinates of vertices
        y = fid.variables['y'][:]+fid.yllcorner   # y-coordinates of vertices
        elevation = fid.variables['elevation'][:]
        time=fid.variables['time'][:]+fid.starttime
        
        
        # Check that stored permutation is as per default
        permutation = list(range(len(x)))
        stored_permutation = fid.variables['permutation'][:]
        msg = 'Permutation was not stored correctly. I got '
        msg += str(stored_permutation)
        assert num.allclose(stored_permutation, permutation), msg        

        # Get quantity data from sts file
        quantity_names=['stage','xmomentum','ymomentum']
        quantities = {}
        for i, name in enumerate(quantity_names):
            quantities[name] = fid.variables[name][:]

        # For each station, compare urs2sts output to known urs output   
        delta_t = 0.1
        
        # Make sure start time from sts file is the minimum starttime 
        # across all stations (for this source)
        starttime = min(time_start_z[:])
        sts_starttime = fid.starttime
        msg = 'sts starttime was %f. Should have been %f'\
            %(sts_starttime, starttime)
        assert num.allclose(sts_starttime, starttime), msg
    
        #stations = [1,2,3]
        #for j in stations: 
        for j in range(len(x)):
            index_start_urs_z = 0
            index_end_urs_z = 0
            index_start_urs_e = 0
            index_end_urs_e = 0
            index_start_urs_n = 0
            index_end_urs_n = 0
            count = 0

            # read in urs test data for stage, e and n velocity
            urs_file_name_z = 'z_combined_'+str(j)+'.csv'
            dict = load_csv_as_array(os.path.join(testdir, urs_file_name_z))
            urs_stage = dict['urs_stage']
            urs_file_name_e = 'e_combined_'+str(j)+'.csv'
            dict = load_csv_as_array(os.path.join(testdir, urs_file_name_e))
            urs_e = dict['urs_e']
            urs_file_name_n = 'n_combined_'+str(j)+'.csv'
            dict = load_csv_as_array(os.path.join(testdir, urs_file_name_n))
            urs_n = dict['urs_n']

            # find start and end time for stage         
            for i in range(len(urs_stage)):
                if urs_stage[i] == 0.0:
                    index_start_urs_z = i+1
                if int(urs_stage[i]) == 99 and count != 1:
                    count +=1
                    index_end_urs_z = i

            if count == 0: index_end_urs_z = len(urs_stage)

            start_times_z = time_start_z[j]

            start_times_e = time_start_e[j]
            index_start_urs_e = index_start_urs_z

            start_times_n = time_start_n[j]
            index_start_urs_n = index_start_urs_z
               
            # Check that actual start time matches header information for stage
            msg = 'stage start time from urs file is not the same as the '
            msg += 'header file at station %i' %(j)
            assert num.allclose(index_start_urs_z,start_times_z/delta_t), msg

            msg = 'e velocity start time from urs file is not the same as the '
            msg += 'header file at station %i' %(j)
            assert num.allclose(index_start_urs_e,start_times_e/delta_t), msg

            msg = 'n velocity start time from urs file is not the same as the '
            msg += 'header file at station %i' %(j)
            assert num.allclose(index_start_urs_n,start_times_n/delta_t), msg
                
            # get index for start and end time for sts quantities
            index_start_stage = 0
            index_end_stage = 0
            index_start_x = 0
            index_end_x = 0
            index_start_y = 0
            index_end_y = 0
            count = 0
            count1 = 0
            sts_stage = quantities['stage'][:,j]
            for i in range(len(sts_stage)):
                if sts_stage[i] != 0.0 and count != 1:
                    count += 1
                    index_start_stage = i
                if int(urs_stage[i]) == 99 and count != 1:
                    count +=1
                    index_end_stage = i
                
            index_end_stage = index_start_stage + len(urs_stage[index_start_urs_z:index_end_urs_z])

            index_start_x = index_start_stage
            index_end_x = index_start_x + len(urs_stage[index_start_urs_e:index_end_urs_e])
            sts_xmom = quantities['ymomentum'][:,j]

            index_start_y = index_start_stage
            index_end_y = index_start_y + len(urs_stage[index_start_urs_n:index_end_urs_n])
            sts_ymom = quantities['ymomentum'][:,j]

            # check that urs stage and sts stage are the same
            msg = 'urs stage is not equal to sts stage for station %i' %j
            #print 'urs stage', urs_stage[index_start_urs_z:index_end_urs_z]
            #print 'sts stage', sts_stage[index_start_stage:index_end_stage]
            #print 'diff', max(urs_stage[index_start_urs_z:index_end_urs_z]-sts_stage[index_start_stage:index_end_stage])
            #print 'index', index_start_stage, index_end_stage, len(sts_stage)
            assert num.allclose(urs_stage[index_start_urs_z:index_end_urs_z],
                            sts_stage[index_start_stage:index_end_stage], 
                                rtol=1.0e-5, atol=1.0e-4 ), msg                                
                                
            # check that urs e velocity and sts xmomentum are the same          
            msg = 'urs e velocity is not equivalent to sts xmomentum for station %i' %j
            assert num.allclose(urs_e[index_start_urs_e:index_end_urs_e]*(urs_stage[index_start_urs_e:index_end_urs_e]-elevation[j]),
                            sts_xmom[index_start_x:index_end_x], 
                            rtol=1.0e-5, atol=1.0e-4 ), msg
                
            # check that urs n velocity and sts ymomentum are the same                            
            msg = 'urs n velocity is not equivalent to sts ymomentum for station %i' %j
            assert num.allclose(urs_n[index_start_urs_n:index_end_urs_n]*(urs_stage[index_start_urs_n:index_end_urs_n]-elevation[j]),
                            sts_ymom[index_start_y:index_end_y], 
                            rtol=1.0e-5, atol=1.0e-4 ), msg

        fid.close()
        
        os.remove(sts_name_out+'.sts')
        
        
        
    def test_urs2sts_ordering(self):
        """Test multiple sources with ordering file
        """
        
        tide = 0.35
        time_step_count = 6 # I made this a little longer (Ole)
        time_step = 2
        lat_long_points =[(-21.5,114.5),(-21,114.5),(-21.5,115), (-21.,115.)]
        n=len(lat_long_points)
        first_tstep=num.ones(n,int)
        first_tstep[0]+=1
        first_tstep[2]+=1
        last_tstep=(time_step_count)*num.ones(n,int)
        last_tstep[0]-=1

        gauge_depth=20*num.ones(n,float)
        ha=2*num.ones((n,time_step_count),float)
        ha[0]=num.arange(0,time_step_count)
        ha[1]=num.arange(time_step_count,2*time_step_count)
        ha[2]=num.arange(2*time_step_count,3*time_step_count)
        ha[3]=num.arange(3*time_step_count,4*time_step_count)
        ua=5*num.ones((n,time_step_count),float)
        va=-10*num.ones((n,time_step_count),float)

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

                                               
        # Create ordering file
        permutation = [3,0,2]

        _, ordering_filename = tempfile.mkstemp('')
        order_fid = open(ordering_filename, 'w')  
        order_fid.write('index, longitude, latitude\n')
        for index in permutation:
            order_fid.write('%d, %f, %f\n' %(index, 
                                             lat_long_points[index][1], 
                                             lat_long_points[index][0]))
        order_fid.close()
        
            

                                               
        # Call urs2sts with multiple mux files
        urs2sts([base_nameI, base_nameII], 
                basename_out=base_nameI, 
                ordering_filename=ordering_filename,
                weights=[1.0, 1.0],
                mean_stage=tide,
                verbose=False)

        # now I want to check the sts file ...
        sts_file = base_nameI + '.sts'

        #Let's interrogate the sts file
        # Note, the sts info is not gridded.  It is point data.
        fid = NetCDFFile(sts_file)
        
        # Check that original indices have been stored
        stored_permutation = fid.variables['permutation'][:]
        msg = 'Permutation was not stored correctly. I got '
        msg += str(stored_permutation)
        assert num.allclose(stored_permutation, permutation), msg
        

        # Make x and y absolute
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]

        geo_reference = Geo_reference(NetCDFObject=fid)
        points = geo_reference.get_absolute(list(zip(x, y)))
        points = ensure_numeric(points)

        x = points[:,0]
        y = points[:,1]

        #print
        #print x
        #print y
        for i, index in enumerate(permutation):
            # Check that STS points are stored in the correct order
            
            # Work out the UTM coordinates sts point i
            zone, e, n = redfearn(lat_long_points[index][0], 
                                  lat_long_points[index][1])             

            #print i, [x[i],y[i]], [e,n]
            assert num.allclose([x[i],y[i]], [e,n])
            
                        
        # Check the time vector
        times = fid.variables['time'][:]

        times_actual = []
        for i in range(time_step_count):
            times_actual.append(time_step * i)

        assert num.allclose(ensure_numeric(times),
                            ensure_numeric(times_actual))
                        

        # Check sts values
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
        
        ha_permutation = num.take(ha, permutation, axis=0) 
        ua_permutation = num.take(ua, permutation, axis=0)
        va_permutation = num.take(va, permutation, axis=0)
        gauge_depth_permutation = num.take(gauge_depth, permutation, axis=0)

        
        assert num.allclose(2.0*num.transpose(ha_permutation)+tide, stage)  # Meters

        #Check the momentums - ua
        #momentum = velocity*(stage-elevation)
        # elevation = - depth
        #momentum = velocity_ua *(stage+depth)

        depth=num.zeros((len(lat_long_points),time_step_count),float)
        for i in range(len(lat_long_points)):
            depth[i]=gauge_depth[i]+tide+2.0*ha[i]
            #2.0*ha necessary because using two files with weights=1 are used
            
        depth_permutation = num.take(depth, permutation, axis=0)                     
        

        # The xmomentum stored in the .sts file should be the sum of the ua
        # in the two mux2 files multiplied by the depth.
        assert num.allclose(2.0*num.transpose(ua_permutation*depth_permutation), xmomentum) 

        #Check the momentums - va
        #momentum = velocity*(stage-elevation)
        # elevation = - depth
        #momentum = velocity_va *(stage+depth)

        # The ymomentum stored in the .sts file should be the sum of the va
        # in the two mux2 files multiplied by the depth.
        assert num.allclose(2.0*num.transpose(va_permutation*depth_permutation), ymomentum)

        # check the elevation values.
        # -ve since urs measures depth, sww meshers height,
        assert num.allclose(-gauge_depth_permutation, elevation)  #Meters

        fid.close()
        self.delete_mux(filesI)
        self.delete_mux(filesII)
        os.remove(sts_file)
        

        
        
        
    def Xtest_urs2sts_ordering_exception(self):
        """Test that inconsistent lats and lons in ordering file are caught.
        """
        
        tide=0
        time_step_count = 3
        time_step = 2
        lat_long_points =[(-21.5,114.5),(-21,114.5),(-21.5,115), (-21.,115.)]
        n=len(lat_long_points)
        first_tstep=num.ones(n,int)
        first_tstep[0]+=1
        first_tstep[2]+=1
        last_tstep=(time_step_count)*num.ones(n,int)
        last_tstep[0]-=1

        gauge_depth=20*num.ones(n,float)
        ha=2*num.ones((n,time_step_count),float)
        ha[0]=num.arange(0,time_step_count)
        ha[1]=num.arange(time_step_count,2*time_step_count)
        ha[2]=num.arange(2*time_step_count,3*time_step_count)
        ha[3]=num.arange(3*time_step_count,4*time_step_count)
        ua=5*num.ones((n,time_step_count),float)
        va=-10*num.ones((n,time_step_count),float)

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

                                               
        # Create ordering file
        permutation = [3,0,2]

        # Do it wrongly and check that exception is being raised
        _, ordering_filename = tempfile.mkstemp('')
        order_fid = open(ordering_filename, 'w')  
        order_fid.write('index, longitude, latitude\n')
        for index in permutation:
            order_fid.write('%d, %f, %f\n' %(index, 
                                             lat_long_points[index][0], 
                                             lat_long_points[index][1]))
        order_fid.close()
        
        try:
            urs2sts([base_nameI, base_nameII], 
                    basename_out=base_nameI, 
                    ordering_filename=ordering_filename,
                    weights=[1.0, 1.0],
                    mean_stage=tide,
                    verbose=False)  
            os.remove(ordering_filename)            
        except:
            pass
        else:
            msg = 'Should have caught wrong lat longs'
            raise(Exception, msg)

        
        self.delete_mux(filesI)
        self.delete_mux(filesII)

        

        
    def test_urs2sts_ordering_different_sources(self):
        """Test multiple sources with ordering file, different source files and weights.
           This test also has more variable values than the previous ones
        """
        
        tide = 1.5
        time_step_count = 10
        time_step = 0.2
        
        times_ref = num.arange(0, time_step_count*time_step, time_step)
        #print 'time vector', times_ref
        
        lat_long_points = [(-21.5,114.5), (-21,114.5), (-21.5,115), (-21.,115.), (-22., 117.)]
        n = len(lat_long_points)
        
        # Create non-trivial weights
        #weights = [0.8, 1.5] # OK
        #weights = [0.8, 10.5] # Fail (up to allclose tolerance)
        #weights = [10.5, 10.5] # OK
        #weights = [0.0, 10.5] # OK
        #weights = [0.8, 0.] # OK                
        #weights = [8, 0.1] # OK                        
        #weights = [0.8, 10.0] # OK                                
        #weights = [0.8, 10.6] # OK           
        weights = [3.8, 7.6] # OK                   
        #weights = [0.5, 0.5] # OK                           
        #weights = [2., 2.] # OK                            
        #weights = [0.0, 0.5] # OK                                          
        #weights = [1.0, 1.0] # OK                                                  
        
        
        # Create different timeseries starting and ending at different times 
        first_tstep=num.ones(n,int)
        first_tstep[0]+=2   # Point 0 starts at 2
        first_tstep[1]+=4   # Point 1 starts at 4        
        first_tstep[2]+=3   # Point 2 starts at 3
        
        last_tstep=(time_step_count)*num.ones(n,int)
        last_tstep[0]-=1    # Point 0 ends 1 step early
        last_tstep[1]-=2    # Point 1 ends 2 steps early                
        last_tstep[4]-=3    # Point 4 ends 3 steps early        
        
        #print
        #print 'time_step_count', time_step_count
        #print 'time_step', time_step
        #print 'first_tstep', first_tstep
        #print 'last_tstep', last_tstep                
        
        
        # Create varying elevation data (positive values for seafloor)
        gauge_depth=20*num.ones(n,float)
        for i in range(n):
            gauge_depth[i] += i**2
            
        #print 'gauge_depth', gauge_depth
        
        # Create data to be written to first mux file        
        ha0=2*num.ones((n,time_step_count),float)
        ha0[0]=num.arange(0,time_step_count)
        ha0[1]=num.arange(time_step_count,2*time_step_count)
        ha0[2]=num.arange(2*time_step_count,3*time_step_count)
        ha0[3]=num.arange(3*time_step_count,4*time_step_count)
        ua0=5*num.ones((n,time_step_count),float)
        va0=-10*num.ones((n,time_step_count),float)

        # Ensure data used to write mux file to be zero when gauges are
        # not recording
        for i in range(n):
             # For each point
             
             for j in list(range(0, first_tstep[i]-1)) + list(range(last_tstep[i], time_step_count)):
                 # For timesteps before and after recording range
                 ha0[i][j] = ua0[i][j] = va0[i][j] = 0.0                                  


                 
        #print 
        #print 'using varying start and end time'
        #print 'ha0', ha0[4]
        #print 'ua0', ua0
        #print 'va0', va0        
        
        # Write first mux file to be combined by urs2sts
        base_nameI, filesI = self.write_mux2(lat_long_points,
                                             time_step_count, time_step,
                                             first_tstep, last_tstep,
                                             depth=gauge_depth,
                                             ha=ha0,
                                             ua=ua0,
                                             va=va0)

                                             
                                             
        # Create data to be written to second mux file        
        ha1=num.ones((n,time_step_count),float)
        ha1[0]=num.sin(times_ref)
        ha1[1]=2*num.sin(times_ref - 3)
        ha1[2]=5*num.sin(4*times_ref)
        ha1[3]=num.sin(times_ref)
        ha1[4]=num.sin(2*times_ref-0.7)
                
        ua1=num.zeros((n,time_step_count),float)
        ua1[0]=3*num.cos(times_ref)        
        ua1[1]=2*num.sin(times_ref-0.7)   
        ua1[2]=num.arange(3*time_step_count,4*time_step_count)
        ua1[4]=2*num.ones(time_step_count)
        
        va1=num.zeros((n,time_step_count),float)
        va1[0]=2*num.cos(times_ref-0.87)        
        va1[1]=3*num.ones(time_step_count, int)       #array default#
        va1[3]=2*num.sin(times_ref-0.71)        
        
        
        # Ensure data used to write mux file to be zero when gauges are
        # not recording
        for i in range(n):
             # For each point
             
             for j in list(range(0, first_tstep[i]-1)) + list(range(last_tstep[i], time_step_count)):
                 # For timesteps before and after recording range
                 ha1[i][j] = ua1[i][j] = va1[i][j] = 0.0                                  


        #print 
        #print 'using varying start and end time'
        #print 'ha1', ha1[4]
        #print 'ua1', ua1
        #print 'va1', va1        
                                             
                                             
        # Write second mux file to be combined by urs2sts                                             
        base_nameII, filesII = self.write_mux2(lat_long_points,
                                               time_step_count, time_step,
                                               first_tstep, last_tstep,
                                               depth=gauge_depth,
                                               ha=ha1,
                                               ua=ua1,
                                               va=va1)

                                               
        # Create ordering file
        permutation = [4,0,2]

        _, ordering_filename = tempfile.mkstemp('')
        order_fid = open(ordering_filename, 'w')  
        order_fid.write('index, longitude, latitude\n')
        for index in permutation:
            order_fid.write('%d, %f, %f\n' %(index, 
                                             lat_long_points[index][1], 
                                             lat_long_points[index][0]))
        order_fid.close()
        
            

        #------------------------------------------------------------
        # Now read the mux files one by one without weights and test
        
        # Call urs2sts with mux file #0
        urs2sts([base_nameI], 
                basename_out=base_nameI, 
                ordering_filename=ordering_filename,
                mean_stage=tide,
                verbose=False)

        # Now read the sts file and check that values have been stored correctly.
        sts_file = base_nameI + '.sts'
        fid = NetCDFFile(sts_file)
        

        # Check that original indices have been stored
        stored_permutation = fid.variables['permutation'][:]
        msg = 'Permutation was not stored correctly. I got '
        msg += str(stored_permutation)
        assert num.allclose(stored_permutation, permutation), msg
        

        
        
        # Make x and y absolute
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]

        geo_reference = Geo_reference(NetCDFObject=fid)
        points = geo_reference.get_absolute(list(zip(x, y)))
        points = ensure_numeric(points)

        x = points[:,0]
        y = points[:,1]

        for i, index in enumerate(permutation):
            # Check that STS points are stored in the correct order
            
            # Work out the UTM coordinates sts point i
            zone, e, n = redfearn(lat_long_points[index][0], 
                                  lat_long_points[index][1])             

            #print i, [x[i],y[i]], [e,n]
            assert num.allclose([x[i],y[i]], [e,n])
            
                        
        # Check the time vector
        times = fid.variables['time'][:]
        assert num.allclose(ensure_numeric(times),
                            ensure_numeric(times_ref))
                        

        # Check sts values for mux #0
        stage0 = fid.variables['stage'][:].copy()
        xmomentum0 = fid.variables['xmomentum'][:].copy()
        ymomentum0 = fid.variables['ymomentum'][:].copy()
        elevation0 = fid.variables['elevation'][:].copy()

        
        #print 'stage', stage0
        #print 'xmomentum', xmomentum0
        #print 'ymomentum', ymomentum0        
        #print 'elevation', elevation0
        
        # The quantities stored in the .sts file should be the weighted sum of the 
        # quantities written to the mux2 files subject to the permutation vector.
        
        ha_ref = num.take(ha0, permutation, axis=0)
        ua_ref = num.take(ua0, permutation, axis=0)        
        va_ref = num.take(va0, permutation, axis=0)                

        gauge_depth_ref = num.take(gauge_depth, permutation, axis=0)                      
        
        assert num.allclose(num.transpose(ha_ref)+tide, stage0)  # Meters
        
        
        
        #Check the momentums - ua
        #momentum = velocity*(stage-elevation)
        # elevation = - depth
        #momentum = velocity_ua *(stage+depth)

        depth_ref = num.zeros((len(permutation), time_step_count), float)
        for i in range(len(permutation)):
            depth_ref[i]=gauge_depth_ref[i]+tide+ha_ref[i]


        # The xmomentum stored in the .sts file should be the sum of the ua
        # in the two mux2 files multiplied by the depth.
        assert num.allclose(num.transpose(ua_ref*depth_ref), xmomentum0) 

        #Check the momentums - va
        #momentum = velocity*(stage-elevation)
        # elevation = - depth
        #momentum = velocity_va *(stage+depth)

        # The ymomentum stored in the .sts file should be the sum of the va
        # in the two mux2 files multiplied by the depth.
        
        
        #print transpose(va_ref*depth_ref)
        #print ymomentum
        assert num.allclose(num.transpose(va_ref*depth_ref), ymomentum0)        

        # check the elevation values.
        # -ve since urs measures depth, sww meshers height,
        assert num.allclose(-gauge_depth_ref, elevation0) 

        fid.close()
        os.remove(sts_file)
        
        

        
        # Call urs2sts with mux file #1
        urs2sts([base_nameII], 
                basename_out=base_nameI, 
                ordering_filename=ordering_filename,
                mean_stage=tide,
                verbose=False)

        # Now read the sts file and check that values have been stored correctly.
        sts_file = base_nameI + '.sts'
        fid = NetCDFFile(sts_file)
        
        
        # Check that original indices have been stored
        stored_permutation = fid.variables['permutation'][:]
        msg = 'Permutation was not stored correctly. I got '
        msg += str(stored_permutation)
        assert num.allclose(stored_permutation, permutation), msg
        
        # Make x and y absolute
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]

        geo_reference = Geo_reference(NetCDFObject=fid)
        points = geo_reference.get_absolute(list(zip(x, y)))
        points = ensure_numeric(points)

        x = points[:,0]
        y = points[:,1]

        for i, index in enumerate(permutation):
            # Check that STS points are stored in the correct order
            
            # Work out the UTM coordinates sts point i
            zone, e, n = redfearn(lat_long_points[index][0], 
                                  lat_long_points[index][1])             

            #print i, [x[i],y[i]], [e,n]
            assert num.allclose([x[i],y[i]], [e,n])
            
                        
        # Check the time vector
        times = fid.variables['time'][:]
        assert num.allclose(ensure_numeric(times),
                            ensure_numeric(times_ref))
                        

        # Check sts values for mux #1 
        stage1 = fid.variables['stage'][:].copy()
        xmomentum1 = fid.variables['xmomentum'][:].copy()
        ymomentum1 = fid.variables['ymomentum'][:].copy()
        elevation1 = fid.variables['elevation'][:].copy()

        
        #print 'stage', stage1
        #print 'xmomentum', xmomentum1
        #print 'ymomentum', ymomentum1       
        #print 'elevation', elevation1
        
        # The quantities stored in the .sts file should be the weighted sum of the 
        # quantities written to the mux2 files subject to the permutation vector.
        
        ha_ref = num.take(ha1, permutation, axis=0)
        ua_ref = num.take(ua1, permutation, axis=0)        
        va_ref = num.take(va1, permutation, axis=0)                

        gauge_depth_ref = num.take(gauge_depth, permutation, axis=0)                         


        #print 
        #print stage1
        #print transpose(ha_ref)+tide - stage1
        

        assert num.allclose(num.transpose(ha_ref)+tide, stage1)  # Meters
        #import sys; sys.exit()

        #Check the momentums - ua
        #momentum = velocity*(stage-elevation)
        # elevation = - depth
        #momentum = velocity_ua *(stage+depth)

        depth_ref = num.zeros((len(permutation), time_step_count), float)
        for i in range(len(permutation)):
            depth_ref[i]=gauge_depth_ref[i]+tide+ha_ref[i]


        # The xmomentum stored in the .sts file should be the sum of the ua
        # in the two mux2 files multiplied by the depth.
        assert num.allclose(num.transpose(ua_ref*depth_ref), xmomentum1) 

        #Check the momentums - va
        #momentum = velocity*(stage-elevation)
        # elevation = - depth
        #momentum = velocity_va *(stage+depth)

        # The ymomentum stored in the .sts file should be the sum of the va
        # in the two mux2 files multiplied by the depth.
        
        
        #print transpose(va_ref*depth_ref)
        #print ymomentum
        assert num.allclose(num.transpose(va_ref*depth_ref), ymomentum1)        

        # check the elevation values.
        # -ve since urs measures depth, sww meshers height,
        assert num.allclose(-gauge_depth_ref, elevation1) 

        fid.close()
        os.remove(sts_file)
        
        #----------------------
        # Then read the mux files together and test
        
                                               
        # Call urs2sts with multiple mux files
        urs2sts([base_nameI, base_nameII], 
                basename_out=base_nameI, 
                ordering_filename=ordering_filename,
                weights=weights,
                mean_stage=tide,
                verbose=False)

        # Now read the sts file and check that values have been stored correctly.
        sts_file = base_nameI + '.sts'
        fid = NetCDFFile(sts_file)

        # Make x and y absolute
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]

        geo_reference = Geo_reference(NetCDFObject=fid)
        points = geo_reference.get_absolute(list(zip(x, y)))
        points = ensure_numeric(points)

        x = points[:,0]
        y = points[:,1]

        for i, index in enumerate(permutation):
            # Check that STS points are stored in the correct order
            
            # Work out the UTM coordinates sts point i
            zone, e, n = redfearn(lat_long_points[index][0], 
                                  lat_long_points[index][1])             

            #print i, [x[i],y[i]], [e,n]
            assert num.allclose([x[i],y[i]], [e,n])
            
                        
        # Check the time vector
        times = fid.variables['time'][:]
        assert num.allclose(ensure_numeric(times),
                            ensure_numeric(times_ref))
                        

        # Check sts values
        stage = fid.variables['stage'][:].copy()
        xmomentum = fid.variables['xmomentum'][:].copy()
        ymomentum = fid.variables['ymomentum'][:].copy()
        elevation = fid.variables['elevation'][:].copy()

        
        #print 'stage', stage
        #print 'elevation', elevation
        
        # The quantities stored in the .sts file should be the weighted sum of the 
        # quantities written to the mux2 files subject to the permutation vector.
        
        ha_ref = (weights[0]*num.take(ha0, permutation, axis=0)
                  + weights[1]*num.take(ha1, permutation, axis=0))
        ua_ref = (weights[0]*num.take(ua0, permutation, axis=0)
                  + weights[1]*num.take(ua1, permutation, axis=0))
        va_ref = (weights[0]*num.take(va0, permutation, axis=0)
                  + weights[1]*num.take(va1, permutation, axis=0))

        gauge_depth_ref = num.take(gauge_depth, permutation, axis=0)                         


        #print 
        #print stage
        #print transpose(ha_ref)+tide - stage

        assert num.allclose(num.transpose(ha_ref)+tide, stage)  # Meters

        #Check the momentums - ua
        #momentum = velocity*(stage-elevation)
        # elevation = - depth
        #momentum = velocity_ua *(stage+depth)

        depth_ref = num.zeros((len(permutation), time_step_count), float)
        for i in range(len(permutation)):
            depth_ref[i]=gauge_depth_ref[i]+tide+ha_ref[i]

            
        

        # The xmomentum stored in the .sts file should be the sum of the ua
        # in the two mux2 files multiplied by the depth.
        assert num.allclose(num.transpose(ua_ref*depth_ref), xmomentum) 

        #Check the momentums - va
        #momentum = velocity*(stage-elevation)
        # elevation = - depth
        #momentum = velocity_va *(stage+depth)

        # The ymomentum stored in the .sts file should be the sum of the va
        # in the two mux2 files multiplied by the depth.
        
        
        #print transpose(va_ref*depth_ref)
        #print ymomentum

        assert num.allclose(num.transpose(va_ref*depth_ref), ymomentum)

        # check the elevation values.
        # -ve since urs measures depth, sww meshers height,
        assert num.allclose(-gauge_depth_ref, elevation)  #Meters

        fid.close()
        self.delete_mux(filesI)
        self.delete_mux(filesII)
        os.remove(sts_file)
        
        #---------------
        # "Manually" add the timeseries up with weights and test
        # Tide is discounted from individual results and added back in       
        #

        stage_man = weights[0]*(stage0-tide) + weights[1]*(stage1-tide) + tide
        assert num.allclose(stage_man, stage)
                
        
    def test_file_boundary_stsI(self):
        """test_file_boundary_stsI(self):
        """
        
        bounding_polygon=[[6.0,97.0],[6.01,97.0],[6.02,97.0],[6.02,97.02],[6.00,97.02]]
        tide = 0.37
        time_step_count = 5
        time_step = 2
        lat_long_points =bounding_polygon[0:3]
        n=len(lat_long_points)
        first_tstep=num.ones(n,int)
        last_tstep=(time_step_count)*num.ones(n,int)

        h = 20        
        w = 2
        u = 10
        v = -10
        gauge_depth=h*num.ones(n,float)
        ha=w*num.ones((n,time_step_count),float)
        ua=u*num.ones((n,time_step_count),float)
        va=v*num.ones((n,time_step_count),float)
        base_name, files = self.write_mux2(lat_long_points,
                                           time_step_count, time_step,
                                           first_tstep, last_tstep,
                                           depth=gauge_depth,
                                           ha=ha,
                                           ua=ua,
                                           va=va)

        sts_file=base_name
        urs2sts(base_name,
                sts_file,
                mean_stage=tide,
                verbose=False)
        self.delete_mux(files)

        #print 'start create mesh from regions'
        for i in range(len(bounding_polygon)):
            zone,bounding_polygon[i][0],bounding_polygon[i][1]=redfearn(bounding_polygon[i][0],bounding_polygon[i][1])
        extent_res=1000000
        meshname = 'urs_test_mesh' + '.tsh'
        interior_regions=None
        boundary_tags={'ocean': [0,1], 'otherocean': [2,3,4]}
        create_mesh_from_regions(bounding_polygon,
                                 boundary_tags=boundary_tags,
                                 maximum_triangle_area=extent_res,
                                 filename=meshname,
                                 interior_regions=interior_regions,
                                 verbose=False)
        
        domain_fbound = Domain(meshname)
        domain_fbound.set_quantity('stage', tide)
        Bf = File_boundary(sts_file+'.sts', 
                           domain_fbound, 
                           boundary_polygon=bounding_polygon)
        Br = Reflective_boundary(domain_fbound)
        Bd = Dirichlet_boundary([w+tide, u*(w+h+tide), v*(w+h+tide)])        

        domain_fbound.set_boundary({'ocean': Bf,'otherocean': Br})
        
        # Check boundary object evaluates as it should
        for i, ((vol_id, edge_id), B) in enumerate(domain_fbound.boundary_objects):
            if B is Bf:
            
                qf = B.evaluate(vol_id, edge_id)  # File boundary
                qd = Bd.evaluate(vol_id, edge_id) # Dirichlet boundary

                assert num.allclose(qf, qd) 
                
        
        # Evolve
        finaltime=time_step*(time_step_count-1)
        yieldstep=time_step
        temp_fbound=num.zeros(int(finaltime/yieldstep)+1,float)

        for i, t in enumerate(domain_fbound.evolve(yieldstep=yieldstep,
                                                   finaltime=finaltime, 
                                                   skip_initial_step=False)):
                                                   
            D = domain_fbound
            temp_fbound[i]=D.quantities['stage'].centroid_values[2]

            # Check that file boundary object has populated 
            # boundary array correctly  
            # FIXME (Ole): Do this for the other tests too!
            for j, val in enumerate(D.get_quantity('stage').boundary_values):
            
                (vol_id, edge_id), B = D.boundary_objects[j]
                if isinstance(B, File_boundary):
                    #print j, val
                    assert num.allclose(val, w + tide)


        
        domain_drchlt = Domain(meshname)
        domain_drchlt.set_starttime(2.0)
        domain_drchlt.set_quantity('stage', tide)
        Br = Reflective_boundary(domain_drchlt)

        domain_drchlt.set_boundary({'ocean': Bd,'otherocean': Br})
        temp_drchlt=num.zeros(int(finaltime/yieldstep)+1,float)

        for i, t in enumerate(domain_drchlt.evolve(yieldstep=yieldstep,
                                                   finaltime=finaltime, 
                                                   skip_initial_step=False)):
            temp_drchlt[i]=domain_drchlt.quantities['stage'].centroid_values[2]

        #print domain_fbound.quantities['stage'].vertex_values
        #print domain_drchlt.quantities['stage'].vertex_values
                    
        assert num.allclose(temp_fbound,temp_drchlt)
        
        assert num.allclose(domain_fbound.quantities['stage'].vertex_values,
                            domain_drchlt.quantities['stage'].vertex_values)
                        
        assert num.allclose(domain_fbound.quantities['xmomentum'].vertex_values,
                            domain_drchlt.quantities['xmomentum'].vertex_values) 
                        
        assert num.allclose(domain_fbound.quantities['ymomentum'].vertex_values,
                            domain_drchlt.quantities['ymomentum'].vertex_values)
        
        
        os.remove(sts_file+'.sts')
        os.remove(meshname)
                
        
    def test_file_boundary_stsI_beyond_model_time(self):
        """test_file_boundary_stsI(self):
        
        Test that file_boundary can work when model time
        exceeds available data.
        This is optional and requires the user to specify a default 
        boundary object.
        """
        
        # Don't do warnings in unit test
        import warnings
        warnings.simplefilter('ignore')

        bounding_polygon=[[6.0,97.0],[6.01,97.0],[6.02,97.0],
                          [6.02,97.02],[6.00,97.02]]
        tide = 0.37
        time_step_count = 5
        time_step = 2
        lat_long_points = bounding_polygon[0:3]
        n=len(lat_long_points)
        first_tstep=num.ones(n,int)
        last_tstep=(time_step_count)*num.ones(n,int)

        h = 20        
        w = 2
        u = 10
        v = -10
        gauge_depth=h*num.ones(n,float)
        ha=w*num.ones((n,time_step_count),float)
        ua=u*num.ones((n,time_step_count),float)
        va=v*num.ones((n,time_step_count),float)
        base_name, files = self.write_mux2(lat_long_points,
                                           time_step_count, time_step,
                                           first_tstep, last_tstep,
                                           depth=gauge_depth,
                                           ha=ha,
                                           ua=ua,
                                           va=va)

        sts_file=base_name
        urs2sts(base_name,
                sts_file,
                mean_stage=tide,
                verbose=False)
        self.delete_mux(files)

        #print 'start create mesh from regions'
        for i in range(len(bounding_polygon)):
            zone,\
            bounding_polygon[i][0],\
            bounding_polygon[i][1]=redfearn(bounding_polygon[i][0],
                                            bounding_polygon[i][1])
                                            
        extent_res=1000000
        meshname = 'urs_test_mesh' + '.tsh'
        interior_regions=None
        boundary_tags={'ocean': [0,1], 'otherocean': [2,3,4]}
        create_mesh_from_regions(bounding_polygon,
                                 boundary_tags=boundary_tags,
                                 maximum_triangle_area=extent_res,
                                 filename=meshname,
                                 interior_regions=interior_regions,
                                 verbose=False)
        
        domain_fbound = Domain(meshname)
        domain_fbound.set_quantity('stage', tide)
        
        Br = Reflective_boundary(domain_fbound)
        Bd = Dirichlet_boundary([w+tide, u*(w+h+tide), v*(w+h+tide)])        
        Bf = File_boundary(sts_file+'.sts', 
                           domain_fbound, 
                           boundary_polygon=bounding_polygon,
                           default_boundary=Bd) # Condition to be used 
                                                # if model time exceeds 
                                                # available data

        domain_fbound.set_boundary({'ocean': Bf,'otherocean': Br})
        
        # Check boundary object evaluates as it should
        for i, ((vol_id, edge_id), B) in enumerate(domain_fbound.boundary_objects):
            if B is Bf:
            
                qf = B.evaluate(vol_id, edge_id)  # File boundary
                qd = Bd.evaluate(vol_id, edge_id) # Dirichlet boundary

                assert num.allclose(qf, qd) 
                
        
        # Evolve
        data_finaltime = time_step*(time_step_count-1)
        finaltime = data_finaltime + 10 # Let model time exceed available data
        yieldstep = time_step
        temp_fbound=num.zeros(int(finaltime/yieldstep)+1, float)

        for i, t in enumerate(domain_fbound.evolve(yieldstep=yieldstep,
                                                   finaltime=finaltime, 
                                                   skip_initial_step=False)):
                                                   
            D = domain_fbound
            temp_fbound[i]=D.quantities['stage'].centroid_values[2]

            # Check that file boundary object has populated 
            # boundary array correctly  
            # FIXME (Ole): Do this for the other tests too!
            for j, val in enumerate(D.get_quantity('stage').boundary_values):
            
                (vol_id, edge_id), B = D.boundary_objects[j]
                if isinstance(B, File_boundary):
                    #print j, val
                    assert num.allclose(val, w + tide)


        domain_drchlt = Domain(meshname)
        domain_drchlt.set_starttime(2.0)
        domain_drchlt.set_quantity('stage', tide)
        Br = Reflective_boundary(domain_drchlt)

        domain_drchlt.set_boundary({'ocean': Bd,'otherocean': Br})
        temp_drchlt=num.zeros(int(finaltime/yieldstep)+1,float)

        for i, t in enumerate(domain_drchlt.evolve(yieldstep=yieldstep,
                                                   finaltime=finaltime, 
                                                   skip_initial_step=False)):
            temp_drchlt[i]=domain_drchlt.quantities['stage'].centroid_values[2]

        #print domain_fbound.quantities['stage'].vertex_values
        #print domain_drchlt.quantities['stage'].vertex_values
                    
        assert num.allclose(temp_fbound,temp_drchlt)
        
        assert num.allclose(domain_fbound.quantities['stage'].vertex_values,
                            domain_drchlt.quantities['stage'].vertex_values)
                        
        assert num.allclose(domain_fbound.quantities['xmomentum'].vertex_values,
                            domain_drchlt.quantities['xmomentum'].vertex_values) 
                        
        assert num.allclose(domain_fbound.quantities['ymomentum'].vertex_values,
                            domain_drchlt.quantities['ymomentum'].vertex_values) 
        
        os.remove(sts_file+'.sts')
        os.remove(meshname)
                
        
    def test_field_boundary_stsI_beyond_model_time(self):
        """test_field_boundary(self):
        
        Test that field_boundary can work when model time
        exceeds available data whilst adjusting mean_stage.
        
        """
        
        # Don't do warnings in unit test
        import warnings
        warnings.simplefilter('ignore')

        bounding_polygon=[[6.0,97.0],[6.01,97.0],[6.02,97.0],
                          [6.02,97.02],[6.00,97.02]]
        tide = 0.37
        time_step_count = 5
        time_step = 2
        lat_long_points = bounding_polygon[0:3]
        n=len(lat_long_points)
        first_tstep=num.ones(n,int)
        last_tstep=(time_step_count)*num.ones(n,int)

        h = 20        
        w = 2
        u = 10
        v = -10
        gauge_depth=h*num.ones(n,float)
        ha=w*num.ones((n,time_step_count),float)
        ua=u*num.ones((n,time_step_count),float)
        va=v*num.ones((n,time_step_count),float)
        base_name, files = self.write_mux2(lat_long_points,
                                           time_step_count, time_step,
                                           first_tstep, last_tstep,
                                           depth=gauge_depth,
                                           ha=ha,
                                           ua=ua,
                                           va=va)

        sts_file=base_name
        urs2sts(base_name,
                sts_file,
                mean_stage=0.0, # Deliberately let Field_boundary do the adjustment
                verbose=False)
        self.delete_mux(files)

        #print 'start create mesh from regions'
        for i in range(len(bounding_polygon)):
            zone,\
            bounding_polygon[i][0],\
            bounding_polygon[i][1]=redfearn(bounding_polygon[i][0],
                                            bounding_polygon[i][1])
                                            
        extent_res=1000000
        meshname = 'urs_test_mesh' + '.tsh'
        interior_regions=None
        boundary_tags={'ocean': [0,1], 'otherocean': [2,3,4]}
        create_mesh_from_regions(bounding_polygon,
                                 boundary_tags=boundary_tags,
                                 maximum_triangle_area=extent_res,
                                 filename=meshname,
                                 interior_regions=interior_regions,
                                 verbose=False)
        
        domain_fbound = Domain(meshname)
        domain_fbound.set_quantity('stage', tide)
        
        Br = Reflective_boundary(domain_fbound)
        Bd = Dirichlet_boundary([w+tide, u*(w+h), v*(w+h)])
        Bdefault = Dirichlet_boundary([w, u*(w+h), v*(w+h)])        
                
        Bf = Field_boundary(sts_file+'.sts', 
                           domain_fbound, 
                           mean_stage=tide, # Field boundary to adjust for tide
                           boundary_polygon=bounding_polygon,
                           default_boundary=Bdefault) # Condition to be used 
                                                      # if model time exceeds 
                                                      # available data

        domain_fbound.set_boundary({'ocean': Bf,'otherocean': Br})
        
        # Check boundary object evaluates as it should
        for i, ((vol_id, edge_id), B) in enumerate(domain_fbound.boundary_objects):
            if B is Bf:
            
                qf = B.evaluate(vol_id, edge_id)  # Field boundary
                qd = Bd.evaluate(vol_id, edge_id) # Dirichlet boundary
                
                msg = 'Got %s, should have been %s' %(qf, qd)
                assert num.allclose(qf, qd), msg 
                
        # Evolve
        data_finaltime = time_step*(time_step_count-1)
        finaltime = data_finaltime + 10 # Let model time exceed available data
        yieldstep = time_step
        temp_fbound=num.zeros(int(finaltime/yieldstep)+1, float)

        for i, t in enumerate(domain_fbound.evolve(yieldstep=yieldstep,
                                                   finaltime=finaltime, 
                                                   skip_initial_step=False)):
                                                   
            D = domain_fbound
            #temp_fbound[i]=D.quantities['stage'].centroid_values[2]
            
            # Check that file boundary object has populated 
            # boundary array correctly  
            # FIXME (Ole): Do this for the other tests too!
            for j, val in enumerate(D.get_quantity('stage').boundary_values):
            
                (vol_id, edge_id), B = D.boundary_objects[j]
                if isinstance(B, Field_boundary):
                    msg = 'Got %f should have been %f' %(val, w+tide)
                    assert num.allclose(val, w + tide), msg


    def test_file_boundary_stsII(self):
        """test_file_boundary_stsII(self):
        
         mux2 file has points not included in boundary
         mux2 gauges are not stored with the same order as they are 
         found in bounding_polygon. This does not matter as long as bounding
         polygon passed to file_function contains the mux2 points needed (in
         the correct order).
         """

        bounding_polygon=[[6.01,97.0],[6.02,97.0],[6.02,97.02],[6.00,97.02],[6.0,97.0]]
        tide = -2.20 
        time_step_count = 20
        time_step = 2
        lat_long_points=bounding_polygon[0:2]
        lat_long_points.insert(0,bounding_polygon[len(bounding_polygon)-1])
        lat_long_points.insert(0,[6.0,97.01])
        n=len(lat_long_points)
        first_tstep=num.ones(n,int)
        last_tstep=(time_step_count)*num.ones(n,int)
        gauge_depth=20*num.ones(n,float)
        ha=2*num.ones((n,time_step_count),float)
        ua=10*num.ones((n,time_step_count),float)
        va=-10*num.ones((n,time_step_count),float)
        base_name, files = self.write_mux2(lat_long_points,
                                           time_step_count,
                                           time_step,
                                           first_tstep,
                                           last_tstep,
                                           depth=gauge_depth,
                                           ha=ha,
                                           ua=ua,
                                           va=va)

        sts_file=base_name
        urs2sts(base_name,sts_file,mean_stage=tide,verbose=False)
        self.delete_mux(files)

        #print 'start create mesh from regions'
        for i in range(len(bounding_polygon)):
            zone,\
            bounding_polygon[i][0],\
            bounding_polygon[i][1]=redfearn(bounding_polygon[i][0],
                                            bounding_polygon[i][1])
            
        extent_res=1000000
        meshname = 'urs_test_mesh' + '.tsh'
        interior_regions=None
        boundary_tags={'ocean': [0,1], 'otherocean': [2,3,4]}
        # have to change boundary tags from last example because now bounding
        # polygon starts in different place.
        create_mesh_from_regions(bounding_polygon,boundary_tags=boundary_tags,
                         maximum_triangle_area=extent_res,filename=meshname,
                         interior_regions=interior_regions,verbose=False)
        
        domain_fbound = Domain(meshname)
        domain_fbound.set_quantity('stage', tide)
        Bf = File_boundary(sts_file+'.sts',
                           domain_fbound,
                           boundary_polygon=bounding_polygon)
        Br = Reflective_boundary(domain_fbound)

        domain_fbound.set_boundary({'ocean': Bf,'otherocean': Br})
        finaltime=time_step*(time_step_count-1)
        yieldstep=time_step
        temp_fbound=num.zeros(int(finaltime/yieldstep)+1,float)
        
        for i, t in enumerate(domain_fbound.evolve(yieldstep=yieldstep,
                                                   finaltime=finaltime, 
                                                   skip_initial_step = False)):
            temp_fbound[i]=domain_fbound.quantities['stage'].centroid_values[2]
        
        domain_drchlt = Domain(meshname)
        domain_drchlt.set_quantity('stage', tide)
        Br = Reflective_boundary(domain_drchlt)
        Bd = Dirichlet_boundary([2.0+tide,220+10*tide,-220-10*tide])
        domain_drchlt.set_boundary({'ocean': Bd,'otherocean': Br})
        temp_drchlt=num.zeros(int(finaltime/yieldstep)+1,float)

        for i, t in enumerate(domain_drchlt.evolve(yieldstep=yieldstep,
                                                   finaltime=finaltime, 
                                                   skip_initial_step = False)):
            temp_drchlt[i]=domain_drchlt.quantities['stage'].centroid_values[2]


        assert num.allclose(temp_fbound,temp_drchlt)            
            
        #print domain_fbound.quantities['stage'].vertex_values
        #print domain_drchlt.quantities['stage'].vertex_values
                    
            
        assert num.allclose(domain_fbound.quantities['stage'].vertex_values,
                            domain_drchlt.quantities['stage'].vertex_values)
                        
        assert num.allclose(domain_fbound.quantities['xmomentum'].vertex_values,
                            domain_drchlt.quantities['xmomentum'].vertex_values)
                        
        assert num.allclose(domain_fbound.quantities['ymomentum'].vertex_values,
                            domain_drchlt.quantities['ymomentum'].vertex_values)
            
            

        os.remove(sts_file+'.sts')
        os.remove(meshname)

        
        
    def test_file_boundary_stsIII_ordering(self):
        """test_file_boundary_stsIII_ordering(self):
        Read correct points from ordering file and apply sts to boundary
        """

        lat_long_points=[[6.01,97.0],[6.02,97.0],[6.05,96.9],[6.0,97.0]]
        bounding_polygon=[[6.0,97.0],[6.01,97.0],[6.02,97.0],
                          [6.02,97.02],[6.00,97.02]]
        tide = 3.0
        time_step_count = 50
        time_step = 2
        n=len(lat_long_points)
        first_tstep=num.ones(n,int)
        last_tstep=(time_step_count)*num.ones(n,int)
        gauge_depth=20*num.ones(n,float)
        ha=2*num.ones((n,time_step_count),float)
        ua=10*num.ones((n,time_step_count),float)
        va=-10*num.ones((n,time_step_count),float)
        base_name, files = self.write_mux2(lat_long_points,
                                           time_step_count,
                                           time_step,
                                           first_tstep,
                                           last_tstep,
                                           depth=gauge_depth,
                                           ha=ha,
                                           ua=ua,
                                           va=va)
        # base name will not exist, but 3 other files are created

        # Write order file
        file_handle, order_base_name = tempfile.mkstemp("")
        os.close(file_handle)
        os.remove(order_base_name)
        d=","
        order_file=order_base_name+'order.txt'
        fid=open(order_file,'w')
        
        # Write Header
        header='index, longitude, latitude\n'
        fid.write(header)
        indices=[3,0,1]
        for i in indices:
            line=str(i)+d+str(lat_long_points[i][1])+d+\
                str(lat_long_points[i][0])+"\n"
            fid.write(line)
        fid.close()

        sts_file=base_name
        urs2sts(base_name,
                basename_out=sts_file,
                ordering_filename=order_file,
                mean_stage=tide,
                verbose=False)
        self.delete_mux(files)

        assert(os.access(sts_file+'.sts', os.F_OK))

        boundary_polygon = create_sts_boundary(base_name)

        os.remove(order_file)

        # Append the remaining part of the boundary polygon to be defined by
        # the user
        bounding_polygon_utm=[]
        for point in bounding_polygon:
            zone,easting,northing=redfearn(point[0],point[1])
            bounding_polygon_utm.append([easting,northing])

        boundary_polygon.append(bounding_polygon_utm[3])
        boundary_polygon.append(bounding_polygon_utm[4])

        plot=False
        if plot:
            from pylab import plot,show,axis
            boundary_polygon=ensure_numeric(boundary_polygon)
            bounding_polygon_utm=ensure_numeric(bounding_polygon_utm)
            #lat_long_points=ensure_numeric(lat_long_points)
            #plot(lat_long_points[:,0],lat_long_points[:,1],'o')
            plot(boundary_polygon[:,0], boundary_polygon[:,1],'d')
            plot(bounding_polygon_utm[:,0],bounding_polygon_utm[:,1],'o')
            show()

        assert num.allclose(bounding_polygon_utm,boundary_polygon)


        extent_res=1000000
        meshname = 'urs_test_mesh' + '.tsh'
        interior_regions=None
        boundary_tags={'ocean': [0,1], 'otherocean': [2,3,4]}
        
        # have to change boundary tags from last example because now bounding
        # polygon starts in different place.
        create_mesh_from_regions(boundary_polygon,boundary_tags=boundary_tags,
                         maximum_triangle_area=extent_res,filename=meshname,
                         interior_regions=interior_regions,verbose=False)
        
        domain_fbound = Domain(meshname)
        domain_fbound.set_quantity('stage', tide)
        Bf = File_boundary(sts_file+'.sts',
                           domain_fbound,
                           boundary_polygon=boundary_polygon)
        Br = Reflective_boundary(domain_fbound)

        domain_fbound.set_boundary({'ocean': Bf,'otherocean': Br})
        finaltime=time_step*(time_step_count-1)
        yieldstep=time_step
        temp_fbound=num.zeros(int(finaltime/yieldstep)+1,float)
    
        for i, t in enumerate(domain_fbound.evolve(yieldstep=yieldstep,
                                                   finaltime=finaltime, 
                                                   skip_initial_step = False)):
            temp_fbound[i]=domain_fbound.quantities['stage'].centroid_values[2]
    
        
        domain_drchlt = Domain(meshname)
        domain_drchlt.set_starttime(2.0)
        domain_drchlt.set_quantity('stage', tide)
        Br = Reflective_boundary(domain_drchlt)
        Bd = Dirichlet_boundary([2.0+tide,220+10*tide,-220-10*tide])
        domain_drchlt.set_boundary({'ocean': Bd,'otherocean': Br})
        temp_drchlt=num.zeros(int(finaltime/yieldstep)+1,float)
        
        for i, t in enumerate(domain_drchlt.evolve(yieldstep=yieldstep,
                                                   finaltime=finaltime, 
                                                   skip_initial_step = False)):
            temp_drchlt[i]=domain_drchlt.quantities['stage'].centroid_values[2]

        
        #print domain_fbound.quantities['stage'].vertex_values
        #print domain_drchlt.quantities['stage'].vertex_values
                    
        assert num.allclose(temp_fbound,temp_drchlt)

        
        assert num.allclose(domain_fbound.quantities['stage'].vertex_values,
                            domain_drchlt.quantities['stage'].vertex_values)
                        
        assert num.allclose(domain_fbound.quantities['xmomentum'].vertex_values,
                            domain_drchlt.quantities['xmomentum'].vertex_values)                        
                        
        assert num.allclose(domain_fbound.quantities['ymomentum'].vertex_values,
                            domain_drchlt.quantities['ymomentum'].vertex_values)
        
        # Use known Dirichlet condition (if sufficient timesteps have been taken)

        #FIXME: What do these assertions test? Also do they assume tide =0
        #print domain_fbound.quantities['stage'].vertex_values
        #assert allclose(domain_drchlt.quantities['stage'].vertex_values[6], 2)        
        #assert allclose(domain_fbound.quantities['stage'].vertex_values[6], 2)

        if not sys.platform == 'win32':
            os.remove(sts_file+'.sts')
        
        os.remove(meshname)
        
        
        
    def test_file_boundary_sts_time_limit(self):
        """test_file_boundary_stsIV_sinewave_ordering(self):
        Read correct points from ordering file and apply sts to boundary
        This one uses a sine wave and compares to time boundary
        
        This one tests that times used can be limited by upper limit
        """
        
        lat_long_points=[[6.01,97.0],[6.02,97.0],[6.05,96.9],[6.0,97.0]]
        bounding_polygon=[[6.0,97.0],[6.01,97.0],[6.02,97.0],[6.02,97.02],[6.00,97.02]]
        tide = 0.35
        time_step_count = 50
        time_step = 0.1
        times_ref = num.arange(0, time_step_count*time_step, time_step)
        
        n=len(lat_long_points)
        first_tstep=num.ones(n,int)
        last_tstep=(time_step_count)*num.ones(n,int)
        
        gauge_depth=20*num.ones(n,float)
        
        ha1=num.ones((n,time_step_count),float)
        ua1=3.*num.ones((n,time_step_count),float)
        va1=2.*num.ones((n,time_step_count),float)
        for i in range(n):
            ha1[i]=num.sin(times_ref)
        
        
        base_name, files = self.write_mux2(lat_long_points,
                                           time_step_count, time_step,
                                           first_tstep, last_tstep,
                                           depth=gauge_depth,
                                           ha=ha1,
                                           ua=ua1,
                                           va=va1)

        # Write order file
        file_handle, order_base_name = tempfile.mkstemp("")
        os.close(file_handle)
        os.remove(order_base_name)
        d=","
        order_file=order_base_name+'order.txt'
        fid=open(order_file,'w')
        
        # Write Header
        header='index, longitude, latitude\n'
        fid.write(header)
        indices=[3,0,1]
        for i in indices:
            line=str(i)+d+str(lat_long_points[i][1])+d+\
                str(lat_long_points[i][0])+"\n"
            fid.write(line)
        fid.close()

        sts_file=base_name
        urs2sts(base_name, basename_out=sts_file,
                ordering_filename=order_file,
                mean_stage=tide,
                verbose=False)
        self.delete_mux(files)
        
        
        
        # Now read the sts file and check that values have been stored correctly.
        fid = NetCDFFile(sts_file + '.sts')

        # Check the time vector
        times = fid.variables['time'][:]
        starttime = fid.starttime
        #print times
        #print starttime

        # Check sts quantities
        stage = fid.variables['stage'][:]
        xmomentum = fid.variables['xmomentum'][:]
        ymomentum = fid.variables['ymomentum'][:]
        elevation = fid.variables['elevation'][:]

        

        # Create beginnings of boundary polygon based on sts_boundary
        boundary_polygon = create_sts_boundary(base_name)
        
        os.remove(order_file)

        # Append the remaining part of the boundary polygon to be defined by
        # the user
        bounding_polygon_utm=[]
        for point in bounding_polygon:
            zone,easting,northing=redfearn(point[0],point[1])
            bounding_polygon_utm.append([easting,northing])

        boundary_polygon.append(bounding_polygon_utm[3])
        boundary_polygon.append(bounding_polygon_utm[4])

        #print 'boundary_polygon', boundary_polygon
        

        assert num.allclose(bounding_polygon_utm,boundary_polygon)


        extent_res=1000000
        meshname = 'urs_test_mesh' + '.tsh'
        interior_regions=None
        boundary_tags={'ocean': [0,1], 'otherocean': [2,3,4]}
        
        # have to change boundary tags from last example because now bounding
        # polygon starts in different place.
        create_mesh_from_regions(boundary_polygon,
                                 boundary_tags=boundary_tags,
                                 maximum_triangle_area=extent_res,
                                 filename=meshname,
                                 interior_regions=interior_regions,
                                 verbose=False)
        
        domain_fbound = Domain(meshname)
        domain_fbound.set_quantity('stage', tide)
        
        
        Bf = File_boundary(sts_file+'.sts', 
                           domain_fbound, 
                           boundary_polygon=boundary_polygon)
        time_vec = Bf.F.get_time()
        assert num.allclose(Bf.F.starttime, starttime)
        assert num.allclose(time_vec, times_ref)                                   
        
        for time_limit in [0.1, 0.2, 0.5, 1.0, 2.2, 3.0, 4.3, 6.0, 10.0]:
            Bf = File_boundary(sts_file+'.sts', 
                               domain_fbound, 
                               time_limit=time_limit+starttime,
                               boundary_polygon=boundary_polygon)
        
            time_vec = Bf.F.get_time()
            assert num.allclose(Bf.F.starttime, starttime)            
            assert num.alltrue(time_vec < time_limit)
            
            
        try:    
            Bf = File_boundary(sts_file+'.sts', 
                               domain_fbound, 
                               time_limit=-1+starttime,
                               boundary_polygon=boundary_polygon)            
            time_vec = Bf.F.get_time()    
            print(time_vec)    
        except AssertionError:
            pass
        else:
            raise Exception('Should have raised Exception here')

#-------------------------------------------------------------

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_Urs2Sts,'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)

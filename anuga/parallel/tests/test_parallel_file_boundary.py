
"""Run parallel shallow water domain.

   run using command like:

   mpirun -np m python run_parallel_sw_merimbula.py

   where m is the number of processors to be used.
   
   Will produce sww files with names domain_Pn_m.sww where m is number of processors and
   n in [0, m-1] refers to specific processor that owned this part of the partitioned mesh.
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------

import os
import sys
import time
import pypar
import numpy as num
import unittest
import tempfile
from struct import pack, unpack
from anuga.file.netcdf import NetCDFFile
import copy

#------------------------
# ANUGA Modules
#------------------------
from anuga.utilities.numerical_tools import ensure_numeric
from anuga.utilities.util_ext        import double_precision
from anuga.utilities.norms           import l1_norm, l2_norm, linf_norm
	
from anuga import Domain
from anuga import Reflective_boundary
from anuga import Dirichlet_boundary
from anuga import Time_boundary
from anuga import Transmissive_boundary
from anuga import File_boundary

from anuga.file.mux import WAVEHEIGHT_MUX2_LABEL, EAST_VELOCITY_MUX2_LABEL, \
                NORTH_VELOCITY_MUX2_LABEL
from anuga.file.mux import read_mux2_py
from anuga.file_conversion.urs2sts import urs2sts
from anuga.file.urs import Read_urs
from anuga.file.sts import create_sts_boundary

from anuga.utilities.numerical_tools import ensure_numeric
from anuga.coordinate_transforms.redfearn import redfearn
from anuga.coordinate_transforms.geo_reference import Geo_reference

from anuga import rectangular_cross_domain
from anuga.pmesh.mesh_interface import create_mesh_from_regions
from anuga import create_domain_from_file

from anuga.parallel import distribute, myid, numprocs, send, receive, barrier, finalize

from anuga.file.tests.test_mux import Test_Mux

verbose = False


class Test_urs2sts_parallel(Test_Mux):
    """ A suite of tests to test urs2sts file conversion functions.
        These tests are quite coarse-grained: converting a file
        and checking that its headers and some of its contents
        are correct.
    """ 

    def sequential_time_varying_file_boundary_sts(self):
        """sequential_ltest_time_varying_file_boundary_sts_sequential(self):
        Read correct points from ordering file and apply sts to boundary. The boundary is time varying. FIXME add to test_urs2sts.
        """
        lat_long_points=[[6.01,97.0],[6.02,97.0],[6.05,96.9],[6.0,97.0]]
        bounding_polygon=[[6.0,97.0],[6.01,97.0],[6.02,97.0],
                          [6.02,97.02],[6.00,97.02]]
        tide = 3.0
        time_step_count = 65
        time_step = 2.
        n=len(lat_long_points)
        first_tstep=num.ones(n,num.int)
        last_tstep=(time_step_count)*num.ones(n,num.int)
        finaltime=num.float(time_step*(time_step_count-1))
        yieldstep=num.float(time_step)
        gauge_depth=20*num.ones(n,num.float)
        ha=2*num.ones((n,time_step_count),num.float)
        ua=10*num.ones((n,time_step_count),num.float)
        va=-10*num.ones((n,time_step_count),num.float)

        times=num.arange(0., num.float(time_step_count*time_step), time_step)
        for i in range(n):
            #ha[i]+=num.sin(times)
            ha[i]+=times/finaltime



        sts_file="test"
        if myid==0:
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

            urs2sts(base_name,
                    basename_out=sts_file,
                    ordering_filename=order_file,
                    mean_stage=tide,
                    verbose=verbose)
            self.delete_mux(files)

            assert(os.access(sts_file+'.sts', os.F_OK))

            os.remove(order_file)

        barrier()
        boundary_polygon = create_sts_boundary(sts_file)

        # Append the remaining part of the boundary polygon to be defined by
        # the user
        bounding_polygon_utm=[]
        for point in bounding_polygon:
            zone,easting,northing=redfearn(point[0],point[1])
            bounding_polygon_utm.append([easting,northing])

        boundary_polygon.append(bounding_polygon_utm[3])
        boundary_polygon.append(bounding_polygon_utm[4])

        assert num.allclose(bounding_polygon_utm,boundary_polygon)


        extent_res=1000000
        meshname = 'urs_test_mesh' + '.tsh'
        interior_regions=None
        boundary_tags={'ocean': [0,1], 'otherocean': [2,3,4]}
        
        # have to change boundary tags from last example because now bounding
        # polygon starts in different place.
        if myid==0:
            create_mesh_from_regions(boundary_polygon,
                                     boundary_tags=boundary_tags,
                                     maximum_triangle_area=extent_res,
                                     filename=meshname,
                                     interior_regions=interior_regions,
                                     verbose=verbose)

        barrier()
        
        domain_fbound = Domain(meshname)
        domain_fbound.set_quantities_to_be_stored(None)
        domain_fbound.set_quantity('stage', tide)
        if verbose: print "Creating file boundary condition"
        Bf = File_boundary(sts_file+'.sts',
                           domain_fbound,
                           boundary_polygon=boundary_polygon)
        Br = Reflective_boundary(domain_fbound)

        domain_fbound.set_boundary({'ocean': Bf,'otherocean': Br})

        temp_fbound=num.zeros(int(finaltime/yieldstep)+1,num.float)
        if verbose: print "Evolving domain with file boundary condition"
        for i, t in enumerate(domain_fbound.evolve(yieldstep=yieldstep,
                                                   finaltime=finaltime, 
                                                   skip_initial_step = False)):
            temp_fbound[i]=domain_fbound.quantities['stage'].centroid_values[2]
            if verbose: domain_fbound.write_time()
            
        
        domain_drchlt = Domain(meshname)
        domain_drchlt.set_quantities_to_be_stored(None)
        domain_drchlt.set_starttime(time_step)
        domain_drchlt.set_quantity('stage', tide)
        Br = Reflective_boundary(domain_drchlt)
        #Bd = Dirichlet_boundary([2.0+tide,220+10*tide,-220-10*tide])
        Bd = Time_boundary(domain=domain_drchlt, f=lambda t: [2.0+t/finaltime+tide,220.+10.*tide+10.*t/finaltime,-220.-10.*tide-10.*t/finaltime])
        #Bd = Time_boundary(domain=domain_drchlt,f=lambda t: [2.0+num.sin(t)+tide,10.*(2+20.+num.sin(t)+tide),-10.*(2+20.+num.sin(t)+tide)])
        domain_drchlt.set_boundary({'ocean': Bd,'otherocean': Br})
        temp_drchlt=num.zeros(int(finaltime/yieldstep)+1,num.float)
        
        for i, t in enumerate(domain_drchlt.evolve(yieldstep=yieldstep,
                                                   finaltime=finaltime, 
                                                   skip_initial_step = False)):
            temp_drchlt[i]=domain_drchlt.quantities['stage'].centroid_values[2]
            #domain_drchlt.write_time()
        
        #print domain_fbound.quantities['stage'].vertex_values
        #print domain_drchlt.quantities['stage'].vertex_values
                    
        assert num.allclose(temp_fbound,temp_drchlt),temp_fbound-temp_drchlt

        
        assert num.allclose(domain_fbound.quantities['stage'].vertex_values,
                            domain_drchlt.quantities['stage'].vertex_values)
                        
        assert num.allclose(domain_fbound.quantities['xmomentum'].vertex_values,
                            domain_drchlt.quantities['xmomentum'].vertex_values)                        
                        
        assert num.allclose(domain_fbound.quantities['ymomentum'].vertex_values,
                            domain_drchlt.quantities['ymomentum'].vertex_values)
        
        if not sys.platform == 'win32':
            if myid==0: os.remove(sts_file+'.sts')
        
        if myid==0: os.remove(meshname)

    def parallel_time_varying_file_boundary_sts(self):
        """ parallel_test_time_varying_file_boundary_sts_sequential(self):
            Read correct points from ordering file and apply sts to boundary. 
            The boundary is time varying. Compares sequential result with 
            distributed result found using anuga_parallel
        """

        #------------------------------------------------------------
        # Define test variables
        #------------------------------------------------------------
        lat_long_points=[[6.01,97.0],[6.02,97.0],[6.05,96.9],[6.0,97.0]]
        bounding_polygon=[[6.0,97.0],[6.01,97.0],[6.02,97.0],
                          [6.02,97.02],[6.00,97.02]]
        tide = 3.0
        time_step_count = 65
        time_step = 2
        n=len(lat_long_points)
        first_tstep=num.ones(n,num.int)
        last_tstep=(time_step_count)*num.ones(n,num.int)
        finaltime=num.float(time_step*(time_step_count-1))
        yieldstep=num.float(time_step)
        gauge_depth=20*num.ones(n,num.float)
        ha=2*num.ones((n,time_step_count),num.float)
        ua=10*num.ones((n,time_step_count),num.float)
        va=-10*num.ones((n,time_step_count),num.float)

        times=num.arange(0, time_step_count*time_step, time_step)
        for i in range(n):
            #ha[i]+=num.sin(times)
            ha[i]+=times/finaltime

        #------------------------------------------------------------
        # Write mux data to file then convert to sts format
        #------------------------------------------------------------
        sts_file="test"
        if myid==0:
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

            urs2sts(base_name,
                    basename_out=sts_file,
                    ordering_filename=order_file,
                    mean_stage=tide,
                    verbose=verbose)
            self.delete_mux(files)

            assert(os.access(sts_file+'.sts', os.F_OK))

            os.remove(order_file)

        barrier()
        #------------------------------------------------------------
        # Define boundary_polygon on each processor. This polygon defines the
        # urs boundary and lies on a portion of the bounding_polygon
        #------------------------------------------------------------
        boundary_polygon = create_sts_boundary(sts_file)

        # Append the remaining part of the boundary polygon to be defined by
        # the user
        bounding_polygon_utm=[]
        for point in bounding_polygon:
            zone,easting,northing=redfearn(point[0],point[1])
            bounding_polygon_utm.append([easting,northing])

        boundary_polygon.append(bounding_polygon_utm[3])
        boundary_polygon.append(bounding_polygon_utm[4])


        assert num.allclose(bounding_polygon_utm,boundary_polygon)

        extent_res=10000
        meshname = 'urs_test_mesh' + '.tsh'
        interior_regions=None
        boundary_tags={'ocean': [0,1], 'otherocean': [2,3,4]}
        
        #------------------------------------------------------------
        # Create mesh on the master processor and store in file. This file
        # is read in by each slave processor when needed
        #------------------------------------------------------------
        if myid==0:
            create_mesh_from_regions(boundary_polygon,
                                     boundary_tags=boundary_tags,
                                     maximum_triangle_area=extent_res,
                                     filename=meshname,
                                     interior_regions=interior_regions,
                                     verbose=verbose)
        

            # barrier()
            domain_fbound = Domain(meshname)
            domain_fbound.set_quantities_to_be_stored(None)
            domain_fbound.set_quantity('stage', tide)
            # print domain_fbound.mesh.get_boundary_polygon()
        else:
            domain_fbound=None

        barrier()
        if ( verbose and myid == 0 ): 
            print 'DISTRIBUTING PARALLEL DOMAIN'
        domain_fbound = distribute(domain_fbound)

        #--------------------------------------------------------------------
        # Find which sub_domain in which the interpolation points are located 
        #
        # Sometimes the interpolation points sit exactly
        # between two centroids, so in the parallel run we
        # reset the interpolation points to the centroids
        # found in the sequential run
        #--------------------------------------------------------------------
        interpolation_points = [[279000,664000], [280250,664130], 
                                    [279280,665400], [280500,665000]]

        interpolation_points=num.array(interpolation_points)

        #if myid==0:
        #    import pylab as P
        #    boundary_polygon=num.array(boundary_polygon)
        #    P.plot(boundary_polygon[:,0],boundary_polygon[:,1])
        #    P.plot(interpolation_points[:,0],interpolation_points[:,1],'ko')
        #    P.show()

        fbound_gauge_values = []
        fbound_proc_tri_ids = []
        for i, point in enumerate(interpolation_points):
            fbound_gauge_values.append([]) # Empty list for timeseries

            try:
                k = domain_fbound.get_triangle_containing_point(point)
                if domain_fbound.tri_full_flag[k] == 1:
                    fbound_proc_tri_ids.append(k)
                else:
                    fbound_proc_tri_ids.append(-1)            
            except:
                fbound_proc_tri_ids.append(-2)


        if verbose: print 'P%d has points = %s' %(myid, fbound_proc_tri_ids)

        #------------------------------------------------------------
        # Set boundary conditions
        #------------------------------------------------------------
        Bf = File_boundary(sts_file+'.sts',
                           domain_fbound,
                           boundary_polygon=boundary_polygon)
        Br = Reflective_boundary(domain_fbound)
    
        domain_fbound.set_boundary({'ocean': Bf,'otherocean': Br})

        #------------------------------------------------------------
        # Evolve the domain on each processor
        #------------------------------------------------------------  
        for i, t in enumerate(domain_fbound.evolve(yieldstep=yieldstep,
                                                   finaltime=finaltime, 
                                                   skip_initial_step = False)):

            stage = domain_fbound.get_quantity('stage')
            for i in range(4):
                if fbound_proc_tri_ids[i] > -1:
                    fbound_gauge_values[i].append(stage.centroid_values[fbound_proc_tri_ids[i]])
        
        #------------------------------------------------------------
        # Create domain to be run sequntially on each processor
        #------------------------------------------------------------
        domain_drchlt = Domain(meshname)
        domain_drchlt.set_quantities_to_be_stored(None)
        domain_drchlt.set_starttime(time_step)
        domain_drchlt.set_quantity('stage', tide)
        Br = Reflective_boundary(domain_drchlt)
        #Bd = Dirichlet_boundary([2.0+tide,220+10*tide,-220-10*tide])
        Bd = Time_boundary(domain=domain_drchlt, function=lambda t: [2.0+t/finaltime+tide,220.+10.*tide+10.*t/finaltime,-220.-10.*tide-10.*t/finaltime])
        #Bd = Time_boundary(domain=domain_drchlt,function=lambda t: [2.0+num.sin(t)+tide,10.*(2+20.+num.sin(t)+tide),-10.*(2+20.+num.sin(t)+tide)])
        domain_drchlt.set_boundary({'ocean': Bd,'otherocean': Br})
       
        drchlt_gauge_values = []
        drchlt_proc_tri_ids = []
        for i, point in enumerate(interpolation_points):
            drchlt_gauge_values.append([]) # Empty list for timeseries

            try:
                k = domain_drchlt.get_triangle_containing_point(point)
                if domain_drchlt.tri_full_flag[k] == 1:
                    drchlt_proc_tri_ids.append(k)
                else:
                    drchlt_proc_tri_ids.append(-1)            
            except:
                drchlt_proc_tri_ids.append(-2)


        if verbose: print 'P%d has points = %s' %(myid, drchlt_proc_tri_ids)

        #------------------------------------------------------------
        # Evolve entire domain on each processor
        #------------------------------------------------------------
        for i, t in enumerate(domain_drchlt.evolve(yieldstep=yieldstep,
                                                   finaltime=finaltime, 
                                                   skip_initial_step = False)):

            stage = domain_drchlt.get_quantity('stage')
            for i in range(4):
                drchlt_gauge_values[i].append(stage.centroid_values[drchlt_proc_tri_ids[i]])

        #------------------------------------------------------------
        # Compare sequential values with parallel values
        #------------------------------------------------------------
        barrier()
        success = True
        for i in range(4):
            if fbound_proc_tri_ids[i] > -1:
                fbound_gauge_values[i]=num.array(fbound_gauge_values[i])
                drchlt_gauge_values[i]=num.array(drchlt_gauge_values[i])
                #print i,fbound_gauge_values[i][4]
                #print i,drchlt_gauge_values[i][4]
                success = success and num.allclose(fbound_gauge_values[i], drchlt_gauge_values[i])
                assert success#, (fbound_gauge_values[i]-drchlt_gauge_values[i])

        #assert_(success)       

        if not sys.platform == 'win32':
            if myid==0: os.remove(sts_file+'.sts')
        
        if myid==0: os.remove(meshname)
    
# Because we are doing assertions outside of the TestCase class
# the PyUnit defined assert_ function can't be used.
def assert_(condition, msg="Assertion Failed"):
    if condition == False:
        #pypar.finalize()
        raise AssertionError, msg


# Test an nprocs-way run of the shallow water equations
# against the sequential code.

if __name__=="__main__":
    #verbose=False
    if myid ==0 and verbose: 
        print 'PARALLEL START'
    suite = unittest.makeSuite(Test_urs2sts_parallel,'parallel_test')
    #suite = unittest.makeSuite(Test_urs2sts_parallel,'sequential_test')
    runner = unittest.TextTestRunner()
    runner.run(suite)

    #------------------------------------------
    # Run the code code and compare sequential
    # results at 4 gauge stations
    #------------------------------------------

    
    finalize()

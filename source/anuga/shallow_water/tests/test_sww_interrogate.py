import unittest
import copy
import os
import numpy as num

from anuga.coordinate_transforms.geo_reference import Geo_reference 
from anuga.geometry.polygon import is_inside_polygon
from anuga.abstract_2d_finite_volumes.util import file_function
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a
from anuga.config import g

from anuga.shallow_water.boundaries import Reflective_boundary, \
            Field_boundary, Transmissive_momentum_set_stage_boundary, \
            Transmissive_stage_zero_momentum_boundary
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions\
     import Transmissive_boundary, Dirichlet_boundary, \
            Time_boundary, File_boundary, AWI_boundary

from anuga.file.sww import get_mesh_and_quantities_from_file
            
from anuga.shallow_water.shallow_water_domain import Domain

from anuga.abstract_2d_finite_volumes.mesh_factory \
        import rectangular_cross, rectangular
        
from anuga.shallow_water.sww_interrogate import get_maximum_inundation_elevation, \
            get_maximum_inundation_location, get_maximum_inundation_data, \
            get_flow_through_cross_section, get_energy_through_cross_section
            
            
                

class Test_sww_Interrogate(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        for file in ['flowtest.sww', 'flowtest_uniquely.sww', 'runup_test_2.sww']:
            try:
                os.remove(file)
            except:
                pass     
    
    
    def test_get_maximum_inundation(self):
        """Test that sww information can be converted correctly to maximum
        runup elevation and location (without and with georeferencing)

        This test creates a slope and a runup which is maximal (~11m) at around 10s
        and levels out to the boundary condition (1m) at about 30s.
        """

        import time, os
        from anuga.file.netcdf import NetCDFFile

        #Setup

        #from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh (100m x 100m)
        points, vertices, boundary = rectangular(20, 5, 100, 50)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.default_order = 2
        domain.set_minimum_storable_height(0.01)

        filename = 'runup_test_3'
        domain.set_name(filename)
        swwfile = domain.get_name() + '.sww'

        domain.set_datadir('.')
        domain.format = 'sww'
        domain.smooth = True

        # FIXME (Ole): Backwards compatibility
        # Look at sww file and see what happens when
        # domain.tight_slope_limiters = 1
        domain.tight_slope_limiters = 0
        domain.use_centroid_velocities = 0 # Backwards compatibility (7/5/8)        
        
        Br = Reflective_boundary(domain)
        Bd = Dirichlet_boundary([1.0,0,0])


        #---------- First run without geo referencing
        
        domain.set_quantity('elevation', lambda x,y: -0.2*x + 14) # Slope
        domain.set_quantity('stage', -6)
        domain.set_boundary( {'left': Br, 'right': Bd, 'top': Br, 'bottom': Br})

        for t in domain.evolve(yieldstep=1, finaltime = 50):
            pass


        # Check maximal runup
        runup = get_maximum_inundation_elevation(swwfile)
        location = get_maximum_inundation_location(swwfile)
        #print 'Runup, location', runup, location
        assert num.allclose(runup, 6.33333333) or \
               num.allclose(runup, 6) or \
               num.allclose(runup, 12) # old limiters
        assert num.allclose(location[0], 38.33333333) or \
               num.allclose(location[0], 40.0) or \
               num.allclose(location[0], 10)

        # Check final runup
        runup = get_maximum_inundation_elevation(swwfile, time_interval=[45,50])
        location = get_maximum_inundation_location(swwfile, time_interval=[45,50])
        #print 'Runup, location:',runup, location

 
        assert num.allclose(runup, 1.666666666)
        assert num.allclose(location[0], 61.666666)

        # Check runup restricted to a polygon
        p = [[50,1], [99,1], [99,49], [50,49]]
        runup = get_maximum_inundation_elevation(swwfile, polygon=p)
        location = get_maximum_inundation_location(swwfile, polygon=p)
        #print runup, location

        assert num.allclose(runup, 3.6666666)
        assert num.allclose(location[0], 51.6666666)                

        # Check that mimimum_storable_height works
        fid = NetCDFFile(swwfile, netcdf_mode_r) # Open existing file
        
        stage = fid.variables['stage'][:]
        z = fid.variables['elevation'][:]
        xmomentum = fid.variables['xmomentum'][:]
        ymomentum = fid.variables['ymomentum'][:]        

        
        
        for i in range(stage.shape[0]):
            h = stage[i]-z # depth vector at time step i
            
            # Check every node location
            for j in range(stage.shape[1]):
                # Depth being either exactly zero implies
                # momentum being zero.
                # Or else depth must be greater than or equal to
                # the minimal storable height
                if h[j] == 0.0:
                    assert xmomentum[i,j] == 0.0
                    assert ymomentum[i,j] == 0.0                
                else:
                    assert h[j] >= domain.minimum_storable_height
        
        fid.close()

        # Cleanup
        os.remove(swwfile)
        


        #------------- Now the same with georeferencing

        domain.time=0.0
        E = 308500
        N = 6189000
        #E = N = 0
        domain.geo_reference = Geo_reference(56, E, N)

        domain.set_quantity('elevation', lambda x,y: -0.2*x + 14) # Slope
        domain.set_quantity('stage', -6)
        domain.set_boundary( {'left': Br, 'right': Bd, 'top': Br, 'bottom': Br})

        for t in domain.evolve(yieldstep=1, finaltime = 50):
            pass

        # Check maximal runup
        runup = get_maximum_inundation_elevation(swwfile)
        location = get_maximum_inundation_location(swwfile)

        #print runup, location

        assert num.allclose(runup,6.33333333) or \
               num.allclose(runup, 6) or \
               num.allclose(runup, 12) # old limiters
        assert num.allclose(location[0], 38.34+E) or \
               num.allclose(location[0], 40+E) or \
               num.allclose(location[0], 10+E)

        # Check final runup
        runup = get_maximum_inundation_elevation(swwfile, time_interval=[45,50])
        location = get_maximum_inundation_location(swwfile, time_interval=[45,50])
        #print runup, location
        #1.66666666667 [308561.66, 6189006.5]

        assert num.allclose(runup, 1.666666666)
        assert num.allclose(location[0], 61.66+E)

        # Check runup restricted to a polygon
        p = num.array([[50,1], [99,1], [99,49], [50,49]], num.int) + num.array([E, N], num.int)      #array default#

        runup = get_maximum_inundation_elevation(swwfile, polygon=p)
        location = get_maximum_inundation_location(swwfile, polygon=p)

        #print runup, location

        assert num.allclose(runup, 3.66666666)
        assert num.allclose(location[0], 51.66+E)                


        # Cleanup
        os.remove(swwfile)



    def test_get_flow_through_cross_section(self):
        """test_get_flow_through_cross_section(self):

        Test that the total flow through a cross section can be
        correctly obtained from an sww file.
        
        This test creates a flat bed with a known flow through it and tests
        that the function correctly returns the expected flow.

        The specifics are
        u = 2 m/s
        h = 1 m
        w = 3 m (width of channel)

        q = u*h*w = 6 m^3/s

        #---------- First run without geo referencing        
        
        """

        import time, os
        from anuga.file.netcdf import NetCDFFile

        # Setup
        #from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh (20m x 3m)
        width = 3
        length = 20
        t_end = 3
        points, vertices, boundary = rectangular(length, width,
                                                 length, width)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.default_order = 2
        domain.set_minimum_storable_height(0.01)

        domain.set_name('flowtest')
        swwfile = domain.get_name() + '.sww'

        domain.set_datadir('.')
        domain.format = 'sww'
        domain.smooth = True

        h = 1.0
        u = 2.0
        uh = u*h

        Br = Reflective_boundary(domain)     # Side walls
        Bd = Dirichlet_boundary([h, uh, 0])  # 2 m/s across the 3 m inlet: 


        
        domain.set_quantity('elevation', 0.0)
        domain.set_quantity('stage', h)
        domain.set_quantity('xmomentum', uh)
        domain.set_boundary( {'left': Bd, 'right': Bd, 'top': Br, 'bottom': Br})

        for t in domain.evolve(yieldstep=1, finaltime = t_end):
            pass

        # Check that momentum is as it should be in the interior

        I = [[0, width/2.],
             [length/2., width/2.],
             [length, width/2.]]
        
        f = file_function(swwfile,
                          quantities=['stage', 'xmomentum', 'ymomentum'],
                          interpolation_points=I,
                          verbose=False)
        for t in range(t_end+1):
            for i in range(3):
                assert num.allclose(f(t, i), [1, 2, 0], atol=1.0e-6)
            

        # Check flows through the middle
        for i in range(5):
            x = length/2. + i*0.23674563 # Arbitrary
            cross_section = [[x, 0], [x, width]]
            time, Q = get_flow_through_cross_section(swwfile,
                                                     cross_section,
                                                     verbose=False)

            assert num.allclose(Q, uh*width)


       
        # Try the same with partial lines
        x = length/2.
        for i in range(5):
            start_point = [length/2., i*width/5.]
            #print start_point
                            
            cross_section = [start_point, [length/2., width]]
            time, Q = get_flow_through_cross_section(swwfile,
                                                     cross_section,
                                                     verbose=False)

            #print i, Q, (width-start_point[1])
            assert num.allclose(Q, uh*(width-start_point[1]))


        # Verify no flow when line is parallel to flow
        cross_section = [[length/2.-10, width/2.], [length/2.+10, width/2.]]
        time, Q = get_flow_through_cross_section(swwfile,
                                                 cross_section,
                                                 verbose=False)

        #print i, Q
        assert num.allclose(Q, 0, atol=1.0e-5)        


        # Try with lines on an angle (all flow still runs through here)
        cross_section = [[length/2., 0], [length/2.+width, width]]
        time, Q = get_flow_through_cross_section(swwfile,
                                                 cross_section,
                                                 verbose=False)

        assert num.allclose(Q, uh*width)        
        


    def test_get_flow_through_cross_section_stored_uniquely(self):
        """test_get_flow_through_cross_section_stored_uniquely(self):

        Test that the total flow through a cross section can be
        correctly obtained from an sww file.
        
        This test creates a flat bed with a known flow through it and tests
        that the function correctly returns the expected flow.

        The specifics are
        u = 2 m/s
        h = 1 m
        w = 3 m (width of channel)

        q = u*h*w = 6 m^3/s
       
        
        """

        import time, os
        from anuga.file.netcdf import NetCDFFile

        # Setup
        #from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh (20m x 3m)
        width = 3
        length = 20
        t_end = 3
        points, vertices, boundary = rectangular(length, width,
                                                 length, width)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.default_order = 2
        domain.set_minimum_storable_height(0.01)

        domain.set_name('flowtest_uniquely')
        swwfile = domain.get_name() + '.sww'

        domain.set_store_vertices_uniquely()
        
        domain.set_datadir('.')
        domain.format = 'sww'
        domain.smooth = True

        h = 1.0
        u = 2.0
        uh = u*h

        Br = Reflective_boundary(domain)     # Side walls
        Bd = Dirichlet_boundary([h, uh, 0])  # 2 m/s across the 3 m inlet: 


        
        domain.set_quantity('elevation', 0.0)
        domain.set_quantity('stage', h)
        domain.set_quantity('xmomentum', uh)
        domain.set_boundary( {'left': Bd, 'right': Bd, 'top': Br, 'bottom': Br})

        for t in domain.evolve(yieldstep=1, finaltime = t_end):
            pass

        # Check that momentum is as it should be in the interior

        I = [[0, width/2.],
             [length/2., width/2.],
             [length, width/2.]]
        
        f = file_function(swwfile,
                          quantities=['stage', 'xmomentum', 'ymomentum'],
                          interpolation_points=I,
                          verbose=False)
        for t in range(t_end+1):
            for i in range(3):
                assert num.allclose(f(t, i), [1, 2, 0], atol=1.0e-6)
            

        # Check flows through the middle
        for i in range(5):
            x = length/2. + i*0.23674563 # Arbitrary
            cross_section = [[x, 0], [x, width]]
            time, Q = get_flow_through_cross_section(swwfile,
                                                     cross_section,
                                                     verbose=False)

            assert num.allclose(Q, uh*width)


       
        # Try the same with partial lines
        x = length/2.
        for i in range(5):
            start_point = [length/2., i*width/5.]
            #print start_point
                            
            cross_section = [start_point, [length/2., width]]
            time, Q = get_flow_through_cross_section(swwfile,
                                                     cross_section,
                                                     verbose=False)

            #print i, Q, (width-start_point[1])
            assert num.allclose(Q, uh*(width-start_point[1]))


        # Verify no flow when line is parallel to flow
        cross_section = [[length/2.-10, width/2.], [length/2.+10, width/2.]]
        time, Q = get_flow_through_cross_section(swwfile,
                                                 cross_section,
                                                 verbose=False)

        #print i, Q
        assert num.allclose(Q, 0, atol=1.0e-5)        


        # Try with lines on an angle (all flow still runs through here)
        cross_section = [[length/2., 0], [length/2.+width, width]]
        time, Q = get_flow_through_cross_section(swwfile,
                                                 cross_section,
                                                 verbose=False)

        assert num.allclose(Q, uh*width)        
        


                                      
    def test_get_flow_through_cross_section_with_geo(self):
        """test_get_flow_through_cross_section(self):

        Test that the total flow through a cross section can be
        correctly obtained from an sww file.
        
        This test creates a flat bed with a known flow through it and tests
        that the function correctly returns the expected flow.

        The specifics are
        u = 2 m/s
        h = 2 m
        w = 3 m (width of channel)

        q = u*h*w = 12 m^3/s


        This run tries it with georeferencing and with elevation = -1
        
        """

        import time, os
        from anuga.file.netcdf import NetCDFFile

        # Setup
        #from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh (20m x 3m)
        width = 3
        length = 20
        t_end = 1
        points, vertices, boundary = rectangular(length, width,
                                                 length, width)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary,
                        geo_reference = Geo_reference(56,308500,6189000))

        domain.default_order = 2
        domain.set_minimum_storable_height(0.01)

        domain.set_name('flowtest')
        swwfile = domain.get_name() + '.sww'

        domain.set_datadir('.')
        domain.format = 'sww'
        domain.smooth = True

        e = -1.0
        w = 1.0
        h = w-e
        u = 2.0
        uh = u*h

        Br = Reflective_boundary(domain)     # Side walls
        Bd = Dirichlet_boundary([w, uh, 0])  # 2 m/s across the 3 m inlet: 



        
        domain.set_quantity('elevation', e)
        domain.set_quantity('stage', w)
        domain.set_quantity('xmomentum', uh)
        domain.set_boundary( {'left': Bd, 'right': Bd, 'top': Br, 'bottom': Br})

        for t in domain.evolve(yieldstep=1, finaltime = t_end):
            pass

        # Check that momentum is as it should be in the interior

        I = [[0, width/2.],
             [length/2., width/2.],
             [length, width/2.]]
        
        I = domain.geo_reference.get_absolute(I)
        f = file_function(swwfile,
                          quantities=['stage', 'xmomentum', 'ymomentum'],
                          interpolation_points=I,
                          verbose=False)

        for t in range(t_end+1):
            for i in range(3):
                #print i, t, f(t, i)            
                assert num.allclose(f(t, i), [w, uh, 0], atol=1.0e-6)
            

        # Check flows through the middle
        for i in range(5):
            x = length/2. + i*0.23674563 # Arbitrary
            cross_section = [[x, 0], [x, width]]

            cross_section = domain.geo_reference.get_absolute(cross_section)            
            time, Q = get_flow_through_cross_section(swwfile,
                                                     cross_section,
                                                     verbose=False)

            assert num.allclose(Q, uh*width)


            
    def test_get_energy_through_cross_section(self):
        """test_get_energy_through_cross_section(self):

        Test that the specific and total energy through a cross section can be
        correctly obtained from an sww file.
        
        This test creates a flat bed with a known flow through it and tests
        that the function correctly returns the expected energies.

        The specifics are
        u = 2 m/s
        h = 1 m
        w = 3 m (width of channel)

        q = u*h*w = 6 m^3/s
        Es = h + 0.5*v*v/g  # Specific energy head [m]
        Et = w + 0.5*v*v/g  # Total energy head [m]        


        This test uses georeferencing
        
        """

        import time, os
        from anuga.file.netcdf import NetCDFFile

        # Setup
        #from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh (20m x 3m)
        width = 3
        length = 20
        t_end = 1
        points, vertices, boundary = rectangular(length, width,
                                                 length, width)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary,
                        geo_reference = Geo_reference(56,308500,6189000))

        domain.default_order = 2
        domain.set_minimum_storable_height(0.01)

        domain.set_name('flowtest')
        swwfile = domain.get_name() + '.sww'

        domain.set_datadir('.')
        domain.format = 'sww'
        domain.smooth = True

        e = -1.0
        w = 1.0
        h = w-e
        u = 2.0
        uh = u*h

        Br = Reflective_boundary(domain)     # Side walls
        Bd = Dirichlet_boundary([w, uh, 0])  # 2 m/s across the 3 m inlet: 

        
        domain.set_quantity('elevation', e)
        domain.set_quantity('stage', w)
        domain.set_quantity('xmomentum', uh)
        domain.set_boundary( {'left': Bd, 'right': Bd, 'top': Br, 'bottom': Br})

        for t in domain.evolve(yieldstep=1, finaltime = t_end):
            pass

        # Check that momentum is as it should be in the interior

        I = [[0, width/2.],
             [length/2., width/2.],
             [length, width/2.]]
        
        I = domain.geo_reference.get_absolute(I)
        f = file_function(swwfile,
                          quantities=['stage', 'xmomentum', 'ymomentum'],
                          interpolation_points=I,
                          verbose=False)

        for t in range(t_end+1):
            for i in range(3):
                #print i, t, f(t, i)
                assert num.allclose(f(t, i), [w, uh, 0], atol=1.0e-6)
            

        # Check energies through the middle
        for i in range(5):
            x = length/2. + i*0.23674563 # Arbitrary
            cross_section = [[x, 0], [x, width]]

            cross_section = domain.geo_reference.get_absolute(cross_section)            
            
            time, Es = get_energy_through_cross_section(swwfile,
                                                       cross_section,
                                                       kind='specific',
                                                       verbose=False)
            assert num.allclose(Es, h + 0.5*u*u/g)
            
            time, Et = get_energy_through_cross_section(swwfile,
                                                        cross_section,
                                                        kind='total',
                                                        verbose=False)
            assert num.allclose(Et, w + 0.5*u*u/g)            






    def test_get_maximum_inundation_from_sww(self):
        """test_get_maximum_inundation_from_sww(self)

        Test of get_maximum_inundation_elevation()
        and get_maximum_inundation_location().
   
        This is based on test_get_maximum_inundation_3(self) but works with the
        stored results instead of with the internal data structure.

        This test uses the underlying get_maximum_inundation_data for tests
        """

        verbose = False
        from anuga.config import minimum_storable_height
        
        initial_runup_height = -0.4
        final_runup_height = -0.3
        filename = 'runup_test_2'

        #--------------------------------------------------------------
        # Setup computational domain
        #--------------------------------------------------------------
        N = 10
        points, vertices, boundary = rectangular_cross(N, N)
        domain = Domain(points, vertices, boundary)
        domain.set_name(filename)
        domain.set_maximum_allowed_speed(1.0)
        #domain.set_minimum_storable_height(1.0e-5)
        domain.set_store_vertices_uniquely()

        # FIXME: This works better with old limiters so far
        domain.tight_slope_limiters = 0

        #--------------------------------------------------------------
        # Setup initial conditions
        #--------------------------------------------------------------
        def topography(x, y):
            return -x/2                             # linear bed slope

        # Use function for elevation
        domain.set_quantity('elevation', topography)
        domain.set_quantity('friction', 0.)                # Zero friction
        # Constant negative initial stage
        domain.set_quantity('stage', initial_runup_height)

        #--------------------------------------------------------------
        # Setup boundary conditions
        #--------------------------------------------------------------
        Br = Reflective_boundary(domain)                       # Reflective wall
        Bd = Dirichlet_boundary([final_runup_height, 0, 0])    # Constant inflow

        # All reflective to begin with (still water)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        #--------------------------------------------------------------
        # Test initial inundation height
        #--------------------------------------------------------------
        indices = domain.get_wet_elements()
        z = domain.get_quantity('elevation').\
                get_values(location='centroids', indices=indices)
        assert num.alltrue(z < initial_runup_height)

        q_ref = domain.get_maximum_inundation_elevation(minimum_height=minimum_storable_height)
        # First order accuracy
        assert num.allclose(q_ref, initial_runup_height, rtol=1.0/N)

        #--------------------------------------------------------------
        # Let triangles adjust
        #--------------------------------------------------------------
        q_max = None
        for t in domain.evolve(yieldstep = 0.1, finaltime = 1.0):
            q = domain.get_maximum_inundation_elevation(minimum_height=minimum_storable_height)

            if verbose:
                domain.write_time()
                print q
                
            if q > q_max:
                q_max = q

        #--------------------------------------------------------------
        # Test inundation height again
        #--------------------------------------------------------------
        #q_ref = domain.get_maximum_inundation_elevation()
        q = get_maximum_inundation_elevation(filename+'.sww')
        msg = 'We got %f, should have been %f' % (q, q_max)
        assert num.allclose(q, q_max, rtol=2.0/N), msg

        msg = 'We got %f, should have been %f' % (q, initial_runup_height)
        assert num.allclose(q, initial_runup_height, rtol = 1.0/N), msg

        # Test error condition if time interval is out
        try:
            q = get_maximum_inundation_elevation(filename+'.sww',
                                                 time_interval=[2.0, 3.0])
        except ValueError:
            pass
        else:
            msg = 'should have caught wrong time interval'
            raise Exception, msg

        # Check correct time interval
        q, loc = get_maximum_inundation_data(filename+'.sww',
                                             time_interval=[0.0, 3.0])
        msg = 'We got %f, should have been %f' % (q, initial_runup_height)
        assert num.allclose(q, initial_runup_height, rtol = 1.0/N), msg
        assert num.allclose(-loc[0]/2, q)    # From topography formula

        #--------------------------------------------------------------
        # Update boundary to allow inflow
        #--------------------------------------------------------------
        domain.set_boundary({'right': Bd})

        #--------------------------------------------------------------
        # Evolve system through time
        #--------------------------------------------------------------
        
        for t in domain.evolve(yieldstep = 0.1, finaltime = 3.0,
                               skip_initial_step = True):
            q = domain.get_maximum_inundation_elevation(minimum_height=minimum_storable_height)

            if verbose:
                domain.write_time()
                print q

            if q > q_max:
                q_max = q
                

        #--------------------------------------------------------------
        # Test inundation height again
        #--------------------------------------------------------------
        indices = domain.get_wet_elements()
        z = domain.get_quantity('elevation').\
                get_values(location='centroids', indices=indices)


        assert num.alltrue(z < final_runup_height+1.0/N)

        q = domain.get_maximum_inundation_elevation()
        # First order accuracy
        assert num.allclose(q, final_runup_height, rtol=1.0/N)

        q, loc = get_maximum_inundation_data(filename+'.sww',
                                             time_interval=[3.0, 3.0])
        msg = 'We got %f, should have been %f' % (q, final_runup_height)
        assert num.allclose(q, final_runup_height, rtol=1.0/N), msg
        assert num.allclose(-loc[0]/2, q)    # From topography formula

        q = get_maximum_inundation_elevation(filename+'.sww',verbose = verbose)
        loc = get_maximum_inundation_location(filename+'.sww')
        
        
        msg = 'We got %f, should have been %f' % (q, q_max)
        assert num.allclose(q, q_max, rtol=1.0/N), msg
        assert num.allclose(-loc[0]/2, q)    # From topography formula

        q = get_maximum_inundation_elevation(filename+'.sww',
                                             time_interval=[0, 3])
        msg = 'We got %f, should have been %f' % (q, q_max)
        assert num.allclose(q, q_max, rtol=1.0/N), msg

        # Check polygon mode
        # Runup region
        polygon = [[0.3, 0.0], [0.9, 0.0], [0.9, 1.0], [0.3, 1.0]]
        q = get_maximum_inundation_elevation(filename+'.sww',
                                             polygon = polygon,
                                             time_interval=[0, 3])
        msg = 'We got %f, should have been %f' % (q, q_max)
        assert num.allclose(q, q_max, rtol=1.0/N), msg

        # Offshore region
        polygon = [[0.9, 0.0], [1.0, 0.0], [1.0, 1.0], [0.9, 1.0]]
        q, loc = get_maximum_inundation_data(filename+'.sww',
                                             polygon = polygon,
                                             time_interval=[0, 3])
        msg = 'We got %f, should have been %f' % (q, -0.475)
        assert num.allclose(q, -0.475, rtol=1.0/N), msg
        assert is_inside_polygon(loc, polygon)
        assert num.allclose(-loc[0]/2, q)    # From topography formula

        # Dry region
        polygon = [[0.0, 0.0], [0.4, 0.0], [0.4, 1.0], [0.0, 1.0]]
        q, loc = get_maximum_inundation_data(filename+'.sww',
                                             polygon = polygon,
                                             time_interval=[0, 3])
        msg = 'We got %s, should have been None' % (q)
        assert q is None, msg
        msg = 'We got %s, should have been None' % (loc)
        assert loc is None, msg

        # Check what happens if no time point is within interval
        try:
            q = get_maximum_inundation_elevation(filename+'.sww',
                                                 time_interval=[2.75, 2.75])
        except AssertionError:
            pass
        else:
            msg = 'Time interval should have raised an exception'
            raise Exception, msg

        # Cleanup
        try:
            pass
            #os.remove(domain.get_name() + '.sww')
        except:
            pass
            #FIXME(Ole): Windows won't allow removal of this

 
 
 
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_sww_Interrogate, 'test')#_get_maximum_inundation_from_sww')
    runner = unittest.TextTestRunner() #verbosity=2)
    runner.run(suite)
               

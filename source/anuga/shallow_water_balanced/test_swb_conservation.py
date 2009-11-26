#!/usr/bin/env python

import unittest, os
import os.path
from math import pi, sqrt
import tempfile

from anuga.config import g, epsilon
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a
from anuga.utilities.numerical_tools import mean
from anuga.utilities.polygon import is_inside_polygon
from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.abstract_2d_finite_volumes.quantity import Quantity
from anuga.geospatial_data.geospatial_data import Geospatial_data
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross

from anuga.utilities.system_tools import get_pathname_from_package
from swb_domain import *

import numpy as num

# Get gateway to C implementation of flux function for direct testing
from shallow_water_ext import flux_function_central as flux_function




class Test_swb_conservation(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass


    def test_conservation_1(self):
        """Test that stage is conserved globally

        This one uses a flat bed, reflective bdries and a suitable
        initial condition
        """

        from mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order = 2

        # IC
        def x_slope(x, y):
            return x/3

        domain.set_quantity('elevation', 0)
        domain.set_quantity('friction', 0)
        domain.set_quantity('stage', x_slope)

        # Boundary conditions (reflective everywhere)
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        domain.check_integrity()

        initial_volume = domain.quantities['stage'].get_integral()
        initial_xmom = domain.quantities['xmomentum'].get_integral()

        # Evolution
        for t in domain.evolve(yieldstep=0.05, finaltime=5.0):
            volume = domain.quantities['stage'].get_integral()
            assert num.allclose(volume, initial_volume)

            #I don't believe that the total momentum should be the same
            #It starts with zero and ends with zero though
            #xmom = domain.quantities['xmomentum'].get_integral()
            #print xmom
            #assert allclose (xmom, initial_xmom)

        os.remove(domain.get_name() + '.sww')

    def test_conservation_2(self):
        """Test that stage is conserved globally

        This one uses a slopy bed, reflective bdries and a suitable
        initial condition
        """

        from mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order = 2

        # IC
        def x_slope(x, y):
            return x/3

        domain.set_quantity('elevation', x_slope)
        domain.set_quantity('friction', 0)
        domain.set_quantity('stage', 0.4)    # Steady

        # Boundary conditions (reflective everywhere)
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        domain.check_integrity()

        initial_volume = domain.quantities['stage'].get_integral()
        initial_xmom = domain.quantities['xmomentum'].get_integral()

        # Evolution
        for t in domain.evolve(yieldstep=0.05, finaltime=5.0):
            volume = domain.quantities['stage'].get_integral()
            assert num.allclose(volume, initial_volume)

            #FIXME: What would we expect from momentum
            #xmom = domain.quantities['xmomentum'].get_integral()
            #print xmom
            #assert allclose (xmom, initial_xmom)

        os.remove(domain.get_name() + '.sww')

    def test_conservation_3(self):
        """Test that stage is conserved globally

        This one uses a larger grid, convoluted bed, reflective boundaries
        and a suitable initial condition
        """

        from mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(2, 1)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False


        # IC
        def x_slope(x, y):
            z = 0*x
            for i in range(len(x)):
                if x[i] < 0.3:
                    z[i] = x[i]/3
                if 0.3 <= x[i] < 0.5:
                    z[i] = -0.5
                if 0.5 <= x[i] < 0.7:
                    z[i] = 0.39
                if 0.7 <= x[i]:
                    z[i] = x[i]/3
            return z

        domain.set_quantity('elevation', x_slope)
        domain.set_quantity('friction', 0)
        domain.set_quantity('stage', 0.4) #Steady

        # Boundary conditions (reflective everywhere)
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        domain.check_integrity()

        initial_volume = domain.quantities['stage'].get_integral()
        initial_xmom = domain.quantities['xmomentum'].get_integral()

        import copy

        ref_centroid_values = copy.copy(domain.quantities['stage'].\
                                            centroid_values)

        domain.distribute_to_vertices_and_edges()

        assert num.allclose(domain.quantities['stage'].centroid_values,
                            ref_centroid_values)

        # Check that initial limiter doesn't violate cons quan
        assert num.allclose(domain.quantities['stage'].get_integral(),
                            initial_volume)

        # Evolution
        for t in domain.evolve(yieldstep=0.05, finaltime=10):
            volume =  domain.quantities['stage'].get_integral()
                        
            assert num.allclose (volume, initial_volume)

        os.remove(domain.get_name() + '.sww')

    def test_conservation_4(self):
        """Test that stage is conserved globally

        This one uses a larger grid, convoluted bed, reflective boundaries
        and a suitable initial condition
        """

        from mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order = 2

        # IC
        def x_slope(x, y):
            z = 0*x
            for i in range(len(x)):
                if x[i] < 0.3:
                    z[i] = x[i]/3
                if 0.3 <= x[i] < 0.5:
                    z[i] = -0.5
                if 0.5 <= x[i] < 0.7:
                    #z[i] = 0.3     # OK with beta == 0.2
                    z[i] = 0.34     # OK with beta == 0.0
                    #z[i] = 0.35    # Fails after 80 timesteps with an error
                                    # of the order 1.0e-5
                if 0.7 <= x[i]:
                    z[i] = x[i]/3
            return z

        domain.set_quantity('elevation', x_slope)
        domain.set_quantity('friction', 0)
        domain.set_quantity('stage', 0.4) #Steady

        # Boundary conditions (reflective everywhere)
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        domain.check_integrity()


        initial_volume = domain.quantities['stage'].get_integral()
        initial_xmom = domain.quantities['xmomentum'].get_integral()

        import copy

        ref_centroid_values = copy.copy(domain.quantities['stage'].\
                                            centroid_values)

        # Test limiter by itself
        domain.distribute_to_vertices_and_edges()

        # Check that initial limiter doesn't violate cons quan
        assert num.allclose(domain.quantities['stage'].get_integral(),
                            initial_volume)
        # NOTE: This would fail if any initial stage was less than the
        # corresponding bed elevation - but that is reasonable.

        #Evolution
        #print domain.get_time(), initial_volume
        for t in domain.evolve(yieldstep=0.05, finaltime=10.0):
            volume =  domain.quantities['stage'].get_integral()

            #print domain.get_time(), volume
            assert num.allclose (volume, initial_volume)

        os.remove(domain.get_name() + '.sww')

    def test_conservation_5(self):
        """Test that momentum is conserved globally in steady state scenario

        This one uses a slopy bed, dirichlet and reflective bdries
        """

        from mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order = 2

        # IC
        def x_slope(x, y):
            return x/3

        domain.set_quantity('elevation', x_slope)
        domain.set_quantity('friction', 0)
        domain.set_quantity('stage', 0.4) # Steady

        # Boundary conditions (reflective everywhere)
        Br = Reflective_boundary(domain)
        Bleft = Dirichlet_boundary([0.5, 0, 0])
        Bright = Dirichlet_boundary([0.1, 0, 0])
        domain.set_boundary({'left': Bleft, 'right': Bright,
                             'top': Br, 'bottom': Br})

        domain.check_integrity()

        initial_volume = domain.quantities['stage'].get_integral()
        initial_xmom = domain.quantities['xmomentum'].get_integral()

        # Evolution
        for t in domain.evolve(yieldstep=0.05, finaltime=15.0):
            stage =  domain.quantities['stage'].get_integral()
            xmom = domain.quantities['xmomentum'].get_integral()
            ymom = domain.quantities['ymomentum'].get_integral()

            if num.allclose(t, 10):    # Steady state reached
                steady_xmom = domain.quantities['xmomentum'].get_integral()
                steady_ymom = domain.quantities['ymomentum'].get_integral()
                steady_stage = domain.quantities['stage'].get_integral()

            if t > 10:
                msg = 'time=%.2f, xmom=%.10f, steady_xmom=%.10f' % (t, xmom, steady_xmom)
                assert num.allclose(xmom, steady_xmom,atol=1.0e-4), msg

                msg = 'time=%.2f, ymom=%.10f, steady_ymom=%.10f' % (t, ymom, steady_ymom)
                assert num.allclose(ymom, steady_ymom,atol=1.0e-4), msg

                msg = 'time=%.2f, stage=%.10f, steady_stage=%.10f' % (t, stage, steady_stage)
                assert num.allclose(stage, steady_stage,atol=1.0e-4)

        os.remove(domain.get_name() + '.sww')

    def test_conservation_real(self):
        """Test that momentum is conserved globally

        Stephen finally made a test that revealed the problem.
        This test failed with code prior to 25 July 2005
        """

        import sys
        import os.path
        sys.path.append(os.path.join('..', 'abstract_2d_finite_volumes'))
        from mesh_factory import rectangular

        yieldstep = 0.01
        finaltime = 0.05
        min_depth = 1.0e-2

        #Create shallow water domain
        points, vertices, boundary = rectangular(10, 10, len1=500, len2=500)
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order = 1
        domain.minimum_allowed_height = min_depth

        # Set initial condition
        class Set_IC:
            """Set an initial condition with a constant value, for x0<x<x1"""

            def __init__(self, x0=0.25, x1=0.5, h=1.0):
                self.x0 = x0
                self.x1 = x1
                self.h  = h

            def __call__(self, x, y):
                return self.h*((x > self.x0) & (x < self.x1))

        domain.set_quantity('stage', Set_IC(200.0, 300.0, 5.0))

        # Boundaries
        R = Reflective_boundary(domain)
        domain.set_boundary({'left': R, 'right': R, 'top':R, 'bottom': R})

        ref = domain.quantities['stage'].get_integral()

        # Evolution
        for t in domain.evolve(yieldstep=yieldstep, finaltime=finaltime):
            pass

        now = domain.quantities['stage'].get_integral()

        msg = 'Stage not conserved: was %f, now %f' % (ref, now)
        assert num.allclose(ref, now), msg

        os.remove(domain.get_name() + '.sww')


    def test_total_volume(self):        
        """test_total_volume
        
        Test that total volume can be computed correctly
        """            

        #----------------------------------------------------------------------
        # Import necessary modules
        #----------------------------------------------------------------------
        from anuga.abstract_2d_finite_volumes.mesh_factory \
                import rectangular_cross
        from anuga.shallow_water import Domain

        #----------------------------------------------------------------------
        # Setup computational domain
        #----------------------------------------------------------------------

        length = 100.
        width  = 20.
        dx = dy = 5       # Resolution: of grid on both axes
        
        points, vertices, boundary = rectangular_cross(int(length/dx),
                                                       int(width/dy),
                                                       len1=length,
                                                       len2=width)
        domain = Domain(points, vertices, boundary)   

        #----------------------------------------------------------------------
        # Simple flat bottom bathtub
        #----------------------------------------------------------------------

        d = 1.0
        domain.set_quantity('elevation', 0.0)
        domain.set_quantity('stage', d)
        
        assert num.allclose(domain.compute_total_volume(), length*width*d)

        #----------------------------------------------------------------------
        # Slope
        #----------------------------------------------------------------------
                
        slope = 1.0/10          # RHS drops to -10m
        def topography(x, y):
            return -x * slope

        domain.set_quantity('elevation', topography)
        domain.set_quantity('stage', 0.0)       # Domain full
        
        V = domain.compute_total_volume()
        assert num.allclose(V, length*width*10/2)

        domain.set_quantity('stage', -5.0)      # Domain 'half' full
        
        # IMPORTANT: Adjust stage to match elevation
        domain.distribute_to_vertices_and_edges()
        
        V = domain.compute_total_volume()
        assert num.allclose(V, width*(length/2)*5.0/2)


    def test_volumetric_balance_computation(self):
        """test_volumetric_balance_computation
        
        Test that total in and out flows are computed correctly 
        in a steady state situation
        """

        # Set to True if volumetric output is sought
        verbose = False

        #----------------------------------------------------------------------
        # Import necessary modules
        #----------------------------------------------------------------------

        from anuga.abstract_2d_finite_volumes.mesh_factory \
                import rectangular_cross
        from anuga.shallow_water import Domain
        from anuga.shallow_water.shallow_water_domain import Reflective_boundary
        from anuga.shallow_water.shallow_water_domain import Dirichlet_boundary
        from anuga.shallow_water.shallow_water_domain import Inflow
        from anuga.shallow_water.data_manager \
                import get_flow_through_cross_section

        #----------------------------------------------------------------------
        # Setup computational domain
        #----------------------------------------------------------------------

        finaltime = 500.0
        length = 300.
        width  = 20.
        dx = dy = 5       # Resolution: of grid on both axes
        
        # Input parameters
        uh = 1.0
        vh = 0.0
        d = 1.0
        
        # 20 m^3/s in the x direction across entire domain
        ref_flow = uh*d*width

        points, vertices, boundary = rectangular_cross(int(length/dx),
                                                       int(width/dy),
                                                       len1=length,
                                                       len2=width)

        domain = Domain(points, vertices, boundary)   
        domain.set_name('Inflow_flowline_test')              # Output name

        #----------------------------------------------------------------------
        # Setup initial conditions
        #----------------------------------------------------------------------

        domain.set_quantity('elevation', 0.0)  # Flat bed
        domain.set_quantity('friction', 0.0)   # Constant zero friction
                
        domain.set_quantity('stage', expression='elevation + %d' % d) 

        #----------------------------------------------------------------------
        # Setup boundary conditions
        #----------------------------------------------------------------------

        Br = Reflective_boundary(domain)      # Solid reflective wall
                
        # Constant flow in and out of domain
        # Depth = 1m, uh=1 m/s, i.e. a flow of 20 m^3/s 
        Bi = Dirichlet_boundary([d, uh, vh]) 
        Bo = Dirichlet_boundary([d, uh, vh])

        domain.set_boundary({'left': Bi, 'right': Bo, 'top': Br, 'bottom': Br})

        #----------------------------------------------------------------------
        # Evolve system through time
        #----------------------------------------------------------------------

        for t in domain.evolve(yieldstep=50.0, finaltime=finaltime):
            S = domain.volumetric_balance_statistics()
            if verbose :
                print domain.timestepping_statistics()
                print S
                
            if t > 300:
                # Steady state reached
                
                # Square on flowline at 200m
                q = domain.get_flow_through_cross_section([[200.0,  0.0],
                                                           [200.0, 20.0]])
                
                assert num.allclose(q, ref_flow)

        os.remove('Inflow_flowline_test.sww') 

    def test_volume_conservation_inflow(self):
        """test_volume_conservation
        
        Test that total volume in domain is as expected, based on questions
        raised by Petar Milevski in May 2009.
        
        This test adds inflow at a known rate and verifies that the total 
        terminal volume is as expected.
        
        """

        verbose = False
        

        #---------------------------------------------------------------------
        # Import necessary modules
        #---------------------------------------------------------------------
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
        from anuga.shallow_water import Domain
        from anuga.shallow_water.shallow_water_domain import Reflective_boundary
        from anuga.shallow_water.shallow_water_domain import Dirichlet_boundary
        from anuga.shallow_water.shallow_water_domain import Inflow
        from anuga.shallow_water.data_manager import get_flow_through_cross_section

        #----------------------------------------------------------------------
        # Setup computational domain
        #----------------------------------------------------------------------
        finaltime = 200.0

        length = 300.
        width  = 20.
        dx = dy = 5       # Resolution: of grid on both axes
        

        points, vertices, boundary = rectangular_cross(int(length/dx), 
                                                       int(width/dy),
                                                       len1=length, len2=width)


        domain = Domain(points, vertices, boundary)   
        domain.set_name('Inflow_volume_test')              # Output name
                

        #----------------------------------------------------------------------
        # Setup initial conditions
        #----------------------------------------------------------------------
        slope = 0.0
        def topography(x, y):
            z=-x * slope
            return z

        domain.set_quantity('elevation', topography) # Use function for elevation
        domain.set_quantity('friction', 0.0)         # Constant friction
                
        domain.set_quantity('stage',
                            expression='elevation')  # Dry initially
                            

        #--------------------------------------------------------------
        # Setup Inflow
        #--------------------------------------------------------------

        # Fixed Flowrate onto Area 
        fixed_inflow = Inflow(domain,
                              center=(10.0, 10.0),
                              radius=5.00,
                              rate=10.00)                               
                            
        domain.forcing_terms.append(fixed_inflow)                            
        
        #----------------------------------------------------------------------
        # Setup boundary conditions
        #----------------------------------------------------------------------

        Br = Reflective_boundary(domain) # Solid reflective wall
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        
        #----------------------------------------------------------------------
        # Evolve system through time
        #----------------------------------------------------------------------
        ref_volume = 0.0
        ys = 10.0  # Yieldstep
        for t in domain.evolve(yieldstep=ys, finaltime=finaltime):
        
            # Check volume
            assert num.allclose(domain.compute_total_volume(), ref_volume)
        
            if verbose :
                print domain.timestepping_statistics()
                print domain.volumetric_balance_statistics()
                print 'reference volume', ref_volume
            
            
            # Update reference volume
            ref_volume += ys * fixed_inflow.rate


        os.remove('Inflow_volume_test.sww')


        
    def test_volume_conservation_rain(self):
        """test_volume_conservation
        
        Test that total volume in domain is as expected, based on questions
        raised by Petar Milevski in May 2009.
        
        This test adds rain at a known rate and verifies that the total 
        terminal volume is as expected.
        
        """

        verbose = False
        

        #---------------------------------------------------------------------
        # Import necessary modules
        #---------------------------------------------------------------------
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
        from anuga.shallow_water import Domain
        from anuga.shallow_water.shallow_water_domain import Reflective_boundary
        from anuga.shallow_water.shallow_water_domain import Dirichlet_boundary
        from anuga.shallow_water.shallow_water_domain import Rainfall
        from anuga.shallow_water.data_manager import get_flow_through_cross_section

        #----------------------------------------------------------------------
        # Setup computational domain
        #----------------------------------------------------------------------
        finaltime = 200.0

        length = 300.
        width  = 20.
        dx = dy = 5       # Resolution: of grid on both axes
        

        points, vertices, boundary = rectangular_cross(int(length/dx), 
                                                       int(width/dy),
                                                       len1=length, len2=width)


        domain = Domain(points, vertices, boundary)   
        domain.set_name('Rain_volume_test')              # Output name
                

        #----------------------------------------------------------------------
        # Setup initial conditions
        #----------------------------------------------------------------------
        slope = 0.0
        def topography(x, y):
            z=-x * slope
            return z

        domain.set_quantity('elevation', topography) # Use function for elevation
        domain.set_quantity('friction', 0.0)         # Constant friction
                
        domain.set_quantity('stage',
                            expression='elevation')  # Dry initially
                            

        #--------------------------------------------------------------
        # Setup rain
        #--------------------------------------------------------------

        # Fixed rain onto small circular area 
        fixed_rain = Rainfall(domain,
                              center=(10.0, 10.0),
                              radius=5.00,
                              rate=10.00)   # 10 mm/s                            
                            
        domain.forcing_terms.append(fixed_rain)                            
        
        #----------------------------------------------------------------------
        # Setup boundary conditions
        #----------------------------------------------------------------------

        Br = Reflective_boundary(domain) # Solid reflective wall
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        
        #----------------------------------------------------------------------
        # Evolve system through time
        #----------------------------------------------------------------------
        ref_volume = 0.0
        ys = 10.0  # Yieldstep
        for t in domain.evolve(yieldstep=ys, finaltime=finaltime):
        
            # Check volume
            V = domain.compute_total_volume()
            msg = 'V = %e, Ref = %e' % (V, ref_volume)
            assert num.allclose(V, ref_volume), msg
        
            if verbose :
                print domain.timestepping_statistics()
                print domain.volumetric_balance_statistics()
                print 'reference volume', ref_volume
                print V
            
            
            # Update reference volume.
            # FIXME: Note that rate has now been redefined
            # as m/s internally. This is a little confusing
            # when it was specfied as mm/s.
            
            delta_V = fixed_rain.rate*fixed_rain.exchange_area
            ref_volume += ys * delta_V

        os.remove('Rain_volume_test.sww')

    def Xtest_rain_conservation_and_runoff(self):
        """test_rain_conservation_and_runoff
        
        Test that total volume in domain is as expected, based on questions
        raised by Petar Milevski in May 2009.
        
        This test adds rain at a known rate and verifies that the total 
        volume and outflows are as expected.
        
        """

        # FIXME (Ole): Does not work yet. Investigate boundary flows
        
        verbose = True #False
        

        #---------------------------------------------------------------------
        # Import necessary modules
        #---------------------------------------------------------------------
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
        from anuga.shallow_water import Domain
        from anuga.shallow_water.shallow_water_domain import Reflective_boundary
        from anuga.shallow_water.shallow_water_domain import Dirichlet_boundary
        from anuga.shallow_water.shallow_water_domain import Rainfall
        from anuga.shallow_water.data_manager import get_flow_through_cross_section

        #----------------------------------------------------------------------
        # Setup computational domain
        #----------------------------------------------------------------------
        finaltime = 500.0

        length = 300.
        width  = 20.
        dx = dy = 5       # Resolution: of grid on both axes
        

        points, vertices, boundary = rectangular_cross(int(length/dx), 
                                                       int(width/dy),
                                                       len1=length, len2=width)


        domain = Domain(points, vertices, boundary)   
        domain.set_name('Rain_volume_runoff_test')         # Output name
                

        #----------------------------------------------------------------------
        # Setup initial conditions
        #----------------------------------------------------------------------
        slope = 0.0
        def topography(x, y):
            z=-x * slope
            return z

        domain.set_quantity('elevation', topography) # Use function for elevation
        domain.set_quantity('friction', 0.0)         # Constant friction
                
        domain.set_quantity('stage',
                            expression='elevation')  # Dry initially
                            

        #--------------------------------------------------------------
        # Setup rain
        #--------------------------------------------------------------

        # Fixed rain onto small circular area 
        fixed_rain = Rainfall(domain,
                              center=(10.0, 10.0),
                              radius=5.00,
                              rate=10.00)   # 10 mm/s                            
                            
        domain.forcing_terms.append(fixed_rain)                            
        
        #----------------------------------------------------------------------
        # Setup boundary conditions
        #----------------------------------------------------------------------

        Br = Reflective_boundary(domain) # Solid reflective wall
        Bt = Transmissive_stage_zero_momentum_boundary(domain)
        Bd = Dirichlet_boundary([-10, 0, 0])
        domain.set_boundary({'left': Bt, 'right': Bd, 'top': Bt, 'bottom': Bt})

        
        #----------------------------------------------------------------------
        # Evolve system through time
        #----------------------------------------------------------------------
        ref_volume = 0.0
        ys = 10.0  # Yieldstep
        for t in domain.evolve(yieldstep=ys, finaltime=finaltime):
        
            # Check volume
            V = domain.compute_total_volume()
            msg = 'V = %e, Ref = %e' % (V, ref_volume)
            #assert num.allclose(V, ref_volume) or V < ref_volume, msg
        
            if verbose:
                print domain.timestepping_statistics()
                print domain.volumetric_balance_statistics()
                print 'reference volume', ref_volume
                print V
            
            
            # Update reference volume.
            # FIXME: Note that rate has now been redefined
            # as m/s internally. This is a little confusing
            # when it was specfied as mm/s.
            
            delta_V = fixed_rain.rate*fixed_rain.exchange_area
            ref_volume += ys * delta_V
        
            # Compute outflow at right hand downstream boundary
            boundary_flows, inflow , outflow = domain.compute_boundary_flows()
            net_outflow = outflow - inflow
        
            outflow = boundary_flows['right']
            if verbose:
                print 'Outflow', outflow
                print 'Net outflow', net_outflow
        
            # Update reference volume
            ref_volume += ys * outflow            



#################################################################################

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_swb_conservation, 'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)

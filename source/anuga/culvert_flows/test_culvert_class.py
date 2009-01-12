#!/usr/bin/env python


import unittest
import os.path
import sys

from anuga.utilities.system_tools import get_pathname_from_package
from anuga.utilities.polygon import Polygon_function
        
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
from anuga.abstract_2d_finite_volumes.quantity import Quantity

from anuga.shallow_water import Domain, Reflective_boundary,\
    Dirichlet_boundary,\
    Transmissive_boundary, Time_boundary

from anuga.culvert_flows.culvert_class import Culvert_flow, Culvert_flow_rating, Culvert_flow_energy
from anuga.culvert_flows.culvert_routines import boyd_generalised_culvert_model
     
from math import pi,pow,sqrt
from Numeric import choose, greater, ones, sin, cos, exp, cosh, allclose



class Test_Culvert(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass


    def test_that_culvert_runs_rating(self):
        """test_that_culvert_runs_rating
        
        This test exercises the culvert and checks values outside rating curve
        are dealt with       
        """

        path = get_pathname_from_package('anuga.culvert_flows')    
        
        length = 40.
        width = 5.

        dx = dy = 1           # Resolution: Length of subdivisions on both axes

        points, vertices, boundary = rectangular_cross(int(length/dx),
                                                       int(width/dy),
                                                       len1=length, 
                                                       len2=width)
        domain = Domain(points, vertices, boundary)   
        domain.set_name('Test_culvert')                 # Output name
        domain.set_default_order(2)


        #----------------------------------------------------------------------
        # Setup initial conditions
        #----------------------------------------------------------------------

        def topography(x, y):
            """Set up a weir
            
            A culvert will connect either side
            """
            # General Slope of Topography
            z=-x/1000
            
            N = len(x)
            for i in range(N):

               # Sloping Embankment Across Channel
                if 5.0 < x[i] < 10.1:
                    # Cut Out Segment for Culvert FACE                
                    if  1.0+(x[i]-5.0)/5.0 <  y[i]  < 4.0 - (x[i]-5.0)/5.0: 
                       z[i]=z[i]
                    else:
                       z[i] +=  0.5*(x[i] -5.0)    # Sloping Segment  U/S Face
                if 10.0 < x[i] < 12.1:
                   z[i] +=  2.5                    # Flat Crest of Embankment
                if 12.0 < x[i] < 14.5:
                    # Cut Out Segment for Culvert FACE                
                    if  2.0-(x[i]-12.0)/2.5 <  y[i]  < 3.0 + (x[i]-12.0)/2.5:
                       z[i]=z[i]
                    else:
                       z[i] +=  2.5-1.0*(x[i] -12.0) # Sloping D/S Face
                		   
        		
            return z


        domain.set_quantity('elevation', topography) 
        domain.set_quantity('friction', 0.01)         # Constant friction 
        domain.set_quantity('stage',
                            expression='elevation')   # Dry initial condition

        filename=os.path.join(path, 'example_rating_curve.csv')
        culvert = Culvert_flow(domain,
                               culvert_description_filename=filename,        
                               end_point0=[9.0, 2.5], 
                               end_point1=[13.0, 2.5],
                               use_velocity_head=True,
                               verbose=False)

        domain.forcing_terms.append(culvert)
        

        #-----------------------------------------------------------------------
        # Setup boundary conditions
        #-----------------------------------------------------------------------

        # Inflow based on Flow Depth and Approaching Momentum
        Bi = Dirichlet_boundary([0.0, 0.0, 0.0])
        Br = Reflective_boundary(domain)              # Solid reflective wall
        Bo = Dirichlet_boundary([-5, 0, 0])           # Outflow
        
        # Upstream and downstream conditions that will exceed the rating curve
        # I.e produce delta_h outside the range [0, 10] specified in the the 
        # file example_rating_curve.csv
        Btus = Time_boundary(domain, lambda t: [100*sin(2*pi*(t-4)/10), 0.0, 0.0])
        Btds = Time_boundary(domain, lambda t: [-5*(cos(2*pi*(t-4)/20)), 0.0, 0.0])
        domain.set_boundary({'left': Btus, 'right': Btds, 'top': Br, 'bottom': Br})


        #-----------------------------------------------------------------------
        # Evolve system through time
        #-----------------------------------------------------------------------

        min_delta_w = sys.maxint 
        max_delta_w = -min_delta_w
        for t in domain.evolve(yieldstep = 1, finaltime = 25):
            delta_w = culvert.inlet.stage - culvert.outlet.stage
            
            if delta_w > max_delta_w: max_delta_w = delta_w
            if delta_w < min_delta_w: min_delta_w = delta_w            
            
            #print domain.timestepping_statistics()
            pass

        # Check that extreme values in rating curve have been exceeded
        # so that we know that condition has been exercised
        assert min_delta_w < 0
        assert max_delta_w > 10        
        

    def test_that_culvert_dry_bed_rating_does_not_produce_flow(self):
        """test_that_culvert_in_dry_bed_does_not_produce_flow(self):
        
        Test that culvert on a sloping dry bed doesn't produce flows
        although there will be a 'pressure' head due to delta_w > 0

        This one is using the rating curve variant
        """

        path = get_pathname_from_package('anuga.culvert_flows')    
        
        length = 40.
        width = 5.

        dx = dy = 1           # Resolution: Length of subdivisions on both axes

        points, vertices, boundary = rectangular_cross(int(length/dx),
                                                       int(width/dy),
                                                       len1=length, 
                                                       len2=width)
        domain = Domain(points, vertices, boundary)   
        domain.set_name('Test_culvert_dry') # Output name
        domain.set_default_order(2)


        #----------------------------------------------------------------------
        # Setup initial conditions
        #----------------------------------------------------------------------

        def topography(x, y):
            """Set up a weir
            
            A culvert will connect either side
            """
            # General Slope of Topography
            z=-x/1000
            
            N = len(x)
            for i in range(N):

               # Sloping Embankment Across Channel
                if 5.0 < x[i] < 10.1:
                    # Cut Out Segment for Culvert FACE                
                    if  1.0+(x[i]-5.0)/5.0 <  y[i]  < 4.0 - (x[i]-5.0)/5.0: 
                       z[i]=z[i]
                    else:
                       z[i] +=  0.5*(x[i] -5.0)    # Sloping Segment  U/S Face
                if 10.0 < x[i] < 12.1:
                   z[i] +=  2.5                    # Flat Crest of Embankment
                if 12.0 < x[i] < 14.5:
                    # Cut Out Segment for Culvert FACE                
                    if  2.0-(x[i]-12.0)/2.5 <  y[i]  < 3.0 + (x[i]-12.0)/2.5:
                       z[i]=z[i]
                    else:
                       z[i] +=  2.5-1.0*(x[i] -12.0) # Sloping D/S Face
                		   
        		
            return z


        domain.set_quantity('elevation', topography) 
        domain.set_quantity('friction', 0.01)         # Constant friction 
        domain.set_quantity('stage',
                            expression='elevation')   # Dry initial condition


        filename = os.path.join(path, 'example_rating_curve.csv')
        culvert = Culvert_flow(domain,
                               culvert_description_filename=filename,
                               end_point0=[9.0, 2.5], 
                               end_point1=[13.0, 2.5],
                               verbose=False)

        domain.forcing_terms.append(culvert)
        

        #-----------------------------------------------------------------------
        # Setup boundary conditions
        #-----------------------------------------------------------------------

        # Inflow based on Flow Depth and Approaching Momentum

        Br = Reflective_boundary(domain)              # Solid reflective wall
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})


        #-----------------------------------------------------------------------
        # Evolve system through time
        #-----------------------------------------------------------------------

        ref_volume = domain.get_quantity('stage').get_integral()
        for t in domain.evolve(yieldstep = 1, finaltime = 25):
            #print domain.timestepping_statistics()
            new_volume = domain.get_quantity('stage').get_integral()
            
            msg = 'Total volume has changed'
            assert allclose(new_volume, ref_volume), msg
            pass
    

    
    def test_that_culvert_rating_limits_flow_in_shallow_inlet_condition(self):
        """test_that_culvert_rating_limits_flow_in_shallow_inlet_condition
        
        Test that culvert on a sloping dry bed limits flows when very little water
        is present at inlet

        This one is using the rating curve variant
        """

        path = get_pathname_from_package('anuga.culvert_flows')    
        
        length = 40.
        width = 5.

        dx = dy = 1           # Resolution: Length of subdivisions on both axes

        points, vertices, boundary = rectangular_cross(int(length/dx),
                                                       int(width/dy),
                                                       len1=length, 
                                                       len2=width)
        domain = Domain(points, vertices, boundary)   
        domain.set_name('Test_culvert_shallow') # Output name
        domain.set_default_order(2)


        #----------------------------------------------------------------------
        # Setup initial conditions
        #----------------------------------------------------------------------

        def topography(x, y):
            """Set up a weir
            
            A culvert will connect either side
            """
            # General Slope of Topography
            z=-x/1000
            
            N = len(x)
            for i in range(N):

               # Sloping Embankment Across Channel
                if 5.0 < x[i] < 10.1:
                    # Cut Out Segment for Culvert FACE                
                    if  1.0+(x[i]-5.0)/5.0 <  y[i]  < 4.0 - (x[i]-5.0)/5.0: 
                       z[i]=z[i]
                    else:
                       z[i] +=  0.5*(x[i] -5.0)    # Sloping Segment  U/S Face
                if 10.0 < x[i] < 12.1:
                   z[i] +=  2.5                    # Flat Crest of Embankment
                if 12.0 < x[i] < 14.5:
                    # Cut Out Segment for Culvert FACE                
                    if  2.0-(x[i]-12.0)/2.5 <  y[i]  < 3.0 + (x[i]-12.0)/2.5:
                       z[i]=z[i]
                    else:
                       z[i] +=  2.5-1.0*(x[i] -12.0) # Sloping D/S Face
                		   
        		
            return z


        domain.set_quantity('elevation', topography) 
        domain.set_quantity('friction', 0.01)         # Constant friction 
        domain.set_quantity('stage',
                            expression='elevation + 0.2') # Shallow initial condition
                            
        # NOTE: Shallow values may cause this test to fail regardless of the
        # culvert due to initial adjustments. A good value is 0.2


        filename = os.path.join(path, 'example_rating_curve.csv')
        culvert = Culvert_flow(domain,
                               culvert_description_filename=filename,
                               end_point0=[9.0, 2.5], 
                               end_point1=[13.0, 2.5],
                               trigger_depth=0.01,
                               verbose=False)

        domain.forcing_terms.append(culvert)
        

        #-----------------------------------------------------------------------
        # Setup boundary conditions
        #-----------------------------------------------------------------------

        # Inflow based on Flow Depth and Approaching Momentum

        Br = Reflective_boundary(domain)              # Solid reflective wall
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})


        #-----------------------------------------------------------------------
        # Evolve system through time
        #-----------------------------------------------------------------------

        ref_volume = domain.get_quantity('stage').get_integral()
        for t in domain.evolve(yieldstep = 0.1, finaltime = 25):
            #print domain.timestepping_statistics()
            new_volume = domain.get_quantity('stage').get_integral()
            
            msg = 'Total volume has changed: Is %.2f m^3 should have been %.2f m^3'\
                % (new_volume, ref_volume)
            if not allclose(new_volume, ref_volume):
                print msg
            assert allclose(new_volume, ref_volume), msg
    
    
    
    def test_that_culvert_dry_bed_boyd_does_not_produce_flow(self):
        """test_that_culvert_in_dry_bed_boyd_does_not_produce_flow(self):
        
        Test that culvert on a sloping dry bed doesn't produce flows
        although there will be a 'pressure' head due to delta_w > 0

        This one is using the 'Boyd' variant        
        """

        path = get_pathname_from_package('anuga.culvert_flows')    
        
        length = 40.
        width = 5.

        dx = dy = 1           # Resolution: Length of subdivisions on both axes

        points, vertices, boundary = rectangular_cross(int(length/dx),
                                                       int(width/dy),
                                                       len1=length, 
                                                       len2=width)
        domain = Domain(points, vertices, boundary)   
        domain.set_name('Test_culvert_dry') # Output name
        domain.set_default_order(2)


        #----------------------------------------------------------------------
        # Setup initial conditions
        #----------------------------------------------------------------------

        def topography(x, y):
            """Set up a weir
            
            A culvert will connect either side
            """
            # General Slope of Topography
            z=-x/1000
            
            N = len(x)
            for i in range(N):

               # Sloping Embankment Across Channel
                if 5.0 < x[i] < 10.1:
                    # Cut Out Segment for Culvert FACE                
                    if  1.0+(x[i]-5.0)/5.0 <  y[i]  < 4.0 - (x[i]-5.0)/5.0: 
                       z[i]=z[i]
                    else:
                       z[i] +=  0.5*(x[i] -5.0)    # Sloping Segment  U/S Face
                if 10.0 < x[i] < 12.1:
                   z[i] +=  2.5                    # Flat Crest of Embankment
                if 12.0 < x[i] < 14.5:
                    # Cut Out Segment for Culvert FACE                
                    if  2.0-(x[i]-12.0)/2.5 <  y[i]  < 3.0 + (x[i]-12.0)/2.5:
                       z[i]=z[i]
                    else:
                       z[i] +=  2.5-1.0*(x[i] -12.0) # Sloping D/S Face
                		   
        		
            return z


        domain.set_quantity('elevation', topography) 
        domain.set_quantity('friction', 0.01)         # Constant friction 
        domain.set_quantity('stage',
                            expression='elevation')   # Dry initial condition


        filename = os.path.join(path, 'example_rating_curve.csv')


        culvert = Culvert_flow(domain,
                               label='Culvert No. 1',
                               description='This culvert is a test unit 1.2m Wide by 0.75m High',   
                               end_point0=[9.0, 2.5], 
                               end_point1=[13.0, 2.5],
                               width=1.20,height=0.75,
                               culvert_routine=boyd_generalised_culvert_model,
                               number_of_barrels=1,
                               update_interval=2,
                               verbose=True)
        
        domain.forcing_terms.append(culvert)
        

        #-----------------------------------------------------------------------
        # Setup boundary conditions
        #-----------------------------------------------------------------------

        # Inflow based on Flow Depth and Approaching Momentum

        Br = Reflective_boundary(domain)              # Solid reflective wall
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})


        #-----------------------------------------------------------------------
        # Evolve system through time
        #-----------------------------------------------------------------------

        ref_volume = domain.get_quantity('stage').get_integral()
        for t in domain.evolve(yieldstep = 1, finaltime = 25):
            #print domain.timestepping_statistics()
            new_volume = domain.get_quantity('stage').get_integral()
            
            msg = 'Total volume has changed'
            assert allclose(new_volume, ref_volume), msg
            pass
    

    
    

    def test_predicted_boyd_flow(self):
        """test_predicted_boyd_flow
        
        Test that flows predicted by the boyd method are consistent with what what
        is calculated in engineering codes.
        The data was supplied by Petar Milevski
        """

        # FIXME(Ole) this is nowhere near finished
        path = get_pathname_from_package('anuga.culvert_flows')    
        
        length = 12.
        width = 5.

        dx = dy = 0.5           # Resolution: Length of subdivisions on both axes

        points, vertices, boundary = rectangular_cross(int(length/dx),
                                                       int(width/dy),
                                                       len1=length, 
                                                       len2=width)
        domain = Domain(points, vertices, boundary)   
        
        domain.set_name('test_culvert')                 # Output name
        domain.set_default_order(2)


        #----------------------------------------------------------------------
        # Setup initial conditions
        #----------------------------------------------------------------------

        def topography(x, y):
            # General Slope of Topography
            z=-x/10
            
            return z


        domain.set_quantity('elevation', topography) 
        domain.set_quantity('friction', 0.01)         # Constant friction 
        domain.set_quantity('stage', expression='elevation')


        Q0 = domain.get_quantity('stage')
        Q1 = Quantity(domain)
        
        # Add depths to stage        
        head_water_depth = 0.169
        tail_water_depth = 0.089
        
        inlet_poly = [[0,0], [6,0], [6,5], [0,5]]
        outlet_poly = [[6,0], [12,0], [12,5], [6,5]]        
        
        Q1.set_values(Polygon_function([(inlet_poly, head_water_depth),
                                        (outlet_poly, tail_water_depth)]))
        
        domain.set_quantity('stage', Q0 + Q1)



        culvert = Culvert_flow(domain,
                               label='Test culvert',
                               description='4 m test culvert',   
                               end_point0=[4.0, 2.5], 
                               end_point1=[8.0, 2.5],
                               width=1.20, 
                               height=0.75,
                               culvert_routine=boyd_generalised_culvert_model,        
                               number_of_barrels=1,
                               verbose=True)
                               
        
        domain.forcing_terms.append(culvert)
        
        # Call
        culvert(domain)
    
        #print 'Inlet flow', culvert.inlet.rate
        #print 'Outlet flow', culvert.outlet.rate        
        
    
               
#-------------------------------------------------------------
if __name__ == "__main__":
    #suite = unittest.makeSuite(Test_Culvert, 'test_that_culvert_rating_limits_flow_in_shallow_inlet_condition')
    suite = unittest.makeSuite(Test_Culvert, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)
        

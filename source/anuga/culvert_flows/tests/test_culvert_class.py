#!/usr/bin/env python


import unittest
import os.path
import sys

from anuga.utilities.system_tools import get_pathname_from_package
from anuga.geometry.polygon_function import Polygon_function
        
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
from anuga.abstract_2d_finite_volumes.quantity import Quantity

import anuga

from anuga.culvert_flows.culvert_class import Culvert_flow, \
                            Culvert_flow_rating, Culvert_flow_energy
from anuga.culvert_flows.culvert_routines import boyd_generalised_culvert_model
     
from math import pi, pow, sqrt

import numpy as num


# Helper functions
def run_culvert_flow_problem(depth):
    """Run flow with culvert given depth
    """

    length = 40.
    width = 5.

    dx = dy = 1           # Resolution: Length of subdivisions on both axes

    points, vertices, boundary = rectangular_cross(int(length/dx),
                                                   int(width/dy),
                                                   len1=length, 
                                                   len2=width)
    domain = anuga.Domain(points, vertices, boundary)   
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
                # Cut Out Segment for Culvert face                
                if  1.0+(x[i]-5.0)/5.0 <  y[i]  < 4.0 - (x[i]-5.0)/5.0: 
                   z[i]=z[i]
                else:
                   z[i] +=  0.5*(x[i] -5.0)    # Sloping Segment  U/S Face
            if 10.0 < x[i] < 12.1:
               z[i] +=  2.5                    # Flat Crest of Embankment
            if 12.0 < x[i] < 14.5:
                # Cut Out Segment for Culvert face                
                if  2.0-(x[i]-12.0)/2.5 <  y[i]  < 3.0 + (x[i]-12.0)/2.5:
                   z[i]=z[i]
                else:
                   z[i] +=  2.5-1.0*(x[i] -12.0) # Sloping D/S Face
            		   
    		
        return z


    domain.set_quantity('elevation', topography) 
    domain.set_quantity('friction', 0.01)         # Constant friction 
    domain.set_quantity('stage',
                        expression='elevation + %f' % depth) # Shallow initial condition

    # Boyd culvert
    culvert = Culvert_flow(domain,
                           label='Culvert No. 1',
                           description='This culvert is a test unit 1.2m Wide by 0.75m High',   
                           end_point0=[9.0, 2.5], 
                           end_point1=[13.0, 2.5],
                           width=1.20, height=0.75,
                           culvert_routine=boyd_generalised_culvert_model,
                           number_of_barrels=1,
                           update_interval=2,
                           verbose=False)
    

    domain.forcing_terms.append(culvert)
    

    #-----------------------------------------------------------------------
    # Setup boundary conditions
    #-----------------------------------------------------------------------

    # Inflow based on Flow Depth and Approaching Momentum

    Br = anuga.Reflective_boundary(domain)              # Solid reflective wall
    domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})


    
    #-----------------------------------------------------------------------
    # Evolve system through time
    #-----------------------------------------------------------------------

    #print 'depth', depth
    ref_volume = domain.get_quantity('stage').get_integral()
    for t in domain.evolve(yieldstep = 0.1, finaltime = 10):
        new_volume = domain.get_quantity('stage').get_integral()
        
        msg = ('Total volume has changed: Is %.8f m^3 should have been %.8f m^3'
               % (new_volume, ref_volume))
        assert num.allclose(new_volume, ref_volume), msg        
        

    os.remove('Test_culvert_shallow.sww')

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
        path = os.path.join(path, 'tests', 'data')
        
        length = 40.
        width = 5.

        dx = dy = 1           # Resolution: Length of subdivisions on both axes

        points, vertices, boundary = rectangular_cross(int(length/dx),
                                                       int(width/dy),
                                                       len1=length, 
                                                       len2=width)
        domain = anuga.Domain(points, vertices, boundary)   
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
                    # Cut Out Segment for Culvert face                
                    if  1.0+(x[i]-5.0)/5.0 <  y[i]  < 4.0 - (x[i]-5.0)/5.0: 
                       z[i]=z[i]
                    else:
                       z[i] +=  0.5*(x[i] -5.0)    # Sloping Segment  U/S Face
                if 10.0 < x[i] < 12.1:
                   z[i] +=  2.5                    # Flat Crest of Embankment
                if 12.0 < x[i] < 14.5:
                    # Cut Out Segment for Culvert face                
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
                               width=1.00,
                               use_velocity_head=True,
                               verbose=False)

        domain.forcing_terms.append(culvert)
        

        #-----------------------------------------------------------------------
        # Setup boundary conditions
        #-----------------------------------------------------------------------

        # Inflow based on Flow Depth and Approaching Momentum
        Bi = anuga.Dirichlet_boundary([0.0, 0.0, 0.0])
        Br = anuga.Reflective_boundary(domain)              # Solid reflective wall
        Bo = anuga.Dirichlet_boundary([-5, 0, 0])           # Outflow
        
        # Upstream and downstream conditions that will exceed the rating curve
        # I.e produce delta_h outside the range [0, 10] specified in the the 
        # file example_rating_curve.csv
        Btus = anuga.Time_boundary(domain, \
                    lambda t: [100*num.sin(2*pi*(t-4)/10), 0.0, 0.0])
        Btds = anuga.Time_boundary(domain, \
                    lambda t: [-5*(num.cos(2*pi*(t-4)/20)), 0.0, 0.0])
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
            
            pass

        # Check that extreme values in rating curve have been exceeded
        # so that we know that condition has been exercised
        assert min_delta_w < 0
        assert max_delta_w > 10        


        os.remove('Test_culvert.sww')

    def test_that_culvert_dry_bed_rating_does_not_produce_flow(self):
        """test_that_culvert_in_dry_bed_does_not_produce_flow(self):
        
        Test that culvert on a sloping dry bed doesn't produce flows
        although there will be a 'pressure' head due to delta_w > 0

        This one is using the rating curve variant
        """

        path = get_pathname_from_package('anuga.culvert_flows')
        path = os.path.join(path, 'tests', 'data')    
        
        length = 40.
        width = 5.

        dx = dy = 1           # Resolution: Length of subdivisions on both axes

        points, vertices, boundary = rectangular_cross(int(length/dx),
                                                       int(width/dy),
                                                       len1=length, 
                                                       len2=width)
        domain = anuga.Domain(points, vertices, boundary)   
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
                    # Cut Out Segment for Culvert face                
                    if  1.0+(x[i]-5.0)/5.0 <  y[i]  < 4.0 - (x[i]-5.0)/5.0: 
                       z[i]=z[i]
                    else:
                       z[i] +=  0.5*(x[i] -5.0)    # Sloping Segment  U/S Face
                if 10.0 < x[i] < 12.1:
                   z[i] +=  2.5                    # Flat Crest of Embankment
                if 12.0 < x[i] < 14.5:
                    # Cut Out Segment for Culvert face                
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
                               height=0.75,
                               verbose=False)

        domain.forcing_terms.append(culvert)
        

        #-----------------------------------------------------------------------
        # Setup boundary conditions
        #-----------------------------------------------------------------------

        # Inflow based on Flow Depth and Approaching Momentum

        Br = anuga.Reflective_boundary(domain)              # Solid reflective wall
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})


        #-----------------------------------------------------------------------
        # Evolve system through time
        #-----------------------------------------------------------------------

        ref_volume = domain.get_quantity('stage').get_integral()
        for t in domain.evolve(yieldstep = 1, finaltime = 25):
            new_volume = domain.get_quantity('stage').get_integral()
            
            msg = 'Total volume has changed'
            assert num.allclose(new_volume, ref_volume, rtol=1.0e-10), msg
            pass
    

    
        os.remove('Test_culvert_dry.sww')
        
    def test_that_culvert_flows_conserves_volume(self):
        """test_that_culvert_flows_conserves_volume

        Test that culvert on a sloping dry bed limits flows when very little water
        is present at inlet.

        Uses helper function: run_culvert_flow_problem(depth):

        """

        # Try this for a range of depths
        for depth in [0.1, 1.0]: #[0.1, 0.2, 0.5, 1.0]:
            run_culvert_flow_problem(depth)


    def OBSOLETE_XXXtest_that_culvert_rating_limits_flow_in_shallow_inlet_condition(self):
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
        domain = anuga.Domain(points, vertices, boundary)   
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
                    # Cut Out Segment for Culvert face                
                    if  1.0+(x[i]-5.0)/5.0 <  y[i]  < 4.0 - (x[i]-5.0)/5.0: 
                       z[i]=z[i]
                    else:
                       z[i] +=  0.5*(x[i] -5.0)    # Sloping Segment  U/S Face
                if 10.0 < x[i] < 12.1:
                   z[i] +=  2.5                    # Flat Crest of Embankment
                if 12.0 < x[i] < 14.5:
                    # Cut Out Segment for Culvert face                
                    if  2.0-(x[i]-12.0)/2.5 <  y[i]  < 3.0 + (x[i]-12.0)/2.5:
                       z[i]=z[i]
                    else:
                       z[i] +=  2.5-1.0*(x[i] -12.0) # Sloping D/S Face
                		   
        		
            return z


        domain.set_quantity('elevation', topography) 
        domain.set_quantity('friction', 0.01)         # Constant friction 
        domain.set_quantity('stage',
                            expression='elevation + 0.1') # Shallow initial condition

        # Boyd culvert
        culvert = Culvert_flow(domain,
                               label='Culvert No. 1',
                               description='This culvert is a test unit 1.2m Wide by 0.75m High',   
                               end_point0=[9.0, 2.5], 
                               end_point1=[13.0, 2.5],
                               width=1.20, height=0.75,
                               culvert_routine=boyd_generalised_culvert_model,
                               number_of_barrels=1,
                               update_interval=2,
                               verbose=False)
        
        # Rating curve
        #filename = os.path.join(path, 'example_rating_curve.csv')
        #culvert = Culvert_flow(domain,
        #                       culvert_description_filename=filename,
        #                       end_point0=[9.0, 2.5], 
        #                       end_point1=[13.0, 2.5],
        #                       trigger_depth=0.01,
        #                       verbose=False)

        domain.forcing_terms.append(culvert)
        

        #-----------------------------------------------------------------------
        # Setup boundary conditions
        #-----------------------------------------------------------------------

        # Inflow based on Flow Depth and Approaching Momentum

        Br = anuga.Reflective_boundary(domain)              # Solid reflective wall
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})


        
        #-----------------------------------------------------------------------
        # Evolve system through time
        #-----------------------------------------------------------------------

        print 'depth', 0.1
        ref_volume = domain.get_quantity('stage').get_integral()
        for t in domain.evolve(yieldstep = 0.1, finaltime = 25):
            new_volume = domain.get_quantity('stage').get_integral()

            msg = ('Total volume has changed: Is %.8f m^3 should have been %.8f m^3'
                   % (new_volume, ref_volume))
            assert num.allclose(new_volume, ref_volume, rtol=1.0e-10), msg        
        
        
        return
        # Now try this again for a depth of 10 cm and for a range of other depths
        for depth in [0.1, 0.2, 0.5, 1.0]:
            print 'depth', depth
            domain.set_time(0.0)
            
            domain.set_quantity('elevation', topography) 
            domain.set_quantity('friction', 0.01)         # Constant friction 
            domain.set_quantity('stage',
                                expression='elevation + %f' % depth)
            
        
            ref_volume = domain.get_quantity('stage').get_integral()
            for t in domain.evolve(yieldstep = 0.1, finaltime = 25):
                new_volume = domain.get_quantity('stage').get_integral()
            
                msg = 'Total volume has changed: Is %.8f m^3 should have been %.8f m^3'\
                    % (new_volume, ref_volume)

                assert num.allclose(new_volume, ref_volume, rtol=1.0e-10), msg
    
    
    
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
        domain = anuga.Domain(points, vertices, boundary)   
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
                    # Cut Out Segment for Culvert face                
                    if  1.0+(x[i]-5.0)/5.0 <  y[i]  < 4.0 - (x[i]-5.0)/5.0: 
                       z[i]=z[i]
                    else:
                       z[i] +=  0.5*(x[i] -5.0)    # Sloping Segment  U/S Face
                if 10.0 < x[i] < 12.1:
                   z[i] +=  2.5                    # Flat Crest of Embankment
                if 12.0 < x[i] < 14.5:
                    # Cut Out Segment for Culvert face                
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
                               width=1.20, height=0.75,
                               culvert_routine=boyd_generalised_culvert_model,
                               number_of_barrels=1,
                               update_interval=2,
                               verbose=False)
        
        domain.forcing_terms.append(culvert)
        

        #-----------------------------------------------------------------------
        # Setup boundary conditions
        #-----------------------------------------------------------------------

        # Inflow based on Flow Depth and Approaching Momentum

        Br = anuga.Reflective_boundary(domain)              # Solid reflective wall
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})


        #-----------------------------------------------------------------------
        # Evolve system through time
        #-----------------------------------------------------------------------

        ref_volume = domain.get_quantity('stage').get_integral()
        for t in domain.evolve(yieldstep = 1, finaltime = 25):
            
            new_volume = domain.get_quantity('stage').get_integral()

            msg = 'Total volume has changed'
            assert num.allclose(new_volume, ref_volume, rtol=1.0e-10), msg
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
        domain = anuga.Domain(points, vertices, boundary)   
        
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
                               verbose=False)
                               
        
        domain.forcing_terms.append(culvert)
        
        # Call
        culvert(domain)
    

    

    def test_momentum_jet(self):
        """test_momentum_jet
        
        Test that culvert_class can accept keyword use_momentum_jet
        This does not yet imply that the values have been tested. FIXME
        """


        length = 40.
        width = 5.

        dx = dy = 1           # Resolution: Length of subdivisions on both axes

        points, vertices, boundary = rectangular_cross(int(length/dx),
                                                       int(width/dy),
                                                       len1=length, 
                                                       len2=width)
        domain = anuga.Domain(points, vertices, boundary)   
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
                    # Cut Out Segment for Culvert face                
                    if  1.0+(x[i]-5.0)/5.0 <  y[i]  < 4.0 - (x[i]-5.0)/5.0: 
                       z[i]=z[i]
                    else:
                       z[i] +=  0.5*(x[i] -5.0)    # Sloping Segment  U/S Face
                if 10.0 < x[i] < 12.1:
                   z[i] +=  2.5                    # Flat Crest of Embankment
                if 12.0 < x[i] < 14.5:
                    # Cut Out Segment for Culvert face                
                    if  2.0-(x[i]-12.0)/2.5 <  y[i]  < 3.0 + (x[i]-12.0)/2.5:
                       z[i]=z[i]
                    else:
                       z[i] +=  2.5-1.0*(x[i] -12.0) # Sloping D/S Face
                		   
        		
            return z


        domain.set_quantity('elevation', topography) 
        domain.set_quantity('friction', 0.01)         # Constant friction 
        domain.set_quantity('stage',
                            expression='elevation + 1.0') # Shallow initial condition

        # Boyd culvert
        culvert = Culvert_flow(domain,
                               label='Culvert No. 1',
                               description='This culvert is a test unit 1.2m Wide by 0.75m High',   
                               end_point0=[9.0, 2.5], 
                               end_point1=[13.0, 2.5],
                               width=1.20, height=0.75,
                               culvert_routine=boyd_generalised_culvert_model,
                               number_of_barrels=1,
                               use_momentum_jet=True,
                               update_interval=2,
                               verbose=False)
    

        domain.forcing_terms.append(culvert)

        
        # Call
        culvert(domain)
    

        #-----------------------------------------------------------------------
        # Setup boundary conditions
        #-----------------------------------------------------------------------


        Br = anuga.Reflective_boundary(domain)              # Solid reflective wall
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        #-----------------------------------------------------------------------
        # Evolve system through time
        #-----------------------------------------------------------------------

        ref_volume = domain.get_quantity('stage').get_integral()
        for t in domain.evolve(yieldstep = 0.1, finaltime = 25):
            new_volume = domain.get_quantity('stage').get_integral()
            
            msg = ('Total volume has changed: Is %.8f m^3 should have been %.8f m^3'
                   % (new_volume, ref_volume))
            assert num.allclose(new_volume, ref_volume, rtol=1.0e-10), msg        
        

        

               
#-------------------------------------------------------------

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_Culvert, 'test')
    runner = unittest.TextTestRunner() #verbosity=2)
    runner.run(suite)
        

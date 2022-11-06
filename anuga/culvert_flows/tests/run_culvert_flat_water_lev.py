""" Testing CULVERT (Changing from Horizontal Abstraction to Vertical Abstraction

This example includes a Model Topography that shows a TYPICAL Headwall Configuration

The aim is to change the Culvert Routine to Model more precisely the abstraction
from a vertical face.

The inflow must include the impact of Approach velocity.
Similarly the Outflow has MOMENTUM Not just Up welling as in the Horizontal Style
abstraction

"""

print('Starting.... Importing Modules...')

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross

from anuga.shallow_water import Domain, Reflective_boundary,\
     Dirichlet_boundary,\
     Transmissive_boundary, Time_boundary

from anuga.culvert_flows.culvert_class import Culvert_flow
from anuga.culvert_flows.culvert_routines import boyd_generalised_culvert_model
     
from math import pi,pow,sqrt

import numpy as num


#------------------------------------------------------------------------------
# Setup computational domain
#------------------------------------------------------------------------------
print('Setting up domain')

length = 40.
width = 5.

dx = dy = 1           # Resolution: Length of subdivisions on both axes
#dx = dy = .5           # Resolution: Length of subdivisions on both axes
#dx = dy = .5           # Resolution: Length of subdivisions on both axes
#dx = dy = .1           # Resolution: Length of subdivisions on both axes

points, vertices, boundary = rectangular_cross(int(length/dx), int(width/dy),
                                               len1=length, len2=width)
domain = Domain(points, vertices, boundary)   
domain.set_name('Test_Culv_Flat_WL')                 # Output name
domain.set_default_order(2)
domain.H0 = 0.01
domain.tight_slope_limiters = 1

print('Size', len(domain))

#------------------------------------------------------------------------------
# Setup initial conditions
#------------------------------------------------------------------------------

def topography(x, y):
    """Set up a weir
    
    A culvert will connect either side
    """
    # General Slope of Topography
    z=-x/1000
    
    #       NOW Add bits and Pieces to topography
    N = len(x)
    for i in range(N):

       # Sloping Embankment Across Channel
        if 5.0 < x[i] < 10.1:
            if  1.0+(x[i]-5.0)/5.0 <  y[i]  < 4.0 - (x[i]-5.0)/5.0: # Cut Out Segment for Culvert FACE
               z[i]=z[i]
            else:
               z[i] +=  0.5*(x[i] -5.0)    # Sloping Segment  U/S Face
        if 10.0 < x[i] < 12.1:
           z[i] +=  2.5              # Flat Crest of Embankment
        if 12.0 < x[i] < 14.5:
            if  2.0-(x[i]-12.0)/2.5 <  y[i]  < 3.0 + (x[i]-12.0)/2.5: # Cut Out Segment for Culvert FACE
               z[i]=z[i]
            else:
               z[i] +=  2.5-1.0*(x[i] -12.0)       # Sloping D/S Face
        		   
        
		
    return z

print('Setting Quantities....')
domain.set_quantity('elevation', topography)  # Use function for elevation
domain.set_quantity('friction', 0.01)         # Constant friction 
domain.set_quantity('stage',
                    expression='elevation')   # Dry initial condition




#------------------------------------------------------------------------------
# Setup specialised forcing terms
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Setup CULVERT INLETS and OUTLETS in Current Topography
#------------------------------------------------------------------------------
print('DEFINING any Structures if Required')

#  DEFINE CULVERT INLET AND OUTLETS


culvert_rating = Culvert_flow(domain,
                       culvert_description_filename='example_rating_curve.csv',
                       end_point0=[9.0, 2.5], 
                       end_point1=[13.0, 2.5],
                       verbose=True)


culvert_energy = Culvert_flow(domain,
                       label='Culvert No. 1',
                       description='This culvert is a test unit 1.2m Wide by 0.75m High',   
                       end_point0=[9.0, 2.5], 
                       end_point1=[13.0, 2.5],
                       width=1.20,height=0.75,
                       culvert_routine=boyd_generalised_culvert_model,        
                       number_of_barrels=1,
                       update_interval=2,
                       log_file=True,
                       discharge_hydrograph=True,
                       verbose=True)

domain.forcing_terms.append(culvert_energy)

#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
print('Setting Boundary Conditions')
Bi = Dirichlet_boundary([0.0, 0.0, 0.0])          # Inflow based on Flow Depth and Approaching Momentum !!!
Br = Reflective_boundary(domain)              # Solid reflective wall
Bo = Dirichlet_boundary([-5, 0, 0])           # Outflow
Btus = Time_boundary(domain, lambda t: [0.0+ 1.25*(1+num.sin(2*pi*(t-4)/10)), 0.0, 0.0])
Btds = Time_boundary(domain, lambda t: [0.0+ 0.75*(1+num.sin(2*pi*(t-4)/20)), 0.0, 0.0])
domain.set_boundary({'left': Btus, 'right': Btds, 'top': Br, 'bottom': Br})


#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------

#for t in domain.evolve(yieldstep = 1, finaltime = 25):
#    print domain.timestepping_statistics()
    

    
    
#import sys; sys.exit() 
# Profiling code
import time
t0 = time.time()
    
s = 'for t in domain.evolve(yieldstep = 1, finaltime = 25): domain.write_time()'

import profile, pstats
FN = 'profile.dat'

profile.run(s, FN)
    
print('That took %.2f seconds' %(time.time()-t0))

S = pstats.Stats(FN)
#S.sort_stats('time').print_stats(20)
s = S.sort_stats('cumulative').print_stats(30)

print(s)

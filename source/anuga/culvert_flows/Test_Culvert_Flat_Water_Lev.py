""" Testing CULVERT (Changing from Horizontal Abstraction to Vertical Abstraction

This example includes a Model Topography that shows a TYPICAL Headwall Configuration

The aim is to change the Culvert Routine to Model more precisely the abstraction
from a vertical face.

The inflow must include the impact of Approach velocity.
Similarly the Outflow has MOMENTUM Not just Up welling as in the Horizontal Style
abstraction

"""
print 'Starting.... Importing Modules...'
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
from Numeric import choose, greater, ones, sin, exp, cosh
#------------------------------------------------------------------------------
# Setup computational domain
#------------------------------------------------------------------------------
print 'Setting up domain'

# Open file for storing some specific results...
fid = open('Culvert_Headwall', 'w')

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

print 'Size', len(domain)

#------------------------------------------------------------------------------
# Setup initial conditions
#------------------------------------------------------------------------------

def topography(x, y):
    """Set up a weir
    
    A culvert will connect either side
    """
    # General Slope of Topography
    z=-x/1000
    
	# Changing Slopes in the X- Direction
    # N = len(x)
    # for i in range(N):
    #     if 0 <x[i] < 5.1:
    #         z[i] -= -x[i]/10
    #     if 5 <x[i] < 10.1:
    #         z[i] -= -x[i]/100                              
    #     if 10 <x[i]: 
    #         z[i] -= -x[i]/20                            
    # return z

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
        		   
        
        # Constriction
        #if 27 < x[i] < 29 and y[i] > 3:
        #    z[i] += 2        
        
        # Pole
        #if (x[i] - 34)**2 + (y[i] - 2)**2 < 0.4**2:
        #    z[i] += 2

        # HOLE For Pit at Opening[0]
        #if (x[i]-4)**2 + (y[i]-1.5)**2<0.75**2:
         #  z[i]-=1		

	   # HOLE For Pit at Opening[1]
        #if (x[i]-20)**2 + (y[i]-3.5)**2<0.5**2:
         #  z[i]-=1		
		
    return z

print 'Setting Quantities....'
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
print 'DEFINING any Structures if Required'

#  DEFINE CULVERT INLET AND OUTLETS


culvert = Culvert_flow(domain,
                       label='Culvert No. 1',
                       description='This culvert is a test unit 1.2m Wide by 0.75m High',   
                       end_point0=[9.0, 2.5], 
                       end_point1=[13.0, 2.5],
                       width=1.20,height=0.75,
                       culvert_routine=boyd_generalised_culvert_model,        
                       number_of_barrels=1,
                       verbose=True)

domain.forcing_terms.append(culvert)

#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
print 'Setting Boundary Conditions'
#Bi = Dirichlet_boundary([0.5, 0.0, 0.0])          # Inflow based on Flow Depth (0.5m) and Approaching Momentum !!!
Bi = Dirichlet_boundary([0.0, 0.0, 0.0])          # Inflow based on Flow Depth and Approaching Momentum !!!
Br = Reflective_boundary(domain)              # Solid reflective wall
Bo = Dirichlet_boundary([-5, 0, 0])           # Outflow
Btus = Time_boundary(domain, lambda t: [0.0+ 1.25*(1+sin(2*pi*(t-4)/10)), 0.0, 0.0])
Btds = Time_boundary(domain, lambda t: [0.0+ 0.75*(1+sin(2*pi*(t-4)/20)), 0.0, 0.0])
domain.set_boundary({'left': Btus, 'right': Btds, 'top': Br, 'bottom': Br})

#------------------------------------------------------------------------------
# Setup Application of  specialised forcing terms
#------------------------------------------------------------------------------
print 'Setting up Forcing Terms if Required'
# This is the new element implemented by Ole to allow direct input of Inflow in m^3/s
#fixed_flow = Inflow(domain, center=(2.1, 2.1), radius=1.261566, flow=1.00)   #   Fixed Flow Value Over Area of 5m2 at 1m/s = 5m^3/s

#	     flow=file_function('Q/QPMF_Rot_Sub13.tms'))        # Read Time Series in  from File
#             flow=lambda t: min(0.01*t, 0.01942))                             # Time Varying Function   Tap turning up  	

#domain.forcing_terms.append(fixed_flow)





#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------
print 'STARTING EVOLVE ======>'


for t in domain.evolve(yieldstep = 0.01, finaltime = 45):
     if int(domain.time*100) % 100 == 0:
	     domain.write_time()
    
    #if domain.get_time() >= 4 and tap.flow != 0.0:
    #    print 'Turning tap off'
    #    tap.flow = 0.0
	
    #if domain.get_time() >= 3 and sink.flow < 0.0:
    #    print 'Turning drain on'
    #    sink.flow = -0.8	
    # Close

fid.close()


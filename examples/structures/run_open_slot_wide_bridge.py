""" Testing CULVERT (Changing from Horizontal Abstraction to Vertical Abstraction

This example includes a Model Topography that shows a TYPICAL Headwall Configuration

The aim is to change the Culvert Routine to Model more precisely the abstraction
from a vertical face.

The inflow must include the impact of Approach velocity.
Similarly the Outflow has MOMENTUM Not just Up welling as in the Horizontal Style
abstraction

"""
from math import pi, pow, sqrt
import numpy as num
from anuga.culvert_flows.culvert_routines import boyd_generalised_culvert_model
from anuga.culvert_flows.culvert_class import Culvert_flow
from anuga.shallow_water.forcing import Rainfall, Inflow
from anuga import Domain
from anuga import rectangular_cross_domain
import anuga
print('Starting.... Importing Modules...')
# ------------------------------------------------------------------------------
# Import necessary modules
# ------------------------------------------------------------------------------

#from anuga.shallow_water.forcing import Reflective_boundary
#from anuga.shallow_water.forcing import Dirichlet_boundary
#from anuga.shallow_water.forcing import Transmissive_boundary, Time_boundary

#from anuga.culvert_flows.culvert_routines import weir_orifice_channel_culvert_model


# ------------------------------------------------------------------------------
# Setup computational domain
# ------------------------------------------------------------------------------
print('Setting up domain')

length = 200.  # x-Dir
width = 200.  # y-dir

dx = dy = 2.5          # Resolution: Length of subdivisions on both axes
# dx = dy = .5           # Resolution: Length of subdivisions on both axes
# dx = dy = .5           # Resolution: Length of subdivisions on both axes
# dx = dy = .1           # Resolution: Length of subdivisions on both axes

domain = rectangular_cross_domain(int(length/dx), int(width/dy),
                                  len1=length, len2=width)

domain.set_name('Test_open_slot_wide_bridge')                 # Output name
# domain.set_default_order(2)
# domain.H0 = 0.01
# domain.tight_slope_limiters = 1

print('Size', len(domain))

# ------------------------------------------------------------------------------
# Setup initial conditions
# ------------------------------------------------------------------------------


def topography(x, y):
    """Set up a weir

    A culvert will connect either side
    """
    # General Slope of Topography
    z = 10.0-x/100.0  # % Longitudinal Slope

    #       NOW Add bits and Pieces to topography
    bank_hgt = 10.0
    bridge_width = 50.0
    bank_width = 10.0

    us_apron_skew = 1.0  # 1.0 = 1 Length: 1 Width, 2.0 = 2 Length : 1 Width
    us_start_x = 10.0
    top_start_y = 50.0
    us_slope = 3.0  # Horiz : Vertic
    ds_slope = 3.0
    ds_apron_skew = 1.0  # 1.0 = 1 Length: 1 Width, 2.0 = 2 Length : 1 Width
    centre_line_y = top_start_y+bridge_width/2.0

    # CALCULATE PARAMETERS TO FORM THE EMBANKMENT
    us_slope_length = bank_hgt*us_slope
    us_end_x = us_start_x + us_slope_length
    us_toe_start_y = top_start_y - us_slope_length / us_apron_skew
    us_toe_end_y = top_start_y + bridge_width + us_slope_length / us_apron_skew

    top_end_y = top_start_y + bridge_width
    ds_slope_length = bank_hgt*ds_slope
    ds_start_x = us_end_x + bank_width
    ds_end_x = ds_start_x + ds_slope_length

    ds_toe_start_y = top_start_y - ds_slope_length / ds_apron_skew
    ds_toe_end_y = top_start_y + bridge_width + ds_slope_length / ds_apron_skew

    N = len(x)
    for i in range(N):

       # Sloping Embankment Across Channel
        if us_start_x < x[i] < us_end_x + 0.1:   # For UPSLOPE on the Upstream FACE
            # if 5.0 < x[i] < 10.1: # For a Range of X, and over a Range of Y based on X adjust Z
            if us_toe_start_y + (x[i] - us_start_x)/us_apron_skew < y[i] < us_toe_end_y - (x[i] - us_start_x)/us_apron_skew:
                # if  49.0+(x[i]-5.0)/5.0 <  y[i]  < 151.0 - (x[i]-5.0)/5.0: # Cut Out Base Segment for Culvert FACE
                z[i] = z[i]  # Flat Apron
            else:
                # Set elevation for Sloping Segment  U/S Face
                z[i] += z[i] + (x[i] - us_start_x)/us_slope
        if us_end_x < x[i] < ds_start_x + 0.1:  # FOR The TOP of BANK Segment
            if top_start_y < y[i] < top_end_y:
                z[i] = z[i]  # Flat Apron
            else:
                z[i] += z[i]+bank_hgt        # Flat Crest of Embankment
        if ds_start_x < x[i] < ds_end_x:  # DOWN SDLOPE Segment on Downstream face
            # Cut Out Segment for Culvert FACE
            if top_start_y-(x[i]-ds_start_x)/ds_apron_skew < y[i] < top_end_y + (x[i]-ds_start_x)/ds_apron_skew:
                z[i] = z[i]  # Flat Apron
            else:
                z[i] += z[i]+bank_hgt-(x[i] - ds_start_x) / \
                    ds_slope       # Sloping D/S Face

    return z


print('Setting Quantities....')
domain.set_quantity('elevation', topography)  # Use function for elevation
domain.set_quantity('friction', 0.01)         # Constant friction
domain.set_quantity('stage',
                    expression='elevation')   # Dry initial condition


# ------------------------------------------------------------------------------
# Setup specialised forcing terms
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Setup CULVERT INLETS and OUTLETS in Current Topography
# ------------------------------------------------------------------------------
print('DEFINING any Structures if Required')

#  DEFINE CULVERT INLET AND OUTLETS

"""
culvert_rating = Culvert_flow(domain,
                       culvert_description_filename='example_rating_curve.csv',
                       end_point0=[40.0, 75.0], 
                       end_point1=[50.0, 75.0],
                       verbose=True)


culvert_energy = Culvert_flow(domain,
                       label='Culvert No. 1',
                       description='This culvert is a test unit 1.2m Wide by 0.75m High',   
                       end_point0=[40.0, 75.0], 
                       end_point1=[50.0, 75.0],
                       width=50.0,depth=5.0,
                       culvert_routine=boyd_generalised_culvert_model,        
                       number_of_barrels=1,
                       #update_interval=0.25,
                       log_file=True,
                       discharge_hydrograph=True,
                       use_velocity_head=False,
                       verbose=True)

domain.forcing_terms.append(culvert_energy)
"""

# ------------------------------------------------------------------------------
# Setup boundary conditions
# ------------------------------------------------------------------------------
print('Setting Boundary Conditions')
Br = anuga.Reflective_boundary(domain)              # Solid reflective wall
# Inflow based on Flow Depth and Approaching Momentum !!!
Bi = anuga.Dirichlet_boundary([0.0, 0.0, 0.0])

Bo = anuga.Dirichlet_boundary([-5.0, 0, 0])           # Outflow water at -5.0
Bd = anuga.Dirichlet_boundary([0, 0, 0])                # Outflow water at 0.0

#Btus = Time_boundary(domain, lambda t: [0.0+ 1.025*(1+num.sin(2*pi*(t-4)/10)), 0.0, 0.0])
#Btds = Time_boundary(domain, lambda t: [0.0+ 0.0075*(1+num.sin(2*pi*(t-4)/20)), 0.0, 0.0])

Btus = anuga.Dirichlet_boundary([18.0, 0, 0])           # Outflow water at 5.0
Btds = anuga.Dirichlet_boundary([0.0, 0, 0])           # Outflow water at 1.0
domain.set_boundary({'left': Btus, 'right': Btds, 'top': Br, 'bottom': Br})


# ------------------------------------------------------------------------------
# Evolve system through time
# ------------------------------------------------------------------------------

for t in domain.evolve(yieldstep=1, finaltime=100):
    domain.print_timestepping_statistics()

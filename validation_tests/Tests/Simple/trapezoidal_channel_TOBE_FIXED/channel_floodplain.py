"""Simple water flow example using ANUGA
Water flowing down a channel with a floodplain
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
# Import standard shallow water domain and standard boundaries.
import anuga
import numpy
from anuga.structures.inlet_operator import Inlet_operator
from anuga import *

#------------------------------------------------------------------------------
# Useful parameters for controlling this case
#------------------------------------------------------------------------------

floodplain_length = 1000.0 # Model domain length
floodplain_width = 14.0 # Model domain width
floodplain_slope = 1./300.
chan_initial_depth = 0.65 # Initial depth of water in the channel
chan_bankfull_depth = 1.0 # Bankfull depth of the channel
chan_width = 10.0 # Bankfull width of the channel
bankwidth = 2. # Width of the bank regions -- note that these protrude into the channel
man_n=0.03 # Manning's n
l0 = 1.000 # Length scale associated with triangle side length in channel (min_triangle area = 0.5*l0^2)

assert chan_width < floodplain_width, \
        ' ERROR: Channel width is greater than floodplain width'

assert bankwidth < chan_width/2., \
        'ERROR: The bank width must be less than half the channel width'

#------------------------------------------------------------------------------
# Setup computational domain
#------------------------------------------------------------------------------

# Define boundary polygon -- in this case clockwise around the proposed boundary
boundary_polygon = [ [0.,0.], 
                     [0., floodplain_length], 
                     [floodplain_width/2. - chan_width/2., floodplain_length], 
                     [floodplain_width/2. + chan_width/2., floodplain_length], 
                     [floodplain_width, floodplain_length], 
                     [floodplain_width, 0.], 
                     [floodplain_width/2. + chan_width/2., 0.], 
                     [floodplain_width/2. - chan_width/2., 0.] 
                     ]

# Define channel polygon, where we can optionally refine the resolution. 
# Note that we set this a distance l0 inside the boundary polygon, so the polygons
# do not overlap. 
channel_polygon = [ [floodplain_width/2. - chan_width/2., +l0],
                    [floodplain_width/2. - chan_width/2., floodplain_length-l0],
                    [floodplain_width/2. + chan_width/2., floodplain_length-l0],
                    [floodplain_width/2. + chan_width/2., +l0]
                    ]


# Define domain with appropriate boundary conditions
domain = create_domain_from_regions( boundary_polygon, 
                                   boundary_tags={'left': [0],
                                                  'top1': [1],
                                                  'chan_out': [2],
                                                  'top2': [3],
                                                  'right': [4],
                                                  'bottom1': [5],
                                                  'chan_in': [6],
                                                  'bottom2': [7] },
                                   maximum_triangle_area = 0.5*l0*l0,
                                   minimum_triangle_angle = 28.0,
                                   mesh_filename = 'channel_floodplain.msh',
                                   interior_regions = [ ],
                                   #interior_regions = [\
                                   #    [channel_polygon, 0.5*l0*l0] ],
                                   use_cache=False,
                                   verbose=True)


domain.set_name('channel_floodplain') # Output name

#------------------------------------------------------------------------------
# Setup Algorithm, either using command line arguments
# or override manually yourself
#------------------------------------------------------------------------------
from anuga.utilities.argparsing import parse_standard_args
alg, cfl = parse_standard_args()
domain.set_flow_algorithm(alg)
domain.set_CFL(cfl)

#------------------------------------------------------------------------------
#
# Setup initial conditions
#
#------------------------------------------------------------------------------

# Function for topography
def topography(x, y):
    # Longitudinally sloping floodplain with channel in centre
	#return -y*floodplain_slope -chan_bankfull_depth*\
    #        (x>(floodplain_width/2. - chan_width/2.) )*\
    #        (x<(floodplain_width/2. + chan_width/2.) ) 

    elev1= -y*floodplain_slope - chan_bankfull_depth*\
            (x>(floodplain_width/2. - chan_width/2.))*\
            (x<(floodplain_width/2. + chan_width/2.)) 
    # Add banks
    if(bankwidth>0.0):
        leftbnk = floodplain_width/2. - chan_width/2.
        rightbnk = floodplain_width/2. + chan_width/2.
        # Left bank
        elev2 = elev1 + (chan_bankfull_depth \
                - chan_bankfull_depth/bankwidth*(x - leftbnk))*\
                (x>leftbnk)*(x < leftbnk + bankwidth)
        # Right bank
        elev2 = elev2 + (chan_bankfull_depth \
                + chan_bankfull_depth/bankwidth*(x - rightbnk))*\
                (x>rightbnk-bankwidth)*(x < rightbnk)
    
    if(bankwidth==0.0):
        elev2 = elev1
    #
    return elev2

#

#Function for stage
def stagetopo(x,y):
    return -y*floodplain_slope -chan_bankfull_depth + chan_initial_depth 

domain.set_quantity('elevation', topography) # Use function for elevation
domain.set_quantity('friction', man_n) # Constant friction
domain.set_quantity('stage', stagetopo) # Use function for stage
#domain.set_minimum_allowed_height(0.01) # 

# Define inlet operator 
flow_in_yval=0.0
if True:
    line1 = [ [floodplain_width/2. - chan_width/2., flow_in_yval],\
              [floodplain_width/2. + chan_width/2., flow_in_yval] \
              ]
    Qin = 0.5*(floodplain_slope*(chan_width*chan_initial_depth)**2.*man_n**(-2.)\
            *chan_initial_depth**(4./3.) )**0.5
    
    Inlet_operator(domain, line1, Qin)
    
    print 'Discharge in = ', Qin 

#------------------------------------------------------------------------------
#
# Setup boundary conditions
#
#------------------------------------------------------------------------------

Br = anuga.Reflective_boundary(domain) # Solid reflective wall
#Bt = anuga.Transmissive_boundary(domain) # Transmissive boundary

def outflow_stage_boundary(t):
    return -floodplain_length*floodplain_slope \
            + chan_initial_depth - chan_bankfull_depth

# Note that the outflow boundary may be slightly incorrect for the trapezoidal channel case, 
# or incorrect more generally if there are numerical problems. But, in the central regions of
# the channel, this shouldn't prevent us reaching steady, uniform flow.
Bout_tmss = anuga.shallow_water.boundaries.Transmissive_momentum_set_stage_boundary(domain, function = outflow_stage_boundary) 

domain.set_boundary({'left': Br, 
                     'right': Br, 
                     'top1': Br, 
                     'top2': Br, 
                     'bottom1': Br, 
                     'bottom2': Br, 
                     'chan_out': Bout_tmss, 
                     'chan_in': Br})

# Set up file to record computations of discharge at several points.
discharge_outfile=open('discharge_outputs.txt', 'w')
discharge_outfile.write('Time (s)'+","+ 'Discharge@10' + ","+ 'Discharge@700'+","+ 'Discharge@1000' + "\n")

#------------------------------------------------------------------------------
#
# Evolve system through time
#
#------------------------------------------------------------------------------
for t in domain.evolve(yieldstep=2.0, finaltime=1200.0):
    print domain.timestepping_statistics()
    xx=domain.quantities['ymomentum'].centroid_values
    dd=(domain.quantities['stage'].centroid_values - domain.quantities['elevation'].centroid_values)
    dd=dd*(dd>0.)

    tmp = xx/(dd+1.0e-06)*(dd>0.0)
    print tmp.max(), tmp.argmax(), tmp.min(),  tmp.argmin()

    # Compute flow through cross-section -- check that the inflow boundary condition is doing its job
    # This also provides another useful steady-state check
    if( numpy.floor(t/100.) == t/100. ):
        print '#### COMPUTING FLOW THROUGH CROSS-SECTIONS########'
        s0 = domain.get_flow_through_cross_section([[0., 10.0], [floodplain_width, 10.0]])
        s1 = domain.get_flow_through_cross_section([[0., floodplain_length-300.0], [floodplain_width, floodplain_length-300.0]])
        s2 = domain.get_flow_through_cross_section([[0., floodplain_length-1.0], [floodplain_width, floodplain_length-1.0]])
        
        print 'Cross sectional flows: ',s0, s1, s2 
        print '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
        discharge_outfile.write(str(t) + "," +str(s0) + ","+ str(s1) +"," + str(s2) + "\n")        


discharge_outfile.close()

"""
Simple tidal example with ANUGA
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
# Import standard shallow water domain and standard boundaries.
import anuga
import numpy
from anuga import Inlet_operator, Boyd_box_operator
from anuga import distribute, myid, numprocs, finalize, barrier



args = anuga.get_args()
alg = args.alg
verbose = args.verbose

#------------------------------------------------------------------------------
# Useful parameters for controlling this case
#------------------------------------------------------------------------------

floodplain_length = 1000.0 # Model domain length
floodplain_width = 40.0 # Model domain width
floodplain_slope = 0. 
chan_initial_depth = 0.25 # Initial depth of water in the channel
chan_bankfull_depth = 1.0 # Bankfull depth of the channel
chan_width = 20.0 # Bankfull width of the channel
bankwidth = 0.01 # Width of the bank regions -- note that these protrude into the channel
man_n=0.03 # Manning's n
l0 = 50.000 # Length scale associated with triangle side length in channel (min_triangle area = 0.5*l0^2)

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

meshBreakLines = { # Put breaklines around edge of channel [think it helps the numerics] 
               'n1': [[floodplain_width/2.0 - chan_width/2., floodplain_length],
                      [floodplain_width/2.0 - chan_width/2., 0.]],
               'n3': [[floodplain_width/2.0 + chan_width/2.0, floodplain_length],
                      [floodplain_width/2.0 + chan_width/2.0,     0.]],
             }

regionPtAreas=[ [0.01, 0.01, 0.5*l0*l0*4],
                [floodplain_width/2., 0.01, 0.5*l0*l0],
                [floodplain_width-0.01, 0.01, 0.5*l0*l0*4], 
              ]

if(myid==0):
    # Define domain with appropriate boundary conditions
    anuga.create_mesh_from_regions( boundary_polygon, 
                                       boundary_tags={'left': [0],
                                                      'top1': [1],
                                                      'chan_out': [2],
                                                      'top2': [3],
                                                      'right': [4],
                                                      'bottom1': [5],
                                                      'chan_in': [6],
                                                      'bottom2': [7] },
                                       maximum_triangle_area = 1.0e+06, #0.5*l0*l0,
                                       minimum_triangle_angle = 28.0,
                                       filename = 'channel_floodplain1.msh',
                                       interior_regions = [ ],
                                       breaklines=meshBreakLines.values(),
                                       regionPtArea=regionPtAreas,
                                       verbose=True)
    #
    domain=anuga.create_domain_from_file('channel_floodplain1.msh')
    outname='channel_floodplain1'#+str(smoothing_timescale)
    domain.set_name(outname) # Output name
    domain.set_flow_algorithm(alg)
else:
    domain=None
#
barrier()
domain=distribute(domain)
barrier()

domain.set_store_vertices_uniquely(True)

#------------------------------------------------------------------------------
#
# Setup initial conditions
#
#------------------------------------------------------------------------------

# Function for topography
def topography(x, y):
    # Longitudinally sloping floodplain with channel in centre
    channel_depth=1. + (floodplain_length-y>100.)*(floodplain_length-y<220.)*((220.-(floodplain_length-y))/120.)*2. + 2.*(floodplain_length-y<=100.)

    elev1= -y*floodplain_slope - channel_depth*\
            (x>(floodplain_width/2. - chan_width/2.))*\
            (x<(floodplain_width/2. + chan_width/2.)) 
    # Add banks
    if(bankwidth>0.0):
        leftbnk = floodplain_width/2. - chan_width/2.
        rightbnk = floodplain_width/2. + chan_width/2.
        # Left bank
        elev2 = elev1 + (channel_depth \
                - channel_depth/bankwidth*(x - leftbnk))*\
                (x>leftbnk)*(x < leftbnk + bankwidth)
        # Right bank
        elev2 = elev2 + (channel_depth \
                + channel_depth/bankwidth*(x - rightbnk))*\
                (x>rightbnk-bankwidth)*(x < rightbnk)
    
    if(bankwidth==0.0):
        elev2 = elev1

    return elev2


domain.set_quantity('elevation', topography, location='centroids') # Use function for elevation
domain.set_quantity('friction', man_n) # Constant friction
domain.set_quantity('stage', 0.7-1.0) # Use function for stage

#------------------------------------------------------------------------------
#
# Setup boundary conditions
#
#------------------------------------------------------------------------------


Br = anuga.Reflective_boundary(domain) # Solid reflective wall

def outflow_stage_boundary(t):
    return 0.7+0.3*numpy.sin(2.*numpy.pi*t/(30.*60.))-1.0 # Note elevation datum differs in hecras and ANUGA

Bout_tmss = anuga.Transmissive_n_momentum_zero_t_momentum_set_stage_boundary(domain, function = outflow_stage_boundary) 

domain.set_boundary({'left': Br, 
                     'right': Br, 
                     'top1': Br, 
                     'top2': Br, 
                     'bottom1': Br, 
                     'bottom2': Br, 
                     'chan_out': Bout_tmss, 
                     'chan_in': Br})

#------------------------------------------------------------------------------
# Produce a documentation of parameters
#------------------------------------------------------------------------------
if myid == 0:
    parameter_file=open('parameters.tex', 'w')
    parameter_file.write('\\begin{verbatim}\n')
    from pprint import pprint
    pprint(domain.get_algorithm_parameters(),parameter_file,indent=4)
    parameter_file.write('\\end{verbatim}\n')
    parameter_file.close()

#------------------------------------------------------------------------------
#
# Evolve system through time
#
#------------------------------------------------------------------------------

barrier()

for t in domain.evolve(yieldstep=10.0, finaltime=9.*60.*60.):
    if(myid==0):
        print(domain.timestepping_statistics())

barrier()

# Run sww merge
if( (myid==0) & (numprocs>1)):
    print('Merging sww files: ', numprocs, myid)
    anuga.utilities.sww_merge.sww_merge_parallel(outname,np=numprocs,verbose=True,delete_old=True)

barrier()
domain=None

# Clean up
finalize()

"""

Water flowing down a channel with a floodplain and a bridge

"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
import anuga
import numpy
from anuga.parallel.parallel_operator_factory import Inlet_operator, Boyd_box_operator, Internal_boundary_operator
from anuga.parallel import distribute, myid, numprocs, finalize, barrier
from anuga.structures.internal_boundary_functions import hecras_internal_boundary_function

args = anuga.get_args()
alg = args.alg
verbose = args.verbose



#------------------------------------------------------------------------------
# Useful parameters for controlling this case
#------------------------------------------------------------------------------

floodplain_length = 1000.0 # Model domain length
floodplain_width = 30.0 # Model domain width
floodplain_slope = 1./200.
chan_initial_depth = 0.25 # Initial depth of water in the channel
chan_bankfull_depth = 1.0 # Bankfull depth of the channel
chan_width = 10.0 # Bankfull width of the channel
bankwidth = 0.01 # Width of the bank regions -- note that these protrude into the channel
man_n = 0.045 # Manning's n
l0 = 5.000 # Length scale associated with triangle side length in channel (min_triangle area = 0.5*l0^2)

bridge_us = 520.
bridge_ds = 480.

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

breakLines = { # Put breaklines around edge of channel [think it helps the numerics] 
               'n1': [[floodplain_width/2.0 - chan_width/2., floodplain_length, -999.+0.*floodplain_length*floodplain_slope],
                      [floodplain_width/2.0 - chan_width/2., 0., -999.0]],
               'n3': [[floodplain_width/2.0 + chan_width/2.0, floodplain_length-l0*0, -999.+0.*floodplain_slope*floodplain_length],
                      [floodplain_width/2.0 + chan_width/2.0,     l0*0, -999.]],
                # Put breaklines for bridge, so the edges are 'sharp'
               'n4': [ [0., bridge_ds, -0.99], [floodplain_width, bridge_ds, -0.99]],
               'n5': [ [0., bridge_us, -0.99], [floodplain_width, bridge_us, -0.99]],
             }

regionPtAreas=[ [0.01, 0.01, 0.5*l0*l0*4], # lower-left
                [floodplain_width/2., 0.01, 0.5*l0*l0], # lower channel
                [floodplain_width-0.01, 0.01, 0.5*l0*l0*4], # lower_right
                [0.01, floodplain_length-0.01, 0.5*l0*l0*4], # upper_left
                [floodplain_width/2., floodplain_length-0.01, 0.5*l0*l0], # Upper_channel
                [floodplain_width-0.01, floodplain_length-0.01, 0.5*l0*l0*4], # upper_right
                # Bridge 
                [floodplain_width/2., floodplain_length/2., 0.5*l0*l0], 
                [floodplain_width-0.01, floodplain_length/2., 0.5*l0*l0], 
                [floodplain_width+0.01, floodplain_length/2., 0.5*l0*l0], 

              ]

if(myid==0):
    # Define domain with appropriate boundary conditions
    anuga.create_mesh_from_regions(boundary_polygon, 
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
                                   breaklines=breakLines.values(),
                                   regionPtArea=regionPtAreas,
                                   verbose=True)
    domain=anuga.create_domain_from_file('channel_floodplain1.msh')
    domain.set_name('channel_floodplain1') # Output name
    domain.set_flow_algorithm(alg)
else:
    domain=None

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

    # Add bridge
    elev2 = elev2*( (y<bridge_ds) + (y>bridge_us)) -1.0*( (y>=bridge_ds)*(y<=bridge_us))
    
    return elev2

#Function for stage
def stagetopo(x,y):
    return -y*floodplain_slope -chan_bankfull_depth + chan_initial_depth 

domain.set_quantity('elevation', topography, location='centroids') # Use function for elevation
domain.set_quantity('friction', man_n) # Constant friction
domain.set_quantity('stage', stagetopo) # Use function for stage

#domain.riverwallData.create_riverwalls(breakLines, 
#    default_riverwallPar = {'Qfactor': 0.65}, # Weir coefficient of 1.1 (0.65*default_value) 
#    output_dir='riverwall_text')

# Define inlet operator 
flow_in_yval = 10.0
if True:
    line1 = [ [floodplain_width/2. - chan_width/2., flow_in_yval],\
              [floodplain_width/2. + chan_width/2., flow_in_yval] \
              ]
    Qdata = [3., 3., 3., 3., 3., 3., 4., 5., 6., 6., 6., 6., 6., 7., 8., 9., 10., 10., 10., 10., 10.,\
           11., 12., 13., 14., 15., 15., 15., 15., 16., 17., 18., 19., 20., 21., 23., 25.,\
           30., 35.,40., 45., 50., 55., 60., 65., 70., 70., 70., 70.]+10*[70.]

    dtQdata=300.
    def Qin(t):
        t_hour=t/dtQdata # Used for time index for Qdata
        indL=numpy.floor(t_hour).astype(int)
        indU=numpy.ceil(t_hour).astype(int)
        w1=(t_hour-1.0*indL)
        w2=(1.0*indU-t_hour)
        Qout=Qdata[indL]*w2+Qdata[indU]*w1
        return Qout
    
    Inlet_operator(domain, line1, Qin)
    
# Set up bridge
bridge_in = [ [floodplain_width/2. - chan_width/2.+0.01, bridge_ds-0.01],\
              [floodplain_width/2. + chan_width/2.-0.01, bridge_ds-0.01] ]
bridge_out = [ [floodplain_width/2. - chan_width/2.+0.01, bridge_us+0.01],\
               [floodplain_width/2. + chan_width/2.-0.01, bridge_us+0.01] ]

hecras_discharge_function = hecras_internal_boundary_function(
    'hecras_bridge_table_high_deck.csv',
    allow_sign_reversal=True)

bridge = Internal_boundary_operator(
    domain,
    hecras_discharge_function,
    exchange_lines=[bridge_in, bridge_out],
    enquiry_gap=0.01,
    use_velocity_head=False,
    smoothing_timescale=30.0,
    logging=verbose)
    
#------------------------------------------------------------------------------
#
# Setup boundary conditions
#
#------------------------------------------------------------------------------

Br = anuga.Reflective_boundary(domain) # Solid reflective wall
Bt = anuga.Transmissive_boundary(domain) # Transmissive boundary


Bout_sub = anuga.Dirichlet_boundary( \
        [-floodplain_length*floodplain_slope - chan_bankfull_depth + \
        chan_initial_depth, 0., 0.]) #An outflow boundary for subcritical steady flow

def outflow_stage_boundary(t):
    return -floodplain_length*floodplain_slope \
            + chan_initial_depth - chan_bankfull_depth

Bout_tmss = anuga.shallow_water.boundaries.Transmissive_momentum_set_stage_boundary(domain, function = outflow_stage_boundary) 

domain.set_boundary({'left': Br, 
                     'right': Br, 
                     'top1': Bout_tmss, 
                     'top2': Bout_tmss, 
                     'bottom1': Br, 
                     'bottom2': Br, 
                     'chan_out': Bout_tmss, 
                     'chan_in': Br})

#------------------------------------------------------------------------------
# Produce a documentation of parameters
#------------------------------------------------------------------------------
from anuga.validation_utilities import save_parameters_tex
save_parameters_tex(domain)


#------------------------------------------------------------------------------
#
# Evolve system through time
#
#------------------------------------------------------------------------------

barrier()

for t in domain.evolve(yieldstep=10.0, finaltime=dtQdata * (len(Qdata) - 2)):
    if(myid == 0 and verbose):
        print(domain.timestepping_statistics())

    vol = domain.report_water_volume_statistics()

barrier()

# Run sww merge
if((myid == 0) & (numprocs > 1)):
    print('Merging sww files: ', numprocs, myid)
    anuga.utilities.sww_merge.sww_merge_parallel('channel_floodplain1',
                                                 np=numprocs,
                                                 verbose=True,
                                                 delete_old=True)

barrier()
finalize()

"""
Simple water flow example using ANUGA
Water flowing down a channel with a floodplain
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
# Import standard shallow water domain and standard boundaries.
import anuga
import numpy
from anuga import myid, finalize, distribute, numprocs, receive, send

#------------------------------------------------------------------------------
# Useful parameters for controlling this case
#------------------------------------------------------------------------------
args = anuga.get_args()
verbose = args.verbose
alg = args.alg

#--------------------------------
# See project file for definitions of these
# variables
#--------------------------------
from project import \
          floodplain_length, floodplain_width, floodplain_slope,\
          chan_initial_depth, chan_bankfull_depth, chan_width, bankwidth, \
          man_n, l0

l0 = 1.0

assert chan_width < floodplain_width, \
        'ERROR: Channel width is greater than floodplain width'

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
                    [floodplain_width/2. - chan_width/2., floodplain_length - l0],
                    [floodplain_width/2. + chan_width/2., floodplain_length - l0],
                    [floodplain_width/2. + chan_width/2., +l0]
                    ]
if myid == 0:
    # Define domain with appropriate boundary conditions
    domain = anuga.create_domain_from_regions( boundary_polygon, 
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
                                       verbose=verbose)

    domain.set_name('channel_floodplain') # Output name
    domain.set_flow_algorithm(alg)

    #------------------------------------------------------------------------------
    # Setup initial conditions
    #------------------------------------------------------------------------------

    # Function for topography
    def elevation(x, y):
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
        return elev2


    #Function for stage
    def stage(x,y):
        return -y*floodplain_slope -chan_bankfull_depth + chan_initial_depth 

    domain.set_quantity('elevation', elevation) # Use function for elevation
    domain.set_quantity('friction', man_n)      # Constant friction
    domain.set_quantity('stage', stage)         # Use function for stage

else:

    domain = None

#=====================================================
# Parallel Domain
#=====================================================
domain = distribute(domain)

#---------------------------------------------------------------------
# Define inlet operator 
#--------------------------------------------------------------------- 
flow_in_yval=5.0
line1 = [ [floodplain_width/2. - chan_width/2., flow_in_yval],\
          [floodplain_width/2. + chan_width/2., flow_in_yval] ]

Qin = 0.5*(floodplain_slope*(chan_width*chan_initial_depth)**2.*man_n**(-2.)\
            *chan_initial_depth**(4./3.) )**0.5
anuga.Inlet_operator(domain, line1, Qin)

if myid == 0 and verbose : print('Discharge in = ', Qin)

#---------------------------------------------------------------------
# Setup boundary conditions
#---------------------------------------------------------------------

Br = anuga.Reflective_boundary(domain) # Solid reflective wall

def outflow_stage_boundary(t):
    return -floodplain_length*floodplain_slope \
            + chan_initial_depth - chan_bankfull_depth

# Note that the outflow boundary may be slightly incorrect for the trapezoidal channel case, 
# or incorrect more generally if there are numerical problems. But, in the central regions of
# the channel, this shouldn't prevent us reaching steady, uniform flow.
Bout_tmss = anuga.Transmissive_momentum_set_stage_boundary(domain, function = outflow_stage_boundary) 

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
# Evolve system through time
#------------------------------------------------------------------------------
for t in domain.evolve(yieldstep=10.0, finaltime=1000.0):
    if myid == 0 and verbose: print(domain.timestepping_statistics())
    
#     if numpy.allclose(t,0.0):
#         exact_volume = domain.get_water_volume()
#     else:
#         exact_volume = exact_volume + Qin*10.0
#          
#     water_volume= domain.get_water_volume()
#     
#     
#     if myid == 0 and verbose: print(anuga.indent,'relative error water volume  exact_volume ', (water_volume - exact_volume)/exact_volume) 

domain.sww_merge(delete_old=True)


finalize()

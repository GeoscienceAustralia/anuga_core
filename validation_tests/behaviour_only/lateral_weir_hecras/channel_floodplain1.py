"""
3 channels connected by riverwalls -- example for comparison with HECRAS
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
floodplain_slope = 3./1000.
chan_initial_depth = 0.001 # Initial depth of water in the channel
chan_bankfull_depth = 0.0 # Bankfull depth of the channel
chan_width = 20.0 # Half- width of the central channel -- width of all others
bankwidth = 0.0 # Width of the bank regions -- note that these protrude into the channel
anuga2ras_friction_conversion = (0.5)**(2./3.)/( ( 0.5*10./(2.*0.5+10.) )**(2./3.))
man_n=0.03*anuga2ras_friction_conversion # Manning's n adjusted to be like hecras ( 
l0 = 10.000 # Length scale associated with triangle side length in channel (min_triangle area = 0.5*l0^2)

assert chan_width < floodplain_width, \
        ' ERROR: Channel width is greater than floodplain width'

assert bankwidth < chan_width/2., \
        'ERROR: The bank width must be less than half the channel width'

#------------------------------------------------------------------------------
# Setup computational domain
#------------------------------------------------------------------------------

# Define boundary polygon
boundary_polygon = [ [0.,0.], 
                     [0., floodplain_length], 
                     [floodplain_width/2. - chan_width/2., floodplain_length], 
                     [floodplain_width/2. + chan_width/2., floodplain_length], 
                     [floodplain_width, floodplain_length], 
                     [floodplain_width, 0.], 
                     [floodplain_width/2. + chan_width/2., 0.], 
                     [floodplain_width/2. - chan_width/2., 0.] 
                     ]

# Define riverwalls
liftWall=0.0 # Raise/lower the wall
HighWall=9.0e+10 # So high that the only gaps in the wall are as described below
riverWalls = { # Wall on the 'right' (facing downstream) can overflow from 201-399m. It is elevation 0.5 above the channel beds
               'n1': [[floodplain_width/2.0 - chan_width/2., floodplain_length, HighWall+liftWall],
                      [floodplain_width/2.0 - chan_width/2., (201.+198.), -0.7+liftWall],
                      [floodplain_width/2.0 - chan_width/2., 201., -0.1+liftWall],
                      [floodplain_width/2.0 - chan_width/2., 0., HighWall+liftWall]],
               # Wall on the 'left' can overflow from 301 - 499m, and is also elevation 0.5 above the channel beds
               'n3': [[floodplain_width/2.0 + chan_width/2.0, floodplain_length, HighWall+liftWall],
                      [floodplain_width/2.0 + chan_width/2.0, (301.+198.), -1.0+liftWall ],
                      [floodplain_width/2.0 + chan_width/2.0, 301., -0.4 +liftWall],
                      [floodplain_width/2.0 + chan_width/2.0, l0*0, HighWall+liftWall]],
             }

RAS_Qfac=1.1/1.7 # 
riverwall_Par={'n1':{'Qfactor':RAS_Qfac}, 'n3':{'Qfactor':RAS_Qfac}}

# RegionPtAreas control mesh resolutions in each region [split by riverWalls]
regionPtAreas=[ [0.01, 0.01, 0.5*l0*l0],
                [floodplain_width/2., 0.01, 0.5*l0*l0],
                [floodplain_width-0.01, floodplain_length-0.01, 0.5*l0*l0]
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
                                       breaklines=riverWalls.values(),
                                       regionPtArea=regionPtAreas,
                                       verbose=verbose)
    domain=anuga.create_domain_from_file('channel_floodplain1.msh')
    domain.set_name('channel_floodplain1') # Output name
    domain.set_flow_algorithm(alg)
else:
    domain=None

barrier()
domain=distribute(domain)
barrier()

# If riverwalls don't store vertices uniquely, ANUGA viewer output will
# look as though water extends 1 triangle behind the riverwall
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

    return elev2

#Function for stage
def stagetopo(x,y):
    return -y*floodplain_slope -chan_bankfull_depth + chan_initial_depth 

domain.set_quantity('elevation', topography, location='centroids') # Use function for elevation
domain.set_quantity('friction', man_n) # Constant friction
domain.set_quantity('stage', stagetopo) # Use function for stage
# Add the riverwalls
domain.riverwallData.create_riverwalls(riverWalls, riverwallPar=riverwall_Par, output_dir='riverwall_points',verbose=verbose)

# Define inlet operator 
flow_in_yval=10.0
if True:
    line1 = [ [floodplain_width/2. - chan_width/2.+0.01, flow_in_yval],\
              [floodplain_width/2. + chan_width/2.-0.01, flow_in_yval] \
              ]
    
    Qdata = [1., 1., 5., 9., 13., 17., 21.]+18*[21.]

    dtQdata=3600.
    def Qin(t):
        t_hour=t/dtQdata # Used for time index for Qdata
        indL=numpy.floor(t_hour).astype(int)
        indU=numpy.ceil(t_hour).astype(int)
        w1=(t_hour-1.0*indL)
        w2=(1.0*indU-t_hour)
        Qout=Qdata[indL]*w2+Qdata[indU]*w1
        return Qout
    
    Inlet_operator(domain, line1, Qin)
   
    if(verbose): 
        print('Discharge in = ', Qin) #,'Velocity at inlet should be ', Qin/(chan_width*chan_initial_depth), \
            #'for rectangular channels, and approximately that for trapezoidal ones'

    # Add 'pilot discharge' for left/right areas.
    # This is needed for comparison with hecras [since hecras channels cannot dry]
    Inlet_operator(domain, [[0.+1.0e-03, flow_in_yval], [10.-1.0e-03, flow_in_yval] ], 0.1)
    Inlet_operator(domain, [[30.+1.0e-03,flow_in_yval], [40.-1.0e-03, flow_in_yval] ], 0.1)


#------------------------------------------------------------------------------
#
# Setup boundary conditions
#
#------------------------------------------------------------------------------


Br = anuga.Reflective_boundary(domain) # Solid reflective wall
Bt = anuga.Transmissive_boundary(domain)
def outflow_stage_boundary(t):
    return -floodplain_length*floodplain_slope \
            + chan_initial_depth - chan_bankfull_depth+0.2

Bout_tmss = anuga.Transmissive_n_momentum_zero_t_momentum_set_stage_boundary(domain, function = outflow_stage_boundary) 

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

for t in domain.evolve(yieldstep=60.0, finaltime=dtQdata * (len(Qdata) - 2)):
    if(myid==0 and verbose):
        print(domain.timestepping_statistics())
        #xx=domain.quantities['ymomentum'].centroid_values
        #dd=(domain.quantities['stage'].centroid_values - domain.quantities['elevation'].centroid_values)
        #dd=dd*(dd>0.)
        #tmp = xx/(dd+1.0e-06)*(dd>0.0)
        #print tmp.max(), tmp.argmax(), tmp.min(),  tmp.argmin()
        #print domain.max_speed.max(), domain.max_speed.argmax()

domain.sww_merge(delete_old=True)

finalize()

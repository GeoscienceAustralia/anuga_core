import os.path
import sys

#from anuga.utilities.system_tools import get_pathname_from_package
#from anuga.geometry.polygon_function import Polygon_function

#from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
#from anuga.abstract_2d_finite_volumes.quantity import Quantity

import anuga
import time

from math import pi
import numpy as num



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



length = 40.
width = 15.

dx = dy = 0.5          # Resolution: Length of subdivisions on both axes


if anuga.myid == 0:

    domain = anuga.rectangular_cross_domain(int(length/dx),
                                            int(width/dy),
                                            len1=length,
                                            len2=width)
    #domain.set_name('anuga_gate_operator')     # Output name
    domain.set_name() 
    domain.set_flow_algorithm('DE0')
    domain.set_store_vertices_uniquely(True)
    domain.set_quantity('elevation', topography)
    domain.set_quantity('friction', 0.01)         # Constant friction
    domain.set_quantity('stage',
                        expression='elevation')   # Dry initial condition
else:

    domain = None


domain = anuga.distribute(domain, verbose = True)

#domain.set_store_vertices_uniquely(False)


gate = anuga.Boyd_box_operator(domain,
                            end_points=[[9.0, 2.5],[13.0, 2.5]],
                            losses=1.5,
                            width=1.5,
                            height = 0.0001,
                            apron=5.0,
                            use_momentum_jet=True,
                            use_velocity_head=False,
                            manning=0.013,
                            verbose=False)

# Close the gate
if gate is not None:
    gate.set_culvert_height(0.0)

line = [[0.0, 5.0], [0.0, 10.0]]
Q = 1.0
in0 = anuga.Inlet_operator(domain, line, Q)




##-----------------------------------------------------------------------
## Setup boundary conditions
##-----------------------------------------------------------------------

## Inflow based on Flow Depth and Approaching Momentum
Br = anuga.Reflective_boundary(domain)              # Solid reflective wall

domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})


##-----------------------------------------------------------------------
## Evolve system through time
##-----------------------------------------------------------------------
anuga.barrier()

t0 = time.time()

for t in domain.evolve(yieldstep = 1.0, finaltime = 50):

    if anuga.myid == 0:
        print(80*'=')
        domain.write_time()

    #================================================
    # gate operator lives on a number of processors
    #
    # Need to run parallel gate operations on all
    # processors where gate lives.
    # But true result only on master processor.
    #================================================
    if gate is not None:
        [s0, s1] = gate.get_enquiry_stages()
        [d0, d1] = gate.get_enquiry_depths()
        [e0, e1] = gate.get_enquiry_elevations()
        [i0, i1] = gate.get_enquiry_invert_elevations()
        [w0, w1] = gate.get_enquiry_water_depths()

        output = gate.discharge_routine()

        if anuga.myid == gate.get_master_proc():
            print('myid, time ', anuga.myid, t, s0,s1)
            print('myid, time ', anuga.myid, t, d0,d1)
            print('myid, time ', anuga.myid, t, e0,e1)
            print('myid, time ', anuga.myid, t, i0,i1)
            print('myid, time ', anuga.myid, t, w0,w1)
            print('myid, time ', anuga.myid, t, output)


            if d0 > 0.2: gate.set_culvert_height(10.0)


    anuga.barrier()


domain.sww_merge(delete_old=True)




anuga.barrier()

for p in range(anuga.numprocs):
    if anuga.myid == p:
        print('Processor %g ' % anuga.myid)
        print('That took %.2f seconds' %(time.time()-t0))
        print('Communication time %.2f seconds'%domain.communication_time)
        print('Reduction Communication time %.2f seconds'%domain.communication_reduce_time)
        print('Broadcast time %.2f seconds'%domain.communication_broadcast_time)
    else:
        pass

    anuga.barrier()


anuga.finalize()

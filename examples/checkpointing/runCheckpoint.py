
"""Run parallel shallow water domain.

with check pointing
"""


# ------------------------------------------------------------------------------
# Import necessary modules
# ------------------------------------------------------------------------------

import os
import sys
import time
import numpy as num

# ------------------------
# ANUGA Modules
# ------------------------
from anuga import Reflective_boundary
from anuga import Transmissive_n_momentum_zero_t_momentum_set_stage_boundary

from anuga import create_domain_from_file

from anuga import distribute, myid, numprocs, finalize, barrier


# --------------------------------------------------------------------------
# Setup parameters
# --------------------------------------------------------------------------

verbose = False
domain_name = 'checkpoint'
checkpoint_dir = 'CHECKPOINTS'
finaltime = 5000.0
useCheckpointing = True

mesh_filename = "merimbula_10785_1.tsh"
x0 = 756000.0
x1 = 756500.0
yieldstep = 50


# --------------------------------------------------------------------------
# Setup procedures
# --------------------------------------------------------------------------
class Set_Stage:
    """Set an initial condition with constant water height, for x0<x<x1
    """

    def __init__(self, x0=0.25, x1=0.5, h=1.0):
        self.x0 = x0
        self.x1 = x1
        self.h = h

    def __call__(self, x, y):
        return self.h*((x > self.x0) & (x < self.x1))


class Set_Elevation:
    """Set an elevation
    """

    def __init__(self, h=1.0):
        self.x0 = x0
        self.x1 = x1
        self.h = h

    def __call__(self, x, y):
        return x/self.h


def wave(t):
    from math import sin
    return sin(t/200)


# ------------------------------------------------------------------------------
# Use a try statement to read in previous checkpoint file and if not possible
# just go ahead as normal and produce domain as usual.
#
# Though in the part of the code where you create the domain as normal,
# remember to turn on checkpointing via
# domain.set_checkpointing(checkpoint_time = 5) (see code below)
# ------------------------------------------------------------------------------
try:
    from anuga import load_checkpoint_file

    domain = load_checkpoint_file(domain_name=domain_name,
                                  checkpoint_dir=checkpoint_dir)
except:
    # --------------------------------------------------------------------------
    # Normal Setup of Domain on processor 0
    # --------------------------------------------------------------------------
    if myid == 0:
        domain = create_domain_from_file(mesh_filename)
        domain.set_quantity('stage', Set_Stage(x0, x1, 0.0))

        domain.set_name(domain_name)
        domain.set_store(True)

        domain.set_store_vertices_smoothly(False)
    else:
        domain = None

    # --------------------------------------------------------------------------
    # Distribute sequential domain on processor 0 to other processors
    # --------------------------------------------------------------------------

    if myid == 0 and verbose:
        print('DISTRIBUTING DOMAIN')
    domain = distribute(domain, verbose=verbose)

    # --------------------------------------------------------------------------
    # On all processors, setup evolve parameters for domains on all processors
    # (all called "domain"
    # --------------------------------------------------------------------------

    domain.set_flow_algorithm('DE1')
    domain.set_CFL(2.0)
    domain.set_store_centroids()

    domain.set_quantities_to_be_stored({'elevation': 1,
                                        'friction': 1,
                                        'stage': 2,
                                        'xmomentum': 2,
                                        'ymomentum': 2})

    # ------------------------------------------------------------------------------
    # Setup boundary conditions
    # This must currently happen *after* domain has been distributed
    # ------------------------------------------------------------------------------
    Br = Reflective_boundary(domain)      # Solid reflective wall

    Bts = Transmissive_n_momentum_zero_t_momentum_set_stage_boundary(
        domain, wave)

    domain.set_boundary({'outflow': Br, 'inflow': Br, 'inner': Br,
                         'exterior': Br, 'open': Bts})

    # -----------------------------------------------------------------------------
    # Turn on checkpointing every 5 sec (just for testing, more reasonable to
    # set to say 15 minutes = 15*60 sec)
    # Only set this on process 0, ie within a if myid == 0: structure
    # -----------------------------------------------------------------------------
    if useCheckpointing:
        domain.set_checkpointing(
            checkpoint_time=5, checkpoint_dir=checkpoint_dir)


# ------------------------------------------------------------------------------
# Evolution
# ------------------------------------------------------------------------------
if myid == 0 and verbose:
    print('EVOLVE')

barrier()
t0 = time.time()


for t in domain.evolve(yieldstep=yieldstep, finaltime=finaltime):
    if myid == 0:
        domain.write_time()


barrier()

for p in range(numprocs):
    if myid == p:
        print('Processor %g ' % myid)
        print('That took %.2f seconds' % (time.time()-t0))
        print('Communication time %.2f seconds' % domain.communication_time)
        print('Reduction Communication time %.2f seconds'
              % domain.communication_reduce_time)
        print('Broadcast time %.2f seconds' %
              domain.communication_broadcast_time)
    else:
        pass

    barrier()


# --------------------------------------------------
# Merge the individual sww files into one file
# But don't delete the sub domain sww files
# --------------------------------------------------
domain.sww_merge()


finalize()

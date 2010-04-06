
"""Run parallel shallow water domain.
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------

import unittest
import os
import sys
import pypar

import numpy as num




from anuga.utilities.numerical_tools import ensure_numeric
from anuga.utilities.util_ext        import double_precision
from anuga.utilities.norms           import l1_norm, l2_norm, linf_norm

from anuga.interface import Domain
from anuga.interface import Reflective_boundary
from anuga.interface import Dirichlet_boundary
from anuga.interface import Time_boundary
from anuga.interface import Transmissive_boundary

from anuga.interface import rectangular_cross
from anuga.interface import create_domain_from_file


from anuga_parallel.interface import distribute, myid, numprocs, finalize


#--------------------------------------------------------------------------
# Setup parameters
#--------------------------------------------------------------------------

mesh_filename = "merimbula_10785_1.tsh"
#mesh_filename = "test-100.tsh"
yieldstep = 1
finaltime = 20
quantity = 'stage'
verbose = True

#--------------------------------------------------------------------------
# Setup procedures
#--------------------------------------------------------------------------
class Set_Stage:
    """Set an initial condition with constant water height, for x<x0
    """

    def __init__(self, x0=0.25, x1=0.5, h=1.0):
        self.x0 = x0
        self.x1 = x1
        self.h  = h

    def __call__(self, x, y):
        return self.h*((x>self.x0)&(x<self.x1))

#--------------------------------------------------------------------------
# Setup Domain only on processor 0
#--------------------------------------------------------------------------
if myid == 0:
    domain = create_domain_from_file(mesh_filename)
    domain.set_quantity('stage', Set_Stage(756000.0, 756500.0, 2.0))
else:
    domain = None

#--------------------------------------------------------------------------
# Distribute sequential domain on processor 0 to other processors
#--------------------------------------------------------------------------

if myid == 0 and verbose: print 'DISTRIBUTING DOMAIN'
domain = distribute(domain)

#------------------------------------------------------------------------------
# Setup boundary conditions
# This must currently happen *after* domain has been distributed
#------------------------------------------------------------------------------
Br = Reflective_boundary(domain)      # Solid reflective wall

domain.set_boundary({'outflow' :Br, 'inflow' :Br, 'inner' :Br, 'exterior' :Br, 'open' :Br})

#------------------------------------------------------------------------------
# Evolution
#------------------------------------------------------------------------------
if myid == 0 and verbose: print 'EVOLVE'
        
for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):
    if myid == 0:
        domain.write_time()

finalize()





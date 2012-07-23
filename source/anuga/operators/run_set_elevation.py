"""Simple test of setting stage in a circular region
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------

import sys
import anuga


from math import cos
from numpy import zeros, float, where
import numpy
from time import localtime, strftime, gmtime



#------------------------------------------------------------------------------
# Setup domain
#------------------------------------------------------------------------------
dx = 100.
dy = dx
L = 10000.
W = L
#W = dx

# structured mesh
domain = anuga.rectangular_cross_domain(int(L/dx), int(W/dy), L, W, (-L/2.0, -W/2.0))


output_file = 'set_bed'
domain.set_name(output_file)                

#------------------------------------------------------------------------------
# Setup Algorithm
#------------------------------------------------------------------------------
domain.set_flow_algorithm(2.0)

#------------------------------------------------------------------------------
# Setup initial conditions
#------------------------------------------------------------------------------
h0 = 1000.0

domain.set_quantity('elevation',0.0)
domain.set_quantity('friction', 0.0)
domain.add_quantity('stage', h0)

#-----------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
Br = anuga.Reflective_boundary(domain)      # Solid reflective wall

# Associate boundary tags with boundary objects
domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

#------------------------------------------------------------------------------
# Setup Operators
#------------------------------------------------------------------------------
from anuga.operators.set_bed_operators import Circular_set_bed_operator

import math

bed1 = lambda t: 20.0 * math.sin(t/3.0)
cop1 = Circular_set_bed_operator(domain, bed=bed1, center=(0.0, 0.0), radius=100.0 )

#stage2 = lambda t: h0 + 30.0 * math.sin(t/6.0)
#cop2 = Circular_set_stage_operator(domain, stage=stage2, center=(2000.0, 1000.0), radius=100.0 )

#print cop1.statistics()
#print cop2.statistics()

#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------

for t in domain.evolve(yieldstep = 1.0, finaltime = 20.0):
    #print domain.timestepping_statistics(track_speeds=True)
    domain.print_timestepping_statistics()
    domain.print_operator_timestepping_statistics()



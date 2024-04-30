"""Simple test of setting stage in a circular region
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------

import sys
import anuga


from math import cos
from numpy import zeros, where
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


domain.set_name()                

#------------------------------------------------------------------------------
# Setup Algorithm (Usual choices DE0 and DE1)
#------------------------------------------------------------------------------
domain.set_flow_algorithm('DE0')

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
from anuga.operators.set_stage_operator import Set_stage_operator

import math

stage1 = lambda t: h0 + 20.0 * math.sin(t/3.0)
#cop1 = Set_stage_operator(domain, stage=stage1, center=(0.0, 0.0), radius=100.0 )

stage2 = lambda t: h0 + 30.0# * math.sin(t/3.0)
cop2 = Set_stage_operator(domain, stage=stage2, center=(2000.0, 1000.0), radius=100.0 )

#print cop1.statistics()
#print cop2.statistics()

#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------

## stage = domain.get_quantity('stage')
## id = domain.get_triangle_containing_point([2000., 1150.])
## times = []
## stages = []

for t in domain.evolve(yieldstep = 1.0, finaltime = 200.0):
    #print domain.timestepping_statistics(track_speeds=True)
    domain.print_timestepping_statistics()
    domain.print_operator_timestepping_statistics()
    #times.append(t)
    #stages.append(stage.centroid_values[id])

## import pylab as pl

## pl.plot(times, stages)
## pl.show()

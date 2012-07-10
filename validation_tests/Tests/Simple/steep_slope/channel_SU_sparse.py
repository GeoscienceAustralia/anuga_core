"""Simple water flow example using ANUGA
Water flowing down a channel
"""
import sys

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
import anuga
from anuga import rectangular_cross as rectangular_cross
from anuga.structures.inlet_operator import Inlet_operator
from anuga import Domain

#------------------------------------------------------------------------------
# Setup computational domain
#------------------------------------------------------------------------------
points, vertices, boundary = rectangular_cross(10, 10,
len1=100.0, len2=100.0) # Mesh
domain = Domain(points, vertices, boundary) # Create domain
domain.set_name('channel_SU_2_v2') # Output name

#------------------------------------------------------------------------------
# Setup Algorithm, either using command line arguments
# or override manually yourself
#------------------------------------------------------------------------------
from anuga.utilities.validations import validation_parse
alg, cfl = validation_parse()
domain.set_flow_algorithm(alg)
domain.set_CFL(cfl)

#------------------------------------------------------------------------------
# Setup initial conditions
#------------------------------------------------------------------------------
def topography(x, y):
	return -x/10. #+abs(y-2.5)/10. # linear bed slope

def stagetopo(x,y):
    stg= -x/10. +0.06309625 # Constant depth
    #topo=topography(x,y)
    return stg#*(stg>topo) + topo*(stg<=topo)

line1=[ [1.,0.], [1., 100.] ]
Qin=20.
Inlet_operator(domain, line1, Qin)

domain.set_quantity('elevation', topography) # Use function for elevation
domain.set_quantity('friction', 0.03) # Constant friction
domain.set_quantity('stage', stagetopo)

#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
#Bi = anuga.Dirichlet_boundary([0.06309625, 0.00, 0.00]) # Inflow for steady uniform flow with a depth integrated velocity of 0.1

Br = anuga.Reflective_boundary(domain) # Solid reflective wall
#Bo= anuga.Dirichlet_boundary([-9.9369037,0.0,0]) # Outflow for steady uniform flow with a depth integrated velocity of 0.1

def constant_water(t):
    return -9.936903

#Bo= anuga.shallow_water.boundaries.Transmissive_n_momentum_zero_t_momentum_set_stage_boundary(domain, constant_water) # Outflow for steady uniform flow with a depth integrated velocity of 0.1
Bo = anuga.Transmissive_boundary(domain)
domain.set_boundary({'left': Br, 'right': Bo, 'top': Br, 'bottom': Br})
#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------
for t in domain.evolve(yieldstep=4.0, finaltime=400.0):
	print domain.timestepping_statistics()

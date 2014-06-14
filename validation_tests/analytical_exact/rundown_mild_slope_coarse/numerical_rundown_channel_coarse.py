"""
Simple water flow example using ANUGA: Water flowing down a channel.
It was called "steep_slope" in an old validation test.
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
points, vertices, boundary = rectangular_cross(14, 10, len1=140.0, len2=100.0)
domain = Domain(points, vertices, boundary) # Create domain
domain.set_name('channel') # Output name

#------------------------------------------------------------------------------
# Setup Algorithm, either using command line arguments
# or override manually yourself
#------------------------------------------------------------------------------
from anuga.utilities.argparsing import parse_standard_args
alg, cfl = parse_standard_args()
domain.set_flow_algorithm(alg)
#domain.set_CFL(cfl)

#------------------------------------------------------------------------------
# Setup initial conditions
#------------------------------------------------------------------------------
Qin=0.1
#fluxin=Qin/100. #The momentum flux at the upstream boundary ( = discharge / width)
#uana= 2.15843634571 # analytical Xvelocity
#dana= 0.0926596702273 # analytical water depth

def topography(x, y):
	return -x/10. # linear bed slope

def init_stage(x,y):
    stg= -x/10.+0.10 # Constant depth: 10 cm.
    return stg
#line0=[ [0.,0.], [0., 100.] ]
#Uin=[uana, 0.0]
#Inlet_operator(domain, line0, Q=Qin, velocity=Uin)

line1=[ [0.,0.], [0., 100.] ]
#Qin=0.1
Inlet_operator(domain, line1,Qin)


domain.set_quantity('elevation', topography) # Use function for elevation
domain.set_quantity('friction', 0.03) # Constant friction
domain.set_quantity('stage', init_stage)

#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
Bt = anuga.Transmissive_boundary(domain)
#BdIN = anuga.Dirichlet_boundary([dana, fluxin, 0.0])
#BdOUT = anuga.Dirichlet_boundary([dana-10., fluxin, 0.0])
Br = anuga.Reflective_boundary(domain) # Solid reflective wall
domain.set_boundary({'left': Br, 'right': Bt, 'top': Br, 'bottom': Br})


#------------------------------------------------------------------------------
# Produce a documentation of parameters
#------------------------------------------------------------------------------
parameter_file=open('parameters.tex', 'w')
parameter_file.write('\\begin{verbatim}\n')
from pprint import pprint
pprint(domain.get_algorithm_parameters(),parameter_file,indent=4)
parameter_file.write('\\end{verbatim}\n')
parameter_file.close()

#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------
for t in domain.evolve(yieldstep=10.0, finaltime=2000.0):
    print domain.timestepping_statistics()
    print (domain.areas*(domain.quantities['stage'].centroid_values - domain.quantities['elevation'].centroid_values)).sum()
    s3 = domain.get_flow_through_cross_section([[30., 0.0], [30., 100.]])
    s4 = domain.get_flow_through_cross_section([[32., 0.0], [32., 100.]])
    s5 = domain.get_flow_through_cross_section([[34., 0.0], [34., 100.]])
    s2 = domain.get_flow_through_cross_section([[45., 0.0], [45., 100.]])
    s1 = domain.get_flow_through_cross_section([[53., 0.0], [53., 100.]])
    s0 = domain.get_flow_through_cross_section([[60., 0.0], [60., 100.]])
    print 'Xsectional flow:', s0, s1, s2, s3, s4, s5

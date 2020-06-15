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
from anuga import myid, finalize, distribute

#===============================================================================
# Setup flow conditions
#===============================================================================
Qin=20.
fluxin=Qin/100. #The momentum flux at the upstream boundary ( = discharge / width)
uana= 2.15843634571 # analytical Xvelocity
dana= 0.0926596702273 # analytical water depth
slope = -0.1
mannings = 0.03

args = anuga.get_args()
alg = args.alg
verbose = args.verbose

#------------------------------------------------------------------------------
# Setup sequential computational domain
#------------------------------------------------------------------------------
if myid == 0:
	points, vertices, boundary = rectangular_cross(50, 50, len1=100.0, len2=100.0)
	domain = Domain(points, vertices, boundary) # Create domain
	domain.set_name('channel') # Output name

	domain.set_flow_algorithm(alg)
	domain.set_store_centroids(True)
	
	#------------------------------------------------------------------------------
	# Setup initial conditions
	#------------------------------------------------------------------------------
	
	def topography(x, y):
		return x*slope # linear bed slope
	
	def init_stage(x,y):
		stg= x*slope+0.01 # Constant depth: 1 cm.
		return stg
	
	domain.set_quantity('elevation', topography) # Use function for elevation
	domain.set_quantity('friction', mannings) # Constant friction
	domain.set_quantity('stage', init_stage)

else:
	
	domain = None

#===============================================================================
# Parallel Domain
#===============================================================================
domain = distribute(domain)

#domain.set_store_centroids(True)
#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
Bt = anuga.Transmissive_boundary(domain)
BdIN = anuga.Dirichlet_boundary([dana, fluxin, 0.0])
#BdOUT = anuga.Dirichlet_boundary([dana-10., fluxin, 0.0])
Br = anuga.Reflective_boundary(domain) # Solid reflective wall
domain.set_boundary({'left': BdIN, 'right': Bt, 'top': Br, 'bottom': Br})


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
for t in domain.evolve(yieldstep=2.0, finaltime=200.0):
	if myid ==0 and verbose: print(domain.timestepping_statistics())


domain.sww_merge(delete_old=True)

finalize()

"""
Simple water flow example using ANUGA: Water flowing down a channel.
It was called "steep_slope" in an old validation test.
"""
import sys

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
import anuga
from anuga import rectangular_cross
from anuga import Inlet_operator
from anuga import Domain
from anuga import myid, finalize, distribute

Qin = 0.1
fluxin=Qin/100. #The momentum flux at the upstream boundary ( = discharge / width)
mann=0.03 # Manning's coef
bedslope=-0.1
uana= ( mann**(-2.)*abs(bedslope)*fluxin**(4./3.) )**(3./10.) # Velocity
dana= fluxin/uana # Depth



args = anuga.get_args()
alg = args.alg
verbose = args.verbose

#------------------------------------------------------------------------------
# Setup sequential computational domain
#------------------------------------------------------------------------------
if myid == 0:
	points, vertices, boundary = rectangular_cross(40, 10, len1=400.0, len2=100.0)
	domain = Domain(points, vertices, boundary) # Create domain
	domain.set_name('channel') # Output name
	
	domain.set_flow_algorithm(alg)

	
	#------------------------------------------------------------------------------
	# Setup initial conditions
	#------------------------------------------------------------------------------

	
	def topography(x, y):
		return -x/10. # linear bed slope
	
	def init_stage(x,y):
		stg= -x/10.+0.004 # Constant depth: 10 cm.
		return stg
	#line0=[ [0.,0.], [0., 100.] ]
	#Uin=[uana, 0.0]
	#Inlet_operator(domain, line0, Q=Qin, velocity=Uin)
	

	
	domain.set_quantity('elevation', topography) # Use function for elevation
	domain.set_quantity('friction', mann) # Constant friction
	domain.set_quantity('stage', init_stage)
	domain.set_quantity('xmomentum', dana*uana)
else:
	
	domain = None
	
#===========================================================================
# Create Parallel Domain
#===========================================================================
domain = distribute(domain)

#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
#
# This one can have outflow boundary issues -- a good approach is just to use a reflective
# boundary, and allow water to 'pool' at the bottom of the domain.
#
#Bt = anuga.Transmissive_boundary(domain)
#Bts = anuga.Transmissive_momentum_set_stage_boundary(domain, dana-160.0)
#Bts = anuga.Transmissive_n_momentum_zero_t_momentum_set_stage_boundary(domain, lambda t: dana-40.0)
##BdIN = anuga.Dirichlet_boundary([dana, fluxin, 0.0])
#BdOUT = anuga.Dirichlet_boundary([dana-40., dana*uana, 0.0])

print(dana-40.)

Br = anuga.Reflective_boundary(domain) # Solid reflective wall
domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})


line1=[ [0.0, 0.], [0.0, 100.] ]
Qin=0.1
inlet = Inlet_operator(domain, line1, Q = Qin)


#if inlet: print inlet.statistics()


stage = domain.quantities['stage']
elev  = domain.quantities['elevation']

print((stage-elev).get_integral())



	


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
for t in domain.evolve(yieldstep=10.0, finaltime=3000.0):
	if myid == 0 and verbose: print(domain.timestepping_statistics())
	#print (stage-elev).get_integral()
    #print (domain.areas*(domain.quantities['stage'].centroid_values - domain.quantities['elevation'].centroid_values)).sum()
    #s3 = domain.get_flow_through_cross_section([[30., 0.0], [30., 100.]])
    #s4 = domain.get_flow_through_cross_section([[32., 0.0], [32., 100.]])
    #s5 = domain.get_flow_through_cross_section([[34., 0.0], [34., 100.]])
    #s2 = domain.get_flow_through_cross_section([[45., 0.0], [45., 100.]])
    #s1 = domain.get_flow_through_cross_section([[53., 0.0], [53., 100.]])
    #s0 = domain.get_flow_through_cross_section([[60., 0.0], [60., 100.]])
    #print 'Xsectional flow:', s0, s1, s2, s3, s4, s5


domain.sww_merge(delete_old=True)

finalize()

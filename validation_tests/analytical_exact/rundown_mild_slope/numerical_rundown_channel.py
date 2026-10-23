"""
Simple water flow example using ANUGA: Water flowing down a channel.
It was called "steep_slope" in an old validation test.

Two cases with the same discharge (q = 0.2 m^2/s) on the same 1:10 slope:

  stable    Manning n = 0.06, Froude 1.21 -- the steady uniform flow is
            stable, and the scheme should hold it (100 m x 100 m).
  rollwave  Manning n = 0.03, Froude 2.26 -- above the roll-wave threshold
            of 1.5 for Manning friction, so the steady uniform flow is
            unstable. Any disturbance (here the one the inflow boundary
            makes) grows down the slope into a train of roll waves. The
            channel is 200 m long so the waves have room to grow and then
            saturate. The inflow discharge carries a 1% sinusoidal
            disturbance (period 3.3 s) so that every solver is given the same
            seed: without one, a scheme that settles to an exact discrete
            steady state grows no roll waves at all, however weakly it damps
            them.

RUNDOWN_CASES (default "stable,rollwave") picks the cases to run, and
RUNDOWN_TAG is appended to the output names; compare_solvers.py uses both to
run the roll-wave case once per solver.
"""
import math
import os
import sys

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
import anuga
from anuga import rectangular_cross as rectangular_cross
from anuga import Domain
from anuga import myid, finalize, distribute, barrier

#===============================================================================
# Setup flow conditions
#===============================================================================
Qin=20.
fluxin=Qin/100. #The momentum flux at the upstream boundary ( = discharge / width)
slope = -0.1

# name: (Manning n, length, width, finaltime)
cases = {
	'stable':   (0.06, 100.0, 100.0, 200.0),
	'rollwave': (0.03, 200.0,  20.0, 400.0),
}
run_cases = os.environ.get('RUNDOWN_CASES', 'stable,rollwave').split(',')
tag = os.environ.get('RUNDOWN_TAG', '')

# Inflow disturbance for the roll-wave case. The period is deliberately not a
# multiple of the 2 s yieldstep, or the stored frames would sample the
# disturbance at a fixed phase and alias the travelling waves.
seed_amplitude = 0.01
seed_period = 3.3
fine_window = 30.0

args = anuga.get_args()
alg = args.alg
verbose = args.verbose

def analytic(mannings):
	"""Steady uniform flow: friction slope = bed slope, depth*velocity = fluxin."""
	uana = ( mannings**(-2.)*abs(slope)*fluxin**(4./3.) )**(3./10.)
	return uana, fluxin/uana

for name, (mannings, length, width, finaltime) in cases.items():
	if name not in run_cases:
		continue
	uana, dana = analytic(mannings)

	#------------------------------------------------------------------------------
	# Setup sequential computational domain
	#------------------------------------------------------------------------------
	if myid == 0:
		points, vertices, boundary = rectangular_cross(int(length/2), int(width/2),
		                                               len1=length, len2=width)
		domain = Domain(points, vertices, boundary) # Create domain
		domain.set_name('channel_' + name + tag) # Output name

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

	#------------------------------------------------------------------------------
	# Setup boundary conditions
	#------------------------------------------------------------------------------
	Bt = anuga.Transmissive_boundary(domain)
	if name == 'rollwave':
		BdIN = anuga.Time_boundary(domain, function=lambda t, dana=dana:
			[dana, fluxin*(1.0 + seed_amplitude*math.sin(2.0*math.pi*t/seed_period)), 0.0])
	else:
		BdIN = anuga.Dirichlet_boundary([dana, fluxin, 0.0])
	Br = anuga.Reflective_boundary(domain) # Solid reflective wall
	domain.set_boundary({'left': BdIN, 'right': Bt, 'top': Br, 'bottom': Br})


	#------------------------------------------------------------------------------
	# Produce a documentation of parameters
	#------------------------------------------------------------------------------
	if myid == 0 and name == 'stable':
		parameter_file=open('parameters.tex', 'w')
		parameter_file.write('\\begin{verbatim}\n')
		from pprint import pprint
		pprint(domain.get_algorithm_parameters(),parameter_file,indent=4)
		parameter_file.write('\\end{verbatim}\n')
		parameter_file.close()

	#------------------------------------------------------------------------------
	# Evolve system through time
	#------------------------------------------------------------------------------
	if myid == 0 and verbose:
		print('Case %s: n = %g, Froude %.2f' % (name, mannings, uana/(9.8*dana)**0.5))
	# The roll-wave case stores its last fine_window seconds every 0.5 s so
	# that individual waves can be followed in space and time (at 2 s they
	# move several wavelengths between frames).
	coarse_end = finaltime - fine_window if name == 'rollwave' else finaltime
	for t in domain.evolve(yieldstep=2.0, finaltime=coarse_end):
		if myid ==0 and verbose: print(domain.timestepping_statistics())
	if coarse_end < finaltime:
		for t in domain.evolve(yieldstep=0.5, finaltime=finaltime):
			if myid ==0 and verbose: print(domain.timestepping_statistics())


	domain.sww_merge(delete_old=True)
	barrier()

finalize()

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
from anuga.shallow_water import Domain
from anuga.shallow_water import Reflective_boundary
from anuga.shallow_water import Dirichlet_boundary


def create_test_sww(sww_name):

	#------------------------------------------------------------------------------
	# Setup computational domain
	#------------------------------------------------------------------------------
	points, vertices, boundary = rectangular_cross(50, 50,
												   len1=50.0, len2=50.0) # Mesh

	domain = Domain(points, vertices, boundary)  # Create domain
	domain.set_name(sww_name)                  # Output name

	#------------------------------------------------------------------------------
	# Setup initial conditions
	#------------------------------------------------------------------------------
	def topography(x, y):
		return -x/20                             # linear bed slope

	domain.set_quantity('elevation', topography) # Use function for elevation
	domain.set_quantity('friction', 0.01)        # Constant friction 
	domain.set_quantity('stage',                 # Dry bed
						expression='elevation')  

	#------------------------------------------------------------------------------
	# Setup boundary conditions
	#------------------------------------------------------------------------------
	Bi = Dirichlet_boundary([0.4, 0, 0])         # Inflow
	Br = Reflective_boundary(domain)             # Solid reflective wall

	domain.set_boundary({'left': Bi, 'right': Br, 'top': Br, 'bottom': Br})

	#------------------------------------------------------------------------------
	# Evolve system through time
	#------------------------------------------------------------------------------
	for t in domain.evolve(yieldstep=10.0, finaltime=100.0):
		print domain.timestepping_statistics()
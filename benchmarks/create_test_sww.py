#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
import anuga

def create_test_sww(sww_name):

	#------------------------------------------------------------------------------
	# Setup computational domain
	#------------------------------------------------------------------------------
	points, vertices, boundary = anuga.rectangular_cross(5, 5,
							     len1=50.0, len2=50.0) # Mesh

	domain = anuga.Domain(points, vertices, boundary)  # Create domain
	domain.set_name(sww_name)                  # Output name

	#------------------------------------------------------------------------------
	# Setup initial conditions
	#------------------------------------------------------------------------------
	def topography(x, y):
		return -x/2                             # linear bed slope

	domain.set_quantity('elevation', topography) # Use function for elevation
	domain.set_quantity('friction', 0.01)        # Constant friction 
	domain.set_quantity('stage',                 # Dry bed
			    expression='elevation')  

	#------------------------------------------------------------------------------
	# Setup boundary conditions
	#------------------------------------------------------------------------------
	Bi = anuga.Dirichlet_boundary([0.4, 0, 0])         # Inflow
	Br = anuga.Reflective_boundary(domain)             # Solid reflective wall

	domain.set_boundary({'left': Bi, 'right': Br, 'top': Br, 'bottom': Br})

	#------------------------------------------------------------------------------
	# Evolve system through time
	#------------------------------------------------------------------------------
	for t in domain.evolve(yieldstep=10.0, finaltime=100.0):
		print domain.timestepping_statistics()

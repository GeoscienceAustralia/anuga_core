"""Simple buildings example using ANUGA

    A simple example using internal boundaries to represent buildings, using new
    functionality found in ANUGA 1.2.
    
    written by James Hudson
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
import anuga

def poly_from_box(W, E, N, S):
    """
        Helper function to create a counter-clockwise polygon
        from the given box.
        
        N, E, S, W are the extents of the box.
    """
    return [[W, N], [W, S], [E, S], [E, N]]

#------------------------------------------------------------------------------
# Setup computational domain - a long skinny channel with obstructing boxes
#------------------------------------------------------------------------------
length = 50
width = 10
resolution = 1.0 # make this number smaller to make the simulation more accurate

# Create the "world" as a long, skinny channel
boundary_poly = poly_from_box(0, length, 0, width)
boundary_tags = {'left': [0], 'bottom': [1], 'right': [2], 'top': [3]}

# Place 3 buildings downstream
building_polys = [  poly_from_box(10, 15, 2.5, 7.5),        # upstream box
                    poly_from_box(22.5, 27.5, 1.5, 6.5),    # middle box
                    poly_from_box(35, 40, 3.5, 8.5)]        # downstream box
					
building_tags = [ { 'upstream'   : [1]},    # Set segment[0], default others
                  None,                     # Use default interior tag
		  { 'downstream' : [0,1]} ] # set segments[0,1],default others


# create a domain mesh, with 3 building holes in it
domain = anuga.create_domain_from_regions(boundary_poly,
                                            boundary_tags=boundary_tags,
                                            maximum_triangle_area = resolution,
                                            mesh_filename = 'building.msh',
                                            interior_holes = building_polys,
                                            hole_tags = building_tags,
                                            use_cache=True, # to speed it up
                                            verbose=True)   # log output on
 
domain.set_name('buildings')                 # Output name


#------------------------------------------------------------------------------
# Setup initial conditions
#------------------------------------------------------------------------------
def topography(x,y):
    return -x/15                             # gentle linear bed slope

domain.set_quantity('elevation', topography) # Use function for elevation
domain.set_quantity('friction', 0.03)        # Constant friction 
domain.set_quantity('stage',
                    expression='elevation')  # Dry initial condition

#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
Bi = anuga.Dirichlet_boundary([1.0, 0, 0])   # Inflow
Br = anuga.Reflective_boundary(domain)       # Solid reflective wall
Bo = anuga.Dirichlet_boundary([-5, 0, 0])    # Outflow

print (domain.get_boundary_tags())

domain.set_boundary({'left': Bi, 'right': Bo, 'top': Br, 'bottom': Br,
                     'interior': Br      ,# default interior boundary tag
		     'downstream': Br    ,# downstream building
		     'upstream' : Br      # upstream building boundary
                    })


#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------
for t in domain.evolve(yieldstep=0.2, finaltime=15.0):
    print (domain.timestepping_statistics())

    
# now turn off the tap
domain.set_boundary({'left': Bo})

for t in domain.evolve(yieldstep=0.1, finaltime=30.0):
    print (domain.timestepping_statistics() )       


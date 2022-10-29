"""Example of shallow water wave equation.

This example models a kitchen sink with inflow from the tap, hydraulic jump and also the drain opening.

Specific methods pertaining to the 2D shallow water equation
are imported from shallow_water
for use with the generic finite volume framework

Conserved quantities are h, uh and vh stored as elements 0, 1 and 2 in the
numerical vector named conserved_quantities.
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------

import anuga

#------------------------------------------------------------------------------
# Setup computational domain
#------------------------------------------------------------------------------

#N = 120
N = 40
points, vertices, boundary = anuga.rectangular(N, N)
domain = anuga.Domain(points, vertices, boundary)   
domain.set_name('kitchen_sink')                 # Output name
print('Size', len(domain))


#------------------------------------------------------------------------------
# Setup initial conditions
#------------------------------------------------------------------------------

domain.set_quantity('elevation', 0.0)  # Flat bed elevation
domain.set_quantity('stage', 0.0)      # Constant stage, dry initially
domain.set_quantity('friction', 0.005) # Constant friction 


#------------------------------------------------------------------------------
# Setup sink and tap via Inlet_operators
#------------------------------------------------------------------------------

sink = anuga.Inlet_operator(domain, region=anuga.Region(domain, center=[0.5, 0.5], radius=0.1), Q=0.0)
tap =  anuga.Inlet_operator(domain, region=anuga.Region(domain, center=(0.2, 0.2), radius=0.02), Q=lambda t: min(0.05*t, 0.05))
 
#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------

Br = anuga.Reflective_boundary(domain)              # Solid reflective wall
domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------
for t in domain.evolve(yieldstep = 0.02, finaltime = 15):
    domain.write_time()
    
    if domain.get_time() >= 10 and tap.get_Q() != 0.0:
        print('    Turning tap off')
        tap.set_Q(0.0)
        
    if domain.get_time() >= 1 and sink.get_Q() == 0.0:
        print('    Turning drain on')
        sink.set_Q(-0.75)       
    

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
# Setup specialised forcing terms
#------------------------------------------------------------------------------

class Inflow:
    """Class Inflow - general 'rain and drain' forcing term.
    
    Useful for implementing flows in and out of the domain.
    
    Inflow(center, radius, flow)
    
    center [m]: Coordinates at center of flow point
    radius [m]: Size of circular area
    flow [m/s]: Rate of change of quantity over the specified area.  
                This parameter can be either a constant or a function of time. 
                Positive values indicate inflow, 
                negative values indicate outflow.
    
    Examples
    
    Inflow((0.7, 0.4), 0.07, -0.2) # Constant drain at 0.2 m/s.
                                   # This corresponds to a flow of 
                                   # 0.07**2*pi*0.2 = 0.00314 m^3/s
                                   
    Inflow((0.5, 0.5), 0.001, lambda t: min(4*t, 5)) # Tap turning up to 
                                                     # a maximum inflow of
                                                     # 5 m/s over the 
                                                     # specified area 
    """
    
    # FIXME (OLE): Add a polygon as an alternative.
    # FIXME (OLE): Generalise to all quantities

    def __init__(self, 
                 center=None, radius=None,
                 flow=0.0,
                 quantity_name = 'stage'):
                 
        if center is not None and radius is not None:
            assert len(center) == 2
        else:
            msg = 'Both center and radius must be specified'
            raise Exception(msg)
    
        self.center = center
        self.radius = radius
        self.flow = flow
        self.quantity = domain.quantities[quantity_name].explicit_update
        
    
    def __call__(self, domain):
    
        # Determine indices in flow area
        if not hasattr(self, 'indices'):
            center = self.center
            radius = self.radius
            
            N = len(domain)    
            self.indices = []
            coordinates = domain.get_centroid_coordinates()     
            for k in range(N):
                x, y = coordinates[k,:] # Centroid
                if ((x - center[0])**2 + (y - center[1])**2) < radius**2:
                    self.indices.append(k)    

        # Update inflow
        if callable(self.flow):
            flow = self.flow(domain.get_time())
        else:
            flow = self.flow
                   
        for k in self.indices:
            self.quantity[k] += flow                    
                    

sink = Inflow(center=[0.7, 0.4], radius=0.0707, flow = 0.0)             
tap = Inflow(center=(0.5, 0.5), radius=0.0316,
             flow=lambda t: min(4*t, 5)) # Tap turning up over the first 1.2s       
        

domain.forcing_terms.append(tap)
domain.forcing_terms.append(sink)

#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------

Br = anuga.Reflective_boundary(domain)              # Solid reflective wall
domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------
for t in domain.evolve(yieldstep = 0.01, finaltime = 15):
    domain.write_time()
    
    if domain.get_time() >= 4 and tap.flow != 0.0:
        print('Turning tap off')
        tap.flow = 0.0
        
    if domain.get_time() >= 3 and sink.flow == 0.0:
        print('Turning drain on')
        sink.flow = -3.5        
    

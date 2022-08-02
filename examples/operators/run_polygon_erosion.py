"""Simple water flow example using ANUGA

Water flowing down a channel with a topography that varies with time
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
from anuga import rectangular_cross
from anuga import Domain
from anuga import Reflective_boundary
from anuga import Dirichlet_boundary
from anuga import Time_boundary



#===============================================================================
# Stick all the class and function definitions here
#===============================================================================

#-------------------------------------------------------------------------------
# Here are the def's for defining quantities
#-------------------------------------------------------------------------------
def topography(x,y):
    """Complex topography defined by a function of vectors x and y."""

    z = -x/100

    # Step
    id = (2 < x) & (x < 4)
    z[id] += 0.4 - 0.05*y[id]

    # Permanent pole
    id = (x - 8)**2 + (y - 2)**2 < 0.4**2
    z[id] += 1

    # Dam
    id = (12 < x) & (x < 13)
    z[id] += 0.4

    return z


#-------------------------------------------------------------------------------
# Inherit from erosion operator
#-------------------------------------------------------------------------------
from anuga.operators.erosion_operators import Polygonal_erosion_operator

class My_polygon_erosion_operator(Polygonal_erosion_operator):
    """
    Local version of erosion confined to a polygon

    """

    def __init__(self, domain,
                 threshold=0.0,
                 base=0.0,
                 polygon=None,
                 verbose=False):


        Polygonal_erosion_operator.__init__(self, domain, threshold, base, polygon, verbose)



    def update_quantities(self):
        """Update the vertex values of the quantities to model erosion
        """
        import numpy as num
        
        t = self.get_time()
        dt = self.get_timestep()

        updated = True

        if self.indices is None:

            #--------------------------------------
            # Update all three vertices for each cell
            #--------------------------------------
            self.elev_v[:] = self.elev_v + 0.0

        else:
            ind = self.indices
            m = num.sqrt(self.xmom_c[ind]**2 + self.ymom_c[ind]**2)
            
            if self.domain.get_using_discontinuous_elevation():
                height = self.stage_c[ind] - self.elev_c[ind]
                
                m = num.where(m>self.threshold, m, 0.0)
                self.elev_c[ind] = num.maximum(self.elev_c[ind] - m*dt, self.base)

                self.stage_c[ind] = self.elev_c[ind] + height
            else:
            #--------------------------------------
            # Update all three vertices for each cell
            #--------------------------------------
            
            
                m = num.vstack((m,m,m)).T
                m = num.where(m>self.threshold, m, 0.0)
                self.elev_v[ind] = num.maximum(self.elev_v[ind] - m*dt, self.base)
                #num.maximum(self.elev_v[ind] - momentum*dt, Z)


        return updated



#===============================================================================
# Now to the standard anuga model setup
#===============================================================================


#-------------------------------------------------------------------------------
# Setup computational domain
#-------------------------------------------------------------------------------
length = 24.
width = 5.
dx = dy = 0.2 #.1           # Resolution: Length of subdivisions on both axes

points, vertices, boundary = rectangular_cross(int(length/dx), int(width/dy),
                                               len1=length, len2=width)
domain = Domain(points, vertices, boundary)
domain.set_name() # Output name based on script
domain.set_flow_algorithm('DE1')
print (domain.statistics())
domain.set_store_vertices_uniquely(True)

domain.set_quantities_to_be_stored({'elevation': 2,
                                    'stage': 2,
                                    'xmomentum': 2,
                                    'ymomentum': 2})

#-------------------------------------------------------------------------------
# Setup initial conditions
#-------------------------------------------------------------------------------


domain.set_quantity('elevation', topography)           # elevation is a function
domain.set_quantity('friction', 0.01)                  # Constant friction
domain.set_quantity('stage', expression='elevation')   # Dry initial condition

#-------------------------------------------------------------------------------
# Setup boundary conditions
#-------------------------------------------------------------------------------
Bi = Dirichlet_boundary([0.5, 0, 0])          # Inflow
Br = Reflective_boundary(domain)              # Solid reflective wall
Bo = Dirichlet_boundary([-5, 0, 0])           # Outflow

domain.set_boundary({'left': Bi, 'right': Bo, 'top': Br, 'bottom': Br})

#-------------------------------------------------------------------------------
# Setup erosion operator in the middle of dam
#-------------------------------------------------------------------------------
polygon1 = [ [12., 0.0], [13., 0.0], [13., 5.0], [12., 5.0] ]
op1 = My_polygon_erosion_operator(domain, threshold=0.0, base=-0.1, polygon=polygon1)




#-------------------------------------------------------------------------------
# Evolve system through time
#-------------------------------------------------------------------------------
for t in domain.evolve(yieldstep=0.2, finaltime=60.0):
    domain.print_timestepping_statistics()
    domain.print_operator_timestepping_statistics()










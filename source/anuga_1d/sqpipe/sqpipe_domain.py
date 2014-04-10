"""Class Domain -
1D interval domains for finite-volume computations of
the shallow water wave equation.

This module contains a specialisation of class Domain from module domain.py
consisting of methods specific to the Shallow Water Wave Equation


U_t + E_x = S

where

U = [w, uh]
E = [uh, u^2h + gh^2/2]
S represents source terms forcing the system
(e.g. gravity, friction, wind stress, ...)

and _t, _x, _y denote the derivative with respect to t, x and y respectiely.

The quantities are

symbol    variable name    explanation
x         x                horizontal distance from origin [m]
z         elevation        elevation of bed on which flow is modelled [m]
h         height           water height above z [m]
w         stage            absolute water level, w = z+h [m]
u                          speed in the x direction [m/s]
uh        xmomentum        momentum in the x direction [m^2/s]

eta                        mannings friction coefficient [to appear]
nu                         wind stress coefficient [to appear]

The conserved quantities are w, uh

For details see e.g.
Christopher Zoppou and Stephen Roberts,
Catastrophic Collapse of Water Supply Reservoirs in Urban Areas,
Journal of Hydraulic Engineering, vol. 127, No. 7 July 1999


John Jakeman, Ole Nielsen, Stephen Roberts, Duncan Gray, Christopher Zoppou
Geoscience Australia, 2006
"""

import numpy

from anuga_1d.base.generic_domain import Generic_domain
from anuga_1d.sqpipe.sqpipe_forcing_terms import *
from anuga_1d.sqpipe.sqpipe_boundary_conditions import *

#Shallow water domain
class Domain(Generic_domain):

    def __init__(self, coordinates, boundary = None, forcing_terms = [], state = None, update_state_flag = True, bulk_modulus = 1.0, tagged_elements = None):
        conserved_quantities = ['area', 'discharge']
        evolved_quantities = ['area', 'discharge', 'elevation', 'height', 'velocity','width','top','stage']
        other_quantities = ['friction']
        Generic_domain.__init__(self,
                                coordinates          = coordinates,
                                boundary             = boundary,
                                conserved_quantities = conserved_quantities,
                                evolved_quantities   = evolved_quantities,
                                other_quantities     = other_quantities,
                                tagged_elements      = tagged_elements)

        self.bulk_modulus = bulk_modulus

        # Allocate space for tracking state (if not provided)
        if (state is None):
            self.state = numpy.zeros(self.number_of_elements, numpy.int)
        else:
            self.state = state            

        self.update_state_flag = update_state_flag

        # Get computational parameters
        from anuga_1d.config import minimum_allowed_height, g, h0, epsilon
        self.minimum_allowed_height = minimum_allowed_height
        self.g = g
        self.h0 = h0
        self.epsilon = epsilon

        #setup manning's friction
        f = manning_friction
        self.forcing_terms.append(f)
        
        #forcing terms
        for f in forcing_terms:
            self.forcing_terms.append(f)

        #Stored output
        self.store = True
        self.format = 'sww'
        self.smooth = True

        #Evolve parametrs
        self.cfl = 1.0
        
        #Reduction operation for get_vertex_values
        from anuga_1d.base.util import mean
        self.reduction = mean
        #self.reduction = min  #Looks better near steep slopes

        self.quantities_to_be_stored = conserved_quantities

	self.__doc__ = 'sqpipe_domain'

        self.check_integrity()


    def set_quantities_to_be_stored(self, q):
        """Specify which quantities will be stored in the sww file.

        q must be either:
          - the name of a quantity
          - a list of quantity names
          - None

        In the two first cases, the named quantities will be stored at each
        yieldstep
        (This is in addition to the quantities elevation and friction)  
        If q is None, storage will be switched off altogether.
        """


        if q is None:
            self.quantities_to_be_stored = []            
            self.store = False
            return

        if isinstance(q, basestring):
            q = [q] # Turn argument into a list

        #Check correctness    
        for quantity_name in q:
            msg = 'Quantity %s is not a valid conserved quantity' %quantity_name
            assert quantity_name in self.conserved_quantities, msg 
        
        self.quantities_to_be_stored = q
        

    def check_integrity(self):
        if (len(self.state) != self.number_of_elements):
                raise Exception("state has invalid length")
        Generic_domain.check_integrity(self)

    def compute_fluxes(self):        
        # Set initial timestep to a large value
        import sys
        timestep = float(sys.maxint)

        # Update state before computing flux (if necessary)
        if self.update_state_flag:
            self.update_state()

        # Import flux method and call it
        from anuga_1d.sqpipe.sqpipe_comp_flux import compute_fluxes_ext

        self.flux_timestep = compute_fluxes_ext(timestep,self)

    def update_state(self):
        h = self.quantities['height'].centroid_values
        t = self.quantities['top'].centroid_values

        new_state = numpy.zeros_like(self.state)

        # Update boundary state
        new_state[0] = update_cell_state(0, [0, 1], h, t, self.state)
        N = self.number_of_elements-1
        new_state[N] = update_cell_state(N, [N-1, N], h, t, self.state)

        # Update interior states
        for i in range(1, N-1):
            new_state[i] = update_cell_state(i, [i-1, i+1], h, t, self.state)

        self.state = new_state
        #self.state = numpy.where(h>=t, 1, 0)

    def distribute_to_vertices_and_edges(self):
        # Shortcuts
        h0 = self.h0
        epsilon = self.epsilon

        area      = self.quantities['area']
        discharge = self.quantities['discharge']
        bed       = self.quantities['elevation']
        height    = self.quantities['height']
        velocity  = self.quantities['velocity']
        width     = self.quantities['width']
        top       = self.quantities['top']
        stage     = self.quantities['stage']

        #Arrays   
        a_C   = area.centroid_values
        d_C   = discharge.centroid_values
        z_C   = bed.centroid_values
        h_C   = height.centroid_values
        u_C   = velocity.centroid_values
        b_C   = width.centroid_values
        t_C   = top.centroid_values
        w_C   = stage.centroid_values

        #Calculate height (and fix negatives)better be non-negative!
        h_C[:] = a_C/(b_C + h0/(b_C + h0)) # Do we need to protect against small b? Make a b0 instead of h0 or use epsilon?
        w_C[:] = z_C + h_C
        u_C[:]  = d_C/(h_C + h0/(h_C + h0))

        for name in ['velocity', 'stage' ]:
            Q = self.quantities[name]
            if self.order == 1:
                Q.extrapolate_first_order()
            elif self.order == 2:
                Q.extrapolate_second_order()
            else:
                raise 'Unknown order'

        # These have been extrapolated
        w_V = self.quantities['stage'].vertex_values
        u_V = self.quantities['velocity'].vertex_values

        # These are the given geometry and remain fixed
        z_V  = bed.vertex_values
        b_V  = width.vertex_values
        t_V  = top.vertex_values

        # These need to be updated
        a_V  = area.vertex_values
        d_V  = discharge.vertex_values
        h_V  = height.vertex_values

        h_V[:,:]    = w_V - z_V

        # Fix up small heights
        h_0 = numpy.where(h_V[:,0] < 0.0, 0.0, h_V[:,0])
        h_1 = numpy.where(h_V[:,0] < 0.0, h_V[:,1]+h_V[:,0], h_V[:,1])

        h_V[:,0] = h_0
        h_V[:,1] = h_1

        h_0 = numpy.where(h_V[:,1] < 0.0, h_V[:,1]+h_V[:,0], h_V[:,0])
        h_1 = numpy.where(h_V[:,1] < 0.0, 0.0, h_V[:,1])

        h_V[:,0] = h_0
        h_V[:,1] = h_1

        # It may still be possible that h is small 
        # If we set h to zero, we should also set u to 0
        h_V[:,:] = numpy.where(h_V < epsilon, 0.0, h_V)
        u_V[:,:] = numpy.where(h_V < epsilon, 0.0, u_V)

        # Reconstruct conserved quantities and make everything consistent
        # with new h values
        w_V[:,:] = z_V + h_V
        a_V[:,:] = h_V * b_V
        d_V[:,:] = u_V * a_V


        return
        
    def evolve(self, yieldstep = None, finaltime = None, duration = None,
               skip_initial_step = False):
        """Specialisation of basic evolve method from parent class
        """

        #Call basic machinery from parent class
        for t in Generic_domain.evolve(self, yieldstep, finaltime,duration,
                                       skip_initial_step):

            #Pass control on to outer loop for more specific actions
            yield(t)


    def get_mass(self):
        """Returns array of mass in cells
        """

        # The equivalent area is the conserved mass quantity
        # It is equal to \rho/\rho_0 A where A is the cross sectional
        # area of the fluid where $\rho_0 is some fixed reference density
        area = self.quantities['area']
        a = area.centroid_values
        
        return a * self.areas
        

# Auxillary methods

# Method to update the state of a cell
# cell is the index of the cell to update
# neighbours is a list of the indices of it's neighbours
# (for boundary cells, this could include the cell itself) 
def update_cell_state(cell, neighbours, h, t, s):
    # If height is bigger than top, then pressurise
    if h[cell] >= t[cell]:
        return 1
    else:
        # If the cell was pressurised, and all it's neighbours
        # are pressurised, it remains pressurised even though
        # height is less than top
        # Otherwise, it's free surface
        if all([s[i] for i in neighbours]):
            return s[cell]
        else:
            return 0

#=============== End of Shallow Water Domain ===============================

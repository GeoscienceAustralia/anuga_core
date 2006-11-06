"""Class Domain -
2D triangular domains for finite-volume computations of
the shallow water wave equation.


$Description:
This module contains a specialisation of class Domain from module domain.py
consisting of methods specific to the Shallow Water Wave Equation


U_t + E_x + G_y = S

where

U = [w, uh, vh]
E = [uh, u^2h + gh^2/2, uvh]
G = [vh, uvh, v^2h + gh^2/2]
S represents source terms forcing the system
(e.g. gravity, friction, wind stress, ...)

and _t, _x, _y denote the derivative with respect to t, x and y respectively.

The quantities are

symbol    variable name    explanation
x         x                horizontal distance from origin [m]
y         y                vertical distance from origin [m]
z         elevation        elevation of bed on which flow is modelled [m]
h         height           water height above z [m]
w         stage            absolute water level, w = z+h [m]
u                          speed in the x direction [m/s]
v                          speed in the y direction [m/s]
uh        xmomentum        momentum in the x direction [m^2/s]
vh        ymomentum        momentum in the y direction [m^2/s]

eta                        mannings friction coefficient [to appear]
nu                         wind stress coefficient [to appear]

The conserved quantities are w, uh, vh

$References
Catastrophic Collapse of Water Supply Reservoirs in Urban Areas,
Christopher Zoppou and Stephen Roberts,
Journal of Hydraulic Engineering, vol. 127, No. 7 July 1999

Hydrodynamic modelling of coastal inundation. 
Nielsen, O., S. Roberts, D. Gray, A. McPherson and A. Hitchman
In Zerger, A. and Argent, R.M. (eds) MODSIM 2005 International Congress on
Modelling and Simulation. Modelling and Simulation Society of Australia and
New Zealand, December 2005, pp. 518-523. ISBN: 0-9758400-2-9.
http://www.mssanz.org.au/modsim05/papers/nielsen.pdf


Ole Nielsen, Stephen Roberts, Duncan Gray, Christopher Zoppou
Geoscience Australia, 2004


$Author$
$Revision$
$Date$
$LastChangedDate$
$LastChangedRevision$
$LastChangedBy$
$HeadURL$
"""

#Subversion keywords:
#
#$LastChangedDate$
#$LastChangedRevision$
#$LastChangedBy$

from Numeric import zeros, ones, Float, array, sum, size
from Numeric import compress, arange


from anuga.abstract_2d_finite_volumes.domain import Domain as Generic_Domain
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions\
     import Boundary
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions\
     import File_boundary
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions\
     import Dirichlet_boundary
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions\
     import Time_boundary
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions\
     import Transmissive_boundary

from anuga.utilities.numerical_tools import gradient, mean
from anuga.config import minimum_storable_height
from anuga.config import minimum_allowed_height, maximum_allowed_speed
from anuga.config import g, beta_h, beta_w, beta_w_dry,\
     beta_uh, beta_uh_dry, beta_vh, beta_vh_dry
from anuga.config import alpha_balance


#Shallow water domain
class Domain(Generic_Domain):

    def __init__(self,
                 coordinates=None,
                 vertices=None,
                 boundary=None,
                 tagged_elements=None,
                 geo_reference=None,
                 use_inscribed_circle=False,
                 mesh_filename=None,
                 use_cache=False,
                 verbose=False,
                 full_send_dict=None,
                 ghost_recv_dict=None,
                 processor=0,
                 numproc=1,
                 number_of_full_nodes=0,
                 number_of_full_triangles=0):


        conserved_quantities = ['stage', 'xmomentum', 'ymomentum']
        other_quantities = ['elevation', 'friction']
        Generic_Domain.__init__(self,
                                coordinates,
                                vertices,
                                boundary,
                                conserved_quantities,
                                other_quantities,
                                tagged_elements,
                                geo_reference,
                                use_inscribed_circle,
                                mesh_filename,
                                use_cache,
                                verbose,
                                full_send_dict,
                                ghost_recv_dict,
                                processor,
                                numproc,
                                number_of_full_nodes=number_of_full_nodes,
                                number_of_full_triangles=number_of_full_triangles) 

        self.minimum_allowed_height = minimum_allowed_height
        self.maximum_allowed_speed = maximum_allowed_speed
        self.g = g
        self.beta_w      = beta_w
        self.beta_w_dry  = beta_w_dry
        self.beta_uh     = beta_uh
        self.beta_uh_dry = beta_uh_dry
        self.beta_vh     = beta_vh
        self.beta_vh_dry = beta_vh_dry
        self.beta_h      = beta_h
        self.alpha_balance = alpha_balance

        self.flux_function = flux_function_central
        #self.flux_function = flux_function_kinetic
        
        self.forcing_terms.append(manning_friction)
        self.forcing_terms.append(gravity)

        #Realtime visualisation
        self.visualiser = None
        self.visualise  = False
        self.visualise_color_stage = False
        self.visualise_stage_range = 1.0
        self.visualise_timer = True
        self.visualise_range_z = None

        #Stored output
        self.store = True
        self.format = 'sww'
        self.set_store_vertices_uniquely(False)
        self.minimum_storable_height = minimum_storable_height
        self.quantities_to_be_stored = ['stage','xmomentum','ymomentum']
                

    def set_all_limiters(self, beta):
        """Shorthand to assign one constant value [0,1[ to all limiters.
        0 Corresponds to first order, where as larger values make use of
        the second order scheme. 
        """

        self.beta_w      = beta
        self.beta_w_dry  = beta
        self.beta_uh     = beta
        self.beta_uh_dry = beta
        self.beta_vh     = beta
        self.beta_vh_dry = beta
        self.beta_h      = beta
        

    def set_store_vertices_uniquely(self, flag, reduction=None):
        """Decide whether vertex values should be stored uniquely as
        computed in the model or whether they should be reduced to one
        value per vertex using self.reduction.
        """
        self.smooth = not flag

        #Reduction operation for get_vertex_values
        if reduction is None:
            self.reduction = mean
            #self.reduction = min  #Looks better near steep slopes


    def set_minimum_storable_height(self, minimum_storable_height):
        """
        Set the minimum depth that will be recognised when writing
        to an sww file. This is useful for removing thin water layers
        that seems to be caused by friction creep.

        The minimum allowed sww depth is in meters.
        """
        self.minimum_storable_height = minimum_storable_height
        

    def set_maximum_allowed_speed(self, maximum_allowed_speed):
        """
        Set the maximum particle speed that is allowed in water
        shallower than minimum_allowed_height. This is useful for
        controlling speeds in very thin layers of water and at the same time
        allow some movement avoiding pooling of water.

        """
        self.maximum_allowed_speed = maximum_allowed_speed


    def set_quantities_to_be_stored(self, q):
        """Specify which quantities will be stored in the sww file.

        q must be either:
          - the name of a quantity
          - a list of quantity names
          - None

        In the two first cases, the named quantities will be stored at
        each yieldstep (This is in addition to the quantities elevation
        and friction)
        
        If q is None, storage will be switched off altogether.
        """


        if q is None:
            self.quantities_to_be_stored = []
            self.store = False
            return

        if isinstance(q, basestring):
            q = [q] # Turn argument into a list

        #Check correcness
        for quantity_name in q:
            msg = 'Quantity %s is not a valid conserved quantity'\
                  %quantity_name
            
            assert quantity_name in self.conserved_quantities, msg

        self.quantities_to_be_stored = q



    def get_wet_elements(self, indices=None):
        """Return indices for elements where h > minimum_allowed_height

        Optional argument:
            indices is the set of element ids that the operation applies to.

        Usage:
            indices = get_wet_elements()

        Note, centroid values are used for this operation            
        """

        # Water depth below which it is considered to be 0 in the model
        # FIXME (Ole): Allow this to be specified as a keyword argument as well
        from anuga.config import minimum_allowed_height

        
        elevation = self.get_quantity('elevation').\
                    get_values(location='centroids', indices=indices)
        stage = self.get_quantity('stage').\
                get_values(location='centroids', indices=indices)       
        depth = stage - elevation

        # Select indices for which depth > 0
        wet_indices = compress(depth > minimum_allowed_height,
                               arange(len(depth)))
        return wet_indices 


    def get_maximum_inundation_elevation(self, indices=None):
        """Return highest elevation where h > 0

        Optional argument:
            indices is the set of element ids that the operation applies to.

        Usage:
            q = get_maximum_inundation_elevation()

        Note, centroid values are used for this operation            
        """

        wet_elements = self.get_wet_elements(indices)
        return self.get_quantity('elevation').\
               get_maximum_value(indices=wet_elements)


    def get_maximum_inundation_location(self, indices=None):
        """Return highest elevation where h > 0

        Optional argument:
            indices is the set of element ids that the operation applies to.

        Usage:
            q = get_maximum_inundation_elevation()

        Note, centroid values are used for this operation            
        """

        wet_elements = self.get_wet_elements(indices)
        return self.get_quantity('elevation').\
               get_maximum_location(indices=wet_elements)    




    def initialise_visualiser(self,scale_z=1.0,rect=None):
        #Realtime visualisation
        if self.visualiser is None:
            from realtime_visualisation_new import Visualiser
            self.visualiser = Visualiser(self,scale_z,rect)
            self.visualiser.setup['elevation']=True
            self.visualiser.updating['stage']=True
        self.visualise = True
        if self.visualise_color_stage == True:
            self.visualiser.coloring['stage'] = True
            self.visualiser.qcolor['stage'] = (0.0, 0.0, 0.8)
        print 'initialise visualiser'
        print self.visualiser.setup
        print self.visualiser.updating

    def check_integrity(self):
        Generic_Domain.check_integrity(self)

        #Check that we are solving the shallow water wave equation

        msg = 'First conserved quantity must be "stage"'
        assert self.conserved_quantities[0] == 'stage', msg
        msg = 'Second conserved quantity must be "xmomentum"'
        assert self.conserved_quantities[1] == 'xmomentum', msg
        msg = 'Third conserved quantity must be "ymomentum"'
        assert self.conserved_quantities[2] == 'ymomentum', msg

    def extrapolate_second_order_sw(self):
        #Call correct module function
        #(either from this module or C-extension)
        extrapolate_second_order_sw(self)

    def compute_fluxes(self):
        #Call correct module function
        #(either from this module or C-extension)
        compute_fluxes(self)

    def distribute_to_vertices_and_edges(self):
        #Call correct module function
        #(either from this module or C-extension)
        distribute_to_vertices_and_edges(self)


    #FIXME: Under construction
#     def set_defaults(self):
#         """Set default values for uninitialised quantities.
#         This is specific to the shallow water wave equation
#         Defaults for 'elevation', 'friction', 'xmomentum' and 'ymomentum'
#         are 0.0. Default for 'stage' is whatever the value of 'elevation'.
#         """

#         for name in self.other_quantities + self.conserved_quantities:
#             print name
#             print self.quantities.keys()
#             if not self.quantities.has_key(name):
#                 if name == 'stage':

#                     if self.quantities.has_key('elevation'):
#                         z = self.quantities['elevation'].vertex_values
#                         self.set_quantity(name, z)
#                     else:
#                         self.set_quantity(name, 0.0)
#                 else:
#                     self.set_quantity(name, 0.0)



#         #Lift negative heights up
#         #z = self.quantities['elevation'].vertex_values
#         #w = self.quantities['stage'].vertex_values

#         #h = w-z

#         #for k in range(h.shape[0]):
#         #    for i in range(3):
#         #        if h[k, i] < 0.0:
#         #            w[k, i] = z[k, i]


#         #self.quantities['stage'].interpolate()


    def evolve(self,
               yieldstep = None,
               finaltime = None,
               duration = None,
               skip_initial_step = False):
        """Specialisation of basic evolve method from parent class
        """

        #Call check integrity here rather than from user scripts
        #self.check_integrity()

        msg = 'Parameter beta_h must be in the interval [0, 1['
        assert 0 <= self.beta_h <= 1.0, msg
        msg = 'Parameter beta_w must be in the interval [0, 1['
        assert 0 <= self.beta_w <= 1.0, msg


        #Initial update of vertex and edge values before any storage
        #and or visualisation
        self.distribute_to_vertices_and_edges()

        #Initialise real time viz if requested
        if self.visualise is True and self.time == 0.0:
            if self.visualiser is None:
                self.initialise_visualiser()

            self.visualiser.update_timer()
            self.visualiser.setup_all()

        #Store model data, e.g. for visualisation
        if self.store is True and self.time == 0.0:
            self.initialise_storage()
            #print 'Storing results in ' + self.writer.filename
        else:
            pass
            #print 'Results will not be stored.'
            #print 'To store results set domain.store = True'
            #FIXME: Diagnostic output should be controlled by
            # a 'verbose' flag living in domain (or in a parent class)

        #Call basic machinery from parent class
        for t in Generic_Domain.evolve(self,
                                       yieldstep=yieldstep,
                                       finaltime=finaltime,
                                       duration=duration,
                                       skip_initial_step=skip_initial_step):
            #Real time viz
            if self.visualise is True:
                self.visualiser.update_all()
                self.visualiser.update_timer()


            #Store model data, e.g. for subsequent visualisation
            if self.store is True:
                self.store_timestep(self.quantities_to_be_stored)

            #FIXME: Could maybe be taken from specified list
            #of 'store every step' quantities

            #Pass control on to outer loop for more specific actions
            yield(t)

    def initialise_storage(self):
        """Create and initialise self.writer object for storing data.
        Also, save x,y and bed elevation
        """

        from anuga.shallow_water.data_manager import get_dataobject

        #Initialise writer
        self.writer = get_dataobject(self, mode = 'w')

        #Store vertices and connectivity
        self.writer.store_connectivity()


    def store_timestep(self, name):
        """Store named quantity and time.

        Precondition:
           self.write has been initialised
        """
        self.writer.store_timestep(name)


#=============== End of Shallow Water Domain ===============================



#Rotation of momentum vector
def rotate(q, normal, direction = 1):
    """Rotate the momentum component q (q[1], q[2])
    from x,y coordinates to coordinates based on normal vector.

    If direction is negative the rotation is inverted.

    Input vector is preserved

    This function is specific to the shallow water wave equation
    """

    assert len(q) == 3,\
           'Vector of conserved quantities must have length 3'\
           'for 2D shallow water equation'

    try:
        l = len(normal)
    except:
        raise 'Normal vector must be an Numeric array'

    assert l == 2, 'Normal vector must have 2 components'


    n1 = normal[0]
    n2 = normal[1]

    r = zeros(len(q), Float) #Rotated quantities
    r[0] = q[0]              #First quantity, height, is not rotated

    if direction == -1:
        n2 = -n2


    r[1] =  n1*q[1] + n2*q[2]
    r[2] = -n2*q[1] + n1*q[2]

    return r



####################################
# Flux computation
def flux_function_central(normal, ql, qr, zl, zr):
    """Compute fluxes between volumes for the shallow water wave equation
    cast in terms of w = h+z using the 'central scheme' as described in

    Kurganov, Noelle, Petrova. 'Semidiscrete Central-Upwind Schemes For
    Hyperbolic Conservation Laws and Hamilton-Jacobi Equations'.
    Siam J. Sci. Comput. Vol. 23, No. 3, pp. 707-740.

    The implemented formula is given in equation (3.15) on page 714

    Conserved quantities w, uh, vh are stored as elements 0, 1 and 2
    in the numerical vectors ql and qr.

    Bed elevations zl and zr.
    """

    from anuga.config import g, epsilon
    from math import sqrt

    #Align momentums with x-axis
    q_left  = rotate(ql, normal, direction = 1)
    q_right = rotate(qr, normal, direction = 1)

    z = (zl+zr)/2 #Take average of field values

    w_left  = q_left[0]   #w=h+z
    h_left  = w_left-z
    uh_left = q_left[1]

    if h_left < epsilon:
        u_left = 0.0  #Could have been negative
        h_left = 0.0
    else:
        u_left  = uh_left/h_left


    w_right  = q_right[0]  #w=h+z
    h_right  = w_right-z
    uh_right = q_right[1]


    if h_right < epsilon:
        u_right = 0.0  #Could have been negative
        h_right = 0.0
    else:
        u_right  = uh_right/h_right

    vh_left  = q_left[2]
    vh_right = q_right[2]

    soundspeed_left  = sqrt(g*h_left)
    soundspeed_right = sqrt(g*h_right)

    #Maximal wave speed
    s_max = max(u_left+soundspeed_left, u_right+soundspeed_right, 0)

    #Minimal wave speed
    s_min = min(u_left-soundspeed_left, u_right-soundspeed_right, 0)

    #Flux computation

    #FIXME(Ole): Why is it again that we don't
    #use uh_left and uh_right directly in the first entries?
    flux_left  = array([u_left*h_left,
                        u_left*uh_left + 0.5*g*h_left**2,
                        u_left*vh_left])
    flux_right = array([u_right*h_right,
                        u_right*uh_right + 0.5*g*h_right**2,
                        u_right*vh_right])

    denom = s_max-s_min
    if denom == 0.0:
        edgeflux = array([0.0, 0.0, 0.0])
        max_speed = 0.0
    else:
        edgeflux = (s_max*flux_left - s_min*flux_right)/denom
        edgeflux += s_max*s_min*(q_right-q_left)/denom

        edgeflux = rotate(edgeflux, normal, direction=-1)
        max_speed = max(abs(s_max), abs(s_min))

    return edgeflux, max_speed

def erfcc(x):

    from math import fabs, exp

    z=fabs(x)
    t=1.0/(1.0+0.5*z)
    result=t*exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+
         t*(.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+
         t*(1.48851587+t*(-.82215223+t*.17087277)))))))))
    if x < 0.0:
        result = 2.0-result

    return result

def flux_function_kinetic(normal, ql, qr, zl, zr):
    """Compute fluxes between volumes for the shallow water wave equation
    cast in terms of w = h+z using the 'central scheme' as described in

    Zhang et. al., Advances in Water Resources, 26(6), 2003, 635-647.


    Conserved quantities w, uh, vh are stored as elements 0, 1 and 2
    in the numerical vectors ql an qr.

    Bed elevations zl and zr.
    """

    from anuga.config import g, epsilon
    from math import sqrt
    from Numeric import array

    #Align momentums with x-axis
    q_left  = rotate(ql, normal, direction = 1)
    q_right = rotate(qr, normal, direction = 1)

    z = (zl+zr)/2 #Take average of field values

    w_left  = q_left[0]   #w=h+z
    h_left  = w_left-z
    uh_left = q_left[1]

    if h_left < epsilon:
        u_left = 0.0  #Could have been negative
        h_left = 0.0
    else:
        u_left  = uh_left/h_left


    w_right  = q_right[0]  #w=h+z
    h_right  = w_right-z
    uh_right = q_right[1]


    if h_right < epsilon:
        u_right = 0.0  #Could have been negative
        h_right = 0.0
    else:
        u_right  = uh_right/h_right

    vh_left  = q_left[2]
    vh_right = q_right[2]

    soundspeed_left  = sqrt(g*h_left)
    soundspeed_right = sqrt(g*h_right)

    #Maximal wave speed
    s_max = max(u_left+soundspeed_left, u_right+soundspeed_right, 0)

    #Minimal wave speed
    s_min = min(u_left-soundspeed_left, u_right-soundspeed_right, 0)

    #Flux computation

    F_left  = 0.0
    F_right = 0.0
    from math import sqrt, pi, exp
    if h_left > 0.0:
        F_left = u_left/sqrt(g*h_left)
    if h_right > 0.0:
        F_right = u_right/sqrt(g*h_right)

    edgeflux = array([0.0, 0.0, 0.0])

    edgeflux[0] = h_left*u_left/2.0*erfcc(-F_left) +  \
          h_left*sqrt(g*h_left)/2.0/sqrt(pi)*exp(-(F_left**2)) + \
          h_right*u_right/2.0*erfcc(F_right) -  \
          h_right*sqrt(g*h_right)/2.0/sqrt(pi)*exp(-(F_right**2))

    edgeflux[1] = (h_left*u_left**2 + g/2.0*h_left**2)/2.0*erfcc(-F_left) + \
          u_left*h_left*sqrt(g*h_left)/2.0/sqrt(pi)*exp(-(F_left**2)) + \
          (h_right*u_right**2 + g/2.0*h_right**2)/2.0*erfcc(F_right) -  \
          u_right*h_right*sqrt(g*h_right)/2.0/sqrt(pi)*exp(-(F_right**2))

    edgeflux[2] = vh_left*u_left/2.0*erfcc(-F_left) + \
          vh_left*sqrt(g*h_left)/2.0/sqrt(pi)*exp(-(F_left**2)) + \
          vh_right*u_right/2.0*erfcc(F_right) - \
          vh_right*sqrt(g*h_right)/2.0/sqrt(pi)*exp(-(F_right**2))


    edgeflux = rotate(edgeflux, normal, direction=-1)
    max_speed = max(abs(s_max), abs(s_min))

    return edgeflux, max_speed



def compute_fluxes(domain):
    """Compute all fluxes and the timestep suitable for all volumes
    in domain.

    Compute total flux for each conserved quantity using "flux_function"

    Fluxes across each edge are scaled by edgelengths and summed up
    Resulting flux is then scaled by area and stored in
    explicit_update for each of the three conserved quantities
    stage, xmomentum and ymomentum

    The maximal allowable speed computed by the flux_function for each volume
    is converted to a timestep that must not be exceeded. The minimum of
    those is computed as the next overall timestep.

    Post conditions:
      domain.explicit_update is reset to computed flux values
      domain.timestep is set to the largest step satisfying all volumes.
    """

    import sys

    N = domain.number_of_elements

    #Shortcuts
    Stage = domain.quantities['stage']
    Xmom = domain.quantities['xmomentum']
    Ymom = domain.quantities['ymomentum']
    Bed = domain.quantities['elevation']

    #Arrays
    stage = Stage.edge_values
    xmom =  Xmom.edge_values
    ymom =  Ymom.edge_values
    bed =   Bed.edge_values

    stage_bdry = Stage.boundary_values
    xmom_bdry =  Xmom.boundary_values
    ymom_bdry =  Ymom.boundary_values

    flux = zeros(3, Float) #Work array for summing up fluxes


    #Loop
    timestep = float(sys.maxint)
    for k in range(N):

        flux[:] = 0.  #Reset work array
        for i in range(3):
            #Quantities inside volume facing neighbour i
            ql = [stage[k, i], xmom[k, i], ymom[k, i]]
            zl = bed[k, i]

            #Quantities at neighbour on nearest face
            n = domain.neighbours[k,i]
            if n < 0:
                m = -n-1 #Convert negative flag to index
                qr = [stage_bdry[m], xmom_bdry[m], ymom_bdry[m]]
                zr = zl #Extend bed elevation to boundary
            else:
                m = domain.neighbour_edges[k,i]
                qr = [stage[n, m], xmom[n, m], ymom[n, m]]
                zr = bed[n, m]


            #Outward pointing normal vector
            normal = domain.normals[k, 2*i:2*i+2]

            #Flux computation using provided function
            edgeflux, max_speed = domain.flux_function(normal, ql, qr, zl, zr)
            flux -= edgeflux * domain.edgelengths[k,i]

            #Update optimal_timestep on full cells
            if  domain.tri_full_flag[k] == 1:
                try:
                    timestep = min(timestep, 0.5*domain.radii[k]/max_speed)
                except ZeroDivisionError:
                    pass

        #Normalise by area and store for when all conserved
        #quantities get updated
        flux /= domain.areas[k]
        Stage.explicit_update[k] = flux[0]
        Xmom.explicit_update[k] = flux[1]
        Ymom.explicit_update[k] = flux[2]


    domain.timestep = timestep

#MH090605 The following method belongs to the shallow_water domain class
#see comments in the corresponding method in shallow_water_ext.c
def extrapolate_second_order_sw_c(domain):
    """Wrapper calling C version of extrapolate_second_order_sw
    """
    import sys

    N = domain.number_of_elements

    #Shortcuts
    Stage = domain.quantities['stage']
    Xmom = domain.quantities['xmomentum']
    Ymom = domain.quantities['ymomentum']
    Elevation = domain.quantities['elevation']
    from shallow_water_ext import extrapolate_second_order_sw
    extrapolate_second_order_sw(domain,
                                domain.surrogate_neighbours,
                                domain.number_of_boundaries,
                                domain.centroid_coordinates,
                                Stage.centroid_values,
                                Xmom.centroid_values,
                                Ymom.centroid_values,
                                Elevation.centroid_values,
                                domain.vertex_coordinates,
                                Stage.vertex_values,
                                Xmom.vertex_values,
                                Ymom.vertex_values,
                                Elevation.vertex_values)

def compute_fluxes_c(domain):
    """Wrapper calling C version of compute fluxes
    """

    import sys

    N = domain.number_of_elements

    #Shortcuts
    Stage = domain.quantities['stage']
    Xmom = domain.quantities['xmomentum']
    Ymom = domain.quantities['ymomentum']
    Bed = domain.quantities['elevation']

    timestep = float(sys.maxint)
    from shallow_water_ext import\
         compute_fluxes_ext_central as compute_fluxes_ext
    
    domain.timestep = compute_fluxes_ext(timestep, domain.epsilon, domain.g,
                                     domain.neighbours,
                                     domain.neighbour_edges,
                                     domain.normals,
                                     domain.edgelengths,
                                     domain.radii,
                                     domain.areas,
                                     domain.tri_full_flag,
                                     Stage.edge_values,
                                     Xmom.edge_values,
                                     Ymom.edge_values,
                                     Bed.edge_values,
                                     Stage.boundary_values,
                                     Xmom.boundary_values,
                                     Ymom.boundary_values,
                                     Stage.explicit_update,
                                     Xmom.explicit_update,
                                     Ymom.explicit_update,
                                     domain.already_computed_flux)


####################################
# Module functions for gradient limiting (distribute_to_vertices_and_edges)

def distribute_to_vertices_and_edges(domain):
    """Distribution from centroids to vertices specific to the
    shallow water wave
    equation.

    It will ensure that h (w-z) is always non-negative even in the
    presence of steep bed-slopes by taking a weighted average between shallow
    and deep cases.

    In addition, all conserved quantities get distributed as per either a
    constant (order==1) or a piecewise linear function (order==2).

    FIXME: more explanation about removal of artificial variability etc

    Precondition:
      All quantities defined at centroids and bed elevation defined at
      vertices.

    Postcondition
      Conserved quantities defined at vertices

    """

    from anuga.config import optimised_gradient_limiter

    #Remove very thin layers of water
    protect_against_infinitesimal_and_negative_heights(domain)

    #Extrapolate all conserved quantities
    if optimised_gradient_limiter:
        #MH090605 if second order,
        #perform the extrapolation and limiting on
        #all of the conserved quantitie

        if (domain._order_ == 1):
            for name in domain.conserved_quantities:
                Q = domain.quantities[name]
                Q.extrapolate_first_order()
        elif domain._order_ == 2:
            domain.extrapolate_second_order_sw()
        else:
            raise 'Unknown order'
    else:
        #old code:
        for name in domain.conserved_quantities:
            Q = domain.quantities[name]
            if domain._order_ == 1:
                Q.extrapolate_first_order()
            elif domain._order_ == 2:
                Q.extrapolate_second_order()
                Q.limit()
            else:
                raise 'Unknown order'


    #Take bed elevation into account when water heights are small
    balance_deep_and_shallow(domain)

    #Compute edge values by interpolation
    for name in domain.conserved_quantities:
        Q = domain.quantities[name]
        Q.interpolate_from_vertices_to_edges()


def protect_against_infinitesimal_and_negative_heights(domain):
    """Protect against infinitesimal heights and associated high velocities
    """

    #Shortcuts
    wc = domain.quantities['stage'].centroid_values
    zc = domain.quantities['elevation'].centroid_values
    xmomc = domain.quantities['xmomentum'].centroid_values
    ymomc = domain.quantities['ymomentum'].centroid_values
    hc = wc - zc  #Water depths at centroids

    #Update
    #FIXME: Modify accroditg to c-version - or discard altogether.
    for k in range(domain.number_of_elements):

        if hc[k] < domain.minimum_allowed_height:
            #Control stage
            if hc[k] < domain.epsilon:
                wc[k] = zc[k] # Contain 'lost mass' error

            #Control momentum
            xmomc[k] = ymomc[k] = 0.0


def protect_against_infinitesimal_and_negative_heights_c(domain):
    """Protect against infinitesimal heights and associated high velocities
    """

    #Shortcuts
    wc = domain.quantities['stage'].centroid_values
    zc = domain.quantities['elevation'].centroid_values
    xmomc = domain.quantities['xmomentum'].centroid_values
    ymomc = domain.quantities['ymomentum'].centroid_values

    from shallow_water_ext import protect

    protect(domain.minimum_allowed_height, domain.maximum_allowed_speed,
            domain.epsilon, wc, zc, xmomc, ymomc)



def h_limiter(domain):
    """Limit slopes for each volume to eliminate artificial variance
    introduced by e.g. second order extrapolator

    limit on h = w-z

    This limiter depends on two quantities (w,z) so it resides within
    this module rather than within quantity.py
    """

    N = domain.number_of_elements
    beta_h = domain.beta_h

    #Shortcuts
    wc = domain.quantities['stage'].centroid_values
    zc = domain.quantities['elevation'].centroid_values
    hc = wc - zc

    wv = domain.quantities['stage'].vertex_values
    zv = domain.quantities['elevation'].vertex_values
    hv = wv-zv

    hvbar = zeros(hv.shape, Float) #h-limited values

    #Find min and max of this and neighbour's centroid values
    hmax = zeros(hc.shape, Float)
    hmin = zeros(hc.shape, Float)

    for k in range(N):
        hmax[k] = hmin[k] = hc[k]
        for i in range(3):
            n = domain.neighbours[k,i]
            if n >= 0:
                hn = hc[n] #Neighbour's centroid value

                hmin[k] = min(hmin[k], hn)
                hmax[k] = max(hmax[k], hn)


    #Diffences between centroids and maxima/minima
    dhmax = hmax - hc
    dhmin = hmin - hc

    #Deltas between vertex and centroid values
    dh = zeros(hv.shape, Float)
    for i in range(3):
        dh[:,i] = hv[:,i] - hc

    #Phi limiter
    for k in range(N):

        #Find the gradient limiter (phi) across vertices
        phi = 1.0
        for i in range(3):
            r = 1.0
            if (dh[k,i] > 0): r = dhmax[k]/dh[k,i]
            if (dh[k,i] < 0): r = dhmin[k]/dh[k,i]

            phi = min( min(r*beta_h, 1), phi )

        #Then update using phi limiter
        for i in range(3):
            hvbar[k,i] = hc[k] + phi*dh[k,i]

    return hvbar



def h_limiter_c(domain):
    """Limit slopes for each volume to eliminate artificial variance
    introduced by e.g. second order extrapolator

    limit on h = w-z

    This limiter depends on two quantities (w,z) so it resides within
    this module rather than within quantity.py

    Wrapper for c-extension
    """

    N = domain.number_of_elements
    beta_h = domain.beta_h

    #Shortcuts
    wc = domain.quantities['stage'].centroid_values
    zc = domain.quantities['elevation'].centroid_values
    hc = wc - zc

    wv = domain.quantities['stage'].vertex_values
    zv = domain.quantities['elevation'].vertex_values
    hv = wv - zv

    #Call C-extension
    from shallow_water_ext import h_limiter_sw as h_limiter
    hvbar = h_limiter(domain, hc, hv)

    return hvbar


def balance_deep_and_shallow(domain):
    """Compute linear combination between stage as computed by
    gradient-limiters limiting using w, and stage computed by
    gradient-limiters limiting using h (h-limiter).
    The former takes precedence when heights are large compared to the
    bed slope while the latter takes precedence when heights are
    relatively small.  Anything in between is computed as a balanced
    linear combination in order to avoid numerical disturbances which
    would otherwise appear as a result of hard switching between
    modes.

    The h-limiter is always applied irrespective of the order.
    """

    #Shortcuts
    wc = domain.quantities['stage'].centroid_values
    zc = domain.quantities['elevation'].centroid_values
    hc = wc - zc

    wv = domain.quantities['stage'].vertex_values
    zv = domain.quantities['elevation'].vertex_values
    hv = wv-zv

    #Limit h
    hvbar = h_limiter(domain)

    for k in range(domain.number_of_elements):
        #Compute maximal variation in bed elevation
        #  This quantitiy is
        #    dz = max_i abs(z_i - z_c)
        #  and it is independent of dimension
        #  In the 1d case zc = (z0+z1)/2
        #  In the 2d case zc = (z0+z1+z2)/3

        dz = max(abs(zv[k,0]-zc[k]),
                 abs(zv[k,1]-zc[k]),
                 abs(zv[k,2]-zc[k]))


        hmin = min( hv[k,:] )

        #Create alpha in [0,1], where alpha==0 means using the h-limited
        #stage and alpha==1 means using the w-limited stage as
        #computed by the gradient limiter (both 1st or 2nd order)

        #If hmin > dz/2 then alpha = 1 and the bed will have no effect
        #If hmin < 0 then alpha = 0 reverting to constant height above bed.

        if dz > 0.0:
            alpha = max( min( 2*hmin/dz, 1.0), 0.0 )
        else:
            #Flat bed
            alpha = 1.0

        #Let
        #
        #  wvi be the w-limited stage (wvi = zvi + hvi)
        #  wvi- be the h-limited state (wvi- = zvi + hvi-)
        #
        #
        #where i=0,1,2 denotes the vertex ids
        #
        #Weighted balance between w-limited and h-limited stage is
        #
        #  wvi := (1-alpha)*(zvi+hvi-) + alpha*(zvi+hvi)
        #
        #It follows that the updated wvi is
        #  wvi := zvi + (1-alpha)*hvi- + alpha*hvi
        #
        # Momentum is balanced between constant and limited


        #for i in range(3):
        #    wv[k,i] = zv[k,i] + hvbar[k,i]

        #return

        if alpha < 1:

            for i in range(3):
                wv[k,i] = zv[k,i] + (1-alpha)*hvbar[k,i] + alpha*hv[k,i]

            #Momentums at centroids
            xmomc = domain.quantities['xmomentum'].centroid_values
            ymomc = domain.quantities['ymomentum'].centroid_values

            #Momentums at vertices
            xmomv = domain.quantities['xmomentum'].vertex_values
            ymomv = domain.quantities['ymomentum'].vertex_values

            # Update momentum as a linear combination of
            # xmomc and ymomc (shallow) and momentum
            # from extrapolator xmomv and ymomv (deep).
            xmomv[k,:] = (1-alpha)*xmomc[k] + alpha*xmomv[k,:]
            ymomv[k,:] = (1-alpha)*ymomc[k] + alpha*ymomv[k,:]


def balance_deep_and_shallow_c(domain):
    """Wrapper for C implementation
    """

    #Shortcuts
    wc = domain.quantities['stage'].centroid_values
    zc = domain.quantities['elevation'].centroid_values
    hc = wc - zc

    wv = domain.quantities['stage'].vertex_values
    zv = domain.quantities['elevation'].vertex_values
    hv = wv - zv

    #Momentums at centroids
    xmomc = domain.quantities['xmomentum'].centroid_values
    ymomc = domain.quantities['ymomentum'].centroid_values

    #Momentums at vertices
    xmomv = domain.quantities['xmomentum'].vertex_values
    ymomv = domain.quantities['ymomentum'].vertex_values

    #Limit h
    hvbar = h_limiter(domain)

    #This is how one would make a first order h_limited value
    #as in the old balancer (pre 17 Feb 2005):
    #from Numeric import zeros, Float
    #hvbar = zeros( (len(hc), 3), Float)
    #for i in range(3):
    #    hvbar[:,i] = hc[:]

    from shallow_water_ext import balance_deep_and_shallow
    balance_deep_and_shallow(domain, wc, zc, hc, wv, zv, hv, hvbar,
                             xmomc, ymomc, xmomv, ymomv)




###############################################
#Boundaries - specific to the shallow water wave equation
class Reflective_boundary(Boundary):
    """Reflective boundary returns same conserved quantities as
    those present in its neighbour volume but reflected.

    This class is specific to the shallow water equation as it
    works with the momentum quantities assumed to be the second
    and third conserved quantities.
    """

    def __init__(self, domain = None):
        Boundary.__init__(self)

        if domain is None:
            msg = 'Domain must be specified for reflective boundary'
            raise msg

        #Handy shorthands
        self.stage   = domain.quantities['stage'].edge_values
        self.xmom    = domain.quantities['xmomentum'].edge_values
        self.ymom    = domain.quantities['ymomentum'].edge_values
        self.normals = domain.normals

        self.conserved_quantities = zeros(3, Float)

    def __repr__(self):
        return 'Reflective_boundary'


    def evaluate(self, vol_id, edge_id):
        """Reflective boundaries reverses the outward momentum
        of the volume they serve.
        """

        q = self.conserved_quantities
        q[0] = self.stage[vol_id, edge_id]
        q[1] = self.xmom[vol_id, edge_id]
        q[2] = self.ymom[vol_id, edge_id]

        normal = self.normals[vol_id, 2*edge_id:2*edge_id+2]


        r = rotate(q, normal, direction = 1)
        r[1] = -r[1]
        q = rotate(r, normal, direction = -1)

        return q



class Transmissive_Momentum_Set_Stage_boundary(Boundary):
    """Returns same momentum conserved quantities as
    those present in its neighbour volume.
    Sets stage by specifying a function f of time which may either be a
    vector function or a scalar function

    Example:

    def waveform(t): 
        return sea_level + normalized_amplitude/cosh(t-25)**2

    Bts = Transmissive_Momentum_Set_Stage_boundary(domain, waveform)
    

    Underlying domain must be specified when boundary is instantiated
    """

    def __init__(self, domain = None, function=None):
        Boundary.__init__(self)

        if domain is None:
            msg = 'Domain must be specified for this type boundary'
            raise msg

        if function is None:
            msg = 'Function must be specified for this type boundary'
            raise msg

        self.domain   = domain
        self.function = function

    def __repr__(self):
        return 'Transmissive_Momentum_Set_Stage_boundary(%s)' %self.domain

    def evaluate(self, vol_id, edge_id):
        """Transmissive Momentum Set Stage boundaries return the edge momentum
        values of the volume they serve.
        """

        q = self.domain.get_conserved_quantities(vol_id, edge = edge_id)
        value = self.function(self.domain.time)

        try:
            x = float(value)
        except:    
            x = float(value[0])
            
        q[0] = x
        return q


        #FIXME: Consider this (taken from File_boundary) to allow
        #spatial variation
        #if vol_id is not None and edge_id is not None:
        #    i = self.boundary_indices[ vol_id, edge_id ]
        #    return self.F(t, point_id = i)
        #else:
        #    return self.F(t)



class Dirichlet_Discharge_boundary(Boundary):
    """
    Sets stage (stage0)
    Sets momentum (wh0) in the inward normal direction.

    Underlying domain must be specified when boundary is instantiated
    """

    def __init__(self, domain = None, stage0=None, wh0=None):
        Boundary.__init__(self)

        if domain is None:
            msg = 'Domain must be specified for this type boundary'
            raise msg

        if stage0 is None:
            raise 'set stage'

        if wh0 is None:
            wh0 = 0.0

        self.domain   = domain
        self.stage0  = stage0
        self.wh0 = wh0

    def __repr__(self):
        return 'Dirichlet_Discharge_boundary(%s)' %self.domain

    def evaluate(self, vol_id, edge_id):
        """Set discharge in the (inward) normal direction
        """

        normal = self.domain.get_normal(vol_id,edge_id)
        q = [self.stage0, -self.wh0*normal[0], -self.wh0*normal[1]]
        return q


        #FIXME: Consider this (taken from File_boundary) to allow
        #spatial variation
        #if vol_id is not None and edge_id is not None:
        #    i = self.boundary_indices[ vol_id, edge_id ]
        #    return self.F(t, point_id = i)
        #else:
        #    return self.F(t)





#########################
#Standard forcing terms:
#
def gravity(domain):
    """Apply gravitational pull in the presence of bed slope
    """

    xmom = domain.quantities['xmomentum'].explicit_update
    ymom = domain.quantities['ymomentum'].explicit_update

    Stage = domain.quantities['stage']
    Elevation = domain.quantities['elevation']
    h = Stage.edge_values - Elevation.edge_values
    v = Elevation.vertex_values

    x = domain.get_vertex_coordinates()
    g = domain.g

    for k in range(domain.number_of_elements):
        avg_h = sum( h[k,:] )/3

        #Compute bed slope
        x0, y0, x1, y1, x2, y2 = x[k,:]
        z0, z1, z2 = v[k,:]

        zx, zy = gradient(x0, y0, x1, y1, x2, y2, z0, z1, z2)

        #Update momentum
        xmom[k] += -g*zx*avg_h
        ymom[k] += -g*zy*avg_h


def gravity_c(domain):
    """Wrapper calling C version
    """

    xmom = domain.quantities['xmomentum'].explicit_update
    ymom = domain.quantities['ymomentum'].explicit_update

    Stage = domain.quantities['stage']
    Elevation = domain.quantities['elevation']
    h = Stage.edge_values - Elevation.edge_values
    v = Elevation.vertex_values

    x = domain.get_vertex_coordinates()
    g = domain.g


    from shallow_water_ext import gravity
    gravity(g, h, v, x, xmom, ymom)


def manning_friction(domain):
    """Apply (Manning) friction to water momentum
    (Python version)
    """

    from math import sqrt

    w = domain.quantities['stage'].centroid_values
    z = domain.quantities['elevation'].centroid_values
    h = w-z

    uh = domain.quantities['xmomentum'].centroid_values
    vh = domain.quantities['ymomentum'].centroid_values
    eta = domain.quantities['friction'].centroid_values

    xmom_update = domain.quantities['xmomentum'].semi_implicit_update
    ymom_update = domain.quantities['ymomentum'].semi_implicit_update

    N = domain.number_of_elements
    eps = domain.minimum_allowed_height
    g = domain.g

    for k in range(N):
        if eta[k] >= eps:
            if h[k] >= eps:
                S = -g * eta[k]**2 * sqrt((uh[k]**2 + vh[k]**2))
                S /= h[k]**(7.0/3)

                #Update momentum
                xmom_update[k] += S*uh[k]
                ymom_update[k] += S*vh[k]


def manning_friction_implicit_c(domain):
    """Wrapper for c version
    """


    #print 'Implicit friction'

    xmom = domain.quantities['xmomentum']
    ymom = domain.quantities['ymomentum']

    w = domain.quantities['stage'].centroid_values
    z = domain.quantities['elevation'].centroid_values

    uh = xmom.centroid_values
    vh = ymom.centroid_values
    eta = domain.quantities['friction'].centroid_values

    xmom_update = xmom.semi_implicit_update
    ymom_update = ymom.semi_implicit_update

    N = domain.number_of_elements
    eps = domain.minimum_allowed_height
    g = domain.g

    from shallow_water_ext import manning_friction
    manning_friction(g, eps, w, z, uh, vh, eta, xmom_update, ymom_update)


def manning_friction_explicit_c(domain):
    """Wrapper for c version
    """

    #print 'Explicit friction'

    xmom = domain.quantities['xmomentum']
    ymom = domain.quantities['ymomentum']

    w = domain.quantities['stage'].centroid_values
    z = domain.quantities['elevation'].centroid_values

    uh = xmom.centroid_values
    vh = ymom.centroid_values
    eta = domain.quantities['friction'].centroid_values

    xmom_update = xmom.explicit_update
    ymom_update = ymom.explicit_update

    N = domain.number_of_elements
    eps = domain.minimum_allowed_height
    g = domain.g

    from shallow_water_ext import manning_friction
    manning_friction(g, eps, w, z, uh, vh, eta, xmom_update, ymom_update)


def linear_friction(domain):
    """Apply linear friction to water momentum

    Assumes quantity: 'linear_friction' to be present
    """

    from math import sqrt

    w = domain.quantities['stage'].centroid_values
    z = domain.quantities['elevation'].centroid_values
    h = w-z

    uh = domain.quantities['xmomentum'].centroid_values
    vh = domain.quantities['ymomentum'].centroid_values
    tau = domain.quantities['linear_friction'].centroid_values

    xmom_update = domain.quantities['xmomentum'].semi_implicit_update
    ymom_update = domain.quantities['ymomentum'].semi_implicit_update

    N = domain.number_of_elements
    eps = domain.minimum_allowed_height
    g = domain.g #Not necessary? Why was this added?

    for k in range(N):
        if tau[k] >= eps:
            if h[k] >= eps:
                S = -tau[k]/h[k]

                #Update momentum
                xmom_update[k] += S*uh[k]
                ymom_update[k] += S*vh[k]



def check_forcefield(f):
    """Check that f is either
    1: a callable object f(t,x,y), where x and y are vectors
       and that it returns an array or a list of same length
       as x and y
    2: a scalar
    """

    if callable(f):
        N = 3
        x = ones(3, Float)
        y = ones(3, Float)
        try:
            q = f(1.0, x=x, y=y)
        except Exception, e:
            msg = 'Function %s could not be executed:\n%s' %(f, e)
            #FIXME: Reconsider this semantics
            raise msg

        try:
            q = array(q).astype(Float)
        except:
            msg = 'Return value from vector function %s could ' %f
            msg += 'not be converted into a Numeric array of floats.\n'
            msg += 'Specified function should return either list or array.'
            raise msg

        #Is this really what we want?
        msg = 'Return vector from function %s ' %f
        msg += 'must have same lenght as input vectors'
        assert len(q) == N, msg

    else:
        try:
            f = float(f)
        except:
            msg = 'Force field %s must be either a scalar' %f
            msg += ' or a vector function'
            raise Exception(msg)
    return f


class Wind_stress:
    """Apply wind stress to water momentum in terms of
    wind speed [m/s] and wind direction [degrees]
    """

    def __init__(self, *args, **kwargs):
        """Initialise windfield from wind speed s [m/s]
        and wind direction phi [degrees]

        Inputs v and phi can be either scalars or Python functions, e.g.

        W = Wind_stress(10, 178)

        #FIXME - 'normal' degrees are assumed for now, i.e. the
        vector (1,0) has zero degrees.
        We may need to convert from 'compass' degrees later on and also
        map from True north to grid north.

        Arguments can also be Python functions of t,x,y as in

        def speed(t,x,y):
            ...
            return s

        def angle(t,x,y):
            ...
            return phi

        where x and y are vectors.

        and then pass the functions in

        W = Wind_stress(speed, angle)

        The instantiated object W can be appended to the list of
        forcing_terms as in

        Alternatively, one vector valued function for (speed, angle)
        can be applied, providing both quantities simultaneously.
        As in
        W = Wind_stress(F), where returns (speed, angle) for each t.

        domain.forcing_terms.append(W)
        """

        from anuga.config import rho_a, rho_w, eta_w
        from Numeric import array, Float

        if len(args) == 2:
            s = args[0]
            phi = args[1]
        elif len(args) == 1:
            #Assume vector function returning (s, phi)(t,x,y)
            vector_function = args[0]
            s = lambda t,x,y: vector_function(t,x=x,y=y)[0]
            phi = lambda t,x,y: vector_function(t,x=x,y=y)[1]
        else:
           #Assume info is in 2 keyword arguments

           if len(kwargs) == 2:
               s = kwargs['s']
               phi = kwargs['phi']
           else:
               raise 'Assumes two keyword arguments: s=..., phi=....'

        self.speed = check_forcefield(s)
        self.phi = check_forcefield(phi)

        self.const = eta_w*rho_a/rho_w


    def __call__(self, domain):
        """Evaluate windfield based on values found in domain
        """

        from math import pi, cos, sin, sqrt
        from Numeric import Float, ones, ArrayType

        xmom_update = domain.quantities['xmomentum'].explicit_update
        ymom_update = domain.quantities['ymomentum'].explicit_update

        N = domain.number_of_elements
        t = domain.time

        if callable(self.speed):
            xc = domain.get_centroid_coordinates()
            s_vec = self.speed(t, xc[:,0], xc[:,1])
        else:
            #Assume s is a scalar

            try:
                s_vec = self.speed * ones(N, Float)
            except:
                msg = 'Speed must be either callable or a scalar: %s' %self.s
                raise msg


        if callable(self.phi):
            xc = domain.get_centroid_coordinates()
            phi_vec = self.phi(t, xc[:,0], xc[:,1])
        else:
            #Assume phi is a scalar

            try:
                phi_vec = self.phi * ones(N, Float)
            except:
                msg = 'Angle must be either callable or a scalar: %s' %self.phi
                raise msg

        assign_windfield_values(xmom_update, ymom_update,
                                s_vec, phi_vec, self.const)


def assign_windfield_values(xmom_update, ymom_update,
                            s_vec, phi_vec, const):
    """Python version of assigning wind field to update vectors.
    A c version also exists (for speed)
    """
    from math import pi, cos, sin, sqrt

    N = len(s_vec)
    for k in range(N):
        s = s_vec[k]
        phi = phi_vec[k]

        #Convert to radians
        phi = phi*pi/180

        #Compute velocity vector (u, v)
        u = s*cos(phi)
        v = s*sin(phi)

        #Compute wind stress
        S = const * sqrt(u**2 + v**2)
        xmom_update[k] += S*u
        ymom_update[k] += S*v



##############################
#OBSOLETE STUFF

def balance_deep_and_shallow_old(domain):
    """Compute linear combination between stage as computed by
    gradient-limiters and stage computed as constant height above bed.
    The former takes precedence when heights are large compared to the
    bed slope while the latter takes precedence when heights are
    relatively small.  Anything in between is computed as a balanced
    linear combination in order to avoid numerical disturbances which
    would otherwise appear as a result of hard switching between
    modes.
    """

    #OBSOLETE

    #Shortcuts
    wc = domain.quantities['stage'].centroid_values
    zc = domain.quantities['elevation'].centroid_values
    hc = wc - zc

    wv = domain.quantities['stage'].vertex_values
    zv = domain.quantities['elevation'].vertex_values
    hv = wv-zv


    #Computed linear combination between constant stages and and
    #stages parallel to the bed elevation.
    for k in range(domain.number_of_elements):
        #Compute maximal variation in bed elevation
        #  This quantitiy is
        #    dz = max_i abs(z_i - z_c)
        #  and it is independent of dimension
        #  In the 1d case zc = (z0+z1)/2
        #  In the 2d case zc = (z0+z1+z2)/3

        dz = max(abs(zv[k,0]-zc[k]),
                 abs(zv[k,1]-zc[k]),
                 abs(zv[k,2]-zc[k]))


        hmin = min( hv[k,:] )

        #Create alpha in [0,1], where alpha==0 means using shallow
        #first order scheme and alpha==1 means using the stage w as
        #computed by the gradient limiter (1st or 2nd order)
        #
        #If hmin > dz/2 then alpha = 1 and the bed will have no effect
        #If hmin < 0 then alpha = 0 reverting to constant height above bed.

        if dz > 0.0:
            alpha = max( min( 2*hmin/dz, 1.0), 0.0 )
        else:
            #Flat bed
            alpha = 1.0

        #Weighted balance between stage parallel to bed elevation
        #(wvi = zvi + hc) and stage as computed by 1st or 2nd
        #order gradient limiter
        #(wvi = zvi + hvi) where i=0,1,2 denotes the vertex ids
        #
        #It follows that the updated wvi is
        #  wvi := (1-alpha)*(zvi+hc) + alpha*(zvi+hvi) =
        #  zvi + hc + alpha*(hvi - hc)
        #
        #Note that hvi = zc+hc-zvi in the first order case (constant).

        if alpha < 1:
            for i in range(3):
                wv[k,i] = zv[k,i] + hc[k] + alpha*(hv[k,i]-hc[k])


            #Momentums at centroids
            xmomc = domain.quantities['xmomentum'].centroid_values
            ymomc = domain.quantities['ymomentum'].centroid_values

            #Momentums at vertices
            xmomv = domain.quantities['xmomentum'].vertex_values
            ymomv = domain.quantities['ymomentum'].vertex_values

            # Update momentum as a linear combination of
            # xmomc and ymomc (shallow) and momentum
            # from extrapolator xmomv and ymomv (deep).
            xmomv[k,:] = (1-alpha)*xmomc[k] + alpha*xmomv[k,:]
            ymomv[k,:] = (1-alpha)*ymomc[k] + alpha*ymomv[k,:]





###########################
###########################
#Geometries


#FIXME: Rethink this way of creating values.


class Weir:
    """Set a bathymetry for weir with a hole and a downstream gutter
    x,y are assumed to be in the unit square
    """

    def __init__(self, stage):
        self.inflow_stage = stage

    def __call__(self, x, y):
        from Numeric import zeros, Float
        from math import sqrt

        N = len(x)
        assert N == len(y)

        z = zeros(N, Float)
        for i in range(N):
            z[i] = -x[i]/2  #General slope

            #Flattish bit to the left
            if x[i] < 0.3:
                z[i] = -x[i]/10

            #Weir
            if x[i] >= 0.3 and x[i] < 0.4:
                z[i] = -x[i]+0.9

            #Dip
            x0 = 0.6
            #depth = -1.3
            depth = -1.0
            #plateaux = -0.9
            plateaux = -0.6
            if y[i] < 0.7:
                if x[i] > x0 and x[i] < 0.9:
                    z[i] = depth

                #RHS plateaux
                if x[i] >= 0.9:
                    z[i] = plateaux


            elif y[i] >= 0.7 and y[i] < 1.5:
                #Restrict and deepen
                if x[i] >= x0 and x[i] < 0.8:
                    z[i] = depth-(y[i]/3-0.3)
                    #z[i] = depth-y[i]/5
                    #z[i] = depth
                elif x[i] >= 0.8:
                    #RHS plateaux
                    z[i] = plateaux

            elif y[i] >= 1.5:
                if x[i] >= x0 and x[i] < 0.8 + (y[i]-1.5)/1.2:
                    #Widen up and stay at constant depth
                    z[i] = depth-1.5/5
                elif x[i] >= 0.8 + (y[i]-1.5)/1.2:
                    #RHS plateaux
                    z[i] = plateaux


            #Hole in weir (slightly higher than inflow condition)
            if x[i] >= 0.3 and x[i] < 0.4 and y[i] > 0.2 and y[i] < 0.4:
                z[i] = -x[i]+self.inflow_stage + 0.02

            #Channel behind weir
            x0 = 0.5
            if x[i] >= 0.4 and x[i] < x0 and y[i] > 0.2 and y[i] < 0.4:
                z[i] = -x[i]+self.inflow_stage + 0.02

            if x[i] >= x0 and x[i] < 0.6 and y[i] > 0.2 and y[i] < 0.4:
                #Flatten it out towards the end
                z[i] = -x0+self.inflow_stage + 0.02 + (x0-x[i])/5

            #Hole to the east
            x0 = 1.1; y0 = 0.35
            #if x[i] < -0.2 and y < 0.5:
            if sqrt((2*(x[i]-x0))**2 + (2*(y[i]-y0))**2) < 0.2:
                z[i] = sqrt(((x[i]-x0))**2 + ((y[i]-y0))**2)-1.0

            #Tiny channel draining hole
            if x[i] >= 1.14 and x[i] < 1.2 and y[i] >= 0.4 and y[i] < 0.6:
                z[i] = -0.9 #North south

            if x[i] >= 0.9 and x[i] < 1.18 and y[i] >= 0.58 and y[i] < 0.65:
                z[i] = -1.0 + (x[i]-0.9)/3 #East west



            #Stuff not in use

            #Upward slope at inlet to the north west
            #if x[i] < 0.0: # and y[i] > 0.5:
            #    #z[i] = -y[i]+0.5  #-x[i]/2
            #    z[i] = x[i]/4 - y[i]**2 + 0.5

            #Hole to the west
            #x0 = -0.4; y0 = 0.35 # center
            #if sqrt((2*(x[i]-x0))**2 + (2*(y[i]-y0))**2) < 0.2:
            #    z[i] = sqrt(((x[i]-x0))**2 + ((y[i]-y0))**2)-0.2





        return z/2

class Weir_simple:
    """Set a bathymetry for weir with a hole and a downstream gutter
    x,y are assumed to be in the unit square
    """

    def __init__(self, stage):
        self.inflow_stage = stage

    def __call__(self, x, y):
        from Numeric import zeros, Float

        N = len(x)
        assert N == len(y)

        z = zeros(N, Float)
        for i in range(N):
            z[i] = -x[i]  #General slope

            #Flat bit to the left
            if x[i] < 0.3:
                z[i] = -x[i]/10  #General slope

            #Weir
            if x[i] > 0.3 and x[i] < 0.4:
                z[i] = -x[i]+0.9

            #Dip
            if x[i] > 0.6 and x[i] < 0.9:
                z[i] = -x[i]-0.5  #-y[i]/5

            #Hole in weir (slightly higher than inflow condition)
            if x[i] > 0.3 and x[i] < 0.4 and y[i] > 0.2 and y[i] < 0.4:
                z[i] = -x[i]+self.inflow_stage + 0.05


        return z/2



class Constant_stage:
    """Set an initial condition with constant stage
    """
    def __init__(self, s):
        self.s = s

    def __call__(self, x, y):
        return self.s

class Constant_height:
    """Set an initial condition with constant water height, e.g
    stage s = z+h
    """

    def __init__(self, W, h):
        self.W = W
        self.h = h

    def __call__(self, x, y):
        if self.W is None:
            from Numeric import ones, Float
            return self.h*ones(len(x), Float)
        else:
            return self.W(x,y) + self.h




def compute_fluxes_python(domain):
    """Compute all fluxes and the timestep suitable for all volumes
    in domain.

    Compute total flux for each conserved quantity using "flux_function"

    Fluxes across each edge are scaled by edgelengths and summed up
    Resulting flux is then scaled by area and stored in
    explicit_update for each of the three conserved quantities
    stage, xmomentum and ymomentum

    The maximal allowable speed computed by the flux_function for each volume
    is converted to a timestep that must not be exceeded. The minimum of
    those is computed as the next overall timestep.

    Post conditions:
      domain.explicit_update is reset to computed flux values
      domain.timestep is set to the largest step satisfying all volumes.
    """

    import sys
    from Numeric import zeros, Float

    N = domain.number_of_elements

    #Shortcuts
    Stage = domain.quantities['stage']
    Xmom = domain.quantities['xmomentum']
    Ymom = domain.quantities['ymomentum']
    Bed = domain.quantities['elevation']

    #Arrays
    stage = Stage.edge_values
    xmom =  Xmom.edge_values
    ymom =  Ymom.edge_values
    bed =   Bed.edge_values

    stage_bdry = Stage.boundary_values
    xmom_bdry =  Xmom.boundary_values
    ymom_bdry =  Ymom.boundary_values

    flux = zeros((N,3), Float) #Work array for summing up fluxes

    #Loop
    timestep = float(sys.maxint)
    for k in range(N):

        for i in range(3):
            #Quantities inside volume facing neighbour i
            ql = [stage[k, i], xmom[k, i], ymom[k, i]]
            zl = bed[k, i]

            #Quantities at neighbour on nearest face
            n = domain.neighbours[k,i]
            if n < 0:
                m = -n-1 #Convert negative flag to index
                qr = [stage_bdry[m], xmom_bdry[m], ymom_bdry[m]]
                zr = zl #Extend bed elevation to boundary
            else:
                m = domain.neighbour_edges[k,i]
                qr = [stage[n, m], xmom[n, m], ymom[n, m]]
                zr = bed[n, m]


            #Outward pointing normal vector
            normal = domain.normals[k, 2*i:2*i+2]

            #Flux computation using provided function
            edgeflux, max_speed = flux_function(normal, ql, qr, zl, zr)

            flux[k,:] = edgeflux

    return flux







##############################################
#Initialise module


from anuga.utilities import compile
if compile.can_use_C_extension('shallow_water_ext.c'):
    #Replace python version with c implementations

    from shallow_water_ext import rotate, assign_windfield_values
    compute_fluxes = compute_fluxes_c
    extrapolate_second_order_sw=extrapolate_second_order_sw_c
    gravity = gravity_c
    manning_friction = manning_friction_implicit_c
    h_limiter = h_limiter_c
    balance_deep_and_shallow = balance_deep_and_shallow_c
    protect_against_infinitesimal_and_negative_heights =\
                    protect_against_infinitesimal_and_negative_heights_c


    #distribute_to_vertices_and_edges =\
    #              distribute_to_vertices_and_edges_c #(like MH's)



#Optimisation with psyco
from anuga.config import use_psyco
if use_psyco:
    try:
        import psyco
    except:
        import os
        if os.name == 'posix' and os.uname()[4] == 'x86_64':
            pass
            #Psyco isn't supported on 64 bit systems, but it doesn't matter
        else:
            msg = 'WARNING: psyco (speedup) could not import'+\
                  ', you may want to consider installing it'
            print msg
    else:
        psyco.bind(Domain.distribute_to_vertices_and_edges)
        psyco.bind(Domain.compute_fluxes)

if __name__ == "__main__":
    pass

# Profiling stuff
import profile
profiler = profile.Profile()

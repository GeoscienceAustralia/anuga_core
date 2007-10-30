"""Finite-volume computations of the shallow water wave equation.

Title: ANGUA shallow_water_domain - 2D triangular domains for finite-volume
       computations of the shallow water wave equation.


Author: Ole Nielsen (Ole.Nielsen@ga.gov.au),
        Stephen Roberts (Stephen.Roberts@anu.edu.au),
        Duncan Gray (Duncan.Gray@ga.gov.au), etc

CreationDate: 2004

Description:
    This module contains a specialisation of class Domain from
    module domain.py consisting of methods specific to the
    Shallow Water Wave Equation

    U_t + E_x + G_y = S

    where

    U = [w, uh, vh]
    E = [uh, u^2h + gh^2/2, uvh]
    G = [vh, uvh, v^2h + gh^2/2]
    S represents source terms forcing the system
    (e.g. gravity, friction, wind stress, ...)

    and _t, _x, _y denote the derivative with respect to t, x and y
    respectively.


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

Reference:
    Catastrophic Collapse of Water Supply Reservoirs in Urban Areas,
    Christopher Zoppou and Stephen Roberts,
    Journal of Hydraulic Engineering, vol. 127, No. 7 July 1999

    Hydrodynamic modelling of coastal inundation. 
    Nielsen, O., S. Roberts, D. Gray, A. McPherson and A. Hitchman
    In Zerger, A. and Argent, R.M. (eds) MODSIM 2005 International Congress on
    Modelling and Simulation. Modelling and Simulation Society of Australia and
    New Zealand, December 2005, pp. 518-523. ISBN: 0-9758400-2-9.
    http://www.mssanz.org.au/modsim05/papers/nielsen.pdf


SeeAlso:
    TRAC administration of ANUGA (User Manuals etc) at
    https://datamining.anu.edu.au/anuga and Subversion repository at
    $HeadURL$

Constraints: See GPL license in the user guide
Version: 1.0 ($Revision$)
ModifiedBy:
    $Author$
    $Date$

"""

# Subversion keywords:
#
# $LastChangedDate$
# $LastChangedRevision$
# $LastChangedBy$

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
from anuga.config import g, epsilon, beta_h, beta_w, beta_w_dry,\
     beta_uh, beta_uh_dry, beta_vh, beta_vh_dry, tight_slope_limiters
from anuga.config import alpha_balance
from anuga.config import optimise_dry_cells

#---------------------
# Shallow water domain
#---------------------
class Domain(Generic_Domain):

    conserved_quantities = ['stage', 'xmomentum', 'ymomentum']
    other_quantities = ['elevation', 'friction']
    
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
                 number_of_full_nodes=None,
                 number_of_full_triangles=None):


        other_quantities = ['elevation', 'friction']
        Generic_Domain.__init__(self,
                                coordinates,
                                vertices,
                                boundary,
                                Domain.conserved_quantities,
                                Domain.other_quantities,
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


        self.set_minimum_allowed_height(minimum_allowed_height)
        
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

        self.tight_slope_limiters = tight_slope_limiters
        self.optimise_dry_cells = optimise_dry_cells

        self.forcing_terms.append(manning_friction_implicit)
        self.forcing_terms.append(gravity)

        # Stored output
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

        # FIXME (Ole): how about using the word continuous vertex values?
        self.smooth = not flag

        # Reduction operation for get_vertex_values
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


    def set_minimum_allowed_height(self, minimum_allowed_height):
        """
        Set the minimum depth that will be recognised in the numerical
        scheme

        The minimum allowed depth is in meters.

        The parameter H0 (Minimal height for flux computation)
        is also set by this function
        """

        #FIXME (Ole): rename H0 to minimum_allowed_height_in_flux_computation

        #FIXME (Ole): Maybe use histogram to identify isolated extreme speeds
        #and deal with them adaptively similarly to how we used to use 1 order
        #steps to recover.
        self.minimum_allowed_height = minimum_allowed_height
        self.H0 = minimum_allowed_height   
        

    def set_maximum_allowed_speed(self, maximum_allowed_speed):
        """
        Set the maximum particle speed that is allowed in water
        shallower than minimum_allowed_height. This is useful for
        controlling speeds in very thin layers of water and at the same time
        allow some movement avoiding pooling of water.

        """
        self.maximum_allowed_speed = maximum_allowed_speed


    def set_points_file_block_line_size(self,points_file_block_line_size):
        """
        Set the minimum depth that will be recognised when writing
        to an sww file. This is useful for removing thin water layers
        that seems to be caused by friction creep.

        The minimum allowed sww depth is in meters.
        """
        self.points_file_block_line_size = points_file_block_line_size
        
        
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

        # Check correcness
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
        """Return location of highest elevation where h > 0

        Optional argument:
            indices is the set of element ids that the operation applies to.

        Usage:
            q = get_maximum_inundation_location()

        Note, centroid values are used for this operation            
        """

        wet_elements = self.get_wet_elements(indices)
        return self.get_quantity('elevation').\
               get_maximum_location(indices=wet_elements)    

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
        # Call correct module function
        # (either from this module or C-extension)
        distribute_to_vertices_and_edges(self)



    def evolve(self,
               yieldstep = None,
               finaltime = None,
               duration = None,
               skip_initial_step = False):
        """Specialisation of basic evolve method from parent class
        """

        # Call check integrity here rather than from user scripts
        # self.check_integrity()

        msg = 'Parameter beta_h must be in the interval [0, 1['
        assert 0 <= self.beta_h <= 1.0, msg
        msg = 'Parameter beta_w must be in the interval [0, 1['
        assert 0 <= self.beta_w <= 1.0, msg


        # Initial update of vertex and edge values before any STORAGE
        # and or visualisation
        # This is done again in the initialisation of the Generic_Domain
        # evolve loop but we do it here to ensure the values are ok for storage
        self.distribute_to_vertices_and_edges()

        if self.store is True and self.time == 0.0:
            self.initialise_storage()
            # print 'Storing results in ' + self.writer.filename
        else:
            pass
            # print 'Results will not be stored.'
            # print 'To store results set domain.store = True'
            # FIXME: Diagnostic output should be controlled by
            # a 'verbose' flag living in domain (or in a parent class)

        # Call basic machinery from parent class
        for t in Generic_Domain.evolve(self,
                                       yieldstep=yieldstep,
                                       finaltime=finaltime,
                                       duration=duration,
                                       skip_initial_step=skip_initial_step):

            # Store model data, e.g. for subsequent visualisation
            if self.store is True:
                self.store_timestep(self.quantities_to_be_stored)

            # FIXME: Could maybe be taken from specified list
            # of 'store every step' quantities

            # Pass control on to outer loop for more specific actions
            yield(t)


    def initialise_storage(self):
        """Create and initialise self.writer object for storing data.
        Also, save x,y and bed elevation
        """

        from anuga.shallow_water.data_manager import get_dataobject

        # Initialise writer
        self.writer = get_dataobject(self, mode = 'w')

        # Store vertices and connectivity
        self.writer.store_connectivity()


    def store_timestep(self, name):
        """Store named quantity and time.

        Precondition:
           self.write has been initialised
        """
        self.writer.store_timestep(name)


#=============== End of class Shallow Water Domain ===============================


#-----------------
# Flux computation
#-----------------

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
    

    This wrapper calls the underlying C version of compute fluxes
    """

    import sys

    N = len(domain) # number_of_triangles

    # Shortcuts
    Stage = domain.quantities['stage']
    Xmom = domain.quantities['xmomentum']
    Ymom = domain.quantities['ymomentum']
    Bed = domain.quantities['elevation']

    timestep = float(sys.maxint)
    from shallow_water_ext import\
         compute_fluxes_ext_central as compute_fluxes_ext


    flux_timestep = compute_fluxes_ext(timestep,
                                       domain.epsilon,
                                       domain.H0,
                                       domain.g,
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
                                       domain.already_computed_flux,
                                       domain.max_speed,
                                       int(domain.optimise_dry_cells))

    domain.flux_timestep = flux_timestep



#---------------------------------------
# Module functions for gradient limiting
#---------------------------------------


# MH090605 The following method belongs to the shallow_water domain class
# see comments in the corresponding method in shallow_water_ext.c
def extrapolate_second_order_sw(domain):
    """Wrapper calling C version of extrapolate_second_order_sw
    """
    import sys

    N = len(domain) # number_of_triangles

    # Shortcuts
    Stage = domain.quantities['stage']
    Xmom = domain.quantities['xmomentum']
    Ymom = domain.quantities['ymomentum']
    Elevation = domain.quantities['elevation']

    from shallow_water_ext import extrapolate_second_order_sw as extrapol2
    extrapol2(domain,
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
              Elevation.vertex_values,
              int(domain.optimise_dry_cells))


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

    # Remove very thin layers of water
    protect_against_infinitesimal_and_negative_heights(domain)

    # Extrapolate all conserved quantities
    if optimised_gradient_limiter:
        # MH090605 if second order,
        # perform the extrapolation and limiting on
        # all of the conserved quantities

        if (domain._order_ == 1):
            for name in domain.conserved_quantities:
                Q = domain.quantities[name]
                Q.extrapolate_first_order()
        elif domain._order_ == 2:
            domain.extrapolate_second_order_sw()
        else:
            raise 'Unknown order'
    else:
        # Old code:
        for name in domain.conserved_quantities:
            Q = domain.quantities[name]

            if domain._order_ == 1:
                Q.extrapolate_first_order()
            elif domain._order_ == 2:

                # Experiment
                #if name == 'stage':
                #    #print name, 'second'
                #    Q.extrapolate_second_order()
                #    Q.limit()
                #else:
                #    #print name, 'first'                
                #    Q.extrapolate_first_order()
                #    #Q.extrapolate_second_order()
                #    #Q.limit()                
                
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

    # Shortcuts
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

    Wrapper for c-extension
    """

    N = len(domain) # number_of_triangles
    beta_h = domain.beta_h

    # Shortcuts
    wc = domain.quantities['stage'].centroid_values
    zc = domain.quantities['elevation'].centroid_values
    hc = wc - zc

    wv = domain.quantities['stage'].vertex_values
    zv = domain.quantities['elevation'].vertex_values
    hv = wv - zv

    #Call C-extension
    from shallow_water_ext import h_limiter_sw
    hvbar = h_limiter_sw(domain, hc, hv)

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

    Wrapper for C implementation
    """

    # FIXME (Ole): I reckon this can be simplified significantly:
    #
    # Always use beta_h == 0, and phase it out.
    # Compute hc and hv in the c-code
    # Omit updating xmomv 
    #
    from shallow_water_ext import balance_deep_and_shallow as balance_deep_and_shallow_c
    
    # Shortcuts
    wc = domain.quantities['stage'].centroid_values
    zc = domain.quantities['elevation'].centroid_values

    wv = domain.quantities['stage'].vertex_values
    zv = domain.quantities['elevation'].vertex_values

    # Momentums at centroids
    xmomc = domain.quantities['xmomentum'].centroid_values
    ymomc = domain.quantities['ymomentum'].centroid_values

    # Momentums at vertices
    xmomv = domain.quantities['xmomentum'].vertex_values
    ymomv = domain.quantities['ymomentum'].vertex_values

    # Limit h
    if domain.beta_h > 0:
        hvbar = h_limiter(domain)
        
        balance_deep_and_shallow_c(domain, domain.beta_h,
                                   wc, zc, wv, zv, hvbar,
                                   xmomc, ymomc, xmomv, ymomv)        
    else:
        # print 'Using first order h-limiter'
        # FIXME: Pass wc in for now - it will be ignored.
        
        # This is how one would make a first order h_limited value
        # as in the old balancer (pre 17 Feb 2005):
        #  If we wish to hard wire this, one should modify the C-code
        # from Numeric import zeros, Float
        # hvbar = zeros( (len(wc), 3), Float)
        # for i in range(3):
        #     hvbar[:,i] = wc[:] - zc[:]

        balance_deep_and_shallow_c(domain, domain.beta_h,
                                   wc, zc, wv, zv, wc, 
                                   xmomc, ymomc, xmomv, ymomv)




#------------------------------------------------------------------
# Boundary conditions - specific to the shallow water wave equation
#------------------------------------------------------------------
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

        # Handy shorthands
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


        # FIXME: Consider this (taken from File_boundary) to allow
        # spatial variation
        # if vol_id is not None and edge_id is not None:
        #     i = self.boundary_indices[ vol_id, edge_id ]
        #     return self.F(t, point_id = i)
        # else:
        #     return self.F(t)



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


        # FIXME: Consider this (taken from File_boundary) to allow
        # spatial variation
        # if vol_id is not None and edge_id is not None:
        #     i = self.boundary_indices[ vol_id, edge_id ]
        #     return self.F(t, point_id = i)
        # else:
        #     return self.F(t)


class Field_boundary(Boundary):
    """Set boundary from given field represented in an sww file containing values
    for stage, xmomentum and ymomentum.
    Optionally, the user can specify mean_stage to offset the stage provided in the
    sww file.

    This function is a thin wrapper around the generic File_boundary. The 
    difference between the file_boundary and field_boundary is only that the
    field_boundary will allow you to change the level of the stage height when
    you read in the boundary condition. This is very useful when running 
    different tide heights in the same area as you need only to convert one 
    boundary condition to a SWW file, ideally for tide height of 0 m 
    (saving disk space). Then you can use field_boundary to read this SWW file
    and change the stage height (tide) on the fly depending on the scenario.
    
    """


    def __init__(self, filename, domain,
                 mean_stage=0.0,
                 time_thinning=1, 
                 use_cache=False,
                 verbose=False):
        """Constructor

        filename: Name of sww file
        domain: pointer to shallow water domain for which the boundary applies
        mean_stage: The mean water level which will be added to stage derived
                    from the sww file
        time_thinning: Will set how many time steps from the sww file read in
                       will be interpolated to the boundary. For example if 
                       the sww file has 1 second time steps and is 24 hours
                       in length it has 86400 time steps. If you set 
                       time_thinning to 1 it will read all these steps. 
                       If you set it to 100 it will read every 100th step eg
                       only 864 step. This parameter is very useful to increase
                       the speed of a model run that you are setting up 
                       and testing.
        use_cache:
        verbose:
        
        """

        # Create generic file_boundary object
        self.file_boundary = File_boundary(filename, domain,
                                           time_thinning=time_thinning,
                                           use_cache=use_cache,
                                           verbose=verbose)
        
        # Record information from File_boundary
        self.F = self.file_boundary.F
        self.domain = self.file_boundary.domain
        
        # Record mean stage
        self.mean_stage = mean_stage


    def __repr__(self):
        return 'Field boundary'


    def evaluate(self, vol_id=None, edge_id=None):
        """Return linearly interpolated values based on domain.time

        vol_id and edge_id are ignored
        """

        # Evaluate file boundary
        q = self.file_boundary.evaluate(vol_id, edge_id)

        # Adjust stage
        for j, name in enumerate(self.domain.conserved_quantities):
            if name == 'stage':
                q[j] += self.mean_stage
        return q

    

#-----------------------
# Standard forcing terms
#-----------------------

def gravity(domain):
    """Apply gravitational pull in the presence of bed slope
    Wrapper calls underlying C implementation
    """

    xmom = domain.quantities['xmomentum'].explicit_update
    ymom = domain.quantities['ymomentum'].explicit_update

    stage = domain.quantities['stage']
    elevation = domain.quantities['elevation']

    h = stage.centroid_values - elevation.centroid_values
    z = elevation.vertex_values

    x = domain.get_vertex_coordinates()
    g = domain.g
    

    from shallow_water_ext import gravity as gravity_c
    gravity_c(g, h, z, x, xmom, ymom) #, 1.0e-6)



def manning_friction_implicit(domain):
    """Apply (Manning) friction to water momentum    
    Wrapper for c version
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

    N = len(domain)
    eps = domain.minimum_allowed_height
    g = domain.g

    from shallow_water_ext import manning_friction as manning_friction_c
    manning_friction_c(g, eps, w, z, uh, vh, eta, xmom_update, ymom_update)


def manning_friction_explicit(domain):
    """Apply (Manning) friction to water momentum    
    Wrapper for c version
    """

    # print 'Explicit friction'

    xmom = domain.quantities['xmomentum']
    ymom = domain.quantities['ymomentum']

    w = domain.quantities['stage'].centroid_values
    z = domain.quantities['elevation'].centroid_values

    uh = xmom.centroid_values
    vh = ymom.centroid_values
    eta = domain.quantities['friction'].centroid_values

    xmom_update = xmom.explicit_update
    ymom_update = ymom.explicit_update

    N = len(domain)
    eps = domain.minimum_allowed_height
    g = domain.g

    from shallow_water_ext import manning_friction as manning_friction_c
    manning_friction_c(g, eps, w, z, uh, vh, eta, xmom_update, ymom_update)


# FIXME (Ole): This was implemented for use with one of the analytical solutions (Sampson?)
# Is it still needed (30 Oct 2007)?
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

    N = len(domain) # number_of_triangles
    eps = domain.minimum_allowed_height
    g = domain.g #Not necessary? Why was this added?

    for k in range(N):
        if tau[k] >= eps:
            if h[k] >= eps:
                S = -tau[k]/h[k]

                #Update momentum
                xmom_update[k] += S*uh[k]
                ymom_update[k] += S*vh[k]



#---------------------------------
# Experimental auxiliary functions
#---------------------------------
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
            # FIXME: Reconsider this semantics
            raise msg

        try:
            q = array(q).astype(Float)
        except:
            msg = 'Return value from vector function %s could ' %f
            msg += 'not be converted into a Numeric array of floats.\n'
            msg += 'Specified function should return either list or array.'
            raise msg

        # Is this really what we want?
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
            # Assume vector function returning (s, phi)(t,x,y)
            vector_function = args[0]
            s = lambda t,x,y: vector_function(t,x=x,y=y)[0]
            phi = lambda t,x,y: vector_function(t,x=x,y=y)[1]
        else:
           # Assume info is in 2 keyword arguments

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

        N = len(domain) # number_of_triangles
        t = domain.time

        if callable(self.speed):
            xc = domain.get_centroid_coordinates()
            s_vec = self.speed(t, xc[:,0], xc[:,1])
        else:
            # Assume s is a scalar

            try:
                s_vec = self.speed * ones(N, Float)
            except:
                msg = 'Speed must be either callable or a scalar: %s' %self.s
                raise msg


        if callable(self.phi):
            xc = domain.get_centroid_coordinates()
            phi_vec = self.phi(t, xc[:,0], xc[:,1])
        else:
            # Assume phi is a scalar

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

        # Convert to radians
        phi = phi*pi/180

        # Compute velocity vector (u, v)
        u = s*cos(phi)
        v = s*sin(phi)

        # Compute wind stress
        S = const * sqrt(u**2 + v**2)
        xmom_update[k] += S*u
        ymom_update[k] += S*v



class Rainfall:
    """Class Rainfall - general 'rain over entire domain' forcing term.
    
    Used for implementing Rainfall over the entire domain.
	
	Current Limited to only One Gauge..
	
	Need to add Spatial Varying Capability 
	(This module came from copying and amending the Inflow Code)
    
    Rainfall(rain)
        
    rain [mm/s]:  Total rain rate over the specified domain.  
	          NOTE: Raingauge Data needs to reflect the time step.
		  IE: if Gauge is mm read at a time step, then the input
                  here is as mm/(timeStep) so 10mm in 5minutes becomes
                  10/(5x60) = 0.0333mm/s.
	
	
                  This parameter can be either a constant or a
                  function of time. Positive values indicate inflow, 
                  negative values indicate outflow.
                  (and be used for Infiltration - Write Seperate Module)
                  The specified flow will be divided by the area of
                  the inflow region and then applied to update the
                  quantity in question. 
    
    Examples
    How to put them in a run File...
	
    #--------------------------------------------------------------------------
    # Setup specialised forcing terms
    #--------------------------------------------------------------------------
    # This is the new element implemented by Ole and Rudy to allow direct
    # input of Inflow in mm/s

    catchmentrainfall = Rainfall(rain=file_function('Q100_2hr_Rain.tms'))  
                        # Note need path to File in String.
                        # Else assumed in same directory

    domain.forcing_terms.append(catchmentrainfall)
    """

    # FIXME (OLE): Add a polygon as an alternative.
    # FIXME (AnyOne) : Add various methods to allow spatial variations
    # FIXME (OLE): Generalise to all quantities

    def __init__(self,
		 rain=0.0,
		 quantity_name='stage'):

        self.rain = rain
	self.quantity_name = quantity_name
    
    def __call__(self, domain):

        # Update rainfall
	if callable(self.rain):
	    rain = self.rain(domain.get_time())
	else:
	    rain = self.rain

        # Now rain is a number
        quantity = domain.quantities[self.quantity_name].explicit_update
        quantity[:] += rain/1000  # Converting mm/s to m/s to apply in ANUGA
		# 1mm of rainfall is equivalent to 1 litre /m2 
		# Flow is expressed as m3/s converted to a stage height in (m)
		
		# Note 1m3 = 1x10^9mm3 (mls)
		# or is that m3 to Litres ??? Check this how is it applied !!!


class Inflow:
    """Class Inflow - general 'rain and drain' forcing term.
    
    Useful for implementing flows in and out of the domain.
    
    Inflow(center, radius, flow)
    
    center [m]: Coordinates at center of flow point
    radius [m]: Size of circular area
    flow [m^3/s]: Total flow rate over the specified area.  
                  This parameter can be either a constant or a
                  function of time. Positive values indicate inflow, 
                  negative values indicate outflow.
                  The specified flow will be divided by the area of
                  the inflow region and then applied to update the
                  quantity in question. 
    
    Examples

    # Constant drain at 0.003 m^3/s.
    # The outflow area is 0.07**2*pi=0.0154 m^2
    # This corresponds to a rate of change of 0.003/0.0154 = 0.2 m/s 
    #				          
    Inflow((0.7, 0.4), 0.07, -0.003)


    # Tap turning up to a maximum inflow of 0.0142 m^3/s.
    # The inflow area is 0.03**2*pi = 0.00283 m^2
    # This corresponds to a rate of change of 0.0142/0.00283 = 5 m/s     
    # over the specified area
    Inflow((0.5, 0.5), 0.03, lambda t: min(0.01*t, 0.0142))

    #--------------------------------------------------------------------------
    # Setup specialised forcing terms
    #--------------------------------------------------------------------------
    # This is the new element implemented by Ole to allow direct input
    # of Inflow in m^3/s

    hydrograph = Inflow(center=(320, 300), radius=10,
                        flow=file_function('Q/QPMF_Rot_Sub13.tms'))

    domain.forcing_terms.append(hydrograph)
    
    """

    # FIXME (OLE): Add a polygon as an alternative.
    # FIXME (OLE): Generalise to all quantities

    def __init__(self,
		 center=None, radius=None,
		 flow=0.0,
		 quantity_name = 'stage'):

        from math import pi

        
		 
        if center is not None and radius is not None:
            assert len(center) == 2
	else:
	    msg = 'Both center and radius must be specified'
	    raise Exception, msg
    
        self.center = center
	self.radius = radius
        self.area = radius**2*pi
        self.flow = flow
	self.quantity_name = quantity_name
    
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
    	 	if ((x-center[0])**2+(y-center[1])**2) < radius**2:
		    self.indices.append(k)    

        # Update inflow
	if callable(self.flow):
	    flow = self.flow(domain.get_time())
	else:
	    flow = self.flow

        # Now flow is a number
        
        quantity = domain.quantities[self.quantity_name].explicit_update
        for k in self.indices:
            quantity[k] += flow/self.area		        


#------------------
# Initialise module
#------------------


from anuga.utilities import compile
if compile.can_use_C_extension('shallow_water_ext.c'):
    # Underlying C implementations can be accessed 

    from shallow_water_ext import rotate, assign_windfield_values
else:
    msg = 'C implementations could not be accessed by %s.\n ' %__file__
    msg += 'Make sure compile_all.py has been run as described in '
    msg += 'the ANUGA installation guide.'
    raise Exception, msg


# Optimisation with psyco
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



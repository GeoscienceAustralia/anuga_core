"""
Finite-volume computations of the shallow water wave equation.

Title: ANGUA shallow_water_domain - 2D triangular domains for finite-volume
       computations of the shallow water wave equation.


Author: Ole Nielsen, Ole.Nielsen@ga.gov.au
        Stephen Roberts, Stephen.Roberts@anu.edu.au
        Duncan Gray, Duncan.Gray@ga.gov.au

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
    $HeadURL: https://datamining.anu.edu.au/svn/anuga/trunk/anuga_core/source/
                anuga/shallow_water/shallow_water_domain.py $

Constraints: See GPL license in the user guide
Version: 1.0 ($Revision$)
ModifiedBy:
    $Author$
    $Date$
"""


import numpy as num

from anuga.abstract_2d_finite_volumes.generic_domain \
                    import Generic_Domain

from anuga.shallow_water.forcing import Cross_section
from anuga.utilities.numerical_tools import mean
from anuga.file.sww import SWW_file  
            
import anuga.utilities.log as log


class Domain(Generic_Domain):
    """ Class for a shallow water domain."""
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
                 conserved_quantities = None,
                 evolved_quantities = None,
                 other_quantities = None,
                 full_send_dict=None,
                 ghost_recv_dict=None,
                 starttime=0,
                 processor=0,
                 numproc=1,
                 number_of_full_nodes=None,
                 number_of_full_triangles=None):
        """
            Instantiate a shallow water domain.
            coordinates - vertex locations for the mesh
            vertices - vertex indices for the mesh
            boundary - boundaries of the mesh
            tagged_elements
            geo_reference
            use_inscribed_circle
            mesh_filename
            use_cache
            verbose
            evolved_quantities
            full_send_dict
            ghost_recv_dict
            processor
            numproc
            number_of_full_nodes
            number_of_full_triangles
        """

        # Define quantities for the shallow_water domain
        if conserved_quantities == None:
            conserved_quantities = ['stage', 'xmomentum', 'ymomentum']

        if evolved_quantities == None:
            evolved_quantities =  ['stage', 'xmomentum', 'ymomentum']
            
        if other_quantities == None:
            other_quantities = ['elevation', 'friction', 'height', 'xvelocity', 'yvelocity', 'x', 'y']
        
        
        Generic_Domain.__init__(self,
                            coordinates,
                            vertices,
                            boundary,
                            conserved_quantities,
                            evolved_quantities,
                            other_quantities,
                            tagged_elements,
                            geo_reference,
                            use_inscribed_circle,
                            mesh_filename,
                            use_cache,
                            verbose,
                            full_send_dict,
                            ghost_recv_dict,
                            starttime,
                            processor,
                            numproc,
                            number_of_full_nodes=number_of_full_nodes,
                            number_of_full_triangles=number_of_full_triangles)

        self.set_defaults()

        #-------------------------------
        # Operators and Forcing Terms
        #-------------------------------
        self.forcing_terms.append(manning_friction_implicit)
        self.forcing_terms.append(gravity)


        self.fractional_step_operators = []
        self.kv_operator = None

        #-------------------------------
        # Stored output
        #-------------------------------
        self.set_store(True)
        self.set_store_vertices_uniquely(False)
        self.quantities_to_be_stored = {'elevation': 1, 
                                        'stage': 2, 
                                        'xmomentum': 2, 
                                        'ymomentum': 2}



        #-------------------------------
        # Useful auxiliary quantity
        #-------------------------------
        n = len(self)
        self.quantities['x'].set_values(self.vertex_coordinates[:,0].reshape(n,3))
        self.quantities['x'].set_boundary_values_from_edges()
        
        self.quantities['y'].set_values(self.vertex_coordinates[:,1].reshape(n,3))
        self.quantities['y'].set_boundary_values_from_edges()



    def set_defaults(self):
        """Set the default values in this routine. That way we can inherit class
        and just over redefine the defaults for the new class
        """

        from anuga.config import minimum_storable_height
        from anuga.config import minimum_allowed_height, maximum_allowed_speed
        from anuga.config import g, beta_w, beta_w_dry, \
             beta_uh, beta_uh_dry, beta_vh, beta_vh_dry, tight_slope_limiters
        from anuga.config import alpha_balance
        from anuga.config import optimise_dry_cells
        from anuga.config import optimised_gradient_limiter
        from anuga.config import use_edge_limiter
        from anuga.config import use_centroid_velocities

        self.set_minimum_allowed_height(minimum_allowed_height)
        self.maximum_allowed_speed = maximum_allowed_speed

        self.g = g
        self.beta_w = beta_w
        self.beta_w_dry = beta_w_dry
        self.beta_uh = beta_uh
        self.beta_uh_dry = beta_uh_dry
        self.beta_vh = beta_vh
        self.beta_vh_dry = beta_vh_dry
        self.alpha_balance = alpha_balance

        self.tight_slope_limiters = tight_slope_limiters
        self.optimise_dry_cells = optimise_dry_cells

        
        self.set_sloped_mannings_function(False)

        self.minimum_storable_height = minimum_storable_height

         # Limiters
        self.use_edge_limiter = use_edge_limiter
        self.optimised_gradient_limiter = optimised_gradient_limiter
        self.use_centroid_velocities = use_centroid_velocities       


    def set_store(self, flag=True):
        """Set whether data saved to sww file.
        """

        self.store = flag

        
    def set_sloped_mannings_function(self, flag=True):
        """Set mannings friction function to use the sloped
        wetted area.
        The flag is tested in the python wrapper
        mannings_friction_implicit
        """
        if flag:
            self.use_sloped_mannings = True
        else:
            self.use_sloped_mannings = False


    def set_use_edge_limiter(self, flag=True):
        """Cludge to allow unit test to pass, but to
        also introduce new edge limiting. The flag is
        tested in distribute_to_vertices_and_edges
        """
        if flag:
            self.use_edge_limiter = True
        else:
            self.use_edge_limiter = False


    def set_use_kinematic_viscosity(self, flag=False):

        from anuga.operators.kinematic_viscosity_operator import Kinematic_Viscosity_Operator

        if flag :
            # Create Operator if necessary
            if self.kv_operator is None:
                self.kv_operator = Kinematic_Viscosity_Operator(self)
        else:
            if self.kv_operator is None:
                return
            else:
                # Remove operator from fractional_step_operators
                self.fractional_step_operators.remove(self.kv_operator)
                self.kv_operator = None
                




    def set_beta(self, beta):
        """Shorthand to assign one constant value [0,2] to all limiters.
        0 Corresponds to first order, where as larger values make use of
        the second order scheme.
        """

        self.beta_w = beta
        self.beta_w_dry = beta
        self.quantities['stage'].beta = beta

        self.beta_uh = beta
        self.beta_uh_dry = beta
        self.quantities['xmomentum'].beta = beta

        self.beta_vh = beta
        self.beta_vh_dry = beta
        self.quantities['ymomentum'].beta = beta


    def set_store_vertices_uniquely(self, flag=True, reduction=None):
        """Decide whether vertex values should be stored uniquely as
        computed in the model (True) or whether they should be reduced to one
        value per vertex using self.reduction (False).
        """

        # FIXME (Ole): how about using the word "continuous vertex values" or 
        # "continuous stage surface"
        self.smooth = not flag

        # Reduction operation for get_vertex_values
        if reduction is None:
            self.reduction = mean
            #self.reduction = min  #Looks better near steep slopes

    def set_store_vertices_smoothly(self, flag=True, reduction=None):
        """Decide whether vertex values should be stored smoothly (one value per vertex)
        or uniquely as
        computed in the model (False)
        """

        # FIXME (Ole): how about using the word "continuous vertex values" or
        # "continuous stage surface"
        self.smooth = flag

        # Reduction operation for get_vertex_values
        if reduction is None:
            self.reduction = mean
            #self.reduction = min  #Looks better near steep slopes

    def set_minimum_storable_height(self, minimum_storable_height):
        """Set the minimum depth that will be written to an SWW file.

        minimum_storable_height  minimum allowed SWW depth is in meters

        This is useful for removing thin water layers that seems to be caused
        by friction creep.
        """

        self.minimum_storable_height = minimum_storable_height

    def set_minimum_allowed_height(self, minimum_allowed_height):
        """Set minimum depth that will be recognised in the numerical scheme.

        minimum_allowed_height  minimum allowed depth in meters

        The parameter H0 (Minimal height for flux computation) is also set by
        this function.
        """

        #FIXME (Ole): rename H0 to minimum_allowed_height_in_flux_computation

        #FIXME (Ole): Maybe use histogram to identify isolated extreme speeds
        #and deal with them adaptively similarly to how we used to use 1 order
        #steps to recover.

        self.minimum_allowed_height = minimum_allowed_height
        self.H0 = minimum_allowed_height

    def set_maximum_allowed_speed(self, maximum_allowed_speed):
        """Set the maximum particle speed that is allowed in water shallower
        than minimum_allowed_height.

        maximum_allowed_speed  

        This is useful for controlling speeds in very thin layers of water and
        at the same time allow some movement avoiding pooling of water.
        """

        self.maximum_allowed_speed = maximum_allowed_speed

    def set_points_file_block_line_size(self, points_file_block_line_size):
        """
        """

        self.points_file_block_line_size = points_file_block_line_size

        
    # FIXME: Probably obsolete in its curren form    
    def set_quantities_to_be_stored(self, q):
        """Specify which quantities will be stored in the SWW file.
        
        q must be either:
          - a dictionary with quantity names
          - a list of quantity names (for backwards compatibility)
          - None

        The format of the dictionary is as follows
        
        quantity_name: flag where flag must be either 1 or 2.
        If flag is 1, the quantity is considered static and will be 
        stored once at the beginning of the simulation in a 1D array.
        
        If flag is 2, the quantity is considered time dependent and 
        it will be stored at each yieldstep by appending it to the 
        appropriate 2D array in the sww file.   
        
        If q is None, storage will be switched off altogether.
        
        Once the simulation has started and thw sww file opened, 
        this function will have no effect.
        
        The format, where q is a list of names is for backwards compatibility
        only.
        It will take the specified quantities to be time dependent and assume 
        'elevation' to be static regardless.
        """

        if q is None:
            self.quantities_to_be_stored = {}
            self.store = False
            return

        # Check correctness
        for quantity_name in q:
            msg = ('Quantity %s is not a valid conserved quantity'
                   % quantity_name)
            assert quantity_name in self.quantities, msg

        assert isinstance(q, dict)
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
        wet_indices = num.compress(depth > minimum_allowed_height,
                                   num.arange(len(depth)))
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


    def get_flow_through_cross_section(self, polyline, verbose=False):
        """Get the total flow through an arbitrary poly line.

        This is a run-time equivalent of the function with same name
        in sww_interrogate.py

        Input:
            polyline: Representation of desired cross section - it may contain
                      multiple sections allowing for complex shapes. Assume
                      absolute UTM coordinates.
                      Format [[x0, y0], [x1, y1], ...]

        Output:
            Q: Total flow [m^3/s] across given segments.
        """


        cross_section = Cross_section(self, polyline, verbose)

        return cross_section.get_flow_through_cross_section()


    def get_energy_through_cross_section(self, polyline,
                                         kind='total',
                                         verbose=False):
        """Obtain average energy head [m] across specified cross section.

        Inputs:
            polyline: Representation of desired cross section - it may contain
                      multiple sections allowing for complex shapes. Assume
                      absolute UTM coordinates.
                      Format [[x0, y0], [x1, y1], ...]
            kind:     Select which energy to compute.
                      Options are 'specific' and 'total' (default)

        Output:
            E: Average energy [m] across given segments for all stored times.

        The average velocity is computed for each triangle intersected by
        the polyline and averaged weighted by segment lengths.

        The typical usage of this function would be to get average energy of
        flow in a channel, and the polyline would then be a cross section
        perpendicular to the flow.

        #FIXME (Ole) - need name for this energy reflecting that its dimension
        is [m].
        """



        cross_section = Cross_section(self, polyline, verbose)

        return cross_section.get_energy_through_cross_section(kind)


    def check_integrity(self):
        """ Run integrity checks on shallow water domain. """
        Generic_Domain.check_integrity(self)

        #Check that we are solving the shallow water wave equation
        msg = 'First conserved quantity must be "stage"'
        assert self.conserved_quantities[0] == 'stage', msg
        msg = 'Second conserved quantity must be "xmomentum"'
        assert self.conserved_quantities[1] == 'xmomentum', msg
        msg = 'Third conserved quantity must be "ymomentum"'
        assert self.conserved_quantities[2] == 'ymomentum', msg

    def extrapolate_second_order_sw(self):
        """Call correct module function
            (either from this module or C-extension)"""
        extrapolate_second_order_sw(self)

    def compute_fluxes(self):
        """Call correct module function
            (either from this module or C-extension)"""
        compute_fluxes(self)


    def distribute_to_vertices_and_edges(self):
        """ Call correct module function """
        if self.use_edge_limiter:
            distribute_using_edge_limiter(self)
        else:
            distribute_using_vertex_limiter(self)

    def protect_against_infinitesimal_and_negative_heights(self):
        """ Clean up the stage and momentum values to ensure non-negative heights
        """

        protect_against_infinitesimal_and_negative_heights(self)
        

    def update_centroids_of_velocities_and_height(self):
        """Calculate the centroid values of velocities and height based
        on the values of the quantities stage and x and y momentum

        Assumes that stage and momentum are up to date

        Useful for kinematic viscosity calculations
        """

        # For shallow water we need to update height xvelocity and yvelocity

        #Shortcuts
        W  = self.quantities['stage']
        UH = self.quantities['xmomentum']
        VH = self.quantities['ymomentum']
        H  = self.quantities['height']
        Z  = self.quantities['elevation']
        U  = self.quantities['xvelocity']
        V  = self.quantities['yvelocity']

        #print num.min(W.centroid_values)

        # Make sure boundary values of conserved quantites
        # are consistent with value of functions at centroids
        self.distribute_to_vertices_and_edges()
        Z.set_boundary_values_from_edges()
        W.set_boundary_values_from_edges()
        UH.set_boundary_values_from_edges()
        VH.set_boundary_values_from_edges()

        # Update height values
        H.set_values(W.centroid_values-Z.centroid_values, location='centroids')
        H.set_boundary_values(W.boundary_values-Z.boundary_values)

        assert num.min(H.centroid_values) >= 0
        assert num.min(H.boundary_values) >= 0

        #Aliases
        uh_C  = UH.centroid_values
        vh_C  = VH.centroid_values
        h_C   = H.centroid_values

        uh_B  = UH.boundary_values
        vh_B  = VH.boundary_values
        h_B   = H.boundary_values

        H0 = 1.0e-8
        
        U.set_values(uh_C/(h_C + H0/h_C), location='centroids')
        V.set_values(vh_C/(h_C + H0/h_C), location='centroids')

        U.set_boundary_values(uh_B/(h_B + H0/h_B))
        V.set_boundary_values(vh_B/(h_B + H0/h_B))
 


    def update_centroids_of_momentum_from_velocity(self):
        """Calculate the centroid value of x and y momentum from height and velocities

        Assumes centroids of height and velocities are up to date

        Useful for kinematic viscosity calculations
        """

        # For shallow water we need to update height xvelocity and yvelocity

        #Shortcuts
        UH = self.quantities['xmomentum']
        VH = self.quantities['ymomentum']
        H  = self.quantities['height']
        Z  = self.quantities['elevation']
        U  = self.quantities['xvelocity']
        V  = self.quantities['yvelocity']


        #Arrays
        u_C  = U.centroid_values
        v_C  = V.centroid_values
        h_C  = H.centroid_values

        u_B  = U.boundary_values
        v_B  = V.boundary_values
        h_B  = H.boundary_values

        UH.set_values(u_C*h_C , location='centroids')
        VH.set_values(v_C*h_C , location='centroids')

        self.distribute_to_vertices_and_edges()


    def evolve(self,
               yieldstep=None,
               finaltime=None,
               duration=None,
               skip_initial_step=False):
        """Specialisation of basic evolve method from parent class.
        
            Evolve the model by 1 step.
        """

        # Call check integrity here rather than from user scripts
        # self.check_integrity()

        msg = 'Attribute self.beta_w must be in the interval [0, 2]'
        assert 0 <= self.beta_w <= 2.0, msg

        # Initial update of vertex and edge values before any STORAGE
        # and or visualisation.
        # This is done again in the initialisation of the Generic_Domain
        # evolve loop but we do it here to ensure the values are ok for storage.
        self.distribute_to_vertices_and_edges()

        if self.store is True and self.get_time() == self.get_starttime():
            self.initialise_storage()

        # Call basic machinery from parent class
        for t in Generic_Domain.evolve(self, yieldstep=yieldstep,
                                       finaltime=finaltime, duration=duration,
                                       skip_initial_step=skip_initial_step):
            # Store model data, e.g. for subsequent visualisation
            if self.store is True:
                self.store_timestep()

            # Pass control on to outer loop for more specific actions
            yield(t)


    def initialise_storage(self):
        """Create and initialise self.writer object for storing data.
        Also, save x,y and bed elevation
        """
        
        # Initialise writer
        self.writer = SWW_file(self)

        # Store vertices and connectivity
        self.writer.store_connectivity()


    def store_timestep(self):
        """Store time dependent quantities and time.

        Precondition:
           self.writer has been initialised
        """

        self.writer.store_timestep()


    def timestepping_statistics(self,
                                track_speeds=False,
                                triangle_id=None):
        """Return string with time stepping statistics for printing or logging

        Optional boolean keyword track_speeds decides whether to report
        location of smallest timestep as well as a histogram and percentile
        report.
        """

        from anuga.config import epsilon, g

        # Call basic machinery from parent class
        msg = Generic_Domain.timestepping_statistics(self, track_speeds,
                                                     triangle_id)

        if track_speeds is True:
            # qwidth determines the text field used for quantities
            qwidth = self.qwidth

            # Selected triangle
            k = self.k

            # Report some derived quantities at vertices, edges and centroid
            # specific to the shallow water wave equation
            z = self.quantities['elevation']
            w = self.quantities['stage']

            Vw = w.get_values(location='vertices', indices=[k])[0]
            Ew = w.get_values(location='edges', indices=[k])[0]
            Cw = w.get_values(location='centroids', indices=[k])

            Vz = z.get_values(location='vertices', indices=[k])[0]
            Ez = z.get_values(location='edges', indices=[k])[0]
            Cz = z.get_values(location='centroids', indices=[k])

            name = 'depth'
            Vh = Vw-Vz
            Eh = Ew-Ez
            Ch = Cw-Cz

            message  = '    %s: vertex_values =  %.4f,\t %.4f,\t %.4f\n'\
                 % (name.ljust(qwidth), Vh[0], Vh[1], Vh[2])

            message += '    %s: edge_values =    %.4f,\t %.4f,\t %.4f\n'\
                 % (name.ljust(qwidth), Eh[0], Eh[1], Eh[2])

            message += '    %s: centroid_value = %.4f\n'\
                 % (name.ljust(qwidth), Ch[0])

            msg += message

            uh = self.quantities['xmomentum']
            vh = self.quantities['ymomentum']

            Vuh = uh.get_values(location='vertices', indices=[k])[0]
            Euh = uh.get_values(location='edges', indices=[k])[0]
            Cuh = uh.get_values(location='centroids', indices=[k])

            Vvh = vh.get_values(location='vertices', indices=[k])[0]
            Evh = vh.get_values(location='edges', indices=[k])[0]
            Cvh = vh.get_values(location='centroids', indices=[k])

            # Speeds in each direction
            Vu = Vuh/(Vh + epsilon)
            Eu = Euh/(Eh + epsilon)
            Cu = Cuh/(Ch + epsilon)
            name = 'U'
            message  = '    %s: vertex_values =  %.4f,\t %.4f,\t %.4f\n' \
                 % (name.ljust(qwidth), Vu[0], Vu[1], Vu[2])

            message += '    %s: edge_values =    %.4f,\t %.4f,\t %.4f\n' \
                 % (name.ljust(qwidth), Eu[0], Eu[1], Eu[2])

            message += '    %s: centroid_value = %.4f\n' \
                 % (name.ljust(qwidth), Cu[0])

            msg += message

            Vv = Vvh/(Vh + epsilon)
            Ev = Evh/(Eh + epsilon)
            Cv = Cvh/(Ch + epsilon)
            name = 'V'
            message  = '    %s: vertex_values =  %.4f,\t %.4f,\t %.4f\n' \
                 % (name.ljust(qwidth), Vv[0], Vv[1], Vv[2])

            message += '    %s: edge_values =    %.4f,\t %.4f,\t %.4f\n' \
                 % (name.ljust(qwidth), Ev[0], Ev[1], Ev[2])

            message += '    %s: centroid_value = %.4f\n'\
                 %(name.ljust(qwidth), Cv[0])

            msg += message

            # Froude number in each direction
            name = 'Froude (x)'
            Vfx = Vu/(num.sqrt(g*Vh) + epsilon)
            Efx = Eu/(num.sqrt(g*Eh) + epsilon)
            Cfx = Cu/(num.sqrt(g*Ch) + epsilon)

            message  = '    %s: vertex_values =  %.4f,\t %.4f,\t %.4f\n'\
                 % (name.ljust(qwidth), Vfx[0], Vfx[1], Vfx[2])

            message += '    %s: edge_values =    %.4f,\t %.4f,\t %.4f\n'\
                 % (name.ljust(qwidth), Efx[0], Efx[1], Efx[2])

            message += '    %s: centroid_value = %.4f\n'\
                 % (name.ljust(qwidth), Cfx[0])

            msg += message

            name = 'Froude (y)'
            Vfy = Vv/(num.sqrt(g*Vh) + epsilon)
            Efy = Ev/(num.sqrt(g*Eh) + epsilon)
            Cfy = Cv/(num.sqrt(g*Ch) + epsilon)

            message  = '    %s: vertex_values =  %.4f,\t %.4f,\t %.4f\n'\
                 % (name.ljust(qwidth), Vfy[0], Vfy[1], Vfy[2])

            message += '    %s: edge_values =    %.4f,\t %.4f,\t %.4f\n'\
                 % (name.ljust(qwidth), Efy[0], Efy[1], Efy[2])

            message += '    %s: centroid_value = %.4f\n'\
                 % (name.ljust(qwidth), Cfy[0])

            msg += message

        return msg
       
    def print_timestepping_statistics(self, *args, **kwargs):

        msg = self.timestepping_statistics(*args, **kwargs)

        print msg


    def compute_boundary_flows(self):
        """Compute boundary flows at current timestep.
			
        Quantities computed are:
           Total inflow across boundary
           Total outflow across boundary
           Flow across each tagged boundary segment
        """
		
        # Run through boundary array and compute for each segment
        # the normal momentum ((uh, vh) dot normal) times segment length.
        # Based on sign accumulate this into boundary_inflow and
        # boundary_outflow.
			
        # Compute flows along boundary
        
        uh = self.get_quantity('xmomentum').get_values(location='edges')
        vh = self.get_quantity('ymomentum').get_values(location='edges')        
        
        # Loop through edges that lie on the boundary and calculate 
        # flows
        boundary_flows = {}
        total_boundary_inflow = 0.0
        total_boundary_outflow = 0.0
        for vol_id, edge_id in self.boundary:
            # Compute normal flow across edge. Since normal vector points
            # away from triangle, a positive sign means that water 
            # flows *out* from this triangle.
             
            momentum = [uh[vol_id, edge_id], vh[vol_id, edge_id]]
            normal = self.mesh.get_normal(vol_id, edge_id)
            length = self.mesh.get_edgelength(vol_id, edge_id)            
            normal_flow = num.dot(momentum, normal)*length
            
            # Reverse sign so that + is taken to mean inflow
            # and - means outflow. This is more intuitive.
            edge_flow = -normal_flow
            
            # Tally up inflows and outflows separately
            if edge_flow > 0:
                # Flow is inflow      
                total_boundary_inflow += edge_flow       
            else:
                # Flow is outflow
                total_boundary_outflow += edge_flow    

            # Tally up flows by boundary tag
            tag = self.boundary[(vol_id, edge_id)]
            
            if tag not in boundary_flows:
                boundary_flows[tag] = 0.0
            boundary_flows[tag] += edge_flow
            
                
        return boundary_flows, total_boundary_inflow, total_boundary_outflow
        

    def compute_forcing_flows(self):
        """Compute flows in and out of domain due to forcing terms.
			
        Quantities computed are:
		
        
           Total inflow through forcing terms
           Total outflow through forcing terms
           Current total volume in domain        

        """

        #FIXME(Ole): We need to separate what part of explicit_update was 
        # due to the normal flux calculations and what is due to forcing terms.
        
        pass
			
        
    def compute_total_volume(self):
        """Compute total volume (m^3) of water in entire domain
        """
        
        area = self.mesh.get_areas()
        
        stage = self.get_quantity('stage').get_values(location='centroids')
        elevation = \
            self.get_quantity('elevation').get_values(location='centroids')
        depth = stage-elevation
        
        return num.sum(depth*area)
        
        
    def volumetric_balance_statistics(self):                
        """Create volumetric balance report suitable for printing or logging.
        """
        
        (boundary_flows, total_boundary_inflow,
         total_boundary_outflow) = self.compute_boundary_flows() 
        
        message = '---------------------------\n'        
        message += 'Volumetric balance report:\n'
        message += '--------------------------\n'
        message += 'Total boundary inflow [m^3/s]: %.2f\n' % total_boundary_inflow
        message += 'Total boundary outflow [m^3/s]: %.2f\n' % total_boundary_outflow
        message += 'Net boundary flow by tags [m^3/s]\n'
        for tag in boundary_flows:
            message += '    %s [m^3/s]: %.2f\n' % (tag, boundary_flows[tag])
        
        message += 'Total net boundary flow [m^3/s]: %.2f\n' % \
                    (total_boundary_inflow + total_boundary_outflow) 
        message += 'Total volume in domain [m^3]: %.2f\n' % \
                    self.compute_total_volume()
        
        # The go through explicit forcing update and record the rate of change
        # for stage and 
        # record into forcing_inflow and forcing_outflow. Finally compute  
        # integral of depth to obtain total volume of domain.
	
        # FIXME(Ole): This part is not yet done.		
        
        return message        
           
################################################################################
# End of class Shallow Water Domain
################################################################################

#-----------------
# Flux computation
#-----------------

def compute_fluxes(domain):
    """Compute fluxes and timestep suitable for all volumes in domain.

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
    from shallow_water_ext import compute_fluxes_ext_central \
                                  as compute_fluxes_ext

    # Shortcuts
    Stage = domain.quantities['stage']
    Xmom = domain.quantities['xmomentum']
    Ymom = domain.quantities['ymomentum']
    Bed = domain.quantities['elevation']

    timestep = float(sys.maxint)

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

################################################################################
# Module functions for gradient limiting
################################################################################

def extrapolate_second_order_sw(domain):
    """Wrapper calling C version of extrapolate_second_order_sw.

    domain  the domain to operate on

    Note MH090605: The following method belongs to the shallow_water domain
    class, see comments in the corresponding method in shallow_water_ext.c
    """

    from shallow_water_ext import extrapolate_second_order_sw as extrapol2

    # Shortcuts
    Stage = domain.quantities['stage']
    Xmom = domain.quantities['xmomentum']
    Ymom = domain.quantities['ymomentum']
    Elevation = domain.quantities['elevation']

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

def distribute_using_vertex_limiter(domain):
    """Distribution from centroids to vertices specific to the SWW equation.

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

    # Remove very thin layers of water
    protect_against_infinitesimal_and_negative_heights(domain)

    # Extrapolate all conserved quantities
    if domain.optimised_gradient_limiter:
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
            raise Exception('Unknown order')
    else:
        # Old code:
        for name in domain.conserved_quantities:
            Q = domain.quantities[name]

            if domain._order_ == 1:
                Q.extrapolate_first_order()
            elif domain._order_ == 2:
                Q.extrapolate_second_order_and_limit_by_vertex()
            else:
                raise Exception('Unknown order')

    # Take bed elevation into account when water heights are small
    balance_deep_and_shallow(domain)

    # Compute edge values by interpolation
    for name in domain.conserved_quantities:
        Q = domain.quantities[name]
        Q.interpolate_from_vertices_to_edges()

def distribute_using_edge_limiter(domain):
    """Distribution from centroids to edges specific to the SWW eqn.

    It will ensure that h (w-z) is always non-negative even in the
    presence of steep bed-slopes by taking a weighted average between shallow
    and deep cases.

    In addition, all conserved quantities get distributed as per either a
    constant (order==1) or a piecewise linear function (order==2).


    Precondition:
      All quantities defined at centroids and bed elevation defined at
      vertices.

    Postcondition
      Conserved quantities defined at vertices
    """

    # Remove very thin layers of water
    protect_against_infinitesimal_and_negative_heights(domain)

    for name in domain.conserved_quantities:
        Q = domain.quantities[name]
        if domain._order_ == 1:
            Q.extrapolate_first_order()
        elif domain._order_ == 2:
            Q.extrapolate_second_order_and_limit_by_edge()
        else:
            raise Exception('Unknown order')

    balance_deep_and_shallow(domain)

    # Compute edge values by interpolation
    for name in domain.conserved_quantities:
        Q = domain.quantities[name]
        Q.interpolate_from_vertices_to_edges()

def protect_against_infinitesimal_and_negative_heights(domain):
    """Protect against infinitesimal heights and associated high velocities"""

    from shallow_water_ext import protect

    # Shortcuts
    wc = domain.quantities['stage'].centroid_values
    zc = domain.quantities['elevation'].centroid_values
    xmomc = domain.quantities['xmomentum'].centroid_values
    ymomc = domain.quantities['ymomentum'].centroid_values

    protect(domain.minimum_allowed_height, domain.maximum_allowed_speed,
            domain.epsilon, wc, zc, xmomc, ymomc)

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

    from shallow_water_ext import balance_deep_and_shallow \
                                  as balance_deep_and_shallow_c

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

    balance_deep_and_shallow_c(domain,
                               wc, zc, wv, zv, wc,
                               xmomc, ymomc, xmomv, ymomv)



################################################################################
# Standard forcing terms
################################################################################

def gravity(domain):
    """Apply gravitational pull in the presence of bed slope
    Wrapper calls underlying C implementation
    """

    from shallow_water_ext import gravity as gravity_c

    xmom_update = domain.quantities['xmomentum'].explicit_update
    ymom_update = domain.quantities['ymomentum'].explicit_update

    stage = domain.quantities['stage']
    elevation = domain.quantities['elevation']

    #FIXME SR Should avoid allocating memory!
    height = stage.centroid_values - elevation.centroid_values
    elevation = elevation.vertex_values

    point = domain.get_vertex_coordinates()

    gravity_c(domain.g, height, elevation, point, xmom_update, ymom_update)

def manning_friction_implicit(domain):
    """Apply (Manning) friction to water momentum
    Wrapper for c version
    """

    from shallow_water_ext import manning_friction_flat
    from shallow_water_ext import manning_friction_sloped

    xmom = domain.quantities['xmomentum']
    ymom = domain.quantities['ymomentum']

    x = domain.get_vertex_coordinates()
    
    w = domain.quantities['stage'].centroid_values
    z = domain.quantities['elevation'].vertex_values

    uh = xmom.centroid_values
    vh = ymom.centroid_values
    eta = domain.quantities['friction'].centroid_values

    xmom_update = xmom.semi_implicit_update
    ymom_update = ymom.semi_implicit_update

    eps = domain.minimum_allowed_height
    g = domain.g

    if domain.use_sloped_mannings:
        manning_friction_sloped(g, eps, x, w, uh, vh, z, eta, xmom_update, \
                                ymom_update)
    else:
        manning_friction_flat(g, eps, w, uh, vh, z, eta, xmom_update, \
                                ymom_update)
    

def manning_friction_explicit(domain):
    """Apply (Manning) friction to water momentum
    Wrapper for c version
    """

    from shallow_water_ext import manning_friction_flat
    from shallow_water_ext import manning_friction_sloped

    xmom = domain.quantities['xmomentum']
    ymom = domain.quantities['ymomentum']

    x = domain.get_vertex_coordinates()
    
    w = domain.quantities['stage'].centroid_values
    z = domain.quantities['elevation'].vertex_values

    uh = xmom.centroid_values
    vh = ymom.centroid_values
    eta = domain.quantities['friction'].centroid_values

    xmom_update = xmom.explicit_update
    ymom_update = ymom.explicit_update

    eps = domain.minimum_allowed_height

    if domain.use_sloped_mannings:
        manning_friction_sloped(domain.g, eps, x, w, uh, vh, z, eta, xmom_update, \
                            ymom_update)
    else:
        manning_friction_flat(domain.g, eps, w, uh, vh, z, eta, xmom_update, \
                            ymom_update)



# FIXME (Ole): This was implemented for use with one of the analytical solutions
def linear_friction(domain):
    """Apply linear friction to water momentum

    Assumes quantity: 'linear_friction' to be present
    """

    w = domain.quantities['stage'].centroid_values
    z = domain.quantities['elevation'].centroid_values
    h = w-z

    uh = domain.quantities['xmomentum'].centroid_values
    vh = domain.quantities['ymomentum'].centroid_values
    tau = domain.quantities['linear_friction'].centroid_values

    xmom_update = domain.quantities['xmomentum'].semi_implicit_update
    ymom_update = domain.quantities['ymomentum'].semi_implicit_update

    num_tris = len(domain)
    eps = domain.minimum_allowed_height

    for k in range(num_tris):
        if tau[k] >= eps:
            if h[k] >= eps:
                S = -tau[k]/h[k]

                #Update momentum
                xmom_update[k] += S*uh[k]
                ymom_update[k] += S*vh[k]

def depth_dependent_friction(domain, default_friction,
                             surface_roughness_data,
                             verbose=False):
    """Returns an array of friction values for each wet element adjusted for
            depth.

    Inputs:
        domain - computational domain object
        default_friction - depth independent bottom friction
        surface_roughness_data - N x 5 array of n0, d1, n1, d2, n2 values 
        for each friction region.

    Outputs:
        wet_friction - Array that can be used directly to update friction as
                        follows:
                       domain.set_quantity('friction', wet_friction)

        
        
    """
    
    default_n0 = 0  # James - this was missing, don't know what it should be
    
    # Create a temp array to store updated depth dependent
    # friction for wet elements
    # EHR this is outwardly inneficient but not obvious how to avoid
    # recreating each call??????

    wet_friction    = num.zeros(len(domain), num.float)
    wet_friction[:] = default_n0  # Initially assign default_n0 to all array so
                                  # sure have no zeros values
    
    # create depth instance for this timestep    
    depth = domain.create_quantity_from_expression('stage - elevation')  
    # Recompute depth as vector  
    d_vals = depth.get_values(location='centroids')
 
    # rebuild the 'friction' values adjusted for depth at this instant
    # loop for each wet element in domain
   
    for i in domain.get_wet_elements():        
        # Get roughness data for each element
        d1 = float(surface_roughness_data[i, 1])
        n1 = float(surface_roughness_data[i, 2])
        d2 = float(surface_roughness_data[i, 3])
        n2 = float(surface_roughness_data[i, 4])
        
        
        # Recompute friction values from depth for this element 
               
        if d_vals[i]   <= d1: 
            ddf = n1
        elif d_vals[i] >= d2:
            ddf = n2
        else:
            ddf = n1 + ((n2-n1)/(d2-d1))*(d_vals[i]-d1)
            
        # check sanity of result
        if (ddf  < 0.010 or \
                            ddf > 9999.0) :
            log.critical('>>>> WARNING: computed depth_dependent friction '
                         'out of range, ddf%f, n1=%f, n2=%f'
                         % (ddf, n1, n2))
        
        # update depth dependent friction  for that wet element
        wet_friction[i] = ddf
        
    # EHR add code to show range of 'friction across domain at this instant as
    # sanity check?????????
    
    if verbose :
        # return array of domain nvals
        nvals = domain.get_quantity('friction').get_values(location='centroids')
        n_min = min(nvals)
        n_max = max(nvals)
        
        log.critical('         ++++ calculate_depth_dependent_friction - '
                     'Updated friction - range  %7.3f to %7.3f'
                     % (n_min, n_max))
    
    return wet_friction



################################################################################
# Initialise module
################################################################################

def _raise_compile_exception():
    """ Raise exception if compiler not available. """
    msg = 'C implementations could not be accessed by %s.\n ' % __file__
    msg += 'Make sure compile_all.py has been run as described in '
    msg += 'the ANUGA installation guide.'
    raise Exception(msg)

from anuga.utilities import compile
if not compile.can_use_C_extension('shallow_water_ext.c'):
    _raise_compile_exception()

if __name__ == "__main__":
    pass

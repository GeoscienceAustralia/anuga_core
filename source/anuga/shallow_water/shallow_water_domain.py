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


import numpy as num

from anuga.abstract_2d_finite_volumes.neighbour_mesh import segment_midpoints
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
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions\
     import AWI_boundary

from anuga.pmesh.mesh_interface import create_mesh_from_regions
from anuga.utilities.numerical_tools import gradient, mean, ensure_numeric
from anuga.geospatial_data.geospatial_data import ensure_geospatial

from anuga.config import minimum_storable_height
from anuga.config import minimum_allowed_height, maximum_allowed_speed
from anuga.config import g, epsilon, beta_w, beta_w_dry,\
     beta_uh, beta_uh_dry, beta_vh, beta_vh_dry, tight_slope_limiters
from anuga.config import alpha_balance
from anuga.config import optimise_dry_cells
from anuga.config import optimised_gradient_limiter
from anuga.config import use_edge_limiter
from anuga.config import use_centroid_velocities
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a

from anuga.fit_interpolate.interpolate import Modeltime_too_late, \
                                              Modeltime_too_early

from anuga.utilities.polygon import inside_polygon, polygon_area, \
                                    is_inside_polygon
import anuga.utilities.log as log

import types
from types import IntType, FloatType
from warnings import warn


################################################################################
# Shallow water domain
################################################################################

##
# @brief Class for a shallow water domain.
class Domain(Generic_Domain):

    #conserved_quantities = ['stage', 'xmomentum', 'ymomentum']
    #other_quantities = ['elevation', 'friction']

    ##
    # @brief Instantiate a shallow water domain.
    # @param coordinates
    # @param vertices
    # @param boundary
    # @param tagged_elements
    # @param geo_reference
    # @param use_inscribed_circle
    # @param mesh_filename
    # @param use_cache
    # @param verbose
    # @param full_send_dict
    # @param ghost_recv_dict
    # @param processor
    # @param numproc
    # @param number_of_full_nodes
    # @param number_of_full_triangles
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

        # Define quantities for the shallow_water domain         
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

        self.forcing_terms.append(manning_friction_implicit)
        self.forcing_terms.append(gravity)

        # Stored output
        self.store = True
        self.set_store_vertices_uniquely(False)
        self.minimum_storable_height = minimum_storable_height
        self.quantities_to_be_stored = {'elevation': 1, 
                                        'stage': 2, 
                                        'xmomentum': 2, 
                                        'ymomentum': 2}

        # Limiters
        self.use_edge_limiter = use_edge_limiter
        self.optimised_gradient_limiter = optimised_gradient_limiter
        self.use_centroid_velocities = use_centroid_velocities

    ##
    # @brief
    # @param beta
    def set_all_limiters(self, beta):
        """Shorthand to assign one constant value [0,1] to all limiters.
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

    ##
    # @brief
    # @param flag
    # @param reduction
    def set_store_vertices_uniquely(self, flag, reduction=None):
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

    ##
    # @brief Set the minimum depth that will be written to an SWW file.
    # @param minimum_storable_height The minimum stored height (in m).
    def set_minimum_storable_height(self, minimum_storable_height):
        """Set the minimum depth that will be recognised when writing
        to an sww file. This is useful for removing thin water layers
        that seems to be caused by friction creep.

        The minimum allowed sww depth is in meters.
        """

        self.minimum_storable_height = minimum_storable_height

    ##
    # @brief
    # @param minimum_allowed_height
    def set_minimum_allowed_height(self, minimum_allowed_height):
        """Set minimum depth that will be recognised in the numerical scheme.

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

    ##
    # @brief
    # @param maximum_allowed_speed
    def set_maximum_allowed_speed(self, maximum_allowed_speed):
        """Set the maximum particle speed that is allowed in water
        shallower than minimum_allowed_height. This is useful for
        controlling speeds in very thin layers of water and at the same time
        allow some movement avoiding pooling of water.
        """

        self.maximum_allowed_speed = maximum_allowed_speed

    ##
    # @brief
    # @param points_file_block_line_size
    def set_points_file_block_line_size(self, points_file_block_line_size):
        """Set the minimum depth that will be recognised when writing
        to an sww file. This is useful for removing thin water layers
        that seems to be caused by friction creep.

        The minimum allowed sww depth is in meters.
        """
        self.points_file_block_line_size = points_file_block_line_size

        
    # FIXME: Probably obsolete in its curren form    
    ##
    # @brief Set the quantities that will be written to an SWW file.
    # @param q The quantities to be written.
    # @note Param 'q' may be None, single quantity or list of quantity strings.
    # @note If 'q' is None, no quantities will be stored in the SWW file.
    def set_quantities_to_be_stored(self, q):
        """Specify which quantities will be stored in the sww file
        
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
        
        The format, where q is a list of names is for backwards compatibility only.
        It will take the specified quantities to be time dependent and assume 
        'elevation' to be static regardless.
        """

        if q is None:
            self.quantities_to_be_stored = {}
            self.store = False
            return

        # Check correcness
        for quantity_name in q:
            msg = ('Quantity %s is not a valid conserved quantity'
                   % quantity_name)
            assert quantity_name in self.quantities, msg

        if type(q) == types.ListType:

            msg = 'List arguments to set_quantities_to_be_stored '
            msg += 'has been deprecated and will be removed in future '
            msg += 'versions of ANUGA.'
            msg += 'Please use dictionary argument instead'
            warn(msg, DeprecationWarning) 

        
        
            # FIXME: Raise deprecation warning
            tmp = {}
            for x in q:
                tmp[x] = 2
            tmp['elevation'] = 1    
            q = tmp     
            
        assert type(q) == types.DictType    
        self.quantities_to_be_stored = q

    ##
    # @brief
    # @param indices
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

    ##
    # @brief
    # @param indices
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

    ##
    # @brief
    # @param indices
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



    ##
    # @brief Get the total flow through an arbitrary poly line.
    # @param polyline Representation of desired cross section.
    # @param verbose True if this method is to be verbose.
    # @note 'polyline' may contain multiple sections allowing complex shapes.
    # @note Assume absolute UTM coordinates.
    def get_flow_through_cross_section(self, polyline, verbose=False):
        """Get the total flow through an arbitrary poly line.

        This is a run-time equivalent of the function with same name
        in data_manager.py

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





    ##
    # @brief Get the total flow through an arbitrary poly line.
    # @param polyline Representation of desired cross section.
    # @param verbose True if this method is to be verbose.
    # @note 'polyline' may contain multiple sections allowing complex shapes.
    # @note Assume absolute UTM coordinates.
    def old_get_flow_through_cross_section(self, polyline, verbose=False):
        """Get the total flow through an arbitrary poly line.

        This is a run-time equivalent of the function with same name
        in data_manager.py

        Input:
            polyline: Representation of desired cross section - it may contain
                      multiple sections allowing for complex shapes. Assume
                      absolute UTM coordinates.
                      Format [[x0, y0], [x1, y1], ...]

        Output:
            Q: Total flow [m^3/s] across given segments.
        """

        # Find all intersections and associated triangles.
        segments = self.get_intersecting_segments(polyline, use_cache=True,
                                                  verbose=verbose)

        # Get midpoints
        midpoints = segment_midpoints(segments)

        # Make midpoints Geospatial instances
        midpoints = ensure_geospatial(midpoints, self.geo_reference)

        # Compute flow
        if verbose:
            log.critical('Computing flow through specified cross section')

        # Get interpolated values
        xmomentum = self.get_quantity('xmomentum')
        ymomentum = self.get_quantity('ymomentum')

        uh = xmomentum.get_values(interpolation_points=midpoints,
                                  use_cache=True)
        vh = ymomentum.get_values(interpolation_points=midpoints,
                                  use_cache=True)

        # Compute and sum flows across each segment
        total_flow = 0
        for i in range(len(uh)):
            # Inner product of momentum vector with segment normal [m^2/s]
            normal = segments[i].normal
            normal_momentum = uh[i]*normal[0] + vh[i]*normal[1]

            # Flow across this segment [m^3/s]
            segment_flow = normal_momentum*segments[i].length

            # Accumulate
            total_flow += segment_flow

        return total_flow

    ##
    # @brief 
    # @param polyline Representation of desired cross section.
    # @param kind Select energy type to compute ('specific' or 'total').
    # @param verbose True if this method is to be verbose.
    # @note 'polyline' may contain multiple sections allowing complex shapes.
    # @note Assume absolute UTM coordinates.
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

        from anuga.config import g, epsilon, velocity_protection as h0

        # Find all intersections and associated triangles.
        segments = self.get_intersecting_segments(polyline, use_cache=True,
                                                  verbose=verbose)

        # Get midpoints
        midpoints = segment_midpoints(segments)

        # Make midpoints Geospatial instances
        midpoints = ensure_geospatial(midpoints, self.geo_reference)

        # Compute energy
        if verbose: log.critical('Computing %s energy' % kind)

        # Get interpolated values
        stage = self.get_quantity('stage')
        elevation = self.get_quantity('elevation')
        xmomentum = self.get_quantity('xmomentum')
        ymomentum = self.get_quantity('ymomentum')

        w = stage.get_values(interpolation_points=midpoints, use_cache=True)
        z = elevation.get_values(interpolation_points=midpoints, use_cache=True)
        uh = xmomentum.get_values(interpolation_points=midpoints,
                                  use_cache=True)
        vh = ymomentum.get_values(interpolation_points=midpoints,
                                  use_cache=True)
        h = w-z                # Depth

        # Compute total length of polyline for use with weighted averages
        total_line_length = 0.0
        for segment in segments:
            total_line_length += segment.length

        # Compute and sum flows across each segment
        average_energy = 0.0
        for i in range(len(w)):
            # Average velocity across this segment
            if h[i] > epsilon:
                # Use protection against degenerate velocities
                u = uh[i]/(h[i] + h0/h[i])
                v = vh[i]/(h[i] + h0/h[i])
            else:
                u = v = 0.0

            speed_squared = u*u + v*v
            kinetic_energy = 0.5*speed_squared/g

            if kind == 'specific':
                segment_energy = h[i] + kinetic_energy
            elif kind == 'total':
                segment_energy = w[i] + kinetic_energy
            else:
                msg = 'Energy kind must be either "specific" or "total".'
                msg += ' I got %s' %kind

            # Add to weighted average
            weigth = segments[i].length/total_line_length
            average_energy += segment_energy*weigth

        return average_energy

    ##
    # @brief Run integrity checks on shallow water domain.
    def check_integrity(self):
        Generic_Domain.check_integrity(self)

        #Check that we are solving the shallow water wave equation
        msg = 'First conserved quantity must be "stage"'
        assert self.conserved_quantities[0] == 'stage', msg
        msg = 'Second conserved quantity must be "xmomentum"'
        assert self.conserved_quantities[1] == 'xmomentum', msg
        msg = 'Third conserved quantity must be "ymomentum"'
        assert self.conserved_quantities[2] == 'ymomentum', msg

    ##
    # @brief 
    def extrapolate_second_order_sw(self):
        #Call correct module function (either from this module or C-extension)
        extrapolate_second_order_sw(self)

    ##
    # @brief 
    def compute_fluxes(self):
        #Call correct module function (either from this module or C-extension)
        compute_fluxes(self)

    ##
    # @brief 
    def distribute_to_vertices_and_edges(self):
        # Call correct module function
        if self.use_edge_limiter:
            distribute_using_edge_limiter(self)
        else:
            distribute_using_vertex_limiter(self)

    ##
    # @brief Evolve the model by one step.
    # @param yieldstep 
    # @param finaltime 
    # @param duration 
    # @param skip_initial_step 
    def evolve(self,
               yieldstep=None,
               finaltime=None,
               duration=None,
               skip_initial_step=False):
        """Specialisation of basic evolve method from parent class"""

        # Call check integrity here rather than from user scripts
        # self.check_integrity()

        msg = 'Attribute self.beta_w must be in the interval [0, 2]'
        assert 0 <= self.beta_w <= 2.0, msg

        # Initial update of vertex and edge values before any STORAGE
        # and or visualisation.
        # This is done again in the initialisation of the Generic_Domain
        # evolve loop but we do it here to ensure the values are ok for storage.
        self.distribute_to_vertices_and_edges()

        if self.store is True and self.time == 0.0:
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

    ##
    # @brief 
    def initialise_storage(self):
        """Create and initialise self.writer object for storing data.
        Also, save x,y and bed elevation
        """

        from anuga.shallow_water.data_manager import SWW_file
        
        # Initialise writer
        self.writer = SWW_file(self)

        # Store vertices and connectivity
        self.writer.store_connectivity()

    ##
    # @brief 
    # @param name 
    def store_timestep(self):
        """Store time dependent quantities and time.

        Precondition:
           self.writer has been initialised
        """

        self.writer.store_timestep()

    ##
    # @brief Get time stepping statistics string for printing.
    # @param track_speeds 
    # @param triangle_id 
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

            s  = '    %s: vertex_values =  %.4f,\t %.4f,\t %.4f\n'\
                 %(name.ljust(qwidth), Vh[0], Vh[1], Vh[2])

            s += '    %s: edge_values =    %.4f,\t %.4f,\t %.4f\n'\
                 %(name.ljust(qwidth), Eh[0], Eh[1], Eh[2])

            s += '    %s: centroid_value = %.4f\n'\
                 %(name.ljust(qwidth), Ch[0])

            msg += s

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
            s  = '    %s: vertex_values =  %.4f,\t %.4f,\t %.4f\n'\
                 %(name.ljust(qwidth), Vu[0], Vu[1], Vu[2])

            s += '    %s: edge_values =    %.4f,\t %.4f,\t %.4f\n'\
                 %(name.ljust(qwidth), Eu[0], Eu[1], Eu[2])

            s += '    %s: centroid_value = %.4f\n'\
                 %(name.ljust(qwidth), Cu[0])

            msg += s

            Vv = Vvh/(Vh + epsilon)
            Ev = Evh/(Eh + epsilon)
            Cv = Cvh/(Ch + epsilon)
            name = 'V'
            s  = '    %s: vertex_values =  %.4f,\t %.4f,\t %.4f\n'\
                 %(name.ljust(qwidth), Vv[0], Vv[1], Vv[2])

            s += '    %s: edge_values =    %.4f,\t %.4f,\t %.4f\n'\
                 %(name.ljust(qwidth), Ev[0], Ev[1], Ev[2])

            s += '    %s: centroid_value = %.4f\n'\
                 %(name.ljust(qwidth), Cv[0])

            msg += s

            # Froude number in each direction
            name = 'Froude (x)'
            Vfx = Vu/(num.sqrt(g*Vh) + epsilon)
            Efx = Eu/(num.sqrt(g*Eh) + epsilon)
            Cfx = Cu/(num.sqrt(g*Ch) + epsilon)

            s  = '    %s: vertex_values =  %.4f,\t %.4f,\t %.4f\n'\
                 %(name.ljust(qwidth), Vfx[0], Vfx[1], Vfx[2])

            s += '    %s: edge_values =    %.4f,\t %.4f,\t %.4f\n'\
                 %(name.ljust(qwidth), Efx[0], Efx[1], Efx[2])

            s += '    %s: centroid_value = %.4f\n'\
                 %(name.ljust(qwidth), Cfx[0])

            msg += s

            name = 'Froude (y)'
            Vfy = Vv/(num.sqrt(g*Vh) + epsilon)
            Efy = Ev/(num.sqrt(g*Eh) + epsilon)
            Cfy = Cv/(num.sqrt(g*Ch) + epsilon)

            s  = '    %s: vertex_values =  %.4f,\t %.4f,\t %.4f\n'\
                 %(name.ljust(qwidth), Vfy[0], Vfy[1], Vfy[2])

            s += '    %s: edge_values =    %.4f,\t %.4f,\t %.4f\n'\
                 %(name.ljust(qwidth), Efy[0], Efy[1], Efy[2])

            s += '    %s: centroid_value = %.4f\n'\
                 %(name.ljust(qwidth), Cfy[0])

            msg += s

        return msg
       
        

    def compute_boundary_flows(self):
        """Compute boundary flows at current timestep.
			
        Quantities computed are:
           Total inflow across boundary
           Total outflow across boundary
           Flow across each tagged boundary segment
        """
		
        # Run through boundary array and compute for each segment
        # the normal momentum ((uh, vh) dot normal) times segment length.
        # Based on sign accumulate this into boundary_inflow and boundary_outflow.
			
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
        volume = 0.0
        
        stage = self.get_quantity('stage').get_values(location='centroids')
        elevation = self.get_quantity('elevation').get_values(location='centroids')        
        depth = stage-elevation
        
        return num.sum(depth*area)
        
        
    def volumetric_balance_statistics(self):                
        """Create volumetric balance report suitable for printing or logging.
        """
        
        (boundary_flows, total_boundary_inflow,
         total_boundary_outflow) = self.compute_boundary_flows() 
        
        s = '---------------------------\n'        
        s += 'Volumetric balance report:\n'
        s += '--------------------------\n'
        s += 'Total boundary inflow [m^3/s]: %.2f\n' % total_boundary_inflow
        s += 'Total boundary outflow [m^3/s]: %.2f\n' % total_boundary_outflow        
        s += 'Net boundary flow by tags [m^3/s]\n'
        for tag in boundary_flows:
            s += '    %s [m^3/s]: %.2f\n' % (tag, boundary_flows[tag])
        
        s += 'Total net boundary flow [m^3/s]: %.2f\n' % (total_boundary_inflow + total_boundary_outflow) 
        s += 'Total volume in domain [m^3]: %.2f\n' % self.compute_total_volume()
        
        # The go through explicit forcing update and record the rate of change for stage and 
        # record into forcing_inflow and forcing_outflow. Finally compute integral 
        # of depth to obtain total volume of domain.
	
        # FIXME(Ole): This part is not yet done.		
        
        return s        
           
################################################################################
# End of class Shallow Water Domain
################################################################################

#-----------------
# Flux computation
#-----------------

## @brief Compute fluxes and timestep suitable for all volumes in domain.
# @param domain The domain to calculate fluxes for.
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

    N = len(domain)    # number_of_triangles

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

##
# @brief Wrapper for C version of extrapolate_second_order_sw.
# @param domain The domain to operate on.
# @note MH090605 The following method belongs to the shallow_water domain class
#       see comments in the corresponding method in shallow_water_ext.c
def extrapolate_second_order_sw(domain):
    """Wrapper calling C version of extrapolate_second_order_sw"""

    import sys
    from shallow_water_ext import extrapolate_second_order_sw as extrapol2

    N = len(domain) # number_of_triangles

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

##
# @brief Distribution from centroids to vertices specific to the SWW eqn.
# @param domain The domain to operate on.
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
            raise 'Unknown order'
    else:
        # Old code:
        for name in domain.conserved_quantities:
            Q = domain.quantities[name]

            if domain._order_ == 1:
                Q.extrapolate_first_order()
            elif domain._order_ == 2:
                Q.extrapolate_second_order_and_limit_by_vertex()
            else:
                raise 'Unknown order'

    # Take bed elevation into account when water heights are small
    balance_deep_and_shallow(domain)

    # Compute edge values by interpolation
    for name in domain.conserved_quantities:
        Q = domain.quantities[name]
        Q.interpolate_from_vertices_to_edges()

##
# @brief Distribution from centroids to edges specific to the SWW eqn.
# @param domain The domain to operate on.
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
            raise 'Unknown order'

    balance_deep_and_shallow(domain)

    # Compute edge values by interpolation
    for name in domain.conserved_quantities:
        Q = domain.quantities[name]
        Q.interpolate_from_vertices_to_edges()

##
# @brief  Protect against infinitesimal heights and associated high velocities.
# @param domain The domain to operate on.
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

##
# @brief Wrapper for C function balance_deep_and_shallow_c().
# @param domain The domain to operate on.
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
# Boundary conditions - specific to the shallow water wave equation
################################################################################

##
# @brief Class for a reflective boundary.
# @note Inherits from Boundary.
class Reflective_boundary(Boundary):
    """Reflective boundary returns same conserved quantities as
    those present in its neighbour volume but reflected.

    This class is specific to the shallow water equation as it
    works with the momentum quantities assumed to be the second
    and third conserved quantities.
    """

    ##
    # @brief Instantiate a Reflective_boundary.
    # @param domain 
    def __init__(self, domain=None):
        Boundary.__init__(self)

        if domain is None:
            msg = 'Domain must be specified for reflective boundary'
            raise Exception, msg

        # Handy shorthands
        self.stage = domain.quantities['stage'].edge_values
        self.xmom = domain.quantities['xmomentum'].edge_values
        self.ymom = domain.quantities['ymomentum'].edge_values
        self.normals = domain.normals

        self.conserved_quantities = num.zeros(3, num.float)

    ##
    # @brief Return a representation of this instance.
    def __repr__(self):
        return 'Reflective_boundary'

    ##
    # @brief Calculate reflections (reverse outward momentum).
    # @param vol_id 
    # @param edge_id 
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


##
# @brief Class for a transmissive boundary.
# @note Inherits from Boundary.
class Transmissive_momentum_set_stage_boundary(Boundary):
    """Returns same momentum conserved quantities as
    those present in its neighbour volume.
    Sets stage by specifying a function f of time which may either be a
    vector function or a scalar function

    Example:

    def waveform(t):
        return sea_level + normalized_amplitude/cosh(t-25)**2

    Bts = Transmissive_momentum_set_stage_boundary(domain, waveform)

    Underlying domain must be specified when boundary is instantiated
    """

    ##
    # @brief Instantiate a Transmissive_momentum_set_stage_boundary.
    # @param domain
    # @param function
    def __init__(self, domain=None, function=None):
        Boundary.__init__(self)

        if domain is None:
            msg = 'Domain must be specified for this type boundary'
            raise Exception, msg

        if function is None:
            msg = 'Function must be specified for this type boundary'
            raise Exception, msg

        self.domain = domain
        self.function = function

    ##
    # @brief Return a representation of this instance.
    def __repr__(self):
        return 'Transmissive_momentum_set_stage_boundary(%s)' %self.domain

    ##
    # @brief Calculate transmissive results.
    # @param vol_id 
    # @param edge_id 
    def evaluate(self, vol_id, edge_id):
        """Transmissive momentum set stage boundaries return the edge momentum
        values of the volume they serve.
        """

        q = self.domain.get_conserved_quantities(vol_id, edge = edge_id)

        t = self.domain.get_time()

        if hasattr(self.function, 'time'):
            # Roll boundary over if time exceeds
            while t > self.function.time[-1]:
                msg = 'WARNING: domain time %.2f has exceeded' % t
                msg += 'time provided in '
                msg += 'transmissive_momentum_set_stage_boundary object.\n'
                msg += 'I will continue, reusing the object from t==0'
                log.critical(msg)
                t -= self.function.time[-1]

        value = self.function(t)
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


##
# @brief Deprecated boundary class.
class Transmissive_Momentum_Set_Stage_boundary(Transmissive_momentum_set_stage_boundary):
    pass


##
# @brief Class for a transmissive boundary.
# @note Inherits from Boundary.
class Transmissive_n_momentum_zero_t_momentum_set_stage_boundary(Boundary):
    """Returns the same normal momentum as that 
    present in neighbour volume edge. Zero out the tangential momentum. 
    Sets stage by specifying a function f of time which may either be a
    vector function or a scalar function

    Example:

    def waveform(t):
        return sea_level + normalized_amplitude/cosh(t-25)**2

    Bts = Transmissive_n_momentum_zero_t_momentum_set_stage_boundary(domain, waveform)

    Underlying domain must be specified when boundary is instantiated
    """

    ##
    # @brief Instantiate a Transmissive_n_momentum_zero_t_momentum_set_stage_boundary.
    # @param domain
    # @param function
    def __init__(self, domain=None, function=None):
        Boundary.__init__(self)

        if domain is None:
            msg = 'Domain must be specified for this type boundary'
            raise Exception, msg

        if function is None:
            msg = 'Function must be specified for this type boundary'
            raise Exception, msg

        self.domain = domain
        self.function = function

    ##
    # @brief Return a representation of this instance.
    def __repr__(self):
        return 'Transmissive_n_momentum_zero_t_momentum_set_stage_boundary(%s)' %self.domain

    ##
    # @brief Calculate transmissive results.
    # @param vol_id 
    # @param edge_id 
    def evaluate(self, vol_id, edge_id):
        """Transmissive_n_momentum_zero_t_momentum_set_stage_boundary
        return the edge momentum values of the volume they serve.
        """

        q = self.domain.get_conserved_quantities(vol_id, edge = edge_id)

        normal = self.domain.get_normal(vol_id, edge_id)


        t = self.domain.get_time()

        if hasattr(self.function, 'time'):
            # Roll boundary over if time exceeds
            while t > self.function.time[-1]:
                msg = 'WARNING: domain time %.2f has exceeded' % t
                msg += 'time provided in '
                msg += 'transmissive_momentum_set_stage_boundary object.\n'
                msg += 'I will continue, reusing the object from t==0'
                log.critical(msg)
                t -= self.function.time[-1]

        value = self.function(t)
        try:
            x = float(value)
        except:
            x = float(value[0])

        ## import math
        ## if vol_id == 9433:
        ##     print 'vol_id = ',vol_id, ' edge_id = ',edge_id, 'q = ', q, ' x = ',x
        ##     print 'normal = ', normal
        ##     print 'n . p = ', (normal[0]*q[1] + normal[1]*q[2])
        ##     print 't . p = ', (normal[1]*q[1] - normal[0]*q[2])


        q[0] = x
        ndotq = (normal[0]*q[1] + normal[1]*q[2])
        q[1] = normal[0]*ndotq
        q[2] = normal[1]*ndotq

            
        return q

##
# @brief A transmissive boundary, momentum set to zero.
# @note Inherits from Bouondary.
class Transmissive_stage_zero_momentum_boundary(Boundary):
    """Return same stage as those present in its neighbour volume.
    Set momentum to zero.

    Underlying domain must be specified when boundary is instantiated
    """

    ##
    # @brief Instantiate a Transmissive (zero momentum) boundary.
    # @param domain
    def __init__(self, domain=None):
        Boundary.__init__(self)

        if domain is None:
            msg = ('Domain must be specified for '
                   'Transmissive_stage_zero_momentum boundary')
            raise Exception, msg

        self.domain = domain

    ##
    # @brief Return a representation of this instance.
    def __repr__(self):
        return 'Transmissive_stage_zero_momentum_boundary(%s)' % self.domain

    ##
    # @brief Calculate transmissive (zero momentum) results.
    # @param vol_id 
    # @param edge_id 
    def evaluate(self, vol_id, edge_id):
        """Transmissive boundaries return the edge values
        of the volume they serve.
        """

        q = self.domain.get_conserved_quantities(vol_id, edge=edge_id)

        q[1] = q[2] = 0.0
        return q


##
# @brief Class for a Dirichlet discharge boundary.
# @note Inherits from Boundary.
class Dirichlet_discharge_boundary(Boundary):
    """
    Sets stage (stage0)
    Sets momentum (wh0) in the inward normal direction.

    Underlying domain must be specified when boundary is instantiated
    """

    ##
    # @brief Instantiate a Dirichlet discharge boundary.
    # @param domain
    # @param stage0
    # @param wh0
    def __init__(self, domain=None, stage0=None, wh0=None):
        Boundary.__init__(self)

        if domain is None:
            msg = 'Domain must be specified for this type of boundary'
            raise Exception, msg

        if stage0 is None:
            raise Exception, 'Stage must be specified for this type of boundary'

        if wh0 is None:
            wh0 = 0.0

        self.domain = domain
        self.stage0 = stage0
        self.wh0 = wh0

    ##
    # @brief Return a representation of this instance.
    def __repr__(self):
        return 'Dirichlet_Discharge_boundary(%s)' % self.domain

    ##
    # @brief Calculate Dirichlet discharge boundary results.
    # @param vol_id 
    # @param edge_id 
    def evaluate(self, vol_id, edge_id):
        """Set discharge in the (inward) normal direction"""

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


# Backward compatibility
# FIXME(Ole): Deprecate
##
# @brief Deprecated
class Dirichlet_Discharge_boundary(Dirichlet_discharge_boundary):
    pass


class Inflow_boundary(Boundary):
    """Apply given flow in m^3/s to boundary segment.
    Depth and momentum is derived using Manning's formula.

    Underlying domain must be specified when boundary is instantiated
    """
    
    # FIXME (Ole): This is work in progress and definitely not finished.
    # The associated test has been disabled

    def __init__(self, domain=None, rate=0.0):
        Boundary.__init__(self)

        if domain is None:
            msg = 'Domain must be specified for '
            msg += 'Inflow boundary'
            raise Exception, msg

        self.domain = domain
        
        # FIXME(Ole): Allow rate to be time dependent as well
        self.rate = rate
        self.tag = None # Placeholder for tag associated with this object.

    def __repr__(self):
        return 'Inflow_boundary(%s)' %self.domain

    def evaluate(self, vol_id, edge_id):
        """Apply inflow rate at each edge of this boundary
        """
        
        # First find all segments having the same tag is vol_id, edge_id
        # This will be done the first time evaluate is called.
        if self.tag is None:
            boundary = self.domain.boundary
            self.tag = boundary[(vol_id, edge_id)]        
            
            # Find total length of boundary with this tag
            length = 0.0
            for v_id, e_id in boundary:
                if self.tag == boundary[(v_id, e_id)]:
                    length += self.domain.mesh.get_edgelength(v_id, e_id)            

            self.length = length
            self.average_momentum = self.rate/length
            
            
        # Average momentum has now been established across this boundary
        # Compute momentum in the inward normal direction 
        
        inward_normal = -self.domain.mesh.get_normal(vol_id, edge_id)       
        xmomentum, ymomentum = self.average_momentum * inward_normal
            
        # Compute depth based on Manning's formula v = 1/n h^{2/3} sqrt(S)
        # Where v is velocity, n is manning's coefficient, h is depth and S is the slope into the domain. 
        # Let mu be the momentum (vh), then this equation becomes: mu = 1/n h^{5/3} sqrt(S) 
        # from which we can isolate depth to get
        # h = (mu n/sqrt(S) )^{3/5} 
        
        slope = 0 # get gradient for this triangle dot normal
        
        # get manning coef from this triangle
        friction = self.domain.get_quantity('friction').get_values(location='edges', 
                                                                   indices=[vol_id])[0]
        mannings_n = friction[edge_id]

        if slope > epsilon and mannings_n > epsilon:
            depth = pow(self.average_momentum * mannings_n/math.sqrt(slope), 3.0/5) 
        else:
            depth = 1.0
            
        # Elevation on this edge    
        
        z = self.domain.get_quantity('elevation').get_values(location='edges', 
                                                             indices=[vol_id])[0]
        elevation = z[edge_id]
            
        # Assign conserved quantities and return
        q = num.array([elevation + depth, xmomentum, ymomentum], num.Float)
        return q


        
    
            
        
class Field_boundary(Boundary):
    """Set boundary from given field represented in an sww file containing
    values for stage, xmomentum and ymomentum.

    Optionally, the user can specify mean_stage to offset the stage provided
    in the sww file.

    This function is a thin wrapper around the generic File_boundary. The
    difference between the file_boundary and field_boundary is only that the
    field_boundary will allow you to change the level of the stage height when
    you read in the boundary condition. This is very useful when running
    different tide heights in the same area as you need only to convert one
    boundary condition to a SWW file, ideally for tide height of 0 m
    (saving disk space). Then you can use field_boundary to read this SWW file
    and change the stage height (tide) on the fly depending on the scenario.
    """

    ##
    # @brief Construct an instance of a 'field' boundary.
    # @param filename Name of SWW file containing stage, x and ymomentum.
    # @param domain Shallow water domain for which the boundary applies.
    # @param mean_stage Mean water level added to stage derived from SWW file.
    # @param time_thinning Time step thinning factor.
    # @param time_limit 
    # @param boundary_polygon 
    # @param default_boundary None or an instance of Boundary.
    # @param use_cache True if caching is to be used.
    # @param verbose True if this method is to be verbose.
    def __init__(self,
                 filename,
                 domain,
                 mean_stage=0.0,
                 time_thinning=1,
                 time_limit=None,
                 boundary_polygon=None,
                 default_boundary=None,
                 use_cache=False,
                 verbose=False):
        """Constructor

        filename: Name of sww file
        domain: pointer to shallow water domain for which the boundary applies
        mean_stage: The mean water level which will be added to stage derived
                    from the boundary condition
        time_thinning: Will set how many time steps from the sww file read in
                       will be interpolated to the boundary. For example if
                       the sww file has 1 second time steps and is 24 hours
                       in length it has 86400 time steps. If you set
                       time_thinning to 1 it will read all these steps.
                       If you set it to 100 it will read every 100th step eg
                       only 864 step. This parameter is very useful to increase
                       the speed of a model run that you are setting up
                       and testing.

        default_boundary: Must be either None or an instance of a
                          class descending from class Boundary.
                          This will be used in case model time exceeds
                          that available in the underlying data.

                          Note that mean_stage will also be added to this.
                                                
        use_cache:
        verbose:
        """

        # Create generic file_boundary object
        self.file_boundary = File_boundary(filename,
                                           domain,
                                           time_thinning=time_thinning,
                                           time_limit=time_limit,
                                           boundary_polygon=boundary_polygon,
                                           default_boundary=default_boundary,
                                           use_cache=use_cache,
                                           verbose=verbose)

        # Record information from File_boundary
        self.F = self.file_boundary.F
        self.domain = self.file_boundary.domain

        # Record mean stage
        self.mean_stage = mean_stage

    ##
    # @note Generate a string representation of this instance.
    def __repr__(self):
        return 'Field boundary'

    ##
    # @brief Calculate 'field' boundary results.
    # @param vol_id 
    # @param edge_id 
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


################################################################################
# Standard forcing terms
################################################################################

##
# @brief Apply gravitational pull in the presence of bed slope.
# @param domain The domain to apply gravity to.
# @note Wrapper for C function gravity_c().
def gravity(domain):
    """Apply gravitational pull in the presence of bed slope
    Wrapper calls underlying C implementation
    """

    from shallow_water_ext import gravity as gravity_c

    xmom = domain.quantities['xmomentum'].explicit_update
    ymom = domain.quantities['ymomentum'].explicit_update

    stage = domain.quantities['stage']
    elevation = domain.quantities['elevation']

    h = stage.centroid_values - elevation.centroid_values
    z = elevation.vertex_values

    x = domain.get_vertex_coordinates()
    g = domain.g

    gravity_c(g, h, z, x, xmom, ymom)    #, 1.0e-6)

##
# @brief Apply friction to a surface (implicit).
# @param domain The domain to apply Manning friction to.
# @note Wrapper for C function manning_friction_c().
def manning_friction_implicit(domain):
    """Apply (Manning) friction to water momentum
    Wrapper for c version
    """

    from shallow_water_ext import manning_friction as manning_friction_c

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

    manning_friction_c(g, eps, w, z, uh, vh, eta, xmom_update, ymom_update)

##
# @brief Apply friction to a surface (explicit).
# @param domain The domain to apply Manning friction to.
# @note Wrapper for C function manning_friction_c().
def manning_friction_explicit(domain):
    """Apply (Manning) friction to water momentum
    Wrapper for c version
    """

    from shallow_water_ext import manning_friction as manning_friction_c

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

    manning_friction_c(g, eps, w, z, uh, vh, eta, xmom_update, ymom_update)


# FIXME (Ole): This was implemented for use with one of the analytical solutions (Sampson?)
##
# @brief Apply linear friction to a surface.
# @param domain The domain to apply Manning friction to.
# @note Is this still used (30 Oct 2007)?
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

def depth_dependent_friction(domain, default_friction,
                             surface_roughness_data,
                             verbose=False):
    """Returns an array of friction values for each wet element adjusted for depth.

    Inputs:
        domain - computational domain object
        default_friction - depth independent bottom friction
        surface_roughness_data - N x 5 array of n0, d1, n1, d2, n2 values for each
        friction region.

    Outputs:
        wet_friction - Array that can be used directly to update friction as follows:
                       domain.set_quantity('friction', wet_friction)

        
        
    """

    import numpy as num
    
    # Create a temp array to store updated depth dependent friction for wet elements
    # EHR this is outwardly inneficient but not obvious how to avoid recreating each call??????
    N=len(domain)
    wet_friction    = num.zeros(N, num.float)
    wet_friction[:] = default_n0   # Initially assign default_n0 to all array so sure have no zeros values
    
    
    depth = domain.create_quantity_from_expression('stage - elevation')  # create depth instance for this timestep
    # Recompute depth as vector  
    d = depth.get_values(location='centroids')
 
    # rebuild the 'friction' values adjusted for depth at this instant
    for i in domain.get_wet_elements():                                  # loop for each wet element in domain
        
        # Get roughness data for each element
        n0 = float(surface_roughness_data[i,0])
        d1 = float(surface_roughness_data[i,1])
        n1 = float(surface_roughness_data[i,2])
        d2 = float(surface_roughness_data[i,3])
        n2 = float(surface_roughness_data[i,4])
        
        
        # Recompute friction values from depth for this element 
               
        if d[i]   <= d1: 
            depth_dependent_friction = n1
        elif d[i] >= d2:
            depth_dependent_friction = n2
        else:
            depth_dependent_friction = n1+((n2-n1)/(d2-d1))*(d[i]-d1)
            
        # check sanity of result
        if (depth_dependent_friction  < 0.010 or depth_dependent_friction > 9999.0) :
            log.critical('%s >>>> WARNING: computed depth_dependent friction '
                         'out of range, ddf%f, n1=%f, n2=%f'
                         % (model_data.basename,
                            depth_dependent_friction, n1, n2))
        
        # update depth dependent friction  for that wet element
        wet_friction[i] = depth_dependent_friction
        
    # EHR add code to show range of 'friction across domain at this instant as sanity check?????????
    
    if verbose :
        nvals=domain.get_quantity('friction').get_values(location='centroids')        # return array of domain nvals
        n_min=min(nvals)
        n_max=max(nvals)
        
        log.critical('         ++++ calculate_depth_dependent_friction - '
                     'Updated friction - range  %7.3f to %7.3f'
                     % (n_min, n_max))
    
    return wet_friction


################################################################################
# Experimental auxiliary functions
################################################################################

##
# @brief Check forcefield parameter.
# @param f Object to check.
# @note 'f' may be a callable object or a scalar value.
def check_forcefield(f):
    """Check that force object is as expected.
    
    Check that f is either:
    1: a callable object f(t,x,y), where x and y are vectors
       and that it returns an array or a list of same length
       as x and y
    2: a scalar
    """

    if callable(f):
        N = 3
        x = num.ones(3, num.float)
        y = num.ones(3, num.float)
        try:
            q = f(1.0, x=x, y=y)
        except Exception, e:
            msg = 'Function %s could not be executed:\n%s' %(f, e)
            # FIXME: Reconsider this semantics
            raise Exception, msg

        try:
            q = num.array(q, num.float)
        except:
            msg = ('Return value from vector function %s could not '
                   'be converted into a numeric array of floats.\nSpecified '
                   'function should return either list or array.' % f)
            raise Exception, msg

        # Is this really what we want?
        # info is "(func name, filename, defining line)"
        func_info = (f.func_name, f.func_code.co_filename,
                     f.func_code.co_firstlineno)
        func_msg = 'Function %s (defined in %s, line %d)' % func_info
        try:
            result_len = len(q)
        except:
            msg = '%s must return vector' % func_msg
            self.fail(msg)
        msg = '%s must return vector of length %d' % (func_msg, N)
        assert result_len == N, msg
    else:
        try:
            f = float(f)
        except:
            msg = ('Force field %s must be a scalar value coercible to float.'
                   % str(f))
            raise Exception, msg

    return f


##
# Class to apply a wind stress to a domain.
class Wind_stress:
    """Apply wind stress to water momentum in terms of
    wind speed [m/s] and wind direction [degrees]
    """

    ##
    # @brief Create an instance of Wind_stress.
    # @param *args 
    # @param **kwargs 
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
               raise Exception, 'Assumes two keyword arguments: s=..., phi=....'

        self.speed = check_forcefield(s)
        self.phi = check_forcefield(phi)

        self.const = eta_w*rho_a/rho_w

    ##
    # @brief 'execute' this class instance.
    # @param domain 
    def __call__(self, domain):
        """Evaluate windfield based on values found in domain"""

        from math import pi, cos, sin, sqrt

        xmom_update = domain.quantities['xmomentum'].explicit_update
        ymom_update = domain.quantities['ymomentum'].explicit_update

        N = len(domain)    # number_of_triangles
        t = domain.time

        if callable(self.speed):
            xc = domain.get_centroid_coordinates()
            s_vec = self.speed(t, xc[:,0], xc[:,1])
        else:
            # Assume s is a scalar
            try:
                s_vec = self.speed * num.ones(N, num.float)
            except:
                msg = 'Speed must be either callable or a scalar: %s' %self.s
                raise msg

        if callable(self.phi):
            xc = domain.get_centroid_coordinates()
            phi_vec = self.phi(t, xc[:,0], xc[:,1])
        else:
            # Assume phi is a scalar

            try:
                phi_vec = self.phi * num.ones(N, num.float)
            except:
                msg = 'Angle must be either callable or a scalar: %s' %self.phi
                raise msg

        assign_windfield_values(xmom_update, ymom_update,
                                s_vec, phi_vec, self.const)


##
# @brief Assign wind field values
# @param xmom_update 
# @param ymom_update
# @param s_vec 
# @param phi_vec 
# @param const 
def assign_windfield_values(xmom_update, ymom_update,
                            s_vec, phi_vec, const):
    """Python version of assigning wind field to update vectors.
    A C version also exists (for speed)
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


##
# @brief A class for a general explicit forcing term.
class General_forcing:
    """General explicit forcing term for update of quantity

    This is used by Inflow and Rainfall for instance


    General_forcing(quantity_name, rate, center, radius, polygon)

    domain:     ANUGA computational domain
    quantity_name: Name of quantity to update.
                   It must be a known conserved quantity.

    rate [?/s]: Total rate of change over the specified area.
                This parameter can be either a constant or a
                function of time. Positive values indicate increases,
                negative values indicate decreases.
                Rate can be None at initialisation but must be specified
                before forcing term is applied (i.e. simulation has started).

    center [m]: Coordinates at center of flow point
    radius [m]: Size of circular area
    polygon:    Arbitrary polygon
    default_rate: Rate to be used if rate fails (e.g. if model time exceeds its data)
                  Admissible types: None, constant number or function of t


    Either center, radius or polygon can be specified but not both.
    If neither are specified the entire domain gets updated.
    All coordinates to be specified in absolute UTM coordinates (x, y) assuming the zone of domain.

    Inflow or Rainfall for examples of use
    """


    # FIXME (AnyOne) : Add various methods to allow spatial variations

    ##
    # @brief Create an instance of this forcing term.
    # @param domain 
    # @param quantity_name 
    # @param rate 
    # @param center 
    # @param radius 
    # @param polygon 
    # @param default_rate 
    # @param verbose 
    def __init__(self,
                 domain,
                 quantity_name,
                 rate=0.0,
                 center=None,
                 radius=None,
                 polygon=None,
                 default_rate=None,
                 verbose=False):

        from math import pi, cos, sin

        if center is None:
            msg = 'I got radius but no center.'
            assert radius is None, msg

        if radius is None:
            msg += 'I got center but no radius.'
            assert center is None, msg

        self.domain = domain
        self.quantity_name = quantity_name
        self.rate = rate
        self.center = ensure_numeric(center)
        self.radius = radius
        self.polygon = polygon
        self.verbose = verbose
        self.value = 0.0    # Can be used to remember value at
                            # previous timestep in order to obtain rate

        # Get boundary (in absolute coordinates)
        bounding_polygon = domain.get_boundary_polygon()

        # Update area if applicable
        if center is not None and radius is not None:
            assert len(center) == 2
            msg = 'Polygon cannot be specified when center and radius are'
            assert polygon is None, msg

            # Check that circle center lies within the mesh.
            msg = 'Center %s specified for forcing term did not' % str(center)
            msg += 'fall within the domain boundary.'
            assert is_inside_polygon(center, bounding_polygon), msg

            # Check that circle periphery lies within the mesh.
            N = 100
            periphery_points = []
            for i in range(N):
                theta = 2*pi*i/100

                x = center[0] + radius*cos(theta)
                y = center[1] + radius*sin(theta)

                periphery_points.append([x,y])

            for point in periphery_points:
                msg = 'Point %s on periphery for forcing term' % str(point)
                msg += ' did not fall within the domain boundary.'
                assert is_inside_polygon(point, bounding_polygon), msg

        if polygon is not None:
            # Check that polygon lies within the mesh.
            for point in self.polygon:
                msg = 'Point %s in polygon for forcing term' % str(point)
                msg += ' did not fall within the domain boundary.'
                assert is_inside_polygon(point, bounding_polygon), msg

        # Pointer to update vector
        self.update = domain.quantities[self.quantity_name].explicit_update

        # Determine indices in flow area
        N = len(domain)
        points = domain.get_centroid_coordinates(absolute=True)

        # Calculate indices in exchange area for this forcing term
        self.exchange_indices = None
        if self.center is not None and self.radius is not None:
            # Inlet is circular
            inlet_region = 'center=%s, radius=%s' % (self.center, self.radius)

            self.exchange_indices = []
            for k in range(N):
                x, y = points[k,:]    # Centroid

                c = self.center
                if ((x-c[0])**2+(y-c[1])**2) < self.radius**2:
                    self.exchange_indices.append(k)

        if self.polygon is not None:
            # Inlet is polygon
            inlet_region = 'polygon=%s' % (self.polygon) 
            self.exchange_indices = inside_polygon(points, self.polygon)

        if self.exchange_indices is None:
            self.exchange_area = polygon_area(bounding_polygon)
        else:    
            if len(self.exchange_indices) == 0:
                msg = 'No triangles have been identified in '
                msg += 'specified region: %s' % inlet_region
                raise Exception, msg

            # Compute exchange area as the sum of areas of triangles identified
            # by circle or polygon
            self.exchange_area = 0.0
            for i in self.exchange_indices:
                self.exchange_area += domain.areas[i]
            

        msg = 'Exchange area in forcing term'
        msg += ' has area = %f' %self.exchange_area
        assert self.exchange_area > 0.0            
            
                

            
        # Check and store default_rate
        msg = ('Keyword argument default_rate must be either None '
               'or a function of time.\nI got %s.' % str(default_rate))
        assert (default_rate is None or
                type(default_rate) in [IntType, FloatType] or
                callable(default_rate)), msg

        if default_rate is not None:
            # If it is a constant, make it a function
            if not callable(default_rate):
                tmp = default_rate
                default_rate = lambda t: tmp

            # Check that default_rate is a function of one argument
            try:
                default_rate(0.0)
            except:
                raise Exception, msg

        self.default_rate = default_rate
        self.default_rate_invoked = False    # Flag

    ##
    # @brief Execute this instance.
    # @param domain 
    def __call__(self, domain):
        """Apply inflow function at time specified in domain, update stage"""

        # Call virtual method allowing local modifications
        t = domain.get_time()
        try:
            rate = self.update_rate(t)
        except Modeltime_too_early, e:
            raise Modeltime_too_early, e
        except Modeltime_too_late, e:
            if self.default_rate is None:
                raise Exception, e    # Reraise exception
            else:
                # Pass control to default rate function
                rate = self.default_rate(t)

                if self.default_rate_invoked is False:
                    # Issue warning the first time
                    msg = ('%s\n'
                           'Instead I will use the default rate: %s\n'
                           'Note: Further warnings will be supressed'
                           % (str(e), str(self.default_rate)))
                    warn(msg)

                    # FIXME (Ole): Replace this crude flag with
                    # Python's ability to print warnings only once.
                    # See http://docs.python.org/lib/warning-filter.html
                    self.default_rate_invoked = True

        if rate is None:
            msg = ('Attribute rate must be specified in General_forcing '
                   'or its descendants before attempting to call it')
            raise Exception, msg

        # Now rate is a number
        if self.verbose is True:
            log.critical('Rate of %s at time = %.2f = %f'
                         % (self.quantity_name, domain.get_time(), rate))

        if self.exchange_indices is None:
            self.update[:] += rate
        else:
            # Brute force assignment of restricted rate
            for k in self.exchange_indices:
                self.update[k] += rate

    ##
    # @brief Update the internal rate.
    # @param t A callable or scalar used to set the rate.
    # @return The new rate.
    def update_rate(self, t):
        """Virtual method allowing local modifications by writing an
        overriding version in descendant
        """

        if callable(self.rate):
            rate = self.rate(t)
        else:
            rate = self.rate

        return rate

    ##
    # @brief Get values for the specified quantity.
    # @param quantity_name Name of the quantity of interest.
    # @return The value(s) of the quantity.
    # @note If 'quantity_name' is None, use self.quantity_name.
    def get_quantity_values(self, quantity_name=None):
        """Return values for specified quantity restricted to opening

        Optionally a quantity name can be specified if values from another
        quantity is sought
        """

        if quantity_name is None:
            quantity_name = self.quantity_name

        q = self.domain.quantities[quantity_name]
        return q.get_values(location='centroids',
                            indices=self.exchange_indices)

    ##
    # @brief Set value for the specified quantity.
    # @param val The value object used to set value.
    # @param quantity_name Name of the quantity of interest.
    # @note If 'quantity_name' is None, use self.quantity_name.
    def set_quantity_values(self, val, quantity_name=None):
        """Set values for specified quantity restricted to opening

        Optionally a quantity name can be specified if values from another
        quantity is sought
        """

        if quantity_name is None:
            quantity_name = self.quantity_name

        q = self.domain.quantities[self.quantity_name]
        q.set_values(val,
                     location='centroids',
                     indices=self.exchange_indices)


##
# @brief A class for rainfall forcing function.
# @note Inherits from General_forcing.
class Rainfall(General_forcing):
    """Class Rainfall - general 'rain over entire domain' forcing term.

    Used for implementing Rainfall over the entire domain.

        Current Limited to only One Gauge..

        Need to add Spatial Varying Capability
        (This module came from copying and amending the Inflow Code)

    Rainfall(rain)

    domain
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
                  stage quantity.

    polygon: Specifies a polygon to restrict the rainfall.

    Examples
    How to put them in a run File...

    #------------------------------------------------------------------------
    # Setup specialised forcing terms
    #------------------------------------------------------------------------
    # This is the new element implemented by Ole and Rudy to allow direct
    # input of Rainfall in mm/s

    catchmentrainfall = Rainfall(rain=file_function('Q100_2hr_Rain.tms'))
                        # Note need path to File in String.
                        # Else assumed in same directory

    domain.forcing_terms.append(catchmentrainfall)
    """

    ##
    # @brief Create an instance of the class.
    # @param domain Domain of interest.
    # @param rate Total rain rate over the specified domain (mm/s).
    # @param center 
    # @param radius 
    # @param polygon Polygon  to restrict rainfall.
    # @param default_rate 
    # @param verbose True if this instance is to be verbose.
    def __init__(self,
                 domain,
                 rate=0.0,
                 center=None,
                 radius=None,
                 polygon=None,
                 default_rate=None,
                 verbose=False):

        # Converting mm/s to m/s to apply in ANUGA)
        if callable(rate):
            rain = lambda t: rate(t)/1000.0
        else:
            rain = rate/1000.0

        if default_rate is not None:
            if callable(default_rate):
                default_rain = lambda t: default_rate(t)/1000.0
            else:
                default_rain = default_rate/1000.0
        else:
            default_rain = None

            
            
        General_forcing.__init__(self,
                                 domain,
                                 'stage',
                                 rate=rain,
                                 center=center,
                                 radius=radius,
                                 polygon=polygon,
                                 default_rate=default_rain,
                                 verbose=verbose)


##
# @brief A class for inflow (rain and drain) forcing function.
# @note Inherits from General_forcing.
class Inflow(General_forcing):
    """Class Inflow - general 'rain and drain' forcing term.

    Useful for implementing flows in and out of the domain.

    Inflow(flow, center, radius, polygon)

    domain
    rate [m^3/s]: Total flow rate over the specified area.
                  This parameter can be either a constant or a
                  function of time. Positive values indicate inflow,
                  negative values indicate outflow.
                  The specified flow will be divided by the area of
                  the inflow region and then applied to update stage.
    center [m]: Coordinates at center of flow point
    radius [m]: Size of circular area
    polygon:    Arbitrary polygon.

    Either center, radius or polygon must be specified

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


    #------------------------------------------------------------------------
    # Setup specialised forcing terms
    #------------------------------------------------------------------------
    # This is the new element implemented by Ole to allow direct input
    # of Inflow in m^3/s

    hydrograph = Inflow(center=(320, 300), radius=10,
                        rate=file_function('Q/QPMF_Rot_Sub13.tms'))

    domain.forcing_terms.append(hydrograph)
    """

    ##
    # @brief Create an instance of the class.
    # @param domain Domain of interest.
    # @param rate Total rain rate over the specified domain (mm/s).
    # @param center 
    # @param radius 
    # @param polygon Polygon  to restrict rainfall.
    # @param default_rate 
    # @param verbose True if this instance is to be verbose.
    def __init__(self,
                 domain,
                 rate=0.0,
                 center=None,
                 radius=None,
                 polygon=None,
                 default_rate=None,
                 verbose=False):
        # Create object first to make area is available
        General_forcing.__init__(self,
                                 domain,
                                 'stage',
                                 rate=rate,
                                 center=center,
                                 radius=radius,
                                 polygon=polygon,
                                 default_rate=default_rate,
                                 verbose=verbose)

    ##
    # @brief Update the instance rate.
    # @param t New rate object.
    def update_rate(self, t):
        """Virtual method allowing local modifications by writing an
        overriding version in descendant

        This one converts m^3/s to m/s which can be added directly
        to 'stage' in ANUGA
        """

        if callable(self.rate):
            _rate = self.rate(t)/self.exchange_area
        else:
            _rate = self.rate/self.exchange_area

        return _rate


##
# @brief A class for creating cross sections.
# @note Inherits from General_forcing.
class Cross_section:
    """Class Cross_section - a class to setup a cross section from
    which you can then calculate flow and energy through cross section


    Cross_section(domain, polyline)

    domain:
    polyline: Representation of desired cross section - it may contain
              multiple sections allowing for complex shapes. Assume
              absolute UTM coordinates.
              Format [[x0, y0], [x1, y1], ...]
    verbose: 
    """

    ##
    # @brief Create an instance of the class.
    # @param domain Domain of interest.
    # @param polyline Polyline defining cross section
    # @param verbose True if this instance is to be verbose.
    def __init__(self,
                 domain,
                 polyline=None,
                 verbose=False):
        
        self.domain = domain
        self.polyline = polyline
        self.verbose = verbose
        
        # Find all intersections and associated triangles.
        self.segments = self.domain.get_intersecting_segments(self.polyline,
                                                              use_cache=True,
                                                              verbose=self.verbose)
        
        # Get midpoints
        self.midpoints = segment_midpoints(self.segments)

        # Make midpoints Geospatial instances
        self.midpoints = ensure_geospatial(self.midpoints, self.domain.geo_reference)


    ##
    # @brief calculate current flow through cross section
    def get_flow_through_cross_section(self):
        """ Output: Total flow [m^3/s] across cross section.
        """

        # Get interpolated values
        xmomentum = self.domain.get_quantity('xmomentum')
        ymomentum = self.domain.get_quantity('ymomentum')

        uh = xmomentum.get_values(interpolation_points=self.midpoints,
                                  use_cache=True)
        vh = ymomentum.get_values(interpolation_points=self.midpoints,
                                  use_cache=True)

        # Compute and sum flows across each segment
        total_flow = 0
        for i in range(len(uh)):
            # Inner product of momentum vector with segment normal [m^2/s]
            normal = self.segments[i].normal
            normal_momentum = uh[i]*normal[0] + vh[i]*normal[1]

            # Flow across this segment [m^3/s]
            segment_flow = normal_momentum*self.segments[i].length

            # Accumulate
            total_flow += segment_flow

        return total_flow
 



################################################################################
# Initialise module
################################################################################

from anuga.utilities import compile
if compile.can_use_C_extension('shallow_water_ext.c'):
    # Underlying C implementations can be accessed
    from shallow_water_ext import rotate, assign_windfield_values
else:
    msg = 'C implementations could not be accessed by %s.\n ' % __file__
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
        if os.name == 'posix' and os.uname()[4] in ['x86_64', 'ia64']:
            pass
            #Psyco isn't supported on 64 bit systems, but it doesn't matter
        else:
            msg = ('WARNING: psyco (speedup) could not be imported, '
                   'you may want to consider installing it')
            log.critical(msg)
    else:
        psyco.bind(Domain.distribute_to_vertices_and_edges)
        psyco.bind(Domain.compute_fluxes)


if __name__ == "__main__":
    pass

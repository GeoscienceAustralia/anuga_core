


from builtins import range
from anuga.shallow_water.shallow_water_domain import *
from anuga.shallow_water.shallow_water_domain import Domain as Sww_domain


##############################################################################
# Shallow Water Balanced Domain -- alternative implementation
# 
# FIXME: Following the methods in CITE MODSIM PAPER
#
##############################################################################

class Domain(Sww_domain):

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
                 starttime=0.0,
                 processor=0,
                 numproc=1,
                 number_of_full_nodes=None,
                 number_of_full_triangles=None):

        conserved_quantities = [ 'stage', 'xmomentum', 'ymomentum']

        evolved_quantities = [ 'stage', 'xmomentum', 'ymomentum']         
        other_quantities   = [ 'elevation', 'friction', 'height', 
                               'xvelocity', 'yvelocity', 'x', 'y' ]

        
        Sww_domain.__init__(self,
                            coordinates = coordinates,
                            vertices = vertices,
                            boundary = boundary,
                            tagged_elements = tagged_elements,
                            geo_reference = geo_reference,
                            use_inscribed_circle = use_inscribed_circle,
                            mesh_filename = mesh_filename,
                            use_cache = use_cache,
                            verbose = verbose,
                            conserved_quantities = conserved_quantities,
                            evolved_quantities = evolved_quantities,
                            other_quantities = other_quantities,
                            full_send_dict = full_send_dict,
                            ghost_recv_dict = ghost_recv_dict,
                            starttime = starttime,
                            processor = processor,
                            numproc = numproc,
                            number_of_full_nodes = number_of_full_nodes,
                            number_of_full_triangles = number_of_full_triangles)
       
        #------------------------------------------------
        # set some defaults
        # Most of these override the options in config.py
        #------------------------------------------------

        self.set_tsunami_defaults()


    def set_tsunami_defaults(self):


        self.set_CFL(1.0)
        self.set_use_kinematic_viscosity(False)
        self.timestepping_method='rk2'
        # The following allows storage of the negative depths associated with this method
        self.minimum_storable_height=-99999999999.0

        self.use_edge_limiter=True
        self.default_order=2
        self.extrapolate_velocity_second_order=True

        # Note that the extrapolation method used in quantity_ext.c (e.g.
        # extrapolate_second_order_and_limit_by_edge) uses a constant value for
        # all the betas.
        self.beta_w=1.0
        self.beta_w_dry=1.0
        self.beta_uh=1.0
        self.beta_uh_dry=1.0
        self.beta_vh=1.0
        self.beta_vh_dry=1.0

        #self.optimise_dry_cells=True
        # "self.optimise_dry_cells=False" presently ensures that the stage is
        # always >= minimum_bed_edge value.  Actually, we should only need to
        # apply 'False' on the very first time step (to deal with stage values
        # that were initialised below the bed by the user). After this, the
        # algorithm should take care of itself, and 'True' should be okay.
        self.optimise_dry_cells=False

        # Because gravity is treated within the flux function,
        # we remove it from the forcing terms.
        #self.forcing_terms.remove(gravity)

        # We need the edge_coordinates for the extrapolation
        self.edge_coordinates=self.get_edge_midpoint_coordinates()

        # We demand that vertex values are stored uniquely
        self.set_store_vertices_smoothly(False)

        self.maximum_allowed_speed=0.0
        #self.forcing_terms.append(manning_friction_explicit)
        #self.forcing_terms.remove(manning_friction_implicit)

        print('##########################################################################')
        print('#')
        print('# Using anuga_tsunami solver in anuga_work/development/gareth/experimental/anuga_tsunami/')
        print("#")
        print("# This solver is experimental. Here are some tips on using it")
        print("#")
        print("# 1) When plotting outputs, I strongly suggest you examine centroid values, not vertex values")
        print("# , as the latter can be completely misleading near strong gradients in the flow. ")
        print("# There is a plot_util.py script in anuga_core/utilities/ which might help you extract")
        print("# quantities at centroid values from sww files.")
        print("# Note that to accuractely compute centroid values from sww files, the files need to store ")
        print("# vertices uniquely. This makes for large sww files (3x), but is the price to pay for the right answer")
        print("# (unless we alter IO to allow centroids to be written to sww files, which would then affect")
        print("# ANUGA viewer as well -- I expect this would be lots of work)")
        print("#")
        print("# 2) In field scale applications (where the depth is typically > 1m), I suggest you set")
        print("# domain.minimum_allowed_height=0.01 (the default is 1.0e-3). ")
        print("#")
        print("# 3) This solver is not expected to perform well in problems with very")
        print("# shallow water flowing down steep slopes (such that the stage_centroid_value ")
        print("# is less than the maximum bed_edge_value on a given triangle). However, analytical tests")
        print("# suggest it can do typical wetting/drying situations very well (parabolic oscillations test case) ")
        print("#")
        print("# 4) This solver allows the stage_centroid_value to drop to slightly below the minimum bed_vertex_value")
        print("# on it's triangle. In other ANUGA versions (e.g. 1.2.1), the limit would be the")
        print("# bed_centroid_value. This means that triangles store slightly more water than they are")
        print("# typically interpreted to store, which might have significance in some applications.")
        print("#")
        print("# You will probably be able to tell this is causing you problems by convergence testing")
        print("#")
        print('# 5) Note that many options in config.py have been overridden by the solver -- I have ')
        print('# deliberately attempted to get the solver to perform well with consistent values of ')
        print('# these parameters -- so I would advise against changing them unless you at least check that ')
        print('# it really does improve things')
        print('#')
        print('##########################################################################')


    #-------------------------------
    # Specialisation of Evolve
    #-------------------------------
    # This alters the 'evolve' method from shallow_water_domain.py so that just
    # before the file-output step, the vertex values are replaced with the
    # centroid values. 
    # This is useful so that we can unambigously know what the centroid values
    # of various parameters are

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

        #print 'Into Evolve'

        msg = 'Attribute self.beta_w must be in the interval [0, 2]'
        assert 0 <= self.beta_w <= 2.0, msg

        # Initial update of vertex and edge values before any STORAGE
        # and or visualisation.
        # This is done again in the initialisation of the Generic_Domain
        # evolve loop but we do it here to ensure the values are ok for storage.
        self.distribute_to_vertices_and_edges()

        if self.store is True and self.get_time() == self.get_starttime():
            self.initialise_storage()
        #print 'Into Generic_Domain Evolve'
        # Call basic machinery from parent class
        for t in Generic_Domain.evolve(self, yieldstep=yieldstep,
                                       finaltime=finaltime, duration=duration,
                                       skip_initial_step=skip_initial_step):
            #print 'Out of Generic_Domain Evolve'
            # Store model data, e.g. for subsequent visualisation
            if self.store is True:
                # FIXME: Extrapolation is done before writing, because I had
                # trouble correctly computing the centroid values from the
                # vertex values (an 'average' did not seem to work correctly in
                # very shallow cells -- there was a discrepency between
                # domain.quantity['blah'].centroid_values and the values
                # computed from the sww using the vertex averge). There should
                # be a much much more disk-efficient way to store the centroid
                # values than this
                self.extrapolate_second_order_edge_sw()

                # Store the timestep
                self.store_timestep()

            # Pass control on to outer loop for more specific actions
            yield(t)





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
        from .swb2_domain_ext import compute_fluxes_ext_central \
                                      as compute_fluxes_ext

        #print "."

        # Shortcuts
        Stage = domain.quantities['stage']
        Xmom = domain.quantities['xmomentum']
        Ymom = domain.quantities['ymomentum']
        Bed = domain.quantities['elevation']

        timestep = float(sys.maxsize)

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
                                           domain.boundary_flux_type,
                                           Stage.explicit_update,
                                           Xmom.explicit_update,
                                           Ymom.explicit_update,
                                           domain.already_computed_flux,
                                           domain.max_speed,
                                           int(domain.optimise_dry_cells),
                                           Stage.centroid_values,
                                           Bed.centroid_values, 
                                           Bed.vertex_values)

        #import pdb
        #pdb.set_trace()

        domain.flux_timestep = flux_timestep


    def protect_against_infinitesimal_and_negative_heights(domain):
        """protect against infinitesimal heights and associated high velocities"""

        from .swb2_domain_ext import protect
        #print'using swb2_protect_against ..'

        # shortcuts
        wc = domain.quantities['stage'].centroid_values
        wv = domain.quantities['stage'].vertex_values
        zc = domain.quantities['elevation'].centroid_values
        zv = domain.quantities['elevation'].vertex_values
        xmomc = domain.quantities['xmomentum'].centroid_values
        ymomc = domain.quantities['ymomentum'].centroid_values
        areas = domain.areas

        protect(domain.minimum_allowed_height, domain.maximum_allowed_speed,
                domain.epsilon, wc, wv, zc,zv, xmomc, ymomc, areas)

    def conserved_values_to_evolved_values(self, q_cons, q_evol):
        """Mapping between conserved quantities and the evolved quantities.
        Used where we have a boundary condition which works with conserved
        quantities and we now want to use them for the new
        code using the evolved quantities

        Typically the old boundary condition will set the values in q_cons,

        q_evol on input will have the values of the evolved quantities at the
        edge point.

        """

        wc  = q_cons[0]
        uhc = q_cons[1]
        vhc = q_cons[2]


        q_evol[0]  = wc
        q_evol[1]  = uhc
        q_evol[2]  = vhc

        return q_evol

    def distribute_to_vertices_and_edges(self):
        """ Call correct module function """



        if self.use_edge_limiter:
            #distribute_using_edge_limiter(self)
            self.distribute_using_edge_limiter()
        else:
            #distribute_using_vertex_limiter(self)
            self.distribute_using_vertex_limiter()


    def distribute_using_vertex_limiter(domain):
        """Distribution from centroids to vertices specific to the SWW equation.

        Precondition:
          All quantities defined at centroids and bed elevation defined at
          vertices.

        Postcondition
          Conserved quantities defined at vertices
        """
        # Remove very thin layers of water
        domain.protect_against_infinitesimal_and_negative_heights()

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
        #balance_deep_and_shallow(domain)

        # Compute edge values by interpolation
        for name in domain.conserved_quantities:
            Q = domain.quantities[name]
            Q.interpolate_from_vertices_to_edges()


    def distribute_using_edge_limiter(domain):
        """Distribution from centroids to edges specific to the SWW eqn.

        Precondition:
          All quantities defined at centroids and bed elevation defined at
          vertices.

        Postcondition
          Conserved quantities defined at vertices
        """
        #from swb2_domain import protect_against_infinitesimal_and_negative_heights

        # Remove very thin layers of water
        #print 'Before protect'
        domain.protect_against_infinitesimal_and_negative_heights()
        #print 'After protect'

        #for name in domain.conserved_quantities:
        #    Q = domain.quantities[name]
        #    if domain._order_ == 1:
        #        Q.extrapolate_first_order()
        #    elif domain._order_ == 2:
        #        #Q.extrapolate_second_order_and_limit_by_edge()
        #        # FIXME: This use of 'break' is hacky
        #        domain.extrapolate_second_order_edge_sw()
        #        break
        #    else:
        #        raise Exception('Unknown order')
        
        #print 'Before extrapolate'
        domain.extrapolate_second_order_edge_sw()
        #print 'After extrapolate'

        #balance_deep_and_shallow(domain)

        # Compute vertex values by interpolation
        #for name in domain.conserved_quantities:
        #    Q = domain.quantities[name]
        #    Q.interpolate_from_edges_to_vertices()
        #    Q.interpolate_from_vertices_to_edges()


    def update_other_quantities(self):
        """ There may be a need to calculates some of the other quantities
        based on the new values of conserved quantities

        But that is not needed in this version of the solver.
        """
        pass
        # The centroid values of height and x and y velocity
        # might not have been setup
        
        #self.update_centroids_of_velocities_and_height()
        #
        #for name in ['height', 'xvelocity', 'yvelocity']:
        #    Q = self.quantities[name]
        #    if self._order_ == 1:
        #        Q.extrapolate_first_order()
        #    elif self._order_ == 2:
        #        Q.extrapolate_second_order_and_limit_by_edge()
        #    else:
        #        raise Exception('Unknown order')


    def update_centroids_of_velocities_and_height(self):
        """Calculate the centroid values of velocities and height based
        on the values of the quantities stage and x and y momentum

        Assumes that stage and momentum are up to date

        Useful for kinematic viscosity calculations
        """
        pass
        ## For shallow water we need to update height xvelocity and yvelocity

        ##Shortcuts
        #W  = self.quantities['stage']
        #UH = self.quantities['xmomentum']
        #VH = self.quantities['ymomentum']
        #H  = self.quantities['height']
        #Z  = self.quantities['elevation']
        #U  = self.quantities['xvelocity']
        #V  = self.quantities['yvelocity']

        ##print num.min(W.centroid_values)

        ## Make sure boundary values of conserved quantites
        ## are consistent with value of functions at centroids
        ##self.distribute_to_vertices_and_edges()
        #Z.set_boundary_values_from_edges()

        ##W.set_boundary_values_from_edges()
        ##UH.set_boundary_values_from_edges()
        ##VH.set_boundary_values_from_edges()

        ## Update height values
        #H.set_values(W.centroid_values-Z.centroid_values, location='centroids')
        #H.set_boundary_values( num.where(W.boundary_values-Z.boundary_values>=0,
        #                                 W.boundary_values-Z.boundary_values, 0.0))

        ##assert num.min(H.centroid_values) >= 0
        ##assert num.min(H.boundary_values) >= 0

        ##Aliases
        #uh_C  = UH.centroid_values
        #vh_C  = VH.centroid_values
        #h_C   = H.centroid_values

        #uh_B  = UH.boundary_values
        #vh_B  = VH.boundary_values
        #h_B   = H.boundary_values

        #H0 = 1.0e-8
        #
        #U.set_values(uh_C/(h_C + H0/h_C), location='centroids')
        #V.set_values(vh_C/(h_C + H0/h_C), location='centroids')

        #U.set_boundary_values(uh_B/(h_B + H0/h_B))
        #V.set_boundary_values(vh_B/(h_B + H0/h_B))



    def extrapolate_second_order_edge_sw(self):
        """Wrapper calling C version of extrapolate_second_order_sw, using
            edges instead of vertices in the extrapolation step. The routine
            then interpolates from edges to vertices. 
            The momentum extrapolation / interpolation is based on either 
            momentum or velocity, depending on the choice of 
            extrapolate_velocity_second_order.

        self  the domain to operate on

        """

        from .swb2_domain_ext import extrapolate_second_order_edge_sw as extrapol2

        # Shortcuts
        Stage = self.quantities['stage']
        Xmom = self.quantities['xmomentum']
        Ymom = self.quantities['ymomentum']
        Elevation = self.quantities['elevation']

        extrapol2(self,
                  self.surrogate_neighbours,
                  self.number_of_boundaries,
                  self.centroid_coordinates,
                  Stage.centroid_values,
                  Xmom.centroid_values,
                  Ymom.centroid_values,
                  Elevation.centroid_values,
                  self.edge_coordinates,
                  Stage.edge_values,
                  Xmom.edge_values,
                  Ymom.edge_values,
                  Elevation.edge_values,
                  Stage.vertex_values,
                  Xmom.vertex_values,
                  Ymom.vertex_values,
                  Elevation.vertex_values,
                  int(self.optimise_dry_cells),
                  int(self.extrapolate_velocity_second_order))

    def set_boundary(self, boundary_map):
        ## Specialisation of set_boundary, which also updates the 'boundary_flux_type'
        Sww_domain.set_boundary(self,boundary_map)
        
        # Add a flag which can be used to distinguish flux boundaries within
        # compute_fluxes_central
        # Initialise to zero (which means 'not a flux_boundary')
        self.boundary_flux_type = self.boundary_edges*0
       
        # HACK to set the values of domain.boundary_flux 
        for i in range(len(self.boundary_objects)):
            # Record the first 10 characters of the name of the boundary.
            # FIXME: There must be a better way than this!
            bndry_name = self.boundary_objects[i][1].__repr__()[0:10]

            if(bndry_name=='zero_mass_'):
                # Create flag 'boundary_flux_type' that can be read in
                # compute_fluxes_central
                self.boundary_flux_type[i]=1



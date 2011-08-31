
from anuga.shallow_water.shallow_water_domain import *
from anuga.shallow_water.shallow_water_domain import Domain as Sww_domain


##############################################################################
# Shallow Water Balanced Domain
#
# Uses extra evolved quantities height, elevation, xvelocity, yvelocity
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

        evolved_quantities = [ 'stage', 'xmomentum', 'ymomentum', \
                               'height', 'elevation', \
                               'xvelocity', 'yvelocity' ]
        
        other_quantities = [ 'friction', 'x', 'y' ]

        
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
        
        #---------------------
        # set some defaults
        #---------------------
        self.set_timestepping_method(2)
        self.set_default_order(2)
        self.set_sloped_mannings_function(True)
        self.set_centroid_transmissive_bc(True)
        self.set_CFL(1.0)
        self.set_beta(1.0)
        self.quantities['height'].set_beta(1.0)


    def check_integrity(self):
        Sww_domain.check_integrity(self)

        #Check that the evolved quantities are correct (order)
        msg = 'First evolved quantity must be "stage"'
        assert self.evolved_quantities[0] == 'stage', msg
        msg = 'Second evolved quantity must be "xmomentum"'
        assert self.evolved_quantities[1] == 'xmomentum', msg
        msg = 'Third evolved quantity must be "ymomentum"'
        assert self.evolved_quantities[2] == 'ymomentum', msg
        msg = 'Fourth evolved quantity must be "height"'
        assert self.evolved_quantities[3] == 'height', msg
        msg = 'Fifth evolved quantity must be "elevation"'
        assert self.evolved_quantities[4] == 'elevation', msg
        msg = 'Sixth evolved quantity must be "xvelocity"'
        assert self.evolved_quantities[5] == 'xvelocity', msg        
        msg = 'Seventh evolved quantity must be "yvelocity"'
        assert self.evolved_quantities[6] == 'yvelocity', msg        

        msg = 'First other quantity must be "friction"'
        assert self.other_quantities[0] == 'friction', msg
        msg = 'Second other quantity must be "x"'
        assert self.other_quantities[1] == 'x', msg
        msg = 'Third other quantity must be "y"'
        assert self.other_quantities[2] == 'y', msg


    def compute_fluxes(self):
        """
        Call correct module function (either from this module or C-extension)
        """

        from swb_domain_ext import compute_fluxes_c

        #Shortcuts
        W  = self.quantities['stage']
        UH = self.quantities['xmomentum']
        VH = self.quantities['ymomentum']
        H  = self.quantities['height']
        Z  = self.quantities['elevation']
        U  = self.quantities['xvelocity']
        V  = self.quantities['yvelocity']

        timestep = self.get_evolve_max_timestep()
        
        self.flux_timestep = \
            compute_fluxes_c(timestep, self, W, UH, VH, H, Z, U, V)


    def distribute_to_vertices_and_edges(self):
        """Distribution from centroids to edges specific to the SWW eqn.

        It will ensure that h (w-z) is always non-negative even in the
        presence of steep bed-slopes by taking a weighted average between shallow
        and deep cases.

        In addition, all conserved quantities get distributed as per either a
        constant (order==1) or a piecewise linear function (order==2).
        

        Precondition:
        All conserved quantities defined at centroids and bed elevation defined at
        edges.
        
        Postcondition
        Evolved quantities defined at vertices and edges
        """

        from anuga.shallow_water.shallow_water_domain import \
                  protect_against_infinitesimal_and_negative_heights as protect
    
        #Shortcuts
        W  = self.quantities['stage']
        UH = self.quantities['xmomentum']
        VH = self.quantities['ymomentum']
        H  = self.quantities['height']
        Z  = self.quantities['elevation']
        U  = self.quantities['xvelocity']
        V  = self.quantities['yvelocity']

        #Arrays   
        w_C   = W.centroid_values    
        uh_C  = UH.centroid_values
        vh_C  = VH.centroid_values    
        z_C   = Z.centroid_values
        h_C   = H.centroid_values
        u_C   = U.centroid_values
        v_C   = V.centroid_values

        num_min = num.min(w_C-z_C)
        if num_min < -1.0e-5:
            print '**** num.min(w_C-z_C)', num_min

    
        w_C[:] = num.maximum(w_C, z_C)
        h_C[:] = w_C - z_C


        assert num.min(h_C) >= 0.0
                
        num.putmask(uh_C, h_C < 1.0e-15, 0.0)
        num.putmask(vh_C, h_C < 1.0e-15, 0.0)
        #num.putmask(h_C, h_C < 1.0e-15, 1.0e-16)

        # Noelle has an alternative method for calculating velcities
        # Check it out in the GPU shallow water paper
        H0 = 1.0e-16
        u_C[:]  = uh_C/(h_C + H0/h_C)
        v_C[:]  = vh_C/(h_C + H0/h_C)

        #num.putmask(h_C, h_C < 1.0e-15, 0.0)
        
        for name in [ 'stage', 'xvelocity', 'yvelocity' ]:
            Q = self.quantities[name]
            if self._order_ == 1:
                Q.extrapolate_first_order()
            elif self._order_ == 2:
                Q.extrapolate_second_order_and_limit_by_edge()
                #Q.extrapolate_second_order_and_limit_by_vertex()
            else:
                raise Exception('Unknown order: %s' % str(self._order_))

        for name in [ 'height' ]:
            Q = self.quantities[name]
            if self._order_ == 1:
                Q.extrapolate_first_order()
            elif self._order_ == 2:
                #Q.extrapolate_second_order_and_limit_by_edge()
                Q.extrapolate_second_order_and_limit_by_vertex()
            else:
                raise Exception('Unknown order: %s' % str(self._order_))


        w_V     = W.vertex_values 
        u_V     = U.vertex_values
        v_V     = V.vertex_values
        z_V     = Z.vertex_values

        h_V     = H.vertex_values
        uh_V    = UH.vertex_values
        vh_V    = VH.vertex_values
        

        # Update other quantities
        #protect(self)

        z_V[:]  = w_V - h_V
        uh_V[:] = u_V * h_V
        vh_V[:] = v_V * h_V

        
        num_min = num.min(h_V)
        if num_min < -1.0e-14:
            print 'num.min(h_V)', num_min

        
        # Compute edge values by interpolation
        for name in ['xmomentum', 'ymomentum', 'elevation']:
            Q = self.quantities[name]
            Q.interpolate_from_vertices_to_edges()





    def distribute_to_vertices_and_edges_h(self):
        """Distribution from centroids to edges specific to the SWW eqn.

        It will ensure that h (w-z) is always non-negative even in the
        presence of steep bed-slopes by taking a weighted average between shallow
        and deep cases.

        In addition, all conserved quantities get distributed as per either a
        constant (order==1) or a piecewise linear function (order==2).
        

        Precondition:
        All conserved quantities defined at centroids and bed elevation defined at
        edges.
        
        Postcondition
        Evolved quantities defined at vertices and edges
        """


        #Shortcuts
        W  = self.quantities['stage']
        UH = self.quantities['xmomentum']
        VH = self.quantities['ymomentum']
        H  = self.quantities['height']
        Z  = self.quantities['elevation']
        U  = self.quantities['xvelocity']
        V  = self.quantities['yvelocity']

        #Arrays   
        w_C   = W.centroid_values    
        uh_C  = UH.centroid_values
        vh_C  = VH.centroid_values    
        z_C   = Z.centroid_values
        h_C   = H.centroid_values
        u_C   = U.centroid_values
        v_C   = V.centroid_values

        w_C[:] = num.maximum(w_C, z_C)
        
        h_C[:] = w_C - z_C


        assert num.min(h_C) >= 0
                
        num.putmask(uh_C, h_C < 1.0e-15, 0.0)
        num.putmask(vh_C, h_C < 1.0e-15, 0.0)
        num.putmask(h_C, h_C < 1.0e-15, 1.0e-16)        
        
        u_C[:]  = uh_C/h_C
        v_C[:]  = vh_C/h_C

        num.putmask(h_C, h_C < 1.0e-15, 0.0)
        
        for name in [ 'stage', 'height', 'xvelocity', 'yvelocity' ]:
            Q = self.quantities[name]
            if self._order_ == 1:
                Q.extrapolate_first_order()
            elif self._order_ == 2:
                Q.extrapolate_second_order_and_limit_by_edge()
                #Q.extrapolate_second_order_and_limit_by_vertex()
            else:
                raise Exception('Unknown order: %s' % str(self._order_))


        w_E     = W.edge_values
        uh_E    = UH.edge_values
        vh_E    = VH.edge_values	
        h_E     = H.edge_values
        z_E     = Z.edge_values	
        u_E     = U.edge_values
        v_E     = V.edge_values		


        #minh_E = num.min(h_E)
        #msg = 'min h_E = %g ' % minh_E
        #assert minh_E >= -1.0e-15, msg

        z_E[:]   = w_E - h_E

        num.putmask(h_E, h_E <= 1.0e-8, 0.0)
        num.putmask(u_E, h_E <= 1.0e-8, 0.0)
        num.putmask(v_E, h_E <= 1.0e-8, 0.0)
        num.putmask(w_E, h_E <= 1.0e-8, z_E)
        #num.putmask(h_E, h_E <= 0.0, 0.0)
        
        uh_E[:] = u_E * h_E
        vh_E[:] = v_E * h_E

        """
        print '=========================================================='
        print 'Time ', self.get_time()
        print h_E
        print uh_E
        print vh_E
        """
        
        # Compute vertex values by interpolation
        for name in self.evolved_quantities:
            Q = self.quantities[name]
            Q.interpolate_from_edges_to_vertices()


        w_V     = W.vertex_values
        uh_V    = UH.vertex_values
        vh_V    = VH.vertex_values	
        z_V     = Z.vertex_values	
        h_V     = H.vertex_values
        u_V     = U.vertex_values
        v_V     = V.vertex_values		


        #w_V[:]    = z_V + h_V

        #num.putmask(u_V, h_V <= 0.0, 0.0)
        #num.putmask(v_V, h_V <= 0.0, 0.0)
        #num.putmask(w_V, h_V <= 0.0, z_V)        
        #num.putmask(h_V, h_V <= 0.0, 0.0)
        
        uh_V[:] = u_V * h_V
        vh_V[:] = v_V * h_V




    def conserved_values_to_evolved_values(self, q_cons, q_evol):
        """Mapping between conserved quantities and the evolved quantities.
        Used where we have a boundary condition which works with conserved
        quantities and we now want to use them for the new well balanced
        code using the evolved quantities

        Typically the old boundary condition will set the values in q_cons,

        q_evol on input will have the values of the evolved quantities at the
        edge point (useful as it provides values for evlevation).
        """

        wc  = q_cons[0]
        uhc = q_cons[1]
        vhc = q_cons[2]

        we  = q_evol[0]
        uhe = q_evol[1]
        vhe = q_evol[2]

        he  = q_evol[3]
        be  = q_evol[4]
        ue  = q_evol[5]
        ve  = q_evol[6]


        hc = wc - be

        if hc <= 0.0:
            hc = 0.0
            uc = 0.0
            vc = 0.0
        else:
            uc = uhc/hc
            vc = vhc/hc

        q_evol[0]  = wc
        q_evol[1]  = uhc
        q_evol[2]  = vhc

        q_evol[3]  = hc
        q_evol[4]  = be
        q_evol[5]  = uc
        q_evol[6]  = vc


        return q_evol

 
################################################################################
# Standard forcing terms
################################################################################

def gravity(domain):
    """Apply gravitational pull in the presence of bed slope
    Wrapper calls underlying C implementation
    """

    from swb_domain_ext import gravity as gravity_c

    xmom_update = domain.quantities['xmomentum'].explicit_update
    ymom_update = domain.quantities['ymomentum'].explicit_update

    stage = domain.quantities['stage']
    elevation = domain.quantities['elevation']


    stage = stage.vertex_values
    elevation = elevation.vertex_values

    points = domain.get_vertex_coordinates()

    gravity_c(domain.g, stage, elevation, points, xmom_update, ymom_update)
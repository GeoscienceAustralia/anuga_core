
from anuga.shallow_water.shallow_water_domain import *
from anuga.shallow_water.shallow_water_domain import Domain as Sww_domain

from swb_boundary_conditions import Transmissive_boundary

##############################################################################
# Shallow Water Balanced Domain
#
# Uses extra evolved quantities height, elevation, xvelocity, yvelocity
##############################################################################

##
# @brief Class for a shallow water balanced domain.
class Domain(Sww_domain):

    ##
    # @brief Instantiate a shallow water balanced domain.
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

        conserved_quantities = [ 'stage', 'xmomentum', 'ymomentum']

        evolved_quantities = [ 'stage', 'xmomentum', 'ymomentum', \
                               'height', 'elevation', 'xvelocity', 'yvelocity']
        
        other_quantities = [ 'friction' ]


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
                            processor = processor,
                            numproc = numproc,
                            number_of_full_nodes = number_of_full_nodes,
                            number_of_full_triangles = number_of_full_triangles)
        
        #---------------------
        # set some defaults
        #---------------------
        self.set_timestepping_method(1)
        self.set_default_order(2)
        self.set_new_mannings_function(True)
        self.set_centroid_transmissive_bc(True)

    ##
    # @brief Run integrity checks on shallow water balanced domain.
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

    ##
    # @brief 
    def compute_fluxes(self):
        #Call correct module function (either from this module or C-extension)
        compute_fluxes(self)

    ##
    # @brief 
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


        #Shortcuts
        Stage  = self.quantities['stage']
        Xmom   = self.quantities['xmomentum']
        Ymom   = self.quantities['ymomentum']
        Elev   = self.quantities['elevation']
        Height = self.quantities['height']
        Xvel   = self.quantities['xvelocity']
        Yvel   = self.quantities['yvelocity']

        #Arrays   
        w_C   = Stage.centroid_values    
        uh_C  = Xmom.centroid_values
        vh_C  = Ymom.centroid_values    
        z_C   = Elev.centroid_values
        h_C   = Height.centroid_values
        u_C   = Xvel.centroid_values
        v_C   = Yvel.centroid_values

        w_C[:] = num.maximum(w_C, z_C)
        
        h_C[:] = w_C - z_C


        assert num.min(h_C) >= 0
                
        num.putmask(uh_C, h_C < 1.0e-15, 0.0)
        num.putmask(vh_C, h_C < 1.0e-15, 0.0)
        num.putmask(h_C, h_C < 1.0e-15, 1.0e-15)        
        
        u_C[:]  = uh_C/h_C
        v_C[:]  = vh_C/h_C
	
        for name in [ 'height', 'xvelocity', 'yvelocity' ]:
            Q = self.quantities[name]
            if self._order_ == 1:
                Q.extrapolate_first_order()
            elif self._order_ == 2:
                Q.extrapolate_second_order_and_limit_by_edge()
            else:
                raise 'Unknown order'


        w_E     = Stage.edge_values
        uh_E    = Xmom.edge_values
        vh_E    = Ymom.edge_values	
        z_E     = Elev.edge_values	
        h_E     = Height.edge_values
        u_E     = Xvel.edge_values
        v_E     = Yvel.edge_values		


        w_E[:]   = z_E + h_E

        #num.putmask(u_E, h_temp <= 0.0, 0.0)
        #num.putmask(v_E, h_temp <= 0.0, 0.0)
        #num.putmask(w_E, h_temp <= 0.0, z_E+h_E)
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


        w_V     = Stage.vertex_values
        uh_V    = Xmom.vertex_values
        vh_V    = Ymom.vertex_values	
        z_V     = Elev.vertex_values	
        h_V     = Height.vertex_values
        u_V     = Xvel.vertex_values
        v_V     = Yvel.vertex_values		


        #w_V[:]    = z_V + h_V

        #num.putmask(u_V, h_V <= 0.0, 0.0)
        #num.putmask(v_V, h_V <= 0.0, 0.0)
        #num.putmask(w_V, h_V <= 0.0, z_V)        
        #num.putmask(h_V, h_V <= 0.0, 0.0)
        
        uh_V[:] = u_V * h_V
        vh_V[:] = v_V * h_V


    ##
    # @brief Code to let us use old shallow water domain BCs
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

 

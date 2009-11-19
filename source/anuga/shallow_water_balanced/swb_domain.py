
from anuga.shallow_water.shallow_water_domain import *
from anuga.shallow_water.shallow_water_domain import Domain as Sww_domain


##############################################################################
# Shallow Water Balanced Domain
#
# Uses extra evolved quantities height, elevation, xvelocity, yvelocity
##############################################################################

##
# @brief Class for a shallow water balanced domain.
class Domain(Sww_domain):

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
        self.set_timestepping_method('euler')
        self.set_default_order(1)
        self.set_new_mannings_function(True)
        self.set_use_edge_limiter(True)



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

        if hc < 0.0:
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

 

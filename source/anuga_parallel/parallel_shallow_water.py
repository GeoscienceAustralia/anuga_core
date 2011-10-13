"""Class Parallel_shallow_water_domain -
2D triangular domains for finite-volume computations of
the shallow water equation, with extra structures to allow
communication between other Parallel_domains and itself

This module contains a specialisation of class Domain
from module shallow_water.py

Ole Nielsen, Stephen Roberts, Duncan Gray, Christopher Zoppou
Geoscience Australia, 2004-2005

"""

from anuga import Domain

from anuga_parallel.parallel_generic_communications import *
from anuga.abstract_2d_finite_volumes.neighbour_mesh import Mesh

import numpy as num



class Parallel_domain(Domain):

    def __init__(self, coordinates, vertices,
                 boundary=None,
                 full_send_dict=None,
                 ghost_recv_dict=None,
                 number_of_full_nodes=None,
                 number_of_full_triangles=None,
                 geo_reference=None): #jj added this

        Domain.__init__(self,
                        coordinates,
                        vertices,
                        boundary,
                        full_send_dict=full_send_dict,
                        ghost_recv_dict=ghost_recv_dict,
                        processor=pypar.rank(),
                        numproc=pypar.size(),
                        number_of_full_nodes=number_of_full_nodes,
                        number_of_full_triangles=number_of_full_triangles,
                        geo_reference=geo_reference) #jj added this
        
        # PETE: Find the number of full nodes and full triangles, this is a temporary fix
        # until the bug with get_number_of_full_[nodes|triangles]() is fixed.

        if number_of_full_nodes is not None:
            self.number_of_full_nodes_tmp = number_of_full_nodes
        else:
            self.number_of_full_nodes_tmp = get_number_of_nodes()

        if number_of_full_triangles is not None:
            self.number_of_full_triangles_tmp = number_of_full_triangles
        else:
            self.number_of_full_triangles_tmp = get_number_of_triangles()

        setup_buffers(self)


    def set_name(self, name):
        """Assign name based on processor number 
        """

        if name.endswith('.sww'):
            name = name[:-4]

        # Call parents method with processor number attached.
        Domain.set_name(self, name + '_P%d_%d' %(self.processor, self.numproc))


    def update_timestep(self, yieldstep, finaltime):
        """Calculate local timestep
        """

        communicate_flux_timestep(self, yieldstep, finaltime)

        Domain.update_timestep(self, yieldstep, finaltime)



    def update_ghosts(self):

        # We must send the information from the full cells and
        # receive the information for the ghost cells


        communicate_ghosts(self)

    def apply_fractional_steps(self):

        for operator in self.fractional_step_operators:
            operator()

        # PETE: Make sure that there are no deadlocks here

        self.update_ghosts()

# =======================================================================
# PETE: NEW METHODS FOR FOR PARALLEL STRUCTURES. Note that we assume the 
# first "number_of_full_[nodes|triangles]" are full [nodes|triangles]
# For full triangles it is possible to enquire self.tri_full_flag == True
# =======================================================================

    def get_number_of_full_triangles(self, *args, **kwargs):
        return self.number_of_full_triangles_tmp

    def get_full_centroid_coordinates(self, *args, **kwargs):
        C = self.mesh.get_centroid_coordinates(*args, **kwargs)
        return C[:self.number_of_full_triangles_tmp, :]

    def get_full_vertex_coordinates(self, *args, **kwargs):
        V = self.mesh.get_vertex_coordinates(*args, **kwargs)
        return V[:3*self.number_of_full_triangles_tmp,:]

    def get_full_triangles(self, *args, **kwargs):
        T = self.mesh.get_triangles(*args, **kwargs)
        return T[:self.number_of_full_triangles_tmp,:]

    def get_full_nodes(self, *args, **kwargs):
        N = self.mesh.get_nodes(*args, **kwargs)
        return N[:self.number_of_full_nodes_tmp,:]


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

import numpy as num



class Parallel_domain(Domain):

    def __init__(self, coordinates, vertices,
                 boundary=None,
                 full_send_dict=None,
                 ghost_recv_dict=None,
                 number_of_full_nodes=None,
                 number_of_full_triangles=None):

        Domain.__init__(self,
                        coordinates,
                        vertices,
                        boundary,
                        full_send_dict=full_send_dict,
                        ghost_recv_dict=ghost_recv_dict,
                        processor=pypar.rank(),
                        numproc=pypar.size(),
                        number_of_full_nodes=number_of_full_nodes,
                        number_of_full_triangles=number_of_full_triangles)

 
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

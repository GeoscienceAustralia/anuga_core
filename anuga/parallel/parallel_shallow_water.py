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

from . import parallel_generic_communications as generic_comms

import anuga.utilities.parallel_abstraction as pypar

#from anuga.abstract_2d_finite_volumes.neighbour_mesh import Mesh

import numpy as num
from os.path import join


#Import matplotlib



class Parallel_domain(Domain):

    def __init__(self, coordinates, vertices,
                 boundary=None,
                 full_send_dict=None,
                 ghost_recv_dict=None,
                 number_of_full_nodes=None,
                 number_of_full_triangles=None,
                 geo_reference=None,
                 processor = None,
                 numproc = None,
                 number_of_global_triangles=None, ## SR added this
                 number_of_global_nodes= None, ## SR added this
                 s2p_map=None,
                 p2s_map=None, #jj added this
                 tri_l2g = None, ## SR added this
                 node_l2g = None, #): ## SR added this
                 ghost_layer_width = 2): ## SR added this



        #-----------------------------------------
        # Sometimes we want to manually
        # create instances of the parallel_domain
        # otherwise ...
        #----------------------------------------
        if processor is None:
            processor = pypar.rank()
        if numproc is None:
            numproc = pypar.size()

        Domain.__init__(self,
                        coordinates,
                        vertices,
                        boundary,
                        full_send_dict=full_send_dict,
                        ghost_recv_dict=ghost_recv_dict,
                        processor=processor,
                        numproc=numproc,
                        number_of_full_nodes=number_of_full_nodes,
                        number_of_full_triangles=number_of_full_triangles,
                        geo_reference=geo_reference, #) #jj added this
                        ghost_layer_width = ghost_layer_width)


        self.parallel = True

        # PETE: Find the number of full nodes and full triangles, this is a temporary fix
        # until the bug with get_number_of_full_[nodes|triangles]() is fixed.

        if number_of_full_nodes is not None:
            self.number_of_full_nodes_tmp = number_of_full_nodes
        else:
            self.number_of_full_nodes_tmp = self.get_number_of_nodes()

        if number_of_full_triangles is not None:
            self.number_of_full_triangles_tmp = number_of_full_triangles
        else:
            self.number_of_full_triangles_tmp = self.get_number_of_triangles()

        generic_comms.setup_buffers(self)

        self.global_name = 'domain'

        self.number_of_global_triangles=number_of_global_triangles
        self.number_of_global_nodes = number_of_global_nodes

        self.s2p_map = s2p_map
        self.p2s_map = p2s_map


        self.s2p_map = None
        self.p2s_map = None

        self.tri_l2g = tri_l2g
        self.node_l2g = node_l2g

        self.ghost_counter = 0


    def set_name(self, name):
        """Assign name based on processor number
        """

        if name.endswith('.sww'):
            name = name[:-4]

        self.global_name = name

        # Call parents method with processor number attached.
        Domain.set_name(self, name + '_P%d_%d' %(self.numproc, self.processor))


    def get_global_name(self):

        return self.global_name


    def update_timestep(self, yieldstep, finaltime):
        """Calculate local timestep
        """

        generic_comms.communicate_flux_timestep(self, yieldstep, finaltime)

        Domain.update_timestep(self, yieldstep, finaltime)



    def update_ghosts(self, quantities=None):
        """We must send the information from the full cells and
        receive the information for the ghost cells
        """

        #generic_comms.communicate_ghosts_asynchronous(self, quantities)
        generic_comms.communicate_ghosts_non_blocking(self, quantities)
        #generic_comms.communicate_ghosts_blocking(self)

    def apply_fractional_steps(self):

        for operator in self.fractional_step_operators:
            operator()

        # PETE: Make sure that there are no deadlocks here

        #self.update_ghosts()



    def sww_merge(self, verbose=False, delete_old=False):
        """Merge all the sub domain sww files into a global sww file

        :param bool verbose: Flag to produce more output
        :param bool delete_old: Flag to delete sub domain sww files after
            creating global sww file

        """

        # make sure all the computations have finished

        pypar.barrier()

        # now on processor 0 pull all the separate sww files together
        if self.processor == 0 and self.numproc > 1 and self.store :
            import anuga.utilities.sww_merge as merge

            global_name = join(self.get_datadir(),self.get_global_name())

            merge.sww_merge_parallel(global_name,self.numproc,verbose,delete_old)

        # make sure all the merge completes on processor 0 before other
        # processors complete (like when finalize is forgotten in main script)

        pypar.barrier()

    def write_time(self):

        if self.processor == 0:
            Domain.write_time(self)


    def dump_triangulation(self, filename="domain.png"):
        """
        Outputs domain triangulation, full triangles are shown in green while ghost triangles are shown in blue.
        The default filename is 'domain.png'
        """

        # Get vertex coordinates, partition full and ghost triangles based on self.tri_full_flag

        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            import matplotlib.tri as tri
        except:
            print("Couldn't import module from matplotlib, probably you need to update matplotlib")
            raise

        vertices = self.get_vertex_coordinates()
        full_mask = num.repeat(self.tri_full_flag == 1, 3)
        ghost_mask = num.repeat(self.tri_full_flag == 0, 3)

        myid = pypar.rank()
        numprocs = pypar.size()

        if myid == 0:

            fig = plt.figure()
            fx = {}
            fy = {}
            gx = {}
            gy = {}

            # Proc 0 gathers full and ghost nodes from self and other processors
            fx[0] = vertices[full_mask,0]
            fy[0] = vertices[full_mask,1]
            gx[0] = vertices[ghost_mask,0]
            gy[0] = vertices[ghost_mask,1]

            for i in range(1,numprocs):
                fx[i] = pypar.receive(i)
                fy[i] = pypar.receive(i)
                gx[i] = pypar.receive(i)
                gy[i] = pypar.receive(i)

            # Plot full triangles
            for i in range(0, numprocs):
                n = int(len(fx[i])//3)

                triang = num.array(list(range(0,3*n)))
                triang.shape = (n, 3)
                plt.triplot(fx[i], fy[i], triang, 'g-', linewidth = 0.5)

            # Plot ghost triangles
            for i in range(0, numprocs):
                n = int(len(gx[i])//3)
                if n > 0:
                    triang = num.array(list(range(0,3*n)))
                    triang.shape = (n, 3)
                    plt.triplot(gx[i], gy[i], triang, 'b--', linewidth = 0.5)

            # Save triangulation to location pointed by filename
            plt.savefig(filename, dpi=600)

        else:
            # Proc 1..numprocs send full and ghost triangles to Proc 0
            pypar.send(vertices[full_mask,0], 0)
            pypar.send(vertices[full_mask,1], 0)
            pypar.send(vertices[ghost_mask,0], 0)
            pypar.send(vertices[ghost_mask,1], 0)


    def dump_local_triangulation(self, filename=None):
        '''
        Outputs domain triangulation, full triangles are shown in green while
        ghost triangles are shown in blue.

        The default filename is self.get_name()+'.png'
        '''

        if filename is None:
            filename = self.get_name() + '.png'

        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            import matplotlib.tri as tri
        except:
            print("Couldn't import module from matplotlib, probably you need to update matplotlib")
            raise

        vertices = self.get_vertex_coordinates()
        full_mask = num.repeat(self.tri_full_flag == 1, 3)
        ghost_mask = num.repeat(self.tri_full_flag == 0, 3)

        plt.figure()

        fx = vertices[full_mask,0]
        fy = vertices[full_mask,1]
        gx = vertices[ghost_mask,0]
        gy = vertices[ghost_mask,1]

        # Plot full triangles
        n = int(len(fx)/3)
        triang = num.array(list(range(0,3*n)))
        triang.shape = (n, 3)
        plt.triplot(fx, fy, triang, 'g-', linewidth = 0.5)

        # Plot ghost triangles
        n = int(len(gx)/3)
        if n > 0:
            triang = num.array(list(range(0,3*n)))
            triang.shape = (n, 3)
            plt.triplot(gx, gy, triang, 'b--', linewidth = 0.5)

        # Save triangulation to location pointed by filename
        plt.savefig(filename, dpi=600)


import sys
from os import sep

import unittest
from math import sqrt, pi

from anuga.config import g, epsilon
from anuga.abstract_2d_finite_volumes.generic_domain import Generic_Domain
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions import \
                Transmissive_boundary, Dirichlet_boundary

from anuga.advection.advection import Advection_Domain

import numpy as num


class Test_Advection(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_init(self):
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe, daf, dae
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Advection_Domain(points, vertices)
        domain.check_integrity()

        assert 'stage' in domain.quantities

        assert domain.get_conserved_quantities(0, edge=1) == 0.


    def test_flux_1_triangle0(self):
        a = [0.0, 0.5]
        b = [0.0, 0.0]
        c = [0.5, 0.5]

        points = [a, b, c]
        vertices = [ [0,1,2] ]
        domain = Advection_Domain(points, vertices)
        domain.check_integrity()


        #Populate boundary array with dirichlet conditions.
        domain.neighbours = num.array([[-1,-2,-3]])
        domain.quantities['stage'].boundary_values[:] = 1.0

        domain.order = 1

        domain.distribute_to_vertices_and_edges() #Use first order default

        domain.check_integrity()

        domain.compute_fluxes()
        U = -domain.quantities['stage'].explicit_update
        R = -0.5/domain.areas[0]

        assert U==R, '%s %s' %(U, R)


    def test_flux_1_triangle1(self):

        a = [0.0, 0.5]
        b = [0.0, 0.0]
        c = [0.5, 0.5]

        points = [a, b, c]
        vertices = [ [0,1,2] ]
        domain = Advection_Domain(points, vertices)
        domain.check_integrity()

        domain.set_quantity('stage', [1.0], location='centroids')

        domain.distribute_to_vertices_and_edges()
        domain.check_integrity()


        domain.compute_fluxes()
        U = -domain.quantities['stage'].explicit_update
        R = 0.5/domain.areas[0]

        assert U==R, '%s %s' %(U, R)



    def test_flux_1_triangle2(self):

        a = [0.0, 0.5]
        b = [0.0, 0.0]
        c = [0.5, 0.5]

        points = [a, b, c]
        vertices = [ [0,1,2] ]
        domain = Advection_Domain(points, vertices)
        domain.check_integrity()


        #Populate boundary array with dirichlet conditions.
        domain.neighbours = num.array([[-1,-2,-3]])
        domain.quantities['stage'].boundary_values[0] = 1.0

        domain.distribute_to_vertices_and_edges() #Use first order default

        domain.check_integrity()

        domain.compute_fluxes()
        U = domain.quantities['stage'].explicit_update
        assert num.allclose(U, 0)




    def test_flux_2_triangles(self):
        """Flow between two triangles
        Check that fluxes have opposite signs
        """

        a = [0.0, 0.5]
        b = [0.0, 0.0]
        c = [0.5, 0.5]
        d = [0.5, 0.0]

        points = [a, b, c, d]
        vertices = [ [0,1,2], [3,2,1] ]
        domain = Advection_Domain(points, vertices)
        domain.check_integrity()


        #Populate boundary array with dirichlet conditions.
        domain.neighbours = num.array([[1,-1,-2], [0,-3,-4]])
        domain.set_quantity('stage', [1.0, 0.0], location='centroids')
        domain.distribute_to_vertices_and_edges()

        domain.compute_fluxes()

        X = domain.quantities['stage'].explicit_update
        assert X[0] == -X[1]


    def test_advection_example(self):
        #Test that system can evolve

        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        points, vertices, boundary = rectangular(6, 6)

        #Create advection domain with direction (1,-1)
        domain = Advection_Domain(points, vertices, boundary,
                                        velocity=[1.0, -1.0])

        # Initial condition is zero by default

        #Boundaries
        T = Transmissive_boundary(domain)
        D = Dirichlet_boundary(num.array([3.1415]))

        domain.set_boundary( {'left': D, 'right': T, 'bottom': T, 'top': T} )
        domain.check_integrity()

        #Check that the boundary value gets propagated to all elements
        for t in domain.evolve(yieldstep = 0.05, finaltime = 10):
            if num.allclose(domain.quantities['stage'].centroid_values, 3.1415):
                break

        assert num.allclose(domain.quantities['stage'].centroid_values, 3.1415)


#-------------------------------------------------------------
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_Advection, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)

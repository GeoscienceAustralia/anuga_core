
import unittest
import anuga
import numpy
import os

boundaryPolygon=[ [0., 0.], [0., 100.], [100.0, 100.0], [100.0, 0.0]]

verbose=False

class Test_boundary_flux_integral_operator(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        try:
            os.remove('test_boundaryfluxintegral.sww')
        except OSError:
            pass

    def create_domain(self, flowalg):
        # Riverwall = list of lists, each with a set of x,y,z (and optional QFactor) values

        # Make the domain
        domain = anuga.create_domain_from_regions(boundaryPolygon,
                                 boundary_tags={'left': [0],
                                                'top': [1],
                                                'right': [2],
                                                'bottom': [3]},
                                   maximum_triangle_area = 200.,
                                   minimum_triangle_angle = 28.0,
                                   use_cache=False,
                                   verbose=verbose)


        # 05/05/2014 -- riverwalls only work with DE0 and DE1
        domain.set_flow_algorithm(flowalg)
        domain.set_name('test_boundaryfluxintegral')

        domain.set_store_vertices_uniquely()

        def topography(x,y):
            return -x/150.

        # NOTE: Setting quantities at centroids is important for exactness of tests
        domain.set_quantity('elevation',topography,location='centroids')
        domain.set_quantity('friction',0.03)
        domain.set_quantity('stage', topography,location='centroids')

        # Boundary conditions
        Br=anuga.Reflective_boundary(domain)
        Bd=anuga.Dirichlet_boundary([0., 0., 0.])
        domain.set_boundary({'left': Br, 'right': Bd, 'top': Br, 'bottom':Br})

        return domain

    def test_boundary_flux_operator_DE0(self):
        """
        A (the) boundary flux operator is instantiated when a domain is created.
        This tests the calculation for euler timestepping
        """

        flowalg = 'DE0'

        domain=self.create_domain(flowalg)

        #domain.print_statistics()
        for t in domain.evolve(yieldstep=1.0,finaltime=5.0):
            if verbose: domain.print_timestepping_statistics()
            if verbose: print(domain.get_water_volume())
            pass
        # The domain was initially dry
        vol=domain.get_water_volume()
        boundaryFluxInt=domain.get_boundary_flux_integral()

        if verbose: print(flowalg, vol, boundaryFluxInt)
        assert(numpy.allclose(vol,boundaryFluxInt))



    def test_boundary_flux_operator_DE1(self):
        """
        A (the) boundary flux operator is instantiated when a domain is created.
        This tests the calculation for rk2 timestepping
        """
        flowalg = 'DE1'

        domain=self.create_domain(flowalg)
        #domain.print_statistics()
        for t in domain.evolve(yieldstep=1.0,finaltime=5.0):
            if verbose: domain.print_timestepping_statistics()
            if verbose: print(domain.get_water_volume())
            pass
        # The domain was initially dry
        vol=domain.get_water_volume()
        boundaryFluxInt=domain.get_boundary_flux_integral()

        if verbose: print(flowalg, vol, boundaryFluxInt)
        assert(numpy.allclose(vol,boundaryFluxInt))



    def test_boundary_flux_operator_DE2(self):
        """
        A (the) boundary flux operator is instantiated when a domain is created.
        This tests the calculation for rk3 timestepping
        """

        flowalg = 'DE2'

        domain=self.create_domain(flowalg)
        #domain.print_statistics()
        for t in domain.evolve(yieldstep=1.0,finaltime=5.0):
            if verbose: domain.print_timestepping_statistics()
            if verbose: print(domain.get_water_volume(), domain.get_boundary_flux_integral())
            pass
        # The domain was initially dry
        vol=domain.get_water_volume()
        boundaryFluxInt=domain.get_boundary_flux_integral()

        if verbose: print(flowalg, vol, boundaryFluxInt)
        assert(numpy.allclose(vol,boundaryFluxInt))

    def test_boundary_flux_operator_DE_ader2(self):
        """
        A (the) boundary flux operator is instantiated when a domain is created.
        This tests the calculation for ader2 timestepping
        """

        flowalg = 'DE_ader2'

        domain=self.create_domain(flowalg)
        for t in domain.evolve(yieldstep=1.0,finaltime=5.0):
            if verbose: domain.print_timestepping_statistics()
            if verbose: print(domain.get_water_volume(), domain.get_boundary_flux_integral())
            pass
        # The domain was initially dry
        vol=domain.get_water_volume()
        boundaryFluxInt=domain.get_boundary_flux_integral()

        if verbose: print(flowalg, vol, boundaryFluxInt)
        assert(numpy.allclose(vol,boundaryFluxInt))


class Test_boundary_flux_integral_operator_extra(unittest.TestCase):
    """Tests for uncovered methods in boundary_flux_integral_operator."""

    def setUp(self):
        from anuga import rectangular_cross_domain, Reflective_boundary
        self.domain = rectangular_cross_domain(2, 2)
        self.domain.set_quantity('elevation', 0.0)
        self.domain.set_quantity('stage', 1.0)
        Br = anuga.Reflective_boundary(self.domain)
        self.domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

    def tearDown(self):
        try:
            import os
            os.remove('domain.sww')
        except OSError:
            pass

    def _make_op(self):
        from anuga.operators.boundary_flux_integral_operator import boundary_flux_integral_operator
        return boundary_flux_integral_operator(self.domain)

    def test_parallel_safe(self):
        op = self._make_op()
        self.assertTrue(op.parallel_safe())

    def test_statistics(self):
        op = self._make_op()
        msg = op.statistics()
        self.assertIsInstance(msg, str)

    def test_timestepping_statistics(self):
        op = self._make_op()
        msg = op.timestepping_statistics()
        self.assertIsInstance(msg, str)

    def test_call_unsupported_method_raises(self):
        """Unsupported timestepping method raises Exception (line 59)."""
        from anuga.operators.boundary_flux_integral_operator import boundary_flux_integral_operator
        op = boundary_flux_integral_operator(self.domain)
        self.domain.timestep = 1.0
        self.domain.timestepping_method = 'unsupported_method'
        with self.assertRaises(Exception):
            op()


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_boundary_flux_integral_operator)
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)


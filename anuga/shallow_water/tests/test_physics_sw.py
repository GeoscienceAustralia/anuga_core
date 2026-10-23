#!/usr/bin/env python
"""Tests for shallow water physics: gravity, friction, bed slope, flat/sloped beds."""
import os
import unittest
from math import pi, sqrt
import numpy as num

from anuga.config import g, epsilon
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross, rectangular
from anuga.shallow_water.shallow_water_domain import Domain
from anuga.shallow_water.boundaries import Reflective_boundary
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions import Dirichlet_boundary

class Test_Shallow_Water(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass


    def test_gravity(self):
        #Assuming no friction

        from anuga.config import g

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #             bac,     bce,     ecf,     dbe
        vertices = [[1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices)

        #Set up for a gradient of (3,0) at mid triangle (bce)
        def slope(x, y):
            return 3*x

        h = 0.1
        def stage(x, y):
            return slope(x,y) + h

        domain.set_quantity('elevation', slope)
        domain.set_quantity('stage', stage)

        for name in domain.conserved_quantities:
            assert num.allclose(domain.quantities[name].explicit_update, 0)
            assert num.allclose(domain.quantities[name].semi_implicit_update, 0)

        #domain.compute_forcing_terms()
        from anuga.shallow_water.sw_domain_openmp_ext import gravity
        gravity(domain)



        assert num.allclose(domain.quantities['stage'].explicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].explicit_update,
                            -g*h*3)
        assert num.allclose(domain.quantities['ymomentum'].explicit_update, 0)


    def test_gravity_2(self):
        #Assuming no friction

        from anuga.config import g

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #             bac,     bce,     ecf,     dbe
        vertices = [[1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices)

        #Set up for a gradient of (3,0) at mid triangle (bce)
        def slope(x, y):
            return 3*x

        h = 15
        def stage(x, y):
            return h

        domain.set_quantity('elevation', slope)
        domain.set_quantity('stage', stage)

        for name in domain.conserved_quantities:
            assert num.allclose(domain.quantities[name].explicit_update, 0)
            assert num.allclose(domain.quantities[name].semi_implicit_update, 0)

        #domain.compute_forcing_terms()
        from anuga.shallow_water.sw_domain_openmp_ext import gravity
        gravity(domain)



        assert num.allclose(domain.quantities['stage'].explicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].explicit_update,
                            [-382.2, -323.4, -205.8, -382.2])
        assert num.allclose(domain.quantities['ymomentum'].explicit_update, 0)


    def test_gravity_wb(self):
        #Assuming no friction

        from anuga.config import g

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #             bac,     bce,     ecf,     dbe
        vertices = [[1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices)

        #Set up for a gradient of (3,0) at mid triangle (bce)
        def slope(x, y):
            return 3*x

        h = 0.1
        def stage(x, y):
            return slope(x,y)+h

        domain.set_quantity('elevation', slope)
        domain.set_quantity('stage', stage)

        for name in domain.conserved_quantities:
            assert num.allclose(domain.quantities[name].explicit_update, 0)
            assert num.allclose(domain.quantities[name].semi_implicit_update, 0)

        from anuga.shallow_water.sw_domain_openmp_ext import gravity_wb
        gravity_wb(domain)




        assert num.allclose(domain.quantities['stage'].explicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].explicit_update,
                            -g*h*3)
        assert num.allclose(domain.quantities['ymomentum'].explicit_update, 0)


    def test_gravity_wb_2(self):
        #Assuming no friction

        from anuga.config import g

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #             bac,     bce,     ecf,     dbe
        vertices = [[1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices)

        #Set up for a gradient of (3,0) at mid triangle (bce)
        def slope(x, y):
            return 3*x

        h = 15.0
        def stage(x, y):
            return h

        domain.set_quantity('elevation', slope)
        domain.set_quantity('stage', stage)

        for name in domain.conserved_quantities:
            assert num.allclose(domain.quantities[name].explicit_update, 0)
            assert num.allclose(domain.quantities[name].semi_implicit_update, 0)

        from anuga.shallow_water.sw_domain_openmp_ext import gravity_wb
        gravity_wb(domain)



        assert num.allclose(domain.quantities['stage'].explicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].explicit_update,
                            [-396.9, -308.7, -220.5, -396.9])
        assert num.allclose(domain.quantities['ymomentum'].explicit_update, 0)


    def test_manning_friction(self):
        """ Assuming flat manning frinction is default
        """
        from anuga.config import g

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #             bac,     bce,     ecf,     dbe
        vertices = [[1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices)
        # Mode-2 ('unified') computes friction on-device and never syncs the host
        # semi_implicit_update array this white-box test reads. Pin legacy to
        # exercise the in-Python forcing machinery it is written for.
        domain.set_compute_mode('legacy')

        domain.set_sloped_mannings_function(False)

        B = Reflective_boundary(domain)
        domain.set_boundary( {'exterior': B})

        #Set up for a gradient of (3,0) at mid triangle (bce)
        def slope(x, y):
            return 3*x

        h = 0.1
        def stage(x, y):
            return slope(x, y) + h

        eta = 0.07
        domain.set_quantity('elevation', slope)
        domain.set_quantity('stage', stage)
        domain.set_quantity('friction', eta)

        for name in domain.conserved_quantities:
            assert num.allclose(domain.quantities[name].explicit_update, 0)
            assert num.allclose(domain.quantities[name].semi_implicit_update, 0)

        domain.compute_forcing_terms()

        assert num.allclose(domain.quantities['stage'].explicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].explicit_update,
                            0)
        assert num.allclose(domain.quantities['ymomentum'].explicit_update, 0)

        assert num.allclose(domain.quantities['stage'].semi_implicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].semi_implicit_update,
                            0)
        assert num.allclose(domain.quantities['ymomentum'].semi_implicit_update,
                            0)

        #Create some momentum for friction to work with
        domain.set_quantity('xmomentum', 1)
        S = -g*eta**2 /  h**(7.0/3)

        domain.compute_forcing_terms()
        assert num.allclose(domain.quantities['stage'].semi_implicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].semi_implicit_update,
                            S)
        assert num.allclose(domain.quantities['ymomentum'].semi_implicit_update,
                            0)

        #A more complex example
        domain.quantities['stage'].semi_implicit_update[:] = 0.0
        domain.quantities['xmomentum'].semi_implicit_update[:] = 0.0
        domain.quantities['ymomentum'].semi_implicit_update[:] = 0.0

        domain.set_quantity('xmomentum', 3)
        domain.set_quantity('ymomentum', 4)

        S = -g*eta**2*5 /  h**(7.0/3)

        domain.compute_forcing_terms()

        assert num.allclose(domain.quantities['stage'].semi_implicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].semi_implicit_update,
                            3*S)
        assert num.allclose(domain.quantities['ymomentum'].semi_implicit_update,
                            4*S)


    def test_flat_manning_friction(self):
        from anuga.config import g

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #             bac,     bce,     ecf,     dbe
        vertices = [[1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices)
        # Mode-2 ('unified') computes friction on-device and never syncs the host
        # semi_implicit_update array this white-box test reads. Pin legacy to
        # exercise the in-Python forcing machinery it is written for.
        domain.set_compute_mode('legacy')
        B = Reflective_boundary(domain)
        domain.set_boundary( {'exterior': B})

        # Use the flat function which doesn't takes into account the extra
        # wetted area due to slope of bed
        domain.set_sloped_mannings_function(False)

        #Set up for a gradient of (3,0) at mid triangle (bce)
        def slope(x, y):
            return 3*x

        h = 0.1
        def stage(x, y):
            return slope(x, y) + h

        eta = 0.07
        domain.set_quantity('elevation', slope)
        domain.set_quantity('stage', stage)
        domain.set_quantity('friction', eta)

        for name in domain.conserved_quantities:
            assert num.allclose(domain.quantities[name].explicit_update, 0)
            assert num.allclose(domain.quantities[name].semi_implicit_update, 0)

        domain.compute_forcing_terms()

        assert num.allclose(domain.quantities['stage'].explicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].explicit_update,
                            0)
        assert num.allclose(domain.quantities['ymomentum'].explicit_update, 0)

        assert num.allclose(domain.quantities['stage'].semi_implicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].semi_implicit_update,
                            0)
        assert num.allclose(domain.quantities['ymomentum'].semi_implicit_update,
                            0)

        #Create some momentum for friction to work with
        domain.set_quantity('xmomentum', 1)
        S = -g*eta**2 /  h**(7.0/3)

        domain.compute_forcing_terms()
        assert num.allclose(domain.quantities['stage'].semi_implicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].semi_implicit_update,
                            S)
        assert num.allclose(domain.quantities['ymomentum'].semi_implicit_update,
                            0)

        #A more complex example
        domain.quantities['stage'].semi_implicit_update[:] = 0.0
        domain.quantities['xmomentum'].semi_implicit_update[:] = 0.0
        domain.quantities['ymomentum'].semi_implicit_update[:] = 0.0

        domain.set_quantity('xmomentum', 3)
        domain.set_quantity('ymomentum', 4)

        S = -g*eta**2*5 /  h**(7.0/3)

        domain.compute_forcing_terms()

        assert num.allclose(domain.quantities['stage'].semi_implicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].semi_implicit_update,
                            3*S)
        assert num.allclose(domain.quantities['ymomentum'].semi_implicit_update,
                            4*S)


    def test_sloped_manning_friction(self):
        from anuga.config import g

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #             bac,     bce,     ecf,     dbe
        vertices = [[1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices)
        # Mode-2 ('unified') computes friction on-device and never syncs the host
        # semi_implicit_update array this white-box test reads. Pin legacy to
        # exercise the in-Python forcing machinery it is written for.
        domain.set_compute_mode('legacy')
        B = Reflective_boundary(domain)
        domain.set_boundary( {'exterior': B})

        # Use the sloped function which takes into account the extra
        # wetted area due to slope of bed
        domain.set_sloped_mannings_function(True)

        #Set up for a gradient of (3,0) at mid triangle (bce)
        def slope(x, y):
            return 3*x

        h = 0.1
        def stage(x, y):
            return slope(x, y) + h

        eta = 0.07
        domain.set_quantity('elevation', slope)
        domain.set_quantity('stage', stage)
        domain.set_quantity('friction', eta)

        for name in domain.conserved_quantities:
            assert num.allclose(domain.quantities[name].explicit_update, 0)
            assert num.allclose(domain.quantities[name].semi_implicit_update, 0)

        domain.compute_forcing_terms()

        assert num.allclose(domain.quantities['stage'].explicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].explicit_update,
                            0)
        assert num.allclose(domain.quantities['ymomentum'].explicit_update, 0)

        assert num.allclose(domain.quantities['stage'].semi_implicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].semi_implicit_update,
                            0)
        assert num.allclose(domain.quantities['ymomentum'].semi_implicit_update,
                            0)

        #Create some momentum for friction to work with
        domain.set_quantity('xmomentum', 1)
        S = -g*eta**2 /  h**(7.0/3) * sqrt(10)

        domain.compute_forcing_terms()
        assert num.allclose(domain.quantities['stage'].semi_implicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].semi_implicit_update,
                            S)
        assert num.allclose(domain.quantities['ymomentum'].semi_implicit_update,
                            0)

        #A more complex example
        domain.quantities['stage'].semi_implicit_update[:] = 0.0
        domain.quantities['xmomentum'].semi_implicit_update[:] = 0.0
        domain.quantities['ymomentum'].semi_implicit_update[:] = 0.0

        domain.set_quantity('xmomentum', 3)
        domain.set_quantity('ymomentum', 4)

        S = -g*eta**2*5 /  h**(7.0/3) * sqrt(10.0)

        domain.compute_forcing_terms()

        assert num.allclose(domain.quantities['stage'].semi_implicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].semi_implicit_update,
                            3*S)
        assert num.allclose(domain.quantities['ymomentum'].semi_implicit_update,
                            4*S)


    def test_second_order_flat_bed_onestep(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        #Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.set_flow_algorithm('DE0')
        domain.smooth = False
        domain.default_order = 2
        domain.beta_w = 0.9
        domain.beta_w_dry = 0.9
        domain.beta_uh = 0.9
        domain.beta_uh_dry = 0.9
        domain.beta_vh = 0.9
        domain.beta_vh_dry = 0.9
        domain.H0 = 1.0e-3

        # Boundary conditions
        Br = Reflective_boundary(domain)
        Bd = Dirichlet_boundary([0.1, 0., 0.])
        domain.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})

        domain.check_integrity()

        # Evolution
        for t in domain.evolve(yieldstep=0.05, finaltime=0.05):
            pass

        # Data from earlier version of abstract_2d_finite_volumes
        assert num.allclose(domain.recorded_min_timestep, 0.03571428571428564)
        assert num.allclose(domain.recorded_max_timestep, 0.03571428571428564)

        msg = ("domain.quantities['stage'].centroid_values[:12]=%s"
               % str(domain.quantities['stage'].centroid_values[:12]))
        assert num.allclose(domain.quantities['stage'].centroid_values[:12],
                            [0.00117244, 0.025897, 0.00200148, 0.025897,
                             0.00200148, 0.025897, 0.00200148, 0.025897,
                             0.00200148, 0.025897, 0.00200148, 0.02672604
                            ]), msg

        domain.distribute_to_vertices_and_edges()

        assert num.allclose(domain.quantities['stage'].vertex_values[:12,0],
                            [-0.00104301, 0.025897, 0.00020015, 0.025897,
                              0.00020015, 0.025897, 0.00020015, 0.025897,
                              0.00020015, 0.025897, 0.00011501, 0.02672604])

        assert num.allclose(domain.quantities['stage'].vertex_values[:12,1],
                            [0.00328283, 0.025897, 0.00020015, 0.025897,
                             0.00020015, 0.025897, 0.00020015, 0.025897,
                             0.00020015, 0.025897, 0.00028528, 0.02672604])

        assert num.allclose(domain.quantities['stage'].vertex_values[:12,2],
                            [0.0012775, 0.025897, 0.00560414, 0.025897,
                             0.00560414, 0.025897, 0.00560414, 0.025897,
                             0.00560414, 0.025897, 0.00560414, 0.02672604])


        assert num.allclose(domain.quantities['xmomentum'].centroid_values[:12],
                            [0.000189, 0.0042, 0.000189, 0.0042,
                             0.000189, 0.0042, 0.000189, 0.0042,
                             0.000189, 0.0042, 0.000189, 0.0042])


        assert num.allclose(domain.quantities['ymomentum'].centroid_values[:12],
                            [-1.89000000e-04, 0.00000000e+00, 0.00000000e+00,
                              0.00000000e+00, -5.57589689e-20, -2.78794844e-20,
                              8.36384533e-20,  0.00000000e+00, -3.06674329e-19,
                              0.00000000e+00,  2.23035875e-19, -1.89000000e-04])


        os.remove(os.path.join(domain.get_datadir(), domain.get_name() + '.sww'))


    def test_second_order_flat_bed_moresteps(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order = 2

        # Boundary conditions
        Br = Reflective_boundary(domain)
        Bd = Dirichlet_boundary([0.1, 0., 0.])
        domain.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})

        domain.check_integrity()

        # Evolution
        for t in domain.evolve(yieldstep=0.05, finaltime=0.1):
            pass

        # Data from earlier version of abstract_2d_finite_volumes
        #assert allclose(domain.recorded_min_timestep, 0.0396825396825)
        #assert allclose(domain.recorded_max_timestep, 0.0396825396825)

        os.remove(os.path.join(domain.get_datadir(), domain.get_name() + '.sww'))


    def test_flatbed_first_order(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        N = 8
        points, vertices, boundary = rectangular(N, N)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.set_flow_algorithm('DE0')
        domain.smooth = False
        domain.default_order = 1
        domain.H0 = 1.0e-3 # As suggested in the manual

        # Boundary conditions
        Br = Reflective_boundary(domain)
        Bd = Dirichlet_boundary([0.2, 0., 0.])

        domain.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})
        domain.check_integrity()

        # Evolution
        for t in domain.evolve(yieldstep=0.02, finaltime=0.5):
            pass


        assert num.allclose(domain.quantities['xmomentum'].\
                            edge_values[:4, 0],
                            [0.07335652, 0.06685681, 0.07071273, 0.06628975])

        assert num.allclose(domain.quantities['xmomentum'].\
                            edge_values[:4, 1],
                            [0.07343497, 0.06685681, 0.07083783, 0.06628975])

        assert num.allclose(domain.quantities['xmomentum'].\
                            edge_values[:4, 2],
                            [0.08124162, 0.06685681, 0.07891946, 0.06628975])
        os.remove(os.path.join(domain.get_datadir(), domain.get_name() + '.sww'))


    def test_flatbed_second_order(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        N = 8
        points, vertices, boundary = rectangular(N, N)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.set_flow_algorithm('DE0')
        domain.smooth = False
        domain.default_order = 2
        domain.beta_w = 0.9
        domain.beta_w_dry = 0.9
        domain.beta_uh = 0.9
        domain.beta_uh_dry = 0.9
        domain.beta_vh = 0.9
        domain.beta_vh_dry = 0.9

        domain.H0 = 1.0e-3 # As suggested in the manual
        domain.use_centroid_velocities = False # Backwards compatibility (8/5/8)
        domain.set_maximum_allowed_speed(1.0)

        # Boundary conditions
        Br = Reflective_boundary(domain)
        Bd = Dirichlet_boundary([0.2, 0., 0.])

        domain.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})
        domain.check_integrity()

        # Evolution
        for t in domain.evolve(yieldstep=0.01, finaltime=0.03):
            pass

        msg = 'Min timestep was %f instead of %f' % (domain.recorded_min_timestep,
                                                     0.018940360)
        assert num.allclose(domain.recorded_min_timestep, 0.018940360), msg

        msg = 'Max timestep was %f instead of %f' % (domain.recorded_max_timestep,
                                                     0.018940360)
        assert num.allclose(domain.recorded_max_timestep, 0.018940360), msg

        W_0 = [-0.00524972, 0.05350326, 0.00077479, 0.05356756]
        UH_0 = [-0.00271431, 0.02744767, 0.00023944, 0.02746294]
        VH_0 = [1.10413075e-03, 2.62134850e-04, -2.72890315e-05, 2.77104392e-04]


        assert num.allclose(domain.quantities['stage'].vertex_values[:4, 0], W_0)
        assert num.allclose(domain.quantities['xmomentum'].vertex_values[:4, 0], UH_0)
        assert num.allclose(domain.quantities['ymomentum'].vertex_values[:4, 0], VH_0)


        os.remove(os.path.join(domain.get_datadir(), domain.get_name() + '.sww'))


    def test_flatbed_second_order_vmax_0(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        N = 8
        points, vertices, boundary = rectangular(N, N)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.set_flow_algorithm('DE0')
        domain.smooth = False
        domain.default_order = 2
        domain.beta_w = 0.9
        domain.beta_w_dry = 0.9
        domain.beta_uh = 0.9
        domain.beta_uh_dry = 0.9
        domain.beta_vh = 0.9
        domain.beta_vh_dry = 0.9
        domain.maximum_allowed_speed = 0.0    # Makes it like the 'oldstyle'
        domain.H0 = 1.0e-3 # As suggested in the manual
        domain.use_centroid_velocities = False # Backwards compatibility (8/5/8)

        # Boundary conditions
        Br = Reflective_boundary(domain)
        Bd = Dirichlet_boundary([0.2, 0., 0.])

        domain.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})
        domain.check_integrity()

        # Evolution
        for t in domain.evolve(yieldstep=0.01, finaltime=0.03):
            pass

        assert num.allclose(domain.recorded_min_timestep, 0.018940360210353942)
        assert num.allclose(domain.recorded_max_timestep, 0.018940360210353942)

        UH_EX = [-0.00271431, 0.02744767,  0.00023944, 0.02746294]
        VH_EX = [1.10413075e-03, 2.62134850e-04, -2.72890315e-05, 2.77104392e-04]

        assert num.allclose(domain.quantities['xmomentum'].vertex_values[:4, 0], UH_EX)
        assert num.allclose(domain.quantities['ymomentum'].vertex_values[:4, 0], VH_EX)

        os.remove(os.path.join(domain.get_datadir(), domain.get_name() + '.sww'))


    def test_flatbed_second_order_distribute(self):
        #Use real data from anuga.abstract_2d_finite_volumes 2
        #painfully setup and extracted.

        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        N = 8
        points, vertices, boundary = rectangular(N, N)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.set_flow_algorithm('DE0')
        domain.smooth = False
        domain.default_order=domain._order_ = 2
        domain.beta_w = 0.9
        domain.beta_w_dry = 0.9
        domain.beta_uh = 0.9
        domain.beta_uh_dry = 0.9
        domain.beta_vh = 0.9
        domain.beta_vh_dry = 0.9

        domain.H0 = 1.0e-3 # ANUGA manual 28/5/9

        # Boundary conditions
        Br = Reflective_boundary(domain)
        Bd = Dirichlet_boundary([0.2, 0., 0.])

        domain.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})
        domain.check_integrity()

        for V in [False, True]:
            if V:
                # Set centroids as if system had been evolved
                L = num.zeros(2*N*N, float)
                L[:32] = [7.21205592e-003, 5.35214298e-002, 1.00910824e-002,
                          5.35439433e-002, 1.00910824e-002, 5.35439433e-002,
                          1.00910824e-002, 5.35439433e-002, 1.00910824e-002,
                          5.35439433e-002, 1.00910824e-002, 5.35439433e-002,
                          1.00910824e-002, 5.35393928e-002, 1.02344264e-002,
                          5.59605058e-002, 0.00000000e+000, 3.31027800e-004,
                          0.00000000e+000, 4.37962142e-005, 0.00000000e+000,
                          4.37962142e-005, 0.00000000e+000, 4.37962142e-005,
                          0.00000000e+000, 4.37962142e-005, 0.00000000e+000,
                          4.37962142e-005, 0.00000000e+000, 4.37962142e-005,
                          0.00000000e+000, 5.57305948e-005]

                X = num.zeros(2*N*N, float)
                X[:32] = [6.48351607e-003, 3.68571894e-002, 8.50733285e-003,
                          3.68731327e-002, 8.50733285e-003, 3.68731327e-002,
                          8.50733285e-003, 3.68731327e-002, 8.50733285e-003,
                          3.68731327e-002, 8.50733285e-003, 3.68731327e-002,
                          8.50733285e-003, 3.68693861e-002, 8.65220973e-003,
                          3.85055387e-002, 0.00000000e+000, 2.86060840e-004,
                          0.00000000e+000, 3.58905503e-005, 0.00000000e+000,
                          3.58905503e-005, 0.00000000e+000, 3.58905503e-005,
                          0.00000000e+000, 3.58905503e-005, 0.00000000e+000,
                          3.58905503e-005, 0.00000000e+000, 3.58905503e-005,
                          0.00000000e+000, 4.57662812e-005]

                Y = num.zeros(2*N*N, float)
                Y[:32] = [-1.39463104e-003, 6.15600298e-004, -6.03637382e-004,
                          6.18272251e-004, -6.03637382e-004, 6.18272251e-004,
                          -6.03637382e-004, 6.18272251e-004, -6.03637382e-004,
                          6.18272251e-004, -6.03637382e-004, 6.18272251e-004,
                          -6.03637382e-004, 6.18599320e-004, -6.74622797e-004,
                          -1.48934756e-004, 0.00000000e+000, -5.35079969e-005,
                          0.00000000e+000, -2.57264987e-005, 0.00000000e+000,
                          -2.57264987e-005, 0.00000000e+000, -2.57264987e-005,
                          0.00000000e+000, -2.57264987e-005, 0.00000000e+000,
                          -2.57264987e-005, 0.00000000e+000, -2.57264987e-005,
                          0.00000000e+000, -2.57635178e-005]

                domain.set_quantity('stage', L, location='centroids')
                domain.set_quantity('xmomentum', X, location='centroids')
                domain.set_quantity('ymomentum', Y, location='centroids')

                domain.check_integrity()
            else:
                # Evolution
                for t in domain.evolve(yieldstep=0.01, finaltime=0.03):
                    pass

                assert num.allclose(domain.recorded_min_timestep, 0.018940360)
                assert num.allclose(domain.recorded_max_timestep, 0.018940360)


            #Centroids were correct but not vertices.
            #Hence the check of distribute below.

            if not V:

                W_EX = [0.00519999, 0.05350326, 0.00786757, 0.05356756]
                UH_EX = [0.0026886, 0.02746638, 0.00346158, 0.02746294]
                VH_EX = [-0.00109367, 0.00026213, -0.0002771, 0.0002771]

                assert num.allclose(domain.quantities['stage'].centroid_values[:4], W_EX)
                assert num.allclose(domain.quantities['xmomentum'].centroid_values[:4], UH_EX)
                assert num.allclose(domain.quantities['ymomentum'].centroid_values[:4], VH_EX)

                assert num.allclose(domain.quantities['xmomentum'].centroid_values[17], 0.0, atol=1.0e-3)
            else:
                assert num.allclose(domain.quantities['xmomentum'].\
                                        centroid_values[17],
                                    0.00028606084)
                return #FIXME - Bailout for V True

            import copy

            XX = copy.copy(domain.quantities['xmomentum'].centroid_values)
            assert num.allclose(domain.quantities['xmomentum'].centroid_values,
                                XX)

            domain.distribute_to_vertices_and_edges()

            assert num.allclose(domain.quantities['xmomentum'].centroid_values[17], 0.0, atol=1.0e-3)


            UH_EX = [-0.00271431,  0.02744767,  0.00023944,  0.02746294]
            VH_EX = [1.10413075e-03,  2.62134850e-04, -2.72890315e-05,  2.77104392e-04]



            assert num.allclose(domain.quantities['xmomentum'].vertex_values[:4,0],
                    UH_EX)

            assert num.allclose(domain.quantities['ymomentum'].vertex_values[:4,0],
                    VH_EX)


        os.remove(os.path.join(domain.get_datadir(), domain.get_name() + '.sww'))


    def test_bedslope_problem_first_order(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order = 1

        # Bed-slope and friction
        def x_slope(x, y):
            return -x / 3

        domain.set_quantity('elevation', x_slope)

        # Boundary conditions
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        # Initial condition
        domain.set_quantity('stage', expression='elevation+0.05')
        domain.check_integrity()

        # Evolution
        for t in domain.evolve(yieldstep=0.05, finaltime=0.05):
            pass

        # FIXME (Ole): Need some other assertion here!
        #assert allclose(domain.recorded_min_timestep, 0.050010003001)
        #assert allclose(domain.recorded_max_timestep, 0.050010003001)

        os.remove(os.path.join(domain.get_datadir(), domain.get_name() + '.sww'))


    def test_bedslope_problem_first_order_moresteps_low_froude_0(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order = 1

        domain.set_low_froude(0)

        # FIXME (Ole): Need tests where these two are commented out
        domain.H0 = 0                         # Backwards compatibility (6/2/7)
        domain.tight_slope_limiters = 0       # Backwards compatibility (14/4/7)
        domain.use_centroid_velocities = 0    # Backwards compatibility (7/5/8)

        # Bed-slope and friction
        def x_slope(x, y):
            return -x / 3

        domain.set_quantity('elevation', x_slope)

        # Boundary conditions
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        # Initial condition
        domain.set_quantity('stage', expression='elevation+0.05')
        domain.check_integrity()

        # Evolution
        for t in domain.evolve(yieldstep=0.05, finaltime=0.5):
            pass


        W_EX = num.array([-0.02933401, -0.01091747, -0.02782667, -0.01057726, -0.02746728,
       -0.01033124, -0.02724838, -0.01019733, -0.02684617, -0.00965007,
       -0.02495054, -0.00706089, -0.06890704, -0.06303447, -0.07568528,
       -0.06161392, -0.07566147, -0.06125286, -0.07537564, -0.06094254,
       -0.07488436, -0.06058909, -0.07465553, -0.06045194, -0.12205956,
       -0.10701804, -0.12318423, -0.10795678, -0.12302672, -0.10763744,
       -0.12279702, -0.10731376, -0.12214772, -0.1071946 , -0.12283814,
       -0.10858945, -0.16389971, -0.15304369, -0.16428248, -0.15298725,
       -0.16430596, -0.15254043, -0.16407102, -0.15236021, -0.16390702,
       -0.15187689, -0.16460111, -0.15549742, -0.2083898 , -0.19912702,
       -0.20845637, -0.19709535, -0.20648985, -0.19663999, -0.20525695,
       -0.19631453, -0.20446298, -0.19578452, -0.20735693, -0.19565282,
       -0.14022868, -0.14262659, -0.13774131, -0.14132395, -0.13707526,
       -0.14041639, -0.13594765, -0.13910709, -0.13533594, -0.1393996 ,
       -0.1328638 , -0.1363085 ])


        #Data from earlier version of abstract_2d_finite_volumes
        assert num.allclose(domain.quantities['stage'].centroid_values, W_EX)


        os.remove(os.path.join(domain.get_datadir(), domain.get_name() + '.sww'))


    def test_bedslope_problem_first_order_moresteps_low_froude_1(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order = 1

        domain.set_low_froude(0)

        # FIXME (Ole): Need tests where these two are commented out
        domain.H0 = 0                         # Backwards compatibility (6/2/7)
        domain.tight_slope_limiters = 0       # Backwards compatibility (14/4/7)
        domain.use_centroid_velocities = 0    # Backwards compatibility (7/5/8)

        # Bed-slope and friction
        def x_slope(x, y):
            return -x / 3

        domain.set_quantity('elevation', x_slope)

        # Boundary conditions
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        # Initial condition
        domain.set_quantity('stage', expression='elevation+0.05')
        domain.check_integrity()

        # Evolution
        for t in domain.evolve(yieldstep=0.05, finaltime=0.5):
            pass


        W_EX = num.array([-0.02933401, -0.01091747, -0.02782667, -0.01057726, -0.02746728,
       -0.01033124, -0.02724838, -0.01019733, -0.02684617, -0.00965007,
       -0.02495054, -0.00706089, -0.06890704, -0.06303447, -0.07568528,
       -0.06161392, -0.07566147, -0.06125286, -0.07537564, -0.06094254,
       -0.07488436, -0.06058909, -0.07465553, -0.06045194, -0.12205956,
       -0.10701804, -0.12318423, -0.10795678, -0.12302672, -0.10763744,
       -0.12279702, -0.10731376, -0.12214772, -0.1071946 , -0.12283814,
       -0.10858945, -0.16389971, -0.15304369, -0.16428248, -0.15298725,
       -0.16430596, -0.15254043, -0.16407102, -0.15236021, -0.16390702,
       -0.15187689, -0.16460111, -0.15549742, -0.2083898 , -0.19912702,
       -0.20845637, -0.19709535, -0.20648985, -0.19663999, -0.20525695,
       -0.19631453, -0.20446298, -0.19578452, -0.20735693, -0.19565282,
       -0.14022868, -0.14262659, -0.13774131, -0.14132395, -0.13707526,
       -0.14041639, -0.13594765, -0.13910709, -0.13533594, -0.1393996 ,
       -0.1328638 , -0.1363085 ])


        #Data from earlier version of abstract_2d_finite_volumes
        assert num.allclose(domain.quantities['stage'].centroid_values, W_EX)


        os.remove(os.path.join(domain.get_datadir(), domain.get_name() + '.sww'))


    def test_bedslope_problem_second_order_one_step(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order = 2
        domain.beta_w = 0.9
        domain.beta_w_dry = 0.9
        domain.beta_uh = 0.9
        domain.beta_uh_dry = 0.9
        domain.beta_vh = 0.9
        domain.beta_vh_dry = 0.9


        # FIXME (Ole): Need tests where this is commented out
        domain.tight_slope_limiters = 0       # Backwards compatibility (14/4/7)
        domain.use_centroid_velocities = 0    # Backwards compatibility (7/5/8)

        # Bed-slope and friction at vertices (and interpolated elsewhere)
        def x_slope(x, y):
            return -x / 3

        domain.set_quantity('elevation', x_slope)

        # Boundary conditions
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        #Initial condition
        domain.set_quantity('stage', expression='elevation+0.05')
        domain.check_integrity()

        assert num.allclose(domain.quantities['stage'].centroid_values,
                            [ 0.01296296,  0.03148148,  0.01296296,
                              0.03148148,  0.01296296,  0.03148148,
                              0.01296296,  0.03148148,  0.01296296,
                              0.03148148,  0.01296296,  0.03148148,
                             -0.04259259, -0.02407407, -0.04259259,
                             -0.02407407, -0.04259259, -0.02407407,
                             -0.04259259, -0.02407407, -0.04259259,
                             -0.02407407, -0.04259259, -0.02407407,
                             -0.09814815, -0.07962963, -0.09814815,
                             -0.07962963, -0.09814815, -0.07962963,
                             -0.09814815, -0.07962963, -0.09814815,
                             -0.07962963, -0.09814815, -0.07962963,
                             -0.1537037,  -0.13518519, -0.1537037,
                             -0.13518519, -0.1537037,  -0.13518519,
                             -0.1537037 , -0.13518519, -0.1537037,
                             -0.13518519, -0.1537037,  -0.13518519,
                             -0.20925926, -0.19074074, -0.20925926,
                             -0.19074074, -0.20925926, -0.19074074,
                             -0.20925926, -0.19074074, -0.20925926,
                             -0.19074074, -0.20925926, -0.19074074,
                             -0.26481481, -0.2462963,  -0.26481481,
                             -0.2462963,  -0.26481481, -0.2462963,
                             -0.26481481, -0.2462963,  -0.26481481,
                             -0.2462963,  -0.26481481, -0.2462963])


        # Evolution
        for t in domain.evolve(yieldstep=0.05, finaltime=0.05):
            pass



        W_EX = [ 0.01571282,  0.02678718,  0.01765727,  0.02678718,  0.01765727,
        0.02678718,  0.01765727,  0.02678718,  0.01765727,  0.02678718,
        0.01765727,  0.02873162, -0.04259259, -0.02407407, -0.04259259,
       -0.02407407, -0.04259259, -0.02407407, -0.04259259, -0.02407407,
       -0.04259259, -0.02407407, -0.04259259, -0.02407407, -0.09814815,
       -0.07962963, -0.09814815, -0.07962963, -0.09814815, -0.07962963,
       -0.09814815, -0.07962963, -0.09814815, -0.07962963, -0.09814815,
       -0.07962963, -0.1537037 , -0.13518519, -0.1537037 , -0.13518519,
       -0.1537037 , -0.13518519, -0.1537037 , -0.13518519, -0.1537037 ,
       -0.13518519, -0.1537037 , -0.13518519, -0.20925926, -0.19074074,
       -0.20925926, -0.19074074, -0.20925926, -0.19074074, -0.20925926,
       -0.19074074, -0.20925926, -0.19074074, -0.20925926, -0.19074074,
       -0.26206496, -0.2509906 , -0.26012051, -0.2509906 , -0.26012051,
       -0.2509906 , -0.26012051, -0.2509906 , -0.26012051, -0.2509906 ,
       -0.26012051, -0.24904616]


        assert num.allclose(domain.quantities['stage'].centroid_values, W_EX)


        os.remove(os.path.join(domain.get_datadir(), domain.get_name() + '.sww'))


    def test_bedslope_problem_second_order_two_steps(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)

        domain.set_low_froude(0)

        domain.smooth = False
        domain.default_order = 2
        domain.beta_w = 0.9
        domain.beta_w_dry = 0.9
        domain.beta_uh = 0.9
        domain.beta_uh_dry = 0.9
        domain.beta_vh = 0.9
        domain.beta_vh_dry = 0.9

        # FIXME (Ole): Need tests where this is commented out
        domain.tight_slope_limiters = 0    # Backwards compatibility (14/4/7)
        domain.H0 = 0    # Backwards compatibility (6/2/7)
        domain.use_centroid_velocities = 0    # Backwards compatibility (7/5/8)

        # Bed-slope and friction at vertices (and interpolated elsewhere)
        def x_slope(x, y):
            return -x / 3

        domain.set_quantity('elevation', x_slope)

        # Boundary conditions
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        # Initial condition
        domain.set_quantity('stage', expression='elevation+0.05')
        domain.check_integrity()

        assert num.allclose(domain.quantities['stage'].centroid_values,
                            [ 0.01296296,  0.03148148,  0.01296296,
                              0.03148148,  0.01296296,  0.03148148,
                              0.01296296,  0.03148148,  0.01296296,
                              0.03148148,  0.01296296,  0.03148148,
                             -0.04259259, -0.02407407, -0.04259259,
                             -0.02407407, -0.04259259, -0.02407407,
                             -0.04259259, -0.02407407, -0.04259259,
                             -0.02407407, -0.04259259, -0.02407407,
                             -0.09814815, -0.07962963, -0.09814815,
                             -0.07962963, -0.09814815, -0.07962963,
                             -0.09814815, -0.07962963, -0.09814815,
                             -0.07962963, -0.09814815, -0.07962963,
                             -0.1537037 , -0.13518519, -0.1537037,
                             -0.13518519, -0.1537037,  -0.13518519,
                             -0.1537037 , -0.13518519, -0.1537037,
                             -0.13518519, -0.1537037,  -0.13518519,
                             -0.20925926, -0.19074074, -0.20925926,
                             -0.19074074, -0.20925926, -0.19074074,
                             -0.20925926, -0.19074074, -0.20925926,
                             -0.19074074, -0.20925926, -0.19074074,
                             -0.26481481, -0.2462963,  -0.26481481,
                             -0.2462963,  -0.26481481, -0.2462963,
                             -0.26481481, -0.2462963,  -0.26481481,
                             -0.2462963,  -0.26481481, -0.2462963])

        # Evolution
        for t in domain.evolve(yieldstep=0.05, finaltime=0.1):
            pass

        # Data from earlier version of abstract_2d_finite_volumes ft=0.1
        assert num.allclose(domain.recorded_min_timestep, 0.0344607459654)
        assert num.allclose(domain.recorded_max_timestep, 0.0391090502542)



        W_EX = [ 0.01308246,  0.02201217,  0.01358687,  0.023575  ,  0.01370201,
        0.0235574 ,  0.01370201,  0.02355753,  0.0136998 ,  0.02361447,
        0.01454146,  0.02507853, -0.04254524, -0.022535  , -0.04260359,
       -0.0225144 , -0.04263851, -0.02252512, -0.04263809, -0.02252512,
       -0.04264303, -0.02249883, -0.04257228, -0.02296247, -0.0981472 ,
       -0.07964098, -0.09814502, -0.07968513, -0.09814337, -0.07968807,
       -0.09814296, -0.07968807, -0.09814296, -0.07968807, -0.09814555,
       -0.07965964, -0.15368138, -0.13518778, -0.15366285, -0.13519037,
       -0.15366285, -0.13519037, -0.15366285, -0.13518996, -0.15366479,
       -0.13518832, -0.15369898, -0.13518613, -0.21012122, -0.19071462,
       -0.21066837, -0.19057453, -0.21064807, -0.19057613, -0.21064807,
       -0.19057519, -0.21067682, -0.19061141, -0.21080033, -0.19066023,
       -0.25755282, -0.24897294, -0.25604086, -0.24820174, -0.25598276,
       -0.24820504, -0.25598272, -0.24820504, -0.2559903 , -0.2480719 ,
       -0.25395141, -0.24799187]


        assert num.allclose(domain.quantities['stage'].centroid_values, W_EX)



        os.remove(os.path.join(domain.get_datadir(), domain.get_name() + '.sww'))


    def test_bedslope_problem_second_order_two_yieldsteps(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        #Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)

        #domain.set_compute_fluxes_method('original')

        domain.set_low_froude(0)

        domain.smooth = False
        domain.default_order = 2
        domain.beta_w = 0.9
        domain.beta_w_dry = 0.9
        domain.beta_uh = 0.9
        domain.beta_uh_dry = 0.9
        domain.beta_vh = 0.9
        domain.beta_vh_dry = 0.9

        # FIXME (Ole): Need tests where this is commented out
        domain.tight_slope_limiters = 0    # Backwards compatibility (14/4/7)
        domain.H0 = 0    # Backwards compatibility (6/2/7)
        domain.use_centroid_velocities = 0    # Backwards compatibility (7/5/8)


        # Bed-slope and friction at vertices (and interpolated elsewhere)
        def x_slope(x, y):
            return -x / 3

        domain.set_quantity('elevation', x_slope)

        # Boundary conditions
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        # Initial condition
        domain.set_quantity('stage', expression='elevation+0.05')
        domain.check_integrity()

        assert num.allclose(domain.quantities['stage'].centroid_values,
                            [ 0.01296296,  0.03148148,  0.01296296,
                              0.03148148,  0.01296296,  0.03148148,
                              0.01296296,  0.03148148,  0.01296296,
                              0.03148148,  0.01296296,  0.03148148,
                             -0.04259259, -0.02407407, -0.04259259,
                             -0.02407407, -0.04259259, -0.02407407,
                             -0.04259259, -0.02407407, -0.04259259,
                             -0.02407407, -0.04259259, -0.02407407,
                             -0.09814815, -0.07962963, -0.09814815,
                             -0.07962963, -0.09814815, -0.07962963,
                             -0.09814815, -0.07962963, -0.09814815,
                             -0.07962963, -0.09814815, -0.07962963,
                             -0.1537037 , -0.13518519, -0.1537037,
                             -0.13518519, -0.1537037,  -0.13518519,
                             -0.1537037 , -0.13518519, -0.1537037,
                             -0.13518519, -0.1537037,  -0.13518519,
                             -0.20925926, -0.19074074, -0.20925926,
                             -0.19074074, -0.20925926, -0.19074074,
                             -0.20925926, -0.19074074, -0.20925926,
                             -0.19074074, -0.20925926, -0.19074074,
                             -0.26481481, -0.2462963,  -0.26481481,
                             -0.2462963,  -0.26481481, -0.2462963,
                             -0.26481481, -0.2462963,  -0.26481481,
                             -0.2462963,  -0.26481481, -0.2462963])


        # Evolution
        for t in domain.evolve(yieldstep=0.05, finaltime=0.1):   #0.05??
            pass




        W_EX = [ 0.01308246,  0.02201217,  0.01358687,  0.023575  ,  0.01370201,
        0.0235574 ,  0.01370201,  0.02355753,  0.0136998 ,  0.02361447,
        0.01454146,  0.02507853, -0.04254524, -0.022535  , -0.04260359,
       -0.0225144 , -0.04263851, -0.02252512, -0.04263809, -0.02252512,
       -0.04264303, -0.02249883, -0.04257228, -0.02296247, -0.0981472 ,
       -0.07964098, -0.09814502, -0.07968513, -0.09814337, -0.07968807,
       -0.09814296, -0.07968807, -0.09814296, -0.07968807, -0.09814555,
       -0.07965964, -0.15368138, -0.13518778, -0.15366285, -0.13519037,
       -0.15366285, -0.13519037, -0.15366285, -0.13518996, -0.15366479,
       -0.13518832, -0.15369898, -0.13518613, -0.21012122, -0.19071462,
       -0.21066837, -0.19057453, -0.21064807, -0.19057613, -0.21064807,
       -0.19057519, -0.21067682, -0.19061141, -0.21080033, -0.19066023,
       -0.25755282, -0.24897294, -0.25604086, -0.24820174, -0.25598276,
       -0.24820504, -0.25598272, -0.24820504, -0.2559903 , -0.2480719 ,
       -0.25395141, -0.24799187]


        assert num.allclose(domain.quantities['stage'].centroid_values, W_EX)


        os.remove(os.path.join(domain.get_datadir(), domain.get_name() + '.sww'))


    def test_bedslope_problem_second_order_more_steps(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)

        #domain.set_compute_fluxes_method('original')

        domain.smooth = False
        domain.default_order = 2
        domain.beta_w = 0.9
        domain.beta_w_dry = 0.9
        domain.beta_uh = 0.9
        domain.beta_uh_dry = 0.9
        domain.beta_vh = 0.9
        domain.beta_vh_dry = 0.9

        domain.set_low_froude(0)

        # FIXME (Ole): Need tests where these two are commented out
        domain.H0 = 0                      # Backwards compatibility (6/2/7)
        domain.tight_slope_limiters = 0    # Backwards compatibility (14/4/7)
        domain.use_centroid_velocities = 0    # Backwards compatibility (7/5/8)

        # Bed-slope and friction at vertices (and interpolated elsewhere)
        def x_slope(x, y):
            return -x / 3

        domain.set_quantity('elevation', x_slope)

        # Boundary conditions
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        # Initial condition
        domain.set_quantity('stage', expression='elevation+0.05')
        domain.check_integrity()



        W_EX = [ 0.01296296,  0.03148148,  0.01296296,  0.03148148,  0.01296296,
        0.03148148,  0.01296296,  0.03148148,  0.01296296,  0.03148148,
        0.01296296,  0.03148148, -0.04259259, -0.02407407, -0.04259259,
       -0.02407407, -0.04259259, -0.02407407, -0.04259259, -0.02407407,
       -0.04259259, -0.02407407, -0.04259259, -0.02407407, -0.09814815,
       -0.07962963, -0.09814815, -0.07962963, -0.09814815, -0.07962963,
       -0.09814815, -0.07962963, -0.09814815, -0.07962963, -0.09814815,
       -0.07962963, -0.1537037 , -0.13518519, -0.1537037 , -0.13518519,
       -0.1537037 , -0.13518519, -0.1537037 , -0.13518519, -0.1537037 ,
       -0.13518519, -0.1537037 , -0.13518519, -0.20925926, -0.19074074,
       -0.20925926, -0.19074074, -0.20925926, -0.19074074, -0.20925926,
       -0.19074074, -0.20925926, -0.19074074, -0.20925926, -0.19074074,
       -0.26481481, -0.2462963 , -0.26481481, -0.2462963 , -0.26481481,
       -0.2462963 , -0.26481481, -0.2462963 , -0.26481481, -0.2462963 ,
       -0.26481481, -0.2462963 ]

        assert num.allclose(domain.quantities['stage'].centroid_values, W_EX)

        # Evolution
        for t in domain.evolve(yieldstep=0.05, finaltime=0.5):
            # Check that diagnostics works
            msg = domain.timestepping_statistics(track_speeds=True)
            #FIXME(Ole): One might check the contents of msg here.







        W_EX = [-0.0301883 , -0.01127593, -0.02834861, -0.0108968 , -0.02806583,
       -0.01074475, -0.02788852, -0.01065176, -0.02753658, -0.01013457,
       -0.02602057, -0.00697983, -0.07333722, -0.06511481, -0.08002821,
       -0.06234795, -0.07884103, -0.06214974, -0.07867692, -0.06199716,
       -0.07848496, -0.06160417, -0.07782372, -0.06214378, -0.12379366,
       -0.10782217, -0.12402788, -0.10970479, -0.12407233, -0.10911802,
       -0.12375966, -0.10870709, -0.1233677 , -0.10897925, -0.1246598 ,
       -0.10913556, -0.16030637, -0.15059851, -0.1619466 , -0.15150383,
       -0.16222821, -0.15126241, -0.16204491, -0.15117909, -0.16170723,
       -0.15047256, -0.16159564, -0.15306895, -0.20650721, -0.19412466,
       -0.20661286, -0.193611  , -0.20611843, -0.1932463 , -0.20476096,
       -0.19336801, -0.20673157, -0.1928249 , -0.20666551, -0.18982871,
       -0.14177628, -0.13867427, -0.138984  , -0.13970989, -0.13837738,
       -0.13785944, -0.13772734, -0.13619726, -0.13708211, -0.13865006,
       -0.13413676, -0.14008116]




        assert num.allclose(domain.quantities['stage'].centroid_values, W_EX)



        UH_EX = [ 0.00439415,  0.00080998,  0.00359639,  0.00086799,  0.00369947,
        0.00089459,  0.00376748,  0.00091612,  0.00394185,  0.00100331,
        0.00423642,  0.00102681,  0.01866626,  0.00841126,  0.01265179,
        0.00989651,  0.01358962,  0.01007976,  0.01375394,  0.01018847,
        0.01388854,  0.01057534,  0.01516146,  0.01155302,  0.03130996,
        0.02572364,  0.03012732,  0.02348833,  0.02999708,  0.02419562,
        0.03046757,  0.02471555,  0.03088795,  0.02415592,  0.02886438,
        0.02425477,  0.06609542,  0.04946308,  0.06236933,  0.04715456,
        0.06176596,  0.04743353,  0.06204203,  0.04760272,  0.06262964,
        0.04858867,  0.06290163,  0.0447266 ,  0.08245285,  0.07715851,
        0.08533206,  0.07817769,  0.08610624,  0.07901377,  0.0879718 ,
        0.07879576,  0.08550203,  0.07963651,  0.08581205,  0.08376399,
        0.01997027,  0.07948531,  0.02703306,  0.08443647,  0.02694374,
        0.08404407,  0.02627307,  0.08400949,  0.02649757,  0.08397135,
        0.03080026,  0.09688469]




        assert num.allclose(domain.quantities['xmomentum'].centroid_values, UH_EX)





        VH_EX = [ -4.55438497e-04,   5.60449890e-05,  -1.05041759e-04,
         1.09474097e-04,   3.37384149e-05,   1.48217622e-04,
         8.15006870e-05,   1.41778277e-04,  -8.95699377e-05,
        -1.50989506e-05,  -1.13309713e-03,  -6.56027069e-04,
         4.08407333e-04,  -5.37824083e-04,  -6.07324848e-04,
        -1.79169050e-04,  -2.26068995e-04,  -9.40081480e-05,
        -1.92103118e-04,  -7.38055580e-05,  -2.61634534e-04,
        -2.18534710e-04,  -3.97532842e-04,  -4.38793389e-04,
         4.02437237e-04,   1.71996550e-04,  -4.71404571e-04,
        -7.07195798e-04,  -2.92606384e-04,  -2.89749829e-04,
        -2.14551908e-04,  -3.47586675e-04,  -4.51408518e-04,
        -2.85326203e-04,   2.60585517e-04,   1.85996351e-04,
         1.76070321e-03,   4.27131620e-04,  -1.36487970e-04,
        -7.48486546e-04,  -2.25542208e-04,  -4.38919415e-04,
        -1.23010639e-05,  -3.49340425e-04,  -2.36098938e-04,
        -7.63395297e-04,   6.66099385e-05,   6.22983544e-04,
         4.52833391e-04,   2.24808947e-03,   4.85972740e-04,
         7.50193061e-04,   4.53207970e-04,   1.05078134e-03,
         6.32114704e-04,   1.24869321e-03,   4.06577106e-04,
         1.01553335e-03,  -5.26250901e-04,   1.78062345e-04,
        -2.04979836e-03,  -2.70644308e-03,  -3.34067897e-03,
        -1.94716541e-03,  -2.28075505e-03,  -1.61494727e-03,
        -2.15457373e-03,  -1.48932625e-03,  -3.24804437e-03,
        -1.09719715e-03,   3.38706650e-03,  -8.98151209e-04]



        assert num.allclose(domain.quantities['ymomentum'].centroid_values, VH_EX)

        os.remove(os.path.join(domain.get_datadir(), domain.get_name() + '.sww'))


    def test_temp_play(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(5, 5)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order = 2
        domain.beta_w = 0.9
        domain.beta_w_dry = 0.9
        domain.beta_uh = 0.9
        domain.beta_uh_dry = 0.9
        domain.beta_vh = 0.9
        domain.beta_vh_dry = 0.9

        # FIXME (Ole): Need tests where these two are commented out
        domain.H0 = 0                         # Backwards compatibility (6/2/7)
        domain.tight_slope_limiters = False   # Backwards compatibility (14/4/7)
        domain.use_centroid_velocities = False # Backwards compatibility (7/5/8)
        domain.low_froude = 0                 # Backwards compatibility (25/6/19)

        # Bed-slope and friction at vertices (and interpolated elsewhere)
        def x_slope(x, y):
            return -x / 3

        domain.set_quantity('elevation', x_slope)

        # Boundary conditions
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        # Initial condition
        domain.set_quantity('stage', expression='elevation+0.05')
        domain.check_integrity()

        # Evolution
        for t in domain.evolve(yieldstep=0.05, finaltime=0.1):
            pass


        W = domain.quantities['stage'].centroid_values[:4]
        UH = domain.quantities['xmomentum'].centroid_values[:4]
        VH = domain.quantities['ymomentum'].centroid_values[:4]

        W_0  = [ 0.001362,    0.01344294,  0.00308829, 0.01470289]
        UH_0 = [ 0.01300239,  0.00537933,  0.01214676,  0.00515825]
        VH_0 = [ -1.13165691e-03,  -6.55330189e-04, -6.62804076e-05,   5.26313051e-05]

        W_1 = [ 0.00707892,  0.01849914,  0.00783274,  0.01997863]
        UH_1 = [ 0.01512518,  0.00354391,  0.01503765,  0.00326075]
        VH_1 = [  5.36531332e-04,  -6.77297008e-04,  -4.58560426e-05, 2.47714988e-05]


        assert num.allclose(W,W_0) or num.allclose(W,W_1)
        assert num.allclose(UH, UH_0) or num.allclose(UH, UH_1)
        assert num.allclose(VH, VH_0) or num.allclose(VH, VH_1)




        # old values pre revision 8402
#        assert num.allclose(domain.quantities['stage'].centroid_values[:4],
#                            [0.00206836, 0.01296714, 0.00363415, 0.01438924])
#        assert num.allclose(domain.quantities['xmomentum'].centroid_values[:4],
#                            [0.01360154, 0.00671133, 0.01264578, 0.00648503])
#        assert num.allclose(domain.quantities['ymomentum'].centroid_values[:4],
#                            [-1.19201077e-003, -7.23647546e-004,
#                             -6.39083123e-005, 6.29815168e-005])

        os.remove(os.path.join(domain.get_datadir(), domain.get_name() + '.sww'))



if __name__ == "__main__":
    unittest.main()

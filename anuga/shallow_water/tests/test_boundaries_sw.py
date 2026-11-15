#!/usr/bin/env python
"""Tests for boundary conditions in the shallow water domain."""
import os
import unittest
import tempfile
import numpy as num

from anuga.config import g, epsilon
from anuga.config import netcdf_mode_r
from anuga.utilities.numerical_tools import mean
from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross, rectangular
from anuga.shallow_water.shallow_water_domain import Domain
from anuga.file.netcdf import NetCDFFile
from anuga import create_domain_from_regions
from anuga.shallow_water.boundaries import (
    Reflective_boundary, Field_boundary,
    Transmissive_momentum_set_stage_boundary,
    Transmissive_n_momentum_zero_t_momentum_set_stage_boundary,
    Transmissive_stage_zero_momentum_boundary,
)
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions import (
    Boundary, Transmissive_boundary, Dirichlet_boundary, Time_boundary,
)

class Test_Shallow_Water(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass


    def test_boundary_conditions(self):
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #             bac,     bce,     ecf,     dbe
        vertices = [[1,0,2], [1,2,4], [4,2,5], [3,1,4]]
        boundary = {(0, 0): 'Third',
                    (0, 2): 'First',
                    (2, 0): 'Second',
                    (2, 1): 'Second',
                    (3, 1): 'Second',
                    (3, 2): 'Third'}

        domain = Domain(points, vertices, boundary)
        domain.check_integrity()

        domain.set_quantity('stage', [[1,2,3], [5,5,5], [0,0,9], [-6,3,3]])

        domain.set_quantity('xmomentum', [[1,1,1], [2,2,2], [3,3,3], [4,4,4]])

        domain.set_quantity('ymomentum', [[10,10,10], [20,20,20],
                                          [30,30,30], [40,40,40]])

        D = Dirichlet_boundary([5, 2, 1])
        T = Transmissive_boundary(domain)
        R = Reflective_boundary(domain)
        domain.set_boundary({'First': D, 'Second': T, 'Third': R})

        domain.update_boundary()

        # Stage
        assert domain.quantities['stage'].boundary_values[0] == 2.5
        # Reflective (2.5)
        assert (domain.quantities['stage'].boundary_values[0] ==
                domain.get_conserved_quantities(0, edge=0)[0])
        # Dirichlet
        assert domain.quantities['stage'].boundary_values[1] == 5.
        # Transmissive (4.5)
        assert (domain.quantities['stage'].boundary_values[2] ==
                domain.get_conserved_quantities(2, edge=0)[0])
        # Transmissive (4.5)
        assert (domain.quantities['stage'].boundary_values[3] ==
                domain.get_conserved_quantities(2, edge=1)[0])
        # Transmissive (-1.5)
        assert (domain.quantities['stage'].boundary_values[4] ==
                domain.get_conserved_quantities(3, edge=1)[0])
        # Reflective (-1.5)
        assert (domain.quantities['stage'].boundary_values[5] ==
                domain.get_conserved_quantities(3, edge=2)[0])

        # Xmomentum

        # Reflective
        assert domain.quantities['xmomentum'].boundary_values[0] == 1.0
        # Dirichlet
        assert domain.quantities['xmomentum'].boundary_values[1] == 2.
        # Transmissive
        assert (domain.quantities['xmomentum'].boundary_values[2] ==
                domain.get_conserved_quantities(2, edge=0)[1])
        # Transmissive
        assert (domain.quantities['xmomentum'].boundary_values[3] ==
                domain.get_conserved_quantities(2, edge=1)[1])
        # Transmissive
        assert (domain.quantities['xmomentum'].boundary_values[4] ==
                domain.get_conserved_quantities(3, edge=1)[1])
        # Reflective
        assert domain.quantities['xmomentum'].boundary_values[5] == -4.0

        # Ymomentum
        # Reflective
        assert domain.quantities['ymomentum'].boundary_values[0] == -10.0
        # Dirichlet
        assert domain.quantities['ymomentum'].boundary_values[1] == 1.
        # Transmissive
        assert domain.quantities['ymomentum'].boundary_values[2] == 30.
        # Transmissive
        assert domain.quantities['ymomentum'].boundary_values[3] == 30.
        # Transmissive
        assert domain.quantities['ymomentum'].boundary_values[4] == 40.
        # Reflective
        assert domain.quantities['ymomentum'].boundary_values[5] == 40.


    def test_boundary_conditionsII(self):
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #             bac,     bce,     ecf,     dbe
        vertices = [[1,0,2], [1,2,4], [4,2,5], [3,1,4]]
        boundary = {(0, 0): 'Third',
                    (0, 2): 'First',
                    (2, 0): 'Second',
                    (2, 1): 'Second',
                    (3, 1): 'Second',
                    (3, 2): 'Third',
                    (0, 1): 'Internal',
                    (1, 2): 'Internal'}

        domain = Domain(points, vertices, boundary)
        domain.check_integrity()

        domain.set_quantity('stage', [[1,2,3], [5,5,5], [0,0,9], [-6,3,3]])

        domain.set_quantity('xmomentum', [[1,1,1], [2,2,2], [3,3,3], [4,4,4]])

        domain.set_quantity('ymomentum', [[10,10,10], [20,20,20],
                                          [30,30,30], [40,40,40]])

        D = Dirichlet_boundary([5, 2, 1])
        T = Transmissive_boundary(domain)
        R = Reflective_boundary(domain)
        domain.set_boundary({'First': D, 'Second': T,
                             'Third': R, 'Internal': None})

        domain.update_boundary()
        domain.check_integrity()


    def test_boundary_conditions_inconsistent_internal(self):
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #             bac,     bce,     ecf,     dbe
        vertices = [[1,0,2], [1,2,4], [4,2,5], [3,1,4]]
        boundary = {(0, 0): 'Third',
                    (0, 2): 'First',
                    (2, 0): 'Second',
                    (2, 1): 'Second',
                    (3, 1): 'Second',
                    (3, 2): 'Third',
                    (0, 1): 'Internal'}

        domain = Domain(points, vertices, boundary)

        try:
            domain.check_integrity()
        except AssertionError:
            pass
        else:
            msg = 'should have thrown assertion exception'
            raise Exception(msg)


    def test_boundary_conditionsIII(self):
        """test_boundary_conditionsIII

        Test Transmissive_stage_zero_momentum_boundary
        """

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #             bac,     bce,     ecf,     dbe
        vertices = [[1,0,2], [1,2,4], [4,2,5], [3,1,4]]
        boundary = {(0, 0): 'Third',
                    (0, 2): 'First',
                    (2, 0): 'Second',
                    (2, 1): 'Second',
                    (3, 1): 'Second',
                    (3, 2): 'Third'}

        domain = Domain(points, vertices, boundary)
        domain.check_integrity()

        domain.set_quantity('stage', [[1,2,3], [5,5,5], [0,0,9], [-6,3,3]])

        domain.set_quantity('xmomentum', [[1,1,1], [2,2,2], [3,3,3], [4,4,4]])

        domain.set_quantity('ymomentum', [[10,10,10], [20,20,20],
                                          [30,30,30], [40,40,40]])

        D = Dirichlet_boundary([5, 2, 1])
        T = Transmissive_stage_zero_momentum_boundary(domain)
        R = Reflective_boundary(domain)
        domain.set_boundary({'First': D, 'Second': T, 'Third': R})

        domain.update_boundary()

        # Stage
        # Reflective (2.5)
        assert domain.quantities['stage'].boundary_values[0] == 2.5
        assert (domain.quantities['stage'].boundary_values[0] ==
                domain.get_conserved_quantities(0, edge=0)[0])
        # Dirichlet
        assert domain.quantities['stage'].boundary_values[1] == 5.
        # Transmissive (4.5)
        assert (domain.quantities['stage'].boundary_values[2] ==
                domain.get_conserved_quantities(2, edge=0)[0])
        # Transmissive (4.5)
        assert (domain.quantities['stage'].boundary_values[3] ==
                domain.get_conserved_quantities(2, edge=1)[0])
        # Transmissive (-1.5)
        assert (domain.quantities['stage'].boundary_values[4] ==
                domain.get_conserved_quantities(3, edge=1)[0])
        # Reflective (-1.5)
        assert (domain.quantities['stage'].boundary_values[5] ==
                domain.get_conserved_quantities(3, edge=2)[0])

        # Xmomentum
        # Reflective
        assert domain.quantities['xmomentum'].boundary_values[0] == 1.0
        # Dirichlet
        assert domain.quantities['xmomentum'].boundary_values[1] == 2.
        assert domain.quantities['xmomentum'].boundary_values[2] == 0.0
        assert domain.quantities['xmomentum'].boundary_values[3] == 0.0
        assert domain.quantities['xmomentum'].boundary_values[4] == 0.0
        # Reflective
        assert domain.quantities['xmomentum'].boundary_values[5] == -4.0

        # Ymomentum
        # Reflective
        assert domain.quantities['ymomentum'].boundary_values[0] == -10.0
        # Dirichlet
        assert domain.quantities['ymomentum'].boundary_values[1] == 1.
        assert domain.quantities['ymomentum'].boundary_values[2] == 0.0
        assert domain.quantities['ymomentum'].boundary_values[3] == 0.0
        assert domain.quantities['ymomentum'].boundary_values[4] == 0.0
        # Reflective
        assert domain.quantities['ymomentum'].boundary_values[5] == 40.


    def test_boundary_condition_time(self):
        """test_boundary_condition_time

        This tests that boundary conditions are evaluated
        using the right time from domain.
        """

        # Setup computational domain
        from anuga.abstract_2d_finite_volumes.mesh_factory \
                import rectangular_cross

        #--------------------------------------------------------------
        # Setup computational domain
        #--------------------------------------------------------------
        N = 5
        points, vertices, boundary = rectangular_cross(N, N)
        domain = Domain(points, vertices, boundary)

        #--------------------------------------------------------------
        # Setup initial conditions
        #--------------------------------------------------------------
        domain.set_quantity('elevation', 0.0)
        domain.set_quantity('friction', 0.0)
        domain.set_quantity('stage', 0.0)

        #--------------------------------------------------------------
        # Setup boundary conditions
        #--------------------------------------------------------------
        # Time dependent boundary
        Bt = Time_boundary(domain=domain, function=lambda t: [t, 0.0, 0.0])

        Br = Reflective_boundary(domain)              # Reflective wall

        domain.set_boundary({'left': Bt, 'right': Br, 'top': Br, 'bottom': Br})

        for t in domain.evolve(yieldstep = 10, finaltime = 20.0):
            q = Bt.evaluate()

            # FIXME (Ole): This test would not have passed in
            # changeset:5846.
            msg = 'Time boundary not evaluated correctly'
            assert num.allclose(t, q[0]), msg


    def test_spatio_temporal_boundary_1(self):
        """Test that boundary values can be read from file and interpolated
        in both time and space.

        Verify that the same steady state solution is arrived at and that
        time interpolation works.

        The full solution history is not exactly the same as
        file boundary must read and interpolate from *smoothed* version
        as stored in sww.
        """

        # Create sww file of simple propagation from left to right
        # through rectangular domain
        verbose = False

        import time
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(3, 3)

        # Create shallow water domain
        domain1 = Domain(points, vertices, boundary)

        domain1.set_low_froude(0)
        domain1.reduction = mean
        domain1.smooth = False    # Exact result

        domain1.default_order = 2
        domain1.store = True
        domain1.set_datadir(tempfile.mkdtemp())
        domain1.set_name('spatio_temporal_boundary_source' + str(time.time()))

        # FIXME: This is extremely important!
        # How can we test if they weren't stored?


        # Bed-slope and friction at vertices (and interpolated elsewhere)
        domain1.set_quantity('elevation', 0)
        domain1.set_quantity('friction', 0)

        # Boundary conditions
        Br = Reflective_boundary(domain1)
        Bd = Dirichlet_boundary([0.3,0,0])
        domain1.set_boundary({'left': Bd, 'top': Bd, 'right': Br, 'bottom': Br})

        # Initial condition
        domain1.set_quantity('stage', 0)
        domain1.check_integrity()

        finaltime = 10

        # Evolution  (full domain - large steps)
        for t in domain1.evolve(yieldstep=0.671, finaltime=finaltime):
            pass

        cv1 = domain1.quantities['stage'].centroid_values

        # Create a triangle shaped domain (reusing coordinates from domain 1),
        # formed from the lower and right hand  boundaries and
        # the sw-ne diagonal
        # from domain 1. Call it domain2

        points = [[0,0], [1.0/3,0], [1.0/3,1.0/3],
                  [2.0/3,0], [2.0/3,1.0/3], [2.0/3,2.0/3],
                  [1,0], [1,1.0/3], [1,2.0/3], [1,1]]

        vertices = [[1,2,0], [3,4,1], [2,1,4], [4,5,2],
                    [6,7,3], [4,3,7], [7,8,4], [5,4,8], [8,9,5]]

        boundary = {(0,1): 'bottom',   (1,1): 'bottom',   (4,1): 'bottom',
                    (4,2): 'right',    (6,2): 'right',    (8,2): 'right',
                    (0,0): 'diagonal', (3,0): 'diagonal', (8,0): 'diagonal'}

        domain2 = Domain(points, vertices, boundary)

        domain2.low_froude = domain1.low_froude
        domain2.reduction = domain1.reduction
        domain2.smooth = False
        domain2.default_order = 2

        # Bed-slope and friction at vertices (and interpolated elsewhere)
        domain2.set_quantity('elevation', 0)
        domain2.set_quantity('friction', 0)
        domain2.set_quantity('stage', 0)

        # Boundary conditions
        Br = Reflective_boundary(domain2)
        Bf = Field_boundary(os.path.join(domain1.get_datadir(), domain1.get_name() + '.sww'), domain2)
        domain2.set_boundary({'right':Br, 'bottom':Br, 'diagonal':Bf})
        domain2.check_integrity()

        # Evolution (small steps)
        for t in domain2.evolve(yieldstep=0.0711, finaltime=finaltime):
            pass

        # Use output from domain1 as spatio-temporal boundary for domain2
        # and verify that results at right hand side are close.
        cv2 = domain2.quantities['stage'].centroid_values

        assert num.allclose(num.take(cv1, (0,8,16), axis=0),
                            num.take(cv2, (0,3,8), axis=0), rtol=1.0e-4)      # Diag
        assert num.allclose(num.take(cv1, (0,6,12), axis=0),
                            num.take(cv2, (0,1,4), axis=0), rtol=1.0e-4)      # Bottom
        assert num.allclose(num.take(cv1, (12,14,16), axis=0),
                            num.take(cv2, (4,6,8), axis=0), rtol=1.0e-4)      # RHS

        # Cleanup
        os.remove(os.path.join(domain1.get_datadir(), domain1.get_name() + '.sww'))
        os.remove(os.path.join(domain2.get_datadir(), domain2.get_name() + '.sww'))


    def test_spatio_temporal_boundary_2(self):
        """Test that boundary values can be read from file and interpolated
        in both time and space.
        This is a more basic test, verifying that boundary object
        produces the expected results
        """

        import time
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create sww file of simple propagation from left to right
        # through rectangular domain

        # Create basic mesh
        points, vertices, boundary = rectangular(3, 3)

        #Create shallow water domain
        domain1 = Domain(points, vertices, boundary)

        domain1.reduction = mean
        domain1.smooth = True    # To mimic MOST output

        domain1.default_order = 2
        domain1.store = True
        domain1.set_datadir(tempfile.mkdtemp())
        domain1.set_name('spatio_temporal_boundary_source' + str(time.time()))

        # FIXME: This is extremely important!
        # How can we test if they weren't stored?
        domain1.quantities_to_be_stored = {'elevation': 1,
                                           'stage': 2,
                                           'xmomentum': 2,
                                           'ymomentum': 2}

        # Bed-slope and friction at vertices (and interpolated elsewhere)
        domain1.set_quantity('elevation', 0)
        domain1.set_quantity('friction', 0)

        # Boundary conditions
        Br = Reflective_boundary(domain1)
        Bd = Dirichlet_boundary([0.3,0,0])
        domain1.set_boundary({'left': Bd, 'top': Bd, 'right': Br, 'bottom': Br})

        # Initial condition
        domain1.set_quantity('stage', 0)
        domain1.check_integrity()

        finaltime = 5

        # Evolution  (full domain - large steps)
        for t in domain1.evolve(yieldstep=1, finaltime=finaltime):
            pass

        # Create an triangle shaped domain (coinciding with some
        # coordinates from domain 1),
        # formed from the lower and right hand  boundaries and
        # the sw-ne diagonal
        # from domain 1. Call it domain2
        points = [[0,0],
                  [1.0/3,0], [1.0/3,1.0/3],
                  [2.0/3,0], [2.0/3,1.0/3], [2.0/3,2.0/3],
                  [1,0],     [1,1.0/3],     [1,2.0/3],     [1,1]]

        vertices = [[1,2,0], [3,4,1], [2,1,4], [4,5,2],
                    [6,7,3], [4,3,7], [7,8,4], [5,4,8], [8,9,5]]

        boundary = {(0,1): 'bottom',   (1,1): 'bottom',   (4,1): 'bottom',
                    (4,2): 'right',    (6,2): 'right',    (8,2): 'right',
                    (0,0): 'diagonal', (3,0): 'diagonal', (8,0): 'diagonal'}

        domain2 = Domain(points, vertices, boundary)

        domain2.reduction = domain1.reduction
        domain2.smooth = False
        domain2.default_order = 2

        # Bed-slope and friction at vertices (and interpolated elsewhere)
        domain2.set_quantity('elevation', 0)
        domain2.set_quantity('friction', 0)
        domain2.set_quantity('stage', 0)

        # Read results for specific timesteps t=1 and t=2
        fid = NetCDFFile(os.path.join(domain1.get_datadir(), domain1.get_name() + '.sww'))

        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        s1 = fid.variables['stage'][1,:]
        s2 = fid.variables['stage'][2,:]
        fid.close()

        shp = (len(x), 1)
        points = num.concatenate((num.reshape(x, shp), num.reshape(y, shp)),
                                 axis=1)

        # The diagonal points of domain 1 are 0, 5, 10, 15
        msg = ('value was\n%s\nshould be\n'
               '[[0,0], [1.0/3, 1.0/3],\n'
               '[2.0/3, 2.0/3], [1,1]]'
               % str(num.take(points, [0,5,10,15], axis=0)))
        assert num.allclose(num.take(points, [0,5,10,15], axis=0),
                            [[0,0], [1.0/3, 1.0/3], [2.0/3, 2.0/3], [1,1]]), msg

        # Boundary conditions
        Br = Reflective_boundary(domain2)
        Bf = Field_boundary(os.path.join(domain1.get_datadir(), domain1.get_name() + '.sww'),
                            domain2, verbose=False)
        domain2.set_boundary({'right':Br, 'bottom':Br, 'diagonal':Bf})
        domain2.check_integrity()

        # Test that interpolation points are the mid points of the all boundary
        # segments
        boundary_midpoints = [[1.0/6, 0], [1.0/2, 0], [5.0/6,0],
                              [1.0, 1.0/6], [1.0, 1.0/2], [1.0, 5.0/6],
                              [1.0/6, 1.0/6], [0.5, 0.5], [5.0/6, 5.0/6]]

        boundary_midpoints.sort()
        R = Bf.F.interpolation_points.tolist()
        R.sort()
        assert num.allclose(boundary_midpoints, R)

        # Check spatially interpolated output at time == 1
        domain2.set_time(1)

        # First diagonal midpoint
        R0 = Bf.evaluate(0, 0)
        assert num.allclose(R0[0], (s1[0] + s1[5]) / 2)

        # Second diagonal midpoint
        R0 = Bf.evaluate(3, 0)
        assert num.allclose(R0[0], (s1[5] + s1[10]) / 2)

        # First diagonal midpoint
        R0 = Bf.evaluate(8, 0)
        assert num.allclose(R0[0], (s1[10] + s1[15]) / 2)

        # Check spatially interpolated output at time == 2
        domain2.set_time(2)

        # First diagonal midpoint
        R0 = Bf.evaluate(0, 0)
        assert num.allclose(R0[0], (s2[0] + s2[5]) / 2)

        # Second diagonal midpoint
        R0 = Bf.evaluate(3, 0)
        assert num.allclose(R0[0], (s2[5] + s2[10]) / 2)

        # First diagonal midpoint
        R0 = Bf.evaluate(8, 0)
        assert num.allclose(R0[0], (s2[10] + s2[15]) / 2)

        # Now check temporal interpolation
        domain2.set_time(1 + 2.0/3)

        # First diagonal midpoint
        R0 = Bf.evaluate(0,0)
        assert num.allclose(R0[0],
                            ((s1[0] + s1[5]) / 2 + 2.0*(s2[0] + s2[5]) / 2) / 3)

        # Second diagonal midpoint
        R0 = Bf.evaluate(3, 0)
        assert num.allclose(R0[0],
                            ((s1[5] + s1[10]) / 2 + 2.0*(s2[5] + s2[10]) / 2) / 3)

        # First diagonal midpoint
        R0 = Bf.evaluate(8, 0)
        assert num.allclose(R0[0],
                            ((s1[10] + s1[15]) / 2 + 2.0*(s2[10] + s2[15]) / 2) / 3)

        # Cleanup
        os.remove(os.path.join(domain1.get_datadir(), domain1.get_name() + '.sww'))


    def test_spatio_temporal_boundary_3(self):
        """Test that boundary values can be read from file and interpolated
        in both time and space.
        This is a more basic test, verifying that boundary object
        produces the expected results

        This tests adjusting using mean_stage
        """

        #Create sww file of simple propagation from left to right
        #through rectangular domain

        import time
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        mean_stage = 5.2    # Adjust stage by this amount in boundary

        # Create basic mesh
        points, vertices, boundary = rectangular(3, 3)

        # Create shallow water domain
        domain1 = Domain(points, vertices, boundary)

        domain1.reduction = mean
        domain1.smooth = True    # To mimic MOST output

        domain1.default_order = 2
        domain1.store = True
        domain1.set_datadir(tempfile.mkdtemp())
        domain1.set_name('spatio_temporal_boundary_source' + str(time.time()))

        # FIXME: This is extremely important!
        # How can we test if they weren't stored?

        # Bed-slope and friction at vertices (and interpolated elsewhere)
        domain1.set_quantity('elevation', 0)
        domain1.set_quantity('friction', 0)

        # Boundary conditions
        Br = Reflective_boundary(domain1)
        Bd = Dirichlet_boundary([0.3, 0, 0])
        domain1.set_boundary({'left': Bd, 'top': Bd, 'right': Br, 'bottom': Br})

        # Initial condition
        domain1.set_quantity('stage', 0)
        domain1.check_integrity()

        finaltime = 5

        # Evolution  (full domain - large steps)
        for t in domain1.evolve(yieldstep=1, finaltime=finaltime):
            pass

        # Create an triangle shaped domain (coinciding with some
        # coordinates from domain 1),
        # formed from the lower and right hand  boundaries and
        # the sw-ne diagonal
        # from domain 1. Call it domain2
        points = [[0,0],
                  [1.0/3,0], [1.0/3,1.0/3],
                  [2.0/3,0], [2.0/3,1.0/3], [2.0/3,2.0/3],
                  [1,0],     [1,1.0/3],     [1,2.0/3],     [1,1]]

        vertices = [[1,2,0],
                    [3,4,1], [2,1,4], [4,5,2],
                    [6,7,3], [4,3,7], [7,8,4], [5,4,8], [8,9,5]]

        boundary = {(0,1): 'bottom',   (1,1): 'bottom',   (4,1): 'bottom',
                    (4,2): 'right',    (6,2): 'right',    (8,2): 'right',
                    (0,0): 'diagonal', (3,0): 'diagonal', (8,0): 'diagonal'}

        domain2 = Domain(points, vertices, boundary)

        domain2.reduction = domain1.reduction
        domain2.smooth = False
        domain2.default_order = 2

        # Bed-slope and friction at vertices (and interpolated elsewhere)
        domain2.set_quantity('elevation', 0)
        domain2.set_quantity('friction', 0)
        domain2.set_quantity('stage', 0)

        # Read results for specific timesteps t=1 and t=2
        fid = NetCDFFile(os.path.join(domain1.get_datadir(), domain1.get_name() + '.sww'))

        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        s1 = fid.variables['stage'][1,:]
        s2 = fid.variables['stage'][2,:]
        fid.close()

        shp = (len(x), 1)
        points = num.concatenate((num.reshape(x, shp), num.reshape(y, shp)),
                                 axis=1)
        #The diagonal points of domain 1 are 0, 5, 10, 15

        msg = ('values was\n%s\nshould be\n'
               '[[0,0], [1.0/3, 1.0/3],\n'
               '[2.0/3, 2.0/3], [1,1]]'
               % str(num.take(points, [0,5,10,15], axis=0)))
        assert num.allclose(num.take(points, [0,5,10,15], axis=0),
                            [[0,0], [1.0/3, 1.0/3], [2.0/3, 2.0/3], [1,1]]), msg

        # Boundary conditions
        Br = Reflective_boundary(domain2)
        Bf = Field_boundary(os.path.join(domain1.get_datadir(), domain1.get_name() + '.sww'),
                            domain2, mean_stage=mean_stage, verbose=False)

        domain2.set_boundary({'right':Br, 'bottom':Br, 'diagonal':Bf})
        domain2.check_integrity()

        # Test that interpolation points are the mid points of the all boundary
        # segments
        boundary_midpoints = [[1.0/6, 0], [1.0/2, 0], [5.0/6,0],
                              [1.0, 1.0/6], [1.0, 1.0/2], [1.0, 5.0/6],
                              [1.0/6, 1.0/6], [0.5, 0.5], [5.0/6, 5.0/6]]

        boundary_midpoints.sort()
        R = Bf.F.interpolation_points.tolist()
        R.sort()
        assert num.allclose(boundary_midpoints, R)

        # Check spatially interpolated output at time == 1
        domain2.set_time(1)

        # First diagonal midpoint
        R0 = Bf.evaluate(0, 0)
        assert num.allclose(R0[0], (s1[0] + s1[5]) / 2 + mean_stage)

        # Second diagonal midpoint
        R0 = Bf.evaluate(3, 0)
        assert num.allclose(R0[0], (s1[5] + s1[10]) / 2 + mean_stage)

        # First diagonal midpoint
        R0 = Bf.evaluate(8, 0)
        assert num.allclose(R0[0], (s1[10] + s1[15]) / 2 + mean_stage)

        # Check spatially interpolated output at time == 2
        domain2.set_time(2)

        # First diagonal midpoint
        R0 = Bf.evaluate(0, 0)
        assert num.allclose(R0[0], (s2[0] + s2[5]) / 2 + mean_stage)

        # Second diagonal midpoint
        R0 = Bf.evaluate(3, 0)
        assert num.allclose(R0[0], (s2[5] + s2[10]) / 2 + mean_stage)

        # First diagonal midpoint
        R0 = Bf.evaluate(8, 0)
        assert num.allclose(R0[0], (s2[10] + s2[15]) / 2 + mean_stage)

        #Now check temporal interpolation
        domain2.set_time(1 + 2.0/3)

        # First diagonal midpoint
        R0 = Bf.evaluate(0, 0)
        assert num.allclose(R0[0],
                            ((s1[0] + s1[5]) / 2 + 2.0*(s2[0] + s2[5]) / 2) / 3 +
                                mean_stage)

        # Second diagonal midpoint
        R0 = Bf.evaluate(3, 0)
        assert num.allclose(R0[0],
                            ((s1[5] + s1[10]) / 2 + 2.0*(s2[5] + s2[10]) / 2) / 3 +
                                mean_stage)

        # First diagonal midpoint
        R0 = Bf.evaluate(8, 0)
        assert num.allclose(R0[0],
                            ((s1[10] + s1[15]) / 2 + 2.0*(s2[10] + s2[15]) / 2) / 3 +
                                mean_stage)

        # Cleanup
        os.remove(os.path.join(domain1.get_datadir(), domain1.get_name() + '.sww'))


    def test_spatio_temporal_boundary_outside(self):
        """Test that field_boundary catches if a point is outside the sww
        that defines it
        """

        import time
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create sww file of simple propagation from left to right
        # through rectangular domain

        # Create basic mesh
        points, vertices, boundary = rectangular(3, 3)

        # Create shallow water domain
        domain1 = Domain(points, vertices, boundary)

        domain1.reduction = mean
        domain1.smooth = True    # To mimic MOST output

        domain1.default_order = 2
        domain1.store = True
        domain1.set_datadir(tempfile.mkdtemp())
        domain1.set_name('spatio_temporal_boundary_source' + str(time.time()))

        # FIXME: This is extremely important!
        # How can we test if they weren't stored?
        domain1.set_quantities_to_be_stored({'elevation': 1,
                                             'stage': 2,
                                             'xmomentum': 2,
                                             'ymomentum': 2})


        # Bed-slope and friction at vertices (and interpolated elsewhere)
        domain1.set_quantity('elevation', 0)
        domain1.set_quantity('friction', 0)

        # Boundary conditions
        Br = Reflective_boundary(domain1)
        Bd = Dirichlet_boundary([0.3, 0, 0])
        domain1.set_boundary({'left': Bd, 'top': Bd, 'right': Br, 'bottom': Br})

        # Initial condition
        domain1.set_quantity('stage', 0)
        domain1.check_integrity()

        finaltime = 5

        # Evolution  (full domain - large steps)
        for t in domain1.evolve(yieldstep=1, finaltime=finaltime):
            pass

        # Create an triangle shaped domain (coinciding with some
        # coordinates from domain 1, but one edge outside!),
        # formed from the lower and right hand  boundaries and
        # the sw-ne diagonal as in the previous test but scaled
        # in the x direction by a factor of 2
        points = [[0,0],
                  [2.0/3,0], [2.0/3,1.0/3],
                  [4.0/3,0], [4.0/3,1.0/3], [4.0/3,2.0/3],
                  [2,0],     [2,1.0/3],     [2,2.0/3],     [2,1]]

        vertices = [[1,2,0],
                    [3,4,1], [2,1,4], [4,5,2],
                    [6,7,3], [4,3,7], [7,8,4], [5,4,8], [8,9,5]]

        boundary = {(0,1): 'bottom',   (1,1): 'bottom',   (4,1): 'bottom',
                    (4,2): 'right',    (6,2): 'right',    (8,2): 'right',
                    (0,0): 'diagonal', (3,0): 'diagonal', (8,0): 'diagonal'}

        domain2 = Domain(points, vertices, boundary)

        domain2.reduction = domain1.reduction
        domain2.smooth = False
        domain2.default_order = 2

        # Bed-slope and friction at vertices (and interpolated elsewhere)
        domain2.set_quantity('elevation', 0)
        domain2.set_quantity('friction', 0)
        domain2.set_quantity('stage', 0)

        # Read results for specific timesteps t=1 and t=2
        fid = NetCDFFile(os.path.join(domain1.get_datadir(), domain1.get_name() + '.sww'))

        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        s1 = fid.variables['stage'][1,:]
        s2 = fid.variables['stage'][2,:]
        fid.close()

        shp = (len(x), 1)
        points = num.concatenate((num.reshape(x, shp), num.reshape(y, shp)),
                                 axis=1)
        #The diagonal points of domain 1 are 0, 5, 10, 15

        assert num.allclose(num.take(points, [0,5,10,15], axis=0),
                            [[0,0], [1.0/3,1.0/3], [2.0/3,2.0/3], [1,1]])

        # Boundary conditions
        Br = Reflective_boundary(domain2)
        Bf = Field_boundary(os.path.join(domain1.get_datadir(), domain1.get_name() + '.sww'),
                            domain2, mean_stage=1, verbose=False)

        domain2.set_boundary({'right':Br, 'bottom':Br, 'diagonal':Bf})
        domain2.check_integrity()

        try:
            for t in domain2.evolve(yieldstep=1, finaltime=finaltime):
                pass
        except Exception:
            pass
        else:
            msg = 'This should have caught NAN at boundary'
            raise Exception(msg)

        #Cleanup
        os.remove(os.path.join(domain1.get_datadir(), domain1.get_name() + '.sww'))



class Test_Prescribed_Boundary_Limiter(unittest.TestCase):
    """The edge limiter on boundary edges.

    A boundary cell has no neighbour beyond its boundary edge, so the end of a
    smooth slope looks like a local extremum. A value the boundary condition
    prescribes stands in for that neighbour: the edge may move towards it but
    not past it. Anything else (walls, transmissive copies, custom classes) has
    no outside information and keeps the plain limiter.
    """

    def _domain(self, make_boundary, nx=10, ny=5, length=10.0, width=5.0):
        points, vertices, boundary = rectangular_cross(nx, ny, len1=length, len2=width)
        domain = Domain(points, vertices, boundary)
        domain.set_flow_algorithm('DE1')
        domain.set_compute_mode('legacy')   # inspects the host edge arrays
        domain.set_low_froude(0)
        domain.set_store(False)
        domain.set_quantity('elevation', 0.0)
        domain.set_quantity('friction', 0.0)
        domain.set_quantity('stage', 1.0)
        wall = Reflective_boundary(domain)
        domain.set_boundary({'left': make_boundary(domain), 'right': wall,
                             'top': wall, 'bottom': wall})
        return domain

    def _end_cell(self, domain):
        """The left-end triangle in the middle row, and its boundary edge."""
        centroids = domain.centroid_coordinates
        for k in range(len(domain)):
            for edge in range(3):
                if (domain.neighbours[k, edge] < 0 and domain.boundary[(k, edge)] == 'left'
                        and abs(centroids[k, 1] - 2.5) < 1.0e-9):
                    return k, edge
        raise AssertionError('no end cell found')

    def _extrapolate(self, domain):
        # The limiter reads the boundary values, and the copying boundaries read
        # the edge values, so extrapolate on both sides of the boundary update.
        domain.distribute_to_vertices_and_edges()
        domain.update_boundary()
        domain.distribute_to_vertices_and_edges()

    def test_linear_profile_reaches_the_prescribed_value(self):
        for name, make in (
            ('Dirichlet', lambda d: Dirichlet_boundary([1.0, 0.0, 0.0])),
            ('Time', lambda d: Time_boundary(d, function=lambda t: [1.0, 0.0, 0.0])),
            ('set stage', lambda d: Transmissive_n_momentum_zero_t_momentum_set_stage_boundary(
                d, function=lambda t: 1.0)),
        ):
            with self.subTest(name):
                domain = self._domain(make)
                domain.set_quantity('stage', lambda x, y: 1.0 + 0.01 * x)
                self._extrapolate(domain)
                k, edge = self._end_cell(domain)
                stage = domain.quantities['stage']
                self.assertGreater(stage.centroid_values[k], 1.001)
                self.assertAlmostEqual(stage.edge_values[k, edge], 1.0, places=6)

    def test_uniform_depth_over_a_sloping_bed_reaches_the_prescribed_depth(self):
        # The depth at the boundary uses the bed at the edge, not at the centroid.
        domain = self._domain(lambda d: Dirichlet_boundary([1.0, 0.0, 0.0]))
        domain.set_quantity('elevation', lambda x, y: 0.1 * x)
        domain.set_quantity('stage', lambda x, y: 1.0 + 0.1 * x)
        self._extrapolate(domain)
        k, edge = self._end_cell(domain)
        self.assertAlmostEqual(domain.quantities['stage'].edge_values[k, edge], 1.0, places=6)
        self.assertAlmostEqual(domain.quantities['height'].edge_values[k, edge], 1.0, places=6)

    def test_edge_stops_at_the_prescribed_value(self):
        # The slope would carry the edge to 1.0, but the boundary says 1.001.
        domain = self._domain(lambda d: Dirichlet_boundary([1.001, 0.0, 0.0]))
        domain.set_quantity('stage', lambda x, y: 1.0 + 0.01 * x)
        self._extrapolate(domain)
        k, edge = self._end_cell(domain)
        self.assertAlmostEqual(domain.quantities['stage'].edge_values[k, edge], 1.001, places=6)

    def test_local_peak_is_not_extrapolated_away_from_the_boundary_value(self):
        # The cell stands above its neighbours and the boundary value. Following
        # the gradient would carry the edge further above them both, which feeds
        # back through the boundary flux and runs away at an outflow.
        for name, make in (
            ('Dirichlet', lambda d: Dirichlet_boundary([1.0, 0.0, 0.0])),
            ('set stage', lambda d: Transmissive_n_momentum_zero_t_momentum_set_stage_boundary(
                d, function=lambda t: 1.0)),
        ):
            with self.subTest(name):
                domain = self._domain(make)
                k, edge = self._end_cell(domain)
                stage = domain.quantities['stage']
                stage.centroid_values[k] = 1.2
                self._extrapolate(domain)
                self.assertTrue(num.allclose(stage.edge_values[k], 1.2, atol=1.0e-5))

    def test_a_boundary_without_outside_information_keeps_the_plain_limiter(self):
        # A wall or a transmissive copy cannot tell a slope from a peak, so both
        # stay flat. Relaxing this overshoots beside walls and runs away at outflows.
        for name, make in (('wall', lambda d: Reflective_boundary(d)),
                           ('transmissive', lambda d: Transmissive_boundary(d))):
            with self.subTest(name + ' peak'):
                domain = self._domain(make)
                k, edge = self._end_cell(domain)
                xmom = domain.quantities['xmomentum']
                xmom.centroid_values[k] = 1.5
                self._extrapolate(domain)
                self.assertTrue(num.allclose(xmom.edge_values[k], 1.5, atol=1.0e-10))
            with self.subTest(name + ' slope'):
                domain = self._domain(make)
                domain.set_quantity('stage', lambda x, y: 1.0 + 0.01 * x)
                self._extrapolate(domain)
                k, edge = self._end_cell(domain)
                stage = domain.quantities['stage']
                self.assertTrue(num.allclose(stage.edge_values[k], stage.centroid_values[k],
                                             atol=1.0e-10))

    def test_round_off_does_not_collapse_the_gradient(self):
        # The y velocity is 0 at the cell and at the boundary edge, with a real
        # gradient along the boundary. Round-off in that zero must not read as an
        # overshoot of a zero-width bound, which would flatten every edge and make
        # results differ between machines.
        def edges(noise):
            domain = self._domain(lambda d: Dirichlet_boundary([1.0, 0.0, 0.0]))
            domain.set_quantity('ymomentum', lambda x, y: 0.1 * (y - 2.5))
            k, edge = self._end_cell(domain)
            domain.quantities['ymomentum'].centroid_values[k] = noise
            self._extrapolate(domain)
            return domain.quantities['ymomentum'].edge_values[k].copy()

        clean = edges(0.0)
        self.assertGreater(num.abs(clean).max(), 0.01)
        for noise in (1.0e-16, -1.0e-15, 1.0e-13):
            self.assertTrue(num.allclose(edges(noise), clean, atol=1.0e-8), noise)

    def test_tiny_bound_does_not_amplify_round_off(self):
        # The y velocity nearly vanishes at the cell and the boundary, so its bound
        # is tiny, and its gradient runs almost along the boundary, so the offset of
        # the boundary edge is tiny as well. Their ratio must not set the gradient:
        # that turned 1e-11 into 3e-6 on the other edges, enough to make results
        # differ visibly between the CPU and GPU paths.
        def edges(value):
            domain = self._domain(lambda d: Dirichlet_boundary([1.0, 0.0, 0.0]))
            domain.set_quantity('ymomentum', lambda x, y: 0.1 * (y - 2.5) + 1.0e-7 * x)
            k, edge = self._end_cell(domain)
            domain.quantities['ymomentum'].centroid_values[k] = value
            self._extrapolate(domain)
            return domain.quantities['ymomentum'].edge_values[k].copy()

        worst = 0.0
        for value in num.linspace(-2.0e-6, 2.0e-6, 41):
            eps = 1.0e-11
            worst = max(worst, num.abs(edges(value + eps) - edges(value)).max() / eps)
        self.assertLess(worst, 10.0)

    def test_limit_mask_follows_the_boundary_class(self):
        class Custom(Boundary):
            def evaluate(self, vol_id=None, edge_id=None):
                return [1.0, 0.0, 0.0]

        for name, make, stage_set, momentum_set in (
            ('Dirichlet', lambda d: Dirichlet_boundary([1.0, 0.0, 0.0]), True, True),
            ('Time', lambda d: Time_boundary(d, function=lambda t: [1.0, 0.0, 0.0]), True, True),
            ('set stage', lambda d: Transmissive_momentum_set_stage_boundary(
                d, function=lambda t: 1.0), True, False),
            ('n set stage', lambda d: Transmissive_n_momentum_zero_t_momentum_set_stage_boundary(
                d, function=lambda t: 1.0), True, False),
            ('wall', lambda d: Reflective_boundary(d), False, False),
            ('transmissive', lambda d: Transmissive_boundary(d), False, False),
            ('custom', lambda d: Custom(), False, False),
        ):
            with self.subTest(name):
                domain = self._domain(make)
                k, edge = self._end_cell(domain)
                expected = (1 << edge if stage_set else 0) | (8 << edge if momentum_set else 0)
                self.assertEqual(domain.boundary_limit_mask[k], expected)

    def test_limit_mask_is_rebuilt_by_set_boundary(self):
        domain = self._domain(lambda d: Dirichlet_boundary([1.0, 0.0, 0.0]))
        self.assertGreater(num.count_nonzero(domain.boundary_limit_mask), 0)
        wall = Reflective_boundary(domain)
        domain.set_boundary({'left': wall, 'right': wall, 'top': wall, 'bottom': wall})
        self.assertEqual(num.count_nonzero(domain.boundary_limit_mask), 0)

    def test_first_step_after_set_boundary_sees_the_new_values(self):
        # A continued evolve() skips the initial boundary update, and the limiter
        # must not bound the first step by the boundary that has just been replaced.
        domain = self._domain(lambda d: Reflective_boundary(d))
        for t in domain.evolve(yieldstep=1.0, finaltime=1.0):
            pass
        wall = Reflective_boundary(domain)
        domain.set_boundary({'left': Dirichlet_boundary([0.2, 0.0, 0.0]), 'right': wall,
                             'top': wall, 'bottom': wall})
        seen = []
        original = domain.evolve_one_rk2_step

        def step(*args, **kwargs):
            seen.append(domain.quantities['stage'].boundary_values.copy())
            return original(*args, **kwargs)

        domain.evolve_one_rk2_step = step
        for t in domain.evolve(yieldstep=0.1, finaltime=1.1):
            pass
        k, edge = self._end_cell(domain)
        self.assertAlmostEqual(seen[0][-domain.neighbours[k, edge] - 1], 0.2)

    def _max_speed(self, domain):
        quantities = domain.quantities
        depth = quantities['stage'].centroid_values - quantities['elevation'].centroid_values
        wet = depth > 0.01
        speed = num.hypot(quantities['xmomentum'].centroid_values,
                          quantities['ymomentum'].centroid_values)[wet] / depth[wet]
        return speed.max()

    def _unstructured(self, polygon, area):
        # Unstructured on purpose: the instabilities below do not appear on the
        # symmetric stencils of a regular mesh.
        domain = create_domain_from_regions(
            polygon, boundary_tags={'left': [0], 'top': [1], 'right': [2], 'bottom': [3]},
            maximum_triangle_area=area, verbose=False)
        domain.set_flow_algorithm('DE1')
        domain.set_compute_mode('legacy')
        domain.set_low_froude(0)
        domain.set_store(False)
        return domain

    def test_outflow_at_a_set_stage_boundary_stays_stable(self):
        # A sloping channel with banks, draining through a set-stage outflow. An
        # edge extrapolated against the boundary value feeds back through the
        # outflow flux: this collapses the timestep within seconds.
        length, width, slope = 40.0, 14.0, 1.0 / 300.0
        domain = self._unstructured([[0, 0], [0, length], [width, length], [width, 0]], 0.2)

        def elevation(x, y):
            return (-slope * y - 1.0 * ((x > 2) & (x < 12))
                    + (1.0 - 0.5 * (x - 2)) * ((x > 2) & (x < 4))
                    + (1.0 + 0.5 * (x - 12)) * ((x > 10) & (x < 12)))

        domain.set_quantity('elevation', elevation)
        domain.set_quantity('friction', 0.03)
        domain.set_quantity('stage', lambda x, y: -slope * y - 0.35)
        wall = Reflective_boundary(domain)
        outflow = Transmissive_n_momentum_zero_t_momentum_set_stage_boundary(
            domain, function=lambda t: -slope * length - 0.35)
        domain.set_boundary({'left': wall, 'right': wall, 'bottom': wall, 'top': outflow})
        for t in domain.evolve(yieldstep=2.0, finaltime=30.0):
            self.assertLess(self._max_speed(domain), 2.0)

    def test_dam_break_beside_walls_does_not_overshoot(self):
        # A closed box has no prescribed value, so nothing may change here. Loosening
        # the bound at walls raises the peak speed from 3.3 to above 4.2 m/s.
        domain = self._unstructured([[0, 0], [0, 5], [20, 5], [20, 0]], 0.02)
        domain.set_quantity('elevation', 0.0)
        domain.set_quantity('friction', 0.0)
        domain.set_quantity('stage', lambda x, y: num.where(x + 0.4 * y < 6.0, 1.0, 0.1))
        wall = Reflective_boundary(domain)
        domain.set_boundary({'left': wall, 'right': wall, 'top': wall, 'bottom': wall})
        for t in domain.evolve(yieldstep=0.25, finaltime=10.0):
            self.assertLess(self._max_speed(domain), 3.6)


if __name__ == "__main__":
    unittest.main()

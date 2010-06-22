#!/usr/bin/env python

import unittest, os
import os.path
from math import pi, sqrt
import tempfile

import anuga

from anuga.config import g, epsilon
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a
from anuga.utilities.numerical_tools import mean
from anuga.geometry.polygon import is_inside_polygon
from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.abstract_2d_finite_volumes.quantity import Quantity
from anuga.geospatial_data.geospatial_data import Geospatial_data
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross

from anuga.utilities.system_tools import get_pathname_from_package
from swb_domain import *

import numpy as num

# Get gateway to C implementation of flux function for direct testing
from shallow_water_ext import flux_function_central as flux_function




class Test_swb_boundary_condition(unittest.TestCase):
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

        D = anuga.Dirichlet_boundary([5, 2, 1])
        T = anuga.Transmissive_boundary(domain)
        R = anuga.Reflective_boundary(domain)
        domain.set_boundary({'First': D, 'Second': T, 'Third': R})

        domain.set_centroid_transmissive_bc(False)
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
                    (0, 1): 'Internal'}

        domain = anuga.Domain(points, vertices, boundary)
        domain.check_integrity()

        domain.set_quantity('stage', [[1,2,3], [5,5,5], [0,0,9], [-6,3,3]])

        domain.set_quantity('xmomentum', [[1,1,1], [2,2,2], [3,3,3], [4,4,4]])

        domain.set_quantity('ymomentum', [[10,10,10], [20,20,20],
                                          [30,30,30], [40,40,40]])

        D = anuga.Dirichlet_boundary([5, 2, 1])
        T = anuga.Transmissive_boundary(domain)
        R = anuga.Reflective_boundary(domain)
        domain.set_boundary({'First': D, 'Second': T,
                             'Third': R, 'Internal': None})

        domain.update_boundary()
        domain.check_integrity()

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

        D = anuga.Dirichlet_boundary([5, 2, 1])
        T = anuga.Transmissive_stage_zero_momentum_boundary(domain)
        R = anuga.Reflective_boundary(domain)
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
        Bt = anuga.Time_boundary(domain=domain, f=lambda t: [t, 0.0, 0.0])

        Br = anuga.Reflective_boundary(domain)              # Reflective wall

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

        import time
        from mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(3, 3)

        # Create shallow water domain
        domain1 = Domain(points, vertices, boundary)

        domain1.reduction = mean
        domain1.smooth = False    # Exact result

        domain1.default_order = 2
        domain1.store = True
        domain1.set_datadir('.')
        domain1.set_name('spatio_temporal_boundary_source' + str(time.time()))

        # FIXME: This is extremely important!
        # How can we test if they weren't stored?


        # Bed-slope and friction at vertices (and interpolated elsewhere)
        domain1.set_quantity('elevation', 0)
        domain1.set_quantity('friction', 0)

        # Boundary conditions
        Br = anuga.Reflective_boundary(domain1)
        Bd = anuga.Dirichlet_boundary([0.3,0,0])
        domain1.set_boundary({'left': Bd, 'top': Bd, 'right': Br, 'bottom': Br})

        # Initial condition
        domain1.set_quantity('stage', 0)
        domain1.check_integrity()

        finaltime = 5

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

        domain2.reduction = domain1.reduction
        domain2.smooth = False
        domain2.default_order = 2

        # Bed-slope and friction at vertices (and interpolated elsewhere)
        domain2.set_quantity('elevation', 0)
        domain2.set_quantity('friction', 0)
        domain2.set_quantity('stage', 0)

        # Boundary conditions
        Br = anuga.Reflective_boundary(domain2)
        Bf = anuga.Field_boundary(domain1.get_name() + '.sww', domain2)
        domain2.set_boundary({'right':Br, 'bottom':Br, 'diagonal':Bf})
        domain2.check_integrity()

        # Evolution (small steps)
        for t in domain2.evolve(yieldstep=0.0711, finaltime=finaltime):
            pass

        # Use output from domain1 as spatio-temporal boundary for domain2
        # and verify that results at right hand side are close.
        cv2 = domain2.quantities['stage'].centroid_values


        assert num.allclose(num.take(cv1, (0,8,16), axis=0),
                            num.take(cv2, (0,3,8), axis=0),atol=1.0e-2)      # Diag
        assert num.allclose(num.take(cv1, (0,6,12), axis=0),
                            num.take(cv2, (0,1,4), axis=0),atol=1.0e-2)      # Bottom
        assert num.allclose(num.take(cv1, (12,14,16), axis=0),
                            num.take(cv2, (4,6,8), axis=0),atol=1.0e-2)      # RHS

        # Cleanup
        os.remove(domain1.get_name() + '.sww')

    def test_spatio_temporal_boundary_2(self):
        """Test that boundary values can be read from file and interpolated
        in both time and space.
        This is a more basic test, verifying that boundary object
        produces the expected results
        """

        import time
        from mesh_factory import rectangular

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
        domain1.set_datadir('.')
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
        Br = anuga.Reflective_boundary(domain1)
        Bd = anuga.Dirichlet_boundary([0.3,0,0])
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
        from Scientific.IO.NetCDF import NetCDFFile
        fid = NetCDFFile(domain1.get_name() + '.sww')

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
        Br = anuga.Reflective_boundary(domain2)
        Bf = anuga.Field_boundary(domain1.get_name() + '.sww',
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
        domain2.time = 1

        # First diagonal midpoint
        R0 = Bf.evaluate(0, 0)
        assert num.allclose(R0[0], (s1[0] + s1[5])/2)

        # Second diagonal midpoint
        R0 = Bf.evaluate(3, 0)
        assert num.allclose(R0[0], (s1[5] + s1[10])/2)

        # First diagonal midpoint
        R0 = Bf.evaluate(8, 0)
        assert num.allclose(R0[0], (s1[10] + s1[15])/2)

        # Check spatially interpolated output at time == 2
        domain2.time = 2

        # First diagonal midpoint
        R0 = Bf.evaluate(0, 0)
        assert num.allclose(R0[0], (s2[0] + s2[5])/2)

        # Second diagonal midpoint
        R0 = Bf.evaluate(3, 0)
        assert num.allclose(R0[0], (s2[5] + s2[10])/2)

        # First diagonal midpoint
        R0 = Bf.evaluate(8, 0)
        assert num.allclose(R0[0], (s2[10] + s2[15])/2)

        # Now check temporal interpolation
        domain2.time = 1 + 2.0/3

        # First diagonal midpoint
        R0 = Bf.evaluate(0,0)
        assert num.allclose(R0[0],
                            ((s1[0] + s1[5])/2 + 2.0*(s2[0] + s2[5])/2)/3)

        # Second diagonal midpoint
        R0 = Bf.evaluate(3, 0)
        assert num.allclose(R0[0],
                            ((s1[5] + s1[10])/2 + 2.0*(s2[5] + s2[10])/2)/3)

        # First diagonal midpoint
        R0 = Bf.evaluate(8, 0)
        assert num.allclose(R0[0],
                            ((s1[10] + s1[15])/2 + 2.0*(s2[10] + s2[15])/2)/3)

        # Cleanup
        os.remove(domain1.get_name() + '.sww')

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
        from mesh_factory import rectangular

        mean_stage = 5.2    # Adjust stage by this amount in boundary

        # Create basic mesh
        points, vertices, boundary = rectangular(3, 3)

        # Create shallow water domain
        domain1 = Domain(points, vertices, boundary)

        domain1.reduction = mean
        domain1.smooth = True    # To mimic MOST output

        domain1.default_order = 2
        domain1.store = True
        domain1.set_datadir('.')
        domain1.set_name('spatio_temporal_boundary_source' + str(time.time()))

        # FIXME: This is extremely important!
        # How can we test if they weren't stored?

        # Bed-slope and friction at vertices (and interpolated elsewhere)
        domain1.set_quantity('elevation', 0)
        domain1.set_quantity('friction', 0)

        # Boundary conditions
        Br = anuga.Reflective_boundary(domain1)
        Bd = anuga.Dirichlet_boundary([0.3, 0, 0])
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
        from Scientific.IO.NetCDF import NetCDFFile
        fid = NetCDFFile(domain1.get_name() + '.sww')

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
        Br = anuga.Reflective_boundary(domain2)
        Bf = anuga.Field_boundary(domain1.get_name() + '.sww',
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
        domain2.time = 1

        # First diagonal midpoint
        R0 = Bf.evaluate(0, 0)
        assert num.allclose(R0[0], (s1[0] + s1[5])/2 + mean_stage)

        # Second diagonal midpoint
        R0 = Bf.evaluate(3, 0)
        assert num.allclose(R0[0], (s1[5] + s1[10])/2 + mean_stage)

        # First diagonal midpoint
        R0 = Bf.evaluate(8, 0)
        assert num.allclose(R0[0], (s1[10] + s1[15])/2 + mean_stage)

        # Check spatially interpolated output at time == 2
        domain2.time = 2

        # First diagonal midpoint
        R0 = Bf.evaluate(0, 0)
        assert num.allclose(R0[0], (s2[0] + s2[5])/2 + mean_stage)

        # Second diagonal midpoint
        R0 = Bf.evaluate(3, 0)
        assert num.allclose(R0[0], (s2[5] + s2[10])/2 + mean_stage)

        # First diagonal midpoint
        R0 = Bf.evaluate(8, 0)
        assert num.allclose(R0[0], (s2[10] + s2[15])/2 + mean_stage)

        #Now check temporal interpolation
        domain2.time = 1 + 2.0/3

        # First diagonal midpoint
        R0 = Bf.evaluate(0, 0)
        assert num.allclose(R0[0],
                            ((s1[0] + s1[5])/2 + 2.0*(s2[0] + s2[5])/2)/3 +
                                mean_stage)

        # Second diagonal midpoint
        R0 = Bf.evaluate(3, 0)
        assert num.allclose(R0[0],
                            ((s1[5] + s1[10])/2 + 2.0*(s2[5] + s2[10])/2)/3 +
                                mean_stage)

        # First diagonal midpoint
        R0 = Bf.evaluate(8, 0)
        assert num.allclose(R0[0],
                            ((s1[10] + s1[15])/2 + 2.0*(s2[10] + s2[15])/2)/3 +
                                mean_stage)

        # Cleanup
        os.remove(domain1.get_name() + '.sww')
        os.remove(domain2.get_name() + '.sww')

        

    def test_spatio_temporal_boundary_outside(self):
        """Test that field_boundary catches if a point is outside the sww
        that defines it
        """

        import time
        from mesh_factory import rectangular

        # Create sww file of simple propagation from left to right
        # through rectangular domain

        # Create basic mesh
        points, vertices, boundary = rectangular(3, 3)

        # Create shallow water domain
        domain1 = anuga.Domain(points, vertices, boundary)

        domain1.reduction = mean
        domain1.smooth = True    # To mimic MOST output

        domain1.default_order = 2
        domain1.store = True
        domain1.set_datadir('.')
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
        Br = anuga.Reflective_boundary(domain1)
        Bd = anuga.Dirichlet_boundary([0.3, 0, 0])
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
        from Scientific.IO.NetCDF import NetCDFFile
        fid = NetCDFFile(domain1.get_name() + '.sww')

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
        Br = anuga.Reflective_boundary(domain2)
        Bf = anuga.Field_boundary(domain1.get_name() + '.sww',
                            domain2, mean_stage=1, verbose=False)

        domain2.set_boundary({'right':Br, 'bottom':Br, 'diagonal':Bf})
        domain2.check_integrity()

        try:
            for t in domain2.evolve(yieldstep=1, finaltime=finaltime):
                pass
        except:
            pass
        else:
            msg = 'This should have caught NAN at boundary'
            raise Exception, msg

        #Cleanup
        os.remove(domain1.get_name() + '.sww')
        os.remove(domain2.get_name() + '.sww')        



#################################################################################

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_swb_boundary_condition, 'test_boundary_condition_time')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)

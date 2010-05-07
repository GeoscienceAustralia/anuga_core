#!/usr/bin/env python

import unittest, os
import os.path
from math import pi, sqrt
import tempfile

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




class Test_swb_clean(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass


    def test_balance_deep_and_shallow(self):
        """Test that balanced limiters preserve conserved quantites.
        This test is using old depth based balanced limiters
        """

        import copy

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #             bac,     bce,     ecf,     dbe
        elements = [[1,0,2], [1,2,4], [4,2,5], [3,1,4] ]

        domain = Domain(points, elements)
        domain.check_integrity()

        # Create a deliberate overshoot
        domain.set_quantity('stage', [[3,0,3], [2,2,6], [5,3,8], [8,3,5]])
        domain.set_quantity('elevation', 0)    # Flat bed
        stage = domain.quantities['stage']

        ref_centroid_values = copy.copy(stage.centroid_values[:])    # Copy

        # Limit
        domain.distribute_to_vertices_and_edges()

        # Assert that quantities are conserved
        for k in range(len(domain)):
            assert num.allclose(ref_centroid_values[k],
                                num.sum(stage.vertex_values[k,:])/3)

        # Now try with a non-flat bed - closely hugging initial stage in places
        domain.set_quantity('stage', [[3,0,3], [2,2,6], [5,3,8], [8,3,5]])
        domain.set_quantity('elevation', [[0,0,0],
                                          [1.8,1.9,5.9],
                                          [4.6,0,0],
                                          [0,2,4]])
        stage = domain.quantities['stage']
        elevation = domain.quantities['elevation']
        height = domain.quantities['height']

        ref_centroid_values = copy.copy(stage.centroid_values[:])    # Copy
        ref_vertex_values = copy.copy(stage.vertex_values[:])        # Copy

        # Limit
        domain.distribute_to_vertices_and_edges()

        # Assert that all vertex quantities have changed
        for k in range(len(domain)):
            assert not num.allclose(ref_vertex_values[k,:],
                                    stage.vertex_values[k,:])
        # and assert that quantities are still conserved
        for k in range(len(domain)):
            assert num.allclose(ref_centroid_values[k],
                                num.sum(stage.vertex_values[k,:])/3)

        # Check actual results
        assert num.allclose(stage.vertex_values,
                            [[ 2.66666667,  0.66666667,  2.66666667],
                             [ 3.33333333,  3.33333333,  3.33333333],
                             [ 3.73333333,  4.93333333,  7.33333333],
                             [ 7.33333333,  4.93333333,  3.73333333]])



    def test_balance_deep_and_shallow_tight_SL(self):
        """Test that balanced limiters preserve conserved quantites.
        This test is using Tight Slope Limiters
        """

        import copy

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #             bac,     bce,     ecf,     dbe
        elements = [[1,0,2], [1,2,4], [4,2,5], [3,1,4] ]

        domain = Domain(points, elements)
        domain.check_integrity()

        # Create a deliberate overshoot
        domain.set_quantity('stage', [[3,0,3], [2,2,6], [5,3,8], [8,3,5]])
        domain.set_quantity('elevation', 0)    # Flat bed
        stage = domain.quantities['stage']

        ref_centroid_values = copy.copy(stage.centroid_values[:])    # Copy

        # Limit
        domain.tight_slope_limiters = 1
        domain.distribute_to_vertices_and_edges()

        # Assert that quantities are conserved
        for k in range(len(domain)):
            assert num.allclose (ref_centroid_values[k],
                                 num.sum(stage.vertex_values[k,:])/3)

        # Now try with a non-flat bed - closely hugging initial stage in places
        # This will create alphas in the range [0, 0.478260, 1]
        domain.set_quantity('stage', [[3,0,3], [2,2,6], [5,3,8], [8,3,5]])
        domain.set_quantity('elevation', [[0,0,0],
                                          [1.8,1.9,5.9],
                                          [4.6,0,0],
                                          [0,2,4]])
        stage = domain.quantities['stage']

        ref_centroid_values = copy.copy(stage.centroid_values[:])    # Copy
        ref_vertex_values = copy.copy(stage.vertex_values[:])        # Copy

        # Limit
        domain.tight_slope_limiters = 1
        domain.distribute_to_vertices_and_edges()

        # Assert that all vertex quantities have changed
        for k in range(len(domain)):
            assert not num.allclose(ref_vertex_values[k,:],
                                    stage.vertex_values[k,:])
        # and assert that quantities are still conserved
        for k in range(len(domain)):
            assert num.allclose(ref_centroid_values[k],
                                num.sum(stage.vertex_values[k,:])/3)

    def test_balance_deep_and_shallow_Froude(self):
        """Test that balanced limiters preserve conserved quantites -
        and also that excessive Froude numbers are dealt with.
        This test is using tight slope limiters.
        """

        import copy

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #             bac,     bce,     ecf,     dbe
        elements = [[1,0,2], [1,2,4], [4,2,5], [3,1,4] ]

        domain = Domain(points, elements)
        domain.check_integrity()
        domain.tight_slope_limiters = True
        domain.use_centroid_velocities = True

        # Create non-flat bed - closely hugging initial stage in places
        # This will create alphas in the range [0, 0.478260, 1]
        domain.set_quantity('stage', [[3,0,3], [2,2,6], [5,3,8], [8,3,5]])
        domain.set_quantity('elevation', [[0,0,0],
                                          [1.8,1.999,5.999],
                                          [4.6,0,0],
                                          [0,2,4]])

        # Create small momenta, that nonetheless will generate large speeds
        # due to shallow depth at isolated vertices
        domain.set_quantity('xmomentum', -0.0058)
        domain.set_quantity('ymomentum', 0.0890)

        stage = domain.quantities['stage']
        elevation = domain.quantities['elevation']
        xmomentum = domain.quantities['xmomentum']
        ymomentum = domain.quantities['ymomentum']

        # Setup triangle #1 to mimick real Froude explosion observed
        # in the Onslow example 13 Nov 2007.
        stage.vertex_values[1,:] = [1.6385, 1.6361, 1.2953]
        elevation.vertex_values[1,:] = [1.6375, 1.6336, 0.4647]
        xmomentum.vertex_values[1,:] = [-0.0058, -0.0050, -0.0066]
        ymomentum.vertex_values[1,:] = [0.0890, 0.0890, 0.0890]

        xmomentum.interpolate()
        ymomentum.interpolate()
        stage.interpolate()
        elevation.interpolate()

        # Verify interpolation
        assert num.allclose(stage.centroid_values[1], 1.5233)
        assert num.allclose(elevation.centroid_values[1], 1.2452667)
        assert num.allclose(xmomentum.centroid_values[1], -0.0058)
        assert num.allclose(ymomentum.centroid_values[1], 0.089)

        # Derived quantities
        depth = stage-elevation
        u = xmomentum/depth
        v = ymomentum/depth

        denom = (depth*g)**0.5
        Fx = u/denom
        Fy = v/denom

        # Verify against Onslow example (14 Nov 2007)
        assert num.allclose(depth.centroid_values[1], 0.278033)
        assert num.allclose(u.centroid_values[1], -0.0208608)
        assert num.allclose(v.centroid_values[1], 0.3201055)

        assert num.allclose(denom.centroid_values[1],
                            num.sqrt(depth.centroid_values[1]*g))

        assert num.allclose(u.centroid_values[1]/denom.centroid_values[1],
                            -0.012637746977)
        assert num.allclose(Fx.centroid_values[1],
                            u.centroid_values[1]/denom.centroid_values[1])

        # Check that Froude numbers are small at centroids.
        assert num.allclose(Fx.centroid_values[1], -0.012637746977)
        assert num.allclose(Fy.centroid_values[1], 0.193924048435)

        # But Froude numbers are huge at some vertices and edges
        assert num.allclose(Fx.vertex_values[1,:], [-5.85888475e+01,
                                                    -1.27775313e+01,
                                                    -2.78511420e-03])

        assert num.allclose(Fx.edge_values[1,:], [-6.89150773e-03,
                                                  -7.38672488e-03,
                                                  -2.35626238e+01])

        assert num.allclose(Fy.vertex_values[1,:], [8.99035764e+02,
                                                    2.27440057e+02,
                                                    3.75568430e-02])

        assert num.allclose(Fy.edge_values[1,:], [1.05748998e-01,
                                                  1.06035244e-01,
                                                  3.88346947e+02])


        # The task is now to arrange the limiters such that Froude numbers
        # remain under control whil at the same time obeying the conservation
        # laws.
        ref_centroid_values = copy.copy(stage.centroid_values[:])    # Copy
        ref_vertex_values = copy.copy(stage.vertex_values[:])        # Copy

        # Limit (and invoke balance_deep_and_shallow)
        domain.tight_slope_limiters = 1
        domain.distribute_to_vertices_and_edges()

        # Redo derived quantities
        depth = stage - elevation
        u = xmomentum/depth
        v = ymomentum/depth

        # Assert that all vertex velocities stay within one
        # order of magnitude of centroid velocities.
        assert num.alltrue(num.absolute(u.vertex_values[1,:]) <=
                           num.absolute(u.centroid_values[1])*10)
        assert num.alltrue(num.absolute(v.vertex_values[1,:]) <=
                           num.absolute(v.centroid_values[1])*10)

        denom = (depth*g)**0.5
        Fx = u/denom
        Fy = v/denom

        # Assert that Froude numbers are less than max value (TBA)
        # at vertices, edges and centroids.
        from anuga.config import maximum_froude_number

        assert num.alltrue(num.absolute(Fx.vertex_values[1,:]) <
                           maximum_froude_number)
        assert num.alltrue(num.absolute(Fy.vertex_values[1,:]) <
                           maximum_froude_number)

        # Assert that all vertex quantities have changed
        for k in range(len(domain)):
            assert not num.allclose(ref_vertex_values[k,:],
                                    stage.vertex_values[k,:])

        # Assert that quantities are still conserved
        for k in range(len(domain)):
            assert num.allclose(ref_centroid_values[k],
                                num.sum(stage.vertex_values[k,:])/3)

        return

        qwidth = 12
        for k in [1]:    # range(len(domain)):
            print 'Triangle %d (C, V, E)' % k

            print ('stage'.ljust(qwidth), stage.centroid_values[k],
                   stage.vertex_values[k,:], stage.edge_values[k,:])
            print ('elevation'.ljust(qwidth), elevation.centroid_values[k],
                   elevation.vertex_values[k,:], elevation.edge_values[k,:])
            print ('depth'.ljust(qwidth), depth.centroid_values[k],
                   depth.vertex_values[k,:], depth.edge_values[k,:])
            print ('xmomentum'.ljust(qwidth), xmomentum.centroid_values[k],
                   xmomentum.vertex_values[k,:], xmomentum.edge_values[k,:])
            print ('ymomentum'.ljust(qwidth), ymomentum.centroid_values[k],
                   ymomentum.vertex_values[k,:], ymomentum.edge_values[k,:])
            print ('u'.ljust(qwidth),u.centroid_values[k],
                   u.vertex_values[k,:], u.edge_values[k,:])
            print ('v'.ljust(qwidth), v.centroid_values[k],
                   v.vertex_values[k,:], v.edge_values[k,:])
            print ('Fx'.ljust(qwidth), Fx.centroid_values[k],
                   Fx.vertex_values[k,:], Fx.edge_values[k,:])
            print ('Fy'.ljust(qwidth), Fy.centroid_values[k],
                   Fy.vertex_values[k,:], Fy.edge_values[k,:])


#################################################################################

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_swb_clean, 'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)

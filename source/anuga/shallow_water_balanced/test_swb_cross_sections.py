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

    def test_get_flow_through_cross_section_with_geo(self):
        """test_get_flow_through_cross_section(self):

        Test that the total flow through a cross section can be
        correctly obtained at run-time from the ANUGA domain.

        This test creates a flat bed with a known flow through it and tests
        that the function correctly returns the expected flow.

        The specifics are
        e = -1 m
        u = 2 m/s
        h = 2 m
        w = 3 m (width of channel)

        q = u*h*w = 12 m^3/s

        This run tries it with georeferencing and with elevation = -1
        """

        import time, os
        from Scientific.IO.NetCDF import NetCDFFile
        from mesh_factory import rectangular

        # Create basic mesh (20m x 3m)
        width = 3
        length = 20
        t_end = 1
        points, vertices, boundary = rectangular(length, width, length, width)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary,
                        geo_reference=Geo_reference(56, 308500, 6189000))

        domain.default_order = 2
        domain.set_quantities_to_be_stored(None)

        e = -1.0
        w = 1.0
        h = w-e
        u = 2.0
        uh = u*h

        Br = anuga.Reflective_boundary(domain)     # Side walls
        Bd = anuga.Dirichlet_boundary([w, uh, 0])  # 2 m/s across the 3 m inlet: 


        # Initial conditions
        domain.set_quantity('elevation', e)
        domain.set_quantity('stage', w)
        domain.set_quantity('xmomentum', uh)
        domain.set_boundary({'left': Bd, 'right': Bd, 'top': Br, 'bottom': Br})

        # Interpolation points down the middle
        I = [[0, width/2.],
             [length/2., width/2.],
             [length, width/2.]]
        interpolation_points = domain.geo_reference.get_absolute(I)

        # Shortcuts to quantites
        stage = domain.get_quantity('stage')
        xmomentum = domain.get_quantity('xmomentum')
        ymomentum = domain.get_quantity('ymomentum')

        for t in domain.evolve(yieldstep=0.1, finaltime=t_end):
            # Check that quantities are they should be in the interior
            w_t = stage.get_values(interpolation_points)
            uh_t = xmomentum.get_values(interpolation_points)
            vh_t = ymomentum.get_values(interpolation_points)

            assert num.allclose(w_t, w)
            assert num.allclose(uh_t, uh)
            assert num.allclose(vh_t, 0.0, atol=1.0e-6)

            # Check flows through the middle
            for i in range(5):
                x = length/2. + i*0.23674563    # Arbitrary
                cross_section = [[x, 0], [x, width]]

                cross_section = domain.geo_reference.get_absolute(cross_section)
                Q = domain.get_flow_through_cross_section(cross_section,
                                                          verbose=False)

                assert num.allclose(Q, uh*width)

    def test_get_energy_through_cross_section_with_geo(self):
        """test_get_energy_through_cross_section(self):

        Test that the total and specific energy through a cross section can be
        correctly obtained at run-time from the ANUGA domain.

        This test creates a flat bed with a known flow through it and tests
        that the function correctly returns the expected energies.

        The specifics are
        e = -1 m
        u = 2 m/s
        h = 2 m
        w = 3 m (width of channel)

        q = u*h*w = 12 m^3/s

        This run tries it with georeferencing and with elevation = -1
        """

        import time, os
        from Scientific.IO.NetCDF import NetCDFFile
        from mesh_factory import rectangular

        # Create basic mesh (20m x 3m)
        width = 3
        length = 20
        t_end = 1
        points, vertices, boundary = rectangular(length, width, length, width)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary,
                        geo_reference=Geo_reference(56, 308500, 6189000))

        domain.default_order = 2
        domain.set_quantities_to_be_stored(None)

        e = -1.0
        w = 1.0
        h = w-e
        u = 2.0
        uh = u*h

        Br = anuga.Reflective_boundary(domain)     # Side walls
        Bd = anuga.Dirichlet_boundary([w, uh, 0])  # 2 m/s across the 3 m inlet:

        # Initial conditions
        domain.set_quantity('elevation', e)
        domain.set_quantity('stage', w)
        domain.set_quantity('xmomentum', uh)
        domain.set_boundary({'left': Bd, 'right': Bd, 'top': Br, 'bottom': Br})

        # Interpolation points down the middle
        I = [[0, width/2.],
             [length/2., width/2.],
             [length, width/2.]]
        interpolation_points = domain.geo_reference.get_absolute(I)

        # Shortcuts to quantites
        stage = domain.get_quantity('stage')
        xmomentum = domain.get_quantity('xmomentum')
        ymomentum = domain.get_quantity('ymomentum')

        for t in domain.evolve(yieldstep=0.1, finaltime=t_end):
            # Check that quantities are they should be in the interior
            w_t = stage.get_values(interpolation_points)
            uh_t = xmomentum.get_values(interpolation_points)
            vh_t = ymomentum.get_values(interpolation_points)

            assert num.allclose(w_t, w)
            assert num.allclose(uh_t, uh)
            assert num.allclose(vh_t, 0.0, atol=1.0e-6)

            # Check energies through the middle
            for i in range(5):
                x = length/2. + i*0.23674563    # Arbitrary
                cross_section = [[x, 0], [x, width]]

                cross_section = domain.geo_reference.get_absolute(cross_section)
                Es = domain.get_energy_through_cross_section(cross_section,
                                                             kind='specific',
                                                             verbose=False)

                assert num.allclose(Es, h + 0.5*u*u/g)

                Et = domain.get_energy_through_cross_section(cross_section,
                                                             kind='total',
                                                             verbose=False)
                assert num.allclose(Et, w + 0.5*u*u/g)


    def test_cross_section_class(self):
        """test_cross_section_class(self):

        Test that the total and specific energy through a cross section can be
        correctly obtained at run-time from the ANUGA cross section class.

        This test creates a flat bed with a known flow through it, creates a cross
        section and tests that the correct flow and energies are calculated

        The specifics are
        e = -1 m
        u = 2 m/s
        h = 2 m
        w = 3 m (width of channel)

        q = u*h*w = 12 m^3/s

        This run tries it with georeferencing and with elevation = -1
        """

        import time, os
        from Scientific.IO.NetCDF import NetCDFFile
        from mesh_factory import rectangular

        # Create basic mesh (20m x 3m)
        width = 3
        length = 20
        t_end = 1
        points, vertices, boundary = rectangular(length, width, length, width)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary,
                        geo_reference=Geo_reference(56, 308500, 6189000))

        domain.default_order = 2
        domain.set_quantities_to_be_stored(None)

        e = -1.0
        w = 1.0
        h = w-e
        u = 2.0
        uh = u*h

        Br = anuga.Reflective_boundary(domain)     # Side walls
        Bd = anuga.Dirichlet_boundary([w, uh, 0])  # 2 m/s across the 3 m inlet:

        # Initial conditions
        domain.set_quantity('elevation', e)
        domain.set_quantity('stage', w)
        domain.set_quantity('xmomentum', uh)
        domain.set_boundary({'left': Bd, 'right': Bd, 'top': Br, 'bottom': Br})

        # Interpolation points down the middle
        I = [[0, width/2.],
             [length/2., width/2.],
             [length, width/2.]]
        interpolation_points = domain.geo_reference.get_absolute(I)

        # Shortcuts to quantites
        stage = domain.get_quantity('stage')
        xmomentum = domain.get_quantity('xmomentum')
        ymomentum = domain.get_quantity('ymomentum')


        # Create some cross sections
        cross_sections = []
        for i in range(5):
            x = length/2. + i*0.23674563    # Arbitrary
            polyline = [[x, 0], [x, width]]

            polyline = domain.geo_reference.get_absolute(polyline)
        
            cross_sections.append(Cross_section(domain,polyline))


 
        for t in domain.evolve(yieldstep=0.1, finaltime=t_end):
            # Check that quantities are they should be in the interior
            w_t = stage.get_values(interpolation_points)
            uh_t = xmomentum.get_values(interpolation_points)
            vh_t = ymomentum.get_values(interpolation_points)

            assert num.allclose(w_t, w)
            assert num.allclose(uh_t, uh)
            assert num.allclose(vh_t, 0.0, atol=1.0e-6)


            # Check flows and energies through the middle
            for cross_section in cross_sections:
                
                Q = cross_section.get_flow_through_cross_section()

                assert num.allclose(Q, uh*width)

                Es = cross_section.get_energy_through_cross_section(kind='specific')

                assert num.allclose(Es, h + 0.5*u*u/g)

                Et = cross_section.get_energy_through_cross_section(kind='total')
                
                assert num.allclose(Et, w + 0.5*u*u/g)



#################################################################################

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_swb_clean, 'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)

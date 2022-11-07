#!/usr/bin/env python


import unittest, os, time
import os.path
from math import pi, sqrt
import tempfile

from anuga.file.netcdf import NetCDFFile
from anuga.file.sww import extent_sww

from anuga.config import g, epsilon
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a
from anuga.utilities.numerical_tools import mean, ensure_numeric
from anuga.geometry.polygon import is_inside_polygon
from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.geospatial_data.geospatial_data import Geospatial_data
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross, \
                                            rectangular
from anuga.abstract_2d_finite_volumes.quantity import Quantity
from anuga.shallow_water.forcing import Inflow, Cross_section
from anuga.geospatial_data.geospatial_data import ensure_geospatial

from anuga.utilities.system_tools import get_pathname_from_package

from anuga.abstract_2d_finite_volumes.generic_boundary_conditions \
        import Dirichlet_boundary
from anuga.shallow_water.forcing import Rainfall, Wind_stress
from anuga.shallow_water.forcing import Inflow, Cross_section
from anuga.shallow_water.sww_interrogate import get_flow_through_cross_section

from anuga.shallow_water.shallow_water_domain import Domain

# boundary functions
from anuga.shallow_water.boundaries import Reflective_boundary, \
            Field_boundary, Transmissive_momentum_set_stage_boundary, \
            Transmissive_stage_zero_momentum_boundary
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions\
     import Transmissive_boundary, Dirichlet_boundary, \
            Time_boundary, File_boundary, AWI_boundary

import numpy as num
from anuga.config import g

# Get gateway to C implementation of flux function for direct testing
from anuga.shallow_water.shallow_water_ext import flux_function_central as flux_function
from anuga.shallow_water.shallow_water_ext import rotate


def set_bottom_friction(tag, elements, domain):
    if tag == "bottom":
        domain.set_quantity('friction', 0.09, indices = elements)

def set_top_friction(tag, elements, domain):
    if tag == "top":
        domain.set_quantity('friction', 1., indices = elements)


def set_all_friction(tag, elements, domain):
    if tag == 'all':
        new_values = domain.get_quantity('friction').get_values(indices = elements) + 10.0

        domain.set_quantity('friction', new_values, indices = elements)


# For test_fitting_using_shallow_water_domain example
def linear_function(point):
    point = num.array(point)
    return point[:,0]+point[:,1]



def scalar_func(t, x, y):
    """Function that returns a scalar.

    Used to test error message when numeric array is expected
    """

    return 17.7



class Test_LoadSave(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        for file in ['domain.sww', 'domain_pickle.pickle']:
            try:
                os.remove(file)
            except:
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

        Br = Reflective_boundary(domain)     # Side walls
        Bd = Dirichlet_boundary([w, uh, 0])  # 2 m/s across the 3 m inlet: 


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

        for t in domain.evolve(yieldstep=0.1, finaltime=0.5):
            # Shortcuts to quantites
            stage = domain.get_quantity('stage')
            xmomentum = domain.get_quantity('xmomentum')
            ymomentum = domain.get_quantity('ymomentum')

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

        import pickle
        fid = open('domain_pickle.pickle', 'wb')
        pickle.dump(domain, fid)
        
        fid = open('domain_pickle.pickle', 'rb')
        domain_restored = pickle.load(fid)

        
        for t in domain_restored.evolve(yieldstep=0.1, finaltime=1.0):
            # Shortcuts to quantites
            stage = domain_restored.get_quantity('stage')
            xmomentum = domain_restored.get_quantity('xmomentum')
            ymomentum = domain_restored.get_quantity('ymomentum')
            
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

                cross_section = domain_restored.geo_reference.get_absolute(cross_section)
                Q = domain_restored.get_flow_through_cross_section(cross_section,
                                                          verbose=False)

                assert num.allclose(Q, uh*width)
                


#-------------------------------------------------------------

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_LoadSave, 'test')
    runner = unittest.TextTestRunner() #verbosity=2)
    runner.run(suite)    

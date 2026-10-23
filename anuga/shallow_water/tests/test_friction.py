"""  Test environmental forcing - rain, wind, etc.
"""


import unittest
import os
import anuga
from anuga import Domain
from anuga import Reflective_boundary

from anuga.shallow_water.friction import *

import numpy as np
import warnings



class Test_Friction(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        for file in ['domain.sww']:
            try:
                os.remove(file)
            except OSError:
                pass



    def test_manning_friction_flat_implicit(self):
        """Test the manning friction implicit forcing term
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

        domain = Domain(points, vertices)
        # Mode-2 ('unified') computes friction on-device and never syncs the host
        # semi_implicit_update array this white-box test reads. Pin legacy to
        # exercise the in-Python forcing machinery it is written for.
        domain.set_compute_mode('legacy')

        #Flat surface with 1m of water
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 1.00)
        domain.set_quantity('xmomentum', 1.00)
        domain.set_quantity('ymomentum', 2.00)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        #Test friction forcing term
        domain.compute_forcing_terms()


        #import pprint
        #pprint.pprint(domain.quantities['xmomentum'].explicit_update)
        #pprint.pprint(domain.quantities['ymomentum'].explicit_update)

        #pprint.pprint(domain.quantities['xmomentum'].semi_implicit_update)
        #pprint.pprint(domain.quantities['ymomentum'].semi_implicit_update)

        xmon_semi_implicit_update = np.array([-21.91346618, -21.91346618, -21.91346618, -21.91346618])
        ymon_semi_implicit_update = np.array([-43.82693236, -43.82693236, -43.82693236, -43.82693236])

        assert num.allclose(domain.quantities['stage'].explicit_update, 0.0)
        assert num.allclose(domain.quantities['xmomentum'].explicit_update, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].explicit_update, 0.0)

        assert num.allclose(domain.quantities['xmomentum'].semi_implicit_update, xmon_semi_implicit_update)
        assert num.allclose(domain.quantities['ymomentum'].semi_implicit_update, ymon_semi_implicit_update)


    def test_manning_friction_sloped_implicit(self):
        """Test the manning friction sloped implicit forcing term
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

        domain = Domain(points, vertices)
        # Mode-2 ('unified') computes friction on-device and never syncs the host
        # semi_implicit_update array this white-box test reads. Pin legacy to
        # exercise the in-Python forcing machinery it is written for.
        domain.set_compute_mode('legacy')

        def slope(x, y):
            """Define a slope for the surface"""
            return 0.5 * x + 0.5 * y

        def stage(x, y):
            """Define a stage that is 1m above the slope"""
            return slope(x, y) + 1.0

        #Flat surface with 1m of water
        domain.set_quantity('elevation', slope)
        domain.set_quantity('stage', stage)
        domain.set_quantity('friction', 1.00)
        domain.set_quantity('xmomentum', 1.00)
        domain.set_quantity('ymomentum', 2.00)

        domain.use_sloped_mannings = True

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        #Test friction forcing term
        domain.compute_forcing_terms()


        #import pprint
        #pprint.pprint(domain.quantities['xmomentum'].explicit_update)
        #pprint.pprint(domain.quantities['ymomentum'].explicit_update)

        #pprint.pprint(domain.quantities['xmomentum'].semi_implicit_update)
        #pprint.pprint(domain.quantities['ymomentum'].semi_implicit_update)

        #xmon_semi_implicit_update = np.array([-21.91346618, -21.91346618, -21.91346618, -21.91346618])
        #ymon_semi_implicit_update = np.array([-43.82693236, -43.82693236, -43.82693236, -43.82693236])

        xmon_semi_implicit_update = np.array([-26.83840532, -26.83840532, -26.83840532, -26.83840532])
        ymon_semi_implicit_update = np.array([-53.67681064, -53.67681064, -53.67681064, -53.67681064])

        assert num.allclose(domain.quantities['stage'].explicit_update, 0.0)
        assert num.allclose(domain.quantities['xmomentum'].explicit_update, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].explicit_update, 0.0)

        assert num.allclose(domain.quantities['xmomentum'].semi_implicit_update, xmon_semi_implicit_update)
        assert num.allclose(domain.quantities['ymomentum'].semi_implicit_update, ymon_semi_implicit_update)



if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_Friction)
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)


# ---------------------------------------------------------------------------
# Plain pytest functions for branches not covered by the TestCase above
# ---------------------------------------------------------------------------

import pytest
from anuga.config import MULTIPROCESSOR_OPENMP, MULTIPROCESSOR_GPU
from anuga.abstract_2d_finite_volumes.quantity import Quantity


def _simple_domain(stage=1.0):
    points = [[0, 0], [0, 2], [2, 0], [0, 4], [2, 2], [4, 0]]
    vertices = [[1, 0, 2], [1, 2, 4], [4, 2, 5], [3, 1, 4]]
    domain = Domain(points, vertices)
    domain.set_quantity('elevation', 0.0)
    domain.set_quantity('stage', stage)
    domain.set_boundary({'exterior': Reflective_boundary(domain)})
    return domain


def test_linear_friction_updates_semi_implicit():
    """linear_friction: tau/h applied correctly to semi_implicit updates."""
    domain = _simple_domain(stage=1.0)   # depth h = 1.0
    domain.set_quantity('xmomentum', 2.0)
    domain.set_quantity('ymomentum', 1.0)
    domain.quantities['linear_friction'] = Quantity(domain)
    domain.set_quantity('linear_friction', 0.1)  # tau = 0.1

    from anuga.shallow_water.friction import linear_friction
    linear_friction(domain)

    # S = -tau/h = -0.1 → xmom_update += -0.1*2 = -0.2, ymom_update += -0.1*1 = -0.1
    np.testing.assert_allclose(
        domain.quantities['xmomentum'].semi_implicit_update, -0.2, rtol=1e-10)
    np.testing.assert_allclose(
        domain.quantities['ymomentum'].semi_implicit_update, -0.1, rtol=1e-10)


def test_linear_friction_dry_no_update():
    """linear_friction skips dry elements (h < eps)."""
    domain = _simple_domain(stage=0.0)   # depth = 0
    domain.set_quantity('xmomentum', 2.0)
    domain.quantities['linear_friction'] = Quantity(domain)
    domain.set_quantity('linear_friction', 0.1)

    from anuga.shallow_water.friction import linear_friction
    linear_friction(domain)

    np.testing.assert_allclose(
        domain.quantities['xmomentum'].semi_implicit_update, 0.0, atol=1e-14)


def test_depth_dependent_friction_below_d1():
    """depth < d1 → returns n1 for all wet elements."""
    domain = _simple_domain(stage=0.3)   # depth 0.3 < d1=0.5
    N = len(domain)
    data = np.zeros((N, 5))
    data[:, 1] = 0.5   # d1
    data[:, 2] = 0.03  # n1
    data[:, 3] = 2.0   # d2
    data[:, 4] = 0.06  # n2

    from anuga.shallow_water.friction import depth_dependent_friction
    result = depth_dependent_friction(domain, 0.03, data)
    wet = domain.get_wet_elements()
    np.testing.assert_allclose(result[wet], 0.03, rtol=1e-10)


def test_depth_dependent_friction_above_d2():
    """depth > d2 → returns n2 for all wet elements."""
    domain = _simple_domain(stage=5.0)   # depth 5 > d2=3
    N = len(domain)
    data = np.zeros((N, 5))
    data[:, 1] = 0.5
    data[:, 2] = 0.03
    data[:, 3] = 3.0   # d2
    data[:, 4] = 0.06  # n2

    from anuga.shallow_water.friction import depth_dependent_friction
    result = depth_dependent_friction(domain, 0.03, data)
    wet = domain.get_wet_elements()
    np.testing.assert_allclose(result[wet], 0.06, rtol=1e-10)


def test_depth_dependent_friction_interpolated():
    """d1 < depth < d2 → linearly interpolated friction."""
    domain = _simple_domain(stage=1.0)   # depth 1.0 between d1=0.5 and d2=2.0
    N = len(domain)
    data = np.zeros((N, 5))
    data[:, 1] = 0.5
    data[:, 2] = 0.03
    data[:, 3] = 2.0
    data[:, 4] = 0.06

    from anuga.shallow_water.friction import depth_dependent_friction
    result = depth_dependent_friction(domain, 0.03, data)
    wet = domain.get_wet_elements()
    expected = 0.03 + (0.06 - 0.03) / (2.0 - 0.5) * (1.0 - 0.5)  # 0.04
    np.testing.assert_allclose(result[wet], expected, rtol=1e-10)


def test_manning_friction_gpu_mode_matches_openmp():
    """GPU mode falls back to the same OpenMP extension and produces identical results."""
    domain1 = _simple_domain(stage=1.0)
    domain1.set_quantity('friction', 0.03)
    domain1.set_quantity('xmomentum', 1.0)
    domain1.multiprocessor_mode = MULTIPROCESSOR_OPENMP
    domain1.use_sloped_mannings = False

    domain2 = _simple_domain(stage=1.0)
    domain2.set_quantity('friction', 0.03)
    domain2.set_quantity('xmomentum', 1.0)
    domain2.multiprocessor_mode = MULTIPROCESSOR_GPU
    domain2.use_sloped_mannings = False

    from anuga.shallow_water.friction import manning_friction_semi_implicit
    manning_friction_semi_implicit(domain1)
    manning_friction_semi_implicit(domain2)

    np.testing.assert_allclose(
        domain1.quantities['xmomentum'].semi_implicit_update,
        domain2.quantities['xmomentum'].semi_implicit_update)


def test_manning_friction_gpu_mode_sloped():
    """GPU mode sloped branch runs without error."""
    domain = _simple_domain(stage=1.0)
    domain.set_quantity('friction', 0.03)
    domain.set_quantity('xmomentum', 1.0)
    domain.multiprocessor_mode = MULTIPROCESSOR_GPU
    domain.use_sloped_mannings = True

    from anuga.shallow_water.friction import manning_friction_semi_implicit
    manning_friction_semi_implicit(domain)   # should not raise


def test_manning_friction_invalid_mode_raises():
    """Unsupported multiprocessor_mode raises ValueError."""
    domain = _simple_domain(stage=1.0)
    domain.multiprocessor_mode = 99

    from anuga.shallow_water.friction import manning_friction_semi_implicit
    with pytest.raises(ValueError):
        manning_friction_semi_implicit(domain)

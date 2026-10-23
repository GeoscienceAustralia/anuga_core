"""
Unit tests for anuga.scenario.setup_sediment: the parsed [sediment] table
reaches the domain's sediment API.
"""
import unittest
from types import SimpleNamespace

import numpy as np

from anuga import rectangular_cross_domain, Reflective_boundary
from anuga.scenario.setup_sediment import setup_sediment


def _make_domain():
    d = rectangular_cross_domain(4, 4, len1=10.0, len2=10.0)
    d.set_flow_algorithm('DE0')
    d.store = False
    d.set_quantity('elevation', 0.0)
    d.set_quantity('stage', 1.0)
    d.set_quantity('friction', 0.03)
    return d


class TestSetupSediment(unittest.TestCase):

    def test_no_table_is_a_no_op(self):
        d = _make_domain()
        self.assertIsNone(setup_sediment(d, SimpleNamespace()))
        self.assertIsNone(setup_sediment(d, SimpleNamespace(sediment_data=None)))
        self.assertEqual(d.n_sediment_classes, 0)

    def test_minimal_adds_the_fraction_with_domain_defaults(self):
        d = _make_domain()
        project = SimpleNamespace(sediment_data={
            'fractions': [{'name': 'sand', 'diameter': 2.0e-4}],
            'erodible_regions': []})
        op = setup_sediment(d, project)
        self.assertIsNotNone(op)
        self.assertEqual(d.get_sediment_names(), ['sand'])
        self.assertAlmostEqual(d.sediment_porosity, 0.30)
        self.assertEqual(d.sediment_shear_closure, 0)
        self.assertEqual(d.sediment_erosion_mode, 0)
        self.assertTrue(d.sediment_bed_evolution)

    def test_every_setting_reaches_the_domain(self):
        d = _make_domain()
        project = SimpleNamespace(sediment_data={
            'porosity': 0.28, 'c_max': 0.25, 'bed_evolution': False, 'rho_w': 1025.0,
            'shear_closure': 'depth_slope', 'max_slope': 0.02, 'freeze_slope': True,
            'bed_material': 'cohesive', 'tau_crit': 0.088,
            'deposition_law': 'threshold', 'tau_d': 0.1, 'near_bed': 'rouse',
            'friction_mode': 'wilson', 'bed': 'gravel', 'grain_size': 0.02,
            'bedload': 'wong_parker_eq24', 'bedload_K': 3.0,
            'bedload_open_boundaries': ['left', 'right'],
            'angle_of_repose': 34.0,
            'erodible_base_depth': 1.5,
            'fractions': [
                {'name': 'sand', 'diameter': 2.0e-4, 'tau_c_star': 0.045,
                 'initial_concentration': 0.001, 'boundary': {'left': 0.01}},
                {'name': 'silt', 'diameter': 2.0e-5}],
            'erodible_regions': [
                {'polygon_points': [[0, 0], [5, 0], [5, 10], [0, 10]], 'erodible': True},
                {'center': [7.5, 5.0], 'radius': 2.0, 'erodible': False}],
        })
        setup_sediment(d, project)
        self.assertEqual(d.get_sediment_names(), ['sand', 'silt'])
        self.assertAlmostEqual(d.sediment_porosity, 0.28)
        self.assertAlmostEqual(d.sediment_c_max, 0.25)
        self.assertFalse(d.sediment_bed_evolution)
        self.assertAlmostEqual(d.sediment_rho_w, 1025.0)
        self.assertEqual(d.sediment_shear_closure, 1)
        self.assertAlmostEqual(d.sediment_max_slope, 0.02)
        self.assertEqual(d.sediment_slope_frozen, 1)
        self.assertEqual(d.sediment_erosion_mode, 1)
        self.assertAlmostEqual(d.sediment_tau_crit, 0.088)
        self.assertEqual(d.sediment_deposition_mode, 1)
        self.assertAlmostEqual(d.sediment_tau_d, 0.1)
        self.assertEqual(d.sediment_d_star_mode, 1)
        self.assertEqual(d.sediment_friction_mode, 2)
        self.assertEqual(d.sediment_wilson_bed, 1)
        self.assertAlmostEqual(d.sediment_bedload_K, 3.0)
        self.assertEqual(d._sediment_bedload_open_tags, ('left', 'right'))
        self.assertEqual(int(d.sediment_bedload_open.sum()),
                         len(d.tag_boundary_cells['left'])
                         + len(d.tag_boundary_cells['right']))
        self.assertGreater(d.sediment_repose_tan, 0.0)
        self.assertTrue(d.sediment_has_z_base)
        # the base is 1.5 m below the bed, except where the locked circle
        # raised it to the bed itself
        xc, yc = d.centroid_coordinates.T
        locked = (xc - 7.5) ** 2 + (yc - 5.0) ** 2 <= 2.0 ** 2
        self.assertTrue(locked.any())
        self.assertTrue(np.allclose(d.sediment_z_base[~locked], -1.5))
        self.assertTrue(np.allclose(d.sediment_z_base[locked], 0.0))
        self.assertAlmostEqual(d.sediment_tau_c_star[0], 0.045)
        c = d.get_tracer('sand')
        self.assertTrue(np.allclose(c, 0.001))
        # the inflow concentration on 'left' and none on the others
        left = d.tag_boundary_cells['left']
        self.assertTrue(np.allclose(d.tracer_boundary_values[0][left], 0.01))
        right = d.tag_boundary_cells['right']
        self.assertTrue(np.allclose(d.tracer_boundary_values[0][right], 0.0))

    def test_de_leeuw_entrainment_reaches_the_domain(self):
        d = _make_domain()
        project = SimpleNamespace(sediment_data={
            'entrainment': 'de_leeuw', 'de_leeuw_fit': 'nghiem_2022',
            'skin_roughness': 0.002,
            'fractions': [{'name': 'mud', 'diameter': 2.0e-5}],
            'erodible_regions': []})
        setup_sediment(d, project)
        self.assertEqual(d.sediment_erosion_mode, 3)
        self.assertAlmostEqual(d.sediment_dl_A, 7.04e-4)
        self.assertAlmostEqual(d.sediment_dl_ks, 0.002)

    def test_configured_domain_evolves(self):
        d = _make_domain()
        Br = Reflective_boundary(d)
        d.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})
        project = SimpleNamespace(sediment_data={
            'shear_closure': 'depth_slope',
            'fractions': [{'name': 'sand', 'diameter': 2.0e-4,
                           'initial_concentration': 0.01}],
            'erodible_regions': []})
        setup_sediment(d, project)
        for _ in d.evolve(yieldstep=0.5, finaltime=1.0):
            pass
        c = d.get_tracer('sand')
        self.assertTrue(np.all(np.isfinite(c)))
        self.assertLess(c.max(), 0.01)          # still water: it settles


if __name__ == '__main__':
    unittest.main()

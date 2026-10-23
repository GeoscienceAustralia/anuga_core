"""A ready-made anuga.Region is accepted wherever a polygon, a circle or
triangle indices are (anuga-community/anuga_core#15): Quantity.set_values
(and so Domain.set_quantity), the erosion operators and Set_w_uh_vh.
"""
import unittest

import numpy as num
import pytest

import anuga
from anuga import Region, Reflective_boundary, rectangular_cross_domain
from anuga.operators.erosion_operators import Erosion_operator
from anuga.operators.sanddune_erosion_operator import Sanddune_erosion_operator
from anuga.operators.set_w_uh_vh_operator import Set_w_uh_vh_operator

POLY = [[2.0, 2.0], [8.0, 2.0], [8.0, 6.0], [2.0, 6.0]]


def make_domain():
    d = rectangular_cross_domain(10, 10, len1=10.0, len2=10.0)
    d.set_flow_algorithm('DE0')
    d.store = False
    d.set_quantity('elevation', 0.5)
    d.set_quantity('stage', 1.0)
    d.set_quantity('friction', 0.0)
    d.set_boundary({t: Reflective_boundary(d) for t in d.get_boundary_tags()})
    return d


class TestSetValuesRegion(unittest.TestCase):

    def test_region_matches_polygon(self):
        d = make_domain()
        d.set_quantity('elevation', 0.0)
        d.set_quantity('elevation', 3.0, polygon=POLY)
        by_polygon = d.quantities['elevation'].centroid_values.copy()
        d.set_quantity('elevation', 0.0)
        d.set_quantity('elevation', 3.0, region=Region(d, polygon=POLY))
        by_region = d.quantities['elevation'].centroid_values
        assert num.allclose(by_region, by_polygon)
        assert (by_region == 3.0).any() and (by_region == 0.0).any()

    def test_circle_and_indices_regions(self):
        d = make_domain()
        circle = Region(d, center=[5.0, 5.0], radius=2.0)
        d.quantities['friction'].set_values(0.07, region=circle)
        f = d.quantities['friction'].centroid_values
        assert num.allclose(f[circle.indices], 0.07)
        outside = num.setdiff1d(num.arange(len(f)), circle.indices)
        assert num.allclose(f[outside], 0.0)
        picked = Region(d, indices=[0, 5, 9])
        d.quantities['friction'].set_values(0.02, region=picked)
        assert num.allclose(f[[0, 5, 9]], 0.02)

    def test_region_with_other_selectors_is_rejected(self):
        d = make_domain()
        r = Region(d, polygon=POLY)
        with pytest.raises(Exception, match='Only one of polygon, region and indices'):
            d.set_quantity('elevation', 1.0, region=r, polygon=POLY)
        with pytest.raises(Exception, match='Only one of polygon, region and indices'):
            d.set_quantity('elevation', 1.0, region=r, indices=[0])
        with pytest.raises(Exception, match='must \\(currently\\) be a constant'):
            d.set_quantity('elevation', lambda x, y: x, region=r)


class TestOperatorRegion(unittest.TestCase):

    def test_erosion_operators_take_a_region(self):
        for cls in (Erosion_operator, Sanddune_erosion_operator):
            d = make_domain()
            by_polygon = cls(d, polygon=POLY)
            by_region = cls(d, region=Region(d, polygon=POLY))
            assert sorted(by_region.indices) == sorted(by_polygon.indices), cls
            assert len(by_region.indices) > 0
            with pytest.raises(ValueError, match='cannot specify both'):
                cls(d, region=Region(d, polygon=POLY), polygon=POLY)

    def test_set_w_uh_vh_takes_a_region(self):
        d = make_domain()
        r = Region(d, center=[5.0, 5.0], radius=2.0)
        op = Set_w_uh_vh_operator(d, w_uh_vh=[3.0, 4.0, 5.0], region=r)
        assert sorted(op.indices) == sorted(r.indices)
        d.timestep = 1.0
        op()
        stage = d.quantities['stage'].centroid_values
        assert num.allclose(stage[r.indices], 3.0)
        outside = num.setdiff1d(num.arange(len(stage)), r.indices)
        assert num.allclose(stage[outside], 1.0)
        with pytest.raises(ValueError, match='cannot specify both'):
            Set_w_uh_vh_operator(d, w_uh_vh=[3.0, 4.0, 5.0], region=r, indices=[0])


if __name__ == '__main__':
    unittest.main()

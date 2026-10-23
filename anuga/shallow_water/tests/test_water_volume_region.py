"""get_water_volume(region=...) / (indices=...) -- issue #22.

A Quantity's get_integral takes a region; the domain's water volume did not,
so summing the water in part of a domain meant rebuilding the depth by hand.
"""
import numpy as np
import pytest

import anuga
from anuga import Domain, Reflective_boundary, Region


def make_domain(depth=0.5):
    points, vertices, boundary = anuga.rectangular_cross(10, 10, 10.0, 10.0)
    domain = Domain(points, vertices, boundary)
    domain.set_flow_algorithm('DE0')
    domain.store = False
    domain.set_quantity('elevation', 0.0)
    domain.set_quantity('friction', 0.0)
    domain.set_quantity('stage', depth)
    Br = Reflective_boundary(domain)
    domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})
    return domain


SQUARE = [[2.0, 2.0], [6.0, 2.0], [6.0, 5.0], [2.0, 5.0]]   # 12 m^2 on mesh lines


def test_region_volume_is_depth_times_region_area():
    domain = make_domain(depth=0.5)
    region = Region(domain, polygon=SQUARE)
    area = domain.get_areas()[region.get_indices()].sum()
    assert np.isclose(area, 12.0)
    assert np.isclose(domain.get_water_volume(region=region), 0.5 * area)
    assert np.isclose(domain.compute_total_volume(region=region), 0.5 * area)


def test_indices_match_region_and_partial_sums_add_up():
    domain = make_domain(depth=0.5)
    region = Region(domain, polygon=SQUARE)
    inside = region.get_indices()
    outside = np.setdiff1d(np.arange(domain.number_of_elements), inside)
    v_in = domain.get_water_volume(indices=inside)
    v_out = domain.get_water_volume(indices=outside)
    assert np.isclose(v_in, domain.get_water_volume(region=region))
    assert np.isclose(v_in + v_out, domain.get_water_volume())
    assert np.isclose(domain.get_water_volume(), 0.5 * 100.0)


def test_region_query_after_evolve_uses_the_evolved_depth():
    """After evolve the depth comes from the height / stage-elevation
    quantities rather than the initial centroid values."""
    domain = make_domain(depth=0.5)
    domain.set_quantity('stage', lambda x, y: 0.5 + 0.2 * (x < 5.0))
    region = Region(domain, polygon=SQUARE)
    for _ in domain.evolve(yieldstep=0.5, finaltime=1.0):
        pass
    h = (domain.quantities['stage'].centroid_values
         - domain.quantities['elevation'].centroid_values)
    idx = region.get_indices()
    expected = (h[idx] * domain.get_areas()[idx]).sum()
    assert np.isclose(domain.get_water_volume(region=region), expected)


def test_regional_queries_do_not_pollute_the_volume_history():
    domain = make_domain()
    region = Region(domain, polygon=SQUARE)
    n0 = len(domain.volume_history)
    domain.get_water_volume(region=region)
    domain.get_water_volume(indices=[0, 1, 2])
    assert len(domain.volume_history) == n0
    domain.get_water_volume()
    assert len(domain.volume_history) == n0 + 1


def test_region_and_indices_together_are_rejected():
    domain = make_domain()
    region = Region(domain, polygon=SQUARE)
    with pytest.raises(AssertionError):
        domain.get_water_volume(region=region, indices=[0])


def test_empty_region_gives_zero():
    domain = make_domain()
    region = Region(domain, polygon=[[20.0, 20.0], [21.0, 20.0], [21.0, 21.0]])
    assert domain.get_water_volume(region=region) == 0.0

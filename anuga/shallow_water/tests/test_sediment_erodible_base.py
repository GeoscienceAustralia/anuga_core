"""The non-erodible base [L-5] and the region restriction (PHYSICS_SPEC 4.5).

Erosion is otherwise bottomless: the bed lowers for as long as the flow can
lift material. Right for a deep alluvial bed, wrong wherever the erodible layer
is finite -- a reach floored by an outcrop, a lined culvert, a dam apron.

The tests that matter are not "does z stop at the base" alone; a clamp on z
would pass that while creating sediment from nothing. They are that the floor
holds per cell including where it varies in space, that the budget still
closes, that bedload stays conservative, that the classes are limited together
so registration order does not matter, and -- the regression gate -- that an
unreachable base reproduces the no-base answer BITWISE.

The mode 1 / mode 2 comparison is in test_sediment_gpu.py.
"""
import numpy as np
import pytest

from anuga import Reflective_boundary, rectangular_cross_domain
from anuga.abstract_2d_finite_volumes.region import Region

POROSITY = 0.3
TWO_CLASSES = (('fine', 1.0e-4), ('coarse', 6.0e-4))
HALF = [[0.0, -1.0], [30.0, -1.0], [30.0, 17.0], [0.0, 17.0]]   # x < 30


def build(base=None, base_depth=None, classes=(('sand', 2.0e-4),),
          bedload=False, slope=0.01, depth=0.6, n_x=30):
    d = rectangular_cross_domain(n_x, 8, len1=60.0, len2=16.0)
    d.set_flow_algorithm('DE0')
    d.set_low_froude(0)
    d.store = False
    d.set_quantity('elevation', lambda x, y: -slope * x)
    d.set_quantity('stage', lambda x, y: -slope * x + depth)
    d.set_quantity('xmomentum', 1.2)
    d.set_quantity('ymomentum', 0.0)
    d.set_quantity('friction', 0.03)
    d.set_boundary({t: Reflective_boundary(d) for t in d.get_boundary_tags()})
    d.set_sediment_parameters(porosity=POROSITY)
    if bedload:
        d.set_bedload('wong_parker_eq24')
    for name, diameter in classes:
        d.add_sediment_fraction(name=name, diameter=diameter)
    if base is not None:
        d.set_erodible_base(elevation=base)
    elif base_depth is not None:
        d.set_erodible_base(depth=base_depth)
    return d


def _mass(d):
    return sum(float((d.tracer_conserved_values[s] * d.areas).sum())
               for s in range(d.n_sediment_classes))


# ---------------------------------------------------------------- the floor

def test_no_cell_erodes_below_its_base():
    d = build(base_depth=0.02)
    d.evolve_to_end(finaltime=40.0)
    z = d.quantities['elevation'].centroid_values
    assert (d.sediment_z_base - z).max() <= 0.0


def test_the_base_is_actually_reached():
    """Without this, the test above passes on a run that never got near it."""
    d = build(base_depth=0.02)
    d.evolve_to_end(finaltime=40.0)
    z = d.quantities['elevation'].centroid_values
    assert int((z - d.sediment_z_base < 1e-12).sum()) > 0


def test_the_unlimited_run_scours_deeper_than_the_layer():
    """Makes the floor test non-vacuous: erosion really did want to go further."""
    d = build()
    z0 = d.quantities['elevation'].centroid_values.copy()
    d.evolve_to_end(finaltime=40.0)
    assert (z0 - d.quantities['elevation'].centroid_values).max() > 0.02


def test_a_base_that_varies_in_space_is_honoured_cell_by_cell():
    """Bedrock is a surface, not a level: a thin layer upstream of a step and a
    thick one below it must be respected separately."""
    d = build()
    zc = d.quantities['elevation'].centroid_values.copy()
    x = d.centroid_coordinates[:, 0]
    d.set_erodible_base(depth=np.where(x < 30.0, 0.01, 0.5))
    base = d.sediment_z_base.copy()
    d.evolve_to_end(finaltime=40.0)
    z = d.quantities['elevation'].centroid_values
    assert (base - z).max() <= 0.0
    thin, thick = x < 30.0, x >= 30.0
    assert (zc - z)[thin].max() <= 0.01 + 1e-12
    assert (zc - z)[thick].max() > 0.01


# ---------------------------------------------------------------- the budget

def test_the_budget_still_closes_when_the_limiter_engages():
    """The limit is on the SOURCE, so sediment not eroded never enters the
    water column. Clamping z instead would leave suspended sediment that came
    from nowhere."""
    d = build(base_depth=0.02)
    z0 = d.quantities['elevation'].centroid_values.copy()
    m0 = _mass(d)
    d.evolve_to_end(finaltime=40.0)
    z1 = d.quantities['elevation'].centroid_values
    m1 = _mass(d)
    bed = float(((1.0 - POROSITY) * (z1 - z0) * d.areas).sum())
    assert abs((m1 + bed) - m0) < 1e-12 * max(abs(m1), 1e-30)
    assert int((z1 - d.sediment_z_base < 1e-12).sum()) > 0, 'limiter never engaged'


def test_bedload_stays_conservative_with_a_base():
    d = build(base_depth=0.01, bedload=True, classes=(('sand', 8.0e-4),))
    z0 = d.quantities['elevation'].centroid_values.copy()
    m0 = _mass(d)
    d.evolve_to_end(finaltime=40.0)
    z1 = d.quantities['elevation'].centroid_values
    m1 = _mass(d)
    bed = float(((1.0 - POROSITY) * (z1 - z0) * d.areas).sum())
    assert abs((m1 + bed) - m0) < 1e-12 * max(abs(m1), 1e-30)


def test_bedload_holds_the_base_to_within_one_step_of_flux():
    """Bedload's floor is NOT exact, and the reason is structural.

    When a cell cannot pay for the step's divergence its removing edges are
    closed -- symmetrically, which is what keeps the scheme conservative -- but
    closing an edge also cancels the INFLOW its neighbour was to receive, so
    the neighbour's own removal grows and the deficit migrates one cell per
    sweep. Driving it to zero needs the exhaustion flag iterated to a fixed
    point with double buffering, so the answer stays independent of thread
    order; that is not implemented.

    Measured 5.1e-6 m on a 1.0e-2 m layer. The bound here is 1% of the layer,
    about twice what is observed.
    """
    d = build(base_depth=0.01, bedload=True, classes=(('sand', 8.0e-4),))
    d.evolve_to_end(finaltime=40.0)
    z = d.quantities['elevation'].centroid_values
    overshoot = (d.sediment_z_base - z).max()
    assert overshoot <= 0.01 * 0.01, 'overshoot %.3e m' % overshoot


# ---------------------------------------------------------------- classes

def test_registration_order_does_not_change_the_answer():
    """The erodible thickness belongs to the CELL, not to a class, so the
    classes are limited together by one shared proportional factor. Serving
    them in order would make the answer depend on the order add_sediment_fraction
    was called, which is not physics."""
    first = build(base_depth=0.004, classes=TWO_CLASSES)
    first.evolve_to_end(finaltime=30.0)
    second = build(base_depth=0.004, classes=TWO_CLASSES[::-1])
    second.evolve_to_end(finaltime=30.0)

    assert np.allclose(first.quantities['elevation'].centroid_values,
                       second.quantities['elevation'].centroid_values,
                       rtol=0, atol=1e-14)
    a_fine = float((first.tracer_conserved_values[0] * first.areas).sum())
    a_coarse = float((first.tracer_conserved_values[1] * first.areas).sum())
    b_coarse = float((second.tracer_conserved_values[0] * second.areas).sum())
    b_fine = float((second.tracer_conserved_values[1] * second.areas).sum())
    assert abs(a_fine - b_fine) < 1e-12 * max(abs(a_fine), 1e-30)
    assert abs(a_coarse - b_coarse) < 1e-12 * max(abs(a_coarse), 1e-30)
    # and the budget was genuinely contested, or the test proves nothing
    assert a_fine > 0.0 and a_coarse > 0.0
    assert int((first.quantities['elevation'].centroid_values
                - first.sediment_z_base < 1e-12).sum()) > 0


# ---------------------------------------------------------------- off by default

def test_an_unreachable_base_reproduces_the_no_base_answer_bitwise():
    """The regression gate: this feature must be invisible when it does not
    bite."""
    off = build()
    off.evolve_to_end(finaltime=25.0)
    deep = build(base_depth=1000.0)
    deep.evolve_to_end(finaltime=25.0)
    assert np.array_equal(off.quantities['elevation'].centroid_values,
                          deep.quantities['elevation'].centroid_values)
    assert np.array_equal(off.tracer_conserved_values[0],
                          deep.tracer_conserved_values[0])


# ---------------------------------------------------------------- interface

def test_a_base_above_the_bed_is_rejected():
    with pytest.raises(ValueError):
        build().set_erodible_base(elevation=1000.0)


def test_elevation_and_depth_together_are_rejected():
    with pytest.raises(ValueError):
        build().set_erodible_base(elevation=0.0, depth=1.0)


def test_the_base_can_be_removed():
    d = build(base_depth=0.1)
    d.set_erodible_base()
    assert d.sediment_has_z_base == 0
    assert d.sediment_z_base is None


def test_erodible_thickness_reports_the_layer():
    d = build(base_depth=0.25)
    assert np.allclose(d.erodible_thickness(), 0.25)


def test_the_summary_reports_the_base():
    d = build(base_depth=0.25)
    assert 'erodible base' in d.sediment_summary()
    assert '[L-5]' in d.sediment_summary()


# ---------------------------------------------------------------- region

def test_only_the_named_region_erodes():
    d = build()
    zc = d.quantities['elevation'].centroid_values.copy()
    x = d.centroid_coordinates[:, 0]
    d.set_erodible_region(polygon=HALF)
    d.evolve_to_end(finaltime=30.0)
    z = d.quantities['elevation'].centroid_values
    assert (zc - z)[x >= 30.0].max() <= 0.0
    assert (zc - z)[x < 30.0].max() > 1e-4


def test_erodible_false_locks_the_named_region_instead():
    d = build()
    zc = d.quantities['elevation'].centroid_values.copy()
    x = d.centroid_coordinates[:, 0]
    d.set_erodible_region(polygon=HALF, erodible=False)
    d.evolve_to_end(finaltime=30.0)
    z = d.quantities['elevation'].centroid_values
    assert (zc - z)[x < 30.0].max() <= 0.0
    assert (zc - z)[x >= 30.0].max() > 1e-4


def test_locked_cells_still_accrete():
    """Locked means UNSCOURABLE, not inert: sediment settles onto a concrete
    apron in the field, and that material is erodible again because it sits
    above the base.

    Needs a depositional case to show. In an erosive channel the limiter
    correctly holds such cells at exactly net zero, scaling erosion back until
    it just cancels deposition -- sand does not pile up on a scoured apron.
    """
    d = rectangular_cross_domain(30, 8, len1=60.0, len2=16.0)
    d.set_flow_algorithm('DE0')
    d.set_low_froude(0)
    d.store = False
    d.set_quantity('elevation', 0.0)
    d.set_quantity('stage', 1.0)
    d.set_quantity('friction', 0.03)
    d.set_boundary({t: Reflective_boundary(d) for t in d.get_boundary_tags()})
    d.set_sediment_parameters(porosity=POROSITY)
    d.add_sediment_fraction(name='silt', diameter=6.0e-5, tau_c_star=1.0e6,
                     initial_concentration=0.01)
    zc = d.quantities['elevation'].centroid_values.copy()
    d.set_erodible_region(polygon=HALF)
    locked = d.centroid_coordinates[:, 0] >= 30.0
    d.evolve_to_end(finaltime=30.0)
    z = d.quantities['elevation'].centroid_values
    assert (z - zc)[locked].min() > 0.0
    assert (d.erodible_thickness()[locked] > 0.0).all()


def test_the_budget_closes_with_a_region_restriction():
    d = build()
    z0 = d.quantities['elevation'].centroid_values.copy()
    m0 = _mass(d)
    d.set_erodible_region(polygon=HALF)
    d.evolve_to_end(finaltime=30.0)
    z1 = d.quantities['elevation'].centroid_values
    bed = float(((1.0 - POROSITY) * (z1 - z0) * d.areas).sum())
    assert abs((_mass(d) + bed) - m0) < 1e-12 * max(abs(_mass(d)), 1e-30)


def test_the_base_and_the_region_compose():
    """The region says WHERE erosion may act, the base says HOW DEEP. Setting
    one must not discard the other."""
    d = build()
    zc = d.quantities['elevation'].centroid_values.copy()
    x = d.centroid_coordinates[:, 0]
    d.set_erodible_base(depth=0.01)
    d.set_erodible_region(polygon=HALF)
    d.evolve_to_end(finaltime=40.0)
    z = d.quantities['elevation'].centroid_values
    assert (zc - z)[x >= 30.0].max() <= 0.0
    assert 1e-4 < (zc - z)[x < 30.0].max() <= 0.01 + 1e-12


def test_the_order_they_are_set_in_does_not_matter():
    d = build()
    d.set_erodible_region(polygon=HALF)
    d.set_erodible_base(depth=0.01)
    assert d._sediment_erodible_mask is not None
    assert d._sediment_user_base is not None
    assert (int((d.erodible_thickness() > 0).sum())
            == int(d._sediment_erodible_mask.sum()))


def test_a_circular_region_works():
    d = build()
    d.set_erodible_region(center=[15.0, 8.0], radius=6.0)
    n = int(d._sediment_erodible_mask.sum())
    assert 0 < n < len(d)


def test_a_region_selecting_nothing_is_rejected():
    """Nearly always a polygon in the wrong coordinates, and either sense of
    the flag then gives a plausible run doing the opposite of what was meant."""
    with pytest.raises(ValueError):
        build().set_erodible_region(polygon=[[900.0, 900.0], [910.0, 900.0],
                                             [910.0, 910.0], [900.0, 910.0]])


def test_a_region_object_selects_the_same_cells_as_its_polygon():
    """Passing a Region is the general form -- it also understands line=,
    poly= and expand_polygon=, which have no keyword on the setter."""
    a = build()
    a.set_erodible_region(Region(a, polygon=HALF))
    b = build()
    b.set_erodible_region(polygon=HALF)
    assert np.array_equal(a._sediment_erodible_mask, b._sediment_erodible_mask)


def test_a_region_built_on_another_domain_is_refused():
    d = build()
    with pytest.raises(ValueError):
        d.set_erodible_region(Region(build(), polygon=HALF))


def test_a_bare_list_of_points_is_refused():
    """Taking it as a region would select every cell and look like it worked."""
    with pytest.raises(TypeError):
        build().set_erodible_region(HALF)


def test_a_region_plus_build_arguments_is_refused():
    d = build()
    with pytest.raises(ValueError):
        d.set_erodible_region(Region(d, polygon=HALF), polygon=HALF)


def test_the_region_restriction_can_be_removed():
    d = build()
    d.set_erodible_region(polygon=HALF)
    d.set_erodible_region()
    assert d._sediment_erodible_mask is None
    assert d.sediment_has_z_base == 0


def test_the_summary_reports_the_region():
    d = build()
    d.set_erodible_region(polygon=HALF)
    assert 'erodible region' in d.sediment_summary()
    assert 'locked' in d.sediment_summary()

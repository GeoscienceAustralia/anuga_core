"""Angle-of-repose relaxation (PHYSICS_SPEC 7), after FG21 §2.2.4.

Where the centroid-to-centroid bed slope exceeds a critical angle, bed material
is diffused downslope until it does not. FG21 are explicit that this is a
NUMERICAL HEURISTIC and not physics -- real slope failures are advective -- and
that it suppresses knickpoint retreat that may be real. It exists to stop the
rest of the model breaking on over-steep slopes, and it is off by default.

Spec 7 sets four requirements; three are tested here. Mass conservation is the
one that separates this from ANUGA's older sanddune_erosion_operator, which
lowers an over-steep cell and lets the material vanish.

The mode 1 / mode 2 comparison is in test_sediment_gpu.py.
"""
import math

import numpy as np
import pytest

from anuga import Reflective_boundary, rectangular_cross_domain

ANGLE = 30.0
# The kernel declares convergence within a 1e-3 relative tolerance on the
# threshold SLOPE, because convergence is asymptotic and a strict test never
# terminates. At 30 degrees that admits at most 30.03.
TOL_DEG = math.degrees(math.atan(math.tan(math.radians(ANGLE)) * 1.001))


def cone(x, y, height=6.0, radius=8.0, cx=30.0, cy=8.0):
    """A cone at about 37 degrees -- steeper than any repose angle here."""
    r = np.sqrt((x - cx) ** 2 + (y - cy) ** 2)
    return np.where(r < radius, height * (1.0 - r / radius), 0.0)


def build(angle=None, relax=1.0, max_sweeps=50, base_depth=None, n_x=40):
    d = rectangular_cross_domain(n_x, 12, len1=60.0, len2=16.0)
    d.set_flow_algorithm('DE0')
    d.set_low_froude(0)
    d.store = False
    x, y = d.centroid_coordinates[:, 0], d.centroid_coordinates[:, 1]
    d.set_quantity('elevation', cone(x, y), location='centroids')
    d.set_quantity('stage', -1.0)          # dry: isolates the bed process
    d.set_quantity('friction', 0.03)
    d.set_boundary({t: Reflective_boundary(d) for t in d.get_boundary_tags()})
    d.set_sediment_parameters(porosity=0.3)
    d.add_sediment_fraction(name='sand', diameter=2.0e-4)
    if angle is not None:
        d.set_angle_of_repose(angle, relax=relax, max_sweeps=max_sweeps)
    if base_depth is not None:
        d.set_erodible_base(depth=base_depth)
    return d


def max_slope_degrees(d):
    z = d.quantities['elevation'].centroid_values
    cc = d.centroid_coordinates
    nb = d.neighbours
    worst = 0.0
    for i in range(3):
        j = nb[:, i]
        ok = j >= 0
        dz = np.abs(z[ok] - z[j[ok]])
        dist = np.sqrt(((cc[ok] - cc[j[ok]]) ** 2).sum(axis=1))
        good = dist > 0
        worst = max(worst, float((dz[good] / dist[good]).max()))
    return math.degrees(math.atan(worst))


def _operator(d):
    return [o for o in d.fractional_step_operators
            if type(o).__name__ == 'Sediment_transport_operator'][0]


def test_the_cone_starts_steeper_than_the_critical_angle():
    assert max_slope_degrees(build()) > ANGLE + 1.0


@pytest.mark.slow
def test_an_over_steep_bed_is_relaxed_to_the_critical_angle():
    """A cold start from a badly over-steep bed is the slow case: this is an
    explicit diffusion solve and it took 793 sweeps to converge. A running
    model needs a handful per step."""
    d = build(angle=ANGLE, max_sweeps=2000)
    d.evolve_to_end(finaltime=1.0)
    assert max_slope_degrees(d) <= TOL_DEG + 1e-3


@pytest.mark.slow
def test_relaxation_lowers_the_peak_without_flattening_everything():
    before = build()
    d = build(angle=ANGLE, max_sweeps=2000)
    d.evolve_to_end(finaltime=1.0)
    z0 = before.quantities['elevation'].centroid_values
    z1 = d.quantities['elevation'].centroid_values
    assert z1.max() < z0.max()
    assert z1.max() > 0.5 * z0.max()


@pytest.mark.slow
def test_bed_volume_is_conserved():
    """Material removed from an over-steep cell is DEPOSITED ON ITS NEIGHBOURS,
    never discarded. Structural rather than incidental: the transfer is
    computed per edge from data both cells share, so both compute the identical
    volume and the pair balances exactly."""
    d = build(angle=ANGLE, max_sweeps=2000)
    v0 = float((d.quantities['elevation'].centroid_values * d.areas).sum())
    d.evolve_to_end(finaltime=1.0)
    v1 = float((d.quantities['elevation'].centroid_values * d.areas).sum())
    assert abs(v1 - v0) < 1e-11 * max(abs(v0), 1.0)


def test_a_cap_too_small_to_converge_is_hit_and_reported():
    """Spec 7 requires the cap be reported rather than swallowed: hitting it
    means the bed may still be over-steep."""
    d = build(angle=ANGLE, max_sweeps=3)
    op = _operator(d)
    d.evolve_to_end(finaltime=0.2)
    assert op.repose_cap_hits > 0
    assert op.repose_sweeps <= 3
    assert max_slope_degrees(d) > ANGLE, 'the cap did not actually bind'


@pytest.mark.slow
def test_with_a_workable_cap_it_converges_and_stops_early():
    d = build(angle=ANGLE, max_sweeps=2000)
    op = _operator(d)
    d.evolve_to_end(finaltime=1.0)
    assert op.repose_cap_hits == 0
    assert op.repose_sweeps < 2000


@pytest.mark.slow
def test_relaxation_cannot_slump_a_cell_below_its_erodible_base():
    """A layer thinner than the relaxation wants, or the test proves nothing:
    unconstrained this cone lowers its peak only 0.3425 m, so a 0.5 m layer
    never binds and a 'constrained' run reproduces the free one exactly."""
    free = build(angle=ANGLE, max_sweeps=2000)
    zc = free.quantities['elevation'].centroid_values.copy()
    free.evolve_to_end(finaltime=1.0)
    drop_free = (zc - free.quantities['elevation'].centroid_values).max()

    d = build(angle=ANGLE, max_sweeps=2000, base_depth=0.1)
    base = d.sediment_z_base.copy()
    d.evolve_to_end(finaltime=1.0)
    z = d.quantities['elevation'].centroid_values
    assert (base - z).max() <= 0.0
    assert (zc - z).max() < drop_free - 1e-6
    assert int((d.erodible_thickness() <= 1e-9).sum()) > 0


def test_off_by_default_changes_nothing():
    off = build()
    off.evolve_to_end(finaltime=1.0)
    none = build(angle=None)
    none.evolve_to_end(finaltime=1.0)
    assert np.array_equal(off.quantities['elevation'].centroid_values,
                          none.quantities['elevation'].centroid_values)
    assert max_slope_degrees(off) > ANGLE, 'the cone should survive untouched'


@pytest.mark.parametrize('angle', [0.0, 90.0, -5.0])
def test_out_of_range_angles(angle):
    d = build()
    if angle == 0.0:
        d.set_angle_of_repose(angle=0.0)      # zero disables
        assert d.sediment_repose_tan == 0.0
    else:
        with pytest.raises(ValueError):
            d.set_angle_of_repose(angle)


def test_relax_and_max_sweeps_are_validated():
    d = build()
    with pytest.raises(ValueError):
        d.set_angle_of_repose(30.0, relax=0.0)
    with pytest.raises(ValueError):
        d.set_angle_of_repose(30.0, max_sweeps=0)


def test_the_angle_round_trips_and_is_reported():
    d = build()
    d.set_angle_of_repose(35.0)
    assert abs(math.degrees(math.atan(d.sediment_repose_tan)) - 35.0) < 1e-9
    assert 'angle of repose' in d.sediment_summary()
    assert '35.0 degrees' in d.sediment_summary()
    d.set_angle_of_repose(None)
    assert d.sediment_repose_tan == 0.0

"""Inlet operators carry tracer mass (inlet_tracers.py).

Before this, an Inlet_operator moved only water: an inflow came in clean, and
an outflow left its tracer behind to concentrate in the outlet cells. Water
added now carries the given inflow concentration (0 by default), and water
removed takes the inlet's mean concentration with it. Both are exact to
roundoff, which is what these tests check, in both compute modes.
"""
import numpy as np
import pytest

import anuga
from anuga import Reflective_boundary, rectangular_cross_domain


def basin(mode='legacy', depth=1.0, c0=0.01):
    """Closed flat basin at rest, one tracer at a uniform concentration."""
    d = rectangular_cross_domain(10, 10, len1=100.0, len2=100.0)
    d.set_flow_algorithm('DE1')
    d.set_compute_mode(mode)
    d.store = False
    d.set_quantity('elevation', 0.0)
    d.set_quantity('stage', depth)
    d.set_quantity('friction', 0.0)
    d.set_boundary({t: Reflective_boundary(d) for t in d.get_boundary_tags()})
    d.add_tracer('dye', initial_value=c0)
    return d


def mass(d):
    return float((d.tracer_conserved_values[0] * d.areas).sum())


def water(d):
    return float(((d.quantities['stage'].centroid_values
                   - d.quantities['elevation'].centroid_values) * d.areas).sum())


REGION = dict(center=(50.0, 50.0), radius=15.0)


@pytest.mark.parametrize('mode', ['legacy', 'unified'])
def test_extraction_takes_the_tracer_with_the_water(mode):
    """Uniform c0 everywhere: pumping water out must remove c0 per unit of
    water and leave c0 in the inlet cells, not concentrate it there."""
    c0 = 0.01
    d = basin(mode, c0=c0)
    region = anuga.Region(d, **REGION)
    op = anuga.Inlet_operator(d, region, Q=-50.0)
    m0, w0 = mass(d), water(d)
    d.evolve_to_end(finaltime=60.0)
    removed_water = w0 - water(d)
    assert removed_water > 1000.0
    assert mass(d) == pytest.approx(m0 - c0 * removed_water, rel=1e-10)
    c = d.tracer_centroid_values[0][region.indices]
    assert np.allclose(c, c0, rtol=1e-9)
    assert op.tracers.total_out['dye'] == pytest.approx(c0 * removed_water, rel=1e-10)


@pytest.mark.parametrize('mode', ['legacy', 'unified'])
def test_inflow_brings_its_concentration(mode):
    """Clean basin, inflow at c_in: the tracer added is c_in times the water
    added, exactly."""
    c_in = 0.02
    d = basin(mode, c0=0.0)
    op = anuga.Inlet_operator(d, anuga.Region(d, **REGION), Q=20.0,
                              tracer_concentrations={'dye': c_in})
    w0 = water(d)
    d.evolve_to_end(finaltime=60.0)
    added = water(d) - w0
    assert added == pytest.approx(20.0 * 60.0, rel=1e-6)
    assert mass(d) == pytest.approx(c_in * added, rel=1e-10)
    assert op.tracers.total_in['dye'] == pytest.approx(c_in * added, rel=1e-10)


@pytest.mark.parametrize('mode', ['legacy', 'unified'])
def test_inflow_without_a_concentration_is_clean_water(mode):
    """The default is the old behaviour for inflow: no tracer is added."""
    d = basin(mode, c0=0.01)
    anuga.Inlet_operator(d, anuga.Region(d, **REGION), Q=20.0)
    m0 = mass(d)
    d.evolve_to_end(finaltime=60.0)
    assert mass(d) == pytest.approx(m0, rel=1e-12)


def test_a_time_varying_concentration():
    """A callable is averaged over each step, like Q; with c_in linear in t
    the added mass is the exact integral."""
    d = basin('legacy', c0=0.0)
    anuga.Inlet_operator(d, anuga.Region(d, **REGION), Q=10.0,
                         tracer_concentrations={'dye': lambda t: 1.0e-3 * t})
    d.evolve_to_end(finaltime=60.0)
    # integral of Q c dt = 10 * 1e-3 * 60^2 / 2
    assert mass(d) == pytest.approx(10.0 * 1.0e-3 * 60.0 ** 2 / 2.0, rel=1e-9)


@pytest.mark.parametrize('mode', ['legacy', 'unified'])
def test_draining_the_inlet_dry_removes_all_its_tracer(mode):
    """Asking for more water than the inlet holds empties it: the tracer in
    those cells goes too, and nothing is left behind in a dry cell."""
    d = basin(mode, depth=0.05, c0=0.01)
    region = anuga.Region(d, **REGION)
    anuga.Inlet_operator(d, region, Q=-1.0e5)
    d.evolve_to_end(finaltime=5.0)
    h = (d.quantities['stage'].centroid_values
         - d.quantities['elevation'].centroid_values)[region.indices]
    m = d.tracer_conserved_values[0][region.indices]
    assert np.all(m <= 0.01 * np.maximum(h, 0.0) + 1e-15)
    assert np.isfinite(d.tracer_centroid_values).all()
    assert d.tracer_centroid_values.max() <= 0.01 * (1 + 1e-9)


def test_an_unknown_tracer_is_rejected():
    d = basin('legacy')
    op = anuga.Inlet_operator(d, anuga.Region(d, **REGION), Q=1.0,
                              tracer_concentrations={'salt': 1.0})
    with pytest.raises(ValueError):
        d.evolve_to_end(finaltime=1.0)


def test_the_two_compute_modes_agree():
    """Extraction then inflow at a concentration, mode 1 against mode 2."""
    out = []
    for mode in ('legacy', 'unified'):
        d = basin(mode, c0=0.01)
        anuga.Inlet_operator(d, anuga.Region(d, center=(30.0, 30.0), radius=12.0), Q=-20.0)
        anuga.Inlet_operator(d, anuga.Region(d, center=(70.0, 70.0), radius=12.0), Q=20.0,
                             tracer_concentrations={'dye': 0.05})
        d.evolve_to_end(finaltime=60.0)
        out.append(d.tracer_conserved_values[0].copy())
    assert np.abs(out[0] - out[1]).max() < 1e-12

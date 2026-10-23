"""Rate operators carry tracers with the water they add and remove.

Positive rates add water at `tracer_concentrations` (0 by default: clean
rain). Negative rates either take each tracer with the water at the cell's
concentration ('carry': pumping, drains) or leave it behind ('retain', the
default: salt under evaporation). Every rule is exact per cell.
"""
import numpy as np
import pytest

import anuga
from anuga import Reflective_boundary, rectangular_cross_domain
from anuga import Rate_operator


def basin(mode='legacy', depth=1.0, tracers=(('dye', 0.01),)):
    d = rectangular_cross_domain(10, 10, len1=100.0, len2=100.0)
    d.set_flow_algorithm('DE1')
    d.set_compute_mode(mode)
    d.store = False
    d.set_quantity('elevation', 0.0)
    d.set_quantity('stage', depth)
    d.set_quantity('friction', 0.0)
    d.set_boundary({t: Reflective_boundary(d) for t in d.get_boundary_tags()})
    for name, c0 in tracers:
        d.add_tracer(name, initial_value=c0)
    return d


def mass(d, s=0):
    return float((d.tracer_conserved_values[s] * d.areas).sum())


def water(d):
    return float(((d.quantities['stage'].centroid_values
                   - d.quantities['elevation'].centroid_values) * d.areas).sum())


MODES = ['legacy', 'unified']


@pytest.mark.parametrize('mode', MODES)
def test_rain_brings_its_concentration(mode):
    d = basin(mode, tracers=(('dye', 0.0),))
    Rate_operator(d, rate=1.0e-4, tracer_concentrations={'dye': 0.03})
    w0 = water(d)
    d.evolve_to_end(finaltime=20.0)
    assert mass(d) == pytest.approx(0.03 * (water(d) - w0), rel=1e-10)


@pytest.mark.parametrize('mode', MODES)
def test_clean_rain_is_the_default(mode):
    d = basin(mode)
    Rate_operator(d, rate=1.0e-4)
    m0 = mass(d)
    d.evolve_to_end(finaltime=20.0)
    assert mass(d) == pytest.approx(m0, rel=1e-12)


@pytest.mark.parametrize('mode', MODES)
def test_extraction_retains_by_default(mode):
    """The old behaviour, and right for salt under evaporation."""
    d = basin(mode)
    Rate_operator(d, rate=-1.0e-4, center=(50.0, 50.0), radius=20.0)
    m0 = mass(d)
    d.evolve_to_end(finaltime=20.0)
    assert mass(d) == pytest.approx(m0, rel=1e-12)


@pytest.mark.parametrize('mode', MODES)
def test_extraction_can_carry_the_tracer_away(mode):
    d = basin(mode, tracers=(('dye', 0.01),))
    op = Rate_operator(d, rate=-1.0e-4, tracer_extraction='carry')
    m0, w0 = mass(d), water(d)
    d.evolve_to_end(finaltime=20.0)
    removed = w0 - water(d)
    assert removed > 1.0
    assert mass(d) == pytest.approx(m0 - 0.01 * removed, rel=1e-10)
    assert np.allclose(d.tracer_centroid_values[0], 0.01, rtol=1e-9)


@pytest.mark.parametrize('mode', MODES)
def test_draining_dry_with_carry_leaves_nothing_behind(mode):
    d = basin(mode, depth=0.01)
    Rate_operator(d, rate=-1.0e-2, tracer_extraction='carry')
    d.evolve_to_end(finaltime=5.0)
    assert mass(d) == pytest.approx(0.0, abs=1e-14)


@pytest.mark.parametrize('mode', MODES)
def test_one_choice_per_tracer(mode):
    d = basin(mode, tracers=(('salt', 0.02), ('dye', 0.01)))
    Rate_operator(d, rate=-1.0e-4, tracer_extraction={'dye': 'carry'})
    m_salt0, m_dye0, w0 = mass(d, 0), mass(d, 1), water(d)
    d.evolve_to_end(finaltime=20.0)
    removed = w0 - water(d)
    assert mass(d, 0) == pytest.approx(m_salt0, rel=1e-12)            # retained
    assert mass(d, 1) == pytest.approx(m_dye0 - 0.01 * removed, rel=1e-10)


def test_a_spatial_rate_of_both_signs():
    """Rain on one half at c_in, pumping from the other with carry: the budget
    is exact cell by cell (host path: spatial rates are not offloaded)."""
    d = basin('legacy', tracers=(('dye', 0.01),))
    Rate_operator(d, rate=lambda x, y: np.where(x < 50.0, 2.0e-4, -1.0e-4),
                  tracer_concentrations={'dye': 0.05}, tracer_extraction='carry')
    x = d.centroid_coordinates[:, 0]
    area = d.areas
    m0 = mass(d)
    h0 = (d.quantities['stage'].centroid_values - d.quantities['elevation'].centroid_values)
    d.evolve_to_end(finaltime=1.0e-3)       # one tiny step: no advection to speak of
    h1 = (d.quantities['stage'].centroid_values - d.quantities['elevation'].centroid_values)
    added = float(((h1 - h0) * area)[x < 50.0].sum())
    removed = -float(((h1 - h0) * area)[x > 50.0].sum())
    assert mass(d) == pytest.approx(m0 + 0.05 * added - 0.01 * removed, rel=1e-8)


@pytest.mark.parametrize('mode', MODES)
def test_the_factories_pass_the_options_through(mode):
    d = basin(mode, tracers=(('dye', 0.0),))
    Rate_operator.rainfall(d, rate=360.0, tracer_concentrations={'dye': 0.02})   # 360 mm/hr
    w0 = water(d)
    d.evolve_to_end(finaltime=10.0)
    assert mass(d) == pytest.approx(0.02 * (water(d) - w0), rel=1e-10)


def test_bad_settings_are_rejected():
    d = basin('legacy')
    with pytest.raises(ValueError):
        Rate_operator(d, rate=-1.0, tracer_extraction='drain')
    op = Rate_operator(d, rate=-1.0e-4, tracer_extraction={'salt': 'carry'})
    with pytest.raises(ValueError):
        d.evolve_to_end(finaltime=1.0)


@pytest.mark.parametrize('kind', ['scalar', 'array'])
def test_the_two_compute_modes_agree(kind):
    out = []
    for mode in MODES:
        d = basin(mode, tracers=(('dye', 0.01),))
        x = d.centroid_coordinates[:, 0]
        rate = -1.0e-4 if kind == 'scalar' else np.where(x < 50.0, 2.0e-4, -1.0e-4)
        Rate_operator(d, rate=rate, tracer_concentrations={'dye': 0.05},
                      tracer_extraction='carry')
        d.evolve_to_end(finaltime=20.0)
        out.append(d.tracer_conserved_values[0].copy())
    assert np.abs(out[0] - out[1]).max() < 1e-10

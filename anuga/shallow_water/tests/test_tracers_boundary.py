"""Tracer boundary values: prescribed on inflow, interior value on outflow.

The flux kernel already picks the upwind value edge by edge from the sign of
the water flux, so the only thing a user can supply -- and the only thing that
is physically well posed -- is the concentration arriving on an INFLOW. An
outflow takes the interior edge value, which is the transmissive condition, and
prescribing there would over-determine the advection.
"""

import numpy as num
import pytest

import anuga


def _domain(n=6):
    d = anuga.rectangular_cross_domain(n, n)
    d.set_quantity('elevation', 0.0)
    d.set_quantity('stage', 1.0)
    b = anuga.Reflective_boundary(d)
    d.set_boundary({'left': b, 'right': b, 'top': b, 'bottom': b})
    d.add_tracer('salinity', initial_value=0.0)
    return d


def test_unset_boundaries_bring_clean_water():
    # Documented default: the array is zero-filled, so an unset inflow carries
    # c = 0. Pinned so a change to it has to be deliberate.
    d = _domain()
    assert num.allclose(d.get_tracer_boundary('salinity', 'left'), 0.0)


def test_a_scalar_applies_to_every_edge_of_the_tag():
    d = _domain()
    d.set_tracer_boundary('salinity', 'left', 0.03)
    assert num.allclose(d.get_tracer_boundary('salinity', 'left'), 0.03)
    # and only that tag
    assert num.allclose(d.get_tracer_boundary('salinity', 'right'), 0.0)


def test_an_array_sets_each_edge():
    d = _domain()
    n_left = len(d.tag_boundary_cells['left'])
    vals = num.linspace(0.0, 1.0, n_left)
    d.set_tracer_boundary('salinity', 'left', vals)
    assert num.allclose(d.get_tracer_boundary('salinity', 'left'), vals)


def test_a_wrong_length_array_is_refused():
    d = _domain()
    n_left = len(d.tag_boundary_cells['left'])
    with pytest.raises(ValueError, match='one per edge'):
        d.set_tracer_boundary('salinity', 'left', num.zeros(n_left + 1))


def test_an_unknown_tag_is_refused():
    d = _domain()
    with pytest.raises(ValueError, match='no boundary tagged'):
        d.set_tracer_boundary('salinity', 'nonexistent', 0.03)


def test_an_unknown_tracer_is_refused():
    d = _domain()
    with pytest.raises(Exception):
        d.set_tracer_boundary('not_a_tracer', 'left', 0.03)


def test_a_callable_is_evaluated_now_and_on_update():
    d = _domain()
    d.set_tracer_boundary('salinity', 'left', lambda t: 0.1 * (t + 1.0))
    # evaluated immediately at t = 0
    assert num.allclose(d.get_tracer_boundary('salinity', 'left'), 0.1)
    d.set_time(4.0)
    d.update_tracer_boundary_values()
    assert num.allclose(d.get_tracer_boundary('salinity', 'left'), 0.5)


def test_a_constant_replaces_a_previous_callable():
    d = _domain()
    d.set_tracer_boundary('salinity', 'left', lambda t: 0.1 * (t + 1.0))
    d.set_tracer_boundary('salinity', 'left', 0.02)
    d.set_time(9.0)
    d.update_tracer_boundary_values()      # must not re-apply the old callable
    assert num.allclose(d.get_tracer_boundary('salinity', 'left'), 0.02)


def test_boundary_values_survive_adding_another_tracer():
    # add_tracer REALLOCATES every tracer array; earlier rows must be kept.
    d = _domain()
    d.set_tracer_boundary('salinity', 'left', 0.03)
    d.add_tracer('dye', initial_value=0.0)
    assert num.allclose(d.get_tracer_boundary('salinity', 'left'), 0.03)
    assert num.allclose(d.get_tracer_boundary('dye', 'left'), 0.0)


def test_outflow_does_not_take_the_boundary_value():
    """The characteristic condition: an outflow ignores whatever is prescribed.

    Interior tracer is 0.5 everywhere and the domain drains outward, so if the
    boundary value were applied on outflow the interior would be dragged toward
    0.9. It must not be.
    """
    d = anuga.rectangular_cross_domain(10, 10)
    d.set_quantity('elevation', lambda x, y: -0.2 * x)
    d.set_quantity('stage', lambda x, y: num.maximum(0.4 - 0.2 * x, 0.02))
    Bt = anuga.Transmissive_boundary(d)
    d.set_boundary({'left': Bt, 'right': Bt, 'top': Bt, 'bottom': Bt})
    d.add_tracer('c', initial_value=0.5)
    d.set_tracer_boundary('c', 'right', 0.9)     # the downstream, outflow side

    for _ in d.evolve(yieldstep=0.5, finaltime=2.0):
        pass

    c = d.get_tracer('c')
    wet = (d.quantities['stage'].centroid_values
           - d.quantities['elevation'].centroid_values) > 1e-3
    assert c[wet].max() <= 0.5 + 1e-6, \
        'an outflow boundary must not inject its prescribed concentration'


def test_inflow_brings_the_prescribed_concentration():
    """The other half: water entering must carry what was set.

    Flat bed, interior initially clean, and a Dirichlet boundary on the left
    holding a higher stage so water flows IN there. The prescribed
    concentration must appear in the domain.
    """
    d = anuga.rectangular_cross_domain(10, 10)
    d.set_quantity('elevation', 0.0)
    d.set_quantity('stage', 0.2)
    Bin = anuga.Dirichlet_boundary([0.6, 0.4, 0.0])   # higher stage, flowing in
    Bout = anuga.Transmissive_boundary(d)
    Br = anuga.Reflective_boundary(d)
    d.set_boundary({'left': Bin, 'right': Bout, 'top': Br, 'bottom': Br})

    d.add_tracer('c', initial_value=0.0)
    d.set_tracer_boundary('c', 'left', 1.0)

    # Stopped part-way: given long enough the inflow floods the whole domain
    # to c = 1, which is correct but says nothing about WHERE it entered.
    for _ in d.evolve(yieldstep=0.1, finaltime=0.4):
        pass

    c = d.get_tracer('c')
    x = d.get_centroid_coordinates()[:, 0]

    assert c.max() > 0.05, \
        'nothing entered: the inflow boundary concentration was not applied'
    assert c.max() <= 1.0 + 1e-6, \
        'concentration exceeded what the boundary supplied'
    assert c[x < 0.3].mean() > c[x > 0.7].mean(), \
        'the tracer should be strongest near the boundary it entered through'

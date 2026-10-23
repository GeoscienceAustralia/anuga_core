"""Running maxima of the tracers, alongside the water ones.

Collect_max_quantities_operator tracks max stage/depth/speed/uh; a tracer --
and so a suspended sediment class -- gets the same treatment, stored as
max_<name>_c.

The subtlety these tests exist for: the operator must DERIVE c = m/h rather
than read domain.tracer_centroid_values. c is only recomputed during
extrapolation, at the START of a step, so by the time a fractional-step
operator runs it is one step behind the conserved m that was just written.
max_stage has no such problem, stage being conserved itself.

These are white-box tests: they read op.max_tracer, which is HOST state. In
mode 2 the maxima live on the device and that array is only filled when the
operator stores to the sww, so every domain here pins itself to legacy -- the
pattern CLAUDE.md describes for tests that inspect internal update arrays. The
mode-2 path is covered in test_tracers_gpu.py, which syncs first.
"""

import numpy as num
import pytest

import anuga
from anuga.operators.collect_max_quantities_operator import (
    Collect_max_quantities_operator)

netCDF4 = pytest.importorskip('netCDF4')


def _domain(n=6):
    """A plume that washes out, so the max and the final value differ."""
    d = anuga.rectangular_cross_domain(n, n, len1=60.0, len2=60.0)
    d.set_quantity('elevation', 0.0)
    d.set_quantity('stage', 1.0)
    d.set_quantity('xmomentum', 2.0)
    Br = anuga.Reflective_boundary(d)
    Bt = anuga.Transmissive_boundary(d)
    d.set_boundary({'left': anuga.Dirichlet_boundary([1.0, 2.0, 0.0]),
                    'top': Br, 'bottom': Br, 'right': Bt})
    d.add_tracer('salinity')
    x = d.centroid_coordinates[:, 0]
    d.set_tracer('salinity', num.where(x < 20.0, 1.0, 0.0))
    d.set_tracer_boundary('salinity', 'left', 0.0)   # clean water follows
    d.store = False
    d.set_compute_mode('legacy')     # reads host-side operator state
    return d


def _run(d, finaltime=5.0):
    for _ in d.evolve(yieldstep=0.5, finaltime=finaltime):
        pass


def test_no_tracers_means_no_tracer_maxima():
    d = anuga.rectangular_cross_domain(4, 4)
    b = anuga.Reflective_boundary(d)
    d.set_boundary({'left': b, 'right': b, 'top': b, 'bottom': b})
    d.set_compute_mode('legacy')
    op = Collect_max_quantities_operator(d)
    assert op.max_tracer.shape == (0, len(d))


def test_the_maximum_is_never_below_the_current_value():
    """The property the whole thing rests on, and the one the stale-c bug broke."""
    d = _domain()
    op = Collect_max_quantities_operator(d)
    _run(d)
    c = d.get_tracer('salinity')
    assert num.all(op.max_tracer[0] >= c - 1e-12), \
        'the running maximum is below the final concentration in %d cells' \
        % int((c > op.max_tracer[0] + 1e-12).sum())


def test_the_maximum_records_a_plume_that_has_passed():
    """Otherwise it is just a copy of the final state."""
    d = _domain()
    op = Collect_max_quantities_operator(d)
    _run(d)
    c = d.get_tracer('salinity')
    assert num.any(op.max_tracer[0] > c + 1e-3), \
        'the maximum never exceeds the final value, so nothing was tracked'


def test_it_beats_sampling_at_yield_steps():
    """The operator runs every timestep, so it cannot be lower than this.

    Sampling starts after the FIRST yield: the operator collects on the
    fractional-step call, which does not run before the initial state, so
    the initial condition is outside what either side is measuring here.
    That exclusion is documented on the operator and is how the existing
    max_stage has always behaved.
    """
    d = _domain()
    op = Collect_max_quantities_operator(d)
    sampled = None
    for _ in d.evolve(yieldstep=0.5, finaltime=5.0):
        if sampled is None:
            sampled = num.zeros(len(d))     # skip the t=0 yield
            continue
        sampled = num.maximum(sampled, d.get_tracer('salinity'))
    assert num.all(op.max_tracer[0] >= sampled - 1e-12)


def test_the_initial_condition_is_not_sampled():
    """Pins the documented limitation, so a change to it is deliberate."""
    d = _domain()
    c0 = d.get_tracer('salinity').copy()
    op = Collect_max_quantities_operator(d)
    _run(d)
    # Some cell started at 1.0 and washed out before the first collection.
    assert num.any(c0 > op.max_tracer[0] + 1e-3), \
        'the initial condition now IS sampled -- update the operator docstring'


def test_a_dry_cell_carries_no_concentration():
    """Same rule the kernel uses when it derives c from m."""
    d = anuga.rectangular_cross_domain(6, 6, len1=60.0, len2=60.0)
    d.set_quantity('elevation', 0.0)
    d.set_quantity('stage', 0.0)               # dry everywhere
    b = anuga.Reflective_boundary(d)
    d.set_boundary({'left': b, 'right': b, 'top': b, 'bottom': b})
    d.add_tracer('salinity', initial_value=1.0)
    d.store = False
    d.set_compute_mode('legacy')
    op = Collect_max_quantities_operator(d)
    _run(d, finaltime=1.0)
    assert num.allclose(op.max_tracer[0], 0.0), \
        'a dry cell reported a concentration'


def test_two_tracers_are_tracked_separately():
    d = _domain()
    d.add_tracer('dye', initial_value=0.25)
    op = Collect_max_quantities_operator(d)
    _run(d)
    assert not num.allclose(op.max_tracer[0], op.max_tracer[1]), \
        'the two tracers share a row'
    # dye starts uniform at 0.25 and is diluted by the clean inflow, so its
    # maximum is at most the starting value and above zero away from the inlet.
    assert op.max_tracer[1].max() <= 0.25 + 1e-9
    assert op.max_tracer[1].max() > 0.2


def test_a_tracer_added_afterwards_is_refused():
    """Its array would be the wrong width and it could get no sww variable."""
    d = _domain()
    op = Collect_max_quantities_operator(d)
    d.add_tracer('late')
    with pytest.raises(RuntimeError, match='tracers changed'):
        _run(d, finaltime=1.0)


# --- and into the file ------------------------------------------------------

def test_the_maxima_reach_the_sww(tmp_path):
    d = _domain()
    op = Collect_max_quantities_operator(d, store_to_sww=True)
    assert d.quantities_to_be_stored['max_salinity'] == 4, \
        'stored as flag 4, i.e. overwritten each yield step'
    d.store = True
    d.set_name('maxtr')
    d.set_datadir(str(tmp_path))
    _run(d)

    with netCDF4.Dataset(str(tmp_path / 'maxtr.sww'), 'r') as fid:
        assert 'max_salinity_c' in fid.variables, \
            'the tracer maximum is missing from the sww'
        stored = num.array(fid.variables['max_salinity_c'][:])
    # Flag 4 is written without a time dimension: one slice, the final maximum.
    assert stored.shape == (len(d),)
    assert num.allclose(stored, op.max_tracer[0], atol=1e-6)

"""Generic passive tracers carried by the shallow water solver.

A tracer is a depth-averaged concentration `c` transported as the conserved
quantity `m = h c`, alongside stage and momentum and through the same flux
kernel. Domains carry none by default; `domain.add_tracer(name)` registers one.

The properties tested here, in the order they matter:

  reconstruction  first order at beta = 0, constants exact, monotone on a step
                  (which is what keeps `c` positive), second order on a smooth
                  monotone field
  flux            consistency (`c = 1` reproduces the stage tendency exactly),
                  linearity, upwinding, and discrete conservation
  striding        with two tracers the tracer-major `(ns, ...)` layout is
                  indexed correctly and the slots do not alias
  integration     mass conserved through a dam break, `c` positive and bounded
                  without being clamped, and actually transported
  API             `add_tracer` shapes, contiguity, cache invalidation, and the
                  errors it raises

The transport tests are parametrised over every flow algorithm, because a
kernel that does not carry tracers drops them silently -- nothing errors, the
concentration just stops moving. `test_every_transport_path_is_covered` fails
if a new path selector appears without the parametrisation growing to match.

Mode 1 only. The mode 1 / mode 2 comparison is in test_tracers_gpu.py, which
has to guard itself against the NVHPC runtime.
"""
import numpy as np
import pytest

import anuga
from anuga import Reflective_boundary, rectangular_cross_domain
from anuga.shallow_water.sw_domain_openmp_ext import (
    compute_fluxes_ext_central, extrapolate_second_order_sw)

LEN = 1000.0


def _flat_domain(nxy, beta=1.0, stage=5.0, names=('c',)):
    """A uniformly wet, flat, still domain -- deep enough that the wet/dry
    hfactor does not throttle the reconstruction."""
    d = rectangular_cross_domain(nxy, nxy, len1=LEN, len2=LEN)
    d.set_flow_algorithm('DE0')
    d.set_low_froude(0)
    d.store = False
    d.set_quantity('elevation', 0.0)
    d.set_quantity('stage', stage)
    d.set_quantity('xmomentum', 0.0)
    d.set_quantity('ymomentum', 0.0)
    d.set_boundary({t: Reflective_boundary(d) for t in d.get_boundary_tags()})
    for i, nm in enumerate(names):
        d.add_tracer(nm, beta=beta if i == 0 else None)
    return d


def _reconstruct(d, cfield, slot=0):
    """Seed m = h c, extrapolate, and return the edge values.

    `c` is derived from the conserved `m` each substep, so the field has to be
    seeded through `m` rather than written to the centroid values directly.
    """
    h = (d.quantities['stage'].centroid_values
         - d.quantities['elevation'].centroid_values)
    d.tracer_conserved_values[slot, :] = np.maximum(h, 0.0) * cfield
    d.tracer_centroid_values[slot, :] = cfield
    d.tracer_edge_values[slot, :] = np.nan      # poison: detect non-writes
    extrapolate_second_order_sw(d, update_domain_c_struct=True)
    ev = d.tracer_edge_values[slot].copy()
    assert not np.isnan(ev).any(), 'some edge values were never written'
    return ev


# ---------------------------------------------------------------- reconstruction

def test_beta_zero_is_first_order():
    """beta = 0 must leave every edge value equal to its centroid value."""
    d = _flat_domain(12, beta=0.0)
    c = 0.5 + 0.4 * np.sin(2 * np.pi * d.centroid_coordinates[:, 0] / LEN)
    ev = _reconstruct(d, c)
    # c is derived as m/h, so the round trip costs about one ulp.
    assert np.max(np.abs(ev - np.repeat(c, 3))) < 1e-15


def test_uniform_field_is_reproduced_exactly():
    """A reconstruction that cannot reproduce a constant is broken outright."""
    d = _flat_domain(12, beta=1.0)
    ev = _reconstruct(d, np.full(len(d), 0.7))
    assert np.allclose(ev, 0.7, rtol=0, atol=1e-15)


def test_step_creates_no_new_extrema_and_stays_positive():
    """Monotonicity on a step is what protects positivity of concentration."""
    d = _flat_domain(24, beta=1.0)
    xc = d.centroid_coordinates[:, 0]
    c = np.where(xc < LEN / 2, 1.0, 0.0)
    ev = _reconstruct(d, c)
    assert ev.max() - c.max() <= 1e-14, 'overshoot above the centroid maximum'
    assert c.min() - ev.min() <= 1e-14, 'undershoot below the centroid minimum'
    assert ev.min() >= 0.0, 'negative edge concentration'


def test_linear_field_is_exact_on_interior_cells():
    """A linear field is the fixed point of a second-order reconstruction.

    Interior cells only: where a cell touches the boundary,
    `surrogate_neighbours` substitutes the cell itself and degrades the
    gradient stencil. That is ANUGA's existing behaviour, and stage shows the
    same error on the same cells.
    """
    d = _flat_domain(16, beta=1.0)
    cc = d.centroid_coordinates
    def f(x, y):
        return 0.5 + 0.3 * (x / LEN) + 0.2 * (y / LEN)
    ev = _reconstruct(d, f(cc[:, 0], cc[:, 1]))
    ec = d.get_edge_midpoint_coordinates()
    err = np.abs(ev - f(ec[:, 0], ec[:, 1]))
    interior = np.repeat(d.number_of_boundaries, 3) == 0
    assert err[interior].max() < 1e-12


@pytest.mark.slow
def test_order_of_accuracy_on_a_smooth_monotone_field():
    """beta = 0 converges at first order, beta = 1 at second.

    The field is deliberately monotone. A limiter cannot be both
    unconditionally monotone and cleanly second order at a smooth extremum
    (Godunov), so it clips at crests and drops locally to first order there --
    which is the price of the monotonicity tested above, and the right trade.
    """
    def f(x, y):
        return 0.5 + 0.4 * (x / LEN) ** 2

    def rms_errors(beta, sizes):
        out = []
        for nxy in sizes:
            d = _flat_domain(nxy, beta=beta)
            cc = d.centroid_coordinates
            ev = _reconstruct(d, f(cc[:, 0], cc[:, 1]))
            ec = d.get_edge_midpoint_coordinates()
            m = np.repeat(d.number_of_boundaries, 3) == 0
            out.append(float(np.sqrt(np.mean((ev[m] - f(ec[m, 0], ec[m, 1])) ** 2))))
        return out

    sizes = (8, 16, 32)
    e1 = rms_errors(0.0, sizes)
    e2 = rms_errors(1.0, sizes)
    o1 = np.mean([np.log2(e1[i] / e1[i + 1]) for i in range(len(e1) - 1)])
    o2 = np.mean([np.log2(e2[i] / e2[i + 1]) for i in range(len(e2) - 1)])
    assert 0.7 < o1 < 1.3, 'beta=0 should be first order, got %.2f' % o1
    assert o2 > 1.95, 'beta=1 should be second order, got %.2f' % o2


# ---------------------------------------------------------------- flux kernel

# Every distinct transport path this build offers. The tracer tests below are
# parametrised over it, which is what makes "a tracer is transported" an
# enforced guarantee rather than an argued one: a flux or update kernel that
# does not carry tracers fails here instead of silently dropping them, and a
# dropped tracer is invisible -- the concentration simply stops moving.
#
# Adding a path is one entry. See test_every_transport_path_is_covered below,
# which fails when a new path selector appears and this list has not grown.
_FLOW_ALGORITHMS = ('DE0', 'DE1', 'DE2')

# Domain attributes, not yet present, that will select an alternative flux or
# update kernel. Named here so their arrival is noticed; see the test.
_FUTURE_PATH_SELECTORS = ('reconstruct_edge_bed',)


def _dam_break_domain(nxy=20, names=('c',), beta=1.0, algorithm='DE0'):
    d = rectangular_cross_domain(nxy, nxy, len1=LEN, len2=LEN)
    d.set_flow_algorithm(algorithm)
    d.set_low_froude(0)
    d.store = False
    d.set_quantity('elevation', 0.0)
    d.set_quantity('stage', lambda x, y: np.where(x < LEN / 2, 2.0, 0.5))
    d.set_quantity('xmomentum', 0.0)
    d.set_quantity('ymomentum', 0.0)
    d.set_boundary({t: Reflective_boundary(d) for t in d.get_boundary_tags()})
    for i, nm in enumerate(names):
        d.add_tracer(nm, beta=beta if i == 0 else None)
    return d


def _host_kernel_domain(*args, **kwargs):
    """A dam-break domain pinned to legacy, for the direct-kernel tests below.

    The tests that use this call the mode-1 flux kernel themselves --
    compute_fluxes_ext_central -- and then assert on the host arrays. Under
    ANUGA_DEFAULT_COMPUTE_MODE=unified that is not a meaningful thing to do:
    the domain evolves on the DEVICE, so the host centroids still hold the
    initial condition, and distribute_to_vertices_and_edges() never fills the
    host tracer EDGE values. The host kernel then computes a real stage flux
    from real stage edges and a ZERO tracer flux from zero tracer edges, and
    the identities below fail -- or worse, hold vacuously.

    They are mode-1 statements about the flux kernel's algebra, so they are
    pinned rather than made mode-aware. Mode 2 is covered by
    test_tracers_gpu.py, which compares the two modes end to end.

    Only the direct-kernel tests use this. The tests that merely evolve a
    domain and check the result keep the unpinned fixture, so they still
    exercise whichever mode the suite is run in.
    """
    d = _dam_break_domain(*args, **kwargs)
    d.set_compute_mode('legacy')
    return d


def _prime(d, finaltime=2.0):
    """Advance far enough that there is real flow through interior edges.

    This also gets the boundary conditions applied. Evaluating fluxes on a
    freshly built domain instead measures nothing useful: the boundary arrays
    are still unset, tracer leaks out through the walls, and the conservation
    check below fails with a residual of 87% rather than 4e-17.
    """
    d.evolve_to_end(finaltime=finaltime)
    return d


def _fluxes(d):
    """Take one flux evaluation with the current tracer state."""
    d.tracer_explicit_update[:] = 0.0
    # Recompute so the stage/momentum edge values match the current centroids.
    d.distribute_to_vertices_and_edges()
    d.update_boundary()
    compute_fluxes_ext_central(d, d.evolve_max_timestep,
                               update_domain_c_struct=True)
    return (d.quantities['stage'].explicit_update.copy(),
            d.tracer_explicit_update.copy())


def _seed_uniform(d, value, slot=0):
    h = (d.quantities['stage'].centroid_values
         - d.quantities['elevation'].centroid_values)
    d.tracer_conserved_values[slot, :] = np.maximum(h, 0.0) * value
    d.tracer_centroid_values[slot, :] = value
    d.tracer_boundary_values[slot, :] = value


def test_flux_consistency_c_equal_one_reproduces_the_stage_tendency():
    """With c = 1 everywhere, m = h, so dm/dt must equal dh/dt exactly."""
    d = _prime(_host_kernel_domain())
    _seed_uniform(d, 1.0)
    stage_eu, tracer_eu = _fluxes(d)
    scale = max(np.abs(stage_eu).max(), 1e-30)
    assert np.abs(tracer_eu[0] - stage_eu).max() < 1e-12 * scale


def test_flux_is_linear_in_concentration():
    d = _prime(_host_kernel_domain())
    K = 0.375
    _seed_uniform(d, K)
    stage_eu, tracer_eu = _fluxes(d)
    scale = max(np.abs(K * stage_eu).max(), 1e-30)
    assert np.abs(tracer_eu[0] - K * stage_eu).max() < 1e-12 * scale


def test_flux_is_conservative_on_a_closed_domain():
    """Reflective everywhere, so the area-weighted tendency must sum to zero."""
    d = _prime(_host_kernel_domain())
    rng = np.random.default_rng(42)
    _seed_uniform(d, 0.0)
    h = (d.quantities['stage'].centroid_values
         - d.quantities['elevation'].centroid_values)
    c = rng.random(len(d))
    d.tracer_conserved_values[0, :] = np.maximum(h, 0.0) * c
    d.tracer_centroid_values[0, :] = c
    d.tracer_boundary_values[0, :] = 0.0
    _, tracer_eu = _fluxes(d)
    total = float((tracer_eu[0] * d.areas).sum())
    scale = float((np.abs(tracer_eu[0]) * d.areas).sum())
    assert abs(total) < 1e-12 * max(scale, 1e-30)


def test_transport_is_upwinded():
    """Cold cells ahead of the front stay cold; cells at the front gain."""
    d = _prime(_host_kernel_domain())
    xc = d.centroid_coordinates[:, 0]
    h = (d.quantities['stage'].centroid_values
         - d.quantities['elevation'].centroid_values)
    c = np.where(xc < LEN / 2, 1.0, 0.0)
    d.tracer_conserved_values[0, :] = np.maximum(h, 0.0) * c
    d.tracer_centroid_values[0, :] = c
    d.tracer_boundary_values[0, :] = 0.0
    _, tracer_eu = _fluxes(d)
    far = xc > 0.75 * LEN
    assert np.abs(tracer_eu[0][far]).max() < 1e-14, 'action at a distance'
    front = (xc > LEN / 2) & (xc < 0.55 * LEN)
    assert tracer_eu[0][front].max() > 0.0, 'no tracer crossed the front'


# ---------------------------------------------------------------- two tracers

def test_two_tracers_are_indexed_and_do_not_alias():
    """The tracer-major (ns, ...) layout is indexed as s*n / s*3n."""
    d = _prime(_host_kernel_domain(names=('a', 'b')))
    _seed_uniform(d, 1.0, slot=0)
    _seed_uniform(d, 0.375, slot=1)
    stage_eu, tracer_eu = _fluxes(d)
    scale = max(np.abs(stage_eu).max(), 1e-30)
    assert np.abs(tracer_eu[0] - stage_eu).max() < 1e-12 * scale
    assert np.abs(tracer_eu[1] - 0.375 * stage_eu).max() < 1e-12 * scale
    assert not np.allclose(tracer_eu[0], tracer_eu[1]), 'the slots alias'


# ---------------------------------------------------------------- integration

@pytest.mark.parametrize('algorithm', _FLOW_ALGORITHMS)
def test_uniform_tracer_survives_a_dam_break_unchanged(algorithm):
    """Free stream: a uniform c must stay uniform however violent the flow."""
    d = _dam_break_domain(nxy=30, algorithm=algorithm)
    d.set_tracer('c', 0.4)
    d.tracer_boundary_values[0, :] = 0.4
    d.evolve_to_end(finaltime=20.0)
    h = (d.quantities['stage'].centroid_values
         - d.quantities['elevation'].centroid_values)
    wet = h > 1e-3
    c = d.tracer_conserved_values[0][wet] / h[wet]
    assert np.abs(c - 0.4).max() < 1e-9


@pytest.mark.parametrize('algorithm', _FLOW_ALGORITHMS)
def test_tracer_mass_is_conserved_positive_and_bounded(algorithm):
    d = _dam_break_domain(nxy=30, algorithm=algorithm)
    xc = d.centroid_coordinates[:, 0]
    d.set_tracer('c', np.where(xc < LEN / 4, 0.8, 0.0))
    d.tracer_boundary_values[0, :] = 0.0
    m0 = float((d.tracer_conserved_values[0] * d.areas).sum())

    cmin, cmax = [], []
    for _ in d.evolve(yieldstep=5.0, finaltime=20.0):
        h = (d.quantities['stage'].centroid_values
             - d.quantities['elevation'].centroid_values)
        wet = h > 1e-3
        if wet.any():
            c = d.tracer_conserved_values[0][wet] / h[wet]
            cmin.append(c.min())
            cmax.append(c.max())

    m1 = float((d.tracer_conserved_values[0] * d.areas).sum())
    assert abs(m1 - m0) < 1e-10 * max(abs(m0), 1.0), 'mass drifted'
    # Positivity and boundedness are EMERGENT here -- nothing clamps c.
    assert min(cmin) > -1e-12, 'concentration went negative'
    assert max(cmax) < 0.8 + 1e-9, 'concentration exceeded its initial maximum'


@pytest.mark.slow
@pytest.mark.parametrize('algorithm', _FLOW_ALGORITHMS)
def test_tracer_is_actually_transported(algorithm):
    """Guards against a suite that would pass on a tracer that never moves.

    Needs the whole upstream reservoir, and long enough to matter: over 20 s
    the centre of mass barely shifts and the test measures nothing.
    """
    d = _dam_break_domain(nxy=30, algorithm=algorithm)
    xc = d.centroid_coordinates[:, 0]
    d.set_tracer('c', np.where(xc < LEN / 2, 1.0, 0.0))
    d.tracer_boundary_values[0, :] = 0.0
    w = d.tracer_conserved_values[0] * d.areas
    x0 = float((w * xc).sum() / max(w.sum(), 1e-30))
    d.evolve_to_end(finaltime=120.0)
    w = d.tracer_conserved_values[0] * d.areas
    x1 = float((w * xc).sum() / max(w.sum(), 1e-30))
    assert x1 > x0 + 1.0, 'the tracer centre of mass did not move downstream'


def test_every_transport_path_is_covered():
    """Fail when a new transport path appears that the tests do not cover.

    A flux or update kernel that does not carry tracers drops them SILENTLY --
    the concentration simply stops moving, with no error and no failing test
    unless something exercises that path. The parametrisation above covers the
    paths that exist today; this exists so that the next one cannot arrive
    unnoticed.

    If this fails, a domain attribute selecting an alternative kernel has been
    added. That is not a defect -- the fix is to decide whether tracers are
    carried on the new path, make sure they are, and add it to _FLOW_ALGORITHMS
    or to a parametrisation of its own so the guarantee stays enforced.
    """
    d = _small_domain()
    present = [name for name in _FUTURE_PATH_SELECTORS if hasattr(d, name)]
    assert not present, (
        'Domain now has %s, which selects an alternative flux or update '
        'kernel. Check that kernel carries tracers, then extend the '
        'parametrisation in this file to cover it and remove the name from '
        '_FUTURE_PATH_SELECTORS.' % ', '.join(repr(n) for n in present))


# ---------------------------------------------------------------- add_tracer

def _small_domain(nxy=10):
    d = rectangular_cross_domain(nxy, nxy, len1=500.0, len2=500.0)
    d.set_flow_algorithm('DE0')
    d.store = False
    d.set_quantity('elevation', 0.0)
    d.set_quantity('stage', 1.0)
    d.set_boundary({t: Reflective_boundary(d) for t in d.get_boundary_tags()})
    return d


def test_a_fresh_domain_has_no_tracers():
    d = _small_domain()
    assert getattr(d, 'number_of_tracers', 0) == 0


def test_add_tracer_registers_and_shapes_the_arrays():
    d = _small_domain()
    assert d.add_tracer('mud') == 0
    assert d.number_of_tracers == 1
    n, bl = len(d), d.boundary_length
    expected = {'tracer_centroid_values': (1, n),
                'tracer_edge_values': (1, 3 * n),
                'tracer_boundary_values': (1, bl),
                'tracer_explicit_update': (1, n),
                'tracer_conserved_values': (1, n),
                'tracer_backup_values': (1, n)}
    for name, shape in expected.items():
        a = getattr(d, name)
        assert a.shape == shape, '%s has shape %s' % (name, a.shape)
        # The kernel indexes these as s*N + k, so the layout is load-bearing.
        assert a.flags['C_CONTIGUOUS'], '%s is not C-contiguous' % name
        assert a.dtype == np.float64, '%s is not float64' % name


def test_add_tracer_invalidates_the_cached_c_struct():
    """Adding a tracer after an evolve must not leave a stale struct behind."""
    d = _small_domain()
    d.add_tracer('a')
    d.evolve_to_end(finaltime=0.2)
    assert d._Domain_C_struct is not None
    d.add_tracer('b')
    assert d._Domain_C_struct is None


def test_growing_to_two_tracers_preserves_the_first():
    d = _small_domain()
    d.add_tracer('a')
    d.set_tracer('a', 0.25)
    d.tracer_explicit_update[0, :] = 7.0        # a distinctive row
    d.add_tracer('b')
    assert np.allclose(d.get_tracer('a'), 0.25), 'tracer 0 was disturbed'
    assert np.allclose(d.tracer_explicit_update[0], 7.0)
    assert np.allclose(d.tracer_explicit_update[1], 0.0), 'new row not zeroed'


def test_set_tracer_accepts_a_scalar_or_a_field():
    d = _small_domain()
    d.add_tracer('s')
    field = np.linspace(0.0, 1.0, len(d))
    d.set_tracer('s', field)
    assert np.allclose(d.get_tracer('s'), field)
    h = (d.quantities['stage'].centroid_values
         - d.quantities['elevation'].centroid_values)
    assert np.allclose(d.tracer_conserved_values[0], h * field)
    d.set_tracer('s', 0.5)
    assert np.all(d.get_tracer('s') == 0.5)


def test_the_api_rejects_mistakes():
    # Deliberately not 'x'/'y': those are quantities on every domain, and a
    # tracer may not take a quantity's name -- see
    # test_a_tracer_cannot_be_named_after_a_quantity.
    d = _small_domain()
    d.add_tracer('a', beta=1.0)
    with pytest.raises(ValueError):
        d.add_tracer('a')                       # duplicate name
    with pytest.raises(ValueError):
        d.set_tracer('a', np.zeros(len(d) + 1))  # wrong length
    with pytest.raises(ValueError):
        d.set_tracer('nope', 1.0)               # unknown name
    with pytest.raises(ValueError):
        # beta is shared by every tracer, so a conflicting value must be
        # refused rather than silently applied to all of them.
        d.add_tracer('b', beta=0.0)

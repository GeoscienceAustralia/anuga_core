"""Tracer fields survive a mesh reorder (issue #277).

`reorder()` renumbers the triangles for cache locality. Quantities are permuted
by the loop over `domain.quantities`, but tracers are deliberately not Quantity
objects -- they live in contiguous (ns, N) / (ns, 3N) blocks so the C kernel can
stride them -- so they need permuting explicitly.

The failure this guards against is silent: cell k's concentration ends up on
whichever triangle used to be at k, with no error and a plausible-looking field.
So the tests check the property that actually matters -- a reordered run gives
the same answer as an unreordered one -- rather than just that the call returns.
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
    return d


def _wedge(d):
    """A spatially varying field, so a permutation cannot hide in it."""
    x = d.centroid_coordinates[:, 0]
    y = d.centroid_coordinates[:, 1]
    return 0.25 + 0.5 * x + 0.25 * y


def _shuffled(n, seed=7):
    rng = num.random.default_rng(seed)
    return rng.permutation(n)


def test_reorder_without_tracers_still_works():
    """The tracer hook must not disturb the ordinary path.

    Kept from test_tracers_unsupported_paths.py, which is retired now that
    both paths it guarded (#277 here, #278 in distribute) are implemented.
    """
    d = _domain()
    n = len(d)
    d.reorder(num.arange(n - 1, -1, -1))       # reverse; must not raise
    assert len(d) == n


def test_reorder_with_tracers_no_longer_refuses():
    d = _domain()
    d.add_tracer('salinity', initial_value=0.02)
    d.reorder(_shuffled(len(d)))           # used to raise NotImplementedError
    assert num.allclose(d.get_tracer('salinity'), 0.02)


def test_the_concentration_follows_its_triangle():
    d = _domain()
    d.add_tracer('c')
    d.set_tracer('c', _wedge(d))
    before = d.get_tracer('c').copy()

    order = _shuffled(len(d))
    d.reorder(order)

    # new[i] == old[order[i]] -- the same convention the quantities use.
    assert num.allclose(d.get_tracer('c'), before[order])


def test_the_concentration_still_matches_its_own_centroid():
    """The strongest statement of alignment: c is a known function of (x, y)."""
    d = _domain()
    d.add_tracer('c')
    d.set_tracer('c', _wedge(d))

    d.reorder(_shuffled(len(d)))

    # centroid_coordinates has moved with the mesh; recompute from it.
    assert num.allclose(d.get_tracer('c'), _wedge(d))


def test_the_conserved_variable_moves_with_the_concentration():
    d = _domain()
    d.add_tracer('c')
    d.set_tracer('c', _wedge(d))
    d.reorder(_shuffled(len(d)))

    h = (d.quantities['stage'].centroid_values
         - d.quantities['elevation'].centroid_values)
    assert num.allclose(d.tracer_conserved_values[0], h * d.get_tracer('c')), \
        'm = h*c is inconsistent after reorder'


def test_two_tracers_do_not_get_swapped():
    d = _domain()
    d.add_tracer('a')
    d.add_tracer('b')
    d.set_tracer('a', _wedge(d))
    d.set_tracer('b', 1.0 - _wedge(d))

    d.reorder(_shuffled(len(d)))

    assert num.allclose(d.get_tracer('a'), _wedge(d))
    assert num.allclose(d.get_tracer('b'), 1.0 - _wedge(d))


def test_the_arrays_stay_contiguous_for_the_kernel():
    """The C kernel indexes these as centroid[s*N + k]; a strided view breaks it."""
    d = _domain()
    d.add_tracer('c')
    d.set_tracer('c', _wedge(d))
    d.reorder(_shuffled(len(d)))

    for attr in d._TRACER_ARRAYS:
        arr = getattr(d, attr)
        assert arr.flags['C_CONTIGUOUS'], '%s is no longer C-contiguous' % attr
        assert arr.dtype == num.float64, '%s changed dtype' % attr


# --------------------------------------------------------------------------
# Boundary values: indexed by boundary edge, which reorder() also renumbers.
# --------------------------------------------------------------------------

def _tagged_domain(n=6):
    d = anuga.rectangular_cross_domain(n, n)
    d.set_quantity('elevation', 0.0)
    d.set_quantity('stage', 1.0)
    Br = anuga.Reflective_boundary(d)
    d.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})
    return d


def _boundary_edge_midpoints(d, tag):
    """(x, y) of each boundary edge of `tag`, in that tag's own order."""
    idx = num.asarray(d.tag_boundary_cells[tag], dtype=int)
    cells = d.boundary_cells[idx]
    edges = d.boundary_edges[idx]
    return num.column_stack((d.edge_coordinates[3 * cells + edges, 0],
                             d.edge_coordinates[3 * cells + edges, 1]))


def test_boundary_values_follow_their_edge():
    """The subtle half: boundary indices are renumbered by reorder() too.

    Quantities do not notice, because update_boundary refills their boundary
    values every timestep from the Boundary objects. tracer_boundary_values is
    written once by set_tracer_boundary and then only read, so it has to be
    carried across explicitly.
    """
    d = _tagged_domain()
    d.add_tracer('c')

    # A distinct value per edge, keyed to position so a mismatch is visible.
    mid = _boundary_edge_midpoints(d, 'left')
    d.set_tracer_boundary('c', 'left', 1.0 + mid[:, 1])
    d.set_tracer_boundary('c', 'right', 5.0)

    d.reorder(_shuffled(len(d)))

    new_mid = _boundary_edge_midpoints(d, 'left')
    assert num.allclose(d.get_tracer_boundary('c', 'left'), 1.0 + new_mid[:, 1]), \
        'boundary concentrations no longer sit on the edges they were set for'
    assert num.allclose(d.get_tracer_boundary('c', 'right'), 5.0)


def test_untouched_boundaries_stay_zero():
    d = _tagged_domain()
    d.add_tracer('c')
    d.set_tracer_boundary('c', 'left', 3.0)
    d.reorder(_shuffled(len(d)))

    assert num.allclose(d.get_tracer_boundary('c', 'left'), 3.0)
    for tag in ('right', 'top', 'bottom'):
        assert num.allclose(d.get_tracer_boundary('c', tag), 0.0), \
            'tag %r picked up values from another boundary' % tag


# --------------------------------------------------------------------------
# The property that actually matters: same answer, reordered or not.
# --------------------------------------------------------------------------

def _run(reorder, finaltime=2.0, n=8, seed=3):
    d = anuga.rectangular_cross_domain(n, n)
    d.set_quantity('elevation', lambda x, y: -0.2 * x)
    d.set_quantity('stage', lambda x, y: num.maximum(0.4 - 0.2 * x, 0.05))
    Br = anuga.Reflective_boundary(d)
    Bt = anuga.Transmissive_boundary(d)
    d.set_boundary({'left': Br, 'top': Br, 'bottom': Br, 'right': Bt})
    d.add_tracer('c')
    d.set_tracer('c', _wedge(d))
    d.set_tracer_boundary('c', 'left', 0.9)

    if reorder:
        d.reorder(_shuffled(len(d), seed=seed))

    for _ in d.evolve(yieldstep=0.5, finaltime=finaltime):
        pass

    # Sort by centroid so the two runs can be compared cell for cell.
    coords = d.centroid_coordinates
    idx = num.lexsort((coords[:, 1], coords[:, 0]))
    return d, d.get_tracer('c')[idx], coords[idx]


def test_a_reordered_run_gives_the_same_tracer_field():
    ref_d, ref_c, ref_xy = _run(reorder=False)
    new_d, new_c, new_xy = _run(reorder=True)

    assert num.allclose(ref_xy, new_xy), 'meshes are not comparable'
    # Not a trivially uniform field, or the comparison proves nothing.
    assert num.ptp(ref_c) > 1e-3, 'the tracer field is too flat to discriminate'
    assert num.allclose(ref_c, new_c, atol=1e-12), \
        'reordering changed the tracer solution'


def test_conservation_still_balances_after_a_reorder():
    d, _, _ = _run(reorder=True)
    change, flux, disc = d.check_tracer_conservation('c')
    assert abs(change) > 0.0
    assert abs(disc) < 1e-10 * abs(change), \
        'the boundary accounting lost track of the mesh: %g vs %g' % (change, flux)


@pytest.mark.parametrize('method', ['hilbert', 'morton'])
def test_reorder_domain_helper_handles_tracers(method):
    """The path that actually calls reorder() in practice."""
    from anuga.parallel.partitioning import reorder_domain

    d = _domain()
    d.add_tracer('c')
    d.set_tracer('c', _wedge(d))
    reorder_domain(d, method=method)

    assert num.allclose(d.get_tracer('c'), _wedge(d))

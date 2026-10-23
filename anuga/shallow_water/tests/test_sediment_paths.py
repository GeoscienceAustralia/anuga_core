"""Guard: the sediment kernels are reached on every path.

The sediment physics runs as a fractional step dispatched by Sediment_transport_operator
to three kernels -- suspended exchange, bedload, angle-of-repose relaxation --
through a chain of conditions: the operator must be registered, the mode 1 /
mode 2 branch taken, and each kernel's enable flag read.

Every link in that chain fails SILENTLY. A sediment class that is registered
but never exchanged simply behaves as an inert tracer: no error, no warning,
the bed does not move. The other test_sediment_* modules cover the physics;
this covers whether the physics is reached at all.
"""
import numpy as np
import pytest

from anuga import Reflective_boundary, rectangular_cross_domain

POROSITY = 0.3
FLOW_ALGORITHMS = ('DE0', 'DE1', 'DE2')

# Domain attributes, not present yet, that would select an alternative flux or
# update kernel. Named so their arrival is noticed rather than discovered.
# anuga-community#241 proposes 'reconstruct_edge_bed' for a scatter flux
# kernel; alternative kernels raise the same question for the sediment
# fractional step as for tracers, and the failure is equally silent.
FUTURE_PATH_SELECTORS = ('reconstruct_edge_bed',)


def build(algorithm='DE0', bedload=False, repose=None, sediment=True):
    d = rectangular_cross_domain(30, 8, len1=60.0, len2=16.0)
    d.set_flow_algorithm(algorithm)
    d.set_low_froude(0)
    d.store = False
    x = d.centroid_coordinates[:, 0]
    d.set_quantity('elevation', -0.01 * x, location='centroids')
    d.set_quantity('stage', -0.01 * x + 0.6, location='centroids')
    d.set_quantity('xmomentum', 1.2)
    d.set_quantity('friction', 0.03)
    d.set_boundary({t: Reflective_boundary(d) for t in d.get_boundary_tags()})
    if sediment:
        d.set_sediment_parameters(porosity=POROSITY)
        if bedload:
            d.set_bedload('wong_parker_eq24')
        d.add_sediment_fraction(name='sand', diameter=3.0e-4)
        if repose is not None:
            d.set_angle_of_repose(repose, max_sweeps=50)
    return d


def evolve(d, finaltime=30.0):
    z0 = d.quantities['elevation'].centroid_values.copy()
    d.evolve_to_end(finaltime=finaltime)
    z1 = d.quantities['elevation'].centroid_values
    ns = getattr(d, 'n_sediment_classes', 0)
    m = sum(float((d.tracer_conserved_values[s] * d.areas).sum())
            for s in range(ns))
    bed = float(((1.0 - POROSITY) * (z1 - z0) * d.areas).sum())
    return dict(z0=z0, z1=z1, scour=float((z0 - z1).max()),
                suspended=m, budget=m + bed)


def cone_domain(repose=None, algorithm='DE0'):
    """A dry cone far steeper than any repose angle. Relaxation needs a bed
    steep enough to relax; on the gently sloping channel above it would
    correctly do nothing, and the comparison would prove only that."""
    d = rectangular_cross_domain(40, 12, len1=60.0, len2=16.0)
    d.set_flow_algorithm(algorithm)
    d.set_low_froude(0)
    d.store = False
    x, y = d.centroid_coordinates[:, 0], d.centroid_coordinates[:, 1]
    r = np.sqrt((x - 30.0) ** 2 + (y - 8.0) ** 2)
    d.set_quantity('elevation', np.where(r < 8.0, 6.0 * (1.0 - r / 8.0), 0.0),
                   location='centroids')
    d.set_quantity('stage', -1.0)
    d.set_quantity('friction', 0.03)
    d.set_boundary({t: Reflective_boundary(d) for t in d.get_boundary_tags()})
    d.set_sediment_parameters(porosity=POROSITY)
    d.add_sediment_fraction(name='sand', diameter=3.0e-4)
    if repose is not None:
        d.set_angle_of_repose(repose, max_sweeps=400)
    return d


def test_the_suspended_exchange_is_reached():
    r = evolve(build())
    assert r['scour'] > 1e-4, 'the bed did not move'
    assert r['suspended'] > 0.0, 'no sediment entered suspension'
    # The closed budget shows it is the exchange doing the work, rather than
    # something else moving the bed.
    assert abs(r['budget']) < 1e-10 * max(r['suspended'], 1.0)


def test_the_bedload_kernel_is_reached():
    """A dispatch that silently skipped bedload would leave these identical."""
    off = evolve(build())
    on = evolve(build(bedload=True))
    assert not np.allclose(off['z1'], on['z1'], rtol=0, atol=1e-14)


def test_the_repose_kernel_is_reached_and_conserves():
    off = cone_domain()
    off.evolve_to_end(finaltime=1.0)
    on = cone_domain(repose=30.0)
    on.evolve_to_end(finaltime=1.0)
    z_off = off.quantities['elevation'].centroid_values
    z_on = on.quantities['elevation'].centroid_values
    assert not np.allclose(z_off, z_on, rtol=0, atol=1e-14)
    moved = float((np.abs(z_on - z_off) * on.areas).sum())
    net = float(((z_on - z_off) * on.areas).sum())
    assert abs(net) < 1e-9 * max(moved, 1.0), 'relaxation destroyed material'


@pytest.mark.parametrize('algorithm', FLOW_ALGORITHMS)
def test_sediment_runs_under_every_flow_algorithm(algorithm):
    """Every other sediment module uses DE0. A kernel reached only from the
    DE0 path would pass all of them and silently do nothing elsewhere."""
    r = evolve(build(algorithm=algorithm, bedload=True, repose=33.0))
    assert r['scour'] > 1e-4
    assert r['suspended'] > 0.0
    assert abs(r['budget']) < 1e-9 * max(r['suspended'], 1.0)


def test_the_flow_algorithms_agree_with_each_other():
    """Not identical -- the algorithms differ -- but the same problem, so a
    wildly different answer would mean sediment saw a different flow."""
    ref = evolve(build(algorithm='DE0', bedload=True, repose=33.0))
    for algorithm in FLOW_ALGORITHMS[1:]:
        r = evolve(build(algorithm=algorithm, bedload=True, repose=33.0))
        rel = abs(r['suspended'] - ref['suspended']) / max(ref['suspended'], 1e-30)
        assert rel < 0.05, '%s differs from DE0 by %.1f%%' % (algorithm, 100 * rel)


def test_the_operator_is_registered_automatically():
    names = [type(o).__name__ for o in build().fractional_step_operators]
    assert 'Sediment_transport_operator' in names, names


def test_a_domain_without_sediment_does_not_get_the_operator():
    names = [type(o).__name__
             for o in build(sediment=False).fractional_step_operators]
    assert 'Sediment_transport_operator' not in names


def test_the_sediment_operator_does_not_force_cpu_only_fractional_steps():
    """A host-writing fractional operator reactivates the sync bracket and
    drops the whole run onto the host path -- a large, silent loss."""
    assert not build()._has_cpu_only_fractional_operators()


def test_no_uncovered_transport_path_has_appeared():
    """Fail when a new path selector arrives, so it cannot arrive unnoticed.

    Not a defect when it fires: the fix is to check the new kernel still
    reaches the sediment fractional step, then extend the parametrisation above
    and remove the name from FUTURE_PATH_SELECTORS.
    """
    d = build(sediment=False)
    present = [n for n in FUTURE_PATH_SELECTORS if hasattr(d, n)]
    assert not present, (
        'Domain now has %s, which selects an alternative flux or update '
        'kernel. Check the sediment fractional step is still reached on that '
        'path -- a kernel that does not call it leaves sediment silently '
        'inert -- then cover it here.' % ', '.join(repr(n) for n in present))

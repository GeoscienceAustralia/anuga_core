"""Vegetation drag, spec 8: Baptist et al. (2007) as a friction closure.

The drag of a stem field enters as a friction slope g q |q| / (Cv^2 h^2) with
the Baptist vegetated Chezy coefficient, semi-implicit alongside Manning. The
checks: a uniform flow down a vegetated slope reaches the Chezy normal state
u = Cv sqrt(h S) (with Cv as the closed form), emergent and submerged stems
give the two branches of Cv, unvegetated cells are untouched, and both compute
modes agree.
"""
import numpy as np
import pytest

import anuga
from anuga import Dirichlet_boundary, Reflective_boundary, rectangular_cross_domain

G = 9.8


def chezy_baptist(h, m, D, hv, Cd=1.68, Cb=65.0):
    Cv = (Cb ** -2 + Cd * m * D * min(h, hv) / (2.0 * G)) ** -0.5
    if h > hv:
        Cv += np.sqrt(G) / 0.4 * np.log(h / hv)
    return Cv


def vegetated_channel(m, D, hv, slope=1.0e-3, depth=1.0, mode='legacy', n_manning=0.0):
    """Uniform flow down a plane held by Dirichlet boundaries at both ends."""
    L, W = 200.0, 10.0
    d = rectangular_cross_domain(40, 2, len1=L, len2=W)
    d.set_flow_algorithm('DE1')
    d.set_compute_mode(mode)
    d.store = False
    d.set_quantity('elevation', lambda x, y: slope * (L - x))
    d.set_quantity('friction', n_manning)
    u = chezy_baptist(depth, m, D, hv) * np.sqrt(depth * slope) if m > 0 else 0.5
    d.set_quantity('stage', lambda x, y: slope * (L - x) + depth)
    d.set_quantity('xmomentum', u * depth)
    Bin = Dirichlet_boundary([slope * L + depth, u * depth, 0.0])
    Bout = Dirichlet_boundary([depth, u * depth, 0.0])
    Br = Reflective_boundary(d)
    d.set_boundary({'left': Bin, 'right': Bout, 'top': Br, 'bottom': Br})
    d.set_vegetation_drag(density=m, diameter=D, height=hv)
    return d, u


def normal_velocity_error(d, u):
    x = d.centroid_coordinates[:, 0]
    interior = (x > 40.0) & (x < 160.0)
    h = (d.quantities['stage'].centroid_values - d.quantities['elevation'].centroid_values)[interior]
    uh = d.quantities['xmomentum'].centroid_values[interior]
    return np.abs(uh / h - u).max() / u, np.abs(h - 1.0).max()


@pytest.mark.parametrize('m,D,hv', [(120.0, 0.01, 10.0),    # emergent
                                    (200.0, 0.01, 0.5)])    # submerged
def test_uniform_flow_reaches_the_baptist_normal_state(m, D, hv):
    """Normal flow on a vegetated plane: u = Cv sqrt(h S). Started at the
    normal state, the flow must stay there: the drag balances gravity."""
    d, u = vegetated_channel(m, D, hv)
    d.evolve_to_end(finaltime=600.0)
    eu, eh = normal_velocity_error(d, u)
    assert eu < 2.0e-3, 'velocity off the Chezy normal state by %.2e' % eu
    assert eh < 2.0e-3, 'depth off the normal depth by %.2e m' % eh


def test_without_stems_nothing_changes():
    """m = 0 everywhere: the kernel must be a no-op, so a uniform flow on a
    flat frictionless bed is unchanged to roundoff."""
    d, u = vegetated_channel(0.0, 0.01, 1.0, slope=0.0)
    uh0 = d.quantities['xmomentum'].centroid_values.copy()
    d.evolve_to_end(finaltime=60.0)
    assert np.abs(d.quantities['xmomentum'].centroid_values - uh0).max() < 1e-9


def test_the_drag_backs_the_flow_up():
    """Stems on the downstream half only, outlet stage held: the vegetated
    reach needs more surface slope to pass the discharge, so the water backs
    up and the upstream depth rises above the bare run's."""
    depths = []
    for m in (0.0, 200.0):
        d, _ = vegetated_channel(0.0, 0.01, 1.0)
        x = d.centroid_coordinates[:, 0]
        d.set_vegetation_drag(density=np.where(x > 100.0, m, 0.0), diameter=0.01, height=10.0)
        d.evolve_to_end(finaltime=600.0)
        h = d.quantities['stage'].centroid_values - d.quantities['elevation'].centroid_values
        depths.append(h[(x > 20) & (x < 80)].mean())
    assert depths[1] > depths[0] + 0.05


def test_unknown_formulation_is_rejected():
    d, _ = vegetated_channel(0.0, 0.01, 1.0)
    with pytest.raises(ValueError):
        d.set_vegetation_drag(1.0, 0.01, 1.0, formulation='nepf')
    with pytest.raises(ValueError):
        d.set_vegetation_drag(density=1.0)


def test_off_removes_the_drag():
    d, u = vegetated_channel(200.0, 0.01, 10.0)
    d.set_vegetation_drag(formulation='off')
    assert d.vegetation_mode == 0
    uh0 = d.quantities['xmomentum'].centroid_values.copy()
    d.evolve_to_end(finaltime=30.0)
    # frictionless now, and the state is the vegetated normal state, so the
    # flow accelerates: the momentum must have grown in the interior
    x = d.centroid_coordinates[:, 0]
    mid = (x > 60) & (x < 140)
    assert d.quantities['xmomentum'].centroid_values[mid].mean() > uh0[mid].mean()


def test_the_two_compute_modes_agree():
    """Same run in mode 1 and mode 2. On a CPU build mode 2 executes the
    same kernel on the host, on a GPU build it offloads; either way the
    two must agree to roundoff."""
    results = []
    for mode in ('legacy', 'unified'):
        d, _ = vegetated_channel(200.0, 0.01, 0.5, mode=mode)
        d.evolve_to_end(finaltime=120.0)
        results.append(d.quantities['xmomentum'].centroid_values.copy())
    assert np.abs(results[0] - results[1]).max() < 1e-10

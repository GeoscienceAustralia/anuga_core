"""A bed change under a thin film keeps the velocity, not the momentum.

Every sediment kernel moves the bed with the stage held, so the depth changes
by -dz. If the momentum is left as it was, a deposit close to the depth of a
thin film leaves the film nearly gone but carrying its full momentum, and the
velocity uh/h explodes. On the Zaragoza pier validation case under DE0 one
bedload step put 0.088 mm of sand under a 0.088 mm film and the speed went
from 1.8 to 972 m/s, which collapsed the timestep (DE0's protect threshold,
1e-12, is far below the depth that was left).

The bedload, repose and suspended kernels share one bed update
(core_shift_bed). It is driven here through the repose kernel, which moves
the bed without the flow and so isolates it.
"""
import numpy as np

from anuga.shallow_water.sw_domain_openmp_ext import apply_repose
from anuga.shallow_water.tests.test_sediment_repose import ANGLE, build

U0, V0 = 0.5, -0.2


def one_sweep(film):
    """The bed change and the state after one repose sweep under a film of
    the given depths, flowing at (U0, V0)."""
    d = build(angle=ANGLE, max_sweeps=1)
    z0 = d.quantities['elevation'].centroid_values.copy()
    d.set_quantity('stage', z0 + film, location='centroids')
    d.set_quantity('xmomentum', U0 * film, location='centroids')
    d.set_quantity('ymomentum', V0 * film, location='centroids')
    apply_repose(d, update_domain_c_struct=True)
    q = d.quantities
    dz = q['elevation'].centroid_values - z0
    h = q['stage'].centroid_values - q['elevation'].centroid_values
    return (dz, h, q['xmomentum'].centroid_values.copy(),
            q['ymomentum'].centroid_values.copy())


def test_a_deposit_under_a_thin_film_keeps_its_velocity():
    # The same cone and sweep give the same bed change whatever the water, so
    # a first pass finds dz and a second sets each cell's film from it: twice
    # the deposit (the film survives, at half its depth), half the deposit
    # (the bed rises through the stage), and a deep film where the bed falls.
    dz, _, _, _ = one_sweep(np.full_like(build().areas, 1.0))
    rising = dz > 0.0
    assert rising.sum() > 20 and (dz < 0.0).sum() > 20

    survives = rising & (np.arange(dz.size) % 2 == 0)
    drowned = rising & ~survives
    film = np.where(survives, 2.0 * dz, np.where(drowned, 0.5 * dz, 1.0))

    dz2, h, uh, vh = one_sweep(film)
    assert np.array_equal(dz2, dz)

    # A surviving film halves its depth and keeps its speed. Before the fix
    # the momentum stayed, so the velocity doubled here.
    assert np.allclose(h[survives], dz[survives], rtol=1e-9)
    assert np.allclose(uh[survives] / h[survives], U0, rtol=1e-9)
    assert np.allclose(vh[survives] / h[survives], V0, rtol=1e-9)

    # A bed that rises through the stage leaves no water to carry momentum.
    assert np.all(uh[drowned] == 0.0) and np.all(vh[drowned] == 0.0)

    # Where the bed falls the momentum is left alone: keeping the velocity
    # would multiply the momentum of a nearly dry cell without bound.
    falling = dz < 0.0
    assert np.array_equal(uh[falling], U0 * film[falling])
    assert np.array_equal(vh[falling], V0 * film[falling])

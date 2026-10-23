"""Bed shear closures [T-1], [T-7] and [T-7e] (PHYSICS_SPEC 3.1, 3.4).

    [T-1]  tau_b = rho f_c |v|^2     quadratic drag (default)
    [T-7]  tau_b = rho g h S         depth-slope

[T-7] is the steady uniform flow approximation, kept for reproducing published
anugaSed results. [T-1] is the default because normal-flow equilibrium fails in
exactly the unsteady floods this targets.

The mode 1 / mode 2 comparison for these closures is in test_sediment_gpu.py.
"""
import numpy as np
import pytest

from anuga import Reflective_boundary, rectangular_cross_domain

LEN = 200.0


def channel(slope=0.01, depth=1.0, n_manning=0.03, dt=1.0):
    d = rectangular_cross_domain(20, 10, len1=LEN, len2=LEN / 2)
    d.set_flow_algorithm('DE0')
    d.set_low_froude(0)
    d.store = False
    d.set_quantity('elevation', lambda x, y: -slope * x)
    d.set_quantity('stage', lambda x, y: -slope * x + depth)
    d.set_quantity('friction', n_manning)
    d.set_boundary({t: Reflective_boundary(d) for t in d.get_boundary_tags()})
    d.evolve_max_timestep = dt
    return d


def test_quadratic_drag_is_the_default():
    assert channel().sediment_shear_closure == 0


def test_depth_slope_can_be_selected():
    d = channel()
    d.set_shear_closure('depth_slope')
    assert d.sediment_shear_closure == 1


def test_energy_slope_can_be_selected():
    d = channel()
    d.set_shear_closure('energy_slope')
    assert d.sediment_shear_closure == 2


def test_energy_slope_differs_from_depth_slope_when_the_surface_is_not_parallel():
    """[T-7e] takes S from the free surface, [T-7] from the bed.

    On a bed that is FLAT while the free surface is sloping, [T-7] sees no
    slope at all and predicts no shear, while [T-7e] sees the surface slope
    that is actually driving the flow. If the two ever agree here, the energy
    closure is reading the bed.
    """
    import numpy as np
    from anuga import Dirichlet_boundary

    def run(closure):
        d = rectangular_cross_domain(20, 10, 100.0, 50.0)
        d.set_flow_algorithm('DE1')
        d.set_quantity('elevation', 0.0)              # FLAT bed: grad z == 0
        d.set_quantity('friction', 0.03)
        d.set_quantity('stage', lambda x, y: 2.0 - 0.01 * x)   # sloping surface
        d.set_boundary({'left': Dirichlet_boundary([2.0, 1.0, 0.0]),
                        'right': Dirichlet_boundary([1.0, 0.0, 0.0]),
                        'top': Reflective_boundary(d),
                        'bottom': Reflective_boundary(d)})
        d.add_sediment_fraction('sand', diameter=2.0e-4)
        d.set_shear_closure(closure)
        d.set_bed_material('cohesive', tau_crit=1e-9, K_e=1.0e-5)
        d.set_deposition(law='threshold', tau_d=0.0)
        d.set_datadir('.')
        d.set_name('t7e_' + closure)
        d.set_quantities_to_be_stored(None)
        z0 = d.quantities['elevation'].centroid_values.copy()
        for t in d.evolve(yieldstep=5.0, finaltime=10.0):
            pass
        return np.abs(d.quantities['elevation'].centroid_values - z0).max()

    bed = run('depth_slope')
    energy = run('energy_slope')
    assert energy > bed, (
        'energy_slope (%g) did not exceed depth_slope (%g) on a flat bed with '
        'a sloping surface -- it is reading the wrong gradient' % (energy, bed))


def test_an_unknown_closure_is_rejected():
    with pytest.raises(ValueError):
        channel().set_shear_closure('depth-slope')


def cell_slopes(d, values):
    """The slope the kernel uses for [T-7]/[T-7e]: the magnitude of the
    least-squares gradient of the centroid `values` over each cell and its
    neighbours, one-sided at boundaries. Mirrored here from
    core_tau_b_over_rho so the tests check the kernel rather than restate it.
    """
    cc = d.centroid_coordinates
    nb = d.neighbours
    S = np.zeros(len(values))
    for k in range(len(values)):
        pts = [n for n in nb[k] if n >= 0]
        dx = cc[pts, 0] - cc[k, 0]
        dy = cc[pts, 1] - cc[k, 1]
        df = values[pts] - values[k]
        sxx, sxy, syy = (dx * dx).sum(), (dx * dy).sum(), (dy * dy).sum()
        sxf, syf = (dx * df).sum(), (dy * df).sum()
        det = sxx * syy - sxy * sxy
        tr = sxx + syy
        if det > 1e-12 * tr * tr:
            gx = (syy * sxf - sxy * syf) / det
            gy = (sxx * syf - sxy * sxf) / det
        elif tr > 0:
            gx, gy = sxf / tr, syf / tr
        else:
            gx = gy = 0.0
        S[k] = np.hypot(gx, gy)
    return S


@pytest.mark.parametrize('imposed', [0.01, 0.05])
def test_the_centroid_bed_slope_matches_the_imposed_one(imposed):
    """[T-7] takes grad z from the bed centroids of the cell and its
    neighbours. On a plane that must give the imposed slope in EVERY cell,
    the ones along the walls and in the corners included."""
    d = channel(slope=imposed)
    S = cell_slopes(d, d.quantities['elevation'].centroid_values)
    assert np.allclose(S, imposed, rtol=1e-9, atol=1e-12)


@pytest.mark.parametrize('bed_evolution', [False, True])
def test_every_cell_erodes_at_its_own_rate_under_depth_slope(bed_evolution):
    """Two regressions in one. Still water on a uniform slope: every cell
    must entrain at the rate its own depth and the one slope give.

    Fixed bed: the slope for [T-7] used to be read from the bed EDGE values,
    which the DE extrapolation rebuilds every step as limited stage minus
    limited height. Along a reflective wall the ghost mirrors the interior,
    so the reconstruction saw no gradient and the wall cells never eroded
    while their neighbours did (the centroid gradient at a wall is
    determined by the internal neighbours).

    Evolving bed: the centroid gradient reads the neighbours, and the source
    loop lowers the bed. Taken inside that loop, a cell saw neighbours a
    preceding iteration had already lowered by ~1e-3 m against a 3e-4 m
    rise between centroids -- a slope tens of percent off, loop-order
    dependent, and different between mode 1 and mode 2. The slope is now a
    pass of its own, before anything moves.
    """
    from anuga import rectangular_cross_domain
    S, R, dia, tau_c, gamma0 = 0.001, 1.65, 5.0e-4, 0.04, 0.0024
    d = rectangular_cross_domain(20, 4, len1=20.0, len2=4.0)
    d.set_flow_algorithm('DE1')
    d.store = False
    d.set_quantity('elevation', lambda x, y: -S * x)
    d.set_quantity('friction', 0.0)
    d.set_quantity('stage', 1.0)
    d.set_boundary({t: Reflective_boundary(d) for t in d.get_boundary_tags()})
    d.initialize_sediment_operator(bed_evolution=bed_evolution)
    d.set_shear_closure('depth_slope')
    d.set_deposition(law='d_star', near_bed='constant')
    d.add_sediment_fraction('sand', diameter=dia, d_star=1.0, rho_s=2650.0)
    v_s = float(d.sediment_settling_velocity[0])
    h = (d.quantities['stage'].centroid_values
         - d.quantities['elevation'].centroid_values)

    T = 0.2                       # << h / v_s, so deposition is negligible
    d.evolve_to_end(finaltime=T)
    c = d.get_tracer('sand')
    if bed_evolution:
        # the bed must have lowered, uniformly (E varies 2% across the tank)
        z = d.quantities['elevation'].centroid_values
        dz = z - (-S * d.centroid_coordinates[:, 0])
        assert dz.max() < -1e-4, dz.max()
        assert np.allclose(dz, dz.mean(), rtol=0.05)

    X = h * S / (R * dia) / tau_c - 1.0
    E = v_s * 0.65 * gamma0 * X / (1.0 + gamma0 * X)
    # first-order in T/(h/v_s): m = E T (1 - T v_s/(2h) + ...)
    expected = E * T * (1.0 - 0.5 * T * v_s / h) / h
    ratio = c / expected
    xc = d.centroid_coordinates[:, 0]
    wall = (xc < 0.5) | (xc > 19.5)
    assert wall.any()
    assert np.allclose(ratio[wall], 1.0, atol=5e-3), (
        'wall cells entrain at %s of the expected rate' % ratio[wall])
    assert np.allclose(ratio, 1.0, atol=5e-3)


def test_the_energy_slope_substituted_into_T7_reproduces_T1():
    """Spec 3.4: S_f = f_c |v|^2 / (g h), so rho g h S_f is rho f_c |v|^2.

    This identity is the whole argument for preferring [T-1], so it is checked
    numerically rather than asserted in a comment.
    """
    f_c, vel2, h, g = 0.00883, 4.0, 1.5, 9.81
    S_f = f_c * vel2 / (g * h)
    assert abs(g * h * S_f - f_c * vel2) < 1e-15


def test_the_two_closures_give_materially_different_erosion():
    """PHYSICS_SPEC divergence D1: they are not interchangeable."""
    means = {}
    for closure in ('quadratic_drag', 'depth_slope'):
        d = channel()
        d.set_shear_closure(closure)
        d.add_sediment_fraction(name='sand', diameter=1e-4, initial_concentration=0.0)
        d.evolve_to_end(finaltime=30.0)
        means[closure] = float(d.get_tracer('sand').mean())

    assert all(np.isfinite(v) and v >= 0.0 for v in means.values()), means
    assert (abs(means['quadratic_drag'] - means['depth_slope'])
            > 0.01 * max(means.values())), (
        'the closures should differ materially, got %r' % means)


def _still_tank(slope=0.001, bed_evolution=True, **closure_kwargs):
    from anuga import rectangular_cross_domain
    d = rectangular_cross_domain(20, 4, len1=20.0, len2=4.0)
    d.set_flow_algorithm('DE1')
    d.store = False
    d.set_quantity('elevation', lambda x, y: -slope * x)
    d.set_quantity('friction', 0.0)
    d.set_quantity('stage', 1.0)
    d.set_boundary({t: Reflective_boundary(d) for t in d.get_boundary_tags()})
    d.initialize_sediment_operator(bed_evolution=bed_evolution)
    d.set_shear_closure('depth_slope', **closure_kwargs)
    d.set_deposition(law='d_star', near_bed='constant')
    d.add_sediment_fraction('sand', diameter=5.0e-4, d_star=1.0)
    d.add_sediment_fraction('control', diameter=5.0e-4, d_star=1.0, tau_c_star=10.0)
    return d


def test_bed_slope_magnitude_matches_the_kernel_mirror():
    d = channel(slope=0.02)
    S = d.bed_slope_magnitude()
    assert np.allclose(S, 0.02, rtol=1e-9, atol=1e-12)
    assert np.allclose(S, cell_slopes(d, d.quantities['elevation'].centroid_values))


def test_max_slope_caps_the_stress():
    """Still water on a 1e-2 slope with max_slope = 1e-3 must entrain at
    the 1e-3 rate: the cap, not the bed, sets the stress."""
    S, R, dia, tau_c, gamma0 = 0.01, 1.65, 5.0e-4, 0.04, 0.0024

    def rate(slope, h):
        X = h * slope / (R * dia) / tau_c - 1.0
        return 0.65 * gamma0 * X / (1.0 + gamma0 * X)

    d = _still_tank(slope=S, bed_evolution=False, max_slope=0.001)
    v_s = float(d.sediment_settling_velocity[0])
    h = (d.quantities['stage'].centroid_values - d.quantities['elevation'].centroid_values)
    T = 0.2
    d.evolve_to_end(finaltime=T)
    c = d.get_tracer('sand')
    expected = v_s * rate(0.001, h) * T * (1.0 - 0.5 * T * v_s / h) / h
    assert np.allclose(c / expected, 1.0, atol=5e-3), (c / expected).min()
    # and nowhere near the uncapped rate
    assert (c < 0.5 * v_s * rate(S, h) * T / h).all()


def test_frozen_slope_stops_the_feedback():
    """Bed evolution on. Unfrozen, the closure runs away (the bed scours by
    the order of a metre in a minute and even the tau_c* = 10 control is
    entrained). Frozen at the setup slope, the control stays clean, the
    kernel's slope stays at the setup value, and the bed lowers gently."""
    d = _still_tank(freeze_slope=True)
    assert d.sediment_slope_frozen == 1
    assert np.allclose(d.sediment_slope_work, 0.001)
    z0 = d.quantities['elevation'].centroid_values.copy()
    d.evolve_to_end(finaltime=60.0)
    assert np.allclose(d.sediment_slope_work, 0.001)       # untouched by the kernel
    dz = d.quantities['elevation'].centroid_values - z0
    assert dz.max() < 0.0 and dz.min() > -0.2, (dz.min(), dz.max())
    assert d.get_tracer('control').max() == 0.0
    assert 0.0 < d.get_tracer('sand').max() < 0.1


def test_freeze_slope_is_taken_when_the_grain_is_registered():
    """set_shear_closure before add_sediment_fraction: the array does not
    exist yet, so the slope is recorded when the first grain allocates it."""
    from anuga import rectangular_cross_domain
    d = rectangular_cross_domain(10, 4, len1=10.0, len2=4.0)
    d.set_quantity('elevation', lambda x, y: -0.003 * x)
    d.set_quantity('stage', 1.0)
    d.set_boundary({t: Reflective_boundary(d) for t in d.get_boundary_tags()})
    d.set_shear_closure('depth_slope', freeze_slope=True)
    d.add_sediment_fraction('sand', diameter=2.0e-4)
    assert np.allclose(d.sediment_slope_work, 0.003)


def test_freeze_slope_needs_the_depth_slope_closure():
    d = channel()
    with pytest.raises(ValueError, match='depth_slope'):
        d.set_shear_closure('energy_slope', freeze_slope=True)
    with pytest.raises(ValueError, match='max_slope'):
        d.set_shear_closure('depth_slope', max_slope=-1.0)
    d.set_shear_closure('depth_slope', max_slope=0.02, freeze_slope=True)
    d.add_sediment_fraction('sand', diameter=2.0e-4)
    assert 'capped at 0.02' in d.sediment_summary()
    assert 'frozen' in d.sediment_summary()


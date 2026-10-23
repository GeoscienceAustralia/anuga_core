"""[D-3] adaptation lag of the near-bed concentration.

Galappatti & Vreugdenhil's depth-integrated model relaxes the load toward
the same equilibrium as the instantaneous exchange, but at the rate
alpha v_s rather than d* v_s:

    E - D = alpha v_s (c_eq - c),
    1/alpha = a/h + (1 - a/h) exp[-1.5 (a/h)^(-1/6) w_s/u*]   (Armanini & Di
                                                             Silvio 1988)

Checked cell by cell in uniform flow held by Dirichlet boundaries, with the
sediment kernel's friction set to a uniform Manning n (larsen_lamb) so that
u* = sqrt(f_c) U is known exactly, and on a fixed bed so the flow stays put.
"""
import numpy as np
import pytest

import anuga
from anuga import Dirichlet_boundary, Reflective_boundary, rectangular_cross_domain

G = 9.8
H0 = 1.0
U0 = 1.0
K_S = 0.05


def uniform_flow(mode='legacy', diameter=2.0e-4, c0=0.01, tau_c_star=0.04,
                 near_bed='constant', **dep):
    d = rectangular_cross_domain(10, 4, len1=50.0, len2=20.0)
    d.set_flow_algorithm('DE1')
    d.set_compute_mode(mode)
    d.store = False
    d.set_quantity('elevation', 0.0)
    d.set_quantity('friction', 0.0)
    d.set_quantity('stage', H0)
    d.set_quantity('xmomentum', U0 * H0)
    Bd = Dirichlet_boundary([H0, U0 * H0, 0.0])
    Br = Reflective_boundary(d)
    d.set_boundary({'left': Bd, 'right': Bd, 'top': Br, 'bottom': Br})
    d.initialize_sediment_operator(bed_evolution=False)
    d.set_sediment_friction('larsen_lamb', k_s=K_S)
    dep.setdefault('adaptation', 'none')      # [D-5] is the default; 'plain' here is not
    d.set_deposition(law='d_star', near_bed=near_bed, **dep)
    d.add_sediment_fraction('sand', diameter=diameter, d_star=1.0,
                            tau_c_star=tau_c_star, initial_concentration=c0)
    d.set_tracer_boundary('sand', 'left', c0)
    d.set_tracer_boundary('sand', 'right', c0)
    return d


def expected_alpha(d, h=H0, U=U0):
    """Armanini & Di Silvio's alpha, written out independently of the kernel."""
    n = d.sediment_manning_ll
    f_c = G * n * n / h ** (1.0 / 3.0)
    ustar = np.sqrt(f_c) * U
    v_s = float(d.sediment_settling_velocity[0])
    a_h = max(float(d.sediment_reference_height[0]) / h, d.sediment_a_h_floor)
    return 1.0 / (a_h + (1.0 - a_h) * np.exp(-1.5 * a_h ** (-1.0 / 6.0) * v_s / ustar))


def source_rate(d, dt=0.05):
    """Rate of change of the depth-averaged concentration in an interior cell
    over one step of `dt`. The flow and the initial concentration are uniform
    and the boundaries carry the same concentration, so advection contributes
    nothing and the change is the bed exchange alone."""
    k = d.number_of_elements // 2
    c0 = float(d.get_tracer('sand')[k])
    d.evolve_max_timestep = dt
    for _ in d.evolve(yieldstep=dt, finaltime=dt):
        pass
    return (float(d.get_tracer('sand')[k]) - c0) / dt


@pytest.mark.parametrize('mode', ['legacy', 'unified'])
def test_the_armanini_lag_scales_deposition_by_alpha(mode):
    """Deposition only (no entrainment), well-mixed d* = 1: the rate must be
    alpha times the instantaneous rate, with alpha from the closed form."""
    plain = uniform_flow(mode, tau_c_star=0.0)
    lagged = uniform_flow(mode, tau_c_star=0.0, adaptation='armanini')
    r0 = source_rate(plain)
    r1 = source_rate(lagged)
    alpha = expected_alpha(lagged)
    assert r0 < 0.0 and r1 < 0.0
    assert alpha > 1.0, alpha
    assert r1 / r0 == pytest.approx(alpha, rel=1e-6)


def test_a_constant_alpha_is_applied_as_given():
    plain = uniform_flow(tau_c_star=0.0)
    lagged = uniform_flow(tau_c_star=0.0, adaptation='constant', adaptation_alpha=0.25)
    assert source_rate(lagged) / source_rate(plain) == pytest.approx(0.25, rel=1e-6)


def test_the_well_mixed_limit_has_no_lag():
    """w_s/u* -> 0 gives alpha -> 1: a washload is exchanged at the
    well-mixed rate whether the lag is on or not."""
    plain = uniform_flow(tau_c_star=0.0, diameter=5.0e-6)
    lagged = uniform_flow(tau_c_star=0.0, diameter=5.0e-6, adaptation='armanini')
    ratio = source_rate(lagged) / source_rate(plain)
    assert ratio == pytest.approx(1.0, abs=2e-3), ratio


def test_with_the_rouse_ratio_the_lag_slows_a_stratified_suspension():
    """With near_bed='rouse' the instantaneous rate is d* v_s (c_eq - c) and
    the lagged one alpha v_s (c_eq - c); for a coarse grain in this flow
    d* exceeds alpha, so the lag slows deposition, by exactly alpha/d*."""
    # a dilute suspension, so that d* c stays below the packing cap [L-4]
    c0 = 1.0e-4
    plain = uniform_flow(tau_c_star=0.0, diameter=5.0e-4, near_bed='rouse', c0=c0)
    lagged = uniform_flow(tau_c_star=0.0, diameter=5.0e-4, near_bed='rouse',
                          adaptation='armanini', c0=c0)
    # a short step: d* v_s / h is of order 1/s here, so a longer one would
    # see the concentration fall within the step and bias the rate
    r0, r1 = source_rate(plain, dt=1.0e-3), source_rate(lagged, dt=1.0e-3)
    ratio = r1 / r0
    assert 0.0 < ratio < 1.0, ratio
    # d* is not exposed per cell, but the instantaneous rate is d* v_s c / h
    d_star = -r0 * H0 / (float(plain.sediment_settling_velocity[0]) * c0)
    assert d_star > 1.0
    assert ratio == pytest.approx(expected_alpha(lagged) / d_star, rel=1e-2)


def test_the_equilibrium_is_unchanged_and_the_approach_is_slower():
    """Entrainment on, clear water in: both runs relax toward the same
    c_eq = E*/d* along the channel, the lagged one over the longer
    adaptation length u h/(alpha v_s) instead of u h/(d* v_s)."""
    prof = {}
    for name, dep in (('plain', {}),
                      ('lagged', dict(adaptation='constant', adaptation_alpha=0.3))):
        d = rectangular_cross_domain(60, 2, len1=600.0, len2=20.0)
        d.set_flow_algorithm('DE1')
        d.store = False
        d.set_quantity('elevation', 0.0)
        d.set_quantity('friction', 0.0)
        d.set_quantity('stage', H0)
        d.set_quantity('xmomentum', U0 * H0)
        Bd = Dirichlet_boundary([H0, U0 * H0, 0.0])
        Br = Reflective_boundary(d)
        d.set_boundary({'left': Bd, 'right': Bd, 'top': Br, 'bottom': Br})
        d.initialize_sediment_operator(bed_evolution=False)
        d.set_sediment_friction('larsen_lamb', k_s=K_S)
        d.set_deposition(law='d_star', near_bed='constant', **dep)
        d.add_sediment_fraction('sand', diameter=5.0e-4, d_star=1.0,
                                initial_concentration=0.0)
        d.set_tracer_boundary('sand', 'left', 0.0)
        for _ in d.evolve(yieldstep=100.0, finaltime=1500.0):
            pass
        x = d.centroid_coordinates[:, 0]
        c = d.get_tracer('sand')
        prof[name] = (float(c[np.abs(x - 30.0) < 5.0].mean()),
                      float(c[np.abs(x - 560.0) < 5.0].mean()))
    near_p, far_p = prof['plain']
    near_l, far_l = prof['lagged']
    assert far_p > 0.0
    # the same equilibrium far downstream: the adaptation length is
    # u h/(alpha v_s) ~ 50 m for the lagged run, 15 m for the plain one
    assert far_l == pytest.approx(far_p, rel=1e-2)
    # and the lagged run is further from it near the inflow
    assert near_l < 0.7 * near_p


def test_bad_settings():
    d = uniform_flow()
    with pytest.raises(ValueError):
        d.set_deposition(adaptation='galappatti')
    with pytest.raises(ValueError):
        d.set_deposition(adaptation='constant')
    with pytest.raises(ValueError):
        d.set_deposition(adaptation='constant', adaptation_alpha=0.0)
    with pytest.raises(ValueError):
        d.set_deposition(adaptation='armanini', adaptation_alpha=2.0)
    d.set_deposition(adaptation='armanini')
    assert 'adaptation [D-3]' in d.sediment_summary()
    d.set_deposition()                       # [D-5] is the default closure
    assert 'adaptation [D-5]' in d.sediment_summary()
    assert 'velocity profile' in d.sediment_summary()
    d.set_deposition(adaptation='none')
    assert 'adaptation' not in d.sediment_summary()
    with pytest.raises(ValueError):
        d.set_deposition(adaptation='none', velocity_profile=True)
    with pytest.raises(ValueError):
        d.set_deposition(layer_fraction=0.7)


# ------------------------------------------------ [D-4] carried near-bed ratio

def _channel(adapt, nx=60, length=300.0, bed=None, mode='legacy', c0=0.0, U=1.0):
    d = rectangular_cross_domain(nx, 2, len1=length, len2=10.0)
    d.set_flow_algorithm('DE1')
    d.set_compute_mode(mode)
    d.store = False
    d.set_quantity('elevation', 0.0 if bed is None else bed)
    d.set_quantity('friction', 0.0)
    d.set_quantity('stage', H0)
    d.set_quantity('xmomentum', U * H0)
    Bd_in = Dirichlet_boundary([H0, U * H0, 0.0])
    Bd_out = Dirichlet_boundary([H0, U * H0, 0.0])
    Br = Reflective_boundary(d)
    d.set_boundary({'left': Bd_in, 'right': Bd_out, 'top': Br, 'bottom': Br})
    d.initialize_sediment_operator(bed_evolution=False)
    d.set_sediment_friction('larsen_lamb', k_s=K_S)
    d.set_deposition(law='d_star', near_bed='rouse', adaptation=adapt,
                     layer_fraction=0.0, velocity_profile=False)   # the bare [D-5]
    d.add_sediment_fraction('sand', diameter=2.0e-4, initial_concentration=c0)
    d.set_tracer_boundary('sand', 'left', c0)
    return d


def test_the_centroid_fit_reproduces_the_quadrature():
    """(z_c - a)/h of the Rouse profile from the kernel's coefficients, via
    the Python port, against direct integration."""
    from scipy.integrate import quad
    from anuga import Domain
    worst = 0.0
    for Z in (0.02, 0.1, 0.5, 1.0, 1.5, 2.4):
        for a_h in (0.002, 0.01, 0.05, 0.12):
            f = lambda z: (((1.0 - z) / (1.0 - a_h)) * (a_h / z)) ** Z
            zc = quad(lambda z: z * f(z), a_h, 1.0, limit=400)[0] / quad(f, a_h, 1.0, limit=400)[0]
            fit = Domain.rouse_centroid(Z, a_h)
            worst = max(worst, abs(fit - (zc - a_h)) / (zc - a_h))
    assert worst < 0.015, worst
    assert abs(Domain.rouse_centroid(0.01, 0.01) - 0.49) < 0.01     # well mixed: h/2


def test_the_python_d_star_port_matches_the_kernel_fit():
    """The port evaluates the same polynomial the kernel does; check it
    against the values the rouse test reads out of the kernel source."""
    from anuga import Domain
    from anuga.shallow_water.tests import test_sediment_rouse as R
    for Z, a_h in ((0.05, 0.01), (0.5, 0.02), (1.0, 0.01), (2.0, 0.1)):
        assert Domain.rouse_d_star(Z, a_h) == pytest.approx(R.dstar_fit(Z, a_h), rel=1e-12)


@pytest.mark.parametrize('mode', ['legacy', 'unified'])
def test_in_uniform_flow_the_carried_ratio_changes_nothing(mode):
    """Clear water loaded from the bed in uniform flow: d* never rises, so
    the one-sided lag never fires and the result is bit-identical."""
    plain = _channel('none', mode=mode)
    carried = _channel('carried', mode=mode)
    for d in (plain, carried):
        for _ in d.evolve(yieldstep=100.0, finaltime=400.0):
            pass
    assert carried.number_of_tracers == 2
    assert 'sand_nearbed_ratio' in carried.get_tracer_names() if hasattr(carried, 'get_tracer_names') else True
    assert np.array_equal(plain.get_tracer('sand'), carried.get_tracer('sand'))
    r = carried.get_tracer('sand_nearbed_ratio')
    assert r.min() > 1.0 and np.ptp(r) / r.mean() < 1e-6


def test_over_a_deepening_the_carried_ratio_lags_and_deposition_is_delayed():
    """A step down in the bed slows the flow; d* rises there. Without the
    lag the near-bed ratio is the new d* at once and deposition is largest
    in the first cells past the step; with the carried ratio it starts from
    the upstream value and rises over the settling time, so the bed
    exchange just past the step is smaller and the load persists further."""
    bed = lambda x, y: np.where(x > 150.0, -1.0, 0.0)   # depth doubles past x = 150
    results = {}
    for adapt in ('none', 'carried'):
        d = _channel(adapt, bed=bed, c0=2.0e-3)
        for _ in d.evolve(yieldstep=100.0, finaltime=500.0):
            pass
        x = d.centroid_coordinates[:, 0]
        c = d.get_tracer('sand')
        results[adapt] = (float(c[(x > 152) & (x < 162)].mean()),
                          float(c[(x > 250) & (x < 280)].mean()), d)
    near_p, far_p, _ = results['none']
    near_c, far_c, d = results['carried']
    # more load survives just past the step with the lag
    assert near_c > near_p * 1.02, (near_c, near_p)
    # and the ratio there is below its downstream equilibrium, rising toward it
    r = d.get_tracer('sand_nearbed_ratio'); x = d.centroid_coordinates[:, 0]
    assert r[(x > 152) & (x < 158)].mean() < r[(x > 280)].mean()


def test_fractions_must_be_added_before_the_first_evolve():
    d = _channel('carried')
    for _ in d.evolve(yieldstep=1.0, finaltime=1.0):
        pass
    with pytest.raises(ValueError):
        d.add_sediment_fraction('silt', diameter=5.0e-5)
    assert 'adaptation [D-4]' in d.sediment_summary()


# --------------------------------------------- [D-5] two-layer suspension

@pytest.mark.parametrize('mode', ['legacy', 'unified'])
def test_two_layer_keeps_the_equilibrium_and_lengthens_the_adaptation(mode):
    """Clear water loaded from the bed: the far-downstream equilibrium of
    the two-layer model must be the single-layer one (the exchange is set
    to reproduce d*), while the load near the inflow is lower because the
    upper layer fills only by exchange."""
    plain = _channel('none', mode=mode)
    two = _channel('two_layer', mode=mode)
    for d in (plain, two):
        d.set_deposition(law='d_star', near_bed='rouse', reference_height_floor=0.1,
                         adaptation='two_layer' if d is two else 'none',
                         layer_fraction=0.0, velocity_profile=False)
        for _ in d.evolve(yieldstep=100.0, finaltime=900.0):
            pass
    x = plain.centroid_coordinates[:, 0]
    c_p, c_t = plain.get_tracer('sand'), two.get_tracer('sand')
    far = x > 270.0
    near = (x > 25.0) & (x < 35.0)
    assert two.number_of_tracers == 2
    assert c_t[far].mean() == pytest.approx(c_p[far].mean(), rel=2e-3)
    assert c_t[near].mean() < 0.85 * c_p[near].mean()
    m2 = two.get_tracer('sand_upper')
    share = (m2 / np.maximum(c_t, 1e-30))[far].mean()
    assert 0.0 < share < 1.0
    assert 'adaptation [D-5]' in two.sediment_summary()


def test_two_layer_partition_matches_the_rouse_ratio_at_equilibrium():
    """At equilibrium the lower-layer concentration over the depth-averaged
    one must be the fitted d* of the local flow, by construction of K."""
    from anuga import Domain
    d = _channel('two_layer')
    d.set_deposition(law='d_star', near_bed='rouse', reference_height_floor=0.1,
                     adaptation='two_layer', layer_fraction=0.0, velocity_profile=False)
    for _ in d.evolve(yieldstep=100.0, finaltime=900.0):
        pass
    x = d.centroid_coordinates[:, 0]
    far = x > 270.0
    c = d.get_tracer('sand')[far].mean()
    m2 = d.get_tracer('sand_upper')[far].mean()          # m2 / h as a tracer value
    a_h = 0.1                                            # the floor, since a = 2 d << 0.1 h
    c1 = (c - m2) / a_h                                  # (m - m2) / h1, per unit h
    n = d.sediment_manning_ll
    f_c = G * n * n / H0 ** (1.0 / 3.0)
    ustar = np.sqrt(f_c) * U0
    Z = float(d.sediment_settling_velocity[0]) / (0.41 * ustar)
    assert c1 / c == pytest.approx(Domain.rouse_d_star(Z, a_h), rel=2e-2)


def _profile_channel(mode, vp, nx=40, length=40.0):
    """Sloping channel with a log-law-consistent Manning n, held by a
    Dirichlet inflow at an equilibrium-ish concentration."""
    d = rectangular_cross_domain(nx, 2, len1=length, len2=2.0)
    d.set_flow_algorithm('DE1')
    d.set_compute_mode(mode)
    d.store = False
    d.set_quantity('elevation', lambda x, y: 0.002 * (length - x))
    d.set_quantity('stage', lambda x, y: 0.002 * (length - x) + 0.4)
    d.set_quantity('friction', 0.02)
    d.set_quantity('xmomentum', 0.2)
    Bi = Dirichlet_boundary([0.48, 0.2, 0.0])
    Bo = Dirichlet_boundary([0.4, 0.2, 0.0])
    Br = Reflective_boundary(d)
    d.set_boundary({'left': Bi, 'right': Bo, 'top': Br, 'bottom': Br})
    d.initialize_sediment_operator(porosity=0.4, bed_evolution=False)
    d.set_sediment_friction('larsen_lamb', k_s=0.025)
    d.set_deposition(law='d_star', near_bed='rouse', reference_height_floor=0.1,
                     adaptation='two_layer', layer_fraction=0.0, velocity_profile=vp)
    d.add_sediment_fraction('sand', diameter=1.4e-4, initial_concentration=1e-4)
    d.set_tracer_boundary('sand', 'left', 1e-4)
    for _ in d.evolve(yieldstep=30.0, finaltime=60.0):
        pass
    return d


def test_the_velocity_profile_slows_the_near_bed_layer_and_leads_the_upper():
    """[D-5v] the speed factors are the log-law layer means: the total (most
    of it in the near-bed layer) below the depth mean, the upper layer just
    above it, and the two-mode results identical."""
    d = _profile_channel('legacy', True)
    sf = d.tracer_speed_factor
    assert sf.shape == (2, d.number_of_elements)
    assert 0.3 < sf[0].min() and sf[0].max() < 1.0            # total lags the flow
    assert 1.0 < sf[1].min() and sf[1].max() < 1.2            # upper layer leads
    # log law at f1 = 0.1 (the floor), L = kappa / sqrt(f_c) with the
    # sediment friction's Manning n: u2/u = 1 - f1 ln f1 / ((1 - f1) L)
    k = d.number_of_elements // 2
    h = d.quantities['stage'].centroid_values[k] - d.quantities['elevation'].centroid_values[k]
    n = d.sediment_manning_ll
    L = 0.41 / np.sqrt(G * n * n / h ** (1.0 / 3.0))
    assert sf[1][k] == pytest.approx(1.0 - 0.1 * np.log(0.1) / (0.9 * L), rel=1e-6)
    assert 'velocity profile' in d.sediment_summary()
    u = _profile_channel('unified', True)
    assert np.abs(u.get_tracer('sand') - d.get_tracer('sand')).max() < 1e-9
    # and it is a real change against the depth-averaged advection
    p = _profile_channel('legacy', False)
    assert p.tracer_speed_factor is None
    assert np.abs(p.get_tracer('sand') - d.get_tracer('sand')).max() > 1e-6


def test_the_velocity_profile_conserves_tracer_mass():
    """A closed basin: the factors move mass between cells but never create
    it, whatever the donor cell's factor."""
    d = rectangular_cross_domain(10, 4, len1=50.0, len2=20.0)
    d.set_flow_algorithm('DE1')
    d.store = False
    d.set_quantity('elevation', 0.0)
    d.set_quantity('friction', 0.02)
    d.set_quantity('stage', lambda x, y: 1.0 + 0.2 * (x > 25.0))    # sloshes
    d.set_boundary({t: Reflective_boundary(d) for t in d.get_boundary_tags()})
    d.initialize_sediment_operator(bed_evolution=False)
    d.set_sediment_friction('larsen_lamb', k_s=0.05)
    d.set_deposition(law='threshold', tau_d=0.0, adaptation='two_layer',
                     velocity_profile=True)                          # no bed exchange
    d.add_sediment_fraction('sand', diameter=2.0e-4, tau_c_star=1.0e9,
                            initial_concentration=1e-3)
    areas = d.areas
    m0 = float((d.tracer_conserved_values[0] * areas).sum())
    for _ in d.evolve(yieldstep=5.0, finaltime=20.0):
        pass
    m1 = float((d.tracer_conserved_values[0] * areas).sum())
    assert m1 == pytest.approx(m0, rel=1e-12)


def test_the_velocity_profile_needs_the_two_layer_model():
    d = uniform_flow()
    with pytest.raises(ValueError):
        d.set_deposition(law='d_star', near_bed='rouse', adaptation='carried',
                         velocity_profile=True)


@pytest.mark.parametrize('mode', ['legacy', 'unified'])
def test_the_rouse_correction_flattens_the_profile(mode):
    """[S-2b] van Rijn's beta divides Z by 1 + 2 (w_s/u*)^2 and the scale
    multiplies it; both lower d* and with it the deposition rate, and the
    kernel's d* must be the fitted one at the corrected Z."""
    from anuga import Domain
    plain = uniform_flow(mode, tau_c_star=0.0, near_bed='rouse', reference_height_floor=0.1)
    beta = uniform_flow(mode, tau_c_star=0.0, near_bed='rouse', reference_height_floor=0.1,
                        rouse_beta='van_rijn')
    scaled = uniform_flow(mode, tau_c_star=0.0, near_bed='rouse', reference_height_floor=0.1,
                          rouse_scale=0.5)
    r0, r1, r2 = source_rate(plain), source_rate(beta), source_rate(scaled)
    n = plain.sediment_manning_ll
    ustar = np.sqrt(G * n * n / H0 ** (1.0 / 3.0)) * U0
    v_s = float(plain.sediment_settling_velocity[0])
    Z = v_s / (0.41 * ustar)
    b = min(2.0, 1.0 + 2.0 * (v_s / ustar) ** 2)
    assert r0 < 0.0 and r1 < 0.0 and r2 < 0.0
    assert r1 / r0 == pytest.approx(Domain.rouse_d_star(Z / b, 0.1) / Domain.rouse_d_star(Z, 0.1), rel=1e-3)
    assert r2 / r0 == pytest.approx(Domain.rouse_d_star(0.5 * Z, 0.1) / Domain.rouse_d_star(Z, 0.1), rel=1e-3)
    assert 'Rouse number' in beta.sediment_summary()
    with pytest.raises(ValueError):
        plain.set_deposition(near_bed='rouse', rouse_beta='other')
    with pytest.raises(ValueError):
        plain.set_deposition(near_bed='rouse', rouse_scale=0.0)


@pytest.mark.parametrize('mode', ['legacy', 'unified'])
def test_a_layered_adaptation_coexists_with_a_passive_tracer(mode):
    """The near-bed tracers go at the END of the tracer list, so an ordinary
    passive tracer added after the fractions is untouched by the sediment
    kernels and carries on being advected. Nothing may be added after the
    first evolve, when they are registered."""
    def build(adapt):
        d = rectangular_cross_domain(6, 6, len1=50.0, len2=50.0)
        d.set_flow_algorithm('DE0')
        d.set_compute_mode(mode)
        d.store = False
        d.set_quantity('elevation', 0.0)
        d.set_quantity('stage', 1.0)
        d.set_quantity('xmomentum', 0.3)
        d.set_boundary({t: Reflective_boundary(d) for t in d.get_boundary_tags()})
        d.initialize_sediment_operator(porosity=0.3, bed_evolution=False)
        d.set_deposition(law='d_star', near_bed='rouse', reference_height_floor=0.1,
                         adaptation=adapt)
        d.add_sediment_fraction('sand', diameter=2.0e-4, tau_c_star=1.0e9,
                                initial_concentration=0.01)
        d.add_tracer('salt', initial_value=0.5)
        for _ in d.evolve(yieldstep=1.0, finaltime=2.0):
            pass
        return d
    plain = build('none')
    two = build('two_layer')
    assert plain.number_of_tracers == 2              # sand, salt
    assert two.number_of_tracers == 3                # sand, salt, sand_upper
    assert two.sediment_nearbed_base == 2            # the near-bed tracer is last
    # the passive tracer is untouched by the near-bed machinery
    assert np.allclose(two.get_tracer('salt'), plain.get_tracer('salt'), rtol=1e-12)
    # and the sediment did something different, as the closure intends
    assert not np.allclose(two.get_tracer('sand'), plain.get_tracer('sand'), rtol=1e-6)
    with pytest.raises(ValueError):
        two.add_tracer('late')


@pytest.mark.parametrize('mode', ['legacy', 'unified'])
def test_the_two_layer_default_reduces_to_the_instantaneous_exchange_at_d_star_one(mode):
    """[D-5] with the well-mixed near-bed ratio d* = 1 the partition has
    rho = 1, so the lower-layer concentration IS the depth-averaged one and
    the closure is exactly [D-1]. That is why the analytical sediment
    validation cases, which all run near_bed='constant' with d* = 1, are
    unaffected by [D-5] becoming the default."""
    def run(adapt):
        d = uniform_flow(mode, tau_c_star=0.0, near_bed='constant',
                         adaptation=adapt)
        for _ in d.evolve(yieldstep=1.0, finaltime=5.0):
            pass
        return d.get_tracer('sand').copy()
    plain = run('none')
    two = run('two_layer')
    assert plain.mean() > 0.0
    assert np.allclose(two, plain, rtol=1e-9, atol=1e-14)


def test_the_summary_reports_the_whole_configuration():
    """Every switch a run can set should be visible in sediment_summary(),
    so a log says what was actually computed."""
    d = uniform_flow(near_bed='rouse', reference_height_floor=0.1)
    d.set_deposition(law='d_star', near_bed='rouse', reference_height_floor=0.1)
    d.set_bedload('wong_parker_eq24', min_depth=4.0e-4)
    d.set_shear_amplification(lambda x, y: 1.0 + 3.0 * (x > 25.0))
    s = d.sediment_summary()
    assert 'adaptation [D-5]' in s                 # the default closure
    assert 'near-bed layer 0.2 h' in s
    assert 'velocity profile' in s
    assert 'c_b carried' in s and '[D-5]' in s     # not the bare [D-1] form
    assert 'bedload h_min [K-7]' in s and '0.0004' in s
    assert 'amplified by 1 to 4 per cell' in s and '[T-16]' in s
    # and the plain closure still reads as itself
    plain = uniform_flow(near_bed='rouse', adaptation='none')
    assert 'D = d* c v_s   [D-1]' in plain.sediment_summary()
    assert 'adaptation' not in plain.sediment_summary()


# ------------------------------------------------- the summary stays complete
#
# A setting that sediment_summary() does not report is a setting a run's log
# hides. Three had crept in that way. Rather than trust the next author to
# remember, every `sediment_*` attribute is inventoried below: either it names
# the closure that has to be switched on for the summary to mention it, and
# the test proves that changing it changes the summary, or it is exempt with
# the reason. A new attribute in neither list fails the test.

def _summary_domain():
    d = rectangular_cross_domain(6, 4, len1=30.0, len2=20.0)
    d.set_flow_algorithm('DE1')
    d.store = False
    d.set_quantity('elevation', 0.0)
    d.set_quantity('stage', 1.0)
    d.set_quantity('xmomentum', 1.0)
    d.set_boundary({t: Reflective_boundary(d) for t in d.get_boundary_tags()})
    d.initialize_sediment_operator(porosity=0.4, bed_evolution=True)
    d.add_sediment_fraction('sand', diameter=2.0e-4, initial_concentration=1.0e-3)
    d.set_bedload('wong_parker_eq24')
    return d


# closure -> what to call so that its own parameters are in play
ACTIVATE = {
    'base': lambda d: None,
    'smith_mclean': lambda d: d.set_bed_material('noncohesive'),
    'de_leeuw': lambda d: d.set_bed_material('noncohesive', entrainment='de_leeuw'),
    'cohesive': lambda d: d.set_bed_material('cohesive'),
    'partheniades': lambda d: d.set_bed_material('partheniades'),
    'larsen_lamb': lambda d: d.set_sediment_friction('larsen_lamb', k_s=0.02),
    'wilson': lambda d: d.set_sediment_friction('wilson', grain_size=0.02),
    'depth_slope': lambda d: d.set_shear_closure('depth_slope', max_slope=0.1),
    'repose': lambda d: d.set_angle_of_repose(32.0, relax=0.5),
    'threshold_dep': lambda d: d.set_deposition(law='threshold', tau_d=0.5),
    'rouse': lambda d: d.set_deposition(law='d_star', near_bed='rouse'),
    'constant_lag': lambda d: d.set_deposition(near_bed='rouse', adaptation='constant',
                                               adaptation_alpha=2.0),
    'rouse_corr': lambda d: d.set_deposition(near_bed='rouse', rouse_beta='van_rijn',
                                             rouse_scale=0.7),
    'bedload_ramp': lambda d: d.set_bedload('wong_parker_eq24', min_depth=4.0e-4),
}

# every scalar sediment_* attribute, and the closure under which the summary
# must mention it
REPORTED = {
    'sediment_porosity': 'base', 'sediment_rho_w': 'base',
    'sediment_c_max': 'base', 'sediment_c_pack': 'base',
    'sediment_morphological_factor': 'base',
    'sediment_bedload_K': 'base', 'sediment_bedload_m': 'base',
    'sediment_bedload_tau_c_star': 'base',
    'sediment_bedload_h_min': 'bedload_ramp',
    'sediment_gamma0': 'smith_mclean',
    'sediment_tau_crit': 'cohesive', 'sediment_K_e': 'cohesive',
    'sediment_K_partheniades': 'partheniades',
    'sediment_dl_A': 'de_leeuw', 'sediment_dl_alpha': 'de_leeuw',
    'sediment_dl_beta': 'de_leeuw', 'sediment_dl_ks': 'de_leeuw',
    'sediment_dl_ks_factor': 'de_leeuw', 'sediment_dl_threshold': 'de_leeuw',
    'sediment_manning_ll': 'larsen_lamb', 'sediment_wilson_D': 'wilson',
    'sediment_max_slope': 'depth_slope',
    'sediment_repose_tan': 'repose', 'sediment_repose_relax': 'repose',
    'sediment_repose_max_sweeps': 'repose',
    'sediment_tau_d': 'threshold_dep',
    'sediment_a_h_floor': 'rouse',
    'sediment_adaptation_alpha': 'constant_lag',
    'sediment_rouse_scale': 'rouse_corr',
    'sediment_layer_fraction': 'base', 'sediment_exchange_factor': 'base',
}

EXEMPT = {
    # mode selectors: the summary names the closure in words instead
    'sediment_erosion_mode', 'sediment_deposition_mode', 'sediment_d_star_mode',
    'sediment_shear_closure', 'sediment_friction_mode', 'sediment_bedload_mode',
    'sediment_adaptation_mode', 'sediment_rouse_beta_mode',
    'sediment_wilson_bed', 'sediment_vegetation_shear',
    'sediment_velocity_profile', 'sediment_bed_evolution',
    'sediment_slope_frozen', 'sediment_has_z_base', 'sediment_bed_composition',
    # per-fraction or per-cell state, reported in their own sections
    'sediment_diameter', 'sediment_R', 'sediment_tau_c_star', 'sediment_d_star',
    'sediment_settling_velocity', 'sediment_reference_height',
    'sediment_z_base', 'sediment_shear_factor', 'sediment_nearbed_base',
    'sediment_bedload_open', 'sediment_bedload_supply',
    # bed composition: its own summary line (active layer and mean fractions)
    # when it is on; with it off these are 0 / None
    'sediment_active_layer', 'sediment_bed_active', 'sediment_bed_substrate',
    # kernel scratch, not settings
    'sediment_qbx', 'sediment_qby', 'sediment_qba', 'sediment_slope_work',
    'sediment_source_limited', 'sediment_bed_exhausted', 'sediment_repose_dz',
}


def test_every_sediment_setting_is_inventoried():
    """A new sediment_* attribute must be classified, so it cannot be added
    without deciding whether the summary should report it."""
    d = _summary_domain()
    found = {k for k in vars(d) if k.startswith('sediment_')}
    missing = sorted(found - set(REPORTED) - EXEMPT)
    assert missing == [], (
        'these sediment settings are in neither REPORTED nor EXEMPT: %s. Add a '
        'line to sediment_summary() and list it in REPORTED, or explain it in '
        'EXEMPT.' % missing)


@pytest.mark.parametrize('name,closure', sorted(REPORTED.items()))
def test_a_setting_that_changes_changes_the_summary(name, closure):
    """Under the closure that uses it, perturbing a setting must show."""
    d = _summary_domain()
    ACTIVATE[closure](d)
    before = d.sediment_summary()
    old = getattr(d, name)
    setattr(d, name, old * 2.0 + 1.0)
    after = d.sediment_summary()
    setattr(d, name, old)
    assert after != before, (
        '%s is not reported by sediment_summary() under %r' % (name, closure))

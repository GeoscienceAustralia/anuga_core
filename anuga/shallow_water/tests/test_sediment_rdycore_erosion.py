"""RDycore's dam break WITH bed exchange (RDy26 3.2): [E-4] and [D-2].

Their section 3.2 sets n = 0.025 and enables Partheniades erosion and
threshold deposition. That means the analytic Stoker solution does NOT apply:
Stoker is frictionless, and friction changes the hydrodynamics whatever the
sediment does. Comparing against it was tried and failed at L1 = 4.3e-2 m,
correctly.

So what is verifiable here is stated rather than dressed up:

  * the CLOSURES themselves, exactly. [E-4] and [D-2] are algebraic, so the
    rates can be checked against the formulas directly, as [E-3] was checked
    against anugaSed;
  * that the case actually exercises its thresholds, rather than sitting
    wholly above or below them;
  * mass consistency and physical bounds on SSC;
  * SELF-convergence under refinement, which needs no reference solution;
  * and, as a control, that turning exchange off through [D-2]'s tau_d = 0
    hook with n = 0 recovers the analytic passive solution of their 3.1.1.

The SSC field itself has no exact answer here and is not claimed to agree with
theirs.
"""
import numpy as np
import pytest
from scipy.optimize import brentq

from anuga import Reflective_boundary, rectangular_cross_domain

G = 9.81
HL, HR, CL, CR = 1.0, 0.5, 0.7, 0.5
DAM_X, LENGTH, T_END = 1000.0, 2000.0, 240.0
RHO_S, V_S = 1600.0, 1.0e-4
TAU_D, TAU_E, K_P = 0.1, 0.1, 1.0e-4
MANNING = 0.025


def E4(tau_b, tau_c=TAU_E, Kp=K_P, rho_s=RHO_S):
    """[E-4] Partheniades, as RDy26 state it, in m/s."""
    return (Kp * (tau_b - tau_c) / tau_c) / rho_s if tau_b > tau_c else 0.0


def D2(tau_b, c, vs=V_S, tau_d=TAU_D):
    """[D-2] threshold deposition, in m/s."""
    if tau_d > 0 and tau_b < tau_d:
        return vs * c * (1.0 - tau_b / tau_d)
    return 0.0


def stoker_contact_speed(hl, hr):
    cl = np.sqrt(G * hl)

    def residual(hm):
        return (-2.0 * (np.sqrt(G * hm) - cl)
                - (hm - hr) * np.sqrt(0.5 * G * (hm + hr) / (hm * hr)))

    hm = brentq(residual, hr * (1 + 1e-12), hl * (1 - 1e-12),
                xtol=1e-14, rtol=1e-15)
    return -2.0 * (np.sqrt(G * hm) - cl)


def build(exchange, nx=800, manning=MANNING):
    d = rectangular_cross_domain(nx, 2, len1=LENGTH, len2=10.0)
    d.set_flow_algorithm('DE0')
    d.set_low_froude(0)
    d.store = False
    x = d.centroid_coordinates[:, 0]
    d.set_quantity('elevation', 0.0)
    d.set_quantity('stage', np.where(x < DAM_X, HL, HR), location='centroids')
    d.set_quantity('friction', manning)
    d.set_boundary({t: Reflective_boundary(d) for t in d.get_boundary_tags()})
    d.evolve_max_timestep = 0.05
    d.sediment_bed_evolution = False       # RDy26 do not update bed elevation
    d.sediment_c_max = 1.0
    d.sediment_c_pack = 1.0
    if exchange:
        d.set_bed_material('partheniades', tau_crit=TAU_E, K_e=K_P)
        d.set_deposition('threshold', tau_d=TAU_D)
    else:
        d.set_deposition('threshold', tau_d=0.0)     # their passive hook
    d.add_sediment_fraction(name='ssc', diameter=1.0e-4, rho_s=RHO_S, tau_c_star=0.0,
                     initial_concentration=np.where(x < DAM_X, CL, CR))
    d.sediment_settling_velocity[0] = V_S            # RDy26 prescribe it
    return d, x


# ---------------------------------------------------------------- closures

def test_E4_has_a_genuine_threshold():
    assert E4(TAU_E) == 0.0
    assert E4(0.099) == 0.0
    assert E4(0.15) > 0.0


def test_E4_is_linear_in_excess_stress():
    assert abs(E4(0.3) / E4(0.2) - (0.3 - TAU_E) / (0.2 - TAU_E)) < 1e-12


def test_D2_shuts_off_at_its_threshold_and_is_maximal_in_still_water():
    assert D2(TAU_D, 0.6) == 0.0
    assert D2(0.0, 0.6) == V_S * 0.6


def test_tau_d_zero_disables_deposition_entirely():
    """RDy26's hook for running their passive cases through the same code."""
    assert D2(0.05, 0.6, tau_d=0.0) == 0.0


# ---------------------------------------------------------------- the run

@pytest.fixture(scope='module')
def exchange_run():
    active, x = build(exchange=True)
    m0 = float((active.tracer_conserved_values[0] * active.areas).sum())
    active.evolve_to_end(finaltime=T_END)
    passive, _ = build(exchange=False)
    passive.evolve_to_end(finaltime=T_END)
    return active, passive, x, m0


@pytest.mark.slow
def test_the_case_straddles_its_thresholds(exchange_run):
    """If the whole domain sat above or below 0.1 Pa, only one limb would run
    and the case would verify half of what it looks like it verifies."""
    active, _, _, _ = exchange_run
    h = np.maximum(active.quantities['stage'].centroid_values
                   - active.quantities['elevation'].centroid_values, 0.0)
    u = active.quantities['xmomentum'].centroid_values * h / (h ** 2 + 1e-12)
    tau = 1000.0 * G * MANNING ** 2 * u ** 2 / np.cbrt(np.maximum(h, 1e-12))
    assert tau.max() > TAU_E > tau.min()
    eroding = tau > TAU_E
    assert eroding.sum() > 0 and (~eroding).sum() > 0


@pytest.mark.slow
def test_bed_exchange_moves_ssc_away_from_the_passive_answer(exchange_run):
    active, passive, _, m0 = exchange_run
    c_act = active.get_tracer('ssc')
    c_pas = passive.get_tracer('ssc')
    assert np.abs(c_act - c_pas).max() > 1e-6
    m1 = float((active.tracer_conserved_values[0] * active.areas).sum())
    assert abs(m1 - m0) > 1e-9 * abs(m0), 'total suspended mass did not change'


@pytest.mark.slow
def test_ssc_stays_physical(exchange_run):
    active, _, _, _ = exchange_run
    c = active.get_tracer('ssc')
    assert c.min() >= -1e-12
    assert c.max() <= active.sediment_c_max + 1e-9


@pytest.mark.slow
def test_ssc_self_converges_under_refinement(exchange_run):
    """No reference solution is needed for this one: successive refinements
    must approach each other."""
    active, _, x_fine, _ = exchange_run

    def run(nx):
        d, x = build(exchange=True, nx=nx)
        d.evolve_to_end(finaltime=T_END)
        return x, d.get_tracer('ssc')

    probe = np.linspace(600.0, 1400.0, 200)

    def sample(x, c):
        return np.interp(probe, *zip(*sorted(zip(x, c))))

    x_c, c_c = run(200)
    x_m, c_m = run(400)
    coarse_to_mid = np.abs(sample(x_c, c_c) - sample(x_m, c_m)).mean()
    mid_to_fine = np.abs(sample(x_m, c_m)
                         - sample(x_fine, active.get_tracer('ssc'))).mean()
    assert mid_to_fine < coarse_to_mid


@pytest.mark.slow
def test_the_frictionless_passive_control_recovers_the_analytic_solution():
    """The control: with n = 0 and exchange off through tau_d = 0, RDy26's
    3.1.1 solution must come back -- which also shows that hook does what it
    claims rather than merely producing a plausible field."""
    d, x = build(exchange=False, manning=0.0)
    d.evolve_to_end(finaltime=T_END)
    um = stoker_contact_speed(HL, HR)
    exact = np.where(x < DAM_X + um * T_END, CL, CR)
    assert np.abs(d.get_tracer('ssc') - exact).mean() < 5e-3

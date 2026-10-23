"""RDycore's passive sediment-transport benchmarks (RDy26 3.1).

Two cases in which the sediment is a passive tracer, so the right answer is
known independently of any sediment closure:

  lake at rest    a discontinuous SSC field over a bump, in still water. The
                  free surface, the momentum and the SSC field must ALL be
                  unchanged. Well-balancedness plus zero transport.
  wet-bed dam     Stoker's analytical solution. The SSC contact rides the flow
                  at the star-state velocity, so the sediment field must be a
                  step at the CONTACT and flat across both the rarefaction and
                  the shock -- sediment rides the flow, not the water waves.

The second case is what verifies the advective half of [G-3] to machine
precision, which is why test_sediment_mms.py can restrict itself to the source
term.
"""
import numpy as np
import pytest
from scipy.optimize import brentq

from anuga import Reflective_boundary, rectangular_cross_domain

G = 9.81
HL, HR, CL, CR = 1.0, 0.5, 0.7, 0.5
DAM_X, DAM_T, LENGTH = 1000.0, 240.0, 2000.0


# ---------------------------------------------------------------- lake at rest

@pytest.fixture(scope='module')
def lake():
    d = rectangular_cross_domain(200, 2, len1=20.0, len2=2.0)
    d.set_flow_algorithm('DE0')
    d.set_low_froude(0)
    d.store = False
    x = d.centroid_coordinates[:, 0]
    d.set_quantity('elevation', 0.2 * np.exp(-((x - 10.0) / 1.0) ** 2),
                   location='centroids')
    d.set_quantity('stage', 1.0)
    d.set_quantity('friction', 0.0)
    d.set_boundary({t: Reflective_boundary(d) for t in d.get_boundary_tags()})
    # 1 g/L is 1e-3 by volume. d* = 0 disables exchange: this is a passive case.
    c0 = np.where((x > 8.0) & (x < 12.0), 1.0e-3, 0.0)
    d.add_sediment_fraction(name='ssc', diameter=1e-4, tau_c_star=0.0, d_star=0.0,
                     initial_concentration=c0)
    w0 = d.quantities['stage'].centroid_values.copy()
    z0 = d.quantities['elevation'].centroid_values.copy()
    d.evolve_to_end(finaltime=100.0)
    return d, c0, w0, z0


def test_the_free_surface_is_preserved_over_the_bump(lake):
    d, _, w0, _ = lake
    assert np.allclose(d.quantities['stage'].centroid_values, w0,
                       rtol=0, atol=1e-12)


def test_no_spurious_momentum_is_generated(lake):
    d = lake[0]
    assert np.abs(d.quantities['xmomentum'].centroid_values).max() < 1e-12


def test_the_discontinuous_ssc_field_does_not_move(lake):
    """u = 0 and there is no diffusion, so transport must be exactly zero. A
    scheme with any numerical diffusion smears this step."""
    d, c0, _, _ = lake
    assert np.allclose(d.get_tracer('ssc'), c0, rtol=0, atol=1e-14)


def test_the_bed_does_not_move(lake):
    d, _, _, z0 = lake
    assert np.allclose(d.quantities['elevation'].centroid_values, z0,
                       rtol=0, atol=1e-14)


# ---------------------------------------------------------------- dam break

def stoker(hl, hr, t, xs):
    """Classical wet-bed dam break: left rarefaction, right shock, contact
    between them. Solved for the star state by matching the rarefaction and
    shock relations."""
    cl = np.sqrt(G * hl)

    def residual(hm):
        u_raref = -2.0 * (np.sqrt(G * hm) - cl)
        u_shock = (hm - hr) * np.sqrt(0.5 * G * (hm + hr) / (hm * hr))
        return u_raref - u_shock

    hm = brentq(residual, hr * (1 + 1e-12), hl * (1 - 1e-12),
                xtol=1e-14, rtol=1e-15)
    um = -2.0 * (np.sqrt(G * hm) - cl)
    cm = np.sqrt(G * hm)
    shock_s = hm * um / (hm - hr)          # Rankine-Hugoniot mass jump

    h = np.empty_like(xs)
    xi = (xs - DAM_X) / t
    for i, s in enumerate(xi):
        if s <= -cl:
            h[i] = hl
        elif s <= um - cm:
            h[i] = ((-s + 2.0 * cl) / 3.0) ** 2 / G      # rarefaction fan
        elif s <= shock_s:
            h[i] = hm
        else:
            h[i] = hr
    return h, um, hm, shock_s


def run_dambreak(nx):
    d = rectangular_cross_domain(nx, 2, len1=LENGTH, len2=10.0)
    d.set_flow_algorithm('DE0')
    d.set_low_froude(0)
    d.store = False
    x = d.centroid_coordinates[:, 0]
    d.set_quantity('elevation', 0.0)
    d.set_quantity('stage', np.where(x < DAM_X, HL, HR), location='centroids')
    d.set_quantity('friction', 0.0)
    d.set_boundary({t: Reflective_boundary(d) for t in d.get_boundary_tags()})
    d.add_sediment_fraction(name='ssc', diameter=1e-4, tau_c_star=0.0, d_star=0.0,
                     initial_concentration=np.where(x < DAM_X, CL, CR))
    d.evolve_to_end(finaltime=DAM_T)
    h = np.maximum(d.quantities['stage'].centroid_values
                   - d.quantities['elevation'].centroid_values, 0.0)
    return x, h, d.get_tracer('ssc'), LENGTH / nx


@pytest.fixture(scope='module')
def dambreak():
    out = {}
    for nx in (400, 800):
        x, h_num, c_num, dx = run_dambreak(nx)
        h_ex, um, hm, shock_s = stoker(HL, HR, DAM_T, x)
        c_ex = np.where(x < DAM_X + um * DAM_T, CL, CR)
        out[nx] = dict(dx=dx, x=x, h=h_num, c=c_num, um=um, hm=hm,
                       shock_s=shock_s,
                       L1h=float(np.abs(h_num - h_ex).mean()),
                       L1c=float(np.abs(c_num - c_ex).mean()))
    return out


@pytest.mark.slow
def test_depth_agrees_with_the_analytical_stoker_solution(dambreak):
    assert dambreak[800]['L1h'] < 5e-3


@pytest.mark.slow
def test_ssc_agrees_with_the_passive_contact_solution(dambreak):
    assert dambreak[800]['L1c'] < 5e-3


@pytest.mark.slow
def test_both_converge_under_refinement(dambreak):
    """L1 and convergence, not a max norm: the solution is discontinuous, so a
    max norm measures only how the scheme smears the shock."""
    assert dambreak[800]['L1h'] < dambreak[400]['L1h']
    assert dambreak[800]['L1c'] < dambreak[400]['L1c']


@pytest.mark.slow
def test_ssc_is_flat_across_the_rarefaction(dambreak):
    """Sediment rides the FLOW, not the water waves: a nonlinear rarefaction in
    h must leave no signature at all in c."""
    r = dambreak[800]
    fan = ((r['x'] > DAM_X - np.sqrt(G * HL) * DAM_T * 0.9)
           & (r['x'] < DAM_X + (r['um'] - np.sqrt(G * r['hm'])) * DAM_T * 0.9))
    assert fan.sum() > 10, 'the fan was not resolved, so this proves nothing'
    assert np.ptp(r['c'][fan]) < 1e-3


@pytest.mark.slow
def test_ssc_is_flat_between_the_contact_and_the_shock(dambreak):
    r = dambreak[800]
    zone = ((r['x'] > DAM_X + r['um'] * DAM_T * 1.1)
            & (r['x'] < DAM_X + r['shock_s'] * DAM_T * 0.9))
    assert zone.sum() > 10
    assert np.ptp(r['c'][zone]) < 1e-3


@pytest.mark.slow
def test_ssc_develops_no_new_extrema(dambreak):
    r = dambreak[800]
    assert r['c'].min() >= CR - 1e-9
    assert r['c'].max() <= CL + 1e-9

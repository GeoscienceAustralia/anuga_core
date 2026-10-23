"""The bed shear the sediment kernel sees in vegetated cells.

Uniform flow down a vegetated plane at the Baptist normal state (exact and
steady), with Smith-McLean entrainment: the concentration gained in a short
time is E T / h, and E follows from the friction factor the option predicts:

  'bed'    f_c = g (Cv_r / Cv)^2 / Cb^2   (the bed's share, Baptist et al.)
  'total'  f_c = g / Cv^2                 (the whole vegetated resistance)
  'ignore' the cell's own closure: n = 0 there, so no shear at all.
"""
import numpy as np
import pytest

from anuga import Dirichlet_boundary, Reflective_boundary, rectangular_cross_domain

G = 9.8
CD, CB = 1.68, 65.0
SLOPE, H0 = 1.0e-3, 1.0


def chezy(h, m, D, hv):
    Cv_r = (CB ** -2 + CD * m * D * min(h, hv) / (2 * G)) ** -0.5
    Cv = Cv_r + (np.sqrt(G) / 0.4 * np.log(h / hv) if h > hv else 0.0)
    return Cv_r, Cv


def f_c(option, h, m, D, hv):
    Cv_r, Cv = chezy(h, m, D, hv)
    return {'bed': G * (Cv_r / Cv) ** 2 / CB ** 2, 'total': G / Cv ** 2, 'ignore': 0.0}[option]


def channel(m, D, hv, option, mode='legacy', diameter=1.0e-4):
    L, W = 200.0, 10.0
    d = rectangular_cross_domain(40, 2, len1=L, len2=W)
    d.set_flow_algorithm('DE1')
    d.set_compute_mode(mode)
    d.store = False
    _, Cv = chezy(H0, m, D, hv)
    u = Cv * np.sqrt(H0 * SLOPE)
    d.set_quantity('elevation', lambda x, y: SLOPE * (L - x))
    d.set_quantity('friction', 0.0)
    d.set_quantity('stage', lambda x, y: SLOPE * (L - x) + H0)
    d.set_quantity('xmomentum', u * H0)
    d.set_boundary({'left': Dirichlet_boundary([SLOPE * L + H0, u * H0, 0.0]),
                    'right': Dirichlet_boundary([H0, u * H0, 0.0]),
                    'top': Reflective_boundary(d), 'bottom': Reflective_boundary(d)})
    d.set_vegetation_drag(density=m, diameter=D, height=hv, sediment_shear=option)
    d.initialize_sediment_operator(bed_evolution=False)
    d.set_deposition(law='d_star', near_bed='constant')
    d.add_sediment_fraction('silt', diameter=diameter, d_star=1.0, tau_c_star=0.04,
                            initial_concentration=0.0)
    return d, u


def smith_mclean(d, fc, u):
    v_s = float(d.sediment_settling_velocity[0])
    R = float(d.sediment_R[0])
    diam = float(d.sediment_diameter[0])
    tau_star = fc * u * u / (R * G * diam)
    X = tau_star / 0.04 - 1.0
    if X <= 0.0:
        return 0.0, v_s
    g0 = float(d.sediment_gamma0)
    return v_s * 0.65 * g0 * X / (1.0 + g0 * X), v_s


def gained(d, T=0.5):
    d.evolve_to_end(finaltime=T)
    x = d.centroid_coordinates[:, 0]
    return d.get_tracer('silt')[(x > 40.0) & (x < 160.0)]


@pytest.mark.parametrize('stems', [(5.0, 0.01, 10.0),      # emergent, sparse
                                   (20.0, 0.01, 0.5)])     # submerged
@pytest.mark.parametrize('option', ['bed', 'total'])
@pytest.mark.parametrize('mode', ['legacy', 'unified'])
def test_the_kernel_sees_the_vegetated_shear(stems, option, mode):
    m, D, hv = stems
    d, u = channel(m, D, hv, option, mode)
    E, v_s = smith_mclean(d, f_c(option, H0, m, D, hv), u)
    assert E > 0.0
    T = 0.5
    c = gained(d, T)
    expected = E * T * (1.0 - 0.5 * v_s * T / H0) / H0
    assert np.allclose(c, expected, rtol=5e-3)


def test_ignore_keeps_the_old_zero_shear():
    """n = 0 in a vegetated class, and the kernel told to ignore the stems:
    no shear, so nothing is entrained (the behaviour before this option)."""
    d, _ = channel(120.0, 0.01, 10.0, 'ignore')
    assert np.all(gained(d) == 0.0)


def test_total_exceeds_the_bed_share():
    for stems in ((120.0, 0.01, 10.0), (20.0, 0.01, 0.5)):
        assert f_c('total', H0, *stems) > f_c('bed', H0, *stems)


def test_emergent_bed_share_is_the_bed_chezy():
    """Emergent stems: the canopy velocity is the mean velocity, so the bed
    feels exactly its own Chezy friction, whatever the stems."""
    assert f_c('bed', H0, 120.0, 0.01, 10.0) == pytest.approx(G / CB ** 2, rel=1e-12)


def test_unvegetated_cells_keep_their_manning_friction():
    """Stems on half the channel: the bare half entrains at the Manning
    friction factor, as before."""
    d, u = channel(0.0, 0.01, 10.0, 'bed')
    d.set_quantity('friction', 0.02)
    x = d.centroid_coordinates[:, 0]
    d.set_vegetation_drag(density=np.where(x > 100.0, 120.0, 0.0), diameter=0.01,
                          height=10.0, sediment_shear='bed')
    E_bare, v_s = smith_mclean(d, G * 0.02 ** 2 / H0 ** (1.0 / 3.0), u)
    d.evolve_to_end(finaltime=0.2)
    c = d.get_tracer('silt')
    bare = (x > 20.0) & (x < 60.0)
    assert np.allclose(c[bare], E_bare * 0.2 * (1 - 0.1 * v_s / H0) / H0, rtol=2e-2)


def test_a_bad_option_is_rejected():
    d = rectangular_cross_domain(2, 2)
    with pytest.raises(ValueError):
        d.set_vegetation_drag(1.0, 0.01, 1.0, sediment_shear='canopy')

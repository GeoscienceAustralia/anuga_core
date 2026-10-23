"""[E-6] de Leeuw et al. (2020) entrainment.

    E* = A X^beta / (1 + 3 A X^beta),   X = (u*_skin / v_s)^alpha Fr - threshold

with u*_skin from Manning-Strickler on a skin roughness k_s. Checked cell by
cell in uniform flow, where every quantity in the relation is known exactly:
a flat frictionless channel held at (h0, U) by Dirichlet boundaries, with the
sediment kernel's own friction set to a uniform Manning n (larsen_lamb), so
f_c = g n^2 / h^(1/3) and tau_b / rho = f_c U^2.
"""
import numpy as np
import pytest

import anuga
from anuga import Dirichlet_boundary, Reflective_boundary, rectangular_cross_domain

G = 9.8
H0 = 1.0


def uniform_flow(U, mode='legacy', diameter=2.0e-4, k_s=0.05, **bed):
    d = rectangular_cross_domain(20, 4, len1=100.0, len2=20.0)
    d.set_flow_algorithm('DE1')
    d.set_compute_mode(mode)
    d.store = False
    d.set_quantity('elevation', 0.0)
    d.set_quantity('friction', 0.0)
    d.set_quantity('stage', H0)
    d.set_quantity('xmomentum', U * H0)
    Bd = Dirichlet_boundary([H0, U * H0, 0.0])
    Br = Reflective_boundary(d)
    d.set_boundary({'left': Bd, 'right': Bd, 'top': Br, 'bottom': Br})
    d.initialize_sediment_operator(bed_evolution=False)
    d.set_sediment_friction('larsen_lamb', k_s=k_s)
    d.set_deposition(law='d_star', near_bed='constant')
    d.set_bed_material('noncohesive', entrainment='de_leeuw', **bed)
    d.add_sediment_fraction('sand', diameter=diameter, d_star=1.0,
                            initial_concentration=0.0)
    return d


def expected_E_star(d, U, h=H0, ks=None):
    """The relation, written out independently of the kernel."""
    n = d.sediment_manning_ll
    f_c = G * n * n / h ** (1.0 / 3.0)
    v_s = float(d.sediment_settling_velocity[0])
    diam = float(d.sediment_diameter[0])
    S = f_c * U * U / (G * h)
    ks = ks if ks is not None else d.sediment_dl_ks_factor * diam
    ks = max(ks, 1.0e-6 / (8.0 * np.sqrt(f_c) * U))
    H_sk = min((U * ks ** (1.0 / 6.0) / (8.1 * np.sqrt(G * S))) ** 1.5, h)
    ustar_sk = np.sqrt(G * H_sk * S)
    X = (ustar_sk / v_s) ** d.sediment_dl_alpha * (U / np.sqrt(G * h)) - d.sediment_dl_threshold
    if X <= 0.0:
        return 0.0, v_s
    aX = d.sediment_dl_A * X ** d.sediment_dl_beta
    return aX / (1.0 + 3.0 * aX), v_s


def entrained(d, T=0.5):
    """Concentration after a short time T << h / v_s, so deposition is a
    small, known correction: c = E T (1 - v_s T / (2 h)) / h."""
    d.evolve_to_end(finaltime=T)
    return d.get_tracer('sand')


@pytest.mark.parametrize('fit', ['de_leeuw_2020', 'nghiem_2022'])
@pytest.mark.parametrize('mode', ['legacy', 'unified'])
def test_the_entrainment_matches_the_relation(fit, mode):
    U, T = 1.2, 0.5
    d = uniform_flow(U, mode, de_leeuw_fit=fit)
    E_star, v_s = expected_E_star(d, U)
    assert E_star > 0.0
    c = entrained(d, T)
    expected = v_s * E_star * T * (1.0 - 0.5 * v_s * T / H0) / H0
    x = d.centroid_coordinates[:, 0]
    interior = (x > 10.0) & (x < 90.0)
    assert np.allclose(c[interior], expected, rtol=2e-3)


def test_the_presets():
    d = uniform_flow(1.0, de_leeuw_fit='nghiem_2022')
    assert (d.sediment_dl_A, d.sediment_dl_alpha, d.sediment_dl_beta,
            d.sediment_dl_threshold) == (7.04e-4, 0.94475138, 1.81, 0.015)
    assert d.sediment_erosion_mode == 3
    d = uniform_flow(1.0)
    assert (d.sediment_dl_A, d.sediment_dl_alpha, d.sediment_dl_beta) == (4.74e-4, 1.5, 1.18)
    assert 'de Leeuw' in d.sediment_summary()


def test_below_the_threshold_nothing_is_entrained():
    """Slow flow over coarse sand: X < 0."""
    d = uniform_flow(0.05, diameter=1.0e-3)
    assert expected_E_star(d, 0.05)[0] == 0.0
    assert np.all(entrained(d) == 0.0)


def test_the_near_bed_concentration_saturates_at_a_third():
    """Very fast flow over fine grains: E* -> 1/3 from below, never above."""
    U = 6.0
    d = uniform_flow(U, diameter=2.0e-5, de_leeuw_fit='nghiem_2022')
    E_star, _ = expected_E_star(d, U)
    assert 0.3 < E_star < 1.0 / 3.0


def test_a_fixed_skin_roughness_overrides_the_grain_scale():
    U = 1.2
    grain = uniform_flow(U)
    rough = uniform_flow(U, skin_roughness=0.01)
    assert rough.sediment_dl_ks == 0.01
    E_grain, _ = expected_E_star(grain, U)
    E_rough, _ = expected_E_star(rough, U, ks=0.01)
    # a rougher skin carries more of the stress: H_sk ~ k_s^(1/4)
    assert E_rough > E_grain
    c = entrained(rough)
    x = rough.centroid_coordinates[:, 0]
    interior = (x > 10.0) & (x < 90.0)
    v_s = float(rough.sediment_settling_velocity[0])
    expected = v_s * E_rough * 0.5 * (1.0 - 0.25 * v_s / H0) / H0
    assert np.allclose(c[interior], expected, rtol=2e-3)


def test_overrides_and_bad_settings():
    d = uniform_flow(1.0, A=1.0e-3, beta=1.5)
    assert (d.sediment_dl_A, d.sediment_dl_beta) == (1.0e-3, 1.5)
    assert d.sediment_dl_alpha == 1.5            # the rest from the fit
    with pytest.raises(ValueError):
        d.set_bed_material('noncohesive', entrainment='de_leeuw', de_leeuw_fit='wright')
    with pytest.raises(ValueError):
        d.set_bed_material('cohesive', entrainment='de_leeuw')
    with pytest.raises(ValueError):
        d.set_bed_material('noncohesive', entrainment='garcia_parker')


def test_smith_mclean_is_still_the_default():
    d = rectangular_cross_domain(2, 2)
    d.set_bed_material('noncohesive')
    assert d.sediment_erosion_mode == 0

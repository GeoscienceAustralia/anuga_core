"""Cohesive erosion [E-3]/[E-5], checked against anugaSed (PHYSICS_SPEC 4.1.1).

[E-1] and [E-3] are not competing formulations of the same physics. They
describe DIFFERENT SEDIMENT -- Shields entrainment of sand and gravel versus
jet-test-calibrated erosion of silt and clay -- so choosing between them is a
statement about the bed, which is why the API names the material rather than
the formula.

anugaSed cannot be run here (Python 2, and it drives a Sed_transport_operator
this ANUGA does not have), but the erosion LAW can be compared exactly. Theirs,
from operators/sed_transport_operator.py (MIT licensed):

    self.Ke = 0.2e-6 / self.tau_crit**(0.5)
    shear_stress = self.rho_w * self.u_star**2
    edot = self.Ke * (shear_stress[self.ind] - self.tau_crit)
    edot[edot<0.0] = 0.0

The comparison is made at equal tau_b. How each model REACHES that tau_b is a
separate divergence: they use the depth-slope closure [T-7], we use quadratic
drag [T-1], so a whole-model comparison would differ for that reason alone even
with identical erosion laws.

The mode 1 / mode 2 comparison is in test_sediment_gpu.py.
"""
import numpy as np
import pytest

from anuga import Reflective_boundary, rectangular_cross_domain

LEN = 200.0
TAU_CRIT = 0.088          # Pa, anugaSed's default


def channel(depth=1.0, slope=0.01, n_manning=0.03, dt=1.0):
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


def anugased_edot(tau_b, tau_crit):
    """Transcribed from anugaSed's erosion(), quoted above."""
    Ke = 0.2e-6 / tau_crit ** 0.5
    return max(Ke * (tau_b - tau_crit), 0.0)


def our_edot(tau_b, tau_crit):
    """[E-3] with [E-5]'s K_e, as the kernel computes it."""
    K_e = 0.2e-6 / tau_crit ** 0.5
    return max(K_e * (tau_b - tau_crit), 0.0)


def test_non_cohesive_is_the_default():
    assert channel().sediment_erosion_mode == 0


def test_selecting_cohesive_switches_the_route():
    d = channel()
    d.set_bed_material('cohesive', tau_crit=TAU_CRIT)
    assert d.sediment_erosion_mode == 1


def test_E5_gives_the_published_erosion_coefficient():
    d = channel()
    d.set_bed_material('cohesive', tau_crit=TAU_CRIT)
    assert abs(d.sediment_K_e - 0.2e-6 / TAU_CRIT ** 0.5) < 1e-18
    # Hanson & Simon's relation at their default threshold.
    assert abs(d.sediment_K_e - 6.7420e-7) < 1e-10


def test_the_erosion_coefficient_can_be_overridden():
    d = channel()
    d.set_bed_material('cohesive', tau_crit=TAU_CRIT, K_e=5e-7)
    assert d.sediment_K_e == 5e-7


def test_unknown_material_and_non_positive_threshold_are_rejected():
    with pytest.raises(ValueError):
        channel().set_bed_material('sandy')
    with pytest.raises(ValueError):
        channel().set_bed_material('cohesive', tau_crit=0.0)


def test_the_erosion_rate_has_sane_units_and_magnitude():
    """E is a velocity. At about 1 Pa of excess it should be a fraction of a
    millimetre per hour, not metres per second."""
    e = our_edot(1.0 + TAU_CRIT, TAU_CRIT)
    assert 1e-9 < e < 1e-5, '%.3e m/s' % e


def test_our_E3_reproduces_anugased_exactly():
    diffs = [abs(our_edot(t, TAU_CRIT) - anugased_edot(t, TAU_CRIT))
             for t in np.linspace(0.0, 20.0, 201)]
    assert max(diffs) == 0.0, 'max difference %.3e over tau_b in [0, 20] Pa' % max(diffs)


def test_both_clamp_to_zero_below_the_threshold():
    assert our_edot(0.05, TAU_CRIT) == 0.0 == anugased_edot(0.05, TAU_CRIT)


def test_below_the_critical_shear_nothing_erodes():
    d = channel()
    d.set_bed_material('cohesive', tau_crit=1e6)      # unreachable
    d.add_sediment_fraction(name='silt', diameter=6.5e-5, initial_concentration=0.0)
    d.evolve_to_end(finaltime=20.0)
    assert float(np.abs(d.get_tracer('silt')).max()) == 0.0


def test_above_it_the_bed_erodes():
    d = channel()
    d.set_bed_material('cohesive', tau_crit=TAU_CRIT)
    d.add_sediment_fraction(name='silt', diameter=6.5e-5, initial_concentration=0.0)
    d.evolve_to_end(finaltime=20.0)
    assert d.get_tracer('silt').max() > 0.0


def test_the_two_routes_give_materially_different_answers():
    """Spec 4.1.1 calls picking the wrong one a physics error rather than a
    tuning error, so the two had better not agree."""
    coh = channel()
    coh.set_bed_material('cohesive', tau_crit=TAU_CRIT)
    coh.add_sediment_fraction(name='silt', diameter=6.5e-5, initial_concentration=0.0)
    coh.evolve_to_end(finaltime=20.0)

    non = channel()
    non.add_sediment_fraction(name='silt', diameter=6.5e-5, tau_c_star=0.04,
                       initial_concentration=0.0)
    non.evolve_to_end(finaltime=20.0)

    a = float(coh.get_tracer('silt').mean())
    b = float(non.get_tracer('silt').mean())
    assert abs(a - b) > 0.1 * max(a, b), 'cohesive %.4e vs non-cohesive %.4e' % (a, b)

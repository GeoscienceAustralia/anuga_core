"""Friction closures for the sediment kernel (PHYSICS_SPEC 3.3).

Three closures feed f_c into the bed shear: 'constant' from the domain's
Manning n [T-6], 'larsen_lamb' [T-14]/[T-15], and 'wilson' [T-8]..[T-10]. They
affect only the sediment source term; the hydrodynamic friction operator is
untouched.

The mode 1 / mode 2 comparison is in test_sediment_gpu.py.
"""
import math

import numpy as np
import pytest

from anuga import Reflective_boundary, rectangular_cross_domain

LEN = 500.0

# Wilson (2004) abstract: n = 0.0545 s m^-1/3 for Martian outflow channels.
G_MARS = 3.71
X_W04 = 8.46 * 1000.0 ** 0.1005


def channel():
    d = rectangular_cross_domain(20, 10, len1=LEN, len2=LEN / 2)
    d.set_flow_algorithm('DE0')
    d.set_low_froude(0)
    d.store = False
    d.set_quantity('elevation', lambda x, y: -0.01 * x)
    d.set_quantity('stage', lambda x, y: -0.01 * x + 1.0)
    d.set_quantity('friction', 0.03)
    d.set_boundary({t: Reflective_boundary(d) for t in d.get_boundary_tags()})
    d.evolve_max_timestep = 1.0
    return d


def _depth_for_n(f_c, n, g):
    """Invert [T-6]: n = sqrt(f_c h^(1/3) / g)."""
    return (n * n * g / f_c) ** 3


def test_wilson_warns_when_grain_size_is_implausible_for_the_bed():
    """`bed` picks the curve, `grain_size` sets h/D -- nothing ties them.

    A boulder curve at a sand-sized D is a well-formed calculation of nothing
    real: f_c drops from 0.031 to 0.0016, a factor of 19 in tau_b, silently.
    """
    import warnings

    for bed, D in (('boulder', 2.0e-4), ('boulder', 0.02),
                   ('gravel', 2.0e-4), ('sand', 0.5)):
        d = channel()
        d.add_sediment_fraction('sand', diameter=2.0e-4)
        with pytest.warns(UserWarning, match='implausible'):
            d.set_sediment_friction('wilson', bed=bed, grain_size=D)


def test_wilson_is_quiet_for_sensible_bed_and_grain_size():
    """The warning must not fire on ordinary choices, or it will be ignored."""
    import warnings

    for bed, D in (('sand', 2.0e-4), ('sand', 1.0e-3),
                   ('gravel', 0.01), ('gravel', 0.1),
                   ('boulder', 0.3), ('boulder', 1.0)):
        d = channel()
        d.add_sediment_fraction('sand', diameter=2.0e-4)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            d.set_sediment_friction('wilson', bed=bed, grain_size=D)
        assert not [x for x in w if 'implausible' in str(x.message)], (
            'warned on a sensible combination: bed=%s grain_size=%g' % (bed, D))


def test_the_wilson_factor_of_eight_conversion():
    """W04 write (8/f_c)^1/2, but THEIR f_c is the Darcy-Weisbach f (their
    Eq 4) while ours is f/8. Taking their f_c literally as ours is an
    eight-fold error, and it is silent -- the code runs and gives plausible
    numbers. Anchored to W04's own headline result rather than to arithmetic.
    """
    f_c_ours = 1.0 / (X_W04 * X_W04)
    f_c_literal = 8.0 / (X_W04 * X_W04)

    h_ok = _depth_for_n(f_c_ours, 0.0545, G_MARS)
    h_bad = _depth_for_n(f_c_literal, 0.0545, G_MARS)

    # Martian outflow channels are tens of metres deep.
    assert 5.0 < h_ok < 200.0, 'correct conversion gives h = %.1f m' % h_ok
    assert h_bad < 1.0, 'the literal reading gives h = %.3f m' % h_bad
    assert abs(f_c_literal / f_c_ours - 8.0) < 1e-12


def test_larsen_lamb_reproduces_the_published_manning_n():
    """LL16 report n = 0.065 at Moses Coulee, with sigma_br about 5 m."""
    d = channel()
    d.set_sediment_friction('larsen_lamb', sigma_br=5.0)
    assert abs(d.sediment_manning_ll - 0.065) < 5e-4


def test_larsen_lamb_accepts_k_s_directly():
    a = channel()
    a.set_sediment_friction('larsen_lamb', sigma_br=5.0)   # k_s = 2*2*5 = 20 m
    b = channel()
    b.set_sediment_friction('larsen_lamb', k_s=20.0)
    assert abs(a.sediment_manning_ll - b.sediment_manning_ll) < 1e-12


def test_larsen_lamb_scales_as_the_sixth_root_of_roughness():
    a = channel()
    a.set_sediment_friction('larsen_lamb', k_s=20.0)
    b = channel()
    b.set_sediment_friction('larsen_lamb', k_s=40.0)
    assert b.sediment_manning_ll > a.sediment_manning_ll
    assert abs(b.sediment_manning_ll / a.sediment_manning_ll
               - 2 ** (1 / 6)) < 1e-12


def test_larsen_lamb_requires_a_length_scale():
    """sigma_br is site-measured and has no universal default, so it is refused
    rather than guessed."""
    with pytest.raises(ValueError):
        channel().set_sediment_friction('larsen_lamb')


def _fc_wilson(rel, bed):
    """[T-8]..[T-10] as compiled into the kernel: f_c = 1/X^2."""
    rel = max(rel, 1.0)
    if bed == 'sand':
        X = 8.46 * rel ** 0.1005
    elif bed == 'gravel':
        X = 5.75 * math.log10(rel) + 3.514
    else:
        X = 5.62 * math.log10(rel) + 4.0
    return 1.0 / (X * X)


def test_wilson_f_c_falls_with_relative_submergence():
    assert _fc_wilson(1000, 'sand') < _fc_wilson(10, 'sand')


def test_wilson_bed_types_are_distinct():
    values = {round(_fc_wilson(100, b), 8)
              for b in ('sand', 'gravel', 'boulder')}
    assert len(values) == 3


def test_wilson_submergence_is_floored_at_one():
    """h < D would otherwise send the logarithm negative and f_c to nonsense."""
    assert np.isfinite(_fc_wilson(0.001, 'gravel'))
    assert _fc_wilson(0.001, 'gravel') == _fc_wilson(1.0, 'gravel')


def test_wilson_requires_a_grain_size():
    with pytest.raises(ValueError):
        channel().set_sediment_friction('wilson', bed='gravel')


def test_unknown_modes_and_beds_are_rejected():
    with pytest.raises(ValueError):
        channel().set_sediment_friction('nope')
    with pytest.raises(ValueError):
        channel().set_sediment_friction('wilson', bed='mud', grain_size=1e-3)


def test_constant_is_the_default():
    assert channel().sediment_friction_mode == 0


def test_all_three_closures_run_and_give_distinct_answers():
    """Distinct matters: a closure that silently fell back to 'constant' would
    otherwise pass every check above."""
    means = {}
    for name, kwargs in (('constant', {}),
                         ('larsen_lamb', dict(sigma_br=5.0)),
                         ('wilson', dict(bed='gravel', grain_size=0.05))):
        d = channel()
        if name != 'constant':
            d.set_sediment_friction(name, **kwargs)
        d.sediment_d_star_mode = 1
        d.add_sediment_fraction(name='s', diameter=1e-4, initial_concentration=0.01)
        d.evolve_to_end(finaltime=20.0)
        means[name] = float(d.get_tracer('s').mean())

    assert all(np.isfinite(v) for v in means.values()), means
    assert len({round(v, 12) for v in means.values()}) == 3, means

"""Bedload transport [K-1]-[K-4] and its bed evolution [G-5] (PHYSICS_SPEC 6).

Bedload is a DIVERGENCE, not a source: it moves sediment ALONG the bed rather
than between bed and water column. So its defining property is that in a closed
domain it redistributes bed material and conserves total bed volume EXACTLY.
That is what separates a correct divergence from a plausible one, and it is why
the flux is computed per shared edge rather than from a reconstructed gradient.

The mode 1 / mode 2 comparison is in test_sediment_gpu.py.
"""
import numpy as np
import pytest

from anuga import Reflective_boundary, rectangular_cross_domain

LEN = 100.0


def channel(slope=0.02, depth=0.5, n_manning=0.025, dt=0.5):
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


def _gravel_channel(formula='wong_parker_eq24', concentration=0.0):
    d = channel()
    d.add_sediment_fraction(name='gravel', diameter=5e-3, tau_c_star=0.0,
                     initial_concentration=concentration)
    d.set_bedload(formula)
    return d


@pytest.mark.parametrize('formula,expected', [
    ('wong_parker_eq24', (3.97, 1.5, 0.0495)),
    ('wong_parker_eq23', (4.93, 1.60, 0.0470)),
])
def test_the_published_parameter_sets(formula, expected):
    """Wong & Parker's corrected Meyer-Peter-Muller fits. Which of the two FG21
    used is not settled, so both are available."""
    d = _gravel_channel(formula)
    assert (d.sediment_bedload_K, d.sediment_bedload_m,
            d.sediment_bedload_tau_c_star) == expected


def test_an_unknown_formula_is_rejected():
    with pytest.raises(ValueError):
        channel().set_bedload('mpm')


@pytest.fixture(scope='module')
def transported():
    d = _gravel_channel()
    z0 = d.quantities['elevation'].centroid_values.copy()
    vol0 = float((z0 * d.areas).sum())
    d.evolve_to_end(finaltime=60.0)
    z1 = d.quantities['elevation'].centroid_values
    return d, z0, z1, vol0, float((z1 * d.areas).sum())


def test_bedload_moves_the_bed(transported):
    _, z0, z1, _, _ = transported
    assert np.abs(z1 - z0).max() > 0.0


def test_total_bed_volume_is_conserved_exactly(transported):
    """The property that defines a divergence. FG21 report ~1% loss in their
    bedload operator, which they could not correct; computing the flux as a
    difference across shared edges makes conservation structural instead."""
    _, _, _, vol0, vol1 = transported
    assert abs(vol1 - vol0) <= 1e-10 * max(abs(vol0), 1.0)


def test_it_redistributes_rather_than_only_eroding(transported):
    _, z0, z1, _, _ = transported
    dz = z1 - z0
    assert dz.min() < 0.0 < dz.max()


def test_it_erodes_upstream_and_builds_downstream(transported):
    d, z0, z1, _, _ = transported
    dz = z1 - z0
    x = d.centroid_coordinates[:, 0]
    up, dn = x < LEN * 0.25, x > LEN * 0.75
    assert dz[up].sum() < 0.0 < dz[dn].sum()


def test_engelund_hansen_disables_the_suspended_source():
    """[K-5] is TOTAL load: it already contains suspension, so running the
    suspended source alongside it double counts (spec 6's usage rule)."""
    d = _gravel_channel('engelund_hansen', concentration=0.01)
    assert d._sediment_suspended_enabled is False
    assert d.sediment_bedload_tau_c_star == 0.0, 'total load has no threshold'


def test_under_engelund_hansen_suspended_mass_is_advected_not_exchanged():
    """With the source off, m is advected but never exchanged, so in a closed
    domain its total is conserved. Were the suspended operator still running
    alongside [K-5] -- the double counting spec 6 forbids -- this would change.
    """
    d = _gravel_channel('engelund_hansen', concentration=0.01)
    m0 = float((d.tracer_conserved_values[0] * d.areas).sum())
    d.evolve_to_end(finaltime=20.0)
    m1 = float((d.tracer_conserved_values[0] * d.areas).sum())
    assert abs(m1 - m0) <= 1e-9 * max(abs(m0), 1.0)
    assert np.isfinite(d.quantities['elevation'].centroid_values).all()


# ---------------------------------------------------------------------------
# Grass [K-6] and open boundaries
# ---------------------------------------------------------------------------

def _uniform_channel(open_boundaries=None, K=0.01):
    """Uniform frictionless flow, 10 m deep at 1 m/s, held by Dirichlet
    boundaries at both ends: under Grass every cell carries the same q_b, so
    the bed must not move -- except at a closed inflow, which exports and
    never imports."""
    from anuga import Dirichlet_boundary
    d = rectangular_cross_domain(20, 2, len1=100.0, len2=10.0)
    d.set_flow_algorithm('DE1')
    d.store = False
    d.set_quantity('elevation', 0.1)
    d.set_quantity('friction', 0.0)
    d.set_quantity('stage', 10.0)
    d.set_quantity('xmomentum', 10.0)
    Bd = Dirichlet_boundary([10.0, 10.0, 0.0])
    Br = Reflective_boundary(d)
    d.set_boundary({'left': Bd, 'right': Bd, 'top': Br, 'bottom': Br})
    d.initialize_sediment_operator(porosity=0.4, bed_evolution=True)
    d.add_sediment_fraction(name='sand', diameter=1e-3,
                            initial_concentration=0.0)
    d.set_bedload('grass', K=K, open_boundaries=open_boundaries)
    return d


def test_grass_needs_its_coefficient():
    """A_g is a calibration; there is no defensible default."""
    with pytest.raises(ValueError):
        channel().set_bedload('grass')


def test_grass_is_total_load_with_no_threshold():
    d = _gravel_channel()
    d.set_bedload('grass', K=0.001)
    assert d.sediment_bedload_mode == 3
    assert d.sediment_bedload_m == 3.0
    assert d.sediment_bedload_tau_c_star == 0.0
    assert d._sediment_suspended_enabled is False
    d.set_bedload('grass', K=0.001, m=2.5, tau_c_star=0.1)
    assert d.sediment_bedload_m == 2.5
    assert d.sediment_bedload_tau_c_star == 0.0, 'Grass has no threshold'


def test_grass_transport_is_A_g_times_speed_cubed():
    """[K-6] q_b = A_g |u|^m along the flow, read from the transport vector
    the kernel leaves behind (mode 1 keeps it on the host)."""
    d = _uniform_channel(open_boundaries=('left', 'right'), K=0.01)
    d.set_compute_mode('legacy')
    d.evolve_to_end(finaltime=1.0)
    h = (d.quantities['stage'].centroid_values
         - d.quantities['elevation'].centroid_values)
    u = d.quantities['xmomentum'].centroid_values / h
    v = d.quantities['ymomentum'].centroid_values / h
    speed = np.hypot(u, v)
    q = 0.01 * speed ** 3
    assert np.allclose(d.sediment_qbx, q * u / speed, rtol=1e-6, atol=1e-12)
    assert np.allclose(d.sediment_qby, q * v / speed, rtol=1e-6, atol=1e-12)


def test_an_unknown_open_boundary_tag_is_rejected():
    with pytest.raises(ValueError):
        _uniform_channel(open_boundaries=('inflow',))


def test_open_boundaries_flag_exactly_the_tagged_edges():
    d = _uniform_channel(open_boundaries=('left',))
    flags = d.sediment_bedload_open
    assert flags.shape == (d.boundary_length,)
    left = np.asarray(d.tag_boundary_cells['left'])
    assert flags.sum() == left.size
    assert (flags[left] == 1).all()
    # Re-setting replaces rather than accumulates
    d.set_bedload('grass', K=0.01, open_boundaries=('right',))
    right = np.asarray(d.tag_boundary_cells['right'])
    assert d.sediment_bedload_open.sum() == right.size
    assert (d.sediment_bedload_open[right] == 1).all()


def test_closed_inflow_digs_a_hole_and_open_inflow_does_not():
    """With every boundary closed the first cells export bedload they never
    receive; declaring the inflow and outflow open makes the uniform bed an
    exact steady state of the divergence."""
    results = {}
    for tags in (None, ('left', 'right')):
        d = _uniform_channel(open_boundaries=tags)
        z0 = d.quantities['elevation'].centroid_values.copy()
        d.evolve_to_end(finaltime=60.0)
        dz = d.quantities['elevation'].centroid_values - z0
        results[tags] = dz
    x = _uniform_channel().centroid_coordinates[:, 0]
    closed = results[None]
    assert closed[x < 10.0].min() < -1e-4, 'the closed inflow should erode'
    assert np.abs(results[('left', 'right')]).max() < 1e-6



# ---------------------------------------------------------------- supply

def test_a_zero_supply_is_a_closed_boundary():
    """supply = 0 across an open edge admits nothing: bit-identical to the
    edge being closed."""
    beds = []
    for kw in (dict(), dict(open_boundaries=['left'], supply={'left': 0.0})):
        d = channel()
        d.initialize_sediment_operator(bed_evolution=True)
        d.set_deposition(tau_d=0.0, law='threshold')
        d.add_sediment_fraction('sand', diameter=1.0e-3, tau_c_star=0.0,
                                initial_concentration=0.0)
        d.set_bedload('wong_parker_eq24', **kw)
        for _ in d.evolve(yieldstep=2.0, finaltime=2.0):
            pass
        beds.append(d.quantities['elevation'].centroid_values.copy())
    assert np.array_equal(beds[0], beds[1])


def test_a_prescribed_supply_enters_at_exactly_that_rate():
    """Still water, so the bed carries no bedload of its own: the cells on
    the supplied boundary must rise by S * edge length * dt / (area (1-lambda))
    per step and nothing else may move."""
    d = rectangular_cross_domain(6, 4, len1=LEN, len2=LEN / 2)
    d.set_flow_algorithm('DE0')
    d.store = False
    d.set_quantity('elevation', 0.0)
    d.set_quantity('stage', 1.0)
    d.set_boundary({t: Reflective_boundary(d) for t in d.get_boundary_tags()})
    lam = 0.3
    d.initialize_sediment_operator(porosity=lam, bed_evolution=True)
    d.set_deposition(tau_d=0.0, law='threshold')
    d.add_sediment_fraction('sand', diameter=1.0e-3, initial_concentration=0.0)
    S = 2.0e-6                                  # m2/s per unit width
    d.set_bedload('wong_parker_eq24', supply={'left': S})
    assert d.sediment_bedload_open[d.tag_boundary_cells['left']].all()
    dt = 0.25
    d.evolve_max_timestep = dt
    z0 = d.quantities['elevation'].centroid_values.copy()
    for _ in d.evolve(yieldstep=dt, finaltime=dt):
        pass
    dz = d.quantities['elevation'].centroid_values - z0
    expected = np.zeros_like(dz)
    for b in d.tag_boundary_cells['left']:
        k, i = int(d.boundary_cells[b]), int(d.boundary_edges[b])
        expected[k] += S * d.edgelengths[k, i] * dt / (d.areas[k] * (1.0 - lam))
    assert expected.max() > 0.0
    assert np.allclose(dz, expected, rtol=1e-12, atol=1e-18)


def test_supply_is_validated_and_opens_its_tag():
    d = channel()
    d.initialize_sediment_operator(bed_evolution=True)
    d.add_sediment_fraction('sand', diameter=1.0e-3)
    with pytest.raises(ValueError):
        d.set_bedload('wong_parker_eq24', supply={'left': -1.0e-6})
    with pytest.raises(ValueError):
        d.set_bedload('wong_parker_eq24', supply={'nowhere': 1.0e-6})
    d.set_bedload('wong_parker_eq24', supply={'left': 1.0e-6})
    assert 'left' in d._sediment_bedload_open_tags
    left = d.tag_boundary_cells['left']
    assert (d.sediment_bedload_supply[left] == 1.0e-6).all()
    others = np.ones(d.boundary_length, bool); others[left] = False
    assert (d.sediment_bedload_supply[others] < 0.0).all()


def test_min_depth_ramps_the_bedload_off_in_thin_flow():
    """[K-7] no bedload below min_depth, full above 2 min_depth, linear in
    between. A closed-inflow channel digs a hole at its first cells whose
    depth after one step scales with the transport, so the ratio of the
    ramped to the plain bed change is the ramp factor for the channel's
    (uniform) depth."""
    def hole(depth, min_depth):
        from anuga import Dirichlet_boundary
        d = rectangular_cross_domain(20, 2, len1=100.0, len2=10.0)
        d.set_flow_algorithm('DE1')
        d.store = False
        d.set_quantity('elevation', 0.0)
        d.set_quantity('friction', 0.0)
        d.set_quantity('stage', depth)
        d.set_quantity('xmomentum', depth)          # 1 m/s
        Bd = Dirichlet_boundary([depth, depth, 0.0])
        Br = Reflective_boundary(d)
        d.set_boundary({'left': Bd, 'right': Bd, 'top': Br, 'bottom': Br})
        d.initialize_sediment_operator(porosity=0.4, bed_evolution=True)
        d.add_sediment_fraction(name='sand', diameter=1e-3, initial_concentration=0.0)
        d.set_deposition(law='threshold', tau_d=0.0)
        d.set_bedload('grass', K=0.01, min_depth=min_depth)
        # one step: the bed change is then linear in q_b (rk2's second
        # stage sees a state moved by O(dt), so the ratio holds to ~1e-3)
        d.evolve_max_timestep = 0.05
        for _ in d.evolve(yieldstep=0.05, finaltime=0.05):
            pass
        return d.quantities['elevation'].centroid_values
    h_min = 4.0
    for depth, factor in [(3.0, 0.0), (6.0, 0.5), (10.0, 1.0)]:
        plain = hole(depth, 0.0)
        ramped = hole(depth, h_min)
        assert plain.min() < -1e-6
        if factor == 0.0:
            assert np.abs(ramped).max() == 0.0
        else:
            assert np.allclose(ramped, factor * plain, rtol=1e-2, atol=1e-10)
    d = channel()
    with pytest.raises(ValueError):
        d.set_bedload('wong_parker_eq24', min_depth=-1.0)


@pytest.mark.parametrize('mode', ['legacy', 'unified'])
def test_shear_amplification_scales_f_c_like_the_square_of_manning_n(mode):
    """[T-16] a factor F on the bed shear is f_c -> F f_c, which under Manning
    is n -> sqrt(F) n: the bedload bed change of one step must match."""
    def hole(n_manning, factor):
        from anuga import Dirichlet_boundary
        d = rectangular_cross_domain(20, 2, len1=100.0, len2=10.0)
        d.set_flow_algorithm('DE1')
        d.set_compute_mode(mode)
        d.store = False
        d.set_quantity('elevation', 0.0)
        d.set_quantity('friction', 0.0)
        d.set_quantity('stage', 5.0)
        d.set_quantity('xmomentum', 5.0)
        Bd = Dirichlet_boundary([5.0, 5.0, 0.0])
        Br = Reflective_boundary(d)
        d.set_boundary({'left': Bd, 'right': Bd, 'top': Br, 'bottom': Br})
        d.initialize_sediment_operator(porosity=0.4, bed_evolution=True)
        d.set_sediment_friction('larsen_lamb', k_s=0.05)
        d.sediment_manning_ll = n_manning
        d.set_deposition(law='threshold', tau_d=0.0)
        d.add_sediment_fraction(name='sand', diameter=1e-3, tau_c_star=1e9,
                                initial_concentration=0.0)
        d.set_bedload('wong_parker_eq24', tau_c_star=0.0)
        if factor is not None:
            d.set_shear_amplification(factor)
        d.evolve_max_timestep = 0.05
        for _ in d.evolve(yieldstep=0.05, finaltime=0.05):
            pass
        return d.quantities['elevation'].centroid_values
    n0 = 0.02
    plain = hole(n0, None)
    amplified = hole(n0, 4.0)
    equivalent = hole(2.0 * n0, None)
    assert plain.min() < -1e-9
    assert np.allclose(amplified, equivalent, rtol=1e-9, atol=1e-14)
    assert not np.allclose(amplified, plain, rtol=1e-3, atol=1e-14)
    # a field, and off again
    d = channel()
    d.add_sediment_fraction(name='gravel', diameter=5e-3, initial_concentration=0.0)
    d.set_shear_amplification(lambda x, y: 1.0 + (x > 50.0))
    assert d.sediment_shear_factor.max() == 2.0 and d.sediment_shear_factor.min() == 1.0
    d.set_shear_amplification(None)
    assert d.sediment_shear_factor is None
    with pytest.raises(ValueError):
        d.set_shear_amplification(-1.0)

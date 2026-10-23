"""Full MMS: manufactured hydrodynamics AND sediment (RDy26 3.3).

test_sediment_mms.py verifies the source term in still water. Here the
hydrodynamics EVOLVE, forced by manufactured mass and momentum sources applied
through an operator, so the sediment is advected by a spatially and temporally
varying flow -- which makes this a measurement of SPATIAL convergence of the
sediment equation, and nothing else here does that.

Source terms are derived symbolically with sympy at import rather than by hand.
As a guard against mis-transcribing the manufactured solution, every derivative
sympy computes is checked against the one RDycore states explicitly; a mismatch
fails rather than silently poisoning the sources.

Two deliberate departures:

  * FRICTIONLESS. RDycore manufacture a Manning field too, but ANUGA applies
    friction semi-implicitly so an exactly cancelling source is not
    straightforward. n = 0 keeps the manufactured solution exact, and the
    convergence RATE being compared does not depend on it.
  * FIXED BED, as their MMS also has it: z is time-independent in the
    manufactured solution.

RDycore's own achieved L1 rates, from their YAML, are the target:
h 0.94, hu/hv 0.91, c0 0.94, c1 0.93.
"""
import numpy as np
import pytest

sp = pytest.importorskip('sympy')

import anuga
from anuga import rectangular_cross_domain
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions import Boundary

Hs, Ts, Us, Vs, Zs, Cs = 0.005, 20.0, 0.025, 0.025, 0.0025, 0.5
Ks = float(np.pi / 5.0)
GRAV = 9.81
STOP, DT = 5.0, 0.01
ALPHAS = (1, -1)
# RDycore's achieved rates, for comparison in the failure messages.
RDYCORE_RATES = {'h': 0.94, 'hu': 0.91, 'hv': 0.91, 'c0': 0.94, 'c1': 0.93}


_x, _y, _t = sp.symbols('x y t', real=True)
_H, _T, _U, _V, _Z, _C, _K, _al, _vs, _g = sp.symbols(
    'H T U V Z C K alpha v_s g', real=True)
_E = sp.exp(_t / _T)
_s = sp.sin(_K * _x) * sp.sin(_K * _y)
_h = _H * (1 + _s) * _E
_u = _U * sp.cos(_K * _x) * sp.sin(_K * _y) * _E
_v = _V * sp.sin(_K * _x) * sp.cos(_K * _y) * _E
_z = _Z * _s
_c = _C * (1 + _s) * sp.exp(_al * _t / _T)

# Guard: sympy's derivatives must equal the ones RDycore writes down.
_stated = [
    (sp.diff(_h, _x), _H * _K * sp.sin(_K * _y) * sp.cos(_K * _x) * _E),
    (sp.diff(_h, _y), _H * _K * sp.sin(_K * _x) * sp.cos(_K * _y) * _E),
    (sp.diff(_h, _t), _H / _T * (1 + _s) * _E),
    (sp.diff(_u, _x), -_U * _K * sp.sin(_K * _x) * sp.sin(_K * _y) * _E),
    (sp.diff(_u, _y), _U * _K * sp.cos(_K * _x) * sp.cos(_K * _y) * _E),
    (sp.diff(_v, _x), _K * _V * sp.cos(_K * _x) * sp.cos(_K * _y) * _E),
    (sp.diff(_v, _y), -_K * _V * sp.sin(_K * _x) * sp.sin(_K * _y) * _E),
    (sp.diff(_z, _x), _Z * _K * sp.cos(_K * _x) * sp.sin(_K * _y)),
    (sp.diff(_z, _y), _Z * _K * sp.sin(_K * _x) * sp.cos(_K * _y)),
]
_DERIVS_OK = all(sp.simplify(a - b) == 0 for a, b in _stated)

_S_h = sp.diff(_h, _t) + sp.diff(_h * _u, _x) + sp.diff(_h * _v, _y)
_S_hu = (sp.diff(_h * _u, _t) + sp.diff(_h * _u * _u + _g * _h * _h / 2, _x)
         + sp.diff(_h * _u * _v, _y) + _g * _h * sp.diff(_z, _x))
_S_hv = (sp.diff(_h * _v, _t) + sp.diff(_h * _u * _v, _x)
         + sp.diff(_h * _v * _v + _g * _h * _h / 2, _y) + _g * _h * sp.diff(_z, _y))
_m = _h * _c
_S_ms = sp.diff(_m, _t) + sp.diff(_u * _m, _x) + sp.diff(_v * _m, _y) + _c * _vs

_sub = {_H: Hs, _T: Ts, _U: Us, _V: Vs, _Z: Zs, _C: Cs, _K: Ks, _g: GRAV}
_L = lambda e: sp.lambdify((_x, _y, _t), e.subs(_sub), 'numpy')
h_ex, u_ex, v_ex = _L(_h), _L(_u), _L(_v)
z_ex = sp.lambdify((_x, _y), _z.subs(_sub), 'numpy')
S_h, S_hu, S_hv = _L(_S_h), _L(_S_hu), _L(_S_hv)
c_ex = {a: _L(_c.subs(_al, a)) for a in (1, -1)}


def make_S_ms(alpha, vs):
    return _L(_S_ms.subs(_al, alpha).subs(_vs, vs))


class MMS_boundary(Boundary):
    """Dirichlet on the exact solution -- what MMS requires at the boundary."""

    def __init__(self, domain):
        Boundary.__init__(self)
        self.domain = domain

    def evaluate_segment(self, domain, segment_edges):
        if segment_edges is None or domain is None:
            return
        ids = segment_edges
        vol = domain.boundary_cells[ids]
        edg = domain.boundary_edges[ids]
        em = domain.get_edge_midpoint_coordinates()
        xb, yb = em[3 * vol + edg, 0], em[3 * vol + edg, 1]
        t = domain.get_time()
        hb, ub, vb = h_ex(xb, yb, t), u_ex(xb, yb, t), v_ex(xb, yb, t)
        domain.quantities['stage'].boundary_values[ids] = z_ex(xb, yb) + hb
        domain.quantities['xmomentum'].boundary_values[ids] = hb * ub
        domain.quantities['ymomentum'].boundary_values[ids] = hb * vb
        for s, a in enumerate(ALPHAS):
            domain.tracer_boundary_values[s][ids] = c_ex[a](xb, yb, t)


ALPHAS = (1, -1)


class MMS_source_operator(anuga.operators.base_operator.Operator):
    """Applies the manufactured mass, momentum and sediment sources."""

    def __init__(self, domain, sms):
        anuga.operators.base_operator.Operator.__init__(self, domain)
        self.x = domain.centroid_coordinates[:, 0]
        self.y = domain.centroid_coordinates[:, 1]
        self.sms = sms

    def __call__(self):
        d = self.domain
        dt = d.get_timestep()
        if dt <= 0.0:
            return
        t = d.get_time()
        x, y = self.x, self.y
        d.quantities['stage'].centroid_values[:] += dt * S_h(x, y, t)
        d.quantities['xmomentum'].centroid_values[:] += dt * S_hu(x, y, t)
        d.quantities['ymomentum'].centroid_values[:] += dt * S_hv(x, y, t)
        for name, f in self.sms:
            d.set_tracer_source(name, f(x, y, t))

    def parallel_safe(self):
        return True

    def statistics(self):
        return 'MMS sources'

    def timestepping_statistics(self):
        return ''


def run(nxy, dt=None):
    d = rectangular_cross_domain(nxy, nxy, len1=5.0, len2=5.0)
    d.set_flow_algorithm('DE0')
    d.set_low_froude(0)
    d.store = False
    x, y = d.centroid_coordinates[:, 0], d.centroid_coordinates[:, 1]
    d.set_quantity('elevation', z_ex(x, y), location='centroids')
    d.set_quantity('stage', z_ex(x, y) + h_ex(x, y, 0.0), location='centroids')
    d.set_quantity('xmomentum', h_ex(x, y, 0.0) * u_ex(x, y, 0.0),
                   location='centroids')
    d.set_quantity('ymomentum', h_ex(x, y, 0.0) * v_ex(x, y, 0.0),
                   location='centroids')
    d.set_quantity('friction', 0.0)                # see the docstring
    d.evolve_max_timestep = DT if dt is None else dt
    d.sediment_bed_evolution = False               # z is time-independent
    d.sediment_c_max = 100.0                       # do not let limiters clip
    d.sediment_c_pack = 100.0

    names = []
    for i, a in enumerate(ALPHAS):
        nm = 'c%d' % i
        # Registered WITHOUT building the operator: fractional-step operators
        # run in registration order, and the source has to be applied before
        # the sediment operator consumes it, so the operator is created below
        # -- after MMS_source_operator. Sediment_transport_operator() would
        # otherwise put itself first.
        d._allocate_sediment_fraction(
            nm, diameter=1.0e-4, tau_c_star=0.0,
            initial_concentration=c_ex[a](x, y, 0.0))
        names.append(nm)
    vs = d.sediment_settling_velocity[0]

    d.set_boundary({tg: MMS_boundary(d) for tg in d.get_boundary_tags()})

    # order matters: sources must be set before the sediment operator consumes
    # them, and operators run in registration order.
    MMS_source_operator(d, [(nm, make_S_ms(a, vs))
                            for nm, a in zip(names, ALPHAS)])
    # No grain size argument: it picks up the ones registered above.
    d.initialize_sediment_operator()

    d.evolve_to_end(finaltime=STOP)

    a_ = d.areas
    out = {}
    hn = np.maximum(d.quantities['stage'].centroid_values
                    - d.quantities['elevation'].centroid_values, 0.0)
    fields = {'h': (hn, h_ex(x, y, STOP)),
              'hu': (d.quantities['xmomentum'].centroid_values,
                     h_ex(x, y, STOP) * u_ex(x, y, STOP)),
              'hv': (d.quantities['ymomentum'].centroid_values,
                     h_ex(x, y, STOP) * v_ex(x, y, STOP))}
    for i, a in enumerate(ALPHAS):
        fields['c%d' % i] = (d.get_tracer('c%d' % i), c_ex[a](x, y, STOP))
    for k, (num, exact) in fields.items():
        e = np.abs(num - exact)
        out[k] = (float((e * a_).sum() / a_.sum()),
                  float(np.sqrt((e**2 * a_).sum() / a_.sum())),
                  float(e.max()))
    return out


# Refine SPACE AND TIME TOGETHER. This is the whole ballgame and the first
# attempt got it wrong: refining dx three times while holding dt = 0.01 fixed
# leaves a CONSTANT temporal error, which comes to dominate as the spatial
# error shrinks, and the apparent rate decays -- 0.802 then 0.509. The error
# field gives it away, being spatially UNIFORM (boundary band 4.46e-4 against
# interior 4.90e-4 at the finest mesh), so it is not a boundary artefact and
# not a spatial one. Halving dt with dx recovers first order.
#
# RDycore's own YAML holds dt fixed across refinements and still reports 0.94,
# so their spatial error presumably still dominates at their resolutions. Ours
# does not, so the refinement path has to be diagonal.
def run(nxy, dt=None):
    d = rectangular_cross_domain(nxy, nxy, len1=5.0, len2=5.0)
    d.set_flow_algorithm('DE0')
    d.set_low_froude(0)
    d.store = False
    x, y = d.centroid_coordinates[:, 0], d.centroid_coordinates[:, 1]
    d.set_quantity('elevation', z_ex(x, y), location='centroids')
    d.set_quantity('stage', z_ex(x, y) + h_ex(x, y, 0.0), location='centroids')
    d.set_quantity('xmomentum', h_ex(x, y, 0.0) * u_ex(x, y, 0.0),
                   location='centroids')
    d.set_quantity('ymomentum', h_ex(x, y, 0.0) * v_ex(x, y, 0.0),
                   location='centroids')
    d.set_quantity('friction', 0.0)                # see the docstring
    d.evolve_max_timestep = DT if dt is None else dt
    d.sediment_bed_evolution = False               # z is time-independent
    d.sediment_c_max = 100.0                       # do not let limiters clip
    d.sediment_c_pack = 100.0

    names = []
    for i, a in enumerate(ALPHAS):
        nm = 'c%d' % i
        # Registered WITHOUT building the operator: fractional-step operators
        # run in registration order, and the source has to be applied before
        # the sediment operator consumes it, so the operator is created below
        # -- after MMS_source_operator. Sediment_transport_operator() would
        # otherwise put itself first.
        d._allocate_sediment_fraction(
            nm, diameter=1.0e-4, tau_c_star=0.0,
            initial_concentration=c_ex[a](x, y, 0.0))
        names.append(nm)
    vs = d.sediment_settling_velocity[0]

    d.set_boundary({tg: MMS_boundary(d) for tg in d.get_boundary_tags()})

    # order matters: sources must be set before the sediment operator consumes
    # them, and operators run in registration order.
    MMS_source_operator(d, [(nm, make_S_ms(a, vs))
                            for nm, a in zip(names, ALPHAS)])
    # No grain size argument: it picks up the ones registered above.
    d.initialize_sediment_operator()

    d.evolve_to_end(finaltime=STOP)

    a_ = d.areas
    out = {}
    hn = np.maximum(d.quantities['stage'].centroid_values
                    - d.quantities['elevation'].centroid_values, 0.0)
    fields = {'h': (hn, h_ex(x, y, STOP)),
              'hu': (d.quantities['xmomentum'].centroid_values,
                     h_ex(x, y, STOP) * u_ex(x, y, STOP)),
              'hv': (d.quantities['ymomentum'].centroid_values,
                     h_ex(x, y, STOP) * v_ex(x, y, STOP))}
    for i, a in enumerate(ALPHAS):
        fields['c%d' % i] = (d.get_tracer('c%d' % i), c_ex[a](x, y, STOP))
    for k, (num, exact) in fields.items():
        e = np.abs(num - exact)
        out[k] = (float((e * a_).sum() / a_.sum()),
                  float(np.sqrt((e**2 * a_).sum() / a_.sum())),
                  float(e.max()))
    return out


# Refine SPACE AND TIME TOGETHER. This is the whole ballgame and the first
# attempt got it wrong: refining dx three times while holding dt = 0.01 fixed
# leaves a CONSTANT temporal error, which comes to dominate as the spatial
# error shrinks, and the apparent rate decays -- 0.802 then 0.509. The error
# field gives it away, being spatially UNIFORM (boundary band 4.46e-4 against
# interior 4.90e-4 at the finest mesh), so it is not a boundary artefact and
# not a spatial one. Halving dt with dx recovers first order.
#
# RDycore's own YAML holds dt fixed across refinements and still reports 0.94,
# so their spatial error presumably still dominates at their resolutions. Ours
# does not, so the refinement path has to be diagonal.


GRIDS = ((10, 0.02), (20, 0.01), (40, 0.005))     # dx = 0.5, 0.25, 0.125
FIELDS = ('h', 'hu', 'hv', 'c0', 'c1')


@pytest.fixture(scope='module')
def convergence():
    """L1 errors on three meshes, and the mean observed rate per field."""
    res = {n: run(n, dt=dt) for n, dt in GRIDS}
    meshes = [n for n, _ in GRIDS]
    rates = {}
    for f in FIELDS:
        L1 = [res[n][f][0] for n in meshes]
        r = [np.log(L1[i] / L1[i + 1]) / np.log(2.0) for i in range(len(L1) - 1)]
        rates[f] = float(np.mean(r))
    return res, meshes, rates


def test_the_symbolic_derivatives_match_the_ones_RDycore_states():
    """Guards against mis-transcribing the manufactured solution: a wrong
    derivative would poison every source term silently."""
    assert _DERIVS_OK


@pytest.mark.slow
def test_every_field_converges_under_refinement(convergence):
    _, _, rates = convergence
    assert all(rates[f] > 0.5 for f in FIELDS), rates


@pytest.mark.slow
def test_sediment_converges_at_first_order(convergence):
    """RDycore report 0.94 and 0.93 for the two classes."""
    _, _, rates = convergence
    assert all(0.8 < rates[f] < 1.3 for f in ('c0', 'c1')), rates


@pytest.mark.slow
def test_sediment_converges_no_worse_than_the_flow_carrying_it(convergence):
    _, _, rates = convergence
    assert (min(rates['c0'], rates['c1'])
            > 0.8 * min(rates['h'], rates['hu'])), rates


@pytest.mark.slow
def test_the_absolute_error_is_small_at_the_finest_mesh(convergence):
    """Order alone can be met by a solution that is uniformly wrong."""
    res, meshes, _ = convergence
    for f in ('c0', 'c1'):
        assert res[meshes[-1]][f][0] < 0.05 * max(abs(Cs), Hs)

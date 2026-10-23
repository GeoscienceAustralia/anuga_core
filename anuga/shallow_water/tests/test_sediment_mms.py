"""Method of Manufactured Solutions for the sediment source terms (RDy26 3.3).

MMS is the rigorous way to verify an implementation: prescribe a smooth exact
solution, derive the source term that makes it exact, then check the code
recovers it at the design order. Unlike the other sediment tests, a failure
here cannot be argued away -- the answer is known.

RDy26's own MMS takes its manufactured hydrodynamics from Bisht et al. (2025),
which is not to hand, so their case cannot be reproduced exactly. Their
manufactured SEDIMENT field is given in full (their B2-B5) and is used here:

    c_s(x,y,t) = C [1 + sin(Kx) sin(Ky)] exp(alpha_s t/T)
    K = pi/L,  L = 5 m,  T = 20 s,  C = 0.5,  alpha_1 = +1, alpha_2 = -1

Run in STILL WATER, so u = v = 0 and the advective term vanishes. That is
deliberate rather than a dodge: advection is already verified to machine
precision by the RDy26 dam-break contact test (test_sediment_rdycore.py, SSC
flat to 3.6e-15 across a rarefaction and a shock), while nothing else verifies
the SOURCE-TERM coupling and its order in time -- which is exactly what RDy26
say their MMS adds over the passive cases. The two together cover both halves
of [G-3].

With u = 0 the manufactured source reduces to

    m = h c,  dm/dt = h dc/dt,  E = 0 (no shear),  D = d* c v_s
    => S_ms = h dc/dt + d* c v_s

and the expected convergence is FIRST ORDER in dt: the source is applied as a
fractional step with the full timestep.
"""
import numpy as np
import pytest

import anuga
from anuga import Reflective_boundary, rectangular_cross_domain

C_SCALE, L_MMS, T_MMS = 0.5, 5.0, 20.0
K = np.pi / L_MMS
DEPTH = 1.0
D50 = 1.0e-4
T_END = 4.0


def c_exact(x, y, t, alpha):
    return (C_SCALE * (1.0 + np.sin(K * x) * np.sin(K * y))
            * np.exp(alpha * t / T_MMS))


def dcdt_exact(x, y, t, alpha):
    return c_exact(x, y, t, alpha) * (alpha / T_MMS)


def run(nxy, dt, alpha=1.0):
    """One MMS run. Returns the L1, L2 and Linf errors in the conserved m."""
    d = rectangular_cross_domain(nxy, nxy, len1=L_MMS, len2=L_MMS)
    d.set_flow_algorithm('DE0')
    d.set_low_froude(0)
    d.store = False
    d.set_quantity('elevation', 0.0)
    d.set_quantity('stage', DEPTH)
    d.set_quantity('friction', 0.0)
    d.set_boundary({t: Reflective_boundary(d) for t in d.get_boundary_tags()})
    d.evolve_max_timestep = dt
    d.sediment_bed_evolution = False       # fixed bed: h stays at DEPTH
    d.sediment_c_max = 10.0                # do not let [L-2] clip the solution
    d.sediment_c_pack = 10.0

    x, y = d.centroid_coordinates[:, 0], d.centroid_coordinates[:, 1]
    d.add_sediment_fraction(name='mms', diameter=D50, tau_c_star=0.0,
                     initial_concentration=c_exact(x, y, 0.0, alpha))
    v_s = d.sediment_settling_velocity[0]

    class MMS_source(anuga.operators.base_operator.Operator):
        """The manufactured source, applied through the ordinary hook.

        Using set_tracer_source rather than reaching into the kernel means the
        test exercises the same path a user would.
        """

        def __call__(self):
            t = self.domain.get_time()
            S = (DEPTH * dcdt_exact(x, y, t, alpha)
                 + c_exact(x, y, t, alpha) * v_s)
            self.domain.set_tracer_source('mms', S)

        def parallel_safe(self):
            return True

        def statistics(self):
            return 'MMS source'

        def timestepping_statistics(self):
            return ''

    MMS_source(d)
    d.evolve_to_end(finaltime=T_END)

    m_num = d.tracer_conserved_values[0]
    m_ex = DEPTH * c_exact(x, y, T_END, alpha)
    e = np.abs(m_num - m_ex)
    a = d.areas
    return (float((e * a).sum() / a.sum()),
            float(np.sqrt((e ** 2 * a).sum() / a.sum())),
            float(e.max()))


@pytest.mark.slow
def test_the_manufactured_solution_is_recovered():
    L1, _, _ = run(20, 0.005)
    assert L1 / (DEPTH * C_SCALE) < 0.01, 'L1 = %.4e' % L1


@pytest.mark.slow
def test_the_decaying_class_is_recovered_too():
    """alpha = -1. A source term that only worked for a growing solution would
    otherwise pass."""
    L1, _, _ = run(20, 0.005, alpha=-1.0)
    assert L1 / (DEPTH * C_SCALE) < 0.01


@pytest.mark.slow
def test_the_source_converges_at_first_order_in_dt():
    """First order is the DESIGN order here: the source is a fractional step
    applied with the full timestep. Second order would mean the test is
    measuring something else.

    Note evolve_max_timestep is a CAP, not a setting -- if the CFL condition
    binds first, every refinement level runs at the same achieved dt and the
    measured rate is meaningless. The still, shallow domain here keeps the cap
    binding.
    """
    errs = []
    for dt in (0.008, 0.004, 0.002, 0.001):
        L1, _, _ = run(20, dt)
        errs.append(L1)
    rates = [np.log(errs[i] / errs[i + 1]) / np.log(2.0)
             for i in range(len(errs) - 1)]
    assert all(r > 0 for r in rates), 'error did not fall monotonically: %r' % rates
    assert 0.8 < np.mean(rates) < 1.3, 'observed order %.3f' % np.mean(rates)


def test_an_external_source_is_delivered_against_a_tight_ceiling():
    """[L-2] bounds BED EXCHANGE by what the water column can hold. An external
    supply is neither bed nor water column, so clipping it would also make a
    manufactured solution impossible to impose exactly."""
    d = rectangular_cross_domain(4, 4, len1=L_MMS, len2=L_MMS)
    d.set_flow_algorithm('DE0')
    d.store = False
    d.set_quantity('elevation', 0.0)
    d.set_quantity('stage', DEPTH)
    d.set_boundary({t: Reflective_boundary(d) for t in d.get_boundary_tags()})
    d.sediment_c_max = 1e-6                # absurdly tight
    d.add_sediment_fraction(name='s', diameter=D50, tau_c_star=0.0, d_star=0.0,
                     initial_concentration=0.0)
    d.set_tracer_source('s', 1.0e-3)
    m0 = float((d.tracer_conserved_values[0] * d.areas).sum())
    d.evolve_to_end(finaltime=2.0)
    m1 = float((d.tracer_conserved_values[0] * d.areas).sum())
    assert m1 > m0

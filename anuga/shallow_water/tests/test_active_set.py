"""Tests for active-set stepping (domain.set_use_active_set).

Active-set stepping skips cells that are provably unchangeable this step
(dry, with an all-dry 2-ring neighbourhood).  These tests verify it against
plain full stepping on a mostly-dry domain, that its statistics report real
engagement, and that water added by a rate operator (rain) on inactive cells
is never lost.
"""

import os
import tempfile
import unittest
import warnings

import numpy as np
import pytest

import anuga
from anuga import Reflective_boundary, rectangular_cross_domain


_gpu_error = None
_gpu_avail = None


def gpu_available():
    global _gpu_error, _gpu_avail
    if _gpu_avail is not None:
        return _gpu_avail
    try:
        from anuga.shallow_water.sw_domain_gpu_ext import init_gpu_domain  # noqa: F401
        _gpu_avail = True
    except Exception as e:
        _gpu_avail = False
        _gpu_error = f"{type(e).__name__}: {e}"
    return _gpu_avail


# Same per-process isolation rule as test_DE_gpu_omp.py: on a GPU-offload
# build the NVHPC runtime aborts after many mode-2 domains in one process, so
# this file only runs in-process on CPU builds or under the isolated runner.
if (gpu_available() and anuga.gpu_offload_supported()
        and not os.environ.get('ANUGA_GPU_TESTS_ISOLATED')):
    _skip_reason = (
        "GPU-offload build: run active-set tests via "
        "anuga/shallow_water/tests/run_gpu_tests_isolated.sh")
    warnings.warn(_skip_reason, stacklevel=1)
    pytest.skip(_skip_reason, allow_module_level=True)


def _make_domain(name):
    """A mostly-dry sloped domain: wet pool in the low corner, dry upslope."""
    domain = rectangular_cross_domain(16, 16, len1=100., len2=100.)
    domain.set_flow_algorithm('DE1')
    domain.set_low_froude(0)
    domain.set_name(name)
    domain.set_datadir(tempfile.mkdtemp())
    domain.store = False
    domain.set_quantity('elevation', lambda x, y: x / 10.0)   # 0..10 m
    domain.set_quantity('friction', 0.0)
    domain.set_quantity('stage', 2.0)                         # wet only x < 20
    Br = Reflective_boundary(domain)
    domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})
    return domain


@pytest.mark.skipif(not gpu_available(),
                    reason=_gpu_error or "GPU OpenMP interface not available")
class Test_active_set(unittest.TestCase):

    def _evolve(self, domain, finaltime=3.0):
        for _ in domain.evolve(yieldstep=1.0, finaltime=finaltime):
            pass

    def test_active_matches_full_stepping(self):
        """Active-set evolution matches plain mode-2 evolution.

        The active-set domain uses scatter fluxes (single Riemann solve per
        edge), which agree with the default cell-based kernel only to
        floating-point roundoff, hence the tolerance.
        """
        ref = _make_domain('as_ref')
        ref.set_multiprocessor_mode(2)
        self._evolve(ref)

        act = _make_domain('as_act')
        act.set_multiprocessor_mode(2)
        act.set_use_active_set(True)
        self._evolve(act)

        for q in ('stage', 'xmomentum', 'ymomentum'):
            a = act.quantities[q].centroid_values
            r = ref.quantities[q].centroid_values
            assert np.allclose(a, r, atol=1e-8, rtol=1e-6), \
                f'{q}: max diff {np.abs(a - r).max()}'

        frac, samples = act.get_active_set_stats()
        assert samples > 0, 'active set never engaged'
        assert frac < 0.9, f'active fraction {frac} suspiciously high'
        # the reference never engaged it
        assert ref.get_active_set_stats() == (1.0, 0)

    def test_rain_on_inactive_cells_is_not_lost(self):
        """A rate operator wetting dry (inactive) cells conserves volume."""
        from anuga.operators.rate_operators import Rate_operator

        domain = _make_domain('as_rain')
        domain.set_multiprocessor_mode(2)
        domain.set_use_active_set(True)

        # Rain on a dry upslope region only (x in [60, 90]) -- cells the
        # active set is guaranteed to be skipping.
        region = [[60., 5.], [90., 5.], [90., 95.], [60., 95.]]
        rate = 0.01  # m/s over the polygon
        op = Rate_operator(domain, rate=rate, polygon=region)

        v0 = domain.get_water_volume()
        t_end = 2.0
        for _ in domain.evolve(yieldstep=1.0, finaltime=t_end):
            pass
        v1 = domain.get_water_volume()

        added = v1 - v0
        # The operator rains on cells whose centroids lie in the polygon, so
        # the covered area is the sum of those cell areas, not the polygon's.
        covered = float(np.sum(domain.areas[op.indices]))
        expected = rate * covered * t_end
        assert added > 0.5 * expected, \
            f'rain lost: added {added}, expected ~{expected}'
        assert abs(added - expected) < 1e-6 * expected, \
            f'rain volume off: added {added}, expected {expected}'

        frac, samples = domain.get_active_set_stats()
        assert samples > 0

    def test_tracers_disable_the_active_set(self):
        """With a passive tracer registered the flag is ignored, not honoured
        with a kernel that cannot advect the tracer.

        The active path drives the scatter flux kernel directly, and only the
        cell-based kernel carries the tracer flux.  So the run must (a) never
        engage the active set and (b) reproduce the plain mode-2 run exactly
        -- same kernels, same launch order, hence bit-identical -- for the
        hydrodynamics AND the tracer.
        """
        def build(name):
            d = _make_domain(name)
            # A step in the pool so there is flow to carry the tracer (the
            # base pool is at rest).
            d.set_quantity('stage', lambda x, y: np.where(x < 10.0, 2.5, 2.0))
            d.add_tracer('c', beta=1.0)
            x = d.centroid_coordinates[:, 0]
            d.set_tracer('c', np.where(x < 10.0, 1.0, 0.0))
            d.set_multiprocessor_mode(2)
            return d

        ref = build('as_tr_ref')
        self._evolve(ref)

        act = build('as_tr_act')
        act.set_use_active_set(True)
        self._evolve(act)

        assert act.get_active_set_stats() == (1.0, 0), \
            'active set engaged on a domain with tracers'
        for q in ('stage', 'xmomentum', 'ymomentum'):
            a = act.quantities[q].centroid_values
            r = ref.quantities[q].centroid_values
            assert np.array_equal(a, r), f'{q}: max diff {np.abs(a - r).max()}'
        a = act.get_tracer('c')
        r = ref.get_tracer('c')
        assert np.array_equal(a, r), f'tracer: max diff {np.abs(a - r).max()}'
        # and the tracer actually moved, so the comparison is not vacuous
        initial = np.where(act.centroid_coordinates[:, 0] < 10.0, 1.0, 0.0)
        assert np.abs(a - initial).max() > 1e-3


if __name__ == '__main__':
    unittest.main()

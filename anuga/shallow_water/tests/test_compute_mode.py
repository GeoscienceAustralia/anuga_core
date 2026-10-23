"""Tests for the compute-backend model.

Two orthogonal knobs:
  * per-domain compute mode  — Domain.set_compute_mode('legacy' | 'unified'),
    selecting the internal multiprocessor_mode (1 / 2);
  * process-global GPU offload — anuga.set_gpu_offload(bool), deciding whether
    'unified' domains run on a GPU or CPU-multicore.

Robust to both CPU-only (gpu_offload=false) and GPU (gpu_offload=true) builds.
"""

import unittest
import warnings

import anuga
from anuga import Domain, Reflective_boundary, rectangular_cross
from anuga.config import MULTIPROCESSOR_OPENMP, MULTIPROCESSOR_GPU
from anuga.shallow_water.sw_domain_gpu_ext import gpu_available


def _make_domain():
    points, vertices, boundary = rectangular_cross(6, 6, len1=6.0, len2=6.0)
    domain = Domain(points, vertices, boundary)
    domain.set_quantity('elevation', 0.0)
    domain.set_quantity('stage', 1.0)
    R = Reflective_boundary(domain)
    domain.set_boundary({'left': R, 'right': R, 'top': R, 'bottom': R})
    return domain


class Test_compute_mode(unittest.TestCase):

    def test_default_is_legacy(self):
        # The default is legacy only when not overridden by the opt-in env var
        # ANUGA_DEFAULT_COMPUTE_MODE (the test suite sets it to 'unified' to
        # exercise the whole suite in mode 2). Skip when that override is active.
        import os
        if os.environ.get('ANUGA_DEFAULT_COMPUTE_MODE', 'legacy').lower() != 'legacy':
            self.skipTest('ANUGA_DEFAULT_COMPUTE_MODE overrides the default mode')
        domain = _make_domain()
        self.assertEqual(domain.get_compute_mode(), 'legacy')
        self.assertEqual(domain.multiprocessor_mode, MULTIPROCESSOR_OPENMP)
        self.assertFalse(domain.use_c_rk_loop)

    def test_legacy_maps_to_mode1(self):
        domain = _make_domain()
        domain.set_compute_mode('legacy')
        self.assertEqual(domain.get_compute_mode(), 'legacy')
        self.assertEqual(domain.multiprocessor_mode, MULTIPROCESSOR_OPENMP)
        self.assertFalse(domain.use_c_rk_loop)

    def test_unified_maps_to_mode2(self):
        domain = _make_domain()
        domain.set_compute_mode('unified')
        self.assertEqual(domain.get_compute_mode(), 'unified')
        self.assertEqual(domain.multiprocessor_mode, MULTIPROCESSOR_GPU)
        self.assertTrue(domain.use_c_rk_loop)
        self.assertIsNotNone(domain.gpu_interface)

    def test_cpu_and_gpu_are_not_modes(self):
        domain = _make_domain()
        for bad in ('cpu', 'gpu', 'bogus'):
            with self.assertRaises(ValueError):
                domain.set_compute_mode(bad)

    def test_legacy_int_api_maps_to_compute_mode(self):
        domain = _make_domain()
        domain.set_multiprocessor_mode(1)
        self.assertEqual(domain.get_compute_mode(), 'legacy')
        self.assertEqual(domain.multiprocessor_mode, MULTIPROCESSOR_OPENMP)

        domain.set_multiprocessor_mode(2)
        self.assertEqual(domain.get_compute_mode(), 'unified')
        self.assertEqual(domain.multiprocessor_mode, MULTIPROCESSOR_GPU)

    def test_invalid_int_mode_raises(self):
        domain = _make_domain()
        with self.assertRaises(ValueError):
            domain.set_multiprocessor_mode(3)

    def test_switch_back_to_legacy(self):
        domain = _make_domain()
        domain.set_compute_mode('unified')
        domain.set_compute_mode('legacy')
        self.assertEqual(domain.get_compute_mode(), 'legacy')
        self.assertEqual(domain.multiprocessor_mode, MULTIPROCESSOR_OPENMP)
        self.assertFalse(domain.use_c_rk_loop)

    def test_compute_capabilities(self):
        domain = _make_domain()
        caps = domain.compute_capabilities()
        self.assertEqual(set(caps), {'gpu_offload', 'num_gpu_devices', 'mpi', 'modes'})
        self.assertIn('legacy', caps['modes'])
        # 'unified' is available whenever the gpu_ext extension imports (it does
        # here, since the domain built a gpu_interface in other tests).
        self.assertIn('unified', caps['modes'])
        self.assertEqual(caps['gpu_offload'], anuga.gpu_offload_enabled())

    def test_parallel_without_mpi_build_falls_back_to_legacy(self):
        # Simulate a multi-rank run on a gpu_ext built WITHOUT MPI: 'unified'
        # C ghost exchange would be a silent no-op, so the selector must fall
        # back to 'legacy' rather than compute wrong parallel results.
        domain = _make_domain()
        domain._mode2_mpi_available = lambda: False
        saved = anuga.numprocs
        try:
            anuga.numprocs = 2
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always')
                domain.set_compute_mode('unified')
                warned = any('Python MPI exchange' in str(w.message) for w in caught)
        finally:
            anuga.numprocs = saved

        self.assertEqual(domain.get_compute_mode(), 'legacy')
        self.assertEqual(domain.multiprocessor_mode, MULTIPROCESSOR_OPENMP)
        self.assertTrue(warned)

    def test_parallel_with_mpi_build_keeps_unified(self):
        domain = _make_domain()
        domain._mode2_mpi_available = lambda: True
        saved = anuga.numprocs
        try:
            anuga.numprocs = 2
            domain.set_compute_mode('unified')
        finally:
            anuga.numprocs = saved
        self.assertEqual(domain.get_compute_mode(), 'unified')
        self.assertEqual(domain.multiprocessor_mode, MULTIPROCESSOR_GPU)


class Test_gpu_offload(unittest.TestCase):
    """Process-global offload toggle. Restores state after each test."""

    def setUp(self):
        from anuga.shallow_water import sw_domain_gpu_ext as gpu_ext
        self._saved = gpu_ext.get_offload_enabled()

    def tearDown(self):
        # Restore the process-global offload state: clear any disable env var we
        # set and restore the C flag directly (avoids the not-supported warning
        # that set_gpu_offload(True) emits on CPU builds).
        import os
        from anuga.shallow_water import sw_domain_gpu_ext as gpu_ext
        os.environ.pop('OMP_TARGET_OFFLOAD', None)
        gpu_ext.set_offload_enabled(self._saved)

    def test_disable_is_process_global(self):
        state = anuga.set_gpu_offload(False, verbose=False)
        self.assertFalse(state)
        self.assertFalse(anuga.gpu_offload_enabled())

    def test_enable_on_cpu_only_build_warns_and_stays_off(self):
        if gpu_available():
            self.skipTest("GPU build: enabling offload is valid here")
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            state = anuga.set_gpu_offload(True, verbose=False)
            warned = any('no GPU offload support' in str(w.message) for w in caught)
        self.assertFalse(state)
        self.assertTrue(warned)

    def test_offload_state_drives_domain_bookkeeping(self):
        # With offload disabled, a 'unified' domain reports CPU (offload inactive)
        # regardless of build.
        anuga.set_gpu_offload(False, verbose=False)
        points, vertices, boundary = rectangular_cross(6, 6, len1=6.0, len2=6.0)
        domain = Domain(points, vertices, boundary)
        domain.set_quantity('elevation', 0.0)
        domain.set_quantity('stage', 1.0)
        R = Reflective_boundary(domain)
        domain.set_boundary({'left': R, 'right': R, 'top': R, 'bottom': R})
        domain.set_compute_mode('unified')
        self.assertFalse(domain.gpu_offload_active)


if __name__ == '__main__':
    unittest.main()

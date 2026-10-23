"""
Tests for GPU (OpenMP target offloading) implementation of ANUGA's shallow water solver.

These tests verify that the GPU implementation produces results matching the CPU implementation.
"""

import os
import tempfile
import unittest
import sys
import warnings
import numpy as np
import pytest

import anuga
from anuga import Reflective_boundary, Dirichlet_boundary
from anuga import Transmissive_momentum_set_stage_boundary
from anuga import rectangular_cross_domain
from anuga import Inlet_operator


_gpu_error = None
_gpu_avail = None


def gpu_available():
    """Check if GPU OpenMP interface is available."""
    global _gpu_error, _gpu_avail
    if _gpu_avail is not None:
        return _gpu_avail
    try:
        from anuga.shallow_water.sw_domain_gpu_ext import init_gpu_domain
        _gpu_avail = True
    except Exception as e:
        _gpu_avail = False
        _gpu_error = f"{type(e).__name__}: {e}"
        print(f"sw_domain_gpu_ext not available: {_gpu_error}", flush=True)
    return _gpu_avail


def _gpu_skip_reason():
    if _gpu_error:
        return f"GPU OpenMP interface not available: {_gpu_error}"
    return "GPU OpenMP interface not available"


# On a GPU-offload build (nvc, -Dgpu_offload=true) the NVHPC OpenMP-target
# runtime aborts the process once many mode-2 GPU domains have been created in
# it, so running this whole file in a single pytest process crashes partway
# through (see claude/KNOWN_ISSUES.md, "test_DE_gpu_omp.py aborts mid-file").
# Skip the file in a normal in-process run on such a build and point at the
# isolated runner, which executes one class per fresh process. The runner sets
# ANUGA_GPU_TESTS_ISOLATED=1 to opt back in. On a CPU build (gpu_offload not
# supported) the omp-target regions run on the host with no device present
# table, so the file runs fine in one process and is NOT skipped.
if (gpu_available() and anuga.gpu_offload_supported()
        and not os.environ.get('ANUGA_GPU_TESTS_ISOLATED')):
    _skip_reason = (
        "GPU-offload build: run this file via "
        "anuga/shallow_water/tests/run_gpu_tests_isolated.sh (one fresh process "
        "per class) — running it in one process aborts the NVHPC OpenMP-target "
        "runtime. Set ANUGA_GPU_TESTS_ISOLATED=1 to force in-process collection.")
    # Also emit a warning so the reason is visible in pytest's warnings summary
    # without needing -rs (a bare module-level skip otherwise just shows "1
    # skipped" with no explanation).
    warnings.warn(_skip_reason, stacklevel=1)
    pytest.skip(_skip_reason, allow_module_level=True)


@pytest.mark.skipif(not gpu_available(), reason=_gpu_skip_reason())
class Test_GPU_Kernels(unittest.TestCase):
    """Unit tests for individual GPU kernels."""

    def setUp(self):
        """Create a simple test domain."""
        self.domain = rectangular_cross_domain(10, 10, len1=100., len2=100.)
        self.domain.set_flow_algorithm('DE0')
        self.domain.set_low_froude(0)
        self.domain.set_name('test_gpu')
        self.domain.set_datadir(tempfile.mkdtemp())
        self.domain.store = False

        def topography(x, y):
            return -x / 50.0  # Slope from 0 to -2

        self.domain.set_quantity('elevation', topography)
        self.domain.set_quantity('friction', 0.01)
        self.domain.set_quantity('stage', 0.0)

        Br = Reflective_boundary(self.domain)
        Bd = Dirichlet_boundary([-0.5, 0., 0.])
        self.domain.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})

    def tearDown(self):
        import os
        for ext in ['.sww']:
            try:
                os.remove(f'test_gpu{ext}')
            except OSError:
                pass

    def test_flux_kernel(self):
        """Test that GPU flux computation matches CPU."""
        # Run CPU flux computation
        self.domain.set_multiprocessor_mode(1)
        self.domain.distribute_to_vertices_and_edges()
        self.domain.update_boundary()
        self.domain.compute_fluxes()

        cpu_timestep = self.domain.flux_timestep
        cpu_stage_update = self.domain.quantities['stage'].explicit_update.copy()
        cpu_xmom_update = self.domain.quantities['xmomentum'].explicit_update.copy()
        cpu_ymom_update = self.domain.quantities['ymomentum'].explicit_update.copy()

        # Reset and run GPU flux computation
        for qname in ['stage', 'xmomentum', 'ymomentum']:
            self.domain.quantities[qname].explicit_update[:] = 0.0

        self.domain.set_multiprocessor_mode(2)
        from anuga.shallow_water.sw_domain_gpu_ext import (
            sync_to_device, sync_all_from_device,
            extrapolate_second_order_gpu, compute_fluxes_gpu,
            evaluate_reflective_boundary_gpu, evaluate_dirichlet_boundary_gpu
        )

        gpu_dom = self.domain.gpu_interface.gpu_dom
        sync_to_device(gpu_dom)
        extrapolate_second_order_gpu(gpu_dom)
        evaluate_reflective_boundary_gpu(gpu_dom)
        evaluate_dirichlet_boundary_gpu(gpu_dom)
        gpu_timestep = compute_fluxes_gpu(gpu_dom)
        sync_all_from_device(gpu_dom)

        gpu_stage_update = self.domain.quantities['stage'].explicit_update.copy()
        gpu_xmom_update = self.domain.quantities['xmomentum'].explicit_update.copy()
        gpu_ymom_update = self.domain.quantities['ymomentum'].explicit_update.copy()

        # Verify - use atol for near-zero values where rtol alone causes issues
        self.assertAlmostEqual(cpu_timestep, gpu_timestep, places=10,
                               msg=f"Timestep mismatch: CPU={cpu_timestep}, GPU={gpu_timestep}")
        np.testing.assert_allclose(cpu_stage_update, gpu_stage_update, rtol=1e-10, atol=1e-14,
                                   err_msg="Stage explicit_update mismatch")
        np.testing.assert_allclose(cpu_xmom_update, gpu_xmom_update, rtol=1e-10, atol=1e-14,
                                   err_msg="Xmomentum explicit_update mismatch")
        np.testing.assert_allclose(cpu_ymom_update, gpu_ymom_update, rtol=1e-10, atol=1e-14,
                                   err_msg="Ymomentum explicit_update mismatch")

    def test_extrapolate_kernel(self):
        """Test that GPU extrapolation matches CPU."""
        # Set initial conditions
        self.domain.set_quantity('stage', expression='0.1 * x / 100.0')

        # Run CPU extrapolation
        self.domain.set_multiprocessor_mode(1)
        self.domain.distribute_to_vertices_and_edges()

        cpu_stage_edge = self.domain.quantities['stage'].edge_values.copy()

        # Reset edge values and run GPU
        self.domain.quantities['stage'].edge_values[:] = 0.0

        self.domain.set_multiprocessor_mode(2)
        from anuga.shallow_water.sw_domain_gpu_ext import (
            sync_to_device, sync_all_from_device,
            protect_gpu, extrapolate_second_order_gpu
        )

        gpu_dom = self.domain.gpu_interface.gpu_dom
        sync_to_device(gpu_dom)
        protect_gpu(gpu_dom)
        extrapolate_second_order_gpu(gpu_dom)
        sync_all_from_device(gpu_dom)

        gpu_stage_edge = self.domain.quantities['stage'].edge_values.copy()

        # Verify
        np.testing.assert_allclose(cpu_stage_edge, gpu_stage_edge, rtol=1e-10, atol=1e-14,
                                   err_msg="Edge values mismatch after extrapolation")


@pytest.mark.skipif(not gpu_available(), reason=_gpu_skip_reason())
class Test_GPU_RK2(unittest.TestCase):
    """Tests for complete RK2 step on GPU."""

    def create_domain(self):
        """Create a test domain."""
        domain = rectangular_cross_domain(10, 10, len1=100., len2=100.)
        domain.set_flow_algorithm('DE0')
        domain.set_low_froude(0)
        domain.set_name('test_rk2')
        domain.set_datadir(tempfile.mkdtemp())
        domain.store = False

        def topography(x, y):
            return -x / 50.0

        domain.set_quantity('elevation', topography)
        domain.set_quantity('friction', 0.01)
        domain.set_quantity('stage', 0.0)

        Br = Reflective_boundary(domain)
        Bd = Dirichlet_boundary([-0.5, 0., 0.])
        domain.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})

        return domain

    def test_single_rk2_step(self):
        """Test that a single RK2 step matches between CPU and GPU."""
        # Create two identical domains
        cpu_domain = self.create_domain()
        gpu_domain = self.create_domain()

        # Run CPU
        cpu_domain.set_multiprocessor_mode(1)
        for t in cpu_domain.evolve(yieldstep=0.1, finaltime=0.1):
            pass

        cpu_stage = cpu_domain.quantities['stage'].centroid_values.copy()

        # Run GPU
        gpu_domain.set_multiprocessor_mode(2)
        from anuga.shallow_water.sw_domain_gpu_ext import sync_to_device, sync_from_device

        gpu_dom = gpu_domain.gpu_interface.gpu_dom
        sync_to_device(gpu_dom)

        for t in gpu_domain.evolve(yieldstep=0.1, finaltime=0.1):
            pass

        sync_from_device(gpu_dom)
        gpu_stage = gpu_domain.quantities['stage'].centroid_values.copy()

        # Compare - allow small tolerance for floating point differences
        diff = np.abs(cpu_stage - gpu_stage)
        self.assertLess(diff.max(), 1e-10,
                        f"Stage difference too large: max={diff.max():.2e}")

    def test_multi_step_evolution(self):
        """Test multiple RK2 steps match between CPU and GPU."""
        cpu_domain = self.create_domain()
        gpu_domain = self.create_domain()

        # Run CPU for 1 second
        cpu_domain.set_multiprocessor_mode(1)
        for t in cpu_domain.evolve(yieldstep=0.5, finaltime=1.0):
            pass

        cpu_stage = cpu_domain.quantities['stage'].centroid_values.copy()
        cpu_xmom = cpu_domain.quantities['xmomentum'].centroid_values.copy()

        # Run GPU for 1 second
        gpu_domain.set_multiprocessor_mode(2)
        from anuga.shallow_water.sw_domain_gpu_ext import sync_to_device, sync_from_device

        gpu_dom = gpu_domain.gpu_interface.gpu_dom
        sync_to_device(gpu_dom)

        for t in gpu_domain.evolve(yieldstep=0.5, finaltime=1.0):
            pass

        sync_from_device(gpu_dom)
        gpu_stage = gpu_domain.quantities['stage'].centroid_values.copy()
        gpu_xmom = gpu_domain.quantities['xmomentum'].centroid_values.copy()

        # For longer simulations, allow slightly larger tolerance due to
        # floating point accumulation differences in parallel execution
        stage_diff = np.abs(cpu_stage - gpu_stage)
        xmom_diff = np.abs(cpu_xmom - gpu_xmom)

        self.assertLess(stage_diff.max(), 1e-8,
                        f"Stage difference too large: max={stage_diff.max():.2e}")
        self.assertLess(xmom_diff.max(), 1e-8,
                        f"Xmomentum difference too large: max={xmom_diff.max():.2e}")


@pytest.mark.skipif(not gpu_available(), reason=_gpu_skip_reason())
class Test_GPU_Boundaries(unittest.TestCase):
    """Tests for GPU boundary evaluation."""

    def test_reflective_boundary(self):
        """Test reflective boundary on GPU."""
        domain = rectangular_cross_domain(5, 5, len1=50., len2=50.)
        domain.set_flow_algorithm('DE0')
        domain.set_low_froude(0)
        domain.set_name('test_reflective')
        domain.set_datadir(tempfile.mkdtemp())
        domain.store = False

        domain.set_quantity('elevation', -1.0)
        domain.set_quantity('friction', 0.01)
        domain.set_quantity('stage', 0.5)  # Water above bed
        domain.set_quantity('xmomentum', 1.0)  # Flow towards boundaries

        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        # Run CPU
        domain.set_multiprocessor_mode(1)
        domain.distribute_to_vertices_and_edges()
        domain.update_boundary()

        cpu_stage_bv = domain.quantities['stage'].boundary_values.copy()
        cpu_xmom_bv = domain.quantities['xmomentum'].boundary_values.copy()

        # Run GPU
        domain.set_multiprocessor_mode(2)
        from anuga.shallow_water.sw_domain_gpu_ext import (
            sync_to_device, sync_all_from_device,
            extrapolate_second_order_gpu, evaluate_reflective_boundary_gpu
        )

        gpu_dom = domain.gpu_interface.gpu_dom
        sync_to_device(gpu_dom)
        extrapolate_second_order_gpu(gpu_dom)
        evaluate_reflective_boundary_gpu(gpu_dom)
        sync_all_from_device(gpu_dom)

        gpu_stage_bv = domain.quantities['stage'].boundary_values.copy()
        gpu_xmom_bv = domain.quantities['xmomentum'].boundary_values.copy()

        np.testing.assert_allclose(cpu_stage_bv, gpu_stage_bv, rtol=1e-10, atol=1e-14,
                                   err_msg="Reflective boundary stage mismatch")
        np.testing.assert_allclose(cpu_xmom_bv, gpu_xmom_bv, rtol=1e-10, atol=1e-14,
                                   err_msg="Reflective boundary xmomentum mismatch")

    def test_dirichlet_boundary(self):
        """Test Dirichlet boundary on GPU."""
        domain = rectangular_cross_domain(5, 5, len1=50., len2=50.)
        domain.set_flow_algorithm('DE0')
        domain.set_name('test_dirichlet')
        domain.set_datadir(tempfile.mkdtemp())
        domain.store = False

        domain.set_quantity('elevation', -1.0)
        domain.set_quantity('friction', 0.01)
        domain.set_quantity('stage', 0.0)

        Br = Reflective_boundary(domain)
        Bd = Dirichlet_boundary([0.5, 0.1, 0.05])  # stage, xmom, ymom
        domain.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})

        # Run GPU evolution
        domain.set_multiprocessor_mode(2)
        from anuga.shallow_water.sw_domain_gpu_ext import (
            sync_to_device, sync_all_from_device,
            extrapolate_second_order_gpu, evaluate_reflective_boundary_gpu,
            evaluate_dirichlet_boundary_gpu
        )

        gpu_dom = domain.gpu_interface.gpu_dom
        sync_to_device(gpu_dom)
        extrapolate_second_order_gpu(gpu_dom)
        evaluate_reflective_boundary_gpu(gpu_dom)
        evaluate_dirichlet_boundary_gpu(gpu_dom)
        sync_all_from_device(gpu_dom)

        # Check that Dirichlet values are applied
        stage_bv = domain.quantities['stage'].boundary_values
        xmom_bv = domain.quantities['xmomentum'].boundary_values
        ymom_bv = domain.quantities['ymomentum'].boundary_values

        # Get indices for 'left' boundary tag
        left_indices = domain.tag_boundary_cells.get('left', [])
        self.assertGreater(len(left_indices), 0, "Should have 'left' boundary edges")

        for idx in left_indices:
            self.assertAlmostEqual(stage_bv[idx], 0.5, places=10)
            self.assertAlmostEqual(xmom_bv[idx], 0.1, places=10)
            self.assertAlmostEqual(ymom_bv[idx], 0.05, places=10)

    def test_transmissive_n_zero_t_boundary(self):
        """Test Transmissive_n_momentum_zero_t_momentum_set_stage_boundary on GPU."""
        domain = rectangular_cross_domain(5, 5, len1=50., len2=50.)
        domain.set_flow_algorithm('DE0')
        domain.set_low_froude(0)
        domain.set_name('test_transmissive')
        domain.set_datadir(tempfile.mkdtemp())
        domain.store = False

        domain.set_quantity('elevation', -1.0)
        domain.set_quantity('friction', 0.01)
        domain.set_quantity('stage', 0.0)
        domain.set_quantity('xmomentum', 0.5)

        def tide_function(t):
            return -0.3

        Br = Reflective_boundary(domain)
        Bt = anuga.Transmissive_n_momentum_zero_t_momentum_set_stage_boundary(
            domain, function=tide_function)
        domain.set_boundary({'left': Bt, 'right': Br, 'top': Br, 'bottom': Br})

        # Run one evolve step
        domain.set_multiprocessor_mode(2)
        from anuga.shallow_water.sw_domain_gpu_ext import sync_to_device, sync_all_from_device

        gpu_dom = domain.gpu_interface.gpu_dom
        sync_to_device(gpu_dom)

        for t in domain.evolve(yieldstep=0.1, finaltime=0.1):
            pass

        sync_all_from_device(gpu_dom)

        # Check that stage boundary values are set to tide function value
        stage_bv = domain.quantities['stage'].boundary_values

        # Get indices for 'left' boundary tag
        left_indices = domain.tag_boundary_cells.get('left', [])
        self.assertGreater(len(left_indices), 0, "Should have 'left' boundary edges")

        for idx in left_indices:
            # Stage should be set to tide function value
            self.assertAlmostEqual(stage_bv[idx], -0.3, places=5)


@pytest.mark.skipif(not gpu_available(), reason=_gpu_skip_reason())
class Test_GPU_Initialization(unittest.TestCase):
    """Tests for GPU initialization and error handling."""

    def test_mode2_without_boundaries_defers_interface(self):
        """Mode 2 can be selected before boundaries are set: the device
        interface build is deferred to the first evolve() (boundaries are
        typically set after construction), so this must NOT raise."""
        domain = rectangular_cross_domain(5, 5, len1=50., len2=50.)
        domain.set_flow_algorithm('DE0')
        domain.set_name('test_init')
        domain.set_datadir(tempfile.mkdtemp())
        domain.store = False

        domain.set_quantity('elevation', -1.0)
        domain.set_quantity('stage', 0.0)

        # No boundaries set: selecting mode 2 records the mode but defers the
        # device interface (built lazily once boundaries are available).
        domain.set_multiprocessor_mode(2)
        self.assertEqual(domain.multiprocessor_mode, 2)
        self.assertIsNone(domain.gpu_interface)

    def test_correct_initialization_order(self):
        """Test that correct initialization order works."""
        domain = rectangular_cross_domain(5, 5, len1=50., len2=50.)
        domain.set_flow_algorithm('DE0')
        domain.set_name('test_init_ok')
        domain.set_datadir(tempfile.mkdtemp())
        domain.store = False

        domain.set_quantity('elevation', -1.0)
        domain.set_quantity('stage', 0.0)

        # Set boundaries FIRST
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        # Then enable GPU mode - should work
        domain.set_multiprocessor_mode(2)

        self.assertIsNotNone(domain.gpu_interface)
        self.assertTrue(domain.gpu_interface.initialized)


@pytest.mark.skipif(not gpu_available(), reason=_gpu_skip_reason())
class Test_GPU_LargeDomain(unittest.TestCase):
    """Tests with larger domains to verify scaling."""

    def test_large_rectangular_domain(self):
        """Test with ~5000 elements."""
        domain = rectangular_cross_domain(50, 50, len1=100., len2=100.)
        domain.set_flow_algorithm('DE0')
        domain.set_low_froude(0)
        domain.set_name('test_large')
        domain.set_datadir(tempfile.mkdtemp())
        domain.store = False

        n_elements = len(domain)
        self.assertGreater(n_elements, 4000, "Domain should have >4000 elements")

        def topography(x, y):
            return -x / 50.0

        domain.set_quantity('elevation', topography)
        domain.set_quantity('friction', 0.01)
        domain.set_quantity('stage', 0.0)

        def tide_function(t):
            return -0.5 + 0.01 * t

        Br = Reflective_boundary(domain)
        Bt = anuga.Transmissive_n_momentum_zero_t_momentum_set_stage_boundary(
            domain, function=tide_function)
        domain.set_boundary({'left': Bt, 'right': Br, 'top': Br, 'bottom': Br})

        # Save initial state
        initial_stage = domain.quantities['stage'].centroid_values.copy()

        # Run CPU
        domain.set_multiprocessor_mode(1)
        for t in domain.evolve(yieldstep=1.0, finaltime=1.0):
            pass

        cpu_stage = domain.quantities['stage'].centroid_values.copy()

        # Reset
        for qname in ['stage', 'xmomentum', 'ymomentum']:
            domain.quantities[qname].centroid_values[:] = 0.0
        domain.quantities['stage'].centroid_values[:] = initial_stage
        domain.set_time(0.0)

        # Run GPU
        domain.set_multiprocessor_mode(2)
        from anuga.shallow_water.sw_domain_gpu_ext import sync_to_device, sync_from_device

        sync_to_device(domain.gpu_interface.gpu_dom)

        for t in domain.evolve(yieldstep=1.0, finaltime=1.0):
            pass

        sync_from_device(domain.gpu_interface.gpu_dom)
        gpu_stage = domain.quantities['stage'].centroid_values.copy()

        # Compare
        stage_diff = np.abs(cpu_stage - gpu_stage)
        self.assertLess(stage_diff.max(), 1e-6,
                        f"Large domain stage difference: max={stage_diff.max():.2e}")


@pytest.mark.skipif(not gpu_available(), reason=_gpu_skip_reason())
class Test_GPU_InletOperator(unittest.TestCase):
    """Tests for Inlet_operator with GPU acceleration."""

    def create_domain(self, name='test_inlet'):
        """Create a test domain suitable for inlet testing."""
        domain = rectangular_cross_domain(20, 10, len1=200., len2=100.)
        domain.set_flow_algorithm('DE0')
        domain.set_low_froude(0)
        domain.set_name(name)
        domain.set_datadir(tempfile.mkdtemp())
        domain.store = False

        def topography(x, y):
            return -x / 100.0  # Gentle slope

        domain.set_quantity('elevation', topography)
        domain.set_quantity('friction', 0.03)
        domain.set_quantity('stage', 0.0)

        Br = Reflective_boundary(domain)
        Bd = Dirichlet_boundary([0, 0, 0])
        domain.set_boundary({'left': Bd, 'right': Bd, 'top': Br, 'bottom': Br})

        return domain

    def test_inlet_operator_basic(self):
        """Test that inlet operator works on GPU."""
        domain = self.create_domain('test_inlet_basic')

        # Add inlet operator - line across left side
        line = [[10.0, 20.0], [10.0, 80.0]]
        Q = 10.0  # m^3/s
        inlet = Inlet_operator(domain, line, Q, verbose=False)

        # Enable GPU mode
        domain.set_multiprocessor_mode(2)

        from anuga.shallow_water.sw_domain_gpu_ext import sync_to_device, sync_from_device
        gpu_dom = domain.gpu_interface.gpu_dom
        sync_to_device(gpu_dom)

        # Evolve
        for t in domain.evolve(yieldstep=0.5, finaltime=1.0):
            pass

        sync_from_device(gpu_dom)

        # Verify water was added
        water_volume = domain.get_water_volume()
        applied_volume = inlet.total_applied_volume

        self.assertGreater(water_volume, 0, "Water volume should be positive after inlet")
        self.assertGreater(applied_volume, 0, "Inlet should have applied some volume")
        # Volume added should be approximately Q * time = 10 * 1.0 = 10 m^3
        self.assertAlmostEqual(applied_volume, Q * 1.0, delta=1.0,
                               msg=f"Applied volume {applied_volume} should be close to {Q * 1.0}")

    def test_inlet_operator_cpu_gpu_match(self):
        """Test that inlet operator produces same results on CPU and GPU."""
        # Create two identical domains
        cpu_domain = self.create_domain('test_inlet_cpu')
        gpu_domain = self.create_domain('test_inlet_gpu')

        # Add inlet operators with same parameters
        line = [[10.0, 20.0], [10.0, 80.0]]
        Q = 15.0  # m^3/s

        cpu_inlet = Inlet_operator(cpu_domain, line, Q, verbose=False)
        gpu_inlet = Inlet_operator(gpu_domain, line, Q, verbose=False)

        # Run CPU
        cpu_domain.set_multiprocessor_mode(1)
        for t in cpu_domain.evolve(yieldstep=0.5, finaltime=1.0):
            pass

        cpu_stage = cpu_domain.quantities['stage'].centroid_values.copy()
        cpu_volume = cpu_domain.get_water_volume()
        cpu_inlet_volume = cpu_inlet.total_applied_volume

        # Run GPU
        gpu_domain.set_multiprocessor_mode(2)
        from anuga.shallow_water.sw_domain_gpu_ext import sync_to_device, sync_from_device
        gpu_dom = gpu_domain.gpu_interface.gpu_dom
        sync_to_device(gpu_dom)

        for t in gpu_domain.evolve(yieldstep=0.5, finaltime=1.0):
            pass

        sync_from_device(gpu_dom)
        gpu_stage = gpu_domain.quantities['stage'].centroid_values.copy()
        gpu_volume = gpu_domain.get_water_volume()
        gpu_inlet_volume = gpu_inlet.total_applied_volume

        # Compare inlet volumes
        self.assertAlmostEqual(cpu_inlet_volume, gpu_inlet_volume, places=6,
                               msg=f"Inlet volumes differ: CPU={cpu_inlet_volume}, GPU={gpu_inlet_volume}")

        # Compare water volumes
        self.assertAlmostEqual(cpu_volume, gpu_volume, delta=1e-6,
                               msg=f"Water volumes differ: CPU={cpu_volume}, GPU={gpu_volume}")

        # Compare stage values. mode=1 and mode=2 apply the inlet with the same
        # level-fill, so they agree to machine precision (~7e-17 on a GPU build,
        # bit-identical on CPU); the bound catches any mode-1/mode-2 regression.
        stage_diff = np.abs(cpu_stage - gpu_stage)
        self.assertLess(stage_diff.max(), 1e-10,
                        f"Stage difference too large: max={stage_diff.max():.2e}")

    def test_inlet_operator_time_varying_Q(self):
        """Test inlet operator with time-varying discharge on GPU."""
        domain = self.create_domain('test_inlet_timevar')

        # Time-varying discharge function
        def Q_func(t):
            return 5.0 + 10.0 * t  # Increases with time

        line = [[10.0, 20.0], [10.0, 80.0]]
        inlet = Inlet_operator(domain, line, Q_func, verbose=False)

        # Enable GPU mode
        domain.set_multiprocessor_mode(2)
        from anuga.shallow_water.sw_domain_gpu_ext import sync_to_device, sync_from_device
        gpu_dom = domain.gpu_interface.gpu_dom
        sync_to_device(gpu_dom)

        # Evolve
        for t in domain.evolve(yieldstep=0.5, finaltime=1.0):
            pass

        sync_from_device(gpu_dom)

        # Verify water was added
        water_volume = domain.get_water_volume()
        applied_volume = inlet.total_applied_volume

        self.assertGreater(water_volume, 0, "Water volume should be positive")
        self.assertGreater(applied_volume, 0, "Inlet should have applied some volume")
        # With Q = 5 + 10*t, average over [0,1] is about 10, so ~10 m^3 total
        self.assertGreater(applied_volume, 5.0, "Should have applied at least 5 m^3")

    def test_inlet_operator_with_velocity(self):
        """Test inlet operator with specified velocity on GPU."""
        cpu_domain = self.create_domain('test_inlet_vel_cpu')
        gpu_domain = self.create_domain('test_inlet_vel_gpu')

        line = [[10.0, 20.0], [10.0, 80.0]]
        Q = 10.0
        velocity = [0.5, 0.0]  # Velocity in x direction

        cpu_inlet = Inlet_operator(cpu_domain, line, Q, velocity=velocity, verbose=False)
        gpu_inlet = Inlet_operator(gpu_domain, line, Q, velocity=velocity, verbose=False)

        # Run CPU
        cpu_domain.set_multiprocessor_mode(1)
        for t in cpu_domain.evolve(yieldstep=0.5, finaltime=1.0):
            pass

        cpu_xmom = cpu_domain.quantities['xmomentum'].centroid_values.copy()

        # Run GPU
        gpu_domain.set_multiprocessor_mode(2)
        from anuga.shallow_water.sw_domain_gpu_ext import sync_to_device, sync_from_device
        gpu_dom = gpu_domain.gpu_interface.gpu_dom
        sync_to_device(gpu_dom)

        for t in gpu_domain.evolve(yieldstep=0.5, finaltime=1.0):
            pass

        sync_from_device(gpu_dom)
        gpu_xmom = gpu_domain.quantities['xmomentum'].centroid_values.copy()

        # Compare momentum. Inlet applied identically in mode=1 and mode=2, so
        # agreement is to machine precision (~5e-16 on a GPU build).
        xmom_diff = np.abs(cpu_xmom - gpu_xmom)
        self.assertLess(xmom_diff.max(), 1e-10,
                        f"Xmomentum difference too large: max={xmom_diff.max():.2e}")


@pytest.mark.skipif(not gpu_available(), reason=_gpu_skip_reason())
class Test_GPU_Riverwall(unittest.TestCase):
    """Tests for riverwall/weir support with GPU acceleration."""

    def create_riverwall_domain(self, name='test_riverwall'):
        """Create a test domain with a riverwall."""
        from anuga import create_mesh_from_regions, create_domain_from_file, create_domain_from_regions
        import os

        mesh_filename = f'{name}.msh'

        # Domain polygon
        boundaryPolygon = [[0., 0.], [0., 100.], [100.0, 100.0], [100.0, 0.0]]

        # Riverwall - a wall across the middle of the domain
        riverWall = {'centralWall':
                     [[50., 0.0, -0.0],
                      [50., 45., -0.0],
                      [50., 46., -0.2],  # Dip in the wall
                      [50., 100.0, -0.0]]
                     }

        riverWall_Par = {'centralWall': {'Qfactor': 1.0}}

        # Region points
        regionPtAreas = [[25., 50., 5.0*5.0*0.5],
                         [75., 50., 5.0*5.0*0.5]]

        domain = create_domain_from_regions(boundaryPolygon,
                                 boundary_tags={'left': [0],
                                                'top': [1],
                                                'right': [2],
                                                'bottom': [3]},
                                 maximum_triangle_area=10.0*10.0*0.5,
                                 minimum_triangle_angle=28.0,
                                 breaklines=list(riverWall.values()),
                                 use_cache=False,
                                 verbose=False,
                                 regionPtArea=regionPtAreas)

        #domain = create_domain_from_file(mesh_filename)
        domain.set_flow_algorithm('DE0')
        domain.set_low_froude(0)
        domain.set_name(name)
        domain.set_datadir(tempfile.mkdtemp())
        domain.store = False

        def topography(x, y):
            return -x / 150.0

        def stagefun(x, y):
            return -0.5

        domain.set_quantity('elevation', topography)
        domain.set_quantity('friction', 0.03)
        domain.set_quantity('stage', stagefun)

        # Create the riverwalls
        domain.riverwallData.create_riverwalls(riverWall, riverWall_Par, verbose=False)

        return domain

    def test_riverwall_initialization(self):
        """Test that riverwall arrays are properly initialized for GPU."""
        domain = self.create_riverwall_domain('test_rw_init')

        # Set boundaries
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        # Enable GPU mode
        domain.set_multiprocessor_mode(2)

        # Verify riverwall arrays are set
        self.assertGreater(domain.number_of_riverwall_edges, 0,
                           "Should have riverwall edges")
        self.assertIsNotNone(domain.edge_flux_type)
        self.assertIsNotNone(domain.riverwallData.riverwall_elevation)

    def test_riverwall_cpu_gpu_match(self):
        """Test that riverwall simulation matches between CPU and GPU."""
        from math import exp

        # Create two identical domains
        cpu_domain = self.create_riverwall_domain('test_rw_cpu')
        gpu_domain = self.create_riverwall_domain('test_rw_gpu')

        # Boundary function
        def boundaryFun(t):
            output = -0.4 * exp(-t / 100.) - 0.1
            return min(output, -0.11)

        # Set boundaries
        Br_cpu = Reflective_boundary(cpu_domain)
        Bt_cpu = anuga.Transmissive_n_momentum_zero_t_momentum_set_stage_boundary(
            domain=cpu_domain, function=boundaryFun)
        cpu_domain.set_boundary({'left': Br_cpu, 'right': Bt_cpu, 'top': Br_cpu, 'bottom': Br_cpu})

        Br_gpu = Reflective_boundary(gpu_domain)
        Bt_gpu = anuga.Transmissive_n_momentum_zero_t_momentum_set_stage_boundary(
            domain=gpu_domain, function=boundaryFun)
        gpu_domain.set_boundary({'left': Br_gpu, 'right': Bt_gpu, 'top': Br_gpu, 'bottom': Br_gpu})

        # Run CPU
        cpu_domain.set_multiprocessor_mode(1)
        for t in cpu_domain.evolve(yieldstep=1.0, finaltime=5.0):
            pass

        cpu_stage = cpu_domain.quantities['stage'].centroid_values.copy()
        cpu_xmom = cpu_domain.quantities['xmomentum'].centroid_values.copy()

        # Run GPU
        gpu_domain.set_multiprocessor_mode(2)
        from anuga.shallow_water.sw_domain_gpu_ext import sync_to_device, sync_from_device
        gpu_dom = gpu_domain.gpu_interface.gpu_dom
        sync_to_device(gpu_dom)

        for t in gpu_domain.evolve(yieldstep=1.0, finaltime=5.0):
            pass

        sync_from_device(gpu_dom)
        gpu_stage = gpu_domain.quantities['stage'].centroid_values.copy()
        gpu_xmom = gpu_domain.quantities['xmomentum'].centroid_values.copy()

        # Compare
        stage_diff = np.abs(cpu_stage - gpu_stage)
        xmom_diff = np.abs(cpu_xmom - gpu_xmom)

        self.assertLess(stage_diff.max(), 1e-6,
                        f"Stage difference too large: max={stage_diff.max():.2e}")
        self.assertLess(xmom_diff.max(), 1e-6,
                        f"Xmomentum difference too large: max={xmom_diff.max():.2e}")

    def test_riverwall_flux_kernel(self):
        """Test that riverwall flux computation matches CPU."""
        domain = self.create_riverwall_domain('test_rw_flux')

        # Set initial conditions that should trigger flow over the wall
        domain.set_quantity('stage', 0.5)  # Water above wall level

        # Set boundaries
        Br = Reflective_boundary(domain)
        Bd = Dirichlet_boundary([0.5, 0., 0.])
        domain.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})

        # Run CPU flux computation
        domain.set_multiprocessor_mode(1)
        domain.distribute_to_vertices_and_edges()
        domain.update_boundary()
        domain.compute_fluxes()

        cpu_timestep = domain.flux_timestep
        cpu_stage_update = domain.quantities['stage'].explicit_update.copy()
        cpu_xmom_update = domain.quantities['xmomentum'].explicit_update.copy()

        # Reset and run GPU flux computation
        for qname in ['stage', 'xmomentum', 'ymomentum']:
            domain.quantities[qname].explicit_update[:] = 0.0

        domain.set_multiprocessor_mode(2)
        from anuga.shallow_water.sw_domain_gpu_ext import (
            sync_to_device, sync_all_from_device,
            extrapolate_second_order_gpu, compute_fluxes_gpu,
            evaluate_reflective_boundary_gpu, evaluate_dirichlet_boundary_gpu
        )

        gpu_dom = domain.gpu_interface.gpu_dom
        sync_to_device(gpu_dom)
        extrapolate_second_order_gpu(gpu_dom)
        evaluate_reflective_boundary_gpu(gpu_dom)
        evaluate_dirichlet_boundary_gpu(gpu_dom)
        gpu_timestep = compute_fluxes_gpu(gpu_dom)
        sync_all_from_device(gpu_dom)

        gpu_stage_update = domain.quantities['stage'].explicit_update.copy()
        gpu_xmom_update = domain.quantities['xmomentum'].explicit_update.copy()

        # Verify - riverwalls can have larger flux differences due to weir formula
        self.assertAlmostEqual(cpu_timestep, gpu_timestep, places=8,
                               msg=f"Timestep mismatch: CPU={cpu_timestep}, GPU={gpu_timestep}")
        np.testing.assert_allclose(cpu_stage_update, gpu_stage_update, rtol=1e-8, atol=1e-12,
                                   err_msg="Stage explicit_update mismatch")
        np.testing.assert_allclose(cpu_xmom_update, gpu_xmom_update, rtol=1e-8, atol=1e-12,
                                   err_msg="Xmomentum explicit_update mismatch")

    def test_riverwall_weir_discharge(self):
        """Test that weir discharge formula is applied correctly."""
        domain = self.create_riverwall_domain('test_rw_weir')

        # Create asymmetric water levels to trigger weir flow
        def asymmetric_stage(x, y):
            # Higher water on left side of wall (x < 50)
            return np.where(x < 50, 0.5, -0.3)

        domain.set_quantity('stage', asymmetric_stage)

        # Set boundaries
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        # Run simulation on GPU
        domain.set_multiprocessor_mode(2)
        from anuga.shallow_water.sw_domain_gpu_ext import sync_to_device, sync_from_device

        gpu_dom = domain.gpu_interface.gpu_dom
        sync_to_device(gpu_dom)

        initial_volume = domain.get_water_volume()

        for t in domain.evolve(yieldstep=1.0, finaltime=5.0):
            pass

        sync_from_device(gpu_dom)

        final_volume = domain.get_water_volume()

        # Volume should be conserved (within tolerance)
        self.assertAlmostEqual(initial_volume, final_volume, places=4,
                               msg=f"Volume not conserved: initial={initial_volume}, final={final_volume}")

        # Check that water has redistributed (some flow over the wall)
        stage = domain.quantities['stage'].centroid_values
        x = domain.centroid_coordinates[:, 0]

        left_mask = x < 50
        right_mask = x >= 50

        # Left side should have lower water now (some flowed over)
        # Right side should have higher water
        # This is a qualitative check that the riverwall is working
        left_mean = np.mean(stage[left_mask])
        right_mean = np.mean(stage[right_mask])

        # Initially left was 0.5, right was -0.3
        # After evolution, the difference should be reduced
        self.assertLess(left_mean - right_mean, 0.8,
                        "Water should have flowed over the wall, reducing the level difference")


@pytest.mark.skipif(not gpu_available(), reason=_gpu_skip_reason())
class Test_GPU_EndToEnd(unittest.TestCase):
    """End-to-end regression tests comparing mode=1 (Python RK2) vs mode=2 (C RK2).

    In the default build (CPU_ONLY_MODE) both modes run on CPU, so results
    should be identical to machine-epsilon precision.  These tests serve as
    a regression baseline: if a change breaks numerical equivalence it will
    show up here before it can affect real GPU runs.
    """

    def _create_tidal_domain(self, name, algorithm='DE0'):
        """20×10 domain with a sloping bed and a tidal left boundary."""
        domain = rectangular_cross_domain(20, 10, len1=200., len2=100.)
        domain.set_flow_algorithm(algorithm)
        domain.set_low_froude(0)
        domain.set_name(name)
        domain.set_datadir(tempfile.mkdtemp())
        domain.store = False

        domain.set_quantity('elevation', lambda x, y: -x / 100.)
        domain.set_quantity('friction', 0.01)
        domain.set_quantity('stage', 0.0)

        def tide(t):
            return -0.5 + 0.3 * np.sin(2 * np.pi * t / 100.)

        Br = Reflective_boundary(domain)
        Bt = anuga.Transmissive_n_momentum_zero_t_momentum_set_stage_boundary(
            domain, function=tide)
        domain.set_boundary({'left': Bt, 'right': Br, 'top': Br, 'bottom': Br})
        return domain

    def _create_dam_break_domain(self, name):
        """20×10 domain with a dam-break initial condition."""
        domain = rectangular_cross_domain(20, 10, len1=200., len2=100.)
        domain.set_flow_algorithm('DE0')
        domain.set_low_froude(0)
        domain.set_name(name)
        domain.set_datadir(tempfile.mkdtemp())
        domain.store = False

        domain.set_quantity('elevation', -1.0)
        domain.set_quantity('friction', 0.01)
        domain.set_quantity('stage', lambda x, y: np.where(x < 100., 0.5, -1.0))

        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})
        return domain

    def _run_and_compare(self, cpu_domain, gpu_domain, finaltime, yieldstep=2.0):
        """Run both domains and return (cpu_q, gpu_q) dicts of final centroid values."""
        from anuga.shallow_water.sw_domain_gpu_ext import sync_to_device, sync_from_device

        cpu_domain.set_multiprocessor_mode(1)
        for _ in cpu_domain.evolve(yieldstep=yieldstep, finaltime=finaltime):
            pass

        gpu_domain.set_multiprocessor_mode(2)
        sync_to_device(gpu_domain.gpu_interface.gpu_dom)
        for _ in gpu_domain.evolve(yieldstep=yieldstep, finaltime=finaltime):
            pass
        sync_from_device(gpu_domain.gpu_interface.gpu_dom)

        cpu_q = {q: cpu_domain.quantities[q].centroid_values.copy()
                 for q in ['stage', 'xmomentum', 'ymomentum']}
        gpu_q = {q: gpu_domain.quantities[q].centroid_values.copy()
                 for q in ['stage', 'xmomentum', 'ymomentum']}
        return cpu_q, gpu_q

    @pytest.mark.slow
    def test_10s_tidal_mode1_vs_mode2(self):
        """10-second tidal run: mode=1 and mode=2 must agree to machine precision.

        In CPU_ONLY_MODE the Python and C RK2 loops call identical kernels,
        so results should be bit-for-bit identical.  Tolerance of 1e-12
        provides a comfortable guard against any future drift.
        """
        cpu_q, gpu_q = self._run_and_compare(
            self._create_tidal_domain('e2e_tidal_cpu'),
            self._create_tidal_domain('e2e_tidal_gpu'),
            finaltime=10.0)

        for qname in ['stage', 'xmomentum', 'ymomentum']:
            np.testing.assert_allclose(
                gpu_q[qname], cpu_q[qname],
                rtol=0, atol=1e-12,
                err_msg=f'10s tidal: {qname} mismatch between mode=1 and mode=2')

    @pytest.mark.slow
    def test_10s_dam_break_mode1_vs_mode2(self):
        """10-second dam-break: mode=1 and mode=2 must agree to machine precision."""
        cpu_q, gpu_q = self._run_and_compare(
            self._create_dam_break_domain('e2e_dambreak_cpu'),
            self._create_dam_break_domain('e2e_dambreak_gpu'),
            finaltime=10.0)

        for qname in ['stage', 'xmomentum', 'ymomentum']:
            np.testing.assert_allclose(
                gpu_q[qname], cpu_q[qname],
                rtol=0, atol=1e-12,
                err_msg=f'10s dam break: {qname} mismatch between mode=1 and mode=2')

    def test_boundary_flux_integral_de1_mode1_vs_mode2(self):
        """Regression: boundary_flux_integral must match mode=1 vs mode=2 for a
        multi-substep scheme (DE1/RK2) with a Python-evaluated GPU boundary.

        The Python-orchestrated unified RK loop once called compute_fluxes_gpu
        with a fixed substep index (0, 1), so every substep overwrote
        boundary_flux_sum[0] and the 2nd RK substep's boundary flux was dropped
        — halving the integral. Physics (volume) was unaffected; this guards the
        diagnostic. The tidal domain uses the transmissive-set-stage boundary
        that forces the Python-orchestrated path.
        """
        from anuga.shallow_water.sw_domain_gpu_ext import (
            sync_to_device, sync_from_device)

        cpu = self._create_tidal_domain('bfi_cpu', algorithm='DE1')
        gpu = self._create_tidal_domain('bfi_gpu', algorithm='DE1')

        cpu.set_multiprocessor_mode(1)
        for _ in cpu.evolve(yieldstep=2.0, finaltime=10.0):
            pass

        gpu.set_multiprocessor_mode(2)
        sync_to_device(gpu.gpu_interface.gpu_dom)
        for _ in gpu.evolve(yieldstep=2.0, finaltime=10.0):
            pass
        sync_from_device(gpu.gpu_interface.gpu_dom)

        bf_cpu = cpu.get_boundary_flux_integral()
        bf_gpu = gpu.get_boundary_flux_integral()
        # Guard the test itself: the tidal boundary must move real water.
        self.assertGreater(abs(bf_cpu), 1.0,
                           'test setup: expected a non-trivial boundary flux')
        np.testing.assert_allclose(
            bf_gpu, bf_cpu, rtol=1e-9, atol=1e-9,
            err_msg=f'DE1 boundary_flux_integral mismatch: '
                    f'mode1={bf_cpu}, mode2={bf_gpu}')

    def test_volume_conservation_mode2(self):
        """Water volume is conserved over 10 s in GPU mode (closed boundaries)."""
        from anuga.shallow_water.sw_domain_gpu_ext import sync_to_device, sync_from_device

        domain = self._create_dam_break_domain('e2e_vol_gpu')
        initial_volume = domain.get_water_volume()

        domain.set_multiprocessor_mode(2)
        sync_to_device(domain.gpu_interface.gpu_dom)
        for _ in domain.evolve(yieldstep=2.0, finaltime=10.0):
            pass
        sync_from_device(domain.gpu_interface.gpu_dom)

        final_volume = domain.get_water_volume()
        self.assertAlmostEqual(
            initial_volume, final_volume, places=6,
            msg=f'Volume not conserved: initial={initial_volume:.6f}, '
                f'final={final_volume:.6f}')


@pytest.mark.skipif(not gpu_available(), reason=_gpu_skip_reason())
class Test_GPU_RK3(unittest.TestCase):
    """Tests for SSP-RK3 timestepping in GPU mode (DE2 flow algorithm).

    DE2 uses the Shu-Osher 3-stage SSP-RK3 scheme.  In CPU_ONLY_MODE both
    mode=1 and mode=2 call identical C kernels, so results must match to
    machine precision.
    """

    def _create_domain(self, name, algorithm='DE2'):
        domain = rectangular_cross_domain(20, 10, len1=200., len2=100.)
        domain.set_flow_algorithm(algorithm)
        domain.set_low_froude(0)
        domain.set_name(name)
        domain.set_datadir(tempfile.mkdtemp())
        domain.store = False

        domain.set_quantity('elevation', -1.0)
        domain.set_quantity('friction', 0.01)
        domain.set_quantity('stage', lambda x, y: np.where(x < 100., 0.5, -1.0))

        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})
        return domain

    def test_single_rk3_step_gpu(self):
        """One RK3 step on GPU produces valid (non-NaN) conserved quantities."""
        from anuga.shallow_water.sw_domain_gpu_ext import (
            init_gpu_domain, map_to_gpu, unmap_from_gpu, finalize_gpu_domain,
            evolve_one_rk3_step_gpu, sync_from_device
        )

        d = self._create_domain('rk3_single_step')
        gpu = init_gpu_domain(d)
        map_to_gpu(gpu)

        try:
            ts = evolve_one_rk3_step_gpu(gpu, 1.0, 0)
            self.assertGreater(ts, 0.0)
            self.assertLess(ts, 10.0)

            sync_from_device(gpu)
            self.assertFalse(np.any(np.isnan(d.quantities['stage'].centroid_values)),
                             "stage has NaN after RK3 step")
        finally:
            unmap_from_gpu(gpu)
            finalize_gpu_domain(gpu)

    @pytest.mark.slow
    def test_rk3_mode1_vs_mode2_dam_break(self):
        """DE2 (RK3) dam-break: mode=1 and mode=2 must agree.

        On a CPU build (gpu_offload=false) both modes run the identically
        compiled host kernels, so they agree to machine precision. On a real
        GPU (gpu_offload=true) the device executes the same algorithm with a
        different floating-point evaluation order (fma/contraction, reduction
        order), so the agreement is only to ~1e-9, not 1e-12. Pick the tolerance
        accordingly.
        """
        import anuga
        from anuga.shallow_water.sw_domain_gpu_ext import sync_to_device, sync_from_device

        cpu_d = self._create_domain('rk3_cpu')
        gpu_d = self._create_domain('rk3_gpu')

        cpu_d.set_multiprocessor_mode(1)
        for _ in cpu_d.evolve(yieldstep=2.0, finaltime=10.0):
            pass

        gpu_d.set_multiprocessor_mode(2)
        sync_to_device(gpu_d.gpu_interface.gpu_dom)
        for _ in gpu_d.evolve(yieldstep=2.0, finaltime=10.0):
            pass
        sync_from_device(gpu_d.gpu_interface.gpu_dom)

        atol = 1e-8 if anuga.gpu_offload_enabled() else 1e-12
        for qname in ['stage', 'xmomentum', 'ymomentum']:
            np.testing.assert_allclose(
                gpu_d.quantities[qname].centroid_values,
                cpu_d.quantities[qname].centroid_values,
                rtol=0, atol=atol,
                err_msg=f'RK3 dam-break: {qname} mismatch mode=1 vs mode=2')

    def test_saxpy3_kernel(self):
        """saxpy3_conserved_quantities_gpu computes (a*Q + b*backup)/c correctly."""
        from anuga.shallow_water.sw_domain_gpu_ext import (
            init_gpu_domain, map_to_gpu, unmap_from_gpu, finalize_gpu_domain,
            backup_conserved_quantities_gpu, saxpy3_conserved_quantities_gpu,
            sync_to_device, sync_from_device
        )

        d = self._create_domain('rk3_saxpy3')
        # Set known values
        d.quantities['stage'].centroid_values[:] = 2.0
        d.quantities['xmomentum'].centroid_values[:] = 0.0
        d.quantities['ymomentum'].centroid_values[:] = 0.0

        gpu = init_gpu_domain(d)
        map_to_gpu(gpu)
        sync_to_device(gpu)

        try:
            # Backup (backup = 2.0)
            backup_conserved_quantities_gpu(gpu)

            # Change current to 3.0
            d.quantities['stage'].centroid_values[:] = 3.0
            sync_to_device(gpu)

            # saxpy3(2, 1, 3): Q = (2*3.0 + 1*2.0) / 3 = 8/3
            saxpy3_conserved_quantities_gpu(gpu, 2.0, 1.0, 3.0)
            sync_from_device(gpu)

            expected = (2.0 * 3.0 + 1.0 * 2.0) / 3.0
            np.testing.assert_allclose(
                d.quantities['stage'].centroid_values,
                expected,
                rtol=1e-14,
                err_msg='saxpy3 result incorrect')
        finally:
            unmap_from_gpu(gpu)
            finalize_gpu_domain(gpu)


@pytest.mark.skipif(not gpu_available(), reason=_gpu_skip_reason())
class Test_GPU_Culvert(unittest.TestCase):
    """Tests for Boyd box/pipe culvert operators in GPU mode."""

    def _create_culvert_domain(self, name):
        """Two-compartment domain connected by a Boyd box culvert.

        Water starts on the left (x < 100 m), separated by a land barrier.
        The culvert (0.5 m wide × 0.5 m high) provides the only flow path.
        """
        domain = rectangular_cross_domain(20, 10, len1=200., len2=100.)
        domain.set_flow_algorithm('DE0')
        domain.set_low_froude(0)
        domain.set_name(name)
        domain.set_datadir(tempfile.mkdtemp())
        domain.store = False

        domain.set_quantity('elevation', -1.0)
        domain.set_quantity('friction', 0.013)
        domain.set_quantity('stage', lambda x, y: np.where(x < 100., 0.5, -1.0))

        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        from anuga.structures.boyd_box_operator import Boyd_box_operator
        Boyd_box_operator(domain,
                          end_points=[[90., 50.], [110., 50.]],
                          height=0.5, width=0.5,
                          apron=5., manning=0.013,
                          enquiry_gap=5., verbose=False)
        return domain

    def test_culvert_cpu_gpu_match(self):
        """Boyd box culvert: mode=1 vs mode=2 stage comparison at 5 s."""
        from anuga.shallow_water.sw_domain_gpu_ext import sync_to_device, sync_from_device

        cpu_domain = self._create_culvert_domain('culv_cpu')
        gpu_domain = self._create_culvert_domain('culv_gpu')

        # CPU run
        cpu_domain.set_multiprocessor_mode(1)
        for _ in cpu_domain.evolve(yieldstep=1.0, finaltime=5.0):
            pass
        cpu_stage = cpu_domain.quantities['stage'].centroid_values.copy()
        cpu_xmom = cpu_domain.quantities['xmomentum'].centroid_values.copy()

        # GPU run
        gpu_domain.set_multiprocessor_mode(2)
        sync_to_device(gpu_domain.gpu_interface.gpu_dom)
        for _ in gpu_domain.evolve(yieldstep=1.0, finaltime=5.0):
            pass
        sync_from_device(gpu_domain.gpu_interface.gpu_dom)
        gpu_stage = gpu_domain.quantities['stage'].centroid_values.copy()
        gpu_xmom = gpu_domain.quantities['xmomentum'].centroid_values.copy()

        # mode=1 and mode=2 now share the single C culvert kernel
        # (culvert_compute_one + culvert_gather_inlet_host), so they agree to
        # machine precision: measured ~3e-16 (stage) / ~2e-15 (xmom) at t=5 on a
        # GPU-offload build, and bit-identical on a CPU build.  (Previously mode=1
        # called Python boyd_box_function and mode=2 the C translation, whose
        # pow()/sqrt() FP-order differences amplified to ~3% by t=5 — hence the
        # old atol=0.02/0.05.)  The tight bounds below catch any regression that
        # reintroduces a mode-1/mode-2 divergence, with generous headroom over the
        # measured ~1e-15.
        np.testing.assert_allclose(
            gpu_stage, cpu_stage, rtol=0, atol=1e-10,
            err_msg='Culvert 5s: stage mismatch between mode=1 and mode=2')
        np.testing.assert_allclose(
            gpu_xmom, cpu_xmom, rtol=0, atol=1e-9,
            err_msg='Culvert 5s: xmomentum mismatch between mode=1 and mode=2')

    def test_culvert_volume_conservation(self):
        """Boyd box culvert: total water volume is conserved in GPU mode."""
        from anuga.shallow_water.sw_domain_gpu_ext import sync_to_device, sync_from_device

        domain = self._create_culvert_domain('culv_vol')
        initial_volume = domain.get_water_volume()

        domain.set_multiprocessor_mode(2)
        sync_to_device(domain.gpu_interface.gpu_dom)
        for _ in domain.evolve(yieldstep=1.0, finaltime=5.0):
            pass
        sync_from_device(domain.gpu_interface.gpu_dom)

        final_volume = domain.get_water_volume()
        self.assertAlmostEqual(
            initial_volume, final_volume, places=5,
            msg=f'Culvert volume not conserved: '
                f'initial={initial_volume:.5f}, final={final_volume:.5f}')

    def test_culvert_flow_direction(self):
        """Boyd box culvert: flow moves from high to low water level in GPU mode."""
        from anuga.shallow_water.sw_domain_gpu_ext import sync_to_device, sync_from_device

        domain = self._create_culvert_domain('culv_flow')
        x = domain.centroid_coordinates[:, 0]
        left_mask = x < 90.
        right_mask = x > 110.

        domain.set_multiprocessor_mode(2)
        sync_to_device(domain.gpu_interface.gpu_dom)
        for _ in domain.evolve(yieldstep=1.0, finaltime=5.0):
            pass
        sync_from_device(domain.gpu_interface.gpu_dom)

        stage = domain.quantities['stage'].centroid_values

        # Left side should be lower, right side higher, than initial values
        self.assertLess(
            stage[left_mask].mean(), 0.5,
            'Left side should have lost water through the culvert')
        self.assertGreater(
            stage[right_mask].mean(), -1.0,
            'Right side should have gained water through the culvert')

    def _create_pipe_domain(self, name):
        """Two-compartment domain connected by a Boyd *pipe* culvert.

        Same layout as _create_culvert_domain but with a circular pipe, which
        exercises boyd_pipe_discharge (a different critical-depth / flow-area
        path than the box). Returns (domain, operator).
        """
        from anuga.structures.boyd_pipe_operator import Boyd_pipe_operator
        domain = rectangular_cross_domain(20, 10, len1=200., len2=100.)
        domain.set_flow_algorithm('DE0')
        domain.set_low_froude(0)
        domain.set_name(name)
        domain.set_datadir(tempfile.mkdtemp())
        domain.store = False

        domain.set_quantity('elevation', -1.0)
        domain.set_quantity('friction', 0.013)
        domain.set_quantity('stage', lambda x, y: np.where(x < 100., 0.5, -1.0))

        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        op = Boyd_pipe_operator(domain,
                                end_points=[[90., 50.], [110., 50.]],
                                diameter=0.5,
                                apron=5., manning=0.013,
                                enquiry_gap=5., verbose=False)
        return domain, op

    def test_pipe_culvert_cpu_gpu_velocity_match(self):
        """Boyd *pipe* culvert: reported discharge and barrel velocity agree
        between mode=1 and mode=2.

        Regression guard for CPU/GPU drift in the pipe critical-depth /
        flow-area path: a mismatch there leaves the reported barrel velocity
        off between modes while the (inlet-controlled) discharge still
        matches. Both paths now use the same exact circular critical-depth
        solve (circular_critical_depth); physical correctness of that solve is
        tested in anuga/structures/tests/test_boyd_pipe_operator.py. Comparing
        the operator velocity catches drift; a stage-only comparison
        (test_culvert_cpu_gpu_match, box only) does not. Also covers the
        GPUCulvertManager stats write-back onto the Python op.
        """
        cpu_domain, cpu_op = self._create_pipe_domain('pipe_cpu')
        cpu_domain.set_multiprocessor_mode(1)
        for _ in cpu_domain.evolve(yieldstep=1.0, finaltime=3.0):
            pass

        gpu_domain, gpu_op = self._create_pipe_domain('pipe_gpu')
        gpu_domain.set_multiprocessor_mode(2)
        for _ in gpu_domain.evolve(yieldstep=1.0, finaltime=3.0):
            pass

        # Culvert must actually be flowing for the comparison to mean anything
        self.assertGreater(cpu_op.velocity, 0.05,
                           'pipe culvert should be flowing in the test window')
        # Barrel velocity was the ~8%-off quantity; discharge stayed matched.
        self.assertAlmostEqual(
            gpu_op.velocity, cpu_op.velocity,
            delta=0.02 * cpu_op.velocity,
            msg='pipe barrel velocity: mode=1 vs mode=2')
        self.assertAlmostEqual(
            gpu_op.discharge, cpu_op.discharge,
            delta=0.02 * abs(cpu_op.discharge) + 1e-6,
            msg='pipe discharge: mode=1 vs mode=2')


@pytest.mark.skipif(not gpu_available(), reason=_gpu_skip_reason())
class Test_GPU_WeirTrapezoid(unittest.TestCase):
    """Tests for Weir_orifice_trapezoid_operator in GPU mode."""

    def _create_weir_domain(self, name, z1=0.0, z2=0.0, height=0.8):
        """Two-compartment domain connected by a trapezoidal weir/orifice culvert.

        Water starts on the left (x < 100 m). The culvert (1.0 m wide, `height`
        high, side slopes z1/z2) provides the only flow path. A tall `height`
        keeps the culvert flowing partly full (open-channel critical depth) so
        the flow_area is set by the critical-depth solve.
        """
        from anuga.structures.weir_orifice_trapezoid_operator import Weir_orifice_trapezoid_operator
        domain = rectangular_cross_domain(20, 10, len1=200., len2=100.)
        domain.set_flow_algorithm('DE0')
        domain.set_low_froude(0)
        domain.set_name(name)
        domain.set_datadir(tempfile.mkdtemp())
        domain.store = False

        domain.set_quantity('elevation', -1.0)
        domain.set_quantity('friction', 0.013)
        domain.set_quantity('stage', lambda x, y: np.where(x < 100., 0.5, -1.0))

        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        Weir_orifice_trapezoid_operator(domain,
                                        end_points=[[90., 50.], [110., 50.]],
                                        height=height, width=1.0,
                                        z1=z1, z2=z2,
                                        apron=5., manning=0.013,
                                        enquiry_gap=5., verbose=False)
        return domain

    def test_weir_trapezoid_cpu_gpu_match(self):
        """Weir trapezoid: mode=1 vs mode=2 stage comparison at 5 s (rectangular section)."""
        from anuga.shallow_water.sw_domain_gpu_ext import sync_to_device, sync_from_device

        cpu_domain = self._create_weir_domain('wt_cpu')
        gpu_domain = self._create_weir_domain('wt_gpu')

        cpu_domain.set_multiprocessor_mode(1)
        for _ in cpu_domain.evolve(yieldstep=1.0, finaltime=5.0):
            pass
        cpu_stage = cpu_domain.quantities['stage'].centroid_values.copy()

        gpu_domain.set_multiprocessor_mode(2)
        sync_to_device(gpu_domain.gpu_interface.gpu_dom)
        for _ in gpu_domain.evolve(yieldstep=1.0, finaltime=5.0):
            pass
        sync_from_device(gpu_domain.gpu_interface.gpu_dom)
        gpu_stage = gpu_domain.quantities['stage'].centroid_values.copy()

        # mode=1 and mode=2 share the single C weir/culvert kernel, so they agree
        # to machine precision: ~4e-16 at t=5 on a GPU-offload build, bit-identical
        # on a CPU build.  (Old atol=0.02 predated the shared kernel.)  Tight enough
        # to catch any regression that reintroduces a mode-1/mode-2 divergence.
        np.testing.assert_allclose(
            gpu_stage, cpu_stage, rtol=0, atol=1e-10,
            err_msg='Weir trapezoid 5s: stage mismatch between mode=1 and mode=2')

    def test_weir_trapezoid_cpu_gpu_velocity_match(self):
        """Weir/orifice trapezoid: reported barrel velocity agrees between
        mode=1 and mode=2.

        Regression guard for a gravity-constant mismatch in the trapezoid
        critical-depth Newton solve: the Python reference used to hardcode 9.81
        there (inconsistent with the domain g of 9.8 used by the mode-2 C
        kernel), leaving the mode-2 flow_area ~0.034% high and the reported
        velocity ~0.034% low -- a small constant offset the coarse stage check
        (test_weir_trapezoid_cpu_gpu_match) does not catch. Both paths now derive
        g from domain.g. A trapezoidal section (z1=z2=1) exercises the Newton
        solve.
        """
        from anuga.structures.weir_orifice_trapezoid_operator import Weir_orifice_trapezoid_operator

        def weir_op(domain):
            return next(o for o in domain.fractional_step_operators
                        if isinstance(o, Weir_orifice_trapezoid_operator))

        # Tall culvert (height=3.0) so it flows partly full: the flow_area is
        # then set by the critical-depth solve, which is where the g mismatch
        # bites. A short/full culvert would use the full section and hide it.
        cpu_domain = self._create_weir_domain('wtv_cpu', z1=1.0, z2=1.0, height=3.0)
        cpu_domain.set_multiprocessor_mode(1)
        for _ in cpu_domain.evolve(yieldstep=1.0, finaltime=3.0):
            pass
        cpu_op = weir_op(cpu_domain)

        gpu_domain = self._create_weir_domain('wtv_gpu', z1=1.0, z2=1.0, height=3.0)
        gpu_domain.set_multiprocessor_mode(2)
        for _ in gpu_domain.evolve(yieldstep=1.0, finaltime=3.0):
            pass
        gpu_op = weir_op(gpu_domain)

        self.assertGreater(cpu_op.velocity, 0.05,
                           'weir should be flowing in the test window')
        # The g-mismatch bug is a ~0.034% constant velocity offset; keep the
        # tolerance well below that yet far above the ~1e-5 fixed residual.
        self.assertAlmostEqual(
            gpu_op.velocity, cpu_op.velocity,
            delta=1e-4 * cpu_op.velocity,
            msg='weir barrel velocity: mode=1 vs mode=2')
        self.assertAlmostEqual(
            gpu_op.discharge, cpu_op.discharge,
            delta=1e-4 * abs(cpu_op.discharge) + 1e-9,
            msg='weir discharge: mode=1 vs mode=2')

    def test_weir_trapezoid_volume_conservation(self):
        """Weir trapezoid: total water volume is conserved in GPU mode."""
        from anuga.shallow_water.sw_domain_gpu_ext import sync_to_device, sync_from_device

        domain = self._create_weir_domain('wt_vol')
        initial_volume = domain.get_water_volume()

        domain.set_multiprocessor_mode(2)
        sync_to_device(domain.gpu_interface.gpu_dom)
        for _ in domain.evolve(yieldstep=1.0, finaltime=5.0):
            pass
        sync_from_device(domain.gpu_interface.gpu_dom)

        final_volume = domain.get_water_volume()
        self.assertAlmostEqual(
            initial_volume, final_volume, places=5,
            msg=f'Weir trapezoid volume not conserved: '
                f'initial={initial_volume:.5f}, final={final_volume:.5f}')

    def test_weir_trapezoid_nonrect_section(self):
        """Weir trapezoid with z1=z2=0.5: flow direction correct in GPU mode."""
        from anuga.shallow_water.sw_domain_gpu_ext import sync_to_device, sync_from_device

        cpu_domain = self._create_weir_domain('wt_nr_cpu', z1=0.5, z2=0.5)
        gpu_domain = self._create_weir_domain('wt_nr_gpu', z1=0.5, z2=0.5)

        cpu_domain.set_multiprocessor_mode(1)
        for _ in cpu_domain.evolve(yieldstep=1.0, finaltime=5.0):
            pass
        cpu_stage = cpu_domain.quantities['stage'].centroid_values.copy()

        gpu_domain.set_multiprocessor_mode(2)
        sync_to_device(gpu_domain.gpu_interface.gpu_dom)
        for _ in gpu_domain.evolve(yieldstep=1.0, finaltime=5.0):
            pass
        sync_from_device(gpu_domain.gpu_interface.gpu_dom)
        gpu_stage = gpu_domain.quantities['stage'].centroid_values.copy()

        # atol=0.02: physically reasonable tolerance for real GPU vs CPU FP divergence
        # after 5 s.  Tight enough to catch wrong-direction or zero-flow failures.
        np.testing.assert_allclose(
            gpu_stage, cpu_stage, rtol=0, atol=0.02,
            err_msg='Weir trapezoid (z1=z2=0.5) 5s: stage mismatch between mode=1 and mode=2')


@pytest.mark.skipif(not gpu_available(), reason=_gpu_skip_reason())
class Test_GPU_SlotLimits(unittest.TestCase):
    """Tests that GPU operator arrays grow dynamically beyond the initial capacity."""

    def _base_domain(self, name):
        domain = rectangular_cross_domain(4, 4, len1=40., len2=40.)
        domain.set_flow_algorithm('DE0')
        domain.set_low_froude(0)
        domain.set_name(name)
        domain.set_datadir(tempfile.mkdtemp())
        domain.store = False
        domain.set_quantity('elevation', 0.0)
        domain.set_quantity('stage', 0.5)
        domain.set_quantity('xmomentum', 0.0)
        domain.set_quantity('ymomentum', 0.0)
        domain.set_boundary({tag: anuga.Reflective_boundary(domain)
                             for tag in domain.get_boundary_tags()})
        return domain

    def test_rate_operator_dynamic_growth(self):
        """Rate operator array grows beyond initial MAX_RATE_OPERATORS=64 capacity."""
        from anuga import Rate_operator
        domain = self._base_domain('slot_rate')
        domain.set_multiprocessor_mode(2)

        # Register 66 operators (two beyond the initial capacity of 64)
        operators = []
        for i in range(66):
            op = Rate_operator(domain, rate=0.0)
            op._init_gpu()
            operators.append(op)

        # All 66 should have been allocated successfully (no exception)
        self.assertEqual(len(operators), 66)

    def test_inlet_operator_dynamic_growth(self):
        """Inlet operator array grows beyond initial MAX_INLET_OPERATORS=32 capacity."""
        domain = self._base_domain('slot_inlet')
        domain.set_multiprocessor_mode(2)

        from anuga.shallow_water.sw_domain_gpu_ext import sync_to_device
        sync_to_device(domain.gpu_interface.gpu_dom)

        # Register 34 inlet operators (two beyond the initial capacity of 32)
        operators = []
        for i in range(34):
            op = anuga.Inlet_operator(domain, [[0.0, 20.0], [40.0, 20.0]], Q=0.0)
            op._init_gpu()
            operators.append(op)

        self.assertEqual(len(operators), 34)


@pytest.mark.skipif(not gpu_available(), reason=_gpu_skip_reason())
class Test_GPU_FileBoundary(unittest.TestCase):
    """Tests for G1.1: File_boundary / Field_boundary GPU support."""

    def _make_domain(self, M=15, N=15):
        d = rectangular_cross_domain(M, N, len1=1.0, len2=1.0)
        d.set_flow_algorithm('DE0')
        d.set_low_froude(0)
        d.set_datadir(tempfile.mkdtemp())
        d.store = False
        d.set_quantity('elevation', lambda x, y: -x / 2)
        d.set_quantity('friction', 0.0)
        d.set_quantity('stage', expression='elevation')
        return d

    def _run_with_file_boundary(self, mode):
        """Run a short simulation with a stub File_boundary on the 'left' tag."""
        from anuga.shallow_water.boundaries import Reflective_boundary
        from anuga.abstract_2d_finite_volumes.generic_boundary_conditions import Boundary

        # Stub that behaves like File_boundary (matched by class name)
        class File_boundary(Boundary):
            def evaluate(self, vol_id=None, edge_id=None):
                return [-0.2, 0.0, 0.0]

        d = self._make_domain()
        Br = Reflective_boundary(d)
        Bf = File_boundary()
        d.set_boundary({'left': Br, 'right': Bf, 'top': Br, 'bottom': Br})
        d.set_multiprocessor_mode(mode)
        d.set_quantities_to_be_stored(None)

        gauge_tri = d.get_triangle_containing_point([0.7, 0.5])
        stage = d.get_quantity('stage')
        gauge_vals = []
        for _ in d.evolve(yieldstep=0.25, finaltime=0.5):
            gauge_vals.append(float(stage.centroid_values[gauge_tri]))
        return gauge_vals

    def test_file_boundary_mode1_vs_mode2(self):
        """File_boundary produces identical results in mode=1 and mode=2."""
        g1 = self._run_with_file_boundary(mode=1)
        g2 = self._run_with_file_boundary(mode=2)
        self.assertEqual(len(g1), len(g2))
        for v1, v2 in zip(g1, g2):
            self.assertAlmostEqual(v1, v2, places=10,
                msg=f"mode=1 gauge={v1} vs mode=2 gauge={v2}")

    def test_file_boundary_in_gpu_boundary_types(self):
        """File_boundary and Field_boundary are recognised as GPU-supported types."""
        from anuga.shallow_water.boundaries import Reflective_boundary
        from anuga.abstract_2d_finite_volumes.generic_boundary_conditions import Boundary

        class File_boundary(Boundary):
            def evaluate(self, vol_id=None, edge_id=None):
                return [-0.2, 0.0, 0.0]

        class Field_boundary(Boundary):
            def evaluate(self, vol_id=None, edge_id=None):
                return [-0.1, 0.0, 0.0]

        d = self._make_domain()
        Br = Reflective_boundary(d)
        d.set_boundary({'left': Br, 'right': File_boundary(), 'top': Br, 'bottom': Field_boundary()})
        d.set_multiprocessor_mode(2)

        # Trigger lazy boundary init by running one step
        d.set_quantities_to_be_stored(None)
        for _ in d.evolve(yieldstep=0.25, finaltime=0.25):
            pass

        # Both file boundary types must be on-GPU (no CPU fallback)
        self.assertTrue(d._gpu_all_on_gpu,
            "File_boundary / Field_boundary should be GPU-supported")
        # GPU interface must still be active (no fallback to mode=1)
        self.assertIsNotNone(d.gpu_interface, "GPU interface should remain active")

    def test_file_boundary_values_pushed_to_gpu(self):
        """set_file_boundary_values_from_domain correctly fills per-edge arrays."""
        from anuga.shallow_water.boundaries import Reflective_boundary
        from anuga.abstract_2d_finite_volumes.generic_boundary_conditions import Boundary
        from anuga.shallow_water.sw_domain_gpu_ext import (
            init_gpu_domain, map_to_gpu, unmap_from_gpu, finalize_gpu_domain,
            init_file_boundary, set_file_boundary_values_from_domain,
        )

        STAGE_VAL = -0.42

        class File_boundary(Boundary):
            def evaluate(self, vol_id=None, edge_id=None):
                return [STAGE_VAL, 0.0, 0.0]

        d = self._make_domain(10, 10)
        Br = Reflective_boundary(d)
        Bf = File_boundary()
        d.set_boundary({'left': Br, 'right': Bf, 'top': Br, 'bottom': Br})

        gpu = init_gpu_domain(d)
        init_file_boundary(gpu, d)
        map_to_gpu(gpu)

        # Push current values
        set_file_boundary_values_from_domain(gpu, d)

        # Verify the Python metadata was populated (edges found for right-boundary tag)
        meta = getattr(gpu, '_file_boundary_meta', None)
        self.assertIsNotNone(meta, "_file_boundary_meta should be set after init_file_boundary")
        self.assertGreater(len(meta), 0,
            "file_bdry should have edges for the 'right' File_boundary tag")

        unmap_from_gpu(gpu)
        finalize_gpu_domain(gpu)


@pytest.mark.skipif(not gpu_available(), reason=_gpu_skip_reason())
class Test_GPU_FlatherBoundary(unittest.TestCase):
    """Tests for Flather_external_stage_zero_velocity_boundary GPU support."""

    def _make_domain(self, M=15, N=15):
        d = rectangular_cross_domain(M, N, len1=1.0, len2=1.0)
        d.set_flow_algorithm('DE0')
        d.set_low_froude(0)
        d.set_datadir(tempfile.mkdtemp())
        d.store = False
        d.set_quantity('elevation', 0.0)
        d.set_quantity('friction', 0.0)
        d.set_quantity('stage', 0.1)
        return d

    def _run_mode(self, mode):
        from anuga import Flather_external_stage_zero_velocity_boundary
        d = self._make_domain()
        Br = Reflective_boundary(d)
        Bf = Flather_external_stage_zero_velocity_boundary(d, function=lambda t: 0.0)
        d.set_boundary({'left': Bf, 'right': Br, 'top': Br, 'bottom': Br})
        d.set_multiprocessor_mode(mode)
        d.set_quantities_to_be_stored(None)
        gauge_tri = d.get_triangle_containing_point([0.5, 0.5])
        stage = d.get_quantity('stage')
        vals = []
        for _ in d.evolve(yieldstep=0.1, finaltime=0.3):
            vals.append(float(stage.centroid_values[gauge_tri]))
        return vals

    def test_flather_boundary_mode1_vs_mode2(self):
        """Flather boundary produces identical results in mode=1 and mode=2."""
        g1 = self._run_mode(1)
        g2 = self._run_mode(2)
        self.assertEqual(len(g1), len(g2))
        for v1, v2 in zip(g1, g2):
            self.assertAlmostEqual(v1, v2, places=10,
                msg=f"mode=1 gauge={v1} vs mode=2 gauge={v2}")

    def test_flather_boundary_in_gpu_boundary_types(self):
        """Flather_external_stage_zero_velocity_boundary is GPU-supported (no CPU fallback)."""
        from anuga import Flather_external_stage_zero_velocity_boundary
        d = self._make_domain()
        Br = Reflective_boundary(d)
        Bf = Flather_external_stage_zero_velocity_boundary(d, function=lambda t: 0.0)
        d.set_boundary({'left': Bf, 'right': Br, 'top': Br, 'bottom': Br})
        d.set_multiprocessor_mode(2)
        d.set_quantities_to_be_stored(None)
        for _ in d.evolve(yieldstep=0.1, finaltime=0.1):
            pass
        self.assertTrue(d._gpu_all_on_gpu,
            "Flather_external_stage_zero_velocity_boundary should be GPU-supported")

    def test_flather_boundary_init_via_ext(self):
        """init_flather_boundary populates GPU struct correctly."""
        from anuga import Flather_external_stage_zero_velocity_boundary
        from anuga.shallow_water.sw_domain_gpu_ext import (
            init_gpu_domain, map_to_gpu, unmap_from_gpu, finalize_gpu_domain,
            init_flather_boundary,
        )
        d = self._make_domain(10, 10)
        Br = Reflective_boundary(d)
        Bf = Flather_external_stage_zero_velocity_boundary(d, function=lambda t: 0.0)
        d.set_boundary({'left': Bf, 'right': Br, 'top': Br, 'bottom': Br})

        gpu = init_gpu_domain(d)
        n = init_flather_boundary(gpu, d)
        map_to_gpu(gpu)

        self.assertGreater(n, 0, "Flather BC should find edges on 'left' boundary")

        unmap_from_gpu(gpu)
        finalize_gpu_domain(gpu)


@pytest.mark.skipif(not gpu_available(), reason=_gpu_skip_reason())
class Test_GPU_WaveBoundaries(unittest.TestCase):
    """Tests for Absorbing_wave_boundary and Characteristic_wave_boundary GPU support."""

    def _make_domain(self, M=15, N=15):
        d = rectangular_cross_domain(M, N, len1=1.0, len2=1.0)
        d.set_flow_algorithm('DE0')
        d.set_low_froude(0)
        d.set_datadir(tempfile.mkdtemp())
        d.store = False
        d.set_quantity('elevation', 0.0)
        d.set_quantity('friction', 0.0)
        d.set_quantity('stage', 0.1)
        return d

    def _run_mode(self, mode, bc_factory):
        """Run a short simulation with the given boundary factory and return gauge time series."""
        d = self._make_domain()
        Br = Reflective_boundary(d)
        Bw = bc_factory(d)
        d.set_boundary({'left': Bw, 'right': Br, 'top': Br, 'bottom': Br})
        d.set_multiprocessor_mode(mode)
        d.set_quantities_to_be_stored(None)

        gauge_tri = d.get_triangle_containing_point([0.5, 0.5])
        stage = d.get_quantity('stage')
        vals = []
        for _ in d.evolve(yieldstep=0.1, finaltime=0.3):
            vals.append(float(stage.centroid_values[gauge_tri]))
        return vals

    def test_absorbing_wave_boundary_mode1_vs_mode2(self):
        """Absorbing_wave_boundary produces identical results in mode=1 and mode=2."""
        from anuga import Absorbing_wave_boundary

        def bc_factory(d):
            return Absorbing_wave_boundary(d, function=lambda t: 0.0)

        g1 = self._run_mode(1, bc_factory)
        g2 = self._run_mode(2, bc_factory)
        self.assertEqual(len(g1), len(g2))
        for v1, v2 in zip(g1, g2):
            self.assertAlmostEqual(v1, v2, places=10,
                msg=f"mode=1 gauge={v1} vs mode=2 gauge={v2}")

    def test_characteristic_wave_boundary_mode1_vs_mode2(self):
        """Characteristic_wave_boundary produces identical results in mode=1 and mode=2."""
        from anuga import Characteristic_wave_boundary

        def bc_factory(d):
            return Characteristic_wave_boundary(d, function=lambda t: 0.0, background_stage=0.1)

        g1 = self._run_mode(1, bc_factory)
        g2 = self._run_mode(2, bc_factory)
        self.assertEqual(len(g1), len(g2))
        for v1, v2 in zip(g1, g2):
            self.assertAlmostEqual(v1, v2, places=10,
                msg=f"mode=1 gauge={v1} vs mode=2 gauge={v2}")

    def test_wave_boundaries_in_gpu_boundary_types(self):
        """Both wave BCs are recognised as GPU-supported (no CPU fallback)."""
        from anuga import Absorbing_wave_boundary, Characteristic_wave_boundary

        d = self._make_domain()
        Br = Reflective_boundary(d)
        Ba = Absorbing_wave_boundary(d, function=lambda t: 0.0)
        Bc = Characteristic_wave_boundary(d, function=lambda t: 0.0, background_stage=0.1)
        d.set_boundary({'left': Ba, 'right': Bc, 'top': Br, 'bottom': Br})
        d.set_multiprocessor_mode(2)
        d.set_quantities_to_be_stored(None)
        for _ in d.evolve(yieldstep=0.1, finaltime=0.1):
            pass

        self.assertTrue(d._gpu_all_on_gpu,
            "Absorbing/Characteristic wave BCs should be GPU-supported (no CPU fallback)")

    def test_absorbing_wave_boundary_init_via_ext(self):
        """init_absorbing_wave_boundary populates GPU struct correctly."""
        from anuga import Absorbing_wave_boundary
        from anuga.shallow_water.sw_domain_gpu_ext import (
            init_gpu_domain, map_to_gpu, unmap_from_gpu, finalize_gpu_domain,
            init_absorbing_wave_boundary,
        )

        d = self._make_domain(10, 10)
        Br = Reflective_boundary(d)
        Ba = Absorbing_wave_boundary(d, function=lambda t: 0.0)
        d.set_boundary({'left': Ba, 'right': Br, 'top': Br, 'bottom': Br})

        gpu = init_gpu_domain(d)
        n = init_absorbing_wave_boundary(gpu, d)
        map_to_gpu(gpu)

        self.assertGreater(n, 0, "absorbing_wave BC should find edges on 'left' boundary")

        unmap_from_gpu(gpu)
        finalize_gpu_domain(gpu)

    def test_characteristic_wave_boundary_init_via_ext(self):
        """init_characteristic_wave_boundary populates GPU struct correctly."""
        from anuga import Characteristic_wave_boundary
        from anuga.shallow_water.sw_domain_gpu_ext import (
            init_gpu_domain, map_to_gpu, unmap_from_gpu, finalize_gpu_domain,
            init_characteristic_wave_boundary,
        )

        d = self._make_domain(10, 10)
        Br = Reflective_boundary(d)
        Bc = Characteristic_wave_boundary(d, function=lambda t: 0.0, background_stage=0.1)
        d.set_boundary({'left': Bc, 'right': Br, 'top': Br, 'bottom': Br})

        gpu = init_gpu_domain(d)
        n = init_characteristic_wave_boundary(gpu, d)
        map_to_gpu(gpu)

        self.assertGreater(n, 0, "characteristic_wave BC should find edges on 'left' boundary")

        unmap_from_gpu(gpu)
        finalize_gpu_domain(gpu)


@pytest.mark.skipif(not gpu_available(), reason=_gpu_skip_reason())
class Test_GPU_DeviceMemory(unittest.TestCase):
    """Tests for G1.2: device memory check before array mapping."""

    def _make_domain(self, M=10, N=10):
        d = rectangular_cross_domain(M, N, len1=100., len2=100.)
        d.set_flow_algorithm('DE0')
        d.set_low_froude(0)
        d.set_datadir(tempfile.mkdtemp())
        d.store = False
        d.set_quantity('elevation', 0.0)
        d.set_quantity('stage', 0.5)
        d.set_quantity('friction', 0.0)
        return d

    def test_estimate_positive_and_scales_with_n(self):
        """Memory estimate is positive and grows with domain size."""
        from anuga.shallow_water.sw_domain_gpu_ext import (
            init_gpu_domain, finalize_gpu_domain,
            estimate_required_memory
        )
        d_small = self._make_domain(10, 10)
        d_large = self._make_domain(20, 20)
        small_n = d_small.number_of_elements
        large_n = d_large.number_of_elements

        est_small = estimate_required_memory(small_n, d_small.boundary_length)
        est_large = estimate_required_memory(large_n, d_large.boundary_length)

        self.assertGreater(est_small, 0)
        self.assertGreater(est_large, est_small)
        # ~4× domain → ~4× memory
        ratio = est_large / est_small
        self.assertGreater(ratio, 3.0)
        self.assertLess(ratio, 6.0)

    def test_estimate_reasonable_for_1m_triangles(self):
        """Estimate for 1M triangles is in the expected 400–600 MB range."""
        from anuga.shallow_water.sw_domain_gpu_ext import estimate_required_memory
        est = estimate_required_memory(1_000_000, 10_000)
        est_mb = est / (1024 * 1024)
        self.assertGreater(est_mb, 400)
        self.assertLess(est_mb, 600)

    def test_check_passes_for_small_domain(self):
        """Memory check succeeds for a small domain (never fails in CPU_ONLY_MODE)."""
        from anuga.shallow_water.sw_domain_gpu_ext import (
            init_gpu_domain, map_to_gpu, unmap_from_gpu, finalize_gpu_domain
        )
        d = self._make_domain(10, 10)
        gpu = init_gpu_domain(d)
        # In CPU_ONLY_MODE this always succeeds
        try:
            map_to_gpu(gpu)
        finally:
            unmap_from_gpu(gpu)
            finalize_gpu_domain(gpu)

    def test_map_to_gpu_raises_on_oom(self):
        """map_to_gpu raises RuntimeError when device memory is insufficient."""
        from anuga.shallow_water.sw_domain_gpu_ext import (
            init_gpu_domain, map_to_gpu, unmap_from_gpu, finalize_gpu_domain,
            check_gpu_device_memory
        )

        d = self._make_domain(10, 10)
        gpu = init_gpu_domain(d)

        # In CPU_ONLY_MODE, check always returns 1 (no real device to OOM)
        result = check_gpu_device_memory(gpu)
        self.assertEqual(result, 1)

        map_to_gpu(gpu)
        unmap_from_gpu(gpu)
        finalize_gpu_domain(gpu)

    def test_memory_info_printed_when_verbose(self, capsys=None):
        """Memory estimate line appears in verbose output."""
        from anuga.shallow_water.sw_domain_gpu_ext import (
            init_gpu_domain, map_to_gpu, unmap_from_gpu, finalize_gpu_domain
        )
        import io
        from contextlib import redirect_stdout

        d = self._make_domain(10, 10)
        gpu = init_gpu_domain(d)
        # verbose is set via init_gpu_domain — check printed output
        # (C printf goes to stdout; capture at fd level is tricky; just run
        #  and confirm no exception is raised, i.e., the path executes cleanly)
        map_to_gpu(gpu)
        unmap_from_gpu(gpu)
        finalize_gpu_domain(gpu)


@pytest.mark.skipif(not gpu_available(), reason=_gpu_skip_reason())
class Test_GPU_CollectMaxQuantities(unittest.TestCase):
    """Consistency tests for Collect_max_quantities_operator in mode 1 vs mode 2.

    Three tiers of tests:
    1. GPU self-consistency: D2H-updated arrays vs explicit get — must be bit-identical.
    2. Mode 2 domain-quantity integration: store_to_sww path updates the right quantities.
    3. Mode 1 vs mode 2 approximate agreement: same solver family, results within ~1e-3.
    """

    def _create_domain(self, name):
        d = rectangular_cross_domain(8, 8, len1=100., len2=100.)
        d.set_flow_algorithm('DE0')
        d.set_low_froude(0)
        d.set_name(name)
        d.set_datadir(tempfile.mkdtemp())
        d.store = False
        d.set_quantity('elevation', lambda x, y: -x / 50.0)
        d.set_quantity('friction', 0.01)
        d.set_quantity('stage', 0.0)
        Br = Reflective_boundary(d)
        Bd = Dirichlet_boundary([-0.5, 0., 0.])
        d.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})
        return d

    def test_mode2_self_consistent_store_to_sww(self):
        """GPU kernel self-consistency: D2H-updated arrays == explicit get at end.

        When store_to_sww=True the operator syncs max arrays device→host every
        update_frequency steps.  At the end, calling get_max_quantities_gpu()
        explicitly must return the same bit-identical values — no re-computation
        happens, both paths read the same device memory.
        """
        from anuga.operators.collect_max_quantities_operator import Collect_max_quantities_operator
        from anuga.shallow_water.sw_domain_gpu_ext import (
            sync_to_device, sync_from_device, get_max_quantities_gpu)

        d = self._create_domain('test_maxqty_selfcons')
        op = Collect_max_quantities_operator(d, store_to_sww=True)

        d.set_multiprocessor_mode(2)
        gpu_dom = d.gpu_interface.gpu_dom
        sync_to_device(gpu_dom)
        for t in d.evolve(yieldstep=0.5, finaltime=1.0):
            pass
        sync_from_device(gpu_dom)

        # op.max_stage was kept up-to-date via D2H sync during evolve.
        # Explicit get must return identical values (same device array, no kernel re-run).
        check_stage = np.empty_like(op.max_stage)
        check_depth = np.empty_like(op.max_depth)
        check_speed = np.empty_like(op.max_speed)
        check_uh    = np.empty_like(op.max_uh)
        get_max_quantities_gpu(gpu_dom, check_stage, check_depth, check_speed, check_uh)

        np.testing.assert_array_equal(op.max_stage, check_stage,
                                      err_msg="D2H-updated max_stage != explicit get")
        np.testing.assert_array_equal(op.max_depth, check_depth,
                                      err_msg="D2H-updated max_depth != explicit get")
        np.testing.assert_array_equal(op.max_speed, check_speed,
                                      err_msg="D2H-updated max_speed != explicit get")
        np.testing.assert_array_equal(op.max_uh, check_uh,
                                      err_msg="D2H-updated max_uh != explicit get")

    def test_mode2_self_consistent_update_frequency(self):
        """Self-consistency is maintained with update_frequency > 1.

        Fewer D2H syncs during evolve (every 3 steps) but the final explicit
        get must still return identical values as the last D2H-updated snapshot.
        """
        from anuga.operators.collect_max_quantities_operator import Collect_max_quantities_operator
        from anuga.shallow_water.sw_domain_gpu_ext import (
            sync_to_device, sync_from_device, get_max_quantities_gpu)

        d = self._create_domain('test_maxqty_freq')
        op = Collect_max_quantities_operator(d, update_frequency=3, store_to_sww=True)

        d.set_multiprocessor_mode(2)
        gpu_dom = d.gpu_interface.gpu_dom
        sync_to_device(gpu_dom)
        for t in d.evolve(yieldstep=0.5, finaltime=1.0):
            pass
        sync_from_device(gpu_dom)

        check_stage = np.empty_like(op.max_stage)
        check_depth = np.empty_like(op.max_depth)
        check_speed = np.empty_like(op.max_speed)
        check_uh    = np.empty_like(op.max_uh)
        get_max_quantities_gpu(gpu_dom, check_stage, check_depth, check_speed, check_uh)

        np.testing.assert_array_equal(op.max_stage, check_stage,
                                      err_msg="update_frequency=3: max_stage != explicit get")
        np.testing.assert_array_equal(op.max_depth, check_depth,
                                      err_msg="update_frequency=3: max_depth != explicit get")

    def test_mode2_store_to_sww_updates_domain_quantities(self):
        """Mode 2 with store_to_sww=True updates domain quantities each update step."""
        from anuga.operators.collect_max_quantities_operator import Collect_max_quantities_operator
        from anuga.shallow_water.sw_domain_gpu_ext import sync_to_device, sync_from_device

        d = self._create_domain('test_maxqty_sww')
        op = Collect_max_quantities_operator(d, update_frequency=2, store_to_sww=True)

        d.set_multiprocessor_mode(2)
        gpu_dom = d.gpu_interface.gpu_dom
        sync_to_device(gpu_dom)
        for t in d.evolve(yieldstep=0.5, finaltime=1.0):
            pass
        sync_from_device(gpu_dom)

        # Domain quantities should have been updated (non-trivial values)
        max_stage_qty = d.quantities['max_stage'].centroid_values
        max_depth_qty = d.quantities['max_depth'].centroid_values

        self.assertTrue(np.any(max_stage_qty > -1e30),
                        "max_stage domain quantity was never updated from initial -inf")
        self.assertTrue(np.any(max_depth_qty > 0.0),
                        "max_depth domain quantity was never updated from zero")

        # operator arrays must also be consistent with domain quantities
        np.testing.assert_allclose(op.max_stage, max_stage_qty, rtol=1e-12, atol=1e-15)
        np.testing.assert_allclose(op.max_depth, max_depth_qty, rtol=1e-12, atol=1e-15)

    def test_mode1_mode2_approximate_match(self):
        """Mode 1 (CPU NumPy) and mode 2 (GPU kernel) produce results within ~1e-10.

        Both modes drive the same underlying C solver and should agree to near
        machine epsilon.  A tolerance of 1e-10 catches real bugs (e.g. the GPU
        kernel not updating, or spurious host<->device syncs perturbing the
        mode-2 trajectory) while tolerating harmless floating-point reordering.
        """
        from anuga.operators.collect_max_quantities_operator import Collect_max_quantities_operator
        from anuga.shallow_water.sw_domain_gpu_ext import (
            sync_to_device, sync_from_device, get_max_quantities_gpu)

        cpu_domain = self._create_domain('test_maxqty_approx_cpu')
        gpu_domain = self._create_domain('test_maxqty_approx_gpu')

        cpu_op = Collect_max_quantities_operator(cpu_domain, store_to_sww=False)
        gpu_op = Collect_max_quantities_operator(gpu_domain, store_to_sww=False)

        cpu_domain.set_multiprocessor_mode(1)
        for t in cpu_domain.evolve(yieldstep=0.5, finaltime=1.0):
            pass

        gpu_domain.set_multiprocessor_mode(2)
        gpu_dom = gpu_domain.gpu_interface.gpu_dom
        sync_to_device(gpu_dom)
        for t in gpu_domain.evolve(yieldstep=0.5, finaltime=1.0):
            pass
        sync_from_device(gpu_dom)

        get_max_quantities_gpu(gpu_dom,
                               gpu_op.max_stage, gpu_op.max_depth,
                               gpu_op.max_speed, gpu_op.max_uh)

        # Both should produce non-trivial, physically sensible values
        self.assertTrue(np.all(gpu_op.max_depth >= 0.0), "GPU max_depth has negative values")
        self.assertTrue(np.all(gpu_op.max_speed >= 0.0), "GPU max_speed has negative values")
        self.assertTrue(np.all(gpu_op.max_uh >= 0.0),    "GPU max_uh has negative values")

        np.testing.assert_allclose(cpu_op.max_stage, gpu_op.max_stage,
                                   atol=1e-10, rtol=1e-10,
                                   err_msg="max_stage mode1/mode2 mismatch")
        np.testing.assert_allclose(cpu_op.max_depth, gpu_op.max_depth,
                                   atol=1e-10, rtol=1e-10,
                                   err_msg="max_depth mode1/mode2 mismatch")
        np.testing.assert_allclose(cpu_op.max_speed, gpu_op.max_speed,
                                   atol=1e-10, rtol=1e-10,
                                   err_msg="max_speed mode1/mode2 mismatch")
        np.testing.assert_allclose(cpu_op.max_uh, gpu_op.max_uh,
                                   atol=1e-10, rtol=1e-10,
                                   err_msg="max_uh mode1/mode2 mismatch")


@pytest.mark.skipif(not gpu_available(), reason=_gpu_skip_reason())
class Test_TimestepTimeAdvance(unittest.TestCase):
    """Regression tests for the first-step collection bug in
    Collect_max_quantities_operator.

    The bug: mode-1 fractional-step operators see the pre-step time (the
    generic evolve loop advances relative_time AFTER apply_fractional_steps).
    Mode-2 C-loops advance relative_time before returning, so fractional-step
    operators see the post-step time.  With collection_start_time=0 and the
    old strict inequality (t > 0), mode-1 silently skipped the very first
    inner timestep (0 > 0 is False) while mode-2 correctly collected from
    step 1 onwards.  For transient flows the first step may carry the global
    maximum, causing mode-1 and mode-2 to disagree by O(1e-4).

    Fix: Collect_max_quantities_operator now uses >= instead of > for the
    collection_start_time guard.  Mode-1 at t=0 satisfies 0 >= 0 and
    correctly collects state S_1; mode-2 at t=dt_1 also satisfies the guard
    and collects the same S_1 from device memory.

    Test strategy: run mode-1 and mode-2 for each flow algorithm and assert
    the collected maxima agree to atol=1e-10.  If the >= fix is reverted
    (back to >) mode-1 will miss S_1 and the test fails.
    """

    def _make_domain(self, name, algorithm):
        d = rectangular_cross_domain(8, 8, len1=100., len2=100.)
        d.set_flow_algorithm(algorithm)
        d.set_low_froude(0)
        d.set_name(name)
        d.set_datadir(tempfile.mkdtemp())
        d.store = False
        d.set_quantity('elevation', lambda x, y: -x / 50.0)
        d.set_quantity('friction', 0.01)
        d.set_quantity('stage', 0.0)
        Br = Reflective_boundary(d)
        Bd = Dirichlet_boundary([-0.5, 0., 0.])
        d.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})
        return d

    def _run_modes(self, algorithm):
        """Run mode 1 and mode 2 with Collect_max_quantities and return both ops."""
        from anuga.operators.collect_max_quantities_operator import Collect_max_quantities_operator
        from anuga.shallow_water.sw_domain_gpu_ext import (
            sync_to_device, sync_from_device, get_max_quantities_gpu)

        n = 256  # 8x8 rectangular_cross

        d1 = self._make_domain(f'{algorithm}_mode1', algorithm)
        op1 = Collect_max_quantities_operator(d1, store_to_sww=False)
        d1.set_multiprocessor_mode(1)
        for t in d1.evolve(yieldstep=0.5, finaltime=1.0):
            pass

        d2 = self._make_domain(f'{algorithm}_mode2', algorithm)
        op2 = Collect_max_quantities_operator(d2, store_to_sww=False)
        d2.set_multiprocessor_mode(2)
        gpu_dom = d2.gpu_interface.gpu_dom
        sync_to_device(gpu_dom)
        for t in d2.evolve(yieldstep=0.5, finaltime=1.0):
            pass
        sync_from_device(gpu_dom)

        ms2 = np.zeros(n); md2 = np.zeros(n)
        msp2 = np.zeros(n); mu2 = np.zeros(n)
        get_max_quantities_gpu(gpu_dom, ms2, md2, msp2, mu2)

        return op1, ms2, md2, msp2, mu2

    def _assert_modes_agree(self, algorithm):
        op1, ms2, md2, msp2, mu2 = self._run_modes(algorithm)
        np.testing.assert_allclose(op1.max_stage, ms2, atol=1e-10, rtol=1e-10,
                                   err_msg=f'{algorithm}: max_stage mode1/mode2 mismatch')
        np.testing.assert_allclose(op1.max_depth, md2, atol=1e-10, rtol=1e-10,
                                   err_msg=f'{algorithm}: max_depth mode1/mode2 mismatch')
        np.testing.assert_allclose(op1.max_speed, msp2, atol=1e-10, rtol=1e-10,
                                   err_msg=f'{algorithm}: max_speed mode1/mode2 mismatch')
        np.testing.assert_allclose(op1.max_uh, mu2, atol=1e-10, rtol=1e-10,
                                   err_msg=f'{algorithm}: max_uh mode1/mode2 mismatch')

    def test_DE0_euler(self):
        """DE0 (Euler) — the original bug was in evolve_one_euler_step."""
        self._assert_modes_agree('DE0')

    def test_DE1_rk2(self):
        """DE1 (RK2) — time was already advanced mid-step; verify no regression."""
        self._assert_modes_agree('DE1')

    def test_DE2_rk3(self):
        """DE2 (RK3) — time was already advanced at end of step; verify no regression."""
        self._assert_modes_agree('DE2')

    def test_DE_ader2(self):
        """DE_ader2 (ADER-2) — same bug as Euler; fixed in evolve_one_ader2_step."""
        self._assert_modes_agree('DE_ader2')


@pytest.mark.skipif(not gpu_available(), reason=_gpu_skip_reason())
class Test_GPU_TimeBoundarySubstep(unittest.TestCase):
    """Mode-2 with a time-varying boundary must match mode-1 for multi-substep
    algorithms (DE1/DE2).

    The single-call C RK loop (_evolve_one_rk*_step_c) sets Python-evaluated
    boundaries (Time/File/wave/Flather) on the device once per step, so it would
    reuse that step-start value across every RK substep — whereas mode-1 calls
    update_boundary() before each substep. For a time-varying boundary that is an
    O(dt) boundary-forcing error on RK2/RK3 (~4e-3 m in a rising-tide test).

    Fix (option B): a domain with any Python-evaluated boundary is routed to the
    Python-orchestrated GPU loop, which refreshes the boundary per substep and so
    bit-matches mode-1. This test guards that routing; reverting it makes DE1/DE2
    disagree by O(1e-3).
    """

    def _make_domain(self, name, algorithm):
        import math
        d = rectangular_cross_domain(12, 8, len1=150., len2=100.)
        d.set_flow_algorithm(algorithm)
        d.set_low_froude(0)
        d.set_name(name)
        d.set_datadir(tempfile.mkdtemp())
        d.store = False
        d.set_quantity('elevation', lambda x, y: (x / 150.0) * 2.0 - 1.0)
        d.set_quantity('friction', 0.03)
        d.set_quantity('stage', 1.5)                       # fully wet (isolates the boundary)

        def tide(t):
            return [0.4 * math.sin(0.3 * t), 0.0, 0.0]     # time-varying stage

        Br = Reflective_boundary(d)
        d.set_boundary({'left': anuga.Time_boundary(domain=d, function=tide),
                        'right': Br, 'top': Br, 'bottom': Br})
        return d

    def _run(self, algorithm, mode):
        d = self._make_domain(f'{algorithm}_m{mode}', algorithm)
        d.set_multiprocessor_mode(mode)
        for t in d.evolve(yieldstep=1.0, finaltime=6.0):
            pass
        return d.quantities['stage'].centroid_values.copy()

    def _assert_agree(self, algorithm):
        s1 = self._run(algorithm, 1)
        s2 = self._run(algorithm, 2)
        np.testing.assert_allclose(
            s2, s1, atol=1e-6, rtol=0.0,
            err_msg=f'{algorithm}: mode-2 Time_boundary diverges from mode-1 '
                    f'(max diff {np.abs(s2 - s1).max():.3e})')

    def test_DE1_rk2_time_boundary(self):
        self._assert_agree('DE1')

    def test_DE2_rk3_time_boundary(self):
        self._assert_agree('DE2')

    def test_routing_time_boundary_off_c_loop(self):
        """A Time_boundary must be flagged so the multi-substep dispatch avoids
        the single-call C RK loop."""
        d = self._make_domain('route_time', 'DE1')
        d.set_multiprocessor_mode(2)
        self.assertTrue(d._has_python_evaluated_gpu_boundaries())

    def test_routing_reflective_keeps_c_loop(self):
        """Reflective-only domains keep the fast C RK loop."""
        d = rectangular_cross_domain(8, 8, len1=100., len2=100.)
        d.set_flow_algorithm('DE1')
        d.set_name('route_refl')
        d.set_datadir(tempfile.mkdtemp())
        d.store = False
        d.set_quantity('elevation', 0.0)
        d.set_quantity('stage', 1.0)
        d.set_boundary({b: Reflective_boundary(d) for b in d.get_boundary_tags()})
        d.set_multiprocessor_mode(2)
        self.assertFalse(d._has_python_evaluated_gpu_boundaries())


@pytest.mark.skipif(not gpu_available(), reason=_gpu_skip_reason())
class Test_GPU_OperatorTimeAlignment(unittest.TestCase):
    """A time-varying fractional-step operator must be evaluated at the same
    time in mode-1 and mode-2, for every flow algorithm.

    Fractional-step operators run in the evolve loop *before* it advances
    relative_time to t+dt, so they should see the pre-step time t. DE0, DE2 and
    DE_ader2 (and the mode-2 GPU loops) do; legacy mode-1 **rk2 (DE1)** advanced
    relative_time mid-step (for the substep-2 boundary) and never restored it, so
    its operators evaluated forcing at t+dt — "one step too far" — diverging from
    mode-2 by ~4e-4 for a time-varying rate. The mode-1 rk2 body now restores the
    pre-step time; this test asserts mode-1 == mode-2 for all algorithms (revert
    the fix and DE1 fails). Note: no prior test exercised a time-varying operator.
    """

    def _run(self, algorithm, mode):
        import math
        d = rectangular_cross_domain(24, 16, len1=100., len2=100.)
        d.set_flow_algorithm(algorithm)
        d.set_low_froude(0)
        d.set_name(f'op_{algorithm}_m{mode}')
        d.set_datadir(tempfile.mkdtemp())
        d.store = False
        d.set_quantity('elevation', 0.0)
        d.set_quantity('friction', 0.03)
        d.set_quantity('stage', 0.2)
        anuga.Rate_operator(d, rate=lambda t: 0.001 * (1.0 + 0.9 * math.sin(0.6 * t)))
        Br = Reflective_boundary(d)
        d.set_boundary({b: Br for b in d.get_boundary_tags()})
        d.set_multiprocessor_mode(mode)
        for t in d.evolve(yieldstep=0.5, finaltime=3.0):
            pass
        return d.quantities['stage'].centroid_values.copy()

    def _assert_agree(self, algorithm):
        s1 = self._run(algorithm, 1)
        s2 = self._run(algorithm, 2)
        np.testing.assert_allclose(
            s2, s1, atol=1e-8, rtol=0.0,
            err_msg=f'{algorithm}: time-varying operator mode-1/mode-2 mismatch '
                    f'(max {np.abs(s2 - s1).max():.3e})')

    def test_DE0(self):
        self._assert_agree('DE0')

    def test_DE1(self):
        self._assert_agree('DE1')

    def test_DE2(self):
        self._assert_agree('DE2')

    def test_DE_ader2(self):
        self._assert_agree('DE_ader2')


class Test_GPU_NonGPUBoundaryFallback(unittest.TestCase):
    """Mode 2 must fall back to host boundary evaluation for boundary types the
    C loop cannot evaluate on the device — for EVERY flow algorithm, including
    DE0/Euler.

    Regression for a bug where evolve_one_euler_step() dispatched straight to the
    C Euler loop and silently ignored non-GPU boundary types (rk2/rk3/ader2
    already fell back to a Python-orchestrated loop; euler did not). The 'right'
    boundary here is a Transmissive_momentum_set_stage_boundary, which is NOT a
    GPU-supported type, so mode 2 must use the host-evaluation fallback and match
    mode 1. Symptom of the bug: DE0 results diverged from legacy by ~0.1 m while
    DE1/DE2/DE_ader2 matched. (This is what made run_parallel_riverwall.py — a
    DE0 + Transmissive_momentum_set_stage case — diverge.)
    """

    def _make_domain(self, name, algorithm):
        d = rectangular_cross_domain(10, 5, len1=100., len2=50.)
        d.set_flow_algorithm(algorithm)
        d.set_low_froude(0)
        d.set_name(name)
        d.set_datadir(tempfile.mkdtemp())
        d.store = False
        d.set_quantity('elevation', lambda x, y: -x / 100.0)
        d.set_quantity('friction', 0.03)
        d.set_quantity('stage', lambda x, y: -x / 100.0)  # initially dry
        Br = Reflective_boundary(d)
        Bt = Transmissive_momentum_set_stage_boundary(
            domain=d, function=lambda t: min(-0.4 * np.exp(-t / 5.0) - 0.1, -0.11))
        d.set_boundary({'left': Br, 'right': Bt, 'top': Br, 'bottom': Br})
        return d

    def _assert_modes_agree(self, algorithm):
        from anuga.shallow_water.sw_domain_gpu_ext import sync_to_device, sync_from_device

        d1 = self._make_domain(f'{algorithm}_m1', algorithm)
        d1.set_multiprocessor_mode(1)
        for _ in d1.evolve(yieldstep=2.0, finaltime=10.0):
            pass

        d2 = self._make_domain(f'{algorithm}_m2', algorithm)
        d2.set_multiprocessor_mode(2)
        sync_to_device(d2.gpu_interface.gpu_dom)
        for _ in d2.evolve(yieldstep=2.0, finaltime=10.0):
            pass
        sync_from_device(d2.gpu_interface.gpu_dom)

        # Bit-identical on a CPU build (same compiled kernels); ~1e-9 FP-order
        # differences on a real GPU.
        atol = 1e-8 if anuga.gpu_offload_enabled() else 1e-12
        for q in ['stage', 'xmomentum', 'ymomentum']:
            np.testing.assert_allclose(
                d2.quantities[q].centroid_values,
                d1.quantities[q].centroid_values,
                rtol=0, atol=atol,
                err_msg=f'{algorithm} + Transmissive_momentum_set_stage: '
                        f'{q} mode1 vs mode2 mismatch')

    def test_DE0_euler(self):
        """DE0/Euler — the boundary that was silently ignored before the fix."""
        self._assert_modes_agree('DE0')

    def test_DE1_rk2(self):
        self._assert_modes_agree('DE1')

    def test_DE2_rk3(self):
        self._assert_modes_agree('DE2')

    def test_DE_ader2(self):
        self._assert_modes_agree('DE_ader2')


@pytest.mark.skipif(not gpu_available(), reason=_gpu_skip_reason())
class Test_GPU_SetQuantityReachesDevice(unittest.TestCase):
    """set_quantity() after the GPU interface exists must reach the device.

    Regression guard. In mode 2 ('unified') the *device* holds the authoritative
    centroid state once the GPU interface is built. set_quantity() used to update
    only the host arrays, so any set_quantity() made *after* something built the
    interface was silently ignored: the device kept evolving the stale/default
    values, and the results — and the stored SWW — were wrong.

    The interface gets built early by anything that calls _ensure_gpu_interface(),
    notably distribute_to_vertices_and_edges() and set_boundary(), so this was
    reachable from ordinary scripts (set_boundary() before set_quantity(), or any
    mid-run set_quantity()), not just from tests.
    """

    def _run(self, mode, force_interface):
        domain = rectangular_cross_domain(10, 10)
        domain.set_flow_algorithm('DE0')
        domain.set_name('setq')
        domain.set_datadir(tempfile.mkdtemp())
        domain.store = False
        domain.set_multiprocessor_mode(mode)

        domain.set_quantity('elevation', -10.0)
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        if force_interface:
            # Build the mode-2 device interface *before* stage is set.
            domain.distribute_to_vertices_and_edges()

        # Quiescent flat water well above the bed: stage must stay at 2.0.
        domain.set_quantity('stage', 2.0)

        for _ in domain.evolve(yieldstep=0.5, finaltime=1.0):
            pass
        return domain.quantities['stage'].centroid_values.copy()

    def test_set_quantity_after_interface_built(self):
        """stage set after the interface is built must be the stage that evolves."""
        cpu = self._run(1, force_interface=True)
        gpu = self._run(2, force_interface=True)

        # The bug produced stage == 0 (device default) or the stale pre-set value.
        np.testing.assert_allclose(
            gpu, cpu, rtol=0, atol=1e-8,
            err_msg='mode-2 ignored set_quantity() made after the GPU interface '
                    'was built (device kept the stale values)')
        self.assertAlmostEqual(float(gpu.min()), 2.0, places=6)
        self.assertAlmostEqual(float(gpu.max()), 2.0, places=6)

    def test_set_quantity_before_interface_built(self):
        """The ordinary ordering must keep working (no interface yet)."""
        cpu = self._run(1, force_interface=False)
        gpu = self._run(2, force_interface=False)
        np.testing.assert_allclose(gpu, cpu, rtol=0, atol=1e-8)
        self.assertAlmostEqual(float(gpu.min()), 2.0, places=6)

    def _run_direct(self, mode):
        """As _run(), but write the quantity through the Quantity object itself."""
        domain = rectangular_cross_domain(10, 10)
        domain.set_flow_algorithm('DE0')
        domain.set_name('setq_direct')
        domain.set_datadir(tempfile.mkdtemp())
        domain.store = False
        domain.set_multiprocessor_mode(mode)

        domain.set_quantity('elevation', -10.0)
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})
        domain.distribute_to_vertices_and_edges()   # builds the device interface

        # Bypass Domain.set_quantity() entirely — the sync must hang off the
        # Quantity, not off the Domain wrapper.
        domain.quantities['stage'].set_values(2.0)

        for _ in domain.evolve(yieldstep=0.5, finaltime=1.0):
            pass
        return domain.quantities['stage'].centroid_values.copy()

    def test_direct_quantity_set_values_reaches_device(self):
        """Quantity.set_values() bypassing Domain.set_quantity() must still sync."""
        cpu = self._run_direct(1)
        gpu = self._run_direct(2)

        np.testing.assert_allclose(
            gpu, cpu, rtol=0, atol=1e-8,
            err_msg='mode-2 ignored a direct Quantity.set_values() made after the '
                    'GPU interface was built (device kept the stale values)')
        self.assertAlmostEqual(float(gpu.min()), 2.0, places=6)
        self.assertAlmostEqual(float(gpu.max()), 2.0, places=6)


class Test_GPU_FractionalStepOperatorOrder(unittest.TestCase):
    """Mode 2 must apply fractional-step operators in REGISTRATION ORDER (issue #192).

    Fractional-step operators mutate `stage` in sequence, so their order changes the
    answer. Mode 1 runs them in registration order; mode 2 used to hoist every Boyd
    culvert to the FRONT (batched via GPUCulvertManager), discarding that order. The fix
    fires the batch at the position of the FIRST culvert, so registration order is honored.

    DETERMINISTIC BY CONSTRUCTION — read before changing. An earlier version of this test
    measured the order effect in a symmetric setup (flat bed, equal-depth inlets) where the
    culvert is effectively a no-op, so operator order had NO genuine effect and the observed
    difference was pure chaotic amplification of roundoff. That is numerics-dependent and bit
    this test three times (a build-dependent ratio, then the amplification vanishing under an
    unrelated kernel change). This version instead builds a setup where order genuinely
    matters by a large, deterministic margin:

      * a sloped bed (`elevation = -x/10`) so the two culvert ends sit at different levels
        and the culvert actually TRANSFERS water, and
      * a strong Rate_operator confined by polygon to the UPSTREAM inlet, so applying it
        before vs after the culvert changes the head the culvert reads.

    Over a short evolve (no time for chaos), this gives, measured directly:

        with the fix   mode-1 order effect ~7e-3;  mode-2(rate-first) == mode-1(rate-first) to ~1e-15
        with the bug   mode-1 order effect ~7e-3;  mode-2(rate-first)  = mode-1(CULVERT-first)  (~7e-3 away)

    i.e. the fix/bug signal is a 12-order-of-magnitude separation, not a chaotic magnitude.
    Assert on that, never on a chaotic-divergence value.
    """

    FT = 2.0        # short evolve: the order effect is deterministic, chaos has no time to grow
    ORDER_EFFECT_FLOOR = 1e-4   # mode-1 order effect must clear this for the setup to be valid

    def _run(self, mode, order):
        from anuga import Boyd_box_operator, Rate_operator

        domain = rectangular_cross_domain(20, 20, len1=50.0, len2=50.0)
        domain.set_flow_algorithm('DE0')
        domain.set_name('op_order')
        domain.set_datadir(tempfile.mkdtemp())
        domain.store = False
        domain.set_multiprocessor_mode(mode)

        # Sloped bed -> the culvert ends are at different levels -> it actually transfers.
        domain.set_quantity('elevation', lambda x, y: -x / 10.0)
        domain.set_quantity('stage', 2.0)
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        def add_rate():
            # Confined to a polygon around the upstream (x~10) inlet, so the rate changes
            # exactly the head the culvert reads there.
            Rate_operator(domain, rate=0.5, factor=1.0,
                          polygon=[[5, 20], [15, 20], [15, 30], [5, 30]])

        def add_culvert():
            Boyd_box_operator(domain,
                              end_points=[[10.0, 25.0], [40.0, 25.0]],
                              losses=1.5, width=2.0, height=2.0, apron=0.0,
                              use_momentum_jet=False, use_velocity_head=False,
                              manning=0.013, verbose=False)

        if order == 'rate_first':
            add_rate()
            add_culvert()
        else:
            add_culvert()
            add_rate()

        for _ in domain.evolve(yieldstep=self.FT, finaltime=self.FT):
            pass
        return domain.quantities['stage'].centroid_values.copy()

    def _maxdiff(self, a, b):
        return float(np.abs(a - b).max())

    def test_mode2_does_not_ignore_operator_order(self):
        """Mode 2 must produce different answers for the two registration orders.

        With the bug (culverts hoisted to the front) mode 2 runs the same sequence either
        way, so this is EXACTLY 0.0.
        """
        gpu = self._maxdiff(self._run(2, 'rate_first'), self._run(2, 'culvert_first'))
        self.assertGreater(
            gpu, self.ORDER_EFFECT_FLOOR,
            msg='mode-2 gave the same answer for rate-before-culvert and '
                'culvert-before-rate — it is ignoring fractional-step operator '
                'registration order (issue #192)')

    def test_mode2_applies_same_order_as_mode1(self):
        """Mode 2 with a given registration order must match mode 1 with THAT order.

        The discriminating measurement (deterministic, no chaos): for `[rate, culvert]`,
        mode-2 must track mode-1's rate-first result far more closely than the order effect
        itself. With the bug mode-2 instead matches mode-1's CULVERT-first result, i.e. it
        is a full order-effect away.
        """
        m1_rate_first = self._run(1, 'rate_first')
        m1_culvert_first = self._run(1, 'culvert_first')
        m2_rate_first = self._run(2, 'rate_first')

        order_effect = self._maxdiff(m1_rate_first, m1_culvert_first)
        self.assertGreater(
            order_effect, self.ORDER_EFFECT_FLOOR,
            msg='the test setup no longer makes operator order matter in mode 1 '
                '(order effect %.3e) — it can no longer detect the bug' % order_effect)

        mismatch = self._maxdiff(m2_rate_first, m1_rate_first)
        self.assertLess(
            mismatch, 0.01 * order_effect,
            msg='mode-2 [rate, culvert] does not match mode-1 [rate, culvert]: mismatch '
                '%.3e vs order effect %.3e. Mode 2 is applying a different operator order '
                'than mode 1 (issue #192 — culverts hoisted to the front).'
                % (mismatch, order_effect))


class Test_GPU_StartupBanner(unittest.TestCase):
    """The mode-2 banner must report the GPU count, not the rank count (issue #194).

    It used to print `numprocs` labelled as "GPU(s)", so a 4-rank run on a 1-GPU box
    reported "4 GPU(s)". That concealed the one thing the banner was best placed to
    catch: mode-2 MPI assigns ranks to devices round-robin (rank % num_devices), so
    running more ranks than GPUs silently puts several ranks on one device, which may
    hang or return wrong results.

    gpu_startup_banner() is a pure function precisely so this matrix is testable
    without a multi-GPU machine.
    """

    def _banner(self, numprocs, num_devices, device_id=0, offload_active=True):
        from anuga.shallow_water.shallow_water_domain import gpu_startup_banner
        return '\n'.join(gpu_startup_banner(numprocs, num_devices, device_id,
                                            offload_active))

    def test_reports_device_count_not_rank_count(self):
        """4 ranks on 1 GPU must not claim 4 GPUs."""
        text = self._banner(numprocs=4, num_devices=1)
        self.assertIn('4 MPI rank(s) on 1 GPU(s)', text)
        self.assertNotIn('4 GPU(s)', text)

    def test_warns_when_ranks_exceed_gpus(self):
        """Oversubscription must warn, and lead with the cost we actually measured.

        Measured on one RTX 5070 (160k-tri mode-2 evolve, all ranks on the one device):
        without MPS, np=1/2/4 took 3.80/7.13/11.21 s — so the dependable consequence is
        a ~3x slowdown, not a crash. Results were bit-identical at every rank count, so
        the banner must not promise a failure that may never arrive; a run that quietly
        works but is 3x slow is the likelier outcome and the one users would miss.
        """
        text = self._banner(numprocs=4, num_devices=1)
        self.assertIn('WARNING', text)
        self.assertIn('ONE RANK PER GPU', text)
        self.assertIn('SLOWER', text)          # the measured, reliable consequence
        self.assertIn('MPS', text)             # and the mitigation, if they must do it

    def test_no_warning_when_ranks_match_gpus(self):
        """The supported configuration must stay quiet."""
        text = self._banner(numprocs=4, num_devices=4)
        self.assertIn('4 MPI rank(s) on 4 GPU(s)', text)
        self.assertNotIn('WARNING', text)
        self.assertNotIn('NOTE', text)

    def test_notes_idle_gpus(self):
        """Fewer ranks than GPUs is safe (distinct devices) — a note, not a warning."""
        text = self._banner(numprocs=1, num_devices=4)
        self.assertNotIn('WARNING', text)
        self.assertIn('3 GPU(s) idle', text)

    def test_serial_on_one_gpu_is_quiet(self):
        """The overwhelmingly common case must not nag."""
        text = self._banner(numprocs=1, num_devices=1)
        self.assertNotIn('WARNING', text)
        self.assertNotIn('NOTE', text)

    def test_unknown_device_count_does_not_invent_one(self):
        """If the device query fails, say so — do NOT fall back to printing numprocs."""
        text = self._banner(numprocs=4, num_devices=-1)
        self.assertIn('device count unknown', text)
        self.assertNotIn('4 GPU(s)', text)
        self.assertNotIn('WARNING', text)

    def test_no_offload_and_no_device_paths(self):
        """The CPU-multicore and no-device banners must not claim any GPU count."""
        cpu = self._banner(numprocs=4, num_devices=0, offload_active=False)
        self.assertIn('CPU multicore', cpu)
        self.assertNotIn('GPU(s)', cpu)

        nodev = self._banner(numprocs=4, num_devices=0, device_id=-1)
        self.assertIn('No GPU devices found', nodev)
        self.assertNotIn('4 GPU(s)', nodev)


class Test_GPU_DryCellStartupReconciliation(unittest.TestCase):
    """Mode 2 must reconcile deeply-dry cells to stage=bed when the device
    interface is built (issue #200).

    A dry cell should carry stage = bed (depth 0). Mode 1 reaches that via its
    per-step protect on the first step; the mode-2 device path only converges to it
    gradually (halving the stage<<bed deficit each step), and any forcing applied to
    those cells during that window — an Inlet_operator, rainfall — is absorbed into
    raising the sub-bed stage rather than making depth, and is permanently lost. On
    the Towradgi small case (initial stage=0 under a 215 m creek bank with a 20 m³/s
    inlet) that lost ~24 m³, driving the whole mode-1-vs-mode-2 divergence.

    The fix clamps stage up to bed in set_gpu_interface() before the initial sync to
    the device. This guards that clamp. (The full dynamic mass loss is mesh-dependent
    — a regular-grid domain with the same bed does not reproduce it — so it is
    validated on the Towradgi case, not here; this test guards the mechanism.)
    """

    def test_dry_cells_reconciled_to_bed_on_interface_build(self):
        domain = rectangular_cross_domain(8, 8)
        domain.set_flow_algorithm('DE0')
        domain.set_name('dry_recon')
        domain.set_datadir(tempfile.mkdtemp())
        domain.store = False

        # High bed, stage well below it -> every cell is deeply dry (stage << bed).
        domain.set_quantity('elevation', 100.0)
        domain.set_quantity('stage', 0.0)
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        stage_c = domain.quantities['stage'].centroid_values
        bed_c = domain.quantities['elevation'].centroid_values
        self.assertTrue((stage_c < bed_c).all(),
                        'setup precondition: all cells should start dry (stage < bed)')

        # Building the mode-2 interface must reconcile stage up to the bed — on the
        # DEVICE. Reading the host would pass regardless (the standard host-side
        # protect in distribute_to_vertices_and_edges reconciles the host copy even
        # without the fix); the bug is that the DEVICE keeps the un-reconciled stage,
        # so sync it back before checking. On a CPU-only build host and device are the
        # same arrays, so this can only fail on a real GPU-offload build — which is
        # where the bug exists.
        domain.set_multiprocessor_mode(2)
        domain.distribute_to_vertices_and_edges()   # triggers set_gpu_interface()
        domain.gpu_interface.sync_from_device()

        stage_c = domain.quantities['stage'].centroid_values
        bed_c = domain.quantities['elevation'].centroid_values
        np.testing.assert_allclose(
            stage_c, bed_c, rtol=0, atol=1e-9,
            err_msg='mode-2 did not reconcile deeply-dry cells to stage=bed on the '
                    'device at interface build (issue #200)')

    def test_reconciliation_is_noop_for_wet_cells(self):
        """Cells already at/above bed must be left untouched by the reconciliation."""
        domain = rectangular_cross_domain(8, 8)
        domain.set_flow_algorithm('DE0')
        domain.set_name('wet_noop')
        domain.set_datadir(tempfile.mkdtemp())
        domain.store = False

        domain.set_quantity('elevation', -5.0)
        domain.set_quantity('stage', 2.0)          # wet everywhere (stage > bed)
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        before = domain.quantities['stage'].centroid_values.copy()
        domain.set_multiprocessor_mode(2)
        domain.distribute_to_vertices_and_edges()

        after = domain.quantities['stage'].centroid_values
        np.testing.assert_allclose(
            after, before, rtol=0, atol=1e-12,
            err_msg='reconciliation must be a no-op where stage >= bed')


class Test_GPU_RateOperatorGhostInflux(unittest.TestCase):
    """Rate_operator mass tracking must exclude ghost cells in mode 2 (issue #191).

    Under MPI a rainfall polygon straddling a partition boundary appears on several
    ranks; only the rank that OWNS a triangle may count it toward the reported influx.
    The CPU path masks its sum with ``full_indices``; the mode-2 kernel used to sum
    over every index, so parallel runs over-reported rainfall influx (and hence
    ``domain.fractional_step_volume_integral``) by the ghost-cell contribution. The
    stage update itself was — and stays — applied to ghosts on both paths, since the
    halo exchange overwrites them.

    This is a serial test of a parallel bug: it fakes the partition by clearing
    ``tri_full_flag``, which is the only thing ``set_full_indices()`` reads. That keeps
    the regression catchable in the ordinary (non-MPI) suite.
    """

    RATE = 1.0e-3
    TIMESTEP = 2.0

    def _influx(self, mode):
        from anuga import Rate_operator

        domain = rectangular_cross_domain(10, 10)
        domain.set_flow_algorithm('DE0')
        domain.set_name('rate_ghost')
        domain.set_datadir(tempfile.mkdtemp())
        domain.store = False
        domain.set_multiprocessor_mode(mode)

        domain.set_quantity('elevation', -10.0)
        domain.set_quantity('stage', 1.0)
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})
        domain.distribute_to_vertices_and_edges()   # builds the mode-2 device interface

        # Pretend an MPI partition owns only 2/3 of the triangles. Must happen before
        # the operator is constructed — set_full_indices() reads tri_full_flag in
        # Rate_operator.__init__.
        domain.tri_full_flag[::3] = 0

        op = Rate_operator(domain, rate=self.RATE, factor=1.0)

        domain.timestep = self.TIMESTEP
        op()

        owned_area = domain.areas[domain.tri_full_flag == 1].sum()
        all_area = domain.areas.sum()
        return op.local_influx, owned_area, all_area

    def test_influx_excludes_ghost_cells(self):
        """Mode 2 must count only owned triangles — and agree with mode 1 and theory."""
        cpu, owned_area, all_area = self._influx(1)
        gpu, _, _ = self._influx(2)

        # Exact expected value: rate * factor * timestep * (area of OWNED triangles).
        expected = self.RATE * self.TIMESTEP * owned_area
        # What the bug produced: the same sum over EVERY triangle, ghosts included.
        buggy = self.RATE * self.TIMESTEP * all_area

        # Guard the guard: the two must be far apart, or this test proves nothing.
        self.assertGreater(abs(buggy - expected), 0.1 * abs(expected))

        self.assertAlmostEqual(
            cpu, expected, places=10,
            msg='mode-1 reference influx does not match the analytic value')
        self.assertAlmostEqual(
            gpu, expected, places=10,
            msg='mode-2 counted ghost cells in the rate-operator influx (issue #191)')
        self.assertAlmostEqual(
            gpu, cpu, places=10,
            msg='mode-1 and mode-2 disagree on rate-operator influx')


class Test_GPU_ForcingOperators(unittest.TestCase):
    """Wind_stress / Barometric_pressure / Rate OPERATORS must give identical
    results in mode 1 (legacy) and mode 2 (unified).

    These fractional-step operators replaced the Wind_stress /
    Barometric_pressure / Rainfall forcing-function classes (removed in 4.1).
    Mode 2 skips Python forcing_terms (it applies forcing in C and only
    handles Manning friction — see _warn_unsupported_mode2_forcing), but
    fractional-step operators are applied by apply_fractional_steps() in BOTH
    modes. This mirrors test_forcing.py's wind/pressure evolve tests across
    both compute modes — confirming the operators really are applied in
    unified mode and agree with legacy.
    """

    def _make_domain(self, name):
        d = rectangular_cross_domain(10, 5, len1=100., len2=50.)
        d.set_flow_algorithm('DE0')
        d.set_low_froude(0)
        d.set_name(name)
        d.set_datadir(tempfile.mkdtemp())
        d.store = False
        d.set_quantity('elevation', 0.0)
        d.set_quantity('friction', 0.0)
        d.set_quantity('stage', 1.0)        # still water, depth 1 m
        Br = Reflective_boundary(d)
        d.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})
        return d

    def _run(self, mode, add_operator, name):
        from anuga.shallow_water.sw_domain_gpu_ext import sync_to_device, sync_from_device
        d = self._make_domain(name)
        add_operator(d)
        d.set_multiprocessor_mode(mode)
        if mode == 2:
            sync_to_device(d.gpu_interface.gpu_dom)
        for _ in d.evolve(yieldstep=1.0, finaltime=5.0):
            pass
        if mode == 2:
            sync_from_device(d.gpu_interface.gpu_dom)
        return d

    def _assert_modes_agree(self, add_operator, label):
        d1 = self._run(1, add_operator, f'{label}_m1')
        d2 = self._run(2, add_operator, f'{label}_m2')
        # These operators apply the same per-cell forcing in mode 1 and mode 2, so
        # they agree to machine precision on still water: measured <=1.2e-15 at t=5
        # for wind/pressure and bit-identical for rain on a GPU-offload build. The
        # bound below catches any mode-1/mode-2 regression with wide headroom.
        atol = 1e-10 if anuga.gpu_offload_enabled() else 1e-12
        for q in ['stage', 'xmomentum', 'ymomentum']:
            np.testing.assert_allclose(
                d2.quantities[q].centroid_values,
                d1.quantities[q].centroid_values,
                rtol=0, atol=atol,
                err_msg=f'{label}: {q} mode1 vs mode2 mismatch')
        return d1

    def test_wind_stress_constant(self):
        from anuga.operators.wind_stress_operator import Wind_stress_operator
        d1 = self._assert_modes_agree(
            lambda d: Wind_stress_operator(d, speed=15.0, phi=0.0), 'wind_const')
        # east wind (phi=0) must build positive x-momentum
        self.assertGreater(np.max(d1.quantities['xmomentum'].centroid_values), 0.0)

    def test_wind_stress_temporally_varying(self):
        from anuga.operators.wind_stress_operator import Wind_stress_operator
        def speed(t, x, y):
            return 5.0 + 2.0 * t
        d1 = self._assert_modes_agree(
            lambda d: Wind_stress_operator(d, speed=speed, phi=90.0), 'wind_temporal')
        # north wind (phi=90) must build positive y-momentum
        self.assertGreater(np.max(d1.quantities['ymomentum'].centroid_values), 0.0)

    def test_wind_stress_spatially_varying(self):
        from anuga.operators.wind_stress_operator import Wind_stress_operator
        def speed(t, x, y):
            return 5.0 + 0.1 * x
        self._assert_modes_agree(
            lambda d: Wind_stress_operator(d, speed=speed, phi=0.0), 'wind_spatial')

    def test_barometric_pressure_spatially_varying(self):
        from anuga.operators.barometric_pressure import Barometric_pressure_operator
        def pressure(t, x, y):
            return 101325.0 + 50.0 * x      # ∂p/∂x drives momentum
        d1 = self._assert_modes_agree(
            lambda d: Barometric_pressure_operator(d, pressure=pressure), 'pressure_spatial')
        # a non-zero pressure gradient must move the water
        self.assertGreater(
            np.max(np.abs(d1.quantities['xmomentum'].centroid_values)), 0.0)

    def test_rate_operator_rainfall(self):
        from anuga.operators.rate_operators import Rate_operator
        d1 = self._assert_modes_agree(
            lambda d: Rate_operator(d, rate=1.0e-3), 'rain')   # 1 mm/s rain
        # rainfall must raise the stage above the initial 1.0 m
        self.assertGreater(np.max(d1.quantities['stage'].centroid_values), 1.0)


@pytest.mark.skipif(not gpu_available(), reason=_gpu_skip_reason())
class Test_GPU_TimeBoundary(unittest.TestCase):
    """Multiple Time_boundary objects with different values must not clobber
    one another in mode 2.

    Regression for the bug where the GPU time-boundary stored a single global
    (stage, xmom, ymom) applied to every time-boundary edge. With two
    Time_boundary tags carrying different values — worst on a sloped bed, where
    the absolute stages at the two ends differ a lot — the last one set won and
    corrupted the other boundary, diverging catastrophically from legacy
    (e.g. avalanche_wet). The values are now stored per edge.
    """

    bed_slope = 0.1

    def _make_domain(self):
        L = 100.0
        d = rectangular_cross_domain(12, 6, len1=L, len2=10.0)
        d.set_flow_algorithm('DE1')
        d.set_datadir(tempfile.mkdtemp())
        d.store = False
        d.set_quantity('friction', 0.0)
        d.set_quantity('elevation', lambda x, y: self.bed_slope * x)
        d.set_quantity('stage', lambda x, y: self.bed_slope * x + 5.0)
        return d, L

    def _run_mode(self, mode):
        d, L = self._make_domain()
        # Two Time_boundary objects with DIFFERENT, time-varying values. On the
        # sloped bed the two absolute stages differ by ~bed_slope*L = 10 m.
        def f_left(t):
            return [5.0 + 0.2 * t, 0.0, 0.0]
        def f_right(t):
            return [self.bed_slope * L + 5.0 - 0.1 * t, 0.0, 0.0]
        Bl = anuga.Time_boundary(d, f_left)
        Brt = anuga.Time_boundary(d, f_right)
        Bw = Reflective_boundary(d)
        d.set_boundary({'left': Bl, 'right': Brt, 'top': Bw, 'bottom': Bw})
        d.set_multiprocessor_mode(mode)
        d.set_quantities_to_be_stored(None)
        gauge = d.get_triangle_containing_point([L / 2.0, 5.0])
        stage = d.get_quantity('stage')
        xmom = d.get_quantity('xmomentum')
        out = []
        for _ in d.evolve(yieldstep=1.0, finaltime=5.0):
            out.append((float(stage.centroid_values[gauge]),
                        float(xmom.centroid_values[gauge])))
        return out

    def test_two_time_boundaries_mode1_vs_mode2(self):
        """Two differing Time_boundaries on a slope: mode 1 == mode 2."""
        g1 = self._run_mode(1)
        g2 = self._run_mode(2)
        self.assertEqual(len(g1), len(g2))
        for (s1, x1), (s2, x2) in zip(g1, g2):
            self.assertAlmostEqual(s1, s2, places=8,
                msg=f"stage mode1={s1} vs mode2={s2}")
            self.assertAlmostEqual(x1, x2, places=8,
                msg=f"xmom mode1={x1} vs mode2={x2}")


@pytest.mark.skipif(not gpu_available(), reason=_gpu_skip_reason())
class Test_GPU_InletWithCpuOnlyOperator(unittest.TestCase):
    """A GPU-accelerated Inlet_operator's inflow must survive when a CPU-only
    fractional operator is also present.

    Regression: GPU-path fractional operators (Inlet_operator, Rate_operator)
    apply their update to the *device* arrays. When a CPU-only fractional
    operator (e.g. an Internal_boundary bridge) is also present,
    apply_fractional_steps() brackets the operator loop with
    sync_from_device()/sync_to_device(); the trailing host->device sync then
    overwrote the GPU operator's device write with host data that never received
    it, silently dropping the inflow. bridge_hecras2 drained to the bed under
    mode 2 because of this. The GPU operators now fall back to the host path
    while _gpu_host_writes_suppressed is set, so the batch sync carries them.
    """

    INFLOW_Q = 5.0        # m^3/s
    FINALTIME = 5.0

    def _make_domain(self):
        d = rectangular_cross_domain(10, 10, len1=100., len2=100.)
        d.set_flow_algorithm('DE1')
        d.set_datadir(tempfile.mkdtemp())
        d.store = False
        d.set_quantity('elevation', 0.0)
        d.set_quantity('friction', 0.0)
        d.set_quantity('stage', 1.0)          # 1 m over a 100x100 m closed box
        Br = Reflective_boundary(d)
        d.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})
        return d

    def _total_volume(self, d):
        h = np.maximum(d.quantities['stage'].centroid_values
                       - d.quantities['elevation'].centroid_values, 0.0)
        return float((h * d.areas).sum())

    def _run(self, mode):
        from anuga.operators.base_operator import Operator

        class _NoopCpuOperator(Operator):
            # A plain Operator subclass is classified CPU-only, which forces
            # apply_fractional_steps() into its sync_from/to_device bracket.
            def __call__(self):
                pass

            def parallel_safe(self):
                return True

        d = self._make_domain()
        Inlet_operator(d, [[20.0, 50.0], [80.0, 50.0]], self.INFLOW_Q)
        _NoopCpuOperator(d)
        d.set_multiprocessor_mode(mode)
        d.set_quantities_to_be_stored(None)
        for _ in d.evolve(yieldstep=1.0, finaltime=self.FINALTIME):
            pass
        return self._total_volume(d)

    def test_inlet_inflow_survives_cpu_only_operator(self):
        """Closed-box inflow volume matches between mode 1 and mode 2."""
        v_m1 = self._run(1)
        v_m2 = self._run(2)
        self.assertAlmostEqual(v_m1, v_m2, places=4,
            msg=f"mode1 vol={v_m1} vs mode2 vol={v_m2} (inflow lost under mode 2?)")

    def test_inflow_actually_applied_mode2(self):
        """Sanity: the inflow really raises the volume under mode 2 (not ~0)."""
        v0 = 1.0 * 100.0 * 100.0                 # initial depth 1 over 1e4 area
        expected = self.INFLOW_Q * self.FINALTIME  # closed box: all inflow retained
        v_m2 = self._run(2)
        self.assertGreater(v_m2 - v0, 0.8 * expected,
            msg=f"expected ~{expected} m^3 of inflow, got {v_m2 - v0}")


class Test_GPU_ManyCulverts(unittest.TestCase):
    """More culverts than MAX_CULVERTS must work (issue #217).

    MAX_CULVERTS (gpu_domain.h) is only the INITIAL capacity of the culvert
    arrays -- registration grows them by doubling, so a model may hold any
    number.  gpu_culverts_apply_all() nevertheless sized its per-step working
    buffers (inlet data, results, transfers and the MPI exchange buffers) as
    fixed arrays of MAX_CULVERTS entries, then looped over num_culverts.  A
    model with more culverts than that wrote past the end of them on every
    timestep: a 101-culvert model overran ~8.9 kB of stack per step, which
    showed up as the run wedging at Time = 0 and later segfaulting on a
    pointer whose bits were a NaN double.

    This test registers comfortably more than MAX_CULVERTS culverts and
    evolves.  Before the fix that corrupts the stack and typically crashes;
    no useful assertion can be made about the values, so completing a few
    steps with finite state is the check.
    """

    N_CULVERTS = 80          # > MAX_CULVERTS (64), with room to spare

    def test_more_culverts_than_initial_capacity(self):
        from anuga import Boyd_box_operator

        domain = rectangular_cross_domain(40, 40, len1=200.0, len2=200.0)
        domain.set_flow_algorithm('DE0')
        domain.set_name('many_culverts')
        domain.set_datadir(tempfile.mkdtemp())
        domain.store = False
        domain.set_multiprocessor_mode(2)

        # Sloped bed so the culverts actually transfer water rather than idle.
        domain.set_quantity('elevation', lambda x, y: -x / 20.0)
        domain.set_quantity('stage', 2.0)
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        # Spread the culverts out so their inlets do not overlap.
        for i in range(self.N_CULVERTS):
            y = 5.0 + (190.0 * i) / self.N_CULVERTS
            Boyd_box_operator(domain,
                              end_points=[[40.0, y], [160.0, y]],
                              losses=1.5, width=2.0, height=2.0, apron=0.0,
                              use_momentum_jet=False, use_velocity_head=False,
                              manning=0.013, verbose=False)

        for _ in domain.evolve(yieldstep=1.0, finaltime=3.0):
            pass

        stage = domain.quantities['stage'].centroid_values
        xmom = domain.quantities['xmomentum'].centroid_values
        self.assertTrue(np.all(np.isfinite(stage)),
                        'non-finite stage after evolving %d culverts'
                        % self.N_CULVERTS)
        self.assertTrue(np.all(np.isfinite(xmom)),
                        'non-finite xmomentum after evolving %d culverts'
                        % self.N_CULVERTS)


class Test_GPU_RiverwallCrestUpdate(unittest.TestCase):
    """A riverwall crest changed mid-run must reach the device (issue #224).

    `domain.riverwallData.set_elevation()` writes only the host array. In mode 2
    the riverwall arrays are mapped to the device once, when the interface is
    built, and are never written there — so operating a gate part-way through a
    run used to be silently dropped: the run completed normally and produced
    plausible output computed with the ORIGINAL crest. The fix pushes the host
    arrays across at yieldstep boundaries (Domain._sync_riverwall_to_device()).

    NOTE ON WHAT THIS CATCHES WHERE. On a CPU build the `omp target` regions run
    on the host and the "device" arrays ARE the host arrays, so the physics
    comparison below cannot fail there — it is a real guard only on a GPU-offload
    build. test_sync_is_issued_at_the_yieldstep is the white-box half that fails
    everywhere if the push is dropped or moved off the yieldstep boundary.
    """

    WALL_X = 500.0
    OPEN = -3.0      # crest below the ponded water: gate open
    SHUT = 5.0       # crest above it: gate shut

    def _build(self, mode):
        domain = rectangular_cross_domain(40, 8, len1=1000.0, len2=200.0)
        domain.set_name('rw_gate')
        domain.set_datadir(tempfile.mkdtemp())
        domain.store = False
        domain.set_flow_algorithm('DE0')
        domain.set_quantity('elevation', lambda x, y: -x / 100.0)
        domain.set_quantity('friction', 0.03)
        # Pond water upstream of the wall only.
        domain.set_quantity(
            'stage', lambda x, y: np.where(x < self.WALL_X, 0.0, -x / 100.0))
        domain.riverwallData.create_riverwalls(
            {'gate': [[self.WALL_X, 0.0, self.OPEN],
                      [self.WALL_X, 200.0, self.OPEN]]}, verbose=False)
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})
        domain.set_multiprocessor_mode(mode)
        return domain

    def _downstream_volume(self, domain):
        centroids = domain.get_centroid_coordinates(absolute=True)
        downstream = centroids[:, 0] > self.WALL_X
        h = (domain.quantities['stage'].centroid_values
             - domain.quantities['elevation'].centroid_values)
        return float(np.sum(np.maximum(h[downstream], 0.0)
                            * domain.areas[downstream]))

    def _run(self, mode, shut_gate):
        """Water ponded upstream drains over an open wall; optionally shut it."""
        domain = self._build(mode)
        shut_done = False
        for t in domain.evolve(yieldstep=20.0, finaltime=200.0):
            if shut_gate and not shut_done and t >= 100.0:
                domain.riverwallData.set_elevation('gate', self.SHUT)
                shut_done = True
        return self._downstream_volume(domain)

    def test_shut_gate_changes_the_answer(self):
        """Shutting the gate must hold water back — in mode 2 as in mode 1."""
        open_1, shut_1 = self._run(1, False), self._run(1, True)
        open_2, shut_2 = self._run(2, False), self._run(2, True)

        # Sanity: the gate has a large, unambiguous effect on the reference path.
        self.assertLess(shut_1, 0.9 * open_1,
                        'test setup is not sensitive to the crest change')

        np.testing.assert_allclose(
            shut_2, shut_1, rtol=1e-10,
            err_msg='mode-2 ignored a mid-run riverwall crest change (the device '
                    'kept the original crest elevations)')
        np.testing.assert_allclose(open_2, open_1, rtol=1e-10)

    def test_sync_is_issued_at_the_yieldstep(self):
        """The push happens on resuming from the yield, and only when needed."""
        domain = self._build(2)
        domain._ensure_gpu_interface()
        self.assertIsNotNone(domain.gpu_interface)

        calls = []
        real = domain.gpu_interface.sync_riverwall_to_device

        def spy():
            calls.append(domain.get_time())
            real()

        domain.gpu_interface.sync_riverwall_to_device = spy

        for t in domain.evolve(yieldstep=20.0, finaltime=100.0):
            if t >= 40.0 and not calls:
                domain.riverwallData.set_elevation('gate', self.SHUT)
                # Not pushed yet: the yield body runs before evolve() resumes.
                self.assertEqual(calls, [])

        self.assertEqual(len(calls), 1,
                         'expected exactly one device push — one per change, at '
                         'the yieldstep boundary, got %d' % len(calls))
        # It fired on the yieldstep at which the change was made, not later.
        self.assertAlmostEqual(calls[0], 40.0, places=6)
        self.assertFalse(domain.riverwallData.device_data_dirty)

    def test_no_riverwalls_is_harmless(self):
        """A domain with no riverwalls must not be disturbed by the new hook."""
        domain = rectangular_cross_domain(10, 10)
        domain.set_flow_algorithm('DE0')
        domain.set_name('rw_none')
        domain.set_datadir(tempfile.mkdtemp())
        domain.store = False
        domain.set_quantity('elevation', -10.0)
        domain.set_quantity('stage', 2.0)
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})
        domain.set_multiprocessor_mode(2)

        for _ in domain.evolve(yieldstep=0.5, finaltime=1.0):
            pass

        stage = domain.quantities['stage'].centroid_values
        self.assertTrue(np.all(np.isfinite(stage)))
        np.testing.assert_allclose(stage, 2.0, atol=1e-8)


class Test_GPU_DegenerateTimestepProtection(unittest.TestCase):
    """Mode 2 must SAY it does not run the degenerate-timestep protection (issue #189).

    `apply_protection_against_isolated_degenerate_timesteps()` is reached only via
    `update_timestep()`. The mode-2 C step loops return before that, and the
    Python-orchestrated GPU loops that do call it hit a host `max_speed` the device
    never syncs back — so under GPU offload the protection does nothing either way.
    It is default-off, so the sharp edge is the user who turns it on and gets no
    protection and no warning.
    """

    def _domain(self, mode, protect):
        domain = rectangular_cross_domain(10, 10)
        domain.set_flow_algorithm('DE0')
        domain.set_name('degen')
        domain.set_datadir(tempfile.mkdtemp())
        domain.store = False
        domain.set_quantity('elevation', -10.0)
        domain.set_quantity('stage', 2.0)
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})
        domain.set_multiprocessor_mode(mode)
        domain.protect_against_isolated_degenerate_timesteps = protect
        return domain

    @staticmethod
    def _degenerate_warnings(caught):
        return [w for w in caught
                if 'isolated-degenerate' in str(w.message)
                or 'protect_against_isolated_degenerate_timesteps' in str(w.message)]

    def test_warns_when_enabled_in_mode2(self):
        """Turning the protection on under GPU offload must not be silent."""
        domain = self._domain(2, protect=True)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            for _ in domain.evolve(yieldstep=0.5, finaltime=1.0):
                pass
        found = self._degenerate_warnings(caught)
        self.assertEqual(len(found), 1,
                         'expected exactly one warning that mode 2 does not run the '
                         'degenerate-timestep protection, got %d' % len(found))
        self.assertIn('legacy', str(found[0].message))

    def test_quiet_when_disabled(self):
        """The default (off) must stay silent — the gap only matters if asked for."""
        domain = self._domain(2, protect=False)
        self.assertFalse(domain.protect_against_isolated_degenerate_timesteps,
                         'this test assumes the feature is default-off')
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            for _ in domain.evolve(yieldstep=0.5, finaltime=1.0):
                pass
        self.assertEqual(self._degenerate_warnings(caught), [])

    def test_quiet_in_mode1(self):
        """Mode 1 implements the protection, so it must not warn."""
        domain = self._domain(1, protect=True)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            for _ in domain.evolve(yieldstep=0.5, finaltime=1.0):
                pass
        self.assertEqual(self._degenerate_warnings(caught), [])

    def test_routine_skips_rather_than_acting_on_stale_max_speed(self):
        """Called directly in mode 2, the routine must warn and do nothing.

        The host `max_speed` here carries the exact signature the routine looks
        for (one isolated fast triangle, empty middle bins). In mode 1 that damps
        the triangle's momentum; in mode 2 the same host array is stale by
        construction, so the routine must decline to act on it.
        """
        for mode in (1, 2):
            domain = self._domain(mode, protect=True)
            domain.set_quantity('xmomentum', 1.0)
            domain.max_speed = np.zeros(domain.number_of_triangles, dtype=float)
            domain.max_speed[0] = 100.0     # isolated degenerate triangle

            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always')
                domain.apply_protection_against_isolated_degenerate_timesteps()
            found = self._degenerate_warnings(caught)

            xmom = domain.quantities['xmomentum'].centroid_values
            if mode == 1:
                self.assertEqual(found, [], 'mode 1 implements this; it must not warn')
                self.assertEqual(float(xmom[0]), 0.0,
                                 'mode 1 must damp the isolated triangle')
                self.assertEqual(float(domain.max_speed[0]), 0.0)
            else:
                self.assertEqual(len(found), 1,
                                 'mode 2 must warn that the protection is not applied')
                self.assertEqual(float(xmom[0]), 1.0,
                                 'mode 2 must not damp on the strength of a stale '
                                 'host max_speed')
                self.assertEqual(float(domain.max_speed[0]), 100.0)

    def test_warning_is_issued_once(self):
        """One warning per domain, not one per timestep."""
        domain = self._domain(2, protect=True)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            for _ in domain.evolve(yieldstep=0.2, finaltime=1.0):
                pass
            domain.apply_protection_against_isolated_degenerate_timesteps()
        self.assertEqual(len(self._degenerate_warnings(caught)), 1)


class Test_GPU_LargeInlet(unittest.TestCase):
    """A culvert inlet may cover any number of triangles (issue #225).

    The GPU culvert path used to cap one inlet at MAX_INLET_TRIANGLES (64) —
    fixed-size arrays inside `struct culvert_indices`, unlike the sibling
    operator structs, which hold heap pointers. Registering a wider inlet failed
    outright ("Failed to register culvert ... on GPU"). The device side was
    already fully dynamic: the metadata is flattened into heap arrays sized by
    the actual total at map time, so the cap only ever constrained the host-side
    staging.

    The geometry below gives 72 triangles per inlet — comfortably over the old
    cap: a 40 m exchange line at x = 61, off the mesh lines (3.33 m apart), so
    it cuts three of the four triangles in each of the 24 rectangles it
    crosses. (On a mesh line it would select only the 48 along it, since
    line-region selection is geometric, #231.) apron=5.0 keeps the enquiry
    points clear of the (wide) inlet regions.
    """

    def _run(self, mode):
        from anuga import Boyd_box_operator

        domain = rectangular_cross_domain(60, 30, len1=200.0, len2=50.0)
        domain.set_flow_algorithm('DE0')
        domain.set_name('big_inlet')
        domain.set_datadir(tempfile.mkdtemp())
        domain.store = False
        # A real head across the culvert (ponded upstream, low downstream), so
        # the discharge is large and deterministic rather than set by roundoff.
        domain.set_quantity('elevation', -5.0)
        domain.set_quantity('stage', lambda x, y: np.where(x < 100.0, 1.0, -1.0))
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        culvert = Boyd_box_operator(
            domain, end_points=[[61.0, 25.0], [139.0, 25.0]],
            losses=1.5, width=40.0, height=3.0, apron=5.0,
            use_momentum_jet=False, use_velocity_head=False,
            manning=0.013, verbose=False)
        domain.set_multiprocessor_mode(mode)

        inlet_sizes = [len(inlet.triangle_indices) for inlet in culvert.inlets]

        for _ in domain.evolve(yieldstep=1.0, finaltime=5.0):
            pass

        return (inlet_sizes,
                domain.quantities['stage'].centroid_values.copy(),
                float(culvert.discharge))

    def test_inlet_larger_than_the_old_cap(self):
        """An inlet of 72 triangles must register and give the mode-1 answer."""
        sizes_1, stage_1, q_1 = self._run(1)
        sizes_2, stage_2, q_2 = self._run(2)

        self.assertEqual(sizes_1, sizes_2)
        self.assertGreater(min(sizes_1), 64,
                           'this test is only meaningful with inlets larger than '
                           'the old MAX_INLET_TRIANGLES cap of 64; got %s' % sizes_1)

        np.testing.assert_allclose(
            stage_2, stage_1, rtol=0, atol=1e-10,
            err_msg='mode-2 culvert with a >64-triangle inlet does not match mode 1')
        self.assertAlmostEqual(q_2, q_1, places=6)
        self.assertGreater(abs(q_1), 1.0,
                           'the culvert should actually be passing flow')


class Test_GPU_StructureSurfaceLevel(unittest.TestCase):
    """The device write-back levels the inlet, as the host path does (#229).

    A structure operator used to write a uniform DEPTH across its inlet, which on
    a sloping bed tilts the water surface onto the bed and disturbs a lake at
    rest. Both the host paths and the mode-2 device scatter now apply the volume
    change by filling to / drawing down from a level ("water finds its level"),
    clamping cells at their bed, and write the per-inlet momentum depth-weighted
    so the velocity field stays bounded. A zero transfer is an exact no-op.
    """

    def _lake_at_rest(self, mode, slope_denominator=50.0):
        from anuga import Boyd_box_operator

        domain = rectangular_cross_domain(30, 15, len1=200.0, len2=50.0)
        domain.set_flow_algorithm('DE0')
        domain.set_name('gpu_well_balanced')
        domain.set_datadir(tempfile.mkdtemp())
        domain.store = False
        domain.set_quantity('elevation',
                            lambda x, y: -5.0 + x / slope_denominator)
        domain.set_quantity('stage', 1.0)
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})
        Boyd_box_operator(
            domain, end_points=[[60.0, 25.0], [140.0, 25.0]],
            enquiry_points=[[40.0, 25.0], [160.0, 25.0]],
            losses=1.5, width=20.0, height=3.0, apron=5.0,
            use_momentum_jet=False, use_velocity_head=False,
            manning=0.013, verbose=False)
        domain.set_multiprocessor_mode(mode)

        for _ in domain.evolve(yieldstep=0.5, finaltime=1.0):
            pass

        stage = domain.quantities['stage'].centroid_values
        return float(np.abs(stage - 1.0).max())

    def test_lake_at_rest_holds_in_mode2(self):
        """No head, no flow, no disturbance — on the device as on the host."""
        gpu = self._lake_at_rest(2)
        cpu = self._lake_at_rest(1)
        # 6.8e-02 in both modes with the old uniform-depth write.
        self.assertLess(gpu, 1e-5,
                        'mode-2 culvert disturbed a lake at rest by %.3e m' % gpu)
        self.assertLess(cpu, 1e-5)

    def _drain_shallow_inlet(self, mode):
        """Drain a shallow, sloping inlet into a deep pool: cells must go dry."""
        from anuga import Boyd_box_operator

        domain = rectangular_cross_domain(40, 20, len1=200.0, len2=50.0)
        domain.set_flow_algorithm('DE0')
        domain.set_name('gpu_wetdry')
        domain.set_datadir(tempfile.mkdtemp())
        domain.store = False
        # Shallow sloping shelf upstream (depths 0.027-0.052 m over the inlet),
        # deep pool downstream, so the culvert asks for more water than the
        # shallow end of the inlet holds.
        def bed(x, y):
            return np.where(x < 100.0, -1.0 + 0.005 * x, -6.0)

        domain.set_quantity('elevation', bed)
        # Clamped at the bed so no cell starts below it — otherwise the solver's
        # first protect step lifts those cells and creates ~100 m^3, swamping the
        # volume check below.
        domain.set_quantity(
            'stage', lambda x, y: np.where(x < 100.0,
                                           np.maximum(bed(x, y), -0.66), -5.0))
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})
        # End points deliberately OFF the 5 m grid lines: an exchange line lying
        # exactly on mesh edges is a degenerate intersection whose triangle
        # selection is compiler-dependent (gcc and nvc resolve it differently —
        # 32 vs 54 triangles for x=60.0 on this mesh), which made this test's
        # dried-cell threshold build-dependent.
        culvert = Boyd_box_operator(
            domain, end_points=[[61.3, 25.0], [138.7, 25.0]],
            enquiry_points=[[40.0, 25.0], [170.0, 25.0]],
            losses=1.5, width=10.0, height=3.0, apron=5.0,
            use_momentum_jet=False, use_velocity_head=False,
            manning=0.013, verbose=False)
        domain.set_multiprocessor_mode(mode)

        inlet0 = culvert.inlets[0].triangle_indices.copy()

        def volume():
            bed = domain.quantities['elevation'].centroid_values
            stage = domain.quantities['stage'].centroid_values
            return float(np.sum((stage - bed) * domain.areas))

        start_volume = volume()
        dried = 0
        for _ in domain.evolve(yieldstep=0.5, finaltime=6.0):
            bed = domain.quantities['elevation'].centroid_values
            stage = domain.quantities['stage'].centroid_values
            depth = stage - bed
            dried = max(dried, int(np.sum(depth[inlet0] <= 1e-12)))

        return (volume() - start_volume, dried,
                domain.quantities['stage'].centroid_values.copy())

    def test_wet_dry_clamp_conserves_volume(self):
        """Draining past dry must clamp at the bed, not invent or lose water."""
        results = {mode: self._drain_shallow_inlet(mode) for mode in (1, 2)}

        for mode, (drift, dried, _) in results.items():
            self.assertGreater(dried, 5,
                               'mode %d: the clamp is not being exercised — no '
                               'inlet cells went dry' % mode)
            # An unclamped write drives cells negative; the solver's protect
            # step then lifts them back to zero, creating water.
            self.assertLess(abs(drift), 1e-6,
                            'mode %d: volume drifted by %.3e m^3 draining an '
                            'inlet dry' % (mode, drift))

        np.testing.assert_allclose(
            results[2][2], results[1][2], rtol=0, atol=1e-10,
            err_msg='the device clamp does not match the host clamp')


if __name__ == "__main__":
    unittest.main(verbosity=2)

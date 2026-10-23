"""
GPU interface using OpenMP target offloading.

This provides a similar interface to sw_domain_cuda.py but uses
OpenMP target offloading via the sw_domain_gpu_ext Cython module.

Key differences from CUDA approach:
- Arrays stay in NumPy (OpenMP handles device mapping)
- No explicit array copies needed
- Uses MPI for multi-GPU (vs single-GPU CUDA approach)
"""

import numpy as np

class GPU_OMP_interface:
    """
    Interface for GPU-accelerated shallow water solver using OpenMP target offloading.

    Usage:
        domain.set_multiprocessor_mode(2)  # This will create the interface
        # ... normal evolve loop ...
    """

    def __init__(self, domain):
        """
        Initialize the GPU interface.

        Parameters
        ----------
        domain : anuga.shallow_water.Domain
            The domain object (should already be partitioned for MPI)
        """
        self.domain = domain
        self.gpu_dom = None
        self.initialized = False
        self.boundaries_initialized = False

    def setup(self):
        """
        Initialize GPU domain and map arrays to device.

        Call this once before starting the evolve loop.
        """
        if self.initialized:
            return

        from anuga.shallow_water.sw_domain_gpu_ext import (
            init_gpu_domain, map_to_gpu,
            init_reflective_boundary, init_dirichlet_boundary, init_transmissive_boundary,
            init_transmissive_n_zero_t_boundary, init_time_boundary,
            init_absorbing_wave_boundary, init_characteristic_wave_boundary,
            init_flather_boundary,
        )

        # Ensure work arrays (including edge_flux_type) are allocated before
        # the C GPU struct tries to dereference them.
        if hasattr(self.domain, '_ensure_work_arrays'):
            self.domain._ensure_work_arrays()

        # Initialize GPU domain structure
        self.gpu_dom = init_gpu_domain(self.domain)

        # Initialize all supported boundary types BEFORE mapping to GPU
        # This allows boundary arrays to be mapped together with main arrays,
        # avoiding OpenMP data region issues
        if self.domain.boundary_map is not None:
            self._init_boundaries()

        # Map arrays to GPU memory (persistent for simulation)
        # This now includes all boundary arrays if initialized above
        map_to_gpu(self.gpu_dom)

        self.initialized = True

    def finalize(self):
        """
        Clean up GPU resources.

        Call this when done with GPU computation.
        """
        if not self.initialized:
            return

        from anuga.shallow_water.sw_domain_gpu_ext import (
            unmap_from_gpu, finalize_gpu_domain
        )

        unmap_from_gpu(self.gpu_dom)
        finalize_gpu_domain(self.gpu_dom)
        self.initialized = False

    def _init_boundaries(self):
        """Initialize GPU boundary structures from domain.boundary_map."""
        if self.boundaries_initialized:
            return
        if self.domain.boundary_map is None:
            return

        from anuga.shallow_water.sw_domain_gpu_ext import (
            init_reflective_boundary, init_dirichlet_boundary, init_transmissive_boundary,
            init_transmissive_n_zero_t_boundary, init_time_boundary, init_file_boundary,
            init_absorbing_wave_boundary, init_characteristic_wave_boundary,
            init_flather_boundary,
        )

        init_reflective_boundary(self.gpu_dom, self.domain)
        init_dirichlet_boundary(self.gpu_dom, self.domain)
        init_transmissive_boundary(self.gpu_dom, self.domain)
        init_transmissive_n_zero_t_boundary(self.gpu_dom, self.domain)
        init_time_boundary(self.gpu_dom, self.domain)
        init_file_boundary(self.gpu_dom, self.domain)
        init_absorbing_wave_boundary(self.gpu_dom, self.domain)
        init_characteristic_wave_boundary(self.gpu_dom, self.domain)
        init_flather_boundary(self.gpu_dom, self.domain)

        self.boundaries_initialized = True

    def ensure_boundaries_initialized(self):
        """
        Verify GPU boundaries are initialized.

        Note: As of 2026-01-18, boundaries MUST be set before calling
        set_multiprocessor_mode(2). This method is kept for backwards
        compatibility but should not be needed.
        """
        if not self.boundaries_initialized:
            raise RuntimeError(
                "GPU boundaries not initialized. This should not happen - "
                "boundaries must be set before calling set_multiprocessor_mode(2)."
            )

    def sync_to_device(self):
        """Sync centroid values from host to device."""
        from anuga.shallow_water.sw_domain_gpu_ext import sync_to_device
        sync_to_device(self.gpu_dom)

    def sync_from_device(self):
        """Sync centroid values from device to host."""
        from anuga.shallow_water.sw_domain_gpu_ext import sync_from_device
        sync_from_device(self.gpu_dom)

    def sync_boundary_values(self):
        """Sync boundary values from host to device."""
        from anuga.shallow_water.sw_domain_gpu_ext import sync_boundary_values
        sync_boundary_values(self.gpu_dom)

    def sync_riverwall_to_device(self):
        """Sync riverwall elevations / hydraulic properties from host to device."""
        from anuga.shallow_water.sw_domain_gpu_ext import sync_riverwall_to_device
        sync_riverwall_to_device(self.gpu_dom)

    def sync_edge_values_from_device(self):
        """Sync edge values from device to host."""
        from anuga.shallow_water.sw_domain_gpu_ext import sync_edge_values_from_device
        sync_edge_values_from_device(self.gpu_dom)

    def exchange_ghosts(self):
        """Exchange ghost cells between MPI ranks."""
        from anuga.shallow_water.sw_domain_gpu_ext import exchange_ghosts
        exchange_ghosts(self.gpu_dom)

    # =========================================================================
    # Kernel wrappers - these match the interface expected by shallow_water_domain.py
    # =========================================================================

    def extrapolate_second_order_edge_sw_kernel(self, domain, distribute_to_vertices=False):
        """
        Extrapolate centroid values to edges.

        Note: distribute_to_vertices is ignored (we only do edge-based extrapolation)
        """
        from anuga.shallow_water.sw_domain_gpu_ext import extrapolate_second_order_gpu
        extrapolate_second_order_gpu(self.gpu_dom)

    def compute_fluxes_ext_central_kernel(self, domain, timestep):
        """
        Compute fluxes and return minimum timestep.

        Returns
        -------
        float
            Local minimum timestep (MPI allreduce needed for global min)
        """
        from anuga.shallow_water.sw_domain_gpu_ext import compute_fluxes_gpu
        return compute_fluxes_gpu(self.gpu_dom)

    def protect_against_infinitesimal_and_negative_heights_kernel(self, domain):
        """Protect against negative water depths."""
        from anuga.shallow_water.sw_domain_gpu_ext import protect_gpu
        return protect_gpu(self.gpu_dom)

    def update_conserved_quantities_kernel(self, domain, timestep):
        """Update conserved quantities with explicit and semi-implicit terms."""
        from anuga.shallow_water.sw_domain_gpu_ext import update_conserved_quantities_gpu
        update_conserved_quantities_gpu(self.gpu_dom, timestep)
        # GPU protect kernel handles negative cells; return (count, volume) = (0, 0.0)
        # to match the openmp wrapper's (num_negative_cells, negative_volume) shape.
        return 0, 0.0

    def backup_conserved_quantities_kernel(self, domain):
        """Backup centroid values for RK2."""
        from anuga.shallow_water.sw_domain_gpu_ext import backup_conserved_quantities_gpu
        backup_conserved_quantities_gpu(self.gpu_dom)

    def saxpy_conserved_quantities_kernel(self, domain, a, b):
        """RK2 combination: Q = a*Q_current + b*Q_backup."""
        from anuga.shallow_water.sw_domain_gpu_ext import saxpy_conserved_quantities_gpu
        saxpy_conserved_quantities_gpu(self.gpu_dom, a, b)

    def saxpy3_conserved_quantities_kernel(self, domain, a, b, c):
        """RK3 final combination: Q = (a*Q_current + b*Q_backup) / c."""
        from anuga.shallow_water.sw_domain_gpu_ext import saxpy3_conserved_quantities_gpu
        saxpy3_conserved_quantities_gpu(self.gpu_dom, a, b, c)

    def manning_friction_kernel(self, domain):
        """Apply Manning friction (semi-implicit)."""
        from anuga.shallow_water.sw_domain_gpu_ext import manning_friction_gpu
        manning_friction_gpu(self.gpu_dom)

    # =========================================================================
    # FLOP Counter API (Gordon Bell Performance Profiling)
    # =========================================================================

    def flop_counters_reset(self):
        """
        Reset all FLOP counters to zero.

        Call this at the start of the profiling period.
        """
        from anuga.shallow_water.sw_domain_gpu_ext import flop_counters_reset
        flop_counters_reset(self.gpu_dom)

    def flop_counters_enable(self, enable=True):
        """
        Enable or disable FLOP counting.

        Parameters
        ----------
        enable : bool
            True to enable counting, False to disable
        """
        from anuga.shallow_water.sw_domain_gpu_ext import flop_counters_enable
        flop_counters_enable(self.gpu_dom, enable)

    def flop_counters_start_timer(self):
        """
        Start the FLOP counter timer.

        Call this at the start of the profiling period.
        """
        from anuga.shallow_water.sw_domain_gpu_ext import flop_counters_start_timer
        flop_counters_start_timer(self.gpu_dom)

    def flop_counters_stop_timer(self):
        """
        Stop the FLOP counter timer.

        Call this at the end of the profiling period.
        """
        from anuga.shallow_water.sw_domain_gpu_ext import flop_counters_stop_timer
        flop_counters_stop_timer(self.gpu_dom)

    def flop_counters_get_total(self):
        """
        Get total FLOP count across all kernels.

        Returns
        -------
        int
            Total FLOPs executed since last reset
        """
        from anuga.shallow_water.sw_domain_gpu_ext import flop_counters_get_total
        return flop_counters_get_total(self.gpu_dom)

    def flop_counters_get_flops(self):
        """
        Get FLOP rate (FLOPs per second).

        Call flop_counters_stop_timer first to record elapsed time.

        Returns
        -------
        float
            FLOP/s rate
        """
        from anuga.shallow_water.sw_domain_gpu_ext import flop_counters_get_flops
        return flop_counters_get_flops(self.gpu_dom)

    def flop_counters_print(self):
        """
        Print detailed FLOP counter summary to stdout.

        Shows per-kernel breakdown and total performance in GFLOP/s.
        """
        from anuga.shallow_water.sw_domain_gpu_ext import flop_counters_print
        flop_counters_print(self.gpu_dom)

    def flop_counters_get_stats(self):
        """
        Get all FLOP counter statistics as a dictionary.

        Returns
        -------
        dict
            Dictionary with per-kernel FLOPs, total, elapsed time, and GFLOP/s
        """
        from anuga.shallow_water.sw_domain_gpu_ext import flop_counters_get_stats
        return flop_counters_get_stats(self.gpu_dom)

    # =========================================================================
    # Global FLOP Counter API (MPI reduction for multi-GPU)
    # =========================================================================

    def flop_counters_get_global_total(self):
        """
        Get global total FLOP count summed across all MPI ranks/GPUs.

        Returns
        -------
        int
            Total FLOPs across all GPUs since last reset
        """
        from anuga.shallow_water.sw_domain_gpu_ext import flop_counters_get_global_total
        return flop_counters_get_global_total(self.gpu_dom)

    def flop_counters_get_global_flops(self):
        """
        Get global FLOP rate (FLOPs per second) across all GPUs.

        Returns
        -------
        float
            Global FLOP/s rate
        """
        from anuga.shallow_water.sw_domain_gpu_ext import flop_counters_get_global_flops
        return flop_counters_get_global_flops(self.gpu_dom)

    def flop_counters_print_global(self):
        """
        Print global FLOP counter summary to stdout (rank 0 only).

        Shows per-kernel breakdown summed across all GPUs with percentages.
        """
        from anuga.shallow_water.sw_domain_gpu_ext import flop_counters_print_global
        flop_counters_print_global(self.gpu_dom)

    def flop_counters_get_global_stats(self):
        """
        Get global FLOP counter statistics as a dictionary (MPI reduced).

        Returns
        -------
        dict
            Dictionary with global total, GFLOP/s, num_gpus, and per-GPU average
        """
        from anuga.shallow_water.sw_domain_gpu_ext import flop_counters_get_global_stats
        return flop_counters_get_global_stats(self.gpu_dom)

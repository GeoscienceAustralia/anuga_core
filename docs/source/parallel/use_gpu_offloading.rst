.. _use_gpu_offloading:

.. currentmodule:: anuga

GPU Acceleration (OpenMP target offloading)
============================================

ANUGA includes an experimental GPU backend that offloads the computationally
intensive parts of the evolve loop — flux computation, extrapolation, friction,
and momentum updates — to a GPU using
`OpenMP target offloading <https://www.openmp.org/spec-html/5.0/openmpsu57.html>`_.
No CUDA or Python GPU libraries are required.

.. include:: /appendices/gpu_acknowledgement.inc

.. note::

   The GPU backend is experimental and under active development.  The API
   (``set_compute_mode``, operator support) may change in future
   releases.  For production runs where result reproducibility matters,
   validate GPU output against ``'legacy'`` before switching.

.. seealso::

   :ref:`install_gpu` — how to build ANUGA with GPU support (NVIDIA HPC SDK /
   ``nvc``). This page covers *using* the backend once it is built.


Hardware and compiler requirements
------------------------------------

+--------------------+-------------------------------------------------------+
| Requirement        | Notes                                                 |
+====================+=======================================================+
| NVIDIA GPU         | Tested on NVIDIA GPUs.  AMD ROCm GPUs supported in    |
|                    | principle via LLVM ``libomptarget``.                  |
+--------------------+-------------------------------------------------------+
| Compiler           | GCC 12+ with offload targets *or* LLVM/Clang with     |
|                    | ``libomptarget``.  Standard GCC on conda-forge does   |
|                    | **not** include GPU offload targets.                  |
+--------------------+-------------------------------------------------------+
| Driver / toolkit   | CUDA toolkit matching your GPU driver version.        |
+--------------------+-------------------------------------------------------+
| OpenMP ≥ 4.5       | Required for ``omp target`` directives.               |
+--------------------+-------------------------------------------------------+

The ANUGA GPU extensions are built only when the compiler supports
``-fopenmp-targets``.  On a standard conda or pip install the GPU extension
(``sw_domain_gpu_ext``) is absent; all GPU code paths fall back to the
CPU-OpenMP path automatically.


CPU-only mode (default)
-----------------------

When ``sw_domain_gpu_ext`` is not available (standard install) or when
``OMP_TARGET_OFFLOAD=disabled`` is set, ``mode=2`` still works correctly — it
runs the same C RK loop on the CPU without any device transfer.  This is called
**CPU_ONLY_MODE** and is the mode used by the test suite in CI.

You can check which mode is active::

   from anuga.shallow_water.sw_domain_gpu_ext import gpu_available
   print(gpu_available())   # True only when a real GPU target is available


Quick start
-----------

The only change required to an existing ANUGA script is one call before the
evolve loop:

.. code-block:: python

   import anuga

   domain = anuga.rectangular_cross_domain(200, 200, len1=10000., len2=5000.)
   domain.set_flow_algorithm('DE0')

   # ... set quantities, boundaries, operators as normal ...

   domain.set_compute_mode('unified')   # enable GPU mode

   for t in domain.evolve(yieldstep=60., finaltime=3600.):
       domain.print_timestepping_statistics()

``set_compute_mode('unified')`` initialises the GPU interface and switches the
evolve loop to a C-side Runge-Kutta loop that keeps all data resident on the
GPU between timesteps.  Data is only transferred back to the host at each
``yieldstep`` for Python I/O.


How it works
-------------

.. code-block:: text

   Python evolve() call
     │
     ├─ map_to_gpu()         host → device (once, at first yieldstep)
     │
     └─ for each yieldstep:
          ├─ C RK loop       all kernels run on device (no CPU round-trips)
          │    ├─ extrapolate_second_order_gpu
          │    ├─ evaluate_boundary_gpu  (all boundary types)
          │    ├─ compute_fluxes_gpu     (riverwall weir/orifice in-kernel)
          │    ├─ update_conserved_quantities_gpu
          │    ├─ protect_against_negatives_gpu
          │    ├─ manning_friction_gpu
          │    └─ [fractional step operators — see below]
          │
          └─ sync_from_device()   device → host (once per yieldstep)

The data path is batched: a single ``H→D`` transfer at startup, a single
``D→H`` transfer at each yieldstep.  The inner timestep loop never touches host
memory, so GPU→CPU bandwidth is not a bottleneck.


Flow algorithm compatibility
-----------------------------

The GPU backend supports the DE (discontinuous elevation) solver family only.

+-------------------+-----------------+
| Algorithm         | GPU support     |
+===================+=================+
| ``'DE0'``         | ✓ Full          |
+-------------------+-----------------+
| ``'DE1'``         | ✓ Full          |
+-------------------+-----------------+
| ``'DE2'`` (RK3)   | ✓ Full (SSP-RK3)|
+-------------------+-----------------+
| ``'1_5'``, etc.   | ✗ CPU only      |
+-------------------+-----------------+

Use ``domain.set_flow_algorithm('DE0')`` (default) or ``'DE1'``/``'DE2'`` for
GPU-accelerated runs.


Supported operators
--------------------

Operators that are registered with the GPU domain execute entirely on the
device.  All other operators execute on the CPU via the standard fractional
step mechanism; the domain is automatically synced to/from host before and
after each CPU operator call.

**Fully GPU-accelerated operators**

+------------------------------------------+--------------------+
| Operator                                 | GPU implementation |
+==========================================+====================+
| ``Rate_operator``                        | Batched; up to     |
|                                          | 64 instances       |
+------------------------------------------+--------------------+
| ``Inlet_operator``                       | Batched; up to     |
|                                          | 32 instances       |
+------------------------------------------+--------------------+
| ``Boyd_box_operator``                    | GPUCulvertManager; |
|                                          | up to 64 total     |
+------------------------------------------+--------------------+
| ``Boyd_pipe_operator``                   | GPUCulvertManager  |
+------------------------------------------+--------------------+
| ``Weir_orifice_trapezoid_operator``      | GPUCulvertManager  |
+------------------------------------------+--------------------+
| Riverwalls (``create_riverwalls``)       | In-kernel (flux    |
|                                          | computation only)  |
+------------------------------------------+--------------------+

**Operators that fall back to CPU**

+--------------------------------------------+----------------------------------+
| Operator                                   | Reason                           |
+============================================+==================================+
| ``Bed_shear_erosion_operator``             | No GPU kernel yet                |
+--------------------------------------------+----------------------------------+
| ``Kinematic_viscosity_operator``           | No GPU kernel yet                |
+--------------------------------------------+----------------------------------+
| ``Weir_orifice_trapezoid_operator``        | GPU via GPUCulvertManager        |
| (parallel)                                 | (cross-boundary MPI supported)   |
+--------------------------------------------+----------------------------------+

When a CPU fractional-step operator runs, ANUGA automatically syncs data from
device before the operator executes and back to device afterwards.  This adds
one round-trip per operator call per timestep, so minimise the number of CPU
operators when GPU performance is important.


Operator limits
---------------

The GPU domain uses dynamically-allocated arrays for operator registration,
so there is no hard limit on the number of ``Rate_operator``, ``Inlet_operator``,
or culvert operators.  The arrays start at their default capacity and double
automatically as more operators are added.

The one remaining fixed limit is **64 triangles per inlet face** — if a
culvert inlet region covers more than 64 mesh triangles, registration will
fail with a ``RuntimeError``.  Coarsen the mesh near the structure or split
the inlet region to stay within this limit.


Supported boundary conditions
------------------------------

All standard boundary types are supported in GPU mode:

- ``Reflective_boundary``
- ``Transmissive_boundary``, ``Transmissive_n_momentum_zero_t_momentum_set_stage_boundary``
- ``Dirichlet_boundary``
- ``File_boundary`` / ``Field_boundary`` (per-edge values pushed to device each sub-step)
- ``Time_boundary``, ``Time_stage_zero_momentum_boundary``
- ``Absorbing_wave_boundary`` (wave scalar pushed from Python; ghost state computed on device)
- ``Characteristic_wave_boundary`` (wave scalar pushed from Python; nonlinear characteristic kernel on device)
- ``Flather_external_stage_zero_velocity_boundary`` (exterior stage scalar pushed from Python; Blayo & Debreu characteristic decomposition kernel on device)

Custom boundary classes not in the above list are evaluated on the CPU with an
automatic device sync.


Parallel (MPI + GPU)
---------------------

The GPU backend is compatible with MPI domain decomposition.  Each MPI rank
uses one GPU device.  Halo exchange between ranks uses CPU-side buffers;
GPU-aware MPI (direct GPU-to-GPU transfer over NVLink/InfiniBand) is planned
for a future release (G2.2).

.. code-block:: bash

   # 4 MPI ranks, each using one GPU
   mpiexec -np 4 python my_parallel_gpu_script.py

In the script each rank calls ``set_compute_mode('unified')`` independently
after domain decomposition.


Checking the build and device
-------------------------------

.. code-block:: python

   import anuga

   # Check whether sw_domain_gpu_ext is present (required for mode=2)
   try:
       from anuga.shallow_water import sw_domain_gpu_ext
       print("GPU extension available")
   except ImportError:
       print("GPU extension not available — running in CPU-only mode")

   domain = anuga.rectangular_cross_domain(100, 100)
   domain.set_flow_algorithm('DE0')
   domain.set_compute_mode('unified')

   # Check device memory (prints estimate and device info if CUDA/HIP available)
   # Raises RuntimeError if estimated memory exceeds device capacity.
   # No-op in CPU_ONLY_MODE.
   from anuga.shallow_water.sw_domain_gpu_ext import check_device_memory
   check_device_memory(domain.gpu_interface.gpu_dom)

   print(f"Mode: {domain.get_compute_mode()}")   # unified


Performance tips
-----------------

1. **Use DE0 or DE1.**  The DE2 (SSP-RK3) scheme does 3× the kernel work per
   timestep in exchange for higher accuracy; only use it if accuracy requires it.

2. **Minimise CPU fractional-step operators.**  Each CPU operator adds a
   full device sync round-trip.  Prefer GPU-registered operators (see table).

3. **Use large meshes.**  GPU occupancy is low for small meshes (< ~50 000
   triangles).  The GPU path is typically faster than OpenMP CPU for meshes
   larger than ~200 000 triangles on modern hardware.

4. **Set OMP_NUM_THREADS=1 for GPU runs.**  The host-side RK loop is
   single-threaded in mode=2; extra OpenMP threads on the host compete for
   memory bandwidth without benefit::

      OMP_NUM_THREADS=1 OMP_TARGET_OFFLOAD=mandatory python my_script.py

   ``OMP_TARGET_OFFLOAD=mandatory`` causes an error if the GPU target is not
   found, rather than silently falling back to CPU.

5. **Yieldstep interval.**  The D→H sync at each yieldstep is fast (~10 ms for
   a 1M-triangle domain), so yieldstep granularity has minimal impact on
   throughput.

6. **Active-set stepping, on mostly-dry domains.**  See below.


Active-set stepping
-------------------

A mostly-dry domain spends most of its time integrating cells that cannot
change: dry, and surrounded by dry.  Active-set stepping skips them.

.. code-block:: python

   domain.set_use_active_set(True)     # off by default

On mostly-dry flood problems -- a dam break, a levee breach, a storm surge
over dry land -- this is the largest single speedup measured for the GPU path.
It costs nothing in accuracy: a cell is skipped only when it is dry and its
whole two-ring neighbourhood is dry, so there is nothing to integrate.  Water
added by rate operators is never lost either; a newly wetted cell activates on
the next step.

.. warning::

   **The flag is a request, not a guarantee.**  Several configurations turn it
   off again, and they do so with a log notice rather than an error -- so a run
   that quietly gets no speedup looks exactly like one that does.  Check the
   log, or ask the domain:

   .. code-block:: python

      frac, rebuilds = domain.get_active_set_stats()
      # frac == 1.0 means the active set never engaged

   ``get_active_set_stats()`` returns the mean active fraction and the number
   of rebuilds, and reports ``1.0`` whenever the active set did not engage.

The flag is ignored, with a notice, when:

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Condition
     - Why
   * - compute mode is ``'legacy'``
     - only the unified (GPU) path implements it
   * - flow algorithm is ``DE2`` (rk3)
     - implemented for ``DE0``/euler, ``DE1``/rk2 and ``DE_ader2``
   * - passive tracers or sediment are registered
     - the active path uses the single-solve scatter flux kernel, which does
       not advect tracers; only the cell-based kernel does
   * - the domain has riverwalls
     - the scatter kernel does not carry them
   * - Manning friction varies across a cell
     - likewise
   * - running under MPI
     - serial domains only for now

Enabling it alongside **rate operators** is allowed and correct, but widespread
rainfall wets the whole mesh, which activates every cell and erodes the speedup
toward zero.  That combination also logs a notice.

.. note::

   Active-set stepping switches the flux kernel to single-solve scatter mode.
   Results differ from the default cell-based kernel only at floating-point
   roundoff.


Troubleshooting
----------------

``ImportError: cannot import name 'init_gpu_domain'``
    The ``sw_domain_gpu_ext`` module was not built.  Rebuild with a compiler
    that supports ``-fopenmp-targets`` (GCC 12+ with offload targets or
    LLVM with ``libomptarget``).

``RuntimeError: GPU inlet face too large``
    A culvert inlet region covers more than 64 mesh triangles.  Coarsen
    the mesh near the structure or split the inlet region.

``RuntimeError: GPU device memory insufficient``
    The domain is too large for the available device memory.  Use a smaller
    mesh, reduce the number of stored quantities, or use a GPU with more VRAM.

Results differ between mode=1 and mode=2
    Differences larger than ~1e-12 (machine epsilon for ``float64``) indicate
    a bug.  Please open an issue at
    `anuga-community/anuga_core <https://github.com/anuga-community/anuga_core/issues>`_
    with the minimal reproducer.  Expected numerical difference is zero for
    CPU_ONLY_MODE and ≤ 1e-10 for real GPU hardware (due to non-associative
    floating-point reduction order).


See Also
---------

- :ref:`use_parallel_openmp` — CPU OpenMP threading (mode=1)
- :ref:`use_parallel_mpi` — MPI domain decomposition
- `ANUGA GPU development plan <https://github.com/anuga-community/anuga_core/blob/develop/claude/GPU_DEVELOPMENT_PLAN.md>`_

.. _profiling_gpu:

.. currentmodule:: anuga

Profiling the GPU build
=======================

.. note::

   This applies to a **GPU (nvc) build** running in mode 2 (``-mpm 2`` /
   ``compute_mode = "unified"`` / ``set_gpu_offload(True)``). See
   :doc:`install_gpu` and :ref:`compute_modes`.

The NVIDIA HPC SDK ships the two NVIDIA profilers, plus a memory checker:

- **Nsight Systems** (``nsys``) — a whole-application **timeline**: kernel
  launches, host↔device transfers, CPU gaps, MPI. Start here to see *where* time
  goes and which kernels dominate.
- **Nsight Compute** (``ncu``) — a **per-kernel** deep dive: occupancy, memory
  throughput, roofline, warp stalls. Use it on the one or two kernels ``nsys``
  identifies as hot.
- **compute-sanitizer** — races / out-of-bounds / uninitialised device memory.

Setup
-----

Put the SDK tools on ``PATH`` (define ``NVC_BASE`` first):

.. code-block:: bash

   export NVC_BASE=/opt/nvidia/hpc_sdk/Linux_x86_64/26.3
   export PATH=$NVC_BASE/compilers/bin:$PATH

Then ``ncu``, ``nsys``, ``ncu-ui``, ``nsys-ui`` and ``compute-sanitizer`` are on
the path. (Adjust the ``26.3`` version to your installed SDK.)

**GPU performance counters.** ``ncu`` reads restricted GPU counters, which by
default require elevated privileges — run it with ``sudo`` (as below), or have an
admin allow non-root profiling once via the driver option
``NVreg_RestrictProfilingToAdminUsers=0``. ``nsys`` timeline profiling does **not**
need ``sudo``.

Whole-run timeline (Nsight Systems)
-----------------------------------

.. code-block:: bash

   nsys profile -o anuga_timeline -f true \
     python run_small_towradgi.py -mpm 2 -ft 200

This writes ``anuga_timeline.nsys-rep``. List the hottest kernels from it with:

.. code-block:: bash

   nsys stats --report cuda_gpu_kern_sum anuga_timeline.nsys-rep

The kernel-summary table gives the exact **kernel names** (see below) to hand to
``ncu -k``.

Per-kernel deep dive (Nsight Compute)
-------------------------------------

Profile a single launch of the flux kernel (the main compute kernel), skipping
warmup, collecting the full metric set:

.. code-block:: bash

   sudo $NVC_BASE/compilers/bin/ncu \
     -k "nvkernel_core_compute_fluxes_central_F1L936_36" \
     -c 1 -s 5 --set=full \
     -o my_profile_anuga_fluxes -f \
     $(which python) run_small_towradgi.py -mpm 2 -ft 200

Flags:

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Flag
     - Meaning
   * - ``-k <regex>``
     - Only profile kernels whose name matches (here, the flux kernel).
   * - ``-s 5``
     - **Skip** the first 5 matching launches (past JIT/warmup, into steady state).
   * - ``-c 1``
     - **Collect** 1 launch after the skip (per-kernel profiling replays the
       launch many times, so keep the count small).
   * - ``--set=full``
     - The full metric set (roofline, memory workload, occupancy, …). Thorough
       but slow; use ``--set=basic`` or ``--set=roofline`` for a quick pass.
   * - ``-o <name>``
     - Report file ``<name>.ncu-rep``.
   * - ``-f``
     - Overwrite an existing report.

Because ``ncu`` replays each profiled kernel many times, always narrow it with
``-k`` **and** a small ``-c`` — never profile a whole ANUGA run unfiltered.

Kernel names
------------

nvc names each OpenMP ``target`` (offload) region
``nvkernel_<function>_F<file>L<line>_<index>``. For example
``nvkernel_core_compute_fluxes_central_F1L936_36`` is the ``omp target`` region
at **line 936** of ``anuga/shallow_water/gpu/core_kernels.c``, inside
``core_compute_fluxes_central`` (the ``F1`` is the compiler's file index and
``_36`` its kernel index). The build prints these ``function: line`` regions as
it compiles; ``nsys stats`` (above) or an unfiltered ``ncu`` also lists them.

Viewing the reports
-------------------

Open interactively (GUI needed):

.. code-block:: bash

   ncu-ui  my_profile_anuga_fluxes.ncu-rep     # Nsight Compute report
   nsys-ui anuga_timeline.nsys-rep             # Nsight Systems timeline

Or copy the ``.ncu-rep`` / ``.nsys-rep`` file to a workstation and open it in the
Nsight GUIs there. For a quick terminal summary of an ``ncu`` report:

.. code-block:: bash

   ncu --import my_profile_anuga_fluxes.ncu-rep --page details | less

.. include:: /appendices/gpu_acknowledgement.inc

.. _compute_modes:

.. currentmodule:: anuga

Compute modes: legacy vs unified (and GPU offload)
==================================================

.. note::

   **Standard users can skip this page.** By default ANUGA runs in the
   **legacy** compute path on the CPU, which is stable and needs no
   configuration. This appendix is for advanced users who want the newer
   unified C kernels or GPU offloading.

ANUGA v4.0 exposes two *orthogonal* knobs:

1. **Per-domain compute path** — *how* each step is computed: the **legacy**
   path (mode 1: the historical ``sw_domain_openmp`` kernels plus Python-level
   operators) or the **unified** path (mode 2: a single set of C kernels shared
   by CPU and GPU).
2. **Process-global GPU offload** — *where* the unified kernels run: on the CPU
   (multicore OpenMP) or offloaded to a GPU. This is a process-wide setting, not
   a per-domain one, because OpenMP target offload is a process-level runtime
   state.

The two combine: ``cpu`` = unified + offload off; ``gpu`` = unified + offload on.
"cpu" and "gpu" are *compositions* of these knobs, not extra modes.


Per-domain compute path
-----------------------

.. code-block:: python

   domain.set_compute_mode('legacy')    # = set_multiprocessor_mode(1)
   domain.set_compute_mode('unified')   # = set_multiprocessor_mode(2)
   domain.get_compute_mode()            # 'legacy' or 'unified'

- **legacy** (default) — the historical CPU path. Widely exercised and the
  reference for correctness.
- **unified** — the mode-2 C kernels. On a CPU-only build these run as multicore
  OpenMP; on a GPU build they can offload (see below). Results match legacy to
  round-off; validate before relying on unified for production.

The process-wide default can be set with the environment variable
``ANUGA_DEFAULT_COMPUTE_MODE=unified`` (default ``legacy``). This is honoured for
**serial** runs only — under MPI (``numprocs > 1``) new domains stay ``legacy``.


Process-global GPU offload
--------------------------

Call these *before* the first ``evolve()`` (and before building the first
unified domain, since its arrays are mapped to the chosen device at init):

.. code-block:: python

   anuga.set_gpu_offload(True)      # offload unified kernels to a GPU (GPU build only)
   anuga.set_gpu_offload(False)     # run unified on the CPU
   anuga.gpu_offload_enabled()      # resolved offload state
   anuga.set_omp_num_threads(16)    # OpenMP thread count for the whole process
   domain.compute_capabilities()    # {gpu_offload, num_gpu_devices, mpi, modes}

On a standard conda/pip install the GPU extension is absent and offload requests
degrade gracefully to the CPU path — ``set_gpu_offload(True)`` warns and returns
``False`` rather than failing. GPU offload requires a GPU-enabled build (see
:ref:`use_gpu_offloading` for hardware/compiler requirements and build steps).


Strict mode: fail instead of falling back
-----------------------------------------

Every fallback above is a warning by default, which is right for interactive
use and wrong for a batch queue: a "GPU job" that quietly runs on the CPU, or
a 'unified' run that quietly became 'legacy', burns its allocation before
anyone reads the log. Strict mode turns each of those into a ``RuntimeError``
raised at setup, before any time is spent:

.. code-block:: python

   anuga.set_gpu_offload(True, strict=True)          # no offload support -> error
   domain.set_compute_mode('unified', strict=True)   # see below

or, for a whole job without touching the script:

.. code-block:: bash

   export ANUGA_STRICT_COMPUTE_MODE=1

With strict set on a domain, three things that would otherwise downgrade
silently are refused:

* the MPI fallback from 'unified' to 'legacy' on a build whose unified
  kernels were compiled without MPI;
* running on the CPU through the OpenMP target host fallback when GPU offload
  is enabled but no device was found;
* boundary types the device path cannot evaluate, which otherwise run on the
  host every step with a host<->device sync each time.

The setting is per domain and persists across later ``set_compute_mode``
calls until you pass ``strict=False``. The environment variable sets the
default for every domain and for ``set_gpu_offload``.


When is parallelism worth it? (mesh-size guidance)
--------------------------------------------------

Neither OpenMP threads nor GPU offload speed up a *small* mesh. The shallow
water solver is memory-bandwidth-bound, and each timestep is a sequence of short
parallel regions. On a small domain the fixed per-region costs — OpenMP
fork/join, and for the GPU the kernel-launch latency — dominate the actual
arithmetic, so extra threads give little and the GPU can be *slower than a single
CPU core*. Parallelism only pays off once the mesh is large enough to amortise
that overhead.

The numbers below are illustrative measurements on **one machine** (32-core CPU
vs an RTX 5070 Laptop GPU, single rank, ``unified`` mode, steady-state timing
with one-time setup excluded). **Treat the crossover sizes, not the exact
ratios, as the takeaway** — both shift with CPU core count, memory bandwidth,
GPU class, and problem.

**OpenMP thread scaling** (speedup over a single core):

.. list-table::
   :header-rows: 1
   :widths: 20 20 20

   * - Triangles
     - 8-thread speedup
     - Notes
   * - ~400
     - 0.4×
     - *slower* — fork/join overhead wins
   * - ~1,600
     - 1.4×
     - barely worthwhile
   * - ~6,400
     - 2.6×
     -
   * - ~25,600
     - 4.4×
     -
   * - ~100,000+
     - ~5–6×
     - near-peak (bandwidth-bound plateau, below the ideal 8×)

**GPU offload crossover** (speedup, same meshes):

.. list-table::
   :header-rows: 1
   :widths: 20 20 20

   * - Triangles
     - vs 1 core
     - vs 8 threads
   * - ~100–400
     - 0.5× (slower)
     - 0.6× (slower)
   * - ~1,600
     - 1.1×
     - 0.7× (slower)
   * - ~6,400
     - 3.4×
     - ~1.0× (tie)
   * - ~25,600
     - 7.9×
     - 1.5×
   * - ~100,000
     - 12×
     - 1.8×
   * - ~400,000
     - 15×
     - 2.9× (still climbing)

**Rule of thumb:**

- Small exploratory meshes (a few hundred to a few thousand triangles, e.g. a
  quick Jupyter notebook) — stay on **legacy, single-threaded**. Threads give
  almost nothing and the GPU is counterproductive.
- From **~10k triangles**, raise ``anuga.set_omp_num_threads()`` for a useful CPU
  speedup.
- Reach for the **GPU above ~25k triangles** versus a multicore CPU (or above
  ~1.5k versus a single core); its lead keeps growing into the hundreds of
  thousands of triangles, which is where these solvers are meant to run.

.. note::

   OpenMP thread count is **process-wide**, so ``anuga.set_omp_num_threads(n)``
   applies to every domain in the session — including ones already constructed
   (each domain's ``omp_num_threads`` reflects the live process-wide value).
   Setting it once at the top of a notebook is enough;
   ``anuga.get_omp_num_threads()`` reports the current value.


Command-line flags
------------------

The standard argument parser (used by the benchmark/example scripts) exposes:

- ``-mpm 1|2`` / ``--multiprocessor_mode`` — legacy vs unified
- ``-nt N`` / ``--num_threads`` — OpenMP thread count
- ``-go`` / ``-ngo`` — GPU offload on / off
- ``-ro <method>`` — mesh reordering for cache locality (e.g. ``metis_rcm``,
  ``hilbert``)


.. seealso::

   :ref:`use_gpu_offloading`
      Hardware/compiler requirements, the GPU build, and CPU-only fallback.

   :ref:`use_parallel_openmp`, :ref:`use_parallel_mpi`
      Multi-threaded and multi-process parallelism (independent of the
      legacy/unified choice).

.. include:: /appendices/gpu_acknowledgement.inc

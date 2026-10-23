.. _install_gpu:

.. currentmodule:: anuga

Installing for GPU (NVIDIA HPC SDK / nvc)
=========================================

.. note::

   This is only needed to **run on a GPU**. A standard conda/pip install runs
   ANUGA on the CPU (including the unified ``mode=2`` kernels as multicore
   OpenMP). See :ref:`compute_modes` for how to *use* GPU offload once built.

What you actually need
----------------------

ANUGA's GPU extension (``sw_domain_gpu_ext``) is written in **standard OpenMP
target offloading** — it is not CUDA code and is not tied to any one vendor.
The real requirement is therefore just:

   **a C compiler with working OpenMP offloading support for your GPU.**

In practice, **``nvc`` from the NVIDIA HPC SDK is the best option at the
moment**, and is the toolchain ANUGA is routinely built and tested with:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Compiler
     - Status
   * - **``nvc``** (NVIDIA HPC SDK)
     - **Recommended.** Mature OpenMP target offloading for NVIDIA GPUs; what
       the rest of this page uses.
   * - GCC (``nvptx`` offload backend)
     - **Does not currently work** — it hits an internal compiler error (ICE)
       on the ANUGA kernels. This is why a stock conda/pip install does **not**
       include the GPU extension.
   * - Others (LLVM/Clang offload, AMD ``AOMP``/ROCm, Intel ``icx``)
     - **Untested** with ANUGA. Because the kernels are plain OpenMP ``target``
       code these are feasible in principle, but expect some porting work.

The rest of this page covers the recommended ``nvc`` route.

Requirements
------------

- A GPU and a matching driver (for the ``nvc`` route: an NVIDIA GPU and a
  matching CUDA driver).
- **A compiler with OpenMP offloading support** — in practice the
  **NVIDIA HPC SDK** (which provides ``nvc``).
- A conda environment created via ``tools/install_miniforge.sh`` (see the
  :doc:`developer install </installation/install_anuga_developers>`).

One-time setup (Ubuntu)
-----------------------

Install the NVIDIA HPC SDK (~5 GB):

.. code-block:: bash

   curl -fsSL https://developer.download.nvidia.com/hpc-sdk/ubuntu/DEB-GPG-KEY-NVIDIA-HPC-SDK \
     | sudo gpg --dearmor -o /usr/share/keyrings/nvidia-hpcsdk-archive-keyring.gpg
   echo 'deb [signed-by=/usr/share/keyrings/nvidia-hpcsdk-archive-keyring.gpg] https://developer.download.nvidia.com/hpc-sdk/ubuntu/amd64 /' \
     | sudo tee /etc/apt/sources.list.d/nvhpc.list
   sudo apt-get update -y && sudo apt-get install -y nvhpc

Build ANUGA for GPU
-------------------

**Easiest — the helper script** (auto-detects ``nvc``, clears any stale build
dir, builds, and runs the isolated GPU tests):

.. code-block:: bash

   bash tools/install_anuga_nvc.sh
   # options via env vars, e.g. build for an A100 with Python 3.13:
   PY=3.13 GPU_ARCH=cc80 bash tools/install_anuga_nvc.sh

``GPU_ARCH`` defaults to the compute capability of the GPU in the machine you
are building on, read from ``nvidia-smi``, so you normally do not set it. Set it
explicitly when building for a *different* card than the one present (e.g. on a
login node for a compute node with other GPUs).

**Manual build:**

.. code-block:: bash

   CC=$(which nvc) pip install --no-build-isolation -e . \
       -Csetup-args=-Dgpu_offload=true \
       -Csetup-args=-Dgpu_arch=cc120

Pick ``gpu_arch`` for your card:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - ``gpu_arch``
     - GPU
   * - ``cc120``
     - RTX 50-series (Blackwell)
   * - ``cc90``
     - H100
   * - ``cc80``
     - A100
   * - ``cc86``
     - RTX 30-series (Ampere)
   * - ``cc89``
     - RTX 40-series / L4 (Ada)
   * - ``cc75``
     - T4 / RTX 20-series (Turing)
   * - ``cc70``
     - V100

Several may be given at once for a binary that runs on any of them, e.g.
``-Dgpu_arch=cc80,cc86,cc90``. Note that the CUDA bundled with recent NVIDIA HPC
SDKs dropped Volta, so ``cc70`` needs an older CUDA selected alongside it:
``-Dgpu_arch=cuda12.9,cc70,cc80``.

.. warning::

   **Build for the right architecture.** The build contains only the
   architectures you ask for. Requesting the wrong one still compiles cleanly
   and then fails at *every kernel launch* — which appears as the GPU tests
   reporting ``CRASH`` after an apparently successful install, not as a build
   error.

   Check what a build actually contains::

       cuobjdump --list-elf $(python -c \
           "import anuga.shallow_water.sw_domain_gpu_ext as m; print(m.__file__)") \
           | grep -o 'sm_[0-9]*' | sort -u

   and compare with your card::

       nvidia-smi --query-gpu=name,compute_cap --format=csv

   ``compute_cap`` 8.6 means you want ``cc86``. ``tools/install_anuga_nvc.sh``
   performs this check for you and reports a mismatch before running the tests.

.. warning::

   meson reads ``CC`` only on the **first** configure of a build directory. When
   switching between a gcc (CPU) and an nvc (GPU) build you **must** remove the
   build dir first, or meson keeps the old compiler and the ``gpu_offload=true``
   build aborts:

   .. code-block:: bash

      rm -rf build/cp*

   (``tools/install_anuga_nvc.sh`` does this automatically.)

Verify
------

Run the GPU test file with the isolated runner (one process per test — required
on a GPU build; see the runner's help for why):

.. code-block:: bash

   anuga_run_isolated_tests    # defaults to test_DE_gpu_omp.py

Switching back to a CPU-only build
----------------------------------

.. code-block:: bash

   rm -rf build/cp*
   pip install --no-build-isolation -e .    # or -Csetup-args=-Dgpu_offload=false

.. seealso::

   :ref:`compute_modes`
      Using the built GPU: ``set_gpu_offload`` / compute-mode selection.

   :ref:`use_gpu_offloading`
      More on the offloading backend and CPU-only fallback.

.. include:: /appendices/gpu_acknowledgement.inc

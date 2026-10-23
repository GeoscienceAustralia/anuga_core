.. _use_parallel_openmp:

.. currentmodule:: anuga

OpenMP Parallelisation
======================

From version 3.2, ANUGA can use multiple CPU threads to speed up the
computationally intensive parts of the evolve loop (flux computation, friction,
and momentum updates).  This is done using
`OpenMP <https://www.openmp.org/>`_, a widely-used API for shared-memory
parallel programming in C, C++, and Fortran.

By default, ANUGA uses **one thread**.  Setting ``OMP_NUM_THREADS`` or calling
``domain.set_omp_num_threads()`` enables multi-threading without any other
changes to your script.


Checking your ANUGA version
-----------------------------

.. code-block:: python

    import anuga
    print(anuga.__version__)   # should be 3.2 or later for OpenMP support


Controlling the number of threads
----------------------------------

**Environment variable (recommended)**

Set ``OMP_NUM_THREADS`` before launching Python:

.. code-block:: bash

    export OMP_NUM_THREADS=4
    python my_anuga_script.py

You can add this line to your shell profile (``~/.bashrc``, ``~/.zshrc``)
to make it permanent, or set it inline for a single run:

.. code-block:: bash

    OMP_NUM_THREADS=8 python my_anuga_script.py

**Inside the script**

Call ``domain.set_omp_num_threads()`` after creating the domain and before
the evolve loop:

.. code-block:: python

    import anuga

    domain = anuga.rectangular_cross_domain(200, 200)
    domain.set_omp_num_threads(4)   # use 4 threads

    for t in domain.evolve(yieldstep=1.0, duration=100.0):
        domain.print_timestepping_statistics()

The method reads ``OMP_NUM_THREADS`` from the environment if no argument is
passed, and defaults to 1 if neither is set.

.. note::

    The optimal thread count depends on your hardware.  A good starting point
    is the number of physical cores on the machine.  Hyper-threaded (logical)
    cores typically give little or no speedup for ANUGA's memory-bound loops.


If the build was compiled without OpenMP
-----------------------------------------

ANUGA's multi-threading is provided by the compiled extension
``sw_domain_openmp_ext``.  If this extension was not built with OpenMP
support (or was not built at all), calling ``domain.set_omp_num_threads()``
will raise an ``ImportError``::

    ImportError: cannot import name 'set_omp_num_threads'
                 from 'anuga.shallow_water.sw_domain_openmp_ext'

In this case the simulation still runs correctly — it simply uses a single
thread.

**How to check**

.. code-block:: python

    import anuga
    domain = anuga.rectangular_cross_domain(10, 10)
    domain.set_omp_num_threads(4, verbose=True)

If OpenMP is active you will see::

    Setting omp_num_threads to 4

If you see an ``ImportError``, the OpenMP extension is missing or was built
without OpenMP.


Installing ANUGA with OpenMP support
--------------------------------------

**conda-forge (recommended)**

Pre-built conda-forge packages are compiled with OpenMP enabled on Linux,
macOS, and Windows:

.. code-block:: bash

    conda install -c conda-forge anuga

**Building from source on Linux**

GCC includes OpenMP by default.  The standard install command is sufficient:

.. code-block:: bash

    conda activate anuga_env_3.12
    pip install --no-build-isolation -v .

**Building from source on macOS**

Apple's Clang does not include OpenMP.  Install ``libomp`` via Homebrew
first:

.. code-block:: bash

    brew install libomp
    conda activate anuga_env_3.12
    pip install --no-build-isolation -v .

**Building from source on Windows**

Use the conda-forge MinGW compilers, which include OpenMP:

.. code-block:: bash

    conda install -c conda-forge m2w64-gcc
    pip install --no-build-isolation -v .


OpenMP with MPI
----------------

OpenMP and MPI can be used together (hybrid parallelism).  Each MPI rank
runs its own OpenMP thread pool independently.  A typical hybrid launch
using 4 MPI ranks with 4 OpenMP threads each (16 cores total):

.. code-block:: bash

    OMP_NUM_THREADS=4 mpiexec -np 4 python my_parallel_script.py

.. note::

    Some MPI launchers do not forward environment variables to child
    processes.  If threads do not activate under ``mpiexec``, set
    ``OMP_NUM_THREADS`` explicitly inside the script before creating the
    domain::

        import os
        os.environ['OMP_NUM_THREADS'] = '4'
        import anuga


GPU acceleration (experimental)
--------------------------------

An experimental GPU backend is included in the ``develop`` branch.  It uses
**OpenMP target offloading** to run flux, friction, and momentum kernels on a
GPU without requiring Python-level changes to your script:

.. code-block:: python

   domain.set_compute_mode('unified')   # enable GPU mode

See :ref:`use_gpu_offloading` for hardware requirements, supported operators,
slot limits, and troubleshooting.


.. seealso::

   :ref:`api_domain` — full API of ``Domain.set_omp_num_threads`` and the other
   ``Domain`` methods.

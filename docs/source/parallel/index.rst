.. currentmodule:: anuga


Parallelisation
===============

ANUGA supports two parallel approaches that require minimal changes to a
standard sequential script:

- **OpenMP** — multi-threaded shared-memory parallelism; enable by setting
  the thread count before creating the domain.
- **MPI (distribute)** — multi-process distributed-memory parallelism;
  create the domain on rank 0 then call :func:`distribute` to partition it
  across ranks.

The two approaches can be combined.

.. toctree::
   :maxdepth: 1

   use_parallel_openmp
   use_parallel_mpi

For long parallel runs, see :ref:`checkpointing` for how to save and restart
simulations from periodic checkpoints.

To diagnose load imbalance caused by wet/dry triangle distributions, see
:ref:`load_balance_statistics`.

.. seealso::

   :ref:`advanced_parallel`
      GPU offloading, mesh-first distribution, offline partitioning, and
      uniform mesh refinement — for large HPC production runs.

   `ANUGA User Manual — Chapter 12: Parallel Simulation
   <https://github.com/anuga-community/anuga_user_manual>`_
   covers the MPI domain decomposition in depth, including ghost cell
   communication, scalability benchmarks, and HPC cluster setup.

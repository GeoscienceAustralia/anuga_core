.. _advanced_parallel:

.. currentmodule:: anuga

Advanced Parallel Methods
=========================

The methods here build on the basic MPI workflow described in
:doc:`use_parallel_mpi`.  They are intended for large production runs, HPC
cluster deployments, or situations where the standard ``distribute`` approach
is not sufficient (e.g. meshes too large to fit in memory on a single node).

.. toctree::
   :maxdepth: 1

   use_distribute_basic_mesh
   use_sequential_domain_io
   use_sequential_mesh_io
   use_mesh_refinement
   use_gpu_offloading

.. seealso::

   :doc:`../parallel/index`
      The simple parallel workflow (OpenMP and basic MPI distribute) that most
      users should start with.

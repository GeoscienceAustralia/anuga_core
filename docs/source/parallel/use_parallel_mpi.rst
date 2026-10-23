.. _use_parallel_mpi:

.. currentmodule:: anuga

MPI Distribute Domain-first parallel workflow (``distribute``)
===============================================================


Choosing an MPI distribution strategy
--------------------------------------

Six functions and function-pairs are available for distributing a domain
across MPI ranks.  The right choice depends on mesh size and whether
initial conditions need to be stored in the partition files.

.. list-table::
   :header-rows: 1
   :widths: 32 18 50

   * - Function
     - Rank-0 peak memory
     - When to use
   * - :func:`distribute`
     - Full domain + quantities
     - Default choice.  Rank 0 builds the complete ``Domain`` (mesh +
       all quantities), then sends submeshes to other ranks.  Simple
       drop-in for serial scripts.
   * - :func:`distribute_collaborative`
     - Full domain + quantities (shared)
     - Same interface as ``distribute`` but rank 0 broadcasts topology
       via MPI shared memory (one copy per node instead of one per rank).
       Use when several ranks share a node and rank-0 memory is tight.
   * - :func:`distribute_basic_mesh`
     - Topology only
     - Rank 0 builds a lightweight :class:`Basic_mesh` (no quantity
       arrays), distributes topology, then every rank sets its own
       quantities.  Best when quantity arrays would exhaust rank-0 memory.
   * - :func:`distribute_basic_mesh_collaborative`
     - Topology only (shared)
     - Like ``distribute_basic_mesh`` but topology is broadcast via shared
       memory.  Use for the largest meshes where even per-rank copies of
       the topology are expensive.
   * - :func:`sequential_distribute_dump` / :func:`sequential_distribute_load`
     - Full domain + quantities (preprocessing only)
     - Offline two-phase workflow: build and partition the full domain once
       (including quantities), write one set of pickle + npy files per rank,
       then load and evolve.  Nothing is held on rank 0 at runtime.
       See :ref:`use_sequential_domain_io`.
   * - :func:`sequential_mesh_dump` / :func:`sequential_mesh_load`
     - Topology only (preprocessing only)
     - Offline two-phase workflow: partition the mesh topology only (no
       quantities), write one NetCDF4 file per rank, then load and set
       quantities per-rank at runtime.  Smallest files; ideal for very
       large meshes or ensemble runs.  See :ref:`use_sequential_mesh_io`.

**Decision guide**

.. code-block:: text

    Does rank 0 have enough RAM to hold the full Domain (mesh + quantities)?
    │
    ├─ Yes ──► Run immediately (no preprocessing)?
    │           ├─ Yes ──► Do multiple ranks share the same node?
    │           │           ├─ No  ──► distribute()
    │           │           └─ Yes ──► distribute_collaborative()
    │           └─ No  ──► Preprocess offline with quantities:
    │                        sequential_distribute_dump / sequential_distribute_load
    │                        (see use_sequential_domain_io)
    │
    └─ No  ──► Does rank 0 have enough RAM for topology only (no quantities)?
                ├─ Yes ──► Do multiple ranks share the same node?
                │           ├─ No  ──► distribute_basic_mesh()
                │           └─ Yes ──► distribute_basic_mesh_collaborative()
                └─ No  ──► Preprocess offline (mesh only):
                            sequential_mesh_dump / sequential_mesh_load
                            (see use_sequential_mesh_io)

As a rough guide, a mesh of N triangles with P quantities requires
approximately ``8 × N × P`` bytes for quantity arrays alone on rank 0
(double precision).  Topology (coordinates + connectivity) is
``~56 × N`` bytes.


All six approaches differ in how much memory and work rank 0 must do
before the evolve loop begins.

.. list-table:: Distribution function comparison
   :header-rows: 1
   :widths: 20 16 16 16 16 16

   * - Feature
     - ``distribute``
     - ``distribute_collaborative``
     - ``distribute_basic_mesh`` / ``_collaborative``
     - ``sequential_distribute_dump`` / ``_load``
     - ``sequential_mesh_dump`` / ``_load``
   * - Mesh built on rank 0
     - Full ``Domain`` with quantities
     - Full ``Domain`` with quantities
     - ``Basic_mesh`` only (no quantities)
     - Preprocessing only (offline, no MPI)
     - Preprocessing only (offline, no MPI)
   * - Quantities set
     - On rank 0, then distributed
     - On rank 0, then distributed
     - Per-rank after distribution
     - Stored in files; loaded per-rank
     - Per-rank after load
   * - Topology broadcast
     - Point-to-point per rank
     - Shared memory per node
     - Point-to-point / shared memory
     - Read from disk per rank
     - Read from disk per rank
   * - Peak rank-0 memory at runtime
     - O(N) mesh + O(N) quantities
     - O(N) mesh + O(N) quantities
     - O(N) mesh only
     - Per-rank partition only
     - Per-rank partition only
   * - Partition file format
     - —
     - —
     - —
     - Pickle + npy arrays
     - NetCDF4
   * - Best for
     - Small–medium meshes, simple scripts
     - Many ranks per node, large meshes
     - Very large meshes; quantities set from functions
     - Offline partitioning with quantities included (fast restart)
     - Very large meshes or ensemble runs; quantities set per-rank


``distribute``
--------------

The simplest workflow.  Rank 0 builds a full ``Domain`` (mesh *and*
quantities), then ``distribute`` partitions it and sends each rank its
submesh.  All ranks can set boundary conditions and operators after the
call.

.. code-block:: python

   import anuga

   if anuga.myid == 0:
       domain = anuga.rectangular_cross_domain(1000, 1000, len1=10.0, len2=10.0)
       domain.set_quantity('elevation', lambda x, y: x / 10.0)
       domain.set_quantity('stage', expression='elevation + 0.2')
   else:
       domain = None

   # Rank 0 partitions; all ranks receive their submesh
   domain = anuga.distribute(domain)

   Br = anuga.Reflective_boundary(domain)
   domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

   rain = anuga.Rate_operator(domain, rate=lambda t: 0.001, factor=1.0)

   for t in domain.evolve(yieldstep=1.0, finaltime=10.0):
       if anuga.myid == 0:
           domain.print_timestepping_statistics()

   domain.sww_merge()
   anuga.finalize()

The main limitation is that rank 0 must hold the complete mesh **and**
all quantity arrays before any MPI communication starts.  For very large
meshes this can exhaust available memory on rank 0.


``distribute_collaborative``
-----------------------------

A drop-in replacement for ``distribute`` that reduces peak memory use and
inter-rank data movement.  The API is identical — rank 0 still builds a
full ``Domain`` with quantities set — but the distribution is done
differently:

1. Rank 0 partitions the mesh (METIS / Morton / Hilbert).
2. The full mesh topology is broadcast using ``MPI.Win.Allocate_shared``:
   one physical copy per compute node, shared by all ranks on that node.
   Between nodes a single ``Bcast`` among node leaders transfers the data.
3. Each rank independently builds its own ghost layer using a Cython BFS
   kernel rather than receiving it from rank 0.
4. Quantities are distributed with ``MPI_Scatterv`` (each rank receives
   only its O(N/P) full-triangle rows) plus targeted ``MPI_Isend`` for
   ghost triangles.

The result is that topology is stored only once per node (not once per
rank), and quantity communication scales as O(N/P) rather than O(N).

.. code-block:: python

   import anuga
   from anuga.parallel.parallel_api import distribute_collaborative

   if anuga.myid == 0:
       domain = anuga.rectangular_cross_domain(1000, 1000, len1=10.0, len2=10.0)
       domain.set_quantity('elevation', lambda x, y: x / 10.0)
       domain.set_quantity('stage', expression='elevation + 0.2')
   else:
       domain = None

   # Shared-memory broadcast + per-rank BFS ghost layer
   domain = distribute_collaborative(domain)

   Br = anuga.Reflective_boundary(domain)
   domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

   for t in domain.evolve(yieldstep=1.0, finaltime=10.0):
       if anuga.myid == 0:
           domain.print_timestepping_statistics()

   domain.sww_merge()
   anuga.finalize()

Use ``distribute_collaborative`` when you have many MPI ranks sharing the
same node (e.g. 32–64 ranks on a large-memory node) or when the mesh is
large enough that topology duplication per rank becomes a concern.

.. note::

   ``distribute_collaborative`` still requires rank 0 to allocate the full
   quantity arrays before distribution.  If that alone exceeds available
   memory, use ``distribute_basic_mesh`` or
   ``distribute_basic_mesh_collaborative`` instead.


Mesh-first workflows
---------------------

For very large meshes where even allocating the full quantity arrays on
rank 0 is impractical, two *mesh-first* strategies are available:

* **Online** (:ref:`use_distribute_basic_mesh`) — rank 0 builds a
  ``Basic_mesh`` at runtime, distributes topology to all ranks, and each
  rank sets its own quantities.  Rank 0 still needs enough RAM for the
  full mesh topology.

* **Offline** (:ref:`use_sequential_mesh_io`) — partition the mesh in a
  separate preprocessing step (no MPI required, can use a large-memory
  workstation or login node), write one NetCDF4 file per rank, then each
  rank reads only its own file at simulation time.  Nothing beyond a
  per-rank partition is ever resident in memory.  This is also the
  recommended approach for ensemble / scenario runs on the same mesh.

The online approach via :func:`distribute_basic_mesh` is described below.

.. code-block:: python

   import anuga
   from anuga.abstract_2d_finite_volumes.basic_mesh import rectangular_cross_basic_mesh
   from anuga.parallel.parallel_api import distribute_basic_mesh

   if anuga.myid == 0:
       bm = rectangular_cross_basic_mesh(1000, 1000, len1=10.0, len2=10.0)
   else:
       bm = None

   # Only mesh topology is distributed; quantities are never on rank 0
   domain = distribute_basic_mesh(bm)

   domain.set_quantity('elevation', lambda x, y: x / 10.0)
   domain.set_quantity('stage', expression='elevation + 0.2')

   Br = anuga.Reflective_boundary(domain)
   domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

   for t in domain.evolve(yieldstep=1.0, finaltime=10.0):
       if anuga.myid == 0:
           domain.print_timestepping_statistics()

   domain.sww_merge()
   anuga.finalize()

See :ref:`use_distribute_basic_mesh` for full details, including the
shared-memory variant ``distribute_basic_mesh_collaborative``.


Running parallel scripts
-------------------------

Use ``mpiexec`` (or ``mpirun``) to launch the script.  For example, to
run on 8 MPI processes::

   mpiexec -np 8 python -u run_model.py

The ``-u`` flag gives unbuffered output so that ``print`` statements from
all ranks appear in real time.

Example scripts are in ``examples/parallel/``:

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - Script
     - Demonstrates
   * - ``run_parallel_rectangular.py``
     - ``distribute`` with a rectangular domain
   * - ``run_parallel_merimbula.py``
     - ``distribute`` with a real DEM and polygon mesh
   * - ``run_parallel_basic_mesh.py``
     - ``distribute_basic_mesh``
   * - ``run_parallel_basic_mesh_collaborative.py``
     - ``distribute_basic_mesh_collaborative``


Merging parallel SWW files
---------------------------

At the end of a parallel run each MPI rank writes its own per-rank SWW
file.  ``domain.sww_merge()`` (called from rank 0) reassembles them into
a single global SWW file.

By default the merge reads all timesteps from every per-rank file at
once, which is the fastest approach but requires enough RAM to hold the
full merged time series for each dynamic quantity:

.. code-block:: python

   domain.sww_merge(delete_old=True)   # default: all timesteps in RAM

For large simulations where the full time series does not fit in RAM, use
the ``chunk_size`` parameter to process only *N* timesteps at a time.
The merge will make multiple passes through the input files — more I/O,
but peak memory is bounded to approximately
``chunk_size × n_global_nodes × 4 bytes`` per dynamic quantity:

.. code-block:: python

   # Process 50 timesteps per pass — reduces peak RAM substantially
   domain.sww_merge(delete_old=True, chunk_size=50)

As a rough sizing guide:

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - Global nodes (smooth) / vertices (non-smooth)
     - chunk_size=50 peak RAM (MB)
     - chunk_size=200 peak RAM (MB)
   * - 100 000
     - ~20 MB
     - ~80 MB
   * - 1 000 000
     - ~200 MB
     - ~800 MB
   * - 10 000 000
     - ~2 GB
     - ~8 GB

*(figures are per dynamic quantity and assume float32 storage)*

You can also call ``sww_merge_parallel`` directly from post-processing
scripts without a running domain:

.. code-block:: python

   from anuga.utilities.sww_merge import sww_merge_parallel

   sww_merge_parallel('cairns', np=8, delete_old=True, chunk_size=100)


Choosing the number of processes
----------------------------------

There is a practical limit to how many MPI ranks benefit a simulation.
As a rough guide, aim for at least 1000–2000 triangles per rank after
partitioning.  Below that threshold, the overhead of ghost-cell
communication and MPI synchronisation outweighs the computational
saving.  Some experimentation for your specific mesh and hardware is
always worthwhile.


.. _load_balance_statistics:

Load balance monitoring
------------------------

METIS partitions the mesh so that each rank owns roughly the same number
of triangles.  However, in inundation simulations the computational cost
per triangle depends on whether it is *wet* or *dry*:

- **Wet triangle**: full flux computation, extrapolation, friction update.
- **Dry triangle**: flux computation is skipped (``optimise_dry_cells``);
  far cheaper.

As the inundation front advances, some ranks end up owning mostly wet
triangles while others own mostly dry coastal-plain triangles.  The dry
ranks finish their compute phase early and then sit idle in the
``MPI_Allreduce`` timestep-synchronisation barrier waiting for the wet
ranks.  This idle time appears in ``domain.communication_reduce_time`` and
is a direct measure of load imbalance caused by the wet/dry distribution.

``load_balance_statistics`` and ``print_load_balance_statistics`` gather
these per-rank numbers via ``MPI_Allgather`` and report them on rank 0.

.. code-block:: python

   for t in domain.evolve(yieldstep=60.0, finaltime=3600.0):
       if anuga.myid == 0:
           domain.print_timestepping_statistics()
       domain.print_load_balance_statistics()   # all ranks participate

Example output for a 2-rank inundation run where rank 0 owns the wet
half and rank 1 owns the dry coastal floodplain::

   Load balance statistics=========================================
   Rank    n_full   n_ghost    wet%  ghost%  compute(s)   comm(s)  Allreduce_wait(s)
   ------------------------------------------------------------------------
      0     19999       145    75.0     0.7       1.42     0.016              0.008
      1     20001       142     0.0     0.7       0.63     0.016              0.792
   ================================================================
     Imbalance ratio (max/mean compute): 1.38
     Pearson r(wet_fraction, compute_time): +1.000
     Interpretation: wetter ranks do more work (|r| = 1.000)

Key fields:

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Field
     - Meaning
   * - ``wet%``
     - Percentage of owned (full) triangles with depth >
       ``minimum_allowed_height``.
   * - ``ghost%``
     - Halo triangles as a fraction of all triangles.  A value above
       ~15% suggests the subdomain boundary is too large relative to
       the interior — consider a different partition or fewer ranks.
   * - ``compute(s)``
     - Wall time minus all communication time.  The limiting rank sets
       the pace of the whole simulation.
   * - ``Allreduce_wait(s)``
     - Time spent idle in the timestep-sync barrier.  Large values on
       a rank indicate it arrives early because it has little wet work
       to do.  This is wasted wall time.
   * - ``Imbalance ratio``
     - ``max(compute_time) / mean(compute_time)`` across ranks.  A
       value of 1.0 is perfect balance; 1.5 means the slowest rank
       does 50% more work than average and the others wait for it.
   * - ``Pearson r``
     - Correlation between wet fraction and compute time across ranks.
       Values near +1 confirm that wet/dry distribution is the dominant
       source of imbalance (as opposed to, for example, ghost-cell
       overhead or structure operator costs).

The method can also be called programmatically:

.. code-block:: python

   stats = domain.load_balance_statistics()
   # stats is a dict of numpy arrays (length numproc) on rank 0; None on others
   if anuga.myid == 0:
       print(f"Imbalance ratio: {stats['imbalance_ratio']:.2f}")
       print(f"Wet fractions:   {stats['wet_fraction']}")

**What to do about imbalance**

Static METIS decomposition cannot know in advance how the inundation
front will evolve, so some wet/dry imbalance is expected for most real
scenarios.  Practical mitigations:

1. **Use fewer ranks for mostly-dry domains.** If most of the domain is
   dry for most of the simulation, MPI overhead may not be worthwhile.

2. **Weight the METIS partition towards expected wet regions.** If you
   know in advance (from a coarse serial run or bathymetry) which regions
   will remain wet, you can provide per-triangle weights to the partitioner
   so that wet triangles are counted more heavily.

3. **Accept the imbalance and use the saved time differently.** If ranks
   are idle, consider adding more output or operator work on those ranks
   rather than spinning in the barrier.

4. **Dynamic repartitioning (future work).** Repartitioning the mesh
   mid-run as the wet front advances is the algorithmically correct
   solution but requires transferring all quantity arrays between ranks
   and is not yet implemented in ANUGA.

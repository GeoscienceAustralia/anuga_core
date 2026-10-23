.. _use_sequential_mesh_io:

.. currentmodule:: anuga

Offline mesh partitioning (``sequential_mesh_dump`` / ``sequential_mesh_load``)
================================================================================

Overview
--------

The ``sequential_mesh_dump`` / ``sequential_mesh_load`` pair implements an
*offline mesh partitioning* workflow:

1. **Preprocessing** — run once, on any machine (workstation, login node, …):

   .. code-block:: bash

      python create_dump.py -np N

   Creates one NetCDF4 file per rank containing mesh topology and the halo
   communication structure.  No quantity data is written.

2. **Simulation** — run as many times as needed, on the HPC cluster:

   .. code-block:: bash

      mpiexec -np N python -u run_evolve.py

   Each rank reads its own file independently (no rank-0 bottleneck), sets
   its own initial conditions from a function or DEM, and proceeds directly
   to ``domain.evolve()``.

This differs from the ``sdpl`` (``sequential_distribute_dump`` /
``sequential_distribute_load``) workflow, which stores the full domain
including quantities.

When to use this approach
--------------------------

.. list-table::
   :header-rows: 1
   :widths: 55 45

   * - Situation
     - Recommendation
   * - Full domain + quantities fits in rank-0 RAM
     - :func:`distribute` or :func:`sequential_distribute_dump`
   * - Mesh fits in rank-0 RAM but quantity arrays do not
     - :func:`distribute_basic_mesh` (live) or this workflow (offline)
   * - Mesh is too large for rank-0 RAM at runtime but a preprocessing
       node has sufficient memory
     - **This workflow**
   * - Many scenario / ensemble variants on the same mesh
     - **This workflow** — dump once, load-and-evolve many times
   * - Initial conditions come from a function or per-rank DEM read
     - **This workflow**

Memory comparison
^^^^^^^^^^^^^^^^^

For a mesh with N triangles and P quantities:

.. list-table::
   :header-rows: 1
   :widths: 40 30 30

   * - Workflow
     - Rank-0 peak RAM
     - Required at runtime
   * - ``distribute``
     - O(N) × (mesh + P qty)
     - MPI job running
   * - ``distribute_basic_mesh``
     - O(N) × mesh only
     - MPI job running
   * - ``sequential_mesh_dump`` / ``load``
     - O(N) × mesh only (preprocessing)
     - Only per-rank partition


API
---

.. autofunction:: sequential_mesh_dump

.. autofunction:: sequential_mesh_load


File format
-----------

Each partition file ``<name>_mesh_P<N>_<rank>.nc`` is a self-contained
NETCDF4 file.  It can be inspected with ``ncdump -h``.

**Global attributes** (scalar metadata)

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Attribute
     - Description
   * - ``rank``, ``numprocs``
     - This partition's rank and total partition count.
   * - ``number_of_full_triangles``
     - Triangles owned by this rank (excludes ghost layer).
   * - ``number_of_full_nodes``
     - Nodes owned by this rank (excludes ghost layer nodes).
   * - ``number_of_global_triangles``
     - Total triangles in the global mesh.
   * - ``number_of_global_nodes``
     - Total nodes in the global mesh.
   * - ``ghost_layer_width``
     - Depth of the halo ghost layer (typically 2).
   * - ``xllcorner``, ``yllcorner``, ``zone``, …
     - :class:`Geo_reference` coordinate origin (same attributes as SWW files).

**Variables**

.. list-table::
   :header-rows: 1
   :widths: 22 20 58

   * - Variable
     - Shape
     - Description
   * - ``points``
     - (node, 2)
     - Node (x, y) coordinates.  Includes ghost nodes.
   * - ``vertices``
     - (tri, 3)
     - Triangle connectivity (local node indices).  Includes ghost triangles.
   * - ``tri_l2g``
     - (tri,)
     - Local-to-global triangle index map.
   * - ``node_l2g``
     - (node,)
     - Local-to-global node index map.
   * - ``boundary_tri``, ``boundary_edge``, ``boundary_tag``
     - (bnd,)
     - Boundary edge records: triangle index, edge index (0–2), and string tag.
   * - ``send_ranks``, ``send_offsets``, ``send_local``, ``send_global``
     - CSR arrays
     - Full-send communication pattern (which local triangles to send to which
       neighbour ranks, with their global IDs).
   * - ``recv_ranks``, ``recv_offsets``, ``recv_local``, ``recv_global``
     - CSR arrays
     - Ghost-receive communication pattern.


Performance for large partition counts
--------------------------------------

For very large meshes split into many partitions, **writing the partition
files dominates** the preprocessing time.  The mesh dump always writes one
self-describing NetCDF file per rank (topology only, no quantities), so — unlike
the domain dump, which also has a single-file layout knob — the lever here is
parallel writing:

``num_workers`` (default ``1``)
   With ``num_workers > 1`` (on a POSIX/fork platform) the per-rank NetCDF files
   are written in parallel by a pool of worker processes that share the
   partitioned mesh copy-on-write, so the write scales toward the filesystem's
   I/O and metadata throughput.  Match ``num_workers`` to the machine doing the
   preprocessing (often a large-memory login/preprocessing node).

.. code-block:: python

   anuga.sequential_mesh_dump(domain, numprocs=18400,
                              partition_dir='Partitions', num_workers=32)

.. note::

   The serial default (``num_workers=1``) releases each rank's memory as it is
   written, keeping rank-0 peak RAM low.  The parallel path keeps the whole
   partitioned mesh resident for the duration of the worker pool — the
   memory-for-speed trade-off — so choose ``num_workers`` with the preprocessing
   node's RAM in mind.


Preprocessing example
---------------------

.. code-block:: python

   # create_partitions.py  — run once on a workstation or login node
   # python create_partitions.py -np 64

   import argparse
   import anuga
   from anuga import rectangular_cross_domain          # or your mesh builder

   parser = argparse.ArgumentParser()
   parser.add_argument('-np', '--numprocs', type=int, default=8)
   args = parser.parse_args()

   domain = rectangular_cross_domain(500, 500, len1=10.0, len2=10.0)
   domain.set_name('flood_mesh')

   anuga.sequential_mesh_dump(
       domain,
       numprocs=args.numprocs,
       partition_dir='Partitions',
       parameters={'partition_scheme': 'metis', 'ghost_layer_width': 2},
       verbose=True,
   )
   # Writes: Partitions/flood_mesh_mesh_P<N>_<rank>.nc

When building a mesh from a DEM or polygon file, replace
``rectangular_cross_domain`` with your usual
:func:`create_domain_from_regions` call.  The ``set_quantity`` calls for
initial conditions are *not* needed here — each rank will set them
independently at runtime.


Parallel load-and-evolve example
---------------------------------

.. code-block:: python

   # run_evolve.py  — run with:  mpiexec -np N python -u run_evolve.py

   import anuga
   from anuga import myid, numprocs, finalize, barrier, Reflective_boundary

   barrier()
   domain = anuga.sequential_mesh_load(name='flood_mesh',
                                       partition_dir='Partitions',
                                       verbose=(myid == 0))
   barrier()

   # --- quantities: each rank evaluates independently ---
   domain.set_quantity('elevation', lambda x, y: 0.1 * x)
   domain.set_quantity('stage',     expression='elevation + 0.5')
   domain.set_quantity('friction',  0.03)

   # --- boundary conditions ---
   Br = Reflective_boundary(domain)
   domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

   # --- evolve ---
   domain.set_name('flood_mesh')
   domain.set_flow_algorithm('DE0')

   for t in domain.evolve(yieldstep=60.0, finaltime=3600.0):
       if myid == 0:
           domain.print_timestepping_statistics()

   domain.sww_merge(delete_old=True)
   finalize()


Using a DEM for elevation
^^^^^^^^^^^^^^^^^^^^^^^^^

Because each rank calls ``set_quantity`` with its own local coordinates,
you can read a DEM file directly without any rank-0 bottleneck:

.. code-block:: python

   from anuga import Quantity
   import anuga

   # file_function reads only the region covered by this rank's triangles
   elev_func = anuga.file_function('topography.asc', domain,
                                    quantities=['elevation'])
   domain.set_quantity('elevation', elev_func)

Alternatively, pass a callable that reads from the DEM independently on
each rank — the call happens only on the local mesh centroids.


Scenario ensembles on the same mesh
-------------------------------------

A common flood risk workflow: one mesh, many hydrological scenarios.

.. code-block:: bash

   # Step 1 — partition once
   python create_partitions.py --mesh flood_mesh --np 32

   # Step 2 — run each scenario in parallel (same partition files)
   for SCENARIO in Q10 Q100 Q500 QPMF; do
       mpiexec -np 32 python run_evolve.py --scenario $SCENARIO
   done

The partition files are reused across all scenarios.  Only the
``set_quantity`` calls (initial water level, rainfall rate, etc.) differ
between runs.


Example scripts
---------------

Ready-to-run examples are in ``examples/parallel/``:

.. list-table::
   :header-rows: 1
   :widths: 52 48

   * - Script
     - Description
   * - ``run_smpl_rectangular_create_dump.py``
     - Creates a rectangular-cross mesh and writes partition files.
       Command line: ``python run_smpl_rectangular_create_dump.py -np N -sn 100``
   * - ``run_smpl_rectangular_load_evolve.py``
     - Loads partition files and runs the evolve loop.
       Command line: ``mpiexec -np N python -u run_smpl_rectangular_load_evolve.py -sn 100``

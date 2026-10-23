.. _use_distribute_basic_mesh:

.. currentmodule:: anuga

MPI Distribute Mesh-first parallel workflow (``distribute_basic_mesh``)
========================================================================

Overview
--------

The ``distribute_basic_mesh`` workflow is the recommended way to start a
parallel ANUGA simulation when your mesh can be described as a structured grid
(rectangular or rectangular-cross) or when you have a ``Pmesh`` object from
``create_pmesh_from_regions``.

The key idea is that **only rank 0 builds the mesh**.  All other ranks pass
``None`` to ``distribute_basic_mesh`` and receive their submesh automatically
via MPI.  Initial conditions (quantities, boundaries, operators) are then set
on every rank's local ``Parallel_domain`` after distribution.

This contrasts with the traditional ``distribute()`` workflow, where rank 0
builds a full ``Domain`` *including quantities*, which means all the elevation,
stage, and momentum arrays are allocated on rank 0 before any MPI communication
takes place.  For large meshes this can be a significant memory and time saving.


Basic_mesh
----------

``Basic_mesh`` is a lightweight mesh container that stores only the topology
needed for partitioning:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Attribute
     - Description
   * - ``nodes``
     - ``(M, 2)`` float array of node x, y coordinates
   * - ``triangles``
     - ``(N, 3)`` int array of triangle vertex indices
   * - ``boundary``
     - dict mapping ``(triangle_id, edge_index)`` → tag string
   * - ``geo_reference``
     - coordinate origin and UTM zone
   * - ``number_of_triangles``
     - ``N``
   * - ``number_of_nodes``
     - ``M``
   * - ``neighbours``
     - ``(N, 3)`` int array, -1 for boundary edges (computed lazily)
   * - ``centroid_coordinates``
     - ``(N, 2)`` float array (computed lazily)

Unlike the full ``Mesh`` class, ``Basic_mesh`` does **not** compute normals,
edge lengths, areas, radii, or vertex coordinates at construction time.  Those
are only needed by the finite-volume solver and are computed later on each rank
when the ``Parallel_domain`` is built.


Factory functions
-----------------

Structured rectangular grids
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from anuga.abstract_2d_finite_volumes.basic_mesh import (
       rectangular_basic_mesh,
       rectangular_cross_basic_mesh,
   )

   # 2*m*n triangles, (m+1)*(n+1) nodes
   bm = rectangular_basic_mesh(m, n, len1=Lx, len2=Ly)

   # 4*m*n triangles, (m+1)*(n+1) + m*n nodes  (recommended for accuracy)
   bm = rectangular_cross_basic_mesh(m, n, len1=Lx, len2=Ly)

Both functions accept an optional ``origin=(x0, y0)`` keyword argument.

Polygon regions
~~~~~~~~~~~~~~~

For meshes created with ``create_pmesh_from_regions``, use
``pmesh_to_basic_mesh`` to obtain a ``Basic_mesh`` without building a full
``Domain``:

.. code-block:: python

   import anuga
   from anuga.abstract_2d_finite_volumes.pmesh2domain import pmesh_to_basic_mesh

   pmesh = anuga.create_pmesh_from_regions(
       bounding_polygon,
       boundary_tags={'left': [0], 'right': [1], 'top': [2], 'bottom': [3]},
       maximum_triangle_area=0.01,
       filename=None,
       use_cache=False,
   )
   bm = pmesh_to_basic_mesh(pmesh)


``distribute_basic_mesh``
-------------------------

.. code-block:: python

   from anuga.parallel.parallel_api import distribute_basic_mesh

   domain = distribute_basic_mesh(basic_mesh, verbose=False, parameters=None)

Partitions the mesh on rank 0 and sends each rank its local submesh.

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Argument
     - Description
   * - ``basic_mesh``
     - A ``Basic_mesh`` on rank 0; ``None`` on all other ranks
   * - ``verbose``
     - Print partitioning progress (default ``False``)
   * - ``parameters``
     - Dict of partitioning options, e.g. ``{'partition_scheme': 'morton'}``

Available partition schemes are ``'metis'`` (default, best quality),
``'morton'``, and ``'hilbert'``.

Returns a ``Parallel_domain`` on every rank, or a plain ``Domain`` when
``numprocs == 1`` so the same script runs serially without modification.

After the call, every rank holds only its local submesh — full and ghost
triangles for that rank's partition.  The caller then sets quantities,
boundaries, operators, and structures on the returned domain exactly as in a
serial script.


General script structure
------------------------

.. code-block:: python

   """Example parallel simulation using distribute_basic_mesh."""

   import anuga
   from anuga import (
       Reflective_boundary,
       Dirichlet_boundary,
       myid, numprocs, barrier, finalize,
   )
   from anuga.abstract_2d_finite_volumes.basic_mesh import rectangular_cross_basic_mesh
   from anuga.parallel.parallel_api import distribute_basic_mesh

   # 1. Build the mesh on rank 0 only
   if myid == 0:
       bm = rectangular_cross_basic_mesh(100, 100, len1=100.0, len2=100.0)
   else:
       bm = None

   # 2. Distribute the mesh to all ranks
   domain = distribute_basic_mesh(bm, parameters={'partition_scheme': 'metis'})

   # Give the simulation a name.
   # Each rank automatically appends _P<nproc>_<rank> to the filename.
   domain.set_name('my_simulation')
   domain.set_datadir('output')

   # 3. Set initial conditions on the local submesh.
   #    Each rank evaluates quantities only over its own triangles.
   def topography(x, y):
       return -x / 100.0

   domain.set_quantity('elevation', topography)
   domain.set_quantity('friction',  0.03)
   domain.set_quantity('stage',     expression='elevation')

   # 4. Set boundary conditions.
   #    ghost and exterior edges are handled automatically by Parallel_domain.
   Br = Reflective_boundary(domain)
   Bd = Dirichlet_boundary([-0.2, 0.0, 0.0])
   domain.set_boundary({'left': Br, 'right': Bd, 'top': Br, 'bottom': Br})

   # 5. Evolve
   for t in domain.evolve(yieldstep=60.0, finaltime=3600.0):
       if myid == 0:
           domain.print_timestepping_statistics()

   # 6. Merge per-rank SWW files into a single file (optional)
   domain.sww_merge(delete_old=False)

   finalize()

Run with::

   mpirun -np 4 python my_simulation.py


Notes
-----

Only rank 0 needs the ``Basic_mesh``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   if myid == 0:
       bm = rectangular_cross_basic_mesh(M, N, len1=Lx, len2=Ly)
   else:
       bm = None                         # other ranks pass None

   domain = distribute_basic_mesh(bm)   # all ranks call this

Allocating the full mesh on every rank defeats the purpose of the workflow and
increases peak memory usage.

Quantities are set after distribution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Unlike ``distribute()``, which copies quantity arrays from rank 0, with
``distribute_basic_mesh`` each rank sets its own quantities independently after
receiving its submesh:

.. code-block:: python

   domain.set_quantity('elevation', topography)  # evaluated per-rank
   domain.set_quantity('stage',     0.0)

This is more memory-efficient and often faster for large initial-condition
functions, since rank 0 never holds the full quantity arrays.

Boundary tags
~~~~~~~~~~~~~

The boundary tags in the ``Basic_mesh`` boundary dict must match the strings
used in ``domain.set_boundary()``.  For rectangular and rectangular-cross
meshes the built-in tags are ``'left'``, ``'right'``, ``'top'``, and
``'bottom'``.

``Parallel_domain.set_boundary`` automatically handles the internal ``'ghost'``
and ``'exterior'`` tags that appear on inter-rank edges — you do not need to
include them in your boundary map.

Serial compatibility
~~~~~~~~~~~~~~~~~~~~~

When ``numprocs == 1``, ``distribute_basic_mesh`` returns a plain ``Domain``
built directly from the ``Basic_mesh`` with no MPI calls.  The same script
therefore runs correctly both serially and in parallel without any
modifications.

SWW output
~~~~~~~~~~

Each rank writes its own SWW file named ``<name>_P<nproc>_<rank>.sww``.  Call
``domain.sww_merge()`` after the evolve loop to combine them into a single
``<name>.sww`` for visualisation.  To suppress SWW output entirely:

.. code-block:: python

   domain.set_quantities_to_be_stored(None)


``distribute_basic_mesh_collaborative``
---------------------------------------

``distribute_basic_mesh_collaborative`` is an alternative to
``distribute_basic_mesh`` that uses shared memory within each compute node
to reduce inter-process data movement.  The API is identical:

.. code-block:: python

   from anuga.parallel.parallel_api import distribute_basic_mesh_collaborative

   domain = distribute_basic_mesh_collaborative(basic_mesh, verbose=False, parameters=None)

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Argument
     - Description
   * - ``basic_mesh``
     - A ``Basic_mesh`` on rank 0; ``None`` on all other ranks
   * - ``verbose``
     - Print partitioning progress (default ``False``)
   * - ``parameters``
     - Dict of partitioning options, e.g. ``{'partition_scheme': 'metis'}``

Returns a ``Parallel_domain`` on every rank, or a plain ``Domain`` when
``numprocs == 1``.


How it differs from ``distribute_basic_mesh``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Both functions share the same mesh-first interface — rank 0 builds a
``Basic_mesh``, all other ranks pass ``None``, and quantities are set
per-rank after distribution.  The difference is in how the mesh topology
reaches each rank:

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - Step
     - ``distribute_basic_mesh``
     - ``distribute_basic_mesh_collaborative``
   * - Partition
     - rank 0 (METIS / Morton / Hilbert)
     - rank 0 (same)
   * - Topology broadcast
     - Point-to-point ``send`` / ``recv`` per rank
     - Shared-memory window on each node; one physical copy per node
   * - Ghost layer
     - Extracted from global submesh data sent by rank 0
     - Per-rank BFS on the local shared-memory copy
   * - Quantities
     - None distributed (set per-rank after return)
     - None distributed (set per-rank after return)

The collaborative variant is advantageous when many MPI ranks share the
same node (e.g. 32–64 ranks on a single large-memory node) because the
full topology is stored only once in shared memory per node rather than
once per rank.


Example script
~~~~~~~~~~~~~~

The script structure is identical to ``distribute_basic_mesh``; simply
swap the function name:

.. code-block:: python

   """Example parallel simulation using distribute_basic_mesh_collaborative."""

   import anuga
   from anuga import (
       Reflective_boundary,
       Dirichlet_boundary,
       myid, numprocs, barrier, finalize,
   )
   from anuga.abstract_2d_finite_volumes.basic_mesh import rectangular_cross_basic_mesh
   from anuga.parallel.parallel_api import distribute_basic_mesh_collaborative

   # 1. Build the mesh on rank 0 only
   if myid == 0:
       bm = rectangular_cross_basic_mesh(100, 100, len1=100.0, len2=100.0)
   else:
       bm = None

   # 2. Distribute the mesh to all ranks (shared-memory broadcast)
   domain = distribute_basic_mesh_collaborative(bm, parameters={'partition_scheme': 'metis'})

   domain.set_name('my_simulation')
   domain.set_datadir('output')

   # 3. Set initial conditions per-rank
   def topography(x, y):
       return -x / 100.0

   domain.set_quantity('elevation', topography)
   domain.set_quantity('friction',  0.03)
   domain.set_quantity('stage',     expression='elevation')

   # 4. Set boundary conditions
   Br = Reflective_boundary(domain)
   Bd = Dirichlet_boundary([-0.2, 0.0, 0.0])
   domain.set_boundary({'left': Br, 'right': Bd, 'top': Br, 'bottom': Br})

   # 5. Evolve
   for t in domain.evolve(yieldstep=60.0, finaltime=3600.0):
       if myid == 0:
           domain.print_timestepping_statistics()

   # 6. Merge per-rank SWW files (optional)
   domain.sww_merge(delete_old=False)

   finalize()

Run with::

   mpirun -np 4 python my_simulation.py


When to prefer ``distribute_basic_mesh_collaborative``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- **Many ranks per node**: shared memory eliminates redundant topology copies.
- **Large meshes**: peak RSS per node is lower because each node holds one
  copy of the full topology instead of one copy per rank.
- **Homogeneous nodes**: all ranks on a node must share memory; mixed
  hardware configurations may require additional tuning.

For small meshes or a small number of ranks the two functions perform
equivalently; prefer ``distribute_basic_mesh`` if ``mpi4py`` shared-memory
windows are not available in your environment.

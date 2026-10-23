.. _use_mesh_refinement:

.. currentmodule:: anuga

Uniform mesh refinement
========================

ANUGA provides three functions for uniformly refining a mesh — replacing
every triangle with four children by connecting edge midpoints (red
refinement):

.. code-block:: text

    parent triangle (v0, v1, v2)
    ─────────────────────────────
        v2
        /\
       /  \
    m02────m12
     / \  / \
    /   \/   \
   v0──m01───v1

    child 0 = (v0,  m01, m02)   corner near v0
    child 1 = (m01, v1,  m12)   corner near v1
    child 2 = (m02, m12, v2 )   corner near v2
    child 3 = (m01, m12, m02)   central


.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Function
     - Purpose
   * - :func:`uniform_refine_domain`
     - Refine a sequential :class:`Domain` in memory; return refined Domain.
   * - :func:`create_parallel_mesh`
     - Partition a coarse mesh and optionally refine it offline; write one
       NetCDF4 file per rank ready for :func:`sequential_mesh_load`.
   * - :func:`sequential_mesh_refine`
     - Refine an existing set of partition files produced by
       :func:`sequential_mesh_dump` or a previous refinement.


Sequential refinement
---------------------

:func:`uniform_refine_domain` takes any sequential :class:`Domain` (or one
built from a :class:`Basic_mesh`) and returns a new quantity-free domain with
4× as many triangles:

.. code-block:: python

    import anuga

    coarse = anuga.rectangular_cross_domain(50, 50, len1=10.0, len2=10.0)

    fine = anuga.uniform_refine_domain(coarse)
    print(fine.number_of_triangles)   # 4 × coarse count

    # Quantities must be set on the refined domain before evolving.
    fine.set_quantity('elevation', lambda x, y: 0.1 * x)
    fine.set_quantity('stage', expression='elevation + 0.5')
    # … boundary conditions, operators, evolve …

Multiple levels:

.. code-block:: python

    r1 = anuga.uniform_refine_domain(coarse)           # 4×
    r2 = anuga.uniform_refine_domain(r1)               # 16×


Parallel workflow: coarse partition + offline refinement
---------------------------------------------------------

The recommended parallel workflow when the *refined* mesh would be very large:

1. **Partition the coarse mesh** on a workstation or login node — this
   requires only coarse-mesh RAM.

2. **Refine each partition** offline (also on the workstation / login node)
   — the refinement touches one partition at a time and needs only
   per-rank RAM.

3. **Load and evolve** in parallel on the HPC cluster — each rank reads its
   refined partition file independently with :func:`sequential_mesh_load`.

:func:`create_parallel_mesh` performs steps 1 and 2 in a single call:

.. code-block:: python

    # create_refined_mesh.py  — run once on any machine
    import anuga
    from anuga import rectangular_cross_domain

    coarse = rectangular_cross_domain(100, 100, len1=10.0, len2=10.0)
    coarse.set_name('flood_mesh')

    # Partition into 32 ranks and refine twice (4² = 16× more triangles).
    # num_workers=4 refines 4 partition files concurrently (Pass 2 only).
    anuga.create_parallel_mesh(
        coarse,
        numprocs=32,
        refinement_levels=2,
        partition_dir='Partitions',
        verbose=True,
        num_workers=4,
    )
    # Writes: Partitions/flood_mesh_mesh_P32_<rank>.nc  (32 files)

.. note::

   ``num_workers`` parallelises **Pass 2** (the per-partition refinement),
   which is embarrassingly parallel.  Pass 1 (global edge collection) is
   always single-process and fast.  Set ``num_workers`` to the number of
   CPU cores available on the preprocessing machine; the wall time for the
   refinement step scales roughly as ``numprocs / num_workers``.

The load-and-evolve step is identical to the non-refined workflow:

.. code-block:: python

    # run_evolve.py  — mpiexec -np 32 python -u run_evolve.py

    import anuga
    from anuga import myid, numprocs, finalize, barrier, Reflective_boundary

    barrier()
    domain = anuga.sequential_mesh_load(name='flood_mesh',
                                        partition_dir='Partitions',
                                        verbose=(myid == 0))
    barrier()

    domain.set_quantity('elevation', lambda x, y: 0.1 * x)
    domain.set_quantity('stage', expression='elevation + 0.5')
    domain.set_quantity('friction', 0.03)

    Br = Reflective_boundary(domain)
    domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

    for t in domain.evolve(yieldstep=60.0, finaltime=3600.0):
        if myid == 0:
            domain.print_timestepping_statistics()

    domain.sww_merge(delete_old=True)
    finalize()


Accepting a ``Basic_mesh``
^^^^^^^^^^^^^^^^^^^^^^^^^^

:func:`create_parallel_mesh` also accepts a :class:`Basic_mesh` directly,
which avoids building a full solver domain just for partitioning:

.. code-block:: python

    from anuga import Basic_mesh

    nodes, triangles, boundary = my_mesh_builder(...)
    mesh = Basic_mesh(nodes, triangles, boundary)

    anuga.create_parallel_mesh(
        mesh,
        numprocs=64,
        refinement_levels=1,
        name='mymesh',
        partition_dir='Partitions',
    )


Refining existing partition files
-----------------------------------

If partition files already exist (written by :func:`sequential_mesh_dump`
or a previous call to :func:`sequential_mesh_refine`),
:func:`sequential_mesh_refine` refines them in place without re-partitioning:

.. code-block:: python

    # Partition the coarse mesh first
    anuga.sequential_mesh_dump(coarse, numprocs=32, partition_dir='Coarse')

    # Refine one level — reads from 'Coarse/', writes to 'Fine/'
    # Use num_workers to refine partitions concurrently.
    anuga.sequential_mesh_refine(
        name='flood_mesh',
        numprocs=32,
        levels=1,
        partition_dir='Coarse',
        output_dir='Fine',
        num_workers=4,
    )

    # Or refine further in a second pass
    anuga.sequential_mesh_refine(
        name='flood_mesh',
        numprocs=32,
        levels=1,
        partition_dir='Fine',
        output_dir='Finest',
        num_workers=4,
    )

This separation is useful when you want to keep coarse partition files for
quick exploratory runs and fine partition files for production runs.


Choosing the number of refinement levels
-----------------------------------------

Each refinement level multiplies the triangle count by 4 and halves the
characteristic edge length, so the CFL timestep also roughly halves.

.. list-table::
   :header-rows: 1
   :widths: 15 18 20 47

   * - Levels
     - Triangle factor
     - Edge-length factor
     - Typical use
   * - 0
     - 1×
     - 1×
     - Production run at coarse resolution, or when METIS gives a good
       balanced partition at target resolution directly.
   * - 1
     - 4×
     - 1/2
     - Single refinement of a well-partitioned coarse mesh; convenient
       when the target mesh is 2–4× too large for direct partitioning RAM.
   * - 2
     - 16×
     - 1/4
     - Common choice: partition a mesh of ~10 K triangles → run at ~160 K.
   * - 3
     - 64×
     - 1/8
     - Large production meshes; partition at ~10 K, run at ~640 K.

.. note::

   The METIS partition quality is determined by the *coarse* mesh topology.
   Very coarse meshes may produce slightly unbalanced partitions that
   persist after refinement.  If load imbalance is a concern, increase the
   coarse mesh resolution or use the ``'hilbert'`` scheme:

   .. code-block:: python

       anuga.create_parallel_mesh(
           coarse, numprocs=64, refinement_levels=2,
           parameters={'partition_scheme': 'hilbert'},
       )


API reference
-------------

.. autofunction:: uniform_refine_domain

.. autofunction:: create_parallel_mesh

.. autofunction:: sequential_mesh_refine


Example scripts
---------------

Ready-to-run examples are in ``examples/parallel/``:

.. list-table::
   :header-rows: 1
   :widths: 52 48

   * - Script
     - Description
   * - ``run_refine_create_mesh.py``
     - Builds a rectangular-cross mesh, calls :func:`create_parallel_mesh`
       with configurable refinement levels, and writes the partition files.
       Command: ``python run_refine_create_mesh.py -np N -rl 2``
   * - ``run_smpl_rectangular_load_evolve.py``
     - Loads partition files (refined or not) and runs the evolve loop.
       Command: ``mpiexec -np N python -u run_smpl_rectangular_load_evolve.py -sn <coarse_sqrt_n>``

.. seealso::

   :doc:`use_sequential_mesh_io`
      The underlying dump / load functions that :func:`create_parallel_mesh`
      builds on.

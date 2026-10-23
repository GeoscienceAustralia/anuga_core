.. currentmodule:: anuga

.. _riverwalls:

Riverwalls
==========

A *riverwall* is an infinitely-thin wall placed along mesh edges.  It acts
like a weir: the flow over the wall is governed by the Villemonte overtopping
formula, blended smoothly with the shallow-water Riemann solution as
submergence increases.  Unlike raised bed elevations, a riverwall can be
higher than the surrounding bed without creating mesh artefacts, and its
crest height can be changed *during* a simulation.

Riverwalls are only available for discontinuous-elevation flow algorithms
(``'DE0'`` ``'DE1'`` ``'DE2'`` ``'DE_ader2'``).

.. contents:: Contents
   :local:
   :depth: 2


How they work
-------------

ANUGA stores one crest elevation value per riverwall *edge*.  During each
flux computation the C solver checks whether an edge belongs to a riverwall
(``edge_flux_type == 1``), and if so it uses the blended weir/shallow-water
discharge formula instead of the plain Riemann flux.

The blending is controlled by two ratios:

* **s** = tailwater head / headwater head  (submergence ratio)
* **h** = tailwater head / weir height     (absolute submergence ratio)

and four parameters *s1, s2, h1, h2*:

.. math::

   w_1 &= \text{clip}\!\left(\frac{s - s_1}{s_2 - s_1},\; 0,\; 1\right) \\
   w_2 &= \text{clip}\!\left(\frac{h - h_1}{h_2 - h_1},\; 0,\; 1\right) \\
   Q   &= \bigl(w_1\,Q_\text{SW} + (1 - w_1)\,Q_\text{ID}\bigr)(1 - w_2)
           + w_2\,Q_\text{SW}

where :math:`Q_\text{ID}` is the ideal Villemonte weir flow and
:math:`Q_\text{SW}` is the shallow-water Riemann flux.  When *s < s1* and
*h < h1* the ideal weir formula dominates; when *s > s2* or *h > h2* the
flow is fully determined by the shallow-water solution.

Optionally, water can also seep *through* the wall body below the crest via
a submerged orifice formula controlled by ``Cd_through`` (default 0 —
fully impermeable):

.. math::

   Q_\text{through} = C_{d,\text{through}}\; h_\text{eff}
                      \sqrt{2 g\,|\Delta\text{stage}|}

where :math:`h_\text{eff}` is the upstream submerged depth (depth below the
crest on the high-stage side).


Quick-start
-----------

A minimal two-step workflow:

1. Pass the wall polylines as ``breaklines`` when building the mesh — this
   forces triangle edges to align exactly with the wall.
2. Call :meth:`Domain.create_riverwalls` (or
   ``domain.riverwallData.create_riverwalls``) after setting initial
   conditions.

.. code-block:: python

    import anuga
    import numpy as np

    bounding_polygon = [[0, 0], [60, 0], [60, 20], [0, 20]]
    boundary_tags    = {'bottom': [0], 'right': [1], 'top': [2], 'left': [3]}

    # Each wall is a list of [x, y, z] points.
    # x, y are plan coordinates; z is the crest elevation (m).
    river_walls = {
        'leveeA': [[20.0,  0.0, 1.0],
                   [20.0, 20.0, 1.0]],
        'leveeB': [[40.0,  0.0, 0.8],
                   [40.0, 20.0, 0.8]],
    }

    # Step 1 — domain/mesh with breaklines aligned to the walls
    anuga.create_domain_from_regions(
        bounding_polygon,
        boundary_tags=boundary_tags,
        maximum_triangle_area=4.0,
        breaklines=list(river_walls.values()),
    )

    domain.set_quantity('elevation', 0.0, location='centroids')
    domain.set_quantity('friction',  0.03, location='centroids')
    domain.set_quantity('stage',     lambda x, y: np.where(x < 20, 1.5, 0.0),
                        location='centroids')

    Br = anuga.Reflective_boundary(domain)
    Bo = anuga.Dirichlet_boundary([0.0, 0.0, 0.0])
    domain.set_boundary({'left': Br, 'top': Br, 'bottom': Br, 'right': Bo})

    # Step 2 — create riverwalls (call AFTER set_boundary)
    domain.create_riverwalls(river_walls, verbose=False)

    for t in domain.evolve(yieldstep=5.0, finaltime=60.0):
        domain.print_timestepping_statistics()


Hydraulic parameters
--------------------

Pass a ``riverwallPar`` dictionary to customise parameters per wall.  Any
parameter not specified uses the default value shown below.

.. list-table::
   :header-rows: 1
   :widths: 15 12 73

   * - Parameter
     - Default
     - Description
   * - ``Qfactor``
     - 1.0
     - Multiplicative calibration factor applied to the ideal weir discharge.
       Values < 1 reduce flow; values > 1 increase it.
   * - ``s1``
     - 0.9
     - Submergence ratio *s = tailwater-head / headwater-head* below which
       the ideal weir formula is used without blending.  Must be < ``s2``.
   * - ``s2``
     - 0.95
     - Submergence ratio above which the shallow-water solution is used
       exclusively.  Must be > ``s1``.
   * - ``h1``
     - 1.0
     - Absolute submergence ratio *h = tailwater-head / weir-height* below
       which the ideal weir formula is used without blending.  Must be
       < ``h2``.
   * - ``h2``
     - 1.5
     - Absolute submergence ratio above which the shallow-water solution is
       used exclusively.  Must be > ``h1``.
   * - ``Cd_through``
     - 0.0
     - Discharge coefficient for flow through the wall body below the crest
       (submerged orifice formula).  Set to 0 (default) for an impermeable
       wall.  Typical values for a culvert-like opening: 0.5–0.8.

.. code-block:: python

    riverwall_par = {
        'leveeA': {'Qfactor': 1.0, 's1': 0.8, 's2': 0.95},
        'leveeB': {'Qfactor': 0.9, 'Cd_through': 0.3},
    }

    domain.create_riverwalls(river_walls, riverwallPar=riverwall_par)


Runtime interface
-----------------

The ``RiverWall`` object (``domain.riverwallData``) exposes methods to
inspect and modify wall state *inside* the evolve loop.  This allows
time-varying crest heights, tide-gate-like logic, and scenario sensitivity
studies without rebuilding the mesh.

.. code-block:: python

    rw = domain.riverwallData

    # List all wall names
    print(rw.get_wall_names())           # ['leveeA', 'leveeB']

    # Read edge coordinates (absolute xy, n × 2 array)
    xy = rw.get_edge_coordinates('leveeA')

    # Read current crest heights (one value per riverwall edge)
    elev = rw.get_elevation('leveeA')    # numpy array, copy

    # Set uniform crest height
    rw.set_elevation('leveeA', 0.5)

    # Set per-edge heights (array must match number of edges for this wall)
    n = len(rw.get_elevation('leveeA'))
    rw.set_elevation('leveeA', numpy.linspace(0.5, 1.0, n))

    # Add an offset to the current height
    rw.set_elevation_offset('leveeB', +0.3)   # raise by 0.3 m

    # Read / write a hydraulic parameter
    q = rw.get_hydraulic_parameter('leveeA', 'Qfactor')
    rw.set_hydraulic_parameter('leveeA', 'Qfactor', 1.5)

Use these calls inside ``domain.evolve`` to create dynamic wall behaviour:

.. code-block:: python

    for t in domain.evolve(yieldstep=5.0, finaltime=120.0):
        domain.print_timestepping_statistics()

        if abs(t - 30.0) < 1e-6:
            # Lower leveeA crest at t = 30 s to release impounded water
            rw.set_elevation('leveeA', 0.3)

        if abs(t - 60.0) < 1e-6:
            # Raise leveeB to prevent downstream flooding
            rw.set_elevation_offset('leveeB', 1.0)

        if abs(t - 90.0) < 1e-6:
            # Reduce Qfactor on leveeA to calibrate peak discharge
            rw.set_hydraulic_parameter('leveeA', 'Qfactor', 0.7)


Variable-height walls
---------------------

A wall does not have to be flat — its crest can vary along its length.
Specify a different *z* value at each point in the polyline and ANUGA
interpolates linearly along each segment:

.. code-block:: python

    river_walls = {
        'spillway': [
            [20.0,  0.0, 2.0],   # 2 m at the southern end
            [20.0, 10.0, 1.0],   # 1 m notch in the middle
            [20.0, 20.0, 2.0],   # 2 m at the northern end
        ],
    }

You can also pass a per-edge array via :meth:`RiverWall.set_elevation`
after creation.


Parallel simulations
--------------------

Riverwalls work in parallel (MPI) but ``create_riverwalls`` must be called
**after** the domain has been decomposed and distributed — i.e. after
``distribute()`` in a parallel script.  Each process only stores the
riverwall edges that belong to its subdomain.

.. code-block:: python

    domain = anuga.distribute(domain)

    domain.create_riverwalls(river_walls, verbose=False)

    for t in domain.evolve(yieldstep=5.0, finaltime=60.0):
        if anuga.myid == 0:
            domain.print_timestepping_statistics()

    anuga.finalize()


Full example
------------

A complete runnable script is available at
``examples/structures/run_riverwall.py`` in the repository.  It builds a
60 m × 20 m domain with two levees, evolves for 60 s, and demonstrates
raising/lowering the crest heights and changing the Qfactor at runtime.


.. seealso::

   * :doc:`notebook_create_domain_with_riverwalls <../examples/notebook_create_domain_with_riverwalls>`
     — Jupyter notebook showing the same concepts interactively with inline
     animations.
   * :doc:`operators` — other structure operators (culverts, weirs, gates).

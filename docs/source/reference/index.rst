.. _api_reference:

.. currentmodule:: anuga

API Reference
=============

This section documents the main classes in the ANUGA public API.  For a
complete alphabetical index of all names exported by the ``anuga`` package
see the :ref:`genindex`.

The three classes below form the foundation of every ANUGA simulation:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Class
     - Role
   * - :class:`Domain`
     - The simulation domain — owns the mesh, quantities, operators, boundary
       conditions, and the evolve loop.  This is the central object in any
       ANUGA script.
   * - :class:`Quantity`
     - Represents a single physical field (elevation, stage, friction,
       xmomentum, ymomentum) on the mesh.  Accessed via
       ``domain.quantities['stage']`` etc.  Provides methods for setting
       values, computing integrals, interpolation, and extrapolation.
   * - :class:`Region`
     - A spatial subset of the mesh defined by a polygon or list of triangle
       indices.  Used to apply operators or set quantities over a specific
       area of the domain.


.. _api_domain:

Domain
------

The ``Domain`` class is the top-level simulation object.  A typical script:

1. Creates a ``Domain`` from a mesh (rectangular grid or unstructured mesh
   from ``create_domain_from_regions``).
2. Sets initial conditions on its ``Quantity`` objects via
   ``domain.set_quantity()``.
3. Assigns boundary conditions to named boundary tags via
   ``domain.set_boundary()``.
4. Optionally attaches operators (rainfall, culverts, inlets).
5. Calls ``domain.evolve()`` to advance the simulation in time.

.. toctree::
   :hidden:

   anuga.Domain

.. autosummary::
   :nosignatures:

   Domain

:doc:`Full Domain API <anuga.Domain>`

Domain creation functions and the domain methods used throughout this guide:

.. autosummary::
   :toctree: generated
   :nosignatures:

   rectangular_cross_domain
   create_domain_from_regions
   size_field_from_regions
   Size_field


Quantity
--------

``Quantity`` objects are created automatically by the ``Domain`` and are
not usually instantiated directly.  They are accessed via
``domain.quantities``:

.. code-block:: python

   elev = domain.quantities['elevation']
   print(elev.centroid_values)     # values at triangle centroids
   print(elev.vertex_values)       # values at triangle vertices
   print(elev.get_integral())      # integral over the domain

Each quantity allocates only the arrays it needs via the ``qty_type``
parameter (``'evolved'``, ``'centroid_only'``, ``'edge_diagnostic'``,
or ``'coordinate'``).  Gradient and phi arrays are lazy for all types.
See :ref:`quantity-memory-layout` for details.

.. toctree::
   :hidden:

   anuga.Quantity

.. autosummary::
   :nosignatures:

   Quantity

:doc:`Full Quantity API <anuga.Quantity>`


Region
------

A ``Region`` identifies a spatial subset of the mesh.  It is typically
created by passing a polygon to an operator, or directly with
``anuga.Region(domain, polygon=...)`` -- a region may also be given as a
centre and radius, a line, or explicit triangle indices.  It provides the
list of triangle indices that fall inside it.  Operators use regions to
apply forcing (rainfall, extraction) only over a defined area.

.. code-block:: python

   import anuga

   polygon = [[0, 0], [5, 0], [5, 5], [0, 5]]
   region = anuga.Region(domain, polygon=polygon)
   print(region.get_indices())    # triangle IDs inside the polygon

.. autosummary::
   :toctree: generated
   :nosignatures:

   Region


Geo_reference
-------------

``Geo_reference`` records the coordinate reference system (CRS) and local
origin of an ANUGA domain.  It supports WGS84 UTM zones (auto-computed EPSG),
national grids such as RD New (EPSG:28992) or British National Grid
(EPSG:27700), and arbitrary local systems for wavetank simulations.

.. code-block:: python

   geo_ref = anuga.Geo_reference(zone=55, hemisphere='southern',
                                 xllcorner=363000.0, yllcorner=8132000.0)
   geo_ref.epsg          # 32755  (auto-computed)
   geo_ref.is_located()  # True

   # Or supply EPSG directly (zone and hemisphere inferred for UTM codes):
   geo_ref = anuga.Geo_reference(epsg=32755, xllcorner=363000.0, yllcorner=8132000.0)

   # National grids work too:
   geo_ref = anuga.Geo_reference(epsg=28992)   # Netherlands RD New

See :doc:`../setup_anuga_script/coordinate_reference` for full usage guidance.

.. _api_geo_reference:

.. autosummary::
   :toctree: generated
   :nosignatures:

   Geo_reference


.. _api_boundaries:

Boundary conditions
-------------------

The boundary classes assigned to named domain edges via
``domain.set_boundary()``. See
:doc:`../setup_anuga_script/boundaries` for a description of each and when to
use it; the full API of each class is below.

.. currentmodule:: anuga

.. autosummary::
   :toctree: generated
   :nosignatures:

   Reflective_boundary
   Dirichlet_boundary
   Time_boundary
   Transmissive_n_momentum_zero_t_momentum_set_stage_boundary
   Flather_external_stage_zero_velocity_boundary
   File_boundary
   Field_boundary
   Absorbing_wave_boundary
   Characteristic_wave_boundary


.. _api_operators:

Operators
---------

Operators modify domain quantities each timestep (rainfall/extraction, inlets,
culverts, setting quantities, viscosity). See
:doc:`../setup_anuga_script/operators` for categories and usage; the full API of
each is below.

.. currentmodule:: anuga

.. autosummary::
   :toctree: generated
   :nosignatures:

   Rate_operator
   Inlet_operator
   Set_quantity_operator
   Set_stage_operator
   Set_elevation_operator
   Boyd_box_operator
   Boyd_pipe_operator
   Weir_orifice_trapezoid_operator
   Internal_boundary_operator
   Kinematic_viscosity_operator


.. _api_logging:

Logging
-------

Configuring ANUGA's log output. See
:doc:`../setup_anuga_script/logging` for usage.

.. currentmodule:: anuga

.. autosummary::
   :toctree: generated
   :nosignatures:

   set_logfile
   utilities.log.TeeStream
   utilities.log.file_only
   utilities.log.verbose
   utilities.log.info
   utilities.log.warning
   utilities.log.debug
   utilities.log.critical


File format reference
---------------------

Reference documentation for all file formats read and written by ANUGA —
SWW, STS, TMS, mesh files (TSH/MSH), DEM formats, and point data.

.. toctree::
   :hidden:

   file_formats

:doc:`File format reference <file_formats>`


Validation test suite
---------------------

Description of the validation tests in ``validation_tests/``, how to run
them, and what physical benchmarks they cover.

.. toctree::
   :hidden:

   validation

:doc:`Validation test suite <validation>`

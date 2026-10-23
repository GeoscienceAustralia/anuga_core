.. _conventions:

.. currentmodule:: anuga

Conventions and units
=====================

ANUGA uses SI units throughout and a small number of conventions that are worth
knowing before you build a model.

Coordinate system
-----------------

The mesh lives in a 2D Cartesian plane: ``x`` points **east**, ``y`` points
**north**, both in **metres**. Models are usually built in a projected
coordinate system (e.g. a UTM zone) so that distances are metric and locally
undistorted — not in latitude/longitude.

To tie the mesh to a real-world CRS, attach a :class:`Geo_reference` (UTM zone
or EPSG code, plus a local origin ``xllcorner``/``yllcorner``). ANUGA works
internally in coordinates **relative to that origin** to keep the numbers small
and preserve floating-point precision; the origin is added back when results are
written or when you request absolute coordinates. See
:doc:`coordinate_reference` for details.

Quantities and their units
--------------------------

Every model tracks these quantities on the mesh:

.. list-table::
   :header-rows: 1
   :widths: 22 20 58

   * - Quantity
     - Units
     - Meaning
   * - ``elevation``
     - m
     - Bed elevation / bathymetry, measured from the vertical datum.
   * - ``stage``
     - m
     - Absolute water-surface elevation, from the **same** datum as elevation.
   * - ``xmomentum``
     - m²/s
     - Depth-integrated momentum in the x-direction (``depth × x-velocity``).
   * - ``ymomentum``
     - m²/s
     - Depth-integrated momentum in the y-direction (``depth × y-velocity``).
   * - ``friction``
     - –
     - Manning roughness coefficient *n* (dimensionless).

**Depth is derived**, not stored: ``depth = stage - elevation``. A cell is *dry*
where ``stage`` equals ``elevation``. Because ``stage`` is measured from the
datum, it is **not** the water depth — a common source of confusion when setting
initial conditions.

Velocity is recovered from momentum as ``velocity = momentum / depth``; in very
shallow water this is regularised by ``velocity_protection`` (see below) to avoid
dividing by a vanishing depth.

Wetting and drying
------------------

ANUGA models moving wet/dry fronts. A cell is treated as dry once its depth
falls below ``minimum_allowed_height`` (default ``1.0e-5`` m); such cells carry
no momentum. Only depths above ``minimum_storable_height`` (default ``1.0e-3``
m) are written to the SWW output, so shallow films do not clutter results.

Physical defaults and constants
-------------------------------

Defined in :mod:`anuga.config` and used unless you override them:

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Constant
     - Default
     - Meaning
   * - ``g``
     - 9.8 m/s²
     - Gravitational acceleration.
   * - ``manning``
     - 0.03
     - Default Manning friction if ``friction`` is not set.
   * - ``minimum_allowed_height``
     - 1.0e-5 m
     - Depth below which a cell is considered dry.
   * - ``velocity_protection``
     - 1.0e-6
     - Regularises velocity = momentum / depth in shallow water.

Time
----

Simulation time is in **seconds**, measured from the Unix epoch —
**00:00:00 UTC on 1 January 1970** — so the default start of ``t = 0``
corresponds to that instant. Set a real-world start with
:meth:`~Domain.set_starttime`, which accepts either a number of seconds or a
Python ``datetime`` object; use the standard library ``datetime`` and
``zoneinfo`` modules to build a timezone-aware start time (and
:meth:`~Domain.set_timezone` to control how times are reported).

``finaltime`` / ``duration`` in :meth:`~Domain.evolve` are likewise in seconds.
The internal timestep is chosen automatically to satisfy the CFL stability
condition and is generally much smaller than the ``yieldstep`` — see
:doc:`evolve`.

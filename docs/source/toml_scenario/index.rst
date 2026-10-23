.. _toml_scenario:

=====================================
Running Scenarios with anuga_toml_run
=====================================

.. currentmodule:: anuga

The ``anuga_toml_run`` script provides a ready-made runner for ANUGA flood and
tsunami scenarios.  All simulation inputs are described in a single
`TOML <https://toml.io>`_ configuration file — no Python coding required for
standard setups.  The same script also accepts legacy Excel (``.xlsx``) files
from the older ``cairns_excel`` interface.

.. contents:: Contents
   :local:
   :depth: 2


Quick Start
-----------

**Serial run:**

.. code-block:: bash

   anuga_toml_run  path/to/scenario.toml

**Parallel run (MPI):**

.. code-block:: bash

   mpirun -np 6  anuga_toml_run  path/to/scenario.toml

All relative paths inside the TOML file are resolved relative to the directory
that contains the TOML file, so the script can be invoked from any working
directory.


Eject to a standalone Python script
-----------------------------------

When you need to go beyond what the TOML supports — a custom operator, a
scripted boundary function, bespoke post-processing — generate an equivalent
standalone script and edit it::

   anuga_toml_run --emit-script my_run.py  path/to/scenario.toml

This writes ``my_run.py`` (and exits without running). The script parses the
same TOML via ``PrepareData`` and drives the standard phases (mesh → initial
conditions → forcing/structures → boundaries → evolve), so it runs identically
to ``anuga_toml_run``::

   python my_run.py
   mpiexec -np 4 python my_run.py        # legacy (CPU) parallel

It is ordinary, editable Python: delete phases you do not need, splice in your
own operators or boundaries between phases, or replace any ``setup_*`` call with
hand-written ANUGA code. Keep the script alongside the TOML (it resolves
TOML-relative paths from its own directory).


Working example — Cairns tsunami scenario
------------------------------------------

A complete, ready-to-run example is provided in
``examples/run_toml/cairns/`` of the repository.  It models a synthetic
tsunami entering Cairns Harbour (Queensland, Australia) and exercises most
TOML features: a shapefile boundary, a DEM elevation raster, Flather open-
ocean boundaries, and one each of a bridge, pumping station, culvert, and weir.
The (≈9 MB) DEM is shared under ``examples/data/cairns/``.

.. code-block:: bash

   cd examples/run_toml/cairns
   anuga_toml_run cairns_example.toml          # serial
   mpirun -np 4 anuga_toml_run cairns_example.toml   # parallel

The same scenario is also available through a legacy Excel front-end under
``examples/cairns_toml_excel/`` (see *Excel Compatibility* below).


Minimal annotated TOML
-----------------------

The snippet below is the smallest valid TOML that will run a real simulation.
Copy it, adjust the file paths, and add extra sections as needed.

.. code-block:: toml

   # ── Project ───────────────────────────────────────────────────────────────
   [project]
   scenario                  = "my_run"      # name for the SWW file & output dir
   output_base_directory     = "OUTPUT/"
   yieldstep                 = 120.0         # seconds between evolve yields
   finaltime                 = 21600.0       # total simulation duration [s]
   projection_information    = "EPSG:32755"  # "EPSG:<code>", UTM zone int (e.g. -55), or proj4 string
   flow_algorithm            = "DE0"         # "DE0" (fast) or "DE1" (accurate)
   output_tif_cellsize       = 50.0          # cell size [m] for output GeoTiff rasters

   # ── Mesh ──────────────────────────────────────────────────────────────────
   [mesh]
   bounding_polygon          = "mesh/boundary.shp"   # shapefile with boundary tags
   default_res               = 1000000.0             # default triangle area [m²]
   # [[mesh.interior_regions]]
   # polygon    = "mesh/fine_zone.csv"
   # resolution = 10000.0   # finer triangles inside this polygon [m²]

   # [[mesh.interior_holes]]
   # polygon = "mesh/buildings.csv"   # cut OUT of the mesh entirely
   # tag     = "building"             # optional; needs a matching boundary condition

   # ── Boundary conditions ───────────────────────────────────────────────────
   [boundary_conditions]
   boundary_tags_attribute_name = "Boundary"   # shapefile attribute holding tag names

   [[boundary_conditions.boundaries]]
   tag        = "Ocean"
   type       = "Flather_Stage"                # weakly-reflecting open boundary
   file       = "boundarycond/tide.csv"        # CSV: columns time_s, stage_m
   start_time = 0.0

   [[boundary_conditions.boundaries]]
   tag  = "Land"
   type = "Reflective"                         # solid wall — no file needed

   # ── Initial conditions ────────────────────────────────────────────────────
   [initial_conditions]

   [[initial_conditions.elevation]]
   polygon  = "Extent"                         # full extent of the raster
   value    = "initialcond/dem.asc"            # GDAL raster (asc, tif, …)
   clip_min = -inf
   clip_max = inf

   [[initial_conditions.elevation]]
   polygon  = "All"                            # catch-all fallback
   value    = 0.0
   clip_min = -inf
   clip_max = inf

   [[initial_conditions.friction]]
   polygon  = "All"
   value    = 0.03                             # Manning's n
   clip_min = 0.01
   clip_max = 0.5

   [[initial_conditions.stage]]
   polygon  = "All"
   value    = 0.0                              # start domain dry
   clip_min = -inf
   clip_max = inf

The full reference for every key is in the `Configuration File Reference`_
section below.


Scenario Directory Layout
--------------------------

A typical scenario directory looks like this::

   my_scenario/
   ├── scenario.toml          ← main configuration file
   ├── user_functions.py      ← optional custom callbacks (see below)
   ├── mesh/
   │   ├── boundary.shp       ← domain boundary shapefile (or .csv)
   │   └── ...
   ├── boundarycond/
   │   ├── tide.csv
   │   └── ...
   └── initialcond/
       └── dem.asc

On completion, outputs are written to a timestamped subdirectory under
``output_base_directory`` (e.g. ``OUTPUT/RUN_20260315_120000_scenario/``).
A copy of all code files — including the TOML config and the runner script
itself — is archived there for reproducibility.


Configuration File Reference
-----------------------------

The TOML file is divided into sections (TOML *tables*).  Every key shown
without a comment is **required**; keys shown with a ``# default:`` comment
are optional.

.. note::

   TOML has no arithmetic, so ``finaltime = 5*60`` is a syntax error. Numeric
   fields therefore also accept a **quoted** arithmetic expression, which the
   parser evaluates — e.g. ``finaltime = "5*60"`` gives ``300``,
   ``finaltime = "7*24*3600"`` gives one week. Only numeric literals and
   ``+ - * / // % **`` with parentheses are permitted (no names or function
   calls). Plain numbers such as ``finaltime = 21600.0`` work as before.

.. _toml-project:

[project]
~~~~~~~~~

Top-level simulation settings.

.. code-block:: toml

   [project]

   # Name used for the output SWW file and the output directory prefix
   scenario = "my_scenario"

   # Base directory where per-run output subdirectories are created
   output_base_directory = "OUTPUT/"

   # Time between successive yields of the evolve loop [seconds].
   # Controls how often timestepping statistics are printed.
   yieldstep = 120.0

   # Total simulation duration [seconds]
   finaltime = 21600.0

   # Model time at the start of the run [seconds]. default: 0.0
   # A negative value runs a burn-in before the forcing data, which start at
   # t = 0: inlet and stage-boundary series hold their first value and rain is
   # zero until t = 0. Model times then match the forcing and observation
   # clocks, with no shift. finaltime is on the same clock, so
   # starttime = "-3*3600" with finaltime = "24*3600" evolves 27 h.
   # starttime = "-3*3600"

   # SWW output interval [seconds].
   # Must be an integer multiple of yieldstep.
   # Omit to write SWW output at every yieldstep (the default).
   # outputstep = 600.0

   # Coordinate reference system. One of:
   #   Integer       — UTM zone (negative = southern hemisphere), e.g. -55
   #   "EPSG:<code>" — an EPSG code, e.g. "EPSG:32755" (= UTM zone 55 south)
   #   String        — a full proj4 string for other projections
   projection_information = "EPSG:32755"

   # Numerical flow algorithm.
   #   "DE0" — first-order (faster)
   #   "DE1" — second-order (more accurate)
   flow_algorithm = "DE0"

   # Enable local extrapolation and flux updating (improves stability
   # near wet/dry fronts)
   # default: false
   use_local_extrapolation_and_flux_updating = false

   # Grid cell size [m] for output GeoTiff rasters
   output_tif_cellsize = 50.0

   # Polygon CSV to clip output GeoTiffs.
   # Leave empty to use the full domain extent.
   output_tif_bounding_polygon = ""

   # How often (every N timesteps) to update tracked maximum quantities
   # default: 1
   max_quantity_update_frequency = 1

   # Time from which maximum-quantity tracking begins [seconds]
   # Must be < finaltime.  default: 0.0 (so a negative starttime burn-in is excluded)
   max_quantity_collection_starttime = 0.0

   # Store each triangle's three vertex values separately.
   # Produces larger files; needed for some visualisation workflows.
   # default: false
   store_vertices_uniquely = false

   # Re-write elevation at every yieldstep.
   # Only needed if elevation changes during the run (erosion operators,
   # sediment transport with bed evolution); those default it to true.
   # default: false
   store_elevation_every_timestep = false

   # Directory for text-format copies of spatial inputs (used for QC)
   # default: "SPATIAL_TEXT"
   spatial_text_output_dir = "SPATIAL_TEXT"

   # --- Per-yieldstep diagnostics ---
   # default: false for all
   report_mass_conservation_statistics      = false
   report_peak_velocity_statistics          = false   # requires user_functions.py
   report_smallest_edge_timestep_statistics = false
   report_operator_statistics               = false   # requires user_functions.py

   # --- Solver threading ---
   # Number of OpenMP threads.  Omit to read OMP_NUM_THREADS env var (default 1).
   # omp_num_threads = 4

   # Compute mode: "legacy" = OpenMP CPU (default), "unified" = shared
   # CPU/GPU kernels, offloading only when GPU offload is enabled
   # process-wide (experimental).  The older integer multiprocessor_mode
   # (1 or 2) is still accepted.
   compute_mode = "legacy"

.. _toml-mesh:

[mesh]
~~~~~~

Mesh geometry and resolution.

.. code-block:: toml

   [mesh]

   # Set true to reuse a previously built partitioned mesh (skips mesh generation)
   # default: false
   use_existing_mesh_pickle = false

   # Domain boundary file — two formats are supported:
   #
   #   Shapefile (.shp): each line feature must carry a tag attribute whose name
   #     is given by boundary_conditions.boundary_tags_attribute_name.
   #
   #   CSV (.csv): a plain x,y polygon.  You must also provide
   #     [[mesh.boundary_tags]] entries (see below).
   bounding_polygon = "mesh/boundary.shp"

   # Default maximum triangle area [m²] across the whole domain
   default_res = 1000000.0

   # Breakline files — polylines that force triangle edges.
   # Glob patterns are accepted (e.g. "mesh/breaklines_*.csv").
   # default: []
   breakline_files = []

   # Riverwall CSV files — like breaklines but each row has three columns:
   #   x [m], y [m], crest_elevation [m]
   # An optional first line may contain hydraulic parameters, e.g.:
   #   Qfactor: 1.5, s1: 0.94, Cd_through: 0.1
   # Glob patterns are accepted (e.g. "mesh/riverwalls_*.csv").
   # default: []
   riverwall_csv_files = []

   # Snapping distance [m] for breakline / bounding-polygon intersections.
   # Set to "ignore" to skip intersection processing.
   # default: 0.1
   breakline_intersection_threshold = 0.1

   # Point-based resolution file — CSV with columns x, y, resolution.
   # Mutually exclusive with [[mesh.interior_regions]].
   region_areas_file = ""

   # [[mesh.interior_holes]] — polygons cut OUT of the mesh entirely, leaving a
   # void rather than a refined region. Use for anything water should neither
   # enter nor flow through: building footprints, tank pads, solid structures.
   #
   # Unlike [[mesh.interior_regions]], which keeps the triangles and only changes
   # their size, so no resolution is given.
   #
   # 'tag' is optional and names the hole's edges so a boundary condition can be
   # bound to them. If set, a matching [[boundary_conditions.boundaries]] entry is
   # REQUIRED, or the run stops with
   #   Tag "..." has not been bound to a boundary object
   # Omit it and the mesh generator applies its own 'interior' default.
   #
   # [[mesh.interior_holes]]
   # polygon = "mesh/buildings.csv"
   # tag     = "building"

   # Interpretation of the resolution column in region_areas_file.
   #   "area"   — maximum triangle area [m²]
   #   "length" — approximate maximum edge length [m]
   # default: "area"
   region_areas_type = "area"

**Interior regions** (polygon-based resolution refinement):

.. code-block:: toml

   [[mesh.interior_regions]]
   polygon    = "mesh/fine_region.csv"
   resolution = 10000.0   # maximum triangle area [m²]

Interior regions and ``region_areas_file`` / breaklines are mutually exclusive.

**Explicit boundary tags** (required when ``bounding_polygon`` is a CSV):

.. code-block:: toml

   [[mesh.boundary_tags]]
   tag   = "ocean"
   edges = [0, 1]

   [[mesh.boundary_tags]]
   tag   = "land"
   edges = [2, 3]

Edge indices are 0-based and refer to the segments between consecutive polygon
vertices.

.. _toml-boundary:

[boundary_conditions]
~~~~~~~~~~~~~~~~~~~~~

One ``[[boundary_conditions.boundaries]]`` entry is required for every tag
that appears in the bounding polygon.

.. code-block:: toml

   [boundary_conditions]

   # Shapefile attribute name that carries the boundary tag.
   # Must match the attribute used in the boundary shapefile.
   boundary_tags_attribute_name = "Boundary"

   # Reflective (solid wall) boundary — no file needed
   [[boundary_conditions.boundaries]]
   tag  = "Land"
   type = "Reflective"

   # Stage boundary — transmissive with prescribed stage timeseries
   [[boundary_conditions.boundaries]]
   tag        = "Ocean"
   type       = "Stage"
   file       = "boundarycond/tide.csv"   # columns: time_s, stage_m
   start_time = 0.0                       # subtract from CSV times [s]

   # Flather radiation boundary with prescribed stage timeseries
   [[boundary_conditions.boundaries]]
   tag        = "Open"
   type       = "Flather_Stage"
   file       = "boundarycond/tide.csv"
   start_time = 0.0

.. _toml-initial:

[initial_conditions]
~~~~~~~~~~~~~~~~~~~~~

Each quantity is built from an ordered list of ``[[initial_conditions.<qty>]]``
entries.  Earlier entries take priority over later ones.

Supported quantities: ``elevation``, ``friction``, ``stage``, ``xmomentum``,
``ymomentum``.

Each quantity also supports an optional spatial-averaging step controlled by
a ``{qty}_spatial_average`` key at the ``[initial_conditions]`` level.
When set, the quantity field is averaged onto a regular grid at the given
spacing (metres) before being applied to the mesh — useful for smoothing
noisy DEMs or variable-resolution friction fields:

.. code-block:: toml

   [initial_conditions]
   elevation_spatial_average  = 100.0   # smooth DEM at 100 m grid
   friction_spatial_average   = 50.0    # smooth friction at  50 m grid
   # stage_spatial_average    = 50.0    # also valid for stage, xmomentum, ymomentum

   [[initial_conditions.elevation]]
   polygon  = "Extent"                    # apply over full raster extent
   value    = "initialcond/dem.asc"       # GDAL-compatible raster
   clip_min = -inf
   clip_max = inf

   [[initial_conditions.elevation]]
   polygon  = "All"                       # catch-all fallback
   value    = 0.0
   clip_min = -inf
   clip_max = inf

   [[initial_conditions.friction]]
   polygon  = "mesh/mannings_zone.csv"    # apply inside polygon
   value    = 0.025
   clip_min = 0.01
   clip_max = 0.5

   [[initial_conditions.friction]]
   polygon  = "All"
   value    = 0.03
   clip_min = 0.01
   clip_max = 0.5

**Polygon options:**

===============  ==============================================================
``"Extent"``     Apply over the full spatial extent of the raster ``value``
``"All"``        Apply everywhere (use as the final catch-all entry)
``"None"``       Skip this entry
path to file     ``.shp`` or ``.csv`` polygon
===============  ==============================================================

**Value options:**

====================  =========================================================
constant number       applied uniformly across the polygon
path to raster        GDAL-compatible raster (``.asc``, ``.tif``, …)
path to CSV/TXT       file with columns ``x``, ``y``, ``z``
====================  =========================================================

.. _toml-additions:

[initial_condition_additions]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Same structure as ``[initial_conditions]``, but values are *added on top of*
whatever was set there.  Useful for local offsets such as a pre-dredged
channel, a localised initial flood pool, or a small friction adjustment.

.. code-block:: toml

   [initial_condition_additions]

   [[initial_condition_additions.elevation]]
   polygon  = "mesh/dredge_channel.csv"
   value    = -2.0
   clip_min = -inf
   clip_max = inf

.. _toml-inlets:

[[inlets]]
~~~~~~~~~~

Line sources / sinks of discharge across a cross-section.

.. code-block:: toml

   [[inlets]]
   name            = "river_inflow"
   line_file       = "mesh/river_line.csv"       # cross-section polyline
   timeseries_file = "boundarycond/river_flow.csv"  # columns: time_s, discharge_m3s
   start_time      = 0.0

Positive discharge adds water; negative removes it.  Multiple ``[[inlets]]``
entries are supported.

.. _toml-rainfall:

[[rainfall]]
~~~~~~~~~~~~

Spatially uniform rainfall rate over a region.

.. code-block:: toml

   [[rainfall]]
   timeseries_file    = "boundarycond/rain.csv"  # columns: time_s, rate_mm_per_hr
   start_time         = 0.0
   interpolation_type = "linear"                 # "linear" or "nearest"
   polygon            = "All"                    # path to CSV polygon, or "All"
   multiplier         = 1.0                      # scale factor (default 1.0)

Multiple ``[[rainfall]]`` entries are supported.

.. _toml-bridges:

[[bridges]]
~~~~~~~~~~~

HEC-RAS rating-curve internal boundary operator.  Set ``enabled = false`` to
disable a bridge without removing its definition.

.. code-block:: toml

   [[bridges]]
   enabled                      = true
   label                        = "bridge_1"
   deck_file                    = "mesh/bridge1_deck.csv"
   deck_elevation               = 4.5         # [m] bed elevation forced inside deck polygon
   exchange_line_0              = "mesh/bridge1_exchange_up.csv"
   exchange_line_1              = "mesh/bridge1_exchange_down.csv"
   enquiry_gap                  = 2.0         # [m] from exchange lines to enquiry points
   internal_boundary_curve_file = "mesh/bridge1_rating_curve.csv"
   vertical_datum_offset        = 0.0         # [m] added to curve elevations
   smoothing_timescale          = 20.0        # [s] exponential smoothing constant

.. _toml-culverts:

[[culverts]]
~~~~~~~~~~~~

Box or pipe culverts using the Boyd (1987) head-discharge algorithm.
Set ``enabled = false`` to disable without removing the definition.

``type`` selects the cross-section shape:

* ``"boyd_box"``  — rectangular barrel (``width`` × ``height``).
* ``"boyd_pipe"`` — circular barrel (``diameter``).

**Geometry** — choose exactly one of:

* ``exchange_line_0`` / ``exchange_line_1``: paths to CSV polyline files
  that define the upstream and downstream exchange zones.  The exchange
  zone widths govern how much of the mesh perimeter couples to the culvert.
  This is the preferred approach for field models.

* ``end_point_0`` / ``end_point_1``: ``[x, y]`` coordinate pairs for the
  two barrel ends.  ANUGA derives short perpendicular exchange lines
  automatically.

**Boyd box example:**

.. code-block:: toml

   [[culverts]]
   enabled              = true
   type                 = "boyd_box"
   label                = "road_culvert_1"
   width                = 0.9              # [m] internal barrel width
   height               = 0.6             # [m] internal barrel height (omit = square)
   exchange_line_0      = "mesh/culvert1_exchange_up.csv"
   exchange_line_1      = "mesh/culvert1_exchange_down.csv"
   enquiry_gap          = 0.2             # [m] from exchange line to enquiry point
   losses               = 0.5             # head-loss coefficient (0.5 = sharp-edge inlet)
   barrels              = 1.0             # number of parallel identical barrels
   blockage             = 0.0             # fractional blockage [0,1]; 0 = fully open
   z1                   = 0.0             # batter slope at end 0 (rise/run)
   z2                   = 0.0             # batter slope at end 1
   apron                = 0.1             # [m] flat apron at each end
   manning              = 0.013           # barrel Manning's n
   smoothing_timescale  = 0.0             # [s] exponential smoothing constant
   use_momentum_jet     = true
   use_velocity_head    = true
   # invert_elevations  = [1.0, 0.8]      # [m] upstream, downstream; sampled from DEM if omitted

**Boyd pipe example:**

.. code-block:: toml

   [[culverts]]
   enabled              = true
   type                 = "boyd_pipe"
   label                = "drain_pipe_1"
   diameter             = 0.6              # [m] internal barrel diameter
   exchange_line_0      = "mesh/pipe1_exchange_up.csv"
   exchange_line_1      = "mesh/pipe1_exchange_down.csv"
   enquiry_gap          = 0.2
   losses               = 0.5
   barrels              = 2.0              # two parallel pipes
   manning              = 0.013
   smoothing_timescale  = 0.0

**Using end points instead of exchange lines:**

.. code-block:: toml

   [[culverts]]
   type        = "boyd_box"
   label       = "simple_culvert"
   width       = 0.9
   end_point_0 = [355420.0, 8132050.0]    # [x, y] upstream barrel end [m]
   end_point_1 = [355435.0, 8132050.0]    # [x, y] downstream barrel end [m]
   losses      = 0.5

Multiple ``[[culverts]]`` entries are supported.


.. _toml-weirs:

[[weirs]]
~~~~~~~~~

Weir / orifice structures with a trapezoidal cross-section, using combined
weir and orifice flow formulae (``Weir_orifice_trapezoid_operator``).
Set ``enabled = false`` to disable without removing the definition.

Geometry is specified with the same ``exchange_line_0`` / ``exchange_line_1``
(file paths) or ``end_point_0`` / ``end_point_1`` (coordinates) choice as for
``[[culverts]]``.

.. code-block:: toml

   [[weirs]]
   enabled              = true
   label                = "outlet_weir"
   width                = 3.0             # [m] bottom width of trapezoidal section
   height               = 1.2             # [m] section height (omit = width)
   exchange_line_0      = "mesh/weir1_exchange_up.csv"
   exchange_line_1      = "mesh/weir1_exchange_down.csv"
   enquiry_gap          = 0.0             # [m]; 0 is typical for weirs
   losses               = 0.5             # head-loss coefficient
   barrels              = 1.0
   blockage             = 0.0
   z1                   = 0.0
   z2                   = 0.0
   apron                = 0.1
   manning              = 0.013
   smoothing_timescale  = 0.0
   use_momentum_jet     = true
   use_velocity_head    = true
   # invert_elevations  = [1.0, 0.8]      # [m]; sampled from DEM if omitted

Multiple ``[[weirs]]`` entries are supported.


.. _toml-pumping:

[[pumping_stations]]
~~~~~~~~~~~~~~~~~~~~~

Pump operator transferring water between a wet-well basin and a discharge
point.  Set ``enabled = false`` to disable without removing the definition.



.. code-block:: toml

   [[pumping_stations]]
   enabled               = true
   label                 = "pump_1"
   pump_capacity         = 2.0     # [m³/s] maximum flow rate
   pump_rate_of_increase = 0.5     # [m³/s per second]
   pump_rate_of_decrease = 0.5     # [m³/s per second]
   hw_to_start_pumping   = 0.3     # headwater depth [m] to switch on
   hw_to_stop_pumping    = 0.1     # headwater depth [m] to switch off
   basin_polygon_file    = "mesh/pump1_basin.csv"
   basin_elevation       = -1.0    # [m] bed elevation inside basin polygon
   exchange_line_0       = "mesh/pump1_wet_well_line.csv"
   exchange_line_1       = "mesh/pump1_discharge_line.csv"
   smoothing_timescale   = 30.0    # [s]

Multiple ``[[pumping_stations]]`` entries are supported.


[[erosion]]
~~~~~~~~~~~

Bed erosion / scour operators. Multiple ``[[erosion]]`` entries are supported.

.. code-block:: toml

   [[erosion]]

   # Which erosion behaviour to apply. One of:
   #   simple      — base erosion operator
   #   bed_shear   — scour driven by bed shear stress (dam-breach style)
   #   flat_slice  — erode down to a flat surface at a given elevation
   #   flat_fill   — as flat_slice, filling rather than only cutting
   #   sand_dune   — sand-dune erosion, limited by an angle of repose
   type = "bed_shear"

   # Region the operator acts on — give EITHER a polygon OR center + radius.
   # (Raw triangle indices are not supported here: they cannot survive a
   # re-mesh, so they have no stable meaning in a configuration file.)
   polygon = "erosion/breach.csv"
   # center = [382000.0, 6354000.0]
   # radius = 50.0

   # Bed shear below which no scour occurs. 0.0 means scour begins immediately.
   # default: 0.0
   threshold = 0.0

   # Minimum bed level [m]; scour does not cut below this.
   # default: 0.0
   base = -5.0

   # Optional label naming this operator's log file when logging = true
   label   = "breach"
   logging = false

   # ── Type-specific parameters ──────────────────────────────────────────────
   # Supplying one on the wrong type is an ERROR, not silently ignored.
   #
   #   shear_factor   bed_shear only   — higher slows the breach (default 75000.0)
   #   elevation      flat_slice / flat_fill only — target surface level [m]
   #   Ra             sand_dune only   — angle of repose [degrees] (default 34.0)
   shear_factor = 75000.0

.. warning::

   Erosion changes bed elevation as the run proceeds, but ANUGA stores elevation
   **once at t=0** by default — so the eroded bed will not appear in the ``.sww``
   file or in any raster derived from it. Set ``store_elevation_every_timestep =
   true`` in ``[project]`` when using erosion operators. A warning is issued if
   ``[[erosion]]`` is present without it.

[sediment]
~~~~~~~~~~

Sediment transport: suspended load carried as tracers, with deposition and
entrainment, and optionally bedload and an evolving bed. The table maps onto
:meth:`~anuga.Domain.initialize_sediment_operator` and its setters, each
``[[sediment.fractions]]`` onto :meth:`~anuga.Domain.add_sediment_fraction`.
See :doc:`/setup_anuga_script/sediment` for what every choice means and its
default; an absent key leaves that default.

.. code-block:: toml

   [sediment]

   # ── Domain-wide ──────────────────────────────────────────────────────────
   porosity      = 0.3      # bed porosity                          default: 0.3
   c_max         = 0.3      # ceiling on depth-averaged concentration  default: 0.3
   bed_evolution = true     # move the bed; false = fixed-bed stage   default: true
   # c_pack = 0.65          # near-bed packing ceiling
   # rho_w  = 1000.0        # water density [kg/m^3]

   # How bed shear stress is formed: quadratic_drag (default; velocity based),
   # depth_slope (bed slope; for reproducing anugaSed) or energy_slope
   # (free-surface slope; what the old Bed_shear_erosion_operator used).
   # For the two slope closures, max_slope caps the slope, and for
   # depth_slope freeze_slope = true takes the slope from the initial bed
   # and keeps it, so the bed can evolve without the closure feeding on the
   # roughness it creates.
   shear_closure = "quadratic_drag"
   # max_slope = 0.05
   # freeze_slope = true

   # Which erosion law, by naming the bed: noncohesive (default; Shields
   # entrainment, threshold per fraction via tau_c_star), cohesive
   # (Hanson & Simon excess shear) or partheniades. tau_crit [Pa] and K_e
   # apply to the cohesive and partheniades routes only.
   bed_material = "noncohesive"
   # tau_crit = 0.088
   # K_e = 1.0e-6
   # Non-cohesive entrainment: smith_mclean (default) or de_leeuw [E-6], with
   # de_leeuw_fit = de_leeuw_2020 (default) or nghiem_2022 and an optional
   # skin_roughness [m] (default 3 x each fraction's diameter).
   # entrainment = "de_leeuw"
   # de_leeuw_fit = "nghiem_2022"
   # skin_roughness = 0.002

   # Deposition: d_star (default) or threshold (with tau_d [Pa]); the
   # near-bed profile factor constant (default) or rouse.
   deposition_law = "d_star"
   near_bed       = "constant"
   # tau_d = 0.0
   # reference_height_floor = 0.01

   # Friction feeding the sediment shear: constant (default; the domain's
   # Manning n), larsen_lamb (needs k_s or sigma_br) or wilson (needs
   # grain_size, and bed = sand | gravel | boulder).
   # friction_mode = "constant"

   # Bedload: off (default), wong_parker_eq24, wong_parker_eq23,
   # engelund_hansen or grass; bedload_K, bedload_m, bedload_tau_c_star
   # override the formula's constants (grass needs bedload_K, its A_g).
   # Boundary tags through which bedload passes (zero gradient); all other
   # boundaries are closed to it.
   # bedload = "wong_parker_eq24"
   # bedload_open_boundaries = ["inflow", "outflow"]

   # Relax bed slopes steeper than this [degrees]; absent = off.
   # angle_of_repose = 35.0
   # repose_relax = 1.0
   # repose_max_sweeps = 50

   # Bedrock: erosion never cuts below it. Give ONE of an elevation [m] or an
   # erodible thickness [m] below the initial bed.
   # erodible_base_elevation = -2.0
   # erodible_base_depth = 1.5

   # ── Fractions (at least one) ─────────────────────────────────────────────
   [[sediment.fractions]]
   name     = "sand"
   diameter = 2.0e-4        # grain diameter [m]
   # rho_s = 2650.0         # particle density [kg/m^3]
   # tau_c_star = 0.04      # critical Shields stress (noncohesive route); 0 = no entrainment
   # d_star = 1.0           # near-bed / depth-averaged concentration ratio
   # initial_concentration = 0.0   # volumetric
   # reference_height = 4.0e-4     # [m], for near_bed = "rouse"
   # nu = 1.0e-6            # settling-velocity constants (Ferguson & Church)
   # C1 = 18.0
   # C2 = 0.4
   # Concentration carried in where water enters across a boundary tag;
   # unset tags bring clean water.
   [sediment.fractions.boundary]
   left = 0.01

   [[sediment.fractions]]
   name     = "silt"
   diameter = 2.0e-5
   tau_c_star = 0.11

   # ── Erodible regions (optional) ──────────────────────────────────────────
   # Restrict erosion to part of the domain, or lock part of it
   # (erodible = false). A polygon file OR center + radius, as for [[erosion]].
   # [[sediment.erodible_regions]]
   # polygon = "sediment/channel.csv"
   # erodible = true

.. note::

   With bed evolution on, elevation is stored every timestep by default
   (``store_elevation_every_timestep`` in ``[project]``), so the evolving bed
   appears in the ``.sww``; a warning is issued if it is explicitly set to
   ``false``. Each fraction's concentration is stored as ``<name>_c``.

.. warning::

   Nothing prevents ``[sediment]`` and ``[[erosion]]`` in the same scenario.
   Both write elevation and the bed changes simply add; one conserves mass and
   the other does not. Pick one.


Input Validation
-----------------

The TOML parser validates all required fields and value ranges before the
simulation starts.  Errors are collected and reported together so you can
fix all problems in a single pass rather than re-running after each one.

A missing required field or an out-of-range value produces a message like::

   TOML configuration errors in 'scenario.toml':
     [project] 'scenario' is required but missing
     [project] 'flow_algorithm' must be one of ('DE0', 'DE1', 'DE2', 'DE0_7', 'DE1_7', 'DE_ader2') — got 'de0'
     culverts['road_culvert_1'] 'width' must be > 0 — got -0.9

The run aborts after reporting all errors.


Custom Callbacks: user_functions.py
-------------------------------------

Place a file called ``user_functions.py`` alongside the TOML file to add
custom per-timestep actions.  The runner imports it automatically; if the file
is absent it is silently skipped.

Two hooks are called when the corresponding TOML flags are ``true``:

``print_velocity_statistics(domain, max_quantities)``
   Called when ``report_peak_velocity_statistics = true``.
   Receives the live domain and the
   :class:`~anuga.operators.collect_max_quantities_operator.Collect_max_quantities_operator`
   instance.

``print_operator_inputs(domain)``
   Called when ``report_operator_statistics = true``.
   Receives the live domain; useful for printing rainfall rates, inlet
   discharges, etc.

A minimal example:

.. code-block:: python

   # user_functions.py
   from anuga.parallel import myid, barrier

   def print_velocity_statistics(domain, max_quantities):
       if myid == 0:
           import numpy as np
           xx = domain.quantities['xmomentum'].centroid_values
           yy = domain.quantities['ymomentum'].centroid_values
           dd = domain.quantities['stage'].centroid_values \
               - domain.quantities['elevation'].centroid_values
           dd = np.maximum(dd, 1.0e-3)
           speed = np.sqrt(xx**2 + yy**2) / dd
           print(f'  Peak speed: {speed.max():.3f} m/s')

   def print_operator_inputs(domain):
       if myid == 0:
           for op in domain.fractional_step_operators:
               if hasattr(op, 'applied_Q'):
                   print(f'  {op.label}: Q = {op.applied_Q:.3f} m³/s')


Output Directory Structure
---------------------------

Each run creates a timestamped directory under ``output_base_directory``::

   OUTPUT/
   └── RUN_20260315_120000_scenario/
       ├── scenario.sww           ← NetCDF output (stage, momentum, elevation)
       ├── Simulation_logfile.log ← full stdout log
       ├── Max_quantities_*.csv   ← peak stage, depth, speed per triangle
       ├── *.tif                  ← GeoTiff rasters of peak quantities
       ├── code/                  ← archived copy of all input files
       │   ├── scenario.toml
       │   ├── anuga_toml_run     ← copy of the runner script
       │   └── user_functions.py
       └── SPATIAL_TEXT/          ← text copies of spatial inputs (for QC)


Excel Compatibility
--------------------

``anuga_toml_run`` also accepts legacy Excel files::

   anuga_toml_run  path/to/ANUGA_setup.xlsx

The Excel format is described in the ``cairns_toml_excel`` example directory.
Attributes that exist only in the TOML interface
(``compute_mode``, ``omp_num_threads``, ``outputstep``,
``report_operator_statistics``) are set to sensible defaults when reading
Excel files.

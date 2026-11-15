.. currentmodule:: anuga

.. _rainfall-operators:

Rainfall Operators
==================

ANUGA provides several classes for applying rainfall to a simulation domain.
All rainfall operators inherit from :class:`~anuga.operators.rate_operators.Rate_operator`
and are GPU-compatible — they precompute centroid rates during initialisation so
that only a fast array lookup occurs each yield step, with no file I/O during
the evolve loop.

Simple scalar or callable rainfall
-----------------------------------

For uniform or time-varying rainfall use :class:`Rate_operator` directly:

.. code-block:: python

   import anuga

   domain = anuga.rectangular_cross_domain(10, 5)

   # Constant 10 mm/hr = 10/3 600 000 m/s over the whole domain
   rain = anuga.Rate_operator(domain, rate=10.0/3.6e6)

   # Or a function of time (e.g. a Gaussian pulse)
   import math
   rain = anuga.Rate_operator(domain,
                              rate=lambda t: math.exp(-((t - 300)**2) / 3600),
                              factor=0.001)


Hyetograph rainfall
-------------------

:class:`Hyetograph_rate_operator` applies rain that is the same everywhere and
changes in steps through time: a gauge record or a design storm. Give the
times that bound the intervals, and the depth (mm) or the intensity (mm/h) in
each interval:

.. code-block:: python

   # 10, 30 and 20 mm in three 5-minute intervals
   rain = anuga.Hyetograph_rate_operator(
       domain,
       times=[0.0, 300.0, 600.0, 900.0],      # s
       depths=[10.0, 30.0, 20.0])             # mm, or intensities=[...] in mm/h

   for t in domain.evolve(yieldstep=300.0, finaltime=3600.0):
       domain.print_timestepping_statistics()

The intervals need not be the same length. Outside them the rate is
``default_rate`` (0 by default). As with the other operators here, end each
``yieldstep`` on a change of rate for an exact total depth. It takes the same
``polygon``, ``region`` and ``losses`` arguments as the operators below.

Raster time-series rainfall
-----------------------------

:class:`Raster_rate_operator` drives rainfall from a
:class:`~anuga.rain.raster_time_slice_data.Raster_time_slice_data` instance
(for example a :class:`~anuga.rain.calibrated_radar_rain.Calibrated_radar_rain`
object loaded from BoM radar NetCDF files).

During ``__init__`` the operator extracts the raster values at **all** domain
centroids for **every** time slice and stores the result as a cache array.
This precomputation means each ``__call__`` involves only an index lookup
and, when the time slice advances, a single :meth:`~Rate_operator.set_rate`
call using a ``centroid_array`` rate type, which the GPU backend handles
natively.

.. code-block:: python

   import anuga
   from anuga.rain import Calibrated_radar_rain, Raster_rate_operator

   # Load BoM radar grids covering the simulation period
   radar = Calibrated_radar_rain(
       radar_dir='/data/bom_radar',
       start_time='20230101_0600',
       final_time='20230101_1800',
   )

   domain = anuga.create_domain_from_regions(...)

   # time_offset aligns raster epoch times with domain t=0
   rain = Raster_rate_operator(domain, radar,
                               time_offset=radar.start_time)

   for t in domain.evolve(yieldstep=rain.time_step, finaltime=43200.0):
       domain.print_timestepping_statistics()


ARR 2016 design-storm rainfall
--------------------------------

The :class:`ARR_rate_operator` combines:

* **Spatial depth** — a per-centroid array of total event depths (mm) derived
  from ARR IFD grid data (:class:`Arr_grd`).
* **Temporal pattern** — a :class:`Single_pattern` from the ARR 2016 temporal
  patterns zip archive.

Like :class:`Raster_rate_operator`, all rate arrays are precomputed at
initialisation.

Loading ARR data
~~~~~~~~~~~~~~~~

**Hub file** — download from the `ARR Data Hub`_ and parse with
:class:`Arr_hub_rain` to get the temporal pattern code for your site:

.. code-block:: python

   from anuga.rain import Arr_hub_rain

   hub = Arr_hub_rain('site_hub_export.txt')
   print(hub.Tpat_code)   # e.g. 'ECsouth' for East Coast South
   print(hub.Loc_Lat, hub.Loc_Lon)

**IFD grid** — download the IFD zip files from the `ARR IFD Data Server`_
(:doc:`../appendices/arr_design_rainfall` shows how) and open the grid for a
specific duration and AEP:

.. code-block:: python

   from anuga.rain import Arr_ifd_rain

   ifd = Arr_ifd_rain('/data/arr_ifd', Lat=-33.0, Lon=151.3)
   grd = ifd.open_grd(Dur=60, aep_percent=1.0)  # 60-min, 1% AEP

   # Extract depths at domain centroid coordinates
   # (centroids are in projected metres — convert to lat/lon first)
   import anuga
   from anuga.coordinate_transforms.redfearn import epsg_to_ll

   domain = anuga.create_domain_from_regions(...)
   xy = domain.get_centroid_coordinates(absolute=True)
   lats, lons = epsg_to_ll(xy[:, 0], xy[:, 1], domain.get_epsg())
   depth_array = grd.get_rain_at_points(lons, lats)  # mm

A point outside the grid raises ``ValueError``. A point on a NODATA cell gives
NaN and a warning, and :class:`ARR_rate_operator` rejects NaN on the cells it
acts on.

**Temporal patterns** — load the patterns zip and pick one pattern:

.. code-block:: python

   from anuga.rain import ARR_point_rainfall_patterns, Single_pattern

   prp = ARR_point_rainfall_patterns('EC_Patterns.zip', hub.Tpat_code)

   # Pattern index 1–720. The depth comes from depth_array: Ev_dep only scales
   # pattern.Rplot and pattern.plot(), and has no effect on ARR_rate_operator.
   pattern = Single_pattern(prp, index=1)

Running the simulation
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from anuga.rain import ARR_rate_operator

   # ifd_duration warns if the IFD grid and the pattern have different durations
   rain = ARR_rate_operator(domain, depth_array, pattern, ifd_duration=grd.dur)

   finaltime = pattern.Ev_dur * 60 + 6 * 3600   # the storm, then 6 hours to drain
   for t in domain.evolve(yieldstep=rain.time_step, finaltime=finaltime):
       domain.print_timestepping_statistics()

The operator switches rate arrays only at pattern time-step boundaries
(every ``pattern.Tstep`` minutes), so the GPU execution path sees a
``centroid_array`` rate with no CPU overhead between those transitions.
Outside the pattern, before its start and after its end, the rate is
``default_rate`` (0 by default). The same applies to
:class:`Raster_rate_operator` outside the time span of its data.

A timestep that starts before a switch applies the old rate for the whole step,
so end each ``yieldstep`` on a change for an exact total depth. ``rain.time_step``
is the time between changes in seconds; use it or a whole fraction of it.
See :doc:`../appendices/arr_design_rainfall` for more on choosing it.

ANUGA does not apply the areal reduction factors that :class:`Arr_hub_rain`
reads from the hub file: adjust ``depth_array`` yourself. Losses are applied
only if you ask for them, as described next.

Rainfall losses
---------------

By default all the rain runs off. To take out the part that soaks in, give a
loss model to :class:`Hyetograph_rate_operator`, :class:`Raster_rate_operator`
or :class:`ARR_rate_operator` with ``losses=``. The operator then adds the
*excess* rain only.

.. code-block:: python

   # US: NRCS (SCS) curve number
   losses = anuga.SCS_curve_number(curve_number=75.0)

   # Australia: ARR initial loss (mm) and continuing loss (mm/h)
   losses = anuga.Initial_continuing_loss(initial=15.0, continuing=2.5)

   # Initial loss (mm), then 30% of the rain lost (a runoff coefficient of 0.7)
   losses = anuga.Initial_proportional_loss(proportional=0.3, initial=5.0)

   rain = anuga.Raster_rate_operator(domain, rain_data, losses=losses)

   for t in domain.evolve(yieldstep=rain.time_step, finaltime=finaltime):
       print(rain.cumulative_influx, rain.cumulative_loss)   # m^3 added, m^3 lost

Each parameter can be a number, an array with one value per triangle, a
``Quantity``, or a function ``f(x, y)`` of the centroid coordinates, so the
losses can follow land use and soil:

.. code-block:: python

   def curve_number(x, y):
       return np.where(x < 500.0, 98.0, 70.0)    # paved west of x = 500 m

   losses = anuga.SCS_curve_number(curve_number=curve_number)

For a map made of polygons and rasters, build the function with
``composite_quantity_setting_function``, as for elevation or friction.
The first matching entry wins:

.. code-block:: python

   from anuga.utilities.quantity_setting_functions import \
       composite_quantity_setting_function

   curve_number = composite_quantity_setting_function(
       [[car_park_polygon, 98.0],          # a polygon and a value
        ['Extent', 'curve_number.tif'],    # a raster, over its own extent
        ['All', 70.0]],                    # everywhere else
       domain)

   losses = anuga.SCS_curve_number(curve_number=curve_number)

:class:`SCS_curve_number`
   The runoff from a cumulative rain depth :math:`P` is
   :math:`Q = (P - I_a)^2 / (P - I_a + S)` once :math:`P` passes the initial
   abstraction :math:`I_a = \lambda S`, with :math:`S = 25400/CN - 254` mm.
   Take ``curve_number`` from the NRCS tables (land use, hydrologic soil
   group, antecedent moisture). ``initial_abstraction_ratio`` is
   :math:`\lambda`, 0.2 by default. A curve number of 100 loses nothing.

:class:`Initial_continuing_loss`
   Rain first fills the initial loss. After that the continuing loss is taken
   from the rain at a steady rate, and the excess is never below zero. The
   ARR Data Hub gives both values for a location; :class:`Arr_hub_rain` reads
   them from the hub file as the strings ``ARR_IL`` and ``ARR_CL``.

:class:`Initial_proportional_loss`
   Rain first fills the initial loss (0 by default). After that the fraction
   ``proportional`` of the rain is lost. This is ARR's initial loss /
   proportional loss model, and with no initial loss it is a fixed runoff
   coefficient or percentage runoff of ``1 - proportional``.

Things to know:

* These are losses on rainfall. Water that runs on to a cell from upslope, or
  ponds on it, does not soak in.
* The rain is taken as uniform within each time slice of the data, so the
  excess is spread evenly over the slice. The total over the storm is exact
  for the curve number method.
* ``factor`` scales the excess, and ``default_rate`` (outside the time span of
  the data) is applied with no losses.
* A loss model holds the state of one storm, so make one for each operator.
* ``cumulative_influx`` is the water added (the excess), and
  ``cumulative_loss`` the volume lost.

API summary
-----------

.. autosummary::
   :nosignatures:

   Hyetograph_rate_operator
   Raster_rate_operator
   ARR_rate_operator
   Arr_hub_rain
   ARR_point_rainfall_patterns
   Single_pattern
   Arr_ifd_rain
   Arr_grd
   SCS_curve_number
   Initial_continuing_loss
   Initial_proportional_loss

.. _ARR Data Hub: https://data.arr-software.org/
.. _ARR IFD Data Server: https://www.bom.gov.au/water/designRainfalls/revised-ifd/

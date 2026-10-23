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

   for t in domain.evolve(yieldstep=360.0, finaltime=43200.0):
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
   print(hub.Tpat_code)   # e.g. 'EC' for East Coast
   print(hub.Loc_Lat, hub.Loc_Lon)

**IFD grid** — download the IFD zip files from the `ARR IFD Data Server`_
and open the grid for a specific duration and AEP:

.. code-block:: python

   from anuga.rain import Arr_ifd_rain

   ifd = Arr_ifd_rain('/data/arr_ifd', Lat=-33.0, Lon=151.3)
   grd = ifd.open_grd(Dur=60, Frq=100)  # 60-min, 1% AEP

   # Extract depths at domain centroid coordinates
   # (centroids are in UTM metres — convert to lon/lat first)
   import anuga
   domain = anuga.create_domain_from_regions(...)
   geo = domain.geo_reference
   lons, lats = geo.get_lonlat_from_xy(
       domain.centroid_coordinates[:, 0],
       domain.centroid_coordinates[:, 1])
   depth_array = grd.get_rain_at_points(lons, lats)  # mm

**Temporal patterns** — load the patterns zip and pick one pattern:

.. code-block:: python

   from anuga.rain import ARR_point_rainfall_patterns, Single_pattern

   prp = ARR_point_rainfall_patterns('EC_Patterns.zip', hub.Tpat_code)

   # Pattern index 1–720; Ev_dep scales the hyetograph to a target depth
   mean_depth = float(depth_array.mean())
   pattern = Single_pattern(prp, index=1, Ev_dep=mean_depth)

Running the simulation
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from anuga.rain import ARR_rate_operator

   rain = ARR_rate_operator(domain, depth_array, pattern)

   for t in domain.evolve(yieldstep=600.0, finaltime=3600.0):
       domain.print_timestepping_statistics()

The operator switches rate arrays only at pattern time-step boundaries
(every ``pattern.Tstep`` minutes), so the GPU execution path sees a
``centroid_array`` rate with no CPU overhead between those transitions.

API summary
-----------

.. autosummary::
   :nosignatures:

   Raster_rate_operator
   ARR_rate_operator
   Arr_hub_rain
   ARR_point_rainfall_patterns
   Single_pattern
   Arr_ifd_rain
   Arr_grd

.. _ARR Data Hub: https://data.arr-software.org/
.. _ARR IFD Data Server: https://www.bom.gov.au/water/designRainfalls/revised-ifd/

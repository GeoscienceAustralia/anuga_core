.. _arr_design_rainfall:

.. currentmodule:: anuga

ARR design rainfall data
========================

This appendix shows how to get the data for an Australian Rainfall and Runoff
(ARR 2016) design storm and load it for :class:`ARR_rate_operator`. See
:doc:`../setup_anuga_script/rainfall` for how the operator is used.

A design storm needs two things from two different places:

* **Rainfall depth at each location (IFD grids)**, from the Bureau of Meteorology's
  `Design Rainfall Data System`_. Read them with :class:`Arr_ifd_rain`, which
  returns an :class:`Arr_grd`.
* **How that depth falls over time (temporal patterns)**, from the `ARR Data Hub`_.
  Read them with :class:`Arr_hub_rain`, :class:`ARR_point_rainfall_patterns` and
  :class:`Single_pattern`.


IFD grids
---------

Download the grids in a browser. The Bureau's
`copyright notice <https://www.bom.gov.au/copyright>`_ does not allow automated
extraction from its site, so don't script this step, and check its terms before
you redistribute the grids, for example in a repository.

1. On the `Design Rainfall Data System`_, choose **Extent** and enter the
   lower-left and upper-right corners of a box around your catchment. The grid
   cells are 0.025° (about 2.5 km), so leave a margin of a cell or two.
2. Select **Submit**. The results panel stays empty.
3. Use the download buttons at the top of the page, and keep the page open until
   the zip file arrives.

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - File
     - AEPs (%)
   * - ``IFD_Extent.zip``
     - 63.212, 50, 39.347, 20, 18.127, 10, 5, 2, 1
   * - ``Rare_Extent.zip``
     - 1, 0.5, 0.2, 0.1, 0.05
   * - ``Very_frequent_Extent.zip``
     - 99.999, 99.752, 98.168, 95.021, 86.466, 63.212, 39.347, 18.127

Each probability comes in the 29 standard durations, from 1 minute to 7 days.

.. tip::

   Download the sets one at a time. The all-patterns button sends all 522 requests
   at once and builds the zip file only when they have all finished. A request
   that fails does not produce an error message: the counter just stops, and the
   grids already received are lost when you reload.

Put the zip files, under any name, in one directory:

.. code-block:: python

   from anuga.rain import Arr_ifd_rain

   ifd = Arr_ifd_rain('/data/lismore_ifd', Lat=-28.81, Lon=153.28)
   grd = ifd.open_grd(Dur=1440, aep_percent=1.0)   # 24 hours, 1% AEP

:meth:`Arr_ifd_rain.open_grd` raises ``ValueError``, listing what is available, if
no grid matches.

A grid holds the depth for one duration and one probability, so a set of grids gives
curves like these at each location.

.. figure:: arr_ifd_curves.*
   :width: 100%
   :alt: Invented rainfall depth against storm duration on log axes, for five AEPs
         from 50% to 0.1%. Depth rises with duration and with rarity. One point,
         24 hours at 1% AEP, is marked at 309 mm.

   Synthetic data, to show the form of the curves. These are not real depths. The
   marked point is the grid that ``open_grd(Dur=1440, aep_percent=1.0)`` would open.

.. warning::

   ``aep_percent`` is the clearer way to choose a grid. If you use ``Frq``, it is
   the label in the file name, which is the AEP in percent times 1000, not a return
   period: ``Frq=1000`` (``01000aep``) is 1% AEP and ``Frq=100`` (``00100aep``) is
   0.1% AEP, or 1 in 1000. Earlier versions of these docs gave ``Frq=100`` for 1%,
   which opens a much rarer storm without any warning if ``Rare_Extent.zip`` is
   present.

Then read the depth at each centroid. The grids are in latitude and longitude
(GDA94), so convert the mesh's coordinates first:

.. code-block:: python

   from anuga.coordinate_transforms.redfearn import epsg_to_ll

   xy = domain.get_centroid_coordinates(absolute=True)
   lats, lons = epsg_to_ll(xy[:, 0], xy[:, 1], domain.get_epsg())
   depth_array = grd.get_rain_at_points(lons, lats)   # mm

Each point takes the value of the grid cell holding it, so the depths are in
steps of about 2.5 km.

.. figure:: arr_ifd_grid.*
   :width: 80%
   :alt: A map of invented depths in square cells, with a black coast line and gray
         ocean cells holding 0. A point inland reads 309 mm. A point on land, next
         to the coast but in a cell whose centre is offshore, reads 0 mm.

   Synthetic grid and coast, read back with :meth:`Arr_grd.get_rain_at_point`. These
   are not real depths. Gray cells hold 0. The point at the coast is on land, but its
   cell is mostly ocean, so it reads 0 mm.

* **Outside the grid:** a point outside the extent raises ``ValueError``.
* **Ocean:** the grids hold 0, not NODATA, over the ocean. A point on such a cell
  gets no rain, and gives a warning. Check this on a coastal model, where a cell
  near the shore could read 0 while the centroid in it is on land.
* **NODATA:** a NODATA cell gives NaN and a warning, and
  :class:`ARR_rate_operator` raises an error if NaN falls on a cell it acts on.


Temporal patterns
-----------------

Enter your catchment on the `ARR Data Hub`_, tick the data you want, including
**Temporal Patterns**, and submit. Save the results with **Download TXT**, and
download the patterns from the **Download (.zip)** link beside them.

.. code-block:: python

   from anuga.rain import Arr_hub_rain, ARR_point_rainfall_patterns, Single_pattern

   hub = Arr_hub_rain('lismore_hub_export.txt')
   prp = ARR_point_rainfall_patterns('Point_ECsouth.zip', hub.Tpat_code)

Each zip file holds 720 patterns. Each of the 24 durations, from 10 minutes to
7 days, has ten patterns in each of three groups: ``frequent``, ``intermediate``
and ``rare``.
``prp.linesInc`` has one row per pattern, after a header row, with the duration in
column 1 and the group in column 4. This picks the ten for a rare 24-hour storm:

.. code-block:: python

   rows = [i for i, row in enumerate(prp.linesInc)
           if i > 0 and int(row[1]) == 1440 and row[4].strip() == 'rare']
   pattern = Single_pattern(prp, index=rows[0])

ARR 2016 runs a storm with each of the ten patterns and works from the set of
results, so a single pattern is not a design answer on its own.

The ten patterns can put the same depth down very differently. Among the ten rare
24-hour patterns for East Coast (South), the wettest hour holds from 7.6% to 17.4% of
the depth, and the first 12 hours from 31.5% to 66.3%.

.. figure:: arr_pattern_ensemble.*
   :width: 100%
   :alt: Ten small bar charts of hourly rainfall as a percentage of the storm depth,
         one per rare 24-hour pattern, sharing axes. The bursts fall at different
         hours in each pattern. Pattern 4661 has the highest hour, at 17.4%.

   The ten rare 24-hour patterns for East Coast (South), labelled by event ID.
   Patterns from the `ARR Data Hub`_.

.. figure:: arr_pattern_cumulative.*
   :width: 100%
   :alt: Cumulative rainfall as a percentage of the storm depth against time for the
         ten patterns. All reach 100% at 24 hours, and they spread from about 32% to
         66% at 12 hours. Pattern 4655, which the example selects, is highlighted.

   The same patterns as cumulative rain. Each delivers the whole depth and differs in
   how it gets there, so the runoff can differ too. Patterns from the `ARR Data Hub`_.


Running the storm
-----------------

.. code-block:: python

   from anuga.rain import ARR_rate_operator

   rain = ARR_rate_operator(domain, depth_array, pattern, ifd_duration=grd.dur)

   finaltime = pattern.Ev_dur * 60 + 6 * 3600   # the storm, then 6 hours to drain
   for t in domain.evolve(yieldstep=rain.time_step, finaltime=finaltime):
       domain.print_timestepping_statistics()

``ifd_duration`` warns if the grid and the pattern have different durations.

``rain.time_step`` is the time between changes of rate, in seconds. It is 3600 for
the 24-hour pattern above, and the patterns range from 5 minutes to 3 hours. A
timestep that starts before a change applies the old rate for the whole step, so the
storm delivers its exact depth only if each ``yieldstep`` ends on a change.

* **More output:** use a whole fraction of it, such as
  ``yieldstep=rain.time_step / 6`` for every 10 minutes here. A fixed value is not
  safe: 600 s fits this pattern, but not one with 15-minute steps.
* **Less output:** keep ``yieldstep`` and add an ``outputstep``, a whole multiple of
  ``yieldstep``, which sets how often the SWW file is written.
* **Starting before the rain:** the pattern starts at model time 0, and yields count
  from the start time. Use a ``starttime`` that is a whole number of yieldsteps
  before 0, for example ``domain.set_starttime(-6 * rain.time_step)``.


What ANUGA does not do for you
------------------------------

ANUGA applies the depths and pattern as they are. These are yours to handle:

* **Climate baseline.** The Bureau's grids are its 2016 values. The ARR Data Hub
  also tabulates the same IFDs adjusted for climate change, and by default shows a
  current climate (2030) table. At Lismore, the 24-hour 1% AEP depth is 353 mm in
  the native table and 388 mm in the 2030 table. To use an adjusted depth, scale
  ``depth_array``.
* **Areal reduction.** A depth from the grids is a point rainfall. The hub's ARF
  parameters are read by :class:`Arr_hub_rain` as strings and not applied.
* **Losses.** The hub's initial and continuing losses are read and not applied.
  :class:`Arr_hub_rain` reads the ``[LOSSES]`` section, which the hub now labels as
  superseded where it has newer values, so check the hub's results for the losses
  it recommends.
* **Jurisdiction guidance.** Some states have their own guidance on losses and
  pre-burst rainfall. The ARR Data Hub has jurisdiction pages, including NSW, and
  this can matter a lot for rare storms.
* **Areal temporal patterns.** ARR 2016 has separate patterns for larger catchments.
  The steps above use the point patterns.

.. _Design Rainfall Data System: https://www.bom.gov.au/water/designRainfalls/revised-ifd/
.. _ARR Data Hub: https://data.arr-software.org/

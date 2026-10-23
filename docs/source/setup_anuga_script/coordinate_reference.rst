.. currentmodule:: anuga

.. _coordinate_reference:

Coordinate Reference Systems
=============================

ANUGA works in a **projected coordinate system** — all domain coordinates are
in metres, not degrees.  A :class:`Geo_reference` object records which
coordinate system is in use and the local origin offset (``xllcorner``,
``yllcorner``) that keeps the mesh coordinates numerically small.

Every SWW output file stores the complete geo-reference so downstream
tools (QGIS, Python post-processing, etc.) can interpret the coordinates
correctly.

---

Why a local origin?
-------------------

Raw UTM or national-grid coordinates for a real site are typically large
numbers (e.g. easting 308 500 m, northing 6 189 000 m for Cairns,
Australia).  Storing every mesh vertex relative to a local origin keeps
the values small and avoids floating-point precision loss during
computation.

Coordinates stored in a SWW file are **relative** to the geo-reference
origin.  To recover absolute coordinates, add ``xllcorner`` / ``yllcorner``
back to the mesh x / y arrays.

---

WGS84 UTM (the common case)
----------------------------

Most real-world ANUGA simulations use a WGS84 UTM zone.  Pass *zone* and
*hemisphere* (and optionally the local origin) to :class:`Geo_reference`:

.. code-block:: python

    import anuga

    # Cairns, Australia — UTM zone 55, southern hemisphere
    geo_ref = anuga.Geo_reference(
        zone       = 55,
        xllcorner  = 363000.0,   # easting of domain origin (m)
        yllcorner  = 8132000.0,  # northing of domain origin (m)
        hemisphere = 'southern',
    )

    domain = anuga.create_domain_from_regions(
        bounding_polygon,
        boundary_tags,
        geo_reference = geo_ref,
    )

The EPSG code is **computed automatically** from the zone and hemisphere and
written to the SWW file:

* Northern hemisphere: ``EPSG = 32600 + zone``  (e.g. zone 31N → EPSG 32631)
* Southern hemisphere: ``EPSG = 32700 + zone``  (e.g. zone 55S → EPSG 32755)

.. code-block:: python

    >>> geo_ref.epsg
    32755

Alternatively, supply the EPSG code directly and let ANUGA infer the zone
and hemisphere:

.. code-block:: python

    geo_ref = anuga.Geo_reference(
        epsg      = 32755,
        xllcorner = 363000.0,
        yllcorner = 8132000.0,
    )
    >>> geo_ref.zone
    55
    >>> geo_ref.hemisphere
    'southern'

---

National and non-UTM coordinate systems
-----------------------------------------

Some countries use national projected CRS that do not follow the UTM zone
structure.  Common examples:

===============================  ==========  ==========================
Country / region                 EPSG code   CRS name
===============================  ==========  ==========================
Netherlands                      28992       Amersfoort / RD New
Great Britain                    27700       OSGB36 / British National Grid
Germany (ETRS89 / UTM zone 32)   25832       ETRS89 / UTM zone 32N
GDA2020 MGA zone 56 (Australia)  7856        GDA2020 / MGA zone 56
WGS84 geographic (lat/lon)       4326        WGS84
===============================  ==========  ==========================

These work fine with ANUGA — pass the EPSG code directly.  If ``pyproj``
is installed, the datum, projection name, and false easting/northing are
populated automatically from the EPSG definition:

.. code-block:: python

    # Dutch simulation in Rijksdriehoekstelsel (RD New)
    geo_ref = anuga.Geo_reference(
        epsg      = 28992,
        xllcorner = 120000.0,   # easting relative to RD New origin (m)
        yllcorner = 480000.0,   # northing relative to RD New origin (m)
    )

    >>> geo_ref.epsg
    28992
    >>> geo_ref.projection
    'Amersfoort / RD New'
    >>> geo_ref.datum
    'Amersfoort'
    >>> geo_ref.false_easting, geo_ref.false_northing
    (155000, 463000)
    >>> repr(geo_ref)
    '(crs=Amersfoort / RD New, easting=120000.000000, northing=480000.000000, epsg=28992)'
    >>> geo_ref.is_located()
    True

.. note::

    ANUGA does **not** perform any reprojection — it treats all coordinates
    as a flat Cartesian grid in metres.  Convert your input data to the
    target projected CRS *before* building the domain.  Tools such as
    ``pyproj``, GDAL, or QGIS can reproject shapefiles, rasters, and point
    clouds as needed.

---

Wavetank and hypothetical simulations
---------------------------------------

When no real geographic location is involved (e.g. a laboratory wave-tank
or an idealised test case) simply omit zone and EPSG:

.. code-block:: python

    geo_ref = anuga.Geo_reference()   # zone=-1, no EPSG
    geo_ref.is_located()              # False

The domain coordinates are then treated as an arbitrary Cartesian system
with no connection to any geographic CRS.

---

Checking what is stored in a SWW file
---------------------------------------

The geo-reference metadata can be read directly from a SWW file:

.. code-block:: python

    from anuga.file.netcdf import NetCDFFile
    from anuga.coordinate_transforms.geo_reference import Geo_reference
    from anuga.config import netcdf_mode_r

    fid = NetCDFFile('my_simulation.sww', netcdf_mode_r)
    geo_ref = Geo_reference(NetCDFObject=fid)
    fid.close()

    print(geo_ref)          # e.g. (zone=55, ..., epsg=32755)
    print(geo_ref.epsg)     # EPSG code stored in file
    print(geo_ref.is_located())

---

.. seealso::

   :ref:`Geo_reference API <api_geo_reference>` in the API Reference.

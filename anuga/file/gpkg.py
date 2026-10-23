"""GeoPackage (GPKG) in and out: polygons for ANUGA regions, buildings and
breaklines, editable in QGIS or any GIS, and the two CSV conventions ANUGA
already reads (anuga-community/anuga_core#39, #40).

Needs ``fiona`` and ``shapely`` (the ``anuga[data]`` extra). They are
imported when a function is called, so the rest of ANUGA does not depend on
them.

Two things to know about the conventions:

* An ANUGA polygon is a list of ``[x, y]`` points, NOT closed (the last point
  is not a repeat of the first). GIS rings are closed. The functions here
  close on the way out and open on the way in.
* Attributes travel with the geometry. ``polygons2gpkg`` takes one dict per
  polygon and ``gpkg2polygons`` returns one; the schema is taken from the
  first dict (int, float or text).
"""
import os
import glob as _glob

import numpy as num


def _gis():
    try:
        import fiona
        from shapely.geometry import shape, mapping, Polygon, LineString
    except ImportError as exc:
        raise ImportError(
            'GeoPackage support needs fiona and shapely: pip install '
            '"anuga[data]" (or conda install fiona shapely)') from exc
    return fiona, shape, mapping, Polygon, LineString


def _crs(crs):
    """fiona accepts an EPSG code, a 'EPSG:nnnn' string, WKT or None."""
    if crs is None:
        return None
    if isinstance(crs, int):
        return 'EPSG:%d' % crs
    return crs


def _schema_properties(attributes):
    """A fiona properties schema from the first attribute dict."""
    if not attributes or not attributes[0]:
        return {}
    props = {}
    for name, value in attributes[0].items():
        if isinstance(value, (bool, num.bool_)):
            props[name] = 'int'
        elif isinstance(value, (int, num.integer)):
            props[name] = 'int'
        elif isinstance(value, (float, num.floating)):
            props[name] = 'float'
        else:
            props[name] = 'str'
    return props


def _native(value):
    """numpy scalars to Python scalars, for fiona."""
    if isinstance(value, num.generic):
        return value.item()
    return value


def polygons2gpkg(polygons, filename, attributes=None, crs=None,
                  layer='polygons', geometry_type='Polygon'):
    """Write ANUGA polygons (or polylines) to a GeoPackage layer.

    Parameters
    ----------
    polygons : sequence of (N, 2) array-like
        One polygon (or polyline) per entry, as ANUGA holds them: a list of
        ``[x, y]`` points, not closed. A closed ring is accepted too.
    filename : str
        The ``.gpkg`` to write. An existing file gains (or replaces) the layer.
    attributes : sequence of dict, optional
        One dict per polygon; the keys become the layer's fields. The types
        are taken from the first dict.
    crs : int or str, optional
        The coordinate reference system, as an EPSG code (``32756``), a
        string (``'EPSG:32756'``) or WKT. Absent means unknown, which QGIS
        will ask about.
    layer : str
        The layer name inside the GeoPackage.
    geometry_type : {'Polygon', 'LineString'}
        Write the entries as polygons (default) or as polylines.

    Examples
    --------
    >>> anuga.polygons2gpkg([house1, house2], 'houses.gpkg',
    ...                     attributes=[{'id': 1, 'floors': 2},
    ...                                 {'id': 2, 'floors': 1}],
    ...                     crs=32756)
    """
    fiona, shape, mapping, Polygon, LineString = _gis()
    polygons = [num.asarray(p, dtype=float) for p in polygons]
    n = len(polygons)
    if attributes is None:
        attributes = [{} for _ in range(n)]
    attributes = list(attributes)
    if len(attributes) != n:
        raise ValueError('polygons2gpkg: %d polygons but %d attribute dicts'
                         % (n, len(attributes)))
    if geometry_type not in ('Polygon', 'LineString'):
        raise ValueError("geometry_type must be 'Polygon' or 'LineString'")

    schema = {'geometry': geometry_type,
              'properties': _schema_properties(attributes)}
    # Writing a named layer to an existing GeoPackage adds that layer and
    # keeps the others ('a' would append features to an existing layer, and
    # is not what "write this layer" means). An existing layer of the same
    # name is replaced.
    if os.path.exists(filename):
        try:
            existing = fiona.listlayers(filename)
        except Exception:
            existing = []
        if layer in existing:
            if len(existing) == 1:
                os.remove(filename)
            else:
                fiona.remove(filename, driver='GPKG', layer=layer)

    with fiona.open(filename, 'w', driver='GPKG', layer=layer,
                    schema=schema, crs=_crs(crs)) as dst:
        for poly, attr in zip(polygons, attributes):
            if poly.ndim != 2 or poly.shape[1] != 2 or len(poly) < 2:
                raise ValueError('each polygon must be an (N, 2) array of '
                                 'points with N >= 2, got shape %s'
                                 % (poly.shape,))
            pts = [(float(x), float(y)) for x, y in poly]
            if geometry_type == 'Polygon':
                geom = Polygon(pts)          # shapely closes the ring
            else:
                geom = LineString(pts)
            props = {k: _native(attr.get(k)) for k in schema['properties']}
            dst.write({'geometry': mapping(geom), 'properties': props})


def gpkg2polygons(filename, layer=None, closed=False):
    """Read a GeoPackage layer as ANUGA polygons with their attributes.

    Parameters
    ----------
    filename : str
        The ``.gpkg`` to read.
    layer : str, optional
        The layer to read; the first layer when not given.
    closed : bool
        Return closed rings (last point repeats the first). Default False,
        ANUGA's convention.

    Returns
    -------
    polygons : list of (N, 2) ndarray
        One per feature; a MultiPolygon or MultiLineString contributes one
        entry per part, with its attributes repeated. Only exterior rings
        are returned: holes are dropped.
    attributes : list of dict
        The feature properties, one dict per returned polygon.

    Examples
    --------
    >>> houses, attrs = anuga.gpkg2polygons('houses.gpkg')
    >>> for poly, a in zip(houses, attrs):
    ...     domain.set_quantity('elevation', a['floors'] * 3.0, polygon=poly)
    """
    fiona, shape, mapping, Polygon, LineString = _gis()
    if layer is None:
        layers = fiona.listlayers(filename)
        if not layers:
            raise ValueError('%s has no layers' % filename)
        layer = layers[0]

    polygons = []
    attributes = []

    def add(coords, props):
        pts = num.asarray(coords, dtype=float)[:, :2]
        if not closed and len(pts) > 1 and num.allclose(pts[0], pts[-1]):
            pts = pts[:-1]
        polygons.append(pts)
        attributes.append(dict(props))

    with fiona.open(filename, 'r', layer=layer) as src:
        for feature in src:
            geom = feature['geometry']
            if geom is None:
                continue
            props = feature['properties']
            g = shape(geom)
            parts = list(g.geoms) if hasattr(g, 'geoms') else [g]
            for part in parts:
                if part.geom_type == 'Polygon':
                    add(part.exterior.coords, props)
                elif part.geom_type == 'LineString':
                    add(part.coords, props)
                else:
                    raise ValueError('gpkg2polygons: unsupported geometry %s '
                                     'in layer %s' % (part.geom_type, layer))
    return polygons, attributes


def gpkg_layers(filename):
    """The layer names in a GeoPackage."""
    fiona = _gis()[0]
    return list(fiona.listlayers(filename))


# ---------------------------------------------------------------------------
# The CSV conventions ANUGA already reads
# ---------------------------------------------------------------------------

def polygon_csv_files2gpkg(csv_files, filename, crs=None, layer='polygons'):
    """One-polygon-per-file CSVs (as ``anuga.read_polygon`` reads them) to a
    GeoPackage layer, one feature per file with a ``name`` field holding the
    file's stem.

    ``csv_files`` may be a directory (all its ``*.csv``), a glob pattern, or
    a list of paths. Each file holds ``x,y`` rows and no header, e.g. the
    Merewether case study's ``houses/house000.csv``.
    """
    from anuga.geometry.polygon import read_polygon
    if isinstance(csv_files, str):
        if os.path.isdir(csv_files):
            files = sorted(_glob.glob(os.path.join(csv_files, '*.csv')))
        else:
            files = sorted(_glob.glob(csv_files))
    else:
        files = list(csv_files)
    if not files:
        raise ValueError('polygon_csv_files2gpkg: no CSV files in %r' % (csv_files,))
    polygons = [read_polygon(f) for f in files]
    attributes = [{'name': os.path.splitext(os.path.basename(f))[0]} for f in files]
    polygons2gpkg(polygons, filename, attributes=attributes, crs=crs, layer=layer)
    return files


def gpkg2polygon_csv_files(filename, directory, layer=None, prefix=None,
                           name_field='name'):
    """A GeoPackage layer to one-polygon-per-file CSVs that
    ``anuga.read_polygon`` reads (``x,y`` rows, no header, not closed).

    Files are named from the ``name`` field when the layer has one, else
    ``<prefix><index>.csv`` with the index zero-padded. Returns the paths.
    """
    from anuga.geometry.polygon import write_polygon
    polygons, attributes = gpkg2polygons(filename, layer=layer)
    os.makedirs(directory, exist_ok=True)
    if prefix is None:
        prefix = (layer or gpkg_layers(filename)[0]) + '_'
    width = max(3, len(str(max(len(polygons) - 1, 0))))
    paths = []
    for i, (poly, attr) in enumerate(zip(polygons, attributes)):
        name = attr.get(name_field)
        stem = str(name) if name not in (None, '') else '%s%0*d' % (prefix, width, i)
        path = os.path.join(directory, stem + '.csv')
        write_polygon(poly, path)
        paths.append(path)
    return paths


def building_csv2gpkg(csv_file, filename, crs=None, layer='buildings',
                      value_name='floors'):
    """The ``easting,northing,id,<value>`` CSV that
    ``anuga.load_csv_as_polygons`` / ``load_csv_as_building_polygons`` read,
    to a GeoPackage layer with ``id`` and ``<value>`` fields.

    The ``id`` and value columns are stored as text and number: ``id`` as an
    integer when every id parses as one, otherwise text; the value as a float.
    """
    from anuga.file.csv_file import load_csv_as_polygons
    polygons, values = load_csv_as_polygons(csv_file, value_name=value_name)
    ids = list(polygons.keys())
    try:
        int_ids = [int(i) for i in ids]
        id_values = int_ids
    except ValueError:
        id_values = [str(i) for i in ids]
    attributes = []
    for key, id_value in zip(ids, id_values):
        attr = {'id': id_value}
        if values is not None:
            attr[value_name] = float(values[key])
        attributes.append(attr)
    polygons2gpkg([polygons[k] for k in ids], filename, attributes=attributes,
                  crs=crs, layer=layer)


def gpkg2building_csv(filename, csv_file, layer=None, id_field='id',
                      value_name='floors', value_field=None):
    """A GeoPackage polygon layer to the ``easting,northing,id,<value>`` CSV
    that ``anuga.load_csv_as_polygons`` / ``load_csv_as_building_polygons``
    read. Rings are written closed (the first point repeated last), as in
    that format's documentation.

    ``id`` comes from ``id_field`` when the layer has it, else from the
    feature index; the value column from ``value_field`` (default: the same
    name as ``value_name``), 0 when absent.
    """
    if value_field is None:
        value_field = value_name
    polygons, attributes = gpkg2polygons(filename, layer=layer, closed=True)
    with open(csv_file, 'w') as f:
        f.write('easting,northing,%s,%s\n' % ('id', value_name))
        for i, (poly, attr) in enumerate(zip(polygons, attributes)):
            pid = attr.get(id_field)
            if pid in (None, ''):
                pid = i
            value = attr.get(value_field)
            if value in (None, ''):
                value = 0
            for x, y in poly:
                f.write('%s,%s,%s,%s\n' % (repr(float(x)), repr(float(y)), pid, value))

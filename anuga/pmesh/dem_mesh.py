"""Meshes sized from a digital elevation model.

Where water goes is written in the DEM. Routing it by D8 (each cell drains
to its steepest downhill neighbour, after the depressions are filled) gives
every cell the area that drains through it, and so the flow paths: streams
in a catchment, streets and overland flow paths in a town. The height of a
cell above the stream it drains to (HAND, height above nearest drainage)
then gives the valley floors those streams spread over when they rise.

:func:`size_field_from_dem` turns both into a :class:`Size_field`: fine
triangles along the flow paths and over the floors, coarse ones on the
slopes between them, graded smoothly, and with no segment added to the
mesh. :func:`create_pmesh_from_dem` meshes with it, and
``anuga.create_domain_from_dem`` also sets the elevation from the DEM.

Averaging a DEM over coarse triangles can fill a gully a cell drained by,
leaving the cell upstream lower than all its neighbours: a pit that holds
water for ever. :func:`remove_mesh_pits` removes such pits on the mesh
itself, by filling or by breaching them.

A DEM is given as a file name (``.asc``, ``.grd``, ``.dem`` or ``.tif``;
``.tif`` needs rasterio) or as a raster tuple ``(x, y, Z)`` as read by
``grd2array``, ``dem2array`` or ``tif2array`` in ``anuga.file_conversion``:
the cell
coordinates ``x`` and ``y`` (absolute, evenly spaced, the same spacing in
both directions) and ``Z[i, j]``, the elevation at ``(x[i], y[j])``, NaN
where there is no data.
"""

import os

import numpy as num

from anuga.pmesh.size_field import (Size_field, _rasterise_polygon,
                                    _refine_edges, _refine_mask,
                                    area_to_edge, edge_to_area)
from anuga.utilities.numerical_tools import ensure_numeric


def read_dem(dem):
    """A DEM as a raster tuple ``(x, y, Z)``; see the module docstring."""
    if isinstance(dem, (str, os.PathLike)):
        ext = os.path.splitext(str(dem))[1].lower()
        if ext in ('.asc', '.grd'):
            from anuga.file_conversion.grd2array import grd2array
            x, y, Z = grd2array(str(dem))
        elif ext == '.dem':
            from anuga.file_conversion.dem2array import dem2array
            x, y, Z = dem2array(str(dem))
        elif ext in ('.tif', '.tiff'):
            from anuga.file_conversion.tif2array import tif2array
            x, y, Z = tif2array(str(dem))
        else:
            raise ValueError('DEM file must be .asc, .grd, .dem or .tif, '
                             'not %s' % dem)
    else:
        x, y, Z = dem
    x = num.asarray(x, float)
    y = num.asarray(y, float)
    Z = num.asarray(Z, float)
    if Z.shape != (len(x), len(y)):
        raise ValueError('DEM Z must have shape (len(x), len(y)) = (%d, %d), '
                         'not %s' % (len(x), len(y), Z.shape))
    return x, y, Z


class _Grid:
    """A DEM as a grid ``z[row, col]``, row 0 to the south, square cells."""

    def __init__(self, dem):
        x, y, Z = read_dem(dem)
        if len(x) < 2 or len(y) < 2:
            raise ValueError('DEM must have at least 2 cells each way')
        dx = num.diff(x)
        dy = num.diff(y)
        cellsize = dx[0]
        if (not num.allclose(dx, cellsize, rtol=1e-6)
                or not num.allclose(dy, cellsize, rtol=1e-6)):
            raise ValueError('DEM cells must be square and evenly spaced')
        self.x = x
        self.y = y
        self.cellsize = float(cellsize)
        # lower left corner of the lower left cell
        self.xll = float(x[0]) - 0.5 * self.cellsize
        self.yll = float(y[0]) - 0.5 * self.cellsize
        self.z = num.ascontiguousarray(Z.T)

    def mask(self, polygon):
        """Cells whose centres are inside ``polygon``."""
        nrow, ncol = self.z.shape
        m = num.zeros((nrow, ncol), bool)
        r0, c0, inside = _rasterise_polygon(polygon, self.xll, self.yll,
                                            self.cellsize, nrow, ncol)
        m[r0:r0 + inside.shape[0], c0:c0 + inside.shape[1]] = inside
        return m

    def extent(self):
        """The DEM's cell centres' bounding rectangle, anticlockwise from
        the lower left: segments 0 bottom, 1 right, 2 top, 3 left."""
        x0, x1 = float(self.x[0]), float(self.x[-1])
        y0, y1 = float(self.y[0]), float(self.y[-1])
        return [[x0, y0], [x1, y0], [x1, y1], [x0, y1]]


def _absolute(polygon, poly_geo_reference):
    polygon = ensure_numeric(polygon, float)
    if poly_geo_reference is not None:
        polygon = ensure_numeric(poly_geo_reference.get_absolute(polygon),
                                 float)
    return polygon


class _Flow:
    """D8 routing of a DEM within a bounding polygon, around holes."""

    def __init__(self, grid, bounding_polygon=None, interior_holes=None,
                 epsilon=1e-3):
        from anuga.pmesh.dem_flow_ext import route

        z = grid.z.copy()
        if bounding_polygon is not None:
            # Water leaves across the domain's boundary
            z[~grid.mask(bounding_polygon)] = num.nan
        for hole in interior_holes or []:
            # but not into a hole
            z[grid.mask(hole)] = num.inf
        self.grid = grid
        self.z = z
        self.receiver, self.order, self.cells = route(z, epsilon)

    @property
    def drainage_area(self):
        """Area draining through each cell, ``[row, col]``."""
        return self.cells * self.grid.cellsize ** 2

    def height_above_drainage(self, drainage_area):
        """HAND, where a stream drains at least ``drainage_area``."""
        from anuga.pmesh.dem_flow_ext import height_above_drainage

        return height_above_drainage(self.z, self.receiver, self.order,
                                     self.cells,
                                     drainage_area / self.grid.cellsize ** 2)


def drainage_area_from_dem(dem, bounding_polygon=None, interior_holes=None,
                           poly_geo_reference=None):
    """Area draining through each cell of a DEM, by D8 routing.

    Depressions are filled first, so that every cell drains. Water leaves
    the DEM across its edge, at nodata and outside ``bounding_polygon``, and
    flows around ``interior_holes``.

    Returns a raster tuple ``(x, y, A)``, with ``A[i, j]`` the area (in the
    DEM's units squared) draining through the cell at ``(x[i], y[j])``,
    itself included, and 0 where nothing is routed.
    """
    grid = _Grid(dem)
    flow = _Flow(grid,
                 None if bounding_polygon is None
                 else _absolute(bounding_polygon, poly_geo_reference),
                 [_absolute(h, poly_geo_reference)
                  for h in interior_holes or []])
    return grid.x, grid.y, flow.drainage_area.T.copy()


def size_field_from_dem(dem,
                        maximum_triangle_area,
                        streams=None,
                        floodplain=None,
                        regions=None,
                        bounding_polygon=None,
                        interior_holes=None,
                        grade=0.5,
                        band=1.5,
                        poly_geo_reference=None,
                        verbose=False):
    """Build a graded :class:`Size_field` from the flow paths of a DEM.

    The DEM is routed by D8 within ``bounding_polygon`` and around
    ``interior_holes``, and the mesh is made fine where water collects:

    - ``streams``: cells that at least ``drainage_area`` drains through get
      triangles of at most ``area``, held for ``band`` triangle edges either
      side, so that the path is a few triangles wide.
    - ``floodplain``: cells no more than ``height`` above the stream they
      drain to, where that stream drains at least ``drainage_area``, get
      triangles of at most ``area``.

    Elsewhere the target is ``maximum_triangle_area``. As in
    :func:`size_field_from_regions`, the target edge grows away from each
    class at ``grade`` per unit distance, and the finest target wins.

    Parameters
    ----------
    dem : str or (x, y, Z)
        The DEM; see the module docstring. Its cells set the raster of the
        size field.
    maximum_triangle_area : float
        Target area away from all streams, floodplains and regions.
    streams : list of (drainage_area, area), optional
        Drainage areas are in the DEM's units squared (m^2 for a DEM in
        metres). Big streams usually want small triangles, so the list
        typically runs from large drainage areas and small triangles to
        small drainage areas and larger ones.
    floodplain : list of (height, drainage_area, area), optional
        Height above the stream, in the DEM's elevation units.
    regions : list of (polygon, area), optional
        Polygons to refine as well, as for
        :func:`size_field_from_regions`: none of their segments enter the
        mesh. Use these for what the DEM cannot show, such as a town or a
        structure.
    bounding_polygon : list of points, optional
        Water leaves the DEM across this polygon. Default: the DEM's edge.
    interior_holes : list of polygons, optional
        Water flows around these (buildings, say), and they are not refined
        themselves.
    grade : float, optional
        Growth of the target edge per unit distance. Default 0.5.
    band : float, optional
        A stream's target holds for this many of its target edge lengths
        either side of the flow path. Default 1.5.
    poly_geo_reference : Geo_reference, optional
        Geo reference of the polygons. If None they are absolute, like the
        DEM.
    verbose : bool, optional
        Print the length of each stream class and the area of each
        floodplain class.

    Returns
    -------
    Size_field
        It also has the D8 drainage area as ``drainage_area``, a raster
        tuple as returned by :func:`drainage_area_from_dem`.
    """
    if not maximum_triangle_area > 0.0:
        raise ValueError('maximum_triangle_area must be positive')
    if grade < 0.0:
        raise ValueError('grade must not be negative')
    streams = sorted(streams or [], reverse=True)
    floodplain = list(floodplain or [])
    for _, area in streams:
        if not area > 0.0:
            raise ValueError('stream triangle areas must be positive')
    for _, _, area in floodplain:
        if not area > 0.0:
            raise ValueError('floodplain triangle areas must be positive')

    grid = _Grid(dem)
    if bounding_polygon is not None:
        bounding_polygon = _absolute(bounding_polygon, poly_geo_reference)
    holes = [_absolute(h, poly_geo_reference) for h in interior_holes or []]
    flow = _Flow(grid, bounding_polygon, holes)
    cs = grid.cellsize
    drained = flow.drainage_area

    h_max = float(area_to_edge(maximum_triangle_area))
    h = num.full(grid.z.shape, h_max)
    for drainage, area in streams:
        edge = min(float(area_to_edge(area)), h_max)
        mask = drained >= drainage
        _refine_mask(h, mask, edge, grade, cs, band=band * edge)
        if verbose:
            print('    streams draining >= %g -> area %g: %.0f long'
                  % (drainage, area, mask.sum() * cs))
    for height, drainage, area in floodplain:
        edge = min(float(area_to_edge(area)), h_max)
        mask = flow.height_above_drainage(drainage) <= height
        _refine_mask(h, mask, edge, grade, cs)
        if verbose:
            print('    within %g above a stream draining >= %g -> area %g: '
                  '%.0f in area' % (height, drainage, area,
                                    mask.sum() * cs * cs))
    if regions:
        _refine_edges(h, grid.xll, grid.yll, cs,
                      [(_absolute(p, poly_geo_reference),
                        min(float(area_to_edge(a)), h_max))
                       for p, a in regions], grade)

    field = Size_field(edge_to_area(h), grid.xll, grid.yll, cs)
    field.drainage_area = (grid.x, grid.y, drained.T.copy())
    return field


def create_pmesh_from_dem(dem,
                          maximum_triangle_area,
                          bounding_polygon=None,
                          boundary_tags=None,
                          streams=None,
                          floodplain=None,
                          regions=None,
                          interior_holes=None,
                          grade=0.5,
                          band=1.5,
                          size_field=None,
                          relax_passes=30,
                          poly_geo_reference=None,
                          verbose=False,
                          **kwargs):
    """Create a Pmesh sized from the flow paths of a DEM.

    Builds a size field with :func:`size_field_from_dem` and meshes with
    :func:`create_pmesh_from_regions`. The elevation is not set: see
    ``anuga.create_domain_from_dem``.

    Parameters
    ----------
    dem : str or (x, y, Z)
        The DEM; see the module docstring.
    maximum_triangle_area : float
        Target area away from all streams, floodplains and regions.
    bounding_polygon : list of points, optional
        Outline of the mesh. Default: the rectangle of the DEM's cell
        centres, with ``boundary_tags`` (if not given) ``bottom``,
        ``right``, ``top`` and ``left``.
    boundary_tags : dict, optional
        As for :func:`create_pmesh_from_regions`.
    streams, floodplain, regions, grade, band :
        As for :func:`size_field_from_dem`.
    interior_holes : list of polygons, optional
        Holes in the mesh (buildings, say); water is routed around them.
    size_field : Size_field or callable, optional
        A further size field; the finer of it and the DEM's wins.
    relax_passes : int, optional
        Passes of worst-cell relaxation after meshing. Default 30.
    poly_geo_reference : Geo_reference, optional
        Geo reference of the polygons. If None they are absolute, like the
        DEM.
    verbose : bool, optional
    **kwargs :
        Passed on to :func:`create_pmesh_from_regions`: for example
        ``breaklines``, ``hole_tags``, ``interior_regions``,
        ``minimum_triangle_angle``, ``mesh_geo_reference`` or ``filename``.

    Returns
    -------
    Pmesh
        Its size field is ``pmesh.size_field``, with the D8 drainage area as
        ``pmesh.size_field.drainage_area``.
    """
    from anuga.pmesh.mesh_interface import create_pmesh_from_regions

    absolute_outline = bounding_polygon is None
    if absolute_outline:
        bounding_polygon = _Grid(dem).extent()
        if boundary_tags is None:
            boundary_tags = {'bottom': [0], 'right': [1], 'top': [2],
                             'left': [3]}

    field = size_field_from_dem(
        dem, maximum_triangle_area, streams=streams, floodplain=floodplain,
        regions=regions,
        bounding_polygon=None if absolute_outline else bounding_polygon,
        interior_holes=interior_holes, grade=grade, band=band,
        poly_geo_reference=poly_geo_reference, verbose=verbose)
    if size_field is not None:
        field = _finer(field, size_field)

    if absolute_outline and poly_geo_reference is not None:
        bounding_polygon = poly_geo_reference.get_relative(bounding_polygon)

    return create_pmesh_from_regions(
        bounding_polygon, boundary_tags,
        maximum_triangle_area=maximum_triangle_area,
        interior_holes=interior_holes,
        poly_geo_reference=poly_geo_reference,
        size_field=field, relax_passes=relax_passes, verbose=verbose,
        **kwargs)


def _finer(field, other):
    """``field``, lowered to ``other`` wherever that is finer.

    ``other`` is sampled at the centres of ``field``'s cells, so the result
    is a raster as fast to look up as ``field``.
    """
    nrow, ncol = field.area.shape
    cs = field.cellsize
    X, Y = num.meshgrid(field.xllcorner + (num.arange(ncol) + 0.5) * cs,
                        field.yllcorner + (num.arange(nrow) + 0.5) * cs)
    try:
        values = num.broadcast_to(num.asarray(other(X, Y), float), X.shape)
    except Exception:
        # a callable that only takes scalars
        values = num.vectorize(lambda x, y: float(other(x, y)))(X, Y)
    finer = Size_field(num.minimum(field.area, values), field.xllcorner,
                       field.yllcorner, cs)
    finer.drainage_area = field.drainage_area
    return finer


def elevation_from_dem(domain, dem):
    """Mean DEM elevation over each triangle of a domain.

    Each triangle's elevation is the mean of the DEM, interpolated
    bilinearly, at four points: its centroid and the points halfway from
    the centroid to each vertex. On triangles much larger than the DEM's
    cells that is closer to the triangle's mean than the centroid value
    alone. Nodata takes the value of the nearest cell with data.

    Returns an array of centroid values, for
    ``domain.set_quantity('elevation', z, location='centroids')``.
    """
    from scipy import ndimage

    grid = _Grid(dem)
    z = grid.z
    bad = ~num.isfinite(z)
    if bad.all():
        raise ValueError('DEM has no data')
    if bad.any():
        idx = ndimage.distance_transform_edt(bad, return_distances=False,
                                             return_indices=True)
        z = z[tuple(idx)]

    v = domain.get_vertex_coordinates(absolute=True).reshape(-1, 3, 2)
    g = v.mean(axis=1)
    pts = [g] + [0.5 * (g + v[:, k]) for k in range(3)]
    pts = num.concatenate(pts)
    col = (pts[:, 0] - grid.x[0]) / grid.cellsize
    row = (pts[:, 1] - grid.y[0]) / grid.cellsize
    s = ndimage.map_coordinates(z, [row, col], order=1, mode='nearest')
    return s.reshape(4, -1).mean(axis=0)


def remove_mesh_pits(domain, outlets, method='fill', max_depth=None,
                     where=None, epsilon=1.0e-3, verbose=False):
    """Remove pits from a domain's elevation, on its cell-neighbour graph.

    A pit is a cell, or a group of cells, from which no path of falling
    cells leads to an outlet: on a coarse mesh, often where averaging the
    DEM over a cell filled the gully it drained by. Such a pit holds rain
    for good. A priority flood from the outlets finds them all.

    ``method='fill'`` raises each pit to its spill level (Barnes et al.
    2014). ``method='breach'`` instead lowers the cells along each pit's
    lowest way out until it drains (Lindsay 2016); a pit whose breach would
    cut deeper than ``max_depth``, or cut a cell outside ``where``, is
    filled instead. Breaching keeps the pit's own bed, and usually changes
    the volume held by the terrain less.

    Only cells in ``where`` are changed. A large cell that averages away a
    gully is a pit to remove; a small one is usually a real depression
    that should keep its water, so restrict ``where`` to the coarse cells
    when the mesh resolves the depressions that matter.

    The domain's elevation is set at the centroids (flat within each
    triangle).

    Parameters
    ----------
    domain : Domain
    outlets : list of str or array of int
        Boundary tags where water leaves (the cells on those boundary
        segments), or the cell ids themselves.
    method : 'fill' or 'breach'
    max_depth : float, optional
        Deepest cut a breach may make. Default no limit.
    where : boolean array, optional
        Cells that may change. Default all.
    epsilon : float, optional
        Fall from cell to cell along a filled or breached path.
    verbose : bool, optional

    Returns
    -------
    dict
        ``raised`` and ``lowered`` (numbers of cells), ``volume_raised``
        and ``volume_lowered``, and ``pits_left`` (cells in ``where``, other
        than outlets, still lower than all their neighbours).
    """
    if method not in ('fill', 'breach'):
        raise ValueError("method must be 'fill' or 'breach'")
    z0 = domain.quantities['elevation'].centroid_values.astype(float)
    n = len(z0)
    neighbours = domain.neighbours
    movable = num.ones(n, bool) if where is None \
        else num.asarray(where, bool)
    if movable.shape != (n,):
        raise ValueError('where must have one value per triangle')

    outlets = list(outlets)
    if outlets and isinstance(outlets[0], str):
        tags = set(outlets)
        missing = tags - set(domain.get_boundary_tags())
        if missing:
            raise ValueError('No boundary tagged %s' % sorted(missing))
        seeds = sorted({vol for (vol, _), tag in domain.boundary.items()
                        if tag in tags})
    else:
        seeds = sorted({int(i) for i in outlets})
    if not seeds:
        raise ValueError('No outlet cells')

    nb = neighbours.tolist()
    if method == 'fill':
        z = _fill(z0, nb, seeds, epsilon)
        z = num.where(movable, z, z0)
    else:
        depth = num.inf if max_depth is None else float(max_depth)
        z = _breach(z0, nb, seeds, movable, depth, epsilon)

    area = domain.areas
    up = z > z0
    down = z < z0
    lower = num.where(neighbours >= 0, z[num.maximum(neighbours, 0)],
                      num.inf).min(axis=1)
    pit = movable & (lower > z)
    pit[seeds] = False                    # outlets drain out of the domain
    stats = {'raised': int(up.sum()),
             'lowered': int(down.sum()),
             'volume_raised': float(((z - z0) * area)[up].sum()),
             'volume_lowered': float(((z0 - z) * area)[down].sum()),
             'pits_left': int(pit.sum())}
    if verbose:
        print('Pits: %d cells raised (%.4g), %d lowered (%.4g); %d left'
              % (stats['raised'], stats['volume_raised'], stats['lowered'],
                 stats['volume_lowered'], stats['pits_left']))
    domain.set_quantity('elevation', z, location='centroids')
    return stats


def _fill(z, neighbours, seeds, epsilon):
    """Priority-flood fill on a graph from ``seeds``."""
    import heapq

    zf = num.full(len(z), num.inf)
    heap = [(float(z[i]), i) for i in seeds]
    for v, i in heap:
        zf[i] = v
    heapq.heapify(heap)
    while heap:
        v, i = heapq.heappop(heap)
        for j in neighbours[i]:
            if j >= 0 and zf[j] == num.inf:
                zf[j] = max(z[j], v + epsilon)
                heapq.heappush(heap, (zf[j], j))
    # Cells no outlet reaches (an island of the mesh) keep their elevation
    return num.where(num.isfinite(zf), zf, z)


def _breach(z, neighbours, seeds, movable, max_depth, epsilon):
    """Priority-flood breaching on a graph from ``seeds``.

    The flood records the cell each cell was reached from, its way out over
    the lowest saddle. Reaching a cell lower than the one it came from (the
    bottom of a pit), it lowers the cells along that way out, each
    ``epsilon`` below the one before, until one is already low enough. A
    breach that would cut a fixed cell, or deeper than ``max_depth``, fills
    the pit instead.
    """
    import heapq

    n = len(z)
    zf = num.array(z, float)
    parent = num.full(n, -1, num.int64)
    seen = num.zeros(n, bool)
    heap = [(float(z[i]), i) for i in seeds]
    seen[seeds] = True
    heapq.heapify(heap)
    while heap:
        _, i = heapq.heappop(heap)
        v = zf[i]                           # may have been lowered by a breach
        for j in neighbours[i]:
            if j < 0 or seen[j]:
                continue
            seen[j] = True
            parent[j] = i
            if z[j] >= v or not movable[j]:
                # not a pit (a flat drains too), or one we may not touch
                zf[j] = z[j]
            else:
                path = []
                level = z[j] - epsilon
                k = i
                ok = True
                while k >= 0 and zf[k] > level:
                    if not movable[k] or z[k] - level > max_depth:
                        ok = False
                        break
                    path.append((k, level))
                    level -= epsilon
                    k = parent[k]
                if ok:
                    for k, lv in path:
                        zf[k] = lv
                    zf[j] = z[j]
                else:
                    zf[j] = v + epsilon
            heapq.heappush(heap, (zf[j], j))
    return zf

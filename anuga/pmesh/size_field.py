"""Mesh resolution given as a field of target triangle area.

``interior_regions`` vary the resolution of a mesh by adding closed polygons
to it, and Triangle honours every segment of those polygons exactly. That is
right when the outline is a physical feature. When it is not -- an elevation
band, a flood extent, a buffer around a town -- the outline's own geometry
(a sharp tooth, two outlines passing close together) forces very small
triangles, and the smallest triangle sets the timestep of the whole run.

A size field says how large a triangle may be at each place without adding
any segment to the mesh. It is applied through Triangle's user test: a
triangle is refined when its area exceeds the field at its centroid.

A size field is either

- any callable ``f(x, y) -> maximum triangle area``, with ``x`` and ``y`` in
  absolute coordinates, or
- a :class:`Size_field`, a raster of target areas, which is much faster to
  look up and is what :func:`size_field_from_regions` builds.

Example
-------
>>> import anuga
>>> field = anuga.size_field_from_regions(
...     bounding_polygon, [(town_polygon, 50.0)], maximum_triangle_area=5000.0)
>>> domain = anuga.create_domain_from_regions(
...     bounding_polygon, boundary_tags, size_field=field, relax_passes=30)
"""

import numpy as num

from anuga.utilities.numerical_tools import ensure_numeric

# Area of an equilateral triangle of unit edge: area = EQUILATERAL * edge**2
EQUILATERAL = 0.4330127018922193


def area_to_edge(area):
    """Edge length of the equilateral triangle with the given area."""
    return num.sqrt(num.asarray(area, float) / EQUILATERAL)


def edge_to_area(edge):
    """Area of the equilateral triangle with the given edge length."""
    return EQUILATERAL * num.asarray(edge, float) ** 2


class Size_field:
    """A raster of target (maximum) triangle area.

    Parameters
    ----------
    area : 2D array
        ``area[row, col]`` is the maximum triangle area in the cell whose
        lower left corner is ``(xllcorner + col*cellsize,
        yllcorner + row*cellsize)``. Row 0 is the SOUTH edge. Values must be
        positive.
    xllcorner, yllcorner : float
        Absolute coordinates of the lower left corner of the raster.
    cellsize : float
        Edge length of a (square) raster cell.

    Notes
    -----
    A point outside the raster takes the value of the nearest cell.

    The field is evaluated only at triangle centroids, so it should vary
    smoothly: if it jumps from coarse to fine over a short distance, a large
    triangle whose centroid lies on the coarse side can straddle the fine
    part unrefined. :func:`size_field_from_regions` grades its field for
    this reason.
    """

    def __init__(self, area, xllcorner, yllcorner, cellsize):
        area = num.ascontiguousarray(area, dtype=float)
        if area.ndim != 2 or area.size == 0:
            raise ValueError('Size_field area must be a non-empty 2D array')
        if not num.all(area > 0.0):          # also catches NaN
            raise ValueError('Size_field areas must all be positive')
        if not cellsize > 0.0:
            raise ValueError('Size_field cellsize must be positive')

        self.area = area
        self.xllcorner = float(xllcorner)
        self.yllcorner = float(yllcorner)
        self.cellsize = float(cellsize)

    def __call__(self, x, y):
        """Target area at absolute ``(x, y)``; scalars or arrays."""
        nrow, ncol = self.area.shape
        col = num.floor((num.asarray(x, float) - self.xllcorner)
                        / self.cellsize).astype(num.int64)
        row = num.floor((num.asarray(y, float) - self.yllcorner)
                        / self.cellsize).astype(num.int64)
        value = self.area[num.clip(row, 0, nrow - 1),
                          num.clip(col, 0, ncol - 1)]
        return float(value) if num.ndim(value) == 0 else value

    def refinement_func(self, xoffset=0.0, yoffset=0.0):
        """Triangle user test for this field.

        Triangle works in coordinates relative to the mesh's geo reference;
        ``xoffset, yoffset`` is that reference's lower left corner.
        """
        item = self.area.item
        inv = 1.0 / self.cellsize
        x0 = (xoffset - self.xllcorner) * inv
        y0 = (yoffset - self.yllcorner) * inv
        maxrow = self.area.shape[0] - 1
        maxcol = self.area.shape[1] - 1
        third = inv / 3.0

        # Called a few times for every triangle in the mesh: keep it lean.
        def too_big(p, tri_area):
            col = int((p[0][0] + p[1][0] + p[2][0]) * third + x0)
            row = int((p[0][1] + p[1][1] + p[2][1]) * third + y0)
            if col < 0:
                col = 0
            elif col > maxcol:
                col = maxcol
            if row < 0:
                row = 0
            elif row > maxrow:
                row = maxrow
            return tri_area > item(row, col)

        return too_big


def make_refinement_func(size_field, xoffset=0.0, yoffset=0.0):
    """Wrap a size field as Triangle's user test.

    Parameters
    ----------
    size_field : Size_field or callable
        A callable is called as ``size_field(x, y)`` with the absolute
        coordinates of a triangle's centroid and must return the maximum
        area allowed there.
    xoffset, yoffset : float
        Lower left corner of the mesh's geo reference, which Triangle's
        coordinates are relative to.
    """
    if hasattr(size_field, 'refinement_func'):
        return size_field.refinement_func(xoffset, yoffset)

    if not callable(size_field):
        raise TypeError('size_field must be a Size_field or a callable '
                        'f(x, y) -> maximum triangle area')

    def too_big(p, tri_area):
        x = (p[0][0] + p[1][0] + p[2][0]) / 3.0 + xoffset
        y = (p[0][1] + p[1][1] + p[2][1]) / 3.0 + yoffset
        target = size_field(x, y)
        if not target > 0.0:       # also catches NaN; would refine for ever
            raise ValueError('size_field(%g, %g) returned %s; the maximum '
                             'triangle area must be positive' % (x, y, target))
        return tri_area > target

    return too_big


def _rasterise_polygon(polygon, xll, yll, cellsize, nrow, ncol):
    """Cells of a raster whose centres lie inside a polygon.

    Scanline fill by the even-odd rule: every polygon edge is intersected
    with the rows of cell centres it spans, and a cell is inside if an odd
    number of those crossings lie to its left. The work is proportional to
    the polygon's perimeter and bounding box, not their product.

    Returns ``(row0, col0, inside)``, with ``inside`` a boolean window of the
    raster whose lower left cell is ``(row0, col0)``.
    """
    empty = (0, 0, num.zeros((0, 0), bool))

    x1 = (polygon[:, 0] - xll) / cellsize - 0.5      # in units of cells,
    y1 = (polygon[:, 1] - yll) / cellsize - 0.5      # relative to cell centres
    x2, y2 = num.roll(x1, -1), num.roll(y1, -1)

    # Rows crossed by each edge: those with  min(y1, y2) <= row < max(y1, y2)
    lo = num.clip(num.ceil(num.minimum(y1, y2)), 0, nrow).astype(num.int64)
    hi = num.clip(num.ceil(num.maximum(y1, y2)), 0, nrow).astype(num.int64)
    count = hi - lo
    total = int(count.sum())
    if total == 0:
        return empty

    edge = num.repeat(num.arange(len(count)), count)
    first = num.cumsum(count) - count
    row = lo[edge] + (num.arange(total) - first[edge])
    x = x1[edge] + (row - y1[edge]) * (x2[edge] - x1[edge]) \
        / (y2[edge] - y1[edge])
    # first column whose centre is right of the crossing
    col = num.clip(num.floor(x) + 1, 0, ncol).astype(num.int64)

    row0, row1 = int(row.min()), int(row.max()) + 1
    col0, col1 = int(col.min()), int(col.max())
    if col1 <= col0:
        return empty

    width = col1 - col0 + 1
    toggles = num.bincount((row - row0) * width + (col - col0),
                           minlength=(row1 - row0) * width)
    inside = num.cumsum(toggles.reshape(row1 - row0, width), axis=1) % 2 == 1
    return row0, col0, inside[:, :-1]


def size_field_from_regions(bounding_polygon,
                            regions,
                            maximum_triangle_area,
                            grade=0.5,
                            cellsize=None,
                            poly_geo_reference=None,
                            max_cells=20000000):
    """Build a graded :class:`Size_field` from polygons and resolutions.

    Takes the same ``(polygon, maximum area)`` list as ``interior_regions``,
    but the polygons only say where the mesh should be fine: none of their
    segments enter the mesh. They may therefore touch, overlap, have sharp
    corners or reach outside the bounding polygon. Where they overlap, the
    finest one wins.

    The target edge length ``h_k`` of each region holds inside it and grows
    away from it at ``grade`` metres per metre,

        h(x) = min( h_max,  min over regions k of  h_k + grade*dist(x, k) )

    with ``h_max`` the edge length of ``maximum_triangle_area``. Edge lengths
    are those of the equilateral triangle of the given area.

    Parameters
    ----------
    bounding_polygon : list of points
        Outline of the domain; the raster covers its bounding box.
    regions : list of (polygon, maximum area) tuples
        Regions to refine, as for ``interior_regions``.
    maximum_triangle_area : float
        Target area away from all regions.
    grade : float, optional
        Growth of the target edge length per unit distance from a region.
        Smaller is smoother and costs more triangles. 0 gives an ungraded
        field, stepping at the region outlines, which is not recommended
        (see :class:`Size_field`). Default 0.5.
    cellsize : float, optional
        Raster cell size. Defaults to the smallest target edge length,
        enlarged if needed to keep the raster under ``max_cells`` cells.
    poly_geo_reference : Geo_reference, optional
        Geo reference of the polygons. If None they are absolute.
    max_cells : int, optional
        Limit on the raster size when ``cellsize`` is not given.

    Returns
    -------
    Size_field
    """
    from scipy import ndimage

    def absolute(polygon):
        polygon = ensure_numeric(polygon, float)
        if poly_geo_reference is not None:
            polygon = ensure_numeric(
                poly_geo_reference.get_absolute(polygon), float)
        return polygon

    if not maximum_triangle_area > 0.0:
        raise ValueError('maximum_triangle_area must be positive')
    if grade < 0.0:
        raise ValueError('grade must not be negative')

    bounding_polygon = absolute(bounding_polygon)
    regions = [(absolute(polygon), float(area)) for polygon, area in regions]
    for _, area in regions:
        if not area > 0.0:
            raise ValueError('Region maximum areas must be positive')

    h_max = float(area_to_edge(maximum_triangle_area))
    edges = [min(float(area_to_edge(area)), h_max) for _, area in regions]

    xmin, ymin = bounding_polygon.min(axis=0)
    xmax, ymax = bounding_polygon.max(axis=0)
    if cellsize is None:
        cellsize = min(edges + [h_max])
        cellsize = max(cellsize,
                       num.sqrt((xmax - xmin) * (ymax - ymin) / max_cells))
    cellsize = float(cellsize)

    # One cell of margin, so that every centroid falls inside the raster
    xll, yll = xmin - cellsize, ymin - cellsize
    ncol = int(num.ceil((xmax - xll) / cellsize)) + 1
    nrow = int(num.ceil((ymax - yll) / cellsize)) + 1

    # Rasterise by cell centre, grouping regions of equal resolution so that
    # each distinct resolution costs one distance transform
    masks = {}
    for (polygon, _), edge in zip(regions, edges):
        mask = masks.setdefault(edge, num.zeros((nrow, ncol), bool))

        r0, c0, inside = _rasterise_polygon(polygon, xll, yll, cellsize,
                                            nrow, ncol)
        mask[r0:r0 + inside.shape[0], c0:c0 + inside.shape[1]] |= inside

        # A region narrower than a cell can miss every cell centre; its own
        # vertices keep it in the field
        col = ((polygon[:, 0] - xll) / cellsize).astype(int)
        row = ((polygon[:, 1] - yll) / cellsize).astype(int)
        ok = (col >= 0) & (col < ncol) & (row >= 0) & (row < nrow)
        mask[row[ok], col[ok]] = True

    h = num.full((nrow, ncol), h_max)
    for edge, mask in masks.items():
        if not mask.any():
            continue
        if grade > 0.0:
            dist = ndimage.distance_transform_edt(~mask, sampling=cellsize)
            num.minimum(h, edge + grade * dist, out=h)
        else:
            h[mask] = num.minimum(h[mask], edge)

    return Size_field(edge_to_area(h), xll, yll, cellsize)

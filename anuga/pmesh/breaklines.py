"""Prepare breaklines for meshing.

Breaklines that cross, touch or run close to each other, or to the outline of
an interior region, make the mesh generator fill the gap between them with
very small triangles. Such triangles can be thousands of times smaller than
the rest of the mesh. They cost nothing while dry, then cut the time step by
an order of magnitude once the flow reaches them.

:func:`clean_breaklines` rearranges a set of breaklines so that this cannot
happen: no two lines come closer than a chosen gap, except where one ends on
another at a shared vertex.

See also ``anuga.utilities.spatialInputUtil.add_intersections_to_domain_features``,
which keeps lines that cross and gives each crossing an exact shared vertex. It
does not deal with lines that nearly touch, run side by side or cross at a
shallow angle, which is where the slivers come from.
"""

import numpy as num


def _thin(coords, spacing):
    """Drop interior vertices closer than spacing to the last one kept or to
    the final vertex."""

    out = [coords[0]]
    last = coords[-1]
    for p in coords[1:-1]:
        if (num.hypot(p[0] - out[-1][0], p[1] - out[-1][1]) >= spacing
                and num.hypot(p[0] - last[0], p[1] - last[1]) >= spacing):
            out.append(p)
    return out + [last]


def _height_along(line):
    """Function z(x, y) giving the height of a 3D polyline at the point of it
    nearest (x, y), by linear interpolation along its length."""

    from shapely.geometry import LineString, Point

    xyz = num.array(line, float)
    path = LineString(xyz[:, :2])
    s = num.concatenate([[0.0], num.cumsum(num.hypot(*num.diff(xyz[:, :2],
                                                               axis=0).T))])

    def z(x, y):
        return float(num.interp(path.project(Point(x, y)), s, xyz[:, 2]))

    return z


def clean_breaklines(breaklines, interior_regions=None, bounding_polygon=None,
                     interior_holes=None, riverwalls=None, riverwallPar=None,
                     gap=10.0, minimum_length=None,
                     minimum_vertex_spacing=None, join=True, verbose=False):
    """Separate breaklines and riverwalls so that they mesh without sliver
    triangles.

    The lines are taken longest first, riverwalls before breaklines. Each
    later line is cut ``gap`` short of the lines already kept and of the
    outlines given, and the pieces shorter than ``minimum_length`` are
    dropped. With ``join``, each cut end is then joined by one straight
    segment to the line it was cut back from, at a vertex the two share, so
    a line that crossed or touched another now ends on it instead. A line
    that crosses an interior region outline is joined to the outline from
    both sides.

    Vertices are inserted in the interior region outlines where lines join
    them, so use the regions returned. The bounding polygon and the holes are
    never changed, because their boundary tags refer to their segments by
    number: lines are only cut short of them.

    All coordinates must be in the same frame (for example all absolute).

    Parameters
    ----------
    breaklines : list of lists of [x, y]
        The lines, each with at least two points. May be empty.
    interior_regions : list of [polygon, resolution], optional
        As passed to :func:`anuga.create_pmesh_from_regions`.
    bounding_polygon : list of [x, y], optional
        Lines are clipped to it and kept ``gap`` inside it.
    interior_holes : list of polygons, optional
        Lines are kept ``gap`` clear of them.
    riverwalls : dict of name -> list of [x, y, z], optional
        As passed to ``domain.riverwallData.create_riverwalls``. Crest
        heights are kept: a vertex added by a cut or a join takes the height
        of the original wall at the point nearest to it.
    riverwallPar : dict of name -> dict, optional
        Parameters of the riverwalls. Returned with an entry for every piece
        of a wall that had one.
    gap : float
        Smallest distance allowed between two lines, or between a line and
        an outline, away from a shared vertex. It should be comparable to the
        edge length of the triangles wanted there.
    minimum_length : float, optional
        Pieces shorter than this are dropped. Default ``4 * gap``.
    minimum_vertex_spacing : float, optional
        Vertices along a line closer together than this are thinned out, and
        a join uses an existing vertex within this distance instead of adding
        one. Default ``0.6 * gap``.
    join : bool
        Join cut ends to the line or region outline they were cut back from.
        If False the lines are left ``gap`` apart.
    verbose : bool

    Returns
    -------
    breaklines : list of lists of [x, y]
    interior_regions : list of [polygon, resolution], or None
        The regions with the join vertices added (None if none were given).
    riverwalls : dict of name -> list of [x, y, z]
        Only returned if ``riverwalls`` was given. A wall cut into several
        pieces gives ``name``, ``name_2``, ``name_3`` ... longest first.
    riverwallPar : dict
        Only returned if ``riverwalls`` was given.

    Notes
    -----
    Where two lines run side by side closer than ``gap``, only the longer is
    kept along that stretch (a riverwall is kept in preference to a
    breakline).

    Riverwalls are taken first, so only another riverwall can cut one. Where
    a wall crosses a longer wall it is joined to it from both sides, and the
    barrier stays closed. A wall that reaches the bounding polygon or a hole
    ends ``gap`` short of it, which leaves an opening there: extend the
    domain or accept the gap.

    The riverwalls must also be given to the mesh generator as breaklines.
    Check the resulting mesh with ``domain.check_mesh_quality``.

    Examples
    --------
    >>> lines, regions, walls, par = anuga.clean_breaklines(
    ...     lines, interior_regions, bounding_polygon, riverwalls=walls,
    ...     riverwallPar=par, gap=10.0)
    >>> wall_lines = [[p[:2] for p in wall] for wall in walls.values()]
    >>> domain = anuga.create_domain_from_regions(
    ...     bounding_polygon, boundary_tags, maximum_triangle_area=1000.0,
    ...     interior_regions=regions, breaklines=lines + wall_lines)
    >>> domain.riverwallData.create_riverwalls(walls, par)
    """

    try:
        from shapely.geometry import LineString, Point, Polygon
        from shapely.ops import linemerge, unary_union
    except ImportError as e:
        raise ImportError('clean_breaklines needs the shapely package') from e

    gap = float(gap)
    if gap <= 0.0:
        raise ValueError('gap must be positive')
    if minimum_length is None:
        minimum_length = 4.0 * gap
    spacing = 0.6 * gap if minimum_vertex_spacing is None \
        else float(minimum_vertex_spacing)

    def ring(polygon):
        c = [[float(p[0]), float(p[1])] for p in polygon]
        if c[0] != c[-1]:
            c.append(list(c[0]))
        return c

    # Lines kept so far, as lists of [x, y]. The outlines come first: region
    # outlines can be joined to, the bounding polygon and holes cannot.
    kept = []
    joinable = []      # may a later line be joined to this one?
    for region in (interior_regions or []):
        kept.append(ring(region[0]))
        joinable.append(join)
    fixed = ([bounding_polygon] if bounding_polygon is not None else [])
    fixed += list(interior_holes or [])
    for polygon in fixed:
        kept.append(ring(polygon))
        joinable.append(False)

    barrier = None     # everything within gap of a kept line
    for c in kept:
        b = LineString(c).buffer(gap)
        barrier = b if barrier is None else barrier.union(b)

    domain = Polygon(bounding_polygon) if bounding_polygon is not None else None

    def candidates(line):
        """The parts of a line inside the domain, split where it crosses
        itself."""
        c = [[float(p[0]), float(p[1])] for p in line]
        c = [p for i, p in enumerate(c) if i == 0 or p != c[i - 1]]
        if len(c) < 2:
            return []
        g = LineString(c)
        if domain is not None:
            g = g.intersection(domain)
        if not g.is_empty and not g.is_simple:
            g = unary_union(g)
            if g.geom_type == 'MultiLineString':
                g = linemerge(g)
        return [p for p in getattr(g, 'geoms', [g])
                if p.geom_type == 'LineString' and p.length >= minimum_length]

    def attach(pt):
        """Vertex on the nearest kept line to join pt to, or None."""
        point = Point(pt)
        dist = [LineString(c).distance(point) for c in kept]
        k = int(num.argmin(dist))
        if not joinable[k]:
            return None
        c = kept[k]
        line = LineString(c)
        s = line.project(point)
        q = line.interpolate(s)
        dv = [num.hypot(q.x - v[0], q.y - v[1]) for v in c]
        m = int(num.argmin(dv))
        if dv[m] < spacing:
            return list(c[m])
        seg = num.hypot(*num.diff(num.array(c), axis=0).T)
        n = int(num.searchsorted(num.concatenate([[0.0], num.cumsum(seg)]), s,
                                 side='right'))
        c.insert(min(max(n, 1), len(c) - 1), [q.x, q.y])
        return [q.x, q.y]

    # (2D line, name of its riverwall or None), riverwalls first, each group
    # longest first
    todo = []
    for name, wall in (riverwalls or {}).items():
        if any(len(p) < 3 for p in wall):
            raise ValueError('riverwall %s needs [x, y, z] points' % str(name))
        todo.append([(p, name) for p in candidates(wall)])
    todo = [sorted(sum(todo, []), key=lambda t: -t[0].length)]
    todo.append(sorted([(p, None) for line in breaklines
                        for p in candidates(line)],
                       key=lambda t: -t[0].length))

    pieces = []        # (index into kept, riverwall name or None)
    ndropped = 0
    for p, name in todo[0] + todo[1]:
        ends = (p.coords[0], p.coords[-1])
        g = p if barrier is None else p.difference(barrier)
        for piece in getattr(g, 'geoms', [g]):
            if piece.geom_type != 'LineString' or piece.length < minimum_length:
                ndropped += piece.geom_type == 'LineString'
                continue
            c = _thin([list(v) for v in piece.coords], spacing)
            if join and kept:
                # only ends made by the cut are joined, not the line's own
                first = None if tuple(c[0]) in ends else attach(c[0])
                last = None if tuple(c[-1]) in ends else attach(c[-1])
                if first is not None and first == last:
                    last = None
                # A join vertex closer than the vertex spacing to the cut end
                # (the other piece of the same line, when only a short stretch
                # was cut out) replaces the end: no segment shorter than the
                # spacing is added.
                if first is not None:
                    if num.hypot(first[0] - c[0][0],
                                 first[1] - c[0][1]) < spacing:
                        c[0] = first
                    else:
                        c = [first] + c
                if last is not None:
                    if num.hypot(last[0] - c[-1][0],
                                 last[1] - c[-1][1]) < spacing:
                        c[-1] = last
                    else:
                        c = c + [last]
            pieces.append((len(kept), name))
            kept.append(c)
            joinable.append(join)
            b = LineString(c).buffer(gap)
            barrier = b if barrier is None else barrier.union(b)

    # Read the lines back from kept: later joins have added vertices to them
    out = [kept[k] for k, name in pieces if name is None]

    regions = None
    if interior_regions is not None:
        regions = [[kept[i][:-1]] + list(region[1:])
                   for i, region in enumerate(interior_regions)]

    if verbose:
        total = sum(LineString(kept[k]).length for k, _ in pieces)
        nin = len(breaklines) + len(riverwalls or {})
        print('clean_breaklines: %d lines in, %d pieces out (%.0f m), '
              '%d short pieces dropped' % (nin, len(pieces), total, ndropped))

    if riverwalls is None:
        return out, regions

    # Riverwall pieces: heights from the original wall, names name, name_2 ...
    walls, par = {}, {}
    height = {name: _height_along(wall) for name, wall in riverwalls.items()}
    count = {}
    for k, name in pieces:
        if name is None:
            continue
        count[name] = count.get(name, 0) + 1
        new = name
        if count[name] > 1:
            new = '%s_%d' % (str(name), count[name])
            while new in riverwalls or new in walls:
                new += '_'
        walls[new] = [[x, y, height[name](x, y)] for x, y in kept[k]]
        if riverwallPar is not None and name in riverwallPar:
            par[new] = dict(riverwallPar[name])

    return out, regions, walls, par

"""Worst-cell relaxation: enlarge the smallest cells of a triangulation.

ANUGA steps at the global minimum of ``radii / wave_speed``, so the smallest
cell in a mesh sets the timestep of the whole run. Triangle guarantees a
minimum angle and a maximum area, not a minimum size, and an ordinary
refinement leaves a thin tail of cells at around a third of the typical size.
Moving a few vertices removes most of that tail at no cost in triangles.
"""

import numpy as num


def segment_sliders(vertices, segments, first_inserted):
    """Vertices that may slide along a segment without changing its shape.

    These are the vertices the mesh generator inserted to split a straight
    input segment: index at least ``first_inserted`` (Triangle numbers the
    input points first), on exactly two segment edges, and in line with the
    two vertices at their other ends. Input vertices -- a hole's corners,
    a breakline's vertices -- are never included, and nor is a vertex where
    segments meet or cross.

    Returns an (K, 3) int array of rows ``(i, a, b)``: vertex ``i`` lies on
    the segment between ``a`` and ``b``.
    """
    v = num.asarray(vertices, dtype=float)
    seg = num.asarray(segments, dtype=num.int64).reshape(-1, 2)
    if len(seg) == 0:
        return num.zeros((0, 3), num.int64)
    ends = num.concatenate([seg, seg[:, ::-1]])            # (vertex, other end)
    ends = ends[num.argsort(ends[:, 0], kind='stable')]
    vid, first, count = num.unique(ends[:, 0], return_index=True,
                                   return_counts=True)
    two = (count == 2) & (vid >= first_inserted)
    i = vid[two]
    a = ends[first[two], 1]
    b = ends[first[two] + 1, 1]
    # In line: the cross product is negligible against the span, and i is
    # between a and b
    ab = v[b] - v[a]
    ai = v[i] - v[a]
    cross = num.abs(ab[:, 0] * ai[:, 1] - ab[:, 1] * ai[:, 0])
    span2 = num.einsum('ij,ij->i', ab, ab)
    proj = num.einsum('ij,ij->i', ai, ab)
    ok = (cross <= 1.0e-9 * span2) & (proj > 0.0) & (proj < span2)
    return num.ascontiguousarray(num.stack([i, a, b], axis=1)[ok],
                                 dtype=num.int64)


def relax_worst_cells(vertices, triangles, fixed=None,
                      passes=30, fraction=0.5, verbose=False, slide=None):
    """Move the vertices of the smallest cells so the smallest cell grows.

    Each pass takes the worst ``fraction`` percent of cells and tries moving
    each of their free vertices: to the centroid of its neighbours, part-way
    there, and a small step in each of eight directions. The move that most
    raises the worst cell around the vertex is kept, and none is made unless
    one raises it. Triangles are never added, removed or reconnected, and a
    move that would invert a triangle is rejected, so the mesh stays valid.

    A cell is scored by ``min(radii, inradius)``, where ``radii`` is ANUGA's
    default (the shortest distance from the centroid to an edge midpoint).
    ``radii`` alone can be satisfied by a thin triangle of almost no area;
    the inradius term rules that out. The two coincide on an equilateral
    triangle.

    Parameters
    ----------
    vertices : (N, 2) float array
        Vertex coordinates.
    triangles : (M, 3) int array
        Vertex indices of each triangle.
    fixed : array of int, optional
        Indices of vertices that must not move. This MUST include every
        vertex on the boundary, and should include every vertex on a
        segment (breakline, region outline).
    passes : int, optional
        Maximum number of passes. Stops earlier once the worst cell has
        stopped improving. Default 30.
    fraction : float, optional
        Percentage of cells treated as "worst" in each pass. Default 0.5.
    verbose : bool, optional
    slide : (K, 3) int array, optional
        Rows ``(i, a, b)`` from :func:`segment_sliders`: vertex ``i`` may
        move along the line through ``a`` and ``b``, strictly between them,
        even if it is in ``fixed``. This frees the vertices the generator
        inserted along a hole's sides or a breakline without changing their
        shape.

    Returns
    -------
    vertices : (N, 2) float array
        The relaxed coordinates (a new array; the input is not modified).
    before, after : float
        The quality of the worst cell before and after.

    Notes
    -----
    The result is in general no longer Delaunay, and a few cells can end up
    slightly larger than the maximum area the mesh was generated with. ANUGA
    requires neither property.
    """
    from anuga.pmesh.mesh_relax_ext import cell_quality, sweep, slide_sweep
    import anuga.utilities.log as log

    v = num.array(vertices, dtype=float, order='C')        # a copy
    t = num.ascontiguousarray(triangles, dtype=num.int64)
    if v.ndim != 2 or v.shape[1] != 2:
        raise ValueError('vertices must be an (N, 2) array')
    if t.ndim != 2 or t.shape[1] != 3:
        raise ValueError('triangles must be an (M, 3) array')
    if len(t) == 0 or passes <= 0:
        return v, num.nan, num.nan

    # Orientation of the triangulation (Triangle's is anticlockwise)
    a, b, c = v[t[:, 0]], v[t[:, 1]], v[t[:, 2]]
    cross = ((b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1])
             - (c[:, 0] - a[:, 0]) * (b[:, 1] - a[:, 1]))
    sign = 1.0 if num.median(cross) > 0.0 else -1.0

    # Triangles around each vertex, in compressed-row form
    flat = t.ravel()
    order = num.argsort(flat, kind='stable')
    vt = num.ascontiguousarray(order // 3, dtype=num.int64)
    start = num.ascontiguousarray(
        num.searchsorted(flat[order], num.arange(len(v) + 1)), dtype=num.int64)

    free = num.ones(len(v), bool)
    if fixed is not None:
        free[num.asarray(fixed, dtype=num.int64).ravel()] = False

    if slide is None:
        slide = num.zeros((0, 3), num.int64)
    slide = num.ascontiguousarray(slide, dtype=num.int64).reshape(-1, 3)
    row_of = num.full(len(v), -1, num.int64)
    row_of[slide[:, 0]] = num.arange(len(slide))
    free[slide[:, 0]] = False                     # they move only by sliding

    q = cell_quality(v, t, sign)
    before = float(q.min())
    history = []
    for _ in range(passes):
        threshold = num.percentile(q, fraction)
        bad = q <= threshold
        touched = num.unique(t[bad])
        rows = row_of[touched]
        rows = rows[rows >= 0]
        cand = touched[free[touched]]
        if len(cand) == 0 and len(rows) == 0:
            break
        # Worst first. A candidate's worst triangle is one of the bad ones.
        worst_at = num.full(len(v), num.inf)
        num.minimum.at(worst_at, t[bad].ravel(), num.repeat(q[bad], 3))
        cand = num.ascontiguousarray(cand[num.argsort(worst_at[cand])],
                                     dtype=num.int64)
        rows = num.ascontiguousarray(
            rows[num.argsort(worst_at[slide[rows, 0]])], dtype=num.int64)
        # strict threshold in the sweep, so step past it for a tiny mesh
        limit = num.nextafter(threshold, num.inf)
        moved = sweep(v, t, q, vt, start, cand, limit, sign)
        moved += slide_sweep(v, t, q, vt, start, slide, rows, limit, sign)
        history.append(float(q.min()))
        if moved == 0:
            break
        if len(history) > 8 and history[-1] < history[-9] * 1.005:
            break

    after = float(q.min())
    if verbose:
        log.critical('Relaxed mesh in %d passes: worst cell %.4g -> %.4g '
                     '(%.2fx)' % (len(history), before, after, after / before))

    return v, before, after

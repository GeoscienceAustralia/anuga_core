"""Worst-cell relaxation: enlarge the smallest cells of a triangulation.

ANUGA steps at the global minimum of ``radii / wave_speed``, so the smallest
cell in a mesh sets the timestep of the whole run. Triangle guarantees a
minimum angle and a maximum area, not a minimum size, and an ordinary
refinement leaves a thin tail of cells at around a third of the typical size.
Moving a few vertices removes most of that tail at no cost in triangles.
"""

import numpy as num


def relax_worst_cells(vertices, triangles, fixed=None,
                      passes=30, fraction=0.5, verbose=False):
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
    from anuga.pmesh.mesh_relax_ext import cell_quality, sweep
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

    q = cell_quality(v, t, sign)
    before = float(q.min())
    history = []
    for _ in range(passes):
        threshold = num.percentile(q, fraction)
        bad = q <= threshold
        cand = num.unique(t[bad])
        cand = cand[free[cand]]
        if len(cand) == 0:
            break
        # Worst first. A candidate's worst triangle is one of the bad ones.
        worst_at = num.full(len(v), num.inf)
        num.minimum.at(worst_at, t[bad].ravel(), num.repeat(q[bad], 3))
        cand = num.ascontiguousarray(cand[num.argsort(worst_at[cand])],
                                     dtype=num.int64)
        # strict threshold in the sweep, so step past it for a tiny mesh
        moved = sweep(v, t, q, vt, start, cand,
                      num.nextafter(threshold, num.inf), sign)
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

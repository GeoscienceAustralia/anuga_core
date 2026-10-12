#cython: wraparound=False, boundscheck=False, cdivision=True, profile=False, nonecheck=False, overflowcheck=False, cdivision_warnings=False, unraisable_tracebacks=False
"""Kernels for anuga.pmesh.mesh_relax.relax_worst_cells."""

import cython
from libc.math cimport hypot, cos, sin, M_PI
from libc.stdint cimport int64_t

import numpy as np
cimport numpy as np

cdef enum:
    NDIRS = 8
    NTRIAL = 19          # 3 Laplacian fractions + 2 step lengths * NDIRS


cdef inline double _quality(double[:, ::1] v, int64_t a, int64_t b, int64_t c,
                            double sign) noexcept nogil:
    """min(ANUGA radii, inradius) of one triangle, or -1 if it is inverted."""
    cdef double ax = v[a, 0], ay = v[a, 1]
    cdef double bx = v[b, 0], by = v[b, 1]
    cdef double cx = v[c, 0], cy = v[c, 1]
    cdef double area2, gx, gy, r, r1, r2, per

    area2 = ((bx - ax) * (cy - ay) - (cx - ax) * (by - ay)) * sign
    if area2 <= 0.0:
        return -1.0

    # ANUGA's default radii: shortest centroid to edge-midpoint distance
    gx = (ax + bx + cx) / 3.0
    gy = (ay + by + cy) / 3.0
    r = hypot(gx - 0.5 * (bx + cx), gy - 0.5 * (by + cy))
    r1 = hypot(gx - 0.5 * (cx + ax), gy - 0.5 * (cy + ay))
    r2 = hypot(gx - 0.5 * (ax + bx), gy - 0.5 * (ay + by))
    if r1 < r:
        r = r1
    if r2 < r:
        r = r2

    # inradius = 2*area/perimeter
    per = (hypot(bx - ax, by - ay) + hypot(cx - bx, cy - by)
           + hypot(ax - cx, ay - cy))
    r1 = area2 / per
    if r1 < r:
        r = r1
    return r


def cell_quality(double[:, ::1] v not None,
                 int64_t[:, ::1] t not None,
                 double sign):
    """Quality of every triangle: min(ANUGA radii, inradius), -1 if inverted.

    ``sign`` is +1 for anticlockwise triangles, -1 for clockwise.
    """
    cdef int64_t k, M = t.shape[0]
    cdef np.ndarray[double, ndim=1] q_arr = np.empty(M, dtype=np.float64)
    cdef double[::1] q = q_arr

    with nogil:
        for k in range(M):
            q[k] = _quality(v, t[k, 0], t[k, 1], t[k, 2], sign)
    return q_arr


def sweep(double[:, ::1] v not None,
          int64_t[:, ::1] t not None,
          double[::1] q not None,
          int64_t[::1] vt not None,
          int64_t[::1] start not None,
          int64_t[::1] cand not None,
          double thr,
          double sign):
    """One relaxation pass over the candidate vertices. Returns the number moved.

    For each vertex in ``cand`` whose worst incident triangle is still below
    ``thr``, try a set of new positions and keep the one that most raises the
    worst incident triangle, if any raises it by at least 1%.

    v     : (N, 2) vertex coordinates, updated in place
    t     : (M, 3) triangles
    q     : (M,) triangle quality, updated in place
    vt    : triangle ids grouped by vertex
    start : (N+1,) offsets into vt; triangles of vertex i are
            vt[start[i]:start[i+1]]
    cand  : candidate vertices. They must be interior vertices.
    """
    cdef int64_t moved = 0
    cdef int64_t n, i, j, k, s, e, m, d, c
    cdef double kx, ky, sx, sy, sl, mx, my, edge, nn
    cdef double old, best, bx, by, worst, w, f
    cdef double trial[NTRIAL][2]
    cdef double dirs[NDIRS][2]
    cdef double lap[3]
    cdef double step[2]

    lap[0] = 1.0
    lap[1] = 0.5
    lap[2] = 0.25
    step[0] = 0.2
    step[1] = 0.07
    for d in range(NDIRS):
        dirs[d][0] = cos(d * M_PI / 4.0)
        dirs[d][1] = sin(d * M_PI / 4.0)

    with nogil:
        for n in range(cand.shape[0]):
            i = cand[n]
            s = start[i]
            e = start[i + 1]
            if e <= s:
                continue

            old = q[vt[s]]
            for j in range(s + 1, e):
                if q[vt[j]] < old:
                    old = q[vt[j]]
            if old >= thr:           # a neighbour's move already rescued it
                continue

            kx = v[i, 0]
            ky = v[i, 1]

            # Centroid of, and mean distance to, the ring of neighbours.
            # Every ring neighbour of an interior vertex is in exactly two
            # incident triangles, so summing over triangles double-counts
            # them all evenly.
            sx = 0.0
            sy = 0.0
            sl = 0.0
            for j in range(s, e):
                k = vt[j]
                for c in range(3):
                    m = t[k, c]
                    if m != i:
                        sx += v[m, 0]
                        sy += v[m, 1]
                        sl += hypot(v[m, 0] - kx, v[m, 1] - ky)
            nn = 2.0 * (e - s)
            mx = sx / nn
            my = sy / nn
            edge = sl / nn

            # Laplacian moves, whole and partial, then small steps in eight
            # directions. The directional steps are what free the cells whose
            # neighbour centroid is no better than where they are.
            for m in range(3):
                trial[m][0] = kx + lap[m] * (mx - kx)
                trial[m][1] = ky + lap[m] * (my - ky)
            for m in range(2):
                f = step[m] * edge
                for d in range(NDIRS):
                    trial[3 + m * NDIRS + d][0] = kx + f * dirs[d][0]
                    trial[3 + m * NDIRS + d][1] = ky + f * dirs[d][1]

            best = old * 1.01
            bx = kx
            by = ky
            for m in range(NTRIAL):
                v[i, 0] = trial[m][0]
                v[i, 1] = trial[m][1]
                worst = 1.0e300
                for j in range(s, e):
                    k = vt[j]
                    w = _quality(v, t[k, 0], t[k, 1], t[k, 2], sign)
                    if w < worst:
                        worst = w
                if worst > best:
                    best = worst
                    bx = trial[m][0]
                    by = trial[m][1]

            v[i, 0] = bx
            v[i, 1] = by
            if bx != kx or by != ky:
                moved += 1
                for j in range(s, e):
                    k = vt[j]
                    q[k] = _quality(v, t[k, 0], t[k, 1], t[k, 2], sign)

    return moved

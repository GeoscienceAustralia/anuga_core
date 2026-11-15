#cython: wraparound=False, boundscheck=False, cdivision=True, profile=False, nonecheck=False, overflowcheck=False, cdivision_warnings=False, unraisable_tracebacks=False
"""Kernels for anuga.pmesh.dem_mesh: D8 flow routing and height above drainage."""

import cython
from libc.math cimport isnan, isinf, INFINITY
from libc.stdint cimport int64_t

import numpy as np
cimport numpy as np


# A binary min-heap of (key, id) pairs held in two arrays

cdef inline void _push(double[::1] key, int64_t[::1] ident, int64_t *size,
                       double k, int64_t i) noexcept nogil:
    cdef int64_t c = size[0], p
    size[0] += 1
    while c > 0:
        p = (c - 1) >> 1
        if key[p] <= k:
            break
        key[c] = key[p]
        ident[c] = ident[p]
        c = p
    key[c] = k
    ident[c] = i


cdef inline int64_t _pop(double[::1] key, int64_t[::1] ident,
                         int64_t *size) noexcept nogil:
    cdef int64_t top = ident[0], c = 0, child, n
    cdef double k
    cdef int64_t i
    size[0] -= 1
    n = size[0]
    k = key[n]
    i = ident[n]
    while True:
        child = 2 * c + 1
        if child >= n:
            break
        if child + 1 < n and key[child + 1] < key[child]:
            child += 1
        if key[child] >= k:
            break
        key[c] = key[child]
        ident[c] = ident[child]
        c = child
    key[c] = k
    ident[c] = i
    return top


def route(double[:, ::1] z not None, double epsilon):
    """D8 routing of a grid, with depressions filled by priority flood.

    ``z[row, col]`` is the elevation on square cells. NaN is nodata, and water
    leaves the grid there and across the grid's edge. +inf is a wall (a
    building, say): water neither enters nor crosses it.

    The surface is filled by priority flood (Barnes et al. 2014) with an
    ``epsilon`` gradient, so every cell has a strictly lower neighbour on the
    way out, and each cell then drains to its steepest-descent neighbour.

    Returns ``(receiver, order, accumulation)``: each cell's downstream flat
    index (-1 where water leaves the grid), the routed cells in ascending
    filled elevation (so a cell's receiver always comes before it), and the
    number of cells draining through each cell, itself included.
    """
    cdef int64_t nr = z.shape[0], nc = z.shape[1]
    cdef int64_t n = nr * nc
    cdef int64_t r, c, rr, cc, i, j, k, m, size = 0
    cdef int dr[8]
    cdef int dc[8]
    cdef double dist[8]
    cdef double v, best, s
    cdef bint edge

    dr[:] = [-1, -1, 0, 1, 1, 1, 0, -1]
    dc[:] = [0, 1, 1, 1, 0, -1, -1, -1]
    dist[:] = [1.0, 1.4142135623730951, 1.0, 1.4142135623730951,
               1.0, 1.4142135623730951, 1.0, 1.4142135623730951]

    zf_a = np.full(n, np.inf)
    done_a = np.zeros(n, np.uint8)
    order_a = np.empty(n, np.int64)
    rec_a = np.full(n, -1, np.int64)
    acc_a = np.zeros(n)
    key_a = np.empty(n, np.float64)
    ident_a = np.empty(n, np.int64)
    cdef double[::1] zf = zf_a
    cdef unsigned char[::1] done = done_a
    cdef int64_t[::1] order = order_a
    cdef int64_t[::1] rec = rec_a
    cdef double[::1] acc = acc_a
    cdef double[::1] key = key_a
    cdef int64_t[::1] ident = ident_a

    with nogil:
        # Seed from cells where water can leave: the grid's edge and cells
        # next to nodata
        for r in range(nr):
            for c in range(nc):
                i = r * nc + c
                v = z[r, c]
                if isnan(v) or isinf(v):
                    done[i] = 1
                    continue
                edge = r == 0 or c == 0 or r == nr - 1 or c == nc - 1
                if not edge:
                    for k in range(8):
                        if isnan(z[r + dr[k], c + dc[k]]):
                            edge = True
                            break
                if edge:
                    zf[i] = v
                    done[i] = 1
                    _push(key, ident, &size, v, i)

        m = 0
        while size > 0:
            i = _pop(key, ident, &size)
            order[m] = i
            m += 1
            v = zf[i]
            r = i // nc
            c = i % nc
            for k in range(8):
                rr = r + dr[k]
                cc = c + dc[k]
                if rr < 0 or cc < 0 or rr >= nr or cc >= nc:
                    continue
                j = rr * nc + cc
                if done[j]:
                    continue
                zf[j] = z[rr, cc] if z[rr, cc] > v + epsilon else v + epsilon
                done[j] = 1
                _push(key, ident, &size, zf[j], j)

        # Steepest descent on the filled surface
        for k in range(m):
            i = order[k]
            r = i // nc
            c = i % nc
            best = 0.0
            for j in range(8):
                rr = r + dr[j]
                cc = c + dc[j]
                if rr < 0 or cc < 0 or rr >= nr or cc >= nc:
                    continue
                if zf[rr * nc + cc] == INFINITY:
                    continue
                s = (zf[i] - zf[rr * nc + cc]) / dist[j]
                if s > best:
                    best = s
                    rec[i] = rr * nc + cc

        # Accumulate from the top down
        for k in range(m - 1, -1, -1):
            i = order[k]
            acc[i] += 1.0
            if rec[i] >= 0:
                acc[rec[i]] += acc[i]

    return rec_a, order_a[:m].copy(), acc_a.reshape(nr, nc)


def height_above_drainage(double[:, ::1] z not None,
                          int64_t[::1] receiver not None,
                          int64_t[::1] order not None,
                          double[:, ::1] accumulation not None,
                          double min_accumulation):
    """Height of each cell above the stream cell it drains to (HAND).

    A stream cell is one with at least ``min_accumulation`` cells draining
    through it. The arguments are ``z`` and the results of :func:`route`.
    Returns a grid, +inf where the flow path leaves the grid without meeting
    a stream, and at nodata and walls.
    """
    cdef int64_t nr = z.shape[0], nc = z.shape[1]
    cdef int64_t k, i
    base_a = np.full(nr * nc, -np.inf)
    out_a = np.full(nr * nc, np.inf)
    cdef double[::1] base = base_a
    cdef double[::1] out = out_a
    cdef double h

    with nogil:
        for k in range(order.shape[0]):          # a receiver comes first
            i = order[k]
            if accumulation[i // nc, i % nc] >= min_accumulation:
                base[i] = z[i // nc, i % nc]
            elif receiver[i] >= 0:
                base[i] = base[receiver[i]]
        for k in range(order.shape[0]):
            i = order[k]
            if base[i] > -INFINITY:
                h = z[i // nc, i % nc] - base[i]
                out[i] = h if h > 0.0 else 0.0

    return out_a.reshape(nr, nc)

#cython: wraparound=False, boundscheck=False, cdivision=True, profile=False, nonecheck=False, overflowcheck=False, cdivision_warnings=False, unraisable_tracebacks=False
import cython
from libc.stdint cimport int64_t

import numpy as np
cimport numpy as np


def build_neighbour_structure(int64_t N,
                               np.ndarray[int64_t, ndim=2, mode="c"] triangles not None,
                               np.ndarray[int64_t, ndim=2, mode="c"] neighbours not None,
                               np.ndarray[int64_t, ndim=2, mode="c"] neighbour_edges not None,
                               np.ndarray[int64_t, ndim=1, mode="c"] number_of_boundaries not None):
    """Build neighbour structure using sort-based directed-edge matching.

    Encodes each directed edge (u, v) as the integer key u*Mmax + v, sorts all
    3M keys, then binary-searches for the reverse key v*Mmax + u to identify
    the neighbouring triangle and shared edge.  O(M log M), cache-friendly, and
    avoids hash-table overhead.

    Parameters
    ----------
    N : int64
        Number of mesh nodes (kept for API compatibility).
    triangles : (M, 3) int64 array
        Vertex indices of each triangle [a, b, c].
    neighbours : (M, 3) int64 array
        Output: neighbouring triangle index per edge (-1 = boundary).
        MUST be pre-filled with -1 by the caller: only matched edges are
        written, and number_of_boundaries is derived by counting entries that
        are still negative. (The C++ implementation this replaced counted down
        from 3 instead, so it did not need the pre-fill; callers in
        basic_mesh.py and neighbour_mesh.py already do it.)
    neighbour_edges : (M, 3) int64 array
        Output: edge index in the neighbouring triangle.
    number_of_boundaries : (M,) int64 array
        Output: number of boundary edges per triangle.
    """

    cdef int64_t M       = triangles.shape[0]
    cdef int64_t Mmax    = np.max(triangles) + 1
    cdef int64_t three_M = 3 * M
    cdef int64_t i, j, lo, hi, mid, rev_key
    cdef int64_t d1, d2

    # ------------------------------------------------------------------
    # Allocate per-edge arrays (3M entries, one per directed edge)
    # ------------------------------------------------------------------
    cdef np.ndarray[int64_t, ndim=1] keys_arr   = np.empty(three_M, dtype=np.int64)
    cdef np.ndarray[int64_t, ndim=1] starts_arr = np.empty(three_M, dtype=np.int64)
    cdef np.ndarray[int64_t, ndim=1] ends_arr   = np.empty(three_M, dtype=np.int64)
    cdef np.ndarray[int64_t, ndim=1] tri_arr    = np.empty(three_M, dtype=np.int64)
    cdef np.ndarray[int64_t, ndim=1] edge_arr   = np.empty(three_M, dtype=np.int64)

    cdef int64_t[::1] keys   = keys_arr
    cdef int64_t[::1] starts = starts_arr
    cdef int64_t[::1] ends   = ends_arr
    cdef int64_t[::1] tris   = tri_arr
    cdef int64_t[::1] edges  = edge_arr

    cdef int64_t[:,::1] tri_mv  = triangles
    cdef int64_t[:,::1] nb_mv   = neighbours
    cdef int64_t[:,::1] nbe_mv  = neighbour_edges
    cdef int64_t[::1]   nbb_mv  = number_of_boundaries

    # ------------------------------------------------------------------
    # Populate edge arrays.
    # Triangle [a, b, c] has three directed edges:
    #   edge 0: b -> c   (opposite vertex a, index 0)
    #   edge 1: c -> a   (opposite vertex b, index 1)
    #   edge 2: a -> b   (opposite vertex c, index 2)
    # ------------------------------------------------------------------
    for i in range(M):
        starts[i]     = tri_mv[i, 1]; ends[i]     = tri_mv[i, 2]
        tris[i]       = i;            edges[i]     = 0
        keys[i]       = tri_mv[i, 1] * Mmax + tri_mv[i, 2]

        starts[M+i]   = tri_mv[i, 2]; ends[M+i]   = tri_mv[i, 0]
        tris[M+i]     = i;            edges[M+i]   = 1
        keys[M+i]     = tri_mv[i, 2] * Mmax + tri_mv[i, 0]

        starts[2*M+i] = tri_mv[i, 0]; ends[2*M+i] = tri_mv[i, 1]
        tris[2*M+i]   = i;            edges[2*M+i] = 2
        keys[2*M+i]   = tri_mv[i, 0] * Mmax + tri_mv[i, 1]

    # ------------------------------------------------------------------
    # Sort edges by key; obtain index array and sorted key array.
    # ------------------------------------------------------------------
    cdef np.ndarray[int64_t, ndim=1] sidx_arr = np.argsort(keys_arr, kind='stable')
    cdef np.ndarray[int64_t, ndim=1] sk_arr   = keys_arr[sidx_arr]
    cdef int64_t[::1] sidx = sidx_arr
    cdef int64_t[::1] sk   = sk_arr

    # ------------------------------------------------------------------
    # Detect duplicate directed edges (invalid mesh topology).
    # In sorted order, duplicates appear as adjacent equal keys.
    # ------------------------------------------------------------------
    for i in range(three_M - 1):
        if sk[i] == sk[i + 1]:
            d1 = sidx[i]
            d2 = sidx[i + 1]
            raise Exception(
                'Edge %d of triangle %d is duplicating edge %d of triangle %d.\n'
                % (edges[d2], tris[d2], edges[d1], tris[d1]))

    # ------------------------------------------------------------------
    # For every directed edge, binary-search for its reverse to find
    # the neighbouring triangle and shared edge index.
    # ------------------------------------------------------------------
    for i in range(three_M):
        rev_key = ends[i] * Mmax + starts[i]

        # bisect_left: find leftmost position where sk[pos] == rev_key
        lo = 0
        hi = three_M
        while lo < hi:
            mid = (lo + hi) >> 1
            if sk[mid] < rev_key:
                lo = mid + 1
            else:
                hi = mid

        if lo < three_M and sk[lo] == rev_key:
            j = sidx[lo]
            nb_mv[tris[i],  edges[i]] = tris[j]
            nbe_mv[tris[i], edges[i]] = edges[j]

    # ------------------------------------------------------------------
    # Count boundary edges (neighbours still == -1) per triangle.
    # ------------------------------------------------------------------
    for i in range(M):
        nbb_mv[i] = 0
        for j in range(3):
            if nb_mv[i, j] < 0:
                nbb_mv[i] += 1

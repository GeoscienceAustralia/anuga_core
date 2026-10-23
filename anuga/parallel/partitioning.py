

import numpy as np

#==============================================================================================================
# Code for partitioning meshes using METIS, Morton codes, and Hilbert codes.
#==============================================================================================================


#==============================================================================================================
# Partitioning using METIS. This minimizes edge cuts and is ideal for parallel computations.
#==============================================================================================================

def metis_partition(domain, n_procs):
    """
    Partition a mesh using METIS partitioning library.

    This function is a wrapper for the METIS partitioning routines in pymetis.
    The specific partitioning routine used depends on the installed METIS version:
    - METIS 4: Uses partMeshNodal
    - METIS 5: Uses part_mesh if available, otherwise part_graph

    Parameters
    ----------
    domain : object
        Domain object containing mesh information including number of triangles,
        nodes, triangle connectivity, and neighbor information.
    n_procs : int
        Number of processors to partition the mesh across.

    Returns
    -------
    epart_order : ndarray
        Integer array of indices that sort elements by partition assignment.
        Can be used to reorder triangles such that all triangles assigned to
        the same processor are contiguous.
    triangles_per_proc : ndarray
        Integer array of length n_procs where triangles_per_proc[i] is the
        number of triangles assigned to processor i.

    Raises
    ------
    ValueError
        If n_procs is greater than the number of triangles in the domain.
    AssertionError
        If partitioning results in any processor owning zero triangles.

    Notes
    -----
    The partitioning minimizes edge cuts to improve communication efficiency
    in parallel computations. For single processor (n_procs=1), returns trivial
    partitioning without invoking METIS.

    The function handles edge cases where METIS returns zero-dimensional arrays
    by flattening them to standard 1D arrays.
    """

    from pymetis import part_graph

    try:
        from pymetis import part_mesh
        metis_version = "5_part_mesh"
    except ImportError:
        metis_version = "5_part_graph"

    n_tri = domain.number_of_triangles

    if n_procs > n_tri:
        raise ValueError("Number of processors must be less than or equal to the number of triangles")

    if n_procs == 1:
        epart_order = np.arange(n_tri, dtype=int)
        triangles_per_proc = [n_tri]
        return epart_order, triangles_per_proc

    # Use metis to partition the mesh.
    # The partitioning routine used depends on the version of metis installed.
    # If metis 4 is installed, partMeshNodal is used. If metis 5 is installed,
    # part_mesh is used if available, otherwise part_graph

    if metis_version == "5_part_mesh":

        objval, epart, npart = part_mesh(n_procs, domain.triangles)


    if metis_version == "5_part_graph":
        # neighbours uses negative entries for boundary edges; pymetis can't
        # handle them, so we build a CSR adjacency that omits them.
        N = domain.neighbours
        mask = N >= 0
        adjncy = N[mask].astype(np.int32)
        counts = mask.sum(axis=1).astype(np.int32)
        xadj = np.empty(N.shape[0] + 1, dtype=np.int32)
        xadj[0] = 0
        np.cumsum(counts, out=xadj[1:])

        cutcount, partvert = part_graph(n_procs, xadj=xadj, adjncy=adjncy)

        epart = partvert

    # Sometimes (usu. on x86_64), partMeshNodal returns an array of zero
    # dimensional arrays. Correct this.
    # TODO: Not sure if this can still happen with metis 5
    if type(epart[0]) == np.ndarray:
        epart_new = np.zeros(len(epart), int)
        epart_new[:] = epart[:][0]
        epart = epart_new
        del epart_new


    triangles_per_proc = np.bincount(epart)

    msg = "Partition created where at least one submesh has no triangles. "
    msg += "Try using a smaller number of mpi processes."
    assert np.all(triangles_per_proc > 0), msg

    #proc_sum = np.zeros(n_procs+1, int)
    #proc_sum[1:] = np.cumsum(triangles_per_proc)

    epart_order = np.argsort(epart, kind='mergesort')

    return epart_order, triangles_per_proc


#==============================================================================================================
# Reverse Cuthill-McKee (RCM) reordering.
# Minimises the bandwidth of the mesh adjacency graph so that neighbours in
# the connectivity graph land close together in memory — directly targeting
# the access pattern of the flux kernel.
#==============================================================================================================

def rcm_partition(domain, n_procs=1):
    """Reverse Cuthill-McKee ordering split into n_procs equal contiguous blocks.

    Applies scipy's RCM algorithm to the full triangle adjacency graph, then
    divides the resulting ordering into n_procs equal-sized blocks (the same
    split strategy as morton_partition and hilbert_partition).  For n_procs=1
    this is a pure reordering; for n_procs>1 it can be used as the partition
    scheme in distribute().

    Parameters
    ----------
    domain : Domain
    n_procs : int
        Number of equal-sized blocks to divide the RCM ordering into.

    Returns
    -------
    order : ndarray of int, shape (N,)
        Permutation that reorders triangles by RCM then splits into blocks.
    triangles_per_proc : ndarray of int, shape (n_procs,)
        Number of triangles in each block.
    """
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import reverse_cuthill_mckee

    N = domain.number_of_triangles
    nbrs = domain.neighbours          # (N, 3), negative = boundary

    mask = (nbrs >= 0).ravel()
    rows = np.repeat(np.arange(N, dtype=np.int32), 3)[mask]
    cols = nbrs.ravel().astype(np.int32)[mask]
    A = csr_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)), shape=(N, N))

    order = reverse_cuthill_mckee(A, symmetric_mode=True)

    triangles_per_proc = np.full(n_procs, N // n_procs, dtype=int)
    triangles_per_proc[:N % n_procs] += 1

    return order, triangles_per_proc


def _rcm_within_partition(part_tris, all_neighbours, N):
    """Apply RCM to the sub-graph induced by part_tris.

    Parameters
    ----------
    part_tris : ndarray of int
        Global triangle indices in this partition.
    all_neighbours : ndarray, shape (N, 3)
        Full domain neighbour array (negative = boundary).
    N : int
        Total number of triangles (for sizing the lookup array).

    Returns
    -------
    ndarray of int
        part_tris reordered by RCM within the sub-graph.
    """
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import reverse_cuthill_mckee

    n = len(part_tris)
    if n <= 1:
        return part_tris

    # Map global index → local 0..n-1 (-1 means not in partition)
    local_of = np.full(N, -1, dtype=np.int32)
    local_of[part_tris] = np.arange(n, dtype=np.int32)

    # Extract neighbour rows for this partition; keep only internal edges
    nbrs = all_neighbours[part_tris]          # (n, 3)
    local_i = np.repeat(np.arange(n, dtype=np.int32), 3)
    global_j = nbrs.ravel()

    valid = global_j >= 0
    gj = global_j[valid]
    li = local_i[valid]
    lj = local_of[gj]
    internal = lj >= 0

    r, c = li[internal], lj[internal]
    if len(r) == 0:
        return part_tris

    A = csr_matrix((np.ones(len(r), dtype=np.int8), (r, c)), shape=(n, n))
    local_order = reverse_cuthill_mckee(A, symmetric_mode=True)
    return part_tris[local_order]


#==============================================================================================================
# Code for computing Morton (Z-order) codes for 2D points,
# used for spatial locality-preserving ordering of mesh elements.
# This is based on the "Bit Twiddling Hacks" by Sean Eron Anderson:
# https://graphics.stanford.edu/~seander/bithacks.html#InterleaveTables
#==============================================================================================================


def morton_partition(domain, n_procs):
    """
    Partition a mesh using Morton (Z-order) codes.

    This function computes Morton codes for the centroids of the triangles
    in the domain, sorts the triangles by these codes, and then divides them
    into contiguous blocks for each processor. This is a simple spatial
    partitioning method that can improve cache locality in parallel computations.

    Parameters
    ----------
    domain : object
        Domain object containing mesh information including number of triangles,
        nodes, triangle connectivity, and neighbor information.
    n_procs : int
        Number of processors to partition the mesh across.

    Returns
    -------
    epart_order : ndarray
        Integer array of indices that sort elements by Morton code order.
        Can be used to reorder triangles such that all triangles assigned to
        the same processor are contiguous.
    triangles_per_proc : ndarray
        Integer array of length n_procs where triangles_per_proc[i] is the
        number of triangles assigned to processor i.

    Raises
    ------
    ValueError
        If n_procs is greater than the number of triangles in the domain.
    """

    n_tri = domain.number_of_triangles

    if n_procs > n_tri:
        raise ValueError("Number of processors must be less than or equal to the number of triangles")


    points = domain.centroid_coordinates
    order = morton_order_from_points(points)

    # Now we have an ordering of the triangles based on their centroids' Morton codes.
    # We can divide this ordering into contiguous blocks for each processor.

    triangles_per_proc = np.full(n_procs, n_tri // n_procs, dtype=int)
    remainder = n_tri % n_procs
    if remainder > 0:
        triangles_per_proc[:remainder] += 1

    assert triangles_per_proc.sum() == n_tri, "Total number of triangles must match after partitioning"

    return order, triangles_per_proc

def morton_encode_2d(x, y):
    """
    Encode 2D coordinates to Morton codes.
    Uses 32 bits per dimension for 64-bit Morton codes.

    Parameters
    ----------
    x: np.ndarray of floats
    y: np.ndarray of floats

    Returns
    -------
    codes: np.ndarray
        Morton codes as uint64
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)

    # Ensure x and y have the same shape
    x, y = np.broadcast_arrays(x, y)

    min_x, max_x = np.min(x), np.max(x)
    min_y, max_y = np.min(y), np.max(y)

    # Add small padding to avoid edge cases
    x_range = max_x - min_x
    y_range = max_y - min_y
    padding = 1e-12  # Smaller padding for maximum resolution

    if x_range == 0:
        x_range = 1.0
        padding = 0.5
    if y_range == 0:
        y_range = 1.0
        padding = 0.5

    min_x -= x_range * padding
    max_x += x_range * padding
    min_y -= y_range * padding
    max_y += y_range * padding

    # Scale coordinates to [0, 2^32 - 1] for maximum resolution
    max_coord = 0xFFFFFFFF  # 2^32 - 1

    x_scaled = ((x - min_x) / (max_x - min_x) * max_coord).astype(np.uint64)
    y_scaled = ((y - min_y) / (max_y - min_y) * max_coord).astype(np.uint64)

    # Clamp to valid range
    x_scaled = np.clip(x_scaled, 0, max_coord)
    y_scaled = np.clip(y_scaled, 0, max_coord)

    # Dilate bits for 32-bit values - maximum resolution bit interleaving
    def dilate_bits_32(vals):
        """Dilate 32-bit values by inserting zeros between bits"""
        vals = vals & 0xFFFFFFFF  # Ensure 32-bit
        vals = (vals | (vals << 32)) & 0x00000000FFFFFFFF
        vals = (vals | (vals << 16)) & 0x0000FFFF0000FFFF
        vals = (vals | (vals << 8))  & 0x00FF00FF00FF00FF
        vals = (vals | (vals << 4))  & 0x0F0F0F0F0F0F0F0F
        vals = (vals | (vals << 2))  & 0x3333333333333333
        vals = (vals | (vals << 1))  & 0x5555555555555555
        return vals

    # Dilate and interleave: x gets odd positions, y gets even positions
    x_dilated = dilate_bits_32(x_scaled)
    y_dilated = dilate_bits_32(y_scaled)

    # Interleave by shifting y left by 1 and OR-ing with x
    morton_codes = x_dilated | (y_dilated << 1)

    return morton_codes


def morton_order_from_points(points):
    """
    Return Morton ordering of triangles from point coordinates.

    The points are first normalized to the unit square based on their
    bounding box, then mapped to an integer grid and a 2D Morton code
    is computed for each point. The resulting permutation is spatially
    locality-preserving.

    Parameters
    ----------
    points : ndarray, shape (N, 2)
        Points as (x, y) coordinates.

    Returns
    -------
    order : ndarray of int, shape (N,)
        Indices that sort the points in Morton (Z-order) order.

    Raises
    ------
    ValueError
        If points is not a 2D array with shape (N, 2).

    Notes
    -----
    This function is well suited for reordering per-element quantities
    in a finite-volume or finite-element code to improve cache locality.
    For example, for an ANUGA domain, you can pass
    ``domain.centroid_coordinates`` as ``points`` and then reorder
    per-triangle arrays with the returned permutation.
    """
    points = np.asarray(points, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 2:
        raise ValueError("points must have shape (N, 2)")

    codes = morton_encode_2d(points[:, 0], points[:, 1])
    order = np.argsort(codes, kind="mergesort")
    return order

#==============================================================================================================
# Code for computing Hilbert (space-filling curve) codes for 2D points,
# used for spatial locality-preserving ordering of mesh elements.
# This is based on the Hilbert curve algorithms described by John Skilling:
# https://www.johndcook.com/blog/2018/10/16/hilbert-curve-algorithms/
#==============================================================================================================

def hilbert_partition(domain, n_procs):
    """
    Partition a mesh using Hilbert (space-filling curve) codes.

    This function computes Hilbert codes for the centroids of the triangles
    in the domain, sorts the triangles by these codes, and then divides them
    into contiguous blocks for each processor. This is a simple spatial
    partitioning method that can improve cache locality in parallel computations.

    Parameters
    ----------
    domain : object
        Domain object containing mesh information including number of triangles,
        nodes, triangle connectivity, and neighbor information.
    n_procs : int
        Number of processors to partition the mesh across.

    Returns
    -------
    epart_order : ndarray
        Integer array of indices that sort elements by Hilbert code order.
        Can be used to reorder triangles such that all triangles assigned to
        the same processor are contiguous.
    triangles_per_proc : ndarray
        Integer array of length n_procs where triangles_per_proc[i] is the
        number of triangles assigned to processor i.

    Raises
    ------
    ValueError
        If n_procs is greater than the number of triangles in the domain.
    """

    n_tri = domain.number_of_triangles

    if n_procs > n_tri:
        raise ValueError("Number of processors must be less than or equal to the number of triangles")


    points = domain.centroid_coordinates
    order = hilbert_order_from_points(points)

    # Now we have an ordering of the triangles based on their centroids' Hilbert codes.
    # We can divide this ordering into contiguous blocks for each processor.

    triangles_per_proc = np.full(n_procs, n_tri // n_procs, dtype=int)
    remainder = n_tri % n_procs
    if remainder > 0:
        triangles_per_proc[:remainder] += 1

    assert triangles_per_proc.sum() == n_tri, "Total number of triangles must match after partitioning"

    return order, triangles_per_proc

def hilbert_index_2d(x, y, p):
    """Integer 2D Hilbert indices for coords in [0, 2**p - 1]."""
    x = np.asarray(x, dtype=np.uint64)
    y = np.asarray(y, dtype=np.uint64)
    x, y = np.broadcast_arrays(x, y)
    x = x.copy()
    y = y.copy()

    M = np.uint64(1) << (p - 1)
    Q = M
    while Q > 1:
        P = Q - 1
        mask = (x & Q) != 0
        y ^= P * mask
        mask = (y & Q) != 0
        tmp = (x ^ y) & (P * mask)
        x ^= tmp
        y ^= tmp
        Q >>= 1

    t = (x ^ y) >> 1
    y ^= x
    x ^= t

    h = np.zeros_like(x, dtype=np.uint64)
    for i in range(p):
        bit = np.uint64(1) << i
        h |= ((x & bit) << i) | ((y & bit) << (i + 1))

    return h.ravel()


def hilbert_codes_from_points(points, p=16):
    """Raw Hilbert codes for 2D float points (no argsort).

    Parameters
    ----------
    points : array_like, shape (N, 2)
    p : int
        Bits per axis; grid is 2**p by 2**p.

    Returns
    -------
    h : ndarray of uint64, shape (N,)
    """
    pts = np.asarray(points, dtype=np.float64)
    if pts.ndim != 2 or pts.shape[1] != 2:
        raise ValueError("points must have shape (N, 2)")
    mins = pts.min(axis=0)
    maxs = pts.max(axis=0)
    span = maxs - mins
    span[span == 0.0] = 1.0
    norm = (pts - mins) / span
    max_int = (1 << p) - 1
    ij = np.clip((norm * max_int).astype(np.uint64), 0, max_int)
    return hilbert_index_2d(ij[:, 0], ij[:, 1], p)


def hilbert_order_from_points(points, p=16):
    """Hilbert order for 2D float points via normalization to 2**p grid.

    Parameters
    ----------
    points : array_like, shape (N, 2)
        2D coordinates (e.g. centroids) as floats.
    p : int
        Bits per axis; grid is 2**p by 2**p.

    Returns
    -------
    order : ndarray of int, shape (N,)
        Indices that sort points along the 2D Hilbert curve.
    """
    h = hilbert_codes_from_points(points, p)
    return np.argsort(h, kind="mergesort")


def reorder_domain(domain, method='hilbert', n_procs=None, verbose=False):
    """Reorder domain triangles for cache locality.

    Computes an ordering using a space-filling curve or Metis partitioning and
    applies it to the domain in-place.  Call after distribute() and before
    create_riverwalls() / operator setup so those indices are built on top of
    the reordered mesh.

    Parameters
    ----------
    domain : Domain
    method : str
        'hilbert' (default), 'morton', or 'metis'
    n_procs : int or None
        Number of partitions for Metis ordering.  If None, reads
        ``OMP_NUM_THREADS`` from the environment (defaulting to 1 when unset).
        For Hilbert/Morton this parameter is ignored.
    verbose : bool

    Returns
    -------
    domain : the same Domain object, reordered in-place
    """
    import os

    if n_procs is None:
        n_procs = int(os.environ.get('OMP_NUM_THREADS', 1))

    if method == 'hilbert':
        epart_order, _ = hilbert_partition(domain, 1)
    elif method == 'morton':
        epart_order, _ = morton_partition(domain, 1)
    elif method == 'rcm':
        epart_order, _ = rcm_partition(domain)
    elif method == 'metis':
        epart_order, _ = metis_partition(domain, max(1, n_procs))
    elif method == 'metis_hilbert':
        metis_order, triangles_per_proc = metis_partition(domain, max(1, n_procs))
        h = hilbert_codes_from_points(domain.centroid_coordinates)
        chunks = []
        offset = 0
        for count in triangles_per_proc:
            part_tris = metis_order[offset:offset + count]
            chunks.append(part_tris[np.argsort(h[part_tris], kind='mergesort')])
            offset += count
        epart_order = np.concatenate(chunks)
    else:  # metis_rcm
        metis_order, triangles_per_proc = metis_partition(domain, max(1, n_procs))
        N = domain.number_of_triangles
        chunks = []
        offset = 0
        for count in triangles_per_proc:
            part_tris = metis_order[offset:offset + count]
            chunks.append(_rcm_within_partition(part_tris, domain.neighbours, N))
            offset += count
        epart_order = np.concatenate(chunks)

    if verbose:
        if method in ('metis', 'metis_hilbert', 'metis_rcm'):
            print(f'reorder_domain: reordering {len(epart_order)} triangles '
                  f'using {method} ordering with {n_procs} partitions')
        else:
            print(f'reorder_domain: reordering {len(epart_order)} triangles '
                  f'using {method} ordering')

    domain.reorder(epart_order)
    return domain

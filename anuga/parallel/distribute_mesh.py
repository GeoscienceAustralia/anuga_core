

#=========================================================================
#
#  Read in a data file and subdivide the triangle list
#
#
#  The final routine, pmesh_divide_metis, does automatic
# grid partitioning. Once testing has finished on this
# routine the others should be removed.
#
#  FIXME SR: This module was built when ANUGA was using
#  vertex values. We now use centroid values. Really only
#  to work with centroid arrays.
#
#  Authors: Linda Stals and Matthew Hardy, June 2005
#  Modified: Linda Stals, Nov 2005
#            Jack Kelly, Nov 2005
#            Steve Roberts, Aug 2009 (updating to numpy)
#
#=========================================================================
import sys
from os import sep
from sys import path
from math import floor

import numpy as num

try:
    import numpy.lib.arraysetops as numset
except ImportError:
    import numpy as numset


from anuga.abstract_2d_finite_volumes.neighbour_mesh import Mesh
from anuga.abstract_2d_finite_volumes.general_mesh import General_mesh
from anuga import indent

try:
    import local_config as config
except ImportError:
    from . import config as config



verbose = False

# The keys of the quantities that are to be saved and reordered.
DEFAULT_DISTRIBUTE_QUANTITY_NAMES = ["stage", "xmomentum", "ymomentum", "elevation", "friction"]

# Tracers are deliberately not Quantity objects (#276), so nothing that walks
# domain.quantities sees them -- which is why distribute() used to hand the
# sub-domains no tracers at all (#278). Rather than build a second partitioning
# path for them, their centroid values ride the existing one under reserved
# keys: the same reordering, ghost layer and MPI transfer, for free.
#
# Keyed by SLOT, not by name. The order is what carries meaning (class s must
# land back in slot s, which add_sediment_fraction relies on), and a
# tracer name is free text that has no business in a dict key. The names travel
# beside the values, in the metadata below.
TRACER_QUANTITY_PREFIX = '__tracer__'

# The names and the shared beta, carried in the kwargs dict that already
# travels with each partition. Popped before the Domain constructor sees it.
# Old partition files simply lack the key, which reads as "no tracers".
TRACER_METADATA_KEY = '__tracer_metadata__'


def tracer_partition_quantities(domain):
    """Reserved-key entries carrying each tracer's concentration `c`.

    `c` rather than the conserved `m = h*c`: the sub-domain reseeds through
    set_tracer, which recomputes m from its OWN depth, so the two stay
    consistent even though h is rebuilt from the partitioned stage/elevation.
    """
    ns = getattr(domain, 'number_of_tracers', 0)
    if not ns:
        return {}
    return {'%s%d' % (TRACER_QUANTITY_PREFIX, s): domain.tracer_centroid_values[s]
            for s in range(ns)}


def tracer_partition_metadata(domain):
    """The tracer description to carry alongside the values, or None."""
    ns = getattr(domain, 'number_of_tracers', 0)
    if not ns:
        return None
    return {'names': list(domain.get_tracer_names()),
            'beta': float(domain.beta_tracer)}


def pop_tracer_quantities(quantities):
    """Take the reserved tracer entries out of a partition's quantity dict.

    Must be called BEFORE the set_quantity loop: left in, each one would be
    created as a Quantity named `__tracer__0`, which is exactly the confusion
    keeping tracers out of domain.quantities is meant to avoid.

    Returns the concentration arrays in slot order, or None if there are none.
    """
    keys = [k for k in quantities if k.startswith(TRACER_QUANTITY_PREFIX)]
    if not keys:
        return None

    slots = sorted(int(k[len(TRACER_QUANTITY_PREFIX):]) for k in keys)
    if slots != list(range(len(slots))):
        raise ValueError(
            'tracer slots in this partition are not 0..n-1: %r -- the tracer '
            'order would not survive the partition' % (slots,))

    return [num.asarray(quantities.pop('%s%d' % (TRACER_QUANTITY_PREFIX, s))
                        ).reshape(-1)
            for s in slots]


def restore_tracers(domain, metadata, values):
    """Re-create the tracers on a sub-domain and seed them (#278).

    Call AFTER the quantities are set: set_tracer derives m = h*c from the
    sub-domain's own depth, so stage and elevation have to be in place first.

    `values` covers full AND ghost triangles, because the tracer arrays rode
    the same partitioning as the quantities -- so the halo starts correct
    rather than waiting for the first communicate_tracer_ghosts().
    """
    if metadata is None and values is None:
        return
    if metadata is None or values is None:
        raise ValueError(
            'this partition carries tracer %s but not the other; it cannot be '
            'reconstructed' % ('values' if metadata is None else 'metadata'))

    names = metadata['names']
    if len(names) != len(values):
        raise ValueError(
            'partition carries %d tracer name(s) but %d value array(s)'
            % (len(names), len(values)))

    beta = metadata['beta']
    for s, name in enumerate(names):
        # initial_value=None: seeded from the partitioned values just below,
        # rather than filled and then immediately overwritten.
        index = domain.add_tracer(name, beta=beta, initial_value=None)
        assert index == s, 'tracer %r landed in slot %d, expected %d' % (
            name, index, s)
        domain.set_tracer(name, values[s])

#=========================================================================
#
# If the triangles list is reordered, the quantities
# assigned to the triangles must also be reorded.
#
# *) quantities contain the quantites in the old ordering
# *) proc_sum[i] contains the number of triangles in
# processor i
# *) tri_index is a map from the old triangle ordering to
# the new ordering, where the new number for triangle
# i is proc_sum[tri_index[i][0]]+tri_index[i][1]
#
# -------------------------------------------------------
#
# *) The quantaties are returned in the new ordering
#
#=========================================================================#


# def reorder(quantities, tri_index, proc_sum):

#     # Find the number triangles

#     N = len(tri_index)

#     # Temporary storage area

#     index = num.zeros(N, int)
#     q_reord = {}

#     # Find the new ordering of the triangles

#     for i in range(N):
#         bin = tri_index[i][0]
#         bin_off_set = tri_index[i][1]
#         index[i] = proc_sum[bin]+bin_off_set

#     # Reorder each quantity according to the new ordering

#     for k in quantities:
#         q_reord[k] = num.zeros((N, 1), float)
#         for i in range(N):
#             q_reord[k][index[i]] = quantities[k].centroid_values[i]
#     del index

#     return q_reord


def reorder_quantities(quantities, epart_order):

    # Find the number triangles

    N = len(epart_order)

    # Temporary storage area

    q_reord = {}

    # Reorder each quantity according to the new ordering

    for k in quantities:
        q = quantities[k]
        # Quantity objects for the real quantities; plain centroid arrays for
        # the reserved tracer entries, which are not Quantity objects (#276).
        values = q.centroid_values if hasattr(q, 'centroid_values') else num.asarray(q)
        q_reord[k] = num.zeros((N, 1), float)
        q_reord[k][:,0] = values[epart_order]

    return q_reord





def partition_mesh(domain, n_procs,
                   parameters=None,
                   verbose=False):
    """Partition a mesh across multiple processors using METIS, Morton, or Hilbert partitioning.

    This function takes a serial mesh and distributes it across multiple processors
    by reordering triangles according to a partitioning scheme. It generates
    mappings between serial and parallel triangle indices and reorders mesh quantities
    accordingly.

    Parameters
    ----------
    domain : Domain
        The serial domain object containing the mesh and quantities to be partitioned.
    n_procs : int
        Number of processors to partition the mesh across.
    parameters : dict, optional
        Additional parameters for the partitioning algorithm. Supported keys:

        - 'partition_scheme' : str, default 'metis'
            Partitioning scheme ('metis', 'morton', or 'hilbert')
        - 'distribute_quantity_names' : list, default DEFAULT_DISTRIBUTE_QUANTITY_NAMES
            Names of quantities to distribute
        - 'in_place' : bool, default False
            If True, reorder mesh in place to save memory
        - 'verbose' : bool, default False
            If True, print verbose output during partitioning

    verbose : bool, optional
        If True, print verbose output during partitioning, by default False.

    Returns
    -------
    new_mesh : Mesh
        The reordered mesh distributed across processors.
    triangles_per_proc : ndarray
        Array of shape (n_procs,) containing the number of triangles assigned
        to each processor.
    new_quantities : dict
        Dictionary of reordered quantities corresponding to the new mesh ordering.
    new_tri_index : ndarray
        Array of shape (n_triangles, 2) mapping serial triangle indices to
        (processor_id, local_index) pairs.
    epart_order : ndarray
        Array of reordering indices used to partition the mesh from serial to
        parallel ordering.

    Notes
    -----
    The function uses METIS, Morton, or Hilbert partitioning to distribute triangles
    across processors. Set `in_place=True` in parameters if the original sequential
    domain will not be used again to save memory.

    Examples
    --------
    >>> new_mesh, tpp, quantities, tri_idx, epart = partition_mesh(
    ...     domain, 4, parameters={'partition_scheme': 'metis'})
    """

    # Set default parameters
    distribute_quantity_names = DEFAULT_DISTRIBUTE_QUANTITY_NAMES
    in_place = False
    partition_scheme = 'metis'

    # Override defaults with user-provided parameters
    if parameters is None:
        parameters = {}

    if 'partition_scheme' in parameters:
        partition_scheme = parameters['partition_scheme']
    if 'distribute_quantity_names' in parameters:
        distribute_quantity_names = parameters['distribute_quantity_names']
    if 'in_place' in parameters:
        in_place = bool(parameters['in_place'])
    if 'verbose' in parameters:
        verbose = bool(parameters['verbose'])

    # Output partitioning parameters if verbose
    if verbose:
        print("partition_mesh: Number of processors: {}".format(n_procs))
        print("partition_mesh: In-place reordering: {}".format(in_place))
        print("partition_mesh: Quantities to distribute: {}".format(distribute_quantity_names))

    # Serial to Parallel and Parallel to Serial Triangle index maps
    tri_index = {}
    r_tri_index = {}  # reverse tri index, parallel to serial triangle index mapping
    n_tri = len(domain.triangles)

    # Support both a full Domain (has .quantities and .mesh) and a
    # Basic_mesh (no quantities; reorder is called on the object itself).
    if hasattr(domain, 'quantities'):
        distribute_quantities = {k: domain.quantities[k]
                                 for k in distribute_quantity_names
                                 if k in domain.quantities}
    else:
        distribute_quantities = {}

    # Tracers are not in domain.quantities, so the comprehension above cannot
    # see them; add them under reserved keys so they ride the same machinery.
    distribute_quantities.update(tracer_partition_quantities(domain))

    from anuga.parallel.partitioning import metis_partition, morton_partition, hilbert_partition, rcm_partition

    if verbose: print("partition_mesh: Computing partitioning using {}...".format(partition_scheme))
    if partition_scheme == 'morton':
        epart_order, triangles_per_proc = morton_partition(domain, n_procs)
    elif partition_scheme == 'hilbert':
        epart_order, triangles_per_proc = hilbert_partition(domain, n_procs)
    elif partition_scheme == 'rcm':
        epart_order, triangles_per_proc = rcm_partition(domain, n_procs)
    else:
        epart_order, triangles_per_proc = metis_partition(domain, n_procs)

    # Build processor IDs and local indices in a fully vectorized way
    proc_ids = num.repeat(num.arange(n_procs, dtype=int), triangles_per_proc)
    local_ids = num.concatenate(
        [num.arange(k, dtype=int) for k in triangles_per_proc]
    )

    new_tri_index = num.empty((n_tri, 2), dtype=int)
    new_tri_index[epart_order] = num.column_stack((proc_ids, local_ids))

    new_quantities = reorder_quantities(distribute_quantities, epart_order)

    # If you are just distributing the sequential domain and will not be using
    # it again, then some memory can be saved by setting in_place = true.
    # For a Basic_mesh there is no separate .mesh attribute -- reorder directly.
    mesh_obj = getattr(domain, 'mesh', domain)
    new_mesh = mesh_obj.reorder(epart_order, in_place=in_place)

    return new_mesh, triangles_per_proc, new_quantities, new_tri_index, epart_order

def partition_mesh_without_map(domain, n_procs):
    # Wrapper of partition_mesh which does not return tri_index or r_tri_index

    mesh, triangles_per_proc, quantities, tri_index, r_tri_index = partition_mesh(
        domain, n_procs)

    return mesh, triangles_per_proc, quantities

#=========================================================================
#
# Subdivide the domain. This module is primarily
# responsible for building the ghost layer and
# communication pattern
#
#
#  Author: Linda Stals, June 2005
#  Modified: Linda Stals, Nov 2005 (optimise python code)
#            Steve Roberts, Aug 2009 (convert to numpy)
#
#
#=========================================================================


#=========================================================================
#
# Subdivide the triangles into non-overlapping domains.
#
#  *)  The subdivision is controlled by triangles_per_proc.
# The first triangles_per_proc[0] triangles are assigned
# to the first processor, the second triangles_per_proc[1]
# are assigned to the second processor etc.
#
#  *) nodes, triangles and boundary contains all of the
# nodes, triangles and boundary tag information for the
# whole domain. The triangles should be orientated in the
# correct way and the nodes number consecutively from 0.
#
# -------------------------------------------------------
#
#  *) A dictionary containing the full_nodes, full_triangles
# and full_boundary information for each processor is
# returned. The node information consists of
# [global_id, x_coord, y_coord].
#
#=========================================================================

def submesh_full(mesh, triangles_per_proc):

    nodes = mesh.nodes
    triangles = mesh.triangles
    boundary = mesh.boundary

    tlower = 0
    nproc = len(triangles_per_proc)
    node_list = []
    triangle_list = []
    boundary_list = []
    submesh = {}

    # Pre-sort boundary keys by triangle id so each processor's slice can be
    # found with a binary search — O(|boundary| log|boundary| + P log|boundary|)
    # instead of the naive O(P x |boundary|).
    sorted_bnd_keys = sorted(boundary.keys(), key=lambda k: k[0])
    sorted_bnd_tri_ids = num.array([k[0] for k in sorted_bnd_keys], dtype=int)

    # Loop over processors
    for p in range(nproc):

        # Find triangles on processor p
        tupper = triangles_per_proc[p] + tlower
        subtriangles = triangles[tlower:tupper]
        triangle_list.append(subtriangles)

        # Find the boundary edges on processor p using binary search
        lo = int(num.searchsorted(sorted_bnd_tri_ids, tlower))
        hi = int(num.searchsorted(sorted_bnd_tri_ids, tupper))
        subboundary = {sorted_bnd_keys[i]: boundary[sorted_bnd_keys[i]]
                       for i in range(lo, hi)}
        boundary_list.append(subboundary)

        # Find nodes in processor p
        ids = num.unique(subtriangles.flat)
        lnodes = nodes[ids]
        x = num.concatenate((num.reshape(ids, (-1, 1)), lnodes), 1)
        node_list.append(x)

        # Move to the next processor
        tlower = tupper

    # Put the results in a dictionary

    submesh["full_nodes"] = node_list
    submesh["full_triangles"] = triangle_list
    submesh["full_boundary"] = boundary_list

    return submesh


#=========================================================================
#
# Build the ghost layer of triangles
#
#  *) Given the triangle subpartion for the processor
# build a ghost layer of triangles. The ghost layer
# consists of two layers of neighbouring triangles.
#
#  *) The vertices in the ghost triangles must also
# be added to the node list for the current processor
#
#
# -------------------------------------------------------
#
#  *) The extra triangles and nodes are returned.
#
#  *)  The node information consists of
# [global_id, x_coord, y_coord].
#
#  *) The triangle information consists of
# [triangle number, t], where t = [v1, v2, v3].
#
#=========================================================================

def ghost_layer(submesh, mesh, p, tupper, tlower, parameters=None):

    layer_width = (parameters or {}).get('ghost_layer_width', 2)

    from .distribute_mesh_ext import ghost_layer_bfs

    # BFS over the global neighbour array — O(ghost_count) instead of
    # O(ghost_count * log(ghost_count)) set operations per layer.
    neighbours = num.ascontiguousarray(mesh.neighbours, dtype=num.int64)
    new_trianglemap = ghost_layer_bfs(neighbours, tlower, tupper, layer_width)

    new_subtriangles = num.concatenate(
        (num.reshape(new_trianglemap, (-1, 1)), mesh.triangles[new_trianglemap]), 1)

    fullnodes = submesh["full_nodes"][p]
    full_nodes_ids = num.array(fullnodes[:, 0], int)

    new_nodes = num.unique(mesh.triangles[new_trianglemap].flat)
    new_nodes = numset.setdiff1d(new_nodes, full_nodes_ids)

    new_subnodes = num.concatenate(
        (num.reshape(new_nodes, (-1, 1)), mesh.nodes[new_nodes]), 1)

    return new_subnodes, new_subtriangles, layer_width

#=========================================================================
#
# Find the edges of the ghost trianlges that do not
# have a neighbour in the current cell. These are
# treated as a special type of boundary edge.
#
#  *) Given the ghost triangles in a particular
# triangle, use the mesh to find its neigbours. If
# the neighbour is not in the processor set it to
# be a boundary edge
#
#  *) The vertices in the ghost triangles must also
# be added to the node list for the current processor
#
#  *) The boundary edges for the ghost triangles are
# ignored.
#
# -------------------------------------------------------
#
#  *) The type assigned to the ghost boundary edges is 'ghost'
#
#  *)  The boundary information is returned as a directorier
# with the key = (triangle id, edge no) and the values
# assigned to the key is 'ghost'
#
#
#=========================================================================


def ghost_bnd_layer(ghosttri, tlower, tupper, mesh, p):

    from .distribute_mesh_ext import ghost_bnd_layer_classify

    boundary = mesh.boundary

    ghost_ids = num.ascontiguousarray(ghosttri[:, 0], dtype=num.int64)
    neighbours = num.ascontiguousarray(mesh.neighbours, dtype=num.int64)

    tri_ids, edge_ids = ghost_bnd_layer_classify(
        neighbours, ghost_ids, int(tlower), int(tupper))

    # All classified edges start as 'ghost'; real boundary tags override.
    subboundary = {(int(t), int(e)): 'ghost'
                   for t, e in zip(tri_ids, edge_ids)}
    subboundary.update((k, boundary[k])
                       for k in subboundary.keys() & boundary.keys())

    return subboundary

#=======================================================================
#
# The ghost triangles on the current processor will need
# to get updated information from the neighbouring
# processor containing the corresponding full triangles.
#
#  *) The tri_per_proc is used to determine which
# processor contains the full node copy.
#
# -------------------------------------------------------
#
#  *) The ghost communication pattern consists of
# [global node number, neighbour processor number].
#
#=========================================================================


def ghost_commun_pattern(subtri, p, tri_per_proc_range):

    global_no = num.reshape(subtri[:, 0], (-1, 1))
    neigh = num.reshape(num.searchsorted(
        tri_per_proc_range, global_no), (-1, 1))

    ghost_commun = num.concatenate((global_no, neigh), axis=1)

    return ghost_commun

#=========================================================================
#
# The full triangles in this processor must communicate
# updated information to neighbouring processor that
# contain ghost triangles
#
#  *) The ghost communication pattern for all of the
# processor must be built before calling this processor.
#
#  *) The full communication pattern is found by looping
# through the ghost communication pattern for all of the
# processors. Recall that this information is stored in
# the form [global node number, neighbour processor number].
# The full communication for the neighbour processor is
# then updated.
#
# -------------------------------------------------------
#
#  *) The full communication pattern consists of
# [global id, [p1, p2, ...]], where p1, p2 etc contain
# a ghost node copy of the triangle global id.
#
#=========================================================================


def full_commun_pattern(submesh, tri_per_proc):
    nproc = len(tri_per_proc)

    # Sparse dicts: only triangles that are ghost-copied somewhere get an entry.
    # The old O(N_triangles) pre-fill loop is not needed because build_local_commun
    # only iterates over keys that have non-empty lists.
    full_commun = [{} for _ in range(nproc)]

    if nproc == 1:
        return full_commun

    # Stack all ghost_commun arrays and record the ghost-side processor.
    # ghost_commun[p] has shape (G_p, 2): col0 = global_id, col1 = owner_proc.
    gc_list = submesh["ghost_commun"]
    parts = []
    for p in range(nproc):
        gc = gc_list[p]
        if len(gc) == 0:
            continue
        ghost_col = num.full(len(gc), p, dtype=int)
        parts.append(num.column_stack([gc, ghost_col]))   # (G_p, 3)

    if not parts:
        return full_commun

    # all_gc columns: [global_id, owner_proc, ghost_proc]
    all_gc      = num.vstack(parts)
    global_ids  = all_gc[:, 0]
    owner_procs = all_gc[:, 1]
    ghost_procs = all_gc[:, 2]

    # Sort by owner_proc (primary) then global_id (secondary) so each
    # unique (owner, global_id) pair forms a contiguous run.
    sort_idx = num.lexsort((global_ids, owner_procs))
    s_global = global_ids[sort_idx]
    s_owner  = owner_procs[sort_idx]
    s_ghost  = ghost_procs[sort_idx]

    # Locate the start of every new (owner, global_id) group.
    pair_changes = num.empty(len(s_global), dtype=bool)
    pair_changes[0] = True
    pair_changes[1:] = (
        (s_global[1:] != s_global[:-1]) | (s_owner[1:] != s_owner[:-1])
    )
    starts = num.where(pair_changes)[0]
    ends   = num.empty_like(starts)
    ends[:-1] = starts[1:]
    ends[-1]  = len(s_global)

    for lo, hi in zip(starts, ends):
        owner = int(s_owner[lo])
        gid   = int(s_global[lo])
        full_commun[owner][gid] = s_ghost[lo:hi].tolist()

    return full_commun


#=========================================================================
#
# Given the non-overlapping grid partition, an extra layer
# of triangles are included to help with the computations.
# The triangles in this extra layer are not updated by
# the processor, their updated values must be sent by the
# processor containing the original, full, copy of the
# triangle. The communication pattern that controls these
# updates must also be built.
#
#  *) Assumes that full triangles, nodes etc have already
# been found and stored in submesh
#
#  *) See the documentation for ghost_layer,
# ghost_commun_pattern and full_commun_pattern
#
# -------------------------------------------------------
#
#  *) The additional information is added to the submesh
# dictionary. See the documentation for ghost_layer,
# ghost_commun_pattern and full_commun_pattern
#
#  *) The ghost_triangles, ghost_nodes, ghost_boundary,
# ghost_commun and full_commun is added to submesh
#=========================================================================

def submesh_ghost(submesh, mesh, triangles_per_proc, parameters=None):

    from .distribute_mesh_ext import ghost_layer_bfs_all, ghost_bnd_layer_classify_all

    nproc = len(triangles_per_proc)
    layer_width = (parameters or {}).get('ghost_layer_width', 2)

    # Build per-processor range arrays for the batch Cython calls.
    cumsum      = num.concatenate([[0], num.cumsum(triangles_per_proc)])
    tlower_arr  = num.ascontiguousarray(cumsum[:-1], dtype=num.int64)
    tupper_arr  = num.ascontiguousarray(cumsum[1:],  dtype=num.int64)
    neighbours  = num.ascontiguousarray(mesh.neighbours, dtype=num.int64)

    # Step 1 — parallel BFS ghost-layer expansion (OpenMP prange over P procs).
    ghost_id_list = ghost_layer_bfs_all(
        neighbours, tlower_arr, tupper_arr, layer_width)

    # Step 2 — parallel ghost-boundary classification (OpenMP prange over P procs).
    bnd_results = ghost_bnd_layer_classify_all(
        neighbours, ghost_id_list, tlower_arr, tupper_arr)

    # Step 3 — assemble submesh entries (serial; all operations are cheap numpy).
    ghost_triangles    = []
    ghost_nodes        = []
    ghost_commun       = []
    ghost_bnd          = []
    ghost_layer_widths = []
    triangles_per_proc_ranges = num.cumsum(triangles_per_proc) - 1
    boundary = mesh.boundary

    for p in range(nproc):
        ghost_ids = ghost_id_list[p]
        tlo = int(tlower_arr[p])
        tup = int(tupper_arr[p])

        # Ghost triangles: [global_id, v0, v1, v2]
        subtri = num.concatenate(
            (num.reshape(ghost_ids, (-1, 1)), mesh.triangles[ghost_ids]), 1)

        # Ghost nodes: unique vertices touched by ghost triangles, minus full nodes.
        full_node_ids = num.array(submesh["full_nodes"][p][:, 0], int)
        new_nodes = numset.setdiff1d(num.unique(mesh.triangles[ghost_ids].flat),
                                     full_node_ids)
        subnodes = num.concatenate(
            (num.reshape(new_nodes, (-1, 1)), mesh.nodes[new_nodes]), 1)

        ghost_triangles.append(subtri)
        ghost_nodes.append(subnodes)
        ghost_layer_widths.append(layer_width)

        # Ghost boundary dict: start with 'ghost' tags, override with real tags.
        tri_ids, edge_ids = bnd_results[p]
        subbnd = {(int(t), int(e)): 'ghost'
                  for t, e in zip(tri_ids, edge_ids)}
        subbnd.update((k, boundary[k])
                      for k in subbnd.keys() & boundary.keys())
        ghost_bnd.append(subbnd)

        # Ghost communication pattern.
        gcommun = ghost_commun_pattern(subtri, p, triangles_per_proc_ranges)
        ghost_commun.append(gcommun)

    submesh["ghost_layer_width"] = ghost_layer_widths
    submesh["ghost_nodes"]       = ghost_nodes
    submesh["ghost_triangles"]   = ghost_triangles
    submesh["ghost_commun"]      = ghost_commun
    submesh["ghost_boundary"]    = ghost_bnd

    full_commun = full_commun_pattern(submesh, triangles_per_proc)
    submesh["full_commun"] = full_commun

    return submesh


#=========================================================================
#
# Certain quantities may be assigned to the triangles,
# these quantities must be subdivided in the same way
# as the triangles
#
#  *) The quantities are ordered in the same way as the
# triangles
#
# -------------------------------------------------------
#
#  *) The quantites attached to the full triangles are
# stored in full_quan
#
#  *) The quantities attached to the ghost triangles are
# stored in ghost_quan
#=========================================================================

def submesh_quantities(submesh, quantities, triangles_per_proc):

    nproc = len(triangles_per_proc)

    lower = 0

    # Build an empty dictionary to hold the quantites

    submesh["full_quan"] = {}
    submesh["ghost_quan"] = {}
    for k in quantities:
        submesh["full_quan"][k] = []
        submesh["ghost_quan"][k] = []

    # Loop through the subdomains

    for p in range(nproc):
        upper = lower + triangles_per_proc[p]

        # Global IDs of ghost triangles for processor p — column 0 of the
        # ghost_triangles array (shape M×4: [global_id, v0, v1, v2]).
        global_id = submesh["ghost_triangles"][p][:, 0]

        # Use the global IDs to extract quantity values for ghost triangles
        for k in quantities:
            submesh["full_quan"][k].append(quantities[k][lower:upper])
            submesh["ghost_quan"][k].append(quantities[k][global_id])

        lower = upper

    return submesh

#=========================================================================
#
# Build the grid partition on the host.
#
#  *) See the documentation for submesh_ghost and
# submesh_full
#
# -------------------------------------------------------
#
#  *) A dictionary containing the full_triangles,
# full_nodes, full_boundary, ghost_triangles, ghost_nodes,
# ghost_boundary, ghost_commun and full_commun and true boundary polygon is returned.
#
#=========================================================================

def _submesh_structure_hash(mesh, triangles_per_proc, parameters):
    """SHA-256 fingerprint of the mesh topology + partition + layer_width.

    Only mesh.triangles is hashed (the neighbour array is derived from it).
    For a 50 M-triangle mesh this costs ~1-2 s — negligible vs the minutes
    saved by avoiding a full rebuild.
    """
    import hashlib
    layer_width = (parameters or {}).get('ghost_layer_width', 2)
    h = hashlib.sha256()
    h.update(mesh.triangles.tobytes())
    h.update(triangles_per_proc.tobytes())
    h.update(str(layer_width).encode())
    return h.hexdigest()[:20]


def build_submesh(mesh, quantities,
                  triangles_per_proc, parameters=None, verbose=False):
    """Build the full submesh data structure for all processors.

    Parameters
    ----------
    mesh : Mesh
    quantities : dict
    triangles_per_proc : ndarray
    parameters : dict, optional
        Recognised keys (in addition to existing ones):

        ``cache_dir`` : str or path-like, optional
            Directory in which to cache the mesh *structure* (ghost layers,
            communication pattern) between runs.  The quantities are never
            cached because they may differ between runs.  Cache files are
            keyed by a SHA-256 hash of the mesh topology, partition, and
            ghost-layer width, so the cache is automatically invalidated
            whenever any of those change.

            Example::

                parameters = {'cache_dir': '.partition_cache'}

    verbose : bool, optional
    """
    import pathlib
    import pickle
    import hashlib

    cache_dir = None
    if parameters is not None:
        cache_dir = parameters.get('cache_dir', None)

    if cache_dir is not None:
        cache_dir = pathlib.Path(cache_dir)
        cache_dir.mkdir(parents=True, exist_ok=True)
        key = _submesh_structure_hash(mesh, triangles_per_proc, parameters)
        P = len(triangles_per_proc)
        cache_file = cache_dir / f'anuga_submesh_{key}_P{P}.pkl'

        if cache_file.exists():
            if verbose:
                print(f'build_submesh: loading cached structure from {cache_file}')
            with open(cache_file, 'rb') as f:
                submeshg = pickle.load(f)
            # Quantities may differ between runs — always recompute.
            return submesh_quantities(submeshg, quantities, triangles_per_proc)

    # Subdivide into non-overlapping partitions
    submeshf = submesh_full(mesh, triangles_per_proc)

    # Add ghost boundary layer and communication pattern
    submeshg = submesh_ghost(submeshf, mesh, triangles_per_proc, parameters)

    if cache_dir is not None:
        if verbose:
            print(f'build_submesh: saving structure to cache {cache_file}')
        with open(cache_file, 'wb') as f:
            pickle.dump(submeshg, f, protocol=pickle.HIGHEST_PROTOCOL)

    return submesh_quantities(submeshg, quantities, triangles_per_proc)

#=========================================================================
#
#  Given the subdivision of the grid assigned to the
# current processor convert it into a form that is
# appropriate for the GA datastructure.
#
#  The main function of these modules is to change the
# node numbering. The GA datastructure assumes they
# are numbered consecutively from 0.
#
#  The module also changes the communication pattern
# datastructure into a form needed by parallel_advection
#
#  Authors: Linda Stals and Matthew Hardy, June 2005
#  Modified: Linda Stals, Nov 2005 (optimise python code)
#            Steve Roberts, Aug 2009 (updating to numpy)
#
#
#=========================================================================


#=========================================================================#
# Convert the format of the data to that used by ANUGA
#
#
# *) Change the nodes global ID's to an integer value,
# starting from 0.
#
# *) The triangles and boundary edges must also be
# updated accordingly.
#
# -------------------------------------------------------
#
# *) The nodes, triangles and boundary edges defined by
# the new numbering scheme are returned
#
#=========================================================================

def build_local_GA(nodes, triangles, boundaries, tri_map):

    Nnodes = len(nodes)
    Ntriangles = len(triangles)

    # Extract the nodes (using the local ID)

    GAnodes = num.take(nodes, (1, 2), 1)

    # Build a global ID to local ID mapping

    NGlobal = 0
    for i in range(Nnodes):
        if nodes[i][0] > NGlobal:
            NGlobal = nodes[i][0]

    node_map = -1*num.ones(int(NGlobal)+1, int)

    num.put(node_map, num.take(nodes, (0,), 1).astype(int),
            num.arange(Nnodes, dtype=int))

    # Change the global IDs in the triangles to the local IDs

    GAtriangles = num.zeros((Ntriangles, 3), int)
    GAtriangles[:, 0] = num.take(node_map, triangles[:, 0])
    GAtriangles[:, 1] = num.take(node_map, triangles[:, 1])
    GAtriangles[:, 2] = num.take(node_map, triangles[:, 2])

    # Change the triangle numbering in the boundaries

    GAboundaries = {}
    for b in boundaries:
        GAboundaries[tri_map[b[0]], b[1]] = boundaries[b]

    return GAnodes, GAtriangles, GAboundaries, node_map


#=========================================================================
# Change the communication format to that needed by the
# parallel advection file.
#
# *) The index contains [global triangle no,
# local triangle no.]
#
# -------------------------------------------------------
#
# *) The ghost_recv and full_send dictionaries are
# returned.
#
# *) ghost_recv dictionary is local id, global id, value
#
# *) full_recv dictionary is local id, global id, value
#
# *) The information is ordered by the global id. This
# means that the communication order is predetermined and
# local and global id do not need to be
# compared when the information is sent/received.
#
#=========================================================================

def build_local_commun(tri_map, ghostc, fullc, nproc):

    # Initialise

    full_send = {}
    ghost_recv = {}

    # Build the ghost_recv dictionary (sort the
    # information by the global numbering)

    ghostc = num.sort(ghostc, 0)

    for c in range(nproc):
        s = ghostc[:, 0]
        d = num.compress(num.equal(ghostc[:, 1], c), s)
        if len(d) > 0:
            ghost_recv[c] = [0, 0]
            ghost_recv[c][1] = d
            ghost_recv[c][0] = num.take(tri_map, d)

    # Build a temporary copy of the full_send dictionary
    # (this version allows the information to be stored
    # by the global numbering)

    tmp_send = {}
    for global_id in fullc:
        for i in range(len(fullc[global_id])):
            neigh = fullc[global_id][i]
            if neigh not in tmp_send:
                tmp_send[neigh] = []
            tmp_send[neigh].append([global_id,
                                    tri_map[global_id]])

    # Extract the full send information and put it in the form
    # required for the full_send dictionary

    for neigh in tmp_send:
        neigh_commun = num.sort(tmp_send[neigh], 0)
        full_send[neigh] = [0, 0]
        full_send[neigh][0] = neigh_commun[:, 1]
        full_send[neigh][1] = neigh_commun[:, 0]

    return ghost_recv, full_send


#=========================================================================
# Convert the format of the data to that used by ANUGA
#
#
# *) Change the nodes global ID's to an integer value,
# starting from 0. The node numbering in the triangles
# must also be updated to take this into account.
#
# *) The triangle number will also change, which affects
# the boundary tag information and the communication
# pattern.
#
# -------------------------------------------------------
#
# *) The nodes, triangles, boundary edges and communication
# pattern defined by the new numbering scheme are returned
#
#=========================================================================

def build_local_mesh(submesh, lower_t, upper_t, nproc):

    # Combine the full nodes and ghost nodes

    nodes = num.concatenate((submesh["full_nodes"],
                             submesh["ghost_nodes"]))

    ghost_layer_width = submesh["ghost_layer_width"]

    # Combine the full triangles and ghost triangles

    gtri = num.take(submesh["ghost_triangles"], (1, 2, 3), 1)
    triangles = num.concatenate((submesh["full_triangles"], gtri))

    # Combine the full boundaries and ghost boundaries

    boundaries = submesh["full_boundary"]
    for b in submesh["ghost_boundary"]:
        boundaries[b] = submesh["ghost_boundary"][b]

    # Make note of the new triangle numbers, including the ghost
    # triangles

    NGlobal = upper_t
    for i in range(len(submesh["ghost_triangles"])):
        id = submesh["ghost_triangles"][i][0]
        if id > NGlobal:
            NGlobal = id
    #index = num.zeros(int(NGlobal)+1, int)
    tri_map = -1*num.ones(int(NGlobal)+1, int)
    tri_map[lower_t:upper_t] = num.arange(upper_t-lower_t)
    for i in range(len(submesh["ghost_triangles"])):
        tri_map[submesh["ghost_triangles"][i][0]] = i+upper_t-lower_t

    # Change the node numbering (and update the numbering in the
    # triangles)

    [GAnodes, GAtriangles, GAboundary, node_map] = \
        build_local_GA(nodes, triangles, boundaries, tri_map)

    # Extract the local quantities

    quantities = {}
    for k in submesh["full_quan"]:
        Nf = len(submesh["full_quan"][k])
        Ng = len(submesh["ghost_quan"][k])
        quantities[k] = num.zeros((Nf+Ng, 1), float)
        quantities[k][0:Nf] = submesh["full_quan"][k]
        quantities[k][Nf:Nf+Ng] = submesh["ghost_quan"][k]

    # Change the communication pattern into a form needed by
    # the parallel_adv

    gcommun = submesh["ghost_commun"]
    fcommun = submesh["full_commun"]
    [ghost_rec, full_send] = \
        build_local_commun(tri_map, gcommun, fcommun, nproc)

    tri_l2g = extract_l2g_map(tri_map)
    node_l2g = extract_l2g_map(node_map)

    return GAnodes, GAtriangles, GAboundary, quantities, ghost_rec, \
        full_send, tri_map, node_map, tri_l2g, node_l2g, ghost_layer_width


#=========================================================================
#
# Handle the communication between the host machine
# (processor 0) and the processors. The host machine is
# responsible for the doing the initial grid partitioning.
#
# The routines given below should be moved to the
# build_submesh.py and build_local.py file to allow
# overlapping of  communication and computation.
# This should be done after more debugging.
#
#
#  Author: Linda Stals, June 2005
#  Modified: Linda Stals, Nov 2005 (optimise python code)
#            Steve Roberts, Aug 2009 (update to numpy)
#
#
#=========================================================================


#=========================================================================
#
# Send the submesh to processor p.
#
# *) The order and form is strongly coupled with
# rec_submesh.
#
# -------------------------------------------------------
#
# *) All of the information has been sent to processor p.
#
#=========================================================================

def send_submesh(submesh, triangles_per_proc, p, verbose=True):

    from anuga.utilities import parallel_abstraction as pypar

    myid = pypar.rank()
    nprocs = pypar.size()

    if verbose:
        print('P%d: Sending submesh to P%d' % (myid, p))

    # build and send the tagmap for the boundary conditions

    tagmap = {}
    counter = 1
    for b in submesh["full_boundary"][p]:
        bkey = submesh["full_boundary"][p][b]
        if bkey not in tagmap:
            tagmap[bkey] = counter
            counter = counter+1
    for b in submesh["ghost_boundary"][p]:
        bkey = submesh["ghost_boundary"][p][b]
        if bkey not in tagmap:
            tagmap[bkey] = counter
            counter = counter+1

    # send boundary tags
    pypar.send(tagmap, p)

    # send the quantities key information
    pypar.send(list(submesh["full_quan"].keys()), p)

    # compress full_commun
    flat_full_commun = []

    for c in submesh["full_commun"][p]:
        for i in range(len(submesh["full_commun"][p][c])):
            flat_full_commun.append([c, submesh["full_commun"][p][c][i]])

    # send the array sizes so memory can be allocated

    setup_array = num.zeros((9,), int)
    setup_array[0] = len(submesh["full_nodes"][p])
    setup_array[1] = len(submesh["ghost_nodes"][p])
    setup_array[2] = len(submesh["full_triangles"][p])
    setup_array[3] = len(submesh["ghost_triangles"][p])
    setup_array[4] = len(submesh["full_boundary"][p])
    setup_array[5] = len(submesh["ghost_boundary"][p])
    setup_array[6] = len(submesh["ghost_commun"][p])
    setup_array[7] = len(flat_full_commun)
    setup_array[8] = len(submesh["full_quan"])

    x = num.array(setup_array, int)
    pypar.send(x, p, bypass=True)

    # ghost layer width
    x = num.array(submesh["ghost_layer_width"][p], int)
    pypar.send(x, p, bypass=True)

    # send the number of triangles per processor
    x = num.array(triangles_per_proc, int)
    pypar.send(x, p, bypass=True)

    # send the nodes
    x = num.array(submesh["full_nodes"][p], float)
    pypar.send(x, p, bypass=True)

    x = num.array(submesh["ghost_nodes"][p], float)
    pypar.send(x, p, bypass=True)

    # send the triangles
    x = num.array(submesh["full_triangles"][p], int)
    pypar.send(x, p, bypass=True)

    # send ghost triangles
    x = num.array(submesh["ghost_triangles"][p], int)
    pypar.send(x, p, bypass=True)

    # send the boundary
    bc = []
    for b in submesh["full_boundary"][p]:
        bc.append([b[0], b[1], tagmap[submesh["full_boundary"][p][b]]])

    x = num.array(bc, int)
    pypar.send(x, p, bypass=True)

    bc = []
    for b in submesh["ghost_boundary"][p]:
        bc.append([b[0], b[1], tagmap[submesh["ghost_boundary"][p][b]]])

    x = num.array(bc, int)
    pypar.send(x, p, bypass=True)

    # send the communication pattern
    x = num.array(submesh["ghost_commun"][p], int)
    pypar.send(x, p, bypass=True)

    x = num.array(flat_full_commun, int)
    pypar.send(x, p, bypass=True)

    # send the quantities
    for k in submesh["full_quan"]:
        x = num.array(submesh["full_quan"][k][p], float)
        pypar.send(x, p, bypass=True)

    for k in submesh["ghost_quan"]:
        x = num.array(submesh["ghost_quan"][k][p], float)
        pypar.send(x, p, bypass=True)


#=====================================================================
#
# Receive the submesh from processor p.
#
# *) The order and form is strongly coupled with
# send_submesh.
#
# -------------------------------------------------------
#
# *) All of the information has been received by the
# processor p and passed into build_local.
#
# *) The information is returned in a form needed by the
# GA datastructure.
#
#=====================================================================

def rec_submesh_flat(p, verbose=True):

    from anuga.utilities import parallel_abstraction as pypar

    numprocs = pypar.size()
    myid = pypar.rank()

    submesh_cell = {}

    if verbose:
        print(indent+'P%d: Receiving submesh from P%d' % (myid, p))

    # receive the tagmap for the boundary conditions

    tagmap = pypar.receive(p)

    itagmap = {}
    for t in tagmap:
        itagmap[tagmap[t]] = t

    # receive the quantities key information
    qkeys = pypar.receive(p)

    # recieve information about the array sizes
    x = num.zeros((9,), int)
    pypar.receive(p, buffer=x,  bypass=True)
    setup_array = x

    no_full_nodes = setup_array[0]
    no_ghost_nodes = setup_array[1]
    no_full_triangles = setup_array[2]
    no_ghost_triangles = setup_array[3]
    no_full_boundary = setup_array[4]
    no_ghost_boundary = setup_array[5]
    no_ghost_commun = setup_array[6]
    no_full_commun = setup_array[7]
    no_quantities = setup_array[8]

    # ghost layer width
    x = num.zeros((1,), int)
    pypar.receive(p, buffer=x,  bypass=True)
    submesh_cell["ghost_layer_width"] = x[0]

    # receive the number of triangles per processor
    x = num.zeros((numprocs,), int)
    pypar.receive(p, buffer=x,  bypass=True)
    triangles_per_proc = x

    # receive the full nodes
    x = num.zeros((no_full_nodes, 3), float)
    pypar.receive(p, buffer=x,  bypass=True)
    submesh_cell["full_nodes"] = x

    # receive the ghost nodes
    x = num.zeros((no_ghost_nodes, 3), float)
    pypar.receive(p, buffer=x,  bypass=True)
    submesh_cell["ghost_nodes"] = x

    # receive the full triangles
    x = num.zeros((no_full_triangles, 3), int)
    pypar.receive(p, buffer=x,  bypass=True)
    submesh_cell["full_triangles"] = x

    # receive the ghost triangles
    x = num.zeros((no_ghost_triangles, 4), int)
    pypar.receive(p, buffer=x,  bypass=True)
    submesh_cell["ghost_triangles"] = x

    # receive the full boundary
    x = num.zeros((no_full_boundary, 3), int)
    pypar.receive(p, buffer=x,  bypass=True)
    bnd_c = x

    submesh_cell["full_boundary"] = {}
    for b in bnd_c:
        submesh_cell["full_boundary"][b[0], b[1]] = itagmap[b[2]]

    # receive the ghost boundary
    x = num.zeros((no_ghost_boundary, 3), int)
    pypar.receive(p, buffer=x,  bypass=True)
    bnd_c = x

    submesh_cell["ghost_boundary"] = {}
    for b in bnd_c:
        submesh_cell["ghost_boundary"][b[0], b[1]] = itagmap[b[2]]

    # receive the ghost communication pattern
    x = num.zeros((no_ghost_commun, 2), int)

    pypar.receive(p, buffer=x,  bypass=True)
    submesh_cell["ghost_commun"] = x

    # receive the full communication pattern
    x = num.zeros((no_full_commun, 2), int)
    pypar.receive(p, buffer=x,  bypass=True)
    full_commun = x

    submesh_cell["full_commun"] = {}
    for c in full_commun:
        submesh_cell["full_commun"][c[0]] = []
    for c in full_commun:
        submesh_cell["full_commun"][c[0]].append(c[1])

    # receive the quantities

    submesh_cell["full_quan"] = {}
    for i in range(no_quantities):
        x = num.zeros((no_full_triangles, 1), float)
        pypar.receive(p, buffer=x, bypass=True)
        submesh_cell["full_quan"][qkeys[i]] = x

    submesh_cell["ghost_quan"] = {}
    for i in range(no_quantities):
        x = num.zeros((no_ghost_triangles, 1), float)
        pypar.receive(p, buffer=x, bypass=True)
        submesh_cell["ghost_quan"][qkeys[i]] = x

    return submesh_cell, triangles_per_proc,\
        no_full_nodes, no_full_triangles


#=========================================================================
#
# Receive the submesh from processor p.
#
# *) The order and form is strongly coupled with
# send_submesh.
#
# -------------------------------------------------------
#
# *) All of the information has been received by the
# processor p and passed into build_local.
#
# *) The information is returned in a form needed by the
# GA datastructure.
#
#=====================================================================

def rec_submesh(p, verbose=True):

    from anuga.utilities import parallel_abstraction as pypar

    numproc = pypar.size()
    myid = pypar.rank()

    [submesh_cell, triangles_per_proc,
     number_of_full_nodes, number_of_full_triangles] = rec_submesh_flat(p, verbose)

    # find the full triangles assigned to this processor

    lower_t = 0
    for i in range(myid):
        lower_t = lower_t+triangles_per_proc[i]
    upper_t = lower_t+triangles_per_proc[myid]

    # convert the information into a form needed by the GA
    # datastructure

    [GAnodes, GAtriangles, boundary, quantities,
     ghost_rec, full_send,
     tri_map, node_map, tri_l2g, node_l2g,
     ghost_layer_width] = \
        build_local_mesh(submesh_cell, lower_t, upper_t, numproc)

    return GAnodes, GAtriangles, boundary, quantities,\
        ghost_rec, full_send,\
        number_of_full_nodes, number_of_full_triangles, tri_map, node_map,\
        tri_l2g, node_l2g, ghost_layer_width


#=========================================================================
#
# Extract the submesh that will belong to the
# processor 0 (i.e. processor zero)
#
#  *) See the documentation for build_submesh
#
# -------------------------------------------------------
#
#  *) A dictionary containing the full_triangles,
# full_nodes, full_boundary, ghost_triangles, ghost_nodes,
# ghost_boundary, ghost_commun and full_commun belonging
# to processor zero are returned.
#
#=========================================================================#
def extract_submesh(submesh, triangles_per_proc, p2s_map=None, p=0):

    submesh_cell = {}
    submesh_cell["ghost_layer_width"] = submesh["ghost_layer_width"][p]
    submesh_cell["full_nodes"] = submesh["full_nodes"][p]
    submesh_cell["ghost_nodes"] = submesh["ghost_nodes"][p]
    submesh_cell["full_triangles"] = submesh["full_triangles"][p]
    submesh_cell["ghost_triangles"] = submesh["ghost_triangles"][p]
    submesh_cell["full_boundary"] = submesh["full_boundary"][p]
    submesh_cell["ghost_boundary"] = submesh["ghost_boundary"][p]
    submesh_cell["ghost_commun"] = submesh["ghost_commun"][p]
    submesh_cell["full_commun"] = submesh["full_commun"][p]
    submesh_cell["full_quan"] = {}
    submesh_cell["ghost_quan"] = {}
    for k in submesh["full_quan"]:
        submesh_cell["full_quan"][k] = submesh["full_quan"][k][p]
        submesh_cell["ghost_quan"][k] = submesh["ghost_quan"][k][p]

    # FIXME SR: I think there is already a structure with this info in the mesh
    lower_t = 0
    for i in range(p):
        lower_t = lower_t+triangles_per_proc[i]
    upper_t = lower_t+triangles_per_proc[p]

    numprocs = len(triangles_per_proc)
    points, vertices, boundary, quantities, ghost_recv_dict, \
        full_send_dict, tri_map, node_map, tri_l2g, node_l2g, \
        ghost_layer_width = \
        build_local_mesh(submesh_cell, lower_t, upper_t, numprocs)

    if p2s_map is None:
        pass
    else:
        try:
            tri_l2g = p2s_map[tri_l2g]
        except (IndexError, KeyError):
            tri_l2g = p2s_map

    return points, vertices, boundary, quantities, ghost_recv_dict, \
        full_send_dict, tri_map, node_map, tri_l2g, node_l2g, ghost_layer_width


def extract_l2g_map(map):
    # Extract l2g data  from corresponding map
    # Maps

    import numpy as num

    b = num.arange(len(map))

    l_ids = num.extract(map > -1, map)
    g_ids = num.extract(map > -1, b)


#    print len(g_ids)
#    print len(l_ids)
#    print l_ids
#    print g_ids

    l2g = num.zeros_like(g_ids)
    l2g[l_ids] = g_ids

    return l2g

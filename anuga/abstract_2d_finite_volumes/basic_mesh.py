"""Basic_mesh -- minimal mesh for partitioning and parallel distribution.

Stores only nodes, triangles, boundary tags, and geo_reference.
All geometric arrays needed by the finite-volume solver (normals,
edge lengths, areas, radii, vertex_coordinates) are NOT computed here --
they are deferred to when each rank constructs its Parallel_domain from
the distributed submesh.

The neighbours array and centroid_coordinates are computed lazily on
first access, since they are needed by the partitioning algorithms.
"""

import numpy as num
from anuga.coordinate_transforms.geo_reference import Geo_reference


class Basic_mesh:
    """Minimal mesh for partitioning and parallel distribution.

    Unlike the full Mesh class (neighbour_mesh.Mesh), Basic_mesh skips
    computation of normals, edge lengths, areas, radii, vertex_coordinates,
    and edge_midpoint_coordinates at construction time.  Those arrays are
    only needed by the finite-volume solver and will be computed later when
    each MPI rank constructs its Parallel_domain from the distributed submesh.

    neighbours and centroid_coordinates are computed lazily on first access
    (needed by the partitioning algorithms).

    If the triangulation was produced by the triangle library (via
    create_mesh_from_regions / pmesh_to_basic_mesh), the pre-computed
    triangle_neighbours array can be passed directly and the neighbour
    structure is not recomputed from scratch.

    Parameters
    ----------
    nodes : array_like, shape (M, 2)
        Node x, y coordinates (relative to geo_reference).
    triangles : array_like, shape (N, 3)
        Triangle vertex indices into nodes.
    boundary : dict, optional
        Boundary tag dict mapping (triangle_id, edge_index) -> tag string.
        If None an empty dict is used.
    geo_reference : Geo_reference, optional
        Coordinate origin and zone.  Defaults to Geo_reference().
    triangle_neighbours : array_like, shape (N, 3), optional
        Pre-computed neighbour array (-1 for boundary edges), as produced
        by the triangle library.  If supplied the neighbour structure is
        assigned directly rather than recomputed.
    """

    def __init__(self, nodes, triangles, boundary=None,
                 geo_reference=None, triangle_neighbours=None):

        self.nodes = num.array(nodes, float)
        self.triangles = num.array(triangles, int)
        self.number_of_triangles = len(self.triangles)
        self.number_of_nodes = len(self.nodes)

        if geo_reference is None:
            self.geo_reference = Geo_reference()
        else:
            self.geo_reference = geo_reference

        self.boundary = {} if boundary is None else boundary

        # Store pre-computed neighbours if supplied by the triangulator.
        self._triangle_neighbours = (
            num.array(triangle_neighbours, int)
            if triangle_neighbours is not None else None)

        # Lazily computed.
        self._neighbours = None
        self._centroid_coordinates = None

    # ------------------------------------------------------------------
    # Container protocol
    # ------------------------------------------------------------------

    def __len__(self):
        return self.number_of_triangles

    def __repr__(self):
        return ('Basic_mesh: %d triangles, %d nodes'
                % (self.number_of_triangles, self.number_of_nodes))

    # ------------------------------------------------------------------
    # Lazy properties
    # ------------------------------------------------------------------

    @property
    def neighbours(self):
        """(N, 3) int array of neighbour triangle indices (-1 = boundary)."""
        if self._neighbours is None:
            self._build_neighbours()
        return self._neighbours

    def _build_neighbours(self):
        N = self.number_of_triangles
        if self._triangle_neighbours is not None:
            self._neighbours = self._triangle_neighbours.copy()
        else:
            from anuga.abstract_2d_finite_volumes import neighbour_table_ext
            nbrs = -num.ones((N, 3), int)
            nbr_edges = -num.ones((N, 3), int)
            n_boundaries = num.zeros(N, int)
            neighbour_table_ext.build_neighbour_structure(
                self.number_of_nodes,
                self.triangles,
                nbrs,
                nbr_edges,
                n_boundaries,
            )
            self._neighbours = nbrs

    @property
    def centroid_coordinates(self):
        """(N, 2) float array of triangle centroid x, y coordinates."""
        if self._centroid_coordinates is None:
            t = self.triangles
            n = self.nodes
            self._centroid_coordinates = (
                (n[t[:, 0]] + n[t[:, 1]] + n[t[:, 2]]) / 3.0)
        return self._centroid_coordinates

    # ------------------------------------------------------------------
    # Methods expected by partition_mesh and the partitioning algorithms
    # ------------------------------------------------------------------

    def get_number_of_nodes(self):
        return self.number_of_nodes

    def reorder(self, new_order, in_place=True):
        """Reorder triangles by a permutation array.

        Called by partition_mesh after computing the partition ordering.

        Parameters
        ----------
        new_order : array_like of int, length N
            Permutation of 0..N-1 giving the new triangle ordering.
        in_place : bool
            If True (default) reorder this mesh in-place and return self.
            If False return a new Basic_mesh sharing the nodes array.

        Returns
        -------
        Basic_mesh
        """
        # Force neighbour computation before any reordering.
        # _triangle_neighbours (if set) contains pre-reorder indices; once
        # the triangle numbering changes those indices become stale.  Building
        # _neighbours now ensures they are remapped correctly below and
        # prevents _build_neighbours() from later reconstructing from the
        # stale _triangle_neighbours cache.
        _ = self.neighbours  # triggers _build_neighbours() if not yet done

        new_order = num.array(new_order, int)
        N = self.number_of_triangles
        inv_order = num.empty_like(new_order)
        inv_order[new_order] = num.arange(N, dtype=int)

        if in_place:
            target = self
        else:
            import copy
            target = copy.copy(self)
            # nodes are not reordered -- share the array
            target.nodes = self.nodes
            target.triangles = self.triangles.copy()
            target.boundary = dict(self.boundary)
            target._neighbours = (
                self._neighbours.copy()
                if self._neighbours is not None else None)
            target._centroid_coordinates = (
                self._centroid_coordinates.copy()
                if self._centroid_coordinates is not None else None)

        target.triangles = target.triangles[new_order]
        target.boundary = {(int(inv_order[i]), j): v
                           for (i, j), v in target.boundary.items()}

        if target._neighbours is not None:
            # Remap neighbour indices then reorder rows.
            flat = target._neighbours.ravel()
            mask = flat >= 0
            flat[mask] = inv_order[flat[mask]]
            target._neighbours = target._neighbours[new_order]

        if target._centroid_coordinates is not None:
            target._centroid_coordinates = (
                target._centroid_coordinates[new_order])

        return target


# ---------------------------------------------------------------------------
# Factory functions
# ---------------------------------------------------------------------------

def rectangular_basic_mesh(m, n, len1=1.0, len2=1.0, origin=(0.0, 0.0)):
    """Create a Basic_mesh from a rectangular grid (2*m*n triangles).

    Parameters
    ----------
    m, n : int
        Number of cells in x and y directions.
    len1, len2 : float
        Physical dimensions in x and y.
    origin : (float, float)
        Bottom-left corner coordinates.

    Returns
    -------
    Basic_mesh with 2*m*n triangles and (m+1)*(n+1) nodes.
    """
    from anuga.abstract_2d_finite_volumes.mesh_factory import \
        rectangular_with_neighbours
    points, elements, boundary, neighbours, neighbour_edges = \
        rectangular_with_neighbours(m, n, len1, len2, origin)
    return Basic_mesh(points, elements, boundary=boundary,
                      triangle_neighbours=neighbours)


def rectangular_cross_basic_mesh(m, n, len1=1.0, len2=1.0, origin=(0.0, 0.0)):
    """Create a Basic_mesh from a rectangular-cross grid (4*m*n triangles).

    Parameters
    ----------
    m, n : int
        Number of cells in x and y directions.
    len1, len2 : float
        Physical dimensions in x and y.
    origin : (float, float)
        Bottom-left corner coordinates.

    Returns
    -------
    Basic_mesh with 4*m*n triangles and (m+1)*(n+1) + m*n nodes.
    """
    from anuga.abstract_2d_finite_volumes.mesh_factory import \
        rectangular_cross_with_neighbours
    points, elements, boundary, neighbours, neighbour_edges = \
        rectangular_cross_with_neighbours(m, n, len1, len2, origin)
    return Basic_mesh(points, elements, boundary=boundary,
                      triangle_neighbours=neighbours)


def basic_mesh_from_mesh_file(mesh_filename, verbose=False):
    """Read a .tsh or .msh mesh file and return a Basic_mesh.

    Vertex attributes and triangle region tags stored in the file are attached
    to the returned object as extra instance attributes so the caller does not
    need to re-open the file.

    Parameters
    ----------
    mesh_filename : str or Path
        Path to the mesh file (.tsh or .msh format).
    verbose : bool, optional
        If True, print progress and summary messages. Default False.

    Returns
    -------
    Basic_mesh
        The mesh object, with the following additional instance attributes:

        ``vertex_attributes`` : numpy.ndarray, shape (M, n) or None
            Per-vertex attribute values in the same row order as
            ``bm.nodes``.  None when the file contains no attributes.
        ``vertex_attribute_titles`` : list of str
            Column names for ``vertex_attributes``; empty list when None.
        ``triangle_tags`` : list of str
            Region tag string for each triangle (length N); empty list when
            the file contains no region tags.

    Examples
    --------
    >>> bm = basic_mesh_from_mesh_file('my_mesh.tsh', verbose=True)
    >>> elev_idx = bm.vertex_attribute_titles.index('elevation')
    >>> elevation = bm.vertex_attributes[:, elev_idx]
    """
    from anuga.load_mesh.loadASCII import import_mesh_file
    from anuga.abstract_2d_finite_volumes.pmesh2domain import pmesh_dict_to_tag_dict

    if verbose:
        print(f'Reading mesh file: {mesh_filename}')

    mesh_dict = import_mesh_file(str(mesh_filename))

    nodes     = mesh_dict['vertices']
    triangles = mesh_dict['triangles']
    geo_ref   = mesh_dict['geo_reference']
    boundary  = pmesh_dict_to_tag_dict(mesh_dict)

    bm = Basic_mesh(nodes, triangles,
                    boundary=boundary,
                    geo_reference=geo_ref)

    # --- vertex attributes ---
    raw_va = mesh_dict.get('vertex_attributes')
    titles = mesh_dict.get('vertex_attribute_titles') or []
    if raw_va is not None and len(raw_va) > 0:
        bm.vertex_attributes = num.array(raw_va, dtype=float)
        bm.vertex_attribute_titles = list(titles)
        if verbose:
            print(f'  vertex attributes ({bm.vertex_attributes.shape[1]}): '
                  f'{bm.vertex_attribute_titles}')
    else:
        bm.vertex_attributes = None
        bm.vertex_attribute_titles = []
        if verbose:
            print('  no vertex attributes')

    # --- triangle tags (region labels) ---
    tri_tags = mesh_dict.get('triangle_tags') or []
    bm.triangle_tags = list(tri_tags)
    if verbose:
        unique_tags = sorted(set(t for t in tri_tags if t))
        if unique_tags:
            print(f'  triangle tags: {unique_tags}')
        else:
            print('  no triangle tags')

    if verbose:
        print(f'  {bm}')

    return bm


from anuga.shallow_water.shallow_water_domain import Domain
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_with_neighbours, rectangular_cross_with_neighbours
from anuga.abstract_2d_finite_volumes.neighbour_mesh import Mesh

from anuga.abstract_2d_finite_volumes.pmesh2domain import pmesh_to_domain_instance, pmesh_to_mesh, pmesh_to_basic_mesh



#-----------------------------
# rectangular domains and meshes
#-----------------------------
def rectangular_cross_domain(*args, **kwargs):
    """Create a rectangular domain.

    The triangular mesh is made up of m by n uniform rectangular cells divided
    into 4 triangles in a cross pattern

    Parameters
    ----------
    m : int
        Number of cells in x direction
    n : int
        Number of cells in y direction
    len1 : float, optional
        Length of domain in x direction (left to right) (default 1.0)
    len2 : float, optional
        Length of domain in y direction (bottom to top) (default 1.0)
    origin : tuple, optional
        Tuple (x, y) specifying location of lower left corner of domain
        (default (0, 0))
    verbose : bool, optional
        Boolean flag to output information (default False)

    Returns
    -------
    Domain
        Shallow water domain instance
    """

    try:
        verbose = kwargs.pop('verbose')
    except KeyError:
        verbose = False


    points, vertices, boundary, neighbours, neighbour_edges = \
        rectangular_cross_with_neighbours(*args, **kwargs)

    mesh = Mesh(points, vertices, boundary,
                triangle_neighbours=neighbours,
                triangle_neighbour_edges=neighbour_edges)

    return Domain(mesh, verbose=verbose)

def rectangular_cross_mesh(*args, **kwargs):
    """Create a rectangular mesh.

    The triangular mesh is made up of m by n uniform rectangular cells divided
    into 4 triangles in a cross pattern

    Parameters
    ----------
    m : int
        Number of cells in x direction
    n : int
        Number of cells in y direction
    len1 : float, optional
        Length of domain in x direction (left to right) (default 1.0)
    len2 : float, optional
        Length of domain in y direction (bottom to top) (default 1.0)
    origin : tuple, optional
        Tuple (x, y) specifying location of lower left corner of domain
        (default (0, 0))
    verbose : bool, optional
        Boolean flag to output information (default False)

    Returns
    -------
    mesh
        Mesh instance
    """

    try:
        verbose = kwargs.pop('verbose')
    except KeyError:
        verbose = False


    points, vertices, boundary, neighbours, neighbour_edges = \
        rectangular_cross_with_neighbours(*args, **kwargs)

    mesh = Mesh(points, vertices, boundary,
                triangle_neighbours=neighbours,
                triangle_neighbour_edges=neighbour_edges)

    return mesh

#---------------------------
# Create basic mesh from regions
#---------------------------

def create_basic_mesh_from_regions(bounding_polygon,
                                   boundary_tags,
                                   maximum_triangle_area=None,
                                   interior_regions=None,
                                   interior_holes=None,
                                   hole_tags=None,
                                   poly_geo_reference=None,
                                   mesh_geo_reference=None,
                                   breaklines=None,
                                   regionPtArea=None,
                                   minimum_triangle_angle=28.0,
                                   fail_if_polygons_outside=True,
                                   use_cache=False,
                                   verbose=False,
                                   size_field=None,
                                   relax_passes=0):
    """Create a Basic_mesh from bounding polygons and resolutions.

    Like :func:`create_domain_from_regions` but returns a
    :class:`~anuga.abstract_2d_finite_volumes.basic_mesh.Basic_mesh` instead
    of a full :class:`~anuga.shallow_water.shallow_water_domain.Domain`.
    Use when you need mesh topology (vertices, triangles, boundary tags,
    geo-reference) but do not yet need normals, edge lengths, areas, or the
    full simulation machinery — for example when constructing a parallel
    partition before building the per-processor domains.

    Parameters
    ----------
    bounding_polygon : list of tuples
        Points in Eastings and Northings, relative to the zone stated in
        poly_geo_reference if specified. Otherwise plain x, y coordinates.
    boundary_tags : dict
        Symbolic tags where each key maps to a list of indices referring to
        segments associated with that tag. Omitted segments are assigned
        the default tag ''.
    maximum_triangle_area : float, optional
        Maximal area per triangle for the bounding polygon, excluding the
        interior regions.
    interior_regions : list of tuples, optional
        List of (polygon, resolution) tuples for each region to be separately
        refined. Polygon lines must not cross or overlap.
    interior_holes : list of polygons, optional
        Polygons for each hole.
    hole_tags : list of dict, optional
        List of tag segment dictionaries for hole boundaries.
    poly_geo_reference : GeoReference, optional
        Geo-reference of the bounding polygon and interior polygons. If None,
        assume absolute coordinates.
    mesh_geo_reference : GeoReference, optional
        Geo-reference of the mesh to be created. If None, one will be
        automatically generated using the lower left corner of
        bounding_polygon.
    breaklines : list of polygons, optional
        Lines to be preserved by the triangulation algorithm.
    regionPtArea : list of 3-tuples, optional
        Points with maximum area for the region containing each point.
    minimum_triangle_angle : float, optional
        Minimum triangle angle in degrees (default: 28.0).
    fail_if_polygons_outside : bool, optional
        If True (default), raise an Exception when interior polygons fall
        outside the bounding polygon. If False, ignore and continue.
    use_cache : bool, optional
        Whether to use caching (default: False).
    verbose : bool, optional
        Output information (default: False).
    size_field : Size_field or callable, optional
        Target triangle area as a function of position, applied as well as
        the maximum areas above: a triangle is refined if its area exceeds
        the field at its centroid. Either an :class:`anuga.Size_field`
        (see :func:`anuga.size_field_from_regions`) or a callable
        ``f(x, y) -> maximum area`` taking absolute coordinates. Unlike
        interior regions it adds no segments to the mesh. Needs meshpy.
    relax_passes : int, optional
        If positive, follow the triangulation with up to this many passes of
        worst-cell relaxation, which moves the vertices of the smallest
        cells to enlarge them and so lengthen the timestep (default: 0).
        See :meth:`anuga.pmesh.mesh.Pmesh.relax_worst_cells`.

    Returns
    -------
    Basic_mesh
        A Basic_mesh instance containing vertices, triangles, boundary tags,
        and geo-reference, but without normals, edge lengths, or areas.

    See Also
    --------
    create_pmesh_from_regions : returns the raw Pmesh object.
    create_domain_from_regions : returns a full shallow-water Domain.
    """

    args = (bounding_polygon, boundary_tags)
    kwargs = {'maximum_triangle_area': maximum_triangle_area,
              'interior_regions': interior_regions,
              'interior_holes': interior_holes,
              'hole_tags': hole_tags,
              'poly_geo_reference': poly_geo_reference,
              'mesh_geo_reference': mesh_geo_reference,
              'breaklines': breaklines,
              'regionPtArea': regionPtArea,
              'minimum_triangle_angle': minimum_triangle_angle,
              'fail_if_polygons_outside': fail_if_polygons_outside,
              'verbose': verbose,
              'size_field': size_field,
              'relax_passes': relax_passes}

    if use_cache is True:
        try:
            from anuga.caching import cache
        except ImportError:
            msg = 'Caching was requested, but caching module ' \
                  'could not be imported'
            raise Exception(msg)

        basic_mesh = cache(_create_basic_mesh_from_regions,
                           args, kwargs,
                           verbose=verbose,
                           compression=False)
    else:
        basic_mesh = _create_basic_mesh_from_regions(*args, **kwargs)

    return basic_mesh


def _create_basic_mesh_from_regions(bounding_polygon,
                                    boundary_tags,
                                    maximum_triangle_area=None,
                                    interior_regions=None,
                                    interior_holes=None,
                                    hole_tags=None,
                                    poly_geo_reference=None,
                                    mesh_geo_reference=None,
                                    breaklines=None,
                                    regionPtArea=None,
                                    minimum_triangle_angle=28.0,
                                    fail_if_polygons_outside=True,
                                    verbose=True,
                                    size_field=None,
                                    relax_passes=0):
    """Internal implementation — see create_basic_mesh_from_regions."""

    from anuga.pmesh.mesh_interface import create_pmesh_from_regions

    pmesh = create_pmesh_from_regions(bounding_polygon,
                                     boundary_tags,
                                     maximum_triangle_area=maximum_triangle_area,
                                     interior_regions=interior_regions,
                                     interior_holes=interior_holes,
                                     hole_tags=hole_tags,
                                     poly_geo_reference=poly_geo_reference,
                                     mesh_geo_reference=mesh_geo_reference,
                                     breaklines=breaklines,
                                     regionPtArea=regionPtArea,
                                     minimum_triangle_angle=minimum_triangle_angle,
                                     fail_if_polygons_outside=fail_if_polygons_outside,
                                     use_cache=False,
                                     verbose=verbose,
                                     size_field=size_field,
                                     relax_passes=relax_passes)

    return pmesh_to_basic_mesh(pmesh, verbose=verbose)


#----------------------------
# Create domain from file
#----------------------------
def create_domain_from_file(filename, DomainClass=Domain):
    """Create a domain from a mesh file.

    Supported formats include .msh (Gmsh format) and other mesh file formats
    compatible with pmesh.

    Parameters
    ----------
    filename : str
        Path to the mesh file to load and convert to a domain instance.
    DomainClass : type, optional
        The domain class to instantiate. Default is Domain.

    Returns
    -------
    domain : DomainClass
        An instance of the specified DomainClass created from the mesh file.

    Notes
    -----
    This function is a wrapper around pmesh_to_domain_instance that converts
    a mesh file into a domain object suitable for simulation.

    Examples
    --------
    >>> domain = create_domain_from_file('mesh.msh')
    >>> domain = create_domain_from_file('mesh.msh', DomainClass=CustomDomain)
    """

    return pmesh_to_domain_instance(filename,DomainClass=DomainClass)

#---------------------------
# Create domain from regions
#---------------------------

def create_domain_from_regions(bounding_polygon,
                               boundary_tags,
                               maximum_triangle_area=None,
                               interior_regions=None,
                               interior_holes=None,
                               hole_tags=None,
                               poly_geo_reference=None,
                               mesh_geo_reference=None,
                               breaklines=None,
                               regionPtArea=None,
                               minimum_triangle_angle=28.0,
                               fail_if_polygons_outside=True,
                               use_cache=False,
                               verbose=False,
                               size_field=None,
                               relax_passes=0):

    """Create domain from bounding polygons and resolutions.

    Parameters
    ----------
    bounding_polygon : list of tuples
        Points in Eastings and Northings, relative to the zone stated in
        poly_geo_reference if specified. Otherwise points are just x, y
        coordinates with no particular association to any location.
    boundary_tags : dict
        Symbolic tags where each key maps to a list of indices referring to
        segments associated with that tag. Omitted segments are assigned
        the default tag ''.
    maximum_triangle_area : float, optional
        Maximal area per triangle for the bounding polygon, excluding the
        interior regions.
    interior_regions : list of tuples, optional
        List of (polygon, resolution) tuples for each region to be separately
        refined. Polygon lines must not cross or overlap, and polygons should
        not be close to each other. Interior regions outside the bounding_polygon
        will raise an error.
    interior_holes : list of polygons, optional
        Polygons for each hole. These polygons do not need to be closed, but
        their points must be specified in counter-clockwise order.
    hole_tags : list of dict, optional
        List of tag segment dictionaries. Segments cannot share points.
    poly_geo_reference : GeoReference, optional
        Geo-reference of the bounding polygon and interior polygons. If None,
        assume absolute coordinates.
    mesh_geo_reference : GeoReference, optional
        Geo-reference of the mesh to be created. If None, one will be
        automatically generated using the lower left corner of bounding_polygon
        as the x and y values.
    breaklines : list of polygons, optional
        Lines to be preserved by the triangulation algorithm (e.g., coastlines,
        walls). Polygons are not closed.
    regionPtArea : list of 3-tuples, optional
        Points with maximum area for the region containing each point.
    minimum_triangle_angle : float, optional
        Minimum triangle angle in degrees (default: 28.0).
    fail_if_polygons_outside : bool, optional
        If True (default), raise an Exception when interior polygons fall
        outside the bounding polygon. If False, ignore and continue.
    use_cache : bool, optional
        Whether to use caching (default: False).
    verbose : bool, optional
        Output information (default: False).
    size_field : Size_field or callable, optional
        Target triangle area as a function of position, applied as well as
        the maximum areas above: a triangle is refined if its area exceeds
        the field at its centroid. Either an :class:`anuga.Size_field`
        (see :func:`anuga.size_field_from_regions`) or a callable
        ``f(x, y) -> maximum area`` taking absolute coordinates. Unlike
        interior regions it adds no segments to the mesh. Needs meshpy.
    relax_passes : int, optional
        If positive, follow the triangulation with up to this many passes of
        worst-cell relaxation, which moves the vertices of the smallest
        cells to enlarge them and so lengthen the timestep (default: 0).
        See :meth:`anuga.pmesh.mesh.Pmesh.relax_worst_cells`.

    Returns
    -------
    Domain
        A shallow water domain instance.

    Notes
    -----
    Interior regions should be fully nested, as overlaps may cause unintended
    resolutions.
    """


    # Build arguments and keyword arguments for use with caching or apply.
    args = (bounding_polygon,
            boundary_tags)

    # if mesh_filename is None:
    #     import tempfile
    #     import time
    #     mesh_filename = 'mesh_%d.msh'%int(time.time())

    kwargs = {'maximum_triangle_area': maximum_triangle_area,
              'interior_regions': interior_regions,
              'interior_holes': interior_holes,
              'hole_tags': hole_tags,
              'poly_geo_reference': poly_geo_reference,
              'mesh_geo_reference': mesh_geo_reference,
              'breaklines' : breaklines,
              'regionPtArea' : regionPtArea,
              'minimum_triangle_angle': minimum_triangle_angle,
              'fail_if_polygons_outside': fail_if_polygons_outside,
              'verbose': verbose,
              'size_field': size_field,
              'relax_passes': relax_passes} #FIXME (Ole): See ticket:14

    # Call underlying engine with or without caching

    # FIXME SR: The _create_domain_from_regions function creates a mesh file which is then read in to
    # create the domain. It would make sense to take the created pmesh Mesh and then convert it to an anuga.Mesh
    # and use that to create the domain.
    if use_cache is True:
        try:
            from anuga.caching import cache
        except ImportError:
            msg = 'Caching was requested, but caching module'+\
                  'could not be imported'
            raise (msg)


        domain = cache(_create_domain_from_regions,
                       args, kwargs,
                       verbose=verbose,
                       compression=False)
    else:
        domain = _create_domain_from_regions(*args, **kwargs)

    return domain


def _create_domain_from_regions(bounding_polygon,
                                boundary_tags,
                                maximum_triangle_area=None,
                                interior_regions=None,
                                interior_holes=None,
                                hole_tags=None,
                                poly_geo_reference=None,
                                mesh_geo_reference=None,
                                breaklines=None,
                                regionPtArea=None,
                                minimum_triangle_angle=28.0,
                                fail_if_polygons_outside=True,
                                verbose=True,
                                size_field=None,
                                relax_passes=0):
    """_create_domain_from_regions - internal function.

    See create_domain_from_regions for documentation.
    """

    #from anuga.shallow_water.shallow_water_domain import Domain
    from anuga.pmesh.mesh_interface import create_pmesh_from_regions

    pmesh = create_pmesh_from_regions(bounding_polygon,
                             boundary_tags,
                             maximum_triangle_area=maximum_triangle_area,
                             interior_regions=interior_regions,
                             interior_holes=interior_holes,
                             hole_tags=hole_tags,
                             poly_geo_reference=poly_geo_reference,
                             mesh_geo_reference=mesh_geo_reference,
                             breaklines=breaklines,
                             regionPtArea=regionPtArea,
                             minimum_triangle_angle=minimum_triangle_angle,
                             fail_if_polygons_outside=fail_if_polygons_outside,
                             use_cache=False,
                             verbose=verbose,
                             size_field=size_field,
                             relax_passes=relax_passes)

    mesh = pmesh_to_mesh(pmesh)

    domain = Domain(mesh, use_cache=False, verbose=verbose)


    return domain


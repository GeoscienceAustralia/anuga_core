
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
        refined. Polygons may be nested or separate. Outlines that cross
        raise PolygonError. For overlapping refinement see
        size_field_from_regions.
    interior_holes : list of polygons, optional
        Polygons for each hole. Outlines must not cross.
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
        The vertices you give stay where they are; a point the generator
        added part-way along a straight segment (a hole's side, a
        breakline) may slide along it, so the segment keeps its shape.
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
# Create domain from a DEM
#----------------------------
def create_domain_from_dem(dem,
                           maximum_triangle_area,
                           bounding_polygon=None,
                           boundary_tags=None,
                           streams=None,
                           floodplain=None,
                           regions=None,
                           interior_holes=None,
                           relax_passes=30,
                           pits=None,
                           pit_outlets=None,
                           pit_max_depth=None,
                           pit_where=None,
                           verbose=False,
                           mesh_filename=None,
                           **kwargs):
    """Create a domain meshed and elevated from a DEM.

    The mesh is sized from the DEM's flow paths, as in
    :func:`anuga.create_pmesh_from_dem`: fine along the streams (or streets
    and overland flow paths) and over the floors of their valleys, coarse
    elsewhere. Each triangle's elevation is then the mean of the DEM over
    it, set at the centroids.

    Parameters
    ----------
    dem : str or (x, y, Z)
        A DEM file (.asc, .grd, .dem or .tif) or a raster tuple
        ``(x, y, Z)``; see :mod:`anuga.pmesh.dem_mesh`.
    maximum_triangle_area : float
        Target area away from all streams, floodplains and regions.
    bounding_polygon : list of points, optional
        Outline of the domain. Default: the DEM's extent, tagged
        ``bottom``, ``right``, ``top`` and ``left``.
    boundary_tags : dict, optional
        As for :func:`create_domain_from_regions`.
    streams : list of (drainage_area, area), optional
        Cells that at least ``drainage_area`` drains through (in the DEM's
        units squared) get triangles of at most ``area``.
    floodplain : list of (height, drainage_area, area), optional
        Cells within ``height`` above the stream they drain to, where it
        drains at least ``drainage_area``, get triangles of at most
        ``area``.
    regions : list of (polygon, area), optional
        Further polygons to refine, without adding their segments.
    interior_holes : list of polygons, optional
        Holes in the mesh (buildings, say); water is routed around them.
    relax_passes : int, optional
        Passes of worst-cell relaxation after meshing. Default 30.
    pits : None, 'fill' or 'breach', optional
        Remove the pits the mesh makes of the DEM with
        :func:`anuga.remove_mesh_pits`. Default None: keep the elevation as
        sampled.
    pit_outlets : list of str, optional
        Boundary tags where water leaves; needed with ``pits``.
    pit_max_depth : float, optional
        Deepest cut a breach may make before it fills instead.
    pit_where : boolean array or callable, optional
        Cells whose elevation pit removal may change, or a function of the
        domain returning them. Default all.
    verbose : bool, optional
    mesh_filename : str, optional
        If given, also write the mesh the domain is built on to this file
        (``.msh`` or ``.tsh``). The file holds the mesh only, not the
        elevation.
    **kwargs :
        Passed on to :func:`anuga.create_pmesh_from_dem` (``grade``,
        ``band``, ``size_field``, ``breaklines``, ``hole_tags``,
        ``minimum_triangle_angle``, ``poly_geo_reference``, ...).

    Returns
    -------
    Domain
        With its elevation set. ``domain.pit_stats`` holds the result of
        pit removal, if any.

    Examples
    --------
    >>> domain = anuga.create_domain_from_dem(
    ...     'topography.asc', maximum_triangle_area=20.0,
    ...     streams=[(5000.0, 1.0), (500.0, 4.0)], interior_holes=houses)
    """
    from anuga.pmesh.dem_mesh import (create_pmesh_from_dem,
                                      elevation_from_dem, remove_mesh_pits)

    if pits is not None and not pit_outlets:
        raise ValueError('pits needs pit_outlets: the boundary tags where '
                         'water leaves')

    pmesh = create_pmesh_from_dem(dem, maximum_triangle_area,
                                  bounding_polygon=bounding_polygon,
                                  boundary_tags=boundary_tags,
                                  streams=streams, floodplain=floodplain,
                                  regions=regions,
                                  interior_holes=interior_holes,
                                  relax_passes=relax_passes,
                                  verbose=verbose, **kwargs)
    mesh = pmesh_to_mesh(pmesh)
    if mesh_filename is not None:
        pmesh.export_mesh_file(mesh_filename)
    domain = Domain(mesh, use_cache=False, verbose=verbose)
    domain.set_quantity('elevation', elevation_from_dem(domain, dem),
                        location='centroids')

    if pits is not None:
        where = pit_where(domain) if callable(pit_where) else pit_where
        domain.pit_stats = remove_mesh_pits(domain, pit_outlets, method=pits,
                                            max_depth=pit_max_depth,
                                            where=where, verbose=verbose)
    return domain


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
                               relax_passes=0,
                               mesh_filename=None):

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
        The vertices you give stay where they are; a point the generator
        added part-way along a straight segment (a hole's side, a
        breakline) may slide along it, so the segment keeps its shape.
        See :meth:`anuga.pmesh.mesh.Pmesh.relax_worst_cells`.
    mesh_filename : str, optional
        If given, also write the mesh the domain is built on to this file
        (``.msh`` or ``.tsh``), to look at or to reuse with
        :func:`anuga.create_domain_from_file`. The domain is built in
        memory either way. With ``use_cache=True`` the file is written only
        when the mesh is generated, not when the domain comes from the cache.

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
              'relax_passes': relax_passes,
              'mesh_filename': mesh_filename} #FIXME (Ole): See ticket:14

    # Call underlying engine with or without caching
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
                                relax_passes=0,
                                mesh_filename=None):
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

    # pmesh_to_mesh has triangulated the pmesh, so the file holds the mesh
    # the domain is built on.
    if mesh_filename is not None:
        pmesh.export_mesh_file(mesh_filename)

    domain = Domain(mesh, use_cache=False, verbose=verbose)


    return domain


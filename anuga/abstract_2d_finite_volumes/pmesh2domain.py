"""Class pmesh2domain - Converting .tsh files to domains


   Copyright 2004
   Ole Nielsen, Stephen Roberts, Duncan Gray, Christopher Zoppou
   Geoscience Australia
"""
import sys
import numpy as num
import anuga.utilities.log as log



def pmesh_to_basic_mesh(pmesh_instance, verbose=False):
    """Convert a Pmesh instance to a Basic_mesh instance.

    Like pmesh_to_mesh but returns a Basic_mesh, which omits computation of
    normals, edge lengths, areas, radii, and vertex_coordinates.  Use when
    only mesh topology is needed (e.g. for parallel distribution without a
    pre-built Domain).

    The triangle_neighbours array produced by the triangle library is passed
    through directly, so the neighbour structure is not recomputed.

    Parameters
    ----------
    pmesh_instance : anuga.pmesh.mesh.Pmesh
    verbose : bool, optional

    Returns
    -------
    Basic_mesh
    """
    from anuga.abstract_2d_finite_volumes.basic_mesh import Basic_mesh

    if pmesh_instance.tri_mesh is None:
        pmesh_instance.generate_mesh(verbose=verbose)

    mesh_dict = pmesh_instance.mesh2io_dict()

    vertex_coordinates = mesh_dict['vertices']
    triangles          = mesh_dict['triangles']
    geo_reference      = mesh_dict['geo_reference']
    triangle_neighbors = mesh_dict.get('triangle_neighbors')
    boundary           = pmesh_dict_to_tag_dict(mesh_dict)

    return Basic_mesh(vertex_coordinates, triangles,
                     boundary=boundary,
                     geo_reference=geo_reference,
                     triangle_neighbours=triangle_neighbors)


def pmesh_to_mesh(pmesh_instance, verbose=False):
    """Convert a Pmesh instance to a neighbour_mesh.Mesh instance.

    Parameters
    ----------
    pmesh_instance : anuga.pmesh.mesh.Pmesh
        A Pmesh instance with a completed triangulation (i.e. generate_mesh
        has been called).
    verbose : bool, optional
        If True, print verbose output. Default is False.

    Returns
    -------
    anuga.abstract_2d_finite_volumes.neighbour_mesh.Mesh
        A Mesh instance built from the Pmesh triangulation.
    """
    from anuga.abstract_2d_finite_volumes.neighbour_mesh import Mesh

    # create_mesh_from_regions only calls generate_mesh when a filename is
    # supplied; if the caller passed the raw pmesh object we may need to
    # triangulate it first.
    if pmesh_instance.tri_mesh is None:
        pmesh_instance.generate_mesh(verbose=verbose)

    mesh_dict = pmesh_instance.mesh2io_dict()

    vertex_coordinates = mesh_dict['vertices']
    triangles = mesh_dict['triangles']
    geo_reference = mesh_dict['geo_reference']
    triangle_neighbors = mesh_dict.get('triangle_neighbors')
    tagged_elements = build_tagged_elements_dictionary(mesh_dict)
    boundary = pmesh_dict_to_tag_dict(mesh_dict)

    return Mesh(vertex_coordinates, triangles,
                boundary=boundary,
                tagged_elements=tagged_elements,
                geo_reference=geo_reference,
                triangle_neighbours=triangle_neighbors,
                verbose=verbose)


def pmesh_to_domain_instance(source, DomainClass, use_cache=False,
                             verbose=False):
    """Convert a mesh file to a Domain instance.

    Parameters
    ----------
    source : str
        The name of the mesh file to convert, including the extension.
        Supported formats: .tsh, .msh
    DomainClass : type
        The Class that will be returned. Must be a subclass of Domain,
        with the same interface as domain.
    use_cache : bool, optional
        If True, caching is attempted for the computed domain.
        Default is False.
    verbose : bool, optional
        If True, print verbose output during processing.
        Default is False.
    Returns
    -------
    Domain
        An instance of DomainClass initialized from the mesh file.
    Notes
    -----
    When caching is enabled, the function uses the anuga.caching module
    to cache the computed domain based on the source file and DomainClass.
    """

    if use_cache is True:
        from anuga.caching import cache
        result = cache(_pmesh_to_domain_instance, (source, DomainClass),
                       dependencies=[source], verbose=verbose)
    else:
        result = _pmesh_to_domain_instance(*(source, DomainClass))

    return result


def _pmesh_to_domain_instance(source, DomainClass):
    """Converts a mesh file(.tsh or .msh), to a Domain instance.

    Internal function. See public interface pmesh_to_domain_instance for details
    """

    from anuga.abstract_2d_finite_volumes.generic_domain import Generic_Domain

    # ensure the required class is a subclass of Domain
    msg = ('The class %s is not a subclass of the generic domain class %s'
           % (DomainClass, Generic_Domain))
    assert issubclass(DomainClass, Generic_Domain), msg

    if type(source).__name__ == 'str':
        parm = {'file_name': source}
    else:
        parm = {'mesh_instance': source}

    (vertex_coordinates, vertices, tag_dict, vertex_quantity_dict,
     tagged_elements_dict, geo_reference) = pmesh_to_domain(**parm)


    domain = DomainClass(coordinates = vertex_coordinates,
                         vertices = vertices,
                         boundary = tag_dict,
                         tagged_elements = tagged_elements_dict,
                         geo_reference = geo_reference )

    # FIXME (Ole): Is this really the right place to apply a default
    # value specific to the shallow water wave equation?
    # The 'assert' above indicates that any subclass of Domain is acceptable.
    # Suggestion - module shallow_water.py will eventually take care of this
    # (when I get around to it) so it should be removed from here.

    # This doesn't work on the domain instance.
    # This is still needed so -ve elevations don't cuase 'lakes'
    # The fixme we discussed was to only create a quantity when its values
    # are set.
    # I think that's the way to go still

    # set the water stage to be the elevation
    if ('elevation' in vertex_quantity_dict and
        'stage' not in vertex_quantity_dict):
        vertex_quantity_dict['stage'] = vertex_quantity_dict['elevation']
    domain.set_quantity_vertices_dict(vertex_quantity_dict)

    return domain


def pmesh_to_domain(file_name=None, mesh_instance=None, use_cache=False,
                    verbose=False):
    """Convert a pmesh file or a pmesh mesh instance to a bunch of lists
    that can be used to instanciate a domain object.

    use_cache: True means that caching is attempted for the computed domain.
    """

    if verbose: log.info('Pmesh_to_Domain: Initialising')

    if use_cache is True:
        from anuga.caching import cache
        result = cache(_pmesh_to_domain, (file_name, mesh_instance),
                       dependencies=[file_name], verbose=verbose)

    else:
        result = _pmesh_to_domain(*(file_name, mesh_instance))

    if verbose: log.info('Pmesh_to_Domain: Done')

    return result


def _pmesh_to_domain(file_name=None, mesh_instance=None, use_cache=False,
                     verbose=False):
    """Convert a pmesh file or a pmesh mesh instance to a bunch of lists
    that can be used to instantiate a domain object.
    """

    from anuga.load_mesh.loadASCII import import_mesh_file

    # get data from mesh instance or file
    if file_name is None:
        mesh_dict = mesh_instance.mesh2io_dict()
    else:
        mesh_dict = import_mesh_file(file_name)

    # extract required data from the mesh dictionary
    vertex_coordinates = mesh_dict['vertices']
    volumes = mesh_dict['triangles']
    vertex_quantity_dict = {}

    # num.transpose(None) gives scalar array of value None
    point_atts = mesh_dict['vertex_attributes']

    point_titles = mesh_dict['vertex_attribute_titles']
    geo_reference = mesh_dict['geo_reference']
    if point_atts is not None:
        point_atts = num.transpose(point_atts)
        for quantity, value_vector in zip(point_titles, point_atts):
            vertex_quantity_dict[quantity] = value_vector
    tag_dict = pmesh_dict_to_tag_dict(mesh_dict)
    tagged_elements_dict = build_tagged_elements_dictionary(mesh_dict)

    return (vertex_coordinates, volumes, tag_dict, vertex_quantity_dict,
            tagged_elements_dict, geo_reference)


def build_tagged_elements_dictionary(mesh_dict):
    """Build the dictionary of element tags.

    tagged_elements is a dictionary of element arrays,
    keyed by tag: { (tag): [e1, e2, e3..] }
    """

    tri_atts = mesh_dict['triangle_tags']
    tagged_elements = {}
    if tri_atts is None:
       tagged_elements[''] = list(range(len(mesh_dict['triangles'])))
    else:
        for tri_att_index in range(len(tri_atts)):
            tagged_elements.setdefault(tri_atts[tri_att_index],
                                       []).append(tri_att_index)

    return tagged_elements


def pmesh_dict_to_tag_dict_old(mesh_dict):
    """ Convert the pmesh dictionary (mesh_dict) description of boundary tags
    to a dictionary of tags, indexed with volume id and face number.
    """

    triangles = mesh_dict['triangles']

    triangles = num.array(triangles,int)

    sides = {}
    for id, triangle in enumerate(triangles):
        a = triangle[0]
        b = triangle[1]
        c = triangle[2]

        sides[a,b] = 3*id+2 #(id, face)
        sides[b,c] = 3*id+0 #(id, face)
        sides[c,a] = 3*id+1 #(id, face)

    tag_dict = {}
    for seg, tag in zip(mesh_dict['segments'], mesh_dict['segment_tags']):
        v1 = int(seg[0])
        v2 = int(seg[1])
        for key in [(v1,v2),(v2,v1)]:
            if key in sides and tag != "":
                #"" represents null.  Don't put these into the dictionary
                #this creates a dict of lists of faces, indexed by tag
                #tagged_edges.setdefault(tag,[]).append(sides[key])
                vol_id = sides[key]//3
                edge_id = sides[key]%3
                tag_dict[vol_id,edge_id] = tag

    return tag_dict

def pmesh_dict_to_tag_dict(mesh_dict):
    """ Convert the pmesh dictionary (mesh_dict) description of boundary tags
    to a dictionary of tags, indexed with volume id and face number.
    """

    triangles = mesh_dict['triangles']
    segments = mesh_dict['segments']
    segment_tags = mesh_dict['segment_tags']

    triangles = num.array(triangles,int)
    segments = num.array(segments,int)
    tag_dict = {}

    #print triangles
    #print segments
    #print segment_tags

    from anuga.abstract_2d_finite_volumes.pmesh2domain_ext import build_boundary_dictionary

    segment_tags = [seg.encode() for seg in segment_tags]  # Convert to binary form
    tag_dict = build_boundary_dictionary(triangles, segments, segment_tags, tag_dict)


    for key in tag_dict.keys():
        x = tag_dict[key]
        tag_dict[key] = x.decode() # Convert back to string form

    return tag_dict


def calc_sides_old(triangles):
    """Build dictionary mapping from sides (2-tuple of points)
    to left hand side neighbouring triangle
    """

    sides = {}
    triangles = num.array(triangles,int)
    for id, triangle in enumerate(triangles):
        a = triangle[0]
        b = triangle[1]
        c = triangle[2]

        sides[a,b] = 3*id+2 #(id, face)
        sides[b,c] = 3*id+0 #(id, face)
        sides[c,a] = 3*id+1 #(id, face)

    return sides



def calc_sides_zip(triangles):
    """ Build dictionary mapping from sides (2-tuple of points)
    to left hand side neighbouring triangle
    """

    sides = {}


    triangles = num.array(triangles,int)


    a = triangles[:,0]
    b = triangles[:,1]
    c = triangles[:,2]

    id = 3*num.arange(len(triangles))

    sides.update(dict(list(zip(list(zip(a,b)),id+2))))
    sides.update(dict(list(zip(list(zip(b,c)),id+0))))
    sides.update(dict(list(zip(list(zip(c,a)),id+1))))

    return sides

def calc_sides_c(triangles):
    """Build dictionary mapping from sides (2-tuple of points)
    to left hand side neighbouring triangle
    """

    sides = {}


    triangles = num.array(triangles,int)
    ntriangles = len(triangles)

#    print 'calc_sides'
#    print type(triangles)

    #print ntriangles

    from anuga.abstract_2d_finite_volumes.pmesh2domain_ext import build_sides_dictionary
    sides = build_sides_dictionary(triangles, sides)


#    old_sides = {}
#
#    for id, triangle in enumerate(triangles):
#        a = int(triangle[0])
#        b = int(triangle[1])
#        c = int(triangle[2])
#
#        old_sides[a,b] = (id, 2) #(id, face)
#        old_sides[b,c] = (id, 0) #(id, face)
#        old_sides[c,a] = (id, 1) #(id, face)


    #print sides
    #print old_sides


    return sides

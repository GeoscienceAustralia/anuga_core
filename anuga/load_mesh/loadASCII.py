"""
The format for a Points dictionary is:

  ['pointlist'] a 2 column array describing points. 1st column x, 2nd column y.
  ['attributelist'], a dictionary of 1D arrays, representing attribute values
  at the point.  The dictionary key is the attribute header.
  ['geo_reference'] a Geo_refernece object. Use if the point information
        is relative. This is optional.
    eg
    dic['pointlist'] = [[1.0,2.0],[3.0,5.0]]
    dic['attributelist']['elevation'] = [[7.0,5.0]


    The dict format for IO with mesh files is;
    (the triangulation)
    vertices: [[x1,y1],[x2,y2],...] (lists of doubles)
    vertex_attributes: [[a11,a12,...],[a21,a22],...] (lists of doubles)
    vertex_attribute_titles:[A1Title, A2Title ...] (A list of strings)
    segments: [[v1,v2],[v3,v4],...] (lists of integers)
    segment_tags : [tag,tag,...] list of strings
    triangles : [(v1,v2,v3), (v4,v5,v6),....] lists of points
    triangle_tags: [s1,s2,...] A list of strings
    triangle_neighbors: [[t1,t2,t3], [t4,t5,t6],..] lists of triangles

    (the outline)
    points: [[x1,y1],[x2,y2],...] (lists of doubles)
    point_attributes: [[a11,a12,...],[a21,a22],...] (lists of doubles)
    outline_segments: [[point1,point2],[p3,p4],...] (lists of integers)
    outline_segment_tags : [tag1,tag2,...] list of strings
    holes : [[x1,y1],...](List of doubles, one inside each hole region)
    regions : [ [x1,y1],...] (List of 4 doubles)
    region_tags : [tag1,tag2,...] (list of strings)
    region_max_areas: [ma1,ma2,...] (A list of doubles)
    {Convension: A -ve max area means no max area}
    geo_reference: a Geo_refernece object. Use if the point information
        is relative. This is optional.

    Points files are .csv for ascii and .pts for NetCDF
    Mesh files are .tsh for ascii and .msh for NetCDF

    Note: point att's have no titles - that's currently ok, since least
    squares adds the point info to vertices, and they have titles

    The format for a .tsh file is (FIXME update this)

    First line:  <# of vertices> <# of attributes>
    Following lines:  <vertex #> <x> <y> [attributes]
    One line:  <# of triangles>
    Following lines: <triangle #> <vertex #>  <vertex #> <vertex #>
                         <neigbouring triangle #> <neigbouring triangle #>
                         <neigbouring triangle #> [attribute of region]
    One line:  <# of segments>
    Following lines:  <segment #> <vertex #>  <vertex #> [boundary tag]
"""

# FIXME (DSG-DSG) Is the dict format mentioned above a list of a numeric array?
#  Needs to be defined


#from string import find
from os.path import splitext



from anuga.coordinate_transforms.geo_reference import Geo_reference, TITLE, \
    TitleError
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a
from anuga.config import netcdf_float, netcdf_char, netcdf_int
from anuga.utilities.system_tools import *
import anuga.utilities.log as log

from anuga.file.netcdf import NetCDFFile
from anuga.utilities.numerical_tools import ensure_numeric

import numpy as num


class TitleAmountError(Exception):
    pass


NOMAXAREA = -999


def import_mesh_file(ofile):
    """Read a mesh file, either .tsh or .msh

    Note: will throw an IOError if it can't load the file.
    """

    try:
        if ofile[-4:] == ".tsh":
            dict = _read_tsh_file(ofile)
        elif ofile[-4:] == ".msh":
            dict = _read_msh_file(ofile)
        else:
            msg = 'Extension .%s is unknown' % ofile[-4:]
            raise IOError(msg)
    # FIXME No test for ValueError
    except (TitleError, SyntaxError, IndexError, ValueError):
        msg = 'File %s could not be opened' % ofile
        raise IOError(msg)

    return dict


def export_mesh_file(ofile, mesh_dict):
    """Write a mesh file given a dictionary.

    First line:       <# of vertices> <# of attributes>
    Following lines:  <vertex #> <x> <y> [attributes]
    One line:         <# of triangles>
    Following lines:  <triangle #> <vertex #>  <vertex #> <vertex #>
                          <neigbouring triangle #> <neigbouring triangle #>
                          <neigbouring triangle #> [attribute of region]
    One line:         <# of segments>
    Following lines:  <segment #> <vertex #>  <vertex #> [boundary tag]
    """

    # Ensure that if required key isn't present, we add it with [] value.
    # .setdefault() for dictionaries was added in python 2.0.
    # The only other neat way to do this is to override the dictionary
    # .__getitem__() method to return [] if key not found.
    reqd_keys = ['points', 'point_attributes', 'outline_segments',
                 'outline_segment_tags', 'holes', 'regions', 'region_tags',
                 'region_max_areas', 'vertices', 'vertex_attributes',
                 'vertex_attribute_titles', 'segments', 'segment_tags',
                 'triangles', 'triangle_tags', 'triangle_neighbors']

    for k in reqd_keys:
        mesh_dict.setdefault(k, [])

    # hand-off to appropriate function
    if (ofile[-4:] == ".tsh"):
        _write_tsh_file(ofile, mesh_dict)
    elif (ofile[-4:] == ".msh"):
        _write_msh_file(ofile, mesh_dict)
    else:
        msg = 'Unknown file type %s ' % ofile
        raise IOError(msg)


def _read_tsh_file(ofile):
    """Read the text file format for meshes"""

    fd = open(ofile, 'r')
    dict = _read_triangulation(fd)
    dict_mesh = _read_outline(fd)
    for element in list(dict_mesh.keys()):
        dict[element] = dict_mesh[element]
    fd.close()

    return dict


def _read_triangulation(fd):
    """Read the generated triangulation, NOT the outline."""

    delimiter = " "

    # loading the point info
    line = fd.readline()
    fragments = line.split()
    if fragments == []:
        NumOfVertices = 0
        NumOfVertAttributes = 0
    else:
        NumOfVertices = fragments[0]
        NumOfVertAttributes = fragments[1]
    points = []
    pointattributes = []
    for index in range(int(NumOfVertices)):
        fragments = fd.readline().split()
        fragments.pop(0)        # pop off the index

        # pop the x & y off so we're left with a list of attributes
        vert = [float(fragments.pop(0)), float(fragments.pop(0))]
        points.append(vert)
        apointattributes = []
        for fragment in fragments:
            apointattributes.append(float(fragment))
        if apointattributes != []:
            pointattributes.append(apointattributes)

    # loading the point title info
    line = fd.readline()
    vertTitle = []
    for index in range(int(NumOfVertAttributes)):
        fragments = fd.readline().strip()
        vertTitle.append(fragments)

    # loading the triangle info
    line = fd.readline()
    fragments = line.split()
    NumOfTriangles = fragments[0]
    triangles = []
    triangleattributes = []
    triangleneighbors = []
    for index in range(int(NumOfTriangles)):
        line = fd.readline()
        line.strip()            # so we can get the region string
        fragments = line.split()
        fragments.pop(0)        # pop off the index

        tri = [int(fragments[0]), int(fragments[1]), int(fragments[2])]
        triangles.append(tri)
        neighbors = [int(fragments[3]), int(fragments[4]), int(fragments[5])]
        triangleneighbors.append(neighbors)
        for x in range(7):  # remove index [<vertex #>] [<neigbouring tri #>]
            line = line[line.find(delimiter):]     # remove index
            line = line.lstrip()
        stringtag = line.strip()
        triangleattributes.append(stringtag)

    # loading the segment info
    line = fd.readline()
    fragments = line.split()
    NumOfSegments = fragments[0]
    segments = []
    segmenttags = []
    for index in range(int(NumOfSegments)):
        line = fd.readline()
        line.strip()            # to get the segment string
        fragments = line.split()
        fragments.pop(0)  # pop off the index
        seg = [int(fragments[0]), int(fragments[1])]
        segments.append(seg)
        line = line[line.find(delimiter):]  # remove index
        line = line.lstrip()
        line = line[line.find(delimiter):]  # remove x
        line = line.lstrip()
        line = line[line.find(delimiter):]  # remove y
        stringtag = line.strip()
        segmenttags.append(stringtag)

    meshDict = {}
    meshDict['vertices'] = points
    if pointattributes == []:
        meshDict['vertex_attributes'] = None
    else:
        meshDict['vertex_attributes'] = pointattributes
    meshDict['triangles'] = triangles
    meshDict['triangle_tags'] = triangleattributes
    meshDict['triangle_neighbors'] = triangleneighbors
    meshDict['segments'] = segments
    meshDict['segment_tags'] = segmenttags
    meshDict['vertex_attribute_titles'] = vertTitle

    return meshDict


def _read_outline(fd):
    """Read a mesh file outline.

    Note, if a file has no mesh info, it can still be read -
    the meshdic returned will be 'empty'.
    """

    delimiter = " "     # warning: split() calls are using default whitespace

    # loading the point info
    line = fd.readline()
    fragments = line.split()
    if fragments == []:
        NumOfVertices = 0
        NumOfVertAttributes = 0
    else:
        NumOfVertices = fragments[0]
        NumOfVertAttributes = fragments[1]
    points = []
    pointattributes = []
    for index in range(int(NumOfVertices)):
        fragments = fd.readline().split()
        fragments.pop(0)                # pop off the index
        # pop the x & y off so we're left with a list of attributes
        vert = [float(fragments.pop(0)), float(fragments.pop(0))]
        points.append(vert)
        apointattributes = []
        for fragment in fragments:
            apointattributes.append(float(fragment))
        pointattributes.append(apointattributes)

    # loading the segment info
    line = fd.readline()
    fragments = line.split()
    if fragments == []:
        NumOfSegments = 0
    else:
        NumOfSegments = fragments[0]
    segments = []
    segmenttags = []
    for index in range(int(NumOfSegments)):
        line = fd.readline()
        fragments = line.split()
        fragments.pop(0)                    # pop off the index
        seg = [int(fragments[0]), int(fragments[1])]
        segments.append(seg)
        line = line[line.find(delimiter):]  # remove index
        line = line.lstrip()
        line = line[line.find(delimiter):]  # remove x
        line = line.lstrip()
        line = line[line.find(delimiter):]  # remove y
        stringtag = line.strip()
        segmenttags.append(stringtag)

    # loading the hole info
    line = fd.readline()
    fragments = line.split()
    if fragments == []:
        numOfHoles = 0
    else:
        numOfHoles = fragments[0]
    holes = []
    for index in range(int(numOfHoles)):
        fragments = fd.readline().split()
        fragments.pop(0)  # pop off the index
        hole = [float(fragments[0]), float(fragments[1])]
        holes.append(hole)

    # loading the region info
    line = fd.readline()
    fragments = line.split()
    if fragments == []:
        numOfRegions = 0
    else:
        numOfRegions = fragments[0]
    regions = []
    regionattributes = []
    for index in range(int(numOfRegions)):
        line = fd.readline()
        fragments = line.split()
        fragments.pop(0)                    # pop off the index
        region = [float(fragments[0]), float(fragments[1])]
        regions.append(region)

        line = line[line.find(delimiter):]  # remove index
        line = line.lstrip()
        line = line[line.find(delimiter):]  # remove x
        line = line.lstrip()
        line = line[line.find(delimiter):]  # remove y
        stringtag = line.strip()
        regionattributes.append(stringtag)

    regionmaxareas = []
    line = fd.readline()
    for index in range(int(numOfRegions)):  # Read in the Max area info
        line = fd.readline()
        fragments = line.split()
        # The try is here for format compatibility
        try:
            fragments.pop(0)                # pop off the index
            if len(fragments) == 0:         # no max area
                regionmaxareas.append(None)
            else:
                regionmaxareas.append(float(fragments[0]))
        except (ValueError, IndexError) as e:
            regionmaxareas.append(None)

    try:
        geo_reference = Geo_reference(ASCIIFile=fd)
    except:
        # geo_ref not compulsory
        geo_reference = None

    meshDict = {}
    meshDict['points'] = points
    meshDict['point_attributes'] = pointattributes
    meshDict['outline_segments'] = segments
    meshDict['outline_segment_tags'] = segmenttags
    meshDict['holes'] = holes
    meshDict['regions'] = regions
    meshDict['region_tags'] = regionattributes
    meshDict['region_max_areas'] = regionmaxareas
    meshDict['geo_reference'] = geo_reference

    return meshDict


def _write_ASCII_triangulation(fd, gen_dict):
    vertices = gen_dict['vertices']
    vertices_attributes = gen_dict['vertex_attributes']

    try:
        vertices_attribute_titles = gen_dict['vertex_attribute_titles']
    except KeyError as e:
        # FIXME is this the best way?
        if vertices_attributes == [] or vertices_attributes[0] == []:
            vertices_attribute_titles = []
        else:
            raise KeyError(e)

    triangles = gen_dict['triangles']
    triangles_attributes = gen_dict['triangle_tags']
    triangle_neighbors = gen_dict['triangle_neighbors']

    segments = gen_dict['segments']
    segment_tags = gen_dict['segment_tags']

    numVert = str(len(vertices))
    # Don't understand why we have to do vertices_attributes[0] is None,
    # but it doesn't work otherwise...
    if (vertices_attributes is None or
        numVert == "0" or
            len(vertices_attributes) == 0):
        numVertAttrib = "0"
    else:
        #numVertAttrib = str(len(vertices_attributes[0]))
        try:
            numVertAttrib = str(len(vertices_attributes[0]))
        except TypeError:
            vertices_attributes = [[entry] for entry in vertices_attributes]
            numVertAttrib = str(len(vertices_attributes[0]))

    fd.write(numVert + " " + numVertAttrib +
             " # <# of verts> <# of vert attributes>, next lines <vertex #> "
             "<x> <y> [attributes] ...Triangulation Vertices...\n")

    # <vertex #> <x> <y> [attributes]
    index = 0
    for vert in vertices:
        attlist = ""

        if vertices_attributes is None or vertices_attributes == []:
            attlist = ""
        else:
            for att in vertices_attributes[index]:
                attlist = attlist + str(att) + " "
        attlist.strip()
        fd.write(str(index) + " " + str(vert[0]) + " " + str(vert[1])
                 + " " + attlist + "\n")
        index += 1

    # write comments for title
    fd.write("# attribute column titles ...Triangulation Vertex Titles...\n")
    for title in vertices_attribute_titles:
        fd.write(title + "\n")

    # <# of triangles>
    n = len(triangles)
    fd.write(str(n)
             + " # <# of triangles>, next lines <triangle #> [<vertex #>] "
               "[<neigbouring triangle #>] [attribute of region] ..."
               "Triangulation Triangles...\n")

    # <triangle #> <vertex #>  <vertex #> <vertex #> <neigbouring triangle #>
    #    <neigbouring triangle #> <neigbouring triangle #> [attribute of region]
    for index in range(n):
        neighbors = ""
        tri = triangles[index]
        if len(triangle_neighbors) == 0:
            neighbors = "-1 -1 -1 "
        else:
            for neighbor in triangle_neighbors[index]:
                if neighbor:
                    neighbors += str(neighbor) + " "
                else:
                    if neighbor == 0:
                        neighbors += "0 "
                    else:
                        neighbors += "-1 "
        # Warning even though a list is passed, only the first value
        # is written.  There's an assumption that the list only
        # contains one item. This assumption is made since the
        # dict that's being passed around is also be used to communicate
        # with triangle, and it seems to have the option of returning
        # more than one value for triangle attributex
        if (triangles_attributes is None or
            triangles_attributes == [] or
                triangles_attributes[index] == ['']):
            att = ""
        else:
            att = str(triangles_attributes[index])
        fd.write(str(index) + " "
                 + str(tri[0]) + " "
                 + str(tri[1]) + " "
                 + str(tri[2]) + " "
                 + neighbors + " "
                 + att + "\n")

    # One line:  <# of segments>
    fd.write(str(len(segments)) +
             " # <# of segments>, next lines <segment #> <vertex #>  "
             "<vertex #> [boundary tag] ...Triangulation Segments...\n")

    # Following lines:  <segment #> <vertex #>  <vertex #> [boundary tag]
    for i in range(len(segments)):
        seg = segments[i]
        fd.write(str(i) + " "
                 + str(seg[0]) + " "
                 + str(seg[1]) + " "
                 + str(segment_tags[i]) + "\n")


def _write_tsh_file(ofile, mesh_dict):
    fd = open(ofile, 'w')
    _write_ASCII_triangulation(fd, mesh_dict)
    _write_ASCII_outline(fd, mesh_dict)
    fd.close()


def _write_ASCII_outline(fd, dict):
    points = dict['points']
    point_attributes = dict['point_attributes']
    segments = dict['outline_segments']
    segment_tags = dict['outline_segment_tags']
    holes = dict['holes']
    regions = dict['regions']
    region_tags = dict['region_tags']
    region_max_areas = dict['region_max_areas']

    num_points = str(len(points))
    if (num_points == '0'):
        num_point_atts = '0'
    else:
        num_point_atts = str(len(point_attributes[0]))

    fd.write(num_points + ' ' + num_point_atts +
             ' # <# of verts> <# of vert attributes>, next lines <vertex #> '
             '<x> <y> [attributes] ...Mesh Vertices...\n')

    # <x> <y> [attributes]
    for i, point in enumerate(points):
        attlist = ''
        for att in point_attributes[i]:
            attlist = attlist + str(att) + ' '
        attlist.strip()
        fd.write(str(i) + ' ' + str(point[0]) + ' ' + str(point[1]) + ' ' +
                 attlist + '\n')

    # One line:  <# of segments>
    fd.write(str(len(segments)) +
             ' # <# of segments>, next lines <segment #> <vertex #>  '
             '<vertex #> [boundary tag] ...Mesh Segments...\n')

    # Following lines:  <vertex #>  <vertex #> [boundary tag]
    for i, seg in enumerate(segments):
        fd.write(str(i) + ' ' + str(seg[0]) + ' ' + str(seg[1]) + ' ' +
                 str(segment_tags[i]) + '\n')

    # One line:  <# of holes>
    fd.write(str(len(holes)) +
             ' # <# of holes>, next lines <Hole #> <x> <y> ...Mesh Holes...\n')
    # <x> <y>
    for i, h in enumerate(holes):
        fd.write(str(i) + ' ' + str(h[0]) + ' ' + str(h[1]) + '\n')

    # One line:  <# of regions>
    fd.write(str(len(regions)) +
             ' # <# of regions>, next lines <Region #> <x> <y> <tag>'
             '...Mesh Regions...\n')

    # <index> <x> <y> <tag>
    for i, r in enumerate(regions):
        fd.write(str(i) + ' ' + str(r[0]) + ' ' + str(r[1]) + ' ' +
                 str(region_tags[i]) + '\n')

    # <index> [<MaxArea>|'']

    # One line:  <# of regions>
    fd.write(str(len(regions)) +
             ' # <# of regions>, next lines <Region #> [Max Area] '
             '...Mesh Regions...\n')
    for i, r in enumerate(regions):
        area = str(region_max_areas[i])

        fd.write(str(i) + ' ' + area + '\n')

    # geo_reference info
    if 'geo_reference' in dict and not dict['geo_reference'] is None:
        dict['geo_reference'].write_ASCII(fd)


def _write_msh_file(file_name, mesh):
    """Write .msh NetCDF file

    WARNING: This function mangles the mesh data structure
    """

    # FIXME(Ole and John): We ran into a problem on Bogong (64 bit)
    # where integers appeared as arrays.  This may be similar to
    # problem seen by Steve in changeset:2778 where he had to wrap
    # them in int.  Now we are trying with the native Integer format
    # (Int == 'l' == Int64). However, that caused casting errors, when
    # 64bit arrays are to be assigned to their NetCDF counterparts. It
    # seems that the NetCDF arrays are 32bit even though they are
    # created with the type Int64. Need to look at the NetCDF library
    # in more detail.

    IntType = num.int32
    #IntType = Int

    # print 'mesh vertices',mesh['vertices'].shape

    # the triangulation
    mesh['vertices'] = num.array(mesh['vertices'], float)
    mesh['vertex_attribute_titles'] = \
        num.array(string_to_char(
            mesh['vertex_attribute_titles']), 'S1')

    num_attributes = len(mesh['vertex_attribute_titles'])
    num_vertices = mesh['vertices'].shape[0]
    # print 'num_attrib ',num_attributes
    if mesh['vertex_attributes'] != None:
        mesh['vertex_attributes'] = \
            num.array(mesh['vertex_attributes'], float)

    if num_attributes > 0:
        mesh['vertex_attributes'] = \
            num.reshape(mesh['vertex_attributes'], (num_vertices, -1))

    mesh['segments'] = num.array(mesh['segments'], IntType)
    mesh['segment_tags'] = num.array(string_to_char(mesh['segment_tags']),
                                     'S1')
    mesh['triangles'] = num.array(mesh['triangles'], IntType)
    mesh['triangle_tags'] = num.array(string_to_char(mesh['triangle_tags']),
                                      'S1')
    mesh['triangle_neighbors'] = \
        num.array(mesh['triangle_neighbors'], IntType)

    # the outline
    mesh['points'] = num.array(mesh['points'], float)
    mesh['point_attributes'] = num.array(mesh['point_attributes'], float)
    mesh['outline_segments'] = num.array(mesh['outline_segments'], IntType)
    mesh['outline_segment_tags'] = \
        num.array(string_to_char(mesh['outline_segment_tags']), 'S1')
    mesh['holes'] = num.array(mesh['holes'], float)
    mesh['regions'] = num.array(mesh['regions'], float)
    mesh['region_tags'] = num.array(
        string_to_char(mesh['region_tags']), 'S1')
    mesh['region_max_areas'] = num.array(mesh['region_max_areas'], float)

    # NetCDF file definition
    try:
        outfile = NetCDFFile(file_name, netcdf_mode_w)
    except IOError:
        msg = 'File %s could not be created' % file_name
        raise Exception(msg)

    # Create new file
    outfile.institution = 'Geoscience Australia'
    outfile.description = 'NetCDF format for compact and portable storage ' + \
                          'of spatial point data'

    # dimension definitions - fixed
    outfile.createDimension('num_of_dimensions', 2)     # This is 2d data
    outfile.createDimension('num_of_segment_ends', 2)   # Segs have two points
    outfile.createDimension('num_of_triangle_vertices', 3)
    outfile.createDimension('num_of_triangle_faces', 3)
    outfile.createDimension('num_of_region_max_area', 1)

    # Create dimensions, variables and set the variables

    # trianglulation
    # vertices
    if (mesh['vertices'].shape[0] > 0):
        outfile.createDimension('num_of_vertices', mesh['vertices'].shape[0])
        outfile.createVariable('vertices', netcdf_float, ('num_of_vertices',
                                                          'num_of_dimensions'))
        outfile.variables['vertices'][:] = mesh['vertices']

        # print 'mesh vertex attributes', mesh['vertex_attributes'].shape

        if (mesh['vertex_attributes'] is not None and
            (mesh['vertex_attributes'].shape[0] > 0 and
             mesh['vertex_attributes'].shape[1] > 0)):
            outfile.createDimension('num_of_vertex_attributes',
                                    mesh['vertex_attributes'].shape[1])
            outfile.createDimension('num_of_vertex_attribute_title_chars',
                                    mesh['vertex_attribute_titles'].shape[1])
            outfile.createVariable('vertex_attributes',
                                   netcdf_float,
                                   ('num_of_vertices',
                                    'num_of_vertex_attributes'))
            outfile.createVariable('vertex_attribute_titles',
                                   netcdf_char,
                                   ('num_of_vertex_attributes',
                                    'num_of_vertex_attribute_title_chars'))
            outfile.variables['vertex_attributes'][:] = \
                mesh['vertex_attributes']
            outfile.variables['vertex_attribute_titles'][:] = \
                mesh['vertex_attribute_titles']

    # segments
    if (mesh['segments'].shape[0] > 0):
        outfile.createDimension('num_of_segments', mesh['segments'].shape[0])
        outfile.createVariable('segments', netcdf_int,
                               ('num_of_segments', 'num_of_segment_ends'))
        outfile.variables['segments'][:] = mesh['segments']
        if (mesh['segment_tags'].shape[1] > 0):
            outfile.createDimension('num_of_segment_tag_chars',
                                    mesh['segment_tags'].shape[1])
            outfile.createVariable('segment_tags',
                                   netcdf_char,
                                   ('num_of_segments',
                                    'num_of_segment_tag_chars'))
            outfile.variables['segment_tags'][:] = mesh['segment_tags']

    # triangles
    if (mesh['triangles'].shape[0] > 0):
        outfile.createDimension('num_of_triangles', mesh['triangles'].shape[0])
        outfile.createVariable('triangles', netcdf_int,
                               ('num_of_triangles', 'num_of_triangle_vertices'))
        outfile.createVariable('triangle_neighbors', netcdf_int,
                               ('num_of_triangles', 'num_of_triangle_faces'))
        outfile.variables['triangles'][:] = mesh['triangles']
        outfile.variables['triangle_neighbors'][:] = mesh['triangle_neighbors']
        if (mesh['triangle_tags'] is not None and
                (mesh['triangle_tags'].shape[1] > 0)):
            outfile.createDimension('num_of_triangle_tag_chars',
                                    mesh['triangle_tags'].shape[1])
            outfile.createVariable('triangle_tags', netcdf_char,
                                   ('num_of_triangles',
                                    'num_of_triangle_tag_chars'))
            outfile.variables['triangle_tags'][:] = mesh['triangle_tags']

    # outline
    # points
    if (mesh['points'].shape[0] > 0):
        outfile.createDimension('num_of_points', mesh['points'].shape[0])
        outfile.createVariable('points', netcdf_float,
                               ('num_of_points', 'num_of_dimensions'))
        outfile.variables['points'][:] = mesh['points']
        if mesh['point_attributes'].shape[0] > 0  \
           and mesh['point_attributes'].shape[1] > 0:
            outfile.createDimension('num_of_point_attributes',
                                    mesh['point_attributes'].shape[1])
            outfile.createVariable('point_attributes', netcdf_float,
                                   ('num_of_points', 'num_of_point_attributes'))
            outfile.variables['point_attributes'][:] = mesh['point_attributes']

    # outline_segments
    if mesh['outline_segments'].shape[0] > 0:
        outfile.createDimension('num_of_outline_segments',
                                mesh['outline_segments'].shape[0])
        outfile.createVariable('outline_segments', netcdf_int,
                               ('num_of_outline_segments',
                                'num_of_segment_ends'))
        outfile.variables['outline_segments'][:] = mesh['outline_segments']
        if mesh['outline_segment_tags'].shape[1] > 0:
            outfile.createDimension('num_of_outline_segment_tag_chars',
                                    mesh['outline_segment_tags'].shape[1])
            outfile.createVariable('outline_segment_tags', netcdf_char,
                                   ('num_of_outline_segments',
                                    'num_of_outline_segment_tag_chars'))
            outfile.variables['outline_segment_tags'][:] = \
                mesh['outline_segment_tags']

    # holes
    if (mesh['holes'].shape[0] > 0):
        outfile.createDimension('num_of_holes', mesh['holes'].shape[0])
        outfile.createVariable('holes', netcdf_float,
                               ('num_of_holes', 'num_of_dimensions'))
        outfile.variables['holes'][:] = mesh['holes']

    # regions
    if (mesh['regions'].shape[0] > 0):
        outfile.createDimension('num_of_regions', mesh['regions'].shape[0])
        outfile.createVariable('regions', netcdf_float,
                               ('num_of_regions', 'num_of_dimensions'))
        outfile.createVariable('region_max_areas', netcdf_float,
                               ('num_of_regions',))
        outfile.variables['regions'][:] = mesh['regions']
        outfile.variables['region_max_areas'][:] = mesh['region_max_areas']
        if (mesh['region_tags'].shape[1] > 0):
            outfile.createDimension('num_of_region_tag_chars',
                                    mesh['region_tags'].shape[1])
            outfile.createVariable('region_tags', netcdf_char,
                                   ('num_of_regions',
                                    'num_of_region_tag_chars'))
            outfile.variables['region_tags'][:] = mesh['region_tags']

    # geo_reference info
    if 'geo_reference' in mesh and not mesh['geo_reference'] is None:
        mesh['geo_reference'].write_NetCDF(outfile)

    outfile.close()


def _read_msh_file(file_name):
    """ Read in an msh file."""

    # Check contents.  Get NetCDF
    fd = open(file_name, 'r')
    fd.close()

    # throws prints to screen if file not present
    fid = NetCDFFile(file_name, netcdf_mode_r)
    mesh = {}

    # Get the variables - the triangulation
    try:
        mesh['vertices'] = fid.variables['vertices'][:]
    except KeyError:
        mesh['vertices'] = num.array([], int)  # array default#

    try:
        mesh['vertex_attributes'] = fid.variables['vertex_attributes'][:]
    except KeyError:
        mesh['vertex_attributes'] = None

    mesh['vertex_attribute_titles'] = []
    try:
        titles = fid.variables['vertex_attribute_titles'][:]
        mesh['vertex_attribute_titles'] = [x.tobytes().decode().strip()
                                           for x in titles]
    except KeyError:
        pass

    try:
        mesh['segments'] = fid.variables['segments'][:]
    except KeyError:
        mesh['segments'] = num.array([], int)  # array default#

    mesh['segment_tags'] = []
    try:
        tags = fid.variables['segment_tags'][:]
        mesh['segment_tags'] = [x.tobytes().decode().strip() for x in tags]
    except KeyError:
        for ob in mesh['segments']:
            mesh['segment_tags'].append('')
            
    try:
        mesh['triangles'] = fid.variables['triangles'][:]
        mesh['triangle_neighbors'] = fid.variables['triangle_neighbors'][:]
    except KeyError:
        mesh['triangles'] = num.array([], int)  # array default#
        mesh['triangle_neighbors'] = num.array([], int)  # array default#


    mesh['triangle_tags'] = []
    try:
        tags = fid.variables['triangle_tags'][:]
        mesh['triangle_tags'] = [x.tobytes().decode().strip() for x in tags]
    except KeyError:
        for ob in mesh['triangles']:
            mesh['triangle_tags'].append('')

    # the outline
    try:
        mesh['points'] = fid.variables['points'][:]
    except KeyError:
        mesh['points'] = []


    try:
        mesh['point_attributes'] = fid.variables['point_attributes'][:]
    except KeyError:
        mesh['point_attributes'] = []
        for point in mesh['points']:
            mesh['point_attributes'].append([])

    try:
        mesh['outline_segments'] = fid.variables['outline_segments'][:]
    except KeyError:
        mesh['outline_segments'] = num.array([], int)  # array default#

    mesh['outline_segment_tags'] = []
    try:
        tags = fid.variables['outline_segment_tags'][:]
        for i, tag in enumerate(tags):
            mesh['outline_segment_tags'].append(tags[i].tobytes().decode().strip())
    except KeyError:
        for ob in mesh['outline_segments']:
            mesh['outline_segment_tags'].append('')

    try:
        mesh['holes'] = fid.variables['holes'][:]
    except KeyError:
        mesh['holes'] = num.array([], int)  # array default#

    try:
        mesh['regions'] = fid.variables['regions'][:]
    except KeyError:
        mesh['regions'] = num.array([], int)  # array default#

    mesh['region_tags'] = []
    try:
        tags = fid.variables['region_tags'][:]
        for i, tag in enumerate(tags):
            mesh['region_tags'].append(tags[i].tobytes().decode().strip())
    except KeyError:
        for ob in mesh['regions']:
            mesh['region_tags'].append('')

    try:
        mesh['region_max_areas'] = fid.variables['region_max_areas'][:]
    except KeyError:
        mesh['region_max_areas'] = num.array([], int)  # array default#
        
    try:
        geo_reference = Geo_reference(NetCDFObject=fid)
        mesh['geo_reference'] = geo_reference
    except AttributeError as e:
        # geo_ref not compulsory
        mesh['geo_reference'] = None

    fid.close()

    return mesh


def export_boundary_file(file_name, points, title, delimiter=','):
    """Export a boundary file.

    Format:
    First line: Title variable
    Following lines:  [point index][delimiter][point index]

    file_name - the name of the new file
    points - List of point index pairs [[p1, p2],[p3, p4]..]
    title - info to write in the first line
    """

    fd = open(file_name, 'w')

    fd.write(title + '\n')

    # [point index][delimiter][point index]
    for point in points:
        fd.write(str(point[0]) + delimiter + str(point[1]) + '\n')

    fd.close()

################################################################################
#  IMPORT/EXPORT POINTS FILES
################################################################################


def extent_point_atts(point_atts):
    """Returns 4 points representing the extent
    This loses attribute info.
    """

    point_atts['pointlist'] = extent(point_atts['pointlist'])
    point_atts['attributelist'] = {}

    return point_atts


def extent(points):
    points = num.array(points, float)

    max_x = min_x = points[0][0]
    max_y = min_y = points[0][1]

    for point in points[1:]:
        x = point[0]
        if x > max_x:
            max_x = x
        if x < min_x:
            min_x = x

        y = point[1]
        if y > max_y:
            max_y = y
        if y < min_y:
            min_y = y

    extent = num.array([[min_x, min_y],
                        [max_x, min_y],
                        [max_x, max_y],
                        [min_x, max_y]])

    return extent


def reduce_pts(infile, outfile, max_points, verbose=False):
    """Reduce a points file until less than given size.

    Reduces a points file by removing every second point until the # of points
    is less than max_points.
    """

    # check out pts2rectangular in least squares, and the use of reduction.
    # Maybe it does the same sort of thing?
    point_atts = _read_pts_file(infile)

    while point_atts['pointlist'].shape[0] > max_points:
        if verbose:
            log.critical("point_atts['pointlist'].shape[0]")
        point_atts = half_pts(point_atts)

    export_points_file(outfile, point_atts)


def produce_half_point_files(infile, max_points, delimiter, verbose=False):
    point_atts = _read_pts_file(infile)
    root, ext = splitext(infile)
    outfiles = []

    if verbose:
        log.critical("# of points", point_atts['pointlist'].shape[0])

    while point_atts['pointlist'].shape[0] > max_points:
        point_atts = half_pts(point_atts)

        if verbose:
            log.critical("# of points = %s"
                         % str(point_atts['pointlist'].shape[0]))

        outfile = root + delimiter + \
            str(point_atts['pointlist'].shape[0]) + ext
        outfiles.append(outfile)
        export_points_file(outfile, point_atts)

    return outfiles


def point_atts2array(point_atts):
    # convert attribute list to array of floats
    point_atts['pointlist'] = num.array(point_atts['pointlist'], float)

    for key in list(point_atts['attributelist'].keys()):
        point_atts['attributelist'][key] = \
            num.array(point_atts['attributelist'][key], float)

    return point_atts


def half_pts(point_atts):
    point_atts2array(point_atts)
    point_atts['pointlist'] = point_atts['pointlist'][::2]

    for key in list(point_atts['attributelist'].keys()):
        point_atts['attributelist'][key] = point_atts['attributelist'][key][::2]

    return point_atts


def concatinate_attributelist(dic):
    """
    giving a dic[attribute title] = attribute
    return list of attribute titles, array of attributes
    """

    point_attributes = num.array([], float)
    keys = list(dic.keys())
    key = keys.pop(0)
    point_attributes = num.reshape(dic[key], (dic[key].shape[0], 1))
    for key in keys:
        reshaped = num.reshape(dic[key], (dic[key].shape[0], 1))
        point_attributes = num.concatenate(
            [point_attributes, reshaped], axis=1)

    return list(dic.keys()), point_attributes


def take_points(dict, indices_to_keep):
    dict = point_atts2array(dict)
    # FIXME maybe the points data structure should become a class?
    dict['pointlist'] = num.take(dict['pointlist'], indices_to_keep, axis=0)

    for key in list(dict['attributelist'].keys()):
        dict['attributelist'][key] = num.take(dict['attributelist'][key],
                                              indices_to_keep, axis=0)

    return dict


def add_point_dictionaries(dict1, dict2):
    """
    """

    dict1 = point_atts2array(dict1)
    dict2 = point_atts2array(dict2)

    combined = {}
    combined['pointlist'] = num.concatenate((dict2['pointlist'],
                                             dict1['pointlist']), axis=0)

    atts = {}
    for key in list(dict2['attributelist'].keys()):
        atts[key] = num.concatenate((dict2['attributelist'][key],
                                     dict1['attributelist'][key]), axis=0)
    combined['attributelist'] = atts
    combined['geo_reference'] = dict1['geo_reference']

    return combined


if __name__ == "__main__":
    pass

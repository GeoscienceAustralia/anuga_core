#!/usr/bin/env python
#
"""General 2D triangular classes for triangular mesh generation.

   Note: A .index attribute is added to objects such as vertices and
   segments, often when reading and writing to files.  This information
   should not be used as persistant information.  It is not the 'index' of
   an element in a list.

   
   Copyright 2003/2004
   Ole Nielsen, Stephen Roberts, Duncan Gray, Christopher Zoppou
   Geoscience Australia
"""

import sys
import math
import re
import os
import pickle

import numpy as num


from anuga.coordinate_transforms.geo_reference import Geo_reference, \
    DEFAULT_ZONE
from anuga.load_mesh.loadASCII import NOMAXAREA, export_mesh_file, \
    import_mesh_file
import anuga.alpha_shape.alpha_shape
from anuga.geospatial_data.geospatial_data import Geospatial_data, \
    ensure_geospatial, ensure_absolute, ensure_numeric
from anuga.mesh_engine.mesh_engine import generate_mesh
import anuga.utilities.log as log

from anuga.file.ungenerate import load_ungenerate

try:
    import kinds
except ImportError:
    # Hand-built mockup of the things we need from the kinds package, since it
    # was recently removed from the standard numeric distro.  Some users may
    # not have it by default.
    class _bunch(object):
        pass

    class _kinds(_bunch):
        default_float_kind = _bunch()
        default_float_kind.MIN = 2.2250738585072014e-308  # smallest +ve number
        default_float_kind.MAX = 1.7976931348623157e+308

    kinds = _kinds()


# 1st and third values must be the same
# FIXME: maybe make this a switch that the user can change? - DSG
initialconversions = ['', 'exterior', '']
SEG_COLOUR = 'blue'


class MeshObject(object):
    """
    An abstract superclass for the basic mesh objects, eg vertex, segment,
    triangle.
    """

    def __init__(self):
        pass


class Point(MeshObject):
    """
    Define a point in a 2D space.
    """

    def __init__(self, X, Y):
        __slots__ = ['x', 'y']
        self.x = X
        self.y = Y

    def DistanceToPoint(self, OtherPoint):
        """
        Returns the distance from this point to another
        """
        SumOfSquares = ((self.x - OtherPoint.x)**2) + \
                       ((self.y - OtherPoint.y)**2)
        return math.sqrt(SumOfSquares)

    def IsInsideCircle(self, Center, Radius):
        """
        Return 1 if this point is inside the circle, 
        0 otherwise
        """

        if (self.DistanceToPoint(Center) < Radius):
            return 1
        else:
            return 0

    def __repr__(self):
        return "(%f,%f)" % (self.x, self.y)

    def cmp_xy(self, point):
        if self.x < point.x:
            return -1
        elif self.x > point.x:
            return 1
        else:
            if self.y < point.y:
                return -1
            elif self.y > point.y:
                return 1
            else:
                return 0

    def same_x_y(self, point):
        if self.x == point.x and self.y == point.y:
            return True
        else:
            return False


class Vertex(Point):
    """A point on the mesh.
    Object attributes based on the Triangle program
    """

    VERTEXSQUARESIDELENGTH = 6

    def __init__(self, X, Y, attributes=None):
        __slots__ = ['x', 'y', 'attributes']

        # we don't care what we get, as long as we can get float *value*
        self.x = float(X)
        self.y = float(Y)

        self.attributes = []
        if not attributes is None:
            self.attributes = attributes

    def setAttributes(self, attributes):
        """attributes is a list. """

        self.attributes = attributes

    def draw(self, canvas, tags, colour='black', scale=1, xoffset=0,
             yoffset=0, ):
        x = scale*(self.x + xoffset)
        y = -1*scale*(self.y + yoffset)  # - since for a canvas - is up
        cornerOffset = self.VERTEXSQUARESIDELENGTH/ 2

        # A hack to see the vert tags
        # note: there will be many tags, since tags will not be removed
        # when zooming
        # canvas.create_text(x+ 2*cornerOffset,
        #                   y+ 2*cornerOffset,
        #                        text=tags)

        # This gives points info.  It is a mess though, since numbers are
        # not deleted when zaooming in.
        # canvas.create_text(x+ 2*cornerOffset,
        #                   y+ 2*cornerOffset,
        #                        text=str(x)+','+str(y))

        return canvas.create_rectangle(x-cornerOffset,
                                       y-cornerOffset,
                                       x+cornerOffset,
                                       y+cornerOffset,
                                       tags=tags,
                                       outline=colour,
                                       fill='white')

        # return tags

    def __repr__(self):
        return "[(%f,%f),%r]" % (self.x, self.y, self.attributes)


class Hole(Point):
    """A region of the mesh were no triangles are generated.
    Defined by a point in the hole enclosed by segments.
    """

    HOLECORNERLENGTH = 6

    def draw(self, canvas, tags, colour='purple', scale=1, xoffset=0,
             yoffset=0, ):
        x = scale*(self.x + xoffset)
        y = -1*scale*(self.y + yoffset)  # - since for a canvas - is up
        cornerOffset = self.HOLECORNERLENGTH/ 2
        return canvas.create_oval(x-cornerOffset,
                                  y-cornerOffset,
                                  x+cornerOffset,
                                  y+cornerOffset,
                                  tags=tags,
                                  outline=colour,
                                  fill='white')


class Region(Point):
    """A region of the mesh.
    Defined by a point in the region enclosed by segments. Used to tag areas.
    """

    CROSSLENGTH = 6
    TAG = 0
    MAXAREA = 1

    def __init__(self, X, Y, tag=None, maxArea=None):
        """Precondition: tag is a string and maxArea is a real"""

        self.x = X
        self.y = Y
        self.attributes = []  # index 0 is the tag string
        # optional index 1 is the max triangle area
        # NOTE the size of this attribute is assumed
        # to be 1 or 2 in regionstrings2int
        if tag is None:
            self.attributes.append("")
        else:
            self.attributes.append(tag)  # this is a string

        if maxArea is not None:
            self.setMaxArea(maxArea)  # maxArea is a number

    def getTag(self,):
        return self.attributes[self.TAG]

    def setTag(self, tag):
        self.attributes[self.TAG] = tag

    def getMaxArea(self):
        """ Returns the Max Area of a Triangle or
        None, if the Max Area has not been set.
        """
        if self.isMaxArea():
            return self.attributes[1]
        else:
            return None

    def setMaxArea(self, MaxArea):
        if MaxArea is not None:
            if self.isMaxArea():
                self.attributes[self.MAXAREA] = float(MaxArea)
            else:
                self.attributes.append(float(MaxArea))

    def deleteMaxArea(self):
        if self.isMaxArea():
            self.attributes.pop(self.MAXAREA)

    def isMaxArea(self):
        return len(self.attributes) > 1

    def draw(self, canvas, tags, scale=1, xoffset=0, yoffset=0,
             colour="black"):
        """
        Draw a black cross, returning the objectID
        """
        x = scale*(self.x + xoffset)
        y = -1*scale*(self.y + yoffset)
        cornerOffset = self.CROSSLENGTH/ 2
        return canvas.create_polygon(x,
                                     y-cornerOffset,
                                     x,
                                     y,
                                     x+cornerOffset,
                                     y,
                                     x,
                                     y,
                                     x,
                                     y+cornerOffset,
                                     x,
                                     y,
                                     x-cornerOffset,
                                     y,
                                     x,
                                     y,
                                     tags=tags,
                                     outline=colour, fill='')

    def __repr__(self):
        if self.isMaxArea():
            area = self.getMaxArea()
            return "(%f,%f,%s,%f)" % (self.x, self.y,
                                      self.getTag(), self.getMaxArea())
        else:
            return "(%f,%f,%s)" % (self.x, self.y,
                                   self.getTag())


class Segment(MeshObject):
    """
    Segments are edges whose presence in the triangulation is enforced.

    """

    def __init__(self, vertex1, vertex2, tag=None):
        """
        Each segment is specified by listing the vertices of its endpoints
        The vertices are Vertex objects
        """
        assert(vertex1 != vertex2)
        self.vertices = [vertex1, vertex2]

        if tag is None:
            self.tag = self.__class__.default
        else:
            self.tag = tag  # this is a string

    def __repr__(self):
        return "[%s,%s]" % (self.vertices, self.tag)

    def draw(self, canvas, tags, scale=1, xoffset=0,
             yoffset=0, colour=SEG_COLOUR):
        x = []
        y = []
        for end in self.vertices:
            # end.draw(canvas,scale, xoffset, yoffset ) # draw the vertices
            x.append(scale*(end.x + xoffset))
            # - since for a canvas - is up
            y.append(-1*scale*(end.y + yoffset))

        return canvas.create_line(x[0], y[0], x[1], y[1],
                                  tags=tags, fill=colour)

    def set_tag(self, tag):
        self.tag = tag

    # Class methods
    def set_default_tag(cls, default):
        cls.default = default

    def get_default_tag(cls):
        return cls.default

    set_default_tag = classmethod(set_default_tag)
    get_default_tag = classmethod(get_default_tag)


Segment.set_default_tag("")


class Rigid_triangulation(object):
    """ 
    This is a triangulation that can't have triangles added or taken away.

    It just represents the triangulation, not the mesh outline needed to
    build the triangulation.

      This is the guts of the data structure;
        generated vertex list: [(x1,y1),(x2,y2),...] (Tuples of doubles)
        generated segment list: [(point1,point2),(p3,p4),...]
            (Tuples of integers) 
        generated segment tag list: [tag,tag,...] list of strings

        generated triangle list: [(p1,p2,p3), (p4,p5,p6),....] tuple of points

        generated triangle attribute list: [s1,s2,...] list of strings

        generated triangle neighbor list: [(t1,t2,t3), (t4,t5,t6),....]
            tuple of triangles

            Should I do the basic structure like triangle, general
            mesh or something else? How about like triangle, since
            that's where the info is from, and it'll fit easier into
            this file..

             Removing the loners is difficult, since all the vert's
             after it must be removed.

             This happens in set_triangulation.

    """

    def __init__(self,
                 triangles,
                 segments,
                 vertices,
                 triangle_tags,
                 triangle_neighbors,
                 segment_tags,
                 vertex_attributes,
                 vertex_attribute_titles=None
                 ):

        self.triangles = ensure_numeric(triangles)
        self.triangle_neighbors = ensure_numeric(triangle_neighbors)
        self.triangle_tags = triangle_tags  # list of strings
        self.segments = ensure_numeric(segments)

        self.segment_tags = segment_tags  # list of strings
        self.vertices = ensure_numeric(vertices)
        # This is needed for __cmp__
        if vertex_attributes is None:
            self.vertex_attributes = []
        else:
            self.vertex_attributes = ensure_numeric(vertex_attributes)
        if vertex_attribute_titles is None:
            self.vertex_attribute_titles = []
        else:
            self.vertex_attribute_titles = vertex_attribute_titles

    def draw_triangulation(self, canvas, scale=1, xoffset=0, yoffset=0,
                           colour="green"):
        """
        Draw a triangle, returning the objectID
        """

        # FIXME(DSG-DSG) This could be a data structure that is
        # remembered and doesn't have any duplicates.
        # but I wouldn't be able to use create_polygon.
        # regard it as drawing a heap of segments

        for tri in self.triangles:
            vertices = []
            for v_index in range(3):
                vertices.append(self.vertices[tri[v_index]])

            objectID = canvas.create_polygon(
                scale*(vertices[1][0] + xoffset),
                scale*-1*(vertices[1][1] + yoffset),
                scale*(vertices[0][0] + xoffset),
                scale*-1*(vertices[0][1] + yoffset),
                scale*(vertices[2][0] + xoffset),
                scale*-1*(vertices[2][1] + yoffset),
                outline=colour, fill=''
            )

    def calc_mesh_area(self):
        area = 0
        for tri in self.triangles:
            vertices = []
            for v_index in range(3):
                vertices.append(self.vertices[tri[v_index]])
            ax = vertices[0][0]
            ay = vertices[0][1]

            bx = vertices[1][0]
            by = vertices[1][1]

            cx = vertices[2][0]
            cy = vertices[2][1]

            area += abs((bx*ay-ax*by)+(cx*by-bx*cy)+(ax*cy-cx*ay))/ 2
        return area


class Mesh(object):
    """
    Representation of a 2D triangular mesh.
    User attributes describe the mesh region/segments/vertices/attributes

    mesh attributes describe the mesh that is produced eg triangles and
    vertices.
    All point information is relative to the geo_reference passed in


    """

    def __repr__(self):
        return """
        mesh Triangles: %s 
        mesh Attribute Titles: %s 
        mesh Segments: %s  
        mesh Vertices: %s 
        user Segments: %s  
        user Vertices: %s  
        holes: %s 
        regions: %s""" % (self.meshTriangles,
                          self.attributeTitles,
                          self.meshSegments,
                          self.meshVertices,
                          self.getUserSegments(),
                          self.userVertices,
                          self.holes,
                          self.regions)

    def __init__(self,
                 userSegments=None,
                 userVertices=None,
                 holes=None,
                 regions=None,
                 geo_reference=None):
        self.meshTriangles = []
        self.attributeTitles = []
        self.meshSegments = []
        self.meshVertices = []

        # Rigid
        self.tri_mesh = None

        self.visualise_graph = True

        if userSegments is None:
            self.userSegments = []
        else:
            self.userSegments = userSegments
        self.alphaUserSegments = []

        if userVertices is None:
            self.userVertices = []
        else:
            self.userVertices = userVertices

        if holes is None:
            self.holes = []
        else:
            self.holes = holes

        if regions is None:
            self.regions = []
        else:
            self.regions = regions

        if geo_reference is None:
            self.geo_reference = Geo_reference(DEFAULT_ZONE, 0, 0)
        else:
            self.geo_reference = geo_reference

        self.shape = None

    def __cmp__(self, other):


        # FIXME(Ole): This should be refactored to use Mesh2MeshDic(self):
        # Like this:
        # self_dict = self.Mesh2MeshDic()
        # other_dict = other.Mesh2MeshDic()        
        
        # A dic for the initial m
        dic = self.Mesh2triangList()
        dic_mesh = self.Mesh2MeshList()
        for element in list(dic_mesh.keys()):
            dic[element] = dic_mesh[element]
        for element in list(dic.keys()):
            dic[element].sort()

        # A dic for the exported/imported m
        dic_other = other.Mesh2triangList()
        dic_mesh = other.Mesh2MeshList()
        for element in list(dic_mesh.keys()):
            dic_other[element] = dic_mesh[element]
        for element in list(dic.keys()):
            dic_other[element].sort()

        # print "dsg************************8"
        # print "dic ",dic
        # print "*******8"
        # print "mesh",dic_other
        # print "dic.__cmp__(dic_o)",dic.__cmp__(dic_other)
        # print "dsg************************8"

        # FIXME (Ole): For backwards compatibility - should just return True or False
        if (dic == dic_other):
           return 0
        else:
           return -1

    def __eq__(self, other):

        self_dict = self.Mesh2MeshDic()
        other_dict = other.Mesh2MeshDic()        
        
        return (self_dict == other_dict and
                self.geo_reference == other.geo_reference)
    
        
    def addUserPoint(self, pointType, x, y):
        if pointType == Vertex:
            point = self.addUserVertex(x, y)
        if pointType == Hole:
            point = self._addHole(x, y)
        if pointType == Region:
            point = self._addRegion(x, y)
        return point

    def addUserVertex(self, x, y):
        v = Vertex(x, y)
        self.userVertices.append(v)
        return v

    def _addHole(self, x, y):
        h = Hole(x, y)
        self.holes.append(h)
        return h

    def add_hole(self, x, y, geo_reference=None):
        """
        adds a point, which represents a hole.

        The point data can have it's own geo_refernece.
        If geo_reference is None the data is asumed to be absolute
        """
        [[x, y]] = self.geo_reference.change_points_geo_ref([x, y],
                                                            points_geo_ref=geo_reference)
        return self._addHole(x, y)

    def _addRegion(self, x, y):
        h = Region(x, y)
        self.regions.append(h)
        return h

    def add_region(self, x, y, geo_reference=None, tag=None):
        """
        adds a point, which represents a region.

        The point data can have it's own geo_refernece.
        If geo_reference is None the data is asumed to be absolute
        """
        # FIXME: have the user set the tag and resolution here,
        # but still return the instance, just in case.
        [[x, y]] = self.geo_reference.change_points_geo_ref([x, y],
                                                            points_geo_ref=geo_reference)
        region = self._addRegion(x, y)
        if tag is not None:
            region.setTag(tag)
        return region

    def build_grid(self,  vert_rows, vert_columns):
        """
        Build a grid with vert_rows number of vertex rows and
        vert_columns number if vertex columns

        Grid spacing of 1, the origin is the lower left hand corner.

        FIXME(DSG-DSG) no test.
        """

        for i in range(vert_rows):
            for j in range(vert_columns):
                self.addUserVertex(j, i)
        self.auto_segment()
        self.generateMesh(mode="Q", minAngle=20.0)

    def add_vertices(self, point_data):
        """
        Add user vertices.

        The point_data can be a list of (x,y) values, a numeric
        array or a geospatial_data instance.
        """
        point_data = ensure_geospatial(point_data)
        # print "point_data",point_data
        # get points relative to the mesh geo_ref
        points = point_data.get_data_points(geo_reference=self.geo_reference)

        for point in points:
            v = Vertex(point[0], point[1])
            self.userVertices.append(v)

    def add_hole_from_polygon(self, polygon, segment_tags=None,
                              geo_reference=None):
        """
        Add a polygon with tags to the current mesh, as a region.
        The maxArea of the region can be specified.

        If a geo_reference of the polygon points is given, this is used.
        If not;
        The x,y info is assumed to be Easting and Northing, absolute,
        for the meshes zone.

        polygon a list of points, in meters that describe the polygon
             (e.g. [[x1,y1],[x2,y2],...]
        tags (e.g.{'wall':[0,1,3],'ocean':[2]})

        This returns the region instance, so if the user whats to modify
        it they can.       
        """
        return self._add_area_from_polygon(polygon,
                                           segment_tags=segment_tags,
                                           hole=True,
                                           geo_reference=geo_reference
                                           )

    def add_region_from_polygon(self, polygon, segment_tags=None,
                                max_triangle_area=None, geo_reference=None,
                                region_tag=None):
        """
        Add a polygon with tags to the current mesh, as a region.
        The maxArea of the region can be specified.

        If a geo_reference of the polygon points is given, this is used.
        If not;
        The x,y info is assumed to be Easting and Northing, absolute,
        for the meshes zone.

        polygon a list of points, in meters that describe the polygon
             (e.g. [[x1,y1],[x2,y2],...]
        segment_tags (e.g.{'wall':[0,1,3],'ocean':[2]}) add tags to the
        segements of the region
        region_tags  - add a tag to all of the triangles in the region.

        This returns the region instance (if a max_triangle_area is given),
        so if the user whats to modify it they can.     
        """
        if max_triangle_area is None and region_tag is None:
            create_region = False
        else:
            create_region = True

        region = self._add_area_from_polygon(polygon,
                                             segment_tags=segment_tags,
                                             geo_reference=geo_reference,
                                             region=create_region)
        if max_triangle_area is not None:
            region.setMaxArea(max_triangle_area)
        if region_tag is not None:
            region.setTag(region_tag)

        return region

    def _add_area_from_polygon(self, polygon, segment_tags=None,
                               geo_reference=None,
                               hole=False,
                               region=False):
        """
        Add a polygon with tags to the current mesh, as a region.
        The maxArea of the region can be specified.

        If a geo_reference of the polygon points is given, this is used.
        If not;
        The x,y info is assumed to be Easting and Northing, absolute,
        for the meshes zone.

        polygon a list of points, in meters that describe the polygon
             (e.g. [[x1,y1],[x2,y2],...]
        segment_tags (e.g.{'wall':[0,1,3],'ocean':[2]})add tags to the
        segements of the region.

        This returns the region instance, so if the user whats to modify
        it they can.

        """
        # Only import this if necessary.
        # Trying to get pmesh working in an uncompiled environment
        from anuga.geometry.polygon import point_in_polygon

        # get absolute values
        if geo_reference is not None:
            polygon = geo_reference.get_absolute(polygon)
        # polygon is now absolute
        # print "polygon  should be absolute",polygon

        # create points, segs and tags
        region_dict = {}
        region_dict['points'] = polygon

        # Create segments
        # E.g. [[0,1], [1,2], [2,3], [3,0]]
        # from polygon
        # [0,1,2,3]
        segments = []
        N = len(polygon)
        for i in range(N):
            lo = i
            hi = (lo + 1) % N
            segments.append([lo, hi])
        region_dict['segments'] = segments
        region_dict['segment_tags'] = self._tag_dict2list(segment_tags, N,
                                                          hole=hole)

        self.addVertsSegs(region_dict)  # this is passing absolute values

        if region is True:
            # get inner point - absolute values
            inner_point = point_in_polygon(polygon)
            inner = self.add_region(inner_point[0], inner_point[1],
                                    geo_reference=None)
        elif hole is True:
            # get inner point - absolute values
            inner_point = point_in_polygon(polygon)
            inner = self.add_hole(inner_point[0], inner_point[1],
                                  geo_reference=None)
        else:
            inner = None

        return inner

    def _tag_dict2list(self, tags, number_of_segs, hole=False):
        """
        Convert a tag dictionary from this sort of format;
        #{'wall':[0,3],'ocean':[2]}

        To a list format
        # ['wall', '', 'ocean', 'wall']

        Note: '' is a default value of nothing
        """
        # FIXME (DSG-DSG): Using '' as a default isn't good.
        # Try None.
        # Due to this default this method is too connected to
        # _add_area_from_polygon

        if hole:
            default_tag = 'interior'
        else:
            default_tag = ''
        segment_tags = [default_tag]*number_of_segs
        if tags is not None:
            for key in tags:
                indices = tags[key]
                for i in indices:
                    segment_tags[i] = key
        return segment_tags

    def add_circle(self, center, radius, segment_count=100,
                   center_geo_reference=None, tag=None,
                   region=False, hole=False):
        """
        center is a point, in absulute or relative to center_geo_ref
        radius is the radius of the circle
        segment_count is the number of segments used to represent the circle
        tag is a string name, given to the segments.
        If region is True a region object is returned.
        If hole is True a hole object is returned.
           (Don't have them both True.)


        """
        # convert center and radius to a polygon
        cuts = []
        factor = 2 * math.pi/ segment_count
        for cut in range(segment_count):
            cuts.append(cut*factor)

        polygon = []
        for cut in cuts:

            x = center[0] + radius * math.cos(cut)
            y = center[1] + radius * math.sin(cut)
            polygon.append([x, y])
        # build the tags
        tags = {tag: list(range(segment_count))}

        return self._add_area_from_polygon(polygon, segment_tags=tags,
                                           region=region, hole=hole,
                                           geo_reference=center_geo_reference)

    def auto_set_geo_reference(self):
        """
        Automatically set the georeference, based in the minimum x and y
        user vertex values.

        Not to be used with the graphical interface

        Not implemented.
        Don't implement now.  using the new georeferenced points class
        will change this?
        """
        # to do
        # find the lower left hand corner
        [xmin, ymin, xmax, ymax] = self.boxsize()

        # set all points to that lower left hand corner.
        # use change_geo_reference
        new_geo = Geo_reference(self.geo_reference.get_zone(), xmin, ymin)
        self.change_geo_reference(new_geo)

    def change_geo_reference(self, new_geo_reference):
        """
        Change from one geo_reference to another.
        Not to be used with the graphical interface
        """
        # FIXME
        # change georeference of;
        # self.userVertices = self.geo_reference.change_points_geo_ref( \
        # self.userVertices)
        #self.holes = self.geo_reference.change_points_geo_ref(self.holes)
        #self.regions = self.geo_reference.change_points_geo_ref(self.regions)
        # The above will not work.
        # since userVertices (etc) is a list of point objects,
        # not a list of lists.
        # add a method to the points class to fix this up.

    def add_segment(self, v1, v2, tag):
        """
        Don't do this function.
        what will the v's be objects?  How is the user suppost to get them?

        Indexes?  If add vertstosegs or add_region_from_polygon is called
        more than once then the actual index is not obvious.  Making this
        function confusing.
        """
        pass

    def add_points_and_segments(self, points,
                                segments=None, segment_tags=None):
        """
        Add an outline of the mesh.
        Points is a list of points a standard representation of points.
        Segments is a list of tuples of integers.  Each tuple defines the
           start and end of the segment by it's point index.
        segment_tags is an optional dictionary which is used to add tags to
           the segments.  The key is the tag name, value is the list of segment
           indexes the tag will apply to.
           eg. {'wall':[0,3],'ocean':[2]}

        """
        # make sure the points are absolute
        # Since addVertsSegs will deal with georeferencing.
        points = ensure_absolute(points)

        if segments is None:
            segments = []
            for i in range(len(points)-1):
                segments.append([i, i+1])

        # create points, segs and tags
        region_dict = {}
        region_dict['points'] = points
        region_dict['segments'] = segments
        region_dict['segment_tags'] = self._tag_dict2list(segment_tags,
                                                          len(segments))
        self.addVertsSegs(region_dict)

    def addVertsSegs(self, outlineDict):
        """
        Add out-line (user Mesh) attributes given a dictionary of the lists
        points: [(x1,y1),(x2,y2),...] (Tuples of doubles)  
        segments: [(point1,point2),(p3,p4),...] (Tuples of integers)
        segment_tags: [S1Tag, S2Tag, ...] (list of strings)

        Assume the values are in Eastings and Northings, with no reference
        point. eg absolute
        """
        if 'segment_tags' not in outlineDict:
            outlineDict['segment_tags'] = []
            for i in range(len(outlineDict['segments'])):
                outlineDict['segment_tags'].append('')
        # print "outlineDict['segment_tags']",outlineDict['segment_tags']
        # print "outlineDict['points']",outlineDict['points']
        # print "outlineDict['segments']",outlineDict['segments']

        i_offset = len(self.userVertices)
        # print "self.userVertices",self.userVertices
        # print "index_offset",index_offset
        for point in outlineDict['points']:
            v = Vertex(point[0]-self.geo_reference.xllcorner,
                       point[1]-self.geo_reference.yllcorner)
            self.userVertices.append(v)

        for seg, seg_tag in zip(outlineDict['segments'],
                                outlineDict['segment_tags']):
            segObject = Segment(self.userVertices[int(seg[0])+i_offset],
                                self.userVertices[int(seg[1])+i_offset])
            if not seg_tag == '':
                segObject.set_tag(seg_tag)
            self.userSegments.append(segObject)

    def get_triangle_count(self):
        return len(self.getTriangulation())

    def getUserVertices(self):
        """
        Note: The x,y values will be relative to the mesh geo_reference
        This returns a list of vertex objects
        """
        return self.userVertices

    def get_user_vertices(self, absolute=True):
        """
        This returns a list of (x,y) values
        (maybe it should return a geospatical object?
        It makes it a bit confusing though.)
        """
        pointlist = []
        for vertex in self.userVertices:
            pointlist.append([vertex.x, vertex.y])
        spat = Geospatial_data(pointlist, geo_reference=self.geo_reference)
        return spat.get_data_points(absolute=absolute)

    def getUserSegments(self):
        allSegments = self.userSegments + self.alphaUserSegments
        # print "self.userSegments",self.userSegments
        # print "self.alphaUserSegments",self.alphaUserSegments
        # print "allSegments",allSegments
        return allSegments

    def deleteUserSegments(self, seg):
        if self.userSegments.count(seg) == 0:
            self.alphaUserSegments.remove(seg)
            pass
        else:
            self.userSegments.remove(seg)

    def clearUserSegments(self):
        self.userSegments = []
        self.alphaUserSegments = []

    # FIXME see where this is used. return an array instead
    def getTriangulation(self):
        # return self.meshTriangles
        return self.tri_mesh.triangles.tolist()

    def getMeshVertices(self):
        # return self.meshVertices
        return self.tri_mesh.vertices

    def getMeshVerticeAttributes(self):
        # return self.meshVertices
        return self.tri_mesh.vertex_attributes

    def getMeshSegments(self):
        # return self.meshSegments
        return self.tri_mesh.segments

    def getMeshSegmentTags(self):
        # return self.meshSegments
        return self.tri_mesh.segment_tags

    def getHoles(self):
        return self.holes

    def getRegions(self):
        return self.regions

    def isTriangulation(self):
        if self.meshVertices == []:
            return False
        else:
            return True

    def addUserSegment(self, v1, v2):
        """
        PRECON: A segment between the two vertices is not already present.
        Check by calling isUserSegmentNew before calling this function.

        """
        s = Segment(v1, v2)
        self.userSegments.append(s)
        return s

    def generate_mesh(self,
                      maximum_triangle_area="",
                      minimum_triangle_angle=28.0,
                      verbose=False):
        if verbose is True:
            silent = ''
        else:
            silent = 'Q'
        self.generateMesh(mode=silent + "pzq"+str(minimum_triangle_angle)
                          + "a"+str(maximum_triangle_area)
                          + "a", verbose=verbose)
        # The last a is so areas for regions will be used

    def generateMesh(self, mode=None, maxArea=None, minAngle=None,
                     isRegionalMaxAreas=True, verbose=False):
        """
        Based on the current user vaules, holes and regions
        generate a new mesh
        mode is a string that sets conditions on the mesh generations
        see triangle_instructions.txt for a definition of the commands

        PreCondition: maxArea is a double between 1e-20 and 1e30 or is a
        string.
        """
        # print "mode ",mode
        if mode is None:
            self.mode = ""
        else:
            self.mode = mode

        if self.mode.find('p') < 0:
            self.mode += 'p'  # p - read a planar straight line graph.
            # there must be segments to use this switch
            # TODO throw an aception if there aren't seg's
            # it's more comlex than this.  eg holes
        if self.mode.find('z') < 0:
            self.mode += 'z'  # z - Number all items starting from zero
            # (rather than one)
        if self.mode.find('n'):
            self.mode += 'n'  # n - output a list of neighboring triangles
        if self.mode.find('A') < 0:
            self.mode += 'A'  # A - output region attribute list for triangles

        if not self.mode.find('V') and not self.mode.find('Q'):
            self.mode += 'V'  # V - output info about what Triangle is doing

        if self.mode.find('q') < 0 and minAngle is not None:
            #   print "**********8minAngle******** ",minAngle
            min_angle = 'q' + str(minAngle)
            self.mode += min_angle  # z - Number all items starting from zero
            # (rather than one)
        if maxArea != None:
            self.mode += 'a' + str(maxArea)
            try:
                self.mode += 'a' + '%20.20f' % maxArea
            except TypeError:
                self.mode += 'a' + str(maxArea)
            # print "self.mode", self.mode
        # FIXME (DSG-DSG) This isn't explained.
        if isRegionalMaxAreas:
            self.mode += 'a'
        # print "mesh#generateMesh# self.mode",self.mode
        meshDict = self.Mesh2triangList()

        # FIXME (DSG-DSG)  move below section into generate_mesh.py
        #                  & 4 functions eg segment_strings2ints
        # Actually, because of region_list.append((1.0,2.0,""))
        # don't move it, without careful thought
        # print "*************************!@!@ This is going to triangle   !@!@"
        # print meshDict
        # print "************************!@!@ This is going to triangle   !@!@"

        # print "meshDict['segmenttaglist']", meshDict['segmenttaglist']
        [meshDict['segmenttaglist'],
         segconverter] = segment_strings2ints(meshDict['segmenttaglist'],
                                              initialconversions)
        # print "regionlist",meshDict['regionlist']
        [meshDict['regionlist'],
         regionconverter] = region_strings2ints(meshDict['regionlist'])
        # print "%%%%%%%%%%%%%%%%%%%%%%%%%%%regionlist",meshDict['regionlist']
        # print "meshDict['segmenttaglist']", meshDict['segmenttaglist'
        # print "self.mode", self.mode
        generatedMesh = generate_mesh(
            meshDict['pointlist'],
            meshDict['segmentlist'],
            meshDict['holelist'],
            meshDict['regionlist'],
            meshDict['pointattributelist'],
            meshDict['segmenttaglist'],
            self.mode,
            meshDict['pointlist'],
            verbose)
        # print "%%%%%%%%%%%%%%%%%%%%%%%%%%%generated",generatedMesh
        generatedMesh['qaa'] = 1
        if generatedMesh['generatedsegmentmarkerlist'] is not None:
            generatedMesh['generatedsegmentmarkerlist'] = \
                segment_ints2strings(generatedMesh['generatedsegmentmarkerlist'],
                                     segconverter)
        # print "processed gen",generatedMesh['generatedsegmentmarkerlist']
        # print "pmesh mesh generatedMesh['generatedtriangleattributelist']", generatedMesh['generatedtriangleattributelist']
        if generatedMesh['generatedtriangleattributelist'] is not None:
            generatedMesh['generatedtriangleattributelist'] = \
                region_ints2strings(generatedMesh['generatedtriangleattributelist'],
                                    regionconverter)

        # print "pmesh mesh generatedMesh['generatedtriangleattributelist']", generatedMesh['generatedtriangleattributelist']
        # FIXME (DSG-DSG)  move above section into generate_mesh.py

        if generatedMesh['generatedpointattributelist'] is None or \
                generatedMesh['generatedpointattributelist'].shape[1] == 0:
            self.attributeTitles = []
        generatedMesh['generatedpointattributetitlelist'] = \
            self.attributeTitles
        # print "################  FROM TRIANGLE"
        # print "generatedMesh",generatedMesh
        # print "##################"
        self.setTriangulation(generatedMesh)

    def clearTriangulation(self):

        # Clear the current generated mesh values
        self.meshTriangles = []
        self.meshSegments = []
        self.meshVertices = []

    def removeDuplicatedUserVertices(self):
        """Pre-condition: There are no user segments
        This function will keep the first duplicate
        """
        assert self.getUserSegments() == []
        self.userVertices, counter = self.removeDuplicatedVertices(
            self.userVertices)
        return counter

    def removeDuplicatedVertices(self, Vertices):
        """
        This function will keep the first duplicate, remove all others
        Precondition: Each vertex has a dupindex, which is the list
        index.

        Note: this removes vertices that have the same x,y values,
        not duplicate instances in the Vertices list.
        """
        remove = []
        index = 0
        for v in Vertices:
            v.dupindex = index
            index += 1
        t = list(Vertices)

        # Using https://docs.python.org/3/howto/sorting.html
        t.sort(key=cmp_to_key(Point.cmp_xy)) # For Python 3.x
        #t.sort(Point.cmp_xy)  # For Python2.7

        length = len(t)
        behind = 0
        ahead = 1
        counter = 0
        while ahead < length:
            b = t[behind]
            ah = t[ahead]
            if (b.y == ah.y and b.x == ah.x):
                remove.append(ah.dupindex)
            behind += 1
            ahead += 1

        # remove the duplicate vertices
        remove.sort()
        remove.reverse()
        for i in remove:
            Vertices.pop(i)
            pass

        # Remove the attribute that this function added
        for v in Vertices:
            del v.dupindex
        return Vertices, counter

    # FIXME (DSG-DSG) Move this to geospatial
    def thinoutVertices(self, delta):
        """Pre-condition: There are no user segments
        This function will keep the first duplicate
        """
        assert self.getUserSegments() == []
        #t = self.userVertices
        #self.userVertices =[]
        boxedVertices = {}
        thinnedUserVertices = []
        delta = round(delta, 1)

        for v in self.userVertices:
            # tag is the center of the boxes
            tag = (round(v.x/delta, 0)*delta,
                   round(v.y/ delta, 0)*delta)
            # this creates a dict of lists of faces, indexed by tag
            boxedVertices.setdefault(tag, []).append(v)

        for [tag, verts] in list(boxedVertices.items()):
            min = delta
            bestVert = None
            tagVert = Vertex(tag[0], tag[1])
            for v in verts:
                dist = v.DistanceToPoint(tagVert)
                if (dist < min):
                    min = dist
                    bestVert = v
            thinnedUserVertices.append(bestVert)
        self.userVertices = thinnedUserVertices

    def isUserSegmentNew(self, v1, v2):
        identicalSegs = [x for x in self.getUserSegments()
                         if (x.vertices[0] == v1 and x.vertices[1] == v2)
                         or (x.vertices[0] == v2 and x.vertices[1] == v1)]

        return len(identicalSegs) == 0

    def deleteSegsOfVertex(self, delVertex):
        """
        Delete this vertex and any segments that connect to it.
        """
        # Find segments that connect to delVertex
        deletedSegments = []
        for seg in self.getUserSegments():
            if (delVertex in seg.vertices):
                deletedSegments.append(seg)
        # Delete segments that connect to delVertex
        for seg in deletedSegments:
            self.deleteUserSegments(seg)
        return deletedSegments

    def deleteMeshObject(self, MeshObject):
        """
        Returns a list of all objects that were removed
        """
        deletedObs = []
        if isinstance(MeshObject, Vertex):
            deletedObs = self.deleteSegsOfVertex(MeshObject)
            deletedObs.append(MeshObject)
            self.userVertices.remove(MeshObject)
        elif isinstance(MeshObject, Segment):
            deletedObs.append(MeshObject)
            self.deleteUserSegments(MeshObject)
        elif isinstance(MeshObject, Hole):
            deletedObs.append(MeshObject)
            self.holes.remove(MeshObject)
        elif isinstance(MeshObject, Region):
            deletedObs.append(MeshObject)
            self.regions.remove(MeshObject)
        return deletedObs

    def Mesh2triangList(self, userVertices=None,
                        userSegments=None,
                        holes=None,
                        regions=None):
        """
        Convert the Mesh to a dictionary of the lists needed for the
        triang module
        points list: [(x1,y1),(x2,y2),...] (Tuples of doubles)
        pointattributelist: [(a11,a12,...),(a21,a22),...] (Tuples of doubles)
        segment list: [(point1,point2),(p3,p4),...] (Tuples of integers) 
        hole list: [(x1,y1),...](Tuples of doubles, one inside each hole)
        regionlist: [ (x1,y1,tag, max area),...] (Tuple of 3-4 doubles)

        Note, this adds an index attribute to the user Vertex objects.

        Used to produce output to triangle
        """
        if userVertices is None:
            userVertices = self.getUserVertices()
        if userSegments is None:
            userSegments = self.getUserSegments()
        if holes is None:
            holes = self.getHoles()
        if regions is None:
            regions = self.getRegions()

        meshDict = {}

        pointlist = []
        pointattributelist = []
        index = 0
        for vertex in userVertices:
            vertex.index = index
            pointlist.append((vertex.x, vertex.y))
            pointattributelist.append((vertex.attributes))

            index += 1
        meshDict['pointlist'] = pointlist
        meshDict['pointattributelist'] = pointattributelist

        segmentlist = []
        segmenttaglist = []
        for seg in userSegments:
            segmentlist.append((seg.vertices[0].index, seg.vertices[1].index))
            segmenttaglist.append(seg.tag)
        meshDict['segmentlist'] = segmentlist
        meshDict['segmenttaglist'] = segmenttaglist

        holelist = []
        for hole in holes:
            holelist.append((hole.x, hole.y))
        meshDict['holelist'] = holelist

        regionlist = []
        for region in regions:
            if (region.getMaxArea() != None):
                regionlist.append((region.x, region.y, region.getTag(),
                                   region.getMaxArea()))
            else:
                regionlist.append((region.x, region.y, region.getTag()))
        meshDict['regionlist'] = regionlist
        # print "*(*("
        # print meshDict
        # print meshDict['regionlist']
        # print "*(*("
        return meshDict

    def Mesh2MeshList(self):
        """
        Convert the Mesh to a dictionary of lists describing the
        triangulation variables;

        This is used by __cmp__
        generated point list: [(x1,y1),(x2,y2),...] (Tuples of doubles)
        generated point attribute list: [(a11,a12,...),(a21,a22),...]
            (Tuples of doubles)
        generated point attribute title list:[A1Title, A2Title ...]
            (list of strings)
        generated segment list: [(point1,point2),(p3,p4),...]
            (Tuples of integers) 
        generated segment tag list: [tag,tag,...] list of strings

        generated triangle list: [(p1,p2,p3), (p4,p5,p6),....] tuple of points

        generated triangle attribute list: [s1,s2,...] list of strings

        generated triangle neighbor list: [(t1,t2,t3), (t4,t5,t6),....]
            tuple of triangles

        Used to produce .tsh file
        """

        meshDict = {}
        # print "old meshDict['generatedpointattributetitlelist']",meshDict['generatedpointattributetitlelist']
        # print "self.tri_mesh", self.tri_mesh
        if self.tri_mesh is not None:
            # print "self.tri_mesh.triangles.tolist()", self.tri_mesh.triangles.tolist()
            meshDict['generatedtrianglelist'] = self.tri_mesh.triangles.tolist()
            meshDict['generatedtriangleattributelist'] = self.tri_mesh.triangle_tags
            meshDict['generatedtriangleneighborlist'] = self.tri_mesh.triangle_neighbors.tolist()
            meshDict['generatedpointlist'] = self.tri_mesh.vertices.tolist()
            if self.tri_mesh.vertex_attributes == []:
                meshDict['generatedpointattributelist'] = []
            #meshDict['generatedpointattributelist'] = self.tri_mesh.vertex_attributes
            meshDict['generatedpointattributetitlelist'] = \
                self.tri_mesh.vertex_attribute_titles
            meshDict['generatedsegmentlist'] = self.tri_mesh.segments.tolist()
            meshDict['generatedsegmenttaglist'] = self.tri_mesh.segment_tags
        else:
            meshDict['generatedtrianglelist'] = []
            meshDict['generatedtriangleattributelist'] = []
            meshDict['generatedtriangleneighborlist'] = []
            meshDict['generatedpointlist'] = []
            meshDict['generatedpointattributelist'] = []
            meshDict['generatedpointattributetitlelist'] = []
            meshDict['generatedsegmentlist'] = []
            meshDict['generatedsegmenttaglist'] = []

        # print "new meshDict['generatedpointattributetitlelist']",meshDict['generatedpointattributetitlelist']
        # print "mesh.Mesh2MeshList*)*)"
        # print meshDict
        # print "mesh.Mesh2MeshList*)*)"

        return meshDict

    def Mesh2MeshDic(self):
        """
        Convert the user and generated info of a mesh to a dictionary
        structure
        """
        dic = self.Mesh2triangList()
        dic_mesh = self.Mesh2MeshList()
        for element in list(dic_mesh.keys()):
            dic[element] = dic_mesh[element]
        return dic

    def setTriangulation(self, genDict):
        """
        Set the mesh attributes given a dictionary of the lists
        returned from the triang module       
        generated point list: [(x1,y1),(x2,y2),...] (Tuples of doubles)  
        generated point attribute list:[(P1att1,P1attt2, ...),
            (P2att1,P2attt2,...),...]- not implemented
        generated point attribute title list:[A1Title, A2Title ...]
            (list of strings) - not implemented
        generated segment list: [(point1,point2),(p3,p4),...]
            (Tuples of integers)
        generated segment marker list: [S1Tag, S2Tag, ...] (list of strings)
        triangle list:  [(point1,point2, point3),(p5,p4, p1),...]
            (Tuples of integers)
        triangle neighbor list: [(triangle1,triangle2, triangle3),
            (t5,t4, t1),...] (Tuples of integers)
            -1 means there's no triangle neighbor
        triangle attribute list: [(T1att), (T2att), ...](list of strings)

        """
        # Setting up the rigid triangulation
        self.tri_mesh = Rigid_triangulation(
            genDict['generatedtrianglelist'], genDict['generatedsegmentlist'], genDict['generatedpointlist'], genDict['generatedtriangleattributelist'], genDict[
                'generatedtriangleneighborlist'], genDict['generatedsegmentmarkerlist'], genDict['generatedpointattributelist'], genDict['generatedpointattributetitlelist']
        )

    def setMesh(self, genDict):
        """
        Set the user Mesh attributes given a dictionary of the lists
        point list: [(x1,y1),(x2,y2),...] (Tuples of doubles)  
        point attribute list:[(P1att1,P1attt2, ...),(P2att1,P2attt2,...),...]
        segment list: [(point1,point2),(p3,p4),...] (Tuples of integers)
        segment tag list: [S1Tag, S2Tag, ...] (list of ints)
        region list: [(x1,y1),(x2,y2),...] (Tuples of doubles)
        region attribute list: ["","reservoir",""] list of strings
        region max area list:[real, None, Real,...] list of None and reals

        mesh is an instance of a mesh object
        """
        # Clear the current user mesh values
        self.clearUserSegments()
        self.userVertices = []
        self.Holes = []
        self.Regions = []

        # print "mesh.setMesh@#@#@#"
        # print genDict
        # print "@#@#@#"

        #index = 0
        for point in genDict['pointlist']:
            v = Vertex(point[0], point[1])
            #v.index = index
            #index +=1
            self.userVertices.append(v)

        #index = 0
        for seg, tag in zip(genDict['segmentlist'],
                            genDict['segmenttaglist']):
            segObject = Segment(self.userVertices[int(seg[0])],
                                self.userVertices[int(seg[1])], tag=tag)
            #segObject.index = index
            #index +=1
            self.userSegments.append(segObject)

# Remove the loading of attribute info.
# Have attribute info added using least_squares in pyvolution
#         index = 0
#         for att in genDict['pointattributelist']:
#             if att is None:
#                 self.userVertices[index].setAttributes([])
#             else:
#                 self.userVertices[index].setAttributes(att)
#            index += 1

        #index = 0
        for point in genDict['holelist']:
            h = Hole(point[0], point[1])
            #h.index = index
            #index +=1
            self.holes.append(h)

        #index = 0
        for reg, att, maxArea in zip(
                genDict['regionlist'],
                genDict['regionattributelist'],
                genDict['regionmaxarealist']):
            Object = Region(reg[0],
                            reg[1],
                            tag=att,
                            maxArea=maxArea)
            #Object.index = index
            #index +=1
            self.regions.append(Object)

    def Testauto_segment(self):
        newsegs = []
        s1 = Segment(self.userVertices[0],
                     self.userVertices[1])
        s2 = Segment(self.userVertices[0],
                     self.userVertices[2])
        s3 = Segment(self.userVertices[2],
                     self.userVertices[1])
        if self.isUserSegmentNew(s1.vertices[0], s1.vertices[1]):
            newsegs.append(s1)
        if self.isUserSegmentNew(s2.vertices[0], s2.vertices[1]):
            newsegs.append(s2)
        if self.isUserSegmentNew(s3.vertices[0], s3.vertices[1]):
            newsegs.append(s3)
        # DSG!!!
        self.userSegments.extend(newsegs)
        return newsegs

    def savePickle(self, currentName):
        fd = open(currentName, 'w')
        pickle.dump(self, fd)
        fd.close()

    def auto_segmentFilter(self, raw_boundary=True,
                           remove_holes=False,
                           smooth_indents=False,
                           expand_pinch=False):
        """
        Change the filters applied on the alpha shape boundary
        """
        if self.shape is None:
            return [], [], 0.0
        return self._boundary2mesh(raw_boundary=raw_boundary,
                                   remove_holes=remove_holes,
                                   smooth_indents=smooth_indents,
                                   expand_pinch=expand_pinch)

    def auto_segment(self, alpha=None,
                     raw_boundary=True,
                     remove_holes=False,
                     smooth_indents=False,
                     expand_pinch=False):
        """
        Precon: There must be 3 or more vertices in the userVertices structure

        This returns alpha_segs_no_user_segs, segs2delete, optimum_alpha
        """
        self._createBoundary(alpha=alpha)
        return self._boundary2mesh(raw_boundary=raw_boundary,
                                   remove_holes=remove_holes,
                                   smooth_indents=smooth_indents,
                                   expand_pinch=expand_pinch)

    def _createBoundary(self, alpha=None):
        """
        """
        points = []
        for vertex in self.getUserVertices():
            points.append((vertex.x, vertex.y))
        self.shape = anuga.alpha_shape.alpha_shape.Alpha_Shape(points,
                                                               alpha=alpha)

    def _boundary2mesh(self, raw_boundary=True,
                       remove_holes=False,
                       smooth_indents=False,
                       expand_pinch=False):
        """
        Precon there must be a shape object.
        """
        self.shape.set_boundary_type(raw_boundary=raw_boundary,
                                     remove_holes=remove_holes,
                                     smooth_indents=smooth_indents,
                                     expand_pinch=expand_pinch)
        boundary_segs = self.shape.get_boundary()
        # print "boundary_segs",boundary_segs
        segs2delete = self.alphaUserSegments
        # FIXME(DSG-DSG) this algorithm needs comments
        new_segs = {}
        #alpha_segs = []
        #user_segs = []
        for seg in boundary_segs:
            v1 = self.userVertices[int(seg[0])]
            v2 = self.userVertices[int(seg[1])]
            boundary_seg = Segment(v1, v2)
            new_segs[(v1, v2)] = boundary_seg

        for user_seg in self.userSegments:
            if (user_seg.vertices[0],
                    user_seg.vertices[1]) in new_segs:
                del new_segs[user_seg.vertices[0],
                             user_seg.vertices[1]]
            elif (user_seg.vertices[1],
                  user_seg.vertices[0]) in new_segs:
                del new_segs[user_seg.vertices[1],
                             user_seg.vertices[0]]

        optimum_alpha = self.shape.get_alpha()
        alpha_segs_no_user_segs = list(new_segs.values())
        self.alphaUserSegments = alpha_segs_no_user_segs
        return alpha_segs_no_user_segs, segs2delete, optimum_alpha

    def representedAlphaUserSegment(self, v1, v2):
        identicalSegs = [x for x in self.alphaUserSegments
                         if (x.vertices[0] == v1 and x.vertices[1] == v2)
                         or (x.vertices[0] == v2 and x.vertices[1] == v1)]

        if identicalSegs == []:
            return None
        else:
            # Only return the first one.
            return identicalSegs[0]

    def representedUserSegment(self, v1, v2):
        identicalSegs = [x for x in self.userSegments
                         if (x.vertices[0] == v1 and x.vertices[1] == v2)
                         or (x.vertices[0] == v2 and x.vertices[1] == v1)]

        if identicalSegs == []:
            return None
        else:
            # Only return the first one.
            return identicalSegs[0]

    def joinVertices(self):
        """
        Return list of segments connecting all userVertices
        in the order they were given

        Precon: There must be 3 or more vertices in the userVertices structure
        """

        newsegs = []

        v1 = self.userVertices[0]
        for v2 in self.userVertices[1:]:
            if self.isUserSegmentNew(v1, v2):
                newseg = Segment(v1, v2)
                newsegs.append(newseg)
            v1 = v2

        # Connect last point to the first
        v2 = self.userVertices[0]
        if self.isUserSegmentNew(v1, v2):
            newseg = Segment(v1, v2)
            newsegs.append(newseg)

        # Update list of user segments
        # DSG!!!
        self.userSegments.extend(newsegs)
        return newsegs

    def normaliseMesh(self, scale, offset, height_scale):
        [xmin, ymin, xmax, ymax] = self.boxsize()
        [attmin0, attmax0] = self.maxMinVertAtt(0)
        # print "[attmin0, attmax0]" ,[attmin0, attmax0]
        [attmin1, attmax1] = self.maxMinVertAtt(1)
        #print [xmin, ymin, xmax, ymax]
        xrange = xmax - xmin
        yrange = ymax - ymin
        if xrange > yrange:
            min, max = xmin, xmax
        else:
            min, max = ymin, ymax

        for obj in self.getUserVertices():
            obj.x = (obj.x - xmin)/ (max - min)*scale + offset
            obj.y = (obj.y - ymin)/ (max - min)*scale + offset
            if len(obj.attributes) > 0 and attmin0 != attmax0:
                obj.attributes[0] = (obj.attributes[0]-attmin0)/(attmax0-attmin0)*height_scale
            if len(obj.attributes) > 1 and attmin1 != attmax1:
                obj.attributes[1] = (obj.attributes[1]-attmin1)/(attmax1-attmin1)*height_scale

        for obj in self.getMeshVertices():
            obj.x = (obj.x - xmin)/ (max - min)*scale + offset
            obj.y = (obj.y - ymin)/ (max - min)*scale + offset
            if len(obj.attributes) > 0 and attmin0 != attmax0:
                obj.attributes[0] = (obj.attributes[0]-attmin0)/(attmax0-attmin0)*height_scale
            if len(obj.attributes) > 1 and attmin1 != attmax1:
                obj.attributes[1] = (obj.attributes[1]-attmin1)/(attmax1-attmin1)*height_scale

        for obj in self.getHoles():
            obj.x = (obj.x - xmin)/ (max - min)*scale + offset
            obj.y = (obj.y - ymin)/ (max - min)*scale + offset
        for obj in self.getRegions():
            obj.x = (obj.x - xmin)/ (max - min)*scale + offset
            obj.y = (obj.y - ymin)/ (max - min)*scale + offset
        [xmin, ymin, xmax, ymax] = self.boxsize()
        #print [xmin, ymin, xmax, ymax]

    def boxsizeVerts(self):
        """
        Returns a list of verts denoting a box or triangle that contains
        verts on the xmin, ymin, xmax and ymax axis.
        Structure: list of verts 
        """

        large = kinds.default_float_kind.MAX
        xmin = large
        xmax = -large
        ymin = large
        ymax = -large
        for vertex in self.userVertices:
            if vertex.x < xmin:
                xmin = vertex.x
                xminVert = vertex
            if vertex.x > xmax:
                xmax = vertex.x
                xmaxVert = vertex

            if vertex.y < ymin:
                ymin = vertex.y
                yminVert = vertex
            if vertex.y > ymax:
                ymax = vertex.y
                ymaxVert = vertex
        verts, count = self.removeDuplicatedVertices([xminVert,
                                                      xmaxVert,
                                                      yminVert,
                                                      ymaxVert])

        return verts

    def boxsize(self):
        """
        Returns a list denoting a box that contains the entire
        structure of vertices
        Structure: [xmin, ymin, xmax, ymax] 
        """

        large = kinds.default_float_kind.MAX
        xmin = large
        xmax = -large
        ymin = large
        ymax = -large
        for vertex in self.userVertices:
            if vertex.x < xmin:
                xmin = vertex.x
            if vertex.x > xmax:
                xmax = vertex.x

            if vertex.y < ymin:
                ymin = vertex.y
            if vertex.y > ymax:
                ymax = vertex.y
        return [xmin, ymin, xmax, ymax]

    def maxMinVertAtt(self, iatt):
        """
        Returns a list denoting a box that contains the entire structure
        of vertices
        Structure: [xmin, ymin, xmax, ymax] 
        """

        large = kinds.default_float_kind.MAX
        min = large
        max = -large
        for vertex in self.userVertices:
            if len(vertex.attributes) > iatt:
                if vertex.attributes[iatt] < min:
                    min = vertex.attributes[iatt]
                if vertex.attributes[iatt] > max:
                    max = vertex.attributes[iatt]
        for vertex in self.meshVertices:
            if len(vertex.attributes) > iatt:
                if vertex.attributes[iatt] < min:
                    min = vertex.attributes[iatt]
                if vertex.attributes[iatt] > max:
                    max = vertex.attributes[iatt]
        return [min, max]

    def scaleoffset(self, WIDTH, HEIGHT):
        """
        Returns a list denoting the scale and offset terms that need to be
        applied when converting  mesh co-ordinates onto grid co-ordinates 
        Structure: [scale, xoffset, yoffset] 
        """
        OFFSET = 0.05*min([WIDTH, HEIGHT])
        [xmin, ymin, xmax, ymax] = self.boxsize()
        SCALE = min([0.9*WIDTH, 0.9*HEIGHT])/max([xmax-xmin, ymax-ymin])

        if SCALE*xmin < OFFSET:
            xoffset = abs(SCALE*xmin) + OFFSET
        if SCALE*xmax > WIDTH - OFFSET:
            xoffset = -(SCALE*xmax - WIDTH + OFFSET)
        if SCALE*ymin < OFFSET:
            b = abs(SCALE*ymin)+OFFSET
        if SCALE*ymax > HEIGHT-OFFSET:
            b = -(SCALE*ymax - HEIGHT + OFFSET)
        yoffset = HEIGHT - b
        return [SCALE, xoffset, yoffset]

    def exportASCIIobj(self, ofile):
        """
        export a file, ofile, with the format
         lines:  v <x> <y> <first attribute>
        f <vertex #>  <vertex #> <vertex #> (of the triangles)
        """
        fd = open(ofile, 'w')
        self.writeASCIIobj(fd)
        fd.close()

    def writeASCIIobj(self, fd):
        fd.write(" # Triangulation as an obj file\n")
        numVert = str(len(self.meshVertices))

        index1 = 1
        for vert in self.meshVertices:
            vert.index1 = index1
            index1 += 1

            fd.write("v "
                     + str(vert.x) + " "
                     + str(vert.y) + " "
                     + str(vert.attributes[0]) + "\n")

        for tri in self.meshTriangles:
            fd.write("f "
                     + str(tri.vertices[0].index1) + " "
                     + str(tri.vertices[1].index1) + " "
                     + str(tri.vertices[2].index1) + "\n")

    def exportASCIIsegmentoutlinefile(self, ofile):
        """Write the boundary user mesh info, eg
         vertices that are connected to segments,
         segments
        """

        verts = {}
        for seg in self.getUserSegments():
            verts[seg.vertices[0]] = seg.vertices[0]
            verts[seg.vertices[1]] = seg.vertices[1]
        meshDict = self.Mesh2IOOutlineDict(userVertices=list(verts.values()))
        export_mesh_file(ofile, meshDict)

        # exportASCIImeshfile   - this function is used
    def export_mesh_file(self, ofile):
        """
        export a file, ofile, with the format
        """

        dict = self.Mesh2IODict()
        export_mesh_file(ofile, dict)

    # FIXME(DSG-DSG):Break this into two functions.
    # One for the outline points.
    # One for the mesh points.
    # Note: this function is not in the gui
    def exportPointsFile(self, ofile):
        """
        export a points file, ofile.

        """

        mesh_dict = self.Mesh2IODict()
        #point_dict = {}
        # point_dict['attributelist'] = {} #this will need to be expanded..
        # if attributes are brought back in.
        #point_dict['geo_reference'] = self.geo_reference
        if len(mesh_dict['vertices']) == 0:
            #point_dict['pointlist'] = mesh_dict['points']
            geo = Geospatial_data(mesh_dict['points'],
                                  geo_reference=self.geo_reference)
        else:
            #point_dict['pointlist'] = mesh_dict['vertices']
            geo = Geospatial_data(mesh_dict['vertices'],
                                  geo_reference=self.geo_reference)

        geo.export_points_file(ofile, absolute=True)

    def import_ungenerate_file(self, ofile, tag=None, region_tag=None):
        """
        Imports an ungenerate file, from arcGIS into mesh.
        This file describes many polygons.

        ofile is the name of the ungenerated file.
        Tag is a string name to be taggged on each segment.

        region_tag is the tag applied to the polygon regions.
        if it is a string the one value will be assigned to all regions
        if it is a list the first value in the list will be applied to the first polygon etc.
        WARNING: size of list and number of polygons isn't checked 

        WARNING values are assumed to be absolute.
        geo-refs are not taken into account..
        """

        dict = load_ungenerate(ofile)
        default_tag = Segment.get_default_tag()
        if tag is not None:
            Segment.set_default_tag(str(tag))

        if region_tag is None:
            self.addVertsSegs(dict)
        else:
            if not isinstance(region_tag, list):
                region_tag = [region_tag]*len(dict['polygons'])
            for a_tag, polygon in zip(region_tag, dict['polygons']):
                segment_tags = {tag: list(range(len(polygon)))}
                self.add_region_from_polygon(polygon, segment_tags=segment_tags,
                                             region_tag=a_tag)

        Segment.set_default_tag(default_tag)

        # change the tag back to it's default


########### IO CONVERTERS ##################
        """
        The dict fromat for IO with .tsh files is;
        (the triangulation)
        vertices: [[x1,y1],[x2,y2],...] (lists of doubles)
        vertex_attributes: [[a11,a12,...],[a21,a22],...] (lists of doubles)
        vertex_attribute_titles:[A1Title, A2Title ...] (A list of strings)
        segments: [[v1,v2],[v3,v4],...] (lists of integers) 
        segment_tags : [tag,tag,...] list of strings
        triangles : [(v1,v2,v3), (v4,v5,v6),....] lists of points
        triangle tags: [s1,s2,...] A list of strings
        triangle neighbors: [[t1,t2,t3], [t4,t5,t6],..] lists of triangles
        
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
        
        """

    def Mesh2IODict(self):
        """
        Convert the triangulation and outline info of a mesh to a dictionary
        structure
        """
        dict = self.Mesh2IOTriangulationDict()
        dict_mesh = self.Mesh2IOOutlineDict()
        for element in list(dict_mesh.keys()):
            dict[element] = dict_mesh[element]

        # add the geo reference
        dict['geo_reference'] = self.geo_reference
        return dict

    def Mesh2IOTriangulationDict(self):
        """
        Convert the Mesh to a dictionary of lists describing the
        triangulation variables;

        Used to produce .tsh file
        """
        meshDict = {}
        if self.tri_mesh is not None:
            meshDict['triangles'] = self.tri_mesh.triangles
            meshDict['triangle_tags'] = self.tri_mesh.triangle_tags
            # print "mesh meshDict['triangle_tags']", meshDict['triangle_tags']
            meshDict['triangle_neighbors'] = self.tri_mesh.triangle_neighbors
            meshDict['vertices'] = self.tri_mesh.vertices
            meshDict['vertex_attributes'] = self.tri_mesh.vertex_attributes
            meshDict['vertex_attribute_titles'] = \
                self.tri_mesh.vertex_attribute_titles
            meshDict['segments'] = self.tri_mesh.segments
            meshDict['segment_tags'] = self.tri_mesh.segment_tags
        else:
            meshDict['triangles'] = []
            meshDict['triangle_tags'] = []
            meshDict['triangle_neighbors'] = []
            meshDict['vertices'] = []
            meshDict['vertex_attributes'] = []
            meshDict['vertex_attribute_titles'] = []
            meshDict['segments'] = []
            meshDict['segment_tags'] = []
        # print "mesh.Mesh2IOTriangulationDict*)*)"
        # print meshDict
        # print "mesh.Mesh2IOTriangulationDict*)*)"

        return meshDict

    def Mesh2IOOutlineDict(self, userVertices=None,
                           userSegments=None,
                           holes=None,
                           regions=None):
        """
        Convert the mesh outline to a dictionary of the lists needed for the
        triang module;

        Note, this adds an index attribute to the user Vertex objects.

        Used to produce .tsh file and output to triangle
        """

        if userVertices is None:
            userVertices = self.getUserVertices()
        if userSegments is None:
            userSegments = self.getUserSegments()
        if holes is None:
            holes = self.getHoles()
        if regions is None:
            regions = self.getRegions()

        meshDict = {}
        # print "userVertices",userVertices
        # print "userSegments",userSegments
        pointlist = []
        pointattributelist = []
        index = 0
        for vertex in userVertices:
            vertex.index = index
            pointlist.append([vertex.x, vertex.y])
            pointattributelist.append(vertex.attributes)

            index += 1
        meshDict['points'] = pointlist
        meshDict['point_attributes'] = pointattributelist

        segmentlist = []
        segmenttaglist = []
        for seg in userSegments:
            segmentlist.append([seg.vertices[0].index, seg.vertices[1].index])
            segmenttaglist.append(seg.tag)
        meshDict['outline_segments'] = segmentlist
        meshDict['outline_segment_tags'] = segmenttaglist

        holelist = []
        for hole in holes:
            holelist.append([hole.x, hole.y])
        meshDict['holes'] = holelist

        regionlist = []
        regiontaglist = []
        regionmaxarealist = []
        for region in regions:
            regionlist.append([region.x, region.y])
            regiontaglist.append(region.getTag())

            if (region.getMaxArea() != None):
                regionmaxarealist.append(region.getMaxArea())
            else:
                regionmaxarealist.append(NOMAXAREA)
        meshDict['regions'] = regionlist
        meshDict['region_tags'] = regiontaglist
        meshDict['region_max_areas'] = regionmaxarealist
        # print "*(*("
        # print meshDict
        # print meshDict['regionlist']
        # print "*(*("
        return meshDict

    def IOTriangulation2Mesh(self, genDict):
        """
        Set the mesh attributes given an tsh IO dictionary
        """
        # Clear the current generated mesh values
        self.tri_mesh = None

        self.tri_mesh = Rigid_triangulation(
            genDict['triangles'], genDict['segments'], genDict['vertices'], genDict['triangle_tags'], genDict[
                'triangle_neighbors'], genDict['segment_tags'], genDict['vertex_attributes'], genDict['vertex_attribute_titles']
        )
        self.attributeTitles = genDict['vertex_attribute_titles']
        self.maxVertexIndex = len(genDict['vertices'])
        # print "self.maxVertexIndex ", self.maxVertexIndex

    def IOOutline2Mesh(self, genDict):
        """
        Set the outline (user Mesh attributes) given a IO tsh dictionary

        mesh is an instance of a mesh object
        """
        # Clear the current user mesh values
        self.clearUserSegments()
        self.userVertices = []
        self.Holes = []
        self.Regions = []

        # print "mesh.IOOutline2Mesh@#@#@#"
        # print "genDict",genDict
        # print "@#@#@#"

        #index = 0
        for point in genDict['points']:
            v = Vertex(point[0], point[1])
            #v.index = index
            #index +=1
            self.userVertices.append(v)

        #index = 0
        for seg, tag in zip(genDict['outline_segments'],
                            genDict['outline_segment_tags']):

            segObject = Segment(self.userVertices[int(seg[0])],
                                self.userVertices[int(seg[1])], tag=tag)
            #segObject.index = index
            #index +=1
            self.userSegments.append(segObject)

# Remove the loading of attribute info.
# Have attribute info added using least_squares in pyvolution
#         index = 0
#         for att in genDict['point_attributes']:
#             if att is None:
#                 self.userVertices[index].setAttributes([])
#             else:
#                 self.userVertices[index].setAttributes(att)
#            index += 1

        #index = 0
        for point in genDict['holes']:
            h = Hole(point[0], point[1])
            #h.index = index
            #index +=1
            self.holes.append(h)

        #index = 0
        for reg, att, maxArea in zip(
                genDict['regions'],
                genDict['region_tags'],
                genDict['region_max_areas']):
            if maxArea > 0:  # maybe I should ref NOMAXAREA? Prob' not though
                Object = Region(reg[0],
                                reg[1],
                                tag=att,
                                maxArea=maxArea)
            else:
                Object = Region(reg[0],
                                reg[1],
                                tag=att)

            #Object.index = index
            #index +=1
            self.regions.append(Object)


def importMeshFromFile(ofile):
    """returns a mesh object, made from a points file or .tsh/.msh file
    Often raises IOError,RuntimeError
    """
    newmesh = None
    if (ofile[-4:] == ".pts" or ofile[-4:] == ".txt" or
            ofile[-4:] == ".csv"):
        geospatial = Geospatial_data(ofile)
        dict = {}
        dict['points'] = geospatial.get_data_points(absolute=False)
        dict['outline_segments'] = []
        dict['outline_segment_tags'] = []
        dict['regions'] = []
        dict['region_tags'] = []
        dict['region_max_areas'] = []
        dict['holes'] = []
        newmesh = Mesh(geo_reference=geospatial.geo_reference)
        newmesh.IOOutline2Mesh(dict)
        counter = newmesh.removeDuplicatedUserVertices()
        if (counter > 0):
            log.critical(
                "%i duplicate vertices removed from dataset" % counter)
    elif (ofile[-4:] == ".tsh" or ofile[-4:] == ".msh"):
        dict = import_mesh_file(ofile)
        # print "********"
        # print "zq mesh.dict",dict
        # print "********"

        newmesh = Mesh()
        newmesh.IOOutline2Mesh(dict)
        newmesh.IOTriangulation2Mesh(dict)
    else:
        raise RuntimeError

    if 'geo_reference' in dict and not dict['geo_reference'] is None:
        newmesh.geo_reference = dict['geo_reference']

    return newmesh


def loadPickle(currentName):
    fd = open(currentName)
    mesh = pickle.load(fd)
    fd.close()
    return mesh


def square_outline(side_length=1, up="top", left="left", right="right",
                   down="bottom", regions=False):

    a = Vertex(0, 0)
    b = Vertex(0, side_length)
    c = Vertex(side_length, 0)
    d = Vertex(side_length, side_length)

    s2 = Segment(b, d, tag=up)
    s3 = Segment(b, a, tag=left)
    s4 = Segment(d, c, tag=right)
    s5 = Segment(a, c, tag=down)

    if regions:
        e = Vertex(side_length/ 2, side_length/ 2)
        s6 = Segment(a, e, tag=down + left)
        s7 = Segment(b, e, tag=up + left)
        s8 = Segment(c, e, tag=down + right)
        s9 = Segment(d, e, tag=up + right)
        r1 = Region(side_length/ 2,
                    3.*side_length/ 4, tag=up)
        r2 = Region(1.*side_length/ 4,
                    side_length/ 2, tag=left)
        r3 = Region(3.*side_length/ 4,
                    side_length/ 2, tag=right)
        r4 = Region(side_length/ 2, 1.*side_length/4, tag=down)
        mesh = Mesh(userVertices=[a, b, c, d, e],
                    userSegments=[s2, s3, s4, s5, s6, s7, s8, s9],
                    regions=[r1, r2, r3, r4])
    else:
        mesh = Mesh(userVertices=[a, b, c, d],
                    userSegments=[s2, s3, s4, s5])

    return mesh


def region_strings2ints(region_list):
    """Given a list of (x_int,y_int,tag_string) lists it returns a list of
    (x_int,y_int,tag_int) and a list to convert the tag_int's back to
    the tag_strings  
    """
    # Make sure "" has an index of 0
    region_list.reverse()
    region_list.append((1.0, 2.0, ""))
    region_list.reverse()
    convertint2string = []
    for i in range(len(region_list)):
        convertint2string.append(region_list[i][2])
        if len(region_list[i]) == 4:  # there's an area value
            region_list[i] = (region_list[i][0],
                              region_list[i][1], i, region_list[i][3])
        elif len(region_list[i]) == 3:  # no area value
            region_list[i] = (region_list[i][0], region_list[i][1], i)
        else:
            log.critical("The region list has a bad size")
            # raise an error ..
            raise Exception

    # remove "" from the region_list
    region_list.pop(0)

    return [region_list, convertint2string]


def region_ints2strings(region_list, convertint2string):
    """Reverses the transformation of region_strings2ints
    """
    # print 'region_ints2strings region_list', region_list

    returned_region_list = []
    # may not need (not region_list[0] == [])
    # or region_list[0] == [0.0]
    if (region_list[0].size > 0):  # or region_list[0] == [0.0]:
        # print "in loop"
        for i in range(len(region_list)):
            temp = region_list[i]
            returned_region_list.append(convertint2string[int(temp[0])])
    return returned_region_list


def segment_ints2strings(intlist, convertint2string):
    """Reverses the transformation of segment_strings2ints """

    stringlist = []
    for x in intlist:
        stringlist.append(convertint2string[int(x)])
    return stringlist


def segment_strings2ints(stringlist, preset):
    """Given a list of strings return a list of 0 to n ints which represent
    the strings and a converting list of the strings, indexed by 0 to n ints.
    Also, input an initial converting list of the strings
    Note, the converting list starts off with
    ["internal boundary", "external boundary", "internal boundary"]
    example input and output 
    input ["a","b","a","c"], ["c"]
    output [[2, 1, 2, 0], ['c', 'b', 'a']]

    the first element in the converting list will be
    overwritten with "".
    ?This will become the third item in the converting list?

    # note, order the initial converting list is important,
     since the index = the int tag
    """
    nodups = unique(stringlist)
    # note, order is important, the index = the int tag
    #preset = ["internal boundary", "external boundary"]
    # Remove the preset tags from the list with no duplicates
    nodups = [x for x in nodups if x not in preset]

    try:
        nodups.remove("")  # this has to go to zero
    except ValueError:
        pass

    # Add the preset list at the beginning of no duplicates
    preset.reverse()
    nodups.extend(preset)
    nodups.reverse()

    convertstring2int = {}
    convertint2string = []
    index = 0
    for x in nodups:
        convertstring2int[x] = index
        convertint2string.append(x)
        index += 1
    convertstring2int[""] = 0

    intlist = []
    for x in stringlist:
        intlist.append(convertstring2int[x])
    return [intlist, convertint2string]


def unique(s):
    """Return a list of the elements in s, but without duplicates.

    For example, unique([1,2,3,1,2,3]) is some permutation of [1,2,3],
    unique("abcabc") some permutation of ["a", "b", "c"], and
    unique(([1, 2], [2, 3], [1, 2])) some permutation of
    [[2, 3], [1, 2]].

    For best speed, all sequence elements should be hashable.  Then
    unique() will usually work in linear time.

    If not possible, the sequence elements should enjoy a total
    ordering, and if list(s).sort() doesn't raise TypeError it's
    assumed that they do enjoy a total ordering.  Then unique() will
    usually work in O(N*log2(N)) time.

    If that's not possible either, the sequence elements must support
    equality-testing.  Then unique() will usually work in quadratic
    time.
    """

    n = len(s)
    if n == 0:
        return []

    # Try using a dict first, as that's the fastest and will usually
    # work.  If it doesn't work, it will usually fail quickly, so it
    # usually doesn't cost much to *try* it.  It requires that all the
    # sequence elements be hashable, and support equality comparison.
    u = {}
    try:
        for x in s:
            u[x] = 1
    except TypeError:
        del u  # move on to the next method
    else:
        return list(u.keys())

    # We can't hash all the elements.  Second fastest is to sort,
    # which brings the equal elements together; then duplicates are
    # easy to weed out in a single pass.
    # NOTE:  Python's list.sort() was designed to be efficient in the
    # presence of many duplicate elements.  This isn't true of all
    # sort functions in all languages or libraries, so this approach
    # is more effective in Python than it may be elsewhere.
    try:
        t = list(s)
        t.sort()
    except TypeError:
        del t  # move on to the next method
    else:
        assert n > 0
        last = t[0]
        lasti = i = 1
        while i < n:
            if t[i] != last:
                t[lasti] = last = t[i]
                lasti += 1
            i += 1
        return t[:lasti]

    # Brute force is all that's left.
    u = []
    for x in s:
        if x not in u:
            u.append(x)
    return u

# From https://docs.python.org/3/howto/sorting.html
# Used with cmp_xy
def cmp_to_key(mycmp):
    'Convert a cmp= function into a key= function'
    class K:
        def __init__(self, obj, *args):
            self.obj = obj
        def __lt__(self, other):
            return mycmp(self.obj, other.obj) < 0
        def __gt__(self, other):
            return mycmp(self.obj, other.obj) > 0
        def __eq__(self, other):
            return mycmp(self.obj, other.obj) == 0
        def __le__(self, other):
            return mycmp(self.obj, other.obj) <= 0
        def __ge__(self, other):
            return mycmp(self.obj, other.obj) >= 0
        def __ne__(self, other):
            return mycmp(self.obj, other.obj) != 0
    return K



# FIXME (DSG-DSG)
# To do-
# Create a clear interface. eg
# have the interface methods more at the top of this file and add comments
# for the interface functions/methods, use function_name (not functionName),

# Currently
# function_name methods assume absolute values.  Geo-refs can be passed in.
#

# instead of functionName
if __name__ == "__main__":
    pass

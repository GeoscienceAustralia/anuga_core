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

import types
import exceptions
from Numeric import array, Float, Int


#class NoTrianglesError(exceptions.Exception): pass
 
#import load_mesh
from anuga.coordinate_transforms.geo_reference import Geo_reference, \
     DEFAULT_ZONE
from  anuga.load_mesh.loadASCII import NOMAXAREA, export_mesh_file, \
     import_mesh_file 
import anuga.alpha_shape.alpha_shape
from anuga.geospatial_data.geospatial_data import Geospatial_data, \
     ensure_geospatial, ensure_absolute, ensure_numeric
from anuga.mesh_engine.mesh_engine import generate_mesh

try:  
    import kinds  
except ImportError:  
    # Hand-built mockup of the things we need from the kinds package, since it
    # was recently removed from the standard Numeric distro.  Some users may  
    # not have it by default.  
    class _bunch:  
        pass  
         
    class _kinds(_bunch):  
        default_float_kind = _bunch()  
        default_float_kind.MIN = 2.2250738585072014e-308  #smallest +ve number
        default_float_kind.MAX = 1.7976931348623157e+308  
     
    kinds = _kinds()
    
SET_COLOUR='red'

#FIXME: this is not tested.
from anuga.utilities.numerical_tools import gradient



# 1st and third values must be the same
# FIXME: maybe make this a switch that the user can change? - DSG
initialconversions = ['', 'exterior', '']

#from os import sep
#sys.path.append('..'+sep+'pmesh')
#print "sys.path",sys.path

class MeshObject:
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
    def __init__(self,X,Y):
        __slots__ = ['x','y']
        self.x=X
        self.y=Y
        
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
        
        if (self.DistanceToPoint(Center)<Radius):
            return 1
        else:
            return 0
        
    def __repr__(self):
        return "(%f,%f)" % (self.x,self.y) 

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
    """
    A point on the mesh.
    Object attributes based on the Triangle program
    """
    def __init__(self,X,Y, attributes = None):
        __slots__ = ['x','y','attributes']
        
        assert (type(X) == types.FloatType or type(X) == types.IntType)
        assert (type(Y) == types.FloatType or type(Y) == types.IntType)
        self.x=X
        self.y=Y        
        self.attributes=[] 
        
        if attributes is None:
            self.attributes=[]
        else:
            self.attributes=attributes
    

    def setAttributes(self,attributes):
        """
        attributes is a list.
        """
        self.attributes = attributes
        
    VERTEXSQUARESIDELENGTH = 6
    def draw(self, canvas, tags, colour = 'black',scale = 1, xoffset = 0,
             yoffset =0, ):
        x =  scale*(self.x + xoffset)
        y = -1*scale*(self.y + yoffset)  # - since for a canvas - is up
        #print "draw x:", x
        #print "draw y:", y
        cornerOffset= self.VERTEXSQUARESIDELENGTH/2

        # A hack to see the vert tags
        # note: there will be many tags, since tags will not be removed
        #when zooming
        #canvas.create_text(x+ 2*cornerOffset,
        #                   y+ 2*cornerOffset,
        #                        text=tags)
        
        return canvas.create_rectangle(x-cornerOffset,
                                       y-cornerOffset,
                                       x+cornerOffset,
                                       y+cornerOffset,
                                       tags = tags,
                                       outline=colour,
                                       fill = 'white')
    
        #return tags
     
    def __repr__(self):
        return "[(%f,%f),%r]" % (self.x,self.y,self.attributes)
    
class Hole(Point):
    """
    A region of the mesh were no triangles are generated.
    Defined by a point in the hole enclosed by segments.
    """

    HOLECORNERLENGTH = 6
    
    def draw(self, canvas, tags, colour = 'purple',scale = 1, xoffset = 0,
             yoffset =0, ):
        x =  scale*(self.x + xoffset)
        y = -1*scale*(self.y + yoffset)  # - since for a canvas - is up
        #print "draw x:", x
        #print "draw y:", y
        cornerOffset= self.HOLECORNERLENGTH/2
        return canvas.create_oval(x-cornerOffset,
                                       y-cornerOffset,
                                       x+cornerOffset,
                                       y+cornerOffset,
                                       tags = tags,
                                       outline=colour,
                                       fill = 'white')
    
class Region(Point):
    """ 
    A region of the mesh, defined by a point in the region
    enclosed by segments. Used to tag areas.
    """
    CROSSLENGTH = 6
    TAG = 0
    MAXAREA = 1
    
    def __init__(self,X,Y, tag = None, maxArea = None):
        """Precondition: tag is a string and maxArea is a real
        """
        # This didn't work.  
        #super(Region,self)._init_(self,X,Y)
        self.x=X
        self.y=Y    
        self.attributes =[] # index 0 is the tag string
                            #optoinal index 1 is the max triangle area
                            #NOTE the size of this attribute is assumed
                            # to be 1 or 2 in regionstrings2int
        if tag is None:
            self.attributes.append("")
        else:
            self.attributes.append(tag) #this is a string
            
        if maxArea is not None:
            self.setMaxArea(maxArea) # maxArea is a number
            
    def getTag(self,):
        return self.attributes[self.TAG]
    
    def setTag(self,tag):
        self.attributes[self.TAG] = tag
        
    def getMaxArea(self):
        """ Returns the Max Area of a Triangle or
        None, if the Max Area has not been set.
        """
        if self.isMaxArea():
            return self.attributes[1]
        else:
            return None
    
    def setMaxArea(self,MaxArea):
        if MaxArea is not None:
            if self.isMaxArea(): 
                self.attributes[self.MAXAREA] = float(MaxArea)
            else:
                self.attributes.append( float(MaxArea) )
    
    def deleteMaxArea(self):
        if self.isMaxArea():
            self.attributes.pop(self.MAXAREA)
            
    def isMaxArea(self):
        return len(self.attributes)> 1
    
    def draw(self, canvas, tags, scale=1, xoffset = 0, yoffset =0,
             colour = "black"):
        """
        Draw a black cross, returning the objectID
        """
        x =  scale*(self.x + xoffset)
        y = -1*scale*(self.y + yoffset) 
        cornerOffset= self.CROSSLENGTH/2
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
                                     tags = tags,
                                     outline = colour,fill = '')
    
    def __repr__(self):
        if self.isMaxArea():
            area = self.getMaxArea() 
            return "(%f,%f,%s,%f)" % (self.x,self.y,
                                      self.getTag(), self.getMaxArea())
        else:
            return "(%f,%f,%s)" % (self.x,self.y,
                                   self.getTag())
        
class Triangle(MeshObject):
    """
    A triangle element, defined by 3 vertices.
    Attributes based on the Triangle program.
    """

    def __init__(self, vertex1, vertex2, vertex3, attribute = None,
                 neighbors = None ):
        """
        Vertices, the initial arguments, are listed in counterclockwise order.
        """
        self.vertices= [vertex1,vertex2, vertex3 ]
        
        if attribute is None:
            self.attribute =""
        else:
            self.attribute = attribute #this is a string
            
        if neighbors is None:
            self.neighbors=[]
        else:
            self.neighbors=neighbors

    def replace(self,new_triangle):
        self = new_triangle

    def longestSideID(self):
        ax = self.vertices[0].x
        ay = self.vertices[0].y
        
        bx = self.vertices[1].x
        by = self.vertices[1].y
        
        cx = self.vertices[2].x
        cy = self.vertices[2].y

        lenA = ((cx-bx)**2+(cy-by)**2)**0.5
        lenB = ((ax-cx)**2+(ay-cy)**2)**0.5
        lenC = ((bx-ax)**2+(by-ay)**2)**0.5
 
        len = [lenA,lenB,lenC]
        return len.index(max(len))

    def rotate(self,offset):
        """
        permute the order of the sides of the triangle
        offset must be 0,1 or 2
        """

        if offset == 0:
            pass
        else:
            if offset == 1:
                self.vertices = [self.vertices[1],self.vertices[2],
                                 self.vertices[0]]
                self.neighbors = [self.neighbors[1],self.neighbors[2],
                                  self.neighbors[0]]
            if offset == 2:
                self.vertices = [self.vertices[2],self.vertices[0],
                                 self.vertices[1]]
                self.neighbors = [self.neighbors[2],self.neighbors[0],
                                  self.neighbors[1]]

    def rotate_longest_side(self):
        self.rotate(self.longestSideID())

    def getVertices(self):
        return self.vertices
    
    def get_vertices(self):
        """
        Return a list of the vertices.  The x and y values will be relative
        Easting and Northings for the zone of the current geo_ref.
        """
        return self.vertices
    
    def calcArea(self):
        ax = self.vertices[0].x
        ay = self.vertices[0].y
        
        bx = self.vertices[1].x
        by = self.vertices[1].y
        
        cx = self.vertices[2].x
        cy = self.vertices[2].y
        
        return abs((bx*ay-ax*by)+(cx*by-bx*cy)+(ax*cy-cx*ay))/2
    
    def calcP(self):
        #calculate the perimeter
        ax = self.vertices[0].x
        ay = self.vertices[0].y
        
        bx = self.vertices[1].x
        by = self.vertices[1].y
        
        cx = self.vertices[2].x
        cy = self.vertices[2].y

        a =  ((cx-bx)**2+(cy-by)**2)**0.5 
        b =  ((ax-cx)**2+(ay-cy)**2)**0.5 
        c =  ((bx-ax)**2+(by-ay)**2)**0.5 

        return a+b+c
            
    def setNeighbors(self,neighbor1=None, neighbor2=None, neighbor3=None):
        """
        neighbor1 is the triangle opposite vertex1 and so on.
        Null represents no neighbor
        """
        self.neighbors = [neighbor1, neighbor2, neighbor3]

    def setAttribute(self,attribute):
        """
        neighbor1 is the triangle opposite vertex1 and so on.
        Null represents no neighbor
        """
        self.attribute = attribute #this is a string
        
    def __repr__(self):
        return "[%s,%s]" % (self.vertices,self.attribute)
        

    def draw(self, canvas, tags, scale=1, xoffset=0, yoffset=0,
             colour="green"):
        """
        Draw a triangle, returning the objectID
        """
        return canvas.create_polygon(scale*(self.vertices[1].x + xoffset),
                                     scale*-1*(self.vertices[1].y + yoffset),
                                     scale*(self.vertices[0].x + xoffset),
                                     scale*-1*(self.vertices[0].y + yoffset),
                                     scale*(self.vertices[2].x + xoffset),
                                     scale*-1*(self.vertices[2].y + yoffset),
                                     tags = tags,
                                     outline = colour,fill = '')

class Segment(MeshObject):
    """
    Segments are edges whose presence in the triangulation is enforced.
    
    """
    def __init__(self, vertex1, vertex2, tag = None ):
        """
        Each segment is specified by listing the vertices of its endpoints
        The vertices are Vertex objects
        """
        assert(vertex1 != vertex2)
        self.vertices = [vertex1,vertex2 ]
        
        if tag is None:
            self.tag = self.__class__.default
        else:
            self.tag = tag #this is a string 
        
    def __repr__(self):
        return "[%s,%s]" % (self.vertices,self.tag)
            
        
    def draw(self, canvas, tags,scale=1, xoffset=0, yoffset=0,colour='blue'):
        x=[]
        y=[]
        for end in self.vertices:
            #end.draw(canvas,scale, xoffset, yoffset ) # draw the vertices
            x.append(scale*(end.x + xoffset))
            y.append(-1*scale*(end.y + yoffset)) # - since for a canvas - is up
        
        return canvas.create_line(x[0], y[0], x[1], y[1],
                                  tags = tags,fill=colour)
    def set_tag(self,tag):
        self.tag = tag
        
    # Class methods
    def set_default_tag(cls, default):
        cls.default = default 
    
    def get_default_tag(cls):
        return cls.default
    
    set_default_tag = classmethod(set_default_tag)  
    get_default_tag = classmethod(get_default_tag)

Segment.set_default_tag("")       

class Rigid_triangulation:
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
        self.triangle_tags = triangle_tags # list of strings
        self.segments = ensure_numeric(segments) 
        
        self.segment_tags = segment_tags # list of strings
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
                outline = colour,fill = ''
                )
 
        
class Mesh:
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
        self.meshTriangles=[] 
        self.attributeTitles=[] 
        self.meshSegments=[]
        self.meshVertices=[]

        # Rigid
        self.tri_mesh=None
        
        self.setID={}
        #a dictionary of names.
        #multiple sets are allowed, but the gui does not yet
        #support this
        
        self.setID['None']=0
        #contains the names of the sets pointing to the indexes
        #in the list. 
        
	self.sets=[[]]
        #Contains the lists of triangles (triangle sets)

       
        self.visualise_graph = True

        if userSegments is None:
            self.userSegments=[]
        else:
            self.userSegments=userSegments
        self.alphaUserSegments=[]
            
        if userVertices is None:
            self.userVertices=[]
        else:
            self.userVertices=userVertices
            
        if holes is None:
            self.holes=[]
        else:
            self.holes=holes
            
        if regions is None:
            self.regions=[]
        else:
            self.regions=regions

        if geo_reference is None:
            self.geo_reference = Geo_reference(DEFAULT_ZONE,0,0) 
        else:
            self.geo_reference = geo_reference

        self.shape = None
            
    def __cmp__(self,other):
        
        # A dic for the initial m
        dic = self.Mesh2triangList()
        dic_mesh = self.Mesh2MeshList()
        for element in dic_mesh.keys():
            dic[element] = dic_mesh[element]
        for element in dic.keys():
            dic[element].sort()
            
        # A dic for the exported/imported m
        dic_other = other.Mesh2triangList()
        dic_mesh = other.Mesh2MeshList()
        for element in dic_mesh.keys():
            dic_other[element] = dic_mesh[element]
        for element in dic.keys():
            dic_other[element].sort()

        #print "dsg************************8"
        #print "dic ",dic
        #print "*******8"
        #print "mesh",dic_other
        #print "dic.__cmp__(dic_o)",dic.__cmp__(dic_other)
        #print "dsg************************8"
        
        return (dic.__cmp__(dic_other))
    
    def addUserPoint(self, pointType, x,y):
        if pointType == Vertex:
            point = self.addUserVertex(x,y)
        if pointType == Hole:
            point = self._addHole(x,y)
        if pointType == Region:
            point = self._addRegion(x,y)
        return point
    
    def addUserVertex(self, x,y):
        v=Vertex(x, y)
        self.userVertices.append(v)
        return v

    def _addHole(self, x,y):
        h=Hole(x, y)
        self.holes.append(h)
        return h
   
    def add_hole(self, x,y, geo_reference=None):
        """
        adds a point, which represents a hole.

        The point data can have it's own geo_refernece.
        If geo_reference is None the data is asumed to be absolute
        """
        [[x,y]] = self.geo_reference.change_points_geo_ref([x,y],
                                                 points_geo_ref=geo_reference)
        return self._addHole(x, y)

    def _addRegion(self, x,y):
        h=Region(x, y)
        self.regions.append(h)
        return h
   
    def add_region(self, x,y, geo_reference=None):
        """
        adds a point, which represents a region.

        The point data can have it's own geo_refernece.
        If geo_reference is None the data is asumed to be absolute
        """
        #FIXME: have the user set the tag and resolution here,
        # but still return the instance, just in case.
        [[x,y]] = self.geo_reference.change_points_geo_ref([x,y],
                                                 points_geo_ref=geo_reference)
        return self._addRegion(x, y)

    def build_grid(self,  vert_rows, vert_columns):
        """
        Build a grid with vert_rows number of vertex rows and
        vert_columns number if vertex columns

        Grid spacing of 1, the origin is the lower left hand corner.

        FIXME(DSG-DSG) no test.
        """

        for i in range(vert_rows):
            for j in range(vert_columns):
                self.addUserVertex(j,i)
        self.auto_segment()
        self.generateMesh(mode = "Q", minAngle=20.0)
        
    # Depreciated
    def addRegionEN(self, x,y):
        print "depreciated, use add_region"
        return self.add_region(x,y)

    
    def add_vertices(self, point_data):
        """
        Add user vertices.

        The point_data can be a list of (x,y) values, a numeric
        array or a geospatial_data instance.
        """
        point_data = ensure_geospatial(point_data)
        #print "point_data",point_data 
        # get points relative to the mesh geo_ref
        points = point_data.get_data_points(geo_reference=self.geo_reference)
    
        for point in points:
            v=Vertex(point[0], point[1])
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
        if max_triangle_area is None:
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
        from anuga.utilities.polygon import point_in_polygon
        
        #get absolute values
        if geo_reference is not None:
            polygon = geo_reference.get_absolute(polygon)
        # polygon is now absolute
        #print "polygon  should be absolute",polygon
        
        #create points, segs and tags
        region_dict = {}
        region_dict['points'] = polygon
        
        #Create segments
        #E.g. [[0,1], [1,2], [2,3], [3,0]]
        #from polygon
        #[0,1,2,3]
        segments = []
        N = len(polygon)
        for i in range(N):
            lo = i
            hi = (lo + 1) % N
            segments.append( [lo, hi] ) 
        region_dict['segments'] = segments
        region_dict['segment_tags'] = self._tag_dict2list(segment_tags, N)
       
    
        self.addVertsSegs(region_dict) #this is passing absolute values

        if region is True:
            #get inner point - absolute values
            inner_point = point_in_polygon(polygon)
            inner = self.add_region(inner_point[0], inner_point[1],
                                    geo_reference=None) 
        elif hole is True:
            #get inner point - absolute values
            inner_point = point_in_polygon(polygon)
            inner = self.add_hole(inner_point[0], inner_point[1],
                                    geo_reference=None)
        else:
            inner = None
            
        return inner

    def _tag_dict2list(self, tags, number_of_segs):
        """
        Convert a tag dictionary from this sort of format;
        #{'wall':[0,3],'ocean':[2]}

        To a list format
        # ['wall', '', 'ocean', 'wall']

        Note: '' is a default value of nothing
        """
        # FIXME (DSG-DSG): Using '' as a default isn't good.
        #Try None.
        # Due to this default this method is too connected to
        # _add_area_from_polygon
        
        segment_tags = ['']*number_of_segs
        if tags is not None:
            for key in tags:
                indices = tags[key]
                for i in indices:
                    segment_tags[i] = key
        return segment_tags
        
    def add_circle(self, center, radius, segment_count=100,
                   center_geo_reference=None, tag = None,
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
        factor = 2* math.pi/segment_count
        for cut in range(segment_count):
             cuts.append(cut*factor)

        polygon = []
        for cut in cuts:
            
            x = center[0] + radius * math.cos(cut)
            y = center[1] + radius * math.sin(cut)
            polygon.append([x,y])
        # build the tags
        tags = {tag:range(segment_count)}
        
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
        #to do
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
        #self.userVertices = self.geo_reference.change_points_geo_ref( \
        #self.userVertices)
        #self.holes = self.geo_reference.change_points_geo_ref(self.holes)
        #self.regions = self.geo_reference.change_points_geo_ref(self.regions)
        # The above will not work.
        # since userVertices (etc) is a list of point objects,
        #not a list of lists.
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
                                  segments, segment_tags = None):
        """
        Add an outline of the mesh.
        Vertices is a list of points/ a standard representation of points.
        Segments is a list of tuples of integers.  Each tuple defines the
           start and end of the segment by it's vertex index, in relation to
           the list of vertices.
        segment_tags is an optional dictionary which is used to add tags to
           the segments.  The key is the tag name, value is the list of segment
           indexes the tag will apply to.
           eg. {'wall':[0,3],'ocean':[2]}
           
        """
        #make sure the points are absolute
        points = ensure_absolute(points)
        
        #create points, segs and tags
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
        if not outlineDict.has_key('segment_tags'):
            outlineDict['segment_tags'] = []
            for i in range(len(outlineDict['segments'])):
                outlineDict['segment_tags'].append('')
        #print "outlineDict['segment_tags']",outlineDict['segment_tags']
        #print "outlineDict['points']",outlineDict['points']
        #print "outlineDict['segments']",outlineDict['segments']
        
        i_offset = len(self.userVertices)
        #print "self.userVertices",self.userVertices 
        #print "index_offset",index_offset 
        for point in outlineDict['points']:
            v=Vertex(point[0]-self.geo_reference.xllcorner,
                     point[1]-self.geo_reference.yllcorner)
            self.userVertices.append(v)
            
        for seg,seg_tag in map(None,outlineDict['segments'],
                       outlineDict['segment_tags']):
            segObject = Segment(self.userVertices[int(seg[0])+i_offset],
                           self.userVertices[int(seg[1])+i_offset] )
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
        pointlist=[]
        for vertex in self.userVertices:
            pointlist.append([vertex.x,vertex.y])
        spat = Geospatial_data(pointlist, geo_reference=self.geo_reference)
        return spat.get_data_points(absolute=absolute)
    
    def getUserSegments(self):
        allSegments = self.userSegments + self.alphaUserSegments
        #print "self.userSegments",self.userSegments
        #print "self.alphaUserSegments",self.alphaUserSegments
        #print "allSegments",allSegments
        return allSegments
    
    def deleteUserSegments(self,seg):
        if self.userSegments.count(seg) == 0:
            self.alphaUserSegments.remove(seg)
            pass
        else:
            self.userSegments.remove(seg)
            
    def clearUserSegments(self):
        self.userSegments = []
        self.alphaUserSegments = []

       #FIXME see where this is used. return an array instead
    def getTriangulation(self):
        #return self.meshTriangles
        return self.tri_mesh.triangles.tolist()
    
    def getMeshVertices(self):
        #return self.meshVertices
        return self.tri_mesh.vertices
 
    def getMeshVerticeAttributes(self):
        #return self.meshVertices
        return self.tri_mesh.vertex_attributes
    
    def getMeshSegments(self):
        #return self.meshSegments
        return self.tri_mesh.segments
    
    def getMeshSegmentTags(self):
        #return self.meshSegments
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
    
    def addUserSegment(self, v1,v2):
        """
        PRECON: A segment between the two vertices is not already present.
        Check by calling isUserSegmentNew before calling this function.
        
        """
        s=Segment( v1,v2)
        self.userSegments.append(s)
        return s
        
    def generate_mesh(self,
                      maximum_triangle_area="",
                      minimum_triangle_angle=28.0,
                      verbose=True):
        if verbose is True:
            silent = ''
        else:
            silent = 'Q'
        self.generateMesh(mode = silent +"pzq"+str(minimum_triangle_angle)
                                  +"a"+str(maximum_triangle_area)
                                  +"a")
        #The last a is so areas for regions will be used
        
    def generateMesh(self, mode = None, maxArea = None, minAngle= None,
                     isRegionalMaxAreas = True):
        """
        Based on the current user vaules, holes and regions
        generate a new mesh
        mode is a string that sets conditions on the mesh generations
        see triangle_instructions.txt for a definition of the commands
        
        PreCondition: maxArea is a double between 1e-20 and 1e30 or is a
        string.
        """
        #print "mode ",mode
        if mode == None:
            self.mode = ""
        else:
            self.mode = mode 
        
        if self.mode.find('p') < 0:
            self.mode += 'p' #p - read a planar straight line graph.
            #there must be segments to use this switch
            # TODO throw an aception if there aren't seg's
            # it's more comlex than this.  eg holes
        if self.mode.find('z') < 0:
            self.mode += 'z' # z - Number all items starting from zero
                             # (rather than one)
        if self.mode.find('n'):
            self.mode += 'n' # n - output a list of neighboring triangles 
        if self.mode.find('A') < 0:
            self.mode += 'A' # A - output region attribute list for triangles

        if not self.mode.find('V') and not self.mode.find('Q'):
            self.mode += 'V' # V - output info about what Triangle is doing
        
        if self.mode.find('q') < 0 and minAngle is not None:
            #   print "**********8minAngle******** ",minAngle
            min_angle = 'q' + str(minAngle)
            self.mode += min_angle # z - Number all items starting from zero
                             # (rather than one)
        if maxArea != None:
            self.mode += 'a' + str(maxArea)
            try:
                self.mode += 'a' + '%20.20f' %maxArea
            except TypeError:
                self.mode += 'a' + str(maxArea)
            #print "self.mode", self.mode
        #FIXME (DSG-DSG) This isn't explained. 
        if isRegionalMaxAreas:
            self.mode += 'a'
        #print "mesh#generateMesh# self.mode",self.mode  
        meshDict = self.Mesh2triangList()

        #FIXME (DSG-DSG)  move below section into generate_mesh.py
        #                  & 4 functions eg segment_strings2ints
        # Actually, because of region_list.append((1.0,2.0,""))
        # don't move it, without careful thought
        #print "*************************!@!@ This is going to triangle   !@!@"
        #print meshDict
        #print "************************!@!@ This is going to triangle   !@!@"

        #print "meshDict['segmenttaglist']", meshDict['segmenttaglist']
        [meshDict['segmenttaglist'],
         segconverter] =  segment_strings2ints(meshDict['segmenttaglist'],
                                             initialconversions)
        #print "regionlist",meshDict['regionlist']
        [meshDict['regionlist'],
         regionconverter] =  region_strings2ints(meshDict['regionlist'])
        #print "%%%%%%%%%%%%%%%%%%%%%%%%%%%regionlist",meshDict['regionlist']
        #print "meshDict['segmenttaglist']", meshDict['segmenttaglist'
        #print "self.mode", self.mode
        generatedMesh = generate_mesh(
                              meshDict['pointlist'],
                              meshDict['segmentlist'],
                              meshDict['holelist'],
                              meshDict['regionlist'],
                              meshDict['pointattributelist'],
                              meshDict['segmenttaglist'],
                              self.mode,
                              meshDict['pointlist'])
        #print "%%%%%%%%%%%%%%%%%%%%%%%%%%%generated",generatedMesh
        generatedMesh['qaa'] = 1
        if generatedMesh['generatedsegmentmarkerlist'] is not None:
            generatedMesh['generatedsegmentmarkerlist'] = \
              segment_ints2strings(generatedMesh['generatedsegmentmarkerlist'],
                                   segconverter)
        #print "processed gen",generatedMesh['generatedsegmentmarkerlist']
        #print "pmesh mesh generatedMesh['generatedtriangleattributelist']", generatedMesh['generatedtriangleattributelist']
        if generatedMesh['generatedtriangleattributelist'] is not None:
            generatedMesh['generatedtriangleattributelist'] = \
            region_ints2strings(generatedMesh['generatedtriangleattributelist'],
                                  regionconverter)

        #print "pmesh mesh generatedMesh['generatedtriangleattributelist']", generatedMesh['generatedtriangleattributelist']
        #FIXME (DSG-DSG)  move above section into generate_mesh.py
       
        if generatedMesh['generatedpointattributelist'] is None or \
               generatedMesh['generatedpointattributelist'].shape[1] ==0:
            self.attributeTitles = []
        generatedMesh['generatedpointattributetitlelist']= \
                                            self.attributeTitles
        #print "################  FROM TRIANGLE"
        #print "generatedMesh",generatedMesh
        #print "##################"
        self.setTriangulation(generatedMesh)
    
    def clearTriangulation(self):

        #Clear the current generated mesh values
        self.meshTriangles=[] 
        self.meshSegments=[]
        self.meshVertices=[]

    def removeDuplicatedUserVertices(self):
        """Pre-condition: There are no user segments
        This function will keep the first duplicate
        """
        assert self.getUserSegments() == []
        self.userVertices, counter =  self.removeDuplicatedVertices(
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
        t.sort(Point.cmp_xy)
    
        length = len(t)
        behind = 0
        ahead  = 1
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

        #Remove the attribute that this function added
        for v in Vertices:
            del v.dupindex
        return Vertices,counter

    # FIXME (DSG-DSG) Move this to geospatial
    def thinoutVertices(self, delta):
        """Pre-condition: There are no user segments
        This function will keep the first duplicate
        """
        assert self.getUserSegments() == []
        #t = self.userVertices
        #self.userVertices =[]
        boxedVertices = {}
        thinnedUserVertices =[]
        delta = round(delta,1)
        
        for v in self.userVertices :
            # tag is the center of the boxes
            tag = (round(v.x/delta,0)*delta,round(v.y/delta,0)*delta)
            #this creates a dict of lists of faces, indexed by tag
            boxedVertices.setdefault(tag,[]).append(v)

        for [tag,verts] in boxedVertices.items():
            min = delta
            bestVert = None
            tagVert = Vertex(tag[0],tag[1])
            for v in verts:
                dist = v.DistanceToPoint(tagVert)
                if (dist<min):
                    min = dist
                    bestVert = v
            thinnedUserVertices.append(bestVert)
        self.userVertices =thinnedUserVertices
        
            
    def isUserSegmentNew(self, v1,v2):
        identicalSegs= [x for x in self.getUserSegments() \
                        if (x.vertices[0] == v1 and x.vertices[1] == v2)
        or (x.vertices[0] == v2 and x.vertices[1] == v1) ]
        
        return len(identicalSegs) == 0

        
    def deleteSegsOfVertex(self, delVertex):
        """
        Delete this vertex and any segments that connect to it.
        """
        #Find segments that connect to delVertex
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
        if isinstance(MeshObject, Vertex ):
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
        
        pointlist=[]
        pointattributelist=[]
        index = 0
        for vertex in userVertices:
            vertex.index = index 
            pointlist.append((vertex.x,vertex.y))
            pointattributelist.append((vertex.attributes))
            
            index += 1
        meshDict['pointlist'] = pointlist
        meshDict['pointattributelist'] = pointattributelist

        segmentlist=[]
        segmenttaglist=[]
        for seg in userSegments:
            segmentlist.append((seg.vertices[0].index,seg.vertices[1].index))
            segmenttaglist.append(seg.tag)
        meshDict['segmentlist'] =segmentlist 
        meshDict['segmenttaglist'] =segmenttaglist
        
        holelist=[]
        for hole in holes:
            holelist.append((hole.x,hole.y)) 
        meshDict['holelist'] = holelist
        
        regionlist=[]
        for region in regions:
            if (region.getMaxArea() != None): 
                regionlist.append((region.x,region.y,region.getTag(),
                               region.getMaxArea()))
            else:
                regionlist.append((region.x,region.y,region.getTag()))
        meshDict['regionlist'] = regionlist
        #print "*(*("
        #print meshDict
        #print meshDict['regionlist']
        #print "*(*("
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
        #print "old meshDict['generatedpointattributetitlelist']",meshDict['generatedpointattributetitlelist']
        #print "self.tri_mesh", self.tri_mesh
        if self.tri_mesh is not None:
            #print "self.tri_mesh.triangles.tolist()", self.tri_mesh.triangles.tolist()
            meshDict['generatedtrianglelist'] = self.tri_mesh.triangles.tolist()
            meshDict['generatedtriangleattributelist'] = self.tri_mesh.triangle_tags
            meshDict['generatedtriangleneighborlist'] = self.tri_mesh.triangle_neighbors.tolist()
            meshDict['generatedpointlist'] = self.tri_mesh.vertices.tolist()
            if  self.tri_mesh.vertex_attributes == []:
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
            
        #print "new meshDict['generatedpointattributetitlelist']",meshDict['generatedpointattributetitlelist']
        #print "mesh.Mesh2MeshList*)*)"
        #print meshDict
        #print "mesh.Mesh2MeshList*)*)"

        return meshDict

                               
    def Mesh2MeshDic(self):
        """
        Convert the user and generated info of a mesh to a dictionary
        structure
        """
        dic = self.Mesh2triangList()
        dic_mesh = self.Mesh2MeshList()
        for element in dic_mesh.keys():
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
            genDict['generatedtrianglelist']
            ,genDict['generatedsegmentlist']
            ,genDict['generatedpointlist']
            ,genDict['generatedtriangleattributelist']
            ,genDict['generatedtriangleneighborlist']
            ,genDict['generatedsegmentmarkerlist']
            ,genDict['generatedpointattributelist']
            ,genDict['generatedpointattributetitlelist']
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
        #Clear the current user mesh values
        self.clearUserSegments()
        self.userVertices=[]
        self.Holes=[]
        self.Regions=[]

        #print "mesh.setMesh@#@#@#"
        #print genDict
        #print "@#@#@#"
        
        #index = 0
        for point in genDict['pointlist']:
            v=Vertex(point[0], point[1])
            #v.index = index
            #index +=1
            self.userVertices.append(v)

        #index = 0
        for seg,tag in map(None,genDict['segmentlist'],
                           genDict['segmenttaglist']):
            segObject = Segment( self.userVertices[int(seg[0])],
                           self.userVertices[int(seg[1])], tag = tag )
            #segObject.index = index
            #index +=1
            self.userSegments.append(segObject)

# Remove the loading of attribute info.
# Have attribute info added using least_squares in pyvolution
#         index = 0
#         for att in genDict['pointattributelist']:
#             if att == None:
#                 self.userVertices[index].setAttributes([])
#             else:
#                 self.userVertices[index].setAttributes(att)
#            index += 1
        
        #index = 0
        for point in genDict['holelist']:
            h=Hole(point[0], point[1])
            #h.index = index
            #index +=1
            self.holes.append(h)

        #index = 0
        for reg,att,maxArea in map(None,
                                   genDict['regionlist'],
                                   genDict['regionattributelist'],
                                   genDict['regionmaxarealist']):
            Object = Region( reg[0],
                             reg[1],
                             tag = att,
                             maxArea = maxArea)
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
        if self.isUserSegmentNew(s1.vertices[0],s1.vertices[1]):
            newsegs.append(s1)
        if self.isUserSegmentNew(s2.vertices[0],s2.vertices[1]):
            newsegs.append(s2)
        if self.isUserSegmentNew(s3.vertices[0],s3.vertices[1]):
            newsegs.append(s3)
        #DSG!!!
        self.userSegments.extend(newsegs)
        return newsegs

    
    def savePickle(self, currentName):
        fd = open(currentName, 'w')
        pickle.dump(self,fd)
        fd.close()

    def auto_segmentFilter(self,raw_boundary=True,
                    remove_holes=False,
                    smooth_indents=False,
                    expand_pinch=False):
        """
        Change the filters applied on the alpha shape boundary
        """
        if self.shape is None:
            return [],[],0.0
        return self._boundary2mesh(raw_boundary=raw_boundary,
                    remove_holes=remove_holes,
                    smooth_indents=smooth_indents,
                    expand_pinch=expand_pinch)
    
        
    
    def auto_segment(self, alpha = None,
                    raw_boundary=True,
                    remove_holes=False,
                    smooth_indents=False,
                    expand_pinch=False): 
        """
        Precon: There must be 3 or more vertices in the userVertices structure
        """
        self._createBoundary(alpha=alpha)
        return self._boundary2mesh(raw_boundary=raw_boundary,
                    remove_holes=remove_holes,
                    smooth_indents=smooth_indents,
                    expand_pinch=expand_pinch)

    def _createBoundary(self,alpha=None):
        """
        """
        points=[]
        for vertex in self.getUserVertices():
            points.append((vertex.x,vertex.y))
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
        #print "boundary_segs",boundary_segs
        segs2delete = self.alphaUserSegments
        #FIXME(DSG-DSG) this algorithm needs comments
        new_segs = {}
        #alpha_segs = []
        #user_segs = []
        for seg in boundary_segs:
            v1 = self.userVertices[int(seg[0])]
            v2 = self.userVertices[int(seg[1])]
            boundary_seg = Segment(v1, v2)
            new_segs[(v1,v2)] = boundary_seg

        for user_seg in self.userSegments:
            if new_segs.has_key((user_seg.vertices[0],
                                user_seg.vertices[1])):
                del new_segs[user_seg.vertices[0],
                                user_seg.vertices[1]]
            elif new_segs.has_key((user_seg.vertices[1],
                                user_seg.vertices[0])):
                del new_segs[user_seg.vertices[1],
                                user_seg.vertices[0]]
                
        optimum_alpha = self.shape.get_alpha()
        alpha_segs_no_user_segs  = new_segs.values()
        self.alphaUserSegments = alpha_segs_no_user_segs
        return alpha_segs_no_user_segs, segs2delete, optimum_alpha
    
    def _boundary2mesh_old(self, raw_boundary=True,
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
        #print "boundary_segs",boundary_segs
        segs2delete = self.alphaUserSegments

        #FIXME(DSG-DSG) this algorithm needs comments
        #FIXME(DSG-DSG) can it be sped up?  It's slow
        new_segs = []
        alpha_segs = []
        user_segs = []
        for seg in boundary_segs:
            v1 = self.userVertices[int(seg[0])]
            v2 = self.userVertices[int(seg[1])]
            alpha_seg = self.representedAlphaUserSegment(v1, v2)
            user_seg = self.representedUserSegment(v1, v2)
            #DSG!!!
            assert not(not (alpha_seg == None) and not (user_seg == None))
            if not alpha_seg == None:
                alpha_segs.append(alpha_seg)
            elif not user_seg  == None:
                user_segs.append(user_seg)
            else:
                unique_seg = Segment(v1, v2)
                new_segs.append(unique_seg)
                
            for seg in alpha_segs:
                try:
                    segs2delete.remove(seg)
                except:
                    pass
        
        self.alphaUserSegments = []
        self.alphaUserSegments.extend(new_segs)
        self.alphaUserSegments.extend(alpha_segs)

        optimum_alpha = self.shape.get_alpha()
        # need to draw newsegs
        return new_segs, segs2delete, optimum_alpha
    
    def representedAlphaUserSegment(self, v1,v2):
        identicalSegs= [x for x in self.alphaUserSegments \
                        if (x.vertices[0] == v1 and x.vertices[1] == v2)
        or (x.vertices[0] == v2 and x.vertices[1] == v1) ]

        if identicalSegs == []:
            return None
        else:
            # Only return the first one.
            return identicalSegs[0]
    
    def representedUserSegment(self, v1,v2):
        identicalSegs= [x for x in self.userSegments \
                        if (x.vertices[0] == v1 and x.vertices[1] == v2)
        or (x.vertices[0] == v2 and x.vertices[1] == v1) ]

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
            if self.isUserSegmentNew(v1,v2):            
                newseg = Segment(v1, v2)
                newsegs.append(newseg)
            v1 = v2

        #Connect last point to the first
        v2 = self.userVertices[0]        
        if self.isUserSegmentNew(v1,v2):            
                newseg = Segment(v1, v2)
                newsegs.append(newseg)
            

        #Update list of user segments        
        #DSG!!!
        self.userSegments.extend(newsegs)                
        return newsegs
    
    def normaliseMesh(self,scale, offset, height_scale):
        [xmin, ymin, xmax, ymax] = self.boxsize()
        [attmin0, attmax0] = self.maxMinVertAtt(0)
        #print "[attmin0, attmax0]" ,[attmin0, attmax0] 
        [attmin1, attmax1] = self.maxMinVertAtt(1)
        #print [xmin, ymin, xmax, ymax]
        xrange = xmax - xmin
        yrange = ymax - ymin
        if xrange > yrange:
            min,max = xmin, xmax
        else:
            min,max = ymin, ymax
            
        for obj in self.getUserVertices():
            obj.x = (obj.x - xmin)/(max- min)*scale + offset
            obj.y = (obj.y - ymin)/(max- min)*scale + offset
            if len(obj.attributes)  > 0 and attmin0 != attmax0:
                obj.attributes[0] = (obj.attributes[0]-attmin0)/ \
                                    (attmax0-attmin0)*height_scale
            if len(obj.attributes) > 1 and attmin1 != attmax1:
                obj.attributes[1] = (obj.attributes[1]-attmin1)/ \
                                    (attmax1-attmin1)*height_scale
            
        for obj in self.getMeshVertices():
            obj.x = (obj.x - xmin)/(max- min)*scale + offset
            obj.y = (obj.y - ymin)/(max- min)*scale + offset
            if len(obj.attributes)  > 0 and attmin0 != attmax0:
                obj.attributes[0] = (obj.attributes[0]-attmin0)/ \
                                    (attmax0-attmin0)*height_scale
            if len(obj.attributes) > 1 and attmin1 != attmax1:
                obj.attributes[1] = (obj.attributes[1]-attmin1)/ \
                                    (attmax1-attmin1)*height_scale
                
        for obj in self.getHoles():
            obj.x = (obj.x - xmin)/(max- min)*scale + offset
            obj.y = (obj.y - ymin)/(max- min)*scale + offset
        for obj in self.getRegions():
            obj.x = (obj.x - xmin)/(max- min)*scale + offset
            obj.y = (obj.y - ymin)/(max- min)*scale + offset
        [xmin, ymin, xmax, ymax] = self.boxsize()
        #print [xmin, ymin, xmax, ymax]
     
    def boxsizeVerts(self):
        """
        Returns a list of verts denoting a box or triangle that contains
        verts on the xmin, ymin, xmax and ymax axis.
        Structure: list of verts 
        """
       
        large = kinds.default_float_kind.MAX
        xmin= large
        xmax=-large
        ymin= large
        ymax=-large
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
        xmin= large
        xmax=-large
        ymin= large
        ymax=-large
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
        min= large
        max=-large
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
            xoffset= -(SCALE*xmax - WIDTH + OFFSET)
        if SCALE*ymin < OFFSET:
            b = abs(SCALE*ymin)+OFFSET
        if SCALE*ymax > HEIGHT-OFFSET:
            b = -(SCALE*ymax - HEIGHT + OFFSET)
        yoffset = HEIGHT - b
        return [SCALE, xoffset, yoffset]

         
    def exportASCIIobj(self,ofile):
        """
        export a file, ofile, with the format
         lines:  v <x> <y> <first attribute>
        f <vertex #>  <vertex #> <vertex #> (of the triangles)
        """
        fd = open(ofile,'w')
        self.writeASCIIobj(fd)   
        fd.close()


    def writeASCIIobj(self,fd):
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
           
    def exportASCIIsegmentoutlinefile(self,ofile):
        """Write the boundary user mesh info, eg
         vertices that are connected to segments,
         segments
        """
        
        verts = {}
        for seg in self.getUserSegments():
            verts[seg.vertices[0]] = seg.vertices[0]
            verts[seg.vertices[1]] = seg.vertices[1]
        meshDict = self.Mesh2IOOutlineDict(userVertices=verts.values())
        export_mesh_file(ofile,meshDict)
        
        # exportASCIImeshfile   - this function is used
    def export_mesh_file(self,ofile):
        """
        export a file, ofile, with the format
        """
        
        dict = self.Mesh2IODict()
        export_mesh_file(ofile,dict)

    # FIXME(DSG-DSG):Break this into two functions.
    #One for the outline points.
    #One for the mesh points.
    # Note: this function is not in the gui
    def exportPointsFile(self,ofile):
        """
        export a points file, ofile.
        
        """
        
        mesh_dict = self.Mesh2IODict()
        #point_dict = {}
        #point_dict['attributelist'] = {} #this will need to be expanded..
                                         # if attributes are brought back in. 
        #point_dict['geo_reference'] = self.geo_reference
        if mesh_dict['vertices'] == []:
            #point_dict['pointlist'] = mesh_dict['points']
            geo = Geospatial_data(mesh_dict['points'],
                                  geo_reference=self.geo_reference)
        else:
            #point_dict['pointlist'] = mesh_dict['vertices']
            geo = Geospatial_data(mesh_dict['vertices'],
                                  geo_reference=self.geo_reference)

        geo.export_points_file(ofile, absolute=True)
        


    def import_ungenerate_file(self,ofile, tag=None):
        """
        Imports an ungenerate file, from arcGIS into mesh.

        ofile is the name of the ungenerated file.
        Tag is a string name to be taggged on each segment. 

        WARNING values are assumed to be absolute.
        geo-refs are not taken into account..
        """
    
        dict = importUngenerateFile(ofile)
        default_tag = Segment.get_default_tag()
        if tag is not None:
            Segment.set_default_tag(str(tag))
        self.addVertsSegs(dict)
        Segment.set_default_tag(default_tag)

        # change the tag back to  it's default
    
        
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
        for element in dict_mesh.keys():
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
            #print "mesh meshDict['triangle_tags']", meshDict['triangle_tags']
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
        #print "mesh.Mesh2IOTriangulationDict*)*)"
        #print meshDict
        #print "mesh.Mesh2IOTriangulationDict*)*)"

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
        #print "userVertices",userVertices
        #print "userSegments",userSegments 
        pointlist=[]
        pointattributelist=[]
        index = 0
        for vertex in userVertices:
            vertex.index = index 
            pointlist.append([vertex.x,vertex.y])
            pointattributelist.append(vertex.attributes)
            
            index += 1
        meshDict['points'] = pointlist
        meshDict['point_attributes'] = pointattributelist

        segmentlist=[]
        segmenttaglist=[]
        for seg in userSegments:
            segmentlist.append([seg.vertices[0].index,seg.vertices[1].index])
            segmenttaglist.append(seg.tag)
        meshDict['outline_segments'] =segmentlist 
        meshDict['outline_segment_tags'] =segmenttaglist
        
        holelist=[]
        for hole in holes:
            holelist.append([hole.x,hole.y]) 
        meshDict['holes'] = holelist
        
        regionlist=[]
        regiontaglist = []
        regionmaxarealist = []
        for region in regions:
            regionlist.append([region.x,region.y])
            regiontaglist.append(region.getTag())
            
            if (region.getMaxArea() != None): 
                regionmaxarealist.append(region.getMaxArea())
            else:
                regionmaxarealist.append(NOMAXAREA)
        meshDict['regions'] = regionlist
        meshDict['region_tags'] = regiontaglist
        meshDict['region_max_areas'] = regionmaxarealist
        #print "*(*("
        #print meshDict
        #print meshDict['regionlist']
        #print "*(*("
        return meshDict

    def IOTriangulation2Mesh(self, genDict):
        """
        Set the mesh attributes given an tsh IO dictionary
        """
        #Clear the current generated mesh values
        self.tri_mesh = None
        
        self.tri_mesh = Rigid_triangulation(
            genDict['triangles']
            ,genDict['segments']
            ,genDict['vertices']
            ,genDict['triangle_tags']
            ,genDict['triangle_neighbors']
            ,genDict['segment_tags']
            ,genDict['vertex_attributes']
            ,genDict['vertex_attribute_titles']
            )
        self.attributeTitles = genDict['vertex_attribute_titles']
        self.maxVertexIndex = len(genDict['vertices'])
        #print "self.maxVertexIndex ", self.maxVertexIndex
        
    def IOOutline2Mesh(self, genDict):
        """
        Set the outline (user Mesh attributes) given a IO tsh dictionary
        
        mesh is an instance of a mesh object
        """
        #Clear the current user mesh values
        self.clearUserSegments()
        self.userVertices=[]
        self.Holes=[]
        self.Regions=[]

        #print "mesh.IOOutline2Mesh@#@#@#"
        #print "genDict",genDict
        #print "@#@#@#"
        
        #index = 0
        for point in genDict['points']:
            v=Vertex(point[0], point[1])
            #v.index = index
            #index +=1
            self.userVertices.append(v)

        #index = 0
        for seg,tag in map(None,genDict['outline_segments'],
                           genDict['outline_segment_tags']):

            segObject = Segment( self.userVertices[int(seg[0])],
                           self.userVertices[int(seg[1])], tag = tag )
            #segObject.index = index
            #index +=1
            self.userSegments.append(segObject)

# Remove the loading of attribute info.
# Have attribute info added using least_squares in pyvolution
#         index = 0
#         for att in genDict['point_attributes']:
#             if att == None:
#                 self.userVertices[index].setAttributes([])
#             else:
#                 self.userVertices[index].setAttributes(att)
#            index += 1
        
        #index = 0
        for point in genDict['holes']:
            h=Hole(point[0], point[1])
            #h.index = index
            #index +=1
            self.holes.append(h)

        #index = 0
        for reg,att,maxArea in map(None,
                                   genDict['regions'],
                                   genDict['region_tags'],
                                   genDict['region_max_areas']):
            if maxArea > 0:  # maybe I should ref NOMAXAREA? Prob' not though
                Object = Region( reg[0],
                                 reg[1],
                                 tag = att,
                                 maxArea = maxArea)
            else:
                Object = Region( reg[0],
                                 reg[1],
                                 tag = att)
                
            #Object.index = index
            #index +=1
            self.regions.append(Object)
  
############################################

        
    def refineSet(self,setName):
        Triangles = self.sets[self.setID[setName]]
        Refine(self,Triangles)

    def selectAllTriangles(self):
        A=[]
        A.extend(self.meshTriangles)
        if not('All' in self.setID.keys()):
            self.setID['All']=len(self.sets)
            self.sets.append(A)
        else:
            self.sets[self.setID['All']]=A
        return 'All'
        # and objectIDs


    def clearSelection(self):
        A = []
        if not('None' in self.setID.keys()):
            self.setID['None']=len(self.sets)
            self.sets.append(A)
        return 'None'

    def drawSet(self,canvas,setName,SCALE,colour=SET_COLOUR):
    #FIXME Draws over previous triangles - may bloat canvas
        Triangles = self.sets[self.setID[setName]]
        for triangle in Triangles:
            triangle.draw(canvas,1,
                          scale = SCALE,
                          colour = colour)
	    
    def undrawSet(self,canvas,setName,SCALE,colour='green'):
    #FIXME Draws over previous lines - may bloat canvas
        Triangles = self.sets[self.setID[setName]]
        for triangle in Triangles:
            triangle.draw(canvas,1,
                          scale = SCALE,
                          colour = colour)

    def weed(self,Vertices,Segments):
        #Depreciated
        #weed out existing duplicates
        print 'len(self.getUserSegments())'
        print len(self.getUserSegments())
        print 'len(self.getUserVertices())'
        print len(self.getUserVertices())

        point_keys = {}
        for vertex in Vertices:
            point = (vertex.x,vertex.y)
            point_keys[point]=vertex
        #inlined would looks very ugly

        line_keys = {}
        for segment in Segments:
            vertex1 = segment.vertices[0]
            vertex2 = segment.vertices[1]
            point1 = (vertex1.x,vertex1.y)
            point2 = (vertex2.x,vertex2.y)
            segment.vertices[0]=point_keys[point1]
            segment.vertices[1]=point_keys[point2]
            vertex1 = segment.vertices[0]
            vertex2 = segment.vertices[1]
            point1 = (vertex1.x,vertex1.y)
            point2 = (vertex2.x,vertex2.y)
            line1 = (point1,point2)
            line2 = (point2,point1)
            if not (line_keys.has_key(line1) \
                 or line_keys.has_key(line2)):
                 line_keys[line1]=segment
        Vertices=point_keys.values()
        Segments=line_keys.values()
        return Vertices,Segments

    def segs_to_dict(self,segments):
        dict={}
        for segment in segments:
            vertex1 = segment.vertices[0]
            vertex2 = segment.vertices[1]
            point1 = (vertex1.x,vertex1.y)
            point2 = (vertex2.x,vertex2.y)
            line = (point1,point2)
            dict[line]=segment
        return dict

    def seg2line(self,s):
        return ((s.vertices[0].x,s.vertices[0].y,)\
                (s.vertices[1].x,s.vertices[1].y))

    def line2seg(self,line,tag=None):
        point0 = self.point2ver(line[0])
        point1 = self.point2ver(line[1])
        return Segment(point0,point1,tag=tag)

    def ver2point(self,vertex):
        return (vertex.x,vertex.y)

    def point2ver(self,point):
        return Vertex(point[0],point[1])

    def smooth_polySet(self,min_radius=0.05):
        #for all pairs of connecting segments:
        #    propose a new segment that replaces the 2

        #    If the difference between the new segment
        #    and the old lines is small: replace the
        #    old lines.

        seg2line = self.seg2line
        ver2point= self.ver2point
        line2seg = self.line2seg
        point2ver= self.point2ver

        #create dictionaries of lines -> segments
        userSegments = self.segs_to_dict(self.userSegments)
        alphaSegments = self.segs_to_dict(self.alphaUserSegments)

        #lump user and alpha segments
        for key in alphaSegments.keys():
            userSegments[key]=alphaSegments[key]

        #point_keys = tuple -> vertex
        #userVertices = vertex -> [line,line] - lines from that node
        point_keys = {}
        userVertices={}
        for vertex in self.getUserVertices():
            point = ver2point(vertex)
            if not point_keys.has_key(point):
                point_keys[point]=vertex
                userVertices[vertex]=[]
        for key in userSegments.keys():
            line = key
            point_0 = key[0]
            point_1 = key[1]
            userVertices[point_keys[point_0]].append(line)
            userVertices[point_keys[point_1]].append(line)

        for point in point_keys.keys():
            try:
            #removed keys can cause keyerrors
                vertex = point_keys[point]
		lines = userVertices[vertex]
    
                #if there are 2 lines on the node
		if len(lines)==2:
		    line_0 = lines[0]
		    line_1 = lines[1]
    
                    #if the tags are the the same on the 2 lines
		    if userSegments[line_0].tag == userSegments[line_1].tag:
			tag = userSegments[line_0].tag 
    
                        #point_a is one of the next nodes, point_b is the other
			if point==line_0[0]:
			    point_a = line_0[1]
			if point==line_0[1]:
			    point_a = line_0[0]
			if point==line_1[0]:
			    point_b = line_1[1]
			if point==line_1[1]:
			    point_b = line_1[0]
    
    
                        #line_2 is proposed
			line_2 = (point_a,point_b)

			#calculate the area of the triangle between
			#the two existing segments and the proposed
			#new segment
			ax = point_a[0]
			ay = point_a[1]
			bx = point_b[0]
			by = point_b[1]
			cx = point[0]
			cy = point[1]
			area=abs((bx*ay-ax*by)+(cx*by-bx*cy)+(ax*cy-cx*ay))/2

			#calculate the perimeter
			len_a =  ((cx-bx)**2+(cy-by)**2)**0.5 
			len_b =  ((ax-cx)**2+(ay-cy)**2)**0.5 
			len_c =  ((bx-ax)**2+(by-ay)**2)**0.5 
			perimeter = len_a+len_b+len_c

			#calculate the radius
			r = area/(2*perimeter)

			#if the radius is small: then replace the existing
			#segments with the new one
			if r < min_radius:
			    if len_c < min_radius: append = False
			    else: append = True
			    #if the new seg is also time, don't add it
			    if append:
				segment = self.line2seg(line_2,tag=tag)

			    list_a=userVertices[point_keys[point_a]]
			    list_b=userVertices[point_keys[point_b]]

			    if line_0 in list_a:
				list_a.remove(line_0)
			    else:
				list_a.remove(line_1)

			    if line_0 in list_b:
				list_b.remove(line_0)
			    else:
				list_b.remove(line_1)

			    if append:
				list_a.append(line_2)
				list_b.append(line_2)
			    else:
				if len(list_a)==0:
				    userVertices.pop(point_keys[point_a])
				    point_keys.pop(point_a)
				if len(list_b)==0:
				    userVertices.pop(point_keys[point_b])
				    point_keys.pop(point_b)

			    userVertices.pop(point_keys[point])
			    point_keys.pop(point)
			    userSegments.pop(line_0)
			    userSegments.pop(line_1)

			    if append:
				userSegments[line_2]=segment
            except:
                pass

        #self.userVerticies = userVertices.keys()
        #self.userSegments = []
        #for key in userSegments.keys():
        #    self.userSegments.append(userSegments[key])
        #self.alphaUserSegments = []

        self.userVerticies = []
        self.userSegments = []
        self.alphaUserSegments = []

        return userVertices,userSegments,alphaSegments

    def triangles_to_polySet(self,setName):
        #self.smooth_polySet()

        seg2line = self.seg2line
        ver2point= self.ver2point
        line2seg = self.line2seg
        point2ver= self.point2ver

        from Numeric import array,allclose
        #turn the triangles into a set
        Triangles = self.sets[self.setID[setName]]
        Triangles_dict = {}
        for triangle in Triangles:
            Triangles_dict[triangle]=None
 

        #create a dict of points to vertexes (tuple -> object)
        #also create a set of vertexes (object -> True)
        point_keys = {}
        userVertices={}
        for vertex in self.getUserVertices():
            point = ver2point(vertex)
            if not point_keys.has_key(point):
                point_keys[point]=vertex
                userVertices[vertex]=True

        #create a dict of lines to segments (tuple -> object)
        userSegments = self.segs_to_dict(self.userSegments)
        #append the userlines in an affine linespace
        affine_lines = Affine_Linespace()
        for line in userSegments.keys():
            affine_lines.append(line)
        alphaSegments = self.segs_to_dict(self.alphaUserSegments)
        for line in alphaSegments.keys():
            affine_lines.append(line)
        
        for triangle in Triangles:
            for i in (0,1,2):
                #for every triangles neighbour:
                if not Triangles_dict.has_key(triangle.neighbors[i]):
                #if the neighbour is not in the set:
                    a = triangle.vertices[i-1]
                    b = triangle.vertices[i-2]
                    #Get possible matches:
                    point_a = ver2point(a)
                    point_b = ver2point(b)
                    midpoint = ((a.x+b.x)/2,(a.y+b.y)/2)
                    line = (point_a,point_b)
                    tag = None


                    #this bit checks for matching lines
                    possible_lines = affine_lines[line] 
                    possible_lines = unique(possible_lines)
                    found = 0                            
                    for user_line in possible_lines:
                        if self.point_on_line(midpoint,user_line):
                            found+=1
                            assert found<2
                            if userSegments.has_key(user_line):
                                parent_segment = userSegments.pop(user_line)
                            if alphaSegments.has_key(user_line):
                                parent_segment = alphaSegments.pop(user_line)
                            tag = parent_segment.tag
                            offspring = [line]
                            offspring.extend(self.subtract_line(user_line,
                                                                line))
                            affine_lines.remove(user_line)
                            for newline in offspring:
                                line_vertices = []
                                for point in newline:
                                    if point_keys.has_key(point):
                                        vert = point_keys[point]
                                    else:
                                        vert = Vertex(point[0],point[1])
                                        userVertices[vert]=True
                                        point_keys[point]=vert
                                    line_vertices.append(vert)
                                segment = Segment(line_vertices[0],
                                                  line_vertices[1],tag)
                                userSegments[newline]=segment
                                affine_lines.append(newline)
                            #break
                    assert found<2



                    #if no matching lines
                    if not found:
                        line_vertices = []
		        for point in line:
			    if point_keys.has_key(point):
			        vert = point_keys[point]
			    else:
				vert = Vertex(point[0],point[1])
				userVertices[vert]=True
				point_keys[point]=vert
			    line_vertices.append(vert)
			segment = Segment(line_vertices[0],
                                          line_vertices[1],tag)
                        userSegments[line]=segment
                        affine_lines.append(line)
        
        self.userVerticies = []
        self.userSegments = []
        self.alphaUserSegments = []

        return userVertices,userSegments,alphaSegments

    def subtract_line(self,parent,child):
        #Subtracts child from parent
        #Requires that the child is a 
        #subline of parent to work.

        from Numeric import allclose,dot,array
        A= parent[0]
        B= parent[1]
        a = child[0]
        b = child[1]

        A_array = array(parent[0])
        B_array = array(parent[1])
        a_array   = array(child[0])
        b_array   = array(child[1])

        assert not A == B
        assert not a == b

        answer = []

        #if the new line does not share a 
        #vertex with the old one
        if not (allclose(A_array,a_array)\
             or allclose(B_array,b_array)\
             or allclose(A_array,b_array)\
             or allclose(a_array,B_array)):
            if dot(A_array-a_array,A_array-a_array) \
            < dot(A_array-b_array,A_array-b_array):
                sibling1 = (A,a)
                sibling2 = (B,b)
                return [sibling1,sibling2]
            else:
                sibling1 = (A,b)
                sibling2 = (B,a)
                return [sibling1,sibling2]

        elif allclose(A_array,a_array):
            if allclose(B_array,b_array):
                return []
	    else:
                sibling = (b,B)
                return [sibling]
    	elif allclose(B_array,b_array):
            sibling = (a,A)
            return [sibling]

        elif allclose(A_array,b_array):
            if allclose(B,a):
                return []
	    else:
                sibling = (a,B)
                return [sibling]
        elif allclose(a_array,B_array):
            sibling = (b,A)
            return [sibling]

    def point_on_line(self,point,line):
        #returns true within a tolerance of 3 degrees
        x=point[0]
        y=point[1]
        x0=line[0][0]
        x1=line[1][0]
        y0=line[0][1]
        y1=line[1][1]
        from Numeric import array, dot, allclose
        from math import sqrt
        tol = 3. #DEGREES
        tol = tol*3.1415/180

        a = array([x - x0, y - y0]) 
        a_normal = array([a[1], -a[0]])      
        len_a_normal = sqrt(sum(a_normal**2)) 

        b = array([x1 - x0, y1 - y0])  	       
        len_b = sqrt(sum(b**2)) 
    
        if abs(dot(a_normal, b)/(len_b*len_a_normal))< tol:
            #Point is somewhere on the infinite extension of the line

            len_a = sqrt(sum(a**2))     	       	                
            if dot(a, b) >= 0 and len_a <= len_b:
               return True
            else:   
               return False
        else:
          return False

    def line_length(self,line):
        x0=line[0][0]
        x1=line[1][0]
        y0=line[0][1]
        y1=line[1][1]
        return ((x1-x0)**2-(y1-y0)**2)**0.5     

    def threshold(self,setName,min=None,max=None,attribute_name='elevation'):
        """
        threshold using  d
        """
        triangles = self.sets[self.setID[setName]]
        A = []

        if attribute_name in self.attributeTitles:
            i = self.attributeTitles.index(attribute_name)
        else: i = -1#no attribute
        if not max == None:
            for t in triangles:
                if (min<self.av_att(t,i)<max):
                    A.append(t)
        else:
            for t in triangles:
                if (min<self.av_att(t,i)):
                    A.append(t)
        self.sets[self.setID[setName]] = A

    def general_threshold(self,setName,min=None,max=None\
              ,attribute_name = 'elevation',function=None):
        """
        Thresholds the triangles 
        """
        from visual.graph import arange,ghistogram,color as colour
        triangles = self.sets[self.setID[setName]]
        A = []
        data=[]
        #data is for the graph

        if attribute_name in self.attributeTitles:
            i = self.attributeTitles.index(attribute_name)
        else: i = -1
        if not max == None:
            for t in triangles:
                value=function(t,i)
                if (min<value<max):
                    A.append(t)
                data.append(value)
        else:
            for t in triangles:
                value=function(t,i)
                if (min<value):
                    A.append(t)
                data.append(value)
        self.sets[self.setID[setName]] = A

        if self.visualise_graph:
            if len(data)>0:
                max=data[0]
                min=data[0]
                for value in data:
                    if value > max:
                        max = value
                    if value < min:
                        min = value

                inc = (max-min)/100

                histogram = ghistogram(bins=arange(min,max,inc),\
                             color = colour.red)
                histogram.plot(data=data)
        
    def av_att(self,triangle,i):
        if i==-1: return 1
        else:
            #evaluates the average attribute of the vertices of a triangle.
            V = triangle.getVertices()
            a0 = (V[0].attributes[i])
            a1 = (V[1].attributes[i])
            a2 = (V[2].attributes[i])
            return (a0+a1+a2)/3

    def Courant_ratio(self,triangle,index):
        """
        Uses the courant threshold
        """
        e = self.av_att(triangle,index)
        A = triangle.calcArea()
        P = triangle.calcP()
        r = A/(2*P)
        e = max(0.1,abs(e))
        return r/e**0.5

    def Gradient(self,triangle,index):
        V = triangle.vertices
        x0, y0, x1, y1, x2, y2, q0, q1, q2 = V[0].x,V[0].y,V[1].x,V[1].y,V[2].x,V[2].y,V[0].attributes[index],V[1].attributes[index],V[2].attributes[index]
        grad_x,grad_y = gradient(x0, y0, x1, y1, x2, y2, q0, q1, q2)
        if ((grad_x**2)+(grad_y**2))**(0.5)<0:
            print ((grad_x**2)+(grad_y**2))**(0.5)
        return ((grad_x**2)+(grad_y**2))**(0.5)
    

    def append_triangle(self,triangle):
        self.meshTriangles.append(triangle)

    def replace_triangle(self,triangle,replacement):
        i = self.meshTriangles.index(triangle)
        self.meshTriangles[i]=replacement
        assert replacement in self.meshTriangles

def importUngenerateFile(ofile):
    """
    import a file, ofile, with the format
    [poly]
    poly format:
    First line:  <# of vertices> <x centroid> <y centroid>
    Following lines: <x> <y> 
    last line:  "END"

    Note: These are clockwise.
    """
    fd = open(ofile,'r')
    Dict = readUngenerateFile(fd)
    fd.close()
    return Dict

def readUngenerateFile(fd):
    """
    import a file, ofile, with the format
    [poly]
    poly format:
    First line:  <# of polynomial> <x centroid> <y centroid>
    Following lines: <x> <y> 
    last line:  "END"
    """
    END_DELIMITER = 'END'
    
    points = []
    segments = []
    
    isEnd = False
    line = fd.readline() #not used <# of polynomial> <x> <y>
    while not isEnd:
        line = fd.readline()
        fragments = line.split()
        vert = [float(fragments.pop(0)),float(fragments.pop(0))]
        points.append(vert)
        PreviousVertIndex = len(points)-1
        firstVertIndex = PreviousVertIndex
        
        line = fd.readline() #Read the next line
        while not line.startswith(END_DELIMITER): 
            #print "line >" + line + "<"
            fragments = line.split()
            vert = [float(fragments.pop(0)),float(fragments.pop(0))]
            points.append(vert)
            thisVertIndex = len(points)-1
            segment = [PreviousVertIndex,thisVertIndex]
            segments.append(segment)
            PreviousVertIndex = thisVertIndex
            line = fd.readline() #Read the next line
            i =+ 1
        # If the last and first segments are the same,
        # Remove the last segment and the last vertex
        # then add a segment from the second last vert to the 1st vert
        thisVertIndex = len(points)-1
        firstVert = points[firstVertIndex]
        thisVert = points[thisVertIndex]
        #print "firstVert",firstVert
        #print "thisVert",thisVert
        if (firstVert[0] == thisVert[0] and firstVert[1] == thisVert[1]):
            points.pop()
            segments.pop()
            thisVertIndex = len(points)-1
            segments.append([thisVertIndex, firstVertIndex])
        
        line = fd.readline() # read <# of polynomial> <x> <y> OR END
        #print "line >>" + line + "<<"
        if line.startswith(END_DELIMITER):
            isEnd = True
    
    #print "points", points       
    #print "segments", segments
    ungenerated_dict = {}
    ungenerated_dict['points'] = points
    ungenerated_dict['segments'] = segments
    return ungenerated_dict

def importMeshFromFile(ofile):
    """returns a mesh object, made from a points file or .tsh/.msh file
    Often raises IOError,RuntimeError
    """
    newmesh = None
    if (ofile[-4:]== ".pts" or ofile[-4:]== ".txt" or \
        ofile[-4:]== ".csv"):
        geospatial = Geospatial_data(ofile)
        dict = {}
        dict['points'] = geospatial.get_data_points(absolute=False)
        dict['outline_segments'] = []
        dict['outline_segment_tags'] = []
        dict['regions'] = []
        dict['region_tags'] = []
        dict['region_max_areas'] = []
        dict['holes'] = [] 
        newmesh= Mesh(geo_reference = geospatial.geo_reference)
        newmesh.IOOutline2Mesh(dict)
        counter = newmesh.removeDuplicatedUserVertices()
        if (counter >0):
            print "%i duplicate vertices removed from dataset" % (counter)
    elif (ofile[-4:]== ".tsh" or ofile[-4:]== ".msh"):
        dict = import_mesh_file(ofile)
        #print "********"
        #print "zq mesh.dict",dict
        #print "********" 
        newmesh= Mesh()
        newmesh.IOOutline2Mesh(dict)
        newmesh.IOTriangulation2Mesh(dict)
    else:
        raise RuntimeError
    
    if dict.has_key('geo_reference') and not dict['geo_reference'] == None:
        newmesh.geo_reference = dict['geo_reference']
    return newmesh

def loadPickle(currentName):
    fd = open(currentName)
    mesh = pickle.load(fd)
    fd.close()
    return mesh
    
def square_outline(side_length = 1,up = "top", left = "left", right = "right",
                   down = "bottom", regions = False):
    
        a = Vertex (0,0)
        b = Vertex (0,side_length)
        c = Vertex (side_length,0)
        d = Vertex (side_length,side_length)
      
        s2 = Segment(b,d, tag = up)
        s3 = Segment(b,a, tag = left)
        s4 = Segment(d,c, tag = right)
        s5 = Segment(a,c, tag = down)

        if regions:
            e =  Vertex (side_length/2,side_length/2)
            s6 = Segment(a,e, tag = down + left)
            s7 = Segment(b,e, tag = up + left)
            s8 = Segment(c,e, tag = down + right)
            s9 = Segment(d,e, tag = up + right)
            r1 = Region(side_length/2,3.*side_length/4, tag = up)
            r2 = Region(1.*side_length/4,side_length/2, tag = left)
            r3 = Region(3.*side_length/4,side_length/2, tag = right)
            r4 = Region(side_length/2,1.*side_length/4, tag = down)
            mesh = Mesh(userVertices=[a,b,c,d,e],
                        userSegments=[s2,s3,s4,s5,s6,s7,s8,s9],
                        regions = [r1,r2,r3,r4])
        else:
            mesh = Mesh(userVertices=[a,b,c,d],
                 userSegments=[s2,s3,s4,s5])
     
        return mesh

    

def region_strings2ints(region_list):
    """Given a list of (x_int,y_int,tag_string) lists it returns a list of
    (x_int,y_int,tag_int) and a list to convert the tag_int's back to
    the tag_strings  
    """
    # Make sure "" has an index of 0 
    region_list.reverse()
    region_list.append((1.0,2.0,""))
    region_list.reverse()
    convertint2string = []
    for i in xrange(len(region_list)):
        convertint2string.append(region_list[i][2])
        if len(region_list[i]) == 4: # there's an area value
            region_list[i] = (region_list[i][0],
                              region_list[i][1],i,region_list[i][3])
        elif len(region_list[i]) == 3: # no area value
            region_list[i] = (region_list[i][0],region_list[i][1],i)
        else:
            print "The region list has a bad size"
            # raise an error ..
            raise Error

    #remove "" from the region_list
    region_list.pop(0)
    
    return [region_list, convertint2string]

def region_ints2strings(region_list,convertint2string):
    """Reverses the transformation of region_strings2ints
    """
    #print 'region_ints2strings region_list', region_list
    
    returned_region_list = []
    # may not need (not region_list[0] == [])
    # or region_list[0] == [0.0]
    if (not region_list[0] == []): # or region_list[0] == [0.0]:
        #print "in loop"
        for i in xrange(len(region_list)):
            temp = region_list[i]
            returned_region_list.append(convertint2string[int(temp[0])])
    return returned_region_list

def segment_ints2strings(intlist, convertint2string):
    """Reverses the transformation of segment_strings2ints """
    stringlist = []
    for x in intlist:
        stringlist.append(convertint2string[x])
    return stringlist

def segment_strings2ints(stringlist, preset):
    """Given a list of strings return a list of 0 to n ints which represent
    the strings and a converting list of the strings, indexed by 0 to n ints.
    Also, input an initial converting list of the strings
    Note, the converting list starts off with
    ["internal boundary", "external boundary", "internal boundary"]
    example input and output 
    input ["a","b","a","c"],["c"]
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
    #Remove the preset tags from the list with no duplicates
    nodups = [x for x in nodups if x not in preset]

    try:
        nodups.remove("") # this has to go to zero
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
        return u.keys()

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

"""Refines triangles

   Implements the #triangular bisection?# algorithm.
 
   
"""

def Refine(mesh, triangles):
    """
    Given a general_mesh, and a triangle number, split
    that triangle in the mesh in half. Then to prevent
    vertices and edges from meeting, keep refining 
    neighbouring triangles until the mesh is clean.
    """
    state = BisectionState(mesh)
    for triangle in triangles:
        if not state.refined_triangles.has_key(triangle):
            triangle.rotate_longest_side()
            state.start(triangle)
            Refine_mesh(mesh, state)

def Refine_mesh(mesh, state):
    """
    """
    state.getState(mesh)
    refine_triangle(mesh,state)
    state.evolve()
    if not state.end:
        Refine_mesh(mesh,state)

def refine_triangle(mesh,state):
    split(mesh,state.current_triangle,state.new_point)
    if state.case == 'one':
        state.r[3]=state.current_triangle#triangle 2

        new_triangle_id = len(mesh.meshTriangles)-1
        new_triangle = mesh.meshTriangles[new_triangle_id]

        split(mesh,new_triangle,state.old_point)
        state.r[2]=new_triangle#triangle 1.2
        state.r[4]=mesh.meshTriangles[len(mesh.meshTriangles)-1]#triangle 1.1
        r = state.r
        state.repairCaseOne()

    if state.case == 'two':
        state.r[2]=mesh.meshTriangles[len(mesh.meshTriangles)-1]#triangle 1

        new_triangle = state.current_triangle

        split(mesh,new_triangle,state.old_point)

        state.r[3]=mesh.meshTriangles[len(mesh.meshTriangles)-1]#triangle 2.1
        state.r[4]=new_triangle#triangle 2.2
        r = state.r

        state.repairCaseTwo()

    if state.case == 'vertex':
        state.r[2]=state.current_triangle#triangle 2
        state.r[3]=mesh.meshTriangles[len(mesh.meshTriangles)-1]#triangle 1
        r = state.r
        state.repairCaseVertex()
        
    if state.case == 'start':
        state.r[2]=mesh.meshTriangles[len(mesh.meshTriangles)-1]#triangle 1
        state.r[3]=state.current_triangle#triangle 2

    if state.next_case == 'boundary':
        state.repairCaseBoundary()


def split(mesh, triangle, new_point):
	"""
	Given a mesh, triangle_id and a new point,
	split the corrosponding triangle into two
	new triangles and update the mesh.
	"""

	new_triangle1 = Triangle(new_point,triangle.vertices[0],
                                 triangle.vertices[1],
                                 attribute = triangle.attribute,
                                 neighbors = None)
	new_triangle2 = Triangle(new_point,triangle.vertices[2],
                                 triangle.vertices[0],
                                 attribute = triangle.attribute,
                                 neighbors = None)

        new_triangle1.setNeighbors(triangle.neighbors[2],None,new_triangle2)
        new_triangle2.setNeighbors(triangle.neighbors[1],new_triangle1,None)

        mesh.meshTriangles.append(new_triangle1)

        triangle.vertices = new_triangle2.vertices
        triangle.neighbors = new_triangle2.neighbors


class State:

    def __init__(self):
        pass

class BisectionState(State):


    def __init__(self,mesh):
        self.len = len(mesh.meshTriangles)
        self.refined_triangles = {}
        self.mesh = mesh
        self.current_triangle = None
        self.case = 'start'
        self.end = False
        self.r = [None,None,None,None,None]

    def start(self, triangle):
        self.current_triangle = triangle
        self.case = 'start'
        self.end = False
        self.r = [None,None,None,None,None]

    def getState(self,mesh):
        if not self.case == 'vertex':
            self.new_point=self.getNewVertex(mesh, self.current_triangle)
            #self.neighbour=self.getNeighbour(mesh, self.current_triangle)
            self.neighbour = self.current_triangle.neighbors[0]
            if not self.neighbour is None:
                self.neighbour.rotate_longest_side()
            self.next_case = self.get_next_case(mesh,self.neighbour,
                                                self.current_triangle)
        if self.case == 'vertex':
            self.new_point=self.old_point


    def evolve(self):
        if self.case == 'vertex':
            self.end = True

        self.last_case = self.case
        self.case = self.next_case

        self.old_point = self.new_point
        self.current_triangle = self.neighbour

        if self.case == 'boundary':
            self.end = True
        self.refined_triangles[self.r[2]]=1
        self.refined_triangles[self.r[3]]=1
        if not self.r[4] is None:
            self.refined_triangles[self.r[4]]=1
        self.r[0]=self.r[2]
        self.r[1]=self.r[3]


    def getNewVertex(self,mesh,triangle):
	coordinate1 = triangle.vertices[1]
	coordinate2 = triangle.vertices[2]
	a = ([coordinate1.x*1.,coordinate1.y*1.])
	b = ([coordinate2.x*1.,coordinate2.y*1.])
        attributes = []
        for i in range(len(coordinate1.attributes)):
            att = (coordinate1.attributes[i]+coordinate2.attributes[i])/2
            attributes.append(att)
	new_coordinate = [((a[0]-b[0])/2+b[0]),((a[1]-b[1])/2+b[1])]
	newVertex = Vertex(new_coordinate[0],new_coordinate[1],
                           attributes = attributes)
        mesh.maxVertexIndex+=1
        newVertex.index = mesh.maxVertexIndex
        mesh.meshVertices.append(newVertex)
	return newVertex

    def get_next_case(self, mesh,neighbour,triangle):
	"""
	Given the locations of two neighbouring triangles, 
	examine the interior indices of their vertices (i.e. 
	0,1 or 2) to determine what how the neighbour needs
	to be refined.
	"""
	if (neighbour is None):
		next_case = 'boundary'
	else:
		if triangle.vertices[1].x==neighbour.vertices[2].x:
		    if triangle.vertices[1].y==neighbour.vertices[2].y:
			next_case = 'vertex'
                if triangle.vertices[1].x==neighbour.vertices[0].x:
		    if triangle.vertices[1].y==neighbour.vertices[0].y:
			next_case = 'two'
		if triangle.vertices[1].x==neighbour.vertices[1].x:
		    if triangle.vertices[1].y==neighbour.vertices[1].y:
			next_case = 'one'
	return next_case



    def repairCaseVertex(self):

        r = self.r


        self.link(r[0],r[2])
        self.repair(r[0])

        self.link(r[1],r[3])
        self.repair(r[1])

        self.repair(r[2])

        self.repair(r[3])


    def repairCaseOne(self):
        r = self.rkey


        self.link(r[0],r[2])
        self.repair(r[0])

        self.link(r[1],r[4])
        self.repair(r[1])

        self.repair(r[4])

    def repairCaseTwo(self):
        r = self.r

        self.link(r[0],r[4])
        self.repair(r[0])

        self.link(r[1],r[3])
        self.repair(r[1])

        self.repair(r[4])

    def repairCaseBoundary(self):
        r = self.r
        self.repair(r[2])
        self.repair(r[3])



    def repair(self,triangle):
        """
        Given a triangle that knows its neighbours, this will
        force the neighbours to comply.

        However, it needs to compare the vertices of triangles
        for this implementation 

        But it doesn't work for invalid neighbour structures
        """
        n=triangle.neighbors
        for i in (0,1,2):
            if not n[i] is None:
                for j in (0,1,2):#to find which side of the list is broken
                    if not (n[i].vertices[j] in triangle.vertices):
                    #ie if j is the side of n that needs fixing
                        k = j
                n[i].neighbors[k]=triangle



    def link(self,triangle1,triangle2):
        """
        make triangle1 neighbors point to t
                #count = 0riangle2
        """
        count = 0
        for i in (0,1,2):#to find which side of the list is broken
            if not (triangle1.vertices[i] in triangle2.vertices):
                j = i
                count+=1
        assert count == 1
        triangle1.neighbors[j]=triangle2

class Discretised_Tuple_Set:
    """
    if a={(0.0):[(0.01),(0.02)],(0.2):[(0.17)]}
    a[(0.01)]=a[(0.0)]=[(0.01),(0.02)]
    a[(10000)]=[] #NOT KEYERROR

    a.append[(0.01)]
    => {0.0:[(0.01),(0.02),(0.01)],0.2:[(0.17)]}

    #NOT IMPLIMENTED
    a.remove[(0.01)]
    => {(0.0):[(0.02),(0.01)],0.2:[(0.17)]}

    a.remove[(0.17)]
    => {(0.0):[(0.02),(0.01)],0.2:[]}
    #NOT IMPLIMENTED
    at a.precision = 2:
        a.round_up_rel[0.0]=
        a.round_flat[0.0]=
        a.round_down_rel[0.0]=

        a.up((0.1,2.04))=

    If t_rel = 0, nothing gets rounded into
    two bins. If t_rel = 0.5, everything does.

    Ideally, precision can be set high, so that multiple
    entries are rarely in the same bin. And t_rel should
    be low (<0.1 for 1 dimension!,<(0.1/n) for small n!!)
    so that it is rare to put items in mutiple bins.

    Ex bins per entry = product(a,b,c...,n)
    a = 1 or 2 s.t. Ex(a) = 1+2*t_rel
    b = 1 or 2 ... 

    BUT!!! to avoid missing the right bin:
    (-10)**(precision+1)*t_rel must be greater than the 
    greatest possible variation that an identical element
    can display.


    Note that if tol = 0.5 (the max allowed) 0.6 will round to .7 and .5
    but not .6 - this looks wrong, but note that *everything* will round,
    so .6 wont be missed as everything close to it will check in .7 and .5.
    """
    def __init__(self,p_rel = 6,t_rel = 0.01):
        self.__p_rel__ = p_rel
        self.__t_rel__ = t_rel

        self.__p_abs__ = p_rel+1
        self.__t_abs__ = t_rel

        assert t_rel <= 0.5
        self.__items__ = {}
        from math import frexp
        self.frexp = frexp
        roundings = [self.round_up_rel,\
        self.round_down_rel,self.round_flat_rel,\
        self.round_down_abs,self.round_up_abs,\
        self.round_flat_abs]#

        self.roundings = roundings

    def __repr__(self):
        return '%s'%self.__items__

    def rounded_keys(self,key):
        key = tuple(key)
        keys = [key]
        keys = self.__rounded_keys__(key)
        return (keys)

    def __rounded_keys__(self,key):
        keys = []
        rounded_key=list(key)
        rounded_values=list(key)

        roundings = list(self.roundings)

        #initialise rounded_values
        round = roundings.pop(0)
        for i in range(len(rounded_values)):
            rounded_key[i]=round(key[i])
            rounded_values[i]={}
            rounded_values[i][rounded_key[i]]=None
        keys.append(tuple(rounded_key))
        
        for round in roundings:
            for i in range(len(rounded_key)):
                rounded_value=round(key[i])
                if not rounded_values[i].has_key(rounded_value):
                    #ie unless round_up_rel = round_down_rel
                    #so the keys stay unique
                    for j in range(len(keys)):
                        rounded_key = list(keys[j])
                        rounded_key[i]=rounded_value
                        keys.append(tuple(rounded_key))
        return keys

    def append(self,item):
        keys = self.rounded_keys(item)
        for key in keys:
            if self.__items__.has_key(key):
                self.__items__[key].append(item)
            else:
                self.__items__[key]=[item]

    def __getitem__(self,key):
        answer = []
        keys = self.rounded_keys(key)
        for key in keys:
            if self.__items__.has_key(key):
                answer.extend(self.__items__[key])
        #if len(answer)==0:
        #    raise KeyError#FIXME or return KeyError
        #                  #FIXME or just return []?
        else:
            return answer #FIXME or unique(answer)?

    def __delete__(self,item):
        keys = self.rounded_keys(item)
        answer = False
        #if any of the possible keys contains
        #a list, return true
        for key in keys:        
            if self.__items__.has_key(key):
                if item in self.__items__[key]:
                    self.__items__[key].remove(item)

    def remove(self,item):
        self.__delete__(item)

    def __contains__(self,item):

        keys = self.rounded_keys(item)
        answer = False
        #if any of the possible keys contains
        #a list, return true
        for key in keys:        
            if self.__items__.has_key(key):
                if item in self.__items__[key]:
                    answer = True
        return answer


    def has_item(self,item):
        return self.__contains__(item)

    def round_up_rel2(self,value):
         t_rel=self.__t_rel__
         #Rounding up the value
         m,e = self.frexp(value)
         m = m/2
         e = e + 1
         #m is the mantissa, e the exponent
         # 0.5 < |m| < 1.0
         m = m+t_rel*(10**-(self.__p_rel__))
         #bump m up
         m = round(m,self.__p_rel__)
         return m*(2.**e)

    def round_down_rel2(self,value):
         t_rel=self.__t_rel__
         #Rounding down the value
         m,e = self.frexp(value)
         m = m/2
         e = e + 1
         #m is the mantissa, e the exponent
         # 0.5 < m < 1.0
         m = m-t_rel*(10**-(self.__p_rel__))
         #bump the |m| down, by 5% or whatever
         #self.p_rel dictates
         m = round(m,self.__p_rel__)
         return m*(2.**e)

    def round_flat_rel2(self,value):
    #redundant
         m,e = self.frexp(value)
         m = m/2
         e = e + 1
         m = round(m,self.__p_rel__)
         return m*(2.**e)

    def round_up_rel(self,value):
         t_rel=self.__t_rel__
         #Rounding up the value
         m,e = self.frexp(value)
         #m is the mantissa, e the exponent
         # 0.5 < |m| < 1.0
         m = m+t_rel*(10**-(self.__p_rel__))
         #bump m up
         m = round(m,self.__p_rel__)
         return m*(2.**e)

    def round_down_rel(self,value):
         t_rel=self.__t_rel__
         #Rounding down the value
         m,e = self.frexp(value)
         #m is the mantissa, e the exponent
         # 0.5 < m < 1.0
         m = m-t_rel*(10**-(self.__p_rel__))
         #bump the |m| down, by 5% or whatever
         #self.p_rel dictates
         m = round(m,self.__p_rel__)
         return m*(2.**e)

    def round_flat_rel(self,value):
    #redundant
         m,e = self.frexp(value)
         m = round(m,self.__p_rel__)
         return m*(2.**e)

    def round_up_abs(self,value):
         t_abs=self.__t_abs__
         #Rounding up the value
         m = value+t_abs*(10**-(self.__p_abs__))
         #bump m up
         m = round(m,self.__p_abs__)
         return m

    def round_down_abs(self,value):
         t_abs=self.__t_abs__
         #Rounding down the value
         m = value-t_abs*(10**-(self.__p_abs__))
         #bump the |m| down, by 5% or whatever
         #self.p_rel dictates
         m = round(m,self.__p_abs__)
         return m

    def round_flat_abs(self,value):
    #redundant?
         m = round(value,self.__p_abs__)
         return m

    def keys(self):
        return self.__items__.keys()


class Mapped_Discretised_Tuple_Set(Discretised_Tuple_Set):
    """
    This is a discretised tuple set, but 
    based on a mapping. The mapping MUST
    return a sequence.

    example: 
    def weight(animal):
        return [animal.weight]
    
    a = Mapped_Discretised_Tuple_Set(weight)
    a.append[cow]
    a.append[fox]
    a.append[horse]

    a[horse]    -> [cow,horse]
    a[dog]      -> [fox]
    a[elephant] -> []
    """
    def __init__(self,mapping,p_rel = 6, t_rel=0.01):
        Discretised_Tuple_Set.__init__\
         (self,p_rel,t_rel = t_rel)
        self.mapping = mapping

    def rounded_keys(self,key):
        mapped_key = tuple(self.mapping(key))
        keys = self.__rounded_keys__(mapped_key)
        return keys

class Affine_Linespace(Mapped_Discretised_Tuple_Set):
    """
    The affine linespace creates a way to record and compare lines.
    Precision is a bit of a hack, but it creates a way to avoid 
    misses caused by round offs (between two lines of different
    lenghts, the short one gets rounded off more). 
    I am starting to think that a quadratic search would be faster.
    Nearly.
    """
    def __init__(self,p_rel=4,t_rel=0.2):
        Mapped_Discretised_Tuple_Set.__init__\
            (self,self.affine_line,\
            p_rel=p_rel,t_rel=t_rel)

        roundings = \
        [self.round_down_rel,self.round_up_rel,self.round_flat_rel]
        self.roundings = roundings
        #roundings = \
        #[self.round_down_abs,self.round_up_abs,self.round_flat_abs]
        #self.roundings = roundings

    def affine_line(self,line):
        point_1 = line[0]
        point_2 = line[1]
        #returns the equation of a line
        #between two points, in the from
        #(a,b,-c), as in ax+by-c=0
        #or line *dot* (x,y,1) = (0,0,0)

        #Note that it normalises the line
        #(a,b,-c) so Line*Line = 1.
        #This does not change the mathematical
        #properties, but it makes comparism
        #easier.

        #There are probably better algorithms.
        x1 = point_1[0]
        x2 = point_2[0]
        y1 = point_1[1]
        y2 = point_2[1]
        dif_x = x1-x2
        dif_y = y1-y2

        if dif_x == dif_y == 0:
            msg = 'points are the same'
            raise msg
        elif abs(dif_x)>=abs(dif_y):
            alpha = (-dif_y)/dif_x
            #a = alpha * b
            b = -1.
            c = (x1*alpha+x2*alpha+y1+y2)/2.
            a = alpha*b
        else:
            beta = dif_x/(-dif_y)
            #b = beta * a
            a = 1.
            c = (x1+x2+y1*beta+y2*beta)/2.
            b = beta*a
        mag = abs(a)+abs(b)
        #This does not change the mathematical
        #properties, but it makes comparism possible.

        #note that the gradient is b/a, or (a/b)**-1.
        #so 

        #if a == 0:
        #    sign_a = 1.
        #else:
        #    sign_a = a/((a**2)**0.5)
        #if b == 0:
        #    sign_b = 1.
        #else:
        #    sign_b = b/((b**2)**0.5)
        #if c == 0:
        #    sign_c = 1.
        #else:
        #    sign_c = c/((c**2)**0.5)
        #a = a/mag*sign_a
        #b = b/mag*sign_b
        #c = c/mag*sign_c
        a = a/mag
        b = b/mag
        c = c/mag
        return a,b,c


# FIXME (DSG-DSG)
# To do-
# Create a clear interface. eg
# have the interface methods more at the top of this file and add comments
# for the interface functions/methods, use function_name (not functionName),

#Currently
#function_name methods assume absolute values.  Geo-refs can be passed in.
#

# instead of functionName
if __name__ == "__main__":
    pass

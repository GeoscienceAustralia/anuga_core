#!/usr/bin/env python
#

from builtins import zip
from builtins import map
from builtins import str
from builtins import range
import tempfile
import unittest

#try:
from anuga.pmesh.mesh import *
#except ImportError:  
#    from mesh import *


from anuga.load_mesh.loadASCII import *
from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.geospatial_data.geospatial_data import Geospatial_data
from anuga.geometry.polygon import  is_inside_polygon ### inside_polygon

import numpy as num

class meshTestCase(unittest.TestCase):
    def setUp(self):
        pass
    
    def tearDown(self):
        pass

    def testPointDistance(self):
        a = Point(0.0, 0.0)
        b = Point(0.0, 10.0)
        
        self.assertTrue( a.DistanceToPoint(b) == 10.0,
                        'Point DistanceToPoint is wrong!')
    
    def testVertexDistance(self):
        a = Vertex (0.0, 0.0)
        b = Vertex (0.0, 10.0)
        
        self.assertTrue( a.DistanceToPoint(b) == 10.0,
                        'Point DistanceToPoint is wrong!')
        
        
    def testSegment(self):
        a = Vertex (0.0, 0.0)
        b = Vertex (0.0, 10.0)
        s = Segment(a,b, tag = 20)     
        
        self.assertTrue( s.vertices[0].DistanceToPoint(s.vertices[1]) == 10.0,
                        'vertices in a segment are wrong')
        
        self.assertTrue( s.tag == 20.0,
                        'tag in a segment are wrong')

    def testdeleteUserVertex(self):

        
        mesh = Mesh()
        a = mesh.addUserVertex(0.0, 0.0)
        b = mesh.addUserVertex (0.0, 2.0)
        c = mesh.addUserVertex (2.0,0.0)
        
        s1 = mesh.addUserSegment(a,b)
        s2 = mesh.addUserSegment(a,c)
        s3 = mesh.addUserSegment(c,b)

        mesh.deleteMeshObject (a) 
        self.assertTrue(mesh.userSegments[0] == s3,
                        'Bad segment. ')        
        self.assertTrue(len(mesh.userSegments) == 1,
                        'Segments not deleted.')
        self.assertTrue(len(mesh.userVertices) == 2,
                        'Vertex not deleted.')
       
  
    # FIXME add test for minAngle    
    def testgenerateMesh(self):
        a = Vertex (0.0, 0.0)
        d = Vertex (0.0, 4.0)
        f = Vertex (4.0,0.0)

        s1 = Segment(a,d)
        s2 = Segment(d,f)
        s3 = Segment(a,f)

        r1 = Region(0.3, 0.3, tag = 1.3, maxArea = .6)
        #print r1
        m = Mesh(userVertices=[a,d,f], userSegments=[s1,s2,s3], regions=[r1] )
        
        m.generateMesh("Q", maxArea = 2.1 )         

        #print m

        #m.plotMeshTriangle()

        result = 1.414214
        delta  = 0.00001
        tri = m.getTriangulation()
        verts = m.getMeshVertices()
        x = verts[tri[1][0]][0]
        #self.assertTrue((m.meshTriangles[1].vertices[0].x < result + delta) or
         #               (m.meshTriangles[1].vertices[0].x > result - delta),
          #              'generated mesh is wrong!')

        self.assertTrue((x < result + delta) or
                        (x > result - delta),
                        'generated mesh is wrong!')
        
    def test_regionalMaxArea(self):
        v0 = Vertex (0.0, 0.0)
        v1 = Vertex (6.0, 0.0)
        v2 = Vertex (6.0,6.0)
        v3 = Vertex (0.0,6.0)

        s1 = Segment(v0,v1)
        s2 = Segment(v1,v2)
        s3 = Segment(v3,v2)
        s4 = Segment(v3,v0)
        s5 = Segment(v2,v0)

        r1 = Region(3, 1,tag = 1.3)
        #print r1
        m = Mesh(userVertices=[v0,v1,v2,v3], userSegments=[s1,s2,s3,s4,s5], regions=[r1] )
        
        m.generateMesh("Q", maxArea = 36 )         

        #m.plotMeshTriangle()
        #print "len(m.meshTriangles)",len(m.meshTriangles)

        self.assertTrue(len(m.getTriangulation()) == 2, 
                        'test_regionalMaxArea 1:generated mesh is wrong!')
        
        ## Another test case
        r1 = Region(3, 1,tag = 1.3)
        r2 = Region(1, 3,tag = 1.3, maxArea = 8)
        m = Mesh(userVertices=[v0,v1,v2,v3], userSegments=[s1,s2,s3,s4,s5],
                 regions=[r1,r2] )
        m.generateMesh("Q", maxArea = 36 )
        
        self.assertTrue(len(m.getTriangulation()) >= 6,
                        'testregion_with_maxarea 2: # of tris is wrong!')    
       
               
        ## Another test case
        r1 = Region(3, 1, tag = 1.3, maxArea = 8)
        r2 = Region(1, 3, tag = 1.3, maxArea = 8)
        m = Mesh(userVertices=[v0,v1,v2,v3], userSegments=[s1,s2,s3,s4,s5],
                 regions=[r1,r2] )
        m.generateMesh("Q", maxArea = 36 )  
        #print "len(m.meshTriangles)",len(m.meshTriangles)
        
        self.assertTrue(len(m.getTriangulation()) >= 8,
                        'testregion_with_maxarea 3: # of tris is wrong!')
                
                
        ## Another test case
        r1 = Region(3, 1, tag = 1.3 )
        r2 = Region(1, 3, tag = 1.3, maxArea = 8)
        m = Mesh(userVertices=[v0,v1,v2,v3], userSegments=[s1,s2,s3,s4,s5],
                 regions=[r1,r2] )
        m.generateMesh("Q", maxArea = 8 ) 
        self.assertTrue(len(m.getTriangulation()) >= 8,
                        'testregion_with_maxarea 4: # of tris is wrong!')    

        
        ## Another test case
        r1 = Region(3, 1,tag = 1.3, maxArea = 8)
        r2 = Region(1, 3,tag = 1.3, maxArea = 8)
        m = Mesh(userVertices=[v0,v1,v2,v3], userSegments=[s1,s2,s3,s4,s5],
                 regions=[r1,r2] )
        m.generateMesh("Q", maxArea = 36,isRegionalMaxAreas = False )      
        self.assertTrue(len(m.getTriangulation()) == 2, 
                        'test_regionalMaxArea 5:generated mesh is wrong!')
        
        ## Another test case
        r1 = Region(3, 1,tag = 1.3, maxArea = 8)
        r2 = Region(1, 3,tag = 1.3, maxArea = 8)
        m = Mesh(userVertices=[v0,v1,v2,v3], userSegments=[s1,s2,s3,s4,s5],
                 regions=[r1,r2] )
        m.generateMesh("Q",isRegionalMaxAreas = False )
        self.assertTrue(len(m.getTriangulation()) == 2, 
                        'test_regionalMaxArea 5:generated mesh is wrong!')
        
    def test_generate_mesh(self):
        v0 = Vertex (0.0, 0.0)
        v1 = Vertex (6.0, 0.0)
        v2 = Vertex (6.0,6.0)
        v3 = Vertex (0.0,6.0)

        s1 = Segment(v0,v1)
        s2 = Segment(v1,v2)
        s3 = Segment(v3,v2)
        s4 = Segment(v3,v0)
        s5 = Segment(v2,v0)

        r1 = Region(3, 1,tag = 1.3)
        #print r1
        m = Mesh(userVertices=[v0,v1,v2,v3], userSegments=[s1,s2,s3,s4,s5],
                 regions=[r1] )
        
        m.generate_mesh(maximum_triangle_area=36,verbose=False)         

        self.assertTrue(len(m.getTriangulation()) == 2, 
                        'test_regionalMaxArea 1:generated mesh is wrong!')
        
        ## Another test case
        r1 = Region(3, 1,tag = 1.3)
        r2 = Region(1, 3,tag = 1.3, maxArea = 8)
        m = Mesh(userVertices=[v0,v1,v2,v3], userSegments=[s1,s2,s3,s4,s5],
                 regions=[r1,r2] )
        m.generate_mesh(maximum_triangle_area=36,verbose=False)  
        
        self.assertTrue(len(m.getTriangulation()) >= 6,
                        'testregion_with_maxarea 2: # of tris is wrong!')    
               
        ## Another test case
        r1 = Region(3, 1, tag = 1.3, maxArea = 8)
        r2 = Region(1, 3, tag = 1.3, maxArea = 8)
        m = Mesh(userVertices=[v0,v1,v2,v3], userSegments=[s1,s2,s3,s4,s5],
                 regions=[r1,r2] )
        m.generate_mesh(maximum_triangle_area=36,verbose=False)         
        #print "len(m.getTriangulation())",len(m.getTriangulation())
        
        self.assertTrue(len(m.getTriangulation()) >= 8,
                        'testregion_with_maxarea 3: # of tris is wrong!')
                         
        ## Another test case
        r1 = Region(3, 1, tag = 1.3 )
        r2 = Region(1, 3, tag = 1.3, maxArea = 8)
        m = Mesh(userVertices=[v0,v1,v2,v3], userSegments=[s1,s2,s3,s4,s5],
                 regions=[r1,r2] )
        m.generate_mesh(maximum_triangle_area=8,verbose=False)    
        self.assertTrue(len(m.getTriangulation()) >= 8,
                        'testregion_with_maxarea 4: # of tris is wrong!')    

        ## Another test case r1 = Region(3, 1,tag = 1.3, maxArea = 8)
        r2 = Region(1, 3,tag = 1.3, maxArea = 8)
        m = Mesh(userVertices=[v0,v1,v2,v3],
        userSegments=[s1,s2,s3,s4,s5], regions=[r1,r2] )
        m.generate_mesh(verbose=False)
        #print "en(m.getTriangulation())", len(m.getTriangulation())
        self.assertTrue(len(m.getTriangulation()) >= 8,
        'You have issues!')
        
    def testdeleteUserVertex(self):
        mesh = Mesh()
        a = mesh.addUserVertex(0.0, 0.0)
        b = mesh.addUserVertex (0.0, 2.0)
        c = mesh.addUserVertex (2.0,0.0)
        
        s1 = mesh.addUserSegment(a,b)
        s2 = mesh.addUserSegment(a,c)
        s3 = mesh.addUserSegment(c,b)

        mesh.deleteMeshObject (s2)
        
        #print ",s2 in mesh.userSegments" ,s2 in mesh.userSegments
        self.assertTrue(not(s2 in mesh.userSegments),
                        'Bad segment. ')        
        self.assertTrue(len(mesh.userSegments) ==2,
                        'Segments not deleted.')
        self.assertTrue(len(mesh.userVertices) == 3,
                        'Vertex deleted, instead of segment.')

    def testisUserSegmentNew (self):
        mesh = Mesh()
        a = mesh.addUserVertex(0.0, 0.0)
        b = mesh.addUserVertex (0.0, 2.0)
        c = mesh.addUserVertex (2.0,0.0)
        d = mesh.addUserVertex (2.0,3.0)
        
        s1 = mesh.addUserSegment(a,b)
        s2 = mesh.addUserSegment(a,c)
        s3 = mesh.addUserSegment(c,b)

        self.assertTrue(mesh.isUserSegmentNew(a,d) ,
                        'Segment should be new. ')
        self.assertTrue(not(mesh.isUserSegmentNew(a,b)) ,
                        'Segment should not be new. ')


    def testisUserSegmentNewII (self):
        mesh = Mesh()
        a = mesh.addUserVertex(0.0, 0.0)
        b = mesh.addUserVertex (0.0, 2.0)
        c = mesh.addUserVertex (2.0,0.0)
        d = mesh.addUserVertex (2.0,3.0)
        
        s1 = mesh.addUserSegment(a,b)
        s2 = mesh.addUserSegment(a,c)
        s3 = mesh.addUserSegment(c,b)

        self.assertTrue(mesh.representedUserSegment(a,d) is None,
                        'Segment should be new. ')
        self.assertTrue(mesh.representedUserSegment(a,b) == s1 ,
                        'Segment should not be new. ')
        
    def testauto_segment(self):
        p0 = Vertex (0.0, 0.0)
        p1 = Vertex (0.0, 4.0)
        p2 = Vertex (4.0,4.0)
        p3 = Vertex (4.0,0.0)

        s1 = Segment(p0,p1)
        
        m = Mesh(userVertices=[p0, p1, p2, p3], userSegments=[s1] ) 
        m.auto_segment()
        
        #print 'Len', len(m.userSegments)
        self.assertTrue(len(m.getUserSegments()) == 4 ,
                        'userSegments is wrong!')
     
        m.auto_segment()
        self.assertTrue(len(m.getUserSegments()) == 4 ,
                        'userSegments is wrong!')
     
    def testauto_segmentII(self):
        p1 = Vertex (3.0, 4.0)
        p2 = Vertex (3.0,2.0)
        p3 = Vertex (3.0,0.0)
        p4 = Vertex (6.0, 4.0)
        p5 = Vertex (6.0,2.0)
        p0 = Vertex (6.0,0.0)


        s1 = Segment(p2,p3)
        s2 = Segment(p4,p5)
        
        m = Mesh(userVertices=[p0, p1, p2, p3, p4, p5],
                 userSegments=[s1, s2])     

        m.auto_segment()
        
        s3 = m.representedAlphaUserSegment(p3,p0)
        self.assertTrue(not (s3 is None) ,
                        'userSegments is wrong!')

        
        s6 = m.representedAlphaUserSegment(p1,p4)       
        self.assertTrue(not (s6 is None) ,
                        'userSegments is wrong!')
        
        # remove a segment, add a point, auto_segment
        m.alphaUserSegments.remove(s3)
        p6 = Vertex (1.0, 2.0)
        m.userVertices.append(p6)
        
        m.auto_segment()
        
        s1_now = m.representedUserSegment(p3,p2)
        self.assertTrue(s1_now == s1 ,
                        'userSegments is wrong!')
        
        s2_now = m.representedUserSegment(p5,p4)       
        self.assertTrue(s2_now == s2 ,
                        'userSegments is wrong!')
        
        s3 = m.representedAlphaUserSegment(p3,p6)       
        self.assertTrue(not (s3 is None) ,
                        'userSegments is wrong!')
        
        s4 = m.representedAlphaUserSegment(p3,p6)       
        self.assertTrue(not (s4 is None) ,
                        'userSegments is wrong!')
        
        s5 = m.representedAlphaUserSegment(p4,p6)       
        self.assertTrue(s5 is None ,
                        'userSegments is wrong!')
        #print m
        
    def testRegions(self):
        a = Vertex (0.0, 0.0)
        b = Vertex (0.0, 2.0)
        c = Vertex (2.0,0.0)
        d = Vertex (0.0, 4.0)
        e = Vertex (2.0, 2.0)
        f = Vertex (4.0,0.0)
        g = Vertex (0.0,-2.0)
        
        Segment.set_default_tag("")
        s1 = Segment(a,b)
        s2 = Segment(b,e)
        s3 = Segment(e,c)
        s4 = Segment(c,a)
        s5 = Segment(b,d)
        s6 = Segment(e,d)
        s7 = Segment(e,f)
        s8 = Segment(c,f)
        s9 = Segment(g,c)
        s10 = Segment(g,a)

        r1 = Region(0.1,0.1,tag="22")
        r2 = Region(0.1,2.1,tag="11")
        r3 = Region(2.1,0.1)
        
        m = Mesh(userVertices=[a,b,c,d,e,f,g], userSegments=[s1,s2,s3,s4,s5,s6,s7,s8,s9,s10], regions=[r1,r2,r3] )
        m.generateMesh("Q", maxArea = 2.1 )

        # FIXME test the region
        #Triangulation =  m.getTriangulation()
        Triangulation = m.tri_mesh.triangle_tags
        #print Triangulation[0].attribute
        #print Triangulation[1].attribute 
        #print Triangulation[2].attribute 
        #print Triangulation[3].attribute 
        #print Triangulation[4].attribute
       
        self.assertTrue(Triangulation[0] == "" and
                        Triangulation[1] == "22" and
                        Triangulation[2] == "" and
                        Triangulation[3] == "11" and
                        Triangulation[4] == "22" ,
                        'region attributes are wrong!')   

    def test_vertexAttribs(self):
        a = Vertex (0.0, 0.0, attributes = [12.0,2.0])
        d = Vertex (0.0, 4.0, attributes = [9.0,7.0])
        f = Vertex (4.0,0.0, attributes = [14.0,3.0])
    
        Segment.set_default_tag("")
        s1 = Segment(a,d)
        s2 = Segment(d,f)
        s3 = Segment(a,f)
     
        r1 = Region(0.3, 0.3, tag = 88.9)
     
        m = Mesh(userVertices=[a,d,f], userSegments=[s1,s2,s3], regions=[r1])

        m.generateMesh("Q", maxArea = 2.1)

        vert = m.getMeshVerticeAttributes()
        
        self.assertTrue(num.all(vert[0] == [12.0, 2.0]) and
                        num.all(vert[1] == [9.0, 7.0]) and
                        num.all(vert[2] == [14.0,3.0]) and
                        num.all(vert[3] == [12.232233047033631,
                                            4.4142135623730949]) and
                        num.all(vert[4] == [13.0, 2.5]) ,
                        'vertex attributes are wrong!')

        
    def test_vertexAttribs2(self):
    
        a = Vertex (0.0, 0.0)
        d = Vertex (0.0, 4.0)
        f = Vertex (4.0,0.0)
    
        Segment.set_default_tag("")
        s1 = Segment(a,d)
        s2 = Segment(d,f)
        s3 = Segment(a,f)
     
        r1 = Region(0.3, 0.3, tag = 88.9)
     
        m = Mesh(userVertices=[a,d,f], userSegments=[s1,s2,s3], regions=[r1])

        m.generateMesh("Q", maxArea = 2.1 )

        vert = m.getMeshVerticeAttributes()
        #print "vert", vert
        self.assertTrue(vert == [],
                        'vertex attributes are wrong!')

    def test_segtag(self):
    
        a = Vertex (0.0, 0.0)
        d = Vertex (0.0, 4.0)
        f = Vertex (4.0,0.0)
    
        s1 = Segment(a,d,tag = 5)
        s2 = Segment(d,f,tag = 7)
        s3 = Segment(a,f,tag = 9)
     
        m = Mesh(userVertices=[a,d,f], userSegments=[s1,s2,s3])

        m.generateMesh("Q", maxArea = 2.1 )

        #m.export_mesh_file("from_test_mesh.tsh")
        seg = m.getMeshSegmentTags()
        #print "seg",seg
        #print "seg[0].tag"
        #print seg[0].tag
        #print "seg[0].tag"
        
        self.assertTrue(seg[0] == 5 and
                        seg[1] == 7 and
                        seg[2] == 9 and
                        seg[3] == 7 and
                        seg[4] == 9,
                        'seg tags are wrong')
            

    def test_segtag2(self):
    
        a = Vertex (0.0, 0.0)
        d = Vertex (0.0, 4.0)
        f = Vertex (4.0,0.0)
        e = Vertex (1.0,1.0)
    
        s1 = Segment(a,d)
        s2 = Segment(d,f)
        s3 = Segment(a,f)
        s4 = Segment(a,e)
     
        m = Mesh(userVertices=[a,d,f,e], userSegments=[s1,s2,s3,s4])

        m.generateMesh("Q", maxArea = 2.1)

        seg = m.getMeshSegmentTags()
        self.assertTrue(seg[0] == "exterior" and
                        seg[1] == "exterior" and
                        seg[2] == "exterior" and
                        seg[3] == "" and
                        seg[4] == "exterior",
                        '2nd seg tags are wrong')

    def test_asciiFile(self):
    
        a = Vertex (0.0, 0.0)  #, attributes = [1.1])
        d = Vertex (0.0, 4.0)  #, attributes = [1.2])
        f = Vertex (4.0,0.0)  #, attributes = [1.3])
        e = Vertex (1.0,1.0)  #, attributes = [1.4])
    
        s1 = Segment(a,d)
        s2 = Segment(d,f)
        s3 = Segment(a,f)
        s4 = Segment(a,e)
     
        m = Mesh(userVertices=[a,d,f,e], userSegments=[s1,s2,s3,s4])

        m.generateMesh("Q", maxArea = 2.1 )
        seg = m.getMeshSegments()
        
        fileName = tempfile.mktemp(".tsh")
        m.export_mesh_file(fileName)
        file = open(fileName)
        lFile = file.read().split('\n')
        file.close()
        os.remove(fileName)
        
        #print "@^@^"
        #for l in lFile:
        #    print l,"<"
        #print "@^@^"

        # no need to check the title again
        #self.assertTrue(lFile[0] == "5 0 # <vertex #> <x> <y> [attributes] ...Triangulation Vertices..."
          #              ,'Ascii file is wrong, vertex title')
        self.assertTrue(lFile[1] == "0 0.0 0.0 " and #1.1 " and
                        lFile[2] == "1 0.0 4.0 " and #1.2 " and
                        lFile[3] == "2 4.0 0.0 " and #1.3 " and
                        lFile[4] == "3 1.0 1.0 " and #1.4 " and
                        lFile[5] == "4 2.0 2.0 "  #1.25 " 
                        ,
                        'Ascii file is wrong, vertex')
        
        #self.assertTrue(lFile[6] == "# attribute column titles ...Triangulation Vertex Titles..."
          #              ,'Ascii file is wrong, attribute column title')
        self.assertTrue(lFile[8] == "0 3 2 4 -1 2 3  " and
                        lFile[9] == "1 1 0 3 3 2 -1  " and
                        lFile[10] == "2 3 4 1 -1 1 0  " and
                        lFile[11] == "3 2 3 0 1 -1 0  "
                        ,
                        'Ascii file is wrong, triangle') 

        self.assertTrue( lFile[13] == "0 0 1 exterior" and
                        lFile[14] == "1 1 4 exterior" and
                        lFile[15] == "2 2 0 exterior" and
                        lFile[16] == "3 0 3 " and
                        lFile[17] == "4 4 2 exterior" ,
                        'Ascii file is wrong, segment')
        
       # self.assertTrue(lFile[18] == '4 0 # <vertex #> <x> <y> [attributes] ...Mesh Vertices...',
        #                'Ascii file is wrong, Mesh Vertices Title')
        
        self.assertTrue(lFile[19] == '0 0.0 0.0 ' and #1.1 ' and
                        lFile[20] == '1 0.0 4.0 ' and #1.2 ' and
                        lFile[21] == '2 4.0 0.0 ' and #1.3 ' and
                        lFile[22] == '3 1.0 1.0 ' #1.4 '
                        ,
                        'Ascii file is wrong, Mesh Vertices II')
        
        self.assertTrue(lFile[24] == '0 0 1 ' and
                        lFile[25] == '1 1 2 ' and
                        lFile[26] == '2 0 2 ' and
                        lFile[27] == '3 0 3 '
                        ,'Ascii file is wrong, Mesh Segments')       

  
    def test_ascii_file(self):
    
        a = Vertex (0.0, 0.0) #, attributes = [1.1])
        d = Vertex (0.0, 4.0) #, attributes = [1.2])
        f = Vertex (4.0,0.0) #, attributes = [1.3])
        e = Vertex (1.0,1.0) #, attributes = [1.4])
    
        s1 = Segment(a,d)
        s2 = Segment(d,f)
        s3 = Segment(a,f)
        s4 = Segment(a,e)
     
        m = Mesh(userVertices=[a,d,f,e], userSegments=[s1,s2,s3,s4])

        m.generateMesh("Q", maxArea = 2.1 )

        seg = m.getMeshSegments()
        
        fileName = tempfile.mktemp(".tsh")
        m.export_mesh_file(fileName)
        file = open(fileName)
        lFile = file.read().split('\n')
        file.close()
        os.remove(fileName)
        
        #print "@^@^"
        #for l in lFile:
        #    print l,"<"
        #print "@^@^"
        self.assertTrue(lFile[0] == "5 0 # <# of verts> <# of vert attributes>, next lines <vertex #> <x> <y> [attributes] ...Triangulation Vertices..."
                        ,
                        'Ascii file is wrong, vertex title')
        self.assertTrue(lFile[1] == "0 0.0 0.0 " and #1.1 " and
                        lFile[2] == "1 0.0 4.0 " and #1.2 " and
                        lFile[3] == "2 4.0 0.0 " and #1.3 " and
                        lFile[4] == "3 1.0 1.0 " and #1.4 " and
                        lFile[5] == "4 2.0 2.0 "  #1.25 " 
                        ,
                        'Ascii file is wrong, vertex')
        
        self.assertTrue(lFile[6] == "# attribute column titles ...Triangulation Vertex Titles..."
                        ,
                        'Ascii file is wrong, attribute column title')
        self.assertTrue(lFile[7] == "4 # <# of triangles>, next lines <triangle #> [<vertex #>] [<neigbouring triangle #>] [attribute of region] ...Triangulation Triangles..." and
                        lFile[8] == "0 3 2 4 -1 2 3  " and
                        lFile[9] == "1 1 0 3 3 2 -1  " and
                        lFile[10] == "2 3 4 1 -1 1 0  " and
                        lFile[11] == "3 2 3 0 1 -1 0  "
                        ,
                        'Ascii file is wrong, triangle') 

        self.assertTrue(lFile[12] == "5 # <# of segments>, next lines <segment #> <vertex #>  <vertex #> [boundary tag] ...Triangulation Segments..." and
                        lFile[13] == "0 0 1 exterior" and
                        lFile[14] == "1 1 4 exterior" and
                        lFile[15] == "2 2 0 exterior" and
                        lFile[16] == "3 0 3 " and
                        lFile[17] == "4 4 2 exterior" ,
                        'Ascii file is wrong, segment')
        
        self.assertTrue(lFile[18] == '4 0 # <# of verts> <# of vert attributes>, next lines <vertex #> <x> <y> [attributes] ...Mesh Vertices...',
                        'Ascii file is wrong, Mesh Vertices Title')
        
        self.assertTrue(lFile[19] == '0 0.0 0.0 ' and #1.1 ' and
                        lFile[20] == '1 0.0 4.0 ' and #1.2 ' and
                        lFile[21] == '2 4.0 0.0 ' and #1.3 ' and
                        lFile[22] == '3 1.0 1.0 ' #1.4 '
                        ,
                        'Ascii file is wrong, Mesh Vertices II')
        
        self.assertTrue(lFile[23] == '4 # <# of segments>, next lines <segment #> <vertex #>  <vertex #> [boundary tag] ...Mesh Segments...' and
                        lFile[24] == '0 0 1 ' and
                        lFile[25] == '1 1 2 ' and
                        lFile[26] == '2 0 2 ' and
                        lFile[27] == '3 0 3 ' and
                        lFile[28] == '0 # <# of holes>, next lines <Hole #> <x> <y> ...Mesh Holes...' and
                        lFile[29] == '0 # <# of regions>, next lines <Region #> <x> <y> <tag>...Mesh Regions...'
                        ,
                        'Ascii file is wrong, Mesh Segments')       
  

    def test_thinoutVertices(self):

        v1 = Vertex(-20,-20)
        v2 = Vertex(-11,-11)
        v3 = Vertex(-10,-10)
        v4 = Vertex(-9,-1)
        v5 = Vertex(-8,2)
        v6 = Vertex(6,3)
        v7 = Vertex(12,9)
        v8 = Vertex(15,3)
        v9 = Vertex(24,3)
        m = Mesh(userVertices = [v1,v2,v3,v4,v5,v6,v7,v8,v9])
        m.thinoutVertices(10)
         
        self.assertTrue(v1 in m.userVertices,
                        'test_thinoutVertices, test 1 failed')
        self.assertTrue(v3 in m.userVertices,
                        'test_thinoutVertices, test 2 failed')
        self.assertTrue(v4 in m.userVertices,
                        'test_thinoutVertices, test 3 failed')
        self.assertTrue(v6 in m.userVertices,
                        'test_thinoutVertices, test 4 failed')
        self.assertTrue(v7 in m.userVertices,
                        'test_thinoutVertices, test 5 failed')
        self.assertTrue(v9 in m.userVertices,
                        'test_thinoutVertices, test 6 failed')
        self.assertTrue(v5 not in m.userVertices,
                        'test_thinoutVertices, test 7 failed')
        self.assertTrue(v2 not in m.userVertices,
                        'test_thinoutVertices, test 8 failed')
        self.assertTrue(v8 not in m.userVertices,
                        'test_thinoutVertices, test 9 failed')

    def test_same_x_y(self):
        v = Point(7,8)
        f = Point(7,8)
        f.same_x_y(v)

        self.assertTrue(f.same_x_y(v),
                        'same_x_y True failed')
        e = Point(7,9)
        self.assertTrue(not f.same_x_y(e),
                        'same_x_y False failed')

    def test_import_tsh(self):
        
        a_att = [5.0,2.0]
        d_att =[4.0,2.0]
        f_att = [3.0,2.0]
        e_att = [2.0,2.0]
        a_xy = [0.0, 0.0]
        a = Vertex ( a_xy[0],a_xy[1]) #, attributes =a_att)
        d = Vertex (0.0, 4.0) #, attributes =d_att)
        f = Vertex (4.0,0.0) #, attributes =f_att)
        e = Vertex (1.0,1.0) #, attributes =e_att)
    
        s1 = Segment(a,d, tag = "50")
        s2 = Segment(d,f, tag = "40")
        s3 = Segment(a,f, tag = "30")
        s4 = Segment(a,e, tag = "20")
     
        r1 = Region(0.3, 0.3,tag = "1.3")
        geo = Geo_reference(55, 8.9,8.9)
        m = Mesh(userVertices=[a,d,f,e],
                 userSegments=[s1,s2,s3,s4],
                 regions=[r1],
                 geo_reference=geo)

        m.generateMesh("Q", maxArea = 2.1)
        fileName = tempfile.mktemp(".tsh")
        #print "dgs!!!"
        #print "****************** fileName", fileName
        m.export_mesh_file(fileName)
        #print "******************"
        #print "m", m
        #print "******************"
        m_returned = importMeshFromFile(fileName)
        #print "m_returned",m_returned
        #print "******************"
        #print "****************** fileName", fileName
        os.remove(fileName)
        self.assertTrue(0 == m.__cmp__(m_returned),
                        'loading and saving of a mesh failed')
        # do this when .msh supports geo_refs
        #self.assertTrue(m.geo_reference == m_returned.geo_reference,
        #                'loading and saving of a mesh geo refs failed')

    def test_import_mesh(self):
        
        a_att = [5.0,2.0]
        d_att =[4.0,2.0]
        f_att = [3.0,2.0]
        e_att = [2.0,2.0]
        a_xy = [0.0, 0.0]
        a = Vertex ( a_xy[0],a_xy[1]) #, attributes =a_att)
        d = Vertex (0.0, 4.0) #, attributes =d_att)
        f = Vertex (4.0,0.0) #, attributes =f_att)
        e = Vertex (1.0,1.0) #, attributes =e_att)
    
        s1 = Segment(a,d, tag = "50")
        s2 = Segment(d,f, tag = "40")
        s3 = Segment(a,f, tag = "30")
        s4 = Segment(a,e, tag = "20")
     
        r1 = Region(0.3, 0.3,tag = "1.3")
        geo = Geo_reference(55,8.9,8.9)
        m = Mesh(userVertices=[a,d,f,e],
                 userSegments=[s1,s2,s3,s4],
                 regions=[r1],
                 geo_reference=geo)

        m.generateMesh("Q", maxArea = 2.1)
        fileName = tempfile.mktemp(".msh")
        #print "dgs!!!"
        #print "****************** fileName", fileName
        m.export_mesh_file(fileName)
        #print "******************"
        #print "m", m
        #print "******************"
        m_returned = importMeshFromFile(fileName)
        #print "m_returned",m_returned
        #print "******************"
        #print "****************** fileName", fileName
        os.remove(fileName)
        #print "m.geo_reference",m.geo_reference
        #print "m_returned.geo_reference,",m_returned.geo_reference
        self.assertTrue(0 == m.__cmp__(m_returned),
                        'loading and saving of a mesh failed')

        self.assertTrue(m.geo_reference == m_returned.geo_reference,
                        'loading and saving of a mesh geo refs failed')

    def DONTtest_normaliseMesh(self):
        
        a_att = [5.0,2.0]
        d_att =[4.0,2.0]
        f_att = [3.0,2.0]
        e_att = [2.0,2.0]
        a_xy = [10.0, 10.0]
        a = Vertex ( a_xy[0],a_xy[1], attributes =a_att)
        d = Vertex (15.0, 10.0, attributes =d_att)
        f = Vertex (10.0,20.0, attributes =f_att)
        e = Vertex (15.0,20.0, attributes =e_att)
    
        s1 = Segment(a,d, tag = 50)
        s2 = Segment(d,e, tag = 40)
        s3 = Segment(e,f, tag = 30)
        s4 = Segment(f,a, tag = 20)
     
        r1 = Region(0.3, 0.3,tag = 1.3)
        m = Mesh(userVertices=[a,d,f,e],
                 userSegments=[s1,s2,s3,s4],
                 regions=[r1])
        m.normaliseMesh(1,0,1)
        [xmin, ymin, xmax, ymax] = m.boxsize()
        [attmin, attmax] = m.maxMinVertAtt(0)
        self.assertTrue(attmin == 0.0 and attmax == 1.0,
                        'normalise failed')
        self.assertTrue(xmin == 0.0 and ymin == 0.0 and xmax == 0.5 and ymax == 1.0,
                        'normalise failed')
        m.normaliseMesh(200,-100,5)
        [xmin, ymin, xmax, ymax] = m.boxsize()
        [attmin, attmax] = m.maxMinVertAtt(0)
        self.assertTrue(attmin == 0.0 and attmax == 5.0,
                        'normalise failed')
        self.assertTrue(xmin == -100.0 and ymin == -100.0 and xmax == 0.0 and ymax == 100.0,
                        'normalise failed')
        
    def test_exportASCIIsegmentoutlinefile1(self):
        a = Vertex (0,0)
        b = Vertex (0,3)
        c = Vertex (3,3)
        d = Vertex (1,2)
        e = Vertex (3,1)
      
        s1 = Segment(a,b, tag = "50")
        s2 = Segment(b,c, tag = "40")
        s3 = Segment(c,a, tag = "30")
     
        r1 = Region(2, 1,tag = "1.3")
        h1 = Hole(1,4)
        m = Mesh(userVertices=[a,b,c,d,e],
                 userSegments=[s1,s2,s3],
                 regions=[r1],
                 holes = [h1])      

        # vertex e is outside of the outline, so
        # it is a loner and it is removed.
        m.generateMesh("Q", maxArea = 2.1)
        #print "mesh ***************dsg*", m

        fileName = tempfile.mktemp(".tsh")
        m.exportASCIIsegmentoutlinefile(fileName)
        
        m_returned = importMeshFromFile(fileName)

        #print "m_returned ****",m_returned
        #print "****************** fileName", fileName
        os.remove(fileName)

        #Trim mesh, so it should like like m_returned
        m.tri_mesh = None
        m.userVertices=[a,b,c]
        #print("mesh ***************dsg*", m)
        #print "(m.__cmp__(m_returned)", m.__cmp__(m_returned)
        #print('result', m.__cmp__(m))
        self.assertTrue(0 == m.__cmp__(m),
                        'test_exportASCIIsegmentoutlinefile:loading and saving of a mesh failed')
        # Having problems with this on linux.
        #The ordering of the dictionary values wasn't the same as the windows
        #returned value (verts.values())
        #self.assertTrue(0 == m.__cmp__(m_returned),
        #                'test_exportASCIIsegmentoutlinefile:loading and saving of a mesh failed')
        
        self.assertTrue(3 == len(m_returned.userVertices),
                        'segmentoutlinefile:IO of a mesh failed')
        self.assertTrue(len(m.userSegments) == len(m_returned.userSegments),
                        'segmentoutlinefile:IO of a mesh failed')
        for i in range(len(m.userSegments)):
            self.assertTrue(m.userSegments[i].vertices[0].x ==
                            m_returned.userSegments[i].vertices[0].x,
                        'loading and saving of a mesh outline fialed')
            self.assertTrue(m.userSegments[i].vertices[0].y ==
                            m_returned.userSegments[i].vertices[0].y,
                        'loading and saving of a mesh outline fialed')
            self.assertTrue(m.userSegments[i].vertices[1].x ==
                            m_returned.userSegments[i].vertices[1].x,
                        'loading and saving of a mesh outline fialed')
            self.assertTrue(m.userSegments[i].vertices[1].y ==
                            m_returned.userSegments[i].vertices[1].y,
                        'loading and saving of a mesh outline fialed')

 
    def test_exportASCIIsegmentoutlinefile2(self):
        a = Vertex (0,0)
        b = Vertex (0,1)
        c = Vertex (1,0)
        d = Vertex (1,1)
        e = Vertex (0.5,0.5)
        f  = Vertex (0.6,0.6)
      
        s1 = Segment(a,e, tag = "50")
        s2 = Segment(b,e, tag = "40")
        s3 = Segment(c,e, tag = "30")
        s4 = Segment(d,e, tag = "30")
     
        r1 = Region(2, 1,tag = "1.3")
        h1 = Hole(1,4)
        m = Mesh(userVertices=[a,b,c,d,e],
                 userSegments=[s1,s2,s3,s4],
                 regions=[r1],
                 holes = [h1])      
        
        fileName = tempfile.mktemp(".tsh")
        m.exportASCIIsegmentoutlinefile(fileName)
        
        m_returned = importMeshFromFile(fileName)
        #print "****************** fileName", fileName
        os.remove(fileName)

        #Trim mesh, so it should look like m_returned
        m.meshVertices = []
        m.meshTriangles = []
        m.meshSegments = []
        m.userVertices=[a,e,d,b,c]
        #print "mesh ***************dsg*", m
        #print "(m.__cmp__(m_returned)", m.__cmp__(m_returned) 
        self.assertTrue(0 == m.__cmp__(m),
                        'loading and saving of a mesh failed')

        self.assertTrue(5 == len(m_returned.userVertices),
                        'segmentoutlinefile:IO of a mesh failed')
        self.assertTrue(len(m.userSegments) == len(m_returned.userSegments),
                        'segmentoutlinefile:IO of a mesh failed')
        for i in range(len(m.userSegments)):
            self.assertTrue(m.userSegments[i].vertices[0].x ==
                            m_returned.userSegments[i].vertices[0].x,
                        'loading and saving of a mesh outline fialed')
            self.assertTrue(m.userSegments[i].vertices[0].y ==
                            m_returned.userSegments[i].vertices[0].y,
                        'loading and saving of a mesh outline fialed')
            self.assertTrue(m.userSegments[i].vertices[1].x ==
                            m_returned.userSegments[i].vertices[1].x,
                        'loading and saving of a mesh outline fialed')
            self.assertTrue(m.userSegments[i].vertices[1].y ==
                            m_returned.userSegments[i].vertices[1].y,
                        'loading and saving of a mesh outline fialed')


    def test_load_csv(self):
        """
        To test the mesh side of loading csv files.
        Not the loading of csv files
        """
        import os
        import tempfile
       
        fileName = tempfile.mktemp(".csv")
        file = open(fileName,"w")
        file.write("x,y,elevation, speed \n\
1.0, 0.0, 10.0, 0.0\n\
0.0, 1.0, 0.0, 10.0\n\
1.0, 0.0, 10.4, 40.0\n")
        file.close()
        #print fileName
        m = importMeshFromFile(fileName)
        os.remove(fileName)
        self.assertTrue(m.userVertices[0].x == 1.0,
                        'loadxy, test 1 failed')
        self.assertTrue(m.userVertices[0].y == 0.0,
                        'loadxy, test 2 failed')
        #self.assertTrue(m.userVertices[0].attributes == [10.0,0.0],
        #                'loadxy, test 2.2 failed')
        self.assertTrue(m.userVertices[1].x == 0.0,
                        'loadxy, test 3 failed')
        self.assertTrue(m.userVertices[1].y == 1.0,
                        'loadxy, test 4 failed')
        #self.assertTrue(m.userVertices[1].attributes == [0.0,10.0],
        #                'loadxy, test 5 failed')
        
    def test_exportPointsFile(self):
        a = Vertex (0,0)
        b = Vertex (0,3)
        c = Vertex (3,3)
        d = Vertex (1,2)
        e = Vertex (3,1)
        #f = Vertex (3,1)
      
        s1 = Segment(a,b, tag = 50)
        s2 = Segment(b,c, tag = 40)
        s3 = Segment(c,a, tag = 30)
     
        r1 = Region(2, 1,tag = 1.3)
        h1 = Hole(1,4)
        # Warning mesh can't produce this type of data structure its self
        m = Mesh(userVertices=[a,b,c,d,e],
                 userSegments=[s1,s2,s3],
                 regions=[r1],
                 holes = [h1])
        
        fileName = tempfile.mktemp(".txt")
        #fileName = 't.csv'
        #os.remove(fileName)
        m.exportPointsFile(fileName)
        file = open(fileName)
        lFile = file.read().split('\n')
        file.close()
        os.remove(fileName)
        self.assertTrue(lFile[0] == "x,y," and
                        lFile[1] == "0.0,0.0" and
                        lFile[2] == "0.0,3.0" and
                        lFile[3] == "3.0,3.0" 
                        ,
                        'exported Ascii csv file is wrong')
        self.assertTrue(lFile[4] == "1.0,2.0" and
                        lFile[5] == "3.0,1.0" 
                        ,
                        'exported Ascii csv file is wrong')
        
        # vertex e is outside of the outline, so
        # it is a loner and it is removed.
        m.generateMesh("Q", maxArea = 2.1)
        fileName = tempfile.mktemp(".txt")
        #fileName = 't.csv'
        #m.export_mesh_file('m.tsh')
        m.exportPointsFile(fileName)
        file = open(fileName)
        lFile = file.read().split('\n')
        file.close()
        os.remove(fileName)
        
        self.assertTrue(lFile[0] == "x,y," and
                        lFile[1] == "0.0,0.0" and
                        lFile[2] == "0.0,3.0" and
                        lFile[3] == "3.0,3.0" and
                        lFile[4] == "1.0,2.0"
                        ,
                        'exported Ascii csv file is wrong')
     
    def to_be_lone_vert_in_mesh_gen_c_layer(self):
        # currently just a copy of the above test
        a = Vertex (0,0)
        b = Vertex (0,3)
        c = Vertex (3,3)
        d = Vertex (1,2)
        e = Vertex (3,1)
        #f = Vertex (3,1)
      
        s1 = Segment(a,b, tag = 50)
        s2 = Segment(b,c, tag = 40)
        s3 = Segment(c,a, tag = 30)
     
        r1 = Region(2, 1,tag = 1.3)
        h1 = Hole(1,4)
        # Warning mesh can't produce this type of data structure its self
        m = Mesh(userVertices=[a,b,c,d,e],
                 userSegments=[s1,s2,s3],
                 regions=[r1],
                 holes = [h1])
        
        fileName = tempfile.mktemp(".csv")
        #fileName = 't.csv'
        #os.remove(fileName)
        m.exportPointsFile(fileName)
        file = open(fileName)
        lFile = file.read().split('\n')
        file.close()

        os.remove(fileName)
        self.assertTrue(lFile[0] == "x,y" and
                        lFile[1] == "0,0" and
                        lFile[2] == "0,3" and
                        lFile[3] == "3,3" 
                        ,
                        'exported Ascii csv file is wrong')
        self.assertTrue(lFile[4] == "1,2" and
                        lFile[5] == "3,1" 
                        ,
                        'exported Ascii csv file is wrong')
        
        # vertex e is outside of the outline, so
        # it is a loner and it is removed.
        m.generateMesh("Q", maxArea = 2.1)
        fileName = tempfile.mktemp(".csv")
        #fileName = 't.csv'
        #m.export_mesh_file('m.tsh')
        m.exportPointsFile(fileName)
        file = open(fileName)
        lFile = file.read().split('\n')
        file.close()
        os.remove(fileName)
        
        self.assertTrue(lFile[0] == "x,y" and
                        lFile[1] == "0.0,0.0" and
                        lFile[2] == "0.0,3.0" and
                        lFile[3] == "3.0,3.0" and
                        lFile[4] == "1.0,2.0"
                        ,
                        'exported Ascii csv file is wrong')
        
    def NOTtest_exportPointsFilefile2(self):
        #geospatial needs at least one point
        m = Mesh()
        
        fileName = tempfile.mktemp(".csv")
        m.exportPointsFile(fileName)
        file = open(fileName)
        lFile = file.read().split('\n')
        file.close()

        os.remove(fileName)
        #print "************* test_mesh exportPointsFilefile"
        #print "lFile",lFile 
        #print "************* test_mesh exportPointsFilefile"
        self.assertTrue(lFile[0] == "" 
                        ,
                        'exported Ascii csv file is wrong')

    def test_segment_strings2ints(self):
        """ Test directly using its own example from its docstring
            input ['a','b','a','c'], ['c']
            output [[2, 1, 2, 0], ['c', 'b', 'a']]
        """

        stringlist = ['a','b','a','c']
        preset = ['c']
        [intlist, converter] = segment_strings2ints(stringlist, preset)

        assert intlist == [2, 1, 2, 0]
        converter == ['c', 'b', 'a']

        
    def test_strings2ints(self):
        # Expected result 
        outlist = ['sea', 'river inlet', 'moat',
                   'sea', 'moat', 'moat']
        
        # Input
        list = ["sea", "river inlet", "", "sea", "", "moat"]
        preset = ["moat", "internal boundary"]
        [intlist, converter] = segment_strings2ints(list, preset)

        # Instead of testing the converter and the intlist separately, test that they
        # evaluate to the correct result
        for i, k in enumerate(intlist):
            assert converter[k] == outlist[i]

        # And double check the inverse funciont
        newlist = segment_ints2strings(intlist, converter)
        for i, name in enumerate(newlist):
            assert name == outlist[i]            
            
        
    def test_ints2strings1(self):
        list = ["internal boundary","sea","river inlet",
            "","sea","","moat","internal boundary"]
        outlist = ['internal boundary', 'sea', 'river inlet', 'moat',
                   'sea', 'moat', 'moat', 'internal boundary']
        preset = ["moat", "internal boundary"]
        [intlist, converter] = segment_strings2ints(list, preset)
        
        newlist = segment_ints2strings(intlist, converter)
        
        self.assertTrue(outlist == newlist,
                        'test_strings2ints produces wrong result')
        #self.assertTrue(converter == ['moat', 'internal boundary',
        #                              'sea', 'river inlet'],
        #                'test_strings2ints produces bad converter')

        # Instead of testing the converter and the intlist separately, test that they
        # evaluate to the correct result
        for i, k in enumerate(intlist):
            #print(i, k, converter[k], outlist[i])
            assert converter[k] == outlist[i]
        
    def test_ints2strings2(self):
        list = ["","",""]
        preset = ["moat", "internal boundary"]
        [intlist, converter] = segment_strings2ints(list,preset )
        newlist = segment_ints2strings(intlist, converter)
        outlist = ['moat', 'moat', 'moat']
        self.assertTrue(outlist == newlist,
                        'test_strings2ints produces bad intlist')
        self.assertTrue(converter == ['moat', 'internal boundary'],
                        'test_strings2ints produces bad converter')

        
    def test_removeDuplicatedVertices(self):
        a = Vertex (0,0)
        a.index = 0
        b = Vertex (0,3)
        b.index = 1
        c = Vertex (3,3)
        c.index = 2
        d = Vertex (1,1)
        d.index = 3
        e = Vertex (3,1)
        e.index = 4
        f = Vertex (1,1)
        f.index = 5
        g = Vertex (1,1)
        g.index = 6
        inputVerts_noDups = [a,b,c,d,e]
        
        m = Mesh(userVertices=[a,b,c,d,e,f,g])
        counter = m.removeDuplicatedUserVertices()
        UserVerts = m.getUserVertices()
        
        self.assertTrue(UserVerts == inputVerts_noDups,
                            'duplicate verts not removed')
        #for userVert, inputVert in map(None, UserVerts, inputVerts_noDups): 
        #    self.assertTrue(userVert.x == inputVert.x,
        #                    'x duplicate verts not removed')
        #    self.assertTrue(userVert.y == inputVert.y,
        #                    'y duplicate verts not removed')

        
    def test_addVertsSegs(self):
        m = Mesh()
        Segment.set_default_tag("food")
        dict = {}
        dict['points'] = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]
        dict['segments'] = [[0, 1], [1, 2]]
        dict['segment_tags'] = ['','do-op']
        m.addVertsSegs(dict)
        # have to reset this , since it's a class attribute
        Segment.set_default_tag("")

        
        self.assertTrue(len(m.userSegments) ==2,
                        'Wrong segment list length.')
        self.assertTrue(len(m.userVertices) == 3,
                        'Wrong vertex list length.')
        self.assertTrue(m.userSegments[0].tag =='food',
                        'Wrong segment tag length.')
        self.assertTrue(m.userSegments[1].tag =='do-op',
                        'Wrong segment tag.')
        
    def test_addVertsSegs2(self):
        geo = Geo_reference(56,5,10)
        m = Mesh(geo_reference=geo)
        dict = {}
        dict['points'] = [[2.0, 1.0], [3.0, 1.0], [2.0, 2.0]]
        dict['segments'] = [[0, 1], [1, 2], [2,0]]
        dict['segment_tags'] = ['','do-op','']
        m.addVertsSegs(dict)

    def test_addVertsSegs_done_twice(self):
        m = Mesh()
        dict = {}
        dict['points'] = [[0.0, 0.0], [5.0, 0.0], [5.0, 5.0]]
        dict['segments'] = [[0, 1], [1, 2], [2,0]]
        dict['segment_tags'] = ['0','1','2']
        m.addVertsSegs(dict)
        
        dict = {}
        dict['points'] = [[2.0, 1.0], [4.0, 1.0], [4.0, 3.0]]
        dict['segments'] = [[0, 1], [1, 2], [2,0]]
        dict['segment_tags'] = ['3','4','5']
        m.addVertsSegs(dict)

        
        self.assertTrue(m.userSegments[5].vertices[0].y == 3,
                        'Wrong vertex connected.')
        self.assertTrue(m.userSegments[5].vertices[1].y == 1,
                        'Wrong vertex connected.')
            
    def test_add_points_and_segments(self):
        m = Mesh()
        Segment.set_default_tag("food")
        dict = {}
        points =  [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]
        segments = [[0, 1], [1, 2]]
        segment_tags = {'hair':[1]}
        m.add_points_and_segments(points,
                                    segments, segment_tags)
        # have to reset this , since it's a class attribute
        Segment.set_default_tag("")

        
        self.assertTrue(len(m.userSegments) ==2,
                        'Wrong segment list length.')
        self.assertTrue(len(m.userVertices) == 3,
                        'Wrong vertex list length.')
        self.assertTrue(m.userSegments[0].tag =='food',
                        'Wrong segment tag length.')
        self.assertTrue(m.userSegments[1].tag =='hair',
                        'Wrong segment tag.')
        
    def test_add_points_and_segmentsII(self):
        m = Mesh()
        Segment.set_default_tag("food")
        dict = {}
        points =  [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]
        segments = None #[[0, 1], [1, 2]]
        segment_tags = {'hair':[1]}
        m.add_points_and_segments(points,
                                    segments, segment_tags)
        # have to reset this , since it's a class attribute
        Segment.set_default_tag("")

        
        self.assertTrue(len(m.userSegments) ==2,
                        'Wrong segment list length.')
        self.assertTrue(len(m.userVertices) == 3,
                        'Wrong vertex list length.')
        self.assertTrue(m.userSegments[0].tag =='food',
                        'Wrong segment tag length.')
        self.assertTrue(m.userSegments[1].tag =='hair',
                        'Wrong segment tag.')
        
    def test_exportASCIImeshfile(self):
    
        #a_att = [5,2]
        #d_att =[4,2]
        #f_att = [3,2]
        #e_att = [2,2]
        a_xy = [0.0, 0.0]
        a = Vertex ( a_xy[0],a_xy[1]) #, attributes =a_att)
        d = Vertex (0.0, 4.0) #, attributes =d_att)
        f = Vertex (4.0,0.0) #, attributes =f_att)
        e = Vertex (1.0,1.0) #, attributes =e_att)
    
        s1 = Segment(a,d, tag = "50")
        s2 = Segment(d,f, tag = "40")
        s3 = Segment(a,f, tag = "30")
        s4 = Segment(a,e, tag = "20")
     
        r1 = Region(0.3, 0.3,tag = "1.3", maxArea = 36)


        h1 = Hole(0.2,0.6)
        
        m = Mesh(userVertices=[a,d,f,e],
                 userSegments=[s1,s2,s3,s4],
                 regions=[r1],
                 holes=[h1])

        seg = m.getUserSegments()
        points = m.getUserVertices()
        holes = m.getHoles()
        regions = m.getRegions()
        fileName = tempfile.mktemp(".tsh")
        m.export_mesh_file(fileName)
        #print "***************************fileName", fileName
        new_m = importMeshFromFile(fileName)
        os.remove(fileName)
        

        #print '**@@@@@******'
        #print "new_m",new_m 
        #print '**@@@@@******'
        #print "m",m 
        #print '**@@@@@******'

        self.assertTrue(new_m == m,
                        'loadASCIITestCase failed. test new 1')
            
    def test_Mesh2MeshList(self):

        a_att = [5,2]
        d_att =[4,2]
        f_att = [3,2]
        e_att = [2,2]
        a_xy = [0.0, 0.0]
        a = Vertex ( a_xy[0],a_xy[1]) #, attributes =a_att)
        d = Vertex (0.0, 4.0) #, attributes =d_att)
        f = Vertex (4.0,0.0) #, attributes =f_att)
        e = Vertex (1.0,1.0) #, attributes =e_att)
    
        s1 = Segment(a,d, tag = "50")
        s2 = Segment(d,f, tag = "40")
        s3 = Segment(a,f, tag = "30")
        s4 = Segment(a,e, tag = "20")
     
        r1 = Region(0.3, 0.3,tag = "1.3", maxArea = 45)
        m = Mesh(userVertices=[a,d,f,e],
                 userSegments=[s1,s2,s3,s4],
                 regions=[r1])

        m.generateMesh("Qa2.1")

        seg = m.getMeshSegments()
        points = m.getMeshVertices()
        dict = m.Mesh2MeshList()
        #print "dict",dict 
        # test not finished...
  
    def test_Mesh2IOTriangulationDict(self):

        a_att = [5,2]
        d_att =[4,2]
        f_att = [3,2]
        e_att = [2,2]
        a_xy = [0.0, 0.0]
        a = Vertex ( a_xy[0],a_xy[1] , attributes =a_att)
        d = Vertex (0.0, 4.0 , attributes =d_att)
        f = Vertex (4.0,0.0 , attributes =f_att)
        e = Vertex (1.0,1.0 , attributes =e_att)
    
        s1 = Segment(a,d, tag = "50")
        s2 = Segment(d,f, tag = "40")
        s3 = Segment(a,f, tag = "30")
        s4 = Segment(a,e, tag = "20")
     
        r1 = Region(0.3, 0.3,tag = "1.3", maxArea = 45)
        m = Mesh(userVertices=[a,d,f,e],
                 userSegments=[s1,s2,s3,s4],
                 regions=[r1])
        titles = ['ele','friction'] #Feed in directly!
        m.attributeTitles = titles
        m.generateMesh("Qa2.1")

        seg = m.getMeshSegments()
        verts = m.getMeshVertices()
        vert_as = m.getMeshVerticeAttributes()
        seg_tags = m.getMeshSegmentTags()
        dict = m.Mesh2IOTriangulationDict()
        #print "dict",dict 
        
        self.assertTrue( dict['vertex_attribute_titles'] == titles,
                         'test_Mesh2IOTriangulationDict failed. test 1')
        answer = [a_xy,[0.0, 4.0],[4.0,0.0], [1.0,1.0], [2.0,2.0]]
        #print "answer",answer
        #print "dict['vertices']",dict['vertices']
        
        self.assertTrue(num.alltrue(dict['vertices'] == answer),
                        'test_Mesh2IOTriangulationDict failed. test 2')

        self.assertTrue(num.alltrue(dict['vertices'].flatten() ==
                                    verts.flatten()),
                         'test_Mesh2IOTriangulationDict failed. test vert')
        self.assertTrue(num.alltrue(dict['vertex_attributes'].flatten() ==
                                    vert_as.flatten()),
                         'test_Mesh2IOTriangulationDict failed. test vert ats')

        self.assertTrue(num.alltrue(dict['segments'][0] == [0,1]),
                        'test_Mesh2IODict failed. test 3')
        
        self.assertTrue( dict['segment_tags'] == seg_tags,
                        'test_Mesh2IODict failed. test 3')
        #print " dict['triangles'][0]", dict['triangles'][0] 
        self.assertTrue(num.alltrue(dict['triangles'][0] == [3,2,4]),
                        'test_Mesh2IODict failed. test 5')
        self.assertTrue(num.alltrue(dict['triangle_neighbors'][0] == [-1,2,3]),
                        'test_Mesh2IODict failed. test 6')
        #print "dict['triangle_tags'][0]", dict['triangle_tags'][0]
        self.assertTrue( dict['triangle_tags'][0] == "1.3",
                         'test_Mesh2IODict failed. test 7')

  
    def test_Mesh2IODict(self):

        a_att = [5,2]
        d_att =[4,2]
        f_att = [3,2]
        e_att = [2,2]
        a_xy = [0.0, 0.0]
        a = Vertex ( a_xy[0],a_xy[1] , attributes =a_att)
        d = Vertex (0.0, 4.0 , attributes =d_att)
        f = Vertex (4.0,0.0 , attributes =f_att)
        e = Vertex (1.0,1.0 , attributes =e_att)
    
        s1 = Segment(a,d, tag = "50")
        s2 = Segment(d,f, tag = "40")
        s3 = Segment(a,f, tag = "30")
        s4 = Segment(a,e, tag = "20")
     
        r1 = Region(0.3, 0.3,tag = "1.3", maxArea = 45)
        m = Mesh(userVertices=[a,d,f,e],
                 userSegments=[s1,s2,s3,s4],
                 regions=[r1])
        titles = ['ele','friction']
        m.attributeTitles = titles
        m.generateMesh("Qa2.1")

        seg = m.getMeshSegments()
        verts = m.getMeshVertices()
        vert_as = m.getMeshVerticeAttributes()
        dict = m.Mesh2IODict()
        seg_tags = m.getMeshSegmentTags()
        #print "dict",dict 
        
        self.assertTrue( dict['vertex_attribute_titles'] == titles,
                         'test_Mesh2IOTriangulationDict failed. test 1')
        answer = [a_xy,[0.0, 4.0],[4.0,0.0], [1.0,1.0], [2.0,2.0]]
        #print "answer",answer
        #print "dict['vertices']",dict['vertices']
        
        self.assertTrue(num.alltrue(dict['vertices'] == answer),
                        'test_Mesh2IOTriangulationDict failed. test 2')

        self.assertTrue(num.alltrue(dict['vertices'] == verts),
                        'test_Mesh2IOTriangulationDict failed. test vert')
        self.assertTrue(num.alltrue(dict['vertex_attributes'] == vert_as),
                        'test_Mesh2IOTriangulationDict failed. test vert ats')

        self.assertTrue(num.alltrue(dict['segments'][0] == [0,1]),
                        'test_Mesh2IODict failed. test 3')
        
        self.assertTrue(dict['segment_tags'] == seg_tags,
                        'test_Mesh2IODict failed. test 3')
        #print " dict['triangles'][0]", dict['triangles'][0] 
        self.assertTrue(num.alltrue(dict['triangles'][0] == [3,2,4]),
                        'test_Mesh2IODict failed. test 5')
        self.assertTrue(num.alltrue(dict['triangle_neighbors'][0] == [-1,2,3]),
                        'test_Mesh2IODict failed. test 6')
        #print "dict['triangle_tags'][0]", dict['triangle_tags'][0]
        self.assertTrue(dict['triangle_tags'][0] == "1.3",
                        'test_Mesh2IODict failed. test 7')

        seg = m.getUserSegments()
        points = m.getUserVertices()
        holes = m.getHoles()
        regions = m.getRegions()
        
        for pimport,pactual,pimpatt in zip(dict['points'],points,dict['point_attributes']):
            self.assertTrue( pimport == [pactual.x,pactual.y],
                        'test_Mesh2IODict failed. test 1')
            self.assertTrue( pimpatt == pactual.attributes,
                        'test_Mesh2IODict failed. test 1.1')
        self.assertTrue( dict['outline_segments'][0] == [0,1],
                        'test_Mesh2IODict failed. test 3')
        for segimp,segactual in zip(dict['outline_segment_tags'],seg):
            self.assertTrue( segimp == segactual.tag,
                        'test_Mesh2IODict failed. test 4')
        for holeimp,holeactual in zip(dict['holes'],holes):
            self.assertTrue( holeimp == [holeactual.x,holeactual.y],
                        'test_Mesh2IODict failed. test 5')
        
        for regimp,regactual,regattimp, regmaxarea in zip(dict['regions'],regions, dict['region_tags'], dict['region_max_areas']):
            self.assertTrue( regimp == [regactual.x,regactual.y],
                        'loadASCIITestCase failed. test 6')
            self.assertTrue( regattimp == regactual.getTag(),
                        'loadASCIITestCase failed. test 7')
            self.assertTrue( regmaxarea == regactual.getMaxArea(),
                        'loadASCIITestCase failed. test 7')
    
            
        
    def test_Mesh2IOOutlineDict(self):

        a_att = [5,2]
        d_att =[4,2]
        f_att = [3,2]
        e_att = [2,2]
        a_xy = [0.0, 0.0]
        a = Vertex ( a_xy[0],a_xy[1] , attributes =a_att)
        d = Vertex (0.0, 4.0 , attributes =d_att)
        f = Vertex (4.0,0.0 , attributes =f_att)
        e = Vertex (1.0,1.0 , attributes =e_att)
    
        s1 = Segment(a,d, tag = "50")
        s2 = Segment(d,f, tag = "40")
        s3 = Segment(a,f, tag = "30")
        s4 = Segment(a,e, tag = "20")
     
        r1 = Region(0.3, 0.3,tag = "1.3", maxArea = 45)
        m = Mesh(userVertices=[a,d,f,e],
                 userSegments=[s1,s2,s3,s4],
                 regions=[r1])
        titles = ['ele','friction']
        m.attributeTitles = titles
        m.generateMesh("Qa2.1")

        seg = m.getMeshSegments()
        verts = m.getMeshVertices()
        dict = m.Mesh2IOOutlineDict()
        
        seg = m.getUserSegments()
        points = m.getUserVertices()
        holes = m.getHoles()
        regions = m.getRegions()
        
        for pimport,pactual,pimpatt in zip(dict['points'],points,dict['point_attributes']):
            self.assertTrue( pimport == [pactual.x,pactual.y],
                        'loadASCIITestCase failed. test 1')
            self.assertTrue( pimpatt == pactual.attributes,
                        'loadASCIITestCase failed. test 1.1')
        self.assertTrue( dict['outline_segments'][0] == [0,1],
                        'loadASCIITestCase failed. test 3')
        for segimp,segactual in zip(dict['outline_segment_tags'],seg):
            self.assertTrue( segimp == segactual.tag,
                        'loadASCIITestCase failed. test 4')
        for holeimp,holeactual in zip(dict['holes'],holes):
            self.assertTrue( holeimp == [holeactual.x,holeactual.y],
                        'loadASCIITestCase failed. test 5')
        #for regimp,regactual in map(None,dict['regions'],regions):
         #   self.assertTrue( [regimp[0],regimp[1]]==[regactual.x,regactual.y],
          #              'loadASCIITestCase failed. test 6')
           # self.assertTrue( regimp[2] == regactual.getTag(),
            #            'loadASCIITestCase failed. test 7')
            #self.assertTrue( regimp[3] == regactual.getMaxArea(),
             #           'loadASCIITestCase failed. test 7')

            
        for regimp,regactual,regattimp, regmaxarea in zip(dict['regions'],regions, dict['region_tags'], dict['region_max_areas']):
            self.assertTrue( regimp == [regactual.x,regactual.y],
                        'loadASCIITestCase failed. test 6')
            self.assertTrue( regattimp == regactual.getTag(),
                        'loadASCIITestCase failed. test 7')
            self.assertTrue( regmaxarea == regactual.getMaxArea(),
                        'loadASCIITestCase failed. test 7')


    def test_add_region_from_polygon(self):
        m=Mesh()
        region = m.add_region_from_polygon([[0,0],[1,0],[0,1]],
                                  max_triangle_area = 88,
                                           region_tag='cassady')
        self.assertTrue(len(m.regions)==1,
                        'FAILED!')
        self.assertTrue(region.getMaxArea()==88,
                        'FAILED!')
        self.assertTrue(len(m.getUserSegments())==3,
                        'FAILED!')
        self.assertTrue(len(m.userVertices)==3,
                        'FAILED!')
        self.assertTrue(region.getTag()=='cassady',
                        'FAILED!')
       
    def test_add_region_from_polygon2(self):
        m=Mesh()
        m.add_region_from_polygon([[0,0],[1,0],[1,1],[0,1]],
                               {'tagin':[0,1],'bom':[2]},
                                  max_triangle_area=10)
        self.assertTrue(len(m.regions)==1,
                        'FAILED!')
        segs = m.getUserSegments()
        self.assertTrue(len(segs)==4,
                        'FAILED!')
        self.assertTrue(len(m.userVertices)==4,
                        'FAILED!') 
        self.assertTrue(segs[0].tag=='tagin',
                        'FAILED!')  
        self.assertTrue(segs[1].tag=='tagin',
                        'FAILED!') 
         
        self.assertTrue(segs[2].tag=='bom',
                        'FAILED!')
        self.assertTrue(segs[3].tag=='',
                        'FAILED!') 
       
    def test_add_region_from_polygon3(self):
        x=-500
        y=-1000
        m=Mesh(geo_reference=Geo_reference(56,x,y))

        # These are the absolute values
        polygon_absolute = [[0,0],[1,0],[1,1],[0,1]]
        
        x_p = -10
        y_p = -40
        geo_ref_poly = Geo_reference(56, x_p, y_p)
        polygon = geo_ref_poly.change_points_geo_ref(polygon_absolute)
        
        poly_point = m.add_region_from_polygon(polygon,
                                               {'tagin':[0,1],'bom':[2]},
                                               geo_reference=geo_ref_poly,
                                               max_triangle_area=10)
        # poly_point values are relative to the mesh geo-ref
        # make them absolute
        self.assertTrue(is_inside_polygon([poly_point.x+x,poly_point.y+y],
                                       polygon_absolute, closed = False),
                        'FAILED!')
               
        self.assertTrue(len(m.regions)==1,
                        'FAILED!')
        segs = m.getUserSegments()
        self.assertTrue(len(segs)==4,
                        'FAILED!')
        self.assertTrue(len(m.userVertices)==4,
                        'FAILED!') 
        self.assertTrue(segs[0].tag=='tagin',
                        'FAILED!')  
        self.assertTrue(segs[1].tag=='tagin',
                        'FAILED!') 
         
        self.assertTrue(segs[2].tag=='bom',
                        'FAILED!') 
        self.assertTrue(segs[3].tag=='',
                        'FAILED!')
        verts = m.getUserVertices()
        #print "User verts",verts
        #print 'polygon',polygon
        #vert values are relative
        for point,new_point in zip(polygon,verts):
            point_x = point[0] + geo_ref_poly.get_xllcorner()
            new_point_x = new_point.x + m.geo_reference.get_xllcorner()
            point_y = point[1] + geo_ref_poly.get_yllcorner()
            #print "new_point.y",new_point.y
            #print "m.geo_ref.get_yllcorner()",m.geo_reference.get_yllcorner() 
            new_point_y = new_point.y + m.geo_reference.get_yllcorner()
            #print "point_y",point_y
            #print "new_point_y",new_point_y 
            
            self.assertTrue(point_x == new_point_x, ' failed')
            self.assertTrue(point_y == new_point_y, ' failed')
            
         
    def test_add_region_from_polygon4(self):
        x=50000
        y=1000
        m=Mesh(geo_reference=Geo_reference(56,x,y))
        polygon = [[0,0],[1,0],[1,1],[0,1]]
        
        m.add_region_from_polygon(polygon,
                               {'tagin':[0,1],'bom':[2]},
                                  max_triangle_area=10)
        self.assertTrue(len(m.regions)==1,
                        'FAILED!')
        segs = m.getUserSegments()
        self.assertTrue(len(segs)==4,
                        'FAILED!')
        self.assertTrue(len(m.userVertices)==4,
                        'FAILED!') 
        self.assertTrue(segs[0].tag=='tagin',
                        'FAILED!')  
        self.assertTrue(segs[1].tag=='tagin',
                        'FAILED!') 
         
        self.assertTrue(segs[2].tag=='bom',
                        'FAILED!') 
        self.assertTrue(segs[3].tag=='',
                        'FAILED!')
        verts = m.getUserVertices()
        #print "User verts",verts
        #print 'polygon',polygon
        #vert values are relative
        for point,new_point in zip(polygon,verts):
            point_x = point[0] 
            new_point_x = new_point.x + m.geo_reference.get_xllcorner()
            #print "point_x",point_x
            #print "new_point_x",new_point_x 
            point_y = point[1] 
            new_point_y = new_point.y + m.geo_reference.get_yllcorner()
            
            self.assertTrue(point_x == new_point_x, ' failed')
            self.assertTrue(point_y == new_point_y, ' failed')


    def test_add_hole_from_polygon(self):
        x=-500
        y=-1000
        m=Mesh(geo_reference=Geo_reference(56,x,y))

        # These are the absolute values
        polygon_absolute = [[0,0],[1,0],[1,1],[0,1]]
        
        x_p = -10
        y_p = -40
        geo_ref_poly = Geo_reference(56, x_p, y_p)
        polygon = geo_ref_poly.change_points_geo_ref(polygon_absolute)
        
        poly_point = m.add_hole_from_polygon(polygon,
                                               {'tagin':[0,1],'bom':[2]},
                                               geo_reference=geo_ref_poly)
        # poly_point values are relative to the mesh geo-ref
        # make them absolute
        #print "poly_point.x+x",poly_point.x+x
        #print "poly_point.y+y",poly_point.y+y
        #print "polygon_absolute", polygon_absolute 
        self.assertTrue(is_inside_polygon([poly_point.x+x,poly_point.y+y],
                                       polygon_absolute, closed = False),
                        'FAILED!')
               
        self.assertTrue(len(m.holes)==1,
                        'FAILED!')
        segs = m.getUserSegments()
        self.assertTrue(len(segs)==4,
                        'FAILED!')
        self.assertTrue(len(m.userVertices)==4,
                        'FAILED!') 
        self.assertTrue(segs[0].tag=='tagin',
                        'FAILED!')  
        self.assertTrue(segs[1].tag=='tagin',
                        'FAILED!') 
         
        self.assertTrue(segs[2].tag=='bom',
                        'FAILED!') 
        self.assertTrue(segs[3].tag=='interior',
                        'FAILED!')
        verts = m.getUserVertices()
        #print "User verts",verts
        #print 'polygon',polygon
        #vert values are relative
        for point,new_point in zip(polygon,verts):
            point_x = point[0] + geo_ref_poly.get_xllcorner()
            new_point_x = new_point.x + m.geo_reference.get_xllcorner()
            point_y = point[1] + geo_ref_poly.get_yllcorner()
            #print "new_point.y",new_point.y
            #print "m.geo_ref.get_yllcorner()",m.geo_reference.get_yllcorner() 
            new_point_y = new_point.y + m.geo_reference.get_yllcorner()
            #print "point_y",point_y
            #print "new_point_y",new_point_y 
            
            self.assertTrue(point_x == new_point_x, ' failed')
            self.assertTrue(point_y == new_point_y, ' failed')




    def test_add_hole_from_polygon_none_tag(self):
        x=-500
        y=-1000
        m=Mesh(geo_reference=Geo_reference(56,x,y))

        # These are the absolute values
        polygon_absolute = [[0,0],[1,0],[1,1],[0,1]]
        
        x_p = -10
        y_p = -40
        geo_ref_poly = Geo_reference(56, x_p, y_p)
        polygon = geo_ref_poly.change_points_geo_ref(polygon_absolute)
        
        poly_point = m.add_hole_from_polygon(polygon,
                                               None,
                                               geo_reference=geo_ref_poly)
        # poly_point values are relative to the mesh geo-ref
        # make them absolute
        #print "poly_point.x+x",poly_point.x+x
        #print "poly_point.y+y",poly_point.y+y
        #print "polygon_absolute", polygon_absolute 
        self.assertTrue(is_inside_polygon([poly_point.x+x,poly_point.y+y],
                                       polygon_absolute, closed = False),
                        'FAILED!')
               
        self.assertTrue(len(m.holes)==1,
                        'FAILED!')
        segs = m.getUserSegments()
        self.assertTrue(len(segs)==4,
                        'FAILED!')
        self.assertTrue(len(m.userVertices)==4,
                        'FAILED!')
        
        self.assertTrue(segs[0].tag=='interior',
                        'FAILED!')  
        self.assertTrue(segs[1].tag=='interior',
                        'FAILED!') 
         
        self.assertTrue(segs[2].tag=='interior',
                        'FAILED!') 
        self.assertTrue(segs[3].tag=='interior',
                        'FAILED!')
        verts = m.getUserVertices()
        #print "User verts",verts
        #print 'polygon',polygon
        #vert values are relative
        for point,new_point in zip(polygon,verts):
            point_x = point[0] + geo_ref_poly.get_xllcorner()
            new_point_x = new_point.x + m.geo_reference.get_xllcorner()
            point_y = point[1] + geo_ref_poly.get_yllcorner()
            #print "new_point.y",new_point.y
            #print "m.geo_ref.get_yllcorner()",m.geo_reference.get_yllcorner() 
            new_point_y = new_point.y + m.geo_reference.get_yllcorner()
            #print "point_y",point_y
            #print "new_point_y",new_point_y 
            
            self.assertTrue(point_x == new_point_x, ' failed')
            self.assertTrue(point_y == new_point_y, ' failed')            

    def test_add_circle(self):
        x=-500
        y=-1000
        m=Mesh(geo_reference=Geo_reference(56,x,y))

        # These are the absolute values
        tag = 'hey'
        segment_count = 104
        radius = 30
        circle_center_absolute = [100,80]        
        x_p = -.666
        y_p = -.777
        geo_ref_poly = Geo_reference(56, x_p, y_p)
        circle_center = \
                geo_ref_poly.change_points_geo_ref(circle_center_absolute)
        circle_center = circle_center[0] #make a list of lists a list
        poly_point = m.add_circle(circle_center, radius, segment_count,
                                  tag=tag,
                                  region=True,
                                  center_geo_reference=geo_ref_poly)
        # poly_point values are relative to the mesh geo-ref
        # make them absolute
        #print "poly_point.x+x",poly_point.x+x
        #print "polygon_absolute", polygon_absolute 
      
        
        #m.export_mesh_file("aaat.msh")
        
        self.assertTrue(len(m.regions)==1,
                        'FAILED!')
        segs = m.getUserSegments()
        self.assertTrue(len(segs)==segment_count,
                        'FAILED!')
        self.assertTrue(len(m.userVertices)==segment_count,
                        'FAILED!') 
        self.assertTrue(segs[0].tag==tag,
                        'FAILED!')  
        self.assertTrue(segs[1].tag==tag,
                        'FAILED!') 
         
        verts = m.getUserVertices()
        
        #m.export_mesh_file("aaat.msh")
        
    def NOTIMPLEMENTEDtest_auto_set_geo_reference(self):
        x=50000
        y=1000
        m=Mesh(geo_reference=Geo_reference(56,x,y))
        polygon = [[0,0],[1,0],[1,1],[0,1]]
        
        m.add_region_from_polygon(polygon,
                               {'tagin':[0,1],'bom':[2]},
                                  max_triangle_area=10)
        m.auto_set_geo_reference()
        
 
    def test_duplicat_verts_are_removed(self):
    
      
        a = Vertex ( 0.0 ,0.0)
        b = Vertex (0.0, 4.0)
        c = Vertex (4.0,4.0)
        d = Vertex (4.0,0.0)
        e = Vertex (4.0,0.0) # duplicate point
    
        s1 = Segment(a,b, tag = "50")
        s2 = Segment(b,c, tag = "40")
        s3 = Segment(c,d, tag = "30")
        s4 = Segment(d,e, tag = "no where seg")
        s5 = Segment(e,a, tag = "20")

        
        m = Mesh(userVertices=[a,b,c,d,e],
                 userSegments=[s1,s2,s3,s4,s5])

        seg = m.getUserSegments()
        points = m.getUserVertices()
        holes = m.getHoles()
        regions = m.getRegions()
        #fileName = tempfile.mktemp(".tsh")
        #fileName = "badmesh.tsh"
        #m.export_mesh_file(fileName)
        #print "***************************fileName", fileName
        #new_m = importMeshFromFile(fileName)
        #os.remove(fileName)
        
        m.generateMesh("Q", maxArea = 2000.1 )

        #m.export_mesh_file("from_test_mesh.tsh")
        seg = m.getMeshSegments()
        self.assertTrue(4==len(seg),
                        'FAILED!') 

        vert = m.getMeshVertices() 
        self.assertTrue(4==len(vert),
                        'FAILED!')
 
    def test_duplicat_verts_are_removedII(self):
    
      
        a = Vertex ( 0.0 ,0.0)
        b = Vertex (0.0, 4.0)
        c = Vertex (4.0,4.0)
        d = Vertex (4.0,0.0)
        e = Vertex (4.0,0.0) # duplicate point
        f = Vertex (49.0,0.0) # unused point
    
        s1 = Segment(a,b, tag = "50")
        s2 = Segment(b,c, tag = "40")
        s3 = Segment(c,d, tag = "30")
        s4 = Segment(d,e, tag = "no where seg")
        s5 = Segment(e,a, tag = "20")

        
        m = Mesh(userVertices=[a,b,c,d,e,f],
                 userSegments=[s1,s2,s3,s4,s5])

        seg = m.getUserSegments()
        points = m.getUserVertices()
        holes = m.getHoles()
        regions = m.getRegions()
        #fileName = tempfile.mktemp(".tsh")
        #fileName = "badmesh.tsh"
        #m.export_mesh_file(fileName)
        #print "***************************fileName", fileName
        #new_m = importMeshFromFile(fileName)
        #os.remove(fileName)
        
        m.generateMesh("Q", maxArea = 2000.1 )

        #m.export_mesh_file("from_test_mesh.tsh")
        seg = m.getMeshSegments()
        self.assertTrue(4==len(seg),
                        'FAILED!') 

        vert = m.getMeshVertices() 
        self.assertTrue(4==len(vert),
                        'FAILED!')
   
    def test_add_vertices(self):
        points_ab = [[0.1,1],[0.4,.2],[7,5],[10,5]]
        geo =  Geo_reference(56,23,21)
        points = geo.change_points_geo_ref(points_ab)
        spat = Geospatial_data(points, geo_reference=geo)
        
        geo_mesh =  Geo_reference(56,100,200)
        m = Mesh(geo_reference=geo_mesh)
        m.add_vertices(spat)

        vert = m.getUserVertices()
        #print "vert",vert 
        self.assertTrue(4==len(vert),
                        'FAILED!')
        vert= m.get_user_vertices(absolute=True)
        
        self.assertTrue(num.allclose(vert, points_ab),
                        'FAILED!')        

    
    def test_add_vertices_more(self):
        points = [[0.1,1],[0.4,.2],[7,5],[10,5]]
        #spat = Geospatial_data(points)
        
        m = Mesh()
        m.add_vertices(points)

        vert = m.getUserVertices()
        #print "vert",vert 
        self.assertTrue(4==len(vert),
                        'FAILED!')
        vert = m.get_user_vertices(absolute=True)
        
        self.assertTrue(num.alltrue(vert.flatten() ==
                                    num.array(points).flatten()),
                        'FAILED!')
    
    def test_add_verticesII(self):
        points_lat_long = [[-33,152],[-35,152],[-35,150],[-33,150]]
       
        spat = Geospatial_data(data_points=points_lat_long,
                               points_are_lats_longs=True)
        points_ab = spat.get_data_points( absolute = True)
        m = Mesh()
        m.add_vertices(spat)

        vert = m.getUserVertices()
        #print "vert",vert 
        self.assertTrue(4==len(vert),
                        'FAILED!')
        vert= m.get_user_vertices(absolute=True)
        
        self.assertTrue(num.allclose(vert, points_ab),
                        'FAILED!')

        spat = Geospatial_data(data_points=points_lat_long,
                               points_are_lats_longs=True)
        points_ab = spat.get_data_points( absolute = True)
        geo =  Geo_reference(56,400000,6000000)
        spat.set_geo_reference(geo)
        m = Mesh()
        m.add_vertices(spat)

        vert = m.getUserVertices()
        #print "vert",vert 
        self.assertTrue(4==len(vert),
                        'FAILED!')
        vert= m.get_user_vertices(absolute=True)
        
        self.assertTrue(num.allclose(vert, points_ab),
                        'FAILED!')

        #geo =  Geo_reference(56,23,21)
        #points = geo.change_points_geo_ref(points_ab)
        
    def test_get_user_vertices(self):
        points_ab = [[0.1,1],[0.4,.2],[7,5],[10,5]]
        geo =  Geo_reference(56,23,21)
        points = geo.change_points_geo_ref(points_ab)
        spat = Geospatial_data(points, geo_reference=geo)
        
        geo_mesh =  Geo_reference(56,100,200)
        m = Mesh(geo_reference=geo_mesh)
        m.add_vertices(spat)

        vert = m.getUserVertices()
        #print "vert",vert 
        self.assertTrue(4==len(vert),
                        'FAILED!')
        vert= m.get_user_vertices(absolute=True)
        self.assertTrue(num.allclose(vert, points_ab),
                        'FAILED!')
        vert= m.get_user_vertices(absolute=False)
        points_new = m.geo_reference.get_absolute(vert)
        
        self.assertTrue(num.allclose(points_ab, points_new),
                        'FAILED!')

    def mode_string_float_problems(self):
        numbers = [0.0000000001,1000000000000.0, 1e-19,1e19, 1e-25,1e30,1e40,
                   1e41,'0.00001','0.000000000000000000000000000000000001']
        numbers = [1e-21,1e-20,1e30,1e35,1e40]
        for n in numbers:
            mode = 'a' + str(n)
            #print " mode += 'a' + str(n)", mode
            
            try:
                mode = 'a' + '%20.20f' %n
            except TypeError:
                mode = 'a' + str(n)
            print("mode += 'a' + '%20.20f' %n", mode)
        #print "", mode
        
  
    
    def testgenerateMesh_calc_mesh_area(self):
        a = Vertex (0.0, 0.0)
        d = Vertex (0.0, 4.0)
        f = Vertex (4.0,0.0)

        s1 = Segment(a,d)
        s2 = Segment(d,f)
        s3 = Segment(a,f)

        r1 = Region(0.3, 0.3,tag = 1.3,maxArea = .6)
        #print r1
        m = Mesh(userVertices=[a,d,f], userSegments=[s1,s2,s3], regions=[r1] )
        
        m.generateMesh("Q", maxArea = 2.1 )
        calc_mesh_area = m.tri_mesh.calc_mesh_area()
        #print "calc_mesh_area", calc_mesh_area
        delta  = 0.0000000001
        self.assertTrue((8.0 < calc_mesh_area + delta) or
                        (8.0 > calc_mesh_area - delta),
                        'generated mesh is wrong!')
        
def list_comp(A,B):
    yes = len(A)==len(B)
    for item in A:
        if not item in B:
            yes = False
    return yes
   
################################################################################

if __name__ == "__main__":
    suite = unittest.makeSuite(meshTestCase, 'test')
    runner = unittest.TextTestRunner() #verbosity=2)
    runner.run(suite)
    

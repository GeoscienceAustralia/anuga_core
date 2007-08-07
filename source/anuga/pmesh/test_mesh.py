#!/usr/bin/env python
#
import tempfile
import unittest

#try:
from anuga.pmesh.mesh import *
#except ImportError:  
#    from mesh import *


from load_mesh.loadASCII import *
from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.geospatial_data.geospatial_data import Geospatial_data
from anuga.utilities.polygon import  is_inside_polygon ### inside_polygon

class meshTestCase(unittest.TestCase):
    def setUp(self):
        pass
    
    def tearDown(self):
        pass

    def testPointDistance(self):
        a = Point(0.0, 0.0)
        b = Point(0.0, 10.0)
        
        self.failUnless( a.DistanceToPoint(b) == 10.0,
                        'Point DistanceToPoint is wrong!')
    
    def testVertexDistance(self):
        a = Vertex (0.0, 0.0)
        b = Vertex (0.0, 10.0)
        
        self.failUnless( a.DistanceToPoint(b) == 10.0,
                        'Point DistanceToPoint is wrong!')
        
    def testTriangle(self):
        a = Vertex (0.0, 0.0)
        b = Vertex (0.0, 2.0)
        c = Vertex (2.0,0.0)
        d = Vertex (0.0, 4.0)
        e = Vertex (2.0, 2.0)
        f = Vertex (4.0,0.0)
        
        t1 = Triangle(b,a,c)       
        t2 = Triangle(b,c,e)      
        t3 = Triangle(e,c,f)      
        t4 = Triangle(d,b,e)
        t2.setNeighbors(t3,t4,t1)
        
        self.failUnless( t2.neighbors[2].vertices[0] == b, 'Triangle initialisation is wrong!')
    
        
    def testSegment(self):
        a = Vertex (0.0, 0.0)
        b = Vertex (0.0, 10.0)
        s = Segment(a,b, tag = 20)     
        
        self.failUnless( s.vertices[0].DistanceToPoint(s.vertices[1]) == 10.0,
                        'vertices in a segment are wrong')
        
        self.failUnless( s.tag == 20.0,
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
        self.failUnless(mesh.userSegments[0] == s3,
                        'Bad segment. ')        
        self.failUnless(len(mesh.userSegments) ==1,
                        'Segments not deleted.')
        self.failUnless(len(mesh.userVertices) == 2,
                        'Vertex not deleted.')
       
  
    # FIXME add test for minAngle    
    def testgenerateMesh(self):
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

        #print m

        #m.plotMeshTriangle()

        result = 1.414214
        delta  = 0.00001
        
        self.failUnless((m.meshTriangles[1].vertices[0].x < result + delta) or
                        (m.meshTriangles[1].vertices[0].x > result - delta),
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

        self.failUnless(len(m.meshTriangles) == 2, 
                        'test_regionalMaxArea 1:generated mesh is wrong!')
        
        ## Another test case
        r1 = Region(3, 1,tag = 1.3)
        r2 = Region(1, 3,tag = 1.3, maxArea = 8)
        m = Mesh(userVertices=[v0,v1,v2,v3], userSegments=[s1,s2,s3,s4,s5],
                 regions=[r1,r2] )
        m.generateMesh("Q", maxArea = 36 )
        
        self.failUnless(len(m.meshTriangles) >= 6,
                        'testregion_with_maxarea 2: # of tris is wrong!')    
       
               
        ## Another test case
        r1 = Region(3, 1, tag = 1.3, maxArea = 8)
        r2 = Region(1, 3, tag = 1.3, maxArea = 8)
        m = Mesh(userVertices=[v0,v1,v2,v3], userSegments=[s1,s2,s3,s4,s5],
                 regions=[r1,r2] )
        m.generateMesh("Q", maxArea = 36 )  
        #print "len(m.meshTriangles)",len(m.meshTriangles)
        
        self.failUnless(len(m.meshTriangles) >= 8,
                        'testregion_with_maxarea 3: # of tris is wrong!')
                
                
        ## Another test case
        r1 = Region(3, 1, tag = 1.3 )
        r2 = Region(1, 3, tag = 1.3, maxArea = 8)
        m = Mesh(userVertices=[v0,v1,v2,v3], userSegments=[s1,s2,s3,s4,s5],
                 regions=[r1,r2] )
        m.generateMesh("Q", maxArea = 8 ) 
        self.failUnless(len(m.meshTriangles) >= 8,
                        'testregion_with_maxarea 4: # of tris is wrong!')    

        
        ## Another test case
        r1 = Region(3, 1,tag = 1.3, maxArea = 8)
        r2 = Region(1, 3,tag = 1.3, maxArea = 8)
        m = Mesh(userVertices=[v0,v1,v2,v3], userSegments=[s1,s2,s3,s4,s5],
                 regions=[r1,r2] )
        m.generateMesh("Q", maxArea = 36,isRegionalMaxAreas = False )      
        self.failUnless(len(m.meshTriangles) == 2, 
                        'test_regionalMaxArea 5:generated mesh is wrong!')
        
        ## Another test case
        r1 = Region(3, 1,tag = 1.3, maxArea = 8)
        r2 = Region(1, 3,tag = 1.3, maxArea = 8)
        m = Mesh(userVertices=[v0,v1,v2,v3], userSegments=[s1,s2,s3,s4,s5],
                 regions=[r1,r2] )
        m.generateMesh("Q",isRegionalMaxAreas = False )
        self.failUnless(len(m.meshTriangles) == 2, 
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

        self.failUnless(len(m.meshTriangles) == 2, 
                        'test_regionalMaxArea 1:generated mesh is wrong!')
        
        ## Another test case
        r1 = Region(3, 1,tag = 1.3)
        r2 = Region(1, 3,tag = 1.3, maxArea = 8)
        m = Mesh(userVertices=[v0,v1,v2,v3], userSegments=[s1,s2,s3,s4,s5],
                 regions=[r1,r2] )
        m.generate_mesh(maximum_triangle_area=36,verbose=False)  
        
        self.failUnless(len(m.meshTriangles) >= 6,
                        'testregion_with_maxarea 2: # of tris is wrong!')    
               
        ## Another test case
        r1 = Region(3, 1, tag = 1.3, maxArea = 8)
        r2 = Region(1, 3, tag = 1.3, maxArea = 8)
        m = Mesh(userVertices=[v0,v1,v2,v3], userSegments=[s1,s2,s3,s4,s5],
                 regions=[r1,r2] )
        m.generate_mesh(maximum_triangle_area=36,verbose=False)         
        #print "len(m.meshTriangles)",len(m.meshTriangles)
        
        self.failUnless(len(m.meshTriangles) >= 8,
                        'testregion_with_maxarea 3: # of tris is wrong!')
                         
        ## Another test case
        r1 = Region(3, 1, tag = 1.3 )
        r2 = Region(1, 3, tag = 1.3, maxArea = 8)
        m = Mesh(userVertices=[v0,v1,v2,v3], userSegments=[s1,s2,s3,s4,s5],
                 regions=[r1,r2] )
        m.generate_mesh(maximum_triangle_area=8,verbose=False)    
        self.failUnless(len(m.meshTriangles) >= 8,
                        'testregion_with_maxarea 4: # of tris is wrong!')    

        ## Another test case r1 = Region(3, 1,tag = 1.3, maxArea = 8)
        r2 = Region(1, 3,tag = 1.3, maxArea = 8)
        m = Mesh(userVertices=[v0,v1,v2,v3],
        userSegments=[s1,s2,s3,s4,s5], regions=[r1,r2] )
        m.generate_mesh(verbose=False)
        #print "en(m.meshTriangles)", len(m.meshTriangles)
        self.failUnless(len(m.meshTriangles) >= 8,
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
        self.failUnless(not(s2 in mesh.userSegments),
                        'Bad segment. ')        
        self.failUnless(len(mesh.userSegments) ==2,
                        'Segments not deleted.')
        self.failUnless(len(mesh.userVertices) == 3,
                        'Vertex deleted, instead of segment.')

    def testTriangleArea(self):
        a = Vertex (10.0, 10.0)
        b = Vertex (10.0, 20.0)
        c = Vertex (20.0,10.0)
        
        d = Vertex (-20.0, 0.0)
        e = Vertex (-20.0, -20.0)
        f = Vertex (20.0,-20.0)
        
        t1 = Triangle(b,a,c)      
        t2 = Triangle(e,d,f)
        
#         print "t1", t1
#         print "t1 area ", t1.calcArea()
#         print "t2", t2
#         print "t2 area ", t2.calcArea()
        self.failUnless( t1.calcArea() == 50 and t2.calcArea() == 400, 'Triangle area is wrong!')
    def testisUserSegmentNew (self):
        mesh = Mesh()
        a = mesh.addUserVertex(0.0, 0.0)
        b = mesh.addUserVertex (0.0, 2.0)
        c = mesh.addUserVertex (2.0,0.0)
        d = mesh.addUserVertex (2.0,3.0)
        
        s1 = mesh.addUserSegment(a,b)
        s2 = mesh.addUserSegment(a,c)
        s3 = mesh.addUserSegment(c,b)

        self.failUnless(mesh.isUserSegmentNew(a,d) ,
                        'Segment should be new. ')
        self.failUnless(not(mesh.isUserSegmentNew(a,b)) ,
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

        self.failUnless(mesh.representedUserSegment(a,d) == None,
                        'Segment should be new. ')
        self.failUnless(mesh.representedUserSegment(a,b) == s1 ,
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
        self.failUnless(len(m.getUserSegments()) == 4 ,
                        'userSegments is wrong!')
     
        m.auto_segment()
        self.failUnless(len(m.getUserSegments()) == 4 ,
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
        self.failUnless(not (s3 == None) ,
                        'userSegments is wrong!')

        
        s6 = m.representedAlphaUserSegment(p1,p4)       
        self.failUnless(not (s6 == None) ,
                        'userSegments is wrong!')
        
        # remove a segment, add a point, auto_segment
        m.alphaUserSegments.remove(s3)
        p6 = Vertex (1.0, 2.0)
        m.userVertices.append(p6)
        
        m.auto_segment()
        
        s1_now = m.representedUserSegment(p3,p2)
        self.failUnless(s1_now == s1 ,
                        'userSegments is wrong!')
        
        s2_now = m.representedUserSegment(p5,p4)       
        self.failUnless(s2_now == s2 ,
                        'userSegments is wrong!')
        
        s3 = m.representedAlphaUserSegment(p3,p6)       
        self.failUnless(not (s3 == None) ,
                        'userSegments is wrong!')
        
        s4 = m.representedAlphaUserSegment(p3,p6)       
        self.failUnless(not (s4 == None) ,
                        'userSegments is wrong!')
        
        s5 = m.representedAlphaUserSegment(p4,p6)       
        self.failUnless(s5 == None ,
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
        #print m
        Triangulation =  m.getTriangulation()
        #print Triangulation[0].attribute
        #print Triangulation[1].attribute 
        #print Triangulation[2].attribute 
        #print Triangulation[3].attribute 
        #print Triangulation[4].attribute
       
        self.failUnless(Triangulation[0].attribute == "" and
                        Triangulation[1].attribute == "22" and
                        Triangulation[2].attribute == "" and
                        Triangulation[3].attribute == "11" and
                        Triangulation[4].attribute == "22" ,
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

        vert = m.getMeshVertices()
        
        self.failUnless(vert[0].attributes == [12.0, 2.0] and
                        vert[1].attributes == [9.0, 7.0] and
                        vert[2].attributes == [14.0,3.0] and
                        vert[3].attributes == [12.232233047033631, 4.4142135623730949] and
                        vert[4].attributes == [13.0, 2.5] ,
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

        vert = m.getMeshVertices()
        self.failUnless(vert[0].attributes == [] and
                        vert[1].attributes == [] and
                        vert[2].attributes == [] and
                        vert[3].attributes == [] and
                        vert[4].attributes == [],
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
        seg = m.getMeshSegments()
        #print "seg",seg
        #print "seg[0].tag"
        #print seg[0].tag
        #print "seg[0].tag"
        
        self.failUnless(seg[0].tag == 5 and
                        seg[1].tag == 7 and
                        seg[2].tag == 9 and
                        seg[3].tag == 7 and
                        seg[4].tag == 9,
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

        seg = m.getMeshSegments()
        self.failUnless(seg[0].tag == "exterior" and
                        seg[1].tag == "exterior" and
                        seg[2].tag == "exterior" and
                        seg[3].tag == "" and
                        seg[4].tag == "exterior",
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
        #self.failUnless(lFile[0] == "5 0 # <vertex #> <x> <y> [attributes] ...Triangulation Vertices..."
          #              ,'Ascii file is wrong, vertex title')
        self.failUnless(lFile[1] == "0 0.0 0.0 " and #1.1 " and
                        lFile[2] == "1 0.0 4.0 " and #1.2 " and
                        lFile[3] == "2 4.0 0.0 " and #1.3 " and
                        lFile[4] == "3 1.0 1.0 " and #1.4 " and
                        lFile[5] == "4 2.0 2.0 "  #1.25 " 
                        ,
                        'Ascii file is wrong, vertex')
        
        #self.failUnless(lFile[6] == "# attribute column titles ...Triangulation Vertex Titles..."
          #              ,'Ascii file is wrong, attribute column title')
        self.failUnless(lFile[8] == "0 3 2 4 -1 2 3  " and
                        lFile[9] == "1 1 0 3 3 2 -1  " and
                        lFile[10] == "2 3 4 1 -1 1 0  " and
                        lFile[11] == "3 2 3 0 1 -1 0  "
                        ,
                        'Ascii file is wrong, triangle') 

        self.failUnless( lFile[13] == "0 0 1 exterior" and
                        lFile[14] == "1 1 4 exterior" and
                        lFile[15] == "2 2 0 exterior" and
                        lFile[16] == "3 0 3 " and
                        lFile[17] == "4 4 2 exterior" ,
                        'Ascii file is wrong, segment')
        
       # self.failUnless(lFile[18] == '4 0 # <vertex #> <x> <y> [attributes] ...Mesh Vertices...',
        #                'Ascii file is wrong, Mesh Vertices Title')
        
        self.failUnless(lFile[19] == '0 0.0 0.0 ' and #1.1 ' and
                        lFile[20] == '1 0.0 4.0 ' and #1.2 ' and
                        lFile[21] == '2 4.0 0.0 ' and #1.3 ' and
                        lFile[22] == '3 1.0 1.0 ' #1.4 '
                        ,
                        'Ascii file is wrong, Mesh Vertices II')
        
        self.failUnless(lFile[24] == '0 0 1 ' and
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
        self.failUnless(lFile[0] == "5 0 # <# of verts> <# of vert attributes>, next lines <vertex #> <x> <y> [attributes] ...Triangulation Vertices..."
                        ,
                        'Ascii file is wrong, vertex title')
        self.failUnless(lFile[1] == "0 0.0 0.0 " and #1.1 " and
                        lFile[2] == "1 0.0 4.0 " and #1.2 " and
                        lFile[3] == "2 4.0 0.0 " and #1.3 " and
                        lFile[4] == "3 1.0 1.0 " and #1.4 " and
                        lFile[5] == "4 2.0 2.0 "  #1.25 " 
                        ,
                        'Ascii file is wrong, vertex')
        
        self.failUnless(lFile[6] == "# attribute column titles ...Triangulation Vertex Titles..."
                        ,
                        'Ascii file is wrong, attribute column title')
        self.failUnless(lFile[7] == "4 # <# of triangles>, next lines <triangle #> [<vertex #>] [<neigbouring triangle #>] [attribute of region] ...Triangulation Triangles..." and
                        lFile[8] == "0 3 2 4 -1 2 3  " and
                        lFile[9] == "1 1 0 3 3 2 -1  " and
                        lFile[10] == "2 3 4 1 -1 1 0  " and
                        lFile[11] == "3 2 3 0 1 -1 0  "
                        ,
                        'Ascii file is wrong, triangle') 

        self.failUnless(lFile[12] == "5 # <# of segments>, next lines <segment #> <vertex #>  <vertex #> [boundary tag] ...Triangulation Segments..." and
                        lFile[13] == "0 0 1 exterior" and
                        lFile[14] == "1 1 4 exterior" and
                        lFile[15] == "2 2 0 exterior" and
                        lFile[16] == "3 0 3 " and
                        lFile[17] == "4 4 2 exterior" ,
                        'Ascii file is wrong, segment')
        
        self.failUnless(lFile[18] == '4 0 # <# of verts> <# of vert attributes>, next lines <vertex #> <x> <y> [attributes] ...Mesh Vertices...',
                        'Ascii file is wrong, Mesh Vertices Title')
        
        self.failUnless(lFile[19] == '0 0.0 0.0 ' and #1.1 ' and
                        lFile[20] == '1 0.0 4.0 ' and #1.2 ' and
                        lFile[21] == '2 4.0 0.0 ' and #1.3 ' and
                        lFile[22] == '3 1.0 1.0 ' #1.4 '
                        ,
                        'Ascii file is wrong, Mesh Vertices II')
        
        self.failUnless(lFile[23] == '4 # <# of segments>, next lines <segment #> <vertex #>  <vertex #> [boundary tag] ...Mesh Segments...' and
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
         
        self.failUnless(v1 in m.userVertices,
                        'test_thinoutVertices, test 1 failed')
        self.failUnless(v3 in m.userVertices,
                        'test_thinoutVertices, test 2 failed')
        self.failUnless(v4 in m.userVertices,
                        'test_thinoutVertices, test 3 failed')
        self.failUnless(v6 in m.userVertices,
                        'test_thinoutVertices, test 4 failed')
        self.failUnless(v7 in m.userVertices,
                        'test_thinoutVertices, test 5 failed')
        self.failUnless(v9 in m.userVertices,
                        'test_thinoutVertices, test 6 failed')
        self.failUnless(v5 not in m.userVertices,
                        'test_thinoutVertices, test 7 failed')
        self.failUnless(v2 not in m.userVertices,
                        'test_thinoutVertices, test 8 failed')
        self.failUnless(v8 not in m.userVertices,
                        'test_thinoutVertices, test 9 failed')

    def test_same_x_y(self):
        v = Point(7,8)
        f = Point(7,8)
        f.same_x_y(v)

        self.failUnless(f.same_x_y(v),
                        'same_x_y True failed')
        e = Point(7,9)
        self.failUnless(not f.same_x_y(e),
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
        geo = Geo_reference(8.9,8.9,65)
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
        self.failUnless(0 == m.__cmp__(m_returned),
                        'loading and saving of a mesh failed')
        # do this when .msh supports geo_refs
        #self.failUnless(m.geo_reference == m_returned.geo_reference,
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
        geo = Geo_reference(65,8.9,8.9)
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
        self.failUnless(0 == m.__cmp__(m_returned),
                        'loading and saving of a mesh failed')
        self.failUnless(m.geo_reference == m_returned.geo_reference,
                        'loading and saving of a mesh geo refs failed')

    def test_normaliseMesh(self):
        
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
        self.failUnless(attmin == 0.0 and attmax == 1.0,
                        'normalise failed')
        self.failUnless(xmin == 0.0 and ymin == 0.0 and xmax == 0.5 and ymax == 1.0,
                        'normalise failed')
        m.normaliseMesh(200,-100,5)
        [xmin, ymin, xmax, ymax] = m.boxsize()
        [attmin, attmax] = m.maxMinVertAtt(0)
        self.failUnless(attmin == 0.0 and attmax == 5.0,
                        'normalise failed')
        self.failUnless(xmin == -100.0 and ymin == -100.0 and xmax == 0.0 and ymax == 100.0,
                        'normalise failed')
        
    def test_exportASCIIsegmentoutlinefile(self):
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
        m.meshVertices = []
        m.meshTriangles = []
        m.meshSegments = []
        m.userVertices=[a,b,c]
        #print "mesh ***************dsg*", m
        #print "(m.__cmp__(m_returned)", m.__cmp__(m_returned) 
        self.failUnless(0 == m.__cmp__(m),
                        'test_exportASCIIsegmentoutlinefile:loading and saving of a mesh failed')
        # Having problems with this on linux.
        #The ordering of the dictionary values wasn't the same as the windows
        #returned value (verts.values())
        #self.failUnless(0 == m.__cmp__(m_returned),
        #                'test_exportASCIIsegmentoutlinefile:loading and saving of a mesh failed')
        
        self.failUnless(3 == len(m_returned.userVertices),
                        'segmentoutlinefile:IO of a mesh failed')
        self.failUnless(len(m.userSegments) == len(m_returned.userSegments),
                        'segmentoutlinefile:IO of a mesh failed')
        for i in range(len(m.userSegments)):
            self.failUnless(m.userSegments[i].vertices[0].x ==
                            m_returned.userSegments[i].vertices[0].x,
                        'loading and saving of a mesh outline fialed')
            self.failUnless(m.userSegments[i].vertices[0].y ==
                            m_returned.userSegments[i].vertices[0].y,
                        'loading and saving of a mesh outline fialed')
            self.failUnless(m.userSegments[i].vertices[1].x ==
                            m_returned.userSegments[i].vertices[1].x,
                        'loading and saving of a mesh outline fialed')
            self.failUnless(m.userSegments[i].vertices[1].y ==
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
        self.failUnless(0 == m.__cmp__(m),
                        'loading and saving of a mesh failed')

        self.failUnless(5 == len(m_returned.userVertices),
                        'segmentoutlinefile:IO of a mesh failed')
        self.failUnless(len(m.userSegments) == len(m_returned.userSegments),
                        'segmentoutlinefile:IO of a mesh failed')
        for i in range(len(m.userSegments)):
            self.failUnless(m.userSegments[i].vertices[0].x ==
                            m_returned.userSegments[i].vertices[0].x,
                        'loading and saving of a mesh outline fialed')
            self.failUnless(m.userSegments[i].vertices[0].y ==
                            m_returned.userSegments[i].vertices[0].y,
                        'loading and saving of a mesh outline fialed')
            self.failUnless(m.userSegments[i].vertices[1].x ==
                            m_returned.userSegments[i].vertices[1].x,
                        'loading and saving of a mesh outline fialed')
            self.failUnless(m.userSegments[i].vertices[1].y ==
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
        self.failUnless(m.userVertices[0].x == 1.0,
                        'loadxy, test 1 failed')
        self.failUnless(m.userVertices[0].y == 0.0,
                        'loadxy, test 2 failed')
        #self.failUnless(m.userVertices[0].attributes == [10.0,0.0],
        #                'loadxy, test 2.2 failed')
        self.failUnless(m.userVertices[1].x == 0.0,
                        'loadxy, test 3 failed')
        self.failUnless(m.userVertices[1].y == 1.0,
                        'loadxy, test 4 failed')
        #self.failUnless(m.userVertices[1].attributes == [0.0,10.0],
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
        self.failUnless(lFile[0] == "x,y," and
                        lFile[1] == "0.0,0.0" and
                        lFile[2] == "0.0,3.0" and
                        lFile[3] == "3.0,3.0" 
                        ,
                        'exported Ascii csv file is wrong')
        self.failUnless(lFile[4] == "1.0,2.0" and
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
        
        self.failUnless(lFile[0] == "x,y," and
                        lFile[1] == "0.0,0.0" and
                        lFile[2] == "0.0,3.0" and
                        lFile[3] == "3.0,3.0" and
                        lFile[4] == "1.0,2.0"
                        ,
                        'exported Ascii csv file is wrong')
     
    def to_be_test_lone_vert_in_mesh_gen_c_layer(self):
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
        self.failUnless(lFile[0] == "x,y" and
                        lFile[1] == "0,0" and
                        lFile[2] == "0,3" and
                        lFile[3] == "3,3" 
                        ,
                        'exported Ascii csv file is wrong')
        self.failUnless(lFile[4] == "1,2" and
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
        
        self.failUnless(lFile[0] == "x,y" and
                        lFile[1] == "0.0,0.0" and
                        lFile[2] == "0.0,3.0" and
                        lFile[3] == "3.0,3.0" and
                        lFile[4] == "1.0,2.0"
                        ,
                        'exported Ascii csv file is wrong')
        
    def NOT_test_exportPointsFilefile2(self):
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
        self.failUnless(lFile[0] == "" 
                        ,
                        'exported Ascii csv file is wrong')
        
    def test_strings2ints(self):
        list = ["sea","river inlet","","sea","","moat"]
        preset = ["moat", "internal boundary"]
        [intlist, converter] = segment_strings2ints(list,preset )
        self.failUnless(intlist == [2,3 ,0 ,2 ,0 ,0 ]
                        ,
                        'test_strings2ints produces bad intlist')
        self.failUnless(converter == ['moat', 'internal boundary',
                                      'sea', 'river inlet']
                        ,
                        'test_strings2ints produces bad converter')
        
    def test_ints2strings(self):
        list = ["internal boundary","sea","river inlet",
            "","sea","","moat","internal boundary"]
        outlist = ['internal boundary', 'sea', 'river inlet', 'moat',
                   'sea', 'moat', 'moat', 'internal boundary']
        preset = ["moat", "internal boundary"]
        [intlist, converter] = segment_strings2ints(list,preset )
        newlist = segment_ints2strings(intlist, converter)
        self.failUnless(outlist == newlist
                        ,
                        'test_strings2ints produces bad intlist')
        self.failUnless(converter == ['moat', 'internal boundary',
                                      'sea', 'river inlet']
                        ,
                        'test_strings2ints produces bad converter')
        
    def test_ints2strings2(self):
        list = ["","",""]
        preset = ["moat", "internal boundary"]
        [intlist, converter] = segment_strings2ints(list,preset )
        newlist = segment_ints2strings(intlist, converter)
        outlist = ['moat', 'moat', 'moat']
        self.failUnless(outlist == newlist
                        ,
                        'test_strings2ints produces bad intlist')
        self.failUnless(converter == ['moat', 'internal boundary']
                        ,
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
        
         
        self.failUnless(UserVerts == inputVerts_noDups,
                            'duplicate verts not removed')
        #for userVert, inputVert in map(None, UserVerts, inputVerts_noDups): 
        #    self.failUnless(userVert.x == inputVert.x,
        #                    'x duplicate verts not removed')
        #    self.failUnless(userVert.y == inputVert.y,
        #                    'y duplicate verts not removed')

        
    def test_ungenerateFileLoading(self):
        
        fileName = tempfile.mktemp(".txt")
        file = open(fileName,"w")
        file.write("         1       ??      ??\n\
       0.0       0.0\n\
       1.0       0.0\n\
       1.0       1.0\n\
       0.0       1.0\n\
       0.0       0.0\n\
END\n\
         2      ?? ??\n\
       10.0       10.0\n\
       10.0       20.0\n\
       20.0       20.0\n\
       10.0       10.0\n\
END\n\
END\n")
        file.close()
        
        
        a = Vertex (0.0, 0.0) #, attributes = [1.1])
        b = Vertex (0.0, 40.0) #, attributes = [1.2])
        c = Vertex (40.0,40.0) #, attributes = [1.3])
        d = Vertex (40.0,0.0) #, attributes = [1.4])
    
        s1 = Segment(a,b)
        s2 = Segment(b,c)
        s3 = Segment(c,d)
        s4 = Segment(d,a)
     
        m = Mesh(userVertices=[a,b,c,d], userSegments=[s1,s2,s3,s4])
        dict = importUngenerateFile(fileName)
        #os.remove(fileName)

        tag = "DSG"
        Segment.set_default_tag(tag)
        m.addVertsSegs(dict)

        # have to reset this , since it's a class attribute
        Segment.set_default_tag("")
            
        self.failUnless(len(m.userSegments) ==11,
                        'Wrong segment list length.')
        self.failUnless(len(m.userVertices) == 11,
                        'Wrong vertex list length.')
        self.failUnless(m.userSegments[10].vertices[0] == m.userVertices[10],
                        'bad vertex on segment.')
        self.failUnless(m.userSegments[10].vertices[1] == m.userVertices[8],
                        'Bad segment.')
        self.failUnless(m.userSegments[10].tag == tag,
                        'wrong tag.')

        ## let's test the method
        a = Vertex (0.0, 0.0) #, attributes = [1.1])
        b = Vertex (0.0, 40.0) #, attributes = [1.2])
        c = Vertex (40.0,40.0) #, attributes = [1.3])
        d = Vertex (40.0,0.0) #, attributes = [1.4])
    
        s1 = Segment(a,b)
        s2 = Segment(b,c)
        s3 = Segment(c,d)
        s4 = Segment(d,a)
     
        m = Mesh(userVertices=[a,b,c,d], userSegments=[s1,s2,s3,s4])

        tag = "DSG"        
        initial_tag = "PIG"
        Segment.set_default_tag(initial_tag)
        m.import_ungenerate_file(fileName, tag=tag)

        os.remove(fileName)

        self.failUnless(Segment.get_default_tag() == initial_tag,
                        'Wrong segment list length.')
        

        # have to reset this , since it's a class attribute
        Segment.set_default_tag("")
            
        self.failUnless(len(m.userSegments) ==11,
                        'Wrong segment list length.')
        self.failUnless(len(m.userVertices) == 11,
                        'Wrong vertex list length.')
        self.failUnless(m.userSegments[10].vertices[0] == m.userVertices[10],
                        'bad vertex on segment.')
        self.failUnless(m.userSegments[10].vertices[1] == m.userVertices[8],
                        'Bad segment.')
        self.failUnless(m.userSegments[10].tag == tag,
                        'wrong tag.')
        
    def test_ungenerateFileLoadingII(self):
        
        fileName = tempfile.mktemp(".txt")
        file = open(fileName,"w")
        file.write("         1       ??      ??\n\
       0.0       0.0\n\
       1.0       0.0\n\
       1.0       1.0\n\
       0.0       1.0\n\
       0.0       0.0\n\
END\n\
         2      ?? ??\n\
       10.0       10.0\n\
       10.0       20.0\n\
       20.0       20.0\n\
END\n\
END\n")
        file.close()
        
        
        a = Vertex (0.0, 0.0) #, attributes = [1.1])
        b = Vertex (0.0, 40.0) #, attributes = [1.2])
        c = Vertex (40.0,40.0) #, attributes = [1.3])
        d = Vertex (40.0,0.0) #, attributes = [1.4])
    
        s1 = Segment(a,b)
        s2 = Segment(b,c)
        s3 = Segment(c,d)
        s4 = Segment(d,a)
     
        m = Mesh(userVertices=[a,b,c,d], userSegments=[s1,s2,s3,s4])
        dict = importUngenerateFile(fileName)
        #os.remove(fileName)

        tag = "DSG"
        Segment.set_default_tag(tag)
        m.addVertsSegs(dict)

        self.failUnless(len(m.userSegments) ==10,
                        'Wrong segment list length.')
        self.failUnless(len(m.userVertices) == 11,
                        'Wrong vertex list length.')

        # Test the method
        a = Vertex (0.0, 0.0) #, attributes = [1.1])
        b = Vertex (0.0, 40.0) #, attributes = [1.2])
        c = Vertex (40.0,40.0) #, attributes = [1.3])
        d = Vertex (40.0,0.0) #, attributes = [1.4])
    
        s1 = Segment(a,b)
        s2 = Segment(b,c)
        s3 = Segment(c,d)
        s4 = Segment(d,a)
     
        m = Mesh(userVertices=[a,b,c,d], userSegments=[s1,s2,s3,s4])
        tag = "DSG"        
        initial_tag = "PIG"
        Segment.set_default_tag(initial_tag)
        m.import_ungenerate_file(fileName, tag=tag)

        os.remove(fileName)

        self.failUnless(Segment.get_default_tag() == initial_tag,
                        'Wrong segment list length.')
        
        
        self.failUnless(len(m.userSegments) ==10,
                        'Wrong segment list length.')
        self.failUnless(len(m.userVertices) == 11,
                        'Wrong vertex list length.')
        
        # have to reset this , since it's a class attribute
        Segment.set_default_tag("")
        
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

        
        self.failUnless(len(m.userSegments) ==2,
                        'Wrong segment list length.')
        self.failUnless(len(m.userVertices) == 3,
                        'Wrong vertex list length.')
        self.failUnless(m.userSegments[0].tag =='food',
                        'Wrong segment tag length.')
        self.failUnless(m.userSegments[1].tag =='do-op',
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

        
        self.failUnless(m.userSegments[5].vertices[0].y == 3,
                        'Wrong vertex connected.')
        self.failUnless(m.userSegments[5].vertices[1].y == 1,
                        'Wrong vertex connected.')
            
    def test_add_points_and_segments(self):
        m = Mesh()
        Segment.set_default_tag("food")
        dict = {}
        points =  [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]
        segments = [[0, 1], [1, 2]]
        segment_tags = {'do-op':[1]}
        m.add_points_and_segments(points,
                                    segments, segment_tags)
        # have to reset this , since it's a class attribute
        Segment.set_default_tag("")

        
        self.failUnless(len(m.userSegments) ==2,
                        'Wrong segment list length.')
        self.failUnless(len(m.userVertices) == 3,
                        'Wrong vertex list length.')
        self.failUnless(m.userSegments[0].tag =='food',
                        'Wrong segment tag length.')
        self.failUnless(m.userSegments[1].tag =='do-op',
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
        
        self.failUnless( new_m == m,
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
        titles = ['ele','friction']
        m.attributeTitles = titles
        m.generateMesh("Qa2.1")

        seg = m.getMeshSegments()
        verts = m.getMeshVertices()
        dict = m.Mesh2IOTriangulationDict()
        #print "dict",dict 
        
        self.failUnless( dict['vertex_attribute_titles'] == titles,
                         'test_Mesh2IOTriangulationDict failed. test 1')
        answer = [a_xy,[0.0, 4.0],[4.0,0.0], [1.0,1.0], [2.0,2.0]]
        #print "answer",answer
        #print "dict['vertices']",dict['vertices']
        
        self.failUnless( dict['vertices'] == answer,
                         'test_Mesh2IOTriangulationDict failed. test 2')

        
        for pimport,pactual,pimpatt in map(None,dict['vertices'],
                                           verts,dict['vertex_attributes']):
            #assert all_close( pimport, (pactual.x,pactual.y))
            self.failUnless( pimport == [pactual.x,pactual.y],
                        'test_Mesh2IOTriangulationDict failed. test 2.1')
            self.failUnless( pimpatt == pactual.attributes,
                        'test_Mesh2IOTriangulationDict failed. test 2.2')
        self.failUnless( dict['segments'][0] == [0,1],
                        'test_Mesh2IOTriangulationDict failed. test 3')
        for segimp,segactual in map(None,dict['segment_tags'],seg):
            self.failUnless( segimp == segactual.tag,
                        'test_Mesh2IOTriangulationDict failed. test 4')
        #print " dict['triangles'][0]", dict['triangles'][0] 
        self.failUnless( dict['triangles'][0] == [3,2,4],
                        'test_Mesh2IOTriangulationDict failed. test 5')
        self.failUnless( dict['triangle_neighbors'][0] == [-1,2,3],
                        'test_Mesh2IOTriangulationDict failed. test 6')
        self.failUnless( dict['triangle_tags'][0] == "1.3",
                         'test_Mesh2IOTriangulationDict failed. test 7')

  
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
        dict = m.Mesh2IODict()
        #print "dict",dict 
        
        self.failUnless( dict['vertex_attribute_titles'] == titles,
                         'test_Mesh2IOTriangulationDict failed. test 1')
        answer = [a_xy,[0.0, 4.0],[4.0,0.0], [1.0,1.0], [2.0,2.0]]
        #print "answer",answer
        #print "dict['vertices']",dict['vertices']
        
        self.failUnless( dict['vertices'] == answer,
                         'test_Mesh2IOTriangulationDict failed. test 2')

        
        for pimport,pactual,pimpatt in map(None,dict['vertices'],
                                           verts,dict['vertex_attributes']):
            #assert all_close( pimport, (pactual.x,pactual.y))
            self.failUnless( pimport == [pactual.x,pactual.y],
                        'test_Mesh2IODict failed. test 2.1')
            self.failUnless( pimpatt == pactual.attributes,
                        'test_Mesh2IODict failed. test 2.2')
        self.failUnless( dict['segments'][0] == [0,1],
                        'test_Mesh2IODict failed. test 3')
        for segimp,segactual in map(None,dict['segment_tags'],seg):
            self.failUnless( segimp == segactual.tag,
                        'test_Mesh2IODict failed. test 4')
        #print " dict['triangles'][0]", dict['triangles'][0] 
        self.failUnless( dict['triangles'][0] == [3,2,4],
                        'test_Mesh2IODict failed. test 5')
        self.failUnless( dict['triangle_neighbors'][0] == [-1,2,3],
                        'test_Mesh2IODict failed. test 6')
        self.failUnless( dict['triangle_tags'][0] == "1.3",
                         'test_Mesh2IODict failed. test 7')

        seg = m.getUserSegments()
        points = m.getUserVertices()
        holes = m.getHoles()
        regions = m.getRegions()
        
        for pimport,pactual,pimpatt in map(None,dict['points'],points,dict['point_attributes']):
            self.failUnless( pimport == [pactual.x,pactual.y],
                        'test_Mesh2IODict failed. test 1')
            self.failUnless( pimpatt == pactual.attributes,
                        'test_Mesh2IODict failed. test 1.1')
        self.failUnless( dict['outline_segments'][0] == [0,1],
                        'test_Mesh2IODict failed. test 3')
        for segimp,segactual in map(None,dict['outline_segment_tags'],seg):
            self.failUnless( segimp == segactual.tag,
                        'test_Mesh2IODict failed. test 4')
        for holeimp,holeactual in map(None,dict['holes'],holes):
            self.failUnless( holeimp == [holeactual.x,holeactual.y],
                        'test_Mesh2IODict failed. test 5')
        
        for regimp,regactual,regattimp, regmaxarea in map(None,dict['regions'],regions, dict['region_tags'], dict['region_max_areas']):
            self.failUnless( regimp == [regactual.x,regactual.y],
                        'loadASCIITestCase failed. test 6')
            self.failUnless( regattimp == regactual.getTag(),
                        'loadASCIITestCase failed. test 7')
            self.failUnless( regmaxarea == regactual.getMaxArea(),
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
        
        for pimport,pactual,pimpatt in map(None,dict['points'],points,dict['point_attributes']):
            self.failUnless( pimport == [pactual.x,pactual.y],
                        'loadASCIITestCase failed. test 1')
            self.failUnless( pimpatt == pactual.attributes,
                        'loadASCIITestCase failed. test 1.1')
        self.failUnless( dict['outline_segments'][0] == [0,1],
                        'loadASCIITestCase failed. test 3')
        for segimp,segactual in map(None,dict['outline_segment_tags'],seg):
            self.failUnless( segimp == segactual.tag,
                        'loadASCIITestCase failed. test 4')
        for holeimp,holeactual in map(None,dict['holes'],holes):
            self.failUnless( holeimp == [holeactual.x,holeactual.y],
                        'loadASCIITestCase failed. test 5')
        #for regimp,regactual in map(None,dict['regions'],regions):
         #   self.failUnless( [regimp[0],regimp[1]]==[regactual.x,regactual.y],
          #              'loadASCIITestCase failed. test 6')
           # self.failUnless( regimp[2] == regactual.getTag(),
            #            'loadASCIITestCase failed. test 7')
            #self.failUnless( regimp[3] == regactual.getMaxArea(),
             #           'loadASCIITestCase failed. test 7')

            
        for regimp,regactual,regattimp, regmaxarea in map(None,dict['regions'],regions, dict['region_tags'], dict['region_max_areas']):
            self.failUnless( regimp == [regactual.x,regactual.y],
                        'loadASCIITestCase failed. test 6')
            self.failUnless( regattimp == regactual.getTag(),
                        'loadASCIITestCase failed. test 7')
            self.failUnless( regmaxarea == regactual.getMaxArea(),
                        'loadASCIITestCase failed. test 7')


#___________beginning of Peters tests

    def test_set_stuff(self):
        """
        Documentation
        """
        #making a test mesh
        p0=[0.,2.]
        p1=[1.,2.]
        p2=[0.,1.]
        p3=[1.,1.]
        p4=[0.,0.]
        p5=[2.,0.]
        p6=[-1.,1.]
        point_list = [p0,p1,p2,p3,p4,p5,p6]

        a0=[0]
        a1=[0]
        a2=[100]
        a3=[0]
        a4=[0]
        a5=[0]
        a6=[0]
        attribute=[a0,a1,a2,a3,a4,a5,a6]
        
        t0=[0,3,1]
        t1=[0,2,3]
        t2=[2,4,3]
        t3=[4,5,3]
        t4=[1,3,5]
        t5=[2,6,4]

        n0=[4,-1,2]
        n1=[2,0,-1]
        n2=[3,1,5]
        n3=[4,2,-1]
        n4=[3,-1,0]
        n5=[-1,2,-1]

        tri_list = [t0,t1,t2,t3,t4,t5]
        n_list = [n0,n1,n2,n3,n4,n5]
        for i in range(6):
            for j in (0,1,2):
                a=attribute[tri_list[i][j]]
                tri_list[i][j]=point_list[tri_list[i][j]]
                tri_list[i][j]=Vertex(tri_list[i][j][0]\
                                      ,tri_list[i][j][1],a)
            neighbours=n_list[i]
            tri_list[i]=Triangle(tri_list[i][0],\
                                 tri_list[i][1],tri_list[i][2]\
                                ,neighbors=neighbours)

        #testing selectAll
        mesh = Mesh()
        mesh.attributeTitles=['attrib']

        mesh.meshTriangles=tri_list

        mesh.selectAllTriangles()
        A=mesh.sets[mesh.setID['All']]
        assert list_comp(tri_list,A)

       #testing threshold
        mesh = Mesh()
        mesh.attributeTitles=['attrib']

        mesh.meshTriangles=tri_list
        mesh.selectAllTriangles()
        mesh.threshold('All',min=30,max=35,attribute_name = 'attrib')
        A = [tri_list[1],tri_list[2],tri_list[5]]
        B = mesh.sets[mesh.setID['All']]
        assert list_comp(A,B)
       

        A = [tri_list[3],tri_list[2],tri_list[5]]
        assert not list_comp(A,B)

        #testing 

    def test_Discretised_Tuple_Set_rounding(self):
        #This is the hardest bit of DST

        tol = 0.1
        a=Discretised_Tuple_Set(p_rel=1,t_rel= tol)
        m = 0.541
        m_up = 0.6
        m_down = 0.5
        assert m_up == a.round_up_rel(m)
        assert m_down == a.round_down_rel(m)

        tol = 0.1
        a=Discretised_Tuple_Set(p_rel=1,t_rel = tol)
        m = 0.539
        m_up = 0.5
        m_down = 0.5
        assert m_up == a.round_up_rel(m)
        assert m_down == a.round_down_rel(m)

        tol = 0.5
        a=Discretised_Tuple_Set(p_rel=1,t_rel = tol)


        m = 0.6
        m_up = 0.7
        m_down = 0.5
        assert m_up == a.round_up_rel(m)
        assert m_down == a.round_down_rel(m)

        m = 0.599
        m_up = 0.6
        m_down = 0.5
        assert m_up == a.round_up_rel(m)
        assert m_down == a.round_down_rel(m)

    def test_Discretised_Tuple_Set_get(self):
        
        tol = 0.25
        a=Discretised_Tuple_Set(p_rel=1,t_rel = tol)
        b = (1.1,1.1)
        a.append(b)
        list = [(1.2,1.),(1.,1.),(1.,1.2),(1.2,1.2)]
        for key in list:
            assert a[key][0]==b
            assert len(a[key])==1
        
        c = (2.1,1.)
        a.append(c)
        assert a[(2.,1.)][0]==c
        assert a[(2.2,1.)][0]==c

    def test_mapped_Discretised_Tuple_Set(self):

        def map(sequence):
            return [len(sequence)]

        tol = 0.5
        a=Mapped_Discretised_Tuple_Set(map,p_rel=1,t_rel = tol)
        b = range(20)
        a.append(b)
        assert b in a[range(17)] 
        assert b in a[range(22)]

        tol = 0.01
        a=Mapped_Discretised_Tuple_Set(map,p_rel=1,t_rel = tol)
        b = range(20)
        a.append(b)
        assert b in a[range(20)] 
        assert b in a[range(19)] 
        assert not range(17) in a

#___________end of Peters tests

    def test_add_region_from_polygon(self):
        m=Mesh()
        region = m.add_region_from_polygon([[0,0],[1,0],[0,1]],
                                  max_triangle_area = 88)
        self.failUnless(len(m.regions)==1,
                        'FAILED!')
        self.failUnless(region.getMaxArea()==88,
                        'FAILED!')
        self.failUnless(len(m.getUserSegments())==3,
                        'FAILED!')
        self.failUnless(len(m.userVertices)==3,
                        'FAILED!')
       
    def test_add_region_from_polygon2(self):
        m=Mesh()
        m.add_region_from_polygon([[0,0],[1,0],[1,1],[0,1]],
                               {'tagin':[0,1],'bom':[2]},
                                  max_triangle_area=10)
        self.failUnless(len(m.regions)==1,
                        'FAILED!')
        segs = m.getUserSegments()
        self.failUnless(len(segs)==4,
                        'FAILED!')
        self.failUnless(len(m.userVertices)==4,
                        'FAILED!') 
        self.failUnless(segs[0].tag=='tagin',
                        'FAILED!')  
        self.failUnless(segs[1].tag=='tagin',
                        'FAILED!') 
         
        self.failUnless(segs[2].tag=='bom',
                        'FAILED!') 
        self.failUnless(segs[3].tag=='',
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
        self.failUnless(is_inside_polygon([poly_point.x+x,poly_point.y+y],
                                       polygon_absolute, closed = False),
                        'FAILED!')
               
        self.failUnless(len(m.regions)==1,
                        'FAILED!')
        segs = m.getUserSegments()
        self.failUnless(len(segs)==4,
                        'FAILED!')
        self.failUnless(len(m.userVertices)==4,
                        'FAILED!') 
        self.failUnless(segs[0].tag=='tagin',
                        'FAILED!')  
        self.failUnless(segs[1].tag=='tagin',
                        'FAILED!') 
         
        self.failUnless(segs[2].tag=='bom',
                        'FAILED!') 
        self.failUnless(segs[3].tag=='',
                        'FAILED!')
        verts = m.getUserVertices()
        #print "User verts",verts
        #print 'polygon',polygon
        #vert values are relative
        for point,new_point in map(None,polygon,verts):
            point_x = point[0] + geo_ref_poly.get_xllcorner()
            new_point_x = new_point.x + m.geo_reference.get_xllcorner()
            point_y = point[1] + geo_ref_poly.get_yllcorner()
            #print "new_point.y",new_point.y
            #print "m.geo_ref.get_yllcorner()",m.geo_reference.get_yllcorner() 
            new_point_y = new_point.y + m.geo_reference.get_yllcorner()
            #print "point_y",point_y
            #print "new_point_y",new_point_y 
            
            self.failUnless(point_x == new_point_x, ' failed')
            self.failUnless(point_y == new_point_y, ' failed')
            
         
    def test_add_region_from_polygon4(self):
        x=50000
        y=1000
        m=Mesh(geo_reference=Geo_reference(56,x,y))
        polygon = [[0,0],[1,0],[1,1],[0,1]]
        
        m.add_region_from_polygon(polygon,
                               {'tagin':[0,1],'bom':[2]},
                                  max_triangle_area=10)
        self.failUnless(len(m.regions)==1,
                        'FAILED!')
        segs = m.getUserSegments()
        self.failUnless(len(segs)==4,
                        'FAILED!')
        self.failUnless(len(m.userVertices)==4,
                        'FAILED!') 
        self.failUnless(segs[0].tag=='tagin',
                        'FAILED!')  
        self.failUnless(segs[1].tag=='tagin',
                        'FAILED!') 
         
        self.failUnless(segs[2].tag=='bom',
                        'FAILED!') 
        self.failUnless(segs[3].tag=='',
                        'FAILED!')
        verts = m.getUserVertices()
        #print "User verts",verts
        #print 'polygon',polygon
        #vert values are relative
        for point,new_point in map(None,polygon,verts):
            point_x = point[0] 
            new_point_x = new_point.x + m.geo_reference.get_xllcorner()
            #print "point_x",point_x
            #print "new_point_x",new_point_x 
            point_y = point[1] 
            new_point_y = new_point.y + m.geo_reference.get_yllcorner()
            
            self.failUnless(point_x == new_point_x, ' failed')
            self.failUnless(point_y == new_point_y, ' failed')


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
        self.failUnless(is_inside_polygon([poly_point.x+x,poly_point.y+y],
                                       polygon_absolute, closed = False),
                        'FAILED!')
               
        self.failUnless(len(m.holes)==1,
                        'FAILED!')
        segs = m.getUserSegments()
        self.failUnless(len(segs)==4,
                        'FAILED!')
        self.failUnless(len(m.userVertices)==4,
                        'FAILED!') 
        self.failUnless(segs[0].tag=='tagin',
                        'FAILED!')  
        self.failUnless(segs[1].tag=='tagin',
                        'FAILED!') 
         
        self.failUnless(segs[2].tag=='bom',
                        'FAILED!') 
        self.failUnless(segs[3].tag=='',
                        'FAILED!')
        verts = m.getUserVertices()
        #print "User verts",verts
        #print 'polygon',polygon
        #vert values are relative
        for point,new_point in map(None,polygon,verts):
            point_x = point[0] + geo_ref_poly.get_xllcorner()
            new_point_x = new_point.x + m.geo_reference.get_xllcorner()
            point_y = point[1] + geo_ref_poly.get_yllcorner()
            #print "new_point.y",new_point.y
            #print "m.geo_ref.get_yllcorner()",m.geo_reference.get_yllcorner() 
            new_point_y = new_point.y + m.geo_reference.get_yllcorner()
            #print "point_y",point_y
            #print "new_point_y",new_point_y 
            
            self.failUnless(point_x == new_point_x, ' failed')
            self.failUnless(point_y == new_point_y, ' failed')

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
        
        self.failUnless(len(m.regions)==1,
                        'FAILED!')
        segs = m.getUserSegments()
        self.failUnless(len(segs)==segment_count,
                        'FAILED!')
        self.failUnless(len(m.userVertices)==segment_count,
                        'FAILED!') 
        self.failUnless(segs[0].tag==tag,
                        'FAILED!')  
        self.failUnless(segs[1].tag==tag,
                        'FAILED!') 
         
        verts = m.getUserVertices()
        
        #m.export_mesh_file("aaat.msh")
        
    def FIXMEtest_auto_set_geo_reference(self):
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
        self.failUnless(4==len(seg),
                        'FAILED!') 

        vert = m.getMeshVertices() 
        self.failUnless(4==len(vert),
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
        self.failUnless(4==len(seg),
                        'FAILED!') 

        vert = m.getMeshVertices() 
        self.failUnless(4==len(vert),
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
        self.failUnless(4==len(vert),
                        'FAILED!')
        vert= m.get_user_vertices(absolute=True)
        self.failUnless(vert==points_ab,
                        'FAILED!')

    
    def test_add_vertices_more(self):
        points = [[0.1,1],[0.4,.2],[7,5],[10,5]]
        #spat = Geospatial_data(points)
        
        m = Mesh()
        m.add_vertices(points)

        vert = m.getUserVertices()
        #print "vert",vert 
        self.failUnless(4==len(vert),
                        'FAILED!')
        vert= m.get_user_vertices(absolute=True)
        self.failUnless(vert==points,
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
        self.failUnless(4==len(vert),
                        'FAILED!')
        vert= m.get_user_vertices(absolute=True)
        self.failUnless(vert==points_ab,
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
        self.failUnless(4==len(vert),
                        'FAILED!')
        vert= m.get_user_vertices(absolute=True)
        self.failUnless(vert==points_ab,
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
        self.failUnless(4==len(vert),
                        'FAILED!')
        vert= m.get_user_vertices(absolute=True)
        self.failUnless(vert==points_ab,
                        'FAILED!')
        vert= m.get_user_vertices(absolute=False)
        points_new = m.geo_reference.get_absolute(vert)
        
        self.failUnless(points_ab==points_new,
                        'FAILED!')
        
        
def list_comp(A,B):
    yes = len(A)==len(B)
    for item in A:
        if not item in B:
            yes = False
    return yes


            
#-------------------------------------------------------------
if __name__ == "__main__":
    suite = unittest.makeSuite(meshTestCase,'test')
    #suite = unittest.makeSuite(meshTestCase,'test_exportoutlinefile')
    #suite = unittest.makeSuite(meshTestCase,'test_exportPointsFile')
    runner = unittest.TextTestRunner() # verbosity=2)
    runner.run(suite)
    

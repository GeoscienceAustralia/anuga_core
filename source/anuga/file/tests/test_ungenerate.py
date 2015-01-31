#!/usr/bin/env python
#

import os
import unittest
import tempfile

from anuga.pmesh.mesh import Vertex, Segment, Mesh
from anuga.file.ungenerate import load_ungenerate

class ungenerateTestCase(unittest.TestCase):
    def setUp(self):
        pass
    
    def tearDown(self):
        for filename in ['swamp.tsh']:
            try:
                os.remove(filename)
            except:
                pass
        
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
        dict = load_ungenerate(fileName)
        #os.remove(fileName)

        tag = "DSG"
        Segment.set_default_tag(tag)
        m.addVertsSegs(dict)

        # have to reset this , since it's a class attribute
        Segment.set_default_tag("")
            
        self.assertTrue(len(m.userSegments) ==11,
                        'Wrong segment list length.')
        self.assertTrue(len(m.userVertices) == 11,
                        'Wrong vertex list length.')
        self.assertTrue(m.userSegments[10].vertices[0] == m.userVertices[10],
                        'bad vertex on segment.')
        self.assertTrue(m.userSegments[10].vertices[1] == m.userVertices[8],
                        'Bad segment.')
        self.assertTrue(m.userSegments[10].tag == tag,
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

        self.assertTrue(Segment.get_default_tag() == initial_tag,
                        'Wrong segment list length.')
        

        # have to reset this , since it's a class attribute
        Segment.set_default_tag("")
            
        self.assertTrue(len(m.userSegments) ==11,
                        'Wrong segment list length.')
        self.assertTrue(len(m.userVertices) == 11,
                        'Wrong vertex list length.')
        self.assertTrue(m.userSegments[10].vertices[0] == m.userVertices[10],
                        'bad vertex on segment.')
        self.assertTrue(m.userSegments[10].vertices[1] == m.userVertices[8],
                        'Bad segment.')
        self.assertTrue(m.userSegments[10].tag == tag,
                        'wrong tag.')
        
    def test_import_ungenerate_file(self):
        
        fileName = tempfile.mktemp(".txt")
        file = open(fileName,"w")
        file.write("         1       ??      ??\n\
       10.0       10.0\n\
       11.0       10.0\n\
       11.0       11.0\n\
       10.0       11.0\n\
       10.0       10.0\n\
END\n\
         2      ?? ??\n\
       20.0       20.0\n\
       20.0       30.0\n\
       30.0       30.0\n\
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
        dict = load_ungenerate(fileName)
        #os.remove(fileName)

        tag = "DSG"
        Segment.set_default_tag(tag)
        m.addVertsSegs(dict)

        self.assertTrue(len(m.userSegments) ==11,
                        'Wrong segment list length.')
        self.assertTrue(len(m.userVertices) == 11,
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
        m.import_ungenerate_file(fileName, tag=tag, region_tag="swamp")

        os.remove(fileName)

        self.assertTrue(Segment.get_default_tag() == initial_tag,
                        'Wrong segment list length.')
        m.export_mesh_file("swamp.tsh")
        #print "m.userSegments",m.userSegments
        self.assertTrue(len(m.userSegments) ==11,
                        'Wrong segment list length.')
        self.assertTrue(len(m.userVertices) == 11,
                        'Wrong vertex list length.')
        self.assertTrue(len(m.regions) == 2,
                        'Wrong regions list length.')
        self.assertTrue(m.regions[0].getTag() == "swamp",
                        'Wrong regions tag.')
        self.assertTrue(m.regions[1].getTag() == "swamp",
                        'Wrong regions 1 tag.')
        
        # have to reset this , since it's a class attribute
        Segment.set_default_tag("")
        
        
    def test_import_ungenerate_file_different_region_tags(self):
        
        fileName = tempfile.mktemp(".txt")
        file = open(fileName,"w")
        file.write("         1       ??      ??\n\
       10.0       10.0\n\
       11.0       10.0\n\
       11.0       11.0\n\
       10.0       11.0\n\
       10.0       10.0\n\
END\n\
         2      ?? ??\n\
       20.0       20.0\n\
       20.0       30.0\n\
       30.0       30.0\n\
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
        dict = load_ungenerate(fileName)
        #os.remove(fileName)

        tag = "DSG"
        Segment.set_default_tag(tag)
        m.addVertsSegs(dict)

        self.assertTrue(len(m.userSegments) ==11,
                        'Wrong segment list length.')
        self.assertTrue(len(m.userVertices) == 11,
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
        m.import_ungenerate_file(fileName, tag=tag, region_tag=["swamp","coastalp"])

        os.remove(fileName)

        self.assertTrue(Segment.get_default_tag() == initial_tag,
                        'Wrong segment list length.')
        m.export_mesh_file("swamp.tsh")
        #print "m.userSegments",m.userSegments
        self.assertTrue(len(m.userSegments) ==11,
                        'Wrong segment list length.')
        self.assertTrue(len(m.userVertices) == 11,
                        'Wrong vertex list length.')
        self.assertTrue(len(m.regions) == 2,
                        'Wrong regions list length.')
        self.assertTrue(m.regions[0].getTag() == "swamp",
                        'Wrong regions tag.')
        self.assertTrue(m.regions[1].getTag() == "coastalp",
                        'Wrong regions 1 tag.')
        
        
        # have to reset this , since it's a class attribute
        Segment.set_default_tag("")



################################################################################

if __name__ == "__main__":
    suite = unittest.makeSuite(ungenerateTestCase,'test')
    runner = unittest.TextTestRunner() #verbosity=2)
    runner.run(suite)
    

#!/usr/bin/env python
#
"""
I removed lone test vert's, since I'm working on removing lone verts at a lower
level of the code, using the -j flag in triangle.
"""


import sys

import unittest
from anuga.mesh_engine.mesh_engine import generate_mesh

import numpy as num

from anuga.utilities.numerical_tools import ensure_numeric

from anuga.anuga_exceptions import ANUGAError

class triangTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testrectangleIIb(self):
        
        points = []
        seglist = []
        holelist = []
        regionlist = []

        points = [(0.0,0.0),(0.0,10.0),(3.0,0.0),(3.0,10.0)]
        pointattlist = []
        regionlist.append( (1.2,1.2,5.0) )
        seglist = [(0,1),(1,3),(3,2),(2,0)]
        segattlist = []
        
        mode = "Qzp"
        
        data = generate_mesh(points,seglist,holelist,regionlist,
                              pointattlist,segattlist, mode, points)

        correct = num.array([(1, 0, 2), (2, 3, 1)])
        self.assertTrue(num.alltrue(data['generatedtrianglelist'].flat == \
                                    correct.flat),
                        'trianglelist is wrong!')
        correct = num.array([(0, 1), (1, 3), (3, 2), (2, 0)])
        self.assertTrue(num.alltrue(data['generatedsegmentlist'].flat == \
                                    correct.flat),
                        'segmentlist is wrong!')

        correct = num.array([(0.0, 0.0), (0.0, 10.0),
                             (3.0, 0.0), (3.0, 10.0)])
        self.assertTrue(num.allclose(data['generatedpointlist'].flat, \
                                     correct.flat),
                        'Failed')
     
    def testrectangleIIc(self):
        points = []
        seglist = []
        holelist = []
        regionlist = []

        points = [(0.0,0.0),(0.0,10.0),(3.0,0.0),(3.0,10.0)]
        pointattlist = None
        regionlist.append( (1.2,1.2,5.0) )
        seglist = [(0,1),(1,3),(3,2),(2,0)]
        segattlist = None
        
        mode = "Qzp"
        
        data = generate_mesh(points,seglist,holelist,regionlist,
                              pointattlist,segattlist, mode, points)
        correct = num.array([(1, 0, 2), (2, 3, 1)])
        self.assertTrue(num.alltrue(data['generatedtrianglelist'].flat == \
                                    correct.flat),
                        'trianglelist is wrong!')
        correct = num.array([(0, 1), (1, 3), (3, 2), (2, 0)])
        self.assertTrue(num.alltrue(data['generatedsegmentlist'].flat == \
                                    correct.flat),
                        'segmentlist is wrong!')

        correct = num.array([(0.0, 0.0), (0.0, 10.0),
                             (3.0, 0.0), (3.0, 10.0)])
        self.assertTrue(num.allclose(data['generatedpointlist'].flat, \
                                     correct.flat),
                        'Failed')

    def test_bad_point_atts(self):

        points = []
        seglist = []
        holelist = []
        regionlist = []

        points = [(0.0,0.0),(0.0,10.0),(3.0,0.0),(3.0,10.0)]
        pointattlist = [[],[],[]]
        regionlist.append( (1.2,1.2,5.0) )
        seglist = [(0,1),(1,3),(3,2),(2,0)]
        segattlist = [0,0,0,0]
        
        mode = "Qzp"
        
        try:
            data = generate_mesh(points,seglist,holelist,
                                  regionlist,pointattlist,segattlist,
                                   mode, points)
        except ANUGAError:
            pass
        else:
            self.assertTrue(0 ==1,
                        'bad list did not raise error!')
            
        
    def test_bad_point_attsII(self):
        points = []
        seglist = []
        holelist = []
        regionlist = []

        points = [(0.0,0.0),(0.0,10.0),(3.0,0.0),(3.0,10.0)]
        pointattlist = [[1],[2],[3],[4,8]]
        regionlist.append( (1.2,1.2,5.0) )
        seglist = [(0,1),(1,3),(3,2),(2,0)]
        segattlist = [0,0,0,0]
        
        mode = "Qzp"
        
        try:
            data = generate_mesh(points,seglist,holelist,
                                  regionlist,pointattlist,segattlist,
                                   mode, points)
        except ANUGAError:
            pass
        else:
            self.assertTrue(0 ==1,
                        'bad list did not raise error!')
            
    def testsegmarker(self):

        holelist = []
        regionlist = []

        points = [(0.0,0.0),(0.0,10.0),(3.0,0.0),(3.0,10.0)]
        pointattlist = [[],[],[],[]]
        regionlist.append( (1.2,1.2,5.0) )
        seglist = [(0,1),(1,3),(3,2),(2,0)]
        segattlist = [1.0,2.0,3.0,4.0]
        
        mode = "Qzp"
        data = generate_mesh(points,seglist,holelist,regionlist,
                              pointattlist,segattlist, mode, points)

        correct = num.array([(1, 0, 2), (2, 3, 1)])
        self.assertTrue(num.alltrue(data['generatedtrianglelist'].flat == \
                                    correct.flat),
                        'trianglelist is wrong!')
        correct = num.array([(0, 1), (1, 3), (3, 2), (2, 0)])
        self.assertTrue(num.alltrue(data['generatedsegmentlist'].flat == \
                                    correct.flat),
                        'segmentlist is wrong!')

        correct = num.array([(0.0, 0.0), (0.0, 10.0),
                             (3.0, 0.0), (3.0, 10.0)])
        self.assertTrue(num.allclose(data['generatedpointlist'].flat, \
                                     correct.flat),
                        'Failed')
        
        self.assertTrue(num.alltrue(data['generatedsegmentmarkerlist'] == \
                                    num.array([1,2,3,4])),
                        'Failed!')
        
    def testbad_region(self):

        points = []
        seglist = []
        holelist = []
        regionlist = [(1.2,1.2)]

        points = [(0.0,0.0),(0.0,10.0),(3.0,0.0),(3.0,10.0)]
        pointattlist = [[],[],[],[]]
        seglist = [(0,1),(1,3),(3,2),(2,0)]
        segattlist = [0,0,0,0]
        
        mode = "Qzpn"
        try:
            data = generate_mesh(points,seglist,holelist,
                                  regionlist,pointattlist,segattlist,
                                   mode, points)
        except ANUGAError:
            pass
        else:
            self.assertTrue(0 ==1,
                        'bad region list did not raise error!')

    def testbad_regionII(self):

        points = []
        seglist = []
        holelist = []
        regionlist = [(1.2,1.2), (1.2,1.25,1.0)]

        points = [(0.0,0.0),(0.0,10.0),(3.0,0.0),(3.0,10.0)]
        pointattlist = [[],[],[],[]]
        seglist = [(0,1),(1,3),(3,2),(2,0)]
        segattlist = [0,0,0,0]
        
        mode = "Qzpn"
        try:
            data = generate_mesh(points,seglist,holelist,
                                  regionlist,pointattlist,segattlist,
                                   mode, points)
        except ANUGAError:
            pass
        else:
            self.assertTrue(0==1,
                        'bad region list did not raise error!')

    def testregion_with_maxarea(self):

        points = []
        seglist = []
        holelist = []
        regionlist = [(3,1,1.0)]

        points = [(0.0,0.0),(6.0,0.0),(6.0,6.0),(0.0,6.0)]
        pointattlist = [[],[],[],[]]
        seglist = [(0,1),(1,2),(3,2),(3,0),(0,2)]
        segattlist = [0,0,0,0,0]
        
        mode = "Qzpna36a"
        data = generate_mesh(points,seglist,holelist,regionlist,
                              pointattlist,segattlist,  mode, points)

        self.assertTrue(len(data['generatedtrianglelist']) == 2,
                        'testregion_with_maxarea 1: # of tris is wrong!')
        ## Another test case
        regionlist = [(3,1,1.0),(1,3,1.0,8.0)]
        mode = "Qzp21na36a"
        data = generate_mesh(points,seglist,holelist,regionlist,
                              pointattlist,segattlist,  mode, points)
        #print "len(data['generatedtrianglelist']",len(data['generatedtrianglelist'])
        # for Duncan On unix this returns a 7 triangle result.
        # for Duncan on Windows returns a 6 triangle result.
        # for Ole on nautilus this returns 6
        # for Duncan on nautilus this returns 7
        # It seems to be the results from triangle that is
        # causing the different results, and we are treating
        # triangle as a back box.

        self.assertTrue(len(data['generatedtrianglelist']) >= 6,
                        'testregion_with_maxarea 2: # of tris is wrong!')
        ## Another test case
        regionlist = [(3,1,1.0,8.0),(1,3,1.0,8.0)]
        mode = "Qzpna36a"
        data = generate_mesh(points,seglist,holelist,regionlist,
                              pointattlist,segattlist,  mode, points)
        #print "len(data['generatedtrianglelist']",len(data['generatedtrianglelist'])
        # On unix this returns a 10 triangle result.
        # Windows returns a 8 triangle result.
        self.assertTrue(len(data['generatedtrianglelist']) >= 8,
                        'testregion_with_maxarea 3: # of tris is wrong!')

        ## Another test case
        regionlist = [(3,1,1.0),(1,3,1.0,8.0)]
        mode = "Qzpna8a"
        data = generate_mesh(points,seglist,holelist,regionlist,
                              pointattlist,segattlist,  mode, points)
        #print "len(data['generatedtrianglelist']",len(data['generatedtrianglelist'])
        # On unix this returns a 10 triangle result.
        # Windows returns a 8 triangle result.

        self.assertTrue(len(data['generatedtrianglelist']) >= 8,
                        'testregion_with_maxarea 4: # of tris is wrong!')


    def FIXME_bad_hole(self):

        holelist = [(9.0)]
        regionlist = []

        points = [(0.0,0.0),(0.0,10.0),(3.0,0.0),(3.0,10.0)]
        pointattlist = [[],[],[],[]]
        seglist = [(0,1),(1,3),(3,2),(2,0)]
        segattlist = [0,0,0,0]
        
        mode = "Qzpn"
        try:
            data = generate_mesh(points,seglist,holelist,regionlist,
                                  pointattlist,segattlist,  mode, points)

        except TypeError:
            pass
        else:
            self.assertTrue(0 ==1,
                        'bad hole list did not raise error!')


    def testbad_segattlist(self):
        holelist = []
        regionlist = []

        points = [(0.0,0.0),(0.0,10.0),(3.0,0.0),(3.0,10.0)]
        pointattlist = [[],[],[],[]]
        seglist = [(0,1),(1,3),(3,2),(2,0)]
        segattlist = [0,0]
        
        mode = "Qzpn"
        try:
            data = generate_mesh(points,seglist,holelist,regionlist,
                                  pointattlist,segattlist, mode, points)
            
        except ANUGAError:
            pass
        else:
            self.assertTrue(0 ==1,
                        'bad segment attribute list did not raise error!')

    def testrectangle_regions(self):

        points = []
        seglist = []
        holelist = []
        regionlist = []

        points = [(0.0,0.0),(0.0,10.0),(3.0,0.0),(3.0,10.0)]
        pointattlist = [[],[],[],[]]
        #it seems that
        #triangle only associates one region with a triangle
        regionlist.append( (1.3,1.3,88.33) )
        regionlist.append( (1.2,1.2,77,55) )
        seglist = [(0,1),(1,3),(3,2),(2,0)]
        segattlist = [0,0,0,0]
        mode = "QAzpq"
        data = generate_mesh(points,seglist,holelist,regionlist,
                              pointattlist,segattlist, mode, points)

        correct = num.array([[77.0], [77.0], [77.0], [77.0]])
        self.assertTrue(num.allclose(data['generatedtriangleattributelist'].flat, 
                                     correct.flat),
                        'Failed')
        
    def testrectangle_regionsII(self):
        points = []
        seglist = []
        holelist = []
        regionlist = []

        points = [(0.0,0.0),(0.0,10.0),(3.0,0.0),(3.0,10.0)]
        pointattlist = [[],[],[],[]]
        #it seems that
        #triangle only associates one region with a triangle
        regionlist.append( [1.3,1.3,88.33] )
        regionlist.append( [1.2,1.2,77,55] )
        seglist = [(0,1),(1,3),(3,2),(2,0)]
        segattlist = [0,0,0,0]
        
        mode = "QAzpq"
        data = generate_mesh(points,seglist,holelist,regionlist,
                              pointattlist,segattlist, mode, points)

        correct = num.array([[77.0], [77.0], [77.0], [77.0]])
        self.assertTrue(num.allclose(data['generatedtriangleattributelist'].flat, 
                                     correct.flat),
                        'Failed')


    def test_numeric_arrays(self):
        points = []
        seglist = []
        holelist = []
        regionlist = []
        points = [(0.0,0.0),(0.0,10.0),(3.0,0.0),(3.0,10.0)]
        pointattlist = []
        # 5.0 is the region tag, 99.0 is the max area
        tri_tag = 123456.0
        regionlist.append( [0.2,0.2, tri_tag,99.0] )
        seglist = [(0,1),(1,3),(3,2),(2,0)]
        segattlist = [21,22,23,24]
         #The 'A' has to be there to get the region marker stuff working
        mode = "QzpnA"
        #mode = "jQpznAa2000.1a"
        data = generate_mesh(points,seglist,holelist,regionlist,
                              pointattlist,segattlist, mode, points)
        #print "data", data

     
        correct = num.array([(1, 0, 2), (2, 3, 1)])
        self.assertTrue(num.alltrue(data['generatedtrianglelist'].flat == \
                                    correct.flat),
                        'trianglelist is wrong!')
        correct = num.array([(0, 1), (1, 3), (3, 2), (2, 0)])
        self.assertTrue(num.alltrue(data['generatedsegmentlist'].flat == \
                                    correct.flat),
                        'segmentlist is wrong!')

        correct = num.array([(0.0, 0.0), (0.0, 10.0),
                             (3.0, 0.0), (3.0, 10.0)])
        self.assertTrue(num.allclose(data['generatedpointlist'].flat, \
                                     correct.flat),
                        'Failed')
        
        correct = num.array([[tri_tag], [tri_tag]])
        self.assertTrue(num.allclose(data['generatedtriangleattributelist'].flat, \
                                     correct.flat),
                        'Failed')
        correct = num.array([(0, 1), (1, 3), (3, 2), (2, 0)])
        self.assertTrue(num.alltrue(data['generatedsegmentlist'].flat == \
                                    correct.flat),
                        'Failed!')
        
        correct = num.array(segattlist, num.int)
        self.assertTrue(num.allclose(data['generatedsegmentmarkerlist'].flat, 
                                     correct.flat),
                        'Failed')
        
        # I copied these answers from the output, so bad test..
        correct = num.array([(-1, 1, -1), (-1, 0, -1)])
        self.assertTrue(num.alltrue(data['generatedtriangleneighborlist'].flat == \
                                    correct.flat),
                        'Failed!')
        
        
    def test_pointattlist(self):
        # segattlist = []
        points = []
        seglist = []
        holelist = []
        regionlist = []

        points = [(0.0,0.0),(0.0,4.0),(4.0,2.0),(2.0,0.0)]
        pointattlist = [0.,0.,10.,10.]
        regionlist.append( [0.2,0.2,2.1, 99.] )
        seglist = [(0,1),(1,2),(2,3),(3,0)]
        segattlist = [11,12,13,14]
        mode = "Qzp"
        data = generate_mesh(points,seglist,holelist,regionlist,
                              pointattlist,segattlist, mode, points)

        correct = num.array([[0.0],[0.0],[10],[10]])
        self.assertTrue(num.allclose(data['generatedpointattributelist'].flat, \
                                     correct.flat),
                        'Failed')
        
        pointattlist = [[0.],[0.],[10.],[10.]]
        mode = "Qzp"        
        data = generate_mesh(points,seglist,holelist,regionlist,
                              pointattlist,segattlist, mode, points)
        correct = num.array([[0.0],[0.0],[10],[10]])
        self.assertTrue(num.allclose(data['generatedpointattributelist'].flat, \
                                     correct.flat),
                        'Failed')
        
        pointattlist = [[0.,1],[0.,1],[10.,20],[10.,20]]
        mode = "Qzp"        
        data = generate_mesh(points,seglist,holelist,regionlist,
                              pointattlist,segattlist, mode, points)
        #print "data", data
        correct = num.array(pointattlist)
        self.assertTrue(num.allclose(data['generatedpointattributelist'].flat, \
                                     correct.flat),
                        'Failed')
    

if __name__ == "__main__":
    suite = unittest.makeSuite(triangTestCase,'test')
    runner = unittest.TextTestRunner()  #verbosity=2)
    runner.run(suite)

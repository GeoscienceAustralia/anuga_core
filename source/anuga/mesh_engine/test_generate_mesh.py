#!/usr/bin/env python
#

import sys

import unittest
from anuga.mesh_engine.mesh_engine import generate_mesh

from Numeric import array, Float, Int

from anuga.utilities.numerical_tools import ensure_numeric

from anuga.utilities.anuga_exceptions import ANUGAError

class triangTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def metestrectangleII(self):

        points = []
        seglist = []
        holelist = []
        regionlist = []

	points = [(0.0,0.0),(0.0,10.0),(3.0,0.0),(3.0,10.0)]
	pointattlist = [[],[],[],[]]
	regionlist.append( (1.2,1.2,5.0) )
	seglist = [(0,1),(1,3),(3,2),(2,0)]
        segattlist = [0,0,0,0]
        
        mode = "Qzp"
        
        points =  ensure_numeric(points, Float)
        seglist = test = ensure_numeric(seglist, Int)
        data = triang.genMesh(points,seglist,holelist,regionlist,
                              pointattlist,segattlist, mode, test)

        self.failUnless(data['generatedtrianglelist'] ==[(1, 0, 2), (2, 3, 1)],
                        'trianglelist is wrong!')
        self.failUnless(data['generatedsegmentlist'] ==[(0, 1), (1, 3),
                                                        (3, 2), (2, 0)],
                        'segmentlist is wrong!')
        self.failUnless(data['generatedpointlist'] ==[(0.0, 0.0), (0.0, 10.0),
                                                      (3.0, 0.0), (3.0, 10.0)],
                        ' is wrong!')


    def testrectangleIIb(self):
        # segattlist = []
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

        self.failUnless(data['generatedtrianglelist'] ==[(1, 0, 2), (2, 3, 1)],
                        'trianglelist is wrong!')
        self.failUnless(data['generatedsegmentlist'] ==[(0, 1), (1, 3),
                                                        (3, 2), (2, 0)],
                        'segmentlist is wrong!')
        self.failUnless(data['generatedpointlist'] ==[(0.0, 0.0), (0.0, 10.0),
                                                      (3.0, 0.0), (3.0, 10.0)],
                        ' is wrong!')
     
    def testrectangleIIc(self):
        # segattlist = None
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

        self.failUnless(data['generatedtrianglelist'] ==[(1, 0, 2), (2, 3, 1)],
                        'trianglelist is wrong!')
        self.failUnless(data['generatedsegmentlist'] ==[(0, 1), (1, 3),
                                                        (3, 2), (2, 0)],
                        'segmentlist is wrong!')
        self.failUnless(data['generatedpointlist'] ==[(0.0, 0.0), (0.0, 10.0),
                                                      (3.0, 0.0), (3.0, 10.0)],
                        ' is wrong!')
        
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
            self.failUnless(0 ==1,
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
            self.failUnless(0 ==1,
                        'bad list did not raise error!')
            
    def testsegmarker(self):

        points = []
        seglist = []
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

        self.failUnless(data['generatedtrianglelist'] ==[(1, 0, 2), (2, 3, 1)],
                        'trianglelist is wrong!')
        self.failUnless(data['generatedsegmentlist'] ==[(0, 1), (1, 3),
                                                        (3, 2), (2, 0)],
                        'segmentlist is wrong!')
        self.failUnless(data['generatedpointlist'] ==[(0.0, 0.0), (0.0, 10.0),
                                                      (3.0, 0.0), (3.0, 10.0)],
                        'generatedpointlist is wrong!')
        self.failUnless(data['generatedsegmentmarkerlist'] ==[1,2,3,4],
                        'generatedsegmentmarkerlist is wrong!')

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
            self.failUnless(0 ==1,
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
            self.failUnless(0==1,
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

        self.failUnless(len(data['generatedtrianglelist']) == 2,
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
        # ??, it seems to be the results from triangle that is
        # causing the different results, and we are treating
        # triangle as a back box.

        self.failUnless(len(data['generatedtrianglelist']) >= 6,
                        'testregion_with_maxarea 2: # of tris is wrong!')
        ## Another test case
        regionlist = [(3,1,1.0,8.0),(1,3,1.0,8.0)]
        mode = "Qzpna36a"
        data = generate_mesh(points,seglist,holelist,regionlist,
                              pointattlist,segattlist,  mode, points)
        #print "len(data['generatedtrianglelist']",len(data['generatedtrianglelist'])
        # On unix this returns a 10 triangle result.
        # Windows returns a 8 triangle result.
        self.failUnless(len(data['generatedtrianglelist']) >= 8,
                        'testregion_with_maxarea 3: # of tris is wrong!')

        ## Another test case
        regionlist = [(3,1,1.0),(1,3,1.0,8.0)]
        mode = "Qzpna8a"
        data = generate_mesh(points,seglist,holelist,regionlist,
                              pointattlist,segattlist,  mode, points)
        #print "len(data['generatedtrianglelist']",len(data['generatedtrianglelist'])
        # On unix this returns a 10 triangle result.
        # Windows returns a 8 triangle result.

        self.failUnless(len(data['generatedtrianglelist']) >= 8,
                        'testregion_with_maxarea 4: # of tris is wrong!')


    def FIXME_testbad_hole(self):

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
            self.failUnless(0 ==1,
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
            print "data",data
            self.failUnless(data['generatedtrianglelist'] ==[(1, 0, 2), (2, 3, 1)],
                        'trianglelist is wrong!')
        except ANUGAError:
            pass
        else:
            self.failUnless(0 ==1,
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

        self.failUnless(data['generatedtriangleattributelist'] ==[[77.0], [77.0], [77.0], [77.0]],
                        'triangleattributelist is wrong!')

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

        self.failUnless(data['generatedtriangleattributelist'] ==[[77.0], [77.0], [77.0], [77.0]],
                        'triangleattributelist is wrong!')


    def BADtest_lone_verts(self):
        print "test_lone_verts"
        points = []
        seglist = []
        holelist = []
        regionlist = []

	points = [(0.0,0.0),(0.0,10.0),(3.0,0.0),(3.0,10.0),(0.0,10.0)]
	pointattlist = [[],[],[],[],[]]
	regionlist.append( (1.2,1.2,5.0) )
	seglist = [(0,4),(4,1),(1,3),(3,2),(2,0)]
        segattlist = [0,0,0,0,0]
        
        mode = "QpznAa2000.1a"
        data = generate_mesh(points,seglist,holelist,regionlist,
                              pointattlist,segattlist, mode, points)
        #print "data['generatedtrianglelist']", data['generatedtrianglelist']
        self.failUnless(data['generatedtrianglelist'] ==[(4, 0, 2), (2, 3, 4)],
                        'trianglelist is wrong!')
        self.failUnless(data['generatedsegmentlist'] ==[(0, 4), (4, 3),
                                                        (3, 2), (2, 0)],
                        'segmentlist is wrong!')
        self.failUnless(data['generatedpointlist'] ==[(0.0, 0.0), (0.0, 10.0),
                                                      (3.0, 0.0), (3.0, 10.0),
                                                      (0.0,10.0)],
                        ' is wrong!')
        #print "data['lonepointlist']", data['lonepointlist']
        self.failUnless(data['lonepointlist'] ==[1],
                        'lonepointlist is wrong!')

    def test_lone_vertsII(self):

        points = []
        seglist = []
        holelist = []
        regionlist = []

	points = [(0.0,0.0),(0.0,0.0),(0.0,10.0),(0.0,10.0),(10.0,10.0),
                  (10.0,10.0),(0.0,10.0),(10.0,0.0)]

	pointattlist = []
        for point in points:
            pointattlist.append([])
	seglist = [(0,1),(1,2),(2,3),(3,4),(4,5),(5,7),(7,0)]
	segattlist = []
        for seg in seglist:
            segattlist.append(0)
        
        mode = "QpznAa2000.1a"
        data = generate_mesh(points,seglist,holelist,regionlist,
                              pointattlist,segattlist, mode, points)
        #print "data['generatedtrianglelist']", data['generatedtrianglelist']
        self.failUnless(data['generatedtrianglelist'] ==[(6, 1, 7), (7, 5, 6)],
                        'trianglelist is wrong!')
        self.failUnless(data['generatedsegmentlist'] ==[(1, 6), (6, 5),
                                                        (5, 7), (7, 1)],
                        'segmentlist is wrong!')
        self.failUnless(data['generatedpointlist'] ==[(0.0,0.0),(0.0,0.0),
                                                      (0.0,10.0),(0.0,10.0),
                                                      (10.0,10.0),
                                                      (10.0,10.0),(0.0,10.0),
                                                      (10.0,0.0)],
                        ' is wrong!')
        self.failUnless(data['lonepointlist'] ==[0,2,3,4],
                        'lonepointlist is wrong!')

    def test_lone_vertsIII(self):

        points = []
        seglist = []
        holelist = []
        regionlist = []

	points = [(-3.,-5.),(0.0,0.0),(0.0,0.0),
                  (-5.,-8.),(0.0,10.0),(0.0,10.0),
                  (-3.,-5.),(10.0,10.0), (10.0,10.0),
                  (0.0,10.0),
                  (10.0,0.0),
                  (-12.,45.)]

	pointattlist = []
        for point in points:
            pointattlist.append([])
	seglist = [(1,2),(2,4),(4,5),(5,7),(7,8),(8,10),(10,1)]
	segattlist = []
        for seg in seglist:
            segattlist.append(0)
            """
            0 1
            1 2
            2 4
            3 5
            4 7
            5 8
            6 9
            7 10
            """
        
        mode = "QpznAa2000.1a"
        data = generate_mesh(points,seglist,holelist,regionlist,
                              pointattlist,segattlist, mode, points)
        #print "data['generatedtrianglelist']", data['generatedtrianglelist']
        #self.failUnless(data['generatedtrianglelist'] ==[(9, 2, 10),
         #                                                (10, 8, 9)],
          #              'trianglelist is wrong!')
        self.failUnless(data['generatedtrianglelist'] ==[(2, 10, 8),
                                                         (8, 9, 2)],
                        'trianglelist is wrong!')
        #print "data['generatedsegmentlist']",data['generatedsegmentlist']
        self.failUnless(data['generatedsegmentlist'] ==[(9, 2), (9, 8),
                                                        (8, 10), (2, 10)],
                        'segmentlist is wrong!')
        self.failUnless(data['generatedpointlist'] ==points,
                        ' is wrong!')
        self.failUnless(data['lonepointlist'] ==[0,1,3,4,5,6,7,11],
                        'lonepointlist is wrong!')
        
    def test_lone_verts4(self):
        points = []
        seglist = []
        holelist = []
        regionlist = []

	points = [(-56.,-78),(0.0,0.0),(0.0,10.0),(10.,0.)]
	pointattlist = []
	seglist = [(1,2),(2,3),(3,1)]
        segattlist = []
        
        mode = "QpznAa2000.1a"
        data = generate_mesh(points,seglist,holelist,regionlist,
                              pointattlist,segattlist, mode, points)
        #print "data['generatedtrianglelist']", data['generatedtrianglelist']
        #self.failUnless(data['generatedtrianglelist'] ==[(4, 0, 2),(2, 3, 4)],
         #               'trianglelist is wrong!')
        self.failUnless(data['generatedsegmentlist'] ==seglist,
                        'segmentlist is wrong!')
        self.failUnless(data['generatedpointlist'] ==points,
                        ' is wrong!')
        #print "data['lonepointlist']", data['lonepointlist']
        self.failUnless(data['lonepointlist'] ==[0],
                        'lonepointlist is wrong!')
if __name__ == "__main__":

    suite = unittest.makeSuite(triangTestCase,'test')
    #suite = unittest.makeSuite(triangTestCase,'test_lone_verts4')
    #suite = unittest.makeSuite(triangTestCase,'testrectangleIIb')
    runner = unittest.TextTestRunner(verbosity=2) #verbosity=2)
    runner.run(suite)

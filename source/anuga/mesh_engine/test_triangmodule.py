#!/usr/bin/env python
#

import sys

import unittest
import anuga.mesh_engine.triang as triang


class triangTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testrectangle(self):

        points = []
        seglist = []
        segattlist = []
        holelist = []
        regionlist = []

	points = [(0.0,0.0),(0.0,10.0),(3.0,0.0),(3.0,10.0)]
	# This will cause a seg fault
        #points = [[0.0,0.0],[0.0,10.0],[3.0,0.0],[3.0,10.0]]

        # automatically generate pointattlist
	pointattlist = []
        for point in points:
            pointattlist.append([])
        #print "pointattlist",pointattlist
	#pointattlist = [[],[],[],[]]
        trilist = []
        mode = "Qzcn"
        data = triang.genMesh(points,seglist,holelist,regionlist,
                              pointattlist,segattlist,trilist, mode)
        self.failUnless(data['generatedtrianglelist'] ==[(1, 0, 2), (2, 3, 1)],
                        'generatedtrianglelist is wrong!')
        #print "data['generatedsegmentlist']",data['generatedsegmentlist']
        self.failUnless(data['generatedsegmentlist'] ==[(2, 0), (3, 2), (1, 3),
                                                        (0, 1)]
              ,          'generatedsegmentlist is wrong!')
        self.failUnless(data['generatedpointlist'] ==[(0.0, 0.0), (0.0, 10.0),
                                                      (3.0, 0.0), (3.0, 10.0)],
                        ' is wrong!')

        #print "data['generatedtriangleneighborlist']",
        #data['generatedtriangleneighborlist']
        self.failUnless(data['generatedtriangleneighborlist'] ==[(-1, 1, -1),
                                                                 (-1, 0, -1)],
                        'generatedtriangleneighborlist is wrong!')


    def testrectangleII(self):

        points = []
        seglist = []
        holelist = []
        regionlist = []

	points = [(0.0,0.0),(0.0,10.0),(3.0,0.0),(3.0,10.0)]
	pointattlist = [[],[],[],[]]
	regionlist.append( (1.2,1.2,5.0) )
	seglist = [(0,1),(1,3),(3,2),(2,0)]
        segattlist = [0,0,0,0]
        trilist = []
        mode = "Qzp"
        data = triang.genMesh(points,seglist,holelist,regionlist,
                              pointattlist,segattlist,trilist, mode)

        self.failUnless(data['generatedtrianglelist'] ==[(1, 0, 2), (2, 3, 1)],
                        'trianglelist is wrong!')
        self.failUnless(data['generatedsegmentlist'] ==[(0, 1), (1, 3),
                                                        (3, 2), (2, 0)],
                        'segmentlist is wrong!')
        self.failUnless(data['generatedpointlist'] ==[(0.0, 0.0), (0.0, 10.0),
                                                      (3.0, 0.0), (3.0, 10.0)],
                        ' is wrong!')

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
        trilist = []
        mode = "Qzp"
        data = triang.genMesh(points,seglist,holelist,regionlist,
                              pointattlist,segattlist,trilist, mode)

        self.failUnless(data['generatedtrianglelist'] ==[(1, 0, 2), (2, 3, 1)],
                        'trianglelist is wrong!')
        self.failUnless(data['generatedsegmentlist'] ==[(0, 1), (1, 3),
                                                        (3, 2), (2, 0)],
                        'segmentlist is wrong!')
        self.failUnless(data['generatedpointlist'] ==[(0.0, 0.0), (0.0, 10.0),
                                                      (3.0, 0.0), (3.0, 10.0)],
                        ' is wrong!')
        self.failUnless(data['generatedsegmentmarkerlist'] ==[1,2,3,4],
                        ' is wrong!')

    def testbad_region(self):

        points = []
        seglist = []
        holelist = []
        regionlist = [(1.2,1.2)]

	points = [(0.0,0.0),(0.0,10.0),(3.0,0.0),(3.0,10.0)]
	pointattlist = [[],[],[],[]]
	seglist = [(0,1),(1,3),(3,2),(2,0)]
        segattlist = [0,0,0,0]
        trilist = []
        mode = "Qzpn"
        try:
            data = triang.genMesh(points,seglist,holelist,
                                  regionlist,pointattlist,segattlist,
                                  trilist, mode)

        except TypeError:
            pass
        else:
            self.failUnless(0 ==1,
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
        trilist = []
        mode = "Qzpna36a"
        data = triang.genMesh(points,seglist,holelist,regionlist,
                              pointattlist,segattlist, trilist, mode)

        self.failUnless(len(data['generatedtrianglelist']) == 2,
                        'testregion_with_maxarea 1: # of tris is wrong!')
        ## Another test case
        regionlist = [(3,1,1.0),(1,3,1.0,8.0)]
        mode = "Qzp21na36a"
        data = triang.genMesh(points,seglist,holelist,regionlist,
                              pointattlist,segattlist, trilist, mode)
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
        data = triang.genMesh(points,seglist,holelist,regionlist,
                              pointattlist,segattlist, trilist, mode)
        #print "len(data['generatedtrianglelist']",len(data['generatedtrianglelist'])
        # On unix this returns a 10 triangle result.
        # Windows returns a 8 triangle result.
        self.failUnless(len(data['generatedtrianglelist']) >= 8,
                        'testregion_with_maxarea 3: # of tris is wrong!')

        ## Another test case
        regionlist = [(3,1,1.0),(1,3,1.0,8.0)]
        mode = "Qzpna8a"
        data = triang.genMesh(points,seglist,holelist,regionlist,
                              pointattlist,segattlist, trilist, mode)
        #print "len(data['generatedtrianglelist']",len(data['generatedtrianglelist'])
        # On unix this returns a 10 triangle result.
        # Windows returns a 8 triangle result.

        self.failUnless(len(data['generatedtrianglelist']) >= 8,
                        'testregion_with_maxarea 4: # of tris is wrong!')

    def testbad_point(self):

        holelist = []
        regionlist = []

	points = [(0.0,0.0),(0.0,10.0),(3.0,0.0),(3.0)]
	pointattlist = [[],[],[],[]]
	seglist = [(0,1),(1,3),(3,2),(2,0)]
        segattlist = [0,0,0,0]
        trilist = []
        mode = "Qzpn"
        try:
            data = triang.genMesh(points,seglist,holelist,regionlist,
                                  pointattlist,segattlist, trilist, mode)

        except TypeError:
            pass
        else:
            self.failUnless(0 ==1,
                        'bad point list did not raise error!')

    def testbad_hole(self):

        holelist = [(9.0)]
        regionlist = []

	points = [(0.0,0.0),(0.0,10.0),(3.0,0.0),(3.0,10.0)]
	pointattlist = [[],[],[],[]]
	seglist = [(0,1),(1,3),(3,2),(2,0)]
        segattlist = [0,0,0,0]
        trilist = []
        mode = "Qzpn"
        try:
            data = triang.genMesh(points,seglist,holelist,regionlist,
                                  pointattlist,segattlist, trilist, mode)

        except TypeError:
            pass
        else:
            self.failUnless(0 ==1,
                        'bad hole list did not raise error!')

    def testbad_segment(self):

        holelist = []
        regionlist = []

	points = [(0.0,0.0),(0.0,10.0),(3.0,0.0),(3.0,10.0)]
	pointattlist = [[],[],[],[]]
	seglist = [(0,1),(1,3),(3,2),(2)]
        segattlist = [0,0,0,0]
        trilist = []
        mode = "Qzpn"
        try:
            data = triang.genMesh(points,seglist,holelist,regionlist,
                                  pointattlist,segattlist, trilist, mode)

        except TypeError:
            pass
        else:
            self.failUnless(0 ==1,
                        'bad segment list did not raise error!')

    def testbad_segattlist(self):

        holelist = []
        regionlist = []

	points = [(0.0,0.0),(0.0,10.0),(3.0,0.0),(3.0,10.0)]
	pointattlist = [[],[],[],[]]
	seglist = [(0,1),(1,3),(3,2),(2,0)]
        segattlist = [0,0]
        trilist = []
        mode = "Qzpn"
        try:
            data = triang.genMesh(points,seglist,holelist,regionlist,
                                  pointattlist,segattlist, trilist, mode)

            self.failUnless(data['trianglelist'] ==[(1, 0, 2), (2, 3, 1)],
                        'trianglelist is wrong!')
        except TypeError:
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
        trilist = []
        mode = "QAzpq"
        data = triang.genMesh(points,seglist,holelist,regionlist,
                              pointattlist,segattlist, trilist, mode)

        self.failUnless(data['generatedtriangleattributelist'] ==[[77.0], [77.0], [77.0], [77.0]],
                        'triangleattributelist is wrong!')


    def BAD_SEG_FAULT_test_lone_verts(self):

        points = []
        seglist = []
        holelist = []
        regionlist = []

	points = [(0.0,0.0),(0.0,10.0),(3.0,0.0),(3.0,10.0),(0.0,10.0)]
	pointattlist = [[],[],[],[],[]]
	regionlist.append( (1.2,1.2,5.0) )
	seglist = [(0,4),(4,1),(1,3),(3,2),(2,0)]
        segattlist = [0,0,0,0,0]
        trilist = []
        mode = "QpznAa2000.1a"
        data = triang.genMesh(points,seglist,holelist,regionlist,
                              pointattlist,segattlist,trilist, mode)
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
        self.failUnless(data['lonepointlist'] ==[1],
                        'lonepointlist is wrong!')
        #print "", data['lonepointlist']

    def BAD_SEG_FAULT_test_lone_vertsII(self):

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
        trilist = []
        mode = "pznAa2000.1a"
        data = triang.genMesh(points,seglist,holelist,regionlist,
                              pointattlist,segattlist,trilist, mode)
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

if __name__ == "__main__":

    suite = unittest.makeSuite(triangTestCase,'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)

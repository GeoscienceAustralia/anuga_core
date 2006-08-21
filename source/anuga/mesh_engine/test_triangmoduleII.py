#!/usr/bin/env python
#

import sys

import unittest
import mesh_engine.triang as triang

"""
Note, this test is seperated from the other tests because if it is added
the triangle code seg faults.  Don't know why
"""
class triangTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def BAD_test_lone_vertsII(self):

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
        mode = "QpznAa2000.1a"
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

    def test_lone_vertsIII(self):

        points = []
        seglist = []
        holelist = []
        regionlist = []

	points = [(0.0,0.0),(0.0,0.0),(0.0,10.0),(0.0,10.0),(10.0,10.0),
                  (10.0,10.0),(0.0,10.0),(10.0,0.0)]

	pointattlist = []
        for point in points:
            pointattlist.append([])
	seglist = []
	segattlist = []
        for seg in seglist:
            segattlist.append(0)
        trilist = []
        mode = "Qzp"
        data = triang.genMesh(points,seglist,holelist,regionlist,
                              pointattlist,segattlist,trilist, mode)

        self.failUnless(data['lonepointlist'] ==[0,1,2,3,4,5,6,7],
                        'lonepointlist is wrong!')
if __name__ == "__main__":

    suite = unittest.makeSuite(triangTestCase,'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)

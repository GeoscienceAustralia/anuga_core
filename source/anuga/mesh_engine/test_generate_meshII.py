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


    def test_lone_verts(self):
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

if __name__ == "__main__":

    suite = unittest.makeSuite(triangTestCase,'test')
    runner = unittest.TextTestRunner(verbosity=2) #verbosity=2)
    runner.run(suite)

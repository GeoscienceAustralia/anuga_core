#!/usr/bin/env python
from os import sep
from sys import path

from Numeric import array

import unittest

# path.append('..' + sep + 'pymetis')

import metis

class TestMetis(unittest.TestCase):
    def setUp(self):
        pass

    def testHexmesh(self):
        # Hexagonal mesh
        #
        #   1---2  
        #  / \ / \
        # 6---0---3
        #  \ / \ /
        #   5---4
        #
        # Divided 3 ways
        # Calling order is: elements, verticies, edge list
        # element type, number parts
        edgecut, epart, npart = metis.partMeshNodal(6, 7,\
                                                    [0, 1, 2,\
                                                     0, 2, 3,\
                                                     0, 3, 4,\
                                                     0, 4, 5,\
                                                     0, 5, 6,\
                                                     0, 6, 1],\
                                                    1,\
                                                    3,)
        #print edgecut
        #print epart
        #print npart
        epart_expected = array([2, 2, 0, 0, 0, 0], 'i')
        npart_expected = array([0, 2, 2, 2, 0, 0, 0], 'i')
        self.assert_(edgecut == 5)
        for i in range(len(epart)):
            self.assert_(epart[i] == epart_expected[i])
        for i in range(len(npart)):
            self.assert_(npart[i] == npart_expected[i])

if __name__ == '__main__':
    unittest.main()


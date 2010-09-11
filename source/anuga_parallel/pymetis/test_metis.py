#!/usr/bin/env python
from os import sep
from sys import path

from numpy import array, allclose

import unittest

# path.append('..' + sep + 'pymetis')

import anuga_parallel.pymetis.metis  as metis

class TestMetis(unittest.TestCase):
    def setUp(self):
        pass

    def test_Hexmesh(self):
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
        assert allclose(epart, epart_expected)
        assert allclose(npart, npart_expected)        


    def test_Hexmesh2(self):
        # Hexagonal mesh
        #
        #   1---2  
        #  / \ / \
        # 6---0---3
        #  \ / \ /
        #   5---4
        #
        # Divided 2 ways
        # Calling order is: elements, verticies, edge list
        # element type, number parts
        edgecut, epart, npart = metis.partMeshNodal(6, 7,\
                                                    [0, 2, 1,\
                                                     0, 3, 2,\
                                                     0, 4, 3,\
                                                     0, 5, 4,\
                                                     0, 6, 5,\
                                                     0, 1, 6],\
                                                    1,\
                                                    2,)
        #print edgecut
        #print epart
        #print npart
        epart_expected = array([1, 0, 0, 0, 1, 1], 'i')
        npart_expected = array([0, 1, 1, 0, 0, 0, 1], 'i')
        self.assert_(edgecut == 5)
        assert allclose(epart, epart_expected)
        assert allclose(npart, npart_expected)


if __name__ == "__main__":
    suite = unittest.makeSuite(TestMetis,'test_')
    runner = unittest.TextTestRunner()
    runner.run(suite)



#!/usr/bin/env python

import tempfile
import unittest
import os
import tempfile
from os.path import splitext

import numpy as num

from anuga.load_mesh.loadASCII import *
from anuga.coordinate_transforms.geo_reference import Geo_reference
import anuga.load_mesh.loadASCII as loadASCII


class loadASCIITestCase(unittest.TestCase):
    def setUp(self):
        self.dict = {}
        self.dict['outline_segments'] = [(0, 1), (1, 2), (0, 2), (0, 3)]
        self.dict['outline_segment_tags'] = ['50', '40', '30', '20']
        self.dict['holes'] = [(0.2, 0.6)]
        self.dict['point_attributes'] = [[5, 2], [4, 2], [3, 2], [2, 2]]
        self.dict['regions'] = [(0.3, 0.3), (0.3, 0.4)]
        self.dict['region_tags'] = ['1.3', 'yeah']
        self.dict['region_max_areas'] = [36.0, -7.1]
        self.dict['points'] = [(0.0, 0.0), (0.0, 4.0), (4.0, 0.0), (1.0, 1.0)]
        self.dict['vertices'] = [(0.0, 0.0), (0.0, 4.0), (4.0, 0.0),
                                 (1.0, 1.0), (2.0, 2.0)]
        self.dict['triangles'] = [(3, 2, 4), (1, 0, 3), (3, 4,1), (2, 3, 0)]
        self.dict['segments'] = [(0, 1), (1, 4), (2, 0), (0, 3), (4, 2)]
        self.dict['triangle_tags'] = ['1.3', '1.3', '1.3', '1.3']
        self.dict['vertex_attributes'] = [[1.2, 2.], [1.2, 2.], [1.2, 2.],
                                          [1.2, 2.], [1.2, 3.]]
        self.dict['triangle_neighbors'] = [[-1, 2, 3], [3, 2, -1],
                                           [-1, 1, 0], [1, -1, 0]]
        self.dict['segment_tags'] = ['50', '40', '30', '20', '40']
        self.dict['vertex_attribute_titles'] = ['bed elevation', 'height']
        self.dict['geo_reference'] = Geo_reference(56, 1.9, 1.9)
        
        
        self.dict_1 = {}
        self.dict_1['outline_segments'] = [(0, 1), (1, 2), (0, 2), (0, 3)]
        self.dict_1['outline_segment_tags'] = ['50', '40', '30', '20']
        self.dict_1['holes'] = [(0.2, 0.6)]
        self.dict_1['point_attributes'] = [[5], [4], [3], [2]]
        self.dict_1['regions'] = [(0.3, 0.3), (0.3, 0.4)]
        self.dict_1['region_tags'] = ['1.3', 'yeah']
        self.dict_1['region_max_areas'] = [36.0, -7.1]
        self.dict_1['points'] = [(0.0, 0.0), (0.0, 4.0), (4.0, 0.0), (1.0, 1.0)]
        self.dict_1['vertices'] = [(0.0, 0.0), (0.0, 4.0), (4.0, 0.0),
                                 (1.0, 1.0), (2.0, 2.0)]
        self.dict_1['triangles'] = [(3, 2, 4), (1, 0, 3), (3, 4,1), (2, 3, 0)]
        self.dict_1['segments'] = [(0, 1), (1, 4), (2, 0), (0, 3), (4, 2)]
        self.dict_1['triangle_tags'] = ['1.3', '1.3', '1.3', '1.3']
        self.dict_1['vertex_attributes'] = [[1.2], [1.2], [1.2],
                                         [1.2], [1.2]]
        self.dict_1['triangle_neighbors'] = [[-1, 2, 3], [3, 2, -1],
                                           [-1, 1, 0], [1, -1, 0]]
        self.dict_1['segment_tags'] = ['50', '40', '30', '20', '40']
        self.dict_1['vertex_attribute_titles'] = ['height']
        self.dict_1['geo_reference'] = Geo_reference(56, 1.9, 1.9)        

        self.sparse_dict = {}
        self.sparse_dict['outline_segments'] = []
        self.sparse_dict['outline_segment_tags'] = []
        self.sparse_dict['holes'] = []
        self.sparse_dict['points'] = [(0.0, 0.0), (9, 8)]
        self.sparse_dict['point_attributes'] = [[], []] # points don't have to
                                                        # have attributes
        self.sparse_dict['regions'] = []
        self.sparse_dict['region_tags'] = []
        self.sparse_dict['region_max_areas'] = []

        self.sparse_dict['vertices'] = []
        self.sparse_dict['triangles'] = []
        self.sparse_dict['segments'] = []
        self.sparse_dict['triangle_tags'] = []
        self.sparse_dict['vertex_attributes'] = []
        self.sparse_dict['triangle_neighbors'] = []
        self.sparse_dict['segment_tags'] = []
        self.sparse_dict['vertex_attribute_titles'] = []

        self.blank_dict = {}
        self.blank_dict['outline_segments'] = []
        self.blank_dict['outline_segment_tags'] = []
        self.blank_dict['holes'] = []
        self.blank_dict['points'] = []
        self.blank_dict['point_attributes'] = []
        self.blank_dict['regions'] = []
        self.blank_dict['region_tags'] = []
        self.blank_dict['region_max_areas'] = []
        self.blank_dict['vertices'] = []
        self.blank_dict['triangles'] = []
        self.blank_dict['segments'] = []
        self.blank_dict['triangle_tags'] = []
        self.blank_dict['vertex_attributes'] = []
        self.blank_dict['triangle_neighbors'] = []
        self.blank_dict['segment_tags'] = []
        self.blank_dict['vertex_attribute_titles'] = []

        self.tri_dict = {}
        self.tri_dict['outline_segments'] = [[0, 1]]
        self.tri_dict['outline_segment_tags'] = ['']
        self.tri_dict['holes'] = []
        self.tri_dict['points'] = [(9, 8), (7, 8)]
        self.tri_dict['point_attributes'] = [[], []]
        self.tri_dict['regions'] = []
        self.tri_dict['region_tags'] = []
        self.tri_dict['region_max_areas'] = []
        self.tri_dict['vertices'] = [[9, 8], [7, 8], [4, 5]]
        self.tri_dict['triangles'] = [[0, 1, 2]]
        self.tri_dict['segments'] = [[0, 1]]
        self.tri_dict['triangle_tags'] = ['']
        self.tri_dict['vertex_attributes'] = None
        self.tri_dict['triangle_neighbors'] = [[0, 0, 0]]
        self.tri_dict['segment_tags'] = ['']
        self.tri_dict['vertex_attribute_titles'] = []

        self.seg_dict = {}
        self.seg_dict['outline_segments'] = [[0, 1]]
        self.seg_dict['outline_segment_tags'] = ['']
        self.seg_dict['holes'] = []
        self.seg_dict['points'] = [(9, 8), (7, 8)]
        self.seg_dict['point_attributes'] = [[], []]
        self.seg_dict['regions'] = [(5, 4)]
        self.seg_dict['region_tags'] = ['']
        self.seg_dict['region_max_areas'] = [-999]
        self.seg_dict['vertices'] = [(9, 8), (7, 8)]
        self.seg_dict['triangles'] = []
        self.seg_dict['segments'] = [[0, 1]]
        self.seg_dict['triangle_tags'] = []
        self.seg_dict['vertex_attributes'] = None
        self.seg_dict['triangle_neighbors'] = []
        self.seg_dict['segment_tags'] = ['']
        self.seg_dict['vertex_attribute_titles'] = []

        self.reg_dict = {}
        self.reg_dict['outline_segments'] = [[0, 1]]
        self.reg_dict['outline_segment_tags'] = ['']
        self.reg_dict['holes'] = []
        self.reg_dict['points'] = [(9, 8), (7, 8)]
        self.reg_dict['point_attributes'] = [[], []]
        self.reg_dict['regions'] = [(5, 4)]
        self.reg_dict['region_tags'] = ['']
        self.reg_dict['region_max_areas'] = []
        self.reg_dict['vertices'] = [(9, 8), (7, 8)]
        self.reg_dict['triangles'] = []
        self.reg_dict['segments'] = [[0, 1]]
        self.reg_dict['triangle_tags'] = []
        self.reg_dict['vertex_attributes'] = [[], []]
        self.reg_dict['triangle_neighbors'] = []
        self.reg_dict['segment_tags'] = ['']
        self.reg_dict['vertex_attribute_titles'] = []

        self.triangle_tags_dict = {}
        self.triangle_tags_dict['outline_segments'] = [(0, 1), (1, 2),
                                                       (0, 2), (0, 3)]
        self.triangle_tags_dict['outline_segment_tags'] = ['50', '40',
                                                           '30', '20']
        self.triangle_tags_dict['holes'] = [(0.2, 0.6)]
        self.triangle_tags_dict['point_attributes'] = [[5, 2], [4, 2],
                                                       [3, 2], [2,2]]
        self.triangle_tags_dict['regions'] = [(0.3, 0.3), (0.3, 0.4)]
        self.triangle_tags_dict['region_tags'] = ['1.3', 'yeah']
        self.triangle_tags_dict['region_max_areas'] = [36.0, -7.1]
        self.triangle_tags_dict['points'] = [(0.0, 0.0), (0.0, 4.0),
                                             (4.0, 0.0), (1.0, 1.0)]
        self.triangle_tags_dict['vertices'] = [(0.0, 0.0), (0.0, 4.0),
                                               (4.0, 0.0), (1.0, 1.0),
                                               (2.0, 2.0)]
        self.triangle_tags_dict['triangles'] = [(3, 2, 4), (1, 0, 3),
                                                (3, 4, 1), (2, 3, 0)]
        self.triangle_tags_dict['segments'] = [(0, 1), (1, 4), (2, 0),
                                               (0, 3), (4, 2)]
        self.triangle_tags_dict['triangle_tags'] = ['yeah', '1.3', '1.3', '']
        self.triangle_tags_dict['vertex_attributes'] = [[1.2,2.], [1.2,2.],
                                                        [1.2,2.], [1.2,2.],
                                                        [1.2,3.]]
        self.triangle_tags_dict['triangle_neighbors'] = [[-1, 2, 3], [3, 2, -1],
                                                         [-1, 1, 0], [1, -1, 0]]
        self.triangle_tags_dict['segment_tags'] = ['50', '40', '30', '20', '40']
        self.triangle_tags_dict['vertex_attribute_titles'] = ['bed elevation',
                                                              'height']
        self.triangle_tags_dict['geo_reference'] = Geo_reference(56, 1.9, 1.9)

    def tearDown(self):
        pass

  ############### .TSH ##########
    def test_export_mesh_file(self):
        meshDict = self.dict
        fileName = tempfile.mktemp('.tsh')
        export_mesh_file(fileName, meshDict)
        loadedDict = import_mesh_file(fileName)

        self.assertTrue(num.alltrue(num.array(meshDict['vertices']) ==
                                    num.array(loadedDict['vertices'])),
                        'test_export_mesh_file failed. Test 1')
        self.assertTrue(num.alltrue(num.array(meshDict['triangles']) ==
                                    num.array(loadedDict['triangles'])),
                        'test_export_mesh_file failed. Test 2')
        self.assertTrue(num.alltrue(num.array(meshDict['segments']) ==
                                    num.array(loadedDict['segments'])),
                        'test_export_mesh_file failed. Test 3')
        self.assertTrue(num.alltrue(num.array(meshDict['triangle_tags']) ==
                                    num.array(loadedDict['triangle_tags'])),
                        'test_export_mesh_file failed. Test 4')

        self.assertTrue(meshDict['vertex_attributes'] ==
                        loadedDict['vertex_attributes'],
                        'test_export_mesh_file failed. Test 5')
        self.assertTrue(num.alltrue(num.array(meshDict['triangle_neighbors']) ==
                                    num.array(loadedDict['triangle_neighbors'])),
                        'test_export_mesh_file failed. Test 6')
        self.assertTrue(num.alltrue(num.array(meshDict['segment_tags']) ==
                                    num.array(loadedDict['segment_tags'])),
                        'test_export_mesh_file failed. Test 7')
        self.assertTrue(num.alltrue(num.array(meshDict['vertex_attribute_titles']) ==
                                    num.array(loadedDict['vertex_attribute_titles'])),
                        'test_export_mesh_file failed. Test 8')
        self.assertTrue(num.alltrue(num.array(meshDict['geo_reference']) ==
                                    num.array(loadedDict['geo_reference'])),
                        'test_export_mesh_file failed. Test 9')

        os.remove(fileName)
        
        
    def test_export_mesh_1_file(self):
        meshDict = self.dict_1
        fileName = tempfile.mktemp('.tsh')
        export_mesh_file(fileName, meshDict)
        loadedDict = import_mesh_file(fileName)
        

        self.assertTrue(num.alltrue(num.array(meshDict['vertices']) ==
                                    num.array(loadedDict['vertices'])),
                        'test_export_mesh_file failed. Test 1')
        self.assertTrue(num.alltrue(num.array(meshDict['triangles']) ==
                                    num.array(loadedDict['triangles'])),
                        'test_export_mesh_file failed. Test 2')
        self.assertTrue(num.alltrue(num.array(meshDict['segments']) ==
                                    num.array(loadedDict['segments'])),
                        'test_export_mesh_file failed. Test 3')
        self.assertTrue(num.alltrue(num.array(meshDict['triangle_tags']) ==
                                    num.array(loadedDict['triangle_tags'])),
                        'test_export_mesh_file failed. Test 4')

        self.assertTrue(meshDict['vertex_attributes'] ==
                        loadedDict['vertex_attributes'],
                        'test_export_mesh_file failed. Test 5')
        self.assertTrue(num.alltrue(num.array(meshDict['triangle_neighbors']) ==
                                    num.array(loadedDict['triangle_neighbors'])),
                        'test_export_mesh_file failed. Test 6')
        self.assertTrue(num.alltrue(num.array(meshDict['segment_tags']) ==
                                    num.array(loadedDict['segment_tags'])),
                        'test_export_mesh_file failed. Test 7')
        self.assertTrue(num.alltrue(num.array(meshDict['vertex_attribute_titles']) ==
                                    num.array(loadedDict['vertex_attribute_titles'])),
                        'test_export_mesh_file failed. Test 8')
        self.assertTrue(num.alltrue(num.array(meshDict['geo_reference']) ==
                                    num.array(loadedDict['geo_reference'])),
                        'test_export_mesh_file failed. Test 9')


        #os.remove(fileName)        

    def test_read_write_tsh_file(self):
        dict = self.dict.copy()
        fileName = tempfile.mktemp('.tsh')
        export_mesh_file(fileName, dict)
        loaded_dict = import_mesh_file(fileName)
        os.remove(fileName)
        dict = self.dict
        self.check_mesh_dicts(loaded_dict, dict, 'test_read_write_tsh_file')

    def test_read_write_tsh_fileII(self):
        dict = self.sparse_dict.copy()
        fileName = tempfile.mktemp('.tsh')
        export_mesh_file(fileName, dict)
        loaded_dict = import_mesh_file(fileName)
        dict = self.sparse_dict
        self.check_mesh_dicts(loaded_dict, dict, 'test_read_write_tsh_fileII')
        os.remove(fileName)

    def test_read_write_tsh_fileIII(self):
        dict = self.blank_dict.copy()
        fileName = tempfile.mktemp('.tsh')
        export_mesh_file(fileName, dict)
        loaded_dict = import_mesh_file(fileName)
        os.remove(fileName)
        dict = self.blank_dict
        self.check_mesh_dicts(loaded_dict, dict, 'test_read_write_tsh_fileIII')

    def test_read_write_tsh_file4(self):
        dict = self.seg_dict.copy()
        fileName = tempfile.mktemp('.tsh')
        export_mesh_file(fileName, dict)
        loaded_dict = import_mesh_file(fileName)
        os.remove(fileName)
        dict = self.seg_dict
        self.check_mesh_dicts(loaded_dict, dict, 'test_read_write_tsh_file4')

    def test_read_write_tsh_file5(self):
        dict = self.triangle_tags_dict.copy()
        fileName = tempfile.mktemp('.tsh')
        export_mesh_file(fileName, dict)
        loaded_dict = import_mesh_file(fileName)
        dict = self.triangle_tags_dict
        self.check_mesh_dicts(loaded_dict, dict, 'test_read_write_tsh_file5')
        os.remove(fileName)

    def test_read_write_tsh_file6(self):
        dict = self.tri_dict.copy()
        fileName = tempfile.mktemp('.tsh')
        export_mesh_file(fileName, dict)
        loaded_dict = import_mesh_file(fileName)
        dict = self.tri_dict
        self.check_mesh_dicts(loaded_dict, dict, 'test_read_write_tsh_file6')
        os.remove(fileName)

########################## BAD .TSH ##########################

    def test_load_bad_no_file_tsh(self):
        fileName = tempfile.mktemp('.tsh')
        # use self.faileUnlessRaises(IOError, import_mesh_file(fileName))
        try:
            dict = import_mesh_file(fileName)
        except IOError:
            pass
        else:
            self.assertTrue(0 == 1, 'imaginary file did not raise error!')

    def test_read_write_tsh_file_bad(self):
        dict = self.tri_dict.copy()
        fileName = tempfile.mktemp('.xxx')
        try:
            export_mesh_file(fileName, dict)
        except IOError:
            pass
        else:
            self.assertTrue(0 == 1, 'bad tsh file did not raise error!')

    def test_import_tsh_bad(self):
        fileName = tempfile.mktemp('.tsh')
        file = open(fileName, 'w')
        #   this is  a bad tsh file
        file.write('elevn\n\
1.0 what \n\
0.0 the \n\
1.0 !!! \n')
        file.close()
        try:
            dict = import_mesh_file(fileName)
        except IOError:
            pass
        else:
            self.fail('bad tsh file did not raise error!')
        os.remove(fileName)

    def test_import_tsh3(self):
        fileName = tempfile.mktemp('.tsh')
        file = open(fileName, 'w')
        file.write('1.0 \n\
showme1.0 0.0 10.0 \n\
0.0 1.0\n\
13.0 \n')
        file.close()
        try:
            dict = import_mesh_file(fileName)
        except IOError:
            pass
        else:
            self.fail('bad tsh file did not raise error!')
        os.remove(fileName)

  ############### .MSH ##########

    def test_read_write_msh_file1(self):
        dict = self.dict.copy()
        fileName = tempfile.mktemp('.msh')
        export_mesh_file(fileName, dict)
        loaded_dict = loadASCII._read_msh_file(fileName)
        os.remove(fileName)
        dict = self.dict
        self.check_mesh_dicts(loaded_dict, dict, 'test_read_write_msh_file1')

    def test_read_write_msh_fileII(self):
        dict = self.sparse_dict.copy()
        fileName = tempfile.mktemp('.msh')
        export_mesh_file(fileName, dict)
        loaded_dict = loadASCII._read_msh_file(fileName)
        os.remove(fileName)
        dict = self.sparse_dict
        self.check_mesh_dicts(loaded_dict, dict, 'test_read_write_msh_fileII')

    def test_read_write_msh_fileIII(self):
        dict = self.blank_dict.copy()
        fileName = tempfile.mktemp('.msh')
        export_mesh_file(fileName, dict)
        loaded_dict = loadASCII._read_msh_file(fileName)
        os.remove(fileName)
        dict = self.blank_dict
        self.check_mesh_dicts(loaded_dict, dict, 'test_read_write_msh_fileIII')

    def test_read_write_msh_file4(self):
        dict = self.seg_dict.copy()
        fileName = tempfile.mktemp('.msh')
        export_mesh_file(fileName, dict)
        loaded_dict = loadASCII._read_msh_file(fileName)
        os.remove(fileName)
        dict = self.seg_dict
        self.check_mesh_dicts(loaded_dict, dict, 'test_read_write_msh_file4')

    def test_read_write_msh_file5(self):
        dict = self.triangle_tags_dict.copy()
        fileName = tempfile.mktemp('.msh')
        export_mesh_file(fileName, dict)
        loaded_dict = loadASCII._read_msh_file(fileName)
        os.remove(fileName)
        dict = self.triangle_tags_dict
        self.check_mesh_dicts(loaded_dict, dict, 'test_read_write_msh_file5')

    def test_read_write_msh_file6(self):
        dict = self.tri_dict.copy()
        fileName = tempfile.mktemp('.msh')
        export_mesh_file(fileName, dict)
        loaded_dict = loadASCII._read_msh_file(fileName)
        os.remove(fileName)
        dict = self.tri_dict
        self.check_mesh_dicts(loaded_dict, dict, 'test_read_write_msh_file6')

    def check_mesh_dicts(self, loaded_dict, dict, fail_string):
        assert num.allclose(num.array(loaded_dict['points']),
                            num.array(dict['points']))
        assert num.allclose(num.array(loaded_dict['point_attributes']),
                            num.array(dict['point_attributes']))
        assert num.allclose(num.array(loaded_dict['outline_segments']),
                            num.array(dict['outline_segments']))

        self.assertTrue(loaded_dict['outline_segment_tags'] ==
                        dict['outline_segment_tags'],
                        fail_string + ' failed!! Test 4')

        assert num.allclose(num.array(loaded_dict['regions']),
                            num.array(dict['regions']))

        self.assertTrue(loaded_dict['region_tags'] == dict['region_tags'],
                        fail_string + ' failed!! Test 5')

        assert num.allclose(num.array(loaded_dict['region_max_areas']),
                            num.array(dict['region_max_areas']))

        assert num.allclose(num.array(loaded_dict['holes']),
                            num.array(dict['holes']))

        assert num.allclose(num.array(dict['vertices']),
                            num.array(loaded_dict['vertices']))

        assert num.allclose(num.array(dict['triangles']),
                            num.array(loaded_dict['triangles']))

        assert num.allclose(num.array(dict['segments']),
                            num.array(loaded_dict['segments']))
        for ob, ldob in map(None, dict['triangle_tags'],
                            loaded_dict['triangle_tags']):
            msg = ('ob=\n%s\nshould be same as ldob=\n%s' % (str(ob), str(ldob)))
            self.assertTrue(ob == ldob,
                            fail_string + ' failed\n' + msg)

        # A bit hacky
        self.assertTrue(num.alltrue(loaded_dict['vertex_attributes'] ==
                                    dict['vertex_attributes']) or
                                    (loaded_dict['vertex_attributes'] ==
                                    None and
                                    dict['vertex_attributes'] == []),
                        fail_string + ' failed!! Test vertex_attributes')

        assert num.allclose(num.array(dict['triangle_neighbors']),
                            num.array(loaded_dict['triangle_neighbors']))

        for seg, ldseg in map(None, dict['segment_tags'],
                              loaded_dict['segment_tags']):
            msg = ('seg=\n"%s"\nshould be same as ldseg=\n"%s"'
                   % (str(seg), str(ldseg)))
            self.assertTrue(seg == ldseg, fail_string + ' failed\n' + msg)

        try:
            assert num.alltrue(num.array(dict['vertex_attribute_titles'])==
                                num.array(loaded_dict['vertex_attribute_titles']))
        except IndexError:
            self.assertTrue(num.alltrue(num.array(loaded_dict['vertex_attribute_titles'])
                            == num.array(dict['vertex_attribute_titles'])),
                            fail_string + ' failed!! Test 8')

        try:
            msg = ("loaded_dict['geo_reference']=\n%s\n"
                   "should be same as dict['geo_reference']=%s"
                  % (str(loaded_dict['geo_reference']),
                     str(dict['geo_reference'])))
            self.assertTrue(loaded_dict['geo_reference'] ==
                            dict['geo_reference'],
                            fail_string + ' failed\n' + msg)
        except KeyError:        #??# 2 lines below??
            msg = ("'dict' has no key 'geo_reference' "
                   "but loaded_dict['geo_reference'] isn't None")
            self.assertTrue(not dict.has_key('geo_reference') and
                                loaded_dict['geo_reference'] == None,
                            fail_string + ' failed\n' + msg)

########################## BAD .MSH ##########################

    def test_load_bad_no_file_msh(self):
        fileName = tempfile.mktemp('.msh')
        try:
            dict = import_mesh_file(fileName)
        except IOError:
            pass
        else:
            self.fail('imaginary file did not raise error!')

    def throws_error_2_screen_import_mesh_bad(self):
        fileName = tempfile.mktemp('.msh')
        file = open(fileName, 'w')
        # this is  a bad tsh file
        file.write('elevn\n\
1.0 what \n\
0.0 the \n\
1.0 !!! \n')
        file.close()
        try:
            dict = import_mesh_file(fileName)
        except IOError:
            pass
        else:
            self.fail('bad msh file did not raise error!')
        os.remove(fileName)

################################################################################

if __name__ == '__main__':
    suite = unittest.makeSuite(loadASCIITestCase,'test')
    runner = unittest.TextTestRunner() #verbosity=2)
    runner.run(suite)


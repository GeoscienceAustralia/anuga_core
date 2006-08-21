#!/usr/bin/env python
#

import tempfile
import unittest

import os
import tempfile

from os.path import splitext

from Numeric import array, allclose 
    
from loadASCII import *
from coordinate_transforms.geo_reference import Geo_reference
import loadASCII

class loadASCIITestCase(unittest.TestCase):
    def setUp(self):
        self.dict ={}
        self.dict['outline_segments'] = [(0, 1), (1, 2), (0, 2), (0, 3)]
        self.dict['outline_segment_tags'] = ['50', '40', '30', '20']
        self.dict['holes'] = [(0.2, 0.6)]
        self.dict['point_attributes'] = [[5, 2], [4, 2], [3, 2], [2,2]]
        self.dict['regions'] = [(0.3, 0.3),(0.3, 0.4)]
        self.dict['region_tags'] = ['1.3', 'yeah']
        self.dict['region_max_areas'] = [36.0,-7.1]
        self.dict['points'] = [(0.0, 0.0), (0.0, 4.0), (4.0, 0.0), (1.0, 1.0)]
        self.dict['vertices'] = [(0.0, 0.0), (0.0, 4.0),
                                 (4.0, 0.0), (1.0, 1.0), (2.0, 2.0)]
        self.dict['triangles'] = [(3, 2, 4), (1, 0, 3),
                                             (3, 4,1), (2, 3, 0)]
        self.dict['segments'] = [(0, 1), (1, 4), (2, 0),
                                            (0, 3), (4, 2)]
        self.dict['triangle_tags'] = ['1.3', '1.3',
                                      '1.3', '1.3']
        self.dict['vertex_attributes'] = [[1.2,2.], [1.2,2.],
                                          [1.2,2.], [1.2,2.], [1.2,3.]]
        self.dict['triangle_neighbors'] = [[-1, 2, 3], [3, 2, -1],
                                           [-1, 1, 0], [1, -1, 0]]
        self.dict['segment_tags'] = ['50', '40', '30', '20', '40']
        self.dict['vertex_attribute_titles'] = ['bed elevation', 'height']
        self.dict['geo_reference'] = Geo_reference(56,1.9,1.9)
        
        self.sparse_dict ={}
        self.sparse_dict['outline_segments'] = []
        self.sparse_dict['outline_segment_tags'] = []
        self.sparse_dict['holes'] = []
        self.sparse_dict['points'] = [(0.0, 0.0),(9,8)]
        self.sparse_dict['point_attributes'] = [[],[]] # points don't have to have
                                                    # attributes
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
        
        self.blank_dict ={}
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
        
        self.tri_dict ={}
        self.tri_dict['outline_segments'] = [[0,1]]
        self.tri_dict['outline_segment_tags'] = ['']
        self.tri_dict['holes'] = []
        self.tri_dict['points'] = [(9,8),(7,8)]
        self.tri_dict['point_attributes'] = [[],[]] 
        self.tri_dict['regions'] = []
        self.tri_dict['region_tags'] = []
        self.tri_dict['region_max_areas'] = []
        self.tri_dict['vertices'] = [[9,8],[7,8], [4,5]]
        self.tri_dict['triangles'] = [[0,1,2]]
        self.tri_dict['segments'] = [[0,1]]
        self.tri_dict['triangle_tags'] = ['']
        self.tri_dict['vertex_attributes'] = [[],[],[]]
        self.tri_dict['triangle_neighbors'] = [[0,0,0]]
        self.tri_dict['segment_tags'] = ['']
        self.tri_dict['vertex_attribute_titles'] = []
        
        self.seg_dict ={}
        self.seg_dict['outline_segments'] = [[0,1]]
        self.seg_dict['outline_segment_tags'] = ['']
        self.seg_dict['holes'] = []
        self.seg_dict['points'] = [(9,8),(7,8)]
        self.seg_dict['point_attributes'] = [[],[]]  
        self.seg_dict['regions'] = [(5,4)]
        self.seg_dict['region_tags'] = [''] 
        self.seg_dict['region_max_areas'] = [-999]
        self.seg_dict['vertices'] = [(9,8),(7,8)]
        self.seg_dict['triangles'] = []
        self.seg_dict['segments'] = [[0,1]]
        self.seg_dict['triangle_tags'] = []
        self.seg_dict['vertex_attributes'] = [[],[]]
        self.seg_dict['triangle_neighbors'] = []
        self.seg_dict['segment_tags'] = ['']
        self.seg_dict['vertex_attribute_titles'] = []
        
        self.reg_dict ={}
        self.reg_dict['outline_segments'] = [[0,1]]
        self.reg_dict['outline_segment_tags'] = ['']
        self.reg_dict['holes'] = []
        self.reg_dict['points'] = [(9,8),(7,8)]
        self.reg_dict['point_attributes'] = [[],[]] 
        self.reg_dict['regions'] = [(5,4)]
        self.reg_dict['region_tags'] = ['']
        self.reg_dict['region_max_areas'] = []
        self.reg_dict['vertices'] = [(9,8),(7,8)]
        self.reg_dict['triangles'] = []
        self.reg_dict['segments'] = [[0,1]]
        self.reg_dict['triangle_tags'] = []
        self.reg_dict['vertex_attributes'] = [[],[]]
        self.reg_dict['triangle_neighbors'] = []
        self.reg_dict['segment_tags'] = ['']
        self.reg_dict['vertex_attribute_titles'] = []
        
        self.triangle_tags_dict ={}
        self.triangle_tags_dict['outline_segments'] = [(0, 1), (1, 2), (0, 2), (0, 3)]
        self.triangle_tags_dict['outline_segment_tags'] = ['50', '40', '30', '20']
        self.triangle_tags_dict['holes'] = [(0.2, 0.6)]
        self.triangle_tags_dict['point_attributes'] = [[5, 2], [4, 2], [3, 2], [2,2]]
        self.triangle_tags_dict['regions'] = [(0.3, 0.3),(0.3, 0.4)]
        self.triangle_tags_dict['region_tags'] = ['1.3', 'yeah']
        self.triangle_tags_dict['region_max_areas'] = [36.0,-7.1]
        self.triangle_tags_dict['points'] = [(0.0, 0.0), (0.0, 4.0), (4.0, 0.0), (1.0, 1.0)]
        self.triangle_tags_dict['vertices'] = [(0.0, 0.0), (0.0, 4.0),
                                 (4.0, 0.0), (1.0, 1.0), (2.0, 2.0)]
        self.triangle_tags_dict['triangles'] = [(3, 2, 4), (1, 0, 3),
                                             (3, 4,1), (2, 3, 0)]
        self.triangle_tags_dict['segments'] = [(0, 1), (1, 4), (2, 0),
                                            (0, 3), (4, 2)]
        self.triangle_tags_dict['triangle_tags'] = ['yeah', '1.3',
                                      '1.3', '']
        self.triangle_tags_dict['vertex_attributes'] = [[1.2,2.], [1.2,2.],
                                          [1.2,2.], [1.2,2.], [1.2,3.]]
        self.triangle_tags_dict['triangle_neighbors'] = [[-1, 2, 3], [3, 2, -1],
                                           [-1, 1, 0], [1, -1, 0]]
        self.triangle_tags_dict['segment_tags'] = ['50', '40', '30', '20', '40']
        self.triangle_tags_dict['vertex_attribute_titles'] = ['bed elevation', 'height']
        self.triangle_tags_dict['geo_reference'] = Geo_reference(56,1.9,1.9)
        
    def tearDown(self):
        pass

  ############### .TSH ##########      
    def test_export_mesh_file(self):
        import os
        import tempfile
        
        meshDict = self.dict
        fileName = tempfile.mktemp(".tsh")
        export_mesh_file(fileName, meshDict)
        loadedDict = import_mesh_file(fileName)
        
        #print "*(*( meshDict"
        #print meshDict
        #print "*(*(  loadedDcit"
        #print loadedDict
        #print "*(*("
        
        self.failUnless(array(meshDict['vertices'])  ==
                        array(loadedDict['vertices']),
                         'test_export_mesh_file failed. Test 1')
        self.failUnless(array(meshDict['triangles'])  ==
                        array(loadedDict['triangles']),
                         'test_export_mesh_file failed. Test 2')
        self.failUnless(array(meshDict['segments'])  ==
                        array(loadedDict['segments']),
                         'test_export_mesh_file failed. Test 3')
        self.failUnless(array(meshDict['triangle_tags'])  ==
                        array(loadedDict['triangle_tags']),
                         'test_export_mesh_file failed. Test 4')
        
        self.failUnless(meshDict['vertex_attributes']  ==
                        loadedDict['vertex_attributes'],
                         'test_export_mesh_file failed. Test 5')
        self.failUnless(array(meshDict['triangle_neighbors'])  ==
                        array(loadedDict['triangle_neighbors']),
                         'test_export_mesh_file failed. Test 6')
        self.failUnless(array(meshDict['segment_tags'])  ==
                        array(loadedDict['segment_tags']),
                         'test_export_mesh_file failed. Test 7')
        self.failUnless(array(meshDict['vertex_attribute_titles'])  ==
                        array(loadedDict['vertex_attribute_titles']),
                         'test_export_mesh_file failed. Test 8')
        self.failUnless(array(meshDict['geo_reference'])  ==
                        array(loadedDict['geo_reference']),
                         'test_export_mesh_file failed. Test 9')
            
        os.remove(fileName)
 
    def test_read_write_tsh_file(self):
        dict = self.dict.copy()
        fileName = tempfile.mktemp(".tsh")
        export_mesh_file(fileName,dict)
        loaded_dict = import_mesh_file(fileName)
        os.remove(fileName)
        dict = self.dict
        #print "*********************"
        #print dict
        #print "**loaded_dict*******************"
        #print loaded_dict
        #print "*********************"       
        self.check_mesh_dicts(loaded_dict,dict, 'test_read_write_msh_file')
        
    def test_read_write_tsh_fileII(self):
        dict = self.sparse_dict.copy()
        fileName = tempfile.mktemp(".tsh")
        export_mesh_file(fileName,dict)
        loaded_dict = import_mesh_file(fileName)
        dict = self.sparse_dict   
        self.check_mesh_dicts(loaded_dict,dict, 'test_read_write_msh_file')
        os.remove(fileName)
        
    def test_read_write_tsh_fileIII(self):
        dict = self.blank_dict.copy()
        fileName = tempfile.mktemp(".tsh")
        export_mesh_file(fileName,dict)
        loaded_dict = import_mesh_file(fileName)
        os.remove(fileName)
        dict = self.blank_dict
        #print "*********************"
        #print dict
        #print "**loaded_dict*******************"
        #print loaded_dict
        #print "*********************"       
        self.check_mesh_dicts(loaded_dict,dict, 'test_read_write_msh_fileIII')

    def test_read_write_tsh_file4(self):
        dict = self.seg_dict.copy()
        fileName = tempfile.mktemp(".tsh")
        export_mesh_file(fileName,dict)
        loaded_dict = import_mesh_file(fileName)
        os.remove(fileName)
        dict = self.seg_dict   
        self.check_mesh_dicts(loaded_dict,dict, 'test_read_write_msh_file4')

    def test_read_write_tsh_file5(self):
        dict = self.triangle_tags_dict.copy()
        fileName = tempfile.mktemp(".tsh")
        export_mesh_file(fileName,dict)
        loaded_dict = import_mesh_file(fileName)
        dict = self.triangle_tags_dict
        #print "*********************"
        #print dict
        #print "**loaded_dict*******************"
        #print loaded_dict
        #print "*********************"       
        self.check_mesh_dicts(loaded_dict,dict, 'test_read_write_msh_file5')  
        os.remove(fileName)

    def test_read_write_tsh_file6(self):
        dict = self.tri_dict.copy()
        fileName = tempfile.mktemp(".tsh")
        export_mesh_file(fileName,dict)
        loaded_dict = import_mesh_file(fileName)
        dict = self.tri_dict    
        self.check_mesh_dicts(loaded_dict,dict, 'test_read_write_msh_file6')  
        os.remove(fileName)
        
########################## BAD .TSH ##########################

    def test_load_bad_no_file_tsh(self):
        import os
        import tempfile
       
        fileName = tempfile.mktemp(".tsh")
        #print fileName
        try:
            dict = import_mesh_file(fileName)
        except IOError:
            pass
        else:
            self.failUnless(0 ==1,
                        'imaginary file did not raise error!')
         
    def test_read_write_tsh_file_bad(self):
        dict = self.tri_dict.copy()
        fileName = tempfile.mktemp(".xxx")
        try:
            export_mesh_file(fileName,dict)
        except IOError:
            pass
        else:
            self.failUnless(0 ==1,
                        'bad tsh file did not raise error!')       
        
    def test_import_tsh_bad(self):
        import os
        import tempfile
       
        fileName = tempfile.mktemp(".tsh")
        file = open(fileName,"w")
        #   this is  a bad tsh file
        file.write("elevn\n\
1.0 what \n\
0.0 the \n\
1.0 !!! \n")
        file.close()
        #print fileName
        try:
            dict = import_points_file(fileName,delimiter = ' ')
        except IOError:
            pass
        else:
            self.failUnless(0 ==1,
                        'bad tsh file did not raise error!')       
        os.remove(fileName)

    def test_import_tsh3(self):
        import os
        import tempfile
       
        fileName = tempfile.mktemp(".tsh")
        file = open(fileName,"w")
        file.write("1.0 \n\
showme1.0 0.0 10.0 \n\
0.0 1.0\n\
13.0 \n")
        file.close()
        #print fileName
        try:
            dict = import_points_file(fileName,delimiter = ' ')
        except IOError:
            pass
        else:
            self.failUnless(0 ==1,
                        'bad tsh file did not raise error!')
        
        os.remove(fileName)         

            
  ############### .MSH ##########
        
    def test_read_write_msh_file(self):
        dict = self.dict.copy()
        fileName = tempfile.mktemp(".msh")
        export_mesh_file(fileName,dict)
        loaded_dict = loadASCII._read_msh_file(fileName)
        os.remove(fileName)
        dict = self.dict
        #print "*********************"
        #print dict
        #print "**loaded_dict*******************"
        #print loaded_dict
        #print "*********************"
        self.check_mesh_dicts(loaded_dict,dict,'test_read_write_msh_file') 

    def test_read_write_msh_fileII(self):
        dict = self.sparse_dict.copy()
        fileName = tempfile.mktemp(".msh")
        export_mesh_file(fileName,dict)
        loaded_dict = loadASCII._read_msh_file(fileName)
        os.remove(fileName)
        dict = self.sparse_dict
        #print "*********************"
        #print dict
        #print "**loaded_dict*******************"
        #print loaded_dict
        #print "*********************"       
        self.check_mesh_dicts(loaded_dict,dict, 'test_read_write_msh_fileII')
              
    def test_read_write_msh_fileIII(self):
        dict = self.blank_dict.copy()
        fileName = tempfile.mktemp(".msh")
        export_mesh_file(fileName,dict)
        loaded_dict = loadASCII._read_msh_file(fileName)
        os.remove(fileName)
        dict = self.blank_dict
        #print "*********************"
        #print dict
        #print "**loaded_dict*******************"
        #print loaded_dict
        #print "*********************"       
        self.check_mesh_dicts(loaded_dict,dict, 'test_read_write_msh_fileIII')
       
    def test_read_write_msh_file4(self):
        dict = self.seg_dict.copy()
        fileName = tempfile.mktemp(".msh")
        export_mesh_file(fileName,dict)
        loaded_dict = loadASCII._read_msh_file(fileName)
        os.remove(fileName)
        dict = self.seg_dict
        #print "*********************"
        #print dict
        #print "**loaded_dict*******************"
        #print loaded_dict
        #print "*********************" 
        self.check_mesh_dicts(loaded_dict,dict, 'test_read_write_msh_fileIII')
        
    def test_read_write_msh_file5(self):
        dict = self.triangle_tags_dict.copy()
        fileName = tempfile.mktemp(".msh")
        export_mesh_file(fileName,dict)
        loaded_dict = loadASCII._read_msh_file(fileName)
        os.remove(fileName)
        dict = self.triangle_tags_dict
        #print "msh_file5*********************"
        #print dict
        #print "**loaded_dict*******************"
        #print loaded_dict
        #print "*********************" 
        self.check_mesh_dicts(loaded_dict,dict, 'test_read_write_msh_fileIII')
                                       
        
    def test_read_write_msh_file6(self):
        dict = self.tri_dict.copy()
        fileName = tempfile.mktemp(".msh")
        export_mesh_file(fileName,dict)
        loaded_dict = loadASCII._read_msh_file(fileName)
        os.remove(fileName)
        dict = self.tri_dict
        #print "*********************"
        #print dict
        #print "**loaded_dict*******************"
        #print loaded_dict
        #print "*********************" 
        self.check_mesh_dicts(loaded_dict,dict, 'test_read_write_msh_fileIII')
                         
    def check_mesh_dicts(self, loaded_dict, dict, fail_string ):
        assert allclose(array(loaded_dict['points']),
                        array(dict['points']))
                        
	assert allclose(array(loaded_dict['point_attributes']),
                        array(dict['point_attributes']))
        assert allclose(array(loaded_dict['outline_segments']),
                        array(dict['outline_segments']))
        
        self.failUnless(loaded_dict['outline_segment_tags']  ==
                        dict['outline_segment_tags'],
                        fail_string + ' failed!! Test 4')
        
            
        assert allclose(array(loaded_dict['regions']),
                        array(dict['regions']))
        self.failUnless(loaded_dict['region_tags']  ==
                        dict['region_tags'],
                        fail_string + ' failed!! Test 5')
        
        assert allclose(array(loaded_dict['region_max_areas']),
                        array(dict['region_max_areas']))
 
        assert allclose(array(loaded_dict['holes']),
                        array(dict['holes']))
 
        assert allclose(array(dict['vertices']),
                        array(loaded_dict['vertices']))
 
        assert allclose(array(dict['triangles']),
                        array(loaded_dict['triangles']))
 
        assert allclose(array(dict['segments']),
                        array(loaded_dict['segments']))
        for ob, ldob in map(None,dict['triangle_tags'],
                              loaded_dict['triangle_tags']):
            self.failUnless(ob  == ldob,
                        fail_string + ' failed!! Test triangle_tags')
        self.failUnless(loaded_dict['vertex_attributes']  ==
                        dict['vertex_attributes'],
                        fail_string + ' failed!! Test vertex_attributes')
 
        assert allclose(array(dict['triangle_neighbors']),
                        array(loaded_dict['triangle_neighbors']))

        for seg, ldseg in map(None,dict['segment_tags'],
                              loaded_dict['segment_tags']):
            self.failUnless(seg  == ldseg,
                        fail_string + ' failed!! Test 8')
        try:
            assert allclose(array(dict['vertex_attribute_titles']),
                            array(loaded_dict['vertex_attribute_titles']))
        except TypeError:
            self.failUnless(array(loaded_dict['vertex_attribute_titles'])  ==
                        array(dict['vertex_attribute_titles']),
                        fail_string + ' failed!! Test 8')
        try:    
            self.failUnless(loaded_dict['geo_reference']  ==
                        dict['geo_reference'] ,
                        fail_string + ' failed!! Test geo_reference')
        except KeyError:          
            self.failUnless(not dict.has_key('geo_reference' and
                         loaded_dict['geo_reference'] == None)   ,
                        fail_string + ' failed!! Test geo_reference')
  
########################## BAD .MSH ##########################          

    def test_load_bad_no_file_msh(self):
        import os
        import tempfile
       
        fileName = tempfile.mktemp(".msh")
        #print fileName
        try:
            dict = import_mesh_file(fileName)
        except IOError:
            pass
        else:
            self.failUnless(0 ==1,
                        'imaginary file did not raise error!')
            
    def test_import_msh_bad(self):
        import os
        import tempfile
       
        fileName = tempfile.mktemp(".msh")
        file = open(fileName,"w")
        #   this is  a bad tsh file
        file.write("elevn\n\
1.0 what \n\
0.0 the \n\
1.0 !!! \n")
        file.close()
        #print fileName
        try:
            dict = import_points_file(fileName,delimiter = ' ')
        except IOError:
            pass
        else:
            self.failUnless(0 ==1,
                        'bad msh file did not raise error!')        
        os.remove(fileName)         
   
  ###################### .XYA ##############################
        
    def test_export_xya_file(self):
        dict = {}
        att_dict = {}
	dict['pointlist'] = array([[1.0, 0.0],[0.0, 1.0],[1.0, 0.0]])
	att_dict['elevation'] = array([10.0, 0.0, 10.4])
	att_dict['brightness'] = array([10.0, 0.0, 10.4])
        dict['attributelist'] = att_dict
        dict['geo_reference'] = Geo_reference(56,1.9,1.9)
        
        
        fileName = tempfile.mktemp(".xya")
        export_points_file(fileName, dict)
        dict2 = import_points_file(fileName)
        #print "fileName",fileName 
        os.remove(fileName)
        #print "dict2",dict2 
        
	assert allclose(dict2['pointlist'],[[1.0, 0.0],[0.0, 1.0],[1.0, 0.0]])
	assert allclose(dict2['attributelist']['elevation'], [10.0, 0.0, 10.4])
        answer = [10.0, 0.0, 10.4]
	assert allclose(dict2['attributelist']['brightness'], answer)
        #print "dict2['geo_reference']",dict2['geo_reference'] 
        self.failUnless(dict['geo_reference'] == dict2['geo_reference'],
                         'test_writepts failed. Test geo_reference')

    def test_export_xya_file2(self):
        dict = {}
        att_dict = {}
	dict['pointlist'] = array([[1.0, 0.0],[0.0, 1.0],[1.0, 0.0]])
	att_dict['elevation'] = array([10.0, 0.0, 10.4])
	att_dict['brightness'] = array([10.0, 0.0, 10.4])
        dict['attributelist'] = att_dict
        
        
        fileName = tempfile.mktemp(".xya")
        export_points_file(fileName, dict)
        dict2 = import_points_file(fileName)
        #print "fileName",fileName 
        os.remove(fileName)
        #print "dict2",dict2 
        
	assert allclose(dict2['pointlist'],[[1.0, 0.0],[0.0, 1.0],[1.0, 0.0]])
	assert allclose(dict2['attributelist']['elevation'], [10.0, 0.0, 10.4])
        answer = [10.0, 0.0, 10.4]
	assert allclose(dict2['attributelist']['brightness'], answer)

        
    def test_loadxya(self):
        """
        comma delimited
        """
       
        fileName = tempfile.mktemp(".xya")
        file = open(fileName,"w")
        file.write("elevation  , speed \n\
1.0, 0.0, 10.0, 0.0\n\
0.0, 1.0, 0.0, 10.0\n\
1.0, 0.0, 10.4, 40.0\n")
        file.close()
        #print fileName
        dict = import_points_file(fileName,delimiter = ',')
        os.remove(fileName)
	assert allclose(dict['pointlist'], [[1.0, 0.0],[0.0, 1.0],[1.0, 0.0]])
	assert allclose(dict['attributelist']['elevation'], [10.0, 0.0, 10.4])
	assert allclose(dict['attributelist']['speed'], [0.0, 10.0, 40.0])

    #FIXME - redundant test?
    def test_loadxy(self):
        """
        To test the mesh side of loading xya files.
        Not the loading of xya files
        """
        import os
        import tempfile
       
        fileName = tempfile.mktemp(".xya")
        file = open(fileName,"w")
        file.write("elevation speed \n\
1.0 0.0 10.0 0.0\n\
0.0 1.0 0.0 10.0\n\
1.0 0.0 10.4 40.0\n")
        file.close()
        #print fileName
        dict = import_points_file(fileName)
        os.remove(fileName)
	assert allclose(dict['pointlist'], [[1.0, 0.0],[0.0, 1.0],[1.0, 0.0]])
	assert allclose(dict['attributelist']['elevation'], [10.0, 0.0, 10.4])
	assert allclose(dict['attributelist']['speed'], [0.0, 10.0, 40.0])

               
    def test_loadxya2(self):
        """
        space delimited
        """
        import os
        import tempfile
       
        fileName = tempfile.mktemp(".xya")
        file = open(fileName,"w")
        file.write("  elevation   speed \n\
1.0 0.0 10.0 0.0\n\
0.0 1.0 0.0 10.0\n\
1.0 0.0 10.4 40.0\n")
        file.close()
        #print fileName
        dict = import_points_file(fileName,delimiter = ' ')
        os.remove(fileName)
	assert allclose(dict['pointlist'], [[1.0, 0.0],[0.0, 1.0],[1.0, 0.0]])
	assert allclose(dict['attributelist']['elevation'], [10.0, 0.0, 10.4])
	assert allclose(dict['attributelist']['speed'], [0.0, 10.0, 40.0])
     
    def test_loadxya3(self):
        """
        space delimited
        """
        import os
        import tempfile
       
        fileName = tempfile.mktemp(".xya")
        file = open(fileName,"w")
        file.write("  elevation   speed \n\
1.0 0.0 10.0 0.0\n\
0.0 1.0 0.0 10.0\n\
1.0 0.0 10.4 40.0\n\
#geocrap\n\
56\n\
56.6\n\
3\n")
        file.close()
        #print fileName
        dict = import_points_file(fileName,delimiter = ' ')
        os.remove(fileName)
	assert allclose(dict['pointlist'], [[1.0, 0.0],[0.0, 1.0],[1.0, 0.0]])
	assert allclose(dict['attributelist']['elevation'], [10.0, 0.0, 10.4])
	assert allclose(dict['attributelist']['speed'], [0.0, 10.0, 40.0])
     
        geo_reference = Geo_reference(56, 56.6, 3.0)
        
        self.failUnless(geo_reference == dict['geo_reference'],
                         'test_writepts failed. Test geo_reference')
        
########################## BAD .XYA ##########################
  
    def test_loadxy_bad_no_file_xya(self):
        import os
        import tempfile
       
        fileName = tempfile.mktemp(".xya")
        #print fileName
        try:
            dict = import_points_file(fileName,delimiter = ' ')
        except IOError:
            pass
        else:
            self.failUnless(0 ==1,
                        'imaginary file did not raise error!')
  
    def test_read_write_points_file_bad(self):
        dict = self.tri_dict.copy()
        fileName = tempfile.mktemp(".xxx")
        try:
            export_points_file(fileName,dict)
        except IOError:
            pass
        else:
            self.failUnless(0 ==1,
                        'bad points file extension did not raise error!')
     
    def test_read_write_points_file_bad2(self):
        dict = {}
        att_dict = {}
	dict['pointlist'] = array([[1.0, 0.0],[0.0, 1.0],[1.0, 0.0]])
	att_dict['elevation'] = array([10.0, 0.0, 10.4])
	att_dict['brightness'] = array([10.0, 0.0, 10.4])
        dict['attributelist'] = att_dict
        dict['geo_reference'] = Geo_reference(56,1.9,1.9)
        try:
            export_points_file("_???/yeah.xya",dict)
        except IOError:
            pass
        else:
            self.failUnless(0 ==1,
                        'bad points file extension did not raise error!')
                   
    def test_loadxy_bad(self):
        import os
        import tempfile
       
        fileName = tempfile.mktemp(".xya")
        file = open(fileName,"w")
        file.write("  elevation   \n\
1.0 0.0 10.0 0.0\n\
0.0 1.0 0.0 10.0\n\
1.0 0.0 10.4 40.0\n")
        file.close()
        #print fileName
        try:
            dict = import_points_file(fileName,delimiter = ' ')
        except IOError:
            pass
        else:
            self.failUnless(0 ==1,
                        'bad xya file did not raise error!')
        os.remove(fileName)
       
    def test_loadxy_bad2(self):
        import os
        import tempfile
       
        fileName = tempfile.mktemp(".xya")
        file = open(fileName,"w")
        file.write("elevation\n\
1.0 0.0 10.0 \n\
0.0 1.0\n\
1.0 \n")
        file.close()
        #print fileName
        try:
            dict = import_points_file(fileName,delimiter = ' ')
        except IOError:
            pass
        else:
            self.failUnless(0 ==1,
                        'bad xya file did not raise error!')    
        os.remove(fileName)
   
    def test_loadxy_bad3(self):
        """
        specifying wrong delimiter
        """
        import os
        import tempfile
       
        fileName = tempfile.mktemp(".xya")
        file = open(fileName,"w")
        file.write("  elevation  , speed \n\
1.0, 0.0, 10.0, 0.0\n\
0.0, 1.0, 0.0, 10.0\n\
1.0, 0.0, 10.4, 40.0\n")
        file.close()
        try:
            dict = import_points_file(fileName,delimiter = ' ')
        except IOError:
            pass
        else:
            self.failUnless(0 ==1,
                        'bad xya file did not raise error!')    
        os.remove(fileName)
     
    def test_loadxy_bad4(self):
        """
        specifying wrong delimiter
        """
        import os
        import tempfile
       
        fileName = tempfile.mktemp(".xya")
        file = open(fileName,"w")
        file.write("  elevation   speed \n\
1.0 0.0 10.0 0.0\n\
0.0 1.0 0.0 10.0\n\
1.0 0.0 10.4 40.0\n\
yeah")
        file.close()
        try:
            dict = import_points_file(fileName,delimiter = ' ')
        except IOError:
            pass
        else:
            self.failUnless(0 ==1,
                        'bad xya file did not raise error!')    
        os.remove(fileName)
  
    def test_loadxy_bad4(self):
        """
        specifying wrong delimiter
        """
        import os
        import tempfile
       
        fileName = tempfile.mktemp(".xya")
        file = open(fileName,"w")
        file.write("  elevation   speed \n\
1.0 0.0 10.0 0.0\n\
0.0 1.0 0.0 10.0\n\
1.0 0.0 10.4 40.0\n\
#geocrap")
        file.close()
        try:
            dict = import_points_file(fileName,delimiter = ' ')
        except IOError:
            pass
        else:
            self.failUnless(0 ==1,
                        'bad xya file did not raise error!')    
        os.remove(fileName)

    def test_loadxy_bad5(self):
        """
        specifying wrong delimiter
        """
        import os
        import tempfile
       
        fileName = tempfile.mktemp(".xya")
        file = open(fileName,"w")
        file.write("  elevation   speed \n\
1.0 0.0 10.0 0.0\n\
0.0 1.0 0.0 10.0\n\
1.0 0.0 10.4 40.0\n\
#geocrap\n\
crap")
        file.close()
        try:
            dict = import_points_file(fileName,delimiter = ' ')
        except IOError:
            pass
        else:
            self.failUnless(0 ==1,
                        'bad xya file did not raise error!')    
        os.remove(fileName)        
  ############### .PTS ##########
            
    def test_loadpts(self):
        
        from Scientific.IO.NetCDF import NetCDFFile

        fileName = tempfile.mktemp(".pts")
        # NetCDF file definition
        outfile = NetCDFFile(fileName, 'w')
        
        # dimension definitions
        outfile.createDimension('number_of_points', 3)    
        outfile.createDimension('number_of_dimensions', 2) #This is 2d data
    
        # variable definitions
        outfile.createVariable('points', Float, ('number_of_points',
                                                 'number_of_dimensions'))
        outfile.createVariable('elevation', Float, ('number_of_points',))
    
        # Get handles to the variables
        points = outfile.variables['points']
        elevation = outfile.variables['elevation']
 
        points[0, :] = [1.0,0.0]
        elevation[0] = 10.0 
        points[1, :] = [0.0,1.0]
        elevation[1] = 0.0  
        points[2, :] = [1.0,0.0]
        elevation[2] = 10.4    

        outfile.close()
        
        dict = import_points_file(fileName)
        os.remove(fileName)
        answer =  [[1.0, 0.0],[0.0, 1.0],[1.0, 0.0]]
        #print "dict['pointlist']",dict['pointlist']
        #print "answer",answer
	assert allclose(dict['pointlist'], [[1.0, 0.0],[0.0, 1.0],[1.0, 0.0]])
	assert allclose(dict['attributelist']['elevation'], [10.0, 0.0, 10.4])
 
    def test_writepts(self):
        dict = {}
        att_dict = {}
	dict['pointlist'] = array([[1.0, 0.0],[0.0, 1.0],[1.0, 0.0]])
	att_dict['elevation'] = array([10.0, 0.0, 10.4])
	att_dict['brightness'] = array([10.0, 0.0, 10.4])
        dict['attributelist'] = att_dict
        dict['geo_reference'] = Geo_reference(56,1.9,1.9)
        
        
        fileName = tempfile.mktemp(".pts")
        export_points_file(fileName, dict)
        dict2 = import_points_file(fileName)
        #print "fileName",fileName 
        os.remove(fileName)
        #print "dict2",dict2 
        
	assert allclose(dict2['pointlist'],[[1.0, 0.0],[0.0, 1.0],[1.0, 0.0]])
	assert allclose(dict2['attributelist']['elevation'], [10.0, 0.0, 10.4])
        answer = [10.0, 0.0, 10.4]
	assert allclose(dict2['attributelist']['brightness'], answer)

        #print "dict['geo_reference'] ",dict['geo_reference'] 
        #print "dict2['geo_reference']",dict2['geo_reference']
        
        self.failUnless(dict['geo_reference'] == dict2['geo_reference'],
                         'test_writepts failed. Test geo_reference')
        
 ########################## BAD .PTS ##########################          

    def test_load_bad_no_file_pts(self):
        import os
        import tempfile
       
        fileName = tempfile.mktemp(".pts")
        #print fileName
        try:
            dict = import_points_file(fileName)
        except IOError:
            pass
        else:
            self.failUnless(0 ==1,
                        'imaginary file did not raise error!')
            
  ############### .PTS OTHER ##########

    def test_concatinate_attributelist(self):
        dic = {}
        dic['one'] = array([1,2])
        dic['2'] = array([2,7])
        dic['three'] = array([3,79])
        dic['4'] = array([4,47])
        dic['five'] = array([5,17])
        titles, block = concatinate_attributelist(dic)
        #print "titles", titles
        #print "array", block
        self.failUnless(titles  == ['4', '2', 'five', 'three', 'one'],
                         'test_concatinate_attributelist failed.')
	assert allclose(block, [[4,2,5,3,1],[47,7,17,79,2]])

    def test_half_pts(self):
        dict = {}
        att_dict = {}
	dict['pointlist'] = array([[1.0, 0.0],[0.0, 1.0],[1.0, 0.0]])
	att_dict['elevation'] = array([10.0, 0.0, 10.4])
	att_dict['brightness'] = array([10.0, 0.0, 10.4])
        dict['attributelist'] = att_dict
        out_dict = half_pts(dict)
	assert allclose(out_dict['pointlist'],[[1.0, 0.0],[1.0, 0.0]])
	assert allclose(out_dict['attributelist']['elevation'], [10.0, 10.4])
	assert allclose(out_dict['attributelist']['brightness'], [10.0, 10.4])
  
    def test_extent_point_atts(self):
        
        dict = {}
        att_dict = {}
	dict['pointlist'] = array([[1.0, 10.0],[0.0, 1.0],[10.0, -10.0]])
	att_dict['elevation'] = array([30.0, 0.0, 10.4])
	att_dict['brightness'] = array([10.0, 0.0, 10.4])
        dict['attributelist'] = att_dict
        out_dict = extent_point_atts(dict)
        
        #print "out_dict['pointlist']",out_dict #['pointlist'] 
	assert allclose(out_dict['pointlist'],[[0.0, -10.0],[10.0, -10.0],
                                               [10.0,10.0],[0.0, 10.0]])

        self.failUnless(dict['attributelist']  == {},
                         'test_extent_point_atts failed. Test 1')
        
    def test_reduce_pts(self):
        dict = {}
        att_dict = {}
	dict['pointlist'] = array([[1.0, 0.0],[0.0, 1.0],[1.0, 0.0]])
	att_dict['elevation'] = array([10.0, 0.0, 10.4])
	att_dict['brightness'] = array([10.0, 0.0, 10.4])
        dict['attributelist'] = att_dict
        
        inFileName = tempfile.mktemp(".pts")
        export_points_file(inFileName, dict)
        
        outFileName = tempfile.mktemp(".pts")
        
        dict2 = reduce_pts(inFileName,outFileName, 1 )
        os.remove(inFileName)

        dict2 = import_points_file(outFileName)
        os.remove(outFileName)
        #print "dict2",dict2 
        
	assert allclose(dict2['pointlist'],[[1.0, 0.0]])
	assert allclose(dict2['attributelist']['elevation'], [10.0])
	assert allclose(dict2['attributelist']['brightness'], [10.0])
       
    def test_produce_half_point_files(self):
        dict = {}
        att_dict = {}
	dict['pointlist'] = array([[1.0, 0.0],[0.0, 1.0],[1.0, 0.0]])
	att_dict['elevation'] = array([10.0, 0.0, 10.4])
	att_dict['brightness'] = array([10.0, 0.0, 10.4])
        dict['attributelist'] = att_dict
        
        inFileName = tempfile.mktemp(".pts")
        export_points_file(inFileName, dict)
        
        outFileName = tempfile.mktemp(".pts")
        delimiter = '_p'
        outfiles = produce_half_point_files(inFileName, 1,
                                            delimiter, verbose = False )
        os.remove(inFileName)
        root, ext = splitext(inFileName)
        outFileName = root + delimiter + ext
        #print "outFileName",outfiles  
        dict2 = import_points_file(outfiles[1])
        for file in outfiles:
            #print "del file",file 
            os.remove(file)
        
	assert allclose(dict2['pointlist'],[[1.0, 0.0]])
	assert allclose(dict2['attributelist']['elevation'], [10.0])
	assert allclose(dict2['attributelist']['brightness'], [10.0])
                      
#-------------------------------------------------------------
if __name__ == "__main__":

    suite = unittest.makeSuite(loadASCIITestCase,'test')
    #suite = unittest.makeSuite(loadASCIITestCase,'test_writepts')
    runner = unittest.TextTestRunner(verbosity=0)
    runner.run(suite)
    

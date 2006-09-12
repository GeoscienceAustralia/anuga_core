#!/usr/bin/env python

#TEST

import unittest
from math import sqrt


#from least_squares import *
from Numeric import allclose, array, transpose

from anuga.coordinate_transforms.geo_reference import Geo_reference
from combine_pts import *
from load_mesh.loadASCII import import_points_file
#from anuga.geospatial_data.geospatial_data import import_points_file

class Test_combine_pts(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass


    def test_combine_rectangular_points_files(self):
        from load_mesh.loadASCII import export_points_file
        import tempfile
        import os


        # create a fine .pts file
        fine_dict = {}
        fine_dict['pointlist']=[[0,1],[0,4],[4,4],[4,1]]
        att_dict = {}
        att_dict['elevation'] = [10,40,80,50]
        att_dict['sonic'] = [1,4,8,5]
        fine_dict['attributelist'] = att_dict
        fine_dict['geo_reference'] = Geo_reference(56,2.0,1.0)
        
        fine_file = tempfile.mktemp('.pts')
        export_points_file(fine_file, fine_dict)


        # create a coarse .pts file
        coarse_dict = {}
        coarse_dict['pointlist']=[[0,1],[2,1],[0,-2]]
        att_dict = {}
        att_dict['elevation'] = [10,30,-20]
        att_dict['sonic'] = [1,3,-2]
        coarse_dict['attributelist'] = att_dict
        coarse_dict['geo_reference'] = Geo_reference(56,4.0,3.0)
        
        coarse_file = tempfile.mktemp('.pts')
        export_points_file(coarse_file, coarse_dict)
        
        out_file = tempfile.mktemp('.pts')

        combine_rectangular_points_files(fine_file,coarse_file,out_file)

        #clean up
        os.remove(fine_file)
        os.remove(coarse_file)

        results = import_points_file(out_file,
                                     delimiter = ',')
        answer = [[2.0, 0.0],
                  [0.0, 1.0],
                  [0.0, 4.0],
                  [4.0, 4.0],
                  [4.0, 1.0]]
        #print "results",results 
        #print "answer",answer
        
        self.failUnless(len(results['pointlist']) == len(answer),
                        'final number of points wrong. failed.')
        
	assert allclose(results['pointlist'], answer)
	assert allclose(results['attributelist']['sonic'], [ -2.,
                                                             1.,  4.,
                                                             8.,  5.])
	assert allclose(results['attributelist']['elevation'],[ -20.,
                                                               10.,  40.,
                                                               80.,  50.])
        
        self.failUnless(results['geo_reference'] == fine_dict['geo_reference'],
                        ' failed.')
        #clean up
        os.remove(out_file)

    def test_combine_rectangular_points_filesII(self):
        from load_mesh.loadASCII import export_points_file
        import tempfile
        import os

        # create a fine .pts file
        fine_dict = {}
        fine_dict['pointlist']=[[0,1],[0,4],[4,4],[4,1],[3,1],[2,2],[1,3],[3,4]]
        att_dict = {}
        fine_dict['attributelist'] = att_dict
        fine_dict['geo_reference'] = Geo_reference(56,2.0,1.0)
        
        fine_file = tempfile.mktemp('.pts')
        export_points_file(fine_file, fine_dict)


        # create a coarse .pts file
        coarse_dict = {}
        coarse_dict['pointlist']=[[0,1],[0,0],[0.5,0.5],[1,1],
                                  [1.5,1.5],[2,1],[0,-2],[100,10],
                                  [-20,4],[-50,5],[60,70]]
        att_dict = {}
        coarse_dict['attributelist'] = att_dict
        coarse_dict['geo_reference'] = Geo_reference(56,4.0,3.0)
        
        coarse_file = tempfile.mktemp('.pts')
        export_points_file(coarse_file, coarse_dict)
        
        out_file = tempfile.mktemp('.pts')

        combine_rectangular_points_files(fine_file,coarse_file,out_file)

        #clean up
        os.remove(fine_file)
        os.remove(coarse_file)

        results = import_points_file(out_file,
                                  delimiter = ',')
        answer = [[2.0, 0.0],
                  [102.,12.],
                  [-18.,6.],
                  [-48.,7.],
                  [62.,72.],
                  [0.0, 1.0],
                  [0.0, 4.0],
                  [4.0, 4.0],
                  [4.0, 1.0],
                  [3.,1.],
                  [2.,2.],
                  [1.,3.],
                  [3.,4.]]
        #print "results",results['pointlist']
        #print "answer",answer
        #print "len(results['pointlist']",len(results['pointlist'])
        #print "len(answer)",len(answer)
        
        self.failUnless(len(results['pointlist']) == len(answer),
                         'final number of points wrong. failed.')
	assert allclose(results['pointlist'], answer)
        
        self.failUnless(results['geo_reference'] == fine_dict['geo_reference'],
                        ' failed.')
        #clean up
        os.remove(out_file)

    def test_combine_rectangular_points_files_errors(self):
        from load_mesh.loadASCII import export_points_file
        import tempfile
        import os

        # create a fine .pts file
        fine_dict = {}
        fine_dict['pointlist']=[[0,1],[0,4],[4,4],[4,1]]
        att_dict = {}
        att_dict['elevation'] = [1,4,8,5]
        att_dict['pneumonic'] = [1,4,8,5]
        fine_dict['attributelist'] = att_dict
        fine_dict['geo_reference'] = Geo_reference(56,2.0,1.0)
        
        fine_file = tempfile.mktemp('.pts')
        export_points_file(fine_file, fine_dict)


        # create a coarse .pts file
        coarse_dict = {}
        coarse_dict['pointlist']=[[0,1],[2,1],[0,-2]]
        att_dict = {}
        att_dict['elevation'] = [1,3,-2]
        att_dict['sonic'] = [1,3,-2]
        coarse_dict['attributelist'] = att_dict
        coarse_dict['geo_reference'] = Geo_reference(56,4.0,3.0)
        
        coarse_file = tempfile.mktemp('.pts')
        export_points_file(coarse_file, coarse_dict)
        
        out_file = tempfile.mktemp('.pts')
        try:
            combine_rectangular_points_files(fine_file,coarse_file,out_file)
        except AttributeError:
            pass
        else:
            self.failUnless(0 == 1,
                            'bad pts files did not raise error!')
        #clean up
        os.remove(fine_file)
        os.remove(coarse_file)
      
    def test_reduce_points_to_mesh_extent(self):
        from load_mesh.loadASCII import export_points_file, export_mesh_file
        import tempfile
        import os
        x_origin = -45435345.
        y_origin = 433432432.
        # create a fine .pts file
        fine_dict = {}
        fine_dict['pointlist']=[[1.,1.],
                                [1.,7.],
                                [7.,1.],
                                [11.,11.],
                                [7.,8.],
                                [8.,8.],
                                [9.,8.]]
        att_dict = {}
        att_dict['elevation'] = [10,40,80,50,78,78,45]
        att_dict['sonic'] = [1,4,8,5,56,34,213]
        fine_dict['attributelist'] = att_dict
        fine_dict['geo_reference'] = Geo_reference(56,x_origin,y_origin)
        
        points_file = tempfile.mktemp('.pts')
        export_points_file(points_file, fine_dict)


        # create a coarse .pts file
        mesh = {}
        mesh['vertices']=[[0,0],
                          [0,3],
                          [3,3],
                          [3,0],
                          [1,2]
                          ]
        mesh['vertex_attributes']=[[],
                          [],
                          [],
                          [],
                          []
                          ]
        mesh['geo_reference'] = Geo_reference(56,x_origin+7.,y_origin+7.)
        
        mesh_file = tempfile.mktemp('.tsh')
        export_mesh_file(mesh_file, mesh)
        
        out_file = tempfile.mktemp('.pts')

        reduce_points_to_mesh_extent(points_file,mesh_file,out_file)

        results = import_points_file(out_file,
                                  delimiter = ',')
        answer = [
            [7.,8.], 
            [8.,8.],
            [9.,8.]]
        #print "results",results['pointlist']
        #print "answer",answer
        
        self.failUnless(len(results['pointlist']) == len(answer),
                         'final number of points wrong. failed.')
	assert allclose(results['pointlist'],answer)

        answer = [78., 78.,45.]
        
        self.failUnless(len(results['attributelist']['elevation']) == len(answer),
                         'final number of points wrong. failed.')
	assert allclose(results['attributelist']['elevation'], answer)
       
        #clean up
        os.remove(points_file)
        os.remove(mesh_file)
        
        #clean up
        os.remove(out_file)

        #FIXME do test for add points files
        
#-------------------------------------------------------------
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_combine_pts,'test')
    #suite = unittest.makeSuite(Test_combine_pts,'test_reduce_points_to_mesh_extent')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)






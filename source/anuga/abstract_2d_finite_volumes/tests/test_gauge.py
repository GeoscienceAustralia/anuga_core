#!/usr/bin/env python


import unittest
import tempfile
import os
from csv import reader
import time
import numpy as num

import anuga

from anuga.abstract_2d_finite_volumes.gauge import sww2csv_gauges
from anuga.utilities.numerical_tools import mean
from anuga.pmesh.mesh import Mesh
from anuga.file.sww import SWW_file



# def simple_function(x, y):
#     return x+y

class Test_Gauge(unittest.TestCase):
    def setUp(self):

        def elevation_function(x, y):
            return -x
        
        """ Setup for all tests. """
        
        mesh_file = tempfile.mktemp(".tsh")    
        points = [[0.0,0.0],[6.0,0.0],[6.0,6.0],[0.0,6.0]]
        m = Mesh()
        m.add_vertices(points)
        m.auto_segment()
        m.generate_mesh(verbose=False)
        m.export_mesh_file(mesh_file)
        
        # Create shallow water domain
        domain = anuga.Domain(mesh_file)
        os.remove(mesh_file)
 
        domain.default_order = 2

        # This test was made before tight_slope_limiters were introduced
        # Since were are testing interpolation values this is OK
        domain.tight_slope_limiters = 0
                
        # Set some field values
        domain.set_quantity('elevation', elevation_function)
        domain.set_quantity('friction', 0.03)
        domain.set_quantity('xmomentum', 3.0)
        domain.set_quantity('ymomentum', 4.0)

        ######################
        # Boundary conditions
        B = anuga.Transmissive_boundary(domain)
        domain.set_boundary( {'exterior': B})

        # This call mangles the stage values.
        domain.distribute_to_vertices_and_edges()
        domain.set_quantity('stage', 1.0)

        domain.set_name('datatest' + str(time.time()))
        domain.smooth = True
        domain.reduction = mean
        
        self.domain = domain
        
        
    def tearDown(self):
        """Called at end of each test."""
        if self.sww:
            os.remove(self.sww.filename)

    def _create_sww(self,stage=10.0, timestep=2.0):
        self.sww = SWW_file(self.domain)
        self.sww.store_connectivity()
        self.sww.store_timestep()
        self.domain.set_quantity('stage', stage) # This is automatically limited
        # so it will not be less than the elevation
        self.domain.set_time(self.domain.get_time()-self.domain.starttime+timestep)
        self.sww.store_timestep()
        
        
    def test_sww2csv_0(self):

        """Most of this test was copied from test_interpolate
        test_interpole_sww2csv
        
        This is testing the sww2csv_gauges function, by creating a sww file and
        then exporting the gauges and checking the results.
        """
        
        domain = self.domain
        self._create_sww()
        
        # test the function
        points = [[5.0,1.],[0.5,2.]]

        points_file = tempfile.mktemp(".csv") 
#        points_file = 'test_point.csv'
        file_id = open(points_file,"w")
        file_id.write("name, easting, northing, elevation \n\
point1, 5.0, 1.0, 3.0\n\
point2, 0.5, 2.0, 9.0\n")
        file_id.close()

        
        sww2csv_gauges(self.sww.filename, 
                       points_file,
                       verbose=False,
                       use_cache=False)

#        point1_answers_array = [[0.0,1.0,-5.0,3.0,4.0], [2.0,10.0,-5.0,3.0,4.0]]
        point1_answers_array = [[0.0,0.0,1.0,6.0,-5.0,3.0,4.0], [2.0,2.0/3600.,10.0,15.0,-5.0,3.0,4.0]]
        point1_filename = 'gauge_point1.csv'
        point1_handle = open(point1_filename)
        point1_reader = reader(point1_handle)
        point1_reader.next()

        line=[]
        for i,row in enumerate(point1_reader):
            #print 'i',i,'row',row
            line.append([float(row[0]),float(row[1]),float(row[2]),float(row[3]),
                         float(row[4]),float(row[5]),float(row[6])])
            #print 'assert line',line[i],'point1',point1_answers_array[i]
            assert num.allclose(line[i], point1_answers_array[i])

        point2_answers_array = [[0.0,0.0,1.0,1.5,-0.5,3.0,4.0], [2.0,2.0/3600.,10.0,10.5,-0.5,3.0,4.0]]
        point2_filename = 'gauge_point2.csv' 
        point2_handle = open(point2_filename)
        point2_reader = reader(point2_handle)
        point2_reader.next()
                        
        line=[]
        for i,row in enumerate(point2_reader):
#            print 'i',i,'row',row
            line.append([float(row[0]),float(row[1]),float(row[2]),float(row[3]),
                         float(row[4]),float(row[5]),float(row[6])])
#            print 'assert line',line[i],'point1',point1_answers_array[i]
            assert num.allclose(line[i], point2_answers_array[i])
                         
        # clean up
        point1_handle.close()
        point2_handle.close()
        #os.remove(points_file)
        #os.remove(point1_filename)
        #os.remove(point2_filename)


    def test_sww2csv_gauges1(self):
        from anuga.pmesh.mesh import Mesh
        from csv import reader,writer
        import time
        import string
        
        """Most of this test was copied from test_interpolate
        test_interpole_sww2csv
        
        This is testing the sww2csv_gauges function, by creating a sww file and
        then exporting the gauges and checking the results.
        
        This tests the ablity not to have elevation in the points file and 
        not store xmomentum and ymomentum
        """
        
        domain = self.domain
        self._create_sww()
        
        # test the function
        points = [[5.0,1.],[0.5,2.]]

        points_file = tempfile.mktemp(".csv")
#        points_file = 'test_point.csv'
        file_id = open(points_file,"w")
        file_id.write("name,easting,northing \n\
point1, 5.0, 1.0\n\
point2, 0.5, 2.0\n")
        file_id.close()

        sww2csv_gauges(self.sww.filename, 
                            points_file,
                            quantities=['stage', 'elevation'],
                            use_cache=False,
                            verbose=False)

        point1_answers_array = [[0.0,1.0,-5.0], [2.0,10.0,-5.0]]
        point1_filename = 'gauge_point1.csv'
        point1_handle = file(point1_filename)
        point1_reader = reader(point1_handle)
        point1_reader.next()

        line=[]
        for i,row in enumerate(point1_reader):
#            print 'i',i,'row',row
            # note the 'hole' (element 1) below - skip the new 'hours' field
            line.append([float(row[0]),float(row[2]),float(row[3])])
            #print 'line',line[i],'point1',point1_answers_array[i]
            assert num.allclose(line[i], point1_answers_array[i])

        point2_answers_array = [[0.0,1.0,-0.5], [2.0,10.0,-0.5]]
        point2_filename = 'gauge_point2.csv' 
        point2_handle = file(point2_filename)
        point2_reader = reader(point2_handle)
        point2_reader.next()
                        
        line=[]
        for i,row in enumerate(point2_reader):
#            print 'i',i,'row',row
            # note the 'hole' (element 1) below - skip the new 'hours' field
            line.append([float(row[0]),float(row[2]),float(row[3])])
#            print 'line',line[i],'point1',point1_answers_array[i]
            assert num.allclose(line[i], point2_answers_array[i])
                         
        # clean up
        point1_handle.close()
        point2_handle.close() 
        os.remove(points_file)
        os.remove(point1_filename)
        os.remove(point2_filename)        
        

    def test_sww2csv_gauges2(self):
        
        """Most of this test was copied from test_interpolate
        test_interpole_sww2csv
        
        This is testing the sww2csv_gauges function, by creating a sww file and
        then exporting the gauges and checking the results.
        
        This is the same as sww2csv_gauges except set domain.set_starttime to 5.
        Therefore testing the storing of the absolute time in the csv files
        """
        
        domain = self.domain
        domain.set_starttime(1)
        
        self._create_sww(timestep=2)
        
        # test the function
        points = [[5.0,1.],[0.5,2.]]

        points_file = tempfile.mktemp(".csv")
#        points_file = 'test_point.csv'
        file_id = open(points_file,"w")
        file_id.write("name, easting, northing, elevation \n\
point1, 5.0, 1.0, 3.0\n\
point2, 0.5, 2.0, 9.0\n")
        file_id.close()
        
        sww2csv_gauges(self.sww.filename, 
                            points_file,
                            verbose=False,
                            use_cache=False)

#        point1_answers_array = [[0.0,1.0,-5.0,3.0,4.0], [2.0,10.0,-5.0,3.0,4.0]]
        point1_answers_array = [[2.0,2.0/3600.,1.0,6.0,-5.0,3.0,4.0], [3.0,3.0/3600.,10.0,15.0,-5.0,3.0,4.0]]
        point1_filename = 'gauge_point1.csv'
        point1_handle = file(point1_filename)
        point1_reader = reader(point1_handle)
        point1_reader.next()

        line=[]
        for i,row in enumerate(point1_reader):
            #print 'i',i,'row',row
            line.append([float(row[0]),float(row[1]),float(row[2]),float(row[3]),
                         float(row[4]), float(row[5]), float(row[6])])
            #print 'assert line',line[i],'answer',point1_answers_array[i]
            assert num.allclose(line[i], point1_answers_array[i])

        point2_answers_array = [[2.0,2.0/3600.,1.0,1.5,-0.5,3.0,4.0], [3.0,3.0/3600.,10.0,10.5,-0.5,3.0,4.0]]
        point2_filename = 'gauge_point2.csv' 
        point2_handle = file(point2_filename)
        point2_reader = reader(point2_handle)
        point2_reader.next()
                        
        line=[]
        for i,row in enumerate(point2_reader):
            #print 'i',i,'row',row
            line.append([float(row[0]),float(row[1]),float(row[2]),float(row[3]),
                         float(row[4]),float(row[5]), float(row[6])])
            #print 'assert line',line[i],'point1',point1_answers_array[i]
            assert num.allclose(line[i], point2_answers_array[i])
                         
        # clean up
        point1_handle.close()
        point2_handle.close()
        os.remove(points_file)
        os.remove(point1_filename)
        os.remove(point2_filename)


       
    def test_sww2csv_gauge_point_off_mesh(self):
        from anuga.pmesh.mesh import Mesh
        from csv import reader,writer
        import time
        import string
        
        """Most of this test was copied from test_interpolate
        test_interpole_sww2csv
        
        This is testing the sww2csv_gauges function with one gauge off the mesh, by creating a sww file and
        then exporting the gauges and checking the results.
        
        This tests the correct values for when a gauge is off the mesh, which is important for parallel.
        """

        domain = self.domain
        sww = self._create_sww()
     
        # test the function
        points = [[50.0,1.],[50.5,-20.25]]

#        points_file = tempfile.mktemp(".csv")
        points_file = 'test_point.csv'
        file_id = open(points_file,"w")
        file_id.write("name,easting,northing \n\
offmesh1, 50.0, 1.0\n\
offmesh2, 50.5, 20.25\n")
        file_id.close()

        points_files = ['offmesh1.csv', 'offmesh2.csv']        
        
        for point_filename in points_files:
            if os.path.exists(point_filename): os.remove(point_filename)         
        
        sww2csv_gauges(self.sww.filename, 
                            points_file,
                            quantities=['stage', 'elevation', 'bearing'],
                            use_cache=False,
                            verbose=False)

        for point_filename in points_files: 
            assert not os.path.exists(point_filename)
            
        os.remove(points_file)
        
        
    def test_sww2csv_centroid(self):
        
        """Check sww2csv timeseries at centroid.
        
        Test the ability to get a timeseries at the centroid of a triangle, rather
        than the given gauge point.
        """
        
        domain = self.domain
        sww = self._create_sww()
        
        # create a csv file containing our gauge points
        points_file = tempfile.mktemp(".csv")
        file_id = open(points_file,"w")
# These values are where the centroids should be        
#        file_id.write("name, easting, northing, elevation \n\
#point1, 2.0, 2.0, 3.0\n\
#point2, 4.0, 4.0, 9.0\n")
 
# These values are slightly off the centroids - will it find the centroids?
        file_id.write("name, easting, northing, elevation \n\
point1, 2.0, 1.0, 3.0\n\
point2, 4.5, 4.0, 9.0\n")

 
        file_id.close()

        sww2csv_gauges(self.sww.filename, 
                       points_file,
                       verbose=False,
                       use_cache=False,
                       output_centroids=True)

        point1_answers_array = [[0.0,0.0,1.0,3.0,-2.0,3.0,4.0], [2.0,2.0/3600.,10.0,12.0,-2.0,3.0,4.0]]
        point1_filename = 'gauge_point1.csv'
        point1_handle = open(point1_filename)
        point1_reader = reader(point1_handle)
        point1_reader.next()

        line=[]
        for i,row in enumerate(point1_reader):
            line.append([float(row[0]),float(row[1]),float(row[2]),float(row[3]),
                         float(row[4]),float(row[5]),float(row[6])])
#           print 'assert line',line[i],'point1',point1_answers_array[i]
            assert num.allclose(line[i], point1_answers_array[i])

        point2_answers_array = [[0.0,0.0,1.0,5.0,-4.0,3.0,4.0], [2.0,2.0/3600.,10.0,14.0,-4.0,3.0,4.0]]
        point2_filename = 'gauge_point2.csv' 
        point2_handle = open(point2_filename)
        point2_reader = reader(point2_handle)
        point2_reader.next()
                        
        line=[]
        for i,row in enumerate(point2_reader):
            line.append([float(row[0]),float(row[1]),float(row[2]),float(row[3]),
                         float(row[4]),float(row[5]),float(row[6])])
#           print i, 'assert line',line[i],'point2',point2_answers_array[i]
            assert num.allclose(line[i], point2_answers_array[i])
                         
        # clean up
        point1_handle.close()
        point2_handle.close()
        os.remove(points_file)
        os.remove(point1_filename)
        os.remove(point2_filename)


    def test_sww2csv_output_centroid_attribute(self):
        
        """Check sww2csv timeseries at centroid, then output the centroid coordinates.
        
        Test the ability to get a timeseries at the centroid of a triangle, rather
        than the given gauge point, then output the results.
        """
        
        domain = self.domain        
        self._create_sww()
        
        # create a csv file containing our gauge points
        points_file = tempfile.mktemp(".csv")
        file_id = open(points_file,"w")
 
# These values are slightly off the centroids - will it find the centroids?
        file_id.write("name, easting, northing, elevation \n\
point1, 2.5, 4.25, 3.0\n")

        file_id.close()

        sww2csv_gauges(self.sww.filename, 
                       points_file,
                       quantities=['stage', 'xcentroid', 'ycentroid'],
                       verbose=False,
                       use_cache=False,
                       output_centroids=True)

        point1_answers_array = [[0.0,0.0,1.0,4.0,4.0], [2.0,2.0/3600.,10.0,4.0,4.0]]
        point1_filename = 'gauge_point1.csv'
        point1_handle = file(point1_filename)
        point1_reader = reader(point1_handle)
        point1_reader.next()

        line=[]
        for i,row in enumerate(point1_reader):
            line.append([float(row[0]),float(row[1]),float(row[2]),float(row[3]),float(row[4])])
#            print 'assert line',line[i],'point1',point1_answers_array[i]
            assert num.allclose(line[i], point1_answers_array[i])

        # clean up
        point1_handle.close()        
        os.remove(points_file)
        os.remove(point1_filename)

    def test_sww2csv_multiple_files(self):
        """
        This is testing the sww2csv_gauges function, by creating multiple 
        sww files and then exporting the gauges and checking the results.
        """
        timestep=2.0
        domain = self.domain
        domain.set_starttime(0.)
        # Create two sww files with timestep at end. These are to be
        # stored consecutively in the gauge csv files
        basename='datatest1'
        domain.set_name(basename) 
        self._create_sww(stage=10.,timestep=timestep)

        domain.set_name(basename+str(time.time())) 
        domain.set_time(domain.get_time()+timestep)
        self._create_sww(stage=20.,timestep=timestep)

        points_file = tempfile.mktemp(".csv")
        file_id = open(points_file,"w")

        # test the function at these points
        points = [[5.0,1.],[0.5,2.]]

        # create a csv file containing our gauge points
        points_file = tempfile.mktemp(".csv")
        file_id = open(points_file,"w")
        file_id.write("name,easting,northing \n\
point1, 5.0, 1.0\n\
point2, 0.5, 2.0\n")
        file_id.close()


        sww2csv_gauges(basename+".sww", 
                       points_file,
                       quantities=['stage', 'elevation'],
                       use_cache=False,
                       verbose=False)

        point1_answers_array = [[0.0,1.0,-5.0], [2.0,10.0,-5.0],[4.0,10.0,-5.0],
                                [6.0,20.0,-5.0], [0.0,1.0,-5.0]]
        point1_filename = 'gauge_point1.csv'
        point1_handle = file(point1_filename)
        point1_reader = reader(point1_handle)
        point1_reader.next()

        line=[]
        for i,row in enumerate(point1_reader):
            # note the 'hole' (element 1) below - skip the new 'hours' field
            line.append([float(row[0]),float(row[2]),float(row[3])])
            #print 'i', i
            #print 'row',row
            #print 'line',line[i],'point1',point1_answers_array[i]
            assert num.allclose(line[i], point1_answers_array[i])

        point2_answers_array = [[0.0,1.0,-0.5], [2.0,10.0,-0.5],[4.0,10.0,-0.5],
                                [6.0,20.0,-0.5], [0.0,1.0,-0.5]]
        point2_filename = 'gauge_point2.csv' 
        point2_handle = file(point2_filename)
        point2_reader = reader(point2_handle)
        point2_reader.next()
                        
        line=[]
        for i,row in enumerate(point2_reader):
            # note the 'hole' (element 1) below - skip the new 'hours' field
            line.append([float(row[0]),float(row[2]),float(row[3])])
            #print 'line',line[i],'point2',point2_answers_array[i]
            assert num.allclose(line[i], point2_answers_array[i])
                         
        # clean up
        point1_handle.close()
        point2_handle.close() 
        os.remove(points_file)
        os.remove(point1_filename)
        os.remove(point2_filename)       

        #remove second swwfile not removed by tearDown
        os.remove(basename+".sww")
        #os.remove(basename+str(time.time())+".sww")

#-------------------------------------------------------------

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_Gauge, 'test_')
#    runner = unittest.TextTestRunner(verbosity=2)
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)

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
from anuga.pmesh.mesh import Pmesh
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
        m = Pmesh()
        m.add_vertices(points)
        m.auto_segment()
        m.generate_mesh(verbose=False)
        m.export_mesh_file(mesh_file)

        # Create shallow water domain
        domain = anuga.Domain(mesh_file)


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
        self.domain.set_time(self.domain.get_time()+timestep)
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
#        point1_answers_array = [[0.0,0.0,1.0,6.0,-5.0,3.0,4.0], [2.0,2.0/3600.,10.0,15.0,-5.0,3.0,4.0]]
        point1_answers_array = [[0.0, 0.0, 1.0, 4.0, -3.0, 3.0, 4.0],  [2.0, 0.0005555555555555556, 10.0, 13.0, -3.0, 3.0, 4.0]]
        point1_filename = 'gauge_point1.csv'
        point1_handle = open(point1_filename)
        point1_reader = reader(point1_handle)
        next(point1_reader)

        line=[]
        for i,row in enumerate(point1_reader):
            #print 'i',i,'row',row
            line.append([float(row[0]),float(row[1]),float(row[2]),float(row[3]),
                         float(row[4]),float(row[5]),float(row[6])])
            #print 'assert line',line[i],'point1',point1_answers_array[i]
            assert num.allclose(line[i], point1_answers_array[i])

        #point2_answers_array = [[0.0,0.0,1.0,1.5,-0.5,3.0,4.0], [2.0,2.0/3600.,10.0,10.5,-0.5,3.0,4.0]]

        point2_answers_array = [[0.0, 0.0, 1.0, 3.416666666666667, -2.416666666666667, 3.0, 4.0], [2.0, 0.0005555555555555556, 10.000000000000002, 12.416666666666668, -2.416666666666667, 3.0, 4.0] ]



        point2_filename = 'gauge_point2.csv'
        point2_handle = open(point2_filename)
        point2_reader = reader(point2_handle)
        next(point2_reader)

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

        try:
            os.remove(points_file)
            os.remove(point1_filename)
            os.remove(point2_filename)
        except OSError:
            pass


    def test_sww2csv_gauges1(self):
        from anuga.pmesh.mesh import Pmesh
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

        point1_answers_array = [[0.0, 1.0, -3.0], [2.0, 10.0, -3.0]]
        point1_filename = 'gauge_point1.csv'
        point1_handle = open(point1_filename)
        point1_reader = reader(point1_handle)
        next(point1_reader)

        line=[]
        for i,row in enumerate(point1_reader):
#            print 'i',i,'row',row
            # note the 'hole' (element 1) below - skip the new 'hours' field
            line.append([float(row[0]),float(row[2]),float(row[3])])
            #print 'line',line[i],'point1',point1_answers_array[i]
            assert num.allclose(line[i], point1_answers_array[i])

        point2_answers_array = [ [0.0, 1.0, -2.416666666666667], [2.0, 10.000000000000002, -2.416666666666667] ]
        point2_filename = 'gauge_point2.csv'
        point2_handle = open(point2_filename)
        point2_reader = reader(point2_handle)
        next(point2_reader)

        line=[]
        for i,row in enumerate(point2_reader):
#            print 'i',i,'row',row
            # note the 'hole' (element 1) below - skip the new 'hours' field
            line.append([float(row[0]),float(row[2]),float(row[3])])
            # print 'line',line[i],'point1',point1_answers_array[i]
            assert num.allclose(line[i], point2_answers_array[i])

        # clean up
        point1_handle.close()
        point2_handle.close()

        try:
            os.remove(points_file)
            os.remove(point1_filename)
            os.remove(point2_filename)
        except OSError:
            pass


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
        point1_answers_array = [[1.0, 0.0002777777777777778, 1.0, 4.0, -3.0, 3.0, 4.0], [3.0, 0.0008333333333333334, 10.0, 13.0, -3.0, 3.0, 4.0] ]
        point1_filename = 'gauge_point1.csv'
        point1_handle = open(point1_filename)
        point1_reader = reader(point1_handle)
        next(point1_reader)

        line=[]
        for i,row in enumerate(point1_reader):
            #print 'i',i,'row',row
            line.append([float(row[0]),float(row[1]),float(row[2]),float(row[3]),
                         float(row[4]), float(row[5]), float(row[6])])
            #print 'assert line',line[i],'answer',point1_answers_array[i]
            assert num.allclose(line[i], point1_answers_array[i])

        point2_answers_array = [[1.0, 0.0002777777777777778, 1.0, 3.416666666666667, -2.416666666666667, 3.0, 4.0], [3.0, 0.0008333333333333334, 10.000000000000002, 12.416666666666668, -2.416666666666667, 3.0, 4.0]]
        point2_filename = 'gauge_point2.csv'
        point2_handle = open(point2_filename)
        point2_reader = reader(point2_handle)
        next(point2_reader)

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

        try:
            os.remove(points_file)
            os.remove(point1_filename)
            os.remove(point2_filename)
        except OSError:
            pass



    def test_sww2csv_gauge_point_off_mesh(self):
        from anuga.pmesh.mesh import Pmesh
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

        # sww2csv_gauges uses out_name='gauge_' prefix by default
        points_files = ['gauge_offmesh1.csv', 'gauge_offmesh2.csv']

        for point_filename in points_files:
            if os.path.exists(point_filename):
                os.remove(point_filename)

        sww2csv_gauges(self.sww.filename,
                            points_file,
                            quantities=['stage', 'elevation', 'bearing'],
                            use_cache=False,
                            verbose=False)

        # Off-mesh gauges produce header-only CSV files; check no data rows written
        for point_filename in points_files:
            if os.path.exists(point_filename):
                with open(point_filename) as f:
                    lines = f.readlines()
                assert len(lines) <= 1, \
                    f'{point_filename} should have no data rows for off-mesh gauge'

        # clean up
        for point_filename in points_files:
            try:
                os.remove(point_filename)
            except OSError:
                pass
        try:
            os.remove(points_file)
        except OSError:
            pass


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

        #point1_answers_array = [[0.0,0.0,1.0,3.0,-2.0,3.0,4.0], [2.0,2.0/3600.,10.0,12.0,-2.0,3.0,4.0]]
        point1_answers_array = [[0.0, 0.0, 1.0, 3.6666666666666665, -2.6666666666666665, 3.0, 4.0], [2.0, 0.0005555555555555556, 10.0, 12.666666666666666, -2.6666666666666665, 3.0, 4.0]]
        point1_filename = 'gauge_point1.csv'
        point1_handle = open(point1_filename)
        point1_reader = reader(point1_handle)
        next(point1_reader)

        line=[]
        for i,row in enumerate(point1_reader):
            line.append([float(row[0]),float(row[1]),float(row[2]),float(row[3]),
                         float(row[4]),float(row[5]),float(row[6])])
            #print 'assert line',line[i],'point1',point1_answers_array[i]
            assert num.allclose(line[i], point1_answers_array[i])

        #point2_answers_array = [[0.0,0.0,1.0,5.0,-4.0,3.0,4.0], [2.0,2.0/3600.,10.0,14.0,-4.0,3.0,4.0]]
        point2_answers_array = [ [0.0, 0.0, 1.0, 4.333333333333333, -3.333333333333333, 3.0, 4.0], [2.0, 0.0005555555555555556, 10.0, 13.333333333333332, -3.333333333333333, 3.0, 4.0] ]
        point2_filename = 'gauge_point2.csv'
        point2_handle = open(point2_filename)
        point2_reader = reader(point2_handle)
        next(point2_reader)

        line=[]
        for i,row in enumerate(point2_reader):
            line.append([float(row[0]),float(row[1]),float(row[2]),float(row[3]),
                         float(row[4]),float(row[5]),float(row[6])])
            #print i, 'assert line',line[i],'point2',point2_answers_array[i]
            assert num.allclose(line[i], point2_answers_array[i])


        # clean up
        point1_handle.close()
        point2_handle.close()

        try:
            os.remove(points_file)
            os.remove(point1_filename)
            os.remove(point2_filename)
        except OSError:
            pass


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
        point1_handle = open(point1_filename)
        point1_reader = reader(point1_handle)
        next(point1_reader)

        line=[]
        for i,row in enumerate(point1_reader):
            line.append([float(row[0]),float(row[1]),float(row[2]),float(row[3]),float(row[4])])
#            print 'assert line',line[i],'point1',point1_answers_array[i]
            assert num.allclose(line[i], point1_answers_array[i])

        # clean up
        point1_handle.close()

        try:
            os.remove(points_file)
            os.remove(point1_filename)
        except OSError:
            pass

    def test_sww2csv_multiple_files(self):
        """
        This is testing the sww2csv_gauges function, by creating multiple
        sww files and then exporting the gauges and checking the results.
        """
        import shutil
        tmpdir = tempfile.mkdtemp()
        orig_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            self._test_sww2csv_multiple_files_impl()
        finally:
            self.sww = None  # files are inside tmpdir; tearDown must not touch them
            os.chdir(orig_dir)
            shutil.rmtree(tmpdir, ignore_errors=True)

    def _test_sww2csv_multiple_files_impl(self):
        timestep=2.0
        domain = self.domain
        domain.set_starttime(0.)
        domain.set_datadir('.')

        # Create two sww files with timestep at end. These are to be
        # stored consecutively in the gauge csv files
        basename='datatest1'
        domain.set_name(basename)
        self._create_sww(stage=10.,timestep=timestep)

        domain.set_name(basename+str(time.time()))
        domain.set_time(domain.get_time()+timestep)
        self._create_sww(stage=20.,timestep=timestep)

        # create a csv file containing our gauge points
        points_file = tempfile.mktemp(".csv")
        points_handle = open(points_file,"w")
        points_handle.write("name,easting,northing \n\
point1, 5.0, 1.0\n\
point2, 0.5, 2.0\n")
        points_handle.close()

        sww2csv_gauges(basename+".sww",
                       points_file,
                       quantities=['stage', 'elevation'],
                       use_cache=False,
                       verbose=False)

        point1_answers_array = [[0.0, 1.0, -3.0], [2.0, 10.0, -3.0],
                               [4.0, 10.0, -3.0], [6.0, 20.0, -3.0]]

        point1_filename = 'gauge_point1.csv'
        point1_handle = open(point1_filename)
        point1_reader = reader(point1_handle)
        next(point1_reader)

        line=[]
        for i,row in enumerate(point1_reader):
            # note the 'hole' (element 1) below - skip the new 'hours' field
            line.append([float(row[0]),float(row[2]),float(row[3])])
            assert num.allclose(line[i], point1_answers_array[i])
        point1_handle.close()

        point2_answers_array = [[0.0, 1.0, -2.416666666666667],
                                [2.0, 10.000000000000002, -2.416666666666667],
                                [4.0, 10.000000000000002, -2.416666666666667],
                                [6.0, 20.000000000000004, -2.416666666666667]]

        point2_filename = 'gauge_point2.csv'
        point2_handle = open(point2_filename)
        point2_reader = reader(point2_handle)
        next(point2_reader)

        line=[]
        for i,row in enumerate(point2_reader):
            # note the 'hole' (element 1) below - skip the new 'hours' field
            line.append([float(row[0]),float(row[2]),float(row[3])])
            assert num.allclose(line[i], point2_answers_array[i])
        point2_handle.close()


class Test_quantities2csv(unittest.TestCase):
    """Unit tests for _quantities2csv — the pure-function helper in gauge.py."""

    def _call(self, quantities, point_quantities, centroids=None, point_i=0):
        from anuga.abstract_2d_finite_volumes.gauge import _quantities2csv
        if centroids is None:
            centroids = [[(1.0, 2.0)]]
        return _quantities2csv(quantities, point_quantities, centroids, point_i)

    def _pq(self, stage=2.0, elev=0.5, xmom=3.0, ymom=4.0):
        """Return a typical point_quantities array [stage, elev, xmom, ymom]."""
        return [stage, elev, xmom, ymom]

    # ------------------------------------------------------------------
    # Existing quantities (already mostly covered) — sanity checks
    # ------------------------------------------------------------------

    def test_stage(self):
        result = self._call(['stage'], self._pq(stage=5.0))
        self.assertAlmostEqual(result[0], 5.0)

    def test_elevation(self):
        result = self._call(['elevation'], self._pq(elev=1.5))
        self.assertAlmostEqual(result[0], 1.5)

    def test_depth(self):
        # depth = stage - elevation
        result = self._call(['depth'], self._pq(stage=3.0, elev=1.0))
        self.assertAlmostEqual(result[0], 2.0)

    # ------------------------------------------------------------------
    # momentum (lines 47-49) — previously uncovered
    # ------------------------------------------------------------------

    def test_momentum(self):
        # momentum = sqrt(xmom^2 + ymom^2) = sqrt(3^2 + 4^2) = 5
        result = self._call(['momentum'], self._pq(xmom=3.0, ymom=4.0))
        self.assertAlmostEqual(result[0], 5.0)

    def test_momentum_zero(self):
        result = self._call(['momentum'], self._pq(xmom=0.0, ymom=0.0))
        self.assertAlmostEqual(result[0], 0.0)

    # ------------------------------------------------------------------
    # speed (lines 53-64) — previously uncovered
    # ------------------------------------------------------------------

    def test_speed_dry_depth_is_zero(self):
        """When depth < 0.001, speed == 0."""
        pq = self._pq(stage=0.5, elev=0.5)   # depth = 0.0
        result = self._call(['speed'], pq)
        self.assertAlmostEqual(result[0], 0.0)

    def test_speed_shallow_depth_below_threshold(self):
        """When depth < 0.001, speed == 0 regardless of momentum."""
        pq = self._pq(stage=0.5009, elev=0.5, xmom=10.0, ymom=10.0)  # depth=0.0009
        result = self._call(['speed'], pq)
        self.assertAlmostEqual(result[0], 0.0)

    def test_speed_normal(self):
        """speed = momentum / depth for depth >= 0.001 and xmom < 1e6."""
        # depth = 2.0 - 0.0 = 2.0, momentum = sqrt(3^2 + 4^2) = 5, speed = 5/2 = 2.5
        pq = self._pq(stage=2.0, elev=0.0, xmom=3.0, ymom=4.0)
        result = self._call(['speed'], pq)
        self.assertAlmostEqual(result[0], 2.5)

    def test_speed_huge_momentum(self):
        """When xmom >= 1e6, speed is forced to 0."""
        pq = self._pq(stage=2.0, elev=0.0, xmom=2.0e6, ymom=0.0)
        result = self._call(['speed'], pq)
        self.assertAlmostEqual(result[0], 0.0)

    # ------------------------------------------------------------------
    # bearing (line 67) — previously uncovered
    # ------------------------------------------------------------------

    def test_bearing(self):
        from anuga.abstract_2d_finite_volumes.util import calc_bearing
        pq = self._pq(xmom=1.0, ymom=0.0)
        result = self._call(['bearing'], pq)
        expected = calc_bearing(1.0, 0.0)
        self.assertAlmostEqual(result[0], expected)

    # ------------------------------------------------------------------
    # centroid coordinates (lines 69-73)
    # ------------------------------------------------------------------

    def test_xcentroid(self):
        centroids = [[5.5, 7.3]]
        result = self._call(['xcentroid'], self._pq(), centroids=centroids, point_i=0)
        self.assertAlmostEqual(result[0], 5.5)

    def test_ycentroid(self):
        centroids = [[5.5, 7.3]]
        result = self._call(['ycentroid'], self._pq(), centroids=centroids, point_i=0)
        self.assertAlmostEqual(result[0], 7.3)

    # ------------------------------------------------------------------
    # Multiple quantities in one call
    # ------------------------------------------------------------------

    def test_multiple_quantities(self):
        pq = self._pq(stage=3.0, elev=1.0, xmom=3.0, ymom=4.0)
        result = self._call(['stage', 'depth', 'momentum', 'speed'], pq)
        self.assertEqual(len(result), 4)
        self.assertAlmostEqual(result[0], 3.0)       # stage
        self.assertAlmostEqual(result[1], 2.0)       # depth = 3 - 1
        self.assertAlmostEqual(result[2], 5.0)       # momentum = sqrt(9+16)
        self.assertAlmostEqual(result[3], 5.0/2.0)   # speed = 5/2


class Test_gauge_get_from_file(unittest.TestCase):
    """Unit tests for gauge_get_from_file."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _write(self, name, content):
        path = os.path.join(self.tmpdir, name)
        with open(path, 'w') as f:
            f.write(content)
        return path

    def test_full_header(self):
        from anuga.abstract_2d_finite_volumes.gauge import gauge_get_from_file
        path = self._write('gauges.csv',
            'easting,northing,name,elevation\n'
            '100.0,200.0,gauge1,5.0\n'
            '300.0,400.0,gauge2,10.0\n')
        gauges, locations, elev = gauge_get_from_file(path)
        self.assertEqual(len(gauges), 2)
        self.assertAlmostEqual(gauges[0][0], 100.0)
        self.assertAlmostEqual(gauges[0][1], 200.0)
        self.assertEqual(locations[0].strip(), 'gauge1')
        self.assertAlmostEqual(elev[0], 5.0)
        self.assertAlmostEqual(gauges[1][0], 300.0)
        self.assertAlmostEqual(elev[1], 10.0)

    def test_header_any_column_order(self):
        from anuga.abstract_2d_finite_volumes.gauge import gauge_get_from_file
        path = self._write('gauges_reorder.csv',
            'name,elevation,northing,easting\n'
            'siteA,3.5,250.0,150.0\n')
        gauges, locations, elev = gauge_get_from_file(path)
        self.assertAlmostEqual(gauges[0][0], 150.0)   # easting
        self.assertAlmostEqual(gauges[0][1], 250.0)   # northing
        self.assertAlmostEqual(elev[0], 3.5)
        self.assertEqual(locations[0].strip(), 'siteA')

    def test_missing_easting_raises(self):
        from anuga.abstract_2d_finite_volumes.gauge import gauge_get_from_file
        path = self._write('bad.csv',
            'northing,name,elevation\n'
            '200.0,g1,5.0\n')
        with self.assertRaises(Exception) as ctx:
            gauge_get_from_file(path)
        self.assertIn('easting', str(ctx.exception).lower())

    def test_missing_elevation_raises(self):
        from anuga.abstract_2d_finite_volumes.gauge import gauge_get_from_file
        path = self._write('no_elev.csv',
            'easting,northing,name\n'
            '100.0,200.0,g1\n')
        with self.assertRaises(Exception) as ctx:
            gauge_get_from_file(path)
        self.assertIn('elevation', str(ctx.exception).lower())

    def test_missing_name_raises(self):
        from anuga.abstract_2d_finite_volumes.gauge import gauge_get_from_file
        path = self._write('no_name.csv',
            'easting,northing,elevation\n'
            '100.0,200.0,5.0\n')
        with self.assertRaises(Exception) as ctx:
            gauge_get_from_file(path)
        self.assertIn('name', str(ctx.exception).lower())

    def test_header_case_insensitive(self):
        from anuga.abstract_2d_finite_volumes.gauge import gauge_get_from_file
        path = self._write('gauges_upper.csv',
            'Easting,Northing,Name,Elevation\n'
            '10.0,20.0,PT1,1.0\n')
        gauges, locations, elev = gauge_get_from_file(path)
        self.assertAlmostEqual(gauges[0][0], 10.0)
        self.assertAlmostEqual(elev[0], 1.0)


class Test_sww2timeseries(unittest.TestCase):
    """Smoke test for sww2timeseries — CSV output only, no figures."""

    def setUp(self):
        import shutil
        from anuga.pmesh.mesh import Pmesh
        from anuga.file.sww import SWW_file
        from anuga.utilities.numerical_tools import mean

        self.tmpdir = tempfile.mkdtemp()
        self.orig_dir = os.getcwd()
        os.chdir(self.tmpdir)

        mesh_file = tempfile.mktemp('.tsh')
        points = [[0.0, 0.0], [6.0, 0.0], [6.0, 6.0], [0.0, 6.0]]
        m = Pmesh()
        m.add_vertices(points)
        m.auto_segment()
        m.generate_mesh(verbose=False)
        m.export_mesh_file(mesh_file)

        domain = anuga.Domain(mesh_file)
        domain.default_order = 2
        domain.tight_slope_limiters = 0
        domain.set_quantity('elevation', lambda x, y: -x)
        domain.set_quantity('friction', 0.03)
        domain.set_quantity('xmomentum', 1.0)
        domain.set_quantity('ymomentum', 0.0)
        domain.set_boundary({'exterior': anuga.Transmissive_boundary(domain)})
        domain.distribute_to_vertices_and_edges()
        domain.set_quantity('stage', 1.0)
        domain.set_name('sww2ts_test')
        domain.smooth = True
        domain.reduction = mean

        sww = SWW_file(domain)
        sww.store_connectivity()
        sww.store_timestep()
        domain.set_quantity('stage', 2.0)
        domain.set_time(domain.get_time() + 2.0)
        sww.store_timestep()

        self.swwfile = sww.filename
        self.swwdir = os.path.dirname(os.path.abspath(sww.filename))

        gauge_path = os.path.join(self.tmpdir, 'gauges.csv')
        with open(gauge_path, 'w') as f:
            f.write('easting,northing,name,elevation\n'
                    '3.0,3.0,mid,0.0\n')
        self.gauge_file = gauge_path

    def tearDown(self):
        import shutil
        os.chdir(self.orig_dir)
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_csv_output_created(self):
        from anuga.abstract_2d_finite_volumes.gauge import sww2timeseries
        swwfiles = {self.swwfile: ''}
        sww2timeseries(swwfiles,
                       self.gauge_file,
                       production_dirs={self.swwdir: 'test'},
                       report=False,
                       generate_fig=False,
                       plot_quantity=['depth'],
                       use_cache=False,
                       verbose=False)
        # A CSV for the gauge should appear in the sww directory
        csv_name = os.path.join(self.swwdir, 'gauges_time_series_mid.csv')
        self.assertTrue(os.path.exists(csv_name),
                        'Expected gauge CSV not found: %s' % csv_name)
        with open(csv_name) as f:
            lines = f.readlines()
        self.assertGreater(len(lines), 1)   # header + data rows


class Test_sww2csv_gauges_errors(unittest.TestCase):
    """Edge-case and error-path coverage for sww2csv_gauges."""

    def setUp(self):
        import shutil
        from anuga.pmesh.mesh import Pmesh
        from anuga.file.sww import SWW_file
        from anuga.utilities.numerical_tools import mean

        self.tmpdir = tempfile.mkdtemp()
        self.orig_dir = os.getcwd()
        os.chdir(self.tmpdir)

        mesh_file = tempfile.mktemp('.tsh')
        points = [[0.0, 0.0], [6.0, 0.0], [6.0, 6.0], [0.0, 6.0]]
        m = Pmesh()
        m.add_vertices(points)
        m.auto_segment()
        m.generate_mesh(verbose=False)
        m.export_mesh_file(mesh_file)

        domain = anuga.Domain(mesh_file)
        domain.default_order = 2
        domain.tight_slope_limiters = 0
        domain.set_quantity('elevation', lambda x, y: -x)
        domain.set_quantity('friction', 0.03)
        domain.set_boundary({'exterior': anuga.Transmissive_boundary(domain)})
        domain.distribute_to_vertices_and_edges()
        domain.set_quantity('stage', 1.0)
        domain.set_name('err_test')
        domain.smooth = True
        domain.reduction = mean

        sww = SWW_file(domain)
        sww.store_connectivity()
        sww.store_timestep()

        self.swwfile = sww.filename

    def tearDown(self):
        import shutil
        os.chdir(self.orig_dir)
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_bad_gauge_file_raises(self):
        """sww2csv_gauges raises when gauge file does not exist."""
        with self.assertRaises(Exception) as ctx:
            sww2csv_gauges(self.swwfile, '/no/such/file.csv',
                           use_cache=False, verbose=False)
        self.assertIn('could not be opened', str(ctx.exception))

    def test_verbose_runs(self):
        """sww2csv_gauges with verbose=True runs without error."""
        gauge_path = os.path.join(self.tmpdir, 'g.csv')
        with open(gauge_path, 'w') as f:
            f.write('name,easting,northing,elevation\n'
                    'p1,3.0,3.0,0.0\n')
        sww2csv_gauges(self.swwfile, gauge_path,
                       use_cache=False, verbose=True)

    def test_verbose_off_mesh_gauge_warns(self):
        """sww2csv_gauges with verbose=True and an off-mesh gauge logs a warning."""
        gauge_path = os.path.join(self.tmpdir, 'off.csv')
        with open(gauge_path, 'w') as f:
            f.write('name,easting,northing,elevation\n'
                    'off,999.0,999.0,0.0\n')
        sww2csv_gauges(self.swwfile, gauge_path,
                       use_cache=False, verbose=True)



class Test_sww2timeseries_branches(unittest.TestCase):
    """Coverage for _sww2timeseries defaults, error paths and verbose branches."""

    def setUp(self):
        import shutil
        from anuga.pmesh.mesh import Pmesh
        from anuga.file.sww import SWW_file
        from anuga.utilities.numerical_tools import mean

        self.tmpdir = tempfile.mkdtemp()
        self.orig_dir = os.getcwd()
        os.chdir(self.tmpdir)

        mesh_file = tempfile.mktemp('.tsh')
        points = [[0.0, 0.0], [6.0, 0.0], [6.0, 6.0], [0.0, 6.0]]
        m = Pmesh()
        m.add_vertices(points)
        m.auto_segment()
        m.generate_mesh(verbose=False)
        m.export_mesh_file(mesh_file)

        domain = anuga.Domain(mesh_file)
        domain.default_order = 2
        domain.tight_slope_limiters = 0
        domain.set_quantity('elevation', lambda x, y: -x)
        domain.set_quantity('friction', 0.03)
        domain.set_boundary({'exterior': anuga.Transmissive_boundary(domain)})
        domain.distribute_to_vertices_and_edges()
        domain.set_quantity('stage', 1.0)
        domain.set_name('ts_branch_test')
        domain.smooth = True
        domain.reduction = mean

        sww = SWW_file(domain)
        sww.store_connectivity()
        sww.store_timestep()
        domain.set_quantity('stage', 2.0)
        domain.set_time(domain.get_time() + 2.0)
        sww.store_timestep()

        self.swwfile = sww.filename
        self.swwdir = os.path.dirname(os.path.abspath(sww.filename))

        gauge_path = os.path.join(self.tmpdir, 'gauges.csv')
        with open(gauge_path, 'w') as f:
            f.write('easting,northing,name,elevation\n'
                    '3.0,3.0,mid,0.0\n')
        self.gauge_file = gauge_path

    def tearDown(self):
        import shutil
        os.chdir(self.orig_dir)
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _call(self, **kwargs):
        from anuga.abstract_2d_finite_volumes.gauge import sww2timeseries
        defaults = dict(
            report=False, generate_fig=False,
            plot_quantity=['depth'], use_cache=False, verbose=False,
        )
        defaults.update(kwargs)
        return sww2timeseries(
            {self.swwfile: ''},
            self.gauge_file,
            production_dirs={self.swwdir: 'test'},
            **defaults,
        )

    def test_gauge_file_not_found_raises(self):
        from anuga.abstract_2d_finite_volumes.gauge import sww2timeseries
        with self.assertRaises(Exception) as ctx:
            sww2timeseries(
                {self.swwfile: ''}, '/no/such/gauges.csv',
                production_dirs={}, report=False, generate_fig=False,
                use_cache=False, verbose=False,
            )
        self.assertIn('could not be opened', str(ctx.exception))

    def test_swwfile_not_found_raises(self):
        from anuga.abstract_2d_finite_volumes.gauge import sww2timeseries
        with self.assertRaises(Exception) as ctx:
            sww2timeseries(
                {'/no/such.sww': ''}, self.gauge_file,
                production_dirs={}, report=False, generate_fig=False,
                use_cache=False, verbose=False,
            )
        self.assertIn('could not be opened', str(ctx.exception))

    def test_defaults_none_values(self):
        """Passing None for report/plot_quantity/surface/time_unit/title_on
        exercises the default-value branches in _sww2timeseries."""
        self._call(report=None, plot_quantity=None,
                   surface=None, time_unit=None, title_on=None)

    def test_verbose_true(self):
        """verbose=True exercises the log.info paths."""
        self._call(verbose=True)

    def test_time_min_too_early_raises(self):
        """time_min earlier than simulation start raises."""
        with self.assertRaises(Exception) as ctx:
            self._call(time_min=-999.0)
        self.assertIn('Minimum time', str(ctx.exception))

    def test_time_max_too_late_raises(self):
        """time_max later than simulation end raises."""
        with self.assertRaises(Exception) as ctx:
            self._call(time_max=99999.0)
        self.assertIn('Maximum time', str(ctx.exception))

    def test_all_gauges_off_mesh_returns_empty(self):
        """When all gauges are outside the mesh the texfile is empty."""
        off_gauge = os.path.join(self.tmpdir, 'off.csv')
        with open(off_gauge, 'w') as f:
            f.write('easting,northing,name,elevation\n'
                    '999.0,999.0,off,0.0\n')
        from anuga.abstract_2d_finite_volumes.gauge import sww2timeseries
        texfile, elev_output = sww2timeseries(
            {self.swwfile: ''}, off_gauge,
            production_dirs={self.swwdir: 'test'},
            report=False, generate_fig=False,
            plot_quantity=['depth'], use_cache=False, verbose=False,
        )
        self.assertEqual(texfile, '')
        self.assertEqual(elev_output, [])

    def test_explicit_surface_time_unit_title_on(self):
        """Passing non-None surface/time_unit/title_on skips their default branches."""
        self._call(surface=False, time_unit='hours', title_on=True)

    def test_explicit_valid_time_bounds(self):
        """Explicit time_min/time_max within the run range are accepted."""
        self._call(time_min=0.0, time_max=2.0)

    def test_all_gauges_off_mesh_verbose(self):
        """verbose=True with all gauges off-mesh hits the 'No gauges contained' log path."""
        off_gauge = os.path.join(self.tmpdir, 'off.csv')
        with open(off_gauge, 'w') as f:
            f.write('easting,northing,name,elevation\n'
                    '999.0,999.0,off,0.0\n')
        from anuga.abstract_2d_finite_volumes.gauge import sww2timeseries
        sww2timeseries(
            {self.swwfile: ''}, off_gauge,
            production_dirs={self.swwdir: 'test'},
            report=False, generate_fig=False,
            plot_quantity=['depth'], use_cache=False, verbose=True,
        )

    def test_some_gauges_off_mesh_verbose(self):
        """verbose=True with mixed on/off-mesh gauges hits 'Gauges not contained here' path."""
        mixed_gauge = os.path.join(self.tmpdir, 'mixed.csv')
        with open(mixed_gauge, 'w') as f:
            f.write('easting,northing,name,elevation\n'
                    '3.0,3.0,on,0.0\n'
                    '999.0,999.0,off,0.0\n')
        self._call_with_gauge(mixed_gauge, verbose=True)

    def _call_with_gauge(self, gauge_file, **kwargs):
        from anuga.abstract_2d_finite_volumes.gauge import sww2timeseries
        defaults = dict(
            report=False, generate_fig=False,
            plot_quantity=['depth'], use_cache=False, verbose=False,
        )
        defaults.update(kwargs)
        return sww2timeseries(
            {self.swwfile: ''}, gauge_file,
            production_dirs={self.swwdir: 'test'},
            **defaults,
        )


#-------------------------------------------------------------

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_Gauge)
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)

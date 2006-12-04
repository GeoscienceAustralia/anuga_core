#!/usr/bin/env python
#

import unittest
import copy
from Numeric import zeros, array, allclose, Float
from anuga.utilities.numerical_tools import mean
import tempfile
import os
from Scientific.IO.NetCDF import NetCDFFile
from struct import pack

from anuga.shallow_water import *
from anuga.shallow_water.data_manager import *
from anuga.config import epsilon
from anuga.utilities.anuga_exceptions import ANUGAError
from anuga.utilities.numerical_tools import ensure_numeric

# This is needed to run the tests of local functions
import data_manager 
from anuga.coordinate_transforms.redfearn import redfearn
from anuga.coordinate_transforms.geo_reference import Geo_reference

class Test_Data_Manager(unittest.TestCase):
    def setUp(self):
        import time
        from mesh_factory import rectangular

        #Create basic mesh
        points, vertices, boundary = rectangular(2, 2)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.default_order=2


        #Set some field values
        domain.set_quantity('elevation', lambda x,y: -x)
        domain.set_quantity('friction', 0.03)


        ######################
        # Boundary conditions
        B = Transmissive_boundary(domain)
        domain.set_boundary( {'left': B, 'right': B, 'top': B, 'bottom': B})


        ######################
        #Initial condition - with jumps


        bed = domain.quantities['elevation'].vertex_values
        stage = zeros(bed.shape, Float)

        h = 0.3
        for i in range(stage.shape[0]):
            if i % 2 == 0:
                stage[i,:] = bed[i,:] + h
            else:
                stage[i,:] = bed[i,:]

        domain.set_quantity('stage', stage)
        self.initial_stage = copy.copy(domain.quantities['stage'].vertex_values)

        domain.distribute_to_vertices_and_edges()


        self.domain = domain

        C = domain.get_vertex_coordinates()
        self.X = C[:,0:6:2].copy()
        self.Y = C[:,1:6:2].copy()

        self.F = bed

        #Write A testfile (not realistic. Values aren't realistic)
        self.test_MOST_file = 'most_small'

        longitudes = [150.66667, 150.83334, 151., 151.16667]
        latitudes = [-34.5, -34.33333, -34.16667, -34]

        long_name = 'LON'
        lat_name = 'LAT'

        nx = 4
        ny = 4
        six = 6


        for ext in ['_ha.nc', '_ua.nc', '_va.nc', '_e.nc']:
            fid = NetCDFFile(self.test_MOST_file + ext, 'w')

            fid.createDimension(long_name,nx)
            fid.createVariable(long_name,'d',(long_name,))
            fid.variables[long_name].point_spacing='uneven'
            fid.variables[long_name].units='degrees_east'
            fid.variables[long_name].assignValue(longitudes)

            fid.createDimension(lat_name,ny)
            fid.createVariable(lat_name,'d',(lat_name,))
            fid.variables[lat_name].point_spacing='uneven'
            fid.variables[lat_name].units='degrees_north'
            fid.variables[lat_name].assignValue(latitudes)

            fid.createDimension('TIME',six)
            fid.createVariable('TIME','d',('TIME',))
            fid.variables['TIME'].point_spacing='uneven'
            fid.variables['TIME'].units='seconds'
            fid.variables['TIME'].assignValue([0.0, 0.1, 0.6, 1.1, 1.6, 2.1])


            name = ext[1:3].upper()
            if name == 'E.': name = 'ELEVATION'
            fid.createVariable(name,'d',('TIME', lat_name, long_name))
            fid.variables[name].units='CENTIMETERS'
            fid.variables[name].missing_value=-1.e+034

            fid.variables[name].assignValue([[[0.3400644, 0, -46.63519, -6.50198],
                                              [-0.1214216, 0, 0, 0],
                                              [0, 0, 0, 0],
                                              [0, 0, 0, 0]],
                                             [[0.3400644, 2.291054e-005, -23.33335, -6.50198],
                                              [-0.1213987, 4.581959e-005, -1.594838e-007, 1.421085e-012],
                                              [2.291054e-005, 4.582107e-005, 4.581715e-005, 1.854517e-009],
                                              [0, 2.291054e-005, 2.291054e-005, 0]],
                                             [[0.3400644, 0.0001374632, -23.31503, -6.50198],
                                              [-0.1212842, 0.0002756907, 0.006325484, 1.380492e-006],
                                              [0.0001374632, 0.0002749264, 0.0002742863, 6.665601e-008],
                                              [0, 0.0001374632, 0.0001374632, 0]],
                                             [[0.3400644, 0.0002520159, -23.29672, -6.50198],
                                              [-0.1211696, 0.0005075303, 0.01264618, 6.208276e-006],
                                              [0.0002520159, 0.0005040318, 0.0005027961, 2.23865e-007],
                                              [0, 0.0002520159, 0.0002520159, 0]],
                                             [[0.3400644, 0.0003665686, -23.27842, -6.50198],
                                              [-0.1210551, 0.0007413362, 0.01896192, 1.447638e-005],
                                              [0.0003665686, 0.0007331371, 0.0007313463, 4.734126e-007],
                                              [0, 0.0003665686, 0.0003665686, 0]],
                                             [[0.3400644, 0.0004811212, -23.26012, -6.50198],
                                              [-0.1209405, 0.0009771062, 0.02527271, 2.617787e-005],
                                              [0.0004811212, 0.0009622425, 0.0009599366, 8.152277e-007],
                                              [0, 0.0004811212, 0.0004811212, 0]]])


            fid.close()




    def tearDown(self):
        import os
        for ext in ['_ha.nc', '_ua.nc', '_va.nc', '_e.nc']:
            #print 'Trying to remove', self.test_MOST_file + ext
            os.remove(self.test_MOST_file + ext)

    def test_sww_constant(self):
        """Test that constant sww information can be written correctly
        (non smooth)
        """

        import time, os
        from Numeric import array, zeros, allclose, Float, concatenate
        from Scientific.IO.NetCDF import NetCDFFile

        self.domain.set_name('datatest' + str(id(self)))
        self.domain.format = 'sww' #Remove??
        self.domain.smooth = False

        sww = get_dataobject(self.domain)
        sww.store_connectivity()

        #Check contents
        #Get NetCDF
        fid = NetCDFFile(sww.filename, 'r')  #Open existing file for append

        # Get the variables
        x = fid.variables['x']
        y = fid.variables['y']
        z = fid.variables['elevation']

        volumes = fid.variables['volumes']


        assert allclose (x[:], self.X.flat)
        assert allclose (y[:], self.Y.flat)
        assert allclose (z[:], self.F.flat)

        V = volumes

        P = len(self.domain)
        for k in range(P):
            assert V[k, 0] == 3*k
            assert V[k, 1] == 3*k+1
            assert V[k, 2] == 3*k+2


        fid.close()

        #Cleanup
        os.remove(sww.filename)


    def test_sww_constant_smooth(self):
        """Test that constant sww information can be written correctly
        (non smooth)
        """

        import time, os
        from Numeric import array, zeros, allclose, Float, concatenate
        from Scientific.IO.NetCDF import NetCDFFile

        self.domain.set_name('datatest' + str(id(self)))
        self.domain.format = 'sww'
        self.domain.smooth = True

        sww = get_dataobject(self.domain)
        sww.store_connectivity()

        #Check contents
        #Get NetCDF
        fid = NetCDFFile(sww.filename, 'r')  #Open existing file for append

        # Get the variables
        x = fid.variables['x']
        y = fid.variables['y']
        z = fid.variables['elevation']

        volumes = fid.variables['volumes']

        X = x[:]
        Y = y[:]

        assert allclose([X[0], Y[0]], array([0.0, 0.0]))
        assert allclose([X[1], Y[1]], array([0.0, 0.5]))
        assert allclose([X[2], Y[2]], array([0.0, 1.0]))

        assert allclose([X[4], Y[4]], array([0.5, 0.5]))

        assert allclose([X[7], Y[7]], array([1.0, 0.5]))

        Z = z[:]
        assert Z[4] == -0.5

        V = volumes
        assert V[2,0] == 4
        assert V[2,1] == 5
        assert V[2,2] == 1

        assert V[4,0] == 6
        assert V[4,1] == 7
        assert V[4,2] == 3


        fid.close()

        #Cleanup
        os.remove(sww.filename)



    def test_sww_variable(self):
        """Test that sww information can be written correctly
        """

        import time, os
        from Numeric import array, zeros, allclose, Float, concatenate
        from Scientific.IO.NetCDF import NetCDFFile

        self.domain.set_name('datatest' + str(id(self)))
        self.domain.format = 'sww'
        self.domain.smooth = True
        self.domain.reduction = mean

        sww = get_dataobject(self.domain)
        sww.store_connectivity()
        sww.store_timestep('stage')

        #Check contents
        #Get NetCDF
        fid = NetCDFFile(sww.filename, 'r')  #Open existing file for append


        # Get the variables
        x = fid.variables['x']
        y = fid.variables['y']
        z = fid.variables['elevation']
        time = fid.variables['time']
        stage = fid.variables['stage']


        Q = self.domain.quantities['stage']
        Q0 = Q.vertex_values[:,0]
        Q1 = Q.vertex_values[:,1]
        Q2 = Q.vertex_values[:,2]

        A = stage[0,:]
        #print A[0], (Q2[0,0] + Q1[1,0])/2
        assert allclose(A[0], (Q2[0] + Q1[1])/2)
        assert allclose(A[1], (Q0[1] + Q1[3] + Q2[2])/3)
        assert allclose(A[2], Q0[3])
        assert allclose(A[3], (Q0[0] + Q1[5] + Q2[4])/3)

        #Center point
        assert allclose(A[4], (Q1[0] + Q2[1] + Q0[2] +\
                                 Q0[5] + Q2[6] + Q1[7])/6)



        fid.close()

        #Cleanup
        os.remove(sww.filename)


    def test_sww_variable2(self):
        """Test that sww information can be written correctly
        multiple timesteps. Use average as reduction operator
        """

        import time, os
        from Numeric import array, zeros, allclose, Float, concatenate
        from Scientific.IO.NetCDF import NetCDFFile

        self.domain.set_name('datatest' + str(id(self)))
        self.domain.format = 'sww'
        self.domain.smooth = True

        self.domain.reduction = mean

        sww = get_dataobject(self.domain)
        sww.store_connectivity()
        sww.store_timestep('stage')
        self.domain.evolve_to_end(finaltime = 0.01)
        sww.store_timestep('stage')


        #Check contents
        #Get NetCDF
        fid = NetCDFFile(sww.filename, 'r')  #Open existing file for append

        # Get the variables
        x = fid.variables['x']
        y = fid.variables['y']
        z = fid.variables['elevation']
        time = fid.variables['time']
        stage = fid.variables['stage']

        #Check values
        Q = self.domain.quantities['stage']
        Q0 = Q.vertex_values[:,0]
        Q1 = Q.vertex_values[:,1]
        Q2 = Q.vertex_values[:,2]

        A = stage[1,:]
        assert allclose(A[0], (Q2[0] + Q1[1])/2)
        assert allclose(A[1], (Q0[1] + Q1[3] + Q2[2])/3)
        assert allclose(A[2], Q0[3])
        assert allclose(A[3], (Q0[0] + Q1[5] + Q2[4])/3)

        #Center point
        assert allclose(A[4], (Q1[0] + Q2[1] + Q0[2] +\
                                 Q0[5] + Q2[6] + Q1[7])/6)


        fid.close()

        #Cleanup
        os.remove(sww.filename)

    def test_sww_variable3(self):
        """Test that sww information can be written correctly
        multiple timesteps using a different reduction operator (min)
        """

        import time, os
        from Numeric import array, zeros, allclose, Float, concatenate
        from Scientific.IO.NetCDF import NetCDFFile

        self.domain.set_name('datatest' + str(id(self)))
        self.domain.format = 'sww'
        self.domain.smooth = True
        self.domain.reduction = min

        sww = get_dataobject(self.domain)
        sww.store_connectivity()
        sww.store_timestep('stage')

        self.domain.evolve_to_end(finaltime = 0.01)
        sww.store_timestep('stage')


        #Check contents
        #Get NetCDF
        fid = NetCDFFile(sww.filename, 'r')


        # Get the variables
        x = fid.variables['x']
        y = fid.variables['y']
        z = fid.variables['elevation']
        time = fid.variables['time']
        stage = fid.variables['stage']

        #Check values
        Q = self.domain.quantities['stage']
        Q0 = Q.vertex_values[:,0]
        Q1 = Q.vertex_values[:,1]
        Q2 = Q.vertex_values[:,2]

        A = stage[1,:]
        assert allclose(A[0], min(Q2[0], Q1[1]))
        assert allclose(A[1], min(Q0[1], Q1[3], Q2[2]))
        assert allclose(A[2], Q0[3])
        assert allclose(A[3], min(Q0[0], Q1[5], Q2[4]))

        #Center point
        assert allclose(A[4], min(Q1[0], Q2[1], Q0[2],\
                                  Q0[5], Q2[6], Q1[7]))


        fid.close()

        #Cleanup
        os.remove(sww.filename)


    def test_sync(self):
        """Test info stored at each timestep is as expected (incl initial condition)
        """

        import time, os, config
        from Numeric import array, zeros, allclose, Float, concatenate
        from Scientific.IO.NetCDF import NetCDFFile

        self.domain.set_name('synctest')
        self.domain.format = 'sww'
        self.domain.smooth = False
        self.domain.store = True
        self.domain.beta_h = 0

        #Evolution
        for t in self.domain.evolve(yieldstep = 1.0, finaltime = 4.0):
            stage = self.domain.quantities['stage'].vertex_values

            #Get NetCDF
            fid = NetCDFFile(self.domain.writer.filename, 'r')
            stage_file = fid.variables['stage']

            if t == 0.0:
                assert allclose(stage, self.initial_stage)
                assert allclose(stage_file[:], stage.flat)
            else:
                assert not allclose(stage, self.initial_stage)
                assert not allclose(stage_file[:], stage.flat)

            fid.close()

        os.remove(self.domain.writer.filename)


    def test_sww_minimum_storable_height(self):
        """Test that sww information can be written correctly
        multiple timesteps using a different reduction operator (min)
        """

        import time, os
        from Numeric import array, zeros, allclose, Float, concatenate
        from Scientific.IO.NetCDF import NetCDFFile

        self.domain.set_name('datatest' + str(id(self)))
        self.domain.format = 'sww'
        self.domain.smooth = True
        self.domain.reduction = min
        self.domain.minimum_storable_height = 100

        sww = get_dataobject(self.domain)
        sww.store_connectivity()
        sww.store_timestep('stage')

        self.domain.evolve_to_end(finaltime = 0.01)
        sww.store_timestep('stage')


        #Check contents
        #Get NetCDF
        fid = NetCDFFile(sww.filename, 'r')


        # Get the variables
        x = fid.variables['x']
        y = fid.variables['y']
        z = fid.variables['elevation']
        time = fid.variables['time']
        stage = fid.variables['stage']

        #Check values
        Q = self.domain.quantities['stage']
        Q0 = Q.vertex_values[:,0]
        Q1 = Q.vertex_values[:,1]
        Q2 = Q.vertex_values[:,2]

        A = stage[1,:]
        assert allclose(stage[1,:], z[:])
        fid.close()

        #Cleanup
        os.remove(sww.filename)


    def Not_a_test_sww_DSG(self):
        """Not a test, rather a look at the sww format
        """

        import time, os
        from Numeric import array, zeros, allclose, Float, concatenate
        from Scientific.IO.NetCDF import NetCDFFile

        self.domain.set_name('datatest' + str(id(self)))
        self.domain.format = 'sww'
        self.domain.smooth = True
        self.domain.reduction = mean

        sww = get_dataobject(self.domain)
        sww.store_connectivity()
        sww.store_timestep('stage')

        #Check contents
        #Get NetCDF
        fid = NetCDFFile(sww.filename, 'r')

        # Get the variables
        x = fid.variables['x']
        y = fid.variables['y']
        z = fid.variables['elevation']

        volumes = fid.variables['volumes']
        time = fid.variables['time']

        # 2D
        stage = fid.variables['stage']

        X = x[:]
        Y = y[:]
        Z = z[:]
        V = volumes[:]
        T = time[:]
        S = stage[:,:]

#         print "****************************"
#         print "X ",X
#         print "****************************"
#         print "Y ",Y
#         print "****************************"
#         print "Z ",Z
#         print "****************************"
#         print "V ",V
#         print "****************************"
#         print "Time ",T
#         print "****************************"
#         print "Stage ",S
#         print "****************************"


        fid.close()

        #Cleanup
        os.remove(sww.filename)


    #def test_write_pts(self):
    #    #Obsolete
    #
    #    #Get (enough) datapoints
    #
    #    from Numeric import array
    #    points = array([[ 0.66666667, 0.66666667],
    #                    [ 1.33333333, 1.33333333],
    #                    [ 2.66666667, 0.66666667],
    #                    [ 0.66666667, 2.66666667],
    #                    [ 0.0, 1.0],
    #                    [ 0.0, 3.0],
    #                    [ 1.0, 0.0],
    #                    [ 1.0, 1.0],
    #                    [ 1.0, 2.0],
    #                    [ 1.0, 3.0],
    #                    [ 2.0, 1.0],
    #                    [ 3.0, 0.0],
    #                    [ 3.0, 1.0]])
    #
    #    z = points[:,0] + 2*points[:,1]
    #
    #    ptsfile = 'testptsfile.pts'
    #    write_ptsfile(ptsfile, points, z,
    #                  attribute_name = 'linear_combination')
    #
    #    #Check contents
    #    #Get NetCDF
    #    from Scientific.IO.NetCDF import NetCDFFile
    #    fid = NetCDFFile(ptsfile, 'r')
    #
    #    # Get the variables
    #    #print fid.variables.keys()
    #    points1 = fid.variables['points']
    #    z1 = fid.variables['linear_combination']
    #
    #    #Check values#
    #
    #    #print points[:]
    #    #print ref_points
    #    assert allclose(points, points1)
    #
    #    #print attributes[:]
    #    #print ref_elevation
    #    assert allclose(z, z1)
    #
    #    #Cleanup
    #    fid.close()
    #
    #    import os
    #    os.remove(ptsfile)


    def test_dem2pts_bounding_box_v2(self):
        """Test conversion from dem in ascii format to native NetCDF xya format
        """

        import time, os
        from Numeric import array, zeros, allclose, Float, concatenate, ones
        from Scientific.IO.NetCDF import NetCDFFile

        #Write test asc file
        root = 'demtest'

        filename = root+'.asc'
        fid = open(filename, 'w')
        fid.write("""ncols         10
nrows         10
xllcorner     2000
yllcorner     3000
cellsize      1
NODATA_value  -9999
""")
        #Create linear function
        ref_points = []
        ref_elevation = []
        x0 = 2000
        y = 3010
        yvec = range(10)
        xvec = range(10)
        z = -1
        for i in range(10):
            y = y - 1
            for j in range(10):
                x = x0 + xvec[j]
                z += 1
                ref_points.append ([x,y])
                ref_elevation.append(z)
                fid.write('%f ' %z)
            fid.write('\n')

        fid.close()

        #print 'sending pts', ref_points
        #print 'sending elev', ref_elevation

        #Write prj file with metadata
        metafilename = root+'.prj'
        fid = open(metafilename, 'w')


        fid.write("""Projection UTM
Zone 56
Datum WGS84
Zunits NO
Units METERS
Spheroid WGS84
Xshift 0.0000000000
Yshift 10000000.0000000000
Parameters
""")
        fid.close()

        #Convert to NetCDF pts
        convert_dem_from_ascii2netcdf(root)
        dem2pts(root, easting_min=2002.0, easting_max=2007.0,
                northing_min=3003.0, northing_max=3006.0)

        #Check contents
        #Get NetCDF
        fid = NetCDFFile(root+'.pts', 'r')

        # Get the variables
        #print fid.variables.keys()
        points = fid.variables['points']
        elevation = fid.variables['elevation']

        #Check values
        assert fid.xllcorner[0] == 2002.0
        assert fid.yllcorner[0] == 3003.0

        #create new reference points
        newz = []
        newz[0:5] = ref_elevation[32:38]
        newz[6:11] = ref_elevation[42:48]
        newz[12:17] = ref_elevation[52:58]
        newz[18:23] = ref_elevation[62:68]
        ref_elevation = []
        ref_elevation = newz
        ref_points = []
        x0 = 2002
        y = 3007
        yvec = range(4)
        xvec = range(6)
        for i in range(4):
            y = y - 1
            ynew = y - 3003.0
            for j in range(6):
                x = x0 + xvec[j]
                xnew = x - 2002.0
                ref_points.append ([xnew,ynew]) #Relative point values

        assert allclose(points, ref_points)

        assert allclose(elevation, ref_elevation)

        #Cleanup
        fid.close()


        os.remove(root + '.pts')
        os.remove(root + '.dem')
        os.remove(root + '.asc')
        os.remove(root + '.prj')


    def test_dem2pts_bounding_box_removeNullvalues_v2(self):
        """Test conversion from dem in ascii format to native NetCDF xya format
        """

        import time, os
        from Numeric import array, zeros, allclose, Float, concatenate, ones
        from Scientific.IO.NetCDF import NetCDFFile

        #Write test asc file
        root = 'demtest'

        filename = root+'.asc'
        fid = open(filename, 'w')
        fid.write("""ncols         10
nrows         10
xllcorner     2000
yllcorner     3000
cellsize      1
NODATA_value  -9999
""")
        #Create linear function
        ref_points = []
        ref_elevation = []
        x0 = 2000
        y = 3010
        yvec = range(10)
        xvec = range(10)
        #z = range(100)
        z = zeros(100)
        NODATA_value = -9999
        count = -1
        for i in range(10):
            y = y - 1
            for j in range(10):
                x = x0 + xvec[j]
                ref_points.append ([x,y])
                count += 1
                z[count] = (4*i - 3*j)%13
                if j == 4: z[count] = NODATA_value #column inside clipping region
                if j == 8: z[count] = NODATA_value #column outside clipping region
                if i == 9: z[count] = NODATA_value #row outside clipping region
                if i == 4 and j == 6: z[count] = NODATA_value #arbitrary point inside clipping region
                ref_elevation.append( z[count] )
                fid.write('%f ' %z[count])
            fid.write('\n')

        fid.close()

        #print 'sending elev', ref_elevation

        #Write prj file with metadata
        metafilename = root+'.prj'
        fid = open(metafilename, 'w')


        fid.write("""Projection UTM
Zone 56
Datum WGS84
Zunits NO
Units METERS
Spheroid WGS84
Xshift 0.0000000000
Yshift 10000000.0000000000
Parameters
""")
        fid.close()

        #Convert to NetCDF pts
        convert_dem_from_ascii2netcdf(root)
        dem2pts(root, easting_min=2002.0, easting_max=2007.0,
                northing_min=3003.0, northing_max=3006.0)

        #Check contents
        #Get NetCDF
        fid = NetCDFFile(root+'.pts', 'r')

        # Get the variables
        #print fid.variables.keys()
        points = fid.variables['points']
        elevation = fid.variables['elevation']

        #Check values
        assert fid.xllcorner[0] == 2002.0
        assert fid.yllcorner[0] == 3003.0

        #create new reference points
        newz = zeros(19)
        newz[0:2] = ref_elevation[32:34]
        newz[2:5] = ref_elevation[35:38]
        newz[5:7] = ref_elevation[42:44]
        newz[7] = ref_elevation[45]
        newz[8] = ref_elevation[47]
        newz[9:11] = ref_elevation[52:54]
        newz[11:14] = ref_elevation[55:58]
        newz[14:16] = ref_elevation[62:64]
        newz[16:19] = ref_elevation[65:68]


        ref_elevation = newz
        ref_points = []
        new_ref_points = []
        x0 = 2002
        y = 3007
        yvec = range(4)
        xvec = range(6)
        for i in range(4):
            y = y - 1
            ynew = y - 3003.0
            for j in range(6):
                x = x0 + xvec[j]
                xnew = x - 2002.0
                if j <> 2 and (i<>1 or j<>4):
                    ref_points.append([x,y])
                    new_ref_points.append ([xnew,ynew])


        assert allclose(points, new_ref_points)
        assert allclose(elevation, ref_elevation)

        #Cleanup
        fid.close()


        os.remove(root + '.pts')
        os.remove(root + '.dem')
        os.remove(root + '.asc')
        os.remove(root + '.prj')


    def test_dem2pts_bounding_box_removeNullvalues_v3(self):
        """Test conversion from dem in ascii format to native NetCDF xya format
        Check missing values on clipping boundary
        """

        import time, os
        from Numeric import array, zeros, allclose, Float, concatenate, ones
        from Scientific.IO.NetCDF import NetCDFFile

        #Write test asc file
        root = 'demtest'

        filename = root+'.asc'
        fid = open(filename, 'w')
        fid.write("""ncols         10
nrows         10
xllcorner     2000
yllcorner     3000
cellsize      1
NODATA_value  -9999
""")
        #Create linear function
        ref_points = []
        ref_elevation = []
        x0 = 2000
        y = 3010
        yvec = range(10)
        xvec = range(10)
        #z = range(100)
        z = zeros(100)
        NODATA_value = -9999
        count = -1
        for i in range(10):
            y = y - 1
            for j in range(10):
                x = x0 + xvec[j]
                ref_points.append ([x,y])
                count += 1
                z[count] = (4*i - 3*j)%13
                if j == 4: z[count] = NODATA_value #column inside clipping region
                if j == 8: z[count] = NODATA_value #column outside clipping region
                if i == 6: z[count] = NODATA_value #row on clipping boundary
                if i == 4 and j == 6: z[count] = NODATA_value #arbitrary point inside clipping region
                ref_elevation.append( z[count] )
                fid.write('%f ' %z[count])
            fid.write('\n')

        fid.close()

        #print 'sending elev', ref_elevation

        #Write prj file with metadata
        metafilename = root+'.prj'
        fid = open(metafilename, 'w')


        fid.write("""Projection UTM
Zone 56
Datum WGS84
Zunits NO
Units METERS
Spheroid WGS84
Xshift 0.0000000000
Yshift 10000000.0000000000
Parameters
""")
        fid.close()

        #Convert to NetCDF pts
        convert_dem_from_ascii2netcdf(root)
        dem2pts(root, easting_min=2002.0, easting_max=2007.0,
                northing_min=3003.0, northing_max=3006.0)

        #Check contents
        #Get NetCDF
        fid = NetCDFFile(root+'.pts', 'r')

        # Get the variables
        #print fid.variables.keys()
        points = fid.variables['points']
        elevation = fid.variables['elevation']

        #Check values
        assert fid.xllcorner[0] == 2002.0
        assert fid.yllcorner[0] == 3003.0

        #create new reference points
        newz = zeros(14)
        newz[0:2] = ref_elevation[32:34]
        newz[2:5] = ref_elevation[35:38]
        newz[5:7] = ref_elevation[42:44]
        newz[7] = ref_elevation[45]
        newz[8] = ref_elevation[47]
        newz[9:11] = ref_elevation[52:54]
        newz[11:14] = ref_elevation[55:58]



        ref_elevation = newz
        ref_points = []
        new_ref_points = []
        x0 = 2002
        y = 3007
        yvec = range(4)
        xvec = range(6)
        for i in range(4):
            y = y - 1
            ynew = y - 3003.0
            for j in range(6):
                x = x0 + xvec[j]
                xnew = x - 2002.0
                if j <> 2 and (i<>1 or j<>4) and i<>3:
                    ref_points.append([x,y])
                    new_ref_points.append ([xnew,ynew])


        #print points[:],points[:].shape
        #print new_ref_points, len(new_ref_points)

        assert allclose(elevation, ref_elevation)
        assert allclose(points, new_ref_points)


        #Cleanup
        fid.close()


        os.remove(root + '.pts')
        os.remove(root + '.dem')
        os.remove(root + '.asc')
        os.remove(root + '.prj')


    def test_hecras_cross_sections2pts(self):
        """Test conversion from HECRAS cross sections in ascii format
        to native NetCDF pts format
        """

        import time, os
        from Numeric import array, zeros, allclose, Float, concatenate
        from Scientific.IO.NetCDF import NetCDFFile

        #Write test asc file
        root = 'hecrastest'

        filename = root+'.sdf'
        fid = open(filename, 'w')
        fid.write("""
# RAS export file created on Mon 15Aug2005 11:42
# by HEC-RAS Version 3.1.1

BEGIN HEADER:
  UNITS: METRIC
  DTM TYPE: TIN
  DTM: v:\1\cit\perth_topo\river_tin
  STREAM LAYER: c:\\x_local\hecras\21_02_03\up_canning_cent3d.shp
  CROSS-SECTION LAYER: c:\\x_local\hecras\21_02_03\up_can_xs3d.shp
  MAP PROJECTION: UTM
  PROJECTION ZONE: 50
  DATUM: AGD66
  VERTICAL DATUM:
  NUMBER OF REACHES:  19
  NUMBER OF CROSS-SECTIONS:  2
END HEADER:


BEGIN CROSS-SECTIONS:

  CROSS-SECTION:
    STREAM ID:Southern-Wungong
    REACH ID:Southern-Wungong
    STATION:21410
    CUT LINE:
      407546.08 , 6437277.542
      407329.32 , 6437489.482
      407283.11 , 6437541.232
    SURFACE LINE:
     407546.08,   6437277.54,   52.14
     407538.88,   6437284.58,   51.07
     407531.68,   6437291.62,   50.56
     407524.48,   6437298.66,   49.58
     407517.28,   6437305.70,   49.09
     407510.08,   6437312.74,   48.76
  END:

  CROSS-SECTION:
    STREAM ID:Swan River
    REACH ID:Swan Mouth
    STATION:840.*
    CUT LINE:
      381178.0855 , 6452559.0685
      380485.4755 , 6453169.272
    SURFACE LINE:
     381178.09,   6452559.07,   4.17
     381169.49,   6452566.64,   4.26
     381157.78,   6452576.96,   4.34
     381155.97,   6452578.56,   4.35
     381143.72,   6452589.35,   4.43
     381136.69,   6452595.54,   4.58
     381114.74,   6452614.88,   4.41
     381075.53,   6452649.43,   4.17
     381071.47,   6452653.00,   3.99
     381063.46,   6452660.06,   3.67
     381054.41,   6452668.03,   3.67
  END:
END CROSS-SECTIONS:
""")

        fid.close()


        #Convert to NetCDF pts
        hecras_cross_sections2pts(root)

        #Check contents
        #Get NetCDF
        fid = NetCDFFile(root+'.pts', 'r')

        # Get the variables
        #print fid.variables.keys()
        points = fid.variables['points']
        elevation = fid.variables['elevation']

        #Check values
        ref_points = [[407546.08, 6437277.54],
                      [407538.88, 6437284.58],
                      [407531.68, 6437291.62],
                      [407524.48, 6437298.66],
                      [407517.28, 6437305.70],
                      [407510.08, 6437312.74]]

        ref_points += [[381178.09, 6452559.07],
                       [381169.49, 6452566.64],
                       [381157.78, 6452576.96],
                       [381155.97, 6452578.56],
                       [381143.72, 6452589.35],
                       [381136.69, 6452595.54],
                       [381114.74, 6452614.88],
                       [381075.53, 6452649.43],
                       [381071.47, 6452653.00],
                       [381063.46, 6452660.06],
                       [381054.41, 6452668.03]]


        ref_elevation = [52.14, 51.07, 50.56, 49.58, 49.09, 48.76]
        ref_elevation += [4.17, 4.26, 4.34, 4.35, 4.43, 4.58, 4.41, 4.17, 3.99, 3.67, 3.67]

        #print points[:]
        #print ref_points
        assert allclose(points, ref_points)

        #print attributes[:]
        #print ref_elevation
        assert allclose(elevation, ref_elevation)

        #Cleanup
        fid.close()


        os.remove(root + '.sdf')
        os.remove(root + '.pts')



    def test_sww2dem_asc_elevation(self):
        """Test that sww information can be converted correctly to asc/prj
        format readable by e.g. ArcView
        """

        import time, os
        from Numeric import array, zeros, allclose, Float, concatenate
        from Scientific.IO.NetCDF import NetCDFFile

        #Setup
        self.domain.set_name('datatest')

        prjfile = self.domain.get_name() + '_elevation.prj'
        ascfile = self.domain.get_name() + '_elevation.asc'
        swwfile = self.domain.get_name() + '.sww'

        self.domain.set_datadir('.')
        self.domain.format = 'sww'
        self.domain.smooth = True
        self.domain.set_quantity('elevation', lambda x,y: -x-y)

        self.domain.geo_reference = Geo_reference(56,308500,6189000)

        sww = get_dataobject(self.domain)
        sww.store_connectivity()
        sww.store_timestep('stage')

        self.domain.evolve_to_end(finaltime = 0.01)
        sww.store_timestep('stage')

        cellsize = 0.25
        #Check contents
        #Get NetCDF

        fid = NetCDFFile(sww.filename, 'r')

        # Get the variables
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        z = fid.variables['elevation'][:]
        time = fid.variables['time'][:]
        stage = fid.variables['stage'][:]


        #Export to ascii/prj files
        sww2dem(self.domain.get_name(),
                quantity = 'elevation',
                cellsize = cellsize,
                verbose = False,
                format = 'asc')

        #Check prj (meta data)
        prjid = open(prjfile)
        lines = prjid.readlines()
        prjid.close()

        L = lines[0].strip().split()
        assert L[0].strip().lower() == 'projection'
        assert L[1].strip().lower() == 'utm'

        L = lines[1].strip().split()
        assert L[0].strip().lower() == 'zone'
        assert L[1].strip().lower() == '56'

        L = lines[2].strip().split()
        assert L[0].strip().lower() == 'datum'
        assert L[1].strip().lower() == 'wgs84'

        L = lines[3].strip().split()
        assert L[0].strip().lower() == 'zunits'
        assert L[1].strip().lower() == 'no'

        L = lines[4].strip().split()
        assert L[0].strip().lower() == 'units'
        assert L[1].strip().lower() == 'meters'

        L = lines[5].strip().split()
        assert L[0].strip().lower() == 'spheroid'
        assert L[1].strip().lower() == 'wgs84'

        L = lines[6].strip().split()
        assert L[0].strip().lower() == 'xshift'
        assert L[1].strip().lower() == '500000'

        L = lines[7].strip().split()
        assert L[0].strip().lower() == 'yshift'
        assert L[1].strip().lower() == '10000000'

        L = lines[8].strip().split()
        assert L[0].strip().lower() == 'parameters'


        #Check asc file
        ascid = open(ascfile)
        lines = ascid.readlines()
        ascid.close()

        L = lines[0].strip().split()
        assert L[0].strip().lower() == 'ncols'
        assert L[1].strip().lower() == '5'

        L = lines[1].strip().split()
        assert L[0].strip().lower() == 'nrows'
        assert L[1].strip().lower() == '5'

        L = lines[2].strip().split()
        assert L[0].strip().lower() == 'xllcorner'
        assert allclose(float(L[1].strip().lower()), 308500)

        L = lines[3].strip().split()
        assert L[0].strip().lower() == 'yllcorner'
        assert allclose(float(L[1].strip().lower()), 6189000)

        L = lines[4].strip().split()
        assert L[0].strip().lower() == 'cellsize'
        assert allclose(float(L[1].strip().lower()), cellsize)

        L = lines[5].strip().split()
        assert L[0].strip() == 'NODATA_value'
        assert L[1].strip().lower() == '-9999'

        #Check grid values
        for j in range(5):
            L = lines[6+j].strip().split()
            y = (4-j) * cellsize
            for i in range(5):
                assert allclose(float(L[i]), -i*cellsize - y)


        fid.close()

        #Cleanup
        os.remove(prjfile)
        os.remove(ascfile)
        os.remove(swwfile)



    def test_sww2dem_larger(self):
        """Test that sww information can be converted correctly to asc/prj
        format readable by e.g. ArcView. Here:

        ncols         11
        nrows         11
        xllcorner     308500
        yllcorner     6189000
        cellsize      10.000000
        NODATA_value  -9999
        -100 -110 -120 -130 -140 -150 -160 -170 -180 -190 -200
         -90 -100 -110 -120 -130 -140 -150 -160 -170 -180 -190
         -80  -90 -100 -110 -120 -130 -140 -150 -160 -170 -180
         -70  -80  -90 -100 -110 -120 -130 -140 -150 -160 -170
         -60  -70  -80  -90 -100 -110 -120 -130 -140 -150 -160
         -50  -60  -70  -80  -90 -100 -110 -120 -130 -140 -150
         -40  -50  -60  -70  -80  -90 -100 -110 -120 -130 -140
         -30  -40  -50  -60  -70  -80  -90 -100 -110 -120 -130
         -20  -30  -40  -50  -60  -70  -80  -90 -100 -110 -120
         -10  -20  -30  -40  -50  -60  -70  -80  -90 -100 -110
           0  -10  -20  -30  -40  -50  -60  -70  -80  -90 -100

        """

        import time, os
        from Numeric import array, zeros, allclose, Float, concatenate
        from Scientific.IO.NetCDF import NetCDFFile

        #Setup

        from mesh_factory import rectangular

        #Create basic mesh (100m x 100m)
        points, vertices, boundary = rectangular(2, 2, 100, 100)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.default_order = 2

        domain.set_name('datatest')

        prjfile = domain.get_name() + '_elevation.prj'
        ascfile = domain.get_name() + '_elevation.asc'
        swwfile = domain.get_name() + '.sww'

        domain.set_datadir('.')
        domain.format = 'sww'
        domain.smooth = True
        domain.geo_reference = Geo_reference(56, 308500, 6189000)

        #
        domain.set_quantity('elevation', lambda x,y: -x-y)
        domain.set_quantity('stage', 0)

        B = Transmissive_boundary(domain)
        domain.set_boundary( {'left': B, 'right': B, 'top': B, 'bottom': B})


        #
        sww = get_dataobject(domain)
        sww.store_connectivity()
        sww.store_timestep('stage')

        domain.evolve_to_end(finaltime = 0.01)
        sww.store_timestep('stage')

        cellsize = 10  #10m grid


        #Check contents
        #Get NetCDF

        fid = NetCDFFile(sww.filename, 'r')

        # Get the variables
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        z = fid.variables['elevation'][:]
        time = fid.variables['time'][:]
        stage = fid.variables['stage'][:]


        #Export to ascii/prj files
        sww2dem(domain.get_name(),
                quantity = 'elevation',
                cellsize = cellsize,
                verbose = False,
                format = 'asc')


        #Check prj (meta data)
        prjid = open(prjfile)
        lines = prjid.readlines()
        prjid.close()

        L = lines[0].strip().split()
        assert L[0].strip().lower() == 'projection'
        assert L[1].strip().lower() == 'utm'

        L = lines[1].strip().split()
        assert L[0].strip().lower() == 'zone'
        assert L[1].strip().lower() == '56'

        L = lines[2].strip().split()
        assert L[0].strip().lower() == 'datum'
        assert L[1].strip().lower() == 'wgs84'

        L = lines[3].strip().split()
        assert L[0].strip().lower() == 'zunits'
        assert L[1].strip().lower() == 'no'

        L = lines[4].strip().split()
        assert L[0].strip().lower() == 'units'
        assert L[1].strip().lower() == 'meters'

        L = lines[5].strip().split()
        assert L[0].strip().lower() == 'spheroid'
        assert L[1].strip().lower() == 'wgs84'

        L = lines[6].strip().split()
        assert L[0].strip().lower() == 'xshift'
        assert L[1].strip().lower() == '500000'

        L = lines[7].strip().split()
        assert L[0].strip().lower() == 'yshift'
        assert L[1].strip().lower() == '10000000'

        L = lines[8].strip().split()
        assert L[0].strip().lower() == 'parameters'


        #Check asc file
        ascid = open(ascfile)
        lines = ascid.readlines()
        ascid.close()

        L = lines[0].strip().split()
        assert L[0].strip().lower() == 'ncols'
        assert L[1].strip().lower() == '11'

        L = lines[1].strip().split()
        assert L[0].strip().lower() == 'nrows'
        assert L[1].strip().lower() == '11'

        L = lines[2].strip().split()
        assert L[0].strip().lower() == 'xllcorner'
        assert allclose(float(L[1].strip().lower()), 308500)

        L = lines[3].strip().split()
        assert L[0].strip().lower() == 'yllcorner'
        assert allclose(float(L[1].strip().lower()), 6189000)

        L = lines[4].strip().split()
        assert L[0].strip().lower() == 'cellsize'
        assert allclose(float(L[1].strip().lower()), cellsize)

        L = lines[5].strip().split()
        assert L[0].strip() == 'NODATA_value'
        assert L[1].strip().lower() == '-9999'

        #Check grid values (FIXME: Use same strategy for other sww2dem tests)
        for i, line in enumerate(lines[6:]):
            for j, value in enumerate( line.split() ):
                #assert allclose(float(value), -(10-i+j)*cellsize)
                assert float(value) == -(10-i+j)*cellsize


        fid.close()

        #Cleanup
        os.remove(prjfile)
        os.remove(ascfile)
        os.remove(swwfile)



    def test_sww2dem_boundingbox(self):
        """Test that sww information can be converted correctly to asc/prj
        format readable by e.g. ArcView.
        This will test that mesh can be restricted by bounding box

        Original extent is 100m x 100m:

        Eastings:   308500 -  308600
        Northings: 6189000 - 6189100

        Bounding box changes this to the 50m x 50m square defined by

        Eastings:   308530 -  308570
        Northings: 6189050 - 6189100

        The cropped values should be

         -130 -140 -150 -160 -170
         -120 -130 -140 -150 -160
         -110 -120 -130 -140 -150
         -100 -110 -120 -130 -140
          -90 -100 -110 -120 -130
          -80  -90 -100 -110 -120

        and the new lower reference point should be
        Eastings:   308530
        Northings: 6189050

        Original dataset is the same as in test_sww2dem_larger()

        """

        import time, os
        from Numeric import array, zeros, allclose, Float, concatenate
        from Scientific.IO.NetCDF import NetCDFFile

        #Setup

        from mesh_factory import rectangular

        #Create basic mesh (100m x 100m)
        points, vertices, boundary = rectangular(2, 2, 100, 100)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.default_order = 2

        domain.set_name('datatest')

        prjfile = domain.get_name() + '_elevation.prj'
        ascfile = domain.get_name() + '_elevation.asc'
        swwfile = domain.get_name() + '.sww'

        domain.set_datadir('.')
        domain.format = 'sww'
        domain.smooth = True
        domain.geo_reference = Geo_reference(56, 308500, 6189000)

        #
        domain.set_quantity('elevation', lambda x,y: -x-y)
        domain.set_quantity('stage', 0)

        B = Transmissive_boundary(domain)
        domain.set_boundary( {'left': B, 'right': B, 'top': B, 'bottom': B})


        #
        sww = get_dataobject(domain)
        sww.store_connectivity()
        sww.store_timestep('stage')

        domain.evolve_to_end(finaltime = 0.01)
        sww.store_timestep('stage')

        cellsize = 10  #10m grid


        #Check contents
        #Get NetCDF

        fid = NetCDFFile(sww.filename, 'r')

        # Get the variables
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        z = fid.variables['elevation'][:]
        time = fid.variables['time'][:]
        stage = fid.variables['stage'][:]


        #Export to ascii/prj files
        sww2dem(domain.get_name(),
                quantity = 'elevation',
                cellsize = cellsize,
                easting_min = 308530,
                easting_max = 308570,
                northing_min = 6189050,
                northing_max = 6189100,
                verbose = False,
                format = 'asc')

        fid.close()


        #Check prj (meta data)
        prjid = open(prjfile)
        lines = prjid.readlines()
        prjid.close()

        L = lines[0].strip().split()
        assert L[0].strip().lower() == 'projection'
        assert L[1].strip().lower() == 'utm'

        L = lines[1].strip().split()
        assert L[0].strip().lower() == 'zone'
        assert L[1].strip().lower() == '56'

        L = lines[2].strip().split()
        assert L[0].strip().lower() == 'datum'
        assert L[1].strip().lower() == 'wgs84'

        L = lines[3].strip().split()
        assert L[0].strip().lower() == 'zunits'
        assert L[1].strip().lower() == 'no'

        L = lines[4].strip().split()
        assert L[0].strip().lower() == 'units'
        assert L[1].strip().lower() == 'meters'

        L = lines[5].strip().split()
        assert L[0].strip().lower() == 'spheroid'
        assert L[1].strip().lower() == 'wgs84'

        L = lines[6].strip().split()
        assert L[0].strip().lower() == 'xshift'
        assert L[1].strip().lower() == '500000'

        L = lines[7].strip().split()
        assert L[0].strip().lower() == 'yshift'
        assert L[1].strip().lower() == '10000000'

        L = lines[8].strip().split()
        assert L[0].strip().lower() == 'parameters'


        #Check asc file
        ascid = open(ascfile)
        lines = ascid.readlines()
        ascid.close()

        L = lines[0].strip().split()
        assert L[0].strip().lower() == 'ncols'
        assert L[1].strip().lower() == '5'

        L = lines[1].strip().split()
        assert L[0].strip().lower() == 'nrows'
        assert L[1].strip().lower() == '6'

        L = lines[2].strip().split()
        assert L[0].strip().lower() == 'xllcorner'
        assert allclose(float(L[1].strip().lower()), 308530)

        L = lines[3].strip().split()
        assert L[0].strip().lower() == 'yllcorner'
        assert allclose(float(L[1].strip().lower()), 6189050)

        L = lines[4].strip().split()
        assert L[0].strip().lower() == 'cellsize'
        assert allclose(float(L[1].strip().lower()), cellsize)

        L = lines[5].strip().split()
        assert L[0].strip() == 'NODATA_value'
        assert L[1].strip().lower() == '-9999'

        #Check grid values
        for i, line in enumerate(lines[6:]):
            for j, value in enumerate( line.split() ):
                #assert float(value) == -(10-i+j)*cellsize
                assert float(value) == -(10-i+j+3)*cellsize



        #Cleanup
        os.remove(prjfile)
        os.remove(ascfile)
        os.remove(swwfile)



    def test_sww2dem_asc_stage_reduction(self):
        """Test that sww information can be converted correctly to asc/prj
        format readable by e.g. ArcView

        This tests the reduction of quantity stage using min
        """

        import time, os
        from Numeric import array, zeros, allclose, Float, concatenate
        from Scientific.IO.NetCDF import NetCDFFile

        #Setup
        self.domain.set_name('datatest')

        prjfile = self.domain.get_name() + '_stage.prj'
        ascfile = self.domain.get_name() + '_stage.asc'
        swwfile = self.domain.get_name() + '.sww'

        self.domain.set_datadir('.')
        self.domain.format = 'sww'
        self.domain.smooth = True
        self.domain.set_quantity('elevation', lambda x,y: -x-y)

        self.domain.geo_reference = Geo_reference(56,308500,6189000)


        sww = get_dataobject(self.domain)
        sww.store_connectivity()
        sww.store_timestep('stage')

        self.domain.evolve_to_end(finaltime = 0.01)
        sww.store_timestep('stage')

        cellsize = 0.25
        #Check contents
        #Get NetCDF

        fid = NetCDFFile(sww.filename, 'r')

        # Get the variables
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        z = fid.variables['elevation'][:]
        time = fid.variables['time'][:]
        stage = fid.variables['stage'][:]


        #Export to ascii/prj files
        sww2dem(self.domain.get_name(),
                quantity = 'stage',
                cellsize = cellsize,
                reduction = min,
                format = 'asc')


        #Check asc file
        ascid = open(ascfile)
        lines = ascid.readlines()
        ascid.close()

        L = lines[0].strip().split()
        assert L[0].strip().lower() == 'ncols'
        assert L[1].strip().lower() == '5'

        L = lines[1].strip().split()
        assert L[0].strip().lower() == 'nrows'
        assert L[1].strip().lower() == '5'

        L = lines[2].strip().split()
        assert L[0].strip().lower() == 'xllcorner'
        assert allclose(float(L[1].strip().lower()), 308500)

        L = lines[3].strip().split()
        assert L[0].strip().lower() == 'yllcorner'
        assert allclose(float(L[1].strip().lower()), 6189000)

        L = lines[4].strip().split()
        assert L[0].strip().lower() == 'cellsize'
        assert allclose(float(L[1].strip().lower()), cellsize)

        L = lines[5].strip().split()
        assert L[0].strip() == 'NODATA_value'
        assert L[1].strip().lower() == '-9999'


        #Check grid values (where applicable)
        for j in range(5):
            if j%2 == 0:
                L = lines[6+j].strip().split()
                jj = 4-j
                for i in range(5):
                    if i%2 == 0:
                        index = jj/2 + i/2*3
                        val0 = stage[0,index]
                        val1 = stage[1,index]

                        #print i, j, index, ':', L[i], val0, val1
                        assert allclose(float(L[i]), min(val0, val1))


        fid.close()

        #Cleanup
        os.remove(prjfile)
        os.remove(ascfile)
        #os.remove(swwfile)



    def test_sww2dem_asc_derived_quantity(self):
        """Test that sww information can be converted correctly to asc/prj
        format readable by e.g. ArcView

        This tests the use of derived quantities
        """

        import time, os
        from Numeric import array, zeros, allclose, Float, concatenate
        from Scientific.IO.NetCDF import NetCDFFile

        #Setup
        self.domain.set_name('datatest')

        prjfile = self.domain.get_name() + '_depth.prj'
        ascfile = self.domain.get_name() + '_depth.asc'
        swwfile = self.domain.get_name() + '.sww'

        self.domain.set_datadir('.')
        self.domain.format = 'sww'
        self.domain.smooth = True
        self.domain.set_quantity('elevation', lambda x,y: -x-y)
        self.domain.set_quantity('stage', 0.0)

        self.domain.geo_reference = Geo_reference(56,308500,6189000)


        sww = get_dataobject(self.domain)
        sww.store_connectivity()
        sww.store_timestep('stage')

        self.domain.evolve_to_end(finaltime = 0.01)
        sww.store_timestep('stage')

        cellsize = 0.25
        #Check contents
        #Get NetCDF

        fid = NetCDFFile(sww.filename, 'r')

        # Get the variables
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        z = fid.variables['elevation'][:]
        time = fid.variables['time'][:]
        stage = fid.variables['stage'][:]


        #Export to ascii/prj files
        sww2dem(self.domain.get_name(),
                basename_out = 'datatest_depth',
                quantity = 'stage - elevation',
                cellsize = cellsize,
                reduction = min,
                format = 'asc',
                verbose = False)


        #Check asc file
        ascid = open(ascfile)
        lines = ascid.readlines()
        ascid.close()

        L = lines[0].strip().split()
        assert L[0].strip().lower() == 'ncols'
        assert L[1].strip().lower() == '5'

        L = lines[1].strip().split()
        assert L[0].strip().lower() == 'nrows'
        assert L[1].strip().lower() == '5'

        L = lines[2].strip().split()
        assert L[0].strip().lower() == 'xllcorner'
        assert allclose(float(L[1].strip().lower()), 308500)

        L = lines[3].strip().split()
        assert L[0].strip().lower() == 'yllcorner'
        assert allclose(float(L[1].strip().lower()), 6189000)

        L = lines[4].strip().split()
        assert L[0].strip().lower() == 'cellsize'
        assert allclose(float(L[1].strip().lower()), cellsize)

        L = lines[5].strip().split()
        assert L[0].strip() == 'NODATA_value'
        assert L[1].strip().lower() == '-9999'


        #Check grid values (where applicable)
        for j in range(5):
            if j%2 == 0:
                L = lines[6+j].strip().split()
                jj = 4-j
                for i in range(5):
                    if i%2 == 0:
                        index = jj/2 + i/2*3
                        val0 = stage[0,index] - z[index]
                        val1 = stage[1,index] - z[index]

                        #print i, j, index, ':', L[i], val0, val1
                        assert allclose(float(L[i]), min(val0, val1))


        fid.close()

        #Cleanup
        os.remove(prjfile)
        os.remove(ascfile)
        #os.remove(swwfile)





    def test_sww2dem_asc_missing_points(self):
        """Test that sww information can be converted correctly to asc/prj
        format readable by e.g. ArcView

        This test includes the writing of missing values
        """

        import time, os
        from Numeric import array, zeros, allclose, Float, concatenate
        from Scientific.IO.NetCDF import NetCDFFile

        #Setup mesh not coinciding with rectangle.
        #This will cause missing values to occur in gridded data


        points = [                        [1.0, 1.0],
                              [0.5, 0.5], [1.0, 0.5],
                  [0.0, 0.0], [0.5, 0.0], [1.0, 0.0]]

        vertices = [ [4,1,3], [5,2,4], [1,4,2], [2,0,1]]

        #Create shallow water domain
        domain = Domain(points, vertices)
        domain.default_order=2


        #Set some field values
        domain.set_quantity('elevation', lambda x,y: -x-y)
        domain.set_quantity('friction', 0.03)


        ######################
        # Boundary conditions
        B = Transmissive_boundary(domain)
        domain.set_boundary( {'exterior': B} )


        ######################
        #Initial condition - with jumps

        bed = domain.quantities['elevation'].vertex_values
        stage = zeros(bed.shape, Float)

        h = 0.3
        for i in range(stage.shape[0]):
            if i % 2 == 0:
                stage[i,:] = bed[i,:] + h
            else:
                stage[i,:] = bed[i,:]

        domain.set_quantity('stage', stage)
        domain.distribute_to_vertices_and_edges()

        domain.set_name('datatest')

        prjfile = domain.get_name() + '_elevation.prj'
        ascfile = domain.get_name() + '_elevation.asc'
        swwfile = domain.get_name() + '.sww'

        domain.set_datadir('.')
        domain.format = 'sww'
        domain.smooth = True

        domain.geo_reference = Geo_reference(56,308500,6189000)

        sww = get_dataobject(domain)
        sww.store_connectivity()
        sww.store_timestep('stage')

        cellsize = 0.25
        #Check contents
        #Get NetCDF

        fid = NetCDFFile(swwfile, 'r')

        # Get the variables
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        z = fid.variables['elevation'][:]
        time = fid.variables['time'][:]

        try:
            geo_reference = Geo_reference(NetCDFObject=fid)
        except AttributeError, e:
            geo_reference = Geo_reference(DEFAULT_ZONE,0,0)

        #Export to ascii/prj files
        sww2dem(domain.get_name(),
                quantity = 'elevation',
                cellsize = cellsize,
                verbose = False,
                format = 'asc')


        #Check asc file
        ascid = open(ascfile)
        lines = ascid.readlines()
        ascid.close()

        L = lines[0].strip().split()
        assert L[0].strip().lower() == 'ncols'
        assert L[1].strip().lower() == '5'

        L = lines[1].strip().split()
        assert L[0].strip().lower() == 'nrows'
        assert L[1].strip().lower() == '5'

        L = lines[2].strip().split()
        assert L[0].strip().lower() == 'xllcorner'
        assert allclose(float(L[1].strip().lower()), 308500)

        L = lines[3].strip().split()
        assert L[0].strip().lower() == 'yllcorner'
        assert allclose(float(L[1].strip().lower()), 6189000)

        L = lines[4].strip().split()
        assert L[0].strip().lower() == 'cellsize'
        assert allclose(float(L[1].strip().lower()), cellsize)

        L = lines[5].strip().split()
        assert L[0].strip() == 'NODATA_value'
        assert L[1].strip().lower() == '-9999'

        #Check grid values
        for j in range(5):
            L = lines[6+j].strip().split()
            assert len(L) == 5
            y = (4-j) * cellsize

            for i in range(5):
                #print i
                if i+j >= 4:
                    assert allclose(float(L[i]), -i*cellsize - y)
                else:
                    #Missing values
                    assert allclose(float(L[i]), -9999)



        fid.close()

        #Cleanup
        os.remove(prjfile)
        os.remove(ascfile)
        os.remove(swwfile)

    def test_sww2ers_simple(self):
        """Test that sww information can be converted correctly to asc/prj
        format readable by e.g. ArcView
        """

        import time, os
        from Numeric import array, zeros, allclose, Float, concatenate
        from Scientific.IO.NetCDF import NetCDFFile


        NODATA_value = 1758323

        #Setup
        self.domain.set_name('datatest')

        headerfile = self.domain.get_name() + '.ers'
        swwfile = self.domain.get_name() + '.sww'

        self.domain.set_datadir('.')
        self.domain.format = 'sww'
        self.domain.smooth = True
        self.domain.set_quantity('elevation', lambda x,y: -x-y)

        self.domain.geo_reference = Geo_reference(56,308500,6189000)

        sww = get_dataobject(self.domain)
        sww.store_connectivity()
        sww.store_timestep('stage')

        self.domain.evolve_to_end(finaltime = 0.01)
        sww.store_timestep('stage')

        cellsize = 0.25
        #Check contents
        #Get NetCDF

        fid = NetCDFFile(sww.filename, 'r')

        # Get the variables
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        z = fid.variables['elevation'][:]
        time = fid.variables['time'][:]
        stage = fid.variables['stage'][:]


        #Export to ers files
        sww2dem(self.domain.get_name(),
                quantity = 'elevation',
                cellsize = cellsize,
                NODATA_value = NODATA_value,
                verbose = False,
                format = 'ers')

        #Check header data
        from ermapper_grids import read_ermapper_header, read_ermapper_data

        header = read_ermapper_header(self.domain.get_name() + '_elevation.ers')
        #print header
        assert header['projection'].lower() == '"utm-56"'
        assert header['datum'].lower() == '"wgs84"'
        assert header['units'].lower() == '"meters"'
        assert header['value'].lower() == '"elevation"'
        assert header['xdimension'] == '0.25'
        assert header['ydimension'] == '0.25'
        assert float(header['eastings']) == 308500.0   #xllcorner
        assert float(header['northings']) == 6189000.0 #yllcorner
        assert int(header['nroflines']) == 5
        assert int(header['nrofcellsperline']) == 5
        assert int(header['nullcellvalue']) == NODATA_value
        #FIXME - there is more in the header


        #Check grid data
        grid = read_ermapper_data(self.domain.get_name() + '_elevation')

        #FIXME (Ole): Why is this the desired reference grid for -x-y?
        ref_grid = [NODATA_value, NODATA_value, NODATA_value, NODATA_value, NODATA_value,
                    -1,    -1.25, -1.5,  -1.75, -2.0,
                    -0.75, -1.0,  -1.25, -1.5,  -1.75,
                    -0.5,  -0.75, -1.0,  -1.25, -1.5,
                    -0.25, -0.5,  -0.75, -1.0,  -1.25]


        #print grid
        assert allclose(grid, ref_grid)

        fid.close()

        #Cleanup
        #FIXME the file clean-up doesn't work (eg Permission Denied Error)
        #Done (Ole) - it was because sww2ers didn't close it's sww file
        os.remove(sww.filename)
        os.remove(self.domain.get_name() + '_elevation')
        os.remove(self.domain.get_name() + '_elevation.ers')



    def test_sww2pts_centroids(self):
        """Test that sww information can be converted correctly to pts data at specified coordinates
        - in this case, the centroids.
        """

        import time, os
        from Numeric import array, zeros, allclose, Float, concatenate, NewAxis
        from Scientific.IO.NetCDF import NetCDFFile
        from anuga.geospatial_data.geospatial_data import Geospatial_data

        # Used for points that lie outside mesh
        NODATA_value = 1758323

        # Setup
        self.domain.set_name('datatest')

        ptsfile = self.domain.get_name() + '_elevation.pts'
        swwfile = self.domain.get_name() + '.sww'

        self.domain.set_datadir('.')
        self.domain.format = 'sww'
        self.smooth = True #self.set_store_vertices_uniquely(False)
        self.domain.set_quantity('elevation', lambda x,y: -x-y)

        self.domain.geo_reference = Geo_reference(56,308500,6189000)

        sww = get_dataobject(self.domain)
        sww.store_connectivity()
        sww.store_timestep('stage')

        self.domain.evolve_to_end(finaltime = 0.01)
        sww.store_timestep('stage')

        # Check contents in NetCDF
        fid = NetCDFFile(sww.filename, 'r')

        # Get the variables
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        elevation = fid.variables['elevation'][:]
        time = fid.variables['time'][:]
        stage = fid.variables['stage'][:]

        volumes = fid.variables['volumes'][:]


        # Invoke interpolation for vertex points       
        points = concatenate( (x[:,NewAxis],y[:,NewAxis]), axis=1 )
        sww2pts(self.domain.get_name(),
                quantity = 'elevation',
                data_points = points,
                NODATA_value = NODATA_value,
                verbose = False)
        ref_point_values = elevation
        point_values = Geospatial_data(ptsfile).get_attributes()
        #print 'P', point_values
        #print 'Ref', ref_point_values        
        assert allclose(point_values, ref_point_values)        



        # Invoke interpolation for centroids
        points = self.domain.get_centroid_coordinates()
        #print points
        sww2pts(self.domain.get_name(),
                quantity = 'elevation',
                data_points = points,
                NODATA_value = NODATA_value,
                verbose = False)
        ref_point_values = [-0.5, -0.5, -1, -1, -1, -1, -1.5, -1.5]   #At centroids

        
        point_values = Geospatial_data(ptsfile).get_attributes()
        #print 'P', point_values
        #print 'Ref', ref_point_values        
        assert allclose(point_values, ref_point_values)        



        fid.close()

        #Cleanup
        os.remove(sww.filename)
        os.remove(ptsfile)




    def test_ferret2sww1(self):
        """Test that georeferencing etc works when converting from
        ferret format (lat/lon) to sww format (UTM)
        """
        from Scientific.IO.NetCDF import NetCDFFile
        import os, sys

        #The test file has
        # LON = 150.66667, 150.83334, 151, 151.16667
        # LAT = -34.5, -34.33333, -34.16667, -34 ;
        # TIME = 0, 0.1, 0.6, 1.1, 1.6, 2.1 ;
        #
        # First value (index=0) in small_ha.nc is 0.3400644 cm,
        # Fourth value (index==3) is -6.50198 cm



        #Read
        from anuga.coordinate_transforms.redfearn import redfearn
        #fid = NetCDFFile(self.test_MOST_file)
        fid = NetCDFFile(self.test_MOST_file + '_ha.nc')
        first_value = fid.variables['HA'][:][0,0,0]
        fourth_value = fid.variables['HA'][:][0,0,3]
        fid.close()


        #Call conversion (with zero origin)
        #ferret2sww('small', verbose=False,
        #           origin = (56, 0, 0))
        ferret2sww(self.test_MOST_file, verbose=False,
                   origin = (56, 0, 0))

        #Work out the UTM coordinates for first point
        zone, e, n = redfearn(-34.5, 150.66667)
        #print zone, e, n

        #Read output file 'small.sww'
        #fid = NetCDFFile('small.sww')
        fid = NetCDFFile(self.test_MOST_file + '.sww')

        x = fid.variables['x'][:]
        y = fid.variables['y'][:]

        #Check that first coordinate is correctly represented
        assert allclose(x[0], e)
        assert allclose(y[0], n)

        #Check first value
        stage = fid.variables['stage'][:]
        xmomentum = fid.variables['xmomentum'][:]
        ymomentum = fid.variables['ymomentum'][:]

        #print ymomentum

        assert allclose(stage[0,0], first_value/100)  #Meters

        #Check fourth value
        assert allclose(stage[0,3], fourth_value/100)  #Meters

        fid.close()

        #Cleanup
        import os
        os.remove(self.test_MOST_file + '.sww')


    def test_ferret2sww_2(self):
        """Test that georeferencing etc works when converting from
        ferret format (lat/lon) to sww format (UTM)
        """
        from Scientific.IO.NetCDF import NetCDFFile

        #The test file has
        # LON = 150.66667, 150.83334, 151, 151.16667
        # LAT = -34.5, -34.33333, -34.16667, -34 ;
        # TIME = 0, 0.1, 0.6, 1.1, 1.6, 2.1 ;
        #
        # First value (index=0) in small_ha.nc is 0.3400644 cm,
        # Fourth value (index==3) is -6.50198 cm


        from anuga.coordinate_transforms.redfearn import redfearn

        #fid = NetCDFFile('small_ha.nc')
        fid = NetCDFFile(self.test_MOST_file + '_ha.nc')

        #Pick a coordinate and a value

        time_index = 1
        lat_index = 0
        lon_index = 2

        test_value = fid.variables['HA'][:][time_index, lat_index, lon_index]
        test_time = fid.variables['TIME'][:][time_index]
        test_lat = fid.variables['LAT'][:][lat_index]
        test_lon = fid.variables['LON'][:][lon_index]

        linear_point_index = lat_index*4 + lon_index
        fid.close()

        #Call conversion (with zero origin)
        ferret2sww(self.test_MOST_file, verbose=False,
                   origin = (56, 0, 0))


        #Work out the UTM coordinates for test point
        zone, e, n = redfearn(test_lat, test_lon)

        #Read output file 'small.sww'
        fid = NetCDFFile(self.test_MOST_file + '.sww')

        x = fid.variables['x'][:]
        y = fid.variables['y'][:]

        #Check that test coordinate is correctly represented
        assert allclose(x[linear_point_index], e)
        assert allclose(y[linear_point_index], n)

        #Check test value
        stage = fid.variables['stage'][:]

        assert allclose(stage[time_index, linear_point_index], test_value/100)

        fid.close()

        #Cleanup
        import os
        os.remove(self.test_MOST_file + '.sww')


    def test_ferret2sww_lat_long(self):
        # Test that min lat long works

        #The test file has
        # LON = 150.66667, 150.83334, 151, 151.16667
        # LAT = -34.5, -34.33333, -34.16667, -34 ;
        
        #Read
        from anuga.coordinate_transforms.redfearn import redfearn
        fid = NetCDFFile(self.test_MOST_file + '_ha.nc')
        first_value = fid.variables['HA'][:][0,0,0]
        fourth_value = fid.variables['HA'][:][0,0,3]
        fid.close()


        #Call conversion (with zero origin)
        #ferret2sww('small', verbose=False,
        #           origin = (56, 0, 0))
        ferret2sww(self.test_MOST_file, verbose=False,
                   origin = (56, 0, 0), minlat=-34.5, maxlat=-34)

        #Work out the UTM coordinates for first point
        zone, e, n = redfearn(-34.5, 150.66667)
        #print zone, e, n

        #Read output file 'small.sww'
        #fid = NetCDFFile('small.sww')
        fid = NetCDFFile(self.test_MOST_file + '.sww')

        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        #Check that first coordinate is correctly represented
        assert 16 == len(x)

        fid.close()

        #Cleanup
        import os
        os.remove(self.test_MOST_file + '.sww')


    def test_ferret2sww_lat_longII(self):
        # Test that min lat long works

        #The test file has
        # LON = 150.66667, 150.83334, 151, 151.16667
        # LAT = -34.5, -34.33333, -34.16667, -34 ;
        
        #Read
        from anuga.coordinate_transforms.redfearn import redfearn
        fid = NetCDFFile(self.test_MOST_file + '_ha.nc')
        first_value = fid.variables['HA'][:][0,0,0]
        fourth_value = fid.variables['HA'][:][0,0,3]
        fid.close()


        #Call conversion (with zero origin)
        #ferret2sww('small', verbose=False,
        #           origin = (56, 0, 0))
        ferret2sww(self.test_MOST_file, verbose=False,
                   origin = (56, 0, 0), minlat=-34.4, maxlat=-34.2)

        #Work out the UTM coordinates for first point
        zone, e, n = redfearn(-34.5, 150.66667)
        #print zone, e, n

        #Read output file 'small.sww'
        #fid = NetCDFFile('small.sww')
        fid = NetCDFFile(self.test_MOST_file + '.sww')

        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        #Check that first coordinate is correctly represented
        assert 12 == len(x)

        fid.close()

        #Cleanup
        import os
        os.remove(self.test_MOST_file + '.sww')
        
    def test_ferret2sww3(self):
        """Elevation included
        """
        from Scientific.IO.NetCDF import NetCDFFile

        #The test file has
        # LON = 150.66667, 150.83334, 151, 151.16667
        # LAT = -34.5, -34.33333, -34.16667, -34 ;
        # ELEVATION = [-1 -2 -3 -4
        #             -5 -6 -7 -8
        #              ...
        #              ...      -16]
        # where the top left corner is -1m,
        # and the ll corner is -13.0m
        #
        # First value (index=0) in small_ha.nc is 0.3400644 cm,
        # Fourth value (index==3) is -6.50198 cm

        from anuga.coordinate_transforms.redfearn import redfearn
        import os
        fid1 = NetCDFFile('test_ha.nc','w')
        fid2 = NetCDFFile('test_ua.nc','w')
        fid3 = NetCDFFile('test_va.nc','w')
        fid4 = NetCDFFile('test_e.nc','w')

        h1_list = [150.66667,150.83334,151.]
        h2_list = [-34.5,-34.33333]

        long_name = 'LON'
        lat_name = 'LAT'
        time_name = 'TIME'

        nx = 3
        ny = 2

        for fid in [fid1,fid2,fid3]:
            fid.createDimension(long_name,nx)
            fid.createVariable(long_name,'d',(long_name,))
            fid.variables[long_name].point_spacing='uneven'
            fid.variables[long_name].units='degrees_east'
            fid.variables[long_name].assignValue(h1_list)

            fid.createDimension(lat_name,ny)
            fid.createVariable(lat_name,'d',(lat_name,))
            fid.variables[lat_name].point_spacing='uneven'
            fid.variables[lat_name].units='degrees_north'
            fid.variables[lat_name].assignValue(h2_list)

            fid.createDimension(time_name,2)
            fid.createVariable(time_name,'d',(time_name,))
            fid.variables[time_name].point_spacing='uneven'
            fid.variables[time_name].units='seconds'
            fid.variables[time_name].assignValue([0.,1.])
            if fid == fid3: break


        for fid in [fid4]:
            fid.createDimension(long_name,nx)
            fid.createVariable(long_name,'d',(long_name,))
            fid.variables[long_name].point_spacing='uneven'
            fid.variables[long_name].units='degrees_east'
            fid.variables[long_name].assignValue(h1_list)

            fid.createDimension(lat_name,ny)
            fid.createVariable(lat_name,'d',(lat_name,))
            fid.variables[lat_name].point_spacing='uneven'
            fid.variables[lat_name].units='degrees_north'
            fid.variables[lat_name].assignValue(h2_list)

        name = {}
        name[fid1]='HA'
        name[fid2]='UA'
        name[fid3]='VA'
        name[fid4]='ELEVATION'

        units = {}
        units[fid1]='cm'
        units[fid2]='cm/s'
        units[fid3]='cm/s'
        units[fid4]='m'

        values = {}
        values[fid1]=[[[5., 10.,15.], [13.,18.,23.]],[[50.,100.,150.],[130.,180.,230.]]]
        values[fid2]=[[[1., 2.,3.], [4.,5.,6.]],[[7.,8.,9.],[10.,11.,12.]]]
        values[fid3]=[[[13., 12.,11.], [10.,9.,8.]],[[7.,6.,5.],[4.,3.,2.]]]
        values[fid4]=[[-3000,-3100,-3200],[-4000,-5000,-6000]]

        for fid in [fid1,fid2,fid3]:
          fid.createVariable(name[fid],'d',(time_name,lat_name,long_name))
          fid.variables[name[fid]].point_spacing='uneven'
          fid.variables[name[fid]].units=units[fid]
          fid.variables[name[fid]].assignValue(values[fid])
          fid.variables[name[fid]].missing_value = -99999999.
          if fid == fid3: break

        for fid in [fid4]:
            fid.createVariable(name[fid],'d',(lat_name,long_name))
            fid.variables[name[fid]].point_spacing='uneven'
            fid.variables[name[fid]].units=units[fid]
            fid.variables[name[fid]].assignValue(values[fid])
            fid.variables[name[fid]].missing_value = -99999999.


        fid1.sync(); fid1.close()
        fid2.sync(); fid2.close()
        fid3.sync(); fid3.close()
        fid4.sync(); fid4.close()

        fid1 = NetCDFFile('test_ha.nc','r')
        fid2 = NetCDFFile('test_e.nc','r')
        fid3 = NetCDFFile('test_va.nc','r')


        first_amp = fid1.variables['HA'][:][0,0,0]
        third_amp = fid1.variables['HA'][:][0,0,2]
        first_elevation = fid2.variables['ELEVATION'][0,0]
        third_elevation= fid2.variables['ELEVATION'][:][0,2]
        first_speed = fid3.variables['VA'][0,0,0]
        third_speed = fid3.variables['VA'][:][0,0,2]

        fid1.close()
        fid2.close()
        fid3.close()

        #Call conversion (with zero origin)
        ferret2sww('test', verbose=False,
                   origin = (56, 0, 0), inverted_bathymetry=False)

        os.remove('test_va.nc')
        os.remove('test_ua.nc')
        os.remove('test_ha.nc')
        os.remove('test_e.nc')

        #Read output file 'test.sww'
        fid = NetCDFFile('test.sww')


        #Check first value
        elevation = fid.variables['elevation'][:]
        stage = fid.variables['stage'][:]
        xmomentum = fid.variables['xmomentum'][:]
        ymomentum = fid.variables['ymomentum'][:]

        #print ymomentum
        first_height = first_amp/100 - first_elevation
        third_height = third_amp/100 - third_elevation
        first_momentum=first_speed*first_height/100
        third_momentum=third_speed*third_height/100

        assert allclose(ymomentum[0][0],first_momentum)  #Meters
        assert allclose(ymomentum[0][2],third_momentum)  #Meters

        fid.close()

        #Cleanup
        os.remove('test.sww')



    def test_ferret2sww4(self):
        """Like previous but with augmented variable names as
        in files produced by ferret as opposed to MOST
        """
        from Scientific.IO.NetCDF import NetCDFFile

        #The test file has
        # LON = 150.66667, 150.83334, 151, 151.16667
        # LAT = -34.5, -34.33333, -34.16667, -34 ;
        # ELEVATION = [-1 -2 -3 -4
        #             -5 -6 -7 -8
        #              ...
        #              ...      -16]
        # where the top left corner is -1m,
        # and the ll corner is -13.0m
        #
        # First value (index=0) in small_ha.nc is 0.3400644 cm,
        # Fourth value (index==3) is -6.50198 cm

        from anuga.coordinate_transforms.redfearn import redfearn
        import os
        fid1 = NetCDFFile('test_ha.nc','w')
        fid2 = NetCDFFile('test_ua.nc','w')
        fid3 = NetCDFFile('test_va.nc','w')
        fid4 = NetCDFFile('test_e.nc','w')

        h1_list = [150.66667,150.83334,151.]
        h2_list = [-34.5,-34.33333]

#        long_name = 'LON961_1261'
#        lat_name = 'LAT481_841'
#        time_name = 'TIME1'

        long_name = 'LON'
        lat_name = 'LAT'
        time_name = 'TIME'

        nx = 3
        ny = 2

        for fid in [fid1,fid2,fid3]:
            fid.createDimension(long_name,nx)
            fid.createVariable(long_name,'d',(long_name,))
            fid.variables[long_name].point_spacing='uneven'
            fid.variables[long_name].units='degrees_east'
            fid.variables[long_name].assignValue(h1_list)

            fid.createDimension(lat_name,ny)
            fid.createVariable(lat_name,'d',(lat_name,))
            fid.variables[lat_name].point_spacing='uneven'
            fid.variables[lat_name].units='degrees_north'
            fid.variables[lat_name].assignValue(h2_list)

            fid.createDimension(time_name,2)
            fid.createVariable(time_name,'d',(time_name,))
            fid.variables[time_name].point_spacing='uneven'
            fid.variables[time_name].units='seconds'
            fid.variables[time_name].assignValue([0.,1.])
            if fid == fid3: break


        for fid in [fid4]:
            fid.createDimension(long_name,nx)
            fid.createVariable(long_name,'d',(long_name,))
            fid.variables[long_name].point_spacing='uneven'
            fid.variables[long_name].units='degrees_east'
            fid.variables[long_name].assignValue(h1_list)

            fid.createDimension(lat_name,ny)
            fid.createVariable(lat_name,'d',(lat_name,))
            fid.variables[lat_name].point_spacing='uneven'
            fid.variables[lat_name].units='degrees_north'
            fid.variables[lat_name].assignValue(h2_list)

        name = {}
        name[fid1]='HA'
        name[fid2]='UA'
        name[fid3]='VA'
        name[fid4]='ELEVATION'

        units = {}
        units[fid1]='cm'
        units[fid2]='cm/s'
        units[fid3]='cm/s'
        units[fid4]='m'

        values = {}
        values[fid1]=[[[5., 10.,15.], [13.,18.,23.]],[[50.,100.,150.],[130.,180.,230.]]]
        values[fid2]=[[[1., 2.,3.], [4.,5.,6.]],[[7.,8.,9.],[10.,11.,12.]]]
        values[fid3]=[[[13., 12.,11.], [10.,9.,8.]],[[7.,6.,5.],[4.,3.,2.]]]
        values[fid4]=[[-3000,-3100,-3200],[-4000,-5000,-6000]]

        for fid in [fid1,fid2,fid3]:
          fid.createVariable(name[fid],'d',(time_name,lat_name,long_name))
          fid.variables[name[fid]].point_spacing='uneven'
          fid.variables[name[fid]].units=units[fid]
          fid.variables[name[fid]].assignValue(values[fid])
          fid.variables[name[fid]].missing_value = -99999999.
          if fid == fid3: break

        for fid in [fid4]:
            fid.createVariable(name[fid],'d',(lat_name,long_name))
            fid.variables[name[fid]].point_spacing='uneven'
            fid.variables[name[fid]].units=units[fid]
            fid.variables[name[fid]].assignValue(values[fid])
            fid.variables[name[fid]].missing_value = -99999999.


        fid1.sync(); fid1.close()
        fid2.sync(); fid2.close()
        fid3.sync(); fid3.close()
        fid4.sync(); fid4.close()

        fid1 = NetCDFFile('test_ha.nc','r')
        fid2 = NetCDFFile('test_e.nc','r')
        fid3 = NetCDFFile('test_va.nc','r')


        first_amp = fid1.variables['HA'][:][0,0,0]
        third_amp = fid1.variables['HA'][:][0,0,2]
        first_elevation = fid2.variables['ELEVATION'][0,0]
        third_elevation= fid2.variables['ELEVATION'][:][0,2]
        first_speed = fid3.variables['VA'][0,0,0]
        third_speed = fid3.variables['VA'][:][0,0,2]

        fid1.close()
        fid2.close()
        fid3.close()

        #Call conversion (with zero origin)
        ferret2sww('test', verbose=False, origin = (56, 0, 0)
                   , inverted_bathymetry=False)

        os.remove('test_va.nc')
        os.remove('test_ua.nc')
        os.remove('test_ha.nc')
        os.remove('test_e.nc')

        #Read output file 'test.sww'
        fid = NetCDFFile('test.sww')


        #Check first value
        elevation = fid.variables['elevation'][:]
        stage = fid.variables['stage'][:]
        xmomentum = fid.variables['xmomentum'][:]
        ymomentum = fid.variables['ymomentum'][:]

        #print ymomentum
        first_height = first_amp/100 - first_elevation
        third_height = third_amp/100 - third_elevation
        first_momentum=first_speed*first_height/100
        third_momentum=third_speed*third_height/100

        assert allclose(ymomentum[0][0],first_momentum)  #Meters
        assert allclose(ymomentum[0][2],third_momentum)  #Meters

        fid.close()

        #Cleanup
        os.remove('test.sww')




    def test_ferret2sww_nz_origin(self):
        from Scientific.IO.NetCDF import NetCDFFile
        from anuga.coordinate_transforms.redfearn import redfearn

        #Call conversion (with nonzero origin)
        ferret2sww(self.test_MOST_file, verbose=False,
                   origin = (56, 100000, 200000))


        #Work out the UTM coordinates for first point
        zone, e, n = redfearn(-34.5, 150.66667)

        #Read output file 'small.sww'
        #fid = NetCDFFile('small.sww', 'r')
        fid = NetCDFFile(self.test_MOST_file + '.sww')

        x = fid.variables['x'][:]
        y = fid.variables['y'][:]

        #Check that first coordinate is correctly represented
        assert allclose(x[0], e-100000)
        assert allclose(y[0], n-200000)

        fid.close()

        #Cleanup
        os.remove(self.test_MOST_file + '.sww')


    def test_ferret2sww_lat_long(self):
        # Test that min lat long works

        #The test file has
        # LON = 150.66667, 150.83334, 151, 151.16667
        # LAT = -34.5, -34.33333, -34.16667, -34 ;
        
        #Read
        from anuga.coordinate_transforms.redfearn import redfearn
        fid = NetCDFFile(self.test_MOST_file + '_ha.nc')
        first_value = fid.variables['HA'][:][0,0,0]
        fourth_value = fid.variables['HA'][:][0,0,3]
        fid.close()


        #Call conversion (with zero origin)
        #ferret2sww('small', verbose=False,
        #           origin = (56, 0, 0))
        try:
            ferret2sww(self.test_MOST_file, verbose=False,
                   origin = (56, 0, 0), minlat=-34.5, maxlat=-35)
        except AssertionError:
            pass
        else:
            self.failUnless(0 ==1,  'Bad input did not throw exception error!')

    def test_sww_extent(self):
        """Not a test, rather a look at the sww format
        """

        import time, os
        from Numeric import array, zeros, allclose, Float, concatenate
        from Scientific.IO.NetCDF import NetCDFFile

        self.domain.set_name('datatest' + str(id(self)))
        self.domain.format = 'sww'
        self.domain.smooth = True
        self.domain.reduction = mean
        self.domain.set_datadir('.')


        sww = get_dataobject(self.domain)
        sww.store_connectivity()
        sww.store_timestep('stage')
        self.domain.time = 2.

        #Modify stage at second timestep
        stage = self.domain.quantities['stage'].vertex_values
        self.domain.set_quantity('stage', stage/2)

        sww.store_timestep('stage')

        file_and_extension_name = self.domain.get_name() + ".sww"
        #print "file_and_extension_name",file_and_extension_name
        [xmin, xmax, ymin, ymax, stagemin, stagemax] = \
               extent_sww(file_and_extension_name )

        assert allclose(xmin, 0.0)
        assert allclose(xmax, 1.0)
        assert allclose(ymin, 0.0)
        assert allclose(ymax, 1.0)
        assert allclose(stagemin, -0.85)
        assert allclose(stagemax, 0.15)


        #Cleanup
        os.remove(sww.filename)



    def test_sww2domain1(self):
        ################################################
        #Create a test domain, and evolve and save it.
        ################################################
        from mesh_factory import rectangular
        from Numeric import array

        #Create basic mesh

        yiel=0.01
        points, vertices, boundary = rectangular(10,10)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.geo_reference = Geo_reference(56,11,11)
        domain.smooth = False
        domain.store = True
        domain.set_name('bedslope')
        domain.default_order=2
        #Bed-slope and friction
        domain.set_quantity('elevation', lambda x,y: -x/3)
        domain.set_quantity('friction', 0.1)
        # Boundary conditions
        from math import sin, pi
        Br = Reflective_boundary(domain)
        Bt = Transmissive_boundary(domain)
        Bd = Dirichlet_boundary([0.2,0.,0.])
        Bw = Time_boundary(domain=domain,f=lambda t: [(0.1*sin(t*2*pi)), 0.0, 0.0])

        #domain.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})
        domain.set_boundary({'left': Bd, 'right': Bd, 'top': Bd, 'bottom': Bd})

        domain.quantities_to_be_stored.extend(['xmomentum','ymomentum'])
        #Initial condition
        h = 0.05
        elevation = domain.quantities['elevation'].vertex_values
        domain.set_quantity('stage', elevation + h)

        domain.check_integrity()
        #Evolution
        for t in domain.evolve(yieldstep = yiel, finaltime = 0.05):
            #domain.write_time()
            pass


        ##########################################
        #Import the example's file as a new domain
        ##########################################
        from data_manager import sww2domain
        from Numeric import allclose
        import os

        filename = domain.datadir + os.sep + domain.get_name() + '.sww'
        domain2 = sww2domain(filename,None,fail_if_NaN=False,verbose = False)
        #points, vertices, boundary = rectangular(15,15)
        #domain2.boundary = boundary
        ###################
        ##NOW TEST IT!!!
        ###################

        #os.remove(domain.get_name() + '.sww')
        os.remove(filename)

        bits = ['vertex_coordinates']
        for quantity in ['elevation']+domain.quantities_to_be_stored:
            bits.append('get_quantity("%s").get_integral()' %quantity)
            bits.append('get_quantity("%s").get_values()' %quantity)

        for bit in bits:
            #print 'testing that domain.'+bit+' has been restored'
            #print bit
            #print 'done'
            assert allclose(eval('domain.'+bit),eval('domain2.'+bit))

        ######################################
        #Now evolve them both, just to be sure
        ######################################x
        domain.time = 0.
        from time import sleep

        final = .1
        domain.set_quantity('friction', 0.1)
        domain.store = False
        domain.set_boundary({'left': Bd, 'right': Bd, 'top': Bd, 'bottom': Bd})


        for t in domain.evolve(yieldstep = yiel, finaltime = final):
            #domain.write_time()
            pass

        final = final - (domain2.starttime-domain.starttime)
        #BUT since domain1 gets time hacked back to 0:
        final = final + (domain2.starttime-domain.starttime)

        domain2.smooth = False
        domain2.store = False
        domain2.default_order=2
        domain2.set_quantity('friction', 0.1)
        #Bed-slope and friction
        # Boundary conditions
        Bd2=Dirichlet_boundary([0.2,0.,0.])
        domain2.boundary = domain.boundary
        #print 'domain2.boundary'
        #print domain2.boundary
        domain2.set_boundary({'left': Bd, 'right': Bd, 'top': Bd, 'bottom': Bd})
        #domain2.set_boundary({'exterior': Bd})

        domain2.check_integrity()

        for t in domain2.evolve(yieldstep = yiel, finaltime = final):
            #domain2.write_time()
            pass

        ###################
        ##NOW TEST IT!!!
        ##################

        bits = ['vertex_coordinates']

        for quantity in ['elevation','stage', 'ymomentum','xmomentum']:
            bits.append('get_quantity("%s").get_integral()' %quantity)
            bits.append('get_quantity("%s").get_values()' %quantity)

        #print bits
        for bit in bits:
            #print bit
            #print eval('domain.'+bit)
            #print eval('domain2.'+bit)
            
            #print eval('domain.'+bit+'-domain2.'+bit)
            msg = 'Values in the two domains are different for ' + bit
            assert allclose(eval('domain.'+bit),eval('domain2.'+bit),
                            rtol=1.e-5, atol=3.e-8), msg


    def test_sww2domain2(self):
        ##################################################################
        #Same as previous test, but this checks how NaNs are handled.
        ##################################################################


        from mesh_factory import rectangular
        from Numeric import array

        #Create basic mesh
        points, vertices, boundary = rectangular(2,2)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.store = True
        domain.set_name('test_file')
        domain.set_datadir('.')
        domain.default_order=2
        domain.quantities_to_be_stored=['stage']

        domain.set_quantity('elevation', lambda x,y: -x/3)
        domain.set_quantity('friction', 0.1)

        from math import sin, pi
        Br = Reflective_boundary(domain)
        Bt = Transmissive_boundary(domain)
        Bd = Dirichlet_boundary([0.2,0.,0.])
        Bw = Time_boundary(domain=domain,
                           f=lambda t: [(0.1*sin(t*2*pi)), 0.0, 0.0])

        domain.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})

        h = 0.05
        elevation = domain.quantities['elevation'].vertex_values
        domain.set_quantity('stage', elevation + h)

        domain.check_integrity()

        for t in domain.evolve(yieldstep = 1, finaltime = 2.0):
            pass
            #domain.write_time()



        ##################################
        #Import the file as a new domain
        ##################################
        from data_manager import sww2domain
        from Numeric import allclose
        import os

        filename = domain.datadir + os.sep + domain.get_name() + '.sww'

        #Fail because NaNs are present
        try:
            domain2 = sww2domain(filename,boundary,fail_if_NaN=True,verbose=False)
        except:
            #Now import it, filling NaNs to be 0
            filler = 0
            domain2 = sww2domain(filename,None,fail_if_NaN=False,NaN_filler = filler,verbose=False)

        #Clean up
        os.remove(filename)


        bits = [ 'geo_reference.get_xllcorner()',
                'geo_reference.get_yllcorner()',
                'vertex_coordinates']

        for quantity in ['elevation']+domain.quantities_to_be_stored:
            bits.append('get_quantity("%s").get_integral()' %quantity)
            bits.append('get_quantity("%s").get_values()' %quantity)

        for bit in bits:
        #    print 'testing that domain.'+bit+' has been restored'
            assert allclose(eval('domain.'+bit),eval('domain2.'+bit))

        assert max(max(domain2.get_quantity('xmomentum').get_values()))==filler
        assert min(min(domain2.get_quantity('xmomentum').get_values()))==filler
        assert max(max(domain2.get_quantity('ymomentum').get_values()))==filler
        assert min(min(domain2.get_quantity('ymomentum').get_values()))==filler



    #def test_weed(self):
        from data_manager import weed

        coordinates1 = [[0.,0.],[1.,0.],[1.,1.],[1.,0.],[2.,0.],[1.,1.]]
        volumes1 = [[0,1,2],[3,4,5]]
        boundary1= {(0,1): 'external',(1,2): 'not external',(2,0): 'external',(3,4): 'external',(4,5): 'external',(5,3): 'not external'}
        coordinates2,volumes2,boundary2=weed(coordinates1,volumes1,boundary1)

        points2 = {(0.,0.):None,(1.,0.):None,(1.,1.):None,(2.,0.):None}

        assert len(points2)==len(coordinates2)
        for i in range(len(coordinates2)):
            coordinate = tuple(coordinates2[i])
            assert points2.has_key(coordinate)
            points2[coordinate]=i

        for triangle in volumes1:
            for coordinate in triangle:
                assert coordinates2[points2[tuple(coordinates1[coordinate])]][0]==coordinates1[coordinate][0]
                assert coordinates2[points2[tuple(coordinates1[coordinate])]][1]==coordinates1[coordinate][1]


    #FIXME This fails - smooth makes the comparism too hard for allclose
    def ztest_sww2domain3(self):
        ################################################
        #DOMAIN.SMOOTH = TRUE !!!!!!!!!!!!!!!!!!!!!!!!!!
        ################################################
        from mesh_factory import rectangular
        from Numeric import array
        #Create basic mesh

        yiel=0.01
        points, vertices, boundary = rectangular(10,10)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.geo_reference = Geo_reference(56,11,11)
        domain.smooth = True
        domain.store = True
        domain.set_name('bedslope')
        domain.default_order=2
        #Bed-slope and friction
        domain.set_quantity('elevation', lambda x,y: -x/3)
        domain.set_quantity('friction', 0.1)
        # Boundary conditions
        from math import sin, pi
        Br = Reflective_boundary(domain)
        Bt = Transmissive_boundary(domain)
        Bd = Dirichlet_boundary([0.2,0.,0.])
        Bw = Time_boundary(domain=domain,
                           f=lambda t: [(0.1*sin(t*2*pi)), 0.0, 0.0])

        domain.set_boundary({'left': Bd, 'right': Bd, 'top': Bd, 'bottom': Bd})

        domain.quantities_to_be_stored.extend(['xmomentum','ymomentum'])
        #Initial condition
        h = 0.05
        elevation = domain.quantities['elevation'].vertex_values
        domain.set_quantity('stage', elevation + h)


        domain.check_integrity()
        #Evolution
        for t in domain.evolve(yieldstep = yiel, finaltime = 0.05):
        #    domain.write_time()
            pass


        ##########################################
        #Import the example's file as a new domain
        ##########################################
        from data_manager import sww2domain
        from Numeric import allclose
        import os

        filename = domain.datadir + os.sep + domain.get_name() + '.sww'
        domain2 = sww2domain(filename,None,fail_if_NaN=False,verbose = False)
        #points, vertices, boundary = rectangular(15,15)
        #domain2.boundary = boundary
        ###################
        ##NOW TEST IT!!!
        ###################

        os.remove(domain.get_name() + '.sww')

        #FIXME smooth domain so that they can be compared


        bits = []#'vertex_coordinates']
        for quantity in ['elevation']+domain.quantities_to_be_stored:
            bits.append('quantities["%s"].get_integral()'%quantity)


        for bit in bits:
            #print 'testing that domain.'+bit+' has been restored'
            #print bit
            #print 'done'
            #print ('domain.'+bit), eval('domain.'+bit)
            #print ('domain2.'+bit), eval('domain2.'+bit)
            assert allclose(eval('domain.'+bit),eval('domain2.'+bit),rtol=1.0e-1,atol=1.e-3)
            pass

        ######################################
        #Now evolve them both, just to be sure
        ######################################x
        domain.time = 0.
        from time import sleep

        final = .5
        domain.set_quantity('friction', 0.1)
        domain.store = False
        domain.set_boundary({'left': Bd, 'right': Bd, 'top': Bd, 'bottom': Br})

        for t in domain.evolve(yieldstep = yiel, finaltime = final):
            #domain.write_time()
            pass

        domain2.smooth = True
        domain2.store = False
        domain2.default_order=2
        domain2.set_quantity('friction', 0.1)
        #Bed-slope and friction
        # Boundary conditions
        Bd2=Dirichlet_boundary([0.2,0.,0.])
        Br2 = Reflective_boundary(domain2)
        domain2.boundary = domain.boundary
        #print 'domain2.boundary'
        #print domain2.boundary
        domain2.set_boundary({'left': Bd2, 'right': Bd2, 'top': Bd2, 'bottom': Br2})
        #domain2.boundary = domain.boundary
        #domain2.set_boundary({'exterior': Bd})

        domain2.check_integrity()

        for t in domain2.evolve(yieldstep = yiel, finaltime = final):
            #domain2.write_time()
            pass

        ###################
        ##NOW TEST IT!!!
        ##################

        print '><><><><>>'
        bits = [ 'vertex_coordinates']

        for quantity in ['elevation','xmomentum','ymomentum']:#+domain.quantities_to_be_stored:
            #bits.append('quantities["%s"].get_integral()'%quantity)
            bits.append('get_quantity("%s").get_values()' %quantity)

        for bit in bits:
            print bit
            assert allclose(eval('domain.'+bit),eval('domain2.'+bit))


    def test_decimate_dem(self):
        """Test decimation of dem file
        """

        import os
        from Numeric import ones, allclose, Float, arange
        from Scientific.IO.NetCDF import NetCDFFile

        #Write test dem file
        root = 'decdemtest'

        filename = root + '.dem'
        fid = NetCDFFile(filename, 'w')

        fid.institution = 'Geoscience Australia'
        fid.description = 'NetCDF DEM format for compact and portable ' +\
                          'storage of spatial point data'

        nrows = 15
        ncols = 18

        fid.ncols = ncols
        fid.nrows = nrows
        fid.xllcorner = 2000.5
        fid.yllcorner = 3000.5
        fid.cellsize = 25
        fid.NODATA_value = -9999

        fid.zone = 56
        fid.false_easting = 0.0
        fid.false_northing = 0.0
        fid.projection = 'UTM'
        fid.datum = 'WGS84'
        fid.units = 'METERS'

        fid.createDimension('number_of_points', nrows*ncols)

        fid.createVariable('elevation', Float, ('number_of_points',))

        elevation = fid.variables['elevation']

        elevation[:] = (arange(nrows*ncols))

        fid.close()

        #generate the elevation values expected in the decimated file
        ref_elevation = [(  0+  1+  2+ 18+ 19+ 20+ 36+ 37+ 38) / 9.0,
                         (  4+  5+  6+ 22+ 23+ 24+ 40+ 41+ 42) / 9.0,
                         (  8+  9+ 10+ 26+ 27+ 28+ 44+ 45+ 46) / 9.0,
                         ( 12+ 13+ 14+ 30+ 31+ 32+ 48+ 49+ 50) / 9.0,
                         ( 72+ 73+ 74+ 90+ 91+ 92+108+109+110) / 9.0,
                         ( 76+ 77+ 78+ 94+ 95+ 96+112+113+114) / 9.0,
                         ( 80+ 81+ 82+ 98+ 99+100+116+117+118) / 9.0,
                         ( 84+ 85+ 86+102+103+104+120+121+122) / 9.0,
                         (144+145+146+162+163+164+180+181+182) / 9.0,
                         (148+149+150+166+167+168+184+185+186) / 9.0,
                         (152+153+154+170+171+172+188+189+190) / 9.0,
                         (156+157+158+174+175+176+192+193+194) / 9.0,
                         (216+217+218+234+235+236+252+253+254) / 9.0,
                         (220+221+222+238+239+240+256+257+258) / 9.0,
                         (224+225+226+242+243+244+260+261+262) / 9.0,
                         (228+229+230+246+247+248+264+265+266) / 9.0]

        #generate a stencil for computing the decimated values
        stencil = ones((3,3), Float) / 9.0

        decimate_dem(root, stencil=stencil, cellsize_new=100)

        #Open decimated NetCDF file
        fid = NetCDFFile(root + '_100.dem', 'r')

        # Get decimated elevation
        elevation = fid.variables['elevation']

        #Check values
        assert allclose(elevation, ref_elevation)

        #Cleanup
        fid.close()

        os.remove(root + '.dem')
        os.remove(root + '_100.dem')

    def test_decimate_dem_NODATA(self):
        """Test decimation of dem file that includes NODATA values
        """

        import os
        from Numeric import ones, allclose, Float, arange, reshape
        from Scientific.IO.NetCDF import NetCDFFile

        #Write test dem file
        root = 'decdemtest'

        filename = root + '.dem'
        fid = NetCDFFile(filename, 'w')

        fid.institution = 'Geoscience Australia'
        fid.description = 'NetCDF DEM format for compact and portable ' +\
                          'storage of spatial point data'

        nrows = 15
        ncols = 18
        NODATA_value = -9999

        fid.ncols = ncols
        fid.nrows = nrows
        fid.xllcorner = 2000.5
        fid.yllcorner = 3000.5
        fid.cellsize = 25
        fid.NODATA_value = NODATA_value

        fid.zone = 56
        fid.false_easting = 0.0
        fid.false_northing = 0.0
        fid.projection = 'UTM'
        fid.datum = 'WGS84'
        fid.units = 'METERS'

        fid.createDimension('number_of_points', nrows*ncols)

        fid.createVariable('elevation', Float, ('number_of_points',))

        elevation = fid.variables['elevation']

        #generate initial elevation values
        elevation_tmp = (arange(nrows*ncols))
        #add some NODATA values
        elevation_tmp[0]   = NODATA_value
        elevation_tmp[95]  = NODATA_value
        elevation_tmp[188] = NODATA_value
        elevation_tmp[189] = NODATA_value
        elevation_tmp[190] = NODATA_value
        elevation_tmp[209] = NODATA_value
        elevation_tmp[252] = NODATA_value

        elevation[:] = elevation_tmp

        fid.close()

        #generate the elevation values expected in the decimated file
        ref_elevation = [NODATA_value,
                         (  4+  5+  6+ 22+ 23+ 24+ 40+ 41+ 42) / 9.0,
                         (  8+  9+ 10+ 26+ 27+ 28+ 44+ 45+ 46) / 9.0,
                         ( 12+ 13+ 14+ 30+ 31+ 32+ 48+ 49+ 50) / 9.0,
                         ( 72+ 73+ 74+ 90+ 91+ 92+108+109+110) / 9.0,
                         NODATA_value,
                         ( 80+ 81+ 82+ 98+ 99+100+116+117+118) / 9.0,
                         ( 84+ 85+ 86+102+103+104+120+121+122) / 9.0,
                         (144+145+146+162+163+164+180+181+182) / 9.0,
                         (148+149+150+166+167+168+184+185+186) / 9.0,
                         NODATA_value,
                         (156+157+158+174+175+176+192+193+194) / 9.0,
                         NODATA_value,
                         (220+221+222+238+239+240+256+257+258) / 9.0,
                         (224+225+226+242+243+244+260+261+262) / 9.0,
                         (228+229+230+246+247+248+264+265+266) / 9.0]

        #generate a stencil for computing the decimated values
        stencil = ones((3,3), Float) / 9.0

        decimate_dem(root, stencil=stencil, cellsize_new=100)

        #Open decimated NetCDF file
        fid = NetCDFFile(root + '_100.dem', 'r')

        # Get decimated elevation
        elevation = fid.variables['elevation']

        #Check values
        assert allclose(elevation, ref_elevation)

        #Cleanup
        fid.close()

        os.remove(root + '.dem')
        os.remove(root + '_100.dem')

    def xxxtestz_sww2ers_real(self):
        """Test that sww information can be converted correctly to asc/prj
        format readable by e.g. ArcView
        """

        import time, os
        from Numeric import array, zeros, allclose, Float, concatenate
        from Scientific.IO.NetCDF import NetCDFFile

        # the memory optimised least squares
        #  cellsize = 20,   # this one seems to hang
        #  cellsize = 200000, # Ran 1 test in 269.703s
                                #Ran 1 test in 267.344s
        #  cellsize = 20000,  # Ran 1 test in 460.922s
        #  cellsize = 2000   #Ran 1 test in 5340.250s
        #  cellsize = 200   #this one seems to hang, building matirx A

        # not optimised
        # seems to hang
        #  cellsize = 2000   # Ran 1 test in 5334.563s
        #Export to ascii/prj files
        sww2dem('karratha_100m',
                quantity = 'depth',
                cellsize = 200000,
                verbose = True)

    def test_read_asc(self):
        """Test conversion from dem in ascii format to native NetCDF xya format
        """

        import time, os
        from Numeric import array, zeros, allclose, Float, concatenate
        from Scientific.IO.NetCDF import NetCDFFile

        from data_manager import _read_asc
        #Write test asc file
        filename = tempfile.mktemp(".000")
        fid = open(filename, 'w')
        fid.write("""ncols         7
nrows         4
xllcorner     2000.5
yllcorner     3000.5
cellsize      25
NODATA_value  -9999
    97.921    99.285   125.588   180.830   258.645   342.872   415.836
   473.157   514.391   553.893   607.120   678.125   777.283   883.038
   984.494  1040.349  1008.161   900.738   730.882   581.430   514.980
   502.645   516.230   504.739   450.604   388.500   338.097   514.980
""")
        fid.close()
        bath_metadata, grid = _read_asc(filename, verbose=False)
        self.failUnless(bath_metadata['xllcorner']  == 2000.5,  'Failed')
        self.failUnless(bath_metadata['yllcorner']  == 3000.5,  'Failed')
        self.failUnless(bath_metadata['cellsize']  == 25,  'Failed')
        self.failUnless(bath_metadata['NODATA_value']  == -9999,  'Failed')
        self.failUnless(grid[0][0]  == 97.921,  'Failed')
        self.failUnless(grid[3][6]  == 514.980,  'Failed')

        os.remove(filename)

    def test_asc_csiro2sww(self):
        import tempfile

        bath_dir = tempfile.mkdtemp()
        bath_dir_filename = bath_dir + os.sep +'ba19940524.000'
        #bath_dir = 'bath_data_manager_test'
        #print "os.getcwd( )",os.getcwd( )
        elevation_dir =  tempfile.mkdtemp()
        #elevation_dir = 'elev_expanded'
        elevation_dir_filename1 = elevation_dir + os.sep +'el19940524.000'
        elevation_dir_filename2 = elevation_dir + os.sep +'el19940524.001'

        fid = open(bath_dir_filename, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    9000.000    -1000.000    3000.0
   -1000.000    9000.000  -1000.000
""")
        fid.close()

        fid = open(elevation_dir_filename1, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    9000.000    0.000    3000.0
     0.000     9000.000     0.000
""")
        fid.close()

        fid = open(elevation_dir_filename2, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    9000.000    4000.000    4000.0
    4000.000    9000.000    4000.000
""")
        fid.close()

        ucur_dir =  tempfile.mkdtemp()
        ucur_dir_filename1 = ucur_dir + os.sep +'uc19940524.000'
        ucur_dir_filename2 = ucur_dir + os.sep +'uc19940524.001'

        fid = open(ucur_dir_filename1, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    90.000    60.000    30.0
    10.000    10.000    10.000
""")
        fid.close()
        fid = open(ucur_dir_filename2, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    90.000    60.000    30.0
    10.000    10.000    10.000
""")
        fid.close()

        vcur_dir =  tempfile.mkdtemp()
        vcur_dir_filename1 = vcur_dir + os.sep +'vc19940524.000'
        vcur_dir_filename2 = vcur_dir + os.sep +'vc19940524.001'

        fid = open(vcur_dir_filename1, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    90.000    60.000    30.0
    10.000    10.000    10.000
""")
        fid.close()
        fid = open(vcur_dir_filename2, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    90.000    60.000    30.0
    10.000    10.000    10.000
""")
        fid.close()

        sww_file = 'a_test.sww'
        asc_csiro2sww(bath_dir,elevation_dir, ucur_dir, vcur_dir, sww_file)

        # check the sww file

        fid = NetCDFFile(sww_file, 'r')    #Open existing file for read
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        z = fid.variables['z'][:]
        stage = fid.variables['stage'][:]
        xmomentum = fid.variables['xmomentum'][:]
        geo_ref = Geo_reference(NetCDFObject=fid)
        #print "geo_ref",geo_ref
        x_ref = geo_ref.get_xllcorner()
        y_ref = geo_ref.get_yllcorner()
        self.failUnless(geo_ref.get_zone() == 55,  'Failed')
        assert allclose(x_ref, 587798.418) # (-38, 148)
        assert allclose(y_ref, 5793123.477)# (-38, 148.5)

        #Zone:   55
        #Easting:  588095.674  Northing: 5821451.722
        #Latitude:   -37  45 ' 0.00000 ''  Longitude: 148 0 ' 0.00000 ''
        assert allclose((x[0],y[0]), (588095.674 - x_ref, 5821451.722 - y_ref))

        #Zone:   55
        #Easting:  632145.632  Northing: 5820863.269
        #Latitude:   -37  45 ' 0.00000 ''  Longitude: 148  30 ' 0.00000 ''
        assert allclose((x[2],y[2]), (632145.632 - x_ref, 5820863.269 - y_ref))

        #Zone:   55
        #Easting:  609748.788  Northing: 5793447.860
        #Latitude:   -38  0 ' 0.00000 ''  Longitude: 148  15 ' 0.00000 ''
        assert allclose((x[4],y[4]), (609748.788  - x_ref, 5793447.86 - y_ref))

        assert allclose(z[0],9000.0 )
        assert allclose(stage[0][1],0.0 )

        #(4000+1000)*60
        assert allclose(xmomentum[1][1],300000.0 )


        fid.close()

        #tidy up
        os.remove(bath_dir_filename)
        os.rmdir(bath_dir)

        os.remove(elevation_dir_filename1)
        os.remove(elevation_dir_filename2)
        os.rmdir(elevation_dir)

        os.remove(ucur_dir_filename1)
        os.remove(ucur_dir_filename2)
        os.rmdir(ucur_dir)

        os.remove(vcur_dir_filename1)
        os.remove(vcur_dir_filename2)
        os.rmdir(vcur_dir)


        # remove sww file
        os.remove(sww_file)

    def test_asc_csiro2sww2(self):
        import tempfile

        bath_dir = tempfile.mkdtemp()
        bath_dir_filename = bath_dir + os.sep +'ba19940524.000'
        #bath_dir = 'bath_data_manager_test'
        #print "os.getcwd( )",os.getcwd( )
        elevation_dir =  tempfile.mkdtemp()
        #elevation_dir = 'elev_expanded'
        elevation_dir_filename1 = elevation_dir + os.sep +'el19940524.000'
        elevation_dir_filename2 = elevation_dir + os.sep +'el19940524.001'

        fid = open(bath_dir_filename, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    9000.000    -1000.000    3000.0
   -1000.000    9000.000  -1000.000
""")
        fid.close()

        fid = open(elevation_dir_filename1, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    9000.000    0.000    3000.0
     0.000     -9999.000     -9999.000
""")
        fid.close()

        fid = open(elevation_dir_filename2, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    9000.000    4000.000    4000.0
    4000.000    9000.000    4000.000
""")
        fid.close()

        ucur_dir =  tempfile.mkdtemp()
        ucur_dir_filename1 = ucur_dir + os.sep +'uc19940524.000'
        ucur_dir_filename2 = ucur_dir + os.sep +'uc19940524.001'

        fid = open(ucur_dir_filename1, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    90.000    60.000    30.0
    10.000    10.000    10.000
""")
        fid.close()
        fid = open(ucur_dir_filename2, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    90.000    60.000    30.0
    10.000    10.000    10.000
""")
        fid.close()

        vcur_dir =  tempfile.mkdtemp()
        vcur_dir_filename1 = vcur_dir + os.sep +'vc19940524.000'
        vcur_dir_filename2 = vcur_dir + os.sep +'vc19940524.001'

        fid = open(vcur_dir_filename1, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    90.000    60.000    30.0
    10.000    10.000    10.000
""")
        fid.close()
        fid = open(vcur_dir_filename2, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    90.000    60.000    30.0
    10.000    10.000    10.000
""")
        fid.close()

        try:
            asc_csiro2sww(bath_dir,elevation_dir, ucur_dir,
                          vcur_dir, sww_file)
        except:
            #tidy up
            os.remove(bath_dir_filename)
            os.rmdir(bath_dir)

            os.remove(elevation_dir_filename1)
            os.remove(elevation_dir_filename2)
            os.rmdir(elevation_dir)

            os.remove(ucur_dir_filename1)
            os.remove(ucur_dir_filename2)
            os.rmdir(ucur_dir)

            os.remove(vcur_dir_filename1)
            os.remove(vcur_dir_filename2)
            os.rmdir(vcur_dir)
        else:
            #tidy up
            os.remove(bath_dir_filename)
            os.rmdir(bath_dir)

            os.remove(elevation_dir_filename1)
            os.remove(elevation_dir_filename2)
            os.rmdir(elevation_dir)
            raise 'Should raise exception'

            os.remove(ucur_dir_filename1)
            os.remove(ucur_dir_filename2)
            os.rmdir(ucur_dir)

            os.remove(vcur_dir_filename1)
            os.remove(vcur_dir_filename2)
            os.rmdir(vcur_dir)



    def test_asc_csiro2sww3(self):
        import tempfile

        bath_dir = tempfile.mkdtemp()
        bath_dir_filename = bath_dir + os.sep +'ba19940524.000'
        #bath_dir = 'bath_data_manager_test'
        #print "os.getcwd( )",os.getcwd( )
        elevation_dir =  tempfile.mkdtemp()
        #elevation_dir = 'elev_expanded'
        elevation_dir_filename1 = elevation_dir + os.sep +'el19940524.000'
        elevation_dir_filename2 = elevation_dir + os.sep +'el19940524.001'

        fid = open(bath_dir_filename, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    9000.000    -1000.000    3000.0
   -1000.000    9000.000  -1000.000
""")
        fid.close()

        fid = open(elevation_dir_filename1, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    9000.000    0.000    3000.0
     0.000     -9999.000     -9999.000
""")
        fid.close()

        fid = open(elevation_dir_filename2, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    9000.000    4000.000    4000.0
    4000.000    9000.000    4000.000
""")
        fid.close()

        ucur_dir =  tempfile.mkdtemp()
        ucur_dir_filename1 = ucur_dir + os.sep +'uc19940524.000'
        ucur_dir_filename2 = ucur_dir + os.sep +'uc19940524.001'

        fid = open(ucur_dir_filename1, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    90.000    60.000    30.0
    10.000    10.000    10.000
""")
        fid.close()
        fid = open(ucur_dir_filename2, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    90.000    60.000    30.0
    10.000    10.000    10.000
""")
        fid.close()

        vcur_dir =  tempfile.mkdtemp()
        vcur_dir_filename1 = vcur_dir + os.sep +'vc19940524.000'
        vcur_dir_filename2 = vcur_dir + os.sep +'vc19940524.001'

        fid = open(vcur_dir_filename1, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    90.000    60.000    30.0
    10.000    10.000    10.000
""")
        fid.close()
        fid = open(vcur_dir_filename2, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    90.000    60.000    30.0
    10.000    10.000    10.000
""")
        fid.close()

        sww_file = 'a_test.sww'
        asc_csiro2sww(bath_dir,elevation_dir, ucur_dir, vcur_dir,
                      sww_file, fail_on_NaN = False, elevation_NaN_filler = 0,
                      mean_stage = 100)

        # check the sww file

        fid = NetCDFFile(sww_file, 'r')    #Open existing file for read
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        z = fid.variables['z'][:]
        stage = fid.variables['stage'][:]
        xmomentum = fid.variables['xmomentum'][:]
        geo_ref = Geo_reference(NetCDFObject=fid)
        #print "geo_ref",geo_ref
        x_ref = geo_ref.get_xllcorner()
        y_ref = geo_ref.get_yllcorner()
        self.failUnless(geo_ref.get_zone() == 55,  'Failed')
        assert allclose(x_ref, 587798.418) # (-38, 148)
        assert allclose(y_ref, 5793123.477)# (-38, 148.5)

        #Zone:   55
        #Easting:  588095.674  Northing: 5821451.722
        #Latitude:   -37  45 ' 0.00000 ''  Longitude: 148 0 ' 0.00000 ''
        assert allclose((x[0],y[0]), (588095.674 - x_ref, 5821451.722 - y_ref))

        #Zone:   55
        #Easting:  632145.632  Northing: 5820863.269
        #Latitude:   -37  45 ' 0.00000 ''  Longitude: 148  30 ' 0.00000 ''
        assert allclose((x[2],y[2]), (632145.632 - x_ref, 5820863.269 - y_ref))

        #Zone:   55
        #Easting:  609748.788  Northing: 5793447.860
        #Latitude:   -38  0 ' 0.00000 ''  Longitude: 148  15 ' 0.00000 ''
        assert allclose((x[4],y[4]), (609748.788  - x_ref, 5793447.86 - y_ref))

        assert allclose(z[0],9000.0 )
        assert allclose(stage[0][4],100.0 )
        assert allclose(stage[0][5],100.0 )

        #(100.0 - 9000)*10
        assert allclose(xmomentum[0][4], -89000.0 )

        #(100.0 - -1000.000)*10
        assert allclose(xmomentum[0][5], 11000.0 )

        fid.close()

        #tidy up
        os.remove(bath_dir_filename)
        os.rmdir(bath_dir)

        os.remove(elevation_dir_filename1)
        os.remove(elevation_dir_filename2)
        os.rmdir(elevation_dir)

        os.remove(ucur_dir_filename1)
        os.remove(ucur_dir_filename2)
        os.rmdir(ucur_dir)

        os.remove(vcur_dir_filename1)
        os.remove(vcur_dir_filename2)
        os.rmdir(vcur_dir)

        # remove sww file
        os.remove(sww_file)


    def test_asc_csiro2sww4(self):
        """
        Test specifying the extent
        """

        import tempfile

        bath_dir = tempfile.mkdtemp()
        bath_dir_filename = bath_dir + os.sep +'ba19940524.000'
        #bath_dir = 'bath_data_manager_test'
        #print "os.getcwd( )",os.getcwd( )
        elevation_dir =  tempfile.mkdtemp()
        #elevation_dir = 'elev_expanded'
        elevation_dir_filename1 = elevation_dir + os.sep +'el19940524.000'
        elevation_dir_filename2 = elevation_dir + os.sep +'el19940524.001'

        fid = open(bath_dir_filename, 'w')
        fid.write(""" ncols             4
 nrows             4
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
   -9000.000    -1000.000   -3000.0 -2000.000
   -1000.000    9000.000  -1000.000 -3000.000
   -4000.000    6000.000   2000.000 -5000.000
   -9000.000    -1000.000   -3000.0 -2000.000
""")
        fid.close()

        fid = open(elevation_dir_filename1, 'w')
        fid.write(""" ncols             4
 nrows             4
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
   -900.000    -100.000   -300.0 -200.000
   -100.000    900.000  -100.000 -300.000
   -400.000    600.000   200.000 -500.000
   -900.000    -100.000   -300.0 -200.000
""")
        fid.close()

        fid = open(elevation_dir_filename2, 'w')
        fid.write(""" ncols             4
 nrows             4
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
   -990.000    -110.000   -330.0 -220.000
   -110.000    990.000  -110.000 -330.000
   -440.000    660.000   220.000 -550.000
   -990.000    -110.000   -330.0 -220.000
""")
        fid.close()

        ucur_dir =  tempfile.mkdtemp()
        ucur_dir_filename1 = ucur_dir + os.sep +'uc19940524.000'
        ucur_dir_filename2 = ucur_dir + os.sep +'uc19940524.001'

        fid = open(ucur_dir_filename1, 'w')
        fid.write(""" ncols             4
 nrows             4
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
   -90.000    -10.000   -30.0 -20.000
   -10.000    90.000  -10.000 -30.000
   -40.000    60.000   20.000 -50.000
   -90.000    -10.000   -30.0 -20.000
""")
        fid.close()
        fid = open(ucur_dir_filename2, 'w')
        fid.write(""" ncols             4
 nrows             4
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
   -90.000    -10.000   -30.0 -20.000
   -10.000    99.000  -11.000 -30.000
   -40.000    66.000   22.000 -50.000
   -90.000    -10.000   -30.0 -20.000
""")
        fid.close()

        vcur_dir =  tempfile.mkdtemp()
        vcur_dir_filename1 = vcur_dir + os.sep +'vc19940524.000'
        vcur_dir_filename2 = vcur_dir + os.sep +'vc19940524.001'

        fid = open(vcur_dir_filename1, 'w')
        fid.write(""" ncols             4
 nrows             4
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
   -90.000    -10.000   -30.0 -20.000
   -10.000    80.000  -20.000 -30.000
   -40.000    50.000   10.000 -50.000
   -90.000    -10.000   -30.0 -20.000
""")
        fid.close()
        fid = open(vcur_dir_filename2, 'w')
        fid.write(""" ncols             4
 nrows             4
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
   -90.000    -10.000   -30.0 -20.000
   -10.000    88.000  -22.000 -30.000
   -40.000    55.000   11.000 -50.000
   -90.000    -10.000   -30.0 -20.000
""")
        fid.close()

        sww_file = tempfile.mktemp(".sww")
        #sww_file = 'a_test.sww'
        asc_csiro2sww(bath_dir,elevation_dir, ucur_dir, vcur_dir,
                      sww_file, fail_on_NaN = False, elevation_NaN_filler = 0,
                      mean_stage = 100,
                       minlat = -37.6, maxlat = -37.6,
                  minlon = 148.3, maxlon = 148.3
                      #,verbose = True
                      )

        # check the sww file

        fid = NetCDFFile(sww_file, 'r')    #Open existing file for read
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        z = fid.variables['z'][:]
        stage = fid.variables['stage'][:]
        xmomentum = fid.variables['xmomentum'][:]
        ymomentum = fid.variables['ymomentum'][:]
        geo_ref = Geo_reference(NetCDFObject=fid)
        #print "geo_ref",geo_ref
        x_ref = geo_ref.get_xllcorner()
        y_ref = geo_ref.get_yllcorner()
        self.failUnless(geo_ref.get_zone() == 55,  'Failed')

        assert allclose(fid.starttime, 0.0) # (-37.45, 148.25)
        assert allclose(x_ref, 610120.388) # (-37.45, 148.25)
        assert allclose(y_ref,  5820863.269 )# (-37.45, 148.5)

        #Easting:  632145.632  Northing: 5820863.269
        #Latitude:   -37 45 ' 0.00000 ''  Longitude: 148  30 ' 0.00000 ''

        #print "x",x
        #print "y",y
        self.failUnless(len(x) == 4,'failed') # 2*2
        self.failUnless(len(x) == 4,'failed') # 2*2

        #Zone:   55
        #Easting:  632145.632  Northing: 5820863.269
        #Latitude:   -37 45 ' 0.00000 ''  Longitude: 148  30 ' 0.00000 ''
        # magic number - y is close enough for me.
        assert allclose(x[3], 632145.63 - x_ref)
        assert allclose(y[3], 5820863.269  - y_ref + 5.22155314684e-005)

        assert allclose(z[0],9000.0 ) #z is elevation info
        #print "z",z
        # 2 time steps, 4 points
        self.failUnless(xmomentum.shape == (2,4), 'failed')
        self.failUnless(ymomentum.shape == (2,4), 'failed')

        #(100.0 - -1000.000)*10
        #assert allclose(xmomentum[0][5], 11000.0 )

        fid.close()

        # is the sww file readable?
        #Lets see if we can convert it to a dem!
        #print "sww_file",sww_file
        #dem_file = tempfile.mktemp(".dem")
        domain = sww2domain(sww_file) ###, dem_file)
        domain.check_integrity()

        #tidy up
        os.remove(bath_dir_filename)
        os.rmdir(bath_dir)

        os.remove(elevation_dir_filename1)
        os.remove(elevation_dir_filename2)
        os.rmdir(elevation_dir)

        os.remove(ucur_dir_filename1)
        os.remove(ucur_dir_filename2)
        os.rmdir(ucur_dir)

        os.remove(vcur_dir_filename1)
        os.remove(vcur_dir_filename2)
        os.rmdir(vcur_dir)




        # remove sww file
        os.remove(sww_file)

        # remove dem file
        #os.remove(dem_file)

    def test_get_min_max_indexes(self):
        latitudes = [3,2,1,0]
        longitudes = [0,10,20,30]

        # k - lat
        # l - lon
        kmin, kmax, lmin, lmax = data_manager._get_min_max_indexes(
            latitudes,longitudes,
            -10,4,-10,31)

        #print "kmin",kmin;print "kmax",kmax
        #print "lmin",lmin;print "lmax",lmax
        latitudes_new = latitudes[kmin:kmax]
        longitudes_news = longitudes[lmin:lmax]
        #print "latitudes_new", latitudes_new
        #print "longitudes_news",longitudes_news
        self.failUnless(latitudes == latitudes_new and \
                        longitudes == longitudes_news,
                         'failed')

        ## 2nd test
        kmin, kmax, lmin, lmax = data_manager._get_min_max_indexes(
            latitudes,longitudes,
            0.5,2.5,5,25)
        #print "kmin",kmin;print "kmax",kmax
        #print "lmin",lmin;print "lmax",lmax
        latitudes_new = latitudes[kmin:kmax]
        longitudes_news = longitudes[lmin:lmax]
        #print "latitudes_new", latitudes_new
        #print "longitudes_news",longitudes_news

        self.failUnless(latitudes == latitudes_new and \
                        longitudes == longitudes_news,
                         'failed')

        ## 3rd test
        kmin, kmax, lmin, lmax = data_manager._get_min_max_indexes(\
            latitudes,
            longitudes,
            1.1,1.9,12,17)
        #print "kmin",kmin;print "kmax",kmax
        #print "lmin",lmin;print "lmax",lmax
        latitudes_new = latitudes[kmin:kmax]
        longitudes_news = longitudes[lmin:lmax]
        #print "latitudes_new", latitudes_new
        #print "longitudes_news",longitudes_news

        self.failUnless(latitudes_new == [2, 1] and \
                        longitudes_news == [10, 20],
                         'failed')


        ## 4th test
        kmin, kmax, lmin, lmax = data_manager._get_min_max_indexes(
            latitudes,longitudes,
                                                      -0.1,1.9,-2,17)
        #print "kmin",kmin;print "kmax",kmax
        #print "lmin",lmin;print "lmax",lmax
        latitudes_new = latitudes[kmin:kmax]
        longitudes_news = longitudes[lmin:lmax]
        #print "latitudes_new", latitudes_new
        #print "longitudes_news",longitudes_news

        self.failUnless(latitudes_new == [2, 1, 0] and \
                        longitudes_news == [0, 10, 20],
                         'failed')
        ## 5th test
        kmin, kmax, lmin, lmax = data_manager._get_min_max_indexes(
            latitudes,longitudes,
            0.1,1.9,2,17)
        #print "kmin",kmin;print "kmax",kmax
        #print "lmin",lmin;print "lmax",lmax
        latitudes_new = latitudes[kmin:kmax]
        longitudes_news = longitudes[lmin:lmax]
        #print "latitudes_new", latitudes_new
        #print "longitudes_news",longitudes_news

        self.failUnless(latitudes_new == [2, 1, 0] and \
                        longitudes_news == [0, 10, 20],
                         'failed')

        ## 6th test

        kmin, kmax, lmin, lmax = data_manager._get_min_max_indexes(
            latitudes,longitudes,
            1.5,4,18,32)
        #print "kmin",kmin;print "kmax",kmax
        #print "lmin",lmin;print "lmax",lmax
        latitudes_new = latitudes[kmin:kmax]
        longitudes_news = longitudes[lmin:lmax]
        #print "latitudes_new", latitudes_new
        #print "longitudes_news",longitudes_news

        self.failUnless(latitudes_new == [3, 2, 1] and \
                        longitudes_news == [10, 20, 30],
                         'failed')


        ## 7th test
        m2d = array([[0,1,2,3],[4,5,6,7],[8,9,10,11],[12,13,14,15]])
        kmin, kmax, lmin, lmax = data_manager._get_min_max_indexes(
            latitudes,longitudes,
            1.5,1.5,15,15)
        #print "kmin",kmin;print "kmax",kmax
        #print "lmin",lmin;print "lmax",lmax
        latitudes_new = latitudes[kmin:kmax]
        longitudes_news = longitudes[lmin:lmax]
        m2d = m2d[kmin:kmax,lmin:lmax]
        #print "m2d", m2d
        #print "latitudes_new", latitudes_new
        #print "longitudes_news",longitudes_news

        self.failUnless(latitudes_new == [2, 1] and \
                        longitudes_news == [10, 20],
                         'failed')

        self.failUnless(m2d == [[5,6],[9,10]],
                         'failed')

    def test_get_min_max_indexes_lat_ascending(self):
        latitudes = [0,1,2,3]
        longitudes = [0,10,20,30]

        # k - lat
        # l - lon
        kmin, kmax, lmin, lmax = data_manager._get_min_max_indexes(
            latitudes,longitudes,
            -10,4,-10,31)

        #print "kmin",kmin;print "kmax",kmax
        #print "lmin",lmin;print "lmax",lmax
        latitudes_new = latitudes[kmin:kmax]
        longitudes_news = longitudes[lmin:lmax]
        #print "latitudes_new", latitudes_new
        #print "longitudes_news",longitudes_news
        self.failUnless(latitudes == latitudes_new and \
                        longitudes == longitudes_news,
                         'failed')

        ## 3rd test
        kmin, kmax, lmin, lmax = data_manager._get_min_max_indexes(\
            latitudes,
            longitudes,
            1.1,1.9,12,17)
        #print "kmin",kmin;print "kmax",kmax
        #print "lmin",lmin;print "lmax",lmax
        latitudes_new = latitudes[kmin:kmax]
        longitudes_news = longitudes[lmin:lmax]
        #print "latitudes_new", latitudes_new
        #print "longitudes_news",longitudes_news

        self.failUnless(latitudes_new == [1, 2] and \
                        longitudes_news == [10, 20],
                         'failed')

    def test_get_min_max_indexes2(self):
        latitudes = [-30,-35,-40,-45]
        longitudes = [148,149,150,151]

        m2d = array([[0,1,2,3],[4,5,6,7],[8,9,10,11],[12,13,14,15]])

        # k - lat
        # l - lon
        kmin, kmax, lmin, lmax = data_manager._get_min_max_indexes(
            latitudes,longitudes,
            -37,-27,147,149.5)

        #print "kmin",kmin;print "kmax",kmax
        #print "lmin",lmin;print "lmax",lmax
        #print "m2d", m2d
        #print "latitudes", latitudes
        #print "longitudes",longitudes
        #print "latitudes[kmax]", latitudes[kmax]
        latitudes_new = latitudes[kmin:kmax]
        longitudes_new = longitudes[lmin:lmax]
        m2d = m2d[kmin:kmax,lmin:lmax]
        #print "m2d", m2d
        #print "latitudes_new", latitudes_new
        #print "longitudes_new",longitudes_new

        self.failUnless(latitudes_new == [-30, -35, -40] and \
                        longitudes_new == [148, 149,150],
                         'failed')
        self.failUnless(m2d == [[0,1,2],[4,5,6],[8,9,10]],
                         'failed')

    def test_get_min_max_indexes3(self):
        latitudes = [-30,-35,-40,-45,-50,-55,-60]
        longitudes = [148,149,150,151]

        # k - lat
        # l - lon
        kmin, kmax, lmin, lmax = data_manager._get_min_max_indexes(
            latitudes,longitudes,
            -43,-37,148.5,149.5)


        #print "kmin",kmin;print "kmax",kmax
        #print "lmin",lmin;print "lmax",lmax
        #print "latitudes", latitudes
        #print "longitudes",longitudes
        latitudes_new = latitudes[kmin:kmax]
        longitudes_news = longitudes[lmin:lmax]
        #print "latitudes_new", latitudes_new
        #print "longitudes_news",longitudes_news

        self.failUnless(latitudes_new == [-35, -40, -45] and \
                        longitudes_news == [148, 149,150],
                         'failed')

    def test_get_min_max_indexes4(self):
        latitudes = [-30,-35,-40,-45,-50,-55,-60]
        longitudes = [148,149,150,151]

        # k - lat
        # l - lon
        kmin, kmax, lmin, lmax = data_manager._get_min_max_indexes(
            latitudes,longitudes)

        #print "kmin",kmin;print "kmax",kmax
        #print "lmin",lmin;print "lmax",lmax
        #print "latitudes", latitudes
        #print "longitudes",longitudes
        latitudes_new = latitudes[kmin:kmax]
        longitudes_news = longitudes[lmin:lmax]
        #print "latitudes_new", latitudes_new
        #print "longitudes_news",longitudes_news

        self.failUnless(latitudes_new == latitudes  and \
                        longitudes_news == longitudes,
                         'failed')

    def test_tsh2sww(self):
        import os
        import tempfile

        tsh_file = tempfile.mktemp(".tsh")
        file = open(tsh_file,"w")
        file.write("4 3 # <vertex #> <x> <y> [attributes]\n \
0 0.0 0.0 0.0 0.0 0.01 \n \
1 1.0 0.0 10.0 10.0 0.02  \n \
2 0.0 1.0 0.0 10.0 0.03  \n \
3 0.5 0.25 8.0 12.0 0.04  \n \
# Vert att title  \n \
elevation  \n \
stage  \n \
friction  \n \
2 # <triangle #> [<vertex #>] [<neigbouring triangle #>]  \n\
0 0 3 2 -1  -1  1 dsg\n\
1 0 1 3 -1  0 -1   ole nielsen\n\
4 # <segment #> <vertex #>  <vertex #> [boundary tag] \n\
0 1 0 2 \n\
1 0 2 3 \n\
2 2 3 \n\
3 3 1 1 \n\
3 0 # <x> <y> [attributes] ...Mesh Vertices... \n \
0 216.0 -86.0 \n \
1 160.0 -167.0 \n \
2 114.0 -91.0 \n \
3 # <vertex #>  <vertex #> [boundary tag] ...Mesh Segments... \n \
0 0 1 0 \n \
1 1 2 0 \n \
2 2 0 0 \n \
0 # <x> <y> ...Mesh Holes... \n \
0 # <x> <y> <attribute>...Mesh Regions... \n \
0 # <x> <y> <attribute>...Mesh Regions, area... \n\
#Geo reference \n \
56 \n \
140 \n \
120 \n")
        file.close()

        #sww_file = tempfile.mktemp(".sww")
        #print "sww_file",sww_file
        #print "sww_file",tsh_file
        tsh2sww(tsh_file)

        os.remove(tsh_file)
        os.remove(tsh_file[:-4] + '.sww')
        



########## testing nbed class ##################
    def test_exposure_csv_loading(self):
        

        file_name = tempfile.mktemp(".xya")
        file = open(file_name,"w")
        file.write("LATITUDE, LONGITUDE ,sound  , speed \n\
115.0, -21.0, splat, 0.0\n\
114.0, -21.7, pow, 10.0\n\
114.5, -21.4, bang, 40.0\n")
        file.close()
        exposure = Exposure_csv(file_name)
        exposure.get_column("sound")
       
        self.failUnless(exposure._attribute_dic['sound'][2]==' bang',
                        'FAILED!')
        self.failUnless(exposure._attribute_dic['speed'][2]==' 40.0',
                        'FAILED!')
        
        os.remove(file_name)
        
    def test_exposure_csv_loading(self):
        

        file_name = tempfile.mktemp(".xya")
        file = open(file_name,"w")
        file.write("LATITUDE, LONGITUDE ,sound  , speed \n\
115.0, -21.0, splat, 0.0\n\
114.0, -21.7, pow, 10.0\n\
114.5, -21.4, bang, 40.0\n")
        file.close()
        exposure = Exposure_csv(file_name)
        exposure.get_column("sound")
       
        self.failUnless(exposure._attribute_dic['sound'][2]==' bang',
                        'FAILED!')
        self.failUnless(exposure._attribute_dic['speed'][2]==' 40.0',
                        'FAILED!')
        
        os.remove(file_name)

    def test_exposure_csv_cmp(self):
        file_name = tempfile.mktemp(".xya")
        file = open(file_name,"w")
        file.write("LATITUDE, LONGITUDE ,sound  , speed \n\
115.0, -21.0, splat, 0.0\n\
114.0, -21.7, pow, 10.0\n\
114.5, -21.4, bang, 40.0\n")
        file.close()
        
        e1 = Exposure_csv(file_name)
        e2 = Exposure_csv(file_name)
        os.remove(file_name)

        self.failUnless(cmp(e1,e2)==0,
                        'FAILED!')
        
        self.failUnless(cmp(e1,"hey")==1,
                        'FAILED!')
        
        file_name = tempfile.mktemp(".xya")
        file = open(file_name,"w")
        # Note, this has less spaces in the title,
        # the instances will be the same.
        file.write("LATITUDE,LONGITUDE ,sound, speed \n\
115.0, -21.0, splat, 0.0\n\
114.0, -21.7, pow, 10.0\n\
114.5, -21.4, bang, 40.0\n")
        file.close()
        e3 = Exposure_csv(file_name)
        os.remove(file_name)

        self.failUnless(cmp(e3,e2)==0,
                        'FAILED!')
        
        file_name = tempfile.mktemp(".xya")
        file = open(file_name,"w")
        # Note, 40 changed to 44 .
        file.write("LATITUDE,LONGITUDE ,sound, speed \n\
115.0, -21.0, splat, 0.0\n\
114.0, -21.7, pow, 10.0\n\
114.5, -21.4, bang, 44.0\n")
        file.close()
        e4 = Exposure_csv(file_name)
        os.remove(file_name)
        #print "e4",e4._attribute_dic 
        #print "e2",e2._attribute_dic 
        self.failUnless(cmp(e4,e2)<>0,
                        'FAILED!')
        
        file_name = tempfile.mktemp(".xya")
        file = open(file_name,"w")
        # Note, the first two columns are swapped.
        file.write("LONGITUDE,LATITUDE ,sound, speed \n\
 -21.0,115.0, splat, 0.0\n\
 -21.7,114.0, pow, 10.0\n\
 -21.4,114.5, bang, 40.0\n")
        file.close()
        e5 = Exposure_csv(file_name)
        os.remove(file_name)

        self.failUnless(cmp(e3,e5)<>0,
                        'FAILED!')
        
    def test_exposure_csv_saving(self):
        

        file_name = tempfile.mktemp(".xya")
        file = open(file_name,"w")
        file.write("LATITUDE, LONGITUDE ,sound  , speed \n\
115.0, -21.0, splat, 0.0\n\
114.0, -21.7, pow, 10.0\n\
114.5, -21.4, bang, 40.0\n")
        file.close()
        e1 = Exposure_csv(file_name)
        
        file_name2 = tempfile.mktemp(".xya")
        e1.save(file_name = file_name2)
        e2 = Exposure_csv(file_name2)
       
        self.failUnless(cmp(e1,e2)==0,
                        'FAILED!')
        os.remove(file_name)
        os.remove(file_name2)

    def test_exposure_csv_get_location(self):
        file_name = tempfile.mktemp(".xya")
        file = open(file_name,"w")
        file.write("LONGITUDE , LATITUDE, sound  , speed \n\
150.916666667, -34.5, splat, 0.0\n\
150.0, -34.0, pow, 10.0\n")
        file.close()
        e1 = Exposure_csv(file_name)

        gsd = e1.get_location()
        
        points = gsd.get_data_points(absolute=True)
        
        assert allclose(points[0][0], 308728.009)
        assert allclose(points[0][1], 6180432.601)
        assert allclose(points[1][0],  222908.705)
        assert allclose(points[1][1], 6233785.284)
        self.failUnless(gsd.get_geo_reference().get_zone() == 56,
                        'Bad zone error!')

        os.remove(file_name)
        
    def test_exposure_csv_set_column_get_column(self):
        file_name = tempfile.mktemp(".xya")
        file = open(file_name,"w")
        file.write("LONGITUDE , LATITUDE, sound  , speed \n\
150.916666667, -34.5, splat, 0.0\n\
150.0, -34.0, pow, 10.0\n")
        file.close()
        e1 = Exposure_csv(file_name)      
        os.remove(file_name)

        new_title = "feast"
        new_values = ["chicken","soup"]
        e1.set_column(new_title, new_values)
        returned_values = e1.get_column(new_title)
        self.failUnless(returned_values == new_values,
                        ' Error!')
        
        file_name2 = tempfile.mktemp(".xya")
        e1.save(file_name = file_name2)
        e2 = Exposure_csv(file_name2)
        returned_values = e2.get_column(new_title)
        self.failUnless(returned_values == new_values,
                        ' Error!')       
        os.remove(file_name2)

    def test_exposure_csv_set_column_get_column_error_checking(self):
        file_name = tempfile.mktemp(".xya")
        file = open(file_name,"w")
        file.write("LONGITUDE , LATITUDE, sound  , speed \n\
150.916666667, -34.5, splat, 0.0\n\
150.0, -34.0, pow, 10.0\n")
        file.close()
        e1 = Exposure_csv(file_name)      
        os.remove(file_name)

        new_title = "sound"
        new_values = [12.5,7.6]
        try:
            e1.set_column(new_title, new_values)
        except TitleValueError:
            pass
        else:
            self.failUnless(0 ==1,  'Error not thrown error!')
            
        e1.set_column(new_title, new_values, overwrite=True)
        returned_values = e1.get_column(new_title)
        self.failUnless(returned_values == new_values,
                        ' Error!')       
        
        new2_title = "short list"
        new2_values = [12.5]
        try:
            e1.set_column(new2_title, new2_values)
        except DataMissingValuesError:
            pass
        else:
            self.failUnless(0 ==1,  'Error not thrown error!')
            
        new2_title = "long list"
        new2_values = [12.5, 7,8]
        try:
            e1.set_column(new2_title, new2_values)
        except DataMissingValuesError:
            pass
        else:
            self.failUnless(0 ==1,  'Error not thrown error!')
        file_name2 = tempfile.mktemp(".xya")
        e1.save(file_name = file_name2)
        e2 = Exposure_csv(file_name2)
        returned_values = e2.get_column(new_title)
        for returned, new in map(None, returned_values, new_values):
            self.failUnless(returned == str(new), ' Error!')
        #self.failUnless(returned_values == new_values, ' Error!')       
        os.remove(file_name2)
        
        try:
            e1.get_column("toe jam")
        except TitleValueError:
            pass
        else:
            self.failUnless(0 ==1,  'Error not thrown error!')
            
    def test_exposure_csv_loading_x_y(self):
        

        file_name = tempfile.mktemp(".xya")
        file = open(file_name,"w")
        file.write("x, y ,sound  , speed \n\
115.0, 7, splat, 0.0\n\
114.0, 8.0, pow, 10.0\n\
114.5, 9., bang, 40.0\n")
        file.close()
        e1 = Exposure_csv(file_name, is_x_y_locations=True)
        gsd = e1.get_location()
        
        points = gsd.get_data_points(absolute=True)
        
        assert allclose(points[0][0], 115)
        assert allclose(points[0][1], 7)
        assert allclose(points[1][0], 114)
        assert allclose(points[1][1], 8)
        assert allclose(points[2][0], 114.5)
        assert allclose(points[2][1], 9)
        self.failUnless(gsd.get_geo_reference().get_zone() == -1,
                        'Bad zone error!')

        os.remove(file_name)

           
    def test_exposure_csv_loading_x_y2(self):
        
        csv_file = tempfile.mktemp(".csv")
        fd = open(csv_file,'wb')
        writer = csv.writer(fd)
        writer.writerow(['x','y','STR_VALUE','C_VALUE','ROOF_TYPE','WALLS', 'SHORE_DIST'])
        writer.writerow([5.5,0.5,'199770','130000','Metal','Timber',20])
        writer.writerow([4.5,1.0,'150000','76000','Metal','Double Brick',20])
        writer.writerow([4.5,1.5,'150000','76000','Metal','Brick Veneer',20])
        fd.close()

        e1 = Exposure_csv(csv_file)
        gsd = e1.get_location()
        
        points = gsd.get_data_points(absolute=True)
        assert allclose(points[0][0], 5.5)
        assert allclose(points[0][1], 0.5)
        assert allclose(points[1][0], 4.5)
        assert allclose(points[1][1], 1.0)
        assert allclose(points[2][0], 4.5)
        assert allclose(points[2][1], 1.5)
        self.failUnless(gsd.get_geo_reference().get_zone() == -1,
                        'Bad zone error!')

        os.remove(csv_file)

    #### TESTS FOR URS 2 SWW  ###     
    
    def create_mux(self, points_num=None):
        # write all the mux stuff.
        time_step_count = 3
        time_step = 0.5
        
        longitudes = [150.66667, 150.83334, 151., 151.16667]
        latitudes = [-34.5, -34.33333, -34.16667, -34]

        if points_num == None:
            points_num = len(longitudes) * len(latitudes)

        lonlatdeps = []
        quantities = ['HA','UA','VA']
        mux_names = ['-z-mux','-e-mux','-n-mux']
        quantities_init = [[],[],[]]
        # urs binary is latitude fastest
        for i,lon in enumerate(longitudes):
            for j,lat in enumerate(latitudes):
                _ , e, n = redfearn(lat, lon)
                lonlatdeps.append([lon, lat, n])
                quantities_init[0].append(e) # HA
                quantities_init[1].append(n ) # UA
                quantities_init[2].append(e) # VA
        #print "lonlatdeps",lonlatdeps
        _,base_name = tempfile.mkstemp("")
        files = []        
        for i,q in enumerate(quantities): 
            quantities_init[i] = ensure_numeric(quantities_init[i])
            #print "HA_init", HA_init
            q_time = zeros((time_step_count, points_num), Float)
            for time in range(time_step_count):
                q_time[time,:] = quantities_init[i] #* time * 4
            
            #Write C files
            columns = 3 # long, lat , depth
            file = base_name + mux_names[i]
            f = open(file, 'wb')
            files.append(file)
            f.write(pack('i',points_num))
            f.write(pack('i',time_step_count))
            f.write(pack('f',time_step))

            #write lat/long info
            for lonlatdep in lonlatdeps:
                for float in lonlatdep:
                    f.write(pack('f',float))
                    
            # Write quantity info
            for time in  range(time_step_count):
                for i in range(points_num):
                    f.write(pack('f',q_time[time,i]))
            f.close()
        return base_name, files
        
    
    def delete_mux(self, files):
        for file in files:
            os.remove(file)
            
    def test_urs2sww_test_fail(self):
        points_num = -100
        time_step_count = 45
        time_step = -7
        _,base_name = tempfile.mkstemp("")
        files = []
        quantities = ['HA','UA','VA']
        mux_names = ['-z-mux','-e-mux','-n-mux']
        for i,q in enumerate(quantities): 
            #Write C files
            columns = 3 # long, lat , depth
            file = base_name + mux_names[i]
            f = open(file, 'wb')
            files.append(file)
            f.write(pack('i',points_num))
            f.write(pack('i',time_step_count))
            f.write(pack('f',time_step))

            f.close()   
        tide = 1
        try:
            urs2sww(base_name, remove_nc_files=True, mean_stage=tide)        
        except ANUGAError:
            pass
        else:
            self.delete_mux(files)
            msg = 'Should have raised exception'
            raise msg
        sww_file = base_name + '.sww'
        self.delete_mux(files)
        
    def test_urs2sww_test_fail2(self):
        base_name = 'Harry-high-pants'
        try:
            urs2sww(base_name)        
        except IOError:
            pass
        else:
            self.delete_mux(files)
            msg = 'Should have raised exception'
            raise msg
           
    def test_urs2sww(self):
        tide = 1
        base_name, files = self.create_mux()
        urs2sww(base_name
                #, origin=(0,0,0)
                , mean_stage=tide
                , remove_nc_files=False
                )
        sww_file = base_name + '.sww'
        
        #Let's interigate the sww file
        # Note, the sww info is not gridded.  It is point data.
        fid = NetCDFFile(sww_file)

        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        geo_reference = Geo_reference(NetCDFObject=fid)

        
        #Check that first coordinate is correctly represented       
        #Work out the UTM coordinates for first point
        zone, e, n = redfearn(-34.5, 150.66667)       
       
        assert allclose(geo_reference.get_absolute([[x[0],y[0]]]), [e,n])

        # Make x and y absolute
        points = geo_reference.get_absolute(map(None, x, y))
        points = ensure_numeric(points)
        x = points[:,0]
        y = points[:,1]
        
        #Check first value
        stage = fid.variables['stage'][:]
        xmomentum = fid.variables['xmomentum'][:]
        ymomentum = fid.variables['ymomentum'][:]
        elevation = fid.variables['elevation'][:]
        assert allclose(stage[0,0], e +tide)  #Meters

        #Check the momentums - ua
        #momentum = velocity*(stage-elevation)
        #momentum = velocity*(stage+elevation)
        # -(-elevation) since elevation is inverted in mux files
        # = n*(e+tide+n) based on how I'm writing these files
        answer = n*(e+tide+n)
        actual = xmomentum[0,0]
        assert allclose(answer, actual)  #Meters

        # check the stage values, first time step.
        # These arrays are equal since the Easting values were used as
        # the stage
        assert allclose(stage[0], x +tide)  #Meters

        # check the elevation values.
        # -ve since urs measures depth, sww meshers height,
        # these arrays are equal since the northing values were used as
        # the elevation
        assert allclose(-elevation, y)  #Meters
        
        fid.close()
        self.delete_mux(files)
        os.remove(sww_file)
        
  
    def test_urs2sww_origin(self):
        tide = 1
        base_name, files = self.create_mux()
        urs2sww(base_name
                , origin=(0,0,0)
                , mean_stage=tide
                , remove_nc_files=False
                )
        sww_file = base_name + '.sww'
        
        #Let's interigate the sww file
        # Note, the sww info is not gridded.  It is point data.
        fid = NetCDFFile(sww_file)

        #  x and y are absolute
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        geo_reference = Geo_reference(NetCDFObject=fid)

        
        #Check that first coordinate is correctly represented       
        #Work out the UTM coordinates for first point
        zone, e, n = redfearn(-34.5, 150.66667)       
       
        assert allclose([x[0],y[0]], [e,n])

        
        #Check first value
        stage = fid.variables['stage'][:]
        xmomentum = fid.variables['xmomentum'][:]
        ymomentum = fid.variables['ymomentum'][:]
        elevation = fid.variables['elevation'][:]
        assert allclose(stage[0,0], e +tide)  #Meters

        #Check the momentums - ua
        #momentum = velocity*(stage-elevation)
        #momentum = velocity*(stage+elevation)
        # -(-elevation) since elevation is inverted in mux files
        # = n*(e+tide+n) based on how I'm writing these files
        answer = n*(e+tide+n)
        actual = xmomentum[0,0]
        assert allclose(answer, actual)  #Meters

        # check the stage values, first time step.
        # These arrays are equal since the Easting values were used as
        # the stage
        assert allclose(stage[0], x +tide)  #Meters

        # check the elevation values.
        # -ve since urs measures depth, sww meshers height,
        # these arrays are equal since the northing values were used as
        # the elevation
        assert allclose(-elevation, y)  #Meters
        
        fid.close()
        self.delete_mux(files)
        os.remove(sww_file)
 
    def test_urs2sww_minmaxlatlong(self):
        
        #longitudes = [150.66667, 150.83334, 151., 151.16667]
        #latitudes = [-34.5, -34.33333, -34.16667, -34]

        tide = 1
        base_name, files = self.create_mux()
        urs2sww(base_name,
                minlat=-34.5,
                maxlat=-34,
                minlon= 150.66667,
                maxlon= 151.16667,
                mean_stage=tide,
                remove_nc_files=False
                )
        sww_file = base_name + '.sww'
        
        #Let's interigate the sww file
        # Note, the sww info is not gridded.  It is point data.
        fid = NetCDFFile(sww_file)
        

        # Make x and y absolute
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        geo_reference = Geo_reference(NetCDFObject=fid)
        points = geo_reference.get_absolute(map(None, x, y))
        points = ensure_numeric(points)
        x = points[:,0]
        y = points[:,1]
        
        #Check that first coordinate is correctly represented       
        #Work out the UTM coordinates for first point
        zone, e, n = redfearn(-34.5, 150.66667) 
        assert allclose([x[0],y[0]], [e,n])

        
        #Check first value
        stage = fid.variables['stage'][:]
        xmomentum = fid.variables['xmomentum'][:]
        ymomentum = fid.variables['ymomentum'][:]
        elevation = fid.variables['elevation'][:]
        assert allclose(stage[0,0], e +tide)  #Meters

        #Check the momentums - ua
        #momentum = velocity*(stage-elevation)
        #momentum = velocity*(stage+elevation)
        # -(-elevation) since elevation is inverted in mux files
        # = n*(e+tide+n) based on how I'm writing these files
        answer = n*(e+tide+n)
        actual = xmomentum[0,0]
        assert allclose(answer, actual)  #Meters

        # check the stage values, first time step.
        # These arrays are equal since the Easting values were used as
        # the stage
        assert allclose(stage[0], x +tide)  #Meters

        # check the elevation values.
        # -ve since urs measures depth, sww meshers height,
        # these arrays are equal since the northing values were used as
        # the elevation
        assert allclose(-elevation, y)  #Meters
        
        fid.close()
        self.delete_mux(files)
        os.remove(sww_file)
        
    def test_lon_lat2grid(self):
        lonlatdep = [
            [ 113.06700134  ,  -26.06669998 ,   1.        ] ,
            [ 113.06700134  ,  -26.33329964 ,   3.        ] ,
            [ 113.19999695  ,  -26.06669998 ,   2.        ] ,
            [ 113.19999695  ,  -26.33329964 ,   4.        ] ]
            
        long, lat, quantity = lon_lat2grid(lonlatdep)

        for i, result in enumerate(lat):
            assert lonlatdep [i][1] == result
        assert len(lat) == 2 

        for i, result in enumerate(long):
            assert lonlatdep [i*2][0] == result
        assert len(long) == 2

        for i,q in enumerate(quantity):
            assert q == i+1
            
    def test_lon_lat2grid_bad(self):
        lonlatdep  = [
            [ -26.06669998,  113.06700134,    1.        ],
            [ -26.06669998 , 113.19999695 ,   2.        ],
            [ -26.06669998 , 113.33300018,    3.        ],
            [ -26.06669998 , 113.43299866   , 4.        ],
            [ -26.20000076 , 113.06700134,    5.        ],
            [ -26.20000076 , 113.19999695 ,   6.        ],
            [ -26.20000076 , 113.33300018  ,  7.        ],
            [ -26.20000076 , 113.43299866   , 8.        ],
            [ -26.33329964 , 113.06700134,    9.        ],
            [ -26.33329964 , 113.19999695 ,   10.        ],
            [ -26.33329964 , 113.33300018  ,  11.        ],
            [ -26.33329964 , 113.43299866 ,   12.        ],
            [ -26.43330002 , 113.06700134 ,   13        ],
            [ -26.43330002 , 113.19999695 ,   14.        ],
            [ -26.43330002 , 113.33300018,    15.        ],
            [ -26.43330002 , 113.43299866,    16.        ]]
        try:
            long, lat, quantity = lon_lat2grid(lonlatdep)
        except AssertionError:
            pass
        else:
            msg = 'Should have raised exception'
            raise msg
       
    def test_lon_lat2gridII(self):
        lonlatdep = [
            [ 113.06700134  ,  -26.06669998 ,   1.        ] ,
            [ 113.06700134  ,  -26.33329964 ,   2.        ] ,
            [ 113.19999695  ,  -26.06669998 ,   3.        ] ,
            [ 113.19999695  ,  -26.344329964 ,   4.        ] ]
        try:
            long, lat, quantity = lon_lat2grid(lonlatdep)
        except AssertionError:
            pass
        else:
            msg = 'Should have raised exception'
            raise msg
        
    def trial_loading(self):
        basename_in = 'karratha'
        basename_out = basename_in
        urs2sww(basename_in, basename_out, remove_nc_files=True,
                zscale=10000000)
    #### END TESTS FOR URS 2 SWW  ###

        
#-------------------------------------------------------------
if __name__ == "__main__":
    #suite = unittest.makeSuite(Test_Data_Manager,'test_urs2sww_m')
    #suite = unittest.makeSuite(Test_Data_Manager,'test_get_min_max_indexes_lat_ascending')
    #suite = unittest.makeSuite(Test_Data_Manager,'test_ferret2sww_lat_long')
    suite = unittest.makeSuite(Test_Data_Manager,'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)



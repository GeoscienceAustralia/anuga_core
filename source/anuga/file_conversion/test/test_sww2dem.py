import unittest, os
import numpy as num

from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular
from anuga.shallow_water.shallow_water_domain import Domain
from anuga.config import netcdf_mode_r
from anuga.file.netcdf import NetCDFFile

from anuga.coordinate_transforms.geo_reference import Geo_reference, \
     DEFAULT_ZONE

from anuga.file.sww import SWW_file
     
# boundary functions
from anuga.shallow_water.boundaries import Reflective_boundary, \
            Field_boundary, Transmissive_momentum_set_stage_boundary, \
            Transmissive_stage_zero_momentum_boundary
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions\
     import Transmissive_boundary, Dirichlet_boundary, \
            Time_boundary, File_boundary, AWI_boundary

# local modules
from anuga.file_conversion.sww2dem import sww2dem, sww2dem_batch

class Test_Sww2Dem(unittest.TestCase):
    def setUp(self):
        import time
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular
        
        # Create basic mesh
        points, vertices, boundary = rectangular(2, 2)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.default_order = 2

        # Set some field values
        domain.set_quantity('elevation', lambda x,y: -x)
        domain.set_quantity('friction', 0.03)


        ######################
        # Boundary conditions
        B = Transmissive_boundary(domain)
        domain.set_boundary( {'left': B, 'right': B, 'top': B, 'bottom': B})


        ######################
        #Initial condition - with jumps
        bed = domain.quantities['elevation'].vertex_values
        stage = num.zeros(bed.shape, num.float)

        h = 0.3
        for i in range(stage.shape[0]):
            if i % 2 == 0:
                stage[i,:] = bed[i,:] + h
            else:
                stage[i,:] = bed[i,:]

        domain.set_quantity('stage', stage)


        domain.distribute_to_vertices_and_edges()               
        self.domain = domain

        C = domain.get_vertex_coordinates()
        self.X = C[:,0:6:2].copy()
        self.Y = C[:,1:6:2].copy()

        self.F = bed
        self.verbose = False


    def tearDown(self):
        pass


    def test_sww2dem_asc_elevation_depth(self):
        """test_sww2dem_asc_elevation_depth
        
        Test that sww information can be converted correctly to asc/prj
        format readable by e.g. ArcView
        
        Also check geo_reference is correct
        """

        import time, os

        # Setup
        self.domain.set_name('datatest')

        prjfile = self.domain.get_name() + '_elevation.prj'
        ascfile = self.domain.get_name() + '_elevation.asc'
        swwfile = self.domain.get_name() + '.sww'

        self.domain.set_datadir('.')
        self.domain.format = 'sww'
        self.domain.smooth = True
        self.domain.set_quantity('elevation', lambda x,y: -x-y)
        self.domain.set_quantity('stage', 1.0)

        self.domain.geo_reference = Geo_reference(56, 308500, 6189000)

        sww = SWW_file(self.domain)
        sww.store_connectivity()
        sww.store_timestep()


        self.domain.evolve_to_end(finaltime = 0.01)
        sww.store_timestep()

        cellsize = 0.25
        #Check contents
        #Get NetCDF

        fid = NetCDFFile(sww.filename, netcdf_mode_r)

        # Get the variables
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        z = fid.variables['elevation'][:]
        time = fid.variables['time'][:]
        stage = fid.variables['stage'][:]
        
        # Check georeferencig: zone, xllcorner and yllcorner
        assert fid.zone == 56
        assert fid.xllcorner == 308500
        assert fid.yllcorner == 6189000
                

        fid.close()

        #Export to ascii/prj files
        sww2dem(self.domain.get_name()+'.sww',
                self.domain.get_name()+'_elevation.asc',
                quantity = 'elevation',
                cellsize = cellsize,
                number_of_decimal_places = 9,
                verbose = self.verbose)

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
        assert num.allclose(float(L[1].strip().lower()), 308500)

        L = lines[3].strip().split()
        assert L[0].strip().lower() == 'yllcorner'
        assert num.allclose(float(L[1].strip().lower()), 6189000)

        L = lines[4].strip().split()
        assert L[0].strip().lower() == 'cellsize'
        assert num.allclose(float(L[1].strip().lower()), cellsize)

        L = lines[5].strip().split()
        assert L[0].strip() == 'NODATA_value'
        assert L[1].strip().lower() == '-9999'

        #Check grid values
        for j in range(5):
            L = lines[6+j].strip().split()
            y = (4-j) * cellsize
            for i in range(5):
                assert num.allclose(float(L[i]), -i*cellsize - y)
                
        #Cleanup
        os.remove(prjfile)
        os.remove(ascfile)

        ascfile = self.domain.get_name() + '_depth.asc'
        prjfile = self.domain.get_name() + '_depth.prj'

        #Export to ascii/prj files
        sww2dem(self.domain.get_name()+'.sww',
                ascfile,
                quantity = 'depth',
                cellsize = cellsize,
                number_of_decimal_places = 9,
                verbose = self.verbose)
        
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
        assert num.allclose(float(L[1].strip().lower()), 308500)

        L = lines[3].strip().split()
        assert L[0].strip().lower() == 'yllcorner'
        assert num.allclose(float(L[1].strip().lower()), 6189000)

        L = lines[4].strip().split()
        assert L[0].strip().lower() == 'cellsize'
        assert num.allclose(float(L[1].strip().lower()), cellsize)

        L = lines[5].strip().split()
        assert L[0].strip() == 'NODATA_value'
        assert L[1].strip().lower() == '-9999'

        #Check grid values
        for j in range(5):
            L = lines[6+j].strip().split()
            y = (4-j) * cellsize
            for i in range(5):
                assert num.allclose(float(L[i]), 1 - (-i*cellsize - y))


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
        sww = SWW_file(domain)
        sww.store_connectivity()
        sww.store_timestep()
        
        domain.tight_slope_limiters = 1
        domain.evolve_to_end(finaltime = 0.01)
        sww.store_timestep()

        cellsize = 10  #10m grid


        #Check contents
        #Get NetCDF

        fid = NetCDFFile(sww.filename, netcdf_mode_r)

        # Get the variables
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        z = fid.variables['elevation'][:]
        time = fid.variables['time'][:]
        stage = fid.variables['stage'][:]


        #Export to ascii/prj files
        sww2dem(domain.get_name() + '.sww',
                domain.get_name() + '_elevation.asc',
                quantity = 'elevation',
                cellsize = cellsize,
                number_of_decimal_places = 9,
                verbose = self.verbose,
                block_size=2)


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
        assert num.allclose(float(L[1].strip().lower()), 308500)

        L = lines[3].strip().split()
        assert L[0].strip().lower() == 'yllcorner'
        assert num.allclose(float(L[1].strip().lower()), 6189000)

        L = lines[4].strip().split()
        assert L[0].strip().lower() == 'cellsize'
        assert num.allclose(float(L[1].strip().lower()), cellsize)

        L = lines[5].strip().split()
        assert L[0].strip() == 'NODATA_value'
        assert L[1].strip().lower() == '-9999'

        #Check grid values (FIXME: Use same strategy for other sww2dem tests)
        for i, line in enumerate(lines[6:]):
            for j, value in enumerate( line.split() ):
                assert num.allclose(float(value), -(10-i+j)*cellsize,
                                    atol=1.0e-12, rtol=1.0e-12)

                # Note: Equality can be obtained in this case,
                # but it is better to use allclose.
                #assert float(value) == -(10-i+j)*cellsize


        fid.close()

        #Cleanup
        os.remove(prjfile)
        os.remove(ascfile)
        os.remove(swwfile)



    def test_sww2dem_larger_zero(self):
        """Test example has rows with a large number of zeros

        ncols         2001
        nrows         2
        xllcorner     308500
        yllcorner     6189000
        cellsize      1.000000
        NODATA_value  -9999
        0.0 ....
        """

        import time, os

        #Setup

        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross

        #Create basic mesh (100m x 100m)
        points, vertices, boundary = rectangular_cross(20, 1, 20.0, 1.0)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.default_order = 1

        domain.set_name('datatest')

        prjfile = domain.get_name() + '_elevation.prj'
        ascfile = domain.get_name() + '_elevation.asc'
        swwfile = domain.get_name() + '.sww'

        domain.set_datadir('.')
        domain.format = 'sww'
        domain.smooth = True
        domain.geo_reference = Geo_reference(56, 308500, 6189000)

        #
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', 0)

        B = Transmissive_boundary(domain)
        domain.set_boundary( {'left': B, 'right': B, 'top': B, 'bottom': B})


        #
        sww = SWW_file(domain)
        sww.store_connectivity()
        sww.store_timestep()
        
        domain.tight_slope_limiters = 1
        domain.evolve_to_end(finaltime = 0.01)
        sww.store_timestep()

        cellsize = 1.0  #0.1 grid


        #Check contents
        #Get NetCDF

        fid = NetCDFFile(sww.filename, netcdf_mode_r)

        # Get the variables
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        z = fid.variables['elevation'][:]
        time = fid.variables['time'][:]
        stage = fid.variables['stage'][:]


        #Export to ascii/prj files
        sww2dem(domain.get_name()+'.sww',
                domain.get_name() + '_elevation.asc',
                quantity = 'elevation',
                cellsize = cellsize,
                number_of_decimal_places = 9,
                verbose = self.verbose,
                block_size=2)


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
        assert L[1].strip().lower() == '21'

        L = lines[1].strip().split()
        assert L[0].strip().lower() == 'nrows'
        assert L[1].strip().lower() == '2'

        L = lines[2].strip().split()
        assert L[0].strip().lower() == 'xllcorner'
        assert num.allclose(float(L[1].strip().lower()), 308500)

        L = lines[3].strip().split()
        assert L[0].strip().lower() == 'yllcorner'
        assert num.allclose(float(L[1].strip().lower()), 6189000)

        L = lines[4].strip().split()
        assert L[0].strip().lower() == 'cellsize'
        assert num.allclose(float(L[1].strip().lower()), cellsize)

        L = lines[5].strip().split()
        assert L[0].strip() == 'NODATA_value'
        assert L[1].strip().lower() == '-9999'

        #Check grid values (FIXME: Use same strategy for other sww2dem tests)
        for i, line in enumerate(lines[6:]):
            for j, value in enumerate( line.split() ):
                #print value
                assert num.allclose(float(value), 0.0,
                                    atol=1.0e-12, rtol=1.0e-12)

                # Note: Equality can be obtained in this case,
                # but it is better to use allclose.
                #assert float(value) == -(10-i+j)*cellsize


        fid.close()


        #Cleanup
        os.remove(prjfile)
        os.remove(ascfile)
        os.remove(swwfile)




    def test_sww2dem_boundingbox(self):
        """Test that mesh can be restricted by bounding box

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

        #Setup

        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

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
        sww = SWW_file(domain)
        sww.store_connectivity()
        sww.store_timestep()

        #domain.tight_slope_limiters = 1
        domain.evolve_to_end(finaltime = 0.01)
        sww.store_timestep()

        cellsize = 10  #10m grid


        #Check contents
        #Get NetCDF

        fid = NetCDFFile(sww.filename, netcdf_mode_r)

        # Get the variables
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        z = fid.variables['elevation'][:]
        time = fid.variables['time'][:]
        stage = fid.variables['stage'][:]


        # Export to ascii/prj files
        sww2dem(domain.get_name() + '.sww',
                domain.get_name() + '_elevation.asc',
                quantity = 'elevation',
                cellsize = cellsize,
                number_of_decimal_places = 9,
                easting_min = 308530,
                easting_max = 308570,
                northing_min = 6189050,
                northing_max = 6189100,
                verbose = self.verbose)

        fid.close()


        # Check prj (meta data)
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
        assert num.allclose(float(L[1].strip().lower()), 308530)

        L = lines[3].strip().split()
        assert L[0].strip().lower() == 'yllcorner'
        assert num.allclose(float(L[1].strip().lower()), 6189050)

        L = lines[4].strip().split()
        assert L[0].strip().lower() == 'cellsize'
        assert num.allclose(float(L[1].strip().lower()), cellsize)

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


        sww = SWW_file(self.domain)
        sww.store_connectivity()
        sww.store_timestep()

        #self.domain.tight_slope_limiters = 1
        self.domain.evolve_to_end(finaltime = 0.01)
        sww.store_timestep()

        cellsize = 0.25
        #Check contents
        #Get NetCDF

        fid = NetCDFFile(sww.filename, netcdf_mode_r)

        # Get the variables
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        z = fid.variables['elevation'][:]
        time = fid.variables['time'][:]
        stage = fid.variables['stage'][:]


        #Export to ascii/prj files
        sww2dem(self.domain.get_name() + '.sww',
                self.domain.get_name() + '_stage.asc',
                quantity = 'stage',
                cellsize = cellsize,
                number_of_decimal_places = 9,
                reduction = min)


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
        assert num.allclose(float(L[1].strip().lower()), 308500)

        L = lines[3].strip().split()
        assert L[0].strip().lower() == 'yllcorner'
        assert num.allclose(float(L[1].strip().lower()), 6189000)

        L = lines[4].strip().split()
        assert L[0].strip().lower() == 'cellsize'
        assert num.allclose(float(L[1].strip().lower()), cellsize)

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
                        assert num.allclose(float(L[i]), min(val0, val1))


        fid.close()

        #Cleanup
        os.remove(prjfile)
        os.remove(ascfile)
        os.remove(swwfile)

    def test_sww2dem_asc_stage_time(self):
        """Test that sww information can be converted correctly to asc/prj
        format readable by e.g. ArcView

        This tests the reduction of quantity stage using min
        """

        import time, os

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

        sww = SWW_file(self.domain)
        sww.store_connectivity()
        sww.store_timestep()

        #self.domain.tight_slope_limiters = 1
        self.domain.evolve_to_end(finaltime = 0.01)
        sww.store_timestep()

        cellsize = 0.25
        #Check contents
        #Get NetCDF

        fid = NetCDFFile(sww.filename, netcdf_mode_r)

        # Get the variables
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        z = fid.variables['elevation'][:]
        time = fid.variables['time'][:]
        stage = fid.variables['stage'][:]

        #Export to ascii/prj files
        sww2dem(self.domain.get_name() + '.sww',
                self.domain.get_name() + '_stage.asc',
                quantity = 'stage',
                cellsize = cellsize,
                number_of_decimal_places = 9,
                reduction = 1,
                verbose=self.verbose)


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
        assert num.allclose(float(L[1].strip().lower()), 308500)

        L = lines[3].strip().split()
        assert L[0].strip().lower() == 'yllcorner'
        assert num.allclose(float(L[1].strip().lower()), 6189000)

        L = lines[4].strip().split()
        assert L[0].strip().lower() == 'cellsize'
        assert num.allclose(float(L[1].strip().lower()), cellsize)

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
                        
                        val = stage[1,index]
                   
                        assert num.allclose(float(L[i]), val)

        fid.close()

        #Cleanup
        os.remove(prjfile)
        os.remove(ascfile)
        os.remove(swwfile)


    def test_sww2dem_asc_derived_quantity(self):
        """Test that sww information can be converted correctly to asc/prj
        format readable by e.g. ArcView

        This tests the use of derived quantities
        """

        import time, os

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


        sww = SWW_file(self.domain)
        sww.store_connectivity()
        sww.store_timestep()

        #self.domain.tight_slope_limiters = 1
        self.domain.evolve_to_end(finaltime = 0.01)
        sww.store_timestep()

        cellsize = 0.25
        #Check contents
        #Get NetCDF

        fid = NetCDFFile(sww.filename, netcdf_mode_r)

        # Get the variables
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        z = fid.variables['elevation'][:]
        time = fid.variables['time'][:]
        stage = fid.variables['stage'][:]


        #Export to ascii/prj files
        sww2dem(self.domain.get_name()+'.sww',
                name_out = 'datatest_depth.asc',
                quantity = 'stage - elevation',
                cellsize = cellsize,
                number_of_decimal_places = 9,
                reduction = min,
                verbose = self.verbose)


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
        assert num.allclose(float(L[1].strip().lower()), 308500)

        L = lines[3].strip().split()
        assert L[0].strip().lower() == 'yllcorner'
        assert num.allclose(float(L[1].strip().lower()), 6189000)

        L = lines[4].strip().split()
        assert L[0].strip().lower() == 'cellsize'
        assert num.allclose(float(L[1].strip().lower()), cellsize)

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
                        assert num.allclose(float(L[i]), min(val0, val1))


        fid.close()

        #Cleanup
        os.remove(prjfile)
        os.remove(ascfile)
        os.remove(swwfile)





    def test_sww2dem_asc_missing_points(self):
        """Test that sww information can be converted correctly to asc/prj
        format readable by e.g. ArcView

        This test includes the writing of missing values
        """

        import time, os

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
        stage = num.zeros(bed.shape, num.float)

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

        sww = SWW_file(domain)
        sww.store_connectivity()
        sww.store_timestep()

        cellsize = 0.25
        #Check contents
        #Get NetCDF

        fid = NetCDFFile(swwfile, netcdf_mode_r)

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
        sww2dem(domain.get_name()+'.sww',
                domain.get_name()+'_elevation.asc',
                quantity = 'elevation',
                cellsize = cellsize,
                number_of_decimal_places = 9,
                verbose = self.verbose)


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
        assert num.allclose(float(L[1].strip().lower()), 308500)

        L = lines[3].strip().split()
        assert L[0].strip().lower() == 'yllcorner'
        assert num.allclose(float(L[1].strip().lower()), 6189000)

        L = lines[4].strip().split()
        assert L[0].strip().lower() == 'cellsize'
        assert num.allclose(float(L[1].strip().lower()), cellsize)

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
                    assert num.allclose(float(L[i]), -i*cellsize - y)
                else:
                    #Missing values
                    assert num.allclose(float(L[i]), -9999)



        fid.close()

        #Cleanup
        os.remove(prjfile)
        os.remove(ascfile)
        os.remove(swwfile)



    def test_sww2ers_simple(self):
        """Test that sww information can be converted correctly to ers
        format
        """

        import time, os


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

        sww = SWW_file(self.domain)
        sww.store_connectivity()
        sww.store_timestep()

        #self.domain.tight_slope_limiters = 1
        self.domain.evolve_to_end(finaltime = 0.01)
        sww.store_timestep()

        cellsize = 0.25
        #Check contents
        #Get NetCDF

        fid = NetCDFFile(sww.filename, netcdf_mode_r)

        # Get the variables
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        z = fid.variables['elevation'][:]
        time = fid.variables['time'][:]
        stage = fid.variables['stage'][:]


        #Export to ers files
        outname = self.domain.get_name() + '_elevation.ers'
        sww2dem(self.domain.get_name() + '.sww',
                outname,
                quantity = 'elevation',
                cellsize = cellsize,
                number_of_decimal_places = 9,
                NODATA_value = NODATA_value,
                verbose = self.verbose)

        #Check header data
        from anuga.abstract_2d_finite_volumes.ermapper_grids import read_ermapper_header, read_ermapper_data

        header = read_ermapper_header(outname)

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


        ref_grid = [-1,    -1.25, -1.5,  -1.75, -2.0,
                    -0.75, -1.0,  -1.25, -1.5,  -1.75,
                    -0.5,  -0.75, -1.0,  -1.25, -1.5,
                    -0.25, -0.5,  -0.75, -1.0,  -1.25,
                    -0.0,  -0.25, -0.5,  -0.75, -1.0]


        #print grid.reshape((5,5))
        assert num.allclose(grid, ref_grid)

        fid.close()

        #Cleanup
        #FIXME the file clean-up doesn't work (eg Permission Denied Error)
        #Done (Ole) - it was because sww2ers didn't close it's sww file
        os.remove(sww.filename)
        os.remove(self.domain.get_name() + '_elevation')
        os.remove(self.domain.get_name() + '_elevation.ers')
        
    def test_export_grid_parallel(self):
        """Test that sww information can be converted correctly to asc/prj
        format readable by e.g. ArcView
        """

        import time, os

        base_name = 'tegp'
        #Setup
        self.domain.set_name(base_name+'_P0_8')
        swwfile = self.domain.get_name() + '.sww'

        self.domain.set_datadir('.')
        self.domain.format = 'sww'
        self.domain.smooth = True
        self.domain.set_quantity('elevation', lambda x,y: -x-y)
        self.domain.set_quantity('stage', 1.0)

        self.domain.geo_reference = Geo_reference(56,308500,6189000)

        sww = SWW_file(self.domain)
        sww.store_connectivity()
        sww.store_timestep()
        self.domain.evolve_to_end(finaltime = 0.0001)
        #Setup
        self.domain.set_name(base_name+'_P1_8')
        swwfile2 = self.domain.get_name() + '.sww'
        sww = SWW_file(self.domain)
        sww.store_connectivity()
        sww.store_timestep()
        self.domain.evolve_to_end(finaltime = 0.0002)
        sww.store_timestep()

        cellsize = 0.25
        #Check contents
        #Get NetCDF

        fid = NetCDFFile(sww.filename, netcdf_mode_r)

        # Get the variables
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        z = fid.variables['elevation'][:]
        time = fid.variables['time'][:]
        stage = fid.variables['stage'][:]

        fid.close()

        #Export to ascii/prj files
        extra_name_out = 'yeah'
        sww2dem_batch(base_name,
                    quantities = ['elevation', 'depth'],
                    extra_name_out = extra_name_out,
                    cellsize = cellsize,
                    verbose = self.verbose,
                    format = 'asc')

        prjfile = base_name + '_P0_8_elevation_yeah.prj'
        ascfile = base_name + '_P0_8_elevation_yeah.asc'       
        #Check asc file
        ascid = open(ascfile)
        lines = ascid.readlines()
        ascid.close()
        #Check grid values
        for j in range(5):
            L = lines[6+j].strip().split()
            y = (4-j) * cellsize
            for i in range(5):
                #print " -i*cellsize - y",  -i*cellsize - y
                #print "float(L[i])", float(L[i])
                assert num.allclose(float(L[i]), -i*cellsize - y)               
        #Cleanup
        os.remove(prjfile)
        os.remove(ascfile)

        prjfile = base_name + '_P1_8_elevation_yeah.prj'
        ascfile = base_name + '_P1_8_elevation_yeah.asc'       
        #Check asc file
        ascid = open(ascfile)
        lines = ascid.readlines()
        ascid.close()
        #Check grid values
        for j in range(5):
            L = lines[6+j].strip().split()
            y = (4-j) * cellsize
            for i in range(5):
                #print " -i*cellsize - y",  -i*cellsize - y
                #print "float(L[i])", float(L[i])
                assert num.allclose(float(L[i]), -i*cellsize - y)               
        #Cleanup
        os.remove(prjfile)
        os.remove(ascfile)
        os.remove(swwfile)

        #Check asc file
        ascfile = base_name + '_P0_8_depth_yeah.asc'
        prjfile = base_name + '_P0_8_depth_yeah.prj'
        ascid = open(ascfile)
        lines = ascid.readlines()
        ascid.close()
        #Check grid values
        for j in range(5):
            L = lines[6+j].strip().split()
            y = (4-j) * cellsize
            for i in range(5):
                assert num.allclose(float(L[i]), 1 - (-i*cellsize - y))
        #Cleanup
        os.remove(prjfile)
        os.remove(ascfile)

        #Check asc file
        ascfile = base_name + '_P1_8_depth_yeah.asc'
        prjfile = base_name + '_P1_8_depth_yeah.prj'
        ascid = open(ascfile)
        lines = ascid.readlines()
        ascid.close()
        #Check grid values
        for j in range(5):
            L = lines[6+j].strip().split()
            y = (4-j) * cellsize
            for i in range(5):
                assert num.allclose(float(L[i]), 1 - (-i*cellsize - y))
        #Cleanup
        os.remove(prjfile)
        os.remove(ascfile)
        os.remove(swwfile2)


    def test_export_grid(self):
        """
        test_export_grid(self):
        Test that sww information can be converted correctly to asc/prj
        format readable by e.g. ArcView
        """

        import time, os

        try:
            os.remove('teg*.sww')
        except:
            pass

        #Setup
        self.domain.set_name('teg')

        prjfile = self.domain.get_name() + '_elevation.prj'
        ascfile = self.domain.get_name() + '_elevation.asc'
        swwfile = self.domain.get_name() + '.sww'

        self.domain.set_datadir('.')
        self.domain.smooth = True
        self.domain.set_quantity('elevation', lambda x,y: -x-y)
        self.domain.set_quantity('stage', 1.0)

        self.domain.geo_reference = Geo_reference(56,308500,6189000)

        sww = SWW_file(self.domain)
        sww.store_connectivity()
        sww.store_timestep()
        self.domain.evolve_to_end(finaltime = 0.01)
        sww.store_timestep()

        cellsize = 0.25
        #Check contents
        #Get NetCDF

        fid = NetCDFFile(sww.filename, netcdf_mode_r)

        # Get the variables
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        z = fid.variables['elevation'][:]
        time = fid.variables['time'][:]
        stage = fid.variables['stage'][:]

        fid.close()

        #Export to ascii/prj files
        sww2dem_batch(self.domain.get_name(),
                quantities = 'elevation',
                cellsize = cellsize,
                verbose = self.verbose,
                format = 'asc')

        #Check asc file
        ascid = open(ascfile)
        lines = ascid.readlines()
        ascid.close()

        L = lines[2].strip().split()
        assert L[0].strip().lower() == 'xllcorner'
        assert num.allclose(float(L[1].strip().lower()), 308500)

        L = lines[3].strip().split()
        assert L[0].strip().lower() == 'yllcorner'
        assert num.allclose(float(L[1].strip().lower()), 6189000)

        #Check grid values
        #print '==='
        for j in range(5):
            L = lines[6+j].strip().split()
            y = (4-j) * cellsize
            for i in range(5):
                #print float(L[i])
                assert num.allclose(float(L[i]), -i*cellsize - y)
                
        #Cleanup
        os.remove(prjfile)
        os.remove(ascfile)
        os.remove(swwfile)

    def test_export_gridII(self):
        """
        test_export_gridII(self):
        Test that sww information can be converted correctly to asc/prj
        format readable by e.g. ArcView
        """

        import time, os

        try:
            os.remove('teg*.sww')
        except:
            pass

        #Setup
        self.domain.set_name('tegII')

        swwfile = self.domain.get_name() + '.sww'

        self.domain.set_datadir('.')
        self.domain.smooth = True
        self.domain.set_quantity('elevation', lambda x,y: -x-y)
        self.domain.set_quantity('stage', 1.0)

        self.domain.geo_reference = Geo_reference(56,308500,6189000)

        sww = SWW_file(self.domain)
        sww.store_connectivity()
        sww.store_timestep()
        self.domain.evolve_to_end(finaltime = 0.01)
        sww.store_timestep()

        cellsize = 0.25
        #Check contents
        #Get NetCDF

        fid = NetCDFFile(sww.filename, netcdf_mode_r)

        # Get the variables
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        z = fid.variables['elevation'][:]
        time = fid.variables['time'][:]
        stage = fid.variables['stage'][:]
        xmomentum = fid.variables['xmomentum'][:]
        ymomentum = fid.variables['ymomentum'][:]        

        #print 'stage', stage
        #print 'xmom', xmomentum
        #print 'ymom', ymomentum        

        fid.close()

        #Export to ascii/prj files
        if True:
            sww2dem_batch(self.domain.get_name(),
                        quantities = ['elevation', 'depth'],
                        cellsize = cellsize,
                        verbose = self.verbose,
                        format = 'asc')

        else:
            sww2dem_batch(self.domain.get_name(),
                quantities = ['depth'],
                cellsize = cellsize,
                verbose = self.verbose,
                format = 'asc')


            export_grid(self.domain.get_name(),
                quantities = ['elevation'],
                cellsize = cellsize,
                verbose = self.verbose,
                format = 'asc')

        prjfile = self.domain.get_name() + '_elevation.prj'
        ascfile = self.domain.get_name() + '_elevation.asc'
        
        #Check asc file
        ascid = open(ascfile)
        lines = ascid.readlines()
        ascid.close()

        L = lines[2].strip().split()
        assert L[0].strip().lower() == 'xllcorner'
        assert num.allclose(float(L[1].strip().lower()), 308500)

        L = lines[3].strip().split()
        assert L[0].strip().lower() == 'yllcorner'
        assert num.allclose(float(L[1].strip().lower()), 6189000)

        #print "ascfile", ascfile
        #Check grid values
        for j in range(5):
            L = lines[6+j].strip().split()
            y = (4-j) * cellsize
            for i in range(5):
                #print " -i*cellsize - y",  -i*cellsize - y
                #print "float(L[i])", float(L[i])
                assert num.allclose(float(L[i]), -i*cellsize - y)

        #Cleanup
        os.remove(prjfile)
        os.remove(ascfile)
        
        #Check asc file
        ascfile = self.domain.get_name() + '_depth.asc'
        prjfile = self.domain.get_name() + '_depth.prj'
        ascid = open(ascfile)
        lines = ascid.readlines()
        ascid.close()

        L = lines[2].strip().split()
        assert L[0].strip().lower() == 'xllcorner'
        assert num.allclose(float(L[1].strip().lower()), 308500)

        L = lines[3].strip().split()
        assert L[0].strip().lower() == 'yllcorner'
        assert num.allclose(float(L[1].strip().lower()), 6189000)

        #Check grid values
        for j in range(5):
            L = lines[6+j].strip().split()
            y = (4-j) * cellsize
            for i in range(5):
                #print " -i*cellsize - y",  -i*cellsize - y
                #print "float(L[i])", float(L[i])                
                assert num.allclose(float(L[i]), 1 - (-i*cellsize - y))

        #Cleanup
        os.remove(prjfile)
        os.remove(ascfile)
        os.remove(swwfile)


    def test_export_gridIII(self):
        """
        test_export_gridIII
        Test that sww information can be converted correctly to asc/prj
        format readable by e.g. ArcView
        """

        import time, os

        try:
            os.remove('teg*.sww')
        except:
            pass

        #Setup
        
        self.domain.set_name('tegIII')

        swwfile = self.domain.get_name() + '.sww'

        self.domain.set_datadir('.')
        self.domain.format = 'sww'
        self.domain.smooth = True
        self.domain.set_quantity('elevation', lambda x,y: -x-y)
        self.domain.set_quantity('stage', 1.0)

        self.domain.geo_reference = Geo_reference(56,308500,6189000)
        
        sww = SWW_file(self.domain)
        sww.store_connectivity()
        sww.store_timestep() #'stage')
        self.domain.evolve_to_end(finaltime = 0.01)
        sww.store_timestep() #'stage')

        cellsize = 0.25
        #Check contents
        #Get NetCDF

        fid = NetCDFFile(sww.filename, netcdf_mode_r)

        # Get the variables
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        z = fid.variables['elevation'][:]
        time = fid.variables['time'][:]
        stage = fid.variables['stage'][:]

        fid.close()

        #Export to ascii/prj files
        extra_name_out = 'yeah'
        if True:
            sww2dem_batch(self.domain.get_name(),
                        quantities = ['elevation', 'depth'],
                        extra_name_out = extra_name_out,
                        cellsize = cellsize,
                        verbose = self.verbose,
                        format = 'asc')

        else:
            sww2dem_batch(self.domain.get_name(),
                quantities = ['depth'],
                cellsize = cellsize,
                verbose = self.verbose,
                format = 'asc')


            sww2dem_batch(self.domain.get_name(),
                quantities = ['elevation'],
                cellsize = cellsize,
                verbose = self.verbose,
                format = 'asc')

        prjfile = self.domain.get_name() + '_elevation_yeah.prj'
        ascfile = self.domain.get_name() + '_elevation_yeah.asc'
        
        #Check asc file
        ascid = open(ascfile)
        lines = ascid.readlines()
        ascid.close()

        L = lines[2].strip().split()
        assert L[0].strip().lower() == 'xllcorner'
        assert num.allclose(float(L[1].strip().lower()), 308500)

        L = lines[3].strip().split()
        assert L[0].strip().lower() == 'yllcorner'
        assert num.allclose(float(L[1].strip().lower()), 6189000)

        #print "ascfile", ascfile
        #Check grid values
        for j in range(5):
            L = lines[6+j].strip().split()
            y = (4-j) * cellsize
            for i in range(5):
                #print " -i*cellsize - y",  -i*cellsize - y
                #print "float(L[i])", float(L[i])
                assert num.allclose(float(L[i]), -i*cellsize - y)
                
        #Cleanup
        os.remove(prjfile)
        os.remove(ascfile)
        
        #Check asc file
        ascfile = self.domain.get_name() + '_depth_yeah.asc'
        prjfile = self.domain.get_name() + '_depth_yeah.prj'
        ascid = open(ascfile)
        lines = ascid.readlines()
        ascid.close()

        L = lines[2].strip().split()
        assert L[0].strip().lower() == 'xllcorner'
        assert num.allclose(float(L[1].strip().lower()), 308500)

        L = lines[3].strip().split()
        assert L[0].strip().lower() == 'yllcorner'
        assert num.allclose(float(L[1].strip().lower()), 6189000)

        #Check grid values
        for j in range(5):
            L = lines[6+j].strip().split()
            y = (4-j) * cellsize
            for i in range(5):
                assert num.allclose(float(L[i]), 1 - (-i*cellsize - y))

        #Cleanup
        os.remove(prjfile)
        os.remove(ascfile)
        os.remove(swwfile)

    def test_export_grid_bad(self):
        """Test that Bad input throws exception error
        """

        try:
            sww2dem_batch('a_small_round-egg',
                        quantities = ['elevation', 'depth'],
                        cellsize = 99,
                        verbose = self.verbose,
                        format = 'asc')
        except IOError:
            pass
        else:
            self.assertTrue(0 ==1,  'Bad input did not throw exception error!')
        

#################################################################################

if __name__ == "__main__":
    #suite = unittest.makeSuite(Test_Shallow_Water, 'test_rainfall_forcing_with_evolve')
    suite = unittest.makeSuite(Test_Sww2Dem, 'test_')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)

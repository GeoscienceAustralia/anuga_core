import unittest
import copy
import os
import numpy as num

from anuga.config import netcdf_float

from anuga.coordinate_transforms.geo_reference import Geo_reference 
from anuga.geometry.polygon import is_inside_polygon
from anuga.abstract_2d_finite_volumes.util import file_function
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a
from anuga.config import g

from anuga.shallow_water.boundaries import Reflective_boundary, \
            Field_boundary, Transmissive_momentum_set_stage_boundary, \
            Transmissive_stage_zero_momentum_boundary
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions\
     import Transmissive_boundary, Dirichlet_boundary, \
            Time_boundary, File_boundary, AWI_boundary

from anuga.file.sww import get_mesh_and_quantities_from_file
            
from anuga.shallow_water.shallow_water_domain import Domain

from anuga.abstract_2d_finite_volumes.mesh_factory \
        import rectangular_cross
from anuga.shallow_water.sww_interrogate import \
            get_maximum_inundation_elevation, \
            get_maximum_inundation_location, get_maximum_inundation_data, \
            get_flow_through_cross_section, get_energy_through_cross_section
            
            
from anuga.file_conversion.dem2dem import dem2dem
                

class Test_Dem2Dem(unittest.TestCase):
    def test_decimate_dem(self):
        """Test decimation of dem file
        """

        import os
        from anuga.file.netcdf import NetCDFFile

        #Write test dem file
        root = 'decdemtest'

        filename = root + '.dem'
        fid = NetCDFFile(filename, netcdf_mode_w)

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

        fid.createVariable('elevation', netcdf_float, ('number_of_points',))

        elevation = fid.variables['elevation']

        elevation[:] = (num.arange(nrows*ncols))

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

        # generate a stencil for computing the decimated values
        stencil = num.ones((3,3), num.float) / 9.0

        dem2dem(filename, stencil=stencil, cellsize_new=100)

        # Open decimated NetCDF file
        fid = NetCDFFile(root + '_100.dem', netcdf_mode_r)

        # Get decimated elevation
        elevation = fid.variables['elevation']

        # Check values
        assert num.allclose(elevation, ref_elevation)

        # Cleanup
        fid.close()

        os.remove(root + '.dem')
        os.remove(root + '_100.dem')

    def test_decimate_dem_NODATA(self):
        """Test decimation of dem file that includes NODATA values
        """

        import os
        from anuga.file.netcdf import NetCDFFile

        # Write test dem file
        root = 'decdemtest'

        filename = root + '.dem'
        fid = NetCDFFile(filename, netcdf_mode_w)

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

        fid.createVariable('elevation', netcdf_float, ('number_of_points',))

        elevation = fid.variables['elevation']

        # Generate initial elevation values
        elevation_tmp = (num.arange(nrows*ncols))

        # Add some NODATA values
        elevation_tmp[0]   = NODATA_value
        elevation_tmp[95]  = NODATA_value
        elevation_tmp[188] = NODATA_value
        elevation_tmp[189] = NODATA_value
        elevation_tmp[190] = NODATA_value
        elevation_tmp[209] = NODATA_value
        elevation_tmp[252] = NODATA_value

        elevation[:] = elevation_tmp

        fid.close()

        # Generate the elevation values expected in the decimated file
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

        # Generate a stencil for computing the decimated values
        stencil = num.ones((3,3), num.float) / 9.0

        dem2dem(filename, stencil=stencil, cellsize_new=100)

        # Open decimated NetCDF file
        fid = NetCDFFile(root + '_100.dem', netcdf_mode_r)

        # Get decimated elevation
        elevation = fid.variables['elevation']

        # Check values
        assert num.allclose(elevation, ref_elevation)

        # Cleanup
        fid.close()

        os.remove(root + '.dem')
        os.remove(root + '_100.dem')      

#################################################################################

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_Dem2Dem, 'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
        

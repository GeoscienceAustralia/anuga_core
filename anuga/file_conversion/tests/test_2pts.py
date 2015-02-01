# external modules
import unittest
import tempfile
import shutil
import numpy as num

# ANUGA modules
from anuga.shallow_water.shallow_water_domain import Domain 
from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.file.sww import Write_sww, SWW_file
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions \
                            import Transmissive_boundary
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a, \
                            netcdf_float
from anuga.geospatial_data.geospatial_data import Geospatial_data

# local modules
from anuga.file_conversion.sdf2pts import sdf2pts
from anuga.file_conversion.sww2pts import sww2pts


class Test_2Pts(unittest.TestCase):
    """ Test files that convert to pts format. """
    
    def test_hecras_cross_sections2pts(self):
        """Test conversion from HECRAS cross sections in ascii format
        to native NetCDF pts format
        """

        import time, os
        from anuga.file.netcdf import NetCDFFile

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
        sdf2pts(root+'.sdf')

        #Check contents
        #Get NetCDF
        fid = NetCDFFile(root+'.pts', netcdf_mode_r)

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
        assert num.allclose(points, ref_points)

        #print attributes[:]
        #print ref_elevation
        assert num.allclose(elevation, ref_elevation)

        #Cleanup
        fid.close()


        os.remove(root + '.sdf')
        os.remove(root + '.pts')




    def test_sww2pts_centroids(self):
        """Test that sww information can be converted correctly to pts data at specified coordinates
        - in this case, the centroids.
        """

        import time, os
        from anuga.file.netcdf import NetCDFFile
        # Used for points that lie outside mesh
        NODATA_value = 1758323

        # Setup
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create shallow water domain
        domain = Domain(*rectangular(2, 2))

        B = Transmissive_boundary(domain)
        domain.set_boundary( {'left': B, 'right': B, 'top': B, 'bottom': B})

        domain.set_name('datatest')

        ptsfile = domain.get_name() + '_elevation.pts'
        swwfile = domain.get_name() + '.sww'

        domain.set_datadir('.')
        domain.format = 'sww'
        domain.set_quantity('elevation', lambda x,y: -x-y)

        domain.geo_reference = Geo_reference(56,308500,6189000)

        sww = SWW_file(domain)
        sww.store_connectivity()
        sww.store_timestep()

        #self.domain.tight_slope_limiters = 1
        domain.evolve_to_end(finaltime = 0.01)
        sww.store_timestep()

        # Check contents in NetCDF
        fid = NetCDFFile(sww.filename, netcdf_mode_r)

        # Get the variables
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        elevation = fid.variables['elevation'][:]
        time = fid.variables['time'][:]
        stage = fid.variables['stage'][:]

        volumes = fid.variables['volumes'][:]


        # Invoke interpolation for vertex points       
        points = num.concatenate( (x[:,num.newaxis],y[:,num.newaxis]), axis=1 )
        points = num.ascontiguousarray(points)
        sww2pts(domain.get_name() + '.sww',
                quantity = 'elevation',
                data_points = points,
                NODATA_value = NODATA_value)
        ref_point_values = elevation
        point_values = Geospatial_data(ptsfile).get_attributes()
        #print 'P', point_values
        #print 'Ref', ref_point_values        
        assert num.allclose(point_values, ref_point_values)        



        # Invoke interpolation for centroids
        points = domain.get_centroid_coordinates()
        #print points
        sww2pts(domain.get_name() + '.sww',
                quantity = 'elevation',
                data_points = points,
                NODATA_value = NODATA_value)
        ref_point_values = [-0.5, -0.5, -1, -1, -1, -1, -1.5, -1.5]   #At centroids

        
        point_values = Geospatial_data(ptsfile).get_attributes()
        #print 'P', point_values
        #print 'Ref', ref_point_values        
        assert num.allclose(point_values, ref_point_values)        

        fid.close()

        #Cleanup
        os.remove(sww.filename)
        os.remove(ptsfile)


#-------------------------------------------------------------

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_2Pts, 'test_sww')
    runner = unittest.TextTestRunner() #verbosity=2)
    runner.run(suite)    

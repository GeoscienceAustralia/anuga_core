from builtins import str
import unittest
import copy
import os
import numpy as num


from anuga.file_conversion.tif2array import tif2array
from anuga.file_conversion.tif2point_values import tif2point_values


def make_a_utm_tif():
    # We need to make a .tif to test some functions
    # This does the job
    #
    from anuga.utilities import plot_utils as util
    import numpy
    #
    # Use Make_Geotif to make tif file
    # Pick a domain that makes sense in EPSG:32756
    # WGS 84 / UTM zone 56 South
    x = numpy.linspace(307000., 307100., 101)
    y = numpy.linspace(6193000., 6193100., 101)
    xG, yG = numpy.meshgrid(x, y)
    xG = xG.flatten()
    yG = yG.flatten()
    # Surface is z=x+y
    fakeZ = xG-min(xG)+yG - min(yG)
    dataToGrid = numpy.vstack([xG, yG, fakeZ]).transpose()
    #
    util.Make_Geotif(dataToGrid, output_quantities=['test_utm'],
                     EPSG_CODE=32756, output_dir='.', CellSize=1.0)


def make_a_ll_tif():
    # We need to make a .tif with ll coord to test some functions
    #
    from anuga.utilities import plot_utils as util
    import numpy
    #
    # Use Make_Geotif to make tif file
    # Pick a domain that makes sense in EPSG:32756
    lat = numpy.linspace(-34.39, -34.37, 101)
    long = numpy.linspace(150.90, 150.92, 101)

    xG, yG = numpy.meshgrid(long, lat)
    xG = xG.flatten()
    yG = yG.flatten()
    # Surface is z=x+y

    fakeZ = (xG-min(xG))/(max(xG)-min(xG))+(yG - min(yG))/(max(yG)-min(yG))
    dataToGrid = numpy.vstack([xG, yG, fakeZ]).transpose()
    #
    # print(dataToGrid.shape)
    util.Make_Geotif(dataToGrid, output_quantities=['test_ll'],
                     EPSG_CODE=4326, output_dir='.', CellSize=0.0001)


class Test_tif2(unittest.TestCase):

    def test_tif2array(self):

        import os

        # makes a file Point_Data_test_utm.tif
        # which contains UTM easting, northing data for 
        # WGS 84 / UTM zone 56 South
        make_a_utm_tif()

        os.remove('PointData_test_utm.tif')


    def test_tif_lat_lon(self):

        import os

        # makes a file Point_Data_test_ll.tif
        # which contains lat lon data
        make_a_ll_tif()

        os.remove('PointData_test_ll.tif')




#################################################################################
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_tif2, 'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)

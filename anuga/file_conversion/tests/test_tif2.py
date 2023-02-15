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
    # Pick a domain that makes sense in EPSG:32756 (zone 56, south=True)
    lat = numpy.linspace(-34.39, -34.37, 101)
    lon = numpy.linspace(150.90, 150.92, 101)

    xG, yG = numpy.meshgrid(lon, lat)
    xG = xG.flatten()
    yG = yG.flatten()
    # Surface is z=x+y

    fakeZ = (xG-min(xG))/(max(xG)-min(xG))+(yG - min(yG))/(max(yG)-min(yG))
    dataToGrid = numpy.vstack([xG, yG, fakeZ]).transpose()
    #
    # Create file PointData_test_ll.tif
    util.Make_Geotif(dataToGrid, output_quantities=['test_ll'],
                     EPSG_CODE=4326, output_dir='.', CellSize=0.0001)


class Test_tif2(unittest.TestCase):

    def test_tif2array(self):

        import os
        import numpy
        from anuga.file_conversion.tif2array import  tif2array

        # makes a file Point_Data_test_utm.tif
        # which contains UTM easting, northing data for 
        # WGS 84 / UTM zone 56 South
        make_a_utm_tif()

        x, y, Z = tif2array('PointData_test_utm.tif')

        x_exact = numpy.array([306999.5, 307000.5, 307001.5, 307002.5, 307003.5, 307004.5,
                               307005.5, 307006.5, 307007.5, 307008.5, 307009.5, 307010.5,
                               307011.5, 307012.5, 307013.5, 307014.5, 307015.5, 307016.5,
                               307017.5, 307018.5, 307019.5, 307020.5, 307021.5, 307022.5,
                               307023.5, 307024.5, 307025.5, 307026.5, 307027.5, 307028.5,
                               307029.5, 307030.5, 307031.5, 307032.5, 307033.5, 307034.5,
                               307035.5, 307036.5, 307037.5, 307038.5, 307039.5, 307040.5,
                               307041.5, 307042.5, 307043.5, 307044.5, 307045.5, 307046.5,
                               307047.5, 307048.5, 307049.5, 307050.5, 307051.5, 307052.5,
                               307053.5, 307054.5, 307055.5, 307056.5, 307057.5, 307058.5,
                               307059.5, 307060.5, 307061.5, 307062.5, 307063.5, 307064.5,
                               307065.5, 307066.5, 307067.5, 307068.5, 307069.5, 307070.5,
                               307071.5, 307072.5, 307073.5, 307074.5, 307075.5, 307076.5,
                               307077.5, 307078.5, 307079.5, 307080.5, 307081.5, 307082.5,
                               307083.5, 307084.5, 307085.5, 307086.5, 307087.5, 307088.5,
                               307089.5, 307090.5, 307091.5, 307092.5, 307093.5, 307094.5,
                               307095.5, 307096.5, 307097.5, 307098.5, 307099.5])

        y_exact = numpy.array([6193000.5, 6193001.5, 6193002.5, 6193003.5, 6193004.5, 6193005.5,
                               6193006.5, 6193007.5, 6193008.5, 6193009.5, 6193010.5, 6193011.5,
                               6193012.5, 6193013.5, 6193014.5, 6193015.5, 6193016.5, 6193017.5,
                               6193018.5, 6193019.5, 6193020.5, 6193021.5, 6193022.5, 6193023.5,
                               6193024.5, 6193025.5, 6193026.5, 6193027.5, 6193028.5, 6193029.5,
                               6193030.5, 6193031.5, 6193032.5, 6193033.5, 6193034.5, 6193035.5,
                               6193036.5, 6193037.5, 6193038.5, 6193039.5, 6193040.5, 6193041.5,
                               6193042.5, 6193043.5, 6193044.5, 6193045.5, 6193046.5, 6193047.5,
                               6193048.5, 6193049.5, 6193050.5, 6193051.5, 6193052.5, 6193053.5,
                               6193054.5, 6193055.5, 6193056.5, 6193057.5, 6193058.5, 6193059.5,
                               6193060.5, 6193061.5, 6193062.5, 6193063.5, 6193064.5, 6193065.5,
                               6193066.5, 6193067.5, 6193068.5, 6193069.5, 6193070.5, 6193071.5,
                               6193072.5, 6193073.5, 6193074.5, 6193075.5, 6193076.5, 6193077.5,
                               6193078.5, 6193079.5, 6193080.5, 6193081.5, 6193082.5, 6193083.5,
                               6193084.5, 6193085.5, 6193086.5, 6193087.5, 6193088.5, 6193089.5,
                               6193090.5, 6193091.5, 6193092.5, 6193093.5, 6193094.5, 6193095.5,
                               6193096.5, 6193097.5, 6193098.5, 6193099.5, 6193100.5])

        Z_row_11 = numpy.array([11.,  12.,  13.,  14.,  15.,  16.,  17.,  18.,  19.,  20.,  21.,
                                22.,  23.,  24.,  25.,  26.,  27.,  28.,  29.,  30.,  31.,  32.,
                                33.,  34.,  35.,  36.,  37.,  38.,  39.,  40.,  41.,  42.,  43.,
                                44.,  45.,  46.,  47.,  48.,  49.,  50.,  51.,  52.,  53.,  54.,
                                55.,  56.,  57.,  58.,  59.,  60.,  61.,  62.,  63.,  64.,  65.,
                                66.,  67.,  68.,  69.,  70.,  71.,  72.,  73.,  74.,  75.,  76.,
                                77.,  78.,  79.,  80.,  81.,  82.,  83.,  84.,  85.,  86.,  87.,
                                88.,  89.,  90.,  91.,  92.,  93.,  94.,  95.,  96.,  97.,  98.,
                                99., 100., 101., 102., 103., 104., 105., 106., 107., 108., 109.,
                                110., 111.])

        assert numpy.allclose(x,x_exact)
        assert numpy.allclose(y,y_exact)

        assert numpy.allclose(Z[11,:],Z_row_11)


        os.remove('PointData_test_utm.tif')

        
    def test_tif2point_values_ll(self):

        import os
        import numpy
        from anuga.file_conversion.tif2point_values import tif2point_values

        # makes a file Point_Data_test_ll.tif
        # which contains lat lon data
        make_a_ll_tif()

        x = numpy.linspace(307000., 308000., 11)
        y = numpy.linspace(6193000., 6194000., 11)
        xG, yG = numpy.meshgrid(x, y)
        xG = xG.flatten()
        yG = yG.flatten()
        points = numpy.vstack((xG,yG)).T
        Z = tif2point_values('PointData_test_ll.tif', zone=56, south=True, points = points)

        Z_exact = numpy.array([0.21774116, 0.27, 0.32333332, 0.3740863, 0.42774117,
                               0.48333332, 0.5340863, 0.5933333, 0.64, 0.69225883,
                               0.75, 0.25666666, 0.3159137, 0.36774117, 0.42,
                               0.47591373, 0.53, 0.58, 0.63774115, 0.68774116,
                               0.7366667, 0.79225886, 0.30408627, 0.36, 0.41,
                               0.46225885, 0.52, 0.5759137, 0.62666667, 0.6766667,
                               0.73, 0.78225887, 0.84, 0.35774115, 0.40774116,
                               0.45774117, 0.51, 0.56591374, 0.62, 0.67225885,
                               0.73, 0.78333336, 0.83, 0.88408625, 0.4,
                               0.45591372, 0.50225884, 0.55225885, 0.6066667, 0.6659137,
                               0.72, 0.7722588, 0.8240863, 0.88, 0.93225884,
                               0.44225883, 0.49666667, 0.55225885, 0.60225886, 0.66,
                               0.71, 0.7622588, 0.82, 0.87, 0.92225885,
                               0.9766667, 0.49, 0.54225886, 0.6, 0.65,
                               0.70774114, 0.76, 0.81333333, 0.8640863, 0.9140863,
                               0.97, 1.0222589, 0.54, 0.58666664, 0.6422588,
                               0.69408625, 0.75, 0.80225885, 0.86, 0.91225886,
                               0.9633333, 1.0140862, 1.07, 0.5840863, 0.64,
                               0.69, 0.7366667, 0.79225886, 0.8466667, 0.9059137,
                               0.95666665, 1.01, 1.0677412, 1.1122588, 0.6333333,
                               0.68774116, 0.74, 0.79, 0.84, 0.8922588,
                               0.95, 1.0040863, 1.0540863, 1.1133333, 1.1677412,
                               0.6740863, 0.73, 0.7840863, 0.83774114, 0.8933333,
                               0.94408625, 0.99774116, 1.0533333, 1.1, 1.1559137,
                               1.21])

        assert numpy.allclose(Z, Z_exact)

        os.remove('PointData_test_ll.tif')

    def test_tif_lat_lon_too_small(self):

        import os
        import numpy
        from anuga.file_conversion.tif2point_values import tif2point_values

        # makes a file Point_Data_test_ll.tif
        # which contains lat lon data
        make_a_ll_tif()

        x = numpy.linspace(307000., 310000., 11)
        y = numpy.linspace(6190000., 6200000., 11)
        xG, yG = numpy.meshgrid(x, y)
        xG = xG.flatten()
        yG = yG.flatten()
        points = numpy.vstack((xG,yG)).T

        try:
            Z = tif2point_values('PointData_test_ll.tif', zone=56, south=True, points = points) 
        except ValueError:
            pass
        else:
            #Expected ValueError
            raise Exception()       




#################################################################################
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_tif2, 'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)

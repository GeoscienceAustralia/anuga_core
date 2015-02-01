import unittest
import numpy as num
from anuga.file.netcdf import NetCDFFile
from  anuga.shallow_water.most2nc import most2nc
import os

FN = 'small___.txt'

class Test_most2nc(unittest.TestCase):
    def setUp(self):
        fid = open(FN, 'w')
        fid.write("""4 4 
150.66667
150.83334
151.
151.16667
-34.
-34.16667
-34.33333
-34.5
-1. -2. -3. -4.
-5. -6. -7. -8.
-9. -10. -11. -12.
-13. -14. -15. -16.
""")
        fid.close()
                  
    def tearDown(self):
        os.remove(FN)

    def test_small_nxn(self):
        most2nc(input_file=FN,output_file='test.nc'\
                        ,inverted_bathymetry = False,verbose = False)

        fid = NetCDFFile('test.nc')
        elevation = fid.variables['ELEVATION'][:]
        fid.close()

        z=[[-13., -14., -15., -16.]\
           ,[-9., -10., -11., -12.]\
           ,[-5.,  -6.,  -7.,  -8.]\
           ,[-1.,  -2.,  -3.,  -4.]]
        z = num.asarray(z)

        assert num.allclose(z,elevation)
        import os
        os.remove('test.nc')
        
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_most2nc,'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)

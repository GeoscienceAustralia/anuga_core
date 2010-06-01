import os
import unittest
import tempfile
import numpy as num

from csv_file import load_csv_as_array, load_csv_as_dict


class Test_csv(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass
        
    def test_get_data_from_file1(self):
        fileName = tempfile.mktemp(".txt")
#        print"filename",fileName
        file = open(fileName,"w")
        file.write("elevation stage\n\
1.3 3  \n\
0.0 4 \n\
4.5 3.5 \n\
1.0 6 \n")
        file.close()

        x = load_csv_as_array(fileName, delimiter=' ')  
        
       # header, x = load_csv_as_array(fileName, delimiter=' ')
        os.remove(fileName)

        assert num.allclose(x['elevation'], [1.3, 0.0,4.5, 1.0])
        assert num.allclose(x['stage'], [3.0, 4.0,3.5, 6.0])        
        

#################################################################################

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_csv, 'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)

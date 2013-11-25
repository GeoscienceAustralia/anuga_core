import unittest
import numpy as num

from anuga.utilities.function_utils import determine_function_type


class Test_Function_Utils(unittest.TestCase):
                
    def test_determine_function_type_scalar(self):

        type = determine_function_type(1.0)

        assert type == 'scalar'


        type = determine_function_type(1)

        assert type == 'scalar'


    def test_determine_function_type_exception(self):

        
        type = determine_function_type([])

        assert type == 'array'

    def test_determine_function_type_time_only(self):


        def myfunc(t):
             return t*2

        type = determine_function_type(myfunc)

        assert type == 't'


    def test_determine_function_type_spatial_only(self):


        def myfunc(x,y):
             return x+y

        type = determine_function_type(myfunc)

        assert type == 'x,y'

    def test_determine_function_type_spatial_only(self):


        def myfunc(x,y,t):
             return x+y+t

        type = determine_function_type(myfunc)

        assert type == 'x,y,t'

        
#-------------------------------------------------------------

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_Function_Utils, 'test')
    runner = unittest.TextTestRunner() #verbosity=2)
    runner.run(suite)    

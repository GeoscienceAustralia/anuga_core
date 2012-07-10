"""Automatic verification that the ANUGA code validates against the okushiri
dataset as expected. See anuga_validation/okushiri_2005 for more details
"""

import unittest
import os


class Test_Okushiri(unittest.TestCase):
    def setUp(self):
        
        # Remove garbage
        for file in os.listdir('.'):
            if file.endswith('.stdout') or\
                    file.endswith('.sww') or\
                    file.endswith('.msh'):
                os.remove(file)
                
                
    def tearDown(self):
        pass

    def test_that_output_is_as_expected(self):
        """Test that ANUGA replicates physics of the Okushiri Island
        wave tank experiment
        """
        
        #print
        s = 'create_okushiri.py'
        #print s
        res = os.system('python %s > create_okushiri.stdout' %s)
        assert res == 0

        s = 'run_okushiri.py'
        #print s        
        res = os.system('python %s > run_okushiri.stdout' %s)        
        assert res == 0
        
        
        s = 'compare_timeseries_with_measures.py'
        #print s
        res = os.system('python %s > compare_timeseries_with_measures.stdout'\
                        %s)
        assert res == 0


    def test_caching_of_set_quantity(self):
        """Test that caching of set_quantity works
        """
        
        s = 'create_okushiri.py'
        res = os.system('python %s > create_okushiri_for_caching.stdout' %s)
        assert res == 0


        s = 'test_caching_of_set_quantity.py'
        res = os.system('python %s > test_caching_of_set_quantity.stdout' %s)
        assert res == 0

        

    def manual_test_that_output_is_as_expected(self):
        from run_okushiri import main
        from create_okushiri import create_mesh

        # Note, this has not been tested when code is failing.
        # Next step, use the parrameter elevation_in_mesh
        # and only run create_okishiri once for speed.
        
        res = create_mesh()
        assert res == None
        
        res = main()
        assert res == None
        
        s = 'compare_timeseries_with_measures.py'
        #print s
        res = os.system('python %s > compare_timeseries_with_measures.stdout'\
                        %s)
        assert res == 0
        

#-------------------------------------------------------------
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_Okushiri, 'test')
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite)

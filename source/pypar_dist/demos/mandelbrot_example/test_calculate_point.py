"""
  Test of calculate_point within mandelbrot.py
"""

import unittest
import Numeric
from mandelbrot import calculate_point
#from mandel_ext import calculate_point


#-------------------------------------------------------------

class TestCase(unittest.TestCase):

    def setUp(self):
	pass

	
    def test_basic(self):
        """Test that c==0 always returns kmax
	"""

	c = complex(0, 0)
	
	for kmax in [0,1,2,3,4,10, 20, 33, 45, 100]:
	    count = calculate_point(c, kmax)	
            assert count == kmax    	

	    
    def test_points_clearly_outside_set(self):	
    
    	for kmax in [1,2,3,5,10,20,100]:
	    for c in [complex(5, 4), complex(2,1), complex(-2,1),
	    	      complex(2,-1), complex(-2,-1)]:
	        count = calculate_point(c, kmax)
	        assert count == 1
	
		
    def test_arbitrary_points(self):
    	for kmax in [10,20,100]:	    
	    c = complex(1, 1.5)
	    count = calculate_point(c, kmax)
	    assert count == 2    
	    
	    c = complex(-1, 1.5)
	    count = calculate_point(c, kmax)
	    assert count == 2    	    
	    
	    c = complex(-1, -1.5)
	    count = calculate_point(c, kmax)
	    assert count == 2    	    	    
	    
	    c = complex(1, -1.5)
	    count = calculate_point(c, kmax)
	    assert count == 2    	    	    	    
	    
	    c = complex(1, 0)
	    count = calculate_point(c, kmax)
	    assert count == 3    	    	    	    	    
	    
	    c = complex(0.2, 1)
	    count = calculate_point(c, kmax)
	    assert count == 4    	    	    	    	    	    
	    
	    c = complex(0.2, 0.7)
	    count = calculate_point(c, kmax)
	    assert count == 6    	    	    	    	    	    	    
	    
	    c = complex(-0.2, 0.9)
	    count = calculate_point(c, kmax)
	    assert count == 9
            

    def test_another_known_point(self):
        c = complex(0.5, 0.5)
        count = calculate_point(c, 256)
        assert count == 5        

#-------------------------------------------------------------
if __name__ == "__main__":

    mysuite = unittest.makeSuite(TestCase,'test')
    runner = unittest.TextTestRunner()
    runner.run(mysuite)

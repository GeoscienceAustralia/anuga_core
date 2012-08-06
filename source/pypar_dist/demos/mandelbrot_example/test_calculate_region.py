"""
  Test of calculate_region within mandelbrot.py
"""

import unittest
import numpy
from mandelbrot import calculate_region


#-------------------------------------------------------------

class TestCase(unittest.TestCase):

    def setUp(self):
	pass

	
    def test1(self):
	kmax = 16
	M = N = 5
	
	real_min = -2.0
	real_max =  1.0
	imag_min = -1.5
	imag_max =  1.5
        A = calculate_region(real_min, real_max, 
		imag_min, imag_max, kmax, M, N)

	assert numpy.allclose(A,  		
	[[ 1,  1,  1,  1,  1],
 	[ 1,  3,  5,  5,  3],
 	[ 2,  4,  9,  9,  4],
 	[ 2,  9, 16, 16,  9],
 	[ 2,  3, 15, 15,  3]])


    def test2(self):
	kmax = 32
	M = N = 5
	
	real_min = -2.0
	real_max =  1.0
	imag_min = -1.5
	imag_max =  1.5
        A = calculate_region(real_min, real_max, 
		imag_min, imag_max, kmax, M, N)

	assert numpy.allclose(A,  		
	[[ 1,  1,  1,  1,  1],
 	[ 1,  3,  5,  5,  3],
 	[ 2,  4,  9,  9,  4],
 	[ 2,  9, 32, 32,  9],
 	[ 2,  3, 15, 15,  3]])
	

    def test3(self):
	kmax = 32
	M = N = 100
	
	real_min = -2.0
	real_max =  1.0
	imag_min = -1.5
	imag_max =  1.5
        A = calculate_region(real_min, real_max, 
		imag_min, imag_max, kmax, M, N)

	assert numpy.allclose(A[:, 50], 
	[32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 
	 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32,
         32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 
         32, 32, 32, 32, 32, 32, 32, 32,
         32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 
         32, 32, 32, 32, 32, 32, 32, 32,
         32, 32, 16, 11,  9,  7,  7,  6,  5,  5,  5,  4,  4,  4,  4,  4,  
         3, 3,  3,  3,  3,  3,  3,  3, 3,  3])

	
#-------------------------------------------------------------
if __name__ == "__main__":

    mysuite = unittest.makeSuite(TestCase,'test')
    runner = unittest.TextTestRunner()
    runner.run(mysuite)

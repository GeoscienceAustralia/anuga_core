#!/usr/bin/env python


import unittest
from Numeric import zeros, array, allclose
from Numeric import ArrayType, Float, Int, array, alltrue

from math import sqrt, pi
from anuga.config import epsilon
from numerical_tools import *

def test_function(x, y):
    return x+y

class Test_Numerical_Tools(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_angle1(self):
        """Test angles between one vector and the x-axis
	"""
        assert allclose(angle([1.0, 0.0])/pi*180, 0.0)	    
        assert allclose(angle([1.0, 1.0])/pi*180, 45.0)
        assert allclose(angle([0.0, 1.0])/pi*180, 90.0)		
        assert allclose(angle([-1.0, 1.0])/pi*180, 135.0)		
        assert allclose(angle([-1.0, 0.0])/pi*180, 180.0)
        assert allclose(angle([-1.0, -1.0])/pi*180, 225.0)
        assert allclose(angle([0.0, -1.0])/pi*180, 270.0)
        assert allclose(angle([1.0, -1.0])/pi*180, 315.0)
		
							  
    def test_angle2(self):
        """Test angles between two arbitrary vectors
	"""    
	
        assert allclose(angle([1.0, 0.0], [1.0, 1.0])/pi*180, 315.0)
        assert allclose(angle([1.0, 1.0], [1.0, 0.0])/pi*180, 45.0)
		
        assert allclose(angle([-1.0, -1.0], [1.0, 1.0])/pi*180, 180)	
        assert allclose(angle([-1.0, -1.0], [-1.0, 1.0])/pi*180, 90.0)	
	
        assert allclose(angle([-1.0, 0.0], [1.0, 1.0])/pi*180, 135.0)
        assert allclose(angle([0.0, -1.0], [1.0, 1.0])/pi*180, 225.0)	
	
        assert allclose(angle([1.0, -1.0], [1.0, 1.0])/pi*180, 270.0)	
        assert allclose(angle([1.0, 0.0], [0.0, 1.0])/pi*180, 270.0)

        #From test_get_boundary_polygon_V
        v_prev = [-0.5, -0.5]
        vc = [ 0.0,  -0.5]
        assert allclose(angle(vc, v_prev)/pi*180, 45.0)

        vc = [ 0.5,  0.0]
        assert allclose(angle(vc, v_prev)/pi*180, 135.0)

        vc = [ -0.5,  0.5]
        assert allclose(angle(vc, v_prev)/pi*180, 270.0)                


    def test_anglediff(self):
        assert allclose(anglediff([0.0, 1.], [1.0, 1.0])/pi*180, 45.0)

	
    def test_ensure_numeric(self):
        A = [1,2,3,4]
        B = ensure_numeric(A)
        assert type(B) == ArrayType
        assert B.typecode() == 'l'
        assert B[0] == 1 and B[1] == 2 and B[2] == 3 and B[3] == 4

        A = [1,2,3.14,4]
        B = ensure_numeric(A)
        assert type(B) == ArrayType
        assert B.typecode() == 'd'
        assert B[0] == 1 and B[1] == 2 and B[2] == 3.14 and B[3] == 4

        A = [1,2,3,4]
        B = ensure_numeric(A, Float)
        assert type(B) == ArrayType
        assert B.typecode() == 'd'
        assert B[0] == 1.0 and B[1] == 2.0 and B[2] == 3.0 and B[3] == 4.0

        A = [1,2,3,4]
        B = ensure_numeric(A, Float)
        assert type(B) == ArrayType
        assert B.typecode() == 'd'
        assert B[0] == 1.0 and B[1] == 2.0 and B[2] == 3.0 and B[3] == 4.0

        A = array([1,2,3,4])
        B = ensure_numeric(A)
        assert type(B) == ArrayType
        assert B.typecode() == 'l'        
        assert alltrue(A == B)    
        assert A is B   #Same object

        A = array([1,2,3,4])
        B = ensure_numeric(A, Float)
        assert type(B) == ArrayType
        assert B.typecode() == 'd'        
        assert alltrue(A == B)    
        assert A is not B   #Not the same object

        # Check scalars
        A = 1
        B = ensure_numeric(A, Float)
        #print A, B[0], len(B), type(B) 
        #print B.shape
        assert alltrue(A == B)

        B = ensure_numeric(A, Int)        
        #print A, B
        #print B.shape
        assert alltrue(A == B)

        # Error situation

        B = ensure_numeric('hello', Int)                
        assert allclose(B, [104, 101, 108, 108, 111])


    def test_gradient(self):
        x0 = 0.0; y0 = 0.0; z0 = 0.0
        x1 = 1.0; y1 = 0.0; z1 = -1.0
        x2 = 0.0; y2 = 1.0; z2 = 0.0

        zx, zy = gradient(x0, y0, x1, y1, x2, y2, z0, z1, z2)

        assert zx == -1.0
        assert zy == 0.0


    def test_gradient_more(self):
        x0 = 2.0/3; y0 = 2.0/3
        x1=  8.0/3; y1 = 2.0/3
        x2 = 2.0/3; y2 = 8.0/3

        q0 = 2.0+2.0/3
        q1 = 8.0+2.0/3
        q2 = 2.0+8.0/3

        #Gradient of fitted pwl surface
        a, b = gradient(x0, y0, x1, y1, x2, y2, q0, q1, q2)

        assert abs(a - 3.0) < epsilon
        assert abs(b - 1.0) < epsilon


    def test_gradient2(self):
        """Test two-point gradient
        """
        
        x0 = 5.0; y0 = 5.0; z0 = 10.0
        x1 = 8.0; y1 = 2.0; z1 = 1.0
        x2 = 8.0; y2 = 8.0; z2 = 10.0

        #Reference
        zx, zy = gradient(x0, y0, x1, y1, x2, y2, z0, z1, z2)
        a, b = gradient2(x0, y0, x1, y1, z0, z1)

        assert zx == a
        assert zy == b

        z2_computed = z0 + a*(x2-x0) + b*(y2-y0)
        assert z2_computed == z2

        
    def test_gradient2_more(self):
        """Test two-point gradient more
        """
        x0 = 2.0; y0 = 2.0
        x1 = 8.0; y1 = 3.0
        x2 = 1.0; y2 = 8.0

        q0 = 2.0
        q1 = 8.0
        q2 = q0

        #Gradient of fitted pwl surface
        a_ref, b_ref = gradient(x0, y0, x1, y1, x2, y2, q0, q1, q2)
        a, b = gradient2(x0, y0, x1, y1, q0, q1)        

        assert a == a_ref
        assert b == b_ref


    def test_machine_precision(self):
        """test_machine_precision(self):
        Test the function that calculates epsilon. As this varies on
        different machines, this is only an indication.
        """

        eps = get_machine_precision()

        assert eps < 1.0e-12, 'Machine precision should be better than 1.0e-12'
        assert eps > 0.0
        assert 1.0+eps/2 == 1.0
        
        
    def test_histogram(self):
        """Test histogram with different bin boundaries
        """
        
        a = [1,1,1,1,1,2,1,3,2,3,1,2,3,4,1]

        #There are four elements greater than or equal to 3
        bins = [3]
        assert allclose(histogram(a, bins), [4])

        bins = [ min(a) ]
        assert allclose(histogram(a, bins), [len(a)])

        bins = [ max(a)+0.00001 ]
        assert allclose(histogram(a, bins), [0])        
        
        bins = [1,2,3,4]
        assert allclose(histogram(a, bins), [8,3,3,1])

        bins = [1.1,2,3.1,4]
        #print histogram(a, bins)
        assert allclose(histogram(a, bins), [0,6,0,1])

        bins = [0,1.5,2,3]
        assert allclose(histogram(a, bins), [8,0,3,4])
        assert allclose(histogram(a, [0,3]), histogram(a, [-0.5,3]))

        # Check situation with #bins >= #datapoints
        a = [1.7]
        bins = [0,1.5,2,3]
        assert allclose(histogram(a, bins), [0,1,0,0])

        a = [1.7]
        bins = [0]
        assert allclose(histogram(a, bins), [1])

        a = [-1.7]
        bins = [0]
        assert allclose(histogram(a, bins), [0])

        a = [-1.7]
        bins = [-1.7]
        assert allclose(histogram(a, bins), [1])        
        

    def test_that_C_extension_compiles(self):
        FN = 'util_ext.c'
        try:
            import util_ext
        except:
            from compile import compile

            try:
                compile(FN)
            except:
                raise 'Could not compile %s' %FN
            else:
                import util_ext


    def test_gradient_C_extension(self):
        from util_ext import gradient as gradient_c

        x0 = 2.0/3; y0 = 2.0/3
        x1=  8.0/3; y1 = 2.0/3
        x2 = 2.0/3; y2 = 8.0/3

        q0 = 2.0+2.0/3
        q1 = 8.0+2.0/3
        q2 = 2.0+8.0/3

        #Gradient of fitted pwl surface
        a, b = gradient_c(x0, y0, x1, y1, x2, y2, q0, q1, q2)

        assert abs(a - 3.0) < epsilon
        assert abs(b - 1.0) < epsilon


    def test_gradient_C_extension3(self):
        from util_ext import gradient as gradient_c

        from RandomArray import uniform, seed
        seed(17, 53)

        x0, x1, x2, y0, y1, y2 = uniform(0.0,3.0,6)

        q0 = uniform(0.0, 10.0, 4)
        q1 = uniform(1.0, 3.0, 4)
        q2 = uniform(7.0, 20.0, 4)

        for i in range(4):
            #Gradient of fitted pwl surface
            a_ref, b_ref = gradient_python(x0, y0, x1, y1, x2, y2,
                                           q0[i], q1[i], q2[i])

            #print a_ref, b_ref
            a, b = gradient_c(x0, y0, x1, y1, x2, y2,
                              q0[i], q1[i], q2[i])

            #print a, a_ref, b, b_ref
            assert abs(a - a_ref) < epsilon
            assert abs(b - b_ref) < epsilon

     
    def test_err(self):
        x = [2,5] # diff at first position = 4, 4^2 = 16
        y = [6,7] # diff at secnd position = 2, 2^2 = 4
        # 16 + 4 = 20
        
        # If there is x and y, n=2 and relative=False, this will calc;
        # sqrt(sum_over_x&y((xi - yi)^2))
        err__1 = err(x,y,2,False)
        assert err__1 == sqrt(20)
        #print "err_", err_
        #rmsd_1 = err__1*sqrt(1./len(x))
        #print "err__1*sqrt(1./len(x))", err__1*sqrt(1./len(x))
        #print "sqrt(10)", sqrt(10)
        
        x = [2,7,100]
        y = [5,10,103]
        err__2 = err(x,y,2,False)
        assert err__2 == sqrt(27)
        #rmsd_2 = err__2*sqrt(1./len(x))
        #print "err__2*sqrt(1./len(x))", err__2*sqrt(1./len(x))

        x = [2,5,2,7,100]
        y = [6,7,5,10,103]
        err_3 = err(x,y,2,False)
        assert err_3 == sqrt(47)
        
        #rmsd_3 = err_3*sqrt(1./len(x))
        #print "err__3*sqrt(1./len(x))", err__3*sqrt(1./len(x))
        #print "rmsd_3", rmsd_3
        #print "sqrt(err_1*err__1+err__2*err__2)/sqrt(5)", \
        # sqrt(err__1*err__1+err__2*err__2)/sqrt(5)
        #print "(rmsd_1 + rmsd_2)/2.", (rmsd_1 + rmsd_2)/2.
        #print "sqrt((rmsd_1*rmsd_1 + rmsd_2*rmsd_2))/2.", \
        #sqrt((rmsd_1*rmsd_1 + rmsd_2*rmsd_2))/2.
        
    def test_norm(self):
        x = norm(ensure_numeric([3,4]))
        assert x == 5.

                                    

#-------------------------------------------------------------
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_Numerical_Tools,'test')
    #suite = unittest.makeSuite(Test_Numerical_Tools,'test_err')
    runner = unittest.TextTestRunner()
    runner.run(suite)

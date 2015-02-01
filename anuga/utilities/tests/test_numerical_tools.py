#!/usr/bin/env python


import unittest
import numpy as num
from numpy.random import uniform, seed

from math import sqrt, pi
from anuga.config import epsilon
from anuga.utilities.numerical_tools import *


class Test_Numerical_Tools(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_angle1(self):
        """Test angles between one vector and the x-axis
	"""
        assert num.allclose(angle([1.0, 0.0])/pi*180, 0.0)	    
        assert num.allclose(angle([1.0, 1.0])/pi*180, 45.0)
        assert num.allclose(angle([0.0, 1.0])/pi*180, 90.0)		
        assert num.allclose(angle([-1.0, 1.0])/pi*180, 135.0)		
        assert num.allclose(angle([-1.0, 0.0])/pi*180, 180.0)
        assert num.allclose(angle([-1.0, -1.0])/pi*180, 225.0)
        assert num.allclose(angle([0.0, -1.0])/pi*180, 270.0)
        assert num.allclose(angle([1.0, -1.0])/pi*180, 315.0)
		
							  
    def test_angle2(self):
        """Test angles between two arbitrary vectors
	"""    
	
        assert num.allclose(angle([1.0, 0.0], [1.0, 1.0])/pi*180, 315.0)
        assert num.allclose(angle([1.0, 1.0], [1.0, 0.0])/pi*180, 45.0)
		
        assert num.allclose(angle([-1.0, -1.0], [1.0, 1.0])/pi*180, 180)	
        assert num.allclose(angle([-1.0, -1.0], [-1.0, 1.0])/pi*180, 90.0)	
	
        assert num.allclose(angle([-1.0, 0.0], [1.0, 1.0])/pi*180, 135.0)
        assert num.allclose(angle([0.0, -1.0], [1.0, 1.0])/pi*180, 225.0)	
	
        assert num.allclose(angle([1.0, -1.0], [1.0, 1.0])/pi*180, 270.0)	
        assert num.allclose(angle([1.0, 0.0], [0.0, 1.0])/pi*180, 270.0)

        #From test_get_boundary_polygon_V
        v_prev = [-0.5, -0.5]
        vc = [ 0.0,  -0.5]
        assert num.allclose(angle(vc, v_prev)/pi*180, 45.0)

        vc = [ 0.5,  0.0]
        assert num.allclose(angle(vc, v_prev)/pi*180, 135.0)

        vc = [ -0.5,  0.5]
        assert num.allclose(angle(vc, v_prev)/pi*180, 270.0)                


    def test_anglediff(self):
        assert num.allclose(anglediff([0.0, 1.], [1.0, 1.0])/pi*180, 45.0)

	
    def test_ensure_numeric(self):
        A = [1,2,3,4]
        B = ensure_numeric(A)
        assert isinstance(B, num.ndarray)
        assert B.dtype.char == 'l'
        assert B[0] == 1 and B[1] == 2 and B[2] == 3 and B[3] == 4

        A = [1,2,3.14,4]
        B = ensure_numeric(A)
        assert isinstance(B, num.ndarray)
        assert B.dtype.char == 'd'
        assert B[0] == 1 and B[1] == 2 and B[2] == 3.14 and B[3] == 4

        A = [1,2,3,4]
        B = ensure_numeric(A, num.float)
        assert isinstance(B, num.ndarray)
        assert B.dtype.char == 'd'
        assert B[0] == 1.0 and B[1] == 2.0 and B[2] == 3.0 and B[3] == 4.0

        A = [1,2,3,4]
        B = ensure_numeric(A, num.float)
        assert isinstance(B, num.ndarray)
        assert B.dtype.char == 'd'
        assert B[0] == 1.0 and B[1] == 2.0 and B[2] == 3.0 and B[3] == 4.0

        A = num.array([1,2,3,4])
        B = ensure_numeric(A)
        assert isinstance(B, num.ndarray)
        assert B.dtype.char == 'l'        
        assert num.alltrue(A == B)    
        assert A is B   #Same object

        # check default num.array type, which is supposed to be num.int32
        A = num.array((1,2,3,4))
        assert isinstance(A, num.ndarray)
        msg = "Expected dtype.char='l', got '%s'" % A.dtype.char
        assert A.dtype.char == 'l', msg

        A = num.array([1,2,3,4])
        B = ensure_numeric(A, num.float)
        assert isinstance(B, num.ndarray)
        assert A.dtype.char == 'l'        
        assert B.dtype.char == 'd'        
        assert num.alltrue(A == B)    
        assert A is not B   # Not the same object

        # Check scalars
        A = 1
        B = ensure_numeric(A, num.float)
        assert num.alltrue(A == B)

        B = ensure_numeric(A, num.int)        
        assert num.alltrue(A == B)

#        # try to simulate getting (x,0) shape
#        data_points = [[ 413634. ],]
#        array_data_points = ensure_numeric(data_points)
#        if not (0,) == array_data_points.shape:
#            assert len(array_data_points.shape) == 2
#            assert array_data_points.shape[1] == 2

        # strings input should raise exception
        self.assertRaises(Exception, ensure_numeric(['abc',]))
        self.assertRaises(Exception, ensure_numeric(('abc',)))
        self.assertRaises(Exception, ensure_numeric(num.array(('abc',))))

    def NO_ensure_numeric_char(self):
        '''numpy can't handle this'''

        # Error situation
        B = ensure_numeric('hello', num.int)                
        assert num.allclose(B, [104, 101, 108, 108, 111])


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
        assert num.allclose(histogram(a, bins), [4])

        bins = [ min(a) ]
        assert num.allclose(histogram(a, bins), [len(a)])

        bins = [ max(a)+0.00001 ]
        assert num.allclose(histogram(a, bins), [0])        
        
        bins = [1,2,3,4]
        assert num.allclose(histogram(a, bins), [8,3,3,1])

        bins = [1.1,2,3.1,4]
        #print histogram(a, bins)
        assert num.allclose(histogram(a, bins), [0,6,0,1])

        bins = [0,1.5,2,3]
        assert num.allclose(histogram(a, bins), [8,0,3,4])
        assert num.allclose(histogram(a, [0,3]), histogram(a, [-0.5,3]))

        # Check situation with #bins >= #datapoints
        a = [1.7]
        bins = [0,1.5,2,3]
        assert num.allclose(histogram(a, bins), [0,1,0,0])

        a = [1.7]
        bins = [0]
        assert num.allclose(histogram(a, bins), [1])

        a = [-1.7]
        bins = [0]
        assert num.allclose(histogram(a, bins), [0])

        a = [-1.7]
        bins = [-1.7]
        assert num.allclose(histogram(a, bins), [1])        
        

    def test_that_C_extension_compiles(self):
        FN = 'util_ext.c'
        try:
            import anuga.utilities.util_ext as util_ext
        except:
            from anuga.utilities.compile import compile

            try:
                compile(FN)
            except:
                raise Exception('Could not compile %s' % FN)
            else:
                import anuga.utilities.util_ext as util_ext


    def test_gradient_C_extension(self):
        from anuga.utilities.util_ext import gradient as gradient_c

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
        from anuga.utilities.util_ext import gradient as gradient_c

        seed((17, 53))

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


################################################################################
# Test the is_num_????() functions.
################################################################################

    def test_is_float(self):
        def t(val, expected):
            if expected == True:
                msg = 'should be num.float?'
            else:
                msg = 'should not be num.float?'
            msg = '%s (%s) %s' % (str(val), type(val), msg)
            assert is_num_float(val) == expected, msg

        t(1, False)
        t(1.0, False)
        t('abc', False)
        t(None, False)
        t(num.array(None), False)
        # can't create array(None, num.int)
#        t(num.array(None, num.int), False)
        t(num.array(None, num.float), True)
        t(num.array(()), True)
        t(num.array((), num.int), False)
        t(num.array((), num.float), True)
        t(num.array((1), num.int), False)
        t(num.array((1), num.float), True)

        t(num.array((1,2)), False)
        t(num.array((1,2), num.int), False)
        t(num.array((1,2), num.float), True)
        t(num.array([1,2]), False)
        t(num.array([1,2], num.int), False)
        t(num.array([1,2], num.float), True)

        t(num.array((1.0,2.0)), True)
        t(num.array((1.0,2.0), num.int), False)
        t(num.array((1.0,2.0), num.float), True)
        t(num.array([1.0,2.0]), True)
        t(num.array([1.0,2.0], num.int), False)
        t(num.array([1.0,2.0], num.float), True)

        t(num.array(((1.0,2.0),(3.0,4.0))), True)
        t(num.array(((1.0,2.0),(3.0,4.0)), num.int), False)
        t(num.array(((1.0,2.0),(3.0,4.0)), num.float), True)
        t(num.array([[1.0,2.0],[3.0,4.0]]), True)
        t(num.array([1.0,2.0], num.int), False)
        t(num.array([1.0,2.0], num.float), True)

        t(num.array('abc'), False)
        t(num.array('abc', num.character), False)
        # can't create array as int from string
#        t(num.array('abc', num.int), False)
        # can't create array as float from string
#        t(num.array('abc', num.float), True)

    def test_is_int(self):
        def t(val, expected):
            if expected == True:
                msg = 'should be num.int?'
            else:
                msg = 'should not be num.int?'
            msg = '%s (%s) %s' % (str(val), type(val), msg)
            assert is_num_int(val) == expected, msg

        t(1, False)
        t(1.0, False)
        t('abc', False)
        t(None, False)
        t(num.array(None), False)
        # can't create array(None, num.int)
#        t(num.array(None, num.int), True)
        t(num.array(None, num.float), False)
        t(num.array((), num.int), True)
        t(num.array(()), False)
        t(num.array((), num.float), False)
        t(num.array((1), num.int), True)
        t(num.array((1), num.float), False)

        t(num.array((1,2)), True)
        t(num.array((1,2), num.int), True)
        t(num.array((1,2), num.float), False)
        t(num.array([1,2]), True)
        t(num.array([1,2], num.int), True)
        t(num.array([1,2], num.float), False)

        t(num.array((1.0,2.0)), False)
        t(num.array((1.0,2.0), num.int), True)
        t(num.array((1.0,2.0), num.float), False)
        t(num.array([1.0,2.0]), False)
        t(num.array([1.0,2.0], num.int), True)
        t(num.array([1.0,2.0], num.float), False)

        t(num.array(((1.0,2.0),(3.0,4.0))), False)
        t(num.array(((1.0,2.0),(3.0,4.0)), num.int), True)
        t(num.array(((1.0,2.0),(3.0,4.0)), num.float), False)
        t(num.array([[1.0,2.0],[3.0,4.0]]), False)
        t(num.array([1.0,2.0], num.int), True)
        t(num.array([1.0,2.0], num.float), False)

        t(num.array('abc'), False)
        t(num.array('abc', num.character), False)
        # can't create array as int from string
#        t(num.array('abc', num.int), True)
        # can't create array as float from string
#        t(num.array('abc', num.float), False)

    def test_ensure_numeric_copy(self):
        """Test to see if ensure_numeric() behaves as we expect.

        Under Numeric ensure_numeric() *always* returned a copy (bug).
        Under numpy it copies only when it has to.
        """

        #####
        # Make 'points' a _list_ of coordinates.
        # Should be changed by ensure_numeric().
        #####
        points = [[1.,2.], [3.,4.], [5.,6.]]
        points_id = id(points)

        points_new = ensure_numeric(points, num.float)
        points_new_id = id(points_new)

        msg = 'ensure_numeric() should return a copy of a list'
        self.assertTrue(points_new_id != points_id, msg)

        # should never change it's input parameter
        msg = "ensure_numeric() changed it's input parameter"
        self.assertTrue(points_id == id(points), msg)

        #####
        # Make 'points' a _tuple_ of coordinates.
        # Should be changed by ensure_numeric().
        #####
        points = ((1.,2.), (3.,4.), (5.,6.))
        points_id = id(points)

        points_new = ensure_numeric(points, num.int)
        points_new_id = id(points_new)

        msg = 'ensure_numeric() should return a copy of a list'
        self.assertTrue(points_new_id != points_id, msg)

        # should never change it's input parameter
        msg = "ensure_numeric() changed it's input parameter"
        self.assertTrue(points_id == id(points), msg)

        #####
        # Make 'points' a numeric array of float coordinates.
        # Should NOT be changed by ensure_numeric().
        #####
        points = num.array([[1.,2.], [3.,4.], [5.,6.]], num.float)
        points_id = id(points)

        points_new = ensure_numeric(points, num.float)
        points_new_id = id(points_new)

        msg = 'ensure_numeric() should return the original input'
        self.assertTrue(points_new_id == points_id, msg)

        # should never change it's input parameter
        msg = "ensure_numeric() changed it's input parameter"
        self.assertTrue(points_id == id(points), msg)

        #####
        # Make 'points' a numeric array of int coordinates.
        # Should be changed by ensure_numeric(, num.float).
        #####
        points = num.array([[1,2], [3,4], [5,6]], num.int)
        points_id = id(points)

        points_new = ensure_numeric(points, num.float)
        points_new_id = id(points_new)

        msg = 'ensure_numeric() should return a copy of the input'
        self.assertTrue(points_new_id != points_id, msg)

        # should never change it's input parameter
        msg = "ensure_numeric() changed it's input parameter"
        self.assertTrue(points_id == id(points), msg)

        #####
        # Make 'points' a numeric array of int coordinates.
        # Should NOT be changed by ensure_numeric(, num.int).
        #####
        points = num.array([[1,2], [3,4], [5,6]], num.int)
        points_id = id(points)

        points_new = ensure_numeric(points, num.int)
        points_new_id = id(points_new)

        msg = 'ensure_numeric() should return the original input'
        self.assertTrue(points_new_id == points_id, msg)

        # should never change it's input parameter
        msg = "ensure_numeric() changed it's input parameter"
        self.assertTrue(points_id == id(points), msg)

        #####
        # Make 'points' a numeric array of float32 coordinates.
        # Should NOT be changed by ensure_numeric(, num.float32).
        #####
        points = num.array([[1.,2.], [3.,4.], [5.,6.]], num.float32)
        points_id = id(points)

        points_new = ensure_numeric(points, num.float32)
        points_new_id = id(points_new)

        msg = 'ensure_numeric() should return the original input'
        self.assertTrue(points_new_id == points_id, msg)

        # should never change it's input parameter
        msg = "ensure_numeric() changed it's input parameter"
        self.assertTrue(points_id == id(points), msg)

        #####
        # Make 'points' a numeric array of float32 coordinates.
        # Should be changed by ensure_numeric(, num.float64).
        #####
        points = num.array([[1.,2.], [3.,4.], [5.,6.]], num.float32)
        points_id = id(points)

        points_new = ensure_numeric(points, num.float64)
        points_new_id = id(points_new)

        msg = 'ensure_numeric() should return a copy of the input'
        self.assertTrue(points_new_id != points_id, msg)

        # should never change it's input parameter
        msg = "ensure_numeric() changed it's input parameter"
        self.assertTrue(points_id == id(points), msg)

        #####
        # Make 'points' a numeric array of float coordinates.
        # Should NOT be changed by ensure_numeric(, num.float64).
        #####
        points = num.array([[1.,2.], [3.,4.], [5.,6.]], num.float)
        points_id = id(points)

        points_new = ensure_numeric(points, num.float64)
        points_new_id = id(points_new)

        msg = 'ensure_numeric() should return the original input'
        self.assertTrue(points_new_id == points_id, msg)
        #msg = 'ensure_numeric() should return a copy of the input'
        #self.assertTrue(points_new_id != points_id, msg)

        # should never change it's input parameter
        msg = "ensure_numeric() changed it's input parameter"
        self.assertTrue(points_id == id(points), msg)

################################################################################

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_Numerical_Tools,'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)

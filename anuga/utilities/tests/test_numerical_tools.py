#!/usr/bin/env python


import unittest
import numpy as num
import numpy as np

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
        assert num.allclose((angle([1.0, 0.0]) / pi)*180, 0.0)
        assert num.allclose((angle([1.0, 1.0]) / pi)*180, 45.0)
        assert num.allclose((angle([0.0, 1.0]) / pi)*180, 90.0)
        assert num.allclose((angle([-1.0, 1.0]) / pi)*180, 135.0)
        assert num.allclose((angle([-1.0, 0.0]) / pi)*180, 180.0)
        assert num.allclose((angle([-1.0, -1.0]) / pi)*180, 225.0)
        assert num.allclose((angle([0.0, -1.0]) / pi)*180, 270.0)
        assert num.allclose((angle([1.0, -1.0]) / pi)*180, 315.0)

    def test_angle2(self):
        """Test angles between two arbitrary vectors
        """

        assert num.allclose(
            (angle([1.0, 0.0], [1.0, 1.0]) / pi)*180, 315.0)
        assert num.allclose(
            (angle([1.0, 1.0], [1.0, 0.0]) / pi)*180, 45.0)

        assert num.allclose(
            (angle([-1.0, -1.0], [1.0, 1.0]) / pi)*180, 180)
        assert num.allclose(
            (angle([-1.0, -1.0], [-1.0, 1.0]) / pi)*180, 90.0)

        assert num.allclose(
            (angle([-1.0, 0.0], [1.0, 1.0]) / pi)*180, 135.0)
        assert num.allclose(
            (angle([0.0, -1.0], [1.0, 1.0]) / pi)*180, 225.0)

        assert num.allclose(
            (angle([1.0, -1.0], [1.0, 1.0]) / pi)*180, 270.0)
        assert num.allclose(
            (angle([1.0, 0.0], [0.0, 1.0]) / pi)*180, 270.0)

        # From test_get_boundary_polygon_V
        v_prev = [-0.5, -0.5]
        vc = [0.0,  -0.5]
        assert num.allclose((angle(vc, v_prev) / pi)*180, 45.0)

        vc = [0.5,  0.0]
        assert num.allclose((angle(vc, v_prev) / pi)*180, 135.0)

        vc = [-0.5,  0.5]
        assert num.allclose((angle(vc, v_prev) / pi)*180, 270.0)

    def test_anglediff(self):
        assert num.allclose(
            (anglediff([0.0, 1.], [1.0, 1.0]) / pi)*180, 45.0)

    def test_ensure_numeric(self):
        A = [1, 2, 3, 4]
        B = ensure_numeric(A)
        assert isinstance(B, num.ndarray)
        # check that dtype is long or long long
        assert B.dtype.char == 'l' or B.dtype.char == 'q'
        assert B.dtype.itemsize == 8
        assert B[0] == 1 and B[1] == 2 and B[2] == 3 and B[3] == 4

        A = [1, 2, 3.14, 4]
        B = ensure_numeric(A)
        assert isinstance(B, num.ndarray)
        assert B.dtype.char == 'd'
        assert B.dtype.itemsize == 8
        assert B[0] == 1 and B[1] == 2 and B[2] == 3.14 and B[3] == 4

        A = [1, 2, 3, 4]
        B = ensure_numeric(A, float)
        assert isinstance(B, num.ndarray)
        assert B.dtype.char == 'd'
        assert B.dtype.itemsize == 8
        assert B[0] == 1.0 and B[1] == 2.0 and B[2] == 3.0 and B[3] == 4.0

        A = [1, 2, 3, 4]
        B = ensure_numeric(A, float)
        assert isinstance(B, num.ndarray)
        assert B.dtype.char == 'd'
        assert B.dtype.itemsize == 8
        assert B[0] == 1.0 and B[1] == 2.0 and B[2] == 3.0 and B[3] == 4.0

        A = num.array([1, 2, 3, 4])
        B = ensure_numeric(A)
        assert isinstance(B, num.ndarray)
        assert B.dtype.char == 'l' or B.dtype.char == 'q'
        assert num.all(A == B)
        assert A is B  # Same object

        # check default num.array type, which is supposed to be num.int64
        A = num.array((1, 2, 3, 4))
        assert isinstance(A, num.ndarray)
        msg = "Expected dtype.char='l' or 'q', got '%s'" % A.dtype.char
        assert A.dtype.char == 'l' or A.dtype.char == 'q', msg

        A = num.array([1, 2, 3, 4])
        B = ensure_numeric(A, float)
        assert isinstance(B, num.ndarray)
        assert A.dtype.char == 'l' or A.dtype.char == 'q'
        assert B.dtype.char == 'd'
        assert num.all(A == B)
        assert A is not B   # Not the same object

        # Check scalars
        A = 1
        B = ensure_numeric(A, float)
        assert num.all(A == B)

        B = ensure_numeric(A, int)
        assert num.all(A == B)

#        # try to simulate getting (x,0) shape
#        data_points = [[ 413634. ],]
#        array_data_points = ensure_numeric(data_points)
#        if not (0,) == array_data_points.shape:
#            assert len(array_data_points.shape) == 2
#            assert array_data_points.shape[1] == 2

        # strings input should raise exception
        self.assertRaises(Exception, ensure_numeric(['abc', ]))
        self.assertRaises(Exception, ensure_numeric(('abc',)))
        self.assertRaises(Exception, ensure_numeric(num.array(('abc',))))

    def NO_ensure_numeric_char(self):
        '''numpy can't handle this'''

        # Error situation
        B = ensure_numeric('hello', int)
        assert num.allclose(B, [104, 101, 108, 108, 111])

    def test_gradient(self):
        x0 = 0.0
        y0 = 0.0
        z0 = 0.0
        x1 = 1.0
        y1 = 0.0
        z1 = -1.0
        x2 = 0.0
        y2 = 1.0
        z2 = 0.0

        zx, zy = gradient(x0, y0, x1, y1, x2, y2, z0, z1, z2)

        assert zx == -1.0
        assert zy == 0.0

    def test_gradient_more(self):
        x0 = 2.0/3
        y0 = 2.0/3
        x1 = 8.0/3
        y1 = 2.0/3
        x2 = 2.0/3
        y2 = 8.0/3

        q0 = 2.0+2.0/3
        q1 = 8.0+2.0/3
        q2 = 2.0+8.0/3

        # Gradient of fitted pwl surface
        a, b = gradient(x0, y0, x1, y1, x2, y2, q0, q1, q2)

        assert abs(a - 3.0) < epsilon
        assert abs(b - 1.0) < epsilon

    def test_gradient2(self):
        """Test two-point gradient
        """

        x0 = 5.0
        y0 = 5.0
        z0 = 10.0
        x1 = 8.0
        y1 = 2.0
        z1 = 1.0
        x2 = 8.0
        y2 = 8.0
        z2 = 10.0

        # Reference
        zx, zy = gradient(x0, y0, x1, y1, x2, y2, z0, z1, z2)
        a, b = gradient2(x0, y0, x1, y1, z0, z1)

        assert zx == a
        assert zy == b

        z2_computed = z0 + a*(x2-x0) + b*(y2-y0)
        assert z2_computed == z2

    def test_gradient2_more(self):
        """Test two-point gradient more
        """
        x0 = 2.0
        y0 = 2.0
        x1 = 8.0
        y1 = 3.0
        x2 = 1.0
        y2 = 8.0

        q0 = 2.0
        q1 = 8.0
        q2 = q0

        # Gradient of fitted pwl surface
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
        assert 1.0 + (eps / 2) == 1.0

    def test_histogram(self):
        """Test histogram with different bin boundaries
        """

        a = [1, 1, 1, 1, 1, 2, 1, 3, 2, 3, 1, 2, 3, 4, 1]

        # There are four elements greater than or equal to 3
        bins = [3]
        assert num.allclose(histogram(a, bins), [4])

        bins = [min(a)]
        assert num.allclose(histogram(a, bins), [len(a)])

        bins = [max(a)+0.00001]
        assert num.allclose(histogram(a, bins), [0])

        bins = [1, 2, 3, 4]
        assert num.allclose(histogram(a, bins), [8, 3, 3, 1])

        bins = [1.1, 2, 3.1, 4]
        # print histogram(a, bins)
        assert num.allclose(histogram(a, bins), [0, 6, 0, 1])

        bins = [0, 1.5, 2, 3]
        assert num.allclose(histogram(a, bins), [8, 0, 3, 4])
        assert num.allclose(histogram(a, [0, 3]), histogram(a, [-0.5, 3]))

        # Check situation with #bins >= #datapoints
        a = [1.7]
        bins = [0, 1.5, 2, 3]
        assert num.allclose(histogram(a, bins), [0, 1, 0, 0])

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
        except ImportError:
            from anuga.utilities.compile import compile

            try:
                compile(FN)
            except Exception:
                raise Exception('Could not compile %s' % FN)
            else:
                import anuga.utilities.util_ext as util_ext

    def test_gradient_C_extension(self):
        from anuga.utilities.util_ext import gradient as gradient_c

        x0 = 2.0/3
        y0 = 2.0/3
        x1 = 8.0/3
        y1 = 2.0/3
        x2 = 2.0/3
        y2 = 8.0/3

        q0 = 2.0+2.0/3
        q1 = 8.0+2.0/3
        q2 = 2.0+8.0/3

        # Gradient of fitted pwl surface
        a, b = gradient_c(x0, y0, x1, y1, x2, y2, q0, q1, q2)

        assert abs(a - 3.0) < epsilon
        assert abs(b - 1.0) < epsilon

    def test_gradient_C_extension_degenerate_points(self):
        """Collinear or coincident points have no plane to fit.

        The C helpers used to divide by a zero determinant and hand back
        NaN, which then propagated into the friction and gravity terms.
        They now fall back to a flat (first-order) gradient.
        """
        from anuga.utilities.util_ext import gradient as gradient_c
        from anuga.utilities.util_ext import gradient2 as gradient2_c

        # Three collinear points along the x axis
        a, b = gradient_c(0.0, 0.0, 1.0, 0.0, 2.0, 0.0, 1.0, 2.0, 3.0)
        assert a == 0.0 and b == 0.0

        # Three coincident points
        a, b = gradient_c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 3.0)
        assert a == 0.0 and b == 0.0

        # Two coincident points
        a, b = gradient2_c(1.0, 1.0, 1.0, 1.0, 0.0, 5.0)
        assert a == 0.0 and b == 0.0

        # A well-posed case is untouched by the guard
        a, b = gradient_c(0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 2.0)
        assert abs(a - 1.0) < epsilon and abs(b - 2.0) < epsilon
        a, b = gradient2_c(0.0, 0.0, 2.0, 0.0, 0.0, 4.0)
        assert abs(a - 2.0) < epsilon and abs(b - 0.0) < epsilon

    def test_gradient_C_extension3(self):
        from anuga.utilities.util_ext import gradient as gradient_c

        rng = np.random.default_rng([17, 53])

        x0, x1, x2, y0, y1, y2 = rng.uniform(0.0, 3.0, 6)

        q0 = rng.uniform(0.0, 10.0, 4)
        q1 = rng.uniform(1.0, 3.0, 4)
        q2 = rng.uniform(7.0, 20.0, 4)

        for i in range(4):
            # Gradient of fitted pwl surface
            a_ref, b_ref = gradient_python(x0, y0, x1, y1, x2, y2,
                                           q0[i], q1[i], q2[i])

            # print a_ref, b_ref
            a, b = gradient_c(x0, y0, x1, y1, x2, y2,
                              q0[i], q1[i], q2[i])

            # print a, a_ref, b, b_ref
            assert abs(a - a_ref) < epsilon
            assert abs(b - b_ref) < epsilon

    def test_err(self):
        x = [2, 5]  # diff at first position = 4, 4^2 = 16
        y = [6, 7]  # diff at secnd position = 2, 2^2 = 4
        # 16 + 4 = 20

        # If there is x and y, n=2 and relative=False, this will calc;
        # sqrt(sum_over_x&y((xi - yi)^2))
        err__1 = err(x, y, 2, False)
        assert err__1 == sqrt(20)
        # print "err_", err_
        #rmsd_1 = err__1*sqrt(1./len(x))
        # print "err__1*sqrt(1./len(x))", err__1*sqrt(1./len(x))
        # print "sqrt(10)", sqrt(10)

        x = [2, 7, 100]
        y = [5, 10, 103]
        err__2 = err(x, y, 2, False)
        assert err__2 == sqrt(27)
        #rmsd_2 = err__2*sqrt(1./len(x))
        # print "err__2*sqrt(1./len(x))", err__2*sqrt(1./len(x))

        x = [2, 5, 2, 7, 100]
        y = [6, 7, 5, 10, 103]
        err_3 = err(x, y, 2, False)
        assert err_3 == sqrt(47)

        #rmsd_3 = err_3*sqrt(1./len(x))
        # print "err__3*sqrt(1./len(x))", err__3*sqrt(1./len(x))
        # print "rmsd_3", rmsd_3
        # print "sqrt(err_1*err__1+err__2*err__2)/sqrt(5)", \
        # sqrt(err__1*err__1+err__2*err__2)/sqrt(5)
        # print "(rmsd_1 + rmsd_2)/2.", (rmsd_1 + rmsd_2)/2.
        # print "sqrt((rmsd_1*rmsd_1 + rmsd_2*rmsd_2))/2.", \
        #sqrt((rmsd_1*rmsd_1 + rmsd_2*rmsd_2))/2.

    def test_norm(self):
        x = norm(ensure_numeric([3, 4]))
        assert x == 5.


################################################################################
# Test the is_num_????() functions.
################################################################################

    def test_is_float(self):
        def t(val, expected):
            if expected == True:
                msg = 'should be float?'
            else:
                msg = 'should not be float?'
            msg = '%s (%s) %s' % (str(val), type(val), msg)
            assert is_num_float(val) == expected, msg

        t(1, False)
        t(1.0, False)
        t('abc', False)
        t(None, False)
        t(num.array(None), False)
        # can't create array(None, int)
#        t(num.array(None, int), False)
        t(num.array(None, float), True)
        t(num.array(()), True)
        t(num.array((), int), False)
        t(num.array((), float), True)
        t(num.array((1), int), False)
        t(num.array((1), float), True)

        t(num.array((1, 2)), False)
        t(num.array((1, 2), int), False)
        t(num.array((1, 2), float), True)
        t(num.array([1, 2]), False)
        t(num.array([1, 2], int), False)
        t(num.array([1, 2], float), True)

        t(num.array((1.0, 2.0)), True)
        t(num.array((1.0, 2.0), int), False)
        t(num.array((1.0, 2.0), float), True)
        t(num.array([1.0, 2.0]), True)
        t(num.array([1.0, 2.0], int), False)
        t(num.array([1.0, 2.0], float), True)

        t(num.array(((1.0, 2.0), (3.0, 4.0))), True)
        t(num.array(((1.0, 2.0), (3.0, 4.0)), int), False)
        t(num.array(((1.0, 2.0), (3.0, 4.0)), float), True)
        t(num.array([[1.0, 2.0], [3.0, 4.0]]), True)
        t(num.array([1.0, 2.0], int), False)
        t(num.array([1.0, 2.0], float), True)

        t(num.array('abc'), False)
        t(num.array('abc', 'S1'), False)
        # can't create array as int from string
#        t(num.array('abc', int), False)
        # can't create array as float from string
#        t(num.array('abc', float), True)

    def test_is_int(self):
        def t(val, expected):
            if expected == True:
                msg = 'should be int?'
            else:
                msg = 'should not be int?'
            msg = '%s (%s) %s' % (str(val), type(val), msg)
            assert is_num_int(val) == expected, msg

        t(1, False)
        t(1.0, False)
        t('abc', False)
        t(None, False)
        t(num.array(None), False)
        # can't create array(None, int)
#        t(num.array(None, int), True)
        t(num.array(None, float), False)
        t(num.array((), int), True)
        t(num.array(()), False)
        t(num.array((), float), False)
        t(num.array((1), int), True)
        t(num.array((1), float), False)

        t(num.array((1, 2)), True)
        t(num.array((1, 2), int), True)
        t(num.array((1, 2), float), False)
        t(num.array([1, 2]), True)
        t(num.array([1, 2], int), True)
        t(num.array([1, 2], float), False)

        t(num.array((1.0, 2.0)), False)
        t(num.array((1.0, 2.0), int), True)
        t(num.array((1.0, 2.0), float), False)
        t(num.array([1.0, 2.0]), False)
        t(num.array([1.0, 2.0], int), True)
        t(num.array([1.0, 2.0], float), False)

        t(num.array(((1.0, 2.0), (3.0, 4.0))), False)
        t(num.array(((1.0, 2.0), (3.0, 4.0)), int), True)
        t(num.array(((1.0, 2.0), (3.0, 4.0)), float), False)
        t(num.array([[1.0, 2.0], [3.0, 4.0]]), False)
        t(num.array([1.0, 2.0], int), True)
        t(num.array([1.0, 2.0], float), False)

        t(num.array('abc'), False)
        t(num.array('abc', 'S1'), False)
        # can't create array as int from string
#        t(num.array('abc', int), True)
        # can't create array as float from string
#        t(num.array('abc', float), False)

    def test_ensure_numeric_copy(self):
        """Test to see if ensure_numeric() behaves as we expect.

        Under Numeric ensure_numeric() *always* returned a copy (bug).
        Under numpy it copies only when it has to.
        """

        #####
        # Make 'points' a _list_ of coordinates.
        # Should be changed by ensure_numeric().
        #####
        points = [[1., 2.], [3., 4.], [5., 6.]]
        points_id = id(points)

        points_new = ensure_numeric(points, float)
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
        points = ((1., 2.), (3., 4.), (5., 6.))
        points_id = id(points)

        points_new = ensure_numeric(points, int)
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
        points = num.array([[1., 2.], [3., 4.], [5., 6.]], float)
        points_id = id(points)

        points_new = ensure_numeric(points, float)
        points_new_id = id(points_new)

        msg = 'ensure_numeric() should return the original input'
        self.assertTrue(points_new_id == points_id, msg)

        # should never change it's input parameter
        msg = "ensure_numeric() changed it's input parameter"
        self.assertTrue(points_id == id(points), msg)

        #####
        # Make 'points' a numeric array of int coordinates.
        # Should be changed by ensure_numeric(, float).
        #####
        points = num.array([[1, 2], [3, 4], [5, 6]], int)
        points_id = id(points)

        points_new = ensure_numeric(points, float)
        points_new_id = id(points_new)

        msg = 'ensure_numeric() should return a copy of the input'
        self.assertTrue(points_new_id != points_id, msg)

        # should never change it's input parameter
        msg = "ensure_numeric() changed it's input parameter"
        self.assertTrue(points_id == id(points), msg)

        #####
        # Make 'points' a numeric array of int coordinates.
        # Should NOT be changed by ensure_numeric(, int).
        #####
        points = num.array([[1, 2], [3, 4], [5, 6]], int)
        points_id = id(points)

        points_new = ensure_numeric(points, int)
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
        points = num.array([[1., 2.], [3., 4.], [5., 6.]], num.float32)
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
        points = num.array([[1., 2.], [3., 4.], [5., 6.]], num.float32)
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
        points = num.array([[1., 2.], [3., 4.], [5., 6.]], float)
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

class Test_norms(unittest.TestCase):
    """Tests for anuga.utilities.norms (l1_norm, l2_norm, linf_norm)."""

    def test_l1_norm(self):
        from anuga.utilities.norms import l1_norm
        self.assertAlmostEqual(l1_norm([1.0, -2.0, 3.0]), 6.0)

    def test_l2_norm(self):
        from anuga.utilities.norms import l2_norm
        self.assertAlmostEqual(l2_norm([3.0, 4.0]), 5.0)

    def test_linf_norm(self):
        from anuga.utilities.norms import linf_norm
        self.assertAlmostEqual(linf_norm([1.0, -5.0, 3.0]), 5.0)


class Test_file_length(unittest.TestCase):
    """Tests for anuga.lib.file_length."""

    def test_file_length(self):
        import tempfile
        import os
        from anuga.lib.file_length import file_length
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write("line1\nline2\nline3\n")
            fname = f.name
        try:
            n = file_length(fname)
            self.assertEqual(n, 3)
        finally:
            os.unlink(fname)


class Test_metadata_imports(unittest.TestCase):
    """Cover module-level statements in __metadata__, rain/__init__."""

    def test_metadata_importable(self):
        import anuga.__metadata__ as meta
        self.assertIsNotNone(meta.__version__)

    def test_rain_importable(self):
        import anuga.rain
        self.assertIsNotNone(anuga.rain)


class Test_numerical_tools_extra(unittest.TestCase):
    """Cover previously uncovered lines in numerical_tools.py."""

    from anuga.utilities.numerical_tools import (
        safe_acos, cov, corr, err, histogram, create_bins, gradient2_python
    )

    def test_safe_acos_clamped_negative(self):
        """safe_acos clamps x just below -1 (lines 36-37)."""
        from anuga.utilities.numerical_tools import safe_acos, get_machine_precision
        eps = get_machine_precision()
        x = -1.0 - eps / 2  # just inside tolerance
        result = safe_acos(x)
        import math
        self.assertAlmostEqual(result, math.pi)

    def test_safe_acos_raises_below_tolerance(self):
        """safe_acos raises ValueError when x << -1 (lines 34-35)."""
        from anuga.utilities.numerical_tools import safe_acos
        with self.assertRaises(ValueError):
            safe_acos(-2.0)

    def test_safe_acos_clamped_positive(self):
        """safe_acos clamps x just above 1 (lines 42-43)."""
        from anuga.utilities.numerical_tools import safe_acos, get_machine_precision
        eps = get_machine_precision()
        x = 1.0 + eps / 2
        result = safe_acos(x)
        self.assertAlmostEqual(result, 0.0)

    def test_cov_single_arg(self):
        """cov(x) computes variance (lines 143-156)."""
        from anuga.utilities.numerical_tools import cov
        x = [1.0, 2.0, 3.0]
        result = cov(x)
        self.assertAlmostEqual(result, 2.0 / 3.0)

    def test_cov_two_args(self):
        """cov(x, y) computes covariance."""
        from anuga.utilities.numerical_tools import cov
        x = [1.0, 2.0, 3.0]
        y = [4.0, 5.0, 6.0]
        result = cov(x, y)
        self.assertAlmostEqual(result, 2.0 / 3.0)

    def test_corr_basic(self):
        """corr(x, y) computes correlation (lines 214-226)."""
        from anuga.utilities.numerical_tools import corr
        x = [1.0, 2.0, 3.0]
        y = [1.0, 2.0, 3.0]
        result = corr(x, y)
        self.assertAlmostEqual(result, 1.0)

    def test_corr_zero_variance(self):
        """corr with zero variance returns 0 (line 222)."""
        from anuga.utilities.numerical_tools import corr
        x = [1.0, 1.0, 1.0]
        result = corr(x, x)
        self.assertEqual(result, 0)

    def test_err_default_y_zero(self):
        """err with y=0 (default) covers 178->181 branch (y stays 0)."""
        from anuga.utilities.numerical_tools import err
        import numpy as np
        x = np.array([3.0, 4.0])
        result = err(x, relative=False)
        self.assertAlmostEqual(result, 5.0)

    def test_err_max_norm(self):
        """err with n=None uses max norm (lines 189-195)."""
        from anuga.utilities.numerical_tools import err
        # Use scalar y to avoid array truth value issue
        result = err([1.0, 2.0, 3.0], y=0, n=None, relative=False)
        self.assertAlmostEqual(result, 3.0)

    def test_err_max_norm_relative(self):
        """err with n=None, scalar y, relative=True (line 192-195)."""
        from anuga.utilities.numerical_tools import err
        # pass y=1.0 (scalar) to avoid array truth value ambiguity
        result = err([2.0, 3.0], y=1.0, n=None, relative=True)
        # max(abs([1,2]))/max(abs(1)) = 2/1 = 2
        self.assertAlmostEqual(result, 2.0)

    def test_histogram_relative(self):
        """histogram with relative=True (line 282)."""
        from anuga.utilities.numerical_tools import histogram
        import numpy as np
        a = [1, 2, 2, 3]
        bins = [1, 2, 3, 4]
        result = histogram(a, bins, relative=True)
        self.assertAlmostEqual(sum(result), 1.0)

    def test_create_bins_constant(self):
        """create_bins with constant data returns one bin (line 296)."""
        from anuga.utilities.numerical_tools import create_bins
        result = create_bins([5.0, 5.0, 5.0])
        self.assertEqual(len(result), 1)

    def test_create_bins_no_count(self):
        """create_bins without number_of_bins defaults to 10 (lines 298-301)."""
        from anuga.utilities.numerical_tools import create_bins
        result = create_bins([1.0, 2.0, 3.0, 10.0])
        self.assertEqual(len(result), 10)

    def test_gradient2_python(self):
        """gradient2_python covers lines 356-361."""
        from anuga.utilities.numerical_tools import gradient2_python
        a, b = gradient2_python(0.0, 0.0, 1.0, 0.0, 0.0, 1.0)
        self.assertAlmostEqual(a, 1.0)
        self.assertAlmostEqual(b, 0.0)


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_Numerical_Tools)
    runner = unittest.TextTestRunner()
    runner.run(suite)

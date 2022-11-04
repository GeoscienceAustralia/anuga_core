"""Test suite for interpolate2d.py

Ole Nielsen, 2011
"""

__author__ = 'Ole Nielsen <ole.moller.nielsen@gmail.com>'
__revision__ = '$Format:%H$'
__date__ = '01/11/2011'
__license__ = 'GPL v3'
__copyright__ = 'Copyright 2012, Australia Indonesia Facility for '
__copyright__ += 'Disaster Reduction'

import unittest
import numpy

from anuga.fit_interpolate.interpolate2d import interpolate2d


# Auxiliary functions for testing
def nanallclose(x, y, rtol=1.0e-5, atol=1.0e-8):
    """Numpy allclose function which allows NaN

    Args:
        * x, y: Either scalars or numpy arrays

    Returns:
        * True or False

    Note:
        Returns True if all non-nan elements pass.
    """

    xn = numpy.isnan(x)
    yn = numpy.isnan(y)
    if numpy.any(xn != yn):
        # Presence of NaNs is not the same in x and y
        return False

    if numpy.all(xn):
        # Everything is NaN.
        # This will also take care of x and y being NaN scalars
        return True

    # Filter NaN's out
    if numpy.any(xn):
        x = x[numpy.logical_not(xn)]
        y = y[numpy.logical_not(yn)]

    # Compare non NaN's and return
    return numpy.allclose(x, y, rtol=rtol, atol=atol)


def axes2points(x, y):
    """Generate all combinations of grid point coordinates from x and y axes

    Args:
        * x: x coordinates (array)
        * y: y coordinates (array)

    Returns:
        * P: Nx2 array consisting of coordinates for all
             grid points defined by x and y axes. The x coordinate
             will vary the fastest to match the way 2D numpy
             arrays are laid out by default ('C' order). That way,
             the x and y coordinates will match a corresponding
             2D array A when flattened (A.flat[:] or A.reshape(-1))

    Note:
        Example

        x = [1, 2, 3]
        y = [10, 20]

        P = [[1, 10],
             [2, 10],
             [3, 10],
             [1, 20],
             [2, 20],
             [3, 20]]
    """

    # Reverse y coordinates to have them start at bottom of array
    y = numpy.flipud(y)

    # Repeat x coordinates for each y (fastest varying)
    X = numpy.kron(numpy.ones(len(y)), x)

    # Repeat y coordinates for each x (slowest varying)
    Y = numpy.kron(y, numpy.ones(len(x)))

    # Check
    N = len(X)
    assert len(Y) == N

    # Create Nx2 array of x and y coordinates
    X = numpy.reshape(X, (N, 1))
    Y = numpy.reshape(Y, (N, 1))
    P = numpy.concatenate((X, Y), axis=1)

    # Return
    return P


def linear_function(x, y):
    """Auxiliary function for use with interpolation test
    """

    return x + y / 2.0


class Test_interpolate(unittest.TestCase):

    def test_linear_interpolation_basic(self):
        """Interpolation library works for linear function - basic test
        """

        # Define pixel centers along each direction
        x = [1.0, 2.0, 4.0]
        y = [5.0, 9.0]

        # Define ny by nx array with corresponding values
        A = numpy.zeros((len(x), len(y)))

        # Define values for each x, y pair as a linear function
        for i in range(len(x)):
            for j in range(len(y)):
                A[i, j] = linear_function(x[i], y[j])

        # Test first that original points are reproduced correctly
        for i, xi in enumerate(x):
            for j, eta in enumerate(y):
                val = interpolate2d(x, y, A, [(xi, eta)], mode='linear')[0]
                ref = linear_function(xi, eta)
                assert numpy.allclose(val, ref, rtol=1e-12, atol=1e-12)

        # Then test that genuinly interpolated points are correct
        xis = numpy.linspace(x[0], x[-1], 10)
        etas = numpy.linspace(y[0], y[-1], 10)
        points = axes2points(xis, etas)

        vals = interpolate2d(x, y, A, points, mode='linear')
        refs = linear_function(points[:, 0], points[:, 1])
        assert numpy.allclose(vals, refs, rtol=1e-12, atol=1e-12)

    def test_constant_interpolation_basic(self):
        """Interpolation library works for piecewise constant function
        """

        # Define pixel centers along each direction
        x = numpy.array([1.0, 2.0, 4.0])
        y = numpy.array([5.0, 9.0])

        # Define ny by nx array with corresponding values
        A = numpy.zeros((len(x), len(y)))

        # Define values for each x, y pair as a linear function
        for i in range(len(x)):
            for j in range(len(y)):
                A[i, j] = linear_function(x[i], y[j])

        # Then test that interpolated points are always assigned value of
        # closest neighbour
        xis = numpy.linspace(x[0], x[-1], 10)
        etas = numpy.linspace(y[0], y[-1], 10)
        points = axes2points(xis, etas)

        vals = interpolate2d(x, y, A, points, mode='constant')

        # Find upper neighbours for each interpolation point
        xi = points[:, 0]
        eta = points[:, 1]
        idx = numpy.searchsorted(x, xi, side='left')
        idy = numpy.searchsorted(y, eta, side='left')

        # Get the four neighbours for each interpolation point
        x0 = x[idx - 1]
        x1 = x[idx]
        y0 = y[idy - 1]
        y1 = y[idy]

        z00 = A[idx - 1, idy - 1]
        z01 = A[idx - 1, idy]
        z10 = A[idx, idy - 1]
        z11 = A[idx, idy]

        # Location coefficients
        alpha = (xi - x0)/ (x1 - x0)
        beta = (eta - y0)/ (y1 - y0)

        refs = numpy.zeros(len(vals))
        for i in range(len(refs)):
            if alpha[i] < 0.5 and beta[i] < 0.5:
                refs[i] = z00[i]

            if alpha[i] >= 0.5 and beta[i] < 0.5:
                refs[i] = z10[i]

            if alpha[i] < 0.5 and beta[i] >= 0.5:
                refs[i] = z01[i]

            if alpha[i] >= 0.5 and beta[i] >= 0.5:
                refs[i] = z11[i]

        assert numpy.allclose(vals, refs, rtol=1e-12, atol=1e-12)

    def test_linear_interpolation_range(self):
        """Interpolation library works for linear function - a range of cases
        """

        for x in [[1.0, 2.0, 4.0], [-20, -19, 0], numpy.arange(200) + 1000]:
            for y in [[5.0, 9.0], [100, 200, 10000]]:

                # Define ny by nx array with corresponding values
                A = numpy.zeros((len(x), len(y)))

                # Define values for each x, y pair as a linear function
                for i in range(len(x)):
                    for j in range(len(y)):
                        A[i, j] = linear_function(x[i], y[j])

                # Test that linearly interpolated points are correct
                xis = numpy.linspace(x[0], x[-1], 100)
                etas = numpy.linspace(y[0], y[-1], 100)
                points = axes2points(xis, etas)

                vals = interpolate2d(x, y, A, points, mode='linear')
                refs = linear_function(points[:, 0], points[:, 1])
                assert numpy.allclose(vals, refs, rtol=1e-12, atol=1e-12)

    def test_linear_interpolation_nan_points(self):
        """Interpolation library works with interpolation points being NaN

        This is was the reason for bug reported in:
        https://github.com/AIFDR/riab/issues/155
        """

        # Define pixel centers along each direction
        x = [1.0, 2.0, 4.0]
        y = [5.0, 9.0]

        # Define ny by nx array with corresponding values
        A = numpy.zeros((len(x), len(y)))

        # Define values for each x, y pair as a linear function
        for i in range(len(x)):
            for j in range(len(y)):
                A[i, j] = linear_function(x[i], y[j])

        # Then test that interpolated points can contain NaN
        xis = numpy.linspace(x[0], x[-1], 10)
        etas = numpy.linspace(y[0], y[-1], 10)
        
        xis[6:7] = numpy.nan
        etas[3] = numpy.nan
        points = axes2points(xis, etas)
        
        vals = interpolate2d(x, y, A, points, mode='linear')
        refs = linear_function(points[:, 0], points[:, 1])
        
        assert nanallclose(vals, refs, rtol=1e-12, atol=1e-12)

    def test_linear_interpolation_nan_array(self):
        """Interpolation library works (linear mode) with grid points being NaN
        """

        # Define pixel centers along each direction
        x = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
        y = [4.0, 5.0, 7.0, 9.0, 11.0, 13.0]

        # Define ny by nx array with corresponding values
        A = numpy.zeros((len(x), len(y)))

        # Define values for each x, y pair as a linear function
        for i in range(len(x)):
            for j in range(len(y)):
                A[i, j] = linear_function(x[i], y[j])
        A[2, 3] = numpy.nan  # (x=2.0, y=9.0): NaN

        # Then test that interpolated points can contain NaN
        xis = numpy.linspace(x[0], x[-1], 12)
        etas = numpy.linspace(y[0], y[-1], 10)
        points = axes2points(xis, etas)


        vals = interpolate2d(x, y, A, points, mode='linear')
        refs = linear_function(points[:, 0], points[:, 1])

        # Set reference result with expected NaNs and compare
        for i, (xi, eta) in enumerate(points):
            if (1.0 < xi <= 3.0) and (7.0 < eta <= 11.0):
                refs[i] = numpy.nan

        assert nanallclose(vals, refs, rtol=1e-12, atol=1e-12)

    def test_interpolation_random_array_and_nan(self):
        """Interpolation library (constant and linear) works with NaN
        """

        # Define pixel centers along each direction
        x = numpy.arange(20) * 1.0
        y = numpy.arange(25) * 1.0

        # Define ny by nx array with corresponding values
        A = numpy.zeros((len(x), len(y)))

        # Define arbitrary values for each x, y pair
        numpy.random.seed(17)
        A = numpy.random.random((len(x), len(y))) * 10

        # Create islands of NaN
        A[5, 13] = numpy.nan
        A[6, 14] = A[6, 18] = numpy.nan
        A[7, 14:18] = numpy.nan
        A[8, 13:18] = numpy.nan
        A[9, 12:19] = numpy.nan
        A[10, 14:17] = numpy.nan
        A[11, 15] = numpy.nan

        A[15, 5:6] = numpy.nan

        # Creat interpolation points
        xis = numpy.linspace(x[0], x[-1], 39)   # Hit all mid points
        etas = numpy.linspace(y[0], y[-1], 73)  # Hit thirds
        points = axes2points(xis, etas)

        for mode in ['linear', 'constant']:
            vals = interpolate2d(x, y, A, points, mode=mode)

            # Calculate reference result with expected NaNs and compare
            i = j = 0
            for k, (xi, eta) in enumerate(points):

                # Find indices of nearest higher value in x and y
                i = numpy.searchsorted(x, xi)
                j = numpy.searchsorted(y, eta)

                if i > 0 and j > 0:

                    # Get four neigbours
                    A00 = A[i - 1, j - 1]
                    A01 = A[i - 1, j]
                    A10 = A[i, j - 1]
                    A11 = A[i, j]

                    if numpy.allclose(xi, x[i]):
                        alpha = 1.0
                    else:
                        alpha = 0.5

                    if numpy.allclose(eta, y[j]):
                        beta = 1.0
                    else:
                        beta = eta - y[j - 1]

                    if mode == 'linear':
                        if numpy.any(numpy.isnan([A00, A01, A10, A11])):
                            ref = numpy.nan
                        else:
                            ref = (A00 * (1 - alpha) * (1 - beta) +
                                   A01 * (1 - alpha) * beta +
                                   A10 * alpha * (1 - beta) +
                                   A11 * alpha * beta)
                    elif mode == 'constant':
                        assert alpha >= 0.5  # Only case in this test

                        if beta < 0.5:
                            ref = A10
                        else:
                            ref = A11
                    else:
                        msg = 'Unknown mode: %s' % mode
                        raise Exception(msg)

                    #print i, j, xi, eta, alpha, beta, vals[k], ref
                    assert nanallclose(vals[k], ref, rtol=1e-12, atol=1e-12)

    def test_linear_interpolation_outside_domain(self):
        """Interpolation library sensibly handles values outside the domain
        """

        # Define pixel centers along each direction
        x = [1.0, 2.0, 4.0]
        y = [5.0, 9.0]

        # Define ny by nx array with corresponding values
        A = numpy.zeros((len(x), len(y)))

        # Define values for each x, y pair as a linear function
        for i in range(len(x)):
            for j in range(len(y)):
                A[i, j] = linear_function(x[i], y[j])

        # Simple example first for debugging
        xis = numpy.linspace(0.9, 4.0, 4)
        etas = numpy.linspace(5, 9.1, 3)
        points = axes2points(xis, etas)
        refs = linear_function(points[:, 0], points[:, 1])

        vals = interpolate2d(x, y, A, points, mode='linear',
                             bounds_error=False)
        msg = ('Length of interpolation points %i differs from length '
               'of interpolated values %i' % (len(points), len(vals)))
        assert len(points) == len(vals), msg
        for i, (xi, eta) in enumerate(points):
            if xi < x[0] or xi > x[-1] or eta < y[0] or eta > y[-1]:
                assert numpy.isnan(vals[i])
            else:
                msg = ('Got %.15f for (%f, %f), expected %.15f'
                       % (vals[i], xi, eta, refs[i]))
                assert numpy.allclose(vals[i], refs[i],
                                      rtol=1.0e-12, atol=1.0e-12), msg

        # Try a range of combinations of points outside domain
        # with error_bounds True
        for lox in [x[0], x[0] - 1]:
            for hix in [x[-1], x[-1] + 1]:
                for loy in [y[0], y[0] - 1]:
                    for hiy in [y[-1], y[-1] + 1]:

                        # Then test that points outside domain can be handled
                        xis = numpy.linspace(lox, hix, 4)
                        etas = numpy.linspace(loy, hiy, 4)
                        points = axes2points(xis, etas)

                        if lox < x[0] or hix > x[-1] or \
                                loy < x[0] or hiy > y[-1]:
                            try:
                                vals = interpolate2d(x, y, A, points,
                                                     mode='linear',
                                                     bounds_error=True)
                            except Exception as e:
                                pass
                            else:
                                msg = 'Should have raise bounds error'
                                raise Exception(msg)

        # Try a range of combinations of points outside domain with
        # error_bounds False
        for lox in [x[0], x[0] - 1, x[0] - 10]:
            for hix in [x[-1], x[-1] + 1, x[-1] + 5]:
                for loy in [y[0], y[0] - 1, y[0] - 10]:
                    for hiy in [y[-1], y[-1] + 1, y[-1] + 10]:

                        # Then test that points outside domain can be handled
                        xis = numpy.linspace(lox, hix, 10)
                        etas = numpy.linspace(loy, hiy, 10)
                        points = axes2points(xis, etas)
                        refs = linear_function(points[:, 0], points[:, 1])
                        vals = interpolate2d(x, y, A, points,
                                             mode='linear', bounds_error=False)

                        assert len(points) == len(vals), msg
                        for i, (xi, eta) in enumerate(points):
                            if xi < x[0] or xi > x[-1] or\
                                    eta < y[0] or eta > y[-1]:
                                msg = 'Expected NaN for %f, %f' % (xi, eta)
                                assert numpy.isnan(vals[i]), msg
                            else:
                                msg = ('Got %.15f for (%f, %f), expected '
                                       '%.15f' % (vals[i], xi, eta, refs[i]))
                                assert numpy.allclose(vals[i], refs[i],
                                                      rtol=1.0e-12,
                                                      atol=1.0e-12), msg


if __name__ == '__main__':
    suite = unittest.makeSuite(Test_interpolate, 'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)

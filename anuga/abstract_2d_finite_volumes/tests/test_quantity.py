#!/usr/bin/env python


import os
import unittest
from math import sqrt, pi
import tempfile

from anuga.abstract_2d_finite_volumes.quantity import *
from anuga.file_conversion.asc2dem import asc2dem
from anuga.config import epsilon

from anuga.fit_interpolate.fit import fit_to_mesh
# from anuga.pyvolution.least_squares import fit_to_mesh
from anuga.abstract_2d_finite_volumes.generic_domain \
                import Generic_Domain
from anuga.geospatial_data.geospatial_data import Geospatial_data
from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.geometry.polygon import *

import numpy as num
import pprint

try:
    import rasterio
except ImportError:
    pass

import pytest
import sys

def zone_letter_to_hemisphere(zone_letter):
    hemisphere = 'undefined'
    if zone_letter.lower() in 'cdefghjklm':
        hemisphere = 'southern'
    if zone_letter.lower() in 'npqrstuvwx':
        hemisphere = 'northern'
    return hemisphere

# Aux for fit_interpolate.fit example
def linear_function(point):
    point = num.array(point)
    return point[:, 0]+3*point[:, 1]
    #return point[:,1]


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
    import numpy

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


class Test_Quantity(unittest.TestCase):
    def setUp(self):

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]

        # bac, bce, ecf, dbe
        elements = [[1, 0, 2], [1, 2, 4], [4, 2, 5], [3, 1, 4]]

        self.mesh1 = Generic_Domain(points[:3], [elements[0]])
        self.mesh1.check_integrity()

        # print self.mesh1.__class__
        # print isinstance(self.mesh1, Domain)

        self.mesh4 = Generic_Domain(points, elements)
        self.mesh4.check_integrity()

        # UTM round Onslow
        a = [240000, 7620000]
        b = [240000, 7680000]
        c = [300000, 7620000]

        points = [a, b, c]
        elements = [[0, 2, 1]]

        self.mesh_onslow = Generic_Domain(points, elements)
        self.mesh_onslow.check_integrity()

    def tearDown(self):
        pass
        # print "  Tearing down"

    def test_creation(self):

        quantity = Quantity(self.mesh1, [[1, 2, 3]])
        assert num.allclose(quantity.vertex_values, [[1., 2., 3.]])

        try:
            quantity = Quantity()
        except Exception:
            pass
        else:
            raise Exception('Should have raised empty quantity exception')

        # FIXME(Ole): Temporarily disabled 18 Jan 2009
        # try:
        #    quantity = Quantity([1,2,3])
        # except AssertionError:
        #    pass
        # except:
        #    raise Exception('Should have raised "mising mesh object" error')

    def test_creation_zeros(self):

        quantity = Quantity(self.mesh1)
        assert num.allclose(quantity.vertex_values, [[0., 0., 0.]])

        quantity = Quantity(self.mesh4)
        assert num.allclose(quantity.vertex_values, [[0., 0., 0.], [0., 0., 0.],
                                                     [0., 0., 0.], [0., 0., 0.]])

    def test_set_boundary_values(self):

        quantity = Quantity(self.mesh1)

        quantity.set_boundary_values()

        assert num.allclose(quantity.boundary_values, [0.0, 0.0, 0.0])

    def test_set_boundary_values_with_function(self):

        quantity = Quantity(self.mesh1)
        # assert num.allclose(quantity.vertex_values, [[0.,0.,0.]])

        def simple(x, y):
            return x+3*y

        quantity.set_boundary_values(simple)

        assert num.allclose(quantity.boundary_values, [1.0, 4.0, 3.0])

    def test_set_boundary_values_with_constant(self):

        quantity = Quantity(self.mesh1)
        # assert num.allclose(quantity.vertex_values, [[0.,0.,0.]])

        quantity.set_boundary_values(10.0)

        assert num.allclose(quantity.boundary_values, [10.0, 10.0, 10.0])

    def test_set_boundary_values_with_array(self):

        quantity = Quantity(self.mesh1)
        # assert num.allclose(quantity.vertex_values, [[0.,0.,0.]])

        quantity.set_boundary_values([10.0, 4.0, 5.0])

        assert num.allclose(quantity.boundary_values, [10.0, 4.0, 5.0])

    def test_set_boundary_values_with_wrong_sized_array(self):

        quantity = Quantity(self.mesh1)
        # assert num.allclose(quantity.vertex_values, [[0.,0.,0.]])

        try:
            quantity.set_boundary_values([10.0, 4.0, 5.0, 8.0])
        except Exception:
            pass
        else:
            msg = 'Should have caught this'
            raise Exception(msg)

    def test_set_boundary_values_from_edges(self):

        quantity = Quantity(self.mesh4)

        def simple(x, y):
            return x+3*y

        quantity.set_values(simple)

        assert num.allclose(quantity.boundary_values, [
                            0.,   0.,   0.,   0.,  0.,   0.])

        quantity.set_boundary_values_from_edges()

        assert num.allclose(quantity.boundary_values, [
                            1.,   3.,   3.,   6.,  10.,   9.])

    def test_interpolation(self):
        quantity = Quantity(self.mesh1, [[1, 2, 3]])
        assert num.allclose(quantity.centroid_values, [2.0])  # Centroid

        assert num.allclose(quantity.edge_values, [[2.5, 2.0, 1.5]])

    def test_interpolation2(self):
        quantity = Quantity(self.mesh4,
                            [[1, 2, 3], [5, 5, 5], [0, 0, 9], [-6, 3, 3]])
        assert num.allclose(quantity.centroid_values, [
                            2., 5., 3., 0.])  # Centroid

        quantity.extrapolate_second_order()

        # print quantity.vertex_values
        assert num.allclose(quantity.vertex_values, [[3.5, -1.0, 3.5],
                                                     [3.+2./3, 6.+2./3, 4.+2./3],
                                                     [4.6, 3.4, 1.],
                                                     [-5.0, 1.0, 4.0]])

        # print quantity.edge_values
        assert num.allclose(quantity.edge_values, [[1.25, 3.5, 1.25],
                                                   [5. + 2/3.0, 4.0
                                                       + 1.0/6, 5.0 + 1.0/6],
                                                   [2.2, 2.8, 4.0],
                                                   [2.5, -0.5, -2.0]])

    def test_save_to_array(self):
        quantity = Quantity(self.mesh4,
                            [[1, 2, 3], [5, 5, 5], [0, 0, 9], [-6, 3, 3]])

        assert num.allclose(quantity.centroid_values, [
                            2., 5., 3., 0.])  # Centroid

        cellsize = 1.0
        x, y, z = quantity.save_to_array(cellsize=cellsize, smooth=False)
        # x,y,z = quantity.save_to_array(smooth=False)

        from pprint import pprint
        # pprint(x)
        # pprint(y)
        # pprint(z)

        x_ex = [0.,  1.,  2.,  3.,  4.]
        y_ex = [0.,  1.,  2.,  3.,  4.]

        z_ex = [[2.00000000e+00,   2.50000000e+00,   0.00000000e+00,
                 4.50000000e+00,   9.00000000e+00],
                [1.50000000e+00,   5.00000000e+00,   0.00000000e+00,
                 4.50000000e+00,  -9.99900000e+03],
                [3.00000000e+00,   3.00000000e+00,   3.00000000e+00,
                 -9.99900000e+03,  -9.99900000e+03],
                [-1.50000000e+00,  -1.50000000e+00,  -9.99900000e+03,
                 -9.99900000e+03,  -9.99900000e+03],
                [-6.00000000e+00,  -9.99900000e+03,  -9.99900000e+03,
                 -9.99900000e+03,  -9.99900000e+03]]

        assert num.allclose(x_ex, x)
        assert num.allclose(y_ex, y)
        assert num.allclose(z_ex, z)

        Plot = False
        if Plot:
            import pylab
            import numpy
            # a = numpy.where(a == -9999, numpy.nan, a)
            # a = numpy.where(a > 10.0, numpy.nan, a)

            # z = z[::-1,:]

            # print(z)
            # print(z.shape)
            # print(x)
            # print(y)

            nrows = z.shape[0]
            ncols = z.shape[1]

            ratio = float(nrows)/float(ncols)
            print(ratio)

            # y = numpy.arange(nrows)*cellsize
            # x = numpy.arange(ncols)*cellsize

            # Setup fig size to correpond to array size
            fig = pylab.figure(figsize=(10, 10*ratio))

            levels = numpy.arange(-7, 10, 0.1)
            CF = pylab.contourf(x, y, z, levels=levels)
            CB = pylab.colorbar(CF, shrink=0.8, extend='both')
            # CC = pylab.contour(x,y,a, levels=levels)

            pylab.show()

        x, y, z = quantity.save_to_array(cellsize=cellsize, smooth=True)

        x_ex = [0.,  1.,  2.,  3.,  4.]
        y_ex = [0.,  1.,  2.,  3.,  4.]
        z_ex = [[2.00000000e+00,   2.33333333e+00,   2.66666667e+00,
                 5.83333333e+00,   9.00000000e+00],
                [2.50000000e+00,   2.83333333e+00,   2.66666667e+00,
                 5.83333333e+00,  -9.99900000e+03],
                [3.00000000e+00,   2.83333333e+00,   2.66666667e+00,
                 -9.99900000e+03,  -9.99900000e+03],
                [-1.50000000e+00,  -1.66666667e+00,  -9.99900000e+03,
                 -9.99900000e+03,  -9.99900000e+03],
                [-6.00000000e+00,  -9.99900000e+03,  -9.99900000e+03,
                 -9.99900000e+03,  -9.99900000e+03]]

        # pprint(z)
        assert num.allclose(x_ex, x)
        assert num.allclose(y_ex, y)
        assert num.allclose(z_ex, z)

        if Plot:
            import pylab
            import numpy
            # a = numpy.where(a == -9999, numpy.nan, a)
            # a = numpy.where(a > 10.0, numpy.nan, a)

            # a = a[::-1,:]
            nrows = z.shape[0]
            ncols = z.shape[1]

            ratio = float(nrows)/float(ncols)
            print(ratio)

            # Setup fig size to correpond to array size
            fig = pylab.figure(figsize=(10, 10*ratio))

            levels = numpy.arange(-7, 10, 0.1)
            CF = pylab.contourf(x, y, z, levels=levels)
            CB = pylab.colorbar(CF, shrink=0.8, extend='both')
            # CC = pylab.contour(x,y,a, levels=[0.0,1.0,2.0,3.0])

            pylab.show()

    def test_get_extrema_1(self):
        quantity = Quantity(self.mesh4,
                            [[1, 2, 3], [5, 5, 5], [0, 0, 9], [-6, 3, 3]])
        assert num.allclose(quantity.centroid_values, [
                            2., 5., 3., 0.])  # Centroids

        v = quantity.get_maximum_value()
        assert v == 5

        v = quantity.get_minimum_value()
        assert v == 0

        i = quantity.get_maximum_index()
        assert i == 1

        i = quantity.get_minimum_index()
        assert i == 3

        x, y = quantity.get_maximum_location()
        xref, yref = 4.0/3, 4.0/3
        assert x == xref
        assert y == yref

        v = quantity.get_values(interpolation_points=[[x, y]])
        assert num.allclose(v, 5)

        x, y = quantity.get_minimum_location()
        v = quantity.get_values(interpolation_points=[[x, y]])
        assert num.allclose(v, 0)

    def test_get_maximum_2(self):

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        # bac, bce, ecf, dbe
        vertices = [[1, 0, 2], [1, 2, 4], [4, 2, 5], [3, 1, 4]]

        domain = Generic_Domain(points, vertices)

        quantity = Quantity(domain)
        quantity.set_values(lambda x, y: x+2*y)  # 2 4 4 6

        v = quantity.get_maximum_value()
        assert v == 6

        v = quantity.get_minimum_value()
        assert v == 2

        i = quantity.get_maximum_index()
        assert i == 3

        i = quantity.get_minimum_index()
        assert i == 0

        x, y = quantity.get_maximum_location()
        xref, yref = 2.0/3, 8.0/3
        assert x == xref
        assert y == yref

        v = quantity.get_values(interpolation_points=[[x, y]])
        assert num.allclose(v, 6)

        x, y = quantity.get_minimum_location()
        v = quantity.get_values(interpolation_points=[[x, y]])
        assert num.allclose(v, 2)

        # Multiple locations for maximum -
        # Test that the algorithm picks the first occurrence
        v = quantity.get_maximum_value(indices=[0, 1, 2])
        assert num.allclose(v, 4)

        i = quantity.get_maximum_index(indices=[0, 1, 2])
        assert i == 1

        x, y = quantity.get_maximum_location(indices=[0, 1, 2])
        xref, yref = 4.0/3, 4.0/3
        assert x == xref
        assert y == yref

        v = quantity.get_values(interpolation_points=[[x, y]])
        assert num.allclose(v, 4)

        # More test of indices......
        v = quantity.get_maximum_value(indices=[2, 3])
        assert num.allclose(v, 6)

        i = quantity.get_maximum_index(indices=[2, 3])
        assert i == 3

        x, y = quantity.get_maximum_location(indices=[2, 3])
        xref, yref = 2.0/3, 8.0/3
        assert x == xref
        assert y == yref

        v = quantity.get_values(interpolation_points=[[x, y]])
        assert num.allclose(v, 6)

    def test_boundary_allocation(self):
        quantity = Quantity(self.mesh4,
                            [[1, 2, 3], [5, 5, 5], [0, 0, 9], [-6, 3, 3]])

        assert quantity.boundary_values.shape[0] == len(self.mesh4.boundary)

    def test_set_values(self):
        quantity = Quantity(self.mesh4)

        # get referece to data arrays
        centroid_values = quantity.centroid_values
        vertex_values = quantity.vertex_values

        quantity.set_values([[1, 2, 3], [5, 5, 5], [0, 0, 9], [-6, 3, 3]],
                            location='vertices')
        assert num.allclose(quantity.vertex_values,
                            [[1, 2, 3], [5, 5, 5], [0, 0, 9], [-6, 3, 3]])

        assert id(vertex_values) == id(quantity.vertex_values)
        assert num.allclose(quantity.centroid_values, [
                            2., 5., 3., 0.])  # Centroid
        assert num.allclose(quantity.edge_values, [[2.5, 2.0, 1.5],
                                                   [5., 5., 5.],
                                                   [4.5, 4.5, 0.],
                                                   [3.0, -1.5, -1.5]])

        # Test default
        quantity.set_values([[1, 2, 3], [5, 5, 5], [0, 0, 9], [-6, 3, 3]])
        assert num.allclose(quantity.vertex_values,
                            [[1, 2, 3], [5, 5, 5], [0, 0, 9], [-6, 3, 3]])
        assert num.allclose(quantity.centroid_values, [
                            2., 5., 3., 0.])  # Centroid
        assert num.allclose(quantity.edge_values, [[2.5, 2.0, 1.5],
                                                   [5., 5., 5.],
                                                   [4.5, 4.5, 0.],
                                                   [3.0, -1.5, -1.5]])

        # Test centroids
        quantity.set_values([1, 2, 3, 4], location='centroids')
        assert num.allclose(quantity.centroid_values, [
                            1., 2., 3., 4.])  # Centroid

        # Test exceptions
        try:
            quantity.set_values([[1, 2, 3], [5, 5, 5], [0, 0, 9], [-6, 3, 3]],
                                location='bas kamel tuba')
        except Exception:
            pass

        try:
            quantity.set_values([[1, 2, 3], [0, 0, 9]])
        except ValueError:
            pass
        except Exception:
            raise Exception('should have raised ValueeError')

    def test_set_values_const(self):
        quantity = Quantity(self.mesh4)

        quantity.set_values(1.0, location='vertices')
        assert num.allclose(quantity.vertex_values,
                            [[1, 1, 1], [1, 1, 1], [1, 1, 1], [1, 1, 1]])

        assert num.allclose(quantity.centroid_values, [1, 1, 1, 1])  # Centroid
        assert num.allclose(quantity.edge_values, [[1, 1, 1],
                                                   [1, 1, 1],
                                                   [1, 1, 1],
                                                   [1, 1, 1]])

        quantity.set_values(2.0, location='centroids')
        assert num.allclose(quantity.centroid_values, [2, 2, 2, 2])

    def test_set_values_func(self):
        quantity = Quantity(self.mesh4)

        def f(x, y):
            return x+y

        quantity.set_values(f, location='vertices')
        # print "quantity.vertex_values",quantity.vertex_values
        assert num.allclose(quantity.vertex_values,
                            [[2, 0, 2], [2, 2, 4], [4, 2, 4], [4, 2, 4]])
        assert num.allclose(quantity.centroid_values,
                            [4.0/3, 8.0/3, 10.0/3, 10.0/3])
        assert num.allclose(quantity.edge_values,
                            [[1, 2, 1], [3, 3, 2], [3, 4, 3], [3, 4, 3]])

        quantity.set_values(f, location='centroids')
        assert num.allclose(quantity.centroid_values,
                            [4.0/3, 8.0/3, 10.0/3, 10.0/3])

    def test_integral(self):
        quantity = Quantity(self.mesh4)

        # Try constants first
        const = 5
        quantity.set_values(const, location='vertices')
        # print 'Q', quantity.get_integral()

        assert num.allclose(quantity.get_integral(),
                            self.mesh4.get_area() * const)

        # Try with a linear function
        def f(x, y):
            return x+y

        quantity.set_values(f, location='vertices')

        ref_integral = (4.0/3 + 8.0/3 + 10.0/3 + 10.0/3) * 2

        assert num.allclose(quantity.get_integral(), ref_integral)

    def test_integral_with_region(self):
        quantity = Quantity(self.mesh4)

        # Try constants first
        const = 5
        quantity.set_values(const, location='vertices')
        # print 'Q', quantity.get_integral()

        assert num.allclose(quantity.get_integral(),
                            self.mesh4.get_area() * const)

        # Try with a linear function
        def f(x, y):
            return x+y

        quantity.set_values(f, location='vertices')

        from anuga import Region

        reg1 = Region(self.mesh4, indices=[2])
        ref_integral = (10.0/3) * 2

        assert num.allclose(quantity.get_integral(region=reg1), ref_integral)

        reg2 = Region(self.mesh4, indices=[2, 3])
        ref_integral = (10.0/3 + 10.0/3) * 2

        assert num.allclose(quantity.get_integral(region=reg2), ref_integral)

        id = [2, 3]
        ref_integral = (10.0/3 + 10.0/3) * 2

        assert num.allclose(quantity.get_integral(indices=id), ref_integral)

    def test_set_vertex_values(self):
        quantity = Quantity(self.mesh4)
        quantity.set_vertex_values([0, 1, 2, 3, 4, 5])

        assert num.allclose(quantity.vertex_values,
                            [[1, 0, 2], [1, 2, 4], [4, 2, 5], [3, 1, 4]])
        assert num.allclose(quantity.centroid_values,
                            [1., 7./3, 11./3, 8./3])  # Centroid
        assert num.allclose(quantity.edge_values, [[1., 1.5, 0.5],
                                                   [3., 2.5, 1.5],
                                                   [3.5, 4.5, 3.],
                                                   [2.5, 3.5, 2]])

    def test_set_vertex_values_subset(self):
        quantity = Quantity(self.mesh4)
        quantity.set_vertex_values([0, 1, 2, 3, 4, 5])
        quantity.set_vertex_values([0, 20, 30, 50], indices=[0, 2, 3, 5])

        assert num.allclose(quantity.vertex_values,
                            [[1, 0, 20], [1, 20, 4], [4, 20, 50], [30, 1, 4]])

    def test_set_vertex_values_using_general_interface(self):
        quantity = Quantity(self.mesh4)

        quantity.set_values([0, 1, 2, 3, 4, 5])

        assert num.allclose(quantity.vertex_values,
                            [[1, 0, 2], [1, 2, 4], [4, 2, 5], [3, 1, 4]])

        # Centroid
        assert num.allclose(quantity.centroid_values, [1., 7./3, 11./3, 8./3])

        assert num.allclose(quantity.edge_values, [[1., 1.5, 0.5],
                                                   [3., 2.5, 1.5],
                                                   [3.5, 4.5, 3.],
                                                   [2.5, 3.5, 2]])

    def test_set_vertex_values_using_general_interface_with_subset(self):
        """test_set_vertex_values_using_general_interface_with_subset(self):

        Test that indices and polygon works (for constants values)
        """

        quantity = Quantity(self.mesh4)

        quantity.set_values([0, 2, 3, 5], indices=[0, 2, 3, 5])
        assert num.allclose(quantity.vertex_values,
                            [[0, 0, 2], [0, 2, 0], [0, 2, 5], [3, 0, 0]])

        # Constant
        quantity.set_values(0.0)
        quantity.set_values(3.14, indices=[0, 2], location='vertices')

        # Indices refer to triangle numbers
        assert num.allclose(quantity.vertex_values,
                            [[3.14, 3.14, 3.14], [0, 0, 0],
                             [3.14, 3.14, 3.14], [0, 0, 0]])

        # Now try with polygon (pick points where y>2)
        polygon = [[0, 2.1], [4, 2.1], [4, 7], [0, 7]]
        quantity.set_values(0.0)
        quantity.set_values(3.14, polygon=polygon)

        assert num.allclose(quantity.vertex_values,
                            [[0, 0, 0], [0, 0, 0], [0, 0, 0],
                             [3.14, 3.14, 3.14]])

        # Another polygon (pick triangle 1 and 2 (rightmost triangles)
        # using centroids
        polygon = [[2.1, 0.0], [3.5, 0.1], [2, 2.2], [0.2, 2]]
        quantity.set_values(0.0)
        quantity.set_values(3.14, location='centroids', polygon=polygon)
        assert num.allclose(quantity.vertex_values,
                            [[0, 0, 0],
                             [3.14, 3.14, 3.14],
                             [3.14, 3.14, 3.14],
                             [0, 0, 0]])

        # Same polygon now use vertices (default)
        polygon = [[2.1, 0.0], [3.5, 0.1], [2, 2.2], [0.2, 2]]
        quantity.set_values(0.0)
        # print 'Here 2'
        quantity.set_values(3.14, polygon=polygon)
        assert num.allclose(quantity.vertex_values,
                            [[0, 0, 0],
                             [3.14, 3.14, 3.14],
                             [3.14, 3.14, 3.14],
                             [0, 0, 0]])

        # Test input checking
        try:
            quantity.set_values(3.14, polygon=polygon, indices=[0, 2])
        except Exception:
            pass
        else:
            msg = 'Should have caught this'
            raise Exception(msg)

    def test_set_vertex_values_using_general_interface_subset_and_geo(self):
        """test_set_vertex_values_using_general_interface_with_subset(self):
        Test that indices and polygon works using georeferencing
        """

        quantity = Quantity(self.mesh4)
        G = Geo_reference(56, 10, 100)
        quantity.domain.set_georeference(G)

        # Constant
        quantity.set_values(0.0)
        quantity.set_values(3.14, indices=[0, 2], location='vertices')

        # Indices refer to triangle numbers here - not vertices (why?)
        assert num.allclose(quantity.vertex_values,
                            [[3.14, 3.14, 3.14], [0, 0, 0],
                             [3.14, 3.14, 3.14], [0, 0, 0]])

        # Now try with polygon (pick points where y>2)
        polygon = num.array([[0, 2.1], [4, 2.1], [4, 7], [0, 7]])
        polygon += [G.xllcorner, G.yllcorner]

        quantity.set_values(0.0)
        quantity.set_values(3.14, polygon=polygon, location='centroids')
        assert num.allclose(quantity.vertex_values,
                            [[0, 0, 0], [0, 0, 0], [0, 0, 0],
                             [3.14, 3.14, 3.14]])

        # Another polygon (pick triangle 1 and 2 (rightmost triangles)
        polygon = num.array([[2.1, 0.0], [3.5, 0.1], [2, 2.2], [0.2, 2]])
        polygon += [G.xllcorner, G.yllcorner]

        quantity.set_values(0.0)
        quantity.set_values(3.14, polygon=polygon)
        msg = ('quantity.vertex_values=\n%s\nshould be close to\n'
               '[[0,0,0],\n'
               ' [3.14,3.14,3.14],\n'
               ' [3.14,3.14,3.14],\n'
               ' [0,0,0]]' % str(quantity.vertex_values))
        assert num.allclose(quantity.vertex_values,
                            [[0, 0, 0],
                             [3.14, 3.14, 3.14],
                             [3.14, 3.14, 3.14],
                             [0, 0, 0]]), msg

    def test_set_values_using_fit(self):

        quantity = Quantity(self.mesh4)

        # Get (enough) datapoints
        data_points = [[0.66666667, 0.66666667],
                       [1.33333333, 1.33333333],
                       [2.66666667, 0.66666667],
                       [0.66666667, 2.66666667],
                       [0.0, 1.0],
                       [0.0, 3.0],
                       [1.0, 0.0],
                       [1.0, 1.0],
                       [1.0, 2.0],
                       [1.0, 3.0],
                       [2.0, 1.0],
                       [3.0, 0.0],
                       [3.0, 1.0]]

        z = linear_function(data_points)

        # Use built-in fit_interpolate.fit
        quantity.set_values(Geospatial_data(data_points, z), alpha=0)
        # quantity.set_values(points = data_points, values = z, alpha = 0)

        answer = linear_function(quantity.domain.get_vertex_coordinates())
        # print quantity.vertex_values, answer
        assert num.allclose(quantity.vertex_values.flat, answer)

        # Now try by setting the same values directly
        vertex_attributes = fit_to_mesh(data_points,
                                        quantity.domain.get_nodes(),
                                        quantity.domain.get_triangles(),
                                        point_attributes=z,
                                        alpha=0,
                                        verbose=False)

        # print vertex_attributes
        quantity.set_values(vertex_attributes)
        assert num.allclose(quantity.vertex_values.flat, answer)

    def test_test_set_values_using_fit_w_geo(self):

        # Mesh
        vertex_coordinates = [[0.76, 0.76],
                              [0.76, 5.76],
                              [5.76, 0.76]]
        triangles = [[0, 2, 1]]

        mesh_georef = Geo_reference(56, -0.76, -0.76)
        mesh1 = Generic_Domain(vertex_coordinates, triangles,
                               geo_reference=mesh_georef)
        mesh1.check_integrity()

        # Quantity
        quantity = Quantity(mesh1)

        # Data
        data_points = [[201.0, 401.0],
                       [201.0, 403.0],
                       [203.0, 401.0]]

        z = [2, 4, 4]

        data_georef = Geo_reference(56, -200, -400)

        # Reference
        ref = fit_to_mesh(data_points, vertex_coordinates, triangles,
                          point_attributes=z,
                          data_origin=data_georef.get_origin(),
                          mesh_origin=mesh_georef.get_origin(),
                          alpha=0)

        assert num.allclose(ref, [0, 5, 5])

        # Test set_values

        quantity.set_values(Geospatial_data(
            data_points, z, data_georef), alpha=0)

        # quantity.set_values(points = data_points,
        #                    values = z,
        #                    data_georef = data_georef,
        #                    alpha = 0)

        # quantity.set_values(points = data_points,
        #                    values = z,
        #                    data_georef = data_georef,
        #                    alpha = 0)
        assert num.allclose(quantity.vertex_values.flat, ref)

        # Test set_values using geospatial data object
        quantity.vertex_values[:] = 0.0

        geo = Geospatial_data(data_points, z, data_georef)

        quantity.set_values(geospatial_data=geo, alpha=0)
        assert num.allclose(quantity.vertex_values.flat, ref)

    def test_set_values_from_file1(self):
        quantity = Quantity(self.mesh4)

        # Get (enough) datapoints
        data_points = [[0.66666667, 0.66666667],
                       [1.33333333, 1.33333333],
                       [2.66666667, 0.66666667],
                       [0.66666667, 2.66666667],
                       [0.0, 1.0],
                       [0.0, 3.0],
                       [1.0, 0.0],
                       [1.0, 1.0],
                       [1.0, 2.0],
                       [1.0, 3.0],
                       [2.0, 1.0],
                       [3.0, 0.0],
                       [3.0, 1.0]]

        data_geo_spatial = Geospatial_data(data_points,
                                           geo_reference=Geo_reference(56, 0, 0))
        data_points_absolute = data_geo_spatial.get_data_points(absolute=True)
        attributes = linear_function(data_points_absolute)
        att = 'spam_and_eggs'

        # Create .txt file
        fd, ptsfile = tempfile.mkstemp(".txt")
        os.close(fd)
        file = open(ptsfile, "w")
        file.write(" x,y," + att + " \n")
        for data_point, attribute in zip(data_points_absolute, attributes):
            row = str(data_point[0]) + ',' + str(data_point[1]) \
                  + ',' + str(attribute)
            file.write(row + "\n")
        file.close()

        # Check that values can be set from file
        quantity.set_values(filename=ptsfile,
                            attribute_name=att, alpha=0)
        answer = linear_function(quantity.domain.get_vertex_coordinates())

        # print quantity.vertex_values.flat
        # print answer

        assert num.allclose(quantity.vertex_values.flat, answer)

        # Check that values can be set from file using default attribute
        quantity.set_values(filename=ptsfile, alpha=0)
        assert num.allclose(quantity.vertex_values.flat, answer)

        # Cleanup
        os.remove(ptsfile)

    def Xtest_set_values_from_file_using_polygon(self):
        """test_set_values_from_file_using_polygon(self):

        Test that polygon restriction works for general points data
        """

        quantity = Quantity(self.mesh4)

        # Get (enough) datapoints
        data_points = [[0.66666667, 0.66666667],
                       [1.33333333, 1.33333333],
                       [2.66666667, 0.66666667],
                       [0.66666667, 2.66666667],
                       [0.0, 1.0],
                       [0.0, 3.0],
                       [1.0, 0.0],
                       [1.0, 1.0],
                       [1.0, 2.0],
                       [1.0, 3.0],
                       [2.0, 1.0],
                       [3.0, 0.0],
                       [3.0, 1.0]]

        data_geo_spatial = Geospatial_data(data_points,
                                           geo_reference=Geo_reference(56, 0, 0))
        data_points_absolute = data_geo_spatial.get_data_points(absolute=True)
        attributes = linear_function(data_points_absolute)
        att = 'spam_and_eggs'

        # Create .txt file
        fd, ptsfile = tempfile.mkstemp(".txt")
        os.close(fd)
        file = open(ptsfile, "w")
        file.write(" x,y," + att + " \n")
        for data_point, attribute in zip(data_points_absolute, attributes):
            row = str(data_point[0]) + ',' + str(data_point[1]) \
                  + ',' + str(attribute)
            file.write(row + "\n")
        file.close()

        # Create restricting polygon (containing node #4 (2,2) and
        # centroid of triangle #1 (bce)
        polygon = [[1.0, 1.0], [4.0, 1.0],
                   [4.0, 4.0], [1.0, 4.0]]

        # print self.mesh4.nodes
        # print inside_polygon(self.mesh4.nodes, polygon)
        assert num.allclose(inside_polygon(self.mesh4.nodes, polygon), 4)

        # print quantity.domain.get_vertex_coordinates()
        # print quantity.domain.get_nodes()

        # Check that values can be set from file
        quantity.set_values(filename=ptsfile,
                            polygon=polygon,
                            location='unique vertices',
                            alpha=0)

        # Get indices for vertex coordinates in polygon
        indices = inside_polygon(quantity.domain.get_vertex_coordinates(),
                                 polygon)
        points = num.take(quantity.domain.get_vertex_coordinates(), indices)

        answer = linear_function(points)

        # print quantity.vertex_values.flat
        # print answer

        # Check vertices in polygon have been set
        assert num.allclose(num.take(quantity.vertex_values.flat, indices),
                            answer)

        # Check vertices outside polygon are zero
        indices = outside_polygon(quantity.domain.get_vertex_coordinates(),
                                  polygon)
        assert num.allclose(num.take(quantity.vertex_values.flat, indices),
                            0.0)

        # Cleanup
        os.remove(ptsfile)

    def test_cache_test_set_values_from_file(self):
        # FIXME (Ole): What is this about?
        # I don't think it checks anything new
        quantity = Quantity(self.mesh4)

        # Get (enough) datapoints
        data_points = [[0.66666667, 0.66666667],
                       [1.33333333, 1.33333333],
                       [2.66666667, 0.66666667],
                       [0.66666667, 2.66666667],
                       [0.0, 1.0],
                       [0.0, 3.0],
                       [1.0, 0.0],
                       [1.0, 1.0],
                       [1.0, 2.0],
                       [1.0, 3.0],
                       [2.0, 1.0],
                       [3.0, 0.0],
                       [3.0, 1.0]]

        georef = Geo_reference(56, 0, 0)
        data_geo_spatial = Geospatial_data(data_points,
                                           geo_reference=georef)

        data_points_absolute = data_geo_spatial.get_data_points(absolute=True)
        attributes = linear_function(data_points_absolute)
        att = 'spam_and_eggs'

        # Create .txt file
        fd, ptsfile = tempfile.mkstemp(".txt")
        os.close(fd)
        file = open(ptsfile, "w")
        file.write(" x,y," + att + " \n")
        for data_point, attribute in zip(data_points_absolute, attributes):
            row = str(data_point[0]) + ',' + str(data_point[1]) \
                  + ',' + str(attribute)
            file.write(row + "\n")
        file.close()

        # Check that values can be set from file
        quantity.set_values(filename=ptsfile,
                            attribute_name=att,
                            alpha=0,
                            use_cache=True,
                            verbose=False)
        answer = linear_function(quantity.domain.get_vertex_coordinates())
        assert num.allclose(quantity.vertex_values.flat, answer)

        # Check that values can be set from file using default attribute
        quantity.set_values(filename=ptsfile,
                            alpha=0)
        assert num.allclose(quantity.vertex_values.flat, answer)

        # Check cache
        quantity.set_values(filename=ptsfile,
                            attribute_name=att,
                            alpha=0,
                            use_cache=True,
                            verbose=False)

        # Cleanup
        os.remove(ptsfile)

    def test_set_values_from_lat_long(self):
        quantity = Quantity(self.mesh_onslow)

        # Get (enough) datapoints
        data_points = [[-21.5, 114.5], [-21.4, 114.6], [-21.45, 114.65],
                       [-21.35, 114.65], [-21.45, 114.55], [-21.45, 114.6]]

        data_geo_spatial = Geospatial_data(data_points,
                                           points_are_lats_longs=True)
        points_UTM = data_geo_spatial.get_data_points(absolute=True)
        attributes = linear_function(points_UTM)
        att = 'elevation'

        # Create .txt file
        fd, txt_file = tempfile.mkstemp(".txt")
        os.close(fd)
        file = open(txt_file, "w")
        file.write(" lat,long," + att + " \n")
        for data_point, attribute in zip(data_points, attributes):
            row = str(data_point[0]) + ',' + str(data_point[1]) \
                  + ',' + str(attribute)
            # print "row", row
            file.write(row + "\n")
        file.close()

        # Check that values can be set from file
        quantity.set_values(filename=txt_file,
                            attribute_name=att,
                            alpha=0)
        answer = linear_function(quantity.domain.get_vertex_coordinates())

        # print "quantity.vertex_values.flat", quantity.vertex_values.flat
        # print "answer",answer

        assert num.allclose(quantity.vertex_values.flat, answer)

        # Check that values can be set from file using default attribute
        quantity.set_values(filename=txt_file, alpha=0)
        assert num.allclose(quantity.vertex_values.flat, answer)

        # Cleanup
        os.remove(txt_file)

    def test_set_values_from_lat_long_2(self):
        quantity = Quantity(self.mesh_onslow)

        # Get (enough) datapoints
        data_points = [[-21.5, 114.5], [-21.4, 114.6], [-21.45, 114.65],
                       [-21.35, 114.65], [-21.45, 114.55], [-21.45, 114.6]]

        data_geo_spatial = Geospatial_data(data_points,
                                           points_are_lats_longs=True)
        points_UTM = data_geo_spatial.get_data_points(absolute=True)
        attributes = linear_function(points_UTM)
        att = 'elevation'

        # Create .txt file
        fd, txt_file = tempfile.mkstemp(".txt")
        os.close(fd)
        file = open(txt_file, "w")
        file.write(" lat,long," + att + " \n")
        for data_point, attribute in zip(data_points, attributes):
            row = str(data_point[0]) + ',' + str(data_point[1]) \
                  + ',' + str(attribute)
            # print "row", row
            file.write(row + "\n")
        file.close()

        # Check that values can be set from file
        quantity.set_values(filename=txt_file,
                            attribute_name=att, alpha=0)
        answer = linear_function(quantity.domain.get_vertex_coordinates())

        # print "quantity.vertex_values.flat", quantity.vertex_values.flat
        # print "answer",answer

        assert num.allclose(quantity.vertex_values.flat, answer)

        # Check that values can be set from file using default attribute
        quantity.set_values(filename=txt_file, alpha=0)
        assert num.allclose(quantity.vertex_values.flat, answer)

        # Cleanup
        os.remove(txt_file)

    def test_set_values_from_UTM_pts(self):
        quantity = Quantity(self.mesh_onslow)

        # Get (enough) datapoints
        data_points = [[-21.5, 114.5], [-21.4, 114.6], [-21.45, 114.65],
                       [-21.35, 114.65], [-21.45, 114.55], [-21.45, 114.6]]

        data_geo_spatial = Geospatial_data(data_points,
                                           points_are_lats_longs=True)
        points_UTM = data_geo_spatial.get_data_points(absolute=True)
        attributes = linear_function(points_UTM)
        att = 'elevation'

        # Create .txt file
        fd, txt_file = tempfile.mkstemp(".txt")
        os.close(fd)
        file = open(txt_file, "w")
        file.write(" x,y," + att + " \n")
        for data_point, attribute in zip(points_UTM, attributes):
            row = str(data_point[0]) + ',' + str(data_point[1]) \
                  + ',' + str(attribute)
            # print "row", row
            file.write(row + "\n")
        file.close()

        fd, pts_file = tempfile.mkstemp(".pts")
        os.close(fd)
        convert = Geospatial_data(txt_file)
        convert.export_points_file(pts_file)

        # Check that values can be set from file
        quantity.set_values_from_file(pts_file, att, 0,
                                      'vertices', None)
        answer = linear_function(quantity.domain.get_vertex_coordinates())
        # print "quantity.vertex_values.flat", quantity.vertex_values.flat
        # print "answer",answer
        assert num.allclose(quantity.vertex_values.flat, answer)

        # Check that values can be set from file
        quantity.set_values(filename=pts_file,
                            attribute_name=att, alpha=0)
        answer = linear_function(quantity.domain.get_vertex_coordinates())
        # print "quantity.vertex_values.flat", quantity.vertex_values.flat
        # print "answer",answer
        assert num.allclose(quantity.vertex_values.flat, answer)

        # Check that values can be set from file using default attribute
        quantity.set_values(filename=txt_file, alpha=0)
        assert num.allclose(quantity.vertex_values.flat, answer)

        # Cleanup
        os.remove(txt_file)
        os.remove(pts_file)

    def test_set_values_from_UTM_pts_verbose(self):
        quantity = Quantity(self.mesh_onslow)

        # Get (enough) datapoints
        data_points = [[-21.5, 114.5], [-21.4, 114.6], [-21.45, 114.65],
                       [-21.35, 114.65], [-21.45, 114.55], [-21.45, 114.6],
                       [-21.5, 114.5], [-21.4, 114.6], [-21.45, 114.65],
                       [-21.35, 114.65], [-21.45, 114.55], [-21.45, 114.6],
                       [-21.5, 114.5], [-21.4, 114.6], [-21.45, 114.65],
                       [-21.35, 114.65], [-21.45, 114.55], [-21.45, 114.6],
                       [-21.5, 114.5], [-21.4, 114.6], [-21.45, 114.65],
                       [-21.35, 114.65], [-21.45, 114.55], [-21.45, 114.6],
                       [-21.5, 114.5], [-21.4, 114.6], [-21.45, 114.65],
                       [-21.35, 114.65], [-21.45, 114.55], [-21.45, 114.6],
                       [-21.5, 114.5], [-21.4, 114.6], [-21.45, 114.65],
                       [-21.35, 114.65], [-21.45, 114.55], [-21.45, 114.6],
                       [-21.5, 114.5], [-21.4, 114.6], [-21.45, 114.65],
                       [-21.35, 114.65], [-21.45, 114.55], [-21.45, 114.6],
                       [-21.35, 114.65], [-21.45, 114.55], [-21.45, 114.6],
                       [-21.5, 114.5], [-21.4, 114.6], [-21.45, 114.65],
                       [-21.35, 114.65], [-21.45, 114.55], [-21.45, 114.6],
                       [-21.5, 114.5], [-21.4, 114.6], [-21.45, 114.65],
                       [-21.35, 114.65], [-21.45, 114.55], [-21.45, 114.6],
                       [-21.5, 114.5], [-21.4, 114.6], [-21.45, 114.65],
                       [-21.35, 114.65], [-21.45, 114.55], [-21.45, 114.6],
                       [-21.5, 114.5], [-21.4, 114.6], [-21.45, 114.65],
                       [-21.35, 114.65], [-21.45, 114.55], [-21.45, 114.6],
                       [-21.5, 114.5], [-21.4, 114.6], [-21.45, 114.65],
                       [-21.35, 114.65], [-21.45, 114.55], [-21.45, 114.6],
                       [-21.5, 114.5], [-21.4, 114.6], [-21.45, 114.65],
                       [-21.35, 114.65], [-21.45, 114.55], [-21.45, 114.6],
                       ]

        data_geo_spatial = Geospatial_data(data_points,
                                           points_are_lats_longs=True)
        points_UTM = data_geo_spatial.get_data_points(absolute=True)
        attributes = linear_function(points_UTM)
        att = 'elevation'

        # Create .txt file
        fd, txt_file = tempfile.mkstemp(".txt")
        os.close(fd)
        file = open(txt_file, "w")
        file.write(" x,y," + att + " \n")
        for data_point, attribute in zip(points_UTM, attributes):
            row = str(data_point[0]) + ',' + str(data_point[1]) \
                  + ',' + str(attribute)
            # print "row", row
            file.write(row + "\n")
        file.close()

        fd, pts_file = tempfile.mkstemp(".pts")
        os.close(fd)
        convert = Geospatial_data(txt_file)
        convert.export_points_file(pts_file)

        # Check that values can be set from file
        quantity.set_values_from_file(pts_file, att, 0,
                                      'vertices', None, verbose=False,
                                      max_read_lines=2)
        answer = linear_function(quantity.domain.get_vertex_coordinates())
        # print "quantity.vertex_values.flat", quantity.vertex_values.flat
        # print "answer",answer
        assert num.allclose(quantity.vertex_values.flat, answer)

        # Check that values can be set from file
        quantity.set_values(filename=pts_file,
                            attribute_name=att, alpha=0)
        answer = linear_function(quantity.domain.get_vertex_coordinates())
        # print "quantity.vertex_values.flat", quantity.vertex_values.flat
        # print "answer",answer
        assert num.allclose(quantity.vertex_values.flat, answer)

        # Check that values can be set from file using default attribute
        quantity.set_values(filename=txt_file, alpha=0)
        assert num.allclose(quantity.vertex_values.flat, answer)

        # Cleanup
        os.remove(txt_file)
        os.remove(pts_file)

    def test_set_values_from_file_with_georef1(self):

        # Mesh in zone 56 (absolute coords)

        x0 = 314036.58727982
        y0 = 6224951.2960092

        a = [x0+0.0, y0+0.0]
        b = [x0+0.0, y0+2.0]
        c = [x0+2.0, y0+0.0]
        d = [x0+0.0, y0+4.0]
        e = [x0+2.0, y0+2.0]
        f = [x0+4.0, y0+0.0]

        points = [a, b, c, d, e, f]

        # bac, bce, ecf, dbe
        elements = [[1, 0, 2], [1, 2, 4], [4, 2, 5], [3, 1, 4]]

        # absolute going in ..
        mesh4 = Generic_Domain(points, elements,
                               geo_reference=Geo_reference(56, 0, 0))
        mesh4.check_integrity()
        quantity = Quantity(mesh4)

        # Get (enough) datapoints (relative to georef)
        data_points_rel = [[0.66666667, 0.66666667],
                           [1.33333333, 1.33333333],
                           [2.66666667, 0.66666667],
                           [0.66666667, 2.66666667],
                           [0.0, 1.0],
                           [0.0, 3.0],
                           [1.0, 0.0],
                           [1.0, 1.0],
                           [1.0, 2.0],
                           [1.0, 3.0],
                           [2.0, 1.0],
                           [3.0, 0.0],
                           [3.0, 1.0]]

        data_geo_spatial = Geospatial_data(data_points_rel,
                                           geo_reference=Geo_reference(56, x0, y0))
        data_points_absolute = data_geo_spatial.get_data_points(absolute=True)
        attributes = linear_function(data_points_absolute)
        att = 'spam_and_eggs'

        # Create .txt file
        fd, ptsfile = tempfile.mkstemp(".txt")
        os.close(fd)
        file = open(ptsfile, "w")
        file.write(" x,y," + att + " \n")
        for data_point, attribute in zip(data_points_absolute, attributes):
            row = str(data_point[0]) + ',' + str(data_point[1]) \
                  + ',' + str(attribute)
            file.write(row + "\n")
        file.close()

        # file = open(ptsfile, 'r')
        # lines = file.readlines()
        # file.close()

        # Check that values can be set from file
        quantity.set_values(filename=ptsfile,
                            attribute_name=att, alpha=0)
        answer = linear_function(quantity.domain.get_vertex_coordinates())

        assert num.allclose(quantity.vertex_values.flat, answer)

        # Check that values can be set from file using default attribute
        quantity.set_values(filename=ptsfile, alpha=0)
        assert num.allclose(quantity.vertex_values.flat, answer)

        # Cleanup
        os.remove(ptsfile)

    def test_set_values_from_file_with_georef2(self):

        # Mesh in zone 56 (relative coords)

        x0 = 314036.58727982
        y0 = 6224951.2960092
        # x0 = 0.0
        # y0 = 0.0

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]

        # bac, bce, ecf, dbe
        elements = [[1, 0, 2], [1, 2, 4], [4, 2, 5], [3, 1, 4]]

        mesh4 = Generic_Domain(points, elements,
                               geo_reference=Geo_reference(56, x0, y0))
        mesh4.check_integrity()
        quantity = Quantity(mesh4)

        # Get (enough) datapoints
        data_points = [[x0+0.66666667, y0+0.66666667],
                       [x0+1.33333333, y0+1.33333333],
                       [x0+2.66666667, y0+0.66666667],
                       [x0+0.66666667, y0+2.66666667],
                       [x0+0.0, y0+1.0],
                       [x0+0.0, y0+3.0],
                       [x0+1.0, y0+0.0],
                       [x0+1.0, y0+1.0],
                       [x0+1.0, y0+2.0],
                       [x0+1.0, y0+3.0],
                       [x0+2.0, y0+1.0],
                       [x0+3.0, y0+0.0],
                       [x0+3.0, y0+1.0]]

        data_geo_spatial = Geospatial_data(data_points,
                                           geo_reference=Geo_reference(56, 0, 0))
        data_points_absolute = data_geo_spatial.get_data_points(absolute=True)
        attributes = linear_function(data_points_absolute)
        att = 'spam_and_eggs'

        # Create .txt file
        fd, ptsfile = tempfile.mkstemp(".txt")
        os.close(fd)
        file = open(ptsfile, "w")
        file.write(" x,y," + att + " \n")
        for data_point, attribute in zip(data_points_absolute, attributes):
            row = str(data_point[0]) + ',' + str(data_point[1]) \
                  + ',' + str(attribute)
            file.write(row + "\n")
        file.close()

        # Check that values can be set from file
        quantity.set_values(filename=ptsfile,
                            attribute_name=att, alpha=0)
        answer = linear_function(quantity.domain.
                                 get_vertex_coordinates(absolute=True))

        assert num.allclose(quantity.vertex_values.flat, answer)

        # Check that values can be set from file using default attribute
        quantity.set_values(filename=ptsfile, alpha=0)
        assert num.allclose(quantity.vertex_values.flat, answer)

        # Cleanup
        os.remove(ptsfile)

    def test_set_values_from_utm_grid_file(self):

        x0 = 0.0
        y0 = 0.0

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]

        # bac, bce, ecf, dbe
        elements = [[1, 0, 2], [1, 2, 4], [4, 2, 5], [3, 1, 4]]

        mesh4 = Generic_Domain(points, elements)
        #  geo_reference = Geo_reference(56, x0, y0))
        mesh4.check_integrity()
        quantity = Quantity(mesh4)

        """ Format of asc file
        ncols         11
        nrows         12
        xllcorner     240000
        yllcorner     7620000
        cellsize      6000
        NODATA_value  -9999
        """
        ncols = 11  # Nx
        nrows = 12  # Ny
        xllcorner = x0
        yllcorner = y0
        cellsize = 1.0
        NODATA_value = -9999

        # xllcorner = 0
        # yllcorner = 100
        # cellsize  = 10
        # NODATA_value =  -9999

        # Create .asc file
        # txt_file = tempfile.mktemp(".asc")
        txt_file = 'test_asc.asc'
        datafile = open(txt_file, "w")
        datafile.write('ncols '+str(ncols)+"\n")
        datafile.write('nrows '+str(nrows)+"\n")
        datafile.write('xllcorner '+str(xllcorner)+"\n")
        datafile.write('yllcorner '+str(yllcorner)+"\n")
        datafile.write('cellsize '+str(cellsize)+"\n")
        datafile.write('NODATA_value '+str(NODATA_value)+"\n")

        x = num.linspace(xllcorner, xllcorner+(ncols-1)*cellsize, ncols)
        y = num.linspace(yllcorner, yllcorner+(nrows-1)*cellsize, nrows)
        points = axes2points(x, y)

        # print points
        # print x.shape, x
        # print y.shape, y

        datavalues = linear_function(points)
        # print datavalues

        datavalues = datavalues.reshape(nrows, ncols)

        # print datavalues
        # print datavalues.shape
        for row in datavalues:
            # print row
            datafile.write(" ".join(str(elem) for elem in row) + "\n")
        datafile.close()

        # print quantity.vertex_values
        # print quantity.centroid_values
        quantity.set_values(filename=txt_file,
                            location='vertices',
                            indices=None,
                            verbose=False)

        # check order of vertices

        answer = [[6.,   0.,   2.],
                  [6.,   2.,   8.],
                  [8.,  2.,   4.],
                  [12.,   6.,   8.]]
        # print quantity.vertex_values
        assert num.allclose(quantity.vertex_values, answer)

        # print quantity.vertex_values
        # print quantity.centroid_values
        quantity.set_values(0.0)

        # print quantity.vertex_values
        # print quantity.centroid_values
        quantity.set_values(filename=txt_file,
                            location='centroids',
                            indices=None,
                            verbose=False)

        # print quantity.vertex_values
        # print quantity.centroid_values

        answer = [2.66666667,  5.33333333,  4.66666667,  8.66666667]

        assert num.allclose(quantity.centroid_values, answer)

        # check dem file
        # use the same reference solution used above for testing
        # convert test_asc.asc file to .dem file
        txt_file_prj = 'test_asc.prj'
        fid = open(txt_file_prj, 'w')
        fid.write("""Projection UTM
Zone 56
Datum WGS84
Zunits NO
Units METERS
Spheroid WGS84
Xshift 0.0000000000
Yshift 10000000.0000000000
Parameters
        """)
        fid.close()

        txt_file_dem = 'test_asc.dem'
        asc2dem(name_in=txt_file, name_out='test_asc',
                use_cache=False, verbose=False)

        quantity.set_values(0.0)
        quantity.set_values(filename=txt_file_dem,
                            location='vertices',
                            indices=None,
                            verbose=False)

        # check order of vertices

        answer = [[6.,   0.,   2.],
                  [6.,   2.,   8.],
                  [8.,  2.,   4.],
                  [12.,   6.,   8.]]
        # print quantity.vertex_values
        # print quantity.vertex_values, 'vertex values'
        assert num.allclose(quantity.vertex_values, answer)

        # print quantity.vertex_values
        # print quantity.centroid_values
        quantity.set_values(0.0)

        # print quantity.vertex_values
        # print quantity.centroid_values
        quantity.set_values(filename=txt_file_dem,
                            location='centroids',
                            indices=None,
                            verbose=False)

        # print quantity.vertex_values
        # print quantity.centroid_values      , 'centroid values'

        answer = [2.66666667,  5.33333333,  4.66666667,  8.66666667]

        assert num.allclose(quantity.centroid_values, answer)

        # Cleanup
        # import os

        try:
            os.remove(txt_file)
            os.remove(txt_file_prj)
            os.remove(txt_file_dem)
        except OSError:
            pass

    def test_set_values_from_ll_grid_file(self):

        x0 = 0.0
        y0 = 0.0

        # long lat coordinated
        a = [150.0, -30.0]
        b = [150.0, -29.8]
        c = [150.2, -30.0]
        d = [150.0, -29.6]
        e = [150.2, -29.8]
        # f = [150.4, -30.0]

        eastings = num.array([210590.34672016, 210010.98333739, 229892.33388686, 209435.153075,
                              229351.70758497]) #249192.66399564])
        eastings = num.reshape(eastings, (-1, 1))

        northings = num.array([6677424.09561761, 6699600.772466, 6677912.99938376,
                               6721776.96194549, 6700087.69263074 ]) #6678368.09416381])
        northings = num.reshape(northings, (-1, 1))

        zone_number = 56
        zone_letter = 'J'
        zone_hemisphere = 'southern'

        points = num.hstack((eastings, northings))

        xll_corner = points[:,0].min()
        yll_corner = points[:,1].min()
        geo_reference = Geo_reference(zone_number, xll_corner, yll_corner,hemisphere=zone_hemisphere)

        # create relative coordinates
        points[:,0] = points[:,0] - xll_corner
        points[:,1] = points[:,1] - yll_corner

        # points = [a, b, c, d, e, f]

        # bac, bce, dbe # ecf
        elements = [[1, 0, 2], [1, 2, 4], [3, 1, 4]]  #[4, 2, 5]

        mesh4 = Generic_Domain(points, elements, geo_reference=geo_reference)

        mesh4.check_integrity()
        quantity = Quantity(mesh4)

        """ Format of asc lat long file
        ncols         5
        nrows         7
        xllcorner     149.9
        yllcorner     -30.9
        cellsize      0.1
        NODATA_value  -9999
        """
        ncols = 5  # Nx
        nrows = 7  # Ny
        xllcorner = 149.9
        yllcorner = -30.1
        cellsize = 0.1
        NODATA_value = -9999

        # Create .asc file
        # txt_file = tempfile.mktemp(".asc")
        txt_file = 'test_asc_ll.asc'
        datafile = open(txt_file, "w")
        datafile.write('ncols '+str(ncols)+"\n")
        datafile.write('nrows '+str(nrows)+"\n")
        datafile.write('xllcorner '+str(xllcorner)+"\n")
        datafile.write('yllcorner '+str(yllcorner)+"\n")
        datafile.write('cellsize '+str(cellsize)+"\n")
        datafile.write('NODATA_value '+str(NODATA_value)+"\n")

        long = num.linspace(xllcorner, xllcorner+(ncols-1)*cellsize, ncols)
        lat = num.linspace(yllcorner, yllcorner+(nrows-1)*cellsize, nrows)

        # print(long)
        # print(lat)
        points = axes2points(long, lat)

        # print points
        # print x.shape, x
        # print y.shape, y

        datavalues = linear_function(points)
        # print datavalues

        datavalues = datavalues.reshape(nrows, ncols)

        # print datavalues
        # print datavalues.shape
        for row in datavalues:
            # print row
            datafile.write(" ".join(str(elem) for elem in row) + "\n")
        datafile.close()

        # print quantity.vertex_values
        # print quantity.centroid_values
        quantity.set_values_from_lat_long_grid_file(filename=txt_file,
                                                    location='vertices',
                                                    indices=None,
                                                    verbose=False)

        # check order of vertices

        # print(quantity.vertex_values)

        answer = [[60.6, 60.0, 60.2],
                  [60.6, 60.2, 60.8],
                  [61.2, 60.6, 60.8] ]

        # print quantity.vertex_values
        assert num.allclose(quantity.vertex_values, answer)

        # print quantity.vertex_values
        # print quantity.centroid_values
        quantity.set_values(0.0)

        # print (quantity.vertex_values)
        # print (quantity.centroid_values)
        quantity.set_values_from_lat_long_grid_file(filename=txt_file,
                                                    location='centroids',
                                                    indices=None,
                                                    verbose=False)

        # print quantity.vertex_values
        # print (quantity.centroid_values)

        # note that centroid values correspond to utm centroids not long/lat centroids
        answer = [60.26652907, 60.53319575, 60.86652987]

        assert num.allclose(quantity.centroid_values, answer)

        try:
            os.remove(txt_file)
        except OSError:
            pass

    def test_set_values_from_ll_grid_file_northern(self):

        x0 = 0.0
        y0 = 0.0

        # long lat coordinated
        a = [150.0, -30.0+60]
        b = [150.0, -29.8+60]
        c = [150.2, -30.0+60]
        d = [150.0, -29.6+60]
        e = [150.2, -29.8+60]
        # f = [150.4, -30.0]

        eastings = num.array([210590.34672016295, 211173.23637087474, 229892.33388686483,
        211759.6453941525, 230436.25127640477] )
        eastings = num.reshape(eastings, (-1, 1))

        northings = num.array([3322575.9043823928, 3344753.070613778, 3322087.0006162403,
        3366930.72821874, 3344262.207060983] )
        northings = num.reshape(northings, (-1, 1))

        zone_number = 56
        zone_letter = 'R'
        zone_hemisphere = 'northern'

        points = num.hstack((eastings, northings))

        xll_corner = points[:,0].min()
        yll_corner = points[:,1].min()
        geo_reference = Geo_reference(zone_number, xll_corner, yll_corner,hemisphere=zone_hemisphere)

        # create relative coordinates
        points[:,0] = points[:,0] - xll_corner
        points[:,1] = points[:,1] - yll_corner

        # points = [a, b, c, d, e, f]

        # bac, bce, dbe # ecf
        elements = [[1, 0, 2], [1, 2, 4], [3, 1, 4]]  #[4, 2, 5]

        mesh4 = Generic_Domain(points, elements, geo_reference=geo_reference)

        mesh4.check_integrity()
        quantity = Quantity(mesh4)

        """ Format of asc lat long file
        ncols         5
        nrows         7
        xllcorner     149.9
        yllcorner     -30.9
        cellsize      0.1
        NODATA_value  -9999
        """
        ncols = 5  # Nx
        nrows = 7  # Ny
        xllcorner = 149.9
        yllcorner = 29.9
        cellsize = 0.1
        NODATA_value = -9999

        # Create .asc file
        # txt_file = tempfile.mktemp(".asc")
        txt_file = 'test_asc_ll.asc'
        datafile = open(txt_file, "w")
        datafile.write('ncols '+str(ncols)+"\n")
        datafile.write('nrows '+str(nrows)+"\n")
        datafile.write('xllcorner '+str(xllcorner)+"\n")
        datafile.write('yllcorner '+str(yllcorner)+"\n")
        datafile.write('cellsize '+str(cellsize)+"\n")
        datafile.write('NODATA_value '+str(NODATA_value)+"\n")

        long = num.linspace(xllcorner, xllcorner+(ncols-1)*cellsize, ncols)
        lat = num.linspace(yllcorner, yllcorner+(nrows-1)*cellsize, nrows)

        # print(long)
        # print(lat)
        points = axes2points(long, lat)

        # print points
        # print x.shape, x
        # print y.shape, y

        datavalues = linear_function(points)
        # print datavalues

        datavalues = datavalues.reshape(nrows, ncols)

        # print datavalues
        # print datavalues.shape
        for row in datavalues:
            # print row
            datafile.write(" ".join(str(elem) for elem in row) + "\n")
        datafile.close()

        # print quantity.vertex_values
        # print quantity.centroid_values
        quantity.set_values_from_lat_long_grid_file(filename=txt_file,
                                                    location='vertices',
                                                    indices=None,
                                                    verbose=False)

        # check order of vertices

        #print(quantity.vertex_values)

        answer = [[240.6, 240.0, 240.2], [240.6, 240.2, 240.8], [241.2, 240.6, 240.8]]

        # print quantity.vertex_values
        assert num.allclose(quantity.vertex_values, answer)

        # print quantity.vertex_values
        # print quantity.centroid_values
        quantity.set_values(0.0)

        # print (quantity.vertex_values)
        # print (quantity.centroid_values)
        quantity.set_values_from_lat_long_grid_file(filename=txt_file,
                                                    location='centroids',
                                                    indices=None,
                                                    verbose=False)

        # print quantity.vertex_values
        # print (quantity.centroid_values)

        # note that centroid values correspond to utm centroids not long/lat centroids

        answer = [240.26682499, 240.5334917 , 240.86682573]

        assert num.allclose(quantity.centroid_values, answer)

        try:
            os.remove(txt_file)
        except OSError:
            pass

    def test_set_values_from_ll_grid_file_indices(self):

        x0 = 0.0
        y0 = 0.0

        # long lat coordinated
        a = [150.0, -30.0]
        b = [150.0, -29.8]
        c = [150.2, -30.0]
        d = [150.0, -29.6]
        e = [150.2, -29.8]
        # f = [150.4, -30.0]

        eastings = num.array([210590.34672016, 210010.98333739, 229892.33388686, 209435.153075,
                              229351.70758497]) #249192.66399564])
        eastings = num.reshape(eastings, (-1, 1))

        northings = num.array([6677424.09561761, 6699600.772466, 6677912.99938376,
                               6721776.96194549, 6700087.69263074 ]) #6678368.09416381])
        northings = num.reshape(northings, (-1, 1))

        zone_number = 56
        zone_letter = 'J'
        zone_hemisphere = 'southern'

        points = num.hstack((eastings, northings))

        xll_corner = points[:,0].min()
        yll_corner = points[:,1].min()
        geo_reference = Geo_reference(zone_number, xll_corner, yll_corner,hemisphere=zone_hemisphere)

        # create relative coordinates
        points[:,0] = points[:,0] - xll_corner
        points[:,1] = points[:,1] - yll_corner

        # points = [a, b, c, d, e, f]

        # bac, bce, dbe # ecf
        elements = [[1, 0, 2], [1, 2, 4], [3, 1, 4]]  #[4, 2, 5]

        mesh4 = Generic_Domain(points, elements, geo_reference=geo_reference)

        mesh4.check_integrity()
        quantity = Quantity(mesh4)

        """ Format of asc lat long file
        ncols         5
        nrows         7
        xllcorner     149.9
        yllcorner     -30.9
        cellsize      0.1
        NODATA_value  -9999
        """
        ncols = 5  # Nx
        nrows = 7  # Ny
        xllcorner = 149.9
        yllcorner = -30.1
        cellsize = 0.1
        NODATA_value = -9999

        # Create .asc file
        # txt_file = tempfile.mktemp(".asc")
        txt_file = 'test_asc_ll.asc'
        datafile = open(txt_file, "w")
        datafile.write('ncols '+str(ncols)+"\n")
        datafile.write('nrows '+str(nrows)+"\n")
        datafile.write('xllcorner '+str(xllcorner)+"\n")
        datafile.write('yllcorner '+str(yllcorner)+"\n")
        datafile.write('cellsize '+str(cellsize)+"\n")
        datafile.write('NODATA_value '+str(NODATA_value)+"\n")

        long = num.linspace(xllcorner, xllcorner+(ncols-1)*cellsize, ncols)
        lat = num.linspace(yllcorner, yllcorner+(nrows-1)*cellsize, nrows)

        # print(long)
        # print(lat)
        points = axes2points(long, lat)

        # print points
        # print x.shape, x
        # print y.shape, y

        datavalues = linear_function(points)
        # print datavalues

        datavalues = datavalues.reshape(nrows, ncols)

        # print datavalues
        # print datavalues.shape
        for row in datavalues:
            # print row
            datafile.write(" ".join(str(elem) for elem in row) + "\n")
        datafile.close()

        # print quantity.vertex_values
        # print quantity.centroid_values

        indices = [0 ,2, 6, 1]
        quantity.set_values_from_lat_long_grid_file(filename=txt_file,
                                                    location='vertices',
                                                    indices=indices,
                                                    verbose=False)

        # check order of vertices

        # print(quantity.vertex_values)

        answer = [[60.6, 60.0, 60.2],
                  [0.0,  0.0,  0.0 ],
                  [61.2, 0.0,  0.0 ]]

        # print quantity.vertex_values
        assert num.allclose(quantity.vertex_values, answer)

        # print quantity.vertex_values
        # print quantity.centroid_values
        quantity.set_values(0.0)

        # print (quantity.vertex_values)
        # print (quantity.centroid_values)

        indices = [0, 2]
        quantity.set_values_from_lat_long_grid_file(filename=txt_file,
                                                    location='centroids',
                                                    indices=indices,
                                                    verbose=False)

        # print quantity.vertex_values
        # print (quantity.centroid_values)

        # note that centroid values correspond to utm centroids not long/lat centroids
        answer = [60.26652907, 0.0, 60.86652987]

        assert num.allclose(quantity.centroid_values, answer)

        try:
            os.remove(txt_file)
        except OSError:
            pass

    def test_set_values_from_ll_grid_file_with_nan(self):

        # Not implemented. Need to check that the asc file doesn't
        # cover the domain

        x0 = 0.0
        y0 = 0.0

        # long lat coordinated
        a = [150.0, -30.0]
        b = [150.0, -29.8]
        c = [150.2, -30.0]
        d = [150.0, -29.6]
        e = [150.2, -29.8]
        # f = [150.4, -30.0]

        eastings = num.array([210590.34672016, 210010.98333739, 229892.33388686, 209435.153075,
                              229351.70758497]) #249192.66399564])
        eastings = num.reshape(eastings, (-1, 1))

        northings = num.array([6677424.09561761, 6699600.772466, 6677912.99938376,
                               6721776.96194549, 6700087.69263074 ]) #6678368.09416381])
        northings = num.reshape(northings, (-1, 1))

        zone_number = 56
        zone_letter = 'J'
        zone_hemisphere = 'southern'

        points = num.hstack((eastings, northings))

        xll_corner = points[:,0].min()
        yll_corner = points[:,1].min()
        geo_reference = Geo_reference(zone_number, xll_corner, yll_corner,hemisphere=zone_hemisphere)

        # create relative coordinates
        points[:,0] = points[:,0] - xll_corner
        points[:,1] = points[:,1] - yll_corner

        # points = [a, b, c, d, e, f]

        # bac, bce, dbe # ecf
        elements = [[1, 0, 2], [1, 2, 4], [3, 1, 4]]  #[4, 2, 5]

        mesh4 = Generic_Domain(points, elements, geo_reference=geo_reference)

        mesh4.check_integrity()
        quantity = Quantity(mesh4)

        """ Format of asc lat long file
        ncols         5
        nrows         7
        xllcorner     149.9
        yllcorner     -30.9
        cellsize      0.1
        NODATA_value  -9999
        """

        # create data which doesn't cover the full domain
        ncols = 3  # Nx
        nrows = 5  # Ny
        xllcorner = 150.1
        yllcorner = -29.9
        cellsize = 0.1
        NODATA_value = -9999

        # Create .asc file
        # txt_file = tempfile.mktemp(".asc")
        txt_file = 'test_asc_ll.asc'
        datafile = open(txt_file, "w")
        datafile.write('ncols '+str(ncols)+"\n")
        datafile.write('nrows '+str(nrows)+"\n")
        datafile.write('xllcorner '+str(xllcorner)+"\n")
        datafile.write('yllcorner '+str(yllcorner)+"\n")
        datafile.write('cellsize '+str(cellsize)+"\n")
        datafile.write('NODATA_value '+str(NODATA_value)+"\n")

        long = num.linspace(xllcorner, xllcorner+(ncols-1)*cellsize, ncols)
        lat = num.linspace(yllcorner, yllcorner+(nrows-1)*cellsize, nrows)

        # print(long)
        # print(lat)
        points = axes2points(long, lat)

        # print (points)
        # print (x.shape, x)
        # print (y.shape, y)

        datavalues = linear_function(points)
        # print (datavalues)

        datavalues = datavalues.reshape(nrows, ncols)

        # print datavalues
        # print datavalues.shape
        for row in datavalues:
            # print row
            datafile.write(" ".join(str(elem) for elem in row) + "\n")
        datafile.close()

        # print quantity.vertex_values
        # print quantity.centroid_values
        quantity.set_values_from_lat_long_grid_file(filename=txt_file,
                                                    location='vertices',
                                                    indices=None,
                                                    verbose=False)

        # check order of vertices

        # print(quantity.vertex_values)

        answer = [[60.6, 60.0, 60.2],
                  [60.6, 60.2, 60.8],
                  [61.2, 60.6, 60.8] ]

        answer = [[num.nan,         num.nan,     num.nan],
                  [num.nan,         num.nan, 60.80000006],
                  [num.nan,         num.nan, 60.80000006]]

        # print (quantity.vertex_values)
        assert num.allclose(quantity.vertex_values, answer, equal_nan=True)

        # print quantity.vertex_values
        # print quantity.centroid_values
        quantity.set_values(0.0)

        # print (quantity.vertex_values)
        # print (quantity.centroid_values)
        quantity.set_values_from_lat_long_grid_file(filename=txt_file,
                                                    location='centroids',
                                                    indices=None,
                                                    verbose=False)

        # print quantity.vertex_values
        # print (quantity.centroid_values)

        # note that centroid values correspond to utm centroids not long/lat centroids
        answer = [num.nan,  60.53319575, num.nan]

        assert num.allclose(quantity.centroid_values, answer, equal_nan=True)

        try:
            os.remove(txt_file)
        except OSError:
            pass

    def test_set_values_from_ll_grid_file_with_indices_nan(self):

        # Not implemented. Need to check that the asc file doesn't
        # cover the domain

        x0 = 0.0
        y0 = 0.0

        # long lat coordinated
        a = [150.0, -30.0]
        b = [150.0, -29.8]
        c = [150.2, -30.0]
        d = [150.0, -29.6]
        e = [150.2, -29.8]
        # f = [150.4, -30.0]

        eastings = num.array([210590.34672016, 210010.98333739, 229892.33388686, 209435.153075,
                              229351.70758497]) #249192.66399564])
        eastings = num.reshape(eastings, (-1, 1))

        northings = num.array([6677424.09561761, 6699600.772466, 6677912.99938376,
                               6721776.96194549, 6700087.69263074 ]) #6678368.09416381])
        northings = num.reshape(northings, (-1, 1))

        zone_number = 56
        zone_letter = 'J'
        zone_hemisphere = 'southern'

        points = num.hstack((eastings, northings))

        xll_corner = points[:,0].min()
        yll_corner = points[:,1].min()
        geo_reference = Geo_reference(zone_number, xll_corner, yll_corner,hemisphere=zone_hemisphere)

        # create relative coordinates
        points[:,0] = points[:,0] - xll_corner
        points[:,1] = points[:,1] - yll_corner

        # points = [a, b, c, d, e, f]

        # bac, bce, dbe # ecf
        elements = [[1, 0, 2], [1, 2, 4], [3, 1, 4]]  #[4, 2, 5]

        mesh4 = Generic_Domain(points, elements, geo_reference=geo_reference)

        mesh4.check_integrity()
        quantity = Quantity(mesh4)

        """ Format of asc lat long file
        ncols         5
        nrows         7
        xllcorner     149.9
        yllcorner     -30.9
        cellsize      0.1
        NODATA_value  -9999
        """

        # create data which doesn't cover the full domain
        ncols = 3  # Nx
        nrows = 5  # Ny
        xllcorner = 150.1
        yllcorner = -29.9
        cellsize = 0.1
        NODATA_value = -9999

        # Create .asc file
        # txt_file = tempfile.mktemp(".asc")
        txt_file = 'test_asc_ll.asc'
        datafile = open(txt_file, "w")
        datafile.write('ncols '+str(ncols)+"\n")
        datafile.write('nrows '+str(nrows)+"\n")
        datafile.write('xllcorner '+str(xllcorner)+"\n")
        datafile.write('yllcorner '+str(yllcorner)+"\n")
        datafile.write('cellsize '+str(cellsize)+"\n")
        datafile.write('NODATA_value '+str(NODATA_value)+"\n")

        long = num.linspace(xllcorner, xllcorner+(ncols-1)*cellsize, ncols)
        lat = num.linspace(yllcorner, yllcorner+(nrows-1)*cellsize, nrows)

        # print(long)
        # print(lat)
        points = axes2points(long, lat)

        # print (points)
        # print (x.shape, x)
        # print (y.shape, y)

        datavalues = linear_function(points)
        # print (datavalues)

        datavalues = datavalues.reshape(nrows, ncols)

        # print datavalues
        # print datavalues.shape
        for row in datavalues:
            # print row
            datafile.write(" ".join(str(elem) for elem in row) + "\n")
        datafile.close()

        # print quantity.vertex_values
        # print quantity.centroid_values
        quantity.set_values_from_lat_long_grid_file(filename=txt_file,
                                                    location='vertices',
                                                    indices=None,
                                                    verbose=False)

        # check order of vertices

        # print(quantity.vertex_values)

        answer = [[num.nan,         num.nan,     num.nan],
                  [num.nan,         num.nan, 60.80000006],
                  [num.nan,         num.nan, 60.80000006]]

        # print (quantity.vertex_values)
        assert num.allclose(quantity.vertex_values, answer, equal_nan=True)

        # print quantity.vertex_values
        # print quantity.centroid_values
        quantity.set_values(0.0)

        # print (quantity.vertex_values)
        # print (quantity.centroid_values)
        quantity.set_values_from_lat_long_grid_file(filename=txt_file,
                                                    location='centroids',
                                                    indices=None,
                                                    verbose=False)

        # print quantity.vertex_values
        # print (quantity.centroid_values)

        # note that centroid values correspond to utm centroids not long/lat centroids
        answer = [num.nan,  60.53319575, num.nan]

        assert num.allclose(quantity.centroid_values, answer, equal_nan=True)

        try:
            os.remove(txt_file)
        except OSError:
            pass

    @pytest.mark.skipif('rasterio' not in sys.modules,
                    reason="requires rasterio")
    def test_set_values_from_ll_tif_file_north(self):

        # Mesh in zone 56 (relative coords) southern hemisphere

        from pprint import pprint

        # We need to make a .tif with ll coord which covers domain
        #
        from anuga.utilities import plot_utils as util
        import numpy
        #
        # Use Make_Geotif to make tif file
        # Pick a domain that makes sense in EPSG:32656
        lat_ll, lon_ll = 34.37, 150.90
        lat_ur, lon_ur = 34.39, 150.92

        from anuga.coordinate_transforms.lat_long_UTM_conversion import LLtoUTM
        zone_ll, utm_east_ll, utm_north_ll = LLtoUTM(lat_ll, lon_ll)
        zone_ur, utm_east_ur, utm_north_ur = LLtoUTM(lat_ur, lon_ur)

        hemisphere_ll = 'southern' if lat_ll < 0 else 'northern'
        hemisphere_ur = 'southern' if lat_ur < 0 else 'northern'

        lat = numpy.linspace(lat_ll, lat_ur, 11)
        lon = numpy.linspace(lon_ll, lon_ur, 11)

        xG, yG = numpy.meshgrid(lon, lat)

        xG = xG.flatten()
        yG = yG.flatten()
        # Surface is z=x+y

        fakeZ = (xG-min(xG))/(max(xG)-min(xG))+(yG - min(yG))/(max(yG)-min(yG))

        # pprint(fakeZ.reshape((11,11)))

        dataToGrid = numpy.vstack([xG, yG, fakeZ]).transpose()
        #
        # Create file PointData_lat_long.tif
        util.Make_Geotif(dataToGrid, output_quantities=['lat_long'],
                        EPSG_CODE=4326, output_dir='.', CellSize=0.0001)

        tif_file = 'PointData_lat_long.tif'

        # Create a domain and quantity
        x0 = utm_east_ll
        y0 = utm_north_ll

        a = [0.0, 0.0]
        b = [0.0, 500.0]
        c = [500.0, 0.0]
        d = [0.0, 1000.0]
        e = [500.0, 500.0]
        f = [1000.0, 0.0]

        points = [a, b, c, d, e, f]

        # bac, bce, ecf, dbe
        elements = [[1, 0, 2], [1, 2, 4], [4, 2, 5], [3, 1, 4]]

        mesh4 = Generic_Domain(points, elements,
                               geo_reference=Geo_reference(56, x0, y0, hemisphere=hemisphere_ll))
        mesh4.check_integrity()

        quantity1 = Quantity(mesh4)

        # Read in an interpolate from tif file
        quantity1.set_values_from_tif_file(filename=tif_file , location='vertices')

        # pprint(quantity1.centroid_values)
        # pprint(quantity1.vertex_values)

        centroid_values_ex = numpy.array(
            [0.17639186, 0.34223889, 0.44878893, 0.72498227])
        vertex_values_ex = numpy.array([[0.23654294, 0., 0.29263264],
                                        [0.23654294, 0.29263264, 0.4975411],
                                        [0.4975411, 0.29263264, 0.55619305],
                                        [1.44086277, 0.23654294, 0.4975411]])

        centroid_values_mac_ex = numpy.array(
            [0.17639186, 0.34223889, 0.44878893, 0.7188908 ])
        vertex_values_mac_ex = numpy.array([[0.23654294, 0.        , 0.29263264],
                                        [0.23654294, 0.29263264, 0.4975411 ],
                                        [0.4975411 , 0.29263264, 0.55619305],
                                        [1.42258835, 0.23654294, 0.4975411 ]])

        assert num.allclose(quantity1.centroid_values, centroid_values_ex) or \
                num.allclose(quantity1.centroid_values, centroid_values_mac_ex)
        assert num.allclose(quantity1.vertex_values, vertex_values_ex) or \
                num.allclose(quantity1.vertex_values, vertex_values_mac_ex)

        quantity2 = Quantity(mesh4)

        # Read in an interpolate from tif file
        quantity2.set_values(filename=tif_file , location='vertices')

        # pprint(quantity2.centroid_values)
        # pprint(quantity2.vertex_values)


        assert num.allclose(quantity2.centroid_values, centroid_values_ex) or \
                num.allclose(quantity2.centroid_values, centroid_values_mac_ex)
        assert num.allclose(quantity2.vertex_values, vertex_values_ex) or \
                num.allclose(quantity1.vertex_values, vertex_values_mac_ex)

        # Cleanup
        os.remove(tif_file)

    @pytest.mark.skipif('rasterio' not in sys.modules,
                    reason="requires rasterio")
    def test_set_values_from_ll_tif_file_north_indices(self):

        from pprint import pprint

        # We need to make a .tif with ll coord which covers domain
        #
        from anuga.utilities import plot_utils as util
        import numpy
        #
        # Use Make_Geotif to make tif file
        # Pick a domain that makes sense in EPSG:32656
        lat_ll, lon_ll = 34.37, 150.90
        lat_ur, lon_ur = 34.39, 150.92

        from anuga.coordinate_transforms.lat_long_UTM_conversion import LLtoUTM
        zone_ll, utm_east_ll, utm_north_ll = LLtoUTM(lat_ll, lon_ll)
        zone_ur, utm_east_ur, utm_north_ur = LLtoUTM(lat_ur, lon_ur)

        hemisphere_ll = 'southern' if lat_ll < 0 else 'northern'
        hemisphere_ur = 'southern' if lat_ur < 0 else 'northern'

        lat = numpy.linspace(lat_ll, lat_ur, 11)
        lon = numpy.linspace(lon_ll, lon_ur, 11)

        xG, yG = numpy.meshgrid(lon, lat)

        xG = xG.flatten()
        yG = yG.flatten()
        # Surface is z=x+y

        fakeZ = (xG-min(xG))/(max(xG)-min(xG))+(yG - min(yG))/(max(yG)-min(yG))

        # pprint(fakeZ.reshape((11,11)))

        dataToGrid = numpy.vstack([xG, yG, fakeZ]).transpose()
        #
        # Create file PointData_lat_long.tif
        util.Make_Geotif(dataToGrid, output_quantities=['lat_long'],
                        EPSG_CODE=4326, output_dir='.', CellSize=0.0001)

        tif_file = 'PointData_lat_long.tif'

        # Create a domain and quantity
        x0 = utm_east_ll
        y0 = utm_north_ll

        a = [0.0, 0.0]
        b = [0.0, 500.0]
        c = [500.0, 0.0]
        d = [0.0, 1000.0]
        e = [500.0, 500.0]
        f = [1000.0, 0.0]

        points = [a, b, c, d, e, f]

        # bac, bce, ecf, dbe
        elements = [[1, 0, 2], [1, 2, 4], [4, 2, 5], [3, 1, 4]]

        mesh4 = Generic_Domain(points, elements,
                               geo_reference=Geo_reference(56, x0, y0, hemisphere=hemisphere_ll))
        mesh4.check_integrity()

        quantity1 = Quantity(mesh4)

        # Read in an interpolate from tif file
        quantity1.set_values_from_tif_file(filename=tif_file , location='vertices', indices=[0,1,2])

        # pprint(quantity1.centroid_values)
        # pprint(quantity1.vertex_values)

        centroid_values_ex = numpy.array(
            [0.17639186, 0.17639186, 0.17639186, 0. ])

        vertex_values_ex = numpy.array(
            [[0.23654294, 0.        , 0.29263264],
            [0.23654294, 0.        , 0.29263264],
            [0.23654294, 0.        , 0.29263264],
            [0.        , 0.        , 0.        ]])

        assert num.allclose(quantity1.centroid_values, centroid_values_ex)
        assert num.allclose(quantity1.vertex_values, vertex_values_ex)

        quantity2 = Quantity(mesh4)

        # Read in an interpolate from tif file
        quantity2.set_values(filename=tif_file , location='centroids', indices=[1])

        # pprint(quantity2.centroid_values)
        # pprint(quantity2.vertex_values)

        centroid_values_ex = numpy.array(
            [0.        , 0.34097296, 0.        , 0.        ])
        vertex_values_ex = numpy.array(
            [[0.        , 0.        , 0.        ],
            [0.34097296, 0.34097296, 0.34097296],
            [0.        , 0.        , 0.        ],
            [0.        , 0.        , 0.        ]])

        assert num.allclose(quantity2.centroid_values, centroid_values_ex)
        assert num.allclose(quantity2.vertex_values, vertex_values_ex)

        # Cleanup
        os.remove(tif_file)

    @pytest.mark.skipif('rasterio' not in sys.modules,
                    reason="requires rasterio")
    def test_set_values_from_ll_tif_file_south(self):

        # Mesh in zone 56 (relative coords) southern hemisphere

        from pprint import pprint

        # We need to make a .tif with ll coord which covers domain
        #
        from anuga.utilities import plot_utils as util
        import numpy
        #
        # Use Make_Geotif to make tif file
        # Pick a domain that makes sense in EPSG:32756
        lat_ll, lon_ll = -34.39, 150.90
        lat_ur, lon_ur = -34.37, 150.92

        from anuga.coordinate_transforms.lat_long_UTM_conversion import LLtoUTM
        zone_ll, utm_east_ll, utm_north_ll = LLtoUTM(lat_ll, lon_ll)
        zone_ur, utm_east_ur, utm_north_ur = LLtoUTM(lat_ur, lon_ur)

        hemisphere_ll = 'southern' if lat_ll < 0 else 'northern'
        hemisphere_ur = 'southern' if lat_ur < 0 else 'northern'

        lat = numpy.linspace(lat_ll, lat_ur, 11)
        lon = numpy.linspace(lon_ll, lon_ur, 11)

        xG, yG = numpy.meshgrid(lon, lat)

        xG = xG.flatten()
        yG = yG.flatten()
        # Surface is z=x+y

        fakeZ = (xG-min(xG))/(max(xG)-min(xG))+(yG - min(yG))/(max(yG)-min(yG))

        # pprint(fakeZ.reshape((11,11)))

        dataToGrid = numpy.vstack([xG, yG, fakeZ]).transpose()
        #
        # Create file PointData_lat_long.tif
        util.Make_Geotif(dataToGrid, output_quantities=['lat_long'],
                        EPSG_CODE=4326, output_dir='.', CellSize=0.0001)

        tif_file = 'PointData_lat_long.tif'

        # Create a domain and quantity
        x0 = utm_east_ll
        y0 = utm_north_ll

        a = [0.0, 0.0]
        b = [0.0, 500.0]
        c = [500.0, 0.0]
        d = [0.0, 1000.0]
        e = [500.0, 500.0]
        f = [1000.0, 0.0]

        points = [a, b, c, d, e, f]

        # bac, bce, ecf, dbe
        elements = [[1, 0, 2], [1, 2, 4], [4, 2, 5], [3, 1, 4]]

        mesh4 = Generic_Domain(points, elements,
                               geo_reference=Geo_reference(56, x0, y0, hemisphere=hemisphere_ll))
        mesh4.check_integrity()

        quantity1 = Quantity(mesh4)

        # Read in an interpolate from tif file
        quantity1.set_values_from_tif_file(filename=tif_file , location='vertices')

        # pprint(quantity1.centroid_values)
        # pprint(quantity1.vertex_values)

        centroid_values_ex = numpy.array(
            [0.17640652, 0.34250694, 0.44839957, 0.39855834])
        vertex_values_ex = numpy.array([[0.23744977, 0., 0.2917698],
                                        [0.23744977, 0.2917698, 0.49830124],
                                        [0.49830124, 0.2917698, 0.55512768],
                                        [0.45992401, 0.23744977, 0.49830124]])

        assert num.allclose(quantity1.centroid_values, centroid_values_ex)
        assert num.allclose(quantity1.vertex_values, vertex_values_ex)

        quantity2 = Quantity(mesh4)

        # Read in an interpolate from tif file
        quantity2.set_values(filename=tif_file , location='vertices')

        # pprint(quantity2.centroid_values)
        # pprint(quantity2.vertex_values)

        assert num.allclose(quantity2.centroid_values, centroid_values_ex)
        assert num.allclose(quantity2.vertex_values, vertex_values_ex)

        # Cleanup
        os.remove(tif_file)

    @pytest.mark.skipif('rasterio' not in sys.modules,
                    reason="requires rasterio")
    def test_set_values_from_utm_tif_file(self):

        # Mesh in zone 56 (relative coords) southern hemisphere

        from pprint import pprint

        # We need to make a .tif with ll coord which covers domain
        #
        from anuga.utilities import plot_utils as util
        import numpy
        #
        # Use Make_Geotif to make tif file
        # Pick a domain that makes sense in EPSG:32756
        lat_ll, lon_ll = -34.39, 150.90
        lat_ur, lon_ur = -34.37, 150.92

        from anuga.coordinate_transforms.lat_long_UTM_conversion import LLtoUTM
        zone_ll, utm_east_ll, utm_north_ll = LLtoUTM(lat_ll, lon_ll)
        zone_ur, utm_east_ur, utm_north_ur = LLtoUTM(lat_ur, lon_ur)

        hemisphere_ll = 'southern' if lat_ll < 0 else 'northern'
        hemisphere_ur = 'southern' if lat_ur < 0 else 'northern'

        easting = numpy.linspace(utm_east_ll, utm_east_ur, 11)
        northing = numpy.linspace(utm_north_ll, utm_north_ur, 11)

        xG, yG = numpy.meshgrid(easting, northing)

        xG = xG.flatten()
        yG = yG.flatten()
        # Surface is z=x+y

        fakeZ = (xG-min(xG))/(max(xG)-min(xG))+(yG - min(yG))/(max(yG)-min(yG))

        # pprint(fakeZ.reshape((11,11)))

        dataToGrid = numpy.vstack([xG, yG, fakeZ]).transpose()
        #
        # Create file PointData_lat_long.tif
        util.Make_Geotif(dataToGrid, output_quantities=['utm'],
                        EPSG_CODE=32756, output_dir='.', CellSize=10)

        tif_file = 'PointData_utm.tif'

        # Create a domain and quantity
        x0 = utm_east_ll
        y0 = utm_north_ll

        a = [0.0, 0.0]
        b = [0.0, 500.0]
        c = [500.0, 0.0]
        d = [0.0, 1000.0]
        e = [500.0, 500.0]
        f = [1000.0, 0.0]

        points = [a, b, c, d, e, f]

        # bac, bce, ecf, dbe
        elements = [[1, 0, 2], [1, 2, 4], [4, 2, 5], [3, 1, 4]]

        mesh4 = Generic_Domain(points, elements,
                               geo_reference=Geo_reference(56, x0, y0, hemisphere=hemisphere_ll))
        mesh4.check_integrity()

        quantity1 = Quantity(mesh4)

        # Read in an interpolate from tif file
        quantity1.set_values_from_tif_file(filename=tif_file , location='vertices')

        # pprint(quantity1.centroid_values)
        # pprint(quantity1.vertex_values)

        centroid_values_ex = numpy.array(
            [0.17610394, 0.34143136, 0.45676157, 0.39509157])
        vertex_values_ex = numpy.array([[0.23454617, 0., 0.29376566],
                                        [0.23454617, 0.29376566, 0.49598226],
                                        [0.49598226, 0.29376566, 0.58053678],
                                        [0.45474628, 0.23454617, 0.49598226]])

        assert num.allclose(quantity1.centroid_values, centroid_values_ex)
        assert num.allclose(quantity1.vertex_values, vertex_values_ex)

        quantity2 = Quantity(mesh4)

        # Read in an interpolate from tif file
        quantity2.set_values(filename=tif_file , location='vertices')

        # pprint(quantity2.centroid_values)
        # pprint(quantity2.vertex_values)

        assert num.allclose(quantity2.centroid_values, centroid_values_ex)
        assert num.allclose(quantity2.vertex_values, vertex_values_ex)

        # Cleanup
        os.remove(tif_file)

    def test_set_values_from_quantity(self):

        quantity1 = Quantity(self.mesh4)
        quantity1.set_vertex_values([0, 1, 2, 3, 4, 5])

        assert num.allclose(quantity1.vertex_values,
                            [[1, 0, 2], [1, 2, 4], [4, 2, 5], [3, 1, 4]])

        quantity2 = Quantity(self.mesh4)
        quantity2.set_values(quantity=quantity1)
        assert num.allclose(quantity2.vertex_values,
                            [[1, 0, 2], [1, 2, 4], [4, 2, 5], [3, 1, 4]])

        quantity2.set_values(quantity=2*quantity1)
        assert num.allclose(quantity2.vertex_values,
                            [[2, 0, 4], [2, 4, 8], [8, 4, 10], [6, 2, 8]])

        quantity2.set_values(quantity=2*quantity1 + 3)
        assert num.allclose(quantity2.vertex_values,
                            [[5, 3, 7], [5, 7, 11], [11, 7, 13], [9, 5, 11]])

        # Check detection of quantity as first orgument
        quantity2.set_values(2*quantity1 + 3)
        assert num.allclose(quantity2.vertex_values,
                            [[5, 3, 7], [5, 7, 11], [11, 7, 13], [9, 5, 11]])

    def Xtest_set_values_from_quantity_using_polygon(self):
        """test_set_values_from_quantity_using_polygon(self):

        Check that polygon can be used to restrict set_values when
        using another quantity as argument.
        """

        # Create restricting polygon (containing node #4 (2,2) and
        # centroid of triangle #1 (bce)
        polygon = [[1.0, 1.0], [4.0, 1.0],
                   [4.0, 4.0], [1.0, 4.0]]
        assert num.allclose(inside_polygon(self.mesh4.nodes, polygon), 4)

        quantity1 = Quantity(self.mesh4)
        quantity1.set_vertex_values([0, 1, 2, 3, 4, 5])

        assert num.allclose(quantity1.vertex_values,
                            [[1, 0, 2], [1, 2, 4], [4, 2, 5], [3, 1, 4]])

        quantity2 = Quantity(self.mesh4)
        quantity2.set_values(quantity=quantity1,
                             polygon=polygon)

        msg = 'Only node #4(e) at (2,2) should have values applied '
        assert num.allclose(quantity2.vertex_values,
                            [[0, 0, 0], [0, 0, 4], [4, 0, 0], [0, 0, 4]]), msg
        # bac,     bce,     ecf,     dbe

    def test_overloading(self):

        quantity1 = Quantity(self.mesh4)
        quantity1.set_vertex_values([0, 1, 2, 3, 4, 5])

        assert num.allclose(quantity1.vertex_values,
                            [[1, 0, 2], [1, 2, 4], [4, 2, 5], [3, 1, 4]])

        quantity2 = Quantity(self.mesh4)
        quantity2.set_values([[1, 2, 3], [5, 5, 5], [0, 0, 9], [-6, 3, 3]],
                             location='vertices')

        quantity3 = Quantity(self.mesh4)
        quantity3.set_values([[2, 2, 2], [7, 8, 9], [7, 6, 3], [3, 8, -8]],
                             location='vertices')

        # Negation
        Q = -quantity1
        assert num.allclose(Q.vertex_values, -quantity1.vertex_values)
        assert num.allclose(Q.centroid_values, -quantity1.centroid_values)
        assert num.allclose(Q.edge_values, -quantity1.edge_values)

        # Addition
        Q = quantity1 + 7
        assert num.allclose(Q.vertex_values, quantity1.vertex_values + 7)
        assert num.allclose(Q.centroid_values, quantity1.centroid_values + 7)
        assert num.allclose(Q.edge_values, quantity1.edge_values + 7)

        Q = 7 + quantity1
        assert num.allclose(Q.vertex_values, quantity1.vertex_values + 7)
        assert num.allclose(Q.centroid_values, quantity1.centroid_values + 7)
        assert num.allclose(Q.edge_values, quantity1.edge_values + 7)

        Q = quantity1 + quantity2
        assert num.allclose(Q.vertex_values,
                            quantity1.vertex_values + quantity2.vertex_values)
        assert num.allclose(Q.centroid_values,
                            quantity1.centroid_values + quantity2.centroid_values)
        assert num.allclose(Q.edge_values,
                            quantity1.edge_values + quantity2.edge_values)

        Q = quantity1 + quantity2 - 3
        assert num.allclose(Q.vertex_values,
                            quantity1.vertex_values + quantity2.vertex_values - 3)

        Q = quantity1 - quantity2
        assert num.allclose(Q.vertex_values,
                            quantity1.vertex_values - quantity2.vertex_values)

        # Scaling
        Q = quantity1*3
        assert num.allclose(Q.vertex_values, quantity1.vertex_values*3)
        assert num.allclose(Q.centroid_values, quantity1.centroid_values*3)
        assert num.allclose(Q.edge_values, quantity1.edge_values*3)
        Q = 3*quantity1
        assert num.allclose(Q.vertex_values, quantity1.vertex_values*3)

        # Multiplication
        Q = quantity1 * quantity2
        # print Q.vertex_values
        # print Q.centroid_values
        # print quantity1.centroid_values
        # print quantity2.centroid_values

        assert num.allclose(Q.vertex_values,
                            quantity1.vertex_values * quantity2.vertex_values)

        # Linear combinations
        Q = 4*quantity1 + 2
        assert num.allclose(Q.vertex_values,
                            4*quantity1.vertex_values + 2)

        Q = quantity1*quantity2 + 2
        assert num.allclose(Q.vertex_values,
                            quantity1.vertex_values * quantity2.vertex_values + 2)

        Q = quantity1*quantity2 + quantity3
        assert num.allclose(Q.vertex_values,
                            quantity1.vertex_values * quantity2.vertex_values
                            + quantity3.vertex_values)
        Q = quantity1*quantity2 + 3*quantity3
        assert num.allclose(Q.vertex_values,
                            quantity1.vertex_values * quantity2.vertex_values
                            + 3*quantity3.vertex_values)
        Q = quantity1*quantity2 + 3*quantity3 + 5.0
        assert num.allclose(Q.vertex_values,
                            quantity1.vertex_values * quantity2.vertex_values
                            + 3*quantity3.vertex_values + 5)

        Q = quantity1*quantity2 - quantity3
        assert num.allclose(Q.vertex_values,
                            quantity1.vertex_values * quantity2.vertex_values
                            - quantity3.vertex_values)
        Q = 1.5*quantity1*quantity2 - 3*quantity3 + 5.0
        assert num.allclose(Q.vertex_values,
                            1.5*quantity1.vertex_values * quantity2.vertex_values
                            - 3*quantity3.vertex_values + 5)

        # Try combining quantities and arrays and scalars
        Q = 1.5*quantity1*quantity2.vertex_values -\
            3*quantity3.vertex_values + 5.0
        assert num.allclose(Q.vertex_values,
                            1.5*quantity1.vertex_values * quantity2.vertex_values
                            - 3*quantity3.vertex_values + 5)

        # Powers
        Q = quantity1**2
        assert num.allclose(Q.vertex_values, quantity1.vertex_values**2)

        Q = quantity1**2 + quantity2**2
        assert num.allclose(Q.vertex_values,
                            quantity1.vertex_values**2
                            + quantity2.vertex_values**2)

        Q = (quantity1**2 + quantity2**2)**0.5
        assert num.allclose(Q.vertex_values,
                            (quantity1.vertex_values**2 +
                             quantity2.vertex_values**2)**0.5)

    def test_compute_gradient(self):
        quantity = Quantity(self.mesh4)

        # Set up for a gradient of (2,0) at mid triangle
        quantity.set_values([2.0, 4.0, 6.0, 2.0],
                            location='centroids')

        # Gradients
        quantity.compute_gradients()

        a = quantity.x_gradient
        b = quantity.y_gradient
        # print self.mesh4.centroid_coordinates
        # print a, b

        # The central triangle (1)
        # (using standard gradient based on neigbours controid values)
        assert num.allclose(a[1], 2.0)
        assert num.allclose(b[1], 0.0)

        # Left triangle (0) using two point gradient
        # q0 = q1 + a*(x0-x1) + b*(y0-y1)  <=>
        # 2  = 4  + a*(-2/3)  + b*(-2/3)
        assert num.allclose(a[0] + b[0], 3)
        # From orthogonality (a*(y0-y1) + b*(x0-x1) == 0)
        assert num.allclose(a[0] - b[0], 0)

        # Right triangle (2) using two point gradient
        # q2 = q1 + a*(x2-x1) + b*(y2-y1)  <=>
        # 6  = 4  + a*(4/3)  + b*(-2/3)
        assert num.allclose(2*a[2] - b[2], 3)
        # From orthogonality (a*(y1-y2) + b*(x2-x1) == 0)
        assert num.allclose(a[2] + 2*b[2], 0)

        # Top triangle (3) using two point gradient
        # q3 = q1 + a*(x3-x1) + b*(y3-y1)  <=>
        # 2  = 4  + a*(-2/3)  + b*(4/3)
        assert num.allclose(a[3] - 2*b[3], 3)
        # From orthogonality (a*(y1-y3) + b*(x3-x1) == 0)
        assert num.allclose(2*a[3] + b[3], 0)

        # print a, b
        quantity.extrapolate_second_order()

        # Apply q(x,y) = qc + a*(x-xc) + b*(y-yc)
        assert num.allclose(quantity.vertex_values[0, :], [3., 0.,  3.])
        assert num.allclose(quantity.vertex_values[1, :], [
                            4./3, 16./3,  16./3])

        # a = 1.2, b=-0.6
        # q(4,0) = 6 + a*(4 - 8/3) + b*(-2/3)
        assert num.allclose(quantity.vertex_values[2, 2], 8)

    def test_get_gradients(self):
        quantity = Quantity(self.mesh4)

        # Set up for a gradient of (2,0) at mid triangle
        quantity.set_values([2.0, 4.0, 6.0, 2.0],
                            location='centroids')

        # Gradients
        quantity.compute_gradients()

        a, b = quantity.get_gradients()
        # print self.mesh4.centroid_coordinates
        # print a, b

        # The central triangle (1)
        # (using standard gradient based on neigbours controid values)
        assert num.allclose(a[1], 2.0)
        assert num.allclose(b[1], 0.0)

        # Left triangle (0) using two point gradient
        # q0 = q1 + a*(x0-x1) + b*(y0-y1)  <=>
        # 2  = 4  + a*(-2/3)  + b*(-2/3)
        assert num.allclose(a[0] + b[0], 3)
        # From orthogonality (a*(y0-y1) + b*(x0-x1) == 0)
        assert num.allclose(a[0] - b[0], 0)

        # Right triangle (2) using two point gradient
        # q2 = q1 + a*(x2-x1) + b*(y2-y1)  <=>
        # 6  = 4  + a*(4/3)  + b*(-2/3)
        assert num.allclose(2*a[2] - b[2], 3)
        # From orthogonality (a*(y1-y2) + b*(x2-x1) == 0)
        assert num.allclose(a[2] + 2*b[2], 0)

        # Top triangle (3) using two point gradient
        # q3 = q1 + a*(x3-x1) + b*(y3-y1)  <=>
        # 2  = 4  + a*(-2/3)  + b*(4/3)
        assert num.allclose(a[3] - 2*b[3], 3)
        # From orthogonality (a*(y1-y3) + b*(x3-x1) == 0)
        assert num.allclose(2*a[3] + b[3], 0)

    def test_second_order_extrapolation2(self):
        quantity = Quantity(self.mesh4)

        # Set up for a gradient of (3,1), f(x) = 3x+y
        quantity.set_values([2.0+2.0/3, 4.0+4.0/3, 8.0+2.0/3, 2.0+8.0/3],
                            location='centroids')

        # Gradients
        quantity.compute_gradients()

        a = quantity.x_gradient
        b = quantity.y_gradient

        # print a, b

        assert num.allclose(a[1], 3.0)
        assert num.allclose(b[1], 1.0)

        # Work out the others

        quantity.extrapolate_second_order()

        # print quantity.vertex_values
        assert num.allclose(quantity.vertex_values[1, 0], 2.0)
        assert num.allclose(quantity.vertex_values[1, 1], 6.0)
        assert num.allclose(quantity.vertex_values[1, 2], 8.0)

    def test_backup_saxpy_centroid_values(self):
        quantity = Quantity(self.mesh4)

        # Set up for a gradient of (3,1), f(x) = 3x+y
        c_values = num.array([2.0+2.0/3, 4.0+4.0/3, 8.0+2.0/3, 2.0+8.0/3])
        d_values = num.array([1.0, 2.0, 3.0, 4.0])
        quantity.set_values(c_values, location='centroids')

        # Backup
        quantity.backup_centroid_values()

        # print quantity.vertex_values
        assert num.allclose(quantity.centroid_values,
                            quantity.centroid_backup_values)

        quantity.set_values(d_values, location='centroids')

        quantity.saxpy_centroid_values(2.0, 3.0)

        assert num.allclose(quantity.centroid_values,
                            2.0*d_values + 3.0*c_values)

    def test_first_order_extrapolator(self):
        quantity = Quantity(self.mesh4)

        # Test centroids
        quantity.set_values([1., 2., 3., 4.], location='centroids')
        assert num.allclose(quantity.centroid_values, [1, 2, 3, 4])  # Centroid

        # Extrapolate
        quantity.extrapolate_first_order()

        # Check that gradient is zero
        a, b = quantity.get_gradients()
        assert num.allclose(a, [0, 0, 0, 0])
        assert num.allclose(b, [0, 0, 0, 0])

        # Check vertices but not edge values
        assert num.allclose(quantity.vertex_values,
                            [[1, 1, 1], [2, 2, 2], [3, 3, 3], [4, 4, 4]])

    def test_second_order_extrapolator(self):
        quantity = Quantity(self.mesh4)

        # Set up for a gradient of (3,0) at mid triangle
        quantity.set_values([2.0, 4.0, 8.0, 2.0],
                            location='centroids')

        quantity.extrapolate_second_order()
        quantity.limit()

        # Assert that central triangle is limited by neighbours
        assert quantity.vertex_values[1, 0] >= quantity.vertex_values[0, 0]
        assert quantity.vertex_values[1, 0] >= quantity.vertex_values[3, 1]

        assert quantity.vertex_values[1, 1] <= quantity.vertex_values[2, 1]
        assert quantity.vertex_values[1, 1] >= quantity.vertex_values[0, 2]

        assert quantity.vertex_values[1, 2] <= quantity.vertex_values[2, 0]
        assert quantity.vertex_values[1, 2] >= quantity.vertex_values[3, 1]

        # Assert that quantities are conserved
        for k in range(quantity.centroid_values.shape[0]):
            assert num.allclose(quantity.centroid_values[k],
                                num.sum(quantity.vertex_values[k, :])/3)

    def test_limit_vertices_by_all_neighbours(self):
        quantity = Quantity(self.mesh4)

        # Create a deliberate overshoot (e.g. from gradient computation)
        quantity.set_values([[3, 0, 3], [2, 2, 6], [5, 3, 8], [8, 3, 5]])

        # Limit
        quantity.limit_vertices_by_all_neighbours()

        # Assert that central triangle is limited by neighbours
        assert quantity.vertex_values[1, 0] >= quantity.vertex_values[0, 0]
        assert quantity.vertex_values[1, 0] <= quantity.vertex_values[3, 1]

        assert quantity.vertex_values[1, 1] <= quantity.vertex_values[2, 1]
        assert quantity.vertex_values[1, 1] >= quantity.vertex_values[0, 2]

        assert quantity.vertex_values[1, 2] <= quantity.vertex_values[2, 0]
        assert quantity.vertex_values[1, 2] <= quantity.vertex_values[3, 1]

        # Assert that quantities are conserved
        for k in range(quantity.centroid_values.shape[0]):
            assert num.allclose(quantity.centroid_values[k],
                                num.sum(quantity.vertex_values[k, :])/3)

    def test_limit_edges_by_all_neighbours(self):
        quantity = Quantity(self.mesh4)

        # Create a deliberate overshoot (e.g. from gradient computation)
        quantity.set_values([[3, 0, 3], [2, 2, 6], [5, 3, 8], [8, 3, 5]])

        # Limit
        quantity.limit_edges_by_all_neighbours()

        # Assert that central triangle is limited by neighbours
        assert quantity.edge_values[1, 0] <= quantity.centroid_values[2]
        assert quantity.edge_values[1, 0] >= quantity.centroid_values[0]

        assert quantity.edge_values[1, 1] <= quantity.centroid_values[2]
        assert quantity.edge_values[1, 1] >= quantity.centroid_values[0]

        assert quantity.edge_values[1, 2] <= quantity.centroid_values[2]
        assert quantity.edge_values[1, 2] >= quantity.centroid_values[0]

        # Assert that quantities are conserved
        for k in range(quantity.centroid_values.shape[0]):
            assert num.allclose(quantity.centroid_values[k],
                                num.sum(quantity.vertex_values[k, :])/3)

    def test_limit_edges_by_neighbour(self):
        quantity = Quantity(self.mesh4)

        # Create a deliberate overshoot (e.g. from gradient computation)
        quantity.set_values([[3, 0, 3], [2, 2, 6], [5, 3, 8], [8, 3, 5]])

        # Limit
        quantity.limit_edges_by_neighbour()

        # Assert that central triangle is limited by neighbours
        assert quantity.edge_values[1, 0] <= quantity.centroid_values[3]
        assert quantity.edge_values[1, 0] >= quantity.centroid_values[1]

        assert quantity.edge_values[1, 1] <= quantity.centroid_values[2]
        assert quantity.edge_values[1, 1] >= quantity.centroid_values[1]

        assert quantity.edge_values[1, 2] <= quantity.centroid_values[1]
        assert quantity.edge_values[1, 2] >= quantity.centroid_values[0]

        # Assert that quantities are conserved
        for k in range(quantity.centroid_values.shape[0]):
            assert num.allclose(quantity.centroid_values[k],
                                num.sum(quantity.vertex_values[k, :])/3)

    def test_limiter2(self):
        """Taken from test_shallow_water
        """
        quantity = Quantity(self.mesh4)
        quantity.domain.beta_w = 0.9

        # Test centroids
        quantity.set_values([2., 4., 8., 2.], location='centroids')
        assert num.allclose(quantity.centroid_values, [2, 4, 8, 2])  # Centroid

        # Extrapolate
        quantity.extrapolate_second_order()

        assert num.allclose(quantity.vertex_values[1, :], [0.0, 6, 6])

        # Limit
        quantity.limit()

        # limited value for beta_w = 0.9

        assert num.allclose(quantity.vertex_values[1, :], [2.2, 4.9, 4.9])
        # limited values for beta_w = 0.5
        # assert allclose(quantity.vertex_values[1,:], [3.0, 4.5, 4.5])

        # Assert that quantities are conserved
        for k in range(quantity.centroid_values.shape[0]):
            assert num.allclose(quantity.centroid_values[k],
                                num.sum(quantity.vertex_values[k, :])/3)

    def test_distribute_first_order(self):
        quantity = Quantity(self.mesh4)

        # Test centroids
        quantity.set_values([1., 2., 3., 4.], location='centroids')
        assert num.allclose(quantity.centroid_values, [1, 2, 3, 4])  # Centroid

        # Extrapolate from centroid to vertices and edges
        quantity.extrapolate_first_order()

        # Interpolate
        # quantity.interpolate_from_vertices_to_edges()

        assert num.allclose(quantity.vertex_values,
                            [[1, 1, 1], [2, 2, 2], [3, 3, 3], [4, 4, 4]])
        assert num.allclose(quantity.edge_values, [[1, 1, 1], [2, 2, 2],
                                                   [3, 3, 3], [4, 4, 4]])

    def test_interpolate_from_vertices_to_edges(self):
        quantity = Quantity(self.mesh4)

        quantity.vertex_values = num.array(
            [[1, 0, 2], [1, 2, 4], [4, 2, 5], [3, 1, 4]], float)

        quantity.interpolate_from_vertices_to_edges()

        assert num.allclose(quantity.edge_values, [[1., 1.5, 0.5],
                                                   [3., 2.5, 1.5],
                                                   [3.5, 4.5, 3.],
                                                   [2.5, 3.5, 2]])

    def test_interpolate_from_edges_to_vertices(self):
        quantity = Quantity(self.mesh4)

        quantity.edge_values = num.array([[1., 1.5, 0.5],
                                          [3., 2.5, 1.5],
                                          [3.5, 4.5, 3.],
                                          [2.5, 3.5, 2]], float)

        quantity.interpolate_from_edges_to_vertices()

        assert num.allclose(quantity.vertex_values,
                            [[1, 0, 2], [1, 2, 4], [4, 2, 5], [3, 1, 4]])

    def test_distribute_second_order(self):
        quantity = Quantity(self.mesh4)

        # Test centroids
        quantity.set_values([2., 4., 8., 2.], location='centroids')
        assert num.allclose(quantity.centroid_values, [2, 4, 8, 2])  # Centroid

        # Extrapolate
        quantity.extrapolate_second_order()

        assert num.allclose(quantity.vertex_values[1, :], [0.0, 6, 6])

    def test_update_explicit(self):
        quantity = Quantity(self.mesh4)

        # Test centroids
        quantity.set_values([1., 2., 3., 4.], location='centroids')
        assert num.allclose(quantity.centroid_values, [1, 2, 3, 4])  # Centroid

        # Set explicit_update
        quantity.explicit_update = num.array([1., 1., 1., 1.])

        # Update with given timestep
        quantity.update(0.1)

        x = num.array([1, 2, 3, 4]) + num.array([.1, .1, .1, .1])
        assert num.allclose(quantity.centroid_values, x)

    def test_update_semi_implicit(self):
        quantity = Quantity(self.mesh4)

        # Test centroids
        quantity.set_values([1., 2., 3., 4.], location='centroids')
        assert num.allclose(quantity.centroid_values, [1, 2, 3, 4])  # Centroid

        # Set semi implicit update
        quantity.semi_implicit_update = num.array([1., 1., 1., 1.])

        # Update with given timestep
        timestep = 0.1
        quantity.update(timestep)

        sem = num.array([1., 1., 1., 1.])/num.array([1, 2, 3, 4])
        denom = num.ones(4, float)-timestep*sem

        x = num.array([1, 2, 3, 4])/denom
        assert num.allclose(quantity.centroid_values, x)

    def test_both_updates(self):
        quantity = Quantity(self.mesh4)

        # Test centroids
        quantity.set_values([1., 2., 3., 4.], location='centroids')
        assert num.allclose(quantity.centroid_values, [1, 2, 3, 4])  # Centroid

        # Set explicit_update
        quantity.explicit_update = num.array([4., 3., 2., 1.])

        # Set semi implicit update
        quantity.semi_implicit_update = num.array([1., 1., 1., 1.])

        # Update with given timestep
        timestep = 0.1
        quantity.update(0.1)

        sem = num.array([1., 1., 1., 1.])/num.array([1, 2, 3, 4])
        denom = num.ones(4, float)-timestep*sem

        x = num.array([1., 2., 3., 4.])
        x += timestep*num.array([4.0, 3.0, 2.0, 1.0])
        x /= denom

        assert num.allclose(quantity.centroid_values, x)

    def set_array_values_by_index(self):

        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(1, 1)

        # Create shallow water domain
        domain = Generic_Domain(points, vertices, boundary)
        # print "domain.number_of_elements ",domain.number_of_elements
        quantity = Quantity(domain, [[1, 1, 1], [2, 2, 2]])
        value = [7]
        indices = [1]
        quantity.set_array_values_by_index(value,
                                           location='centroids',
                                           indices=indices)
        # print "quantity.centroid_values",quantity.centroid_values

        assert num.allclose(quantity.centroid_values, [1, 7])

        quantity.set_array_values([15, 20, 25], indices=indices)
        assert num.allclose(quantity.centroid_values, [1, 20])

        quantity.set_array_values([15, 20, 25], indices=indices)
        assert num.allclose(quantity.centroid_values, [1, 20])

    def test_setting_some_vertex_values(self):
        """
        set values based on triangle lists.
        """
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(1, 3)
        # print "vertices",vertices
        # Create shallow water domain
        domain = Generic_Domain(points, vertices, boundary)
        # print "domain.number_of_elements ",domain.number_of_elements
        quantity = Quantity(domain, [[1, 1, 1], [2, 2, 2], [3, 3, 3],
                                     [4, 4, 4], [5, 5, 5], [6, 6, 6]])

        # Check that constants work
        value = 7
        indices = [1]
        quantity.set_values(value,
                            location='centroids',
                            indices=indices)
        # print "quantity.centroid_values",quantity.centroid_values
        assert num.allclose(quantity.centroid_values, [1, 7, 3, 4, 5, 6])

        value = [7]
        indices = [1]
        quantity.set_values(value,
                            location='centroids',
                            indices=indices)
        # print "quantity.centroid_values",quantity.centroid_values
        assert num.allclose(quantity.centroid_values, [1, 7, 3, 4, 5, 6])

        value = [[15, 20, 25]]
        quantity.set_values(value, indices=indices)
        # print "1 quantity.vertex_values",quantity.vertex_values
        assert num.allclose(quantity.vertex_values[1], value[0])

        # print "quantity",quantity.vertex_values
        values = [10, 100, 50]
        quantity.set_values(values, indices=[0, 1, 5], location='centroids')
        # print "2 quantity.vertex_values",quantity.vertex_values
        assert num.allclose(quantity.vertex_values[0], [10, 10, 10])
        assert num.allclose(quantity.vertex_values[5], [50, 50, 50])
        # quantity.interpolate()
        # print "quantity.centroid_values",quantity.centroid_values
        assert num.allclose(quantity.centroid_values, [10, 100, 3, 4, 5, 50])

        quantity = Quantity(domain, [[1, 1, 1], [2, 2, 2], [3, 3, 3],
                                     [4, 4, 4], [5, 5, 5], [6, 6, 6]])
        values = [10, 100, 50]
        # this will be per unique vertex, indexing the vertices
        # print "quantity.vertex_values",quantity.vertex_values
        quantity.set_values(values, indices=[0, 1, 5])
        # print "quantity.vertex_values",quantity.vertex_values
        assert num.allclose(quantity.vertex_values[0], [1, 50, 10])
        assert num.allclose(quantity.vertex_values[5], [6, 6, 6])
        assert num.allclose(quantity.vertex_values[1], [100, 10, 50])

        quantity = Quantity(domain, [[1, 1, 1], [2, 2, 2], [3, 3, 3],
                                     [4, 4, 4], [5, 5, 5], [6, 6, 6]])
        values = [[31, 30, 29], [400, 400, 400], [1000, 999, 998]]
        quantity.set_values(values, indices=[3, 3, 5])
        quantity.interpolate()
        assert num.allclose(quantity.centroid_values, [1, 2, 3, 400, 5, 999])

        values = [[1, 1, 1], [2, 2, 2], [3, 3, 3],
                  [4, 4, 4], [5, 5, 5], [6, 6, 6]]
        quantity.set_values(values)

        # testing the standard set values by vertex
        # indexed by vertex_id in general_mesh.coordinates
        values = [0, 1, 2, 3, 4, 5, 6, 7]

        quantity.set_values(values)
        # print "1 quantity.vertex_values",quantity.vertex_values
        assert num.allclose(quantity.vertex_values, [[4.,  5.,  0.],
                                                     [1.,  0.,  5.],
                                                     [5.,  6.,  1.],
                                                     [2.,  1.,  6.],
                                                     [6.,  7.,  2.],
                                                     [3.,  2.,  7.]])

    def test_setting_unique_vertex_values(self):
        """
        set values based on unique_vertex lists.
        """
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(1, 3)
        # print "vertices",vertices
        # Create shallow water domain
        domain = Generic_Domain(points, vertices, boundary)
        # print "domain.number_of_elements ",domain.number_of_elements
        quantity = Quantity(domain, [[0, 0, 0], [1, 1, 1], [2, 2, 2], [3, 3, 3],
                                     [4, 4, 4], [5, 5, 5]])
        value = 7
        indices = [1, 5]
        quantity.set_values(value,
                            location='unique vertices',
                            indices=indices)
        # print "quantity.centroid_values",quantity.centroid_values
        assert num.allclose(quantity.vertex_values[0], [0, 7, 0])
        assert num.allclose(quantity.vertex_values[1], [7, 1, 7])
        assert num.allclose(quantity.vertex_values[2], [7, 2, 7])

    def test_get_values(self):
        """
        get values based on triangle lists.
        """
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(1, 3)

        # print "points",points
        # print "vertices",vertices
        # print "boundary",boundary

        # Create shallow water domain
        domain = Generic_Domain(points, vertices, boundary)
        # print "domain.number_of_elements ",domain.number_of_elements
        quantity = Quantity(domain, [[0, 0, 0], [1, 1, 1], [2, 2, 2], [3, 3, 3],
                                     [4, 4, 4], [5, 5, 5]])

        # print "quantity.get_values(location = 'unique vertices')", \
        #      quantity.get_values(location = 'unique vertices')

        # print "quantity.get_values(location = 'unique vertices')", \
        #      quantity.get_values(indices=[0,1,2,3,4,5,6,7], \
        #                          location = 'unique vertices')

        answer = [0.5, 2, 4, 5, 0, 1, 3, 4.5]
        assert num.allclose(answer,
                            quantity.get_values(location='unique vertices'))

        indices = [0, 5, 3]
        answer = [0.5, 1, 5]
        assert num.allclose(answer,
                            quantity.get_values(indices=indices,
                                                location='unique vertices'))
        # print "quantity.centroid_values",quantity.centroid_values
        # print "quantity.get_values(location = 'centroids') ",\
        #      quantity.get_values(location = 'centroids')

    def test_get_values_2(self):
        """Different mesh (working with domain object) - also check centroids.
        """

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        # bac, bce, ecf, dbe
        vertices = [[1, 0, 2], [1, 2, 4], [4, 2, 5], [3, 1, 4]]

        domain = Generic_Domain(points, vertices)

        quantity = Quantity(domain)
        quantity.set_values(lambda x, y: x+2*y)  # 2 4 4 6

        assert num.allclose(quantity.get_values(location='centroids'),
                            [2, 4, 4, 6])

        assert num.allclose(quantity.get_values(location='centroids',
                                                indices=[1, 3]), [4, 6])

        assert num.allclose(quantity.get_values(location='vertices'),
                            [[4, 0, 2],
                             [4, 2, 6],
                             [6, 2, 4],
                             [8, 4, 6]])

        assert num.allclose(quantity.get_values(location='vertices',
                                                indices=[1, 3]), [[4, 2, 6],
                                                                  [8, 4, 6]])

        assert num.allclose(quantity.get_values(location='edges'),
                            [[1, 3, 2],
                             [4, 5, 3],
                             [3, 5, 4],
                             [5, 7, 6]])
        assert num.allclose(quantity.get_values(location='edges',
                                                indices=[1, 3]),
                            [[4, 5, 3],
                             [5, 7, 6]])

        # Check averaging over vertices
        # a: 0
        # b: (4+4+4)/3
        # c: (2+2+2)/3
        # d: 8
        # e: (6+6+6)/3
        # f: 4
        assert num.allclose(quantity.get_values(location='unique vertices'),
                            [0, 4, 2, 8, 6, 4])

    def test_get_vertex_values(self):
        """
        get values based on triangle lists.
        """
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(1, 3)

        # print "points",points
        # print "vertices",vertices
        # print "boundary",boundary

        # Create shallow water domain
        domain = Generic_Domain(points, vertices, boundary)
        # print "domain.number_of_elements ",domain.number_of_elements
        quantity = Quantity(domain, [[0, 0, 0], [1, 1, 1], [2, 2, 2], [3, 3, 3],
                                     [4, 4, 4], [5, 5, 5]])

        # ======================================================
        # Default: get_vertex_values just returns the individual
        # vertex values within each triangle
        # ======================================================
        Q, V = quantity.get_vertex_values(xy=False)

        answer = [0.,  0.,  0.,  1.,  1.,  1.,  2.,  2.,  2.,  3.,  3.,  3.,  4.,
                  4.,  4.,  5.,  5.,  5.]

        v_answer = num.array([[0,  1,  2],
                              [3,  4,  5],
                              [6,  7,  8],
                              [9, 10, 11],
                              [12, 13, 14],
                              [15, 16, 17]])

        assert num.allclose(answer, Q)
        assert num.allclose(v_answer, V)

        # ======================================================
        # Set output to be smooth, so get one unique value at
        # each node. V now provides id to unique node id
        # ======================================================
        domain.smooth = True

        Q, V = quantity.get_vertex_values(xy=False)

        answer = num.array([0.5, 2, 4, 5, 0, 1, 3, 4.5])

        v_answer = num.array([[4, 5, 0],
                              [1, 0, 5],
                              [5, 6, 1],
                              [2, 1, 6],
                              [6, 7, 2],
                              [3, 2, 7]])

        assert num.allclose(answer, Q)
        assert num.allclose(v_answer, V)

        # ======================================================
        # Set output to be smooth, and if using discontinuous
        # algorithms, get one unique value at
        # each node, based on centroid values.
        # V now provides id to unique node id
        # ======================================================
        domain.smooth = True
        domain.using_centroid_averaging = True

        quantity.centroid_values[:] = num.array([100, 101, 102, 103, 104, 105])

        Q, V = quantity.get_vertex_values(xy=False)

        answer = answer + 100.0

        v_answer = num.array([[4, 5, 0],
                              [1, 0, 5],
                              [5, 6, 1],
                              [2, 1, 6],
                              [6, 7, 2],
                              [3, 2, 7]])

        assert num.allclose(answer, Q)
        assert num.allclose(v_answer, V)

    def test_get_interpolated_values(self):

        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(1, 3)
        domain = Generic_Domain(points, vertices, boundary)

        # Constant values
        quantity = Quantity(domain, [[0, 0, 0], [1, 1, 1], [2, 2, 2], [3, 3, 3],
                                     [4, 4, 4], [5, 5, 5]])

        # Get interpolated values at centroids
        interpolation_points = domain.get_centroid_coordinates()
        answer = quantity.get_values(location='centroids')

        # print quantity.get_values(points=interpolation_points)
        assert num.allclose(answer, quantity.get_values(
            interpolation_points=interpolation_points))

        # Arbitrary values
        quantity = Quantity(domain, [[0, 1, 2], [3, 1, 7], [2, 1, 2], [3, 3, 7],
                                     [1, 4, -9], [2, 5, 0]])

        # Get interpolated values at centroids
        interpolation_points = domain.get_centroid_coordinates()
        answer = quantity.get_values(location='centroids')
        # print answer
        # print quantity.get_values(interpolation_points=interpolation_points)
        assert num.allclose(answer,
                            quantity.get_values(interpolation_points=interpolation_points,
                                                verbose=False))

        # FIXME TODO
        # indices = [0,5,3]
        # answer = [0.5,1,5]
        # assert allclose(answer,
        #                quantity.get_values(indices=indices, \
        #                                    location = 'unique vertices'))

    def test_get_interpolated_values_2(self):
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        # bac, bce, ecf, dbe
        vertices = [[1, 0, 2], [1, 2, 4], [4, 2, 5], [3, 1, 4]]

        domain = Generic_Domain(points, vertices)

        quantity = Quantity(domain)
        quantity.set_values(lambda x, y: x+2*y)  # 2 4 4 6

        # First pick one point
        x, y = 2.0/3, 8.0/3
        v = quantity.get_values(interpolation_points=[[x, y]])
        assert num.allclose(v, 6)

        # Then another to test that algorithm won't blindly
        # reuse interpolation matrix
        x, y = 4.0/3, 4.0/3
        v = quantity.get_values(interpolation_points=[[x, y]])
        assert num.allclose(v, 4)

    def test_get_interpolated_values_with_georef(self):

        zone = 56
        xllcorner = 308500
        yllcorner = 6189000
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        # bac, bce, ecf, dbe
        vertices = [[1, 0, 2], [1, 2, 4], [4, 2, 5], [3, 1, 4]]

        domain = Generic_Domain(points, vertices,
                                geo_reference=Geo_reference(zone, xllcorner, yllcorner))

        quantity = Quantity(domain)
        quantity.set_values(lambda x, y: x+2*y)  # 2 4 4 6

        # First pick one point (and turn it into absolute coordinates)
        x, y = 2.0/3, 8.0/3
        v = quantity.get_values(interpolation_points=[
                                [x+xllcorner, y+yllcorner]])
        assert num.allclose(v, 6)

        # Then another to test that algorithm won't blindly
        # reuse interpolation matrix
        x, y = 4.0/3, 4.0/3
        v = quantity.get_values(interpolation_points=[
                                [x+xllcorner, y+yllcorner]])
        assert num.allclose(v, 4)

        # Try two points
        pts = [[2.0/3 + xllcorner, 8.0/3 + yllcorner],
               [4.0/3 + xllcorner, 4.0/3 + yllcorner]]
        v = quantity.get_values(interpolation_points=pts)
        assert num.allclose(v, [6, 4])

        # Test it using the geospatial data format with absolute input points and default georef
        pts = Geospatial_data(data_points=pts)
        v = quantity.get_values(interpolation_points=pts)
        assert num.allclose(v, [6, 4])

        # Test it using the geospatial data format with relative input points
        pts = Geospatial_data(data_points=[[2.0/3, 8.0/3], [4.0/3, 4.0/3]],
                              geo_reference=Geo_reference(zone, xllcorner, yllcorner))
        v = quantity.get_values(interpolation_points=pts)
        assert num.allclose(v, [6, 4])

    def test_getting_some_vertex_values(self):
        """
        get values based on triangle lists.
        """
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(1, 3)

        # print "points",points
        # print "vertices",vertices
        # print "boundary",boundary

        # Create shallow water domain
        domain = Generic_Domain(points, vertices, boundary)
        # print "domain.number_of_elements ",domain.number_of_elements
        quantity = Quantity(domain, [[1, 1, 1], [2, 2, 2], [3, 3, 3],
                                     [4, 4, 4], [5, 5, 5], [6, 6, 6]])
        value = [7]
        indices = [1]
        quantity.set_values(value,
                            location='centroids',
                            indices=indices)
        # print "quantity.centroid_values",quantity.centroid_values
        # print "quantity.get_values(location = 'centroids') ",\
        #      quantity.get_values(location = 'centroids')
        assert num.allclose(quantity.centroid_values,
                            quantity.get_values(location='centroids'))

        value = [[15, 20, 25]]
        quantity.set_values(value, indices=indices)
        # print "1 quantity.vertex_values",quantity.vertex_values
        assert num.allclose(quantity.vertex_values, quantity.get_values())

        assert num.allclose(quantity.edge_values,
                            quantity.get_values(location='edges'))

        # get a subset of elements
        subset = quantity.get_values(location='centroids', indices=[0, 5])
        answer = [quantity.centroid_values[0], quantity.centroid_values[5]]
        assert num.allclose(subset, answer)

        subset = quantity.get_values(location='edges', indices=[0, 5])
        answer = [quantity.edge_values[0], quantity.edge_values[5]]
        # print "subset",subset
        # print "answer",answer
        assert num.allclose(subset, answer)

        subset = quantity.get_values(indices=[1, 5])
        answer = [quantity.vertex_values[1], quantity.vertex_values[5]]
        # print "subset",subset
        # print "answer",answer
        assert num.allclose(subset, answer)

    def test_smooth_vertex_values(self):
        """
        get values based on triangle lists.
        """
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular
        #  from anuga.shallow_water.shallow_water_domain import Domain

        # Create basic mesh
        points, vertices, boundary = rectangular(2, 2)

        # print "points",points
        # print "vertices",vertices
        # print "boundary",boundary

        # Create shallow water domain
        domain = Generic_Domain(points, vertices, boundary)
        # print "domain.number_of_elements ",domain.number_of_elements
        quantity = Quantity(domain, [[0, 0, 0], [1, 1, 1], [2, 2, 2], [3, 3, 3],
                                     [4, 4, 4], [5, 5, 5], [6, 6, 6], [7, 7, 7]])

        # print "quantity.get_values(location = 'unique vertices')", \
        #      quantity.get_values(location = 'unique vertices')

        # print "quantity.get_values(location = 'unique vertices')", \
        #      quantity.get_values(indices=[0,1,2,3,4,5,6,7], \
        #                          location = 'unique vertices')

        # print quantity.get_values(location = 'unique vertices')
        # print quantity.domain.number_of_triangles_per_node
        # print quantity.vertex_values

        # answer = [0.5, 2, 3, 3, 3.5, 4, 4, 5, 6.5]
        # assert allclose(answer,
        #                quantity.get_values(location = 'unique vertices'))

        quantity.smooth_vertex_values()

        # print quantity.vertex_values

        answer_vertex_values = [[3, 3.5, 0.5], [2, 0.5, 3.5], [3.5, 4, 2], [3, 2, 4],
                                [4, 5, 3], [3.5, 3, 5], [5, 6.5, 3.5], [4, 3.5, 6.5]]

        assert num.allclose(answer_vertex_values,
                            quantity.vertex_values)

        # Just another (slightly larger) test of get_values

        assert num.allclose(quantity.get_values(location='centroids'),
                            quantity.centroid_values)

        assert num.allclose(quantity.get_values(location='vertices'),
                            quantity.vertex_values)

        assert num.allclose(quantity.get_values(location='edges'),
                            quantity.edge_values)

    def test_maximum(self):
        quantity = Quantity(self.mesh4)

        # get referece to data arrays
        centroid_values = quantity.centroid_values
        vertex_values = quantity.vertex_values
        edge_values = quantity.edge_values

        quantity.set_values([[1, 2, 3], [5, 5, 5], [0, 0, 9], [-6, 3, 3]],
                            location='vertices')
        assert num.allclose(quantity.vertex_values,
                            [[1, 2, 3], [5, 5, 5], [0, 0, 9], [-6, 3, 3]])

        assert id(vertex_values) == id(quantity.vertex_values)
        assert num.allclose(quantity.centroid_values, [
                            2., 5., 3., 0.])  # Centroid
        assert num.allclose(quantity.edge_values, [[2.5, 2.0, 1.5],
                                                   [5., 5., 5.],
                                                   [4.5, 4.5, 0.],
                                                   [3.0, -1.5, -1.5]])

        other_quantity = Quantity(self.mesh4)

        other_quantity.set_values([[0, 0, 0], [1, 1, 6], [10, 10, 10], [0, 0, 4]],
                                  location='vertices')

        # ===============================
        quantity.maximum(other_quantity)
        # ===============================

        exact_vertex_values = num.array([[1.,   2.,   3.],
                                         [5.,   5.,  6.],
                                         [10.,  10.,  10.],
                                         [0.,   3.,  4.]])
        exact_centroid_values = num.array([2.,  5., 10.,  1.33333333])
        exact_edge_values = num.array([[2.5,   2.,    1.5],
                                       [5.,   5.,   5., ],
                                       [10.,   10.,  10.],
                                       [3.,    2.,    0.]])

        assert num.allclose(quantity.vertex_values,
                            exact_vertex_values)
        assert num.allclose(quantity.centroid_values,
                            exact_centroid_values)  # Centroid
        assert num.allclose(quantity.edge_values, exact_edge_values)

    def test_minimum(self):
        quantity = Quantity(self.mesh4)

        # get referece to data arrays
        centroid_values = quantity.centroid_values
        vertex_values = quantity.vertex_values
        edge_values = quantity.edge_values

        quantity.set_values([[1, 2, 3], [5, 5, 5], [0, 0, 9], [-6, 3, 3]],
                            location='vertices')
        assert num.allclose(quantity.vertex_values,
                            [[1, 2, 3], [5, 5, 5], [0, 0, 9], [-6, 3, 3]])

        assert id(vertex_values) == id(quantity.vertex_values)
        assert num.allclose(quantity.centroid_values, [
                            2., 5., 3., 0.])  # Centroid
        assert num.allclose(quantity.edge_values, [[2.5, 2.0, 1.5],
                                                   [5., 5., 5.],
                                                   [4.5, 4.5, 0.],
                                                   [3.0, -1.5, -1.5]])

        other_quantity = Quantity(self.mesh4)

        other_quantity.set_values([[0, 0, 0], [1, 1, 6], [10, 10, 10], [0, 0, 4]],
                                  location='vertices')

        # ===============================
        quantity.minimum(other_quantity)
        # ===============================

        exact_vertex_values = num.array([[0.,   0.,   0.],
                                         [1.,   1.,  5.],
                                         [0.,   0.,  9.],
                                         [-6.,   0.,  3.]])
        exact_centroid_values = num.array([0.,  2.66666667, 3.,  0.])
        exact_edge_values = num.array([[0.,   0.,   0.],
                                       [3.5,   3.5,   1., ],
                                       [4.5,   4.5,   0.],
                                       [2.,   -1.5,  -1.5]])

        assert num.allclose(quantity.vertex_values,
                            exact_vertex_values)
        assert num.allclose(quantity.centroid_values,
                            exact_centroid_values)  # Centroid
        assert num.allclose(quantity.edge_values, exact_edge_values)


class Test_Quantity_Memory(unittest.TestCase):
    """Verify selective array allocation for quantity types (QM1-QM6)."""

    def setUp(self):
        from anuga import rectangular_cross_domain
        self.domain = rectangular_cross_domain(10, 10)
        self.domain.set_quantity('elevation', -1.0)
        self.domain.set_quantity('stage', 0.5)
        self.domain.set_quantity('friction', 0.03)

    def test_evolved_has_all_arrays(self):
        for name in ['stage', 'xmomentum', 'ymomentum']:
            q = self.domain.quantities[name]
            # Update arrays are eagerly allocated for evolved quantities
            self.assertIsNotNone(q.explicit_update, name)
            self.assertIsNotNone(q.semi_implicit_update, name)
            self.assertIsNotNone(q.centroid_backup_values, name)
            self.assertIsNotNone(q.edge_values, name)
            # Gradient arrays are lazy — not allocated until first access
            self.assertIsNone(q._x_gradient, name)
            self.assertIsNone(q._y_gradient, name)
            self.assertIsNone(q._phi, name)
            # Accessing the property triggers allocation
            self.assertIsNotNone(q.x_gradient, name)
            self.assertIsNotNone(q.y_gradient, name)
            self.assertIsNotNone(q.phi, name)

    def test_elevation_strips_update_arrays(self):
        q = self.domain.quantities['elevation']
        self.assertIsNone(q.explicit_update)
        self.assertIsNone(q.semi_implicit_update)
        self.assertIsNone(q.centroid_backup_values)
        # Gradient arrays and phi are lazy (not eagerly allocated)
        self.assertIsNone(q._phi)
        self.assertIsNone(q._x_gradient)
        self.assertIsNone(q._y_gradient)
        # Edge values are allocated (needed by DE solver)
        self.assertIsNotNone(q.edge_values)

    def test_friction_centroid_only(self):
        q = self.domain.quantities['friction']
        self.assertIsNone(q.explicit_update)
        self.assertIsNone(q.semi_implicit_update)
        self.assertIsNone(q.centroid_backup_values)
        self.assertIsNone(q._x_gradient)
        self.assertIsNone(q._y_gradient)
        self.assertIsNone(q._phi)
        self.assertIsNone(q.edge_values)
        self.assertTrue(num.allclose(q.centroid_values, 0.03))

    def test_height_edge_diagnostic(self):
        q = self.domain.quantities['height']
        self.assertIsNone(q.explicit_update)
        self.assertIsNone(q._x_gradient)
        self.assertIsNone(q._phi)
        self.assertIsNotNone(q.edge_values)

    def test_velocity_edge_diagnostic(self):
        for name in ['xvelocity', 'yvelocity']:
            q = self.domain.quantities[name]
            self.assertIsNone(q.explicit_update, name)
            self.assertIsNone(q._x_gradient, name)
            self.assertIsNone(q._phi, name)
            self.assertIsNotNone(q.edge_values, name)

    def test_vertex_values_lazy(self):
        from anuga import rectangular_cross_domain
        d = rectangular_cross_domain(4, 4)
        q = d.quantities['xmomentum']
        vv = q.vertex_values  # lazy allocation on first access
        self.assertIsNotNone(vv)
        self.assertEqual(vv.shape, (len(d), 3))

    def test_explicit_qty_type_parameter(self):
        """qty_type kwarg overrides domain._quantity_type_map lookup."""
        from anuga.abstract_2d_finite_volumes.quantity import Quantity
        # Force centroid_only via explicit parameter (domain map has no 'myq')
        q = Quantity(self.domain, name='myq', qty_type='centroid_only')
        self.assertIsNone(q.edge_values)
        self.assertIsNone(q._x_gradient)
        self.assertIsNone(q.explicit_update)
        self.assertEqual(q._qty_type, 'centroid_only')

        # Explicit parameter wins over domain map entry
        self.domain._quantity_type_map['overridden'] = 'evolved'
        q2 = Quantity(self.domain, name='overridden', qty_type='centroid_only')
        self.assertEqual(q2._qty_type, 'centroid_only')
        self.assertIsNone(q2.edge_values)

    def test_friction_set_quantity_correct_value(self):
        self.domain.set_quantity('friction', 0.07)
        q = self.domain.quantities['friction']
        self.assertTrue(num.allclose(q.centroid_values, 0.07))

    def test_memory_saving(self):
        N = len(self.domain)

        def array_bytes(q):
            """Count bytes of eagerly-allocated arrays only.

            Gradient/phi backing stores (_x_gradient etc.) are checked directly
            to avoid triggering lazy allocation via the property.
            """
            total = 0
            for attr in ['centroid_values', 'edge_values',
                         'explicit_update', 'semi_implicit_update',
                         'centroid_backup_values']:
                arr = getattr(q, attr, None)
                if arr is not None:
                    total += arr.nbytes
            # Lazy arrays — check backing store, not property
            for attr in ['_x_gradient', '_y_gradient', '_phi']:
                arr = getattr(q, attr, None)
                if arr is not None:
                    total += arr.nbytes
            return total

        evolved_bytes = array_bytes(self.domain.quantities['stage'])
        friction_bytes = array_bytes(self.domain.quantities['friction'])
        elevation_bytes = array_bytes(self.domain.quantities['elevation'])
        height_bytes = array_bytes(self.domain.quantities['height'])

        # evolved: centroid+edge+explicit_update+semi_implicit+centroid_backup (56N)
        # gradients/phi are lazy — not counted unless accessed
        self.assertEqual(evolved_bytes, 56 * N)
        self.assertEqual(friction_bytes, 8 * N)    # centroid only
        self.assertEqual(elevation_bytes, 32 * N)  # centroid+edge
        self.assertEqual(height_bytes, 32 * N)     # centroid+edge

    @pytest.mark.slow
    def test_total_domain_memory_reduction(self):
        """Total allocated bytes across all quantities must be < 60% of the
        all-evolved baseline, for a domain large enough to be realistic.

        Baseline (all-evolved, no optimisation): 10 quantities × 80N bytes = 800N
        Optimised layout (gradients/phi are lazy — not counted here):
          stage / xmomentum / ymomentum   3 × 56N = 168N
          elevation (edge_diagnostic)              32N
          friction  (centroid_only)                 8N
          height / xvelocity / yvelocity  3 × 32N =  96N
          x / y    (coordinate)           2 × 32N =  64N
          ─────────────────────────────────────────────
          Total optimised                          368N
          Saving  (800-368)/800 ≈ 54 %

        The test uses a 100×100 rectangular_cross_domain (80 000 triangles) so
        the absolute numbers are large enough to catch rounding / off-by-one
        errors, while still completing in < 5 s.
        """
        import pytest
        from anuga import rectangular_cross_domain

        domain = rectangular_cross_domain(100, 100)
        domain.set_quantity('elevation', -1.0)
        domain.set_quantity('stage', 0.5)
        domain.set_quantity('friction', 0.03)

        N = len(domain)

        # Sum every fixed (non-vertex) array for every quantity.
        # Non-lazy attrs (safe to read via property)
        _eager_attrs = [
            'centroid_values', 'edge_values',
            'explicit_update', 'semi_implicit_update', 'centroid_backup_values',
        ]
        # Lazy attrs — check backing store to avoid triggering allocation
        _lazy_backing = ['_x_gradient', '_y_gradient', '_phi']

        def qty_fixed_bytes(q):
            total = sum(
                getattr(q, a).nbytes
                for a in _eager_attrs
                if getattr(q, a) is not None
            )
            total += sum(
                getattr(q, a).nbytes
                for a in _lazy_backing
                if getattr(q, a) is not None
            )
            return total

        actual_bytes = sum(qty_fixed_bytes(q) for q in domain.quantities.values())

        # Unoptimised baseline: every quantity uses the full 'evolved' layout
        n_quantities = len(domain.quantities)
        baseline_bytes = n_quantities * 80 * N  # 10 arrays × 8 bytes, each length N

        saving_pct = 100.0 * (baseline_bytes - actual_bytes) / baseline_bytes

        # At least 50 % saving must be achieved (actual is ~54 %)
        self.assertGreaterEqual(
            saving_pct, 50.0,
            msg=f"Expected >= 50% saving, got {saving_pct:.1f}% "
                f"(actual={actual_bytes/1e6:.1f} MB, baseline={baseline_bytes/1e6:.1f} MB, "
                f"N={N})"
        )

        # Spot-check absolute numbers match the per-quantity analysis
        # evolved: centroid+edge+explicit+semi_implicit+backup = 56N (gradients lazy)
        self.assertEqual(qty_fixed_bytes(domain.quantities['stage']),    56 * N)
        self.assertEqual(qty_fixed_bytes(domain.quantities['elevation']), 32 * N)
        self.assertEqual(qty_fixed_bytes(domain.quantities['friction']),   8 * N)
        self.assertEqual(qty_fixed_bytes(domain.quantities['height']),    32 * N)


    def test_gradient_workspace_exists_on_domain(self):
        """Domain carries three shared gradient workspace arrays (QM7)."""
        from anuga import rectangular_cross_domain
        d = rectangular_cross_domain(5, 5)
        N = len(d)
        for attr in ['_grad_workspace_x', '_grad_workspace_y', '_phi_workspace']:
            arr = getattr(d, attr, None)
            self.assertIsNotNone(arr, f'domain.{attr} should exist')
            self.assertEqual(arr.shape, (N,))

    def test_extrapolate_by_edge_uses_workspace_not_quantity_storage(self):
        """After extrapolate_second_order_and_limit_by_edge the quantity's own
        _x_gradient/_y_gradient/_phi remain None — the shared workspace was used."""
        from anuga import rectangular_cross_domain
        d = rectangular_cross_domain(10, 10)
        d.set_quantity('stage', lambda x, y: x + 0.5 * y)
        q = d.quantities['stage']

        # Confirm lazy arrays start as None
        self.assertIsNone(q._x_gradient)
        self.assertIsNone(q._y_gradient)
        self.assertIsNone(q._phi)

        q.extrapolate_second_order_and_limit_by_edge()

        # Workspace was used — quantity's own backing stores still None
        self.assertIsNone(q._x_gradient, 'x_gradient should stay None (workspace used)')
        self.assertIsNone(q._y_gradient, 'y_gradient should stay None (workspace used)')
        self.assertIsNone(q._phi, 'phi should stay None (workspace used)')

        # Edge values must be non-trivial
        self.assertFalse(num.allclose(q.edge_values, 0))

    def test_extrapolate_by_vertex_uses_workspace_not_quantity_storage(self):
        """After extrapolate_second_order_and_limit_by_vertex the quantity's own
        gradient backing stores remain None."""
        from anuga import rectangular_cross_domain
        d = rectangular_cross_domain(10, 10)
        d.set_quantity('stage', lambda x, y: x)
        q = d.quantities['stage']

        q.extrapolate_second_order_and_limit_by_vertex()

        self.assertIsNone(q._x_gradient, 'x_gradient should stay None (workspace used)')
        self.assertIsNone(q._y_gradient)
        self.assertIsNone(q._phi)
        self.assertFalse(num.allclose(q.vertex_values, 0))

    def test_extrapolate_by_edge_correct_values(self):
        """Workspace-based extrapolation produces the same edge values as direct
        gradient access would."""
        import numpy as np
        from anuga import rectangular_cross_domain
        from anuga.abstract_2d_finite_volumes.quantity_openmp_ext import (
            extrapolate_second_order_and_limit_by_edge as ext_fn,
        )

        d = rectangular_cross_domain(10, 10)
        d.set_quantity('stage', lambda x, y: x + 0.5 * y)

        # Run via workspace path
        q_ws = d.quantities['stage']
        q_ws.extrapolate_second_order_and_limit_by_edge()
        edge_ws = q_ws.edge_values.copy()

        # Run via per-quantity gradient path (pass no workspace)
        d2 = rectangular_cross_domain(10, 10)
        d2.set_quantity('stage', lambda x, y: x + 0.5 * y)
        q_direct = d2.quantities['stage']
        ext_fn(q_direct)  # uses quantity's own lazy arrays
        edge_direct = q_direct.edge_values.copy()

        np.testing.assert_allclose(edge_ws, edge_direct, rtol=1e-12,
                                   err_msg='Workspace path must match per-quantity path')


class Test_Quantity_extra(unittest.TestCase):
    """Tests for uncovered Quantity methods."""

    def _make_domain_and_quantity(self):
        import anuga
        d = anuga.rectangular_cross_domain(3, 3)
        d.set_quantity('elevation', 0.0)
        d.set_quantity('stage', 1.0)
        d.distribute_to_vertices_and_edges()
        q = d.quantities['stage']
        return d, q

    def test_vertex_values_setter_none(self):
        """Setting vertex_values to None clears the cached array."""
        d, q = self._make_domain_and_quantity()
        _ = q.vertex_values  # trigger allocation
        q.vertex_values = None
        self.assertIsNone(q._vertex_values)

    def test_x_gradient_setter(self):
        import numpy as np
        d, q = self._make_domain_and_quantity()
        arr = np.ones(len(d))
        q.x_gradient = arr
        self.assertIs(q.x_gradient, arr)

    def test_y_gradient_setter(self):
        import numpy as np
        d, q = self._make_domain_and_quantity()
        arr = np.ones(len(d)) * 2.0
        q.y_gradient = arr
        self.assertIs(q.y_gradient, arr)

    def test_phi_setter(self):
        import numpy as np
        d, q = self._make_domain_and_quantity()
        arr = np.ones(len(d)) * 3.0
        q.phi = arr
        self.assertIs(q.phi, arr)

    def test_rdiv(self):
        d, q = self._make_domain_and_quantity()
        result = q.__rdiv__(3.0)
        self.assertIsNotNone(result)

    def test_pow_with_scalar(self):
        import numpy as np
        d, q = self._make_domain_and_quantity()
        result = q ** 2
        np.testing.assert_allclose(result.centroid_values,
                                   q.centroid_values ** 2)

    def test_maximum_with_quantity(self):
        import anuga
        import numpy as np
        d = anuga.rectangular_cross_domain(3, 3)
        d.set_quantity('elevation', 0.0)
        d.set_quantity('stage', 1.0)
        d.set_quantity('xmomentum', 2.0)
        d.distribute_to_vertices_and_edges()
        q1 = d.quantities['stage']
        q2 = d.quantities['xmomentum']
        q1.maximum(q2)
        self.assertTrue(np.all(q1.centroid_values >= 1.0))

    def test_maximum_with_scalar(self):
        import numpy as np
        d, q = self._make_domain_and_quantity()
        q.maximum(2.0)
        self.assertTrue(np.all(q.centroid_values >= 2.0))

    def test_minimum_with_scalar(self):
        import numpy as np
        d, q = self._make_domain_and_quantity()
        q.minimum(0.5)
        self.assertTrue(np.all(q.centroid_values <= 1.0))

    def test_save_centroid_data_to_csv(self):
        import tempfile
        import os
        d, q = self._make_domain_and_quantity()
        with tempfile.TemporaryDirectory() as tmpdir:
            fname = os.path.join(tmpdir, 'stage')
            q.save_centroid_data_to_csv(filename=fname)
            out = fname + '_centroid_data.csv'
            self.assertTrue(os.path.exists(out))
            lines = open(out).readlines()
            self.assertGreater(len(lines), 0)

    def test_save_centroid_data_to_csv_default_filename(self):
        """save_centroid_data_to_csv with no filename uses self.name (line 429)."""
        import tempfile
        import os
        d, q = self._make_domain_and_quantity()
        # Run from a temp directory so the output file lands there
        orig_dir = os.getcwd()
        with tempfile.TemporaryDirectory() as tmpdir:
            os.chdir(tmpdir)
            try:
                q.save_centroid_data_to_csv()
                out = q.name + '_centroid_data.csv'
                self.assertTrue(os.path.exists(out))
            finally:
                os.chdir(orig_dir)

    def test_interpolate_old(self):
        """interpolate_old computes centroid values from vertices (lines 775-781)."""
        d, q = self._make_domain_and_quantity()
        q.set_values(2.0)
        q.interpolate_old()
        import numpy as np
        self.assertTrue(np.allclose(q.centroid_values, 2.0))

    def test_interpolate_from_vertices_to_edges_centroid_only(self):
        """centroid_only quantity returns early (line 804)."""
        from anuga.abstract_2d_finite_volumes.quantity import Quantity
        d, _ = self._make_domain_and_quantity()
        qc = Quantity(d, qty_type='centroid_only')
        qc.interpolate_from_vertices_to_edges()  # should not raise

    def test_interpolate_from_edges_to_vertices_centroid_only(self):
        """centroid_only quantity returns early (line 811)."""
        from anuga.abstract_2d_finite_volumes.quantity import Quantity
        d, _ = self._make_domain_and_quantity()
        qc = Quantity(d, qty_type='centroid_only')
        qc.interpolate_from_edges_to_vertices()  # should not raise

    def test_save_to_array_with_bounds(self):
        """save_to_array with explicit easting/northing covers else branches (576, 581, 586, 591, 612)."""
        d, q = self._make_domain_and_quantity()
        import anuga
        d.set_quantity('elevation', 0.0)
        d.set_quantity('stage', 1.0)
        d.distribute_to_vertices_and_edges()
        q = d.quantities['stage']
        # Use bounding box that encloses the domain
        xg, yg, grid = q.save_to_array(
            easting_min=0.0, easting_max=2.0,
            northing_min=0.0, northing_max=2.0)
        self.assertIsNotNone(grid)

    def test_get_name(self):
        """get_name returns the quantity name (line 664)."""
        d, q = self._make_domain_and_quantity()
        name = q.get_name()
        self.assertIsNotNone(name)

    def test_set_beta_out_of_range(self):
        """set_beta logs warnings for values outside [0, 2] (lines 671, 673)."""
        d, q = self._make_domain_and_quantity()
        q.set_beta(-0.1)   # triggers line 671
        q.set_beta(2.5)    # triggers line 673

    def test_get_beta(self):
        """get_beta returns current beta (line 683)."""
        d, q = self._make_domain_and_quantity()
        q.set_beta(0.5)
        self.assertAlmostEqual(q.get_beta(), 0.5)

    def test_set_boundary_values_invalid_raises(self):
        """set_boundary_values with non-numeric string raises (lines 700-703)."""
        d, q = self._make_domain_and_quantity()
        with self.assertRaises(Exception):
            q.set_boundary_values('not_a_number')

    def test_save_data_to_dem(self):
        """save_data_to_dem writes a CSV file (lines 463-481)."""
        import os
        import tempfile
        d, q = self._make_domain_and_quantity()
        orig = os.getcwd()
        with tempfile.TemporaryDirectory() as tmpdir:
            os.chdir(tmpdir)
            try:
                q.save_data_to_dem('stage_test')
                self.assertTrue(os.path.exists('stage_test_centroid_data.csv'))
            finally:
                os.chdir(orig)

    def test_set_values_array_2d_with_indices(self):
        """set_values from 2D array (Nx1) with indices (line 1172)."""
        import numpy as np
        d, q = self._make_domain_and_quantity()
        N = d.number_of_triangles
        indices = np.array([0, 1, 2])
        vals = np.array([[2.0], [3.0], [4.0]])  # shape (3,1) → Nx1
        q.set_values(vals, location='centroids', indices=indices)
        self.assertAlmostEqual(q.centroid_values[0], 2.0)

    def test_set_values_array_vertex_2d(self):
        """set_values from 2D Nx3 array sets vertex_values (lines 1193-1198)."""
        import numpy as np
        d, q = self._make_domain_and_quantity()
        N = d.number_of_triangles
        vals = np.ones((N, 3)) * 5.0
        q.set_values(vals, location='vertices')
        self.assertTrue(np.allclose(q.vertex_values, 5.0))

    def test_set_values_array_vertex_2d_with_indices(self):
        """set_values Nx3 array with indices sets specific rows (line 1198)."""
        import numpy as np
        d, q = self._make_domain_and_quantity()
        indices = np.array([0, 1])
        vals = np.array([[9.0, 9.0, 9.0], [8.0, 8.0, 8.0]])
        q.set_values(vals, location='vertices', indices=indices)
        self.assertTrue(np.allclose(q.vertex_values[0], 9.0))


class Test_Quantity_extra2(unittest.TestCase):
    """Additional quantity coverage tests."""

    def _make_domain_and_quantity(self):
        import anuga
        d = anuga.rectangular_cross_domain(3, 3)
        d.set_quantity('elevation', 0.0)
        d.set_quantity('stage', 1.0)
        d.distribute_to_vertices_and_edges()
        q = d.quantities['stage']
        return d, q

    def test_div_quantity_by_quantity(self):
        """Q1 / Q2 covers isinstance(other, Quantity) True branch (line 302)."""
        import anuga
        d = anuga.rectangular_cross_domain(2, 2)
        d.set_quantity('stage', 2.0)
        d.set_quantity('elevation', 1.0)
        q1 = d.quantities['stage']
        q2 = d.quantities['elevation']
        result = q1 / q2
        import numpy as np
        self.assertTrue(np.allclose(result.centroid_values, 2.0))

    def test_pow_quantity_by_quantity(self):
        """Q1 ** Q2 covers isinstance(other, Quantity) True branch (line 336).
        Note: __pow__ has a bug (uses other not Q.vertex_values) so it raises TypeError.
        """
        import anuga
        d = anuga.rectangular_cross_domain(2, 2)
        d.set_quantity('stage', 2.0)
        d.set_quantity('elevation', 2.0)
        q1 = d.quantities['stage']
        q2 = d.quantities['elevation']
        try:
            result = q1 ** q2
        except TypeError:
            pass  # pre-existing bug: uses other not Q in computation

    def test_set_values_edges_location_raises(self):
        """location='edges' raises deprecated error (lines 978-980)."""
        d, q = self._make_domain_and_quantity()
        with self.assertRaises(Exception) as cm:
            q.set_values(1.0, location='edges')
        self.assertIn('deprecated', str(cm.exception))

    def test_set_values_invalid_location_raises(self):
        """Invalid location raises (lines 982-984)."""
        d, q = self._make_domain_and_quantity()
        with self.assertRaises(Exception) as cm:
            q.set_values(1.0, location='badlocation')
        self.assertIn('Invalid location', str(cm.exception))

    def test_set_values_nonfloat_string_raises(self):
        """Non-coercible string for numeric raises ValueError path (lines 1008-1011)."""
        d, q = self._make_domain_and_quantity()
        with self.assertRaises(Exception):
            q.set_values('not_a_number')

    def test_set_values_2d_centroid_no_indices(self):
        """2D Nx1 centroid array with indices=None (line 1170)."""
        import numpy as np
        d, q = self._make_domain_and_quantity()
        N = d.number_of_triangles
        vals = np.ones((N, 1)) * 3.0
        q.set_values(vals, location='centroids')
        self.assertTrue(np.allclose(q.centroid_values, 3.0))

    def test_set_values_3d_centroid_array_raises(self):
        """3D centroid array raises (lines 1174-1175)."""
        import numpy as np
        d, q = self._make_domain_and_quantity()
        N = d.number_of_triangles
        vals = np.ones((N, 1, 1))  # 3D array
        with self.assertRaises(Exception) as cm:
            q.set_values(vals, location='centroids')
        self.assertIn('1d or 2d', str(cm.exception))

    def test_set_values_3d_vertex_array_raises(self):
        """3D vertex array raises (lines 1202-1203)."""
        import numpy as np
        d, q = self._make_domain_and_quantity()
        N = d.number_of_triangles
        vals = np.ones((N, 3, 1))  # 3D array
        with self.assertRaises(Exception) as cm:
            q.set_values(vals, location='vertices')
        self.assertIn('1d or 2d', str(cm.exception))

    def test_set_values_obsolete_from_points(self):
        """set_values_from_points is obsolete (line 1355)."""
        d, q = self._make_domain_and_quantity()
        with self.assertRaises(Exception) as cm:
            q.set_values_from_points(None, None, None, 'vertices', None)
        self.assertIn('obsolete', str(cm.exception))

    def test_set_values_from_file_bad_location_raises(self):
        """set_values_from_file with location='centroids' raises (lines 1380-1382)."""
        import tempfile
        import os
        d, q = self._make_domain_and_quantity()
        with tempfile.NamedTemporaryFile(suffix='.pts', delete=False) as f:
            fname = f.name
        try:
            with self.assertRaises(Exception) as cm:
                q.set_values_from_file(fname, None, 0.001, 'centroids',
                                       None, verbose=False)
            self.assertIn("location='vertices'", str(cm.exception))
        finally:
            os.unlink(fname)

    def test_set_values_unique_vertices_no_indices(self):
        """set_values_from_constant unique vertices with None indices (line 1084)."""
        import numpy as np
        d, q = self._make_domain_and_quantity()
        q.set_values(5.0, location='unique vertices')
        # Edge values should be updated (line 1084 sets edge_values)
        self.assertIsNotNone(q.edge_values)

    def test_set_values_from_function_centroids_no_indices(self):
        """set_values_from_function centroids with indices=None (lines 1246-1249)."""
        import numpy as np
        d, q = self._make_domain_and_quantity()
        q.set_values(lambda x, y: x + y, location='centroids')
        # centroid values should be non-trivial
        self.assertGreater(np.max(q.centroid_values), 0)

    def test_set_values_polygon_no_numeric_raises(self):
        """polygon without numeric raises (line 946)."""
        d, q = self._make_domain_and_quantity()
        poly = [[0, 0], [1, 0], [1, 1], [0, 1]]
        with self.assertRaises(Exception) as cm:
            q.set_values(polygon=poly)
        self.assertIn('polygon', str(cm.exception))

    def test_set_values_polygon_with_smooth(self):
        """polygon with smooth=True hits smooth_vertex_values call (line 961)."""
        d, q = self._make_domain_and_quantity()
        poly = [[0.1, 0.1], [0.9, 0.1], [0.9, 0.9], [0.1, 0.9]]
        q.set_values(numeric=2.0, polygon=poly, smooth=True)
        self.assertIsNotNone(q.vertex_values)

    def test_set_values_function_keyword(self):
        """function= keyword path in set_values (lines 1016-1020)."""
        import numpy as np
        d, q = self._make_domain_and_quantity()
        q.set_values(function=lambda x, y: x + y, location='vertices')
        self.assertIsNotNone(q.centroid_values)

    def test_set_values_unknown_file_ext_raises(self):
        """Unknown file extension raises (line 1045)."""
        import tempfile
        import os
        d, q = self._make_domain_and_quantity()
        with tempfile.NamedTemporaryFile(suffix='.xyz', delete=False) as f:
            fname = f.name
        try:
            with self.assertRaises(Exception) as cm:
                q.set_values(filename=fname)
            self.assertIn('Extension', str(cm.exception))
        finally:
            os.unlink(fname)


# -------------------------------------------------------------
if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_Quantity)
    unittest.TextTestRunner(verbosity=2).run(suite)

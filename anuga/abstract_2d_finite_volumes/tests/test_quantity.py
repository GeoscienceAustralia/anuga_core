#!/usr/bin/env python

import unittest
from math import sqrt, pi
import tempfile

from anuga.abstract_2d_finite_volumes.quantity import *
from anuga.config import epsilon

from anuga.fit_interpolate.fit import fit_to_mesh
#from anuga.pyvolution.least_squares import fit_to_mesh         
from anuga.abstract_2d_finite_volumes.generic_domain \
                import Generic_Domain
from anuga.geospatial_data.geospatial_data import Geospatial_data
from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.geometry.polygon import *

import numpy as num


#Aux for fit_interpolate.fit example
def linear_function(point):
    point = num.array(point)
    return point[:,0]+3*point[:,1]
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

        #bac, bce, ecf, dbe
        elements = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4] ]

        self.mesh1 = Generic_Domain(points[:3], [elements[0]])
        self.mesh1.check_integrity()
        
        #print self.mesh1.__class__
        #print isinstance(self.mesh1, Domain)

        self.mesh4 = Generic_Domain(points, elements)
        self.mesh4.check_integrity()

        # UTM round Onslow
        a = [240000, 7620000]
        b = [240000, 7680000]
        c = [300000, 7620000]

        points = [a, b, c]
        elements = [[0,2,1]]
        
        self.mesh_onslow = Generic_Domain(points, elements)
        self.mesh_onslow.check_integrity()
        
    def tearDown(self):
        pass
        #print "  Tearing down"


    def test_creation(self):

        quantity = Quantity(self.mesh1, [[1,2,3]])
        assert num.allclose(quantity.vertex_values, [[1.,2.,3.]])

        try:
            quantity = Quantity()
        except:
            pass
        else:
            raise Exception('Should have raised empty quantity exception')


        # FIXME(Ole): Temporarily disabled 18 Jan 2009
        #try:
        #    quantity = Quantity([1,2,3])
        #except AssertionError:
        #    pass
        #except:
        #    raise Exception('Should have raised "mising mesh object" error')


    def test_creation_zeros(self):

        quantity = Quantity(self.mesh1)
        assert num.allclose(quantity.vertex_values, [[0.,0.,0.]])


        quantity = Quantity(self.mesh4)
        assert num.allclose(quantity.vertex_values, [[0.,0.,0.], [0.,0.,0.],
                                                     [0.,0.,0.], [0.,0.,0.]])

    def test_set_boundary_values(self):

        quantity = Quantity(self.mesh1)

        quantity.set_boundary_values()

        assert  num.allclose(quantity.boundary_values, [0.0, 0.0, 0.0])


    def test_set_boundary_values_with_function(self):

        quantity = Quantity(self.mesh1)
        #assert num.allclose(quantity.vertex_values, [[0.,0.,0.]])

        def simple(x,y):
            return x+3*y

        quantity.set_boundary_values(simple)

        assert  num.allclose(quantity.boundary_values, [1.0, 4.0, 3.0])

    def test_set_boundary_values_with_constant(self):

        quantity = Quantity(self.mesh1)
        #assert num.allclose(quantity.vertex_values, [[0.,0.,0.]])

        quantity.set_boundary_values(10.0)

        assert  num.allclose(quantity.boundary_values, [10.0, 10.0, 10.0])


    def test_set_boundary_values_with_array(self):

        quantity = Quantity(self.mesh1)
        #assert num.allclose(quantity.vertex_values, [[0.,0.,0.]])


        quantity.set_boundary_values([10.0, 4.0, 5.0])

        assert  num.allclose(quantity.boundary_values, [10.0, 4.0, 5.0])

    def test_set_boundary_values_with_wrong_sized_array(self):

        quantity = Quantity(self.mesh1)
        #assert num.allclose(quantity.vertex_values, [[0.,0.,0.]])

        try:
            quantity.set_boundary_values([10.0, 4.0, 5.0, 8.0])
        except:
            pass
        else:
            msg = 'Should have caught this'
            raise Exception(msg)

    def test_set_boundary_values_from_edges(self):

        quantity = Quantity(self.mesh4)

        def simple(x,y):
            return x+3*y

        quantity.set_values(simple)

        assert  num.allclose(quantity.boundary_values, [  0.,   0.,   0.,   0.,  0.,   0.])

        quantity.set_boundary_values_from_edges()

        assert  num.allclose(quantity.boundary_values, [  1.,   3.,   3.,   6.,  10.,   9.])

        
    def test_interpolation(self):
        quantity = Quantity(self.mesh1, [[1,2,3]])
        assert num.allclose(quantity.centroid_values, [2.0]) #Centroid

        assert num.allclose(quantity.edge_values, [[2.5, 2.0, 1.5]])


    def test_interpolation2(self):
        quantity = Quantity(self.mesh4,
                            [[1,2,3], [5,5,5], [0,0,9], [-6, 3, 3]])
        assert num.allclose(quantity.centroid_values, [2., 5., 3., 0.]) #Centroid


        quantity.extrapolate_second_order()

        #print quantity.vertex_values
        assert num.allclose(quantity.vertex_values, [[3.5, -1.0, 3.5],
                                                     [3.+2./3, 6.+2./3, 4.+2./3],
                                                 [4.6, 3.4, 1.],
                                                 [-5.0, 1.0, 4.0]])

        #print quantity.edge_values
        assert num.allclose(quantity.edge_values, [[1.25, 3.5, 1.25],
                                                   [5. + 2/3.0, 4.0 + 1.0/6, 5.0 + 1.0/6],
                                                   [2.2, 2.8, 4.0],
                                                   [2.5, -0.5, -2.0]])
        
    def test_save_to_array(self):
        quantity = Quantity(self.mesh4,
                            [[1,2,3], [5,5,5], [0,0,9], [-6, 3, 3]])
        
        
        assert num.allclose(quantity.centroid_values, [2., 5., 3., 0.]) #Centroid

        cellsize = 1.0
        x,y,z = quantity.save_to_array(cellsize=cellsize, smooth=False)
        #x,y,z = quantity.save_to_array(smooth=False)
        
        from pprint import pprint
        #pprint(x)
        #pprint(y)
        #pprint(z)
        
        
        x_ex = [ 0.,  1.,  2.,  3.,  4.]
        y_ex = [ 0.,  1.,  2.,  3.,  4.]
        
        z_ex = [[  2.00000000e+00,   2.50000000e+00,   0.00000000e+00,
          4.50000000e+00,   9.00000000e+00],
       [  1.50000000e+00,   5.00000000e+00,   0.00000000e+00,
          4.50000000e+00,  -9.99900000e+03],
       [  3.00000000e+00,   3.00000000e+00,   3.00000000e+00,
         -9.99900000e+03,  -9.99900000e+03],
       [ -1.50000000e+00,  -1.50000000e+00,  -9.99900000e+03,
         -9.99900000e+03,  -9.99900000e+03],
       [ -6.00000000e+00,  -9.99900000e+03,  -9.99900000e+03,
         -9.99900000e+03,  -9.99900000e+03]]

        assert num.allclose(x_ex,x)
        assert num.allclose(y_ex,y)
        assert num.allclose(z_ex,z)
                
               
        Plot = False 
        if Plot: 
            import pylab
            import numpy
            #a = numpy.where(a == -9999, numpy.nan, a)
            #a = numpy.where(a > 10.0, numpy.nan, a)
        
            #z = z[::-1,:]

            
            print z
            print z.shape
            print x
            print y

            nrows = z.shape[0]
            ncols = z.shape[1]            
        
            ratio = float(nrows)/float(ncols)
            print ratio
            
            #y = numpy.arange(nrows)*cellsize
            #x = numpy.arange(ncols)*cellsize
        
            #Setup fig size to correpond to array size
            fig = pylab.figure(figsize=(10, 10*ratio))
        
            levels = numpy.arange(-7, 10, 0.1)
            CF = pylab.contourf(x,y,z, levels=levels)
            CB = pylab.colorbar(CF, shrink=0.8, extend='both')
            #CC = pylab.contour(x,y,a, levels=levels)
            
            pylab.show()
        


        x,y,z = quantity.save_to_array(cellsize=cellsize, smooth=True)
        
        
        x_ex = [ 0.,  1.,  2.,  3.,  4.]
        y_ex = [ 0.,  1.,  2.,  3.,  4.]
        z_ex = [[  2.00000000e+00,   2.33333333e+00,   2.66666667e+00,
          5.83333333e+00,   9.00000000e+00],
       [  2.50000000e+00,   2.83333333e+00,   2.66666667e+00,
          5.83333333e+00,  -9.99900000e+03],
       [  3.00000000e+00,   2.83333333e+00,   2.66666667e+00,
         -9.99900000e+03,  -9.99900000e+03],
       [ -1.50000000e+00,  -1.66666667e+00,  -9.99900000e+03,
         -9.99900000e+03,  -9.99900000e+03],
       [ -6.00000000e+00,  -9.99900000e+03,  -9.99900000e+03,
         -9.99900000e+03,  -9.99900000e+03]]

        #pprint(z)
        assert num.allclose(x_ex,x)
        assert num.allclose(y_ex,y)
        assert num.allclose(z_ex,z)
        
        if Plot:
            import pylab
            import numpy
            #a = numpy.where(a == -9999, numpy.nan, a)
            #a = numpy.where(a > 10.0, numpy.nan, a)
        
            #a = a[::-1,:]
            nrows = z.shape[0]
            ncols = z.shape[1]
        
            ratio = float(nrows)/float(ncols)
            print ratio
        
            #Setup fig size to correpond to array size
            fig = pylab.figure(figsize=(10, 10*ratio))
        
            levels = numpy.arange(-7, 10, 0.1)
            CF = pylab.contourf(x,y,z, levels=levels)
            CB = pylab.colorbar(CF, shrink=0.8, extend='both')
            #CC = pylab.contour(x,y,a, levels=[0.0,1.0,2.0,3.0])
            
            pylab.show()        
        




    def test_get_extrema_1(self):
        quantity = Quantity(self.mesh4,
                                      [[1,2,3], [5,5,5], [0,0,9], [-6, 3, 3]])
        assert num.allclose(quantity.centroid_values, [2., 5., 3., 0.]) #Centroids

        v = quantity.get_maximum_value()
        assert v == 5

        v = quantity.get_minimum_value()
        assert v == 0        

        i = quantity.get_maximum_index()
        assert i == 1

        i = quantity.get_minimum_index()
        assert i == 3        
        
        x,y = quantity.get_maximum_location()
        xref, yref = 4.0/3, 4.0/3
        assert x == xref
        assert y == yref

        v = quantity.get_values(interpolation_points = [[x,y]])
        assert num.allclose(v, 5)


        x,y = quantity.get_minimum_location()
        v = quantity.get_values(interpolation_points = [[x,y]])
        assert num.allclose(v, 0)


    def test_get_maximum_2(self):

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [[1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Generic_Domain(points, vertices)

        quantity = Quantity(domain)
        quantity.set_values(lambda x, y: x+2*y) #2 4 4 6
        
        v = quantity.get_maximum_value()
        assert v == 6

        v = quantity.get_minimum_value()
        assert v == 2        

        i = quantity.get_maximum_index()
        assert i == 3

        i = quantity.get_minimum_index()
        assert i == 0        
        
        x,y = quantity.get_maximum_location()
        xref, yref = 2.0/3, 8.0/3
        assert x == xref
        assert y == yref

        v = quantity.get_values(interpolation_points = [[x,y]])
        assert num.allclose(v, 6)

        x,y = quantity.get_minimum_location()        
        v = quantity.get_values(interpolation_points = [[x,y]])
        assert num.allclose(v, 2)

        #Multiple locations for maximum -
        #Test that the algorithm picks the first occurrence        
        v = quantity.get_maximum_value(indices=[0,1,2])
        assert num.allclose(v, 4)

        i = quantity.get_maximum_index(indices=[0,1,2])
        assert i == 1
        
        x,y = quantity.get_maximum_location(indices=[0,1,2])
        xref, yref = 4.0/3, 4.0/3
        assert x == xref
        assert y == yref

        v = quantity.get_values(interpolation_points = [[x,y]])
        assert num.allclose(v, 4)        

        # More test of indices......
        v = quantity.get_maximum_value(indices=[2,3])
        assert num.allclose(v, 6)

        i = quantity.get_maximum_index(indices=[2,3])
        assert i == 3
        
        x,y = quantity.get_maximum_location(indices=[2,3])
        xref, yref = 2.0/3, 8.0/3
        assert x == xref
        assert y == yref

        v = quantity.get_values(interpolation_points = [[x,y]])
        assert num.allclose(v, 6)        

        

    def test_boundary_allocation(self):
        quantity = Quantity(self.mesh4,
                            [[1,2,3], [5,5,5], [0,0,9], [-6, 3, 3]])

        assert quantity.boundary_values.shape[0] == len(self.mesh4.boundary)


    def test_set_values(self):
        quantity = Quantity(self.mesh4)

        # get referece to data arrays
        centroid_values = quantity.centroid_values
        vertex_values = quantity.vertex_values

        quantity.set_values([[1,2,3], [5,5,5], [0,0,9], [-6, 3, 3]],
                            location = 'vertices')
        assert num.allclose(quantity.vertex_values,
                            [[1,2,3], [5,5,5], [0,0,9], [-6, 3, 3]])

        assert id(vertex_values) == id(quantity.vertex_values)
        assert num.allclose(quantity.centroid_values, [2., 5., 3., 0.]) #Centroid
        assert num.allclose(quantity.edge_values, [[2.5, 2.0, 1.5],
                                                   [5., 5., 5.],
                                                   [4.5, 4.5, 0.],
                                                   [3.0, -1.5, -1.5]])


        # Test default
        quantity.set_values([[1,2,3], [5,5,5], [0,0,9], [-6, 3, 3]])
        assert num.allclose(quantity.vertex_values,
                            [[1,2,3], [5,5,5], [0,0,9], [-6, 3, 3]])
        assert num.allclose(quantity.centroid_values, [2., 5., 3., 0.]) #Centroid
        assert num.allclose(quantity.edge_values, [[2.5, 2.0, 1.5],
                                                   [5., 5., 5.],
                                                   [4.5, 4.5, 0.],
                                                   [3.0, -1.5, -1.5]])

        # Test centroids
        quantity.set_values([1,2,3,4], location = 'centroids')
        assert num.allclose(quantity.centroid_values, [1., 2., 3., 4.]) #Centroid

        # Test exceptions
        try:
            quantity.set_values([[1,2,3], [5,5,5], [0,0,9], [-6, 3, 3]],
                                location = 'bas kamel tuba')
        except:
            pass



        try:
            quantity.set_values([[1,2,3], [0,0,9]])
        except ValueError:
            pass
        except:
            raise Exception('should have raised ValueeError')



    def test_set_values_const(self):
        quantity = Quantity(self.mesh4)

        quantity.set_values(1.0, location = 'vertices')
        assert num.allclose(quantity.vertex_values,
                            [[1,1,1], [1,1,1], [1,1,1], [1, 1, 1]])

        assert num.allclose(quantity.centroid_values, [1, 1, 1, 1]) #Centroid
        assert num.allclose(quantity.edge_values, [[1, 1, 1],
                                                   [1, 1, 1],
                                                   [1, 1, 1],
                                                   [1, 1, 1]])


        quantity.set_values(2.0, location = 'centroids')
        assert num.allclose(quantity.centroid_values, [2, 2, 2, 2])


    def test_set_values_func(self):
        quantity = Quantity(self.mesh4)

        def f(x, y):
            return x+y

        quantity.set_values(f, location = 'vertices')
        #print "quantity.vertex_values",quantity.vertex_values
        assert num.allclose(quantity.vertex_values,
                            [[2,0,2], [2,2,4], [4,2,4], [4,2,4]])
        assert num.allclose(quantity.centroid_values,
                            [4.0/3, 8.0/3, 10.0/3, 10.0/3])
        assert num.allclose(quantity.edge_values,
                            [[1,2,1], [3,3,2], [3,4,3], [3,4,3]])


        quantity.set_values(f, location = 'centroids')
        assert num.allclose(quantity.centroid_values,
                            [4.0/3, 8.0/3, 10.0/3, 10.0/3])


    def test_integral(self):
        quantity = Quantity(self.mesh4)

        # Try constants first
        const = 5
        quantity.set_values(const, location = 'vertices')
        #print 'Q', quantity.get_integral()

        assert num.allclose(quantity.get_integral(), self.mesh4.get_area() * const)

        # Try with a linear function
        def f(x, y):
            return x+y

        quantity.set_values(f, location = 'vertices')


        ref_integral = (4.0/3 + 8.0/3 + 10.0/3 + 10.0/3) * 2

        assert num.allclose (quantity.get_integral(), ref_integral)

    def test_integral_with_region(self):
        quantity = Quantity(self.mesh4)

        # Try constants first
        const = 5
        quantity.set_values(const, location = 'vertices')
        #print 'Q', quantity.get_integral()

        assert num.allclose(quantity.get_integral(), self.mesh4.get_area() * const)

        # Try with a linear function
        def f(x, y):
            return x+y

        quantity.set_values(f, location = 'vertices')


        from anuga import Region

        reg1 = Region(self.mesh4, indices=[2])
        ref_integral = (10.0/3) * 2

        assert num.allclose (quantity.get_integral(region=reg1), ref_integral)

        reg2 = Region(self.mesh4, indices=[2,3])
        ref_integral = (10.0/3 + 10.0/3) * 2

        assert num.allclose (quantity.get_integral(region=reg2), ref_integral)

        id=[2,3]
        ref_integral = (10.0/3 + 10.0/3) * 2

        assert num.allclose (quantity.get_integral(indices=id), ref_integral)

    def test_set_vertex_values(self):
        quantity = Quantity(self.mesh4)
        quantity.set_vertex_values([0,1,2,3,4,5])

        assert num.allclose(quantity.vertex_values,
                            [[1,0,2], [1,2,4], [4,2,5], [3,1,4]])
        assert num.allclose(quantity.centroid_values,
                            [1., 7./3, 11./3, 8./3]) #Centroid
        assert num.allclose(quantity.edge_values, [[1., 1.5, 0.5],
                                                   [3., 2.5, 1.5],
                                                   [3.5, 4.5, 3.],
                                                   [2.5, 3.5, 2]])


    def test_set_vertex_values_subset(self):
        quantity = Quantity(self.mesh4)
        quantity.set_vertex_values([0,1,2,3,4,5])
        quantity.set_vertex_values([0,20,30,50], indices = [0,2,3,5])

        assert num.allclose(quantity.vertex_values,
                            [[1,0,20], [1,20,4], [4,20,50], [30,1,4]])


    def test_set_vertex_values_using_general_interface(self):
        quantity = Quantity(self.mesh4)


        quantity.set_values([0,1,2,3,4,5])


        assert num.allclose(quantity.vertex_values,
                            [[1,0,2], [1,2,4], [4,2,5], [3,1,4]])

        #Centroid
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


        quantity.set_values([0,2,3,5], indices=[0,2,3,5])
        assert num.allclose(quantity.vertex_values,
                            [[0,0,2], [0,2,0], [0,2,5], [3,0,0]])


        # Constant
        quantity.set_values(0.0)
        quantity.set_values(3.14, indices=[0,2], location='vertices')

        # Indices refer to triangle numbers
        assert num.allclose(quantity.vertex_values,
                            [[3.14,3.14,3.14], [0,0,0],
                             [3.14,3.14,3.14], [0,0,0]])        
        


        # Now try with polygon (pick points where y>2)
        polygon = [[0,2.1], [4,2.1], [4,7], [0,7]]
        quantity.set_values(0.0)
        quantity.set_values(3.14, polygon=polygon)
        
        assert num.allclose(quantity.vertex_values,
                            [[0,0,0], [0,0,0], [0,0,0],
                             [3.14,3.14,3.14]])                


        # Another polygon (pick triangle 1 and 2 (rightmost triangles) 
        # using centroids
        polygon = [[2.1, 0.0], [3.5,0.1], [2,2.2], [0.2,2]]
        quantity.set_values(0.0)
        quantity.set_values(3.14, location='centroids', polygon=polygon)
        assert num.allclose(quantity.vertex_values,
                            [[0,0,0],
                             [3.14,3.14,3.14],
                             [3.14,3.14,3.14],                         
                             [0,0,0]])                


        # Same polygon now use vertices (default)
        polygon = [[2.1, 0.0], [3.5,0.1], [2,2.2], [0.2,2]]
        quantity.set_values(0.0)
        #print 'Here 2'
        quantity.set_values(3.14, polygon=polygon)
        assert num.allclose(quantity.vertex_values,
                            [[0,0,0],
                             [3.14,3.14,3.14],
                             [3.14,3.14,3.14],                         
                             [0,0,0]])                
        

        # Test input checking
        try:
            quantity.set_values(3.14, polygon=polygon, indices = [0,2])
        except:
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
        quantity.set_values(3.14, indices=[0,2], location='vertices')

        # Indices refer to triangle numbers here - not vertices (why?)
        assert num.allclose(quantity.vertex_values,
                            [[3.14,3.14,3.14], [0,0,0],
                             [3.14,3.14,3.14], [0,0,0]])        
        


        # Now try with polygon (pick points where y>2)
        polygon = num.array([[0,2.1], [4,2.1], [4,7], [0,7]])
        polygon += [G.xllcorner, G.yllcorner]
        
        quantity.set_values(0.0)
        quantity.set_values(3.14, polygon=polygon, location='centroids')
        assert num.allclose(quantity.vertex_values,
                            [[0,0,0], [0,0,0], [0,0,0],
                             [3.14,3.14,3.14]])                


        # Another polygon (pick triangle 1 and 2 (rightmost triangles)
        polygon = num.array([[2.1, 0.0], [3.5,0.1], [2,2.2], [0.2,2]])
        polygon += [G.xllcorner, G.yllcorner]
        
        quantity.set_values(0.0)
        quantity.set_values(3.14, polygon=polygon)
        msg = ('quantity.vertex_values=\n%s\nshould be close to\n'
               '[[0,0,0],\n'
               ' [3.14,3.14,3.14],\n'
               ' [3.14,3.14,3.14],\n'
               ' [0,0,0]]' % str(quantity.vertex_values))
        assert num.allclose(quantity.vertex_values,
                            [[0,0,0],
                             [3.14,3.14,3.14],
                             [3.14,3.14,3.14],                         
                             [0,0,0]]), msg



    def test_set_values_using_fit(self):


        quantity = Quantity(self.mesh4)

        #Get (enough) datapoints
        data_points = [[ 0.66666667, 0.66666667],
                       [ 1.33333333, 1.33333333],
                       [ 2.66666667, 0.66666667],
                       [ 0.66666667, 2.66666667],
                       [ 0.0, 1.0],
                       [ 0.0, 3.0],
                       [ 1.0, 0.0],
                       [ 1.0, 1.0],
                       [ 1.0, 2.0],
                       [ 1.0, 3.0],
                       [ 2.0, 1.0],
                       [ 3.0, 0.0],
                       [ 3.0, 1.0]]

        z = linear_function(data_points)

        #Use built-in fit_interpolate.fit
        quantity.set_values( Geospatial_data(data_points, z), alpha = 0 )
        #quantity.set_values(points = data_points, values = z, alpha = 0)


        answer = linear_function(quantity.domain.get_vertex_coordinates())
        #print quantity.vertex_values, answer
        assert num.allclose(quantity.vertex_values.flat, answer)


        #Now try by setting the same values directly
        vertex_attributes = fit_to_mesh(data_points,
                                        quantity.domain.get_nodes(),
                                        quantity.domain.get_triangles(),
                                        point_attributes=z,
                                        alpha = 0,
                                        verbose=False)

        #print vertex_attributes
        quantity.set_values(vertex_attributes)
        assert num.allclose(quantity.vertex_values.flat, answer)





    def test_test_set_values_using_fit_w_geo(self):


        #Mesh
        vertex_coordinates = [[0.76, 0.76],
                              [0.76, 5.76],
                              [5.76, 0.76]]
        triangles = [[0,2,1]]

        mesh_georef = Geo_reference(56,-0.76,-0.76)
        mesh1 = Generic_Domain(vertex_coordinates, triangles,
                       geo_reference = mesh_georef)
        mesh1.check_integrity()

        #Quantity
        quantity = Quantity(mesh1)

        #Data
        data_points = [[ 201.0, 401.0],
                       [ 201.0, 403.0],
                       [ 203.0, 401.0]]

        z = [2, 4, 4]

        data_georef = Geo_reference(56,-200,-400)


        #Reference
        ref = fit_to_mesh(data_points, vertex_coordinates, triangles,
                          point_attributes=z,
                          data_origin = data_georef.get_origin(),
                          mesh_origin = mesh_georef.get_origin(),
                          alpha = 0)

        assert num.allclose( ref, [0,5,5] )


        #Test set_values

        quantity.set_values( Geospatial_data(data_points, z, data_georef), alpha = 0 )

        #quantity.set_values(points = data_points,
        #                    values = z,
        #                    data_georef = data_georef,
        #                    alpha = 0)


        #quantity.set_values(points = data_points,
        #                    values = z,
        #                    data_georef = data_georef,
        #                    alpha = 0)
        assert num.allclose(quantity.vertex_values.flat, ref)



        #Test set_values using geospatial data object
        quantity.vertex_values[:] = 0.0

        geo = Geospatial_data(data_points, z, data_georef)


        quantity.set_values(geospatial_data = geo, alpha = 0)
        assert num.allclose(quantity.vertex_values.flat, ref)



    def test_set_values_from_file1(self):
        quantity = Quantity(self.mesh4)

        #Get (enough) datapoints
        data_points = [[ 0.66666667, 0.66666667],
                       [ 1.33333333, 1.33333333],
                       [ 2.66666667, 0.66666667],
                       [ 0.66666667, 2.66666667],
                       [ 0.0, 1.0],
                       [ 0.0, 3.0],
                       [ 1.0, 0.0],
                       [ 1.0, 1.0],
                       [ 1.0, 2.0],
                       [ 1.0, 3.0],
                       [ 2.0, 1.0],
                       [ 3.0, 0.0],
                       [ 3.0, 1.0]]

        data_geo_spatial = Geospatial_data(data_points,
                         geo_reference = Geo_reference(56, 0, 0))
        data_points_absolute = data_geo_spatial.get_data_points(absolute=True)
        attributes = linear_function(data_points_absolute)
        att = 'spam_and_eggs'
        
        #Create .txt file
        ptsfile = tempfile.mktemp(".txt")
        file = open(ptsfile,"w")
        file.write(" x,y," + att + " \n")
        for data_point, attribute in map(None, data_points_absolute
                                         ,attributes):
            row = str(data_point[0]) + ',' + str(data_point[1]) \
                  + ',' + str(attribute)
            file.write(row + "\n")
        file.close()


        #Check that values can be set from file
        quantity.set_values(filename = ptsfile,
                            attribute_name = att, alpha = 0)
        answer = linear_function(quantity.domain.get_vertex_coordinates())

        #print quantity.vertex_values.flat
        #print answer


        assert num.allclose(quantity.vertex_values.flat, answer)


        #Check that values can be set from file using default attribute
        quantity.set_values(filename = ptsfile, alpha = 0)
        assert num.allclose(quantity.vertex_values.flat, answer)

        #Cleanup
        import os
        os.remove(ptsfile)



    def Xtest_set_values_from_file_using_polygon(self):
        """test_set_values_from_file_using_polygon(self):
        
        Test that polygon restriction works for general points data
        """
        
        quantity = Quantity(self.mesh4)

        #Get (enough) datapoints
        data_points = [[ 0.66666667, 0.66666667],
                       [ 1.33333333, 1.33333333],
                       [ 2.66666667, 0.66666667],
                       [ 0.66666667, 2.66666667],
                       [ 0.0, 1.0],
                       [ 0.0, 3.0],
                       [ 1.0, 0.0],
                       [ 1.0, 1.0],
                       [ 1.0, 2.0],
                       [ 1.0, 3.0],
                       [ 2.0, 1.0],
                       [ 3.0, 0.0],
                       [ 3.0, 1.0]]

        data_geo_spatial = Geospatial_data(data_points,
                         geo_reference = Geo_reference(56, 0, 0))
        data_points_absolute = data_geo_spatial.get_data_points(absolute=True)
        attributes = linear_function(data_points_absolute)
        att = 'spam_and_eggs'
        
        #Create .txt file
        ptsfile = tempfile.mktemp(".txt")
        file = open(ptsfile,"w")
        file.write(" x,y," + att + " \n")
        for data_point, attribute in map(None, data_points_absolute
                                         ,attributes):
            row = str(data_point[0]) + ',' + str(data_point[1]) \
                  + ',' + str(attribute)
            file.write(row + "\n")
        file.close()

        # Create restricting polygon (containing node #4 (2,2) and 
        # centroid of triangle #1 (bce)
        polygon = [[1.0, 1.0], [4.0, 1.0],
                   [4.0, 4.0], [1.0, 4.0]]

        #print self.mesh4.nodes
        #print inside_polygon(self.mesh4.nodes, polygon)
        assert num.allclose(inside_polygon(self.mesh4.nodes, polygon), 4)

        #print quantity.domain.get_vertex_coordinates()
        #print quantity.domain.get_nodes()        
        
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

        #print quantity.vertex_values.flat
        #print answer

        # Check vertices in polygon have been set
        assert num.allclose(num.take(quantity.vertex_values.flat, indices),
                            answer)

        # Check vertices outside polygon are zero
        indices = outside_polygon(quantity.domain.get_vertex_coordinates(), 
                                  polygon)        
        assert num.allclose(num.take(quantity.vertex_values.flat, indices),
                            0.0)        

        #Cleanup
        import os
        os.remove(ptsfile)


        

    def test_cache_test_set_values_from_file(self):
        # FIXME (Ole): What is this about?
        # I don't think it checks anything new
        quantity = Quantity(self.mesh4)

        #Get (enough) datapoints
        data_points = [[ 0.66666667, 0.66666667],
                       [ 1.33333333, 1.33333333],
                       [ 2.66666667, 0.66666667],
                       [ 0.66666667, 2.66666667],
                       [ 0.0, 1.0],
                       [ 0.0, 3.0],
                       [ 1.0, 0.0],
                       [ 1.0, 1.0],
                       [ 1.0, 2.0],
                       [ 1.0, 3.0],
                       [ 2.0, 1.0],
                       [ 3.0, 0.0],
                       [ 3.0, 1.0]]

        georef = Geo_reference(56, 0, 0)
        data_geo_spatial = Geospatial_data(data_points,
                                           geo_reference=georef)
                                           
        data_points_absolute = data_geo_spatial.get_data_points(absolute=True)
        attributes = linear_function(data_points_absolute)
        att = 'spam_and_eggs'
        
        # Create .txt file
        ptsfile = tempfile.mktemp(".txt")
        file = open(ptsfile,"w")
        file.write(" x,y," + att + " \n")
        for data_point, attribute in map(None, data_points_absolute
                                         ,attributes):
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
        
        
        #Cleanup
        import os
        os.remove(ptsfile)

    def test_set_values_from_lat_long(self):
        quantity = Quantity(self.mesh_onslow)

        #Get (enough) datapoints
        data_points = [[-21.5, 114.5],[-21.4, 114.6],[-21.45,114.65],
                       [-21.35, 114.65],[-21.45, 114.55],[-21.45,114.6]]

        data_geo_spatial = Geospatial_data(data_points,
                                           points_are_lats_longs=True)
        points_UTM = data_geo_spatial.get_data_points(absolute=True)
        attributes = linear_function(points_UTM)
        att = 'elevation'
        
        #Create .txt file
        txt_file = tempfile.mktemp(".txt")
        file = open(txt_file,"w")
        file.write(" lat,long," + att + " \n")
        for data_point, attribute in map(None, data_points, attributes):
            row = str(data_point[0]) + ',' + str(data_point[1]) \
                  + ',' + str(attribute)
            #print "row", row 
            file.write(row + "\n")
        file.close()


        #Check that values can be set from file
        quantity.set_values(filename=txt_file,
                            attribute_name=att, 
                            alpha=0)
        answer = linear_function(quantity.domain.get_vertex_coordinates())

        #print "quantity.vertex_values.flat", quantity.vertex_values.flat
        #print "answer",answer

        assert num.allclose(quantity.vertex_values.flat, answer)


        #Check that values can be set from file using default attribute
        quantity.set_values(filename=txt_file, alpha=0)
        assert num.allclose(quantity.vertex_values.flat, answer)

        #Cleanup
        import os
        os.remove(txt_file)
         
    def test_set_values_from_lat_long_2(self):
        quantity = Quantity(self.mesh_onslow)

        #Get (enough) datapoints
        data_points = [[-21.5, 114.5],[-21.4, 114.6],[-21.45,114.65],
                       [-21.35, 114.65],[-21.45, 114.55],[-21.45,114.6]]

        data_geo_spatial = Geospatial_data(data_points,
                                           points_are_lats_longs=True)
        points_UTM = data_geo_spatial.get_data_points(absolute=True)
        attributes = linear_function(points_UTM)
        att = 'elevation'
        
        #Create .txt file
        txt_file = tempfile.mktemp(".txt")
        file = open(txt_file,"w")
        file.write(" lat,long," + att + " \n")
        for data_point, attribute in map(None, data_points, attributes):
            row = str(data_point[0]) + ',' + str(data_point[1]) \
                  + ',' + str(attribute)
            #print "row", row 
            file.write(row + "\n")
        file.close()


        #Check that values can be set from file
        quantity.set_values(filename=txt_file,
                            attribute_name=att, alpha=0)
        answer = linear_function(quantity.domain.get_vertex_coordinates())

        #print "quantity.vertex_values.flat", quantity.vertex_values.flat
        #print "answer",answer

        assert num.allclose(quantity.vertex_values.flat, answer)


        #Check that values can be set from file using default attribute
        quantity.set_values(filename=txt_file, alpha=0)
        assert num.allclose(quantity.vertex_values.flat, answer)

        #Cleanup
        import os
        os.remove(txt_file)
        
    def test_set_values_from_UTM_pts(self):
        quantity = Quantity(self.mesh_onslow)

        #Get (enough) datapoints
        data_points = [[-21.5, 114.5],[-21.4, 114.6],[-21.45,114.65],
                       [-21.35, 114.65],[-21.45, 114.55],[-21.45,114.6]]

        data_geo_spatial = Geospatial_data(data_points,
                                           points_are_lats_longs=True)
        points_UTM = data_geo_spatial.get_data_points(absolute=True)
        attributes = linear_function(points_UTM)
        att = 'elevation'
        
        #Create .txt file
        txt_file = tempfile.mktemp(".txt")
        file = open(txt_file,"w")
        file.write(" x,y," + att + " \n")
        for data_point, attribute in map(None, points_UTM, attributes):
            row = str(data_point[0]) + ',' + str(data_point[1]) \
                  + ',' + str(attribute)
            #print "row", row 
            file.write(row + "\n")
        file.close()


        pts_file = tempfile.mktemp(".pts")        
        convert = Geospatial_data(txt_file)
        convert.export_points_file(pts_file)

        #Check that values can be set from file
        quantity.set_values_from_file(pts_file, att, 0,
                                      'vertices', None)
        answer = linear_function(quantity.domain.get_vertex_coordinates())
        #print "quantity.vertex_values.flat", quantity.vertex_values.flat
        #print "answer",answer
        assert num.allclose(quantity.vertex_values.flat, answer)

        #Check that values can be set from file
        quantity.set_values(filename=pts_file,
                            attribute_name=att, alpha=0)
        answer = linear_function(quantity.domain.get_vertex_coordinates())
        #print "quantity.vertex_values.flat", quantity.vertex_values.flat
        #print "answer",answer
        assert num.allclose(quantity.vertex_values.flat, answer)


        #Check that values can be set from file using default attribute
        quantity.set_values(filename=txt_file, alpha=0)
        assert num.allclose(quantity.vertex_values.flat, answer)

        #Cleanup
        import os
        os.remove(txt_file)
        os.remove(pts_file)
        
    def test_set_values_from_UTM_pts_verbose(self):
        quantity = Quantity(self.mesh_onslow)

        #Get (enough) datapoints
        data_points = [[-21.5, 114.5],[-21.4, 114.6],[-21.45,114.65],
                       [-21.35, 114.65],[-21.45, 114.55],[-21.45,114.6],
                       [-21.5, 114.5],[-21.4, 114.6],[-21.45,114.65],
                       [-21.35, 114.65],[-21.45, 114.55],[-21.45,114.6],
                       [-21.5, 114.5],[-21.4, 114.6],[-21.45,114.65],
                       [-21.35, 114.65],[-21.45, 114.55],[-21.45,114.6],
                       [-21.5, 114.5],[-21.4, 114.6],[-21.45,114.65],
                       [-21.35, 114.65],[-21.45, 114.55],[-21.45,114.6],
                       [-21.5, 114.5],[-21.4, 114.6],[-21.45,114.65],
                       [-21.35, 114.65],[-21.45, 114.55],[-21.45,114.6],
                       [-21.5, 114.5],[-21.4, 114.6],[-21.45,114.65],
                       [-21.35, 114.65],[-21.45, 114.55],[-21.45,114.6],
                       [-21.5, 114.5],[-21.4, 114.6],[-21.45,114.65],
                       [-21.35, 114.65],[-21.45, 114.55],[-21.45,114.6],
                       [-21.35, 114.65],[-21.45, 114.55],[-21.45,114.6],
                       [-21.5, 114.5],[-21.4, 114.6],[-21.45,114.65],
                       [-21.35, 114.65],[-21.45, 114.55],[-21.45,114.6],
                       [-21.5, 114.5],[-21.4, 114.6],[-21.45,114.65],
                       [-21.35, 114.65],[-21.45, 114.55],[-21.45,114.6],
                       [-21.5, 114.5],[-21.4, 114.6],[-21.45,114.65],
                       [-21.35, 114.65],[-21.45, 114.55],[-21.45,114.6],
                       [-21.5, 114.5],[-21.4, 114.6],[-21.45,114.65],
                       [-21.35, 114.65],[-21.45, 114.55],[-21.45,114.6],
                       [-21.5, 114.5],[-21.4, 114.6],[-21.45,114.65],
                       [-21.35, 114.65],[-21.45, 114.55],[-21.45,114.6],
                       [-21.5, 114.5],[-21.4, 114.6],[-21.45,114.65],
                       [-21.35, 114.65],[-21.45, 114.55],[-21.45,114.6],
                       ]

        data_geo_spatial = Geospatial_data(data_points,
                                           points_are_lats_longs=True)
        points_UTM = data_geo_spatial.get_data_points(absolute=True)
        attributes = linear_function(points_UTM)
        att = 'elevation'
        
        #Create .txt file
        txt_file = tempfile.mktemp(".txt")
        file = open(txt_file,"w")
        file.write(" x,y," + att + " \n")
        for data_point, attribute in map(None, points_UTM, attributes):
            row = str(data_point[0]) + ',' + str(data_point[1]) \
                  + ',' + str(attribute)
            #print "row", row 
            file.write(row + "\n")
        file.close()


        pts_file = tempfile.mktemp(".pts")        
        convert = Geospatial_data(txt_file)
        convert.export_points_file(pts_file)

        #Check that values can be set from file
        quantity.set_values_from_file(pts_file, att, 0,
                                      'vertices', None, verbose = False,
                                      max_read_lines=2)
        answer = linear_function(quantity.domain.get_vertex_coordinates())
        #print "quantity.vertex_values.flat", quantity.vertex_values.flat
        #print "answer",answer
        assert num.allclose(quantity.vertex_values.flat, answer)

        #Check that values can be set from file
        quantity.set_values(filename=pts_file,
                            attribute_name=att, alpha=0)
        answer = linear_function(quantity.domain.get_vertex_coordinates())
        #print "quantity.vertex_values.flat", quantity.vertex_values.flat
        #print "answer",answer
        assert num.allclose(quantity.vertex_values.flat, answer)


        #Check that values can be set from file using default attribute
        quantity.set_values(filename=txt_file, alpha=0)
        assert num.allclose(quantity.vertex_values.flat, answer)

        #Cleanup
        import os
        os.remove(txt_file)
        os.remove(pts_file)
        
    def test_set_values_from_file_with_georef1(self):

        #Mesh in zone 56 (absolute coords)

        x0 = 314036.58727982
        y0 = 6224951.2960092

        a = [x0+0.0, y0+0.0]
        b = [x0+0.0, y0+2.0]
        c = [x0+2.0, y0+0.0]
        d = [x0+0.0, y0+4.0]
        e = [x0+2.0, y0+2.0]
        f = [x0+4.0, y0+0.0]

        points = [a, b, c, d, e, f]

        #bac, bce, ecf, dbe
        elements = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4] ]

        #absolute going in ..
        mesh4 = Generic_Domain(points, elements,
                       geo_reference = Geo_reference(56, 0, 0))
        mesh4.check_integrity()
        quantity = Quantity(mesh4)

        #Get (enough) datapoints (relative to georef)
        data_points_rel = [[ 0.66666667, 0.66666667],
                       [ 1.33333333, 1.33333333],
                       [ 2.66666667, 0.66666667],
                       [ 0.66666667, 2.66666667],
                       [ 0.0, 1.0],
                       [ 0.0, 3.0],
                       [ 1.0, 0.0],
                       [ 1.0, 1.0],
                       [ 1.0, 2.0],
                       [ 1.0, 3.0],
                       [ 2.0, 1.0],
                       [ 3.0, 0.0],
                       [ 3.0, 1.0]]

        data_geo_spatial = Geospatial_data(data_points_rel,
                         geo_reference = Geo_reference(56, x0, y0))
        data_points_absolute = data_geo_spatial.get_data_points(absolute=True)
        attributes = linear_function(data_points_absolute)
        att = 'spam_and_eggs'
        
        #Create .txt file
        ptsfile = tempfile.mktemp(".txt")
        file = open(ptsfile,"w")
        file.write(" x,y," + att + " \n")
        for data_point, attribute in map(None, data_points_absolute
                                         ,attributes):
            row = str(data_point[0]) + ',' + str(data_point[1]) \
                  + ',' + str(attribute)
            file.write(row + "\n")
        file.close()

        #file = open(ptsfile, 'r')
        #lines = file.readlines()
        #file.close()
     

        #Check that values can be set from file
        quantity.set_values(filename=ptsfile,
                            attribute_name=att, alpha=0)
        answer = linear_function(quantity.domain.get_vertex_coordinates())

        assert num.allclose(quantity.vertex_values.flat, answer)


        #Check that values can be set from file using default attribute
        quantity.set_values(filename=ptsfile, alpha=0)
        assert num.allclose(quantity.vertex_values.flat, answer)

        #Cleanup
        import os
        os.remove(ptsfile)


    def test_set_values_from_file_with_georef2(self):

        #Mesh in zone 56 (relative coords)

        x0 = 314036.58727982
        y0 = 6224951.2960092
        #x0 = 0.0
        #y0 = 0.0

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]

        #bac, bce, ecf, dbe
        elements = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4] ]

        mesh4 = Generic_Domain(points, elements,
                       geo_reference = Geo_reference(56, x0, y0))
        mesh4.check_integrity()
        quantity = Quantity(mesh4)

        #Get (enough) datapoints
        data_points = [[ x0+0.66666667, y0+0.66666667],
                       [ x0+1.33333333, y0+1.33333333],
                       [ x0+2.66666667, y0+0.66666667],
                       [ x0+0.66666667, y0+2.66666667],
                       [ x0+0.0, y0+1.0],
                       [ x0+0.0, y0+3.0],
                       [ x0+1.0, y0+0.0],
                       [ x0+1.0, y0+1.0],
                       [ x0+1.0, y0+2.0],
                       [ x0+1.0, y0+3.0],
                       [ x0+2.0, y0+1.0],
                       [ x0+3.0, y0+0.0],
                       [ x0+3.0, y0+1.0]]


        data_geo_spatial = Geospatial_data(data_points,
                         geo_reference = Geo_reference(56, 0, 0))
        data_points_absolute = data_geo_spatial.get_data_points(absolute=True)
        attributes = linear_function(data_points_absolute)
        att = 'spam_and_eggs'
        
        #Create .txt file
        ptsfile = tempfile.mktemp(".txt")
        file = open(ptsfile,"w")
        file.write(" x,y," + att + " \n")
        for data_point, attribute in map(None, data_points_absolute
                                         ,attributes):
            row = str(data_point[0]) + ',' + str(data_point[1]) \
                  + ',' + str(attribute)
            file.write(row + "\n")
        file.close()


        #Check that values can be set from file
        quantity.set_values(filename=ptsfile,
                            attribute_name=att, alpha=0)
        answer = linear_function(quantity.domain. \
                                 get_vertex_coordinates(absolute=True))


        assert num.allclose(quantity.vertex_values.flat, answer)


        #Check that values can be set from file using default attribute
        quantity.set_values(filename=ptsfile, alpha=0)
        assert num.allclose(quantity.vertex_values.flat, answer)

        #Cleanup
        import os
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

        #bac, bce, ecf, dbe
        elements = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4] ]

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
        cellsize  = 1.0
        NODATA_value =  -9999
        
        #xllcorner = 0
        #yllcorner = 100
        #cellsize  = 10
        #NODATA_value =  -9999
        
        #Create .asc file
        #txt_file = tempfile.mktemp(".asc")
        txt_file = 'test_asc.asc'
        datafile = open(txt_file,"w")
        datafile.write('ncols '+str(ncols)+"\n")
        datafile.write('nrows '+str(nrows)+"\n")
        datafile.write('xllcorner '+str(xllcorner)+"\n")
        datafile.write('yllcorner '+str(yllcorner)+"\n")
        datafile.write('cellsize '+str(cellsize)+"\n")
        datafile.write('NODATA_value '+str(NODATA_value)+"\n")
        
        x = num.linspace(xllcorner, xllcorner+(ncols-1)*cellsize, ncols)
        y = num.linspace(yllcorner, yllcorner+(nrows-1)*cellsize, nrows)
        points = axes2points(x, y)
        
        #print points
        #print x.shape, x
        #print y.shape, y
        
        datavalues = linear_function(points)
        #print datavalues 
        
        datavalues = datavalues.reshape(nrows,ncols)

        #print datavalues
        #print datavalues.shape
        for row in datavalues:
            #print row
            datafile.write(" ".join(str(elem) for elem in row) + "\n")         
        datafile.close()

        #print quantity.vertex_values
        #print quantity.centroid_values 

        quantity.set_values_from_utm_grid_file(txt_file,
                             location='vertices',
                             indices=None,
                             verbose=False)

        # check order of vertices

        answer = [[  6.,   0.,   2.],
                  [  6.,   2.,   8.],
                  [  8.,  2.,   4.],
                  [ 12.,   6.,   8.]]
        #print quantity.vertex_values
        assert num.allclose(quantity.vertex_values, answer)
        
        
        #print quantity.vertex_values
        #print quantity.centroid_values 
        quantity.set_values(0.0)
        
        
        #print quantity.vertex_values
        #print quantity.centroid_values 
        quantity.set_values_from_utm_grid_file(txt_file,
                             location='centroids',
                             indices=None,
                             verbose=False)
        
        #print quantity.vertex_values
        #print quantity.centroid_values      
        

        answer = [ 2.66666667,  5.33333333,  4.66666667,  8.66666667]
        
        assert num.allclose(quantity.centroid_values, answer)
        #Cleanup
        #import os
        os.remove(txt_file)        
        


    def test_set_values_from_quantity(self):

        quantity1 = Quantity(self.mesh4)
        quantity1.set_vertex_values([0,1,2,3,4,5])

        assert num.allclose(quantity1.vertex_values,
                            [[1,0,2], [1,2,4], [4,2,5], [3,1,4]])


        quantity2 = Quantity(self.mesh4)
        quantity2.set_values(quantity=quantity1)
        assert num.allclose(quantity2.vertex_values,
                            [[1,0,2], [1,2,4], [4,2,5], [3,1,4]])

        quantity2.set_values(quantity = 2*quantity1)
        assert num.allclose(quantity2.vertex_values,
                            [[2,0,4], [2,4,8], [8,4,10], [6,2,8]])

        quantity2.set_values(quantity = 2*quantity1 + 3)
        assert num.allclose(quantity2.vertex_values,
                            [[5,3,7], [5,7,11], [11,7,13], [9,5,11]])


        #Check detection of quantity as first orgument
        quantity2.set_values(2*quantity1 + 3)
        assert num.allclose(quantity2.vertex_values,
                            [[5,3,7], [5,7,11], [11,7,13], [9,5,11]])



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
        quantity1.set_vertex_values([0,1,2,3,4,5])

        assert num.allclose(quantity1.vertex_values,
                            [[1,0,2], [1,2,4], [4,2,5], [3,1,4]])


        quantity2 = Quantity(self.mesh4)
        quantity2.set_values(quantity=quantity1,
                             polygon=polygon)
                             
        msg = 'Only node #4(e) at (2,2) should have values applied '
        assert num.allclose(quantity2.vertex_values,
                            [[0,0,0], [0,0,4], [4,0,0], [0,0,4]]), msg        
                            #bac,     bce,     ecf,     dbe
                        


    def test_overloading(self):

        quantity1 = Quantity(self.mesh4)
        quantity1.set_vertex_values([0,1,2,3,4,5])

        assert num.allclose(quantity1.vertex_values,
                            [[1,0,2], [1,2,4], [4,2,5], [3,1,4]])


        quantity2 = Quantity(self.mesh4)
        quantity2.set_values([[1,2,3], [5,5,5], [0,0,9], [-6, 3, 3]],
                             location = 'vertices')



        quantity3 = Quantity(self.mesh4)
        quantity3.set_values([[2,2,2], [7,8,9], [7,6,3], [3, 8, -8]],
                             location = 'vertices')


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

        #Scaling
        Q = quantity1*3
        assert num.allclose(Q.vertex_values, quantity1.vertex_values*3)
        assert num.allclose(Q.centroid_values, quantity1.centroid_values*3)
        assert num.allclose(Q.edge_values, quantity1.edge_values*3)
        Q = 3*quantity1
        assert num.allclose(Q.vertex_values, quantity1.vertex_values*3)

        #Multiplication
        Q = quantity1 * quantity2
        #print Q.vertex_values
        #print Q.centroid_values
        #print quantity1.centroid_values
        #print quantity2.centroid_values

        assert num.allclose(Q.vertex_values,
                            quantity1.vertex_values * quantity2.vertex_values)

        #Linear combinations
        Q = 4*quantity1 + 2
        assert num.allclose(Q.vertex_values,
                            4*quantity1.vertex_values + 2)

        Q = quantity1*quantity2 + 2
        assert num.allclose(Q.vertex_values,
                            quantity1.vertex_values * quantity2.vertex_values + 2)

        Q = quantity1*quantity2 + quantity3
        assert num.allclose(Q.vertex_values,
                            quantity1.vertex_values * quantity2.vertex_values +
                        quantity3.vertex_values)
        Q = quantity1*quantity2 + 3*quantity3
        assert num.allclose(Q.vertex_values,
                            quantity1.vertex_values * quantity2.vertex_values +
                            3*quantity3.vertex_values)
        Q = quantity1*quantity2 + 3*quantity3 + 5.0
        assert num.allclose(Q.vertex_values,
                            quantity1.vertex_values * quantity2.vertex_values +
                            3*quantity3.vertex_values + 5)

        Q = quantity1*quantity2 - quantity3
        assert num.allclose(Q.vertex_values,
                            quantity1.vertex_values * quantity2.vertex_values -
                            quantity3.vertex_values)
        Q = 1.5*quantity1*quantity2 - 3*quantity3 + 5.0
        assert num.allclose(Q.vertex_values,
                            1.5*quantity1.vertex_values * quantity2.vertex_values -
                            3*quantity3.vertex_values + 5)

        #Try combining quantities and arrays and scalars
        Q = 1.5*quantity1*quantity2.vertex_values -\
            3*quantity3.vertex_values + 5.0
        assert num.allclose(Q.vertex_values,
                            1.5*quantity1.vertex_values * quantity2.vertex_values -
                            3*quantity3.vertex_values + 5)


        #Powers
        Q = quantity1**2
        assert num.allclose(Q.vertex_values, quantity1.vertex_values**2)

        Q = quantity1**2 +quantity2**2
        assert num.allclose(Q.vertex_values,
                            quantity1.vertex_values**2 + \
                            quantity2.vertex_values**2)

        Q = (quantity1**2 +quantity2**2)**0.5
        assert num.allclose(Q.vertex_values,
                            (quantity1.vertex_values**2 + \
                            quantity2.vertex_values**2)**0.5)







    def test_compute_gradient(self):
        quantity = Quantity(self.mesh4)

        #Set up for a gradient of (2,0) at mid triangle
        quantity.set_values([2.0, 4.0, 6.0, 2.0],
                            location = 'centroids')


        #Gradients
        quantity.compute_gradients()

        a = quantity.x_gradient
        b = quantity.y_gradient
        #print self.mesh4.centroid_coordinates
        #print a, b

        #The central triangle (1)
        #(using standard gradient based on neigbours controid values)
        assert num.allclose(a[1], 2.0)
        assert num.allclose(b[1], 0.0)


        #Left triangle (0) using two point gradient
        #q0 = q1 + a*(x0-x1) + b*(y0-y1)  <=>
        #2  = 4  + a*(-2/3)  + b*(-2/3)
        assert num.allclose(a[0] + b[0], 3)
        #From orthogonality (a*(y0-y1) + b*(x0-x1) == 0)
        assert num.allclose(a[0] - b[0], 0)


        #Right triangle (2) using two point gradient
        #q2 = q1 + a*(x2-x1) + b*(y2-y1)  <=>
        #6  = 4  + a*(4/3)  + b*(-2/3)
        assert num.allclose(2*a[2] - b[2], 3)
        #From orthogonality (a*(y1-y2) + b*(x2-x1) == 0)
        assert num.allclose(a[2] + 2*b[2], 0)


        #Top triangle (3) using two point gradient
        #q3 = q1 + a*(x3-x1) + b*(y3-y1)  <=>
        #2  = 4  + a*(-2/3)  + b*(4/3)
        assert num.allclose(a[3] - 2*b[3], 3)
        #From orthogonality (a*(y1-y3) + b*(x3-x1) == 0)
        assert num.allclose(2*a[3] + b[3], 0)



        #print a, b
        quantity.extrapolate_second_order()

        #Apply q(x,y) = qc + a*(x-xc) + b*(y-yc)
        assert num.allclose(quantity.vertex_values[0,:], [3., 0.,  3.])
        assert num.allclose(quantity.vertex_values[1,:], [4./3, 16./3,  16./3])


        #a = 1.2, b=-0.6
        #q(4,0) = 6 + a*(4 - 8/3) + b*(-2/3)
        assert num.allclose(quantity.vertex_values[2,2], 8)

    def test_get_gradients(self):
        quantity = Quantity(self.mesh4)

        #Set up for a gradient of (2,0) at mid triangle
        quantity.set_values([2.0, 4.0, 6.0, 2.0],
                            location = 'centroids')


        #Gradients
        quantity.compute_gradients()

        a, b = quantity.get_gradients()
        #print self.mesh4.centroid_coordinates
        #print a, b

        #The central triangle (1)
        #(using standard gradient based on neigbours controid values)
        assert num.allclose(a[1], 2.0)
        assert num.allclose(b[1], 0.0)


        #Left triangle (0) using two point gradient
        #q0 = q1 + a*(x0-x1) + b*(y0-y1)  <=>
        #2  = 4  + a*(-2/3)  + b*(-2/3)
        assert num.allclose(a[0] + b[0], 3)
        #From orthogonality (a*(y0-y1) + b*(x0-x1) == 0)
        assert num.allclose(a[0] - b[0], 0)


        #Right triangle (2) using two point gradient
        #q2 = q1 + a*(x2-x1) + b*(y2-y1)  <=>
        #6  = 4  + a*(4/3)  + b*(-2/3)
        assert num.allclose(2*a[2] - b[2], 3)
        #From orthogonality (a*(y1-y2) + b*(x2-x1) == 0)
        assert num.allclose(a[2] + 2*b[2], 0)


        #Top triangle (3) using two point gradient
        #q3 = q1 + a*(x3-x1) + b*(y3-y1)  <=>
        #2  = 4  + a*(-2/3)  + b*(4/3)
        assert num.allclose(a[3] - 2*b[3], 3)
        #From orthogonality (a*(y1-y3) + b*(x3-x1) == 0)
        assert num.allclose(2*a[3] + b[3], 0)


    def test_second_order_extrapolation2(self):
        quantity = Quantity(self.mesh4)

        #Set up for a gradient of (3,1), f(x) = 3x+y
        quantity.set_values([2.0+2.0/3, 4.0+4.0/3, 8.0+2.0/3, 2.0+8.0/3],
                            location = 'centroids')

        #Gradients
        quantity.compute_gradients()

        a = quantity.x_gradient
        b = quantity.y_gradient
        
        #print a, b

        assert num.allclose(a[1], 3.0)
        assert num.allclose(b[1], 1.0)

        #Work out the others

        quantity.extrapolate_second_order()

        #print quantity.vertex_values
        assert num.allclose(quantity.vertex_values[1,0], 2.0)
        assert num.allclose(quantity.vertex_values[1,1], 6.0)
        assert num.allclose(quantity.vertex_values[1,2], 8.0)



    def test_backup_saxpy_centroid_values(self):
        quantity = Quantity(self.mesh4)

        #Set up for a gradient of (3,1), f(x) = 3x+y
        c_values = num.array([2.0+2.0/3, 4.0+4.0/3, 8.0+2.0/3, 2.0+8.0/3])
        d_values = num.array([1.0, 2.0, 3.0, 4.0])
        quantity.set_values(c_values, location = 'centroids')

        #Backup
        quantity.backup_centroid_values()

        #print quantity.vertex_values
        assert num.allclose(quantity.centroid_values, quantity.centroid_backup_values)


        quantity.set_values(d_values, location = 'centroids')

        quantity.saxpy_centroid_values(2.0, 3.0)

        assert num.allclose(quantity.centroid_values, 2.0*d_values + 3.0*c_values)



    def test_first_order_extrapolator(self):
        quantity = Quantity(self.mesh4)

        #Test centroids
        quantity.set_values([1.,2.,3.,4.], location = 'centroids')
        assert num.allclose(quantity.centroid_values, [1, 2, 3, 4]) #Centroid

        #Extrapolate
        quantity.extrapolate_first_order()

        #Check that gradient is zero
        a,b = quantity.get_gradients()
        assert num.allclose(a, [0,0,0,0])
        assert num.allclose(b, [0,0,0,0])

        #Check vertices but not edge values
        assert num.allclose(quantity.vertex_values,
                            [[1,1,1], [2,2,2], [3,3,3], [4, 4, 4]])


    def test_second_order_extrapolator(self):
        quantity = Quantity(self.mesh4)

        #Set up for a gradient of (3,0) at mid triangle
        quantity.set_values([2.0, 4.0, 8.0, 2.0],
                            location = 'centroids')



        quantity.extrapolate_second_order()
        quantity.limit()


        #Assert that central triangle is limited by neighbours
        assert quantity.vertex_values[1,0] >= quantity.vertex_values[0,0]
        assert quantity.vertex_values[1,0] >= quantity.vertex_values[3,1]

        assert quantity.vertex_values[1,1] <= quantity.vertex_values[2,1]
        assert quantity.vertex_values[1,1] >= quantity.vertex_values[0,2]

        assert quantity.vertex_values[1,2] <= quantity.vertex_values[2,0]
        assert quantity.vertex_values[1,2] >= quantity.vertex_values[3,1]


        #Assert that quantities are conserved
        for k in range(quantity.centroid_values.shape[0]):
            assert num.allclose (quantity.centroid_values[k],
                             num.sum(quantity.vertex_values[k,:])/3)





    def test_limit_vertices_by_all_neighbours(self):
        quantity = Quantity(self.mesh4)

        #Create a deliberate overshoot (e.g. from gradient computation)
        quantity.set_values([[3,0,3], [2,2,6], [5,3,8], [8,3,5]])


        #Limit
        quantity.limit_vertices_by_all_neighbours()

        #Assert that central triangle is limited by neighbours
        assert quantity.vertex_values[1,0] >= quantity.vertex_values[0,0]
        assert quantity.vertex_values[1,0] <= quantity.vertex_values[3,1]

        assert quantity.vertex_values[1,1] <= quantity.vertex_values[2,1]
        assert quantity.vertex_values[1,1] >= quantity.vertex_values[0,2]

        assert quantity.vertex_values[1,2] <= quantity.vertex_values[2,0]
        assert quantity.vertex_values[1,2] <= quantity.vertex_values[3,1]



        #Assert that quantities are conserved
        for k in range(quantity.centroid_values.shape[0]):
            assert num.allclose (quantity.centroid_values[k],
                                 num.sum(quantity.vertex_values[k,:])/3)



    def test_limit_edges_by_all_neighbours(self):
        quantity = Quantity(self.mesh4)

        #Create a deliberate overshoot (e.g. from gradient computation)
        quantity.set_values([[3,0,3], [2,2,6], [5,3,8], [8,3,5]])


        #Limit
        quantity.limit_edges_by_all_neighbours()

        #Assert that central triangle is limited by neighbours
        assert quantity.edge_values[1,0] <= quantity.centroid_values[2]
        assert quantity.edge_values[1,0] >= quantity.centroid_values[0]

        assert quantity.edge_values[1,1] <= quantity.centroid_values[2]
        assert quantity.edge_values[1,1] >= quantity.centroid_values[0]

        assert quantity.edge_values[1,2] <= quantity.centroid_values[2]
        assert quantity.edge_values[1,2] >= quantity.centroid_values[0]



        #Assert that quantities are conserved
        for k in range(quantity.centroid_values.shape[0]):
            assert num.allclose (quantity.centroid_values[k],
                                 num.sum(quantity.vertex_values[k,:])/3)


    def test_limit_edges_by_neighbour(self):
        quantity = Quantity(self.mesh4)

        #Create a deliberate overshoot (e.g. from gradient computation)
        quantity.set_values([[3,0,3], [2,2,6], [5,3,8], [8,3,5]])


        #Limit
        quantity.limit_edges_by_neighbour()

        #Assert that central triangle is limited by neighbours
        assert quantity.edge_values[1,0] <= quantity.centroid_values[3]
        assert quantity.edge_values[1,0] >= quantity.centroid_values[1]

        assert quantity.edge_values[1,1] <= quantity.centroid_values[2]
        assert quantity.edge_values[1,1] >= quantity.centroid_values[1]

        assert quantity.edge_values[1,2] <= quantity.centroid_values[1]
        assert quantity.edge_values[1,2] >= quantity.centroid_values[0]



        #Assert that quantities are conserved
        for k in range(quantity.centroid_values.shape[0]):
            assert num.allclose (quantity.centroid_values[k],
                                 num.sum(quantity.vertex_values[k,:])/3)

    def test_limiter2(self):
        """Taken from test_shallow_water
        """
        quantity = Quantity(self.mesh4)
        quantity.domain.beta_w = 0.9
        
        #Test centroids
        quantity.set_values([2.,4.,8.,2.], location = 'centroids')
        assert num.allclose(quantity.centroid_values, [2, 4, 8, 2]) #Centroid


        #Extrapolate
        quantity.extrapolate_second_order()

        assert num.allclose(quantity.vertex_values[1,:], [0.0, 6, 6])

        #Limit
        quantity.limit()

        # limited value for beta_w = 0.9
        
        assert num.allclose(quantity.vertex_values[1,:], [2.2, 4.9, 4.9])
        # limited values for beta_w = 0.5
        #assert allclose(quantity.vertex_values[1,:], [3.0, 4.5, 4.5])


        #Assert that quantities are conserved
        for k in range(quantity.centroid_values.shape[0]):
            assert num.allclose (quantity.centroid_values[k],
                                 num.sum(quantity.vertex_values[k,:])/3)





    def test_distribute_first_order(self):
        quantity = Quantity(self.mesh4)

        #Test centroids
        quantity.set_values([1.,2.,3.,4.], location = 'centroids')
        assert num.allclose(quantity.centroid_values, [1, 2, 3, 4]) #Centroid


        #Extrapolate from centroid to vertices and edges
        quantity.extrapolate_first_order()

        #Interpolate
        #quantity.interpolate_from_vertices_to_edges()

        assert num.allclose(quantity.vertex_values,
                            [[1,1,1], [2,2,2], [3,3,3], [4, 4, 4]])
        assert num.allclose(quantity.edge_values, [[1,1,1], [2,2,2],
                                                   [3,3,3], [4, 4, 4]])


    def test_interpolate_from_vertices_to_edges(self):
        quantity = Quantity(self.mesh4)

        quantity.vertex_values = num.array([[1,0,2], [1,2,4], [4,2,5], [3,1,4]], num.float)

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
                                          [2.5, 3.5, 2]], num.float)

        quantity.interpolate_from_edges_to_vertices()

        assert num.allclose(quantity.vertex_values,
                            [[1,0,2], [1,2,4], [4,2,5], [3,1,4]])



    def test_distribute_second_order(self):
        quantity = Quantity(self.mesh4)

        #Test centroids
        quantity.set_values([2.,4.,8.,2.], location = 'centroids')
        assert num.allclose(quantity.centroid_values, [2, 4, 8, 2]) #Centroid


        #Extrapolate
        quantity.extrapolate_second_order()

        assert num.allclose(quantity.vertex_values[1,:], [0.0, 6, 6])


    def test_update_explicit(self):
        quantity = Quantity(self.mesh4)

        #Test centroids
        quantity.set_values([1.,2.,3.,4.], location = 'centroids')
        assert num.allclose(quantity.centroid_values, [1, 2, 3, 4]) #Centroid

        #Set explicit_update
        quantity.explicit_update = num.array( [1.,1.,1.,1.] )

        #Update with given timestep
        quantity.update(0.1)

        x = num.array([1, 2, 3, 4]) + num.array( [.1,.1,.1,.1] )
        assert num.allclose( quantity.centroid_values, x)

    def test_update_semi_implicit(self):
        quantity = Quantity(self.mesh4)

        #Test centroids
        quantity.set_values([1.,2.,3.,4.], location = 'centroids')
        assert num.allclose(quantity.centroid_values, [1, 2, 3, 4]) #Centroid

        #Set semi implicit update
        quantity.semi_implicit_update = num.array([1.,1.,1.,1.])

        #Update with given timestep
        timestep = 0.1
        quantity.update(timestep)

        sem = num.array([1.,1.,1.,1.])/num.array([1, 2, 3, 4])
        denom = num.ones(4, num.float)-timestep*sem

        x = num.array([1, 2, 3, 4])/denom
        assert num.allclose( quantity.centroid_values, x)


    def test_both_updates(self):
        quantity = Quantity(self.mesh4)

        #Test centroids
        quantity.set_values([1.,2.,3.,4.], location = 'centroids')
        assert num.allclose(quantity.centroid_values, [1, 2, 3, 4]) #Centroid

        #Set explicit_update
        quantity.explicit_update = num.array( [4.,3.,2.,1.] )

        #Set semi implicit update
        quantity.semi_implicit_update = num.array( [1.,1.,1.,1.] )

        #Update with given timestep
        timestep = 0.1
        quantity.update(0.1)

        sem = num.array([1.,1.,1.,1.])/num.array([1, 2, 3, 4])
        denom = num.ones(4, num.float)-timestep*sem

        x = num.array([1., 2., 3., 4.])
        x += timestep*num.array( [4.0, 3.0, 2.0, 1.0] )        
        x /= denom


        assert num.allclose( quantity.centroid_values, x)



    def set_array_values_by_index(self):

        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        #Create basic mesh
        points, vertices, boundary = rectangular(1, 1)

        #Create shallow water domain
        domain = Generic_Domain(points, vertices, boundary)
        #print "domain.number_of_elements ",domain.number_of_elements
        quantity = Quantity(domain,[[1,1,1],[2,2,2]])
        value = [7]
        indices = [1]
        quantity.set_array_values_by_index(value,
                                           location = 'centroids',
                                           indices = indices)
        #print "quantity.centroid_values",quantity.centroid_values

        assert num.allclose(quantity.centroid_values, [1,7])

        quantity.set_array_values([15,20,25], indices = indices)
        assert num.allclose(quantity.centroid_values, [1,20])

        quantity.set_array_values([15,20,25], indices = indices)
        assert num.allclose(quantity.centroid_values, [1,20])

    def test_setting_some_vertex_values(self):
        """
        set values based on triangle lists.
        """
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        #Create basic mesh
        points, vertices, boundary = rectangular(1, 3)
        #print "vertices",vertices
        #Create shallow water domain
        domain = Generic_Domain(points, vertices, boundary)
        #print "domain.number_of_elements ",domain.number_of_elements
        quantity = Quantity(domain,[[1,1,1],[2,2,2],[3,3,3],
                                    [4,4,4],[5,5,5],[6,6,6]])


        # Check that constants work
        value = 7
        indices = [1]
        quantity.set_values(value,
                            location = 'centroids',
                            indices = indices)
        #print "quantity.centroid_values",quantity.centroid_values
        assert num.allclose(quantity.centroid_values, [1,7,3,4,5,6])
        
        value = [7]
        indices = [1]
        quantity.set_values(value,
                            location = 'centroids',
                            indices = indices)
        #print "quantity.centroid_values",quantity.centroid_values
        assert num.allclose(quantity.centroid_values, [1,7,3,4,5,6])

        value = [[15,20,25]]
        quantity.set_values(value, indices = indices)
        #print "1 quantity.vertex_values",quantity.vertex_values
        assert num.allclose(quantity.vertex_values[1], value[0])


        #print "quantity",quantity.vertex_values
        values = [10,100,50]
        quantity.set_values(values, indices = [0,1,5], location = 'centroids')
        #print "2 quantity.vertex_values",quantity.vertex_values
        assert num.allclose(quantity.vertex_values[0], [10,10,10])
        assert num.allclose(quantity.vertex_values[5], [50,50,50])
        #quantity.interpolate()
        #print "quantity.centroid_values",quantity.centroid_values
        assert num.allclose(quantity.centroid_values, [10,100,3,4,5,50])


        quantity = Quantity(domain,[[1,1,1],[2,2,2],[3,3,3],
                                    [4,4,4],[5,5,5],[6,6,6]])
        values = [10,100,50]
        #this will be per unique vertex, indexing the vertices
        #print "quantity.vertex_values",quantity.vertex_values
        quantity.set_values(values, indices = [0,1,5])
        #print "quantity.vertex_values",quantity.vertex_values
        assert num.allclose(quantity.vertex_values[0], [1,50,10])
        assert num.allclose(quantity.vertex_values[5], [6,6,6])
        assert num.allclose(quantity.vertex_values[1], [100,10,50])

        quantity = Quantity(domain,[[1,1,1],[2,2,2],[3,3,3],
                                    [4,4,4],[5,5,5],[6,6,6]])
        values = [[31,30,29],[400,400,400],[1000,999,998]]
        quantity.set_values(values, indices = [3,3,5])
        quantity.interpolate()
        assert num.allclose(quantity.centroid_values, [1,2,3,400,5,999])

        values = [[1,1,1],[2,2,2],[3,3,3],
                                    [4,4,4],[5,5,5],[6,6,6]]
        quantity.set_values(values)

        # testing the standard set values by vertex
        # indexed by vertex_id in general_mesh.coordinates
        values = [0,1,2,3,4,5,6,7]

        quantity.set_values(values)
        #print "1 quantity.vertex_values",quantity.vertex_values
        assert num.allclose(quantity.vertex_values,[[ 4.,  5.,  0.],
                                                    [ 1.,  0.,  5.],
                                                    [ 5.,  6.,  1.],
                                                    [ 2.,  1.,  6.],
                                                    [ 6.,  7.,  2.],
                                                    [ 3.,  2.,  7.]])

    def test_setting_unique_vertex_values(self):
        """
        set values based on unique_vertex lists.
        """
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        #Create basic mesh
        points, vertices, boundary = rectangular(1, 3)
        #print "vertices",vertices
        #Create shallow water domain
        domain = Generic_Domain(points, vertices, boundary)
        #print "domain.number_of_elements ",domain.number_of_elements
        quantity = Quantity(domain,[[0,0,0],[1,1,1],[2,2,2],[3,3,3],
                                    [4,4,4],[5,5,5]])
        value = 7
        indices = [1,5]
        quantity.set_values(value,
                            location='unique vertices',
                            indices=indices)
        #print "quantity.centroid_values",quantity.centroid_values
        assert num.allclose(quantity.vertex_values[0], [0,7,0])
        assert num.allclose(quantity.vertex_values[1], [7,1,7])
        assert num.allclose(quantity.vertex_values[2], [7,2,7])


    def test_get_values(self):
        """
        get values based on triangle lists.
        """
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        #Create basic mesh
        points, vertices, boundary = rectangular(1, 3)

        #print "points",points
        #print "vertices",vertices
        #print "boundary",boundary

        #Create shallow water domain
        domain = Generic_Domain(points, vertices, boundary)
        #print "domain.number_of_elements ",domain.number_of_elements
        quantity = Quantity(domain,[[0,0,0],[1,1,1],[2,2,2],[3,3,3],
                                    [4,4,4],[5,5,5]])

        #print "quantity.get_values(location = 'unique vertices')", \
        #      quantity.get_values(location = 'unique vertices')

        #print "quantity.get_values(location = 'unique vertices')", \
        #      quantity.get_values(indices=[0,1,2,3,4,5,6,7], \
        #                          location = 'unique vertices')

        answer = [0.5,2,4,5,0,1,3,4.5]
        assert num.allclose(answer,
                            quantity.get_values(location='unique vertices'))

        indices = [0,5,3]
        answer = [0.5,1,5]
        assert num.allclose(answer,
                            quantity.get_values(indices=indices,
                                                location='unique vertices'))
        #print "quantity.centroid_values",quantity.centroid_values
        #print "quantity.get_values(location = 'centroids') ",\
        #      quantity.get_values(location = 'centroids')




    def test_get_values_2(self):
        """Different mesh (working with domain object) - also check centroids.
        """

        
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Generic_Domain(points, vertices)

        quantity = Quantity(domain)
        quantity.set_values(lambda x, y: x+2*y) #2 4 4 6
        
        assert num.allclose(quantity.get_values(location='centroids'), 
                            [2,4,4,6])
                            
        assert num.allclose(quantity.get_values(location='centroids', 
                                                indices=[1,3]), [4,6])


        assert num.allclose(quantity.get_values(location='vertices'), 
                            [[4,0,2],
                             [4,2,6],
                             [6,2,4],
                             [8,4,6]])
        
        assert num.allclose(quantity.get_values(location='vertices', 
                                                indices=[1,3]), [[4,2,6],
                                                                 [8,4,6]])


        assert num.allclose(quantity.get_values(location='edges'), 
                            [[1,3,2],
                             [4,5,3],
                             [3,5,4],
                             [5,7,6]])
        assert num.allclose(quantity.get_values(location='edges', 
                                                indices=[1,3]),
                            [[4,5,3],
                             [5,7,6]])        

        # Check averaging over vertices
        #a: 0
        #b: (4+4+4)/3
        #c: (2+2+2)/3
        #d: 8
        #e: (6+6+6)/3        
        #f: 4
        assert num.allclose(quantity.get_values(location='unique vertices'), 
                            [0, 4, 2, 8, 6, 4]) 



    def test_get_interpolated_values(self):

        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        #Create basic mesh
        points, vertices, boundary = rectangular(1, 3)
        domain = Generic_Domain(points, vertices, boundary)

        #Constant values
        quantity = Quantity(domain,[[0,0,0],[1,1,1],[2,2,2],[3,3,3],
                                    [4,4,4],[5,5,5]])

        

        # Get interpolated values at centroids
        interpolation_points = domain.get_centroid_coordinates()
        answer = quantity.get_values(location='centroids')

        
        #print quantity.get_values(points=interpolation_points)
        assert num.allclose(answer, quantity.get_values(interpolation_points=interpolation_points))


        #Arbitrary values
        quantity = Quantity(domain,[[0,1,2],[3,1,7],[2,1,2],[3,3,7],
                                    [1,4,-9],[2,5,0]])


        # Get interpolated values at centroids
        interpolation_points = domain.get_centroid_coordinates()
        answer = quantity.get_values(location='centroids')
        #print answer
        #print quantity.get_values(interpolation_points=interpolation_points)
        assert num.allclose(answer, 
                            quantity.get_values(interpolation_points=interpolation_points,
                                                verbose=False))        
                        

        #FIXME TODO
        #indices = [0,5,3]
        #answer = [0.5,1,5]
        #assert allclose(answer,
        #                quantity.get_values(indices=indices, \
        #                                    location = 'unique vertices'))




    def test_get_interpolated_values_2(self):
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [[1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Generic_Domain(points, vertices)

        quantity = Quantity(domain)
        quantity.set_values(lambda x, y: x+2*y) #2 4 4 6

        #First pick one point
        x, y = 2.0/3, 8.0/3
        v = quantity.get_values(interpolation_points = [[x,y]])
        assert num.allclose(v, 6)        

        # Then another to test that algorithm won't blindly
        # reuse interpolation matrix
        x, y = 4.0/3, 4.0/3
        v = quantity.get_values(interpolation_points = [[x,y]])
        assert num.allclose(v, 4)        



    def test_get_interpolated_values_with_georef(self):
    
        zone = 56
        xllcorner = 308500
        yllcorner = 6189000
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [[1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Generic_Domain(points, vertices,
                        geo_reference=Geo_reference(zone,xllcorner,yllcorner))

        quantity = Quantity(domain)
        quantity.set_values(lambda x, y: x+2*y) #2 4 4 6

        #First pick one point (and turn it into absolute coordinates)
        x, y = 2.0/3, 8.0/3
        v = quantity.get_values(interpolation_points = [[x+xllcorner,y+yllcorner]])
        assert num.allclose(v, 6)
        

        # Then another to test that algorithm won't blindly
        # reuse interpolation matrix 
        x, y = 4.0/3, 4.0/3
        v = quantity.get_values(interpolation_points = [[x+xllcorner,y+yllcorner]])
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
                              geo_reference=Geo_reference(zone,xllcorner,yllcorner))
        v = quantity.get_values(interpolation_points=pts)
        assert num.allclose(v, [6, 4])                        
        
        
        

    def test_getting_some_vertex_values(self):
        """
        get values based on triangle lists.
        """
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        #Create basic mesh
        points, vertices, boundary = rectangular(1, 3)

        #print "points",points
        #print "vertices",vertices
        #print "boundary",boundary

        #Create shallow water domain
        domain = Generic_Domain(points, vertices, boundary)
        #print "domain.number_of_elements ",domain.number_of_elements
        quantity = Quantity(domain,[[1,1,1],[2,2,2],[3,3,3],
                                    [4,4,4],[5,5,5],[6,6,6]])
        value = [7]
        indices = [1]
        quantity.set_values(value,
                            location = 'centroids',
                            indices = indices)
        #print "quantity.centroid_values",quantity.centroid_values
        #print "quantity.get_values(location = 'centroids') ",\
        #      quantity.get_values(location = 'centroids')
        assert num.allclose(quantity.centroid_values,
                            quantity.get_values(location = 'centroids'))


        value = [[15,20,25]]
        quantity.set_values(value, indices = indices)
        #print "1 quantity.vertex_values",quantity.vertex_values
        assert num.allclose(quantity.vertex_values, quantity.get_values())

        assert num.allclose(quantity.edge_values,
                            quantity.get_values(location = 'edges'))

        # get a subset of elements
        subset = quantity.get_values(location='centroids', indices=[0,5])
        answer = [quantity.centroid_values[0],quantity.centroid_values[5]]
        assert num.allclose(subset, answer)


        subset = quantity.get_values(location='edges', indices=[0,5])
        answer = [quantity.edge_values[0],quantity.edge_values[5]]
        #print "subset",subset
        #print "answer",answer
        assert num.allclose(subset, answer)

        subset = quantity.get_values( indices=[1,5])
        answer = [quantity.vertex_values[1],quantity.vertex_values[5]]
        #print "subset",subset
        #print "answer",answer
        assert num.allclose(subset, answer)

    def test_smooth_vertex_values(self):
        """
        get values based on triangle lists.
        """
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular
      #  from anuga.shallow_water.shallow_water_domain import Domain

        #Create basic mesh
        points, vertices, boundary = rectangular(2, 2)

        #print "points",points
        #print "vertices",vertices
        #print "boundary",boundary

        #Create shallow water domain
        domain = Generic_Domain(points, vertices, boundary)
        #print "domain.number_of_elements ",domain.number_of_elements
        quantity = Quantity(domain,[[0,0,0],[1,1,1],[2,2,2],[3,3,3],
                                    [4,4,4],[5,5,5],[6,6,6],[7,7,7]])

        #print "quantity.get_values(location = 'unique vertices')", \
        #      quantity.get_values(location = 'unique vertices')

        #print "quantity.get_values(location = 'unique vertices')", \
        #      quantity.get_values(indices=[0,1,2,3,4,5,6,7], \
        #                          location = 'unique vertices')

        #print quantity.get_values(location = 'unique vertices')
        #print quantity.domain.number_of_triangles_per_node
        #print quantity.vertex_values
        
        #answer = [0.5, 2, 3, 3, 3.5, 4, 4, 5, 6.5]
        #assert allclose(answer,
        #                quantity.get_values(location = 'unique vertices'))

        quantity.smooth_vertex_values()

        #print quantity.vertex_values


        answer_vertex_values = [[3,3.5,0.5],[2,0.5,3.5],[3.5,4,2],[3,2,4],
                                [4,5,3],[3.5,3,5],[5,6.5,3.5],[4,3.5,6.5]]
        
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

        quantity.set_values([[1,2,3], [5,5,5], [0,0,9], [-6, 3, 3]],
                            location = 'vertices')
        assert num.allclose(quantity.vertex_values,
                            [[1,2,3], [5,5,5], [0,0,9], [-6, 3, 3]])

        assert id(vertex_values) == id(quantity.vertex_values)
        assert num.allclose(quantity.centroid_values, [2., 5., 3., 0.]) #Centroid
        assert num.allclose(quantity.edge_values, [[2.5, 2.0, 1.5],
                                                   [5., 5., 5.],
                                                   [4.5, 4.5, 0.],
                                                   [3.0, -1.5, -1.5]])



        other_quantity = Quantity(self.mesh4)

        other_quantity.set_values([[0,0,0], [1,1,6], [10,10,10], [0,0,4]],
                            location = 'vertices')

 
        #===============================
        quantity.maximum(other_quantity)
        #===============================

        exact_vertex_values = num.array([[  1.,   2.,   3.],
                                         [  5.,   5.,  6.],
                                         [ 10.,  10.,  10.],
                                         [  0.,   3.,  4.]])
        exact_centroid_values = num.array([  2.,  5., 10.,  1.33333333])
        exact_edge_values = num.array([[  2.5,   2.,    1.5],
                                       [  5.,   5.,   5., ],
                                       [ 10.,   10.,  10. ],
                                       [  3.,    2.,    0. ]])


        
        assert num.allclose(quantity.vertex_values,
                            exact_vertex_values)
        assert num.allclose(quantity.centroid_values, exact_centroid_values) #Centroid
        assert num.allclose(quantity.edge_values, exact_edge_values)

    def test_minimum(self):
        quantity = Quantity(self.mesh4)

        # get referece to data arrays
        centroid_values = quantity.centroid_values
        vertex_values = quantity.vertex_values
        edge_values = quantity.edge_values

        quantity.set_values([[1,2,3], [5,5,5], [0,0,9], [-6, 3, 3]],
                            location = 'vertices')
        assert num.allclose(quantity.vertex_values,
                            [[1,2,3], [5,5,5], [0,0,9], [-6, 3, 3]])

        assert id(vertex_values) == id(quantity.vertex_values)
        assert num.allclose(quantity.centroid_values, [2., 5., 3., 0.]) #Centroid
        assert num.allclose(quantity.edge_values, [[2.5, 2.0, 1.5],
                                                   [5., 5., 5.],
                                                   [4.5, 4.5, 0.],
                                                   [3.0, -1.5, -1.5]])



        other_quantity = Quantity(self.mesh4)

        other_quantity.set_values([[0,0,0], [1,1,6], [10,10,10], [0,0,4]],
                            location = 'vertices')


        #===============================
        quantity.minimum(other_quantity)
        #===============================

        exact_vertex_values = num.array([[  0.,   0.,   0.],
                                         [  1.,   1.,  5.],
                                         [  0.,   0.,  9.],
                                         [ -6.,   0.,  3.]])
        exact_centroid_values = num.array([  0.,  2.66666667, 3.,  0.])
        exact_edge_values = num.array([[  0.,   0.,   0.],
                                       [  3.5,   3.5,   1., ],
                                       [  4.5,   4.5,   0. ],
                                       [  2.,   -1.5,  -1.5 ]])


        assert num.allclose(quantity.vertex_values,
                            exact_vertex_values)
        assert num.allclose(quantity.centroid_values, exact_centroid_values) #Centroid
        assert num.allclose(quantity.edge_values, exact_edge_values)


#-------------------------------------------------------------

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_Quantity, 'test') #_set_values_from_asc')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)

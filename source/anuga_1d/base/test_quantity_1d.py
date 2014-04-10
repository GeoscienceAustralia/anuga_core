import quantity
#!/usr/bin/env python

import unittest
from math import sqrt, pi


from anuga_1d.base.generic_domain import Generic_domain as Domain
#from shallow_water_domain import flux_function as domain_flux_function

from anuga_1d.base.quantity import *



from numpy import allclose, array, ones, zeros
import numpy


class Test_Quantity(unittest.TestCase):
    def setUp(self):
        self.points3        = [0.0, 1.0, 2.0, 3.0]
        self.vertex_values3 = [[1.0,2.0],[4.0,5.0],[-1.0,2.0]]
        self.domain3        = Domain(self.points3)


        
        self.points4          = [0.0, 1.0, 2.5, 3.0, 5.0]
        self.vertex_values4   = [[1.0,2.0],[4.0,5.0],[-1.0,2.0],[3.0,6.0]]
        self.centroid_values4 = [1.5, 4.5, 0.5, 4.5]
        self.boundary4        = {(0, 0): 'left', (3, 1): 'right'}
        self.domain4          = Domain(self.points4,self.boundary4)

        self.points10 = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
        self.domain10 = Domain(self.points10)

        self.points6 = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        self.domain6 = Domain(self.points6)

    def tearDown(self):
        pass
        #print "  Tearing down"


    def test_creat_with_boundary(self):

        assert self.domain4.boundary  == {(0, 0): 'left', (3, 1): 'right'}

        
    def test_creation(self):

        quantity = Quantity(self.domain3)
        assert allclose(quantity.vertex_values, [[0.0,0.0],[0.0,0.0],[0.0,0.0]])

        
        try:
            quantity = Quantity()
        except:
            pass
        else:
            raise 'Should have raised empty quantity exception'


        try:
            quantity = Quantity([1,2,3])
        except AssertionError:
            pass
        else:
            raise 'Should have raised "mising domain object" error'


    def test_creation_zeros(self):

        quantity = Quantity(self.domain3)
        assert allclose(quantity.centroid_values, [[0.,0.,0.]])


        quantity = Quantity(self.domain4)
        assert allclose(quantity.vertex_values, [[0.,0.], [0.,0.],
                                                 [0.,0.], [0.,0.]])


    def test_interpolation(self):
        quantity = Quantity(self.domain4, self.vertex_values4)
        assert allclose(quantity.centroid_values, self.centroid_values4) #Centroid



    def test_interpolation2(self):
        quantity = Quantity(self.domain4, self.vertex_values4)
        assert allclose(quantity.centroid_values, self.centroid_values4) #Centroid

        quantity.extrapolate_second_order()

        #print quantity.vertex_values
        assert allclose(quantity.vertex_values,[[ 0.3,  2.7],
                                                [ 4.5,  4.5],
                                                [ 0.5,  0.5],
                                                [ 1.3,  7.7]])

 

 
    def test_boundary_allocation(self):
        quantity = Quantity(self.domain4,
                            [[1,2], [5,5], [0,9], [-6, 3]])


        assert quantity.boundary_values.shape[0] == len(self.domain4.boundary)


    def test_set_values(self):
        quantity = Quantity(self.domain4)

        # Test vertices
        quantity.set_values([[1,2], [5,5], [0,9], [-6, 3]], location = 'vertices')
        assert allclose(quantity.vertex_values,
                        [[1,2], [5,5], [0,9], [-6, 3]])
        assert allclose(quantity.centroid_values, [1.5, 5., 4.5, -1.5]) #Centroid



        # Test unique_vertices
        quantity.set_values([1,2,4,-5,6], location='unique_vertices')
        assert allclose(quantity.vertex_values,
                        [[1,2], [2,4], [4,-5], [-5,6]])
        assert allclose(quantity.centroid_values, [1.5, 3., -.5, .5]) #Centroid


        # Test centroids
        quantity.set_values([1,2,3,4], location = 'centroids')
        assert allclose(quantity.centroid_values, [1., 2., 3., 4.]) #Centroid

        # Test exceptions
        try:
            quantity.set_values([[1,3], [5,5], [0,9], [-6, 3]],
                                location = 'bas kamel tuba')
        except:
            pass


        try:
            quantity.set_values([[1,3], [0,9]])
        except AssertionError:
            pass
        except:
            raise 'should have raised Assertionerror'



    def test_set_values_const(self):
        quantity = Quantity(self.domain4)

        quantity.set_values(1.0, location = 'vertices')

        assert allclose(quantity.vertex_values, [[1,1], [1,1], [1,1], [1, 1]])
        assert allclose(quantity.centroid_values, [1, 1, 1, 1]) #Centroid


        quantity.set_values(2.0, location = 'centroids')

        assert allclose(quantity.centroid_values, [2, 2, 2, 2])


    def test_set_values_func(self):
        quantity = Quantity(self.domain4)

        def f(x):
            return x*x



        quantity.set_values(f, location = 'vertices')


        assert allclose(quantity.vertex_values,
                        [[0,1], [1,6.25], [6.25,9], [9,25]])

        assert allclose(quantity.centroid_values,
                        [0.5, 3.625, 7.625, 34*0.5])


        quantity.set_values(f, location = 'centroids')

        
        assert allclose(quantity.centroid_values,
                        [0.25, 3.0625, 7.5625, 16.0])


    def test_integral(self):
        quantity = Quantity(self.domain4)

        #Try constants first
        const = 5
        quantity.set_values(const, location = 'vertices')
        #print 'Q', quantity.get_integral()

        assert allclose(quantity.get_integral(), self.domain4.get_area() * const)

        #Try with a linear function
        def f(x):
            return x

        quantity.set_values(f, location = 'vertices')

 
        assert allclose (quantity.centroid_values,
        [ 0.5,   1.75,  2.75,  4.  ])

        assert allclose (quantity.vertex_values, [[ 0.,   1. ],
                                                  [ 1.,   2.5],
                                                  [ 2.5,  3. ],
                                                  [ 3.,   5. ]])


        ref_integral = 0.5 + 1.5*1.75 + 0.5*2.75 + 2.0*4.0

        assert allclose (quantity.get_integral(), ref_integral)





    def test_set_values_from_quantity(self):

        quantity1 = Quantity(self.domain4)
        quantity1.set_values([0,1,2,3,4], location='unique_vertices')

        assert allclose(quantity1.vertex_values,
                        [[0,1], [1,2], [2,3], [3,4]])


        quantity2 = Quantity(self.domain4)
        quantity2.set_values(quantity1)
        assert allclose(quantity2.vertex_values,
                        [[0,1], [1,2], [2,3], [3,4]])

        quantity2.set_values(2*quantity1)

        assert allclose(quantity2.vertex_values,
                        [[0,2], [2,4], [4,6], [6,8]])

        quantity2.set_values(2*quantity1 + 3)
        assert allclose(quantity2.vertex_values,
                        [[3,5], [5,7], [7,9], [9,11]])



    def test_overloading(self):

        quantity1 = Quantity(self.domain4)
        quantity1.set_values( [[0,1],[1,2],[2,3],[3,4]],
                            location = 'vertices')

        assert allclose(quantity1.vertex_values,
                        [[0,1], [1,2], [2,3], [3,4]])


        quantity2 = Quantity(self.domain4)
        quantity2.set_values([[1,2], [5,5], [0,9], [-6, 3]],
                             location = 'vertices')



        quantity3 = Quantity(self.domain4)
        quantity3.set_values([[2,2], [7,8], [7,6], [3, -8]],
                             location = 'vertices')


        # Negation
        Q = -quantity1
        assert allclose(Q.vertex_values, -quantity1.vertex_values)
        assert allclose(Q.centroid_values, -quantity1.centroid_values)


        # Addition
        Q = quantity1 + 7
        assert allclose(Q.vertex_values, quantity1.vertex_values + 7)
        assert allclose(Q.centroid_values, quantity1.centroid_values + 7)

        Q = 7 + quantity1
        assert allclose(Q.vertex_values, quantity1.vertex_values + 7)
        assert allclose(Q.centroid_values, quantity1.centroid_values + 7)

        Q = quantity1 + quantity2
        assert allclose(Q.vertex_values,
                        quantity1.vertex_values + quantity2.vertex_values)
        assert allclose(Q.centroid_values,
                        quantity1.centroid_values + quantity2.centroid_values)

        Q = quantity1 + quantity2 - 3
        assert allclose(Q.vertex_values,
                        quantity1.vertex_values + quantity2.vertex_values - 3)

        Q = quantity1 - quantity2
        assert allclose(Q.vertex_values,
                        quantity1.vertex_values - quantity2.vertex_values)

        #Scaling
        Q = quantity1*3
        assert allclose(Q.vertex_values, quantity1.vertex_values*3)
        assert allclose(Q.centroid_values, quantity1.centroid_values*3)

        Q = 3*quantity1
        assert allclose(Q.vertex_values, quantity1.vertex_values*3)

        #Multiplication
        Q = quantity1 * quantity2
        assert allclose(Q.vertex_values,
                        quantity1.vertex_values * quantity2.vertex_values)

        #Linear combinations
        Q = 4*quantity1 + 2
        assert allclose(Q.vertex_values,
                        4*quantity1.vertex_values + 2)

        Q = quantity1*quantity2 + 2
        assert allclose(Q.vertex_values,
                        quantity1.vertex_values * quantity2.vertex_values + 2)

        Q = quantity1*quantity2 + quantity3
        assert allclose(Q.vertex_values,
                        quantity1.vertex_values * quantity2.vertex_values +
                        quantity3.vertex_values)
        Q = quantity1*quantity2 + 3*quantity3
        assert allclose(Q.vertex_values,
                        quantity1.vertex_values * quantity2.vertex_values +
                        3*quantity3.vertex_values)
        Q = quantity1*quantity2 + 3*quantity3 + 5.0
        assert allclose(Q.vertex_values,
                        quantity1.vertex_values * quantity2.vertex_values +
                        3*quantity3.vertex_values + 5)

        Q = quantity1*quantity2 - quantity3
        assert allclose(Q.vertex_values,
                        quantity1.vertex_values * quantity2.vertex_values -
                        quantity3.vertex_values)
        Q = 1.5*quantity1*quantity2 - 3*quantity3 + 5.0
        assert allclose(Q.vertex_values,
                        1.5*quantity1.vertex_values * quantity2.vertex_values -
                        3*quantity3.vertex_values + 5)

        #Try combining quantities and arrays and scalars
        Q = 1.5*quantity1*quantity2.vertex_values -\
            3*quantity3.vertex_values + 5.0
        assert allclose(Q.vertex_values,
                        1.5*quantity1.vertex_values * quantity2.vertex_values -
                        3*quantity3.vertex_values + 5)


        #Powers
        Q = quantity1**2
        assert allclose(Q.vertex_values, quantity1.vertex_values**2)

        Q = quantity1**2 +quantity2**2
        assert allclose(Q.vertex_values,
                        quantity1.vertex_values**2 + \
                        quantity2.vertex_values**2)

        Q = (quantity1**2 +quantity2**2)**0.5
        assert allclose(Q.vertex_values,
                        (quantity1.vertex_values**2 + \
                        quantity2.vertex_values**2)**0.5)

    def test_compute_gradient(self):
        quantity = Quantity(self.domain6)

        #Set up for a gradient of (2,0) at mid triangle
        quantity.set_values([2.0, 4.0, 4.0, 5.0, 10.0, 12.0],
                            location = 'centroids')


        #Gradients
        quantity.compute_gradients()

        a = quantity.gradients

        assert allclose(a, [ 3., 1., 0.5, 3., 3.5, 0.5])

        quantity.beta = 2.0
        quantity.limiter = "minmod_kurganov"
        quantity.extrapolate_second_order()


        assert allclose(quantity.vertex_values, [[ 1., 3. ],
                                                 [ 4., 4. ],
                                                 [ 4., 4. ],
                                                 [ 4., 6.],
                                                 [ 8.25, 11.75],
                                                 [ 11.,  13. ]])

                                                 

    def test_second_order_extrapolation2(self):
        quantity = Quantity(self.domain4)

        #Set up for a gradient of (3,1), f(x) = 3x+y
        quantity.set_values([2.0+2.0/3, 4.0+4.0/3, 8.0+2.0/3, 2.0+8.0/3],
                            location = 'centroids')

        #Gradients
        quantity.compute_gradients()

        a = quantity.gradients
        

        assert allclose(a[1], 2.8)

        #Work out the others
        quantity.beta = 2.0
        quantity.limiter = "minmod_kurganov"
        quantity.extrapolate_second_order()
        

        assert allclose(quantity.vertex_values[1,0], 3.33333333)
        assert allclose(quantity.vertex_values[1,1], 7.33333333)

        assert allclose(quantity.centroid_values[1], 0.5*(7.33333333+3.33333333) )




    def test_backup_saxpy_centroid_values(self):
        quantity = Quantity(self.domain4)

        #Set up for a gradient of (3,1), f(x) = 3x+y
        c_values = array([2.0+2.0/3, 4.0+4.0/3, 8.0+2.0/3, 2.0+8.0/3])
        d_values = array([1.0, 2.0, 3.0, 4.0])
        quantity.set_values(c_values, location = 'centroids')

        #Backup
        quantity.backup_centroid_values()

        #print quantity.vertex_values
        assert allclose(quantity.centroid_values, quantity.centroid_backup_values)


        quantity.set_values(d_values, location = 'centroids')

        quantity.saxpy_centroid_values(2.0, 3.0)

        assert allclose(quantity.centroid_values, 2.0*d_values + 3.0*c_values)



    def test_first_order_extrapolator(self):
        quantity = Quantity(self.domain4)

        centroid_values = array([1.,2.,3.,4.])
        quantity.set_values(centroid_values, location = 'centroids')
        assert allclose(quantity.centroid_values, centroid_values) #Centroid

        #Extrapolate
        quantity.extrapolate_first_order()

        #Check that gradient is zero
        a = quantity.gradients
        assert allclose(a, [0,0,0,0])
       

        #Check vertices but not edge values
        assert allclose(quantity.vertex_values,
                        [[1,1], [2,2], [3,3], [4,4]])


    def test_second_order_extrapolator(self):
        quantity = Quantity(self.domain4)

        #Set up for a gradient of (3,0) at mid triangle
        quantity.set_values([2.0, 4.0, 8.0, 2.0],
                            location = 'centroids')



        quantity.extrapolate_second_order()


        #Assert that quantities are conserved
        for k in range(quantity.centroid_values.shape[0]):
            assert allclose (quantity.centroid_values[k],
                             numpy.sum(quantity.vertex_values[k,:])/2.0)


    def test_limit(self):
        quantity = Quantity(self.domain4)

        #Create a deliberate overshoot (e.g. from gradient computation)
        quantity.set_values([[0,0], [2,20], [-20,3], [8,3]])

        #Limit
        quantity.limit_minmod()

       
        #cells 1 and 2 are local max and min
        assert quantity.vertex_values[1][0] == quantity.centroid_values[1]
        assert quantity.vertex_values[1][1] == quantity.centroid_values[1]

        assert quantity.vertex_values[2][0] == quantity.centroid_values[2]
        assert quantity.vertex_values[2][1] == quantity.centroid_values[2]




    def test_limit_minmod(self):
        quantity = Quantity(self.domain4)

        #Create a deliberate overshoot (e.g. from gradient computation)
        quantity.set_values([[0,0], [2,10], [9,13], [12,14]])

        #Limit
        quantity.limiter = 'minmod'
        quantity.limit_range()




        assert allclose(quantity.vertex_values, [[ -2.4,   2.4],
                                                 [  2.4,   9.6],
                                                 [ 10.6,  11.4],
                                                 [ 11.4,  14.6]] )

        from anuga_1d.base.quantity_ext import limit_minmod_ext

        quantity.set_values([[0,0], [2,10], [9,13], [12,14]])

        limit_minmod_ext(quantity)




        assert allclose(quantity.vertex_values, [[ -2.4,   2.4],
                                                 [  2.4,   9.6],
                                                 [ 10.6,  11.4],
                                                 [ 11.4,  14.6]] )



    def test_limit_minmod_kurganov(self):
        quantity = Quantity(self.domain4)

        #Create a deliberate overshoot (e.g. from gradient computation)
        quantity.set_values([[0,0], [2,10], [9,13], [12,14]])

        #Limit
        quantity.limiter = 'minmod_kurganov'
        quantity.beta = 0.5
        quantity.limit_range()


        #print quantity.vertex_values

        assert allclose(quantity.vertex_values, [[ -2.4,   2.4],
                                                 [  4.2,   7.8],
                                                 [ 10.8,  11.2],
                                                 [ 11.4,  14.6]])

        from anuga_1d.base.quantity_ext import limit_minmod_kurganov_ext

        quantity.set_values([[0,0], [2,10], [9,13], [12,14]])

        limit_minmod_kurganov_ext(quantity)

        assert allclose(quantity.vertex_values, [[ -2.4,   2.4],
                                                 [  4.2,   7.8],
                                                 [ 10.8,  11.2],
                                                 [ 11.4,  14.6]])

        #print quantity.vertex_values


    def test_limit_vanleer(self):
        quantity = Quantity(self.domain4)

        #Create a deliberate overshoot (e.g. from gradient computation)
        quantity.set_values([[0,0], [2,10], [9,13], [12,14]])

        #Limit
        quantity.set_limiter('vanleer')
        quantity.limit_range()


        #print quantity.vertex_values

        assert allclose(quantity.vertex_values, [[ -2.4,          2.4       ],
                                                 [  2.32653061,   9.67346939],
                                                 [ 10.39393939,  11.60606061],
                                                 [ 11.4,         14.6       ]] )

        from anuga_1d.base.quantity_ext import limit_vanleer_ext

        quantity.set_values([[0,0], [2,10], [9,13], [12,14]])

        limit_vanleer_ext(quantity)


        #print quantity.vertex_values

        assert allclose(quantity.vertex_values, [[ -2.4,          2.4       ],
                                                 [  2.32653061,   9.67346939],
                                                 [ 10.39393939,  11.60606061],
                                                 [ 11.4,         14.6       ]] )

    def test_limit_vanalbada(self):
        quantity = Quantity(self.domain4)

        #Create a deliberate overshoot (e.g. from gradient computation)
        quantity.set_values([[0,0], [2,10], [9,13], [12,14]])

        #Limit
        quantity.limiter = 'vanalbada'
        quantity.limit_range()


        #print quantity.vertex_values

        assert allclose(quantity.vertex_values, [[ -2.4,          2.4       ],
                                                 [  2.32805995,   9.67194005],
                                                 [ 10.52104499,  11.47895501],
                                                 [ 11.4,         14.6       ]])

        from anuga_1d.base.quantity_ext import limit_vanalbada_ext

        quantity.set_values([[0,0], [2,10], [9,13], [12,14]])

        limit_vanalbada_ext(quantity)

        #print quantity.vertex_values




        assert allclose(quantity.vertex_values, [[ -2.4,          2.4       ],
                                                 [  2.32805995,   9.67194005],
                                                 [ 10.52104499,  11.47895501],
                                                 [ 11.4,         14.6       ]] )




    def test_find_qmax_qmin(self):
        quantity = Quantity(self.domain4)


        quantity.set_values([1.0, 4.0, 8.0, 2.0],
                            location = 'centroids')



        quantity.find_qmax_qmin()


        assert allclose(quantity.qmax, [4.0, 8.0, 8.0, 8.0])
        assert allclose(quantity.qmin, [1.0, 1.0, 2.0, 2.0])


    def test_distribute_first_order(self):
        quantity = Quantity(self.domain4)

        #Test centroids
        centroid_values = array([1.,2.,3.,4.])
        quantity.set_values(centroid_values, location = 'centroids')
        assert allclose(quantity.centroid_values, centroid_values) #Centroid


        #Extrapolate from centroid to vertices and edges
        quantity.extrapolate_first_order()
        
        assert allclose(quantity.vertex_values,[[ 1.,  1.],
                                                [ 2.,  2.],
                                                [ 3.,  3.],
                                                [ 4.,  4.]])
 


    def test_distribute_second_order(self):
        quantity = Quantity(self.domain4)

        #Test centroids
        centroid_values = array([2.,4.,8.,2.])
        quantity.set_values(centroid_values, location = 'centroids')
        assert allclose(quantity.centroid_values, centroid_values) #Centroid


        #Extrapolate
        quantity.beta = 2.0
        quantity.limiter = "minmod_kurganov"
        quantity.extrapolate_second_order()

        assert allclose(quantity.vertex_values, [[ 1.2,  2.8],
                                                 [ 2.,  6. ],
                                                 [ 8.,   8. ],
                                                 [ 6.8, -2.8]])


    def test_update_explicit(self):
        quantity = Quantity(self.domain4)

        #Test centroids
        quantity.set_values([1.,2.,3.,4.], location = 'centroids')
        assert allclose(quantity.centroid_values, [1, 2, 3, 4]) #Centroid

        #Set explicit_update
        quantity.explicit_update = array( [1.,1.,1.,1.] )

        #Update with given timestep
        quantity.update(0.1)

        x = array([1, 2, 3, 4]) + array( [.1,.1,.1,.1] )
        assert allclose( quantity.centroid_values, x)

    def test_update_semi_implicit(self):
        quantity = Quantity(self.domain4)

        #Test centroids
        quantity.set_values([1.,2.,3.,4.], location = 'centroids')
        assert allclose(quantity.centroid_values, [1, 2, 3, 4]) #Centroid

        #Set semi implicit update
        quantity.semi_implicit_update = array([1.,1.,1.,1.])

        #Update with given timestep
        timestep = 0.1
        quantity.update(timestep)

        sem = array([1.,1.,1.,1.])
        denom = ones(4, numpy.float)-timestep*sem

        x = array([1, 2, 3, 4])/denom
        assert allclose( quantity.centroid_values, x)


    def test_both_updates(self):
        quantity = Quantity(self.domain4)

        #Test centroids
        centroid_values = array( [1, 2, 3, 4] )
        quantity.set_values(centroid_values, location = 'centroids')
        assert allclose(quantity.centroid_values, centroid_values) #Centroid

        #Set explicit_update
        explicit_update = array( [4.,3.,2.,1.] )
        quantity.explicit_update[:,] = explicit_update

        #Set semi implicit update
        semi_implicit_update = array( [1.,1.,1.,1.] )
        quantity.semi_implicit_update[:,] = semi_implicit_update

        #Update with given timestep
        timestep = 0.1
        quantity.update(0.1)

        denom = 1.0-timestep*semi_implicit_update
        x = (centroid_values + timestep*explicit_update)/denom
 
        assert allclose( quantity.centroid_values, x)


 
#-------------------------------------------------------------
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_Quantity, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)

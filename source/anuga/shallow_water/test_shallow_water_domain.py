#!/usr/bin/env python

import unittest, os
from math import sqrt, pi
import tempfile

from anuga.config import g, epsilon
from Numeric import allclose, alltrue, array, zeros, ones, Float, take
from anuga.utilities.numerical_tools import mean
from anuga.utilities.polygon import is_inside_polygon
from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.abstract_2d_finite_volumes.quantity import Quantity
from anuga.geospatial_data.geospatial_data import Geospatial_data

from shallow_water_domain import *

# Get gateway to C implementation of flux function for direct testing
from shallow_water_ext import flux_function_central as flux_function

# For test_fitting_using_shallow_water_domain example
def linear_function(point):
    point = array(point)
    return point[:,0]+point[:,1]

class Weir:
    """Set a bathymetry for weir with a hole and a downstream gutter
    x,y are assumed to be in the unit square
    """

    def __init__(self, stage):
        self.inflow_stage = stage

    def __call__(self, x, y):
        from Numeric import zeros, Float

        N = len(x)
        assert N == len(y)

        z = zeros(N, Float)
        for i in range(N):
            z[i] = -x[i]/2  #General slope

            #Flattish bit to the left
            if x[i] < 0.3:
                z[i] = -x[i]/10

            #Weir
            if x[i] >= 0.3 and x[i] < 0.4:
                z[i] = -x[i]+0.9

            #Dip
            x0 = 0.6
            #depth = -1.3
            depth = -1.0
            #plateaux = -0.9
            plateaux = -0.6
            if y[i] < 0.7:
                if x[i] > x0 and x[i] < 0.9:
                    z[i] = depth

                #RHS plateaux
                if x[i] >= 0.9:
                    z[i] = plateaux


            elif y[i] >= 0.7 and y[i] < 1.5:
                #Restrict and deepen
                if x[i] >= x0 and x[i] < 0.8:
                    z[i] = depth-(y[i]/3-0.3)
                    #z[i] = depth-y[i]/5
                    #z[i] = depth
                elif x[i] >= 0.8:
                    #RHS plateaux
                    z[i] = plateaux

            elif y[i] >= 1.5:
                if x[i] >= x0 and x[i] < 0.8 + (y[i]-1.5)/1.2:
                    #Widen up and stay at constant depth
                    z[i] = depth-1.5/5
                elif x[i] >= 0.8 + (y[i]-1.5)/1.2:
                    #RHS plateaux
                    z[i] = plateaux


            #Hole in weir (slightly higher than inflow condition)
            if x[i] >= 0.3 and x[i] < 0.4 and y[i] > 0.2 and y[i] < 0.4:
                z[i] = -x[i]+self.inflow_stage + 0.02

            #Channel behind weir
            x0 = 0.5
            if x[i] >= 0.4 and x[i] < x0 and y[i] > 0.2 and y[i] < 0.4:
                z[i] = -x[i]+self.inflow_stage + 0.02

            if x[i] >= x0 and x[i] < 0.6 and y[i] > 0.2 and y[i] < 0.4:
                #Flatten it out towards the end
                z[i] = -x0+self.inflow_stage + 0.02 + (x0-x[i])/5

            #Hole to the east
            x0 = 1.1; y0 = 0.35
            #if x[i] < -0.2 and y < 0.5:
            if sqrt((2*(x[i]-x0))**2 + (2*(y[i]-y0))**2) < 0.2:
                z[i] = sqrt(((x[i]-x0))**2 + ((y[i]-y0))**2)-1.0

            #Tiny channel draining hole
            if x[i] >= 1.14 and x[i] < 1.2 and y[i] >= 0.4 and y[i] < 0.6:
                z[i] = -0.9 #North south

            if x[i] >= 0.9 and x[i] < 1.18 and y[i] >= 0.58 and y[i] < 0.65:
                z[i] = -1.0 + (x[i]-0.9)/3 #East west



            #Stuff not in use

            #Upward slope at inlet to the north west
            #if x[i] < 0.0: # and y[i] > 0.5:
            #    #z[i] = -y[i]+0.5  #-x[i]/2
            #    z[i] = x[i]/4 - y[i]**2 + 0.5

            #Hole to the west
            #x0 = -0.4; y0 = 0.35 # center
            #if sqrt((2*(x[i]-x0))**2 + (2*(y[i]-y0))**2) < 0.2:
            #    z[i] = sqrt(((x[i]-x0))**2 + ((y[i]-y0))**2)-0.2





        return z/2

class Weir_simple:
    """Set a bathymetry for weir with a hole and a downstream gutter
    x,y are assumed to be in the unit square
    """

    def __init__(self, stage):
        self.inflow_stage = stage

    def __call__(self, x, y):
        from Numeric import zeros, Float

        N = len(x)
        assert N == len(y)

        z = zeros(N, Float)
        for i in range(N):
            z[i] = -x[i]  #General slope

            #Flat bit to the left
            if x[i] < 0.3:
                z[i] = -x[i]/10  #General slope

            #Weir
            if x[i] > 0.3 and x[i] < 0.4:
                z[i] = -x[i]+0.9

            #Dip
            if x[i] > 0.6 and x[i] < 0.9:
                z[i] = -x[i]-0.5  #-y[i]/5

            #Hole in weir (slightly higher than inflow condition)
            if x[i] > 0.3 and x[i] < 0.4 and y[i] > 0.2 and y[i] < 0.4:
                z[i] = -x[i]+self.inflow_stage + 0.05


        return z/2




#Variable windfield implemented using functions
def speed(t,x,y):
    """Large speeds halfway between center and edges
    Low speeds at center and edges
    """

    from math import exp, cos, pi

    x = array(x)
    y = array(y)

    N = len(x)
    s = 0*x  #New array

    for k in range(N):

        r = sqrt(x[k]**2 + y[k]**2)

        factor = exp( -(r-0.15)**2 )

        s[k] = 4000 * factor * (cos(t*2*pi/150) + 2)

    return s


def scalar_func(t,x,y):
    """Function that returns a scalar.
    Used to test error message when Numeric array is expected
    """

    return 17.7


def angle(t,x,y):
    """Rotating field
    """
    from math import atan, pi

    x = array(x)
    y = array(y)

    N = len(x)
    a = 0*x  #New array

    for k in range(N):
        r = sqrt(x[k]**2 + y[k]**2)

        angle = atan(y[k]/x[k])

        if x[k] < 0:
            angle+=pi #atan in ]-pi/2; pi/2[

        #Take normal direction
        angle -= pi/2

        #Ensure positive radians
        if angle < 0:
            angle += 2*pi

        a[k] = angle/pi*180

    return a


class Test_Shallow_Water(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_rotate(self):
        normal = array([0.0,-1.0])

        q = array([1.0,2.0,3.0])

        r = rotate(q, normal, direction = 1)
        assert r[0] == 1
        assert r[1] == -3
        assert r[2] == 2

        w = rotate(r, normal, direction = -1)
        assert allclose(w, q)

        #Check error check
        try:
            rotate(r, array([1,1,1]) )
        except:
            pass
        else:
            raise 'Should have raised an exception'


    # Individual flux tests
    def test_flux_zero_case(self):
        ql = zeros( 3, Float )
        qr = zeros( 3, Float )
        normal = zeros( 2, Float )
        edgeflux = zeros( 3, Float )
        zl = zr = 0.
        H0 = 0.0
        
        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux, epsilon, g, H0)

        assert allclose(edgeflux, [0,0,0])
        assert max_speed == 0.

    def test_flux_constants(self):
        w = 2.0

        normal = array([1.,0])
        ql = array([w, 0, 0])
        qr = array([w, 0, 0])
        edgeflux = zeros(3, Float)        
        zl = zr = 0.
        h = w - (zl+zr)/2
        H0 = 0.0

        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux, epsilon, g, H0)        
        assert allclose(edgeflux, [0., 0.5*g*h**2, 0.])
        assert max_speed == sqrt(g*h)

    #def test_flux_slope(self):
    #    #FIXME: TODO
    #    w = 2.0
    #
    #    normal = array([1.,0])
    #    ql = array([w, 0, 0])
    #    qr = array([w, 0, 0])
    #    zl = zr = 0.
    #    h = w - (zl+zr)/2
    #
    #    flux, max_speed = flux_function(normal, ql, qr, zl, zr)
    #
    #    assert allclose(flux, [0., 0.5*g*h**2, 0.])
    #    assert max_speed == sqrt(g*h)


    def test_flux1(self):
        #Use data from previous version of abstract_2d_finite_volumes
        normal = array([1.,0])
        ql = array([-0.2, 2, 3])
        qr = array([-0.2, 2, 3])
        zl = zr = -0.5
        edgeflux = zeros(3, Float)                

        H0 = 0.0

        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux, epsilon, g, H0)        

        assert allclose(edgeflux, [2.,13.77433333, 20.])
        assert allclose(max_speed, 8.38130948661)


    def test_flux2(self):
        #Use data from previous version of abstract_2d_finite_volumes
        normal = array([0., -1.])
        ql = array([-0.075, 2, 3])
        qr = array([-0.075, 2, 3])
        zl = zr = -0.375

        edgeflux = zeros(3, Float)                
        H0 = 0.0
        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux, epsilon, g, H0)        

        assert allclose(edgeflux, [-3.,-20.0, -30.441])
        assert allclose(max_speed, 11.7146428199)

    def test_flux3(self):
        #Use data from previous version of abstract_2d_finite_volumes
        normal = array([-sqrt(2)/2, sqrt(2)/2])
        ql = array([-0.075, 2, 3])
        qr = array([-0.075, 2, 3])
        zl = zr = -0.375

        edgeflux = zeros(3, Float)                
        H0 = 0.0
        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux, epsilon, g, H0)        

        assert allclose(edgeflux, [sqrt(2)/2, 4.40221112, 7.3829019])
        assert allclose(max_speed, 4.0716654239)

    def test_flux4(self):
        #Use data from previous version of abstract_2d_finite_volumes
        normal = array([-sqrt(2)/2, sqrt(2)/2])
        ql = array([-0.34319278, 0.10254161, 0.07273855])
        qr = array([-0.30683287, 0.1071986, 0.05930515])
        zl = zr = -0.375

        edgeflux = zeros(3, Float)                
        H0 = 0.0
        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux, epsilon, g, H0)                

        assert allclose(edgeflux, [-0.04072676, -0.07096636, -0.01604364])
        assert allclose(max_speed, 1.31414103233)

    def test_flux_computation(self):    
        """test_flux_computation - test flux calculation (actual C implementation)
        This one tests the constant case where only the pressure term contributes to each edge and cancels out 
        once the total flux has been summed up.
        """
                
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe, daf, dae
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices)
        domain.check_integrity()

        # The constant case             
        domain.set_quantity('elevation', -1)
        domain.set_quantity('stage', 1) 
        
        domain.compute_fluxes()
        assert allclose(domain.get_quantity('stage').explicit_update[1], 0) # Central triangle
        

        # The more general case                 
        def surface(x,y):
            return -x/2                    
        
        domain.set_quantity('elevation', -10)
        domain.set_quantity('stage', surface)   
        domain.set_quantity('xmomentum', 1)             
        
        domain.compute_fluxes()
        
        #print domain.get_quantity('stage').explicit_update
        # FIXME (Ole): TODO the general case
        #assert allclose(domain.get_quantity('stage').explicit_update[1], ........??)
                
        
                
    def test_sw_domain_simple(self):
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe, daf, dae
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]


        #from anuga.abstract_2d_finite_volumes.domain import Domain as Generic_domain
        #msg = 'The class %s is not a subclass of the generic domain class %s'\
        #      %(DomainClass, Domain)
        #assert issubclass(DomainClass, Domain), msg

        domain = Domain(points, vertices)
        domain.check_integrity()

        for name in ['stage', 'xmomentum', 'ymomentum',
                     'elevation', 'friction']:
            assert domain.quantities.has_key(name)


        assert domain.get_conserved_quantities(0, edge=1) == 0.


    def test_boundary_conditions(self):

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4] ]
        boundary = { (0, 0): 'Third',
                     (0, 2): 'First',
                     (2, 0): 'Second',
                     (2, 1): 'Second',
                     (3, 1): 'Second',
                     (3, 2): 'Third'}


        domain = Domain(points, vertices, boundary)
        domain.check_integrity()


        domain.set_quantity('stage', [[1,2,3], [5,5,5],
                                      [0,0,9], [-6, 3, 3]])

        domain.set_quantity('xmomentum', [[1,1,1], [2,2,2],
                                          [3,3,3], [4, 4, 4]])

        domain.set_quantity('ymomentum', [[10,10,10], [20,20,20],
                                          [30,30,30], [40, 40, 40]])


        D = Dirichlet_boundary([5,2,1])
        T = Transmissive_boundary(domain)
        R = Reflective_boundary(domain)
        domain.set_boundary( {'First': D, 'Second': T, 'Third': R})

        domain.update_boundary()

        #Stage
        assert domain.quantities['stage'].boundary_values[0] == 2.5
        assert domain.quantities['stage'].boundary_values[0] ==\
               domain.get_conserved_quantities(0, edge=0)[0] #Reflective (2.5)
        assert domain.quantities['stage'].boundary_values[1] == 5. #Dirichlet
        assert domain.quantities['stage'].boundary_values[2] ==\
               domain.get_conserved_quantities(2, edge=0)[0] #Transmissive (4.5)
        assert domain.quantities['stage'].boundary_values[3] ==\
               domain.get_conserved_quantities(2, edge=1)[0] #Transmissive (4.5)
        assert domain.quantities['stage'].boundary_values[4] ==\
               domain.get_conserved_quantities(3, edge=1)[0] #Transmissive (-1.5)
        assert domain.quantities['stage'].boundary_values[5] ==\
               domain.get_conserved_quantities(3, edge=2)[0] #Reflective (-1.5)

        #Xmomentum
        assert domain.quantities['xmomentum'].boundary_values[0] == 1.0 #Reflective
        assert domain.quantities['xmomentum'].boundary_values[1] == 2. #Dirichlet
        assert domain.quantities['xmomentum'].boundary_values[2] ==\
               domain.get_conserved_quantities(2, edge=0)[1] #Transmissive
        assert domain.quantities['xmomentum'].boundary_values[3] ==\
               domain.get_conserved_quantities(2, edge=1)[1] #Transmissive
        assert domain.quantities['xmomentum'].boundary_values[4] ==\
               domain.get_conserved_quantities(3, edge=1)[1] #Transmissive
        assert domain.quantities['xmomentum'].boundary_values[5] == -4.0  #Reflective

        #Ymomentum
        assert domain.quantities['ymomentum'].boundary_values[0] == -10.0 #Reflective
        assert domain.quantities['ymomentum'].boundary_values[1] == 1.  #Dirichlet
        assert domain.quantities['ymomentum'].boundary_values[2] == 30. #Transmissive
        assert domain.quantities['ymomentum'].boundary_values[3] == 30. #Transmissive
        assert domain.quantities['ymomentum'].boundary_values[4] == 40. #Transmissive
        assert domain.quantities['ymomentum'].boundary_values[5] == 40. #Reflective


    def test_boundary_conditionsII(self):

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4] ]
        boundary = { (0, 0): 'Third',
                     (0, 2): 'First',
                     (2, 0): 'Second',
                     (2, 1): 'Second',
                     (3, 1): 'Second',
                     (3, 2): 'Third',
                     (0, 1): 'Internal'}


        domain = Domain(points, vertices, boundary)
        domain.check_integrity()


        domain.set_quantity('stage', [[1,2,3], [5,5,5],
                                      [0,0,9], [-6, 3, 3]])

        domain.set_quantity('xmomentum', [[1,1,1], [2,2,2],
                                          [3,3,3], [4, 4, 4]])

        domain.set_quantity('ymomentum', [[10,10,10], [20,20,20],
                                          [30,30,30], [40, 40, 40]])


        D = Dirichlet_boundary([5,2,1])
        T = Transmissive_boundary(domain)
        R = Reflective_boundary(domain)
        domain.set_boundary( {'First': D, 'Second': T,
                              'Third': R, 'Internal': None})

        domain.update_boundary()
        domain.check_integrity()


    def test_compute_fluxes0(self):
        # Do a full triangle and check that fluxes cancel out for
        # the constant stage case

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices)
        domain.set_quantity('stage', [[2,2,2], [2,2,2],
                                      [2,2,2], [2,2,2]])
        domain.check_integrity()

        assert allclose(domain.neighbours, [[-1,1,-1], [2,3,0], [-1,-1,1],[1,-1,-1]])
        assert allclose(domain.neighbour_edges, [[-1,2,-1], [2,0,1], [-1,-1,0],[1,-1,-1]])

        zl=zr=0. # Assume flat bed

        edgeflux = zeros(3, Float)        
        edgeflux0 = zeros(3, Float)
        edgeflux1 = zeros(3, Float)
        edgeflux2 = zeros(3, Float)                                
        H0 = 0.0        

        # Flux across right edge of volume 1
        normal = domain.get_normal(1,0)
        ql = domain.get_conserved_quantities(vol_id=1, edge=0)
        qr = domain.get_conserved_quantities(vol_id=2, edge=2)
        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux0, epsilon, g, H0)                        

        # Check that flux seen from other triangles is inverse
        tmp = qr; qr=ql; ql=tmp
        normal = domain.get_normal(2,2)
        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux, epsilon, g, H0)                                

        assert allclose(edgeflux0 + edgeflux, 0.)

        # Flux across upper edge of volume 1
        normal = domain.get_normal(1,1)
        ql = domain.get_conserved_quantities(vol_id=1, edge=1)
        qr = domain.get_conserved_quantities(vol_id=3, edge=0)
        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux1, epsilon, g, H0)                                        

        # Check that flux seen from other triangles is inverse
        tmp = qr; qr=ql; ql=tmp
        normal = domain.get_normal(3,0)
        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux, epsilon, g, H0)                                               

        assert allclose(edgeflux1 + edgeflux, 0.)        
        

        # Flux across lower left hypotenuse of volume 1
        normal = domain.get_normal(1,2)
        ql = domain.get_conserved_quantities(vol_id=1, edge=2)
        qr = domain.get_conserved_quantities(vol_id=0, edge=1)
        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux2, epsilon, g, H0)                                                               

        # Check that flux seen from other triangles is inverse
        tmp = qr; qr=ql; ql=tmp
        normal = domain.get_normal(0,1)
        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux, epsilon, g, H0)                                                       
        assert allclose(edgeflux2 + edgeflux, 0.)


        # Scale by edgelengths, add up anc check that total flux is zero
        e0 = domain.edgelengths[1, 0]
        e1 = domain.edgelengths[1, 1]
        e2 = domain.edgelengths[1, 2]

        assert allclose(e0*edgeflux0+e1*edgeflux1+e2*edgeflux2, 0.)

        # Now check that compute_flux yields zeros as well
        domain.compute_fluxes()

        for name in ['stage', 'xmomentum', 'ymomentum']:
            #print name, domain.quantities[name].explicit_update
            assert allclose(domain.quantities[name].explicit_update[1], 0)



    def test_compute_fluxes1(self):
        #Use values from previous version

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices)
        val0 = 2.+2.0/3
        val1 = 4.+4.0/3
        val2 = 8.+2.0/3
        val3 = 2.+8.0/3

        domain.set_quantity('stage', [[val0, val0, val0], [val1, val1, val1],
                                      [val2, val2, val2], [val3, val3, val3]])
        domain.check_integrity()

        zl=zr=0. #Assume flat bed

        edgeflux = zeros(3, Float)        
        edgeflux0 = zeros(3, Float)
        edgeflux1 = zeros(3, Float)
        edgeflux2 = zeros(3, Float)                                
        H0 = 0.0        
        

        # Flux across right edge of volume 1
        normal = domain.get_normal(1,0) #Get normal 0 of triangle 1
        assert allclose(normal, [1, 0])
        
        ql = domain.get_conserved_quantities(vol_id=1, edge=0)
        assert allclose(ql, [val1, 0, 0])
        
        qr = domain.get_conserved_quantities(vol_id=2, edge=2)
        assert allclose(qr, [val2, 0, 0])

        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux0, epsilon, g, H0)                                                       

        # Flux across edge in the east direction (as per normal vector)
        assert allclose(edgeflux0, [-15.3598804, 253.71111111, 0.])
        assert allclose(max_speed, 9.21592824046)


        #Flux across edge in the west direction (opposite sign for xmomentum)
        normal_opposite = domain.get_normal(2,2) #Get normal 2 of triangle 2
        assert allclose(normal_opposite, [-1, 0])

        max_speed = flux_function(normal_opposite, ql, qr, zl, zr, edgeflux, epsilon, g, H0)                                             
        #flux_opposite, max_speed = flux_function([-1, 0], ql, qr, zl, zr)
        assert allclose(edgeflux, [-15.3598804, -253.71111111, 0.])
        

        #Flux across upper edge of volume 1
        normal = domain.get_normal(1,1)
        ql = domain.get_conserved_quantities(vol_id=1, edge=1)
        qr = domain.get_conserved_quantities(vol_id=3, edge=0)
        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux1, epsilon, g, H0)                                                               

        assert allclose(edgeflux1, [2.4098563, 0., 123.04444444])
        assert allclose(max_speed, 7.22956891292)

        #Flux across lower left hypotenuse of volume 1
        normal = domain.get_normal(1,2)
        ql = domain.get_conserved_quantities(vol_id=1, edge=2)
        qr = domain.get_conserved_quantities(vol_id=0, edge=1)
        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux2, epsilon, g, H0)        

        assert allclose(edgeflux2, [9.63942522, -61.59685738, -61.59685738])
        assert allclose(max_speed, 7.22956891292)

        #Scale, add up and check that compute_fluxes is correct for vol 1
        e0 = domain.edgelengths[1, 0]
        e1 = domain.edgelengths[1, 1]
        e2 = domain.edgelengths[1, 2]

        total_flux = -(e0*edgeflux0+e1*edgeflux1+e2*edgeflux2)/domain.areas[1]
        assert allclose(total_flux, [-0.68218178, -166.6, -35.93333333])


        domain.compute_fluxes()

        #assert allclose(total_flux, domain.explicit_update[1,:])
        for i, name in enumerate(['stage', 'xmomentum', 'ymomentum']):
            assert allclose(total_flux[i],
                            domain.quantities[name].explicit_update[1])

        #assert allclose(domain.explicit_update, [
        #    [0., -69.68888889, -69.68888889],
        #    [-0.68218178, -166.6, -35.93333333],
        #    [-111.77316251, 69.68888889, 0.],
        #    [-35.68522449, 0., 69.68888889]])

        assert allclose(domain.quantities['stage'].explicit_update,
                        [0., -0.68218178, -111.77316251, -35.68522449])
        assert allclose(domain.quantities['xmomentum'].explicit_update,
                        [-69.68888889, -166.6, 69.68888889, 0])
        assert allclose(domain.quantities['ymomentum'].explicit_update,
                        [-69.68888889, -35.93333333, 0., 69.68888889])


        #assert allclose(domain.quantities[name].explicit_update





    def test_compute_fluxes2(self):
        #Random values, incl momentum

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices)
        val0 = 2.+2.0/3
        val1 = 4.+4.0/3
        val2 = 8.+2.0/3
        val3 = 2.+8.0/3

        zl=zr=0 #Assume flat zero bed
        edgeflux = zeros(3, Float)        
        edgeflux0 = zeros(3, Float)
        edgeflux1 = zeros(3, Float)
        edgeflux2 = zeros(3, Float)                                
        H0 = 0.0        
        

        domain.set_quantity('elevation', zl*ones( (4,3) ))


        domain.set_quantity('stage', [[val0, val0-1, val0-2],
                                      [val1, val1+1, val1],
                                      [val2, val2-2, val2],
                                      [val3-0.5, val3, val3]])

        domain.set_quantity('xmomentum',
                            [[1, 2, 3], [3, 4, 5],
                             [1, -1, 0], [0, -2, 2]])

        domain.set_quantity('ymomentum',
                            [[1, -1, 0], [0, -3, 2],
                             [0, 1, 0], [-1, 2, 2]])


        domain.check_integrity()



        #Flux across right edge of volume 1
        normal = domain.get_normal(1,0)
        ql = domain.get_conserved_quantities(vol_id=1, edge=0)
        qr = domain.get_conserved_quantities(vol_id=2, edge=2)
        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux0, epsilon, g, H0)                

        #Flux across upper edge of volume 1
        normal = domain.get_normal(1,1)
        ql = domain.get_conserved_quantities(vol_id=1, edge=1)
        qr = domain.get_conserved_quantities(vol_id=3, edge=0)
        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux1, epsilon, g, H0)                        

        #Flux across lower left hypotenuse of volume 1
        normal = domain.get_normal(1,2)
        ql = domain.get_conserved_quantities(vol_id=1, edge=2)
        qr = domain.get_conserved_quantities(vol_id=0, edge=1)
        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux2, epsilon, g, H0)                

        #Scale, add up and check that compute_fluxes is correct for vol 1
        e0 = domain.edgelengths[1, 0]
        e1 = domain.edgelengths[1, 1]
        e2 = domain.edgelengths[1, 2]

        total_flux = -(e0*edgeflux0+e1*edgeflux1+e2*edgeflux2)/domain.areas[1]


        domain.compute_fluxes()
        for i, name in enumerate(['stage', 'xmomentum', 'ymomentum']):
            assert allclose(total_flux[i],
                            domain.quantities[name].explicit_update[1])
        #assert allclose(total_flux, domain.explicit_update[1,:])


    # FIXME (Ole): Need test like this for fluxes in very shallow water.    
    def test_compute_fluxes3(self):
        #Random values, incl momentum

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices)
        val0 = 2.+2.0/3
        val1 = 4.+4.0/3
        val2 = 8.+2.0/3
        val3 = 2.+8.0/3

        zl=zr=-3.75 #Assume constant bed (must be less than stage)
        domain.set_quantity('elevation', zl*ones( (4,3) ))


        edgeflux = zeros(3, Float)        
        edgeflux0 = zeros(3, Float)
        edgeflux1 = zeros(3, Float)
        edgeflux2 = zeros(3, Float)                                
        H0 = 0.0        
        


        domain.set_quantity('stage', [[val0, val0-1, val0-2],
                                      [val1, val1+1, val1],
                                      [val2, val2-2, val2],
                                      [val3-0.5, val3, val3]])

        domain.set_quantity('xmomentum',
                            [[1, 2, 3], [3, 4, 5],
                             [1, -1, 0], [0, -2, 2]])

        domain.set_quantity('ymomentum',
                            [[1, -1, 0], [0, -3, 2],
                             [0, 1, 0], [-1, 2, 2]])


        domain.check_integrity()



        #Flux across right edge of volume 1
        normal = domain.get_normal(1,0)
        ql = domain.get_conserved_quantities(vol_id=1, edge=0)
        qr = domain.get_conserved_quantities(vol_id=2, edge=2)
        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux0, epsilon, g, H0)

        #Flux across upper edge of volume 1
        normal = domain.get_normal(1,1)
        ql = domain.get_conserved_quantities(vol_id=1, edge=1)
        qr = domain.get_conserved_quantities(vol_id=3, edge=0)
        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux1, epsilon, g, H0)        

        #Flux across lower left hypotenuse of volume 1
        normal = domain.get_normal(1,2)
        ql = domain.get_conserved_quantities(vol_id=1, edge=2)
        qr = domain.get_conserved_quantities(vol_id=0, edge=1)
        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux2, epsilon, g, H0)        

        #Scale, add up and check that compute_fluxes is correct for vol 1
        e0 = domain.edgelengths[1, 0]
        e1 = domain.edgelengths[1, 1]
        e2 = domain.edgelengths[1, 2]

        total_flux = -(e0*edgeflux0+e1*edgeflux1+e2*edgeflux2)/domain.areas[1]

        domain.compute_fluxes()
        for i, name in enumerate(['stage', 'xmomentum', 'ymomentum']):
            assert allclose(total_flux[i],
                            domain.quantities[name].explicit_update[1])



    def xtest_catching_negative_heights(self):

        #OBSOLETE
        
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices)
        val0 = 2.+2.0/3
        val1 = 4.+4.0/3
        val2 = 8.+2.0/3
        val3 = 2.+8.0/3

        zl=zr=4  #Too large
        domain.set_quantity('elevation', zl*ones( (4,3) ))
        domain.set_quantity('stage', [[val0, val0-1, val0-2],
                                      [val1, val1+1, val1],
                                      [val2, val2-2, val2],
                                      [val3-0.5, val3, val3]])

        #Should fail
        try:
            domain.check_integrity()
        except:
            pass



    def test_get_wet_elements(self):

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices)
        val0 = 2.+2.0/3
        val1 = 4.+4.0/3
        val2 = 8.+2.0/3
        val3 = 2.+8.0/3

        zl=zr=5
        domain.set_quantity('elevation', zl*ones( (4,3) ))
        domain.set_quantity('stage', [[val0, val0-1, val0-2],
                                      [val1, val1+1, val1],
                                      [val2, val2-2, val2],
                                      [val3-0.5, val3, val3]])



        #print domain.get_quantity('elevation').get_values(location='centroids')
        #print domain.get_quantity('stage').get_values(location='centroids')        
        domain.check_integrity()

        indices = domain.get_wet_elements()
        assert allclose(indices, [1,2])

        indices = domain.get_wet_elements(indices=[0,1,3])
        assert allclose(indices, [1])
        


    def test_get_maximum_inundation_1(self):

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices)

        domain.set_quantity('elevation', lambda x, y: x+2*y) #2 4 4 6
        domain.set_quantity('stage', 3)

        domain.check_integrity()

        indices = domain.get_wet_elements()
        assert allclose(indices, [0])

        q = domain.get_maximum_inundation_elevation()
        assert allclose(q, domain.get_quantity('elevation').get_values(location='centroids')[0])

        x, y = domain.get_maximum_inundation_location()
        assert allclose([x, y], domain.get_centroid_coordinates()[0])


    def test_get_maximum_inundation_2(self):
        """test_get_maximum_inundation_2(self)

        Test multiple wet cells with same elevation
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

        domain = Domain(points, vertices)

        domain.set_quantity('elevation', lambda x, y: x+2*y) #2 4 4 6
        domain.set_quantity('stage', 4.1)

        domain.check_integrity()

        indices = domain.get_wet_elements()
        assert allclose(indices, [0,1,2])

        q = domain.get_maximum_inundation_elevation()
        assert allclose(q, 4)        

        x, y = domain.get_maximum_inundation_location()
        assert allclose([x, y], domain.get_centroid_coordinates()[1])        
        

    def test_get_maximum_inundation_3(self):
        """test_get_maximum_inundation_3(self)

        Test of real runup example:
        """

        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross

        initial_runup_height = -0.4
        final_runup_height = -0.3


        #--------------------------------------------------------------
        # Setup computational domain
        #--------------------------------------------------------------
        N = 5
        points, vertices, boundary = rectangular_cross(N, N) 
        domain = Domain(points, vertices, boundary)
        domain.set_maximum_allowed_speed(1.0)        

        #--------------------------------------------------------------
        # Setup initial conditions
        #--------------------------------------------------------------
        def topography(x,y):
            return -x/2                             # linear bed slope
            

        domain.set_quantity('elevation', topography)       # Use function for elevation
        domain.set_quantity('friction', 0.)                # Zero friction 
        domain.set_quantity('stage', initial_runup_height) # Constant negative initial stage


        #--------------------------------------------------------------
        # Setup boundary conditions
        #--------------------------------------------------------------
        Br = Reflective_boundary(domain)              # Reflective wall
        Bd = Dirichlet_boundary([final_runup_height,  # Constant inflow
                                 0,
                                 0])

        # All reflective to begin with (still water) 
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})


        #--------------------------------------------------------------
        # Test initial inundation height
        #--------------------------------------------------------------

        indices = domain.get_wet_elements()
        z = domain.get_quantity('elevation').\
            get_values(location='centroids', indices=indices)
        assert alltrue(z<initial_runup_height)

        q = domain.get_maximum_inundation_elevation()
        assert allclose(q, initial_runup_height, rtol = 1.0/N) # First order accuracy 

        x, y = domain.get_maximum_inundation_location()

        qref = domain.get_quantity('elevation').get_values(interpolation_points = [[x, y]])
        assert allclose(q, qref)


        wet_elements = domain.get_wet_elements()
        wet_elevations = domain.get_quantity('elevation').get_values(location='centroids',
                                                                     indices=wet_elements)
        assert alltrue(wet_elevations<initial_runup_height)
        assert allclose(wet_elevations, z)        


        #print domain.get_quantity('elevation').get_maximum_value(indices=wet_elements)
        #print domain.get_quantity('elevation').get_maximum_location(indices=wet_elements)
        #print domain.get_quantity('elevation').get_maximum_index(indices=wet_elements)

        
        #--------------------------------------------------------------
        # Let triangles adjust
        #--------------------------------------------------------------
        for t in domain.evolve(yieldstep = 0.1, finaltime = 1.0):
            pass


        #--------------------------------------------------------------
        # Test inundation height again
        #--------------------------------------------------------------

        indices = domain.get_wet_elements()
        z = domain.get_quantity('elevation').\
            get_values(location='centroids', indices=indices)

        assert alltrue(z<initial_runup_height)

        q = domain.get_maximum_inundation_elevation()
        assert allclose(q, initial_runup_height, rtol = 1.0/N) # First order accuracy
        
        x, y = domain.get_maximum_inundation_location()
        qref = domain.get_quantity('elevation').get_values(interpolation_points = [[x, y]])
        assert allclose(q, qref)        


        #--------------------------------------------------------------
        # Update boundary to allow inflow
        #--------------------------------------------------------------
        domain.set_boundary({'right': Bd})

        
        #--------------------------------------------------------------
        # Evolve system through time
        #--------------------------------------------------------------
        for t in domain.evolve(yieldstep = 0.1, finaltime = 3.0):
            #print domain.timestepping_statistics(track_speeds=True)
            #domain.write_time()
            pass
    
        #--------------------------------------------------------------
        # Test inundation height again
        #--------------------------------------------------------------

        indices = domain.get_wet_elements()
        z = domain.get_quantity('elevation').\
            get_values(location='centroids', indices=indices)

        assert alltrue(z<final_runup_height)

        q = domain.get_maximum_inundation_elevation()
        assert allclose(q, final_runup_height, rtol = 1.0/N) # First order accuracy 

        x, y = domain.get_maximum_inundation_location()
        qref = domain.get_quantity('elevation').get_values(interpolation_points = [[x, y]])
        assert allclose(q, qref)        


        wet_elements = domain.get_wet_elements()
        wet_elevations = domain.get_quantity('elevation').get_values(location='centroids',
                                                                     indices=wet_elements)
        assert alltrue(wet_elevations<final_runup_height)
        assert allclose(wet_elevations, z)        
        


    def test_get_maximum_inundation_from_sww(self):
        """test_get_maximum_inundation_from_sww(self)

        Test of get_maximum_inundation_elevation()
        and get_maximum_inundation_location() from data_manager.py
        
        This is based on test_get_maximum_inundation_3(self) but works with the
        stored results instead of with the internal data structure.

        This test uses the underlying get_maximum_inundation_data for tests
        """

        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
        from data_manager import get_maximum_inundation_elevation
        from data_manager import get_maximum_inundation_location
        from data_manager import get_maximum_inundation_data
        

        initial_runup_height = -0.4
        final_runup_height = -0.3


        #--------------------------------------------------------------
        # Setup computational domain
        #--------------------------------------------------------------
        N = 10
        points, vertices, boundary = rectangular_cross(N, N) 
        domain = Domain(points, vertices, boundary)
        domain.set_name('runup_test')
        domain.set_maximum_allowed_speed(1.0)

        domain.tight_slope_limiters = 0 # FIXME: This works better with old limiters so far

        #--------------------------------------------------------------
        # Setup initial conditions
        #--------------------------------------------------------------
        def topography(x,y):
            return -x/2                             # linear bed slope
            

        domain.set_quantity('elevation', topography)       # Use function for elevation
        domain.set_quantity('friction', 0.)                # Zero friction 
        domain.set_quantity('stage', initial_runup_height) # Constant negative initial stage


        #--------------------------------------------------------------
        # Setup boundary conditions
        #--------------------------------------------------------------
        Br = Reflective_boundary(domain)              # Reflective wall
        Bd = Dirichlet_boundary([final_runup_height,  # Constant inflow
                                 0,
                                 0])

        # All reflective to begin with (still water) 
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})


        #--------------------------------------------------------------
        # Test initial inundation height
        #--------------------------------------------------------------

        indices = domain.get_wet_elements()
        z = domain.get_quantity('elevation').\
            get_values(location='centroids', indices=indices)
        assert alltrue(z<initial_runup_height)

        q_ref = domain.get_maximum_inundation_elevation()
        assert allclose(q_ref, initial_runup_height, rtol = 1.0/N) # First order accuracy 

        
        #--------------------------------------------------------------
        # Let triangles adjust
        #--------------------------------------------------------------
        for t in domain.evolve(yieldstep = 0.1, finaltime = 1.0):
            pass


        #--------------------------------------------------------------
        # Test inundation height again
        #--------------------------------------------------------------
       
        q_ref = domain.get_maximum_inundation_elevation()
        q = get_maximum_inundation_elevation('runup_test.sww')
        msg = 'We got %f, should have been %f' %(q, q_ref)
        assert allclose(q, q_ref, rtol=1.0/N), msg


        q = get_maximum_inundation_elevation('runup_test.sww')
        msg = 'We got %f, should have been %f' %(q, initial_runup_height)
        assert allclose(q, initial_runup_height, rtol = 1.0/N), msg 


        # Test error condition if time interval is out
        try:
            q = get_maximum_inundation_elevation('runup_test.sww',
                                                 time_interval=[2.0, 3.0])
        except ValueError:
            pass
        else:
            msg = 'should have caught wrong time interval'
            raise Exception, msg

        # Check correct time interval
        q, loc = get_maximum_inundation_data('runup_test.sww',
                                             time_interval=[0.0, 3.0])        
        msg = 'We got %f, should have been %f' %(q, initial_runup_height)
        assert allclose(q, initial_runup_height, rtol = 1.0/N), msg
        assert allclose(-loc[0]/2, q) # From topography formula         
        

        #--------------------------------------------------------------
        # Update boundary to allow inflow
        #--------------------------------------------------------------
        domain.set_boundary({'right': Bd})

        
        #--------------------------------------------------------------
        # Evolve system through time
        #--------------------------------------------------------------
        q_max = None
        for t in domain.evolve(yieldstep = 0.1, finaltime = 3.0,
                               skip_initial_step = True): 
            q = domain.get_maximum_inundation_elevation()
            if q > q_max: q_max = q

    
        #--------------------------------------------------------------
        # Test inundation height again
        #--------------------------------------------------------------

        indices = domain.get_wet_elements()
        z = domain.get_quantity('elevation').\
            get_values(location='centroids', indices=indices)

        assert alltrue(z<final_runup_height)

        q = domain.get_maximum_inundation_elevation()
        assert allclose(q, final_runup_height, rtol = 1.0/N) # First order accuracy

        q, loc = get_maximum_inundation_data('runup_test.sww', time_interval=[3.0, 3.0])
        msg = 'We got %f, should have been %f' %(q, final_runup_height)
        assert allclose(q, final_runup_height, rtol=1.0/N), msg
        #print 'loc', loc, q        
        assert allclose(-loc[0]/2, q) # From topography formula         

        q = get_maximum_inundation_elevation('runup_test.sww')
        loc = get_maximum_inundation_location('runup_test.sww')        
        msg = 'We got %f, should have been %f' %(q, q_max)
        assert allclose(q, q_max, rtol=1.0/N), msg
        #print 'loc', loc, q
        assert allclose(-loc[0]/2, q) # From topography formula 

        

        q = get_maximum_inundation_elevation('runup_test.sww', time_interval=[0, 3])
        msg = 'We got %f, should have been %f' %(q, q_max)
        assert allclose(q, q_max, rtol=1.0/N), msg


        # Check polygon mode
        polygon = [[0.3, 0.0], [0.9, 0.0], [0.9, 1.0], [0.3, 1.0]] # Runup region
        q = get_maximum_inundation_elevation('runup_test.sww',
                                             polygon = polygon,
                                             time_interval=[0, 3])
        msg = 'We got %f, should have been %f' %(q, q_max)
        assert allclose(q, q_max, rtol=1.0/N), msg

        
        polygon = [[0.9, 0.0], [1.0, 0.0], [1.0, 1.0], [0.9, 1.0]] # Offshore region
        q, loc = get_maximum_inundation_data('runup_test.sww',
                                             polygon = polygon,
                                             time_interval=[0, 3])
        msg = 'We got %f, should have been %f' %(q, -0.475)
        assert allclose(q, -0.475, rtol=1.0/N), msg
        assert is_inside_polygon(loc, polygon)
        assert allclose(-loc[0]/2, q) # From topography formula         


        polygon = [[0.0, 0.0], [0.4, 0.0], [0.4, 1.0], [0.0, 1.0]] # Dry region
        q, loc = get_maximum_inundation_data('runup_test.sww',
                                             polygon = polygon,
                                             time_interval=[0, 3])
        msg = 'We got %s, should have been None' %(q)
        assert q is None, msg
        msg = 'We got %s, should have been None' %(loc)
        assert loc is None, msg        

        # Check what happens if no time point is within interval
        try:
            q = get_maximum_inundation_elevation('runup_test.sww', time_interval=[2.75, 2.75])
        except AssertionError:
            pass
        else:
            msg = 'Time interval should have raised an exception'
            raise msg

        # Cleanup
        try:
            os.remove(domain.get_name() + '.' + domain.format)
        except:
            pass
            #FIXME(Ole): Windows won't allow removal of this
            
        

    def test_another_runup_example(self):
        """test_another_runup_example

        Test runup example where actual timeseries at interpolated
        points are tested.
        """

        #-----------------------------------------------------------------
        # Import necessary modules
        #-----------------------------------------------------------------

        from anuga.pmesh.mesh_interface import create_mesh_from_regions
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
        from anuga.shallow_water import Domain
        from anuga.shallow_water import Reflective_boundary
        from anuga.shallow_water import Dirichlet_boundary


        #-----------------------------------------------------------------
        # Setup computational domain
        #-----------------------------------------------------------------
        points, vertices, boundary = rectangular_cross(10, 10) # Basic mesh
        domain = Domain(points, vertices, boundary) # Create domain
        domain.set_quantities_to_be_stored(None)
        domain.set_maximum_allowed_speed(100) #
        
        # FIXME (Ole): Need tests where this is commented out
        domain.tight_slope_limiters = 0 # Backwards compatibility (14/4/7)                 
        domain.H0 = 0 # Backwards compatibility (6/2/7)
        domain.beta_h = 0.2 # Backwards compatibility (14/2/7)
        domain.use_centroid_velocities = 0 # Backwards compatibility (7/5/8)
        

        #-----------------------------------------------------------------
        # Setup initial conditions
        #-----------------------------------------------------------------

        def topography(x,y): 
            return -x/2                              # linear bed slope

        domain.set_quantity('elevation', topography) 
        domain.set_quantity('friction', 0.0)         
        domain.set_quantity('stage', expression='elevation')            


        #----------------------------------------------------------------
        # Setup boundary conditions
        #----------------------------------------------------------------

        Br = Reflective_boundary(domain)      # Solid reflective wall
        Bd = Dirichlet_boundary([-0.2,0.,0.]) # Constant boundary values
        domain.set_boundary({'left': Br, 'right': Bd, 'top': Br, 'bottom': Br})


        #----------------------------------------------------------------
        # Evolve system through time
        #----------------------------------------------------------------

        interpolation_points = [[0.4,0.5], [0.6,0.5], [0.8,0.5], [0.9,0.5]]
        gauge_values = []
        for _ in interpolation_points:
            gauge_values.append([])

        time = []
        for t in domain.evolve(yieldstep = 0.1, finaltime = 5.0):
            # Record time series at known points
            time.append(domain.get_time())
            
            stage = domain.get_quantity('stage')
            w = stage.get_values(interpolation_points=interpolation_points)
            
            for i, _ in enumerate(interpolation_points):
                gauge_values[i].append(w[i])


        #print
        #print time
        #print
        #for i, (x,y) in enumerate(interpolation_points):
        #    print i, gauge_values[i]
        #    print 

        #Reference (nautilus 13/10/2006)

        G0 = [-0.20000000000000001, -0.19999681443389281, -0.1986192343695303, -0.19147413648863046, -0.19132688908678019, -0.17642317476621105, -0.167376262630034, -0.16192452887426961, -0.15609171725778803, -0.15127107084302249, -0.14048864340360018, -0.19296484125327093, -0.19997006390580363, -0.19999999999937063, -0.19999999999937063, -0.19999999999938772, -0.19999999999938772, -0.19999999999938772, -0.19999999999938772, -0.19974288463035494, -0.19951636867991712, -0.19966301435195755, -0.19981082259800226, -0.19978575003960128, -0.19992942471933109, -0.19999999931029933, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989]


        G1 = [-0.29999999999999993, -0.29988962537199199, -0.29293904425532025, -0.28329367722887888, -0.25999146407696289, -0.22613875068011896, -0.21190705052094994, -0.19900707995208217, -0.18876305176191882, -0.18132447501091936, -0.17395459512711151, -0.15562414200985644, -0.16212999953643359, -0.18964422820514618, -0.20871181844346975, -0.21672207791083464, -0.21774940291862779, -0.21482868050219833, -0.21057786776704043, -0.20649663432591045, -0.20294932949211578, -0.19974459897911329, -0.19733648772704043, -0.19641404599824669, -0.19654095699184146, -0.19709942852191994, -0.19780873983410741, -0.19853259125123518, -0.19916495938961168, -0.19965391267799168, -0.19993539587158982, -0.2001383705551133, -0.20029344332295113, -0.20035349748150011, -0.20029886541561631, -0.20015541958920294, -0.19997273066429103, -0.19979879448668514, -0.19966016997024041, -0.19957558009501869, -0.19955725674938532, -0.19958083002853366, -0.19961752462568647, -0.19965296611330258, -0.19968998132634594, -0.19972532942208607, -0.19975372922008239, -0.19977196116929855, -0.19977951443660594, -0.19977792107284789, -0.19976991595502003]

        G2 = [-0.40000000000000002, -0.39011996186687281, -0.33359026016903887, -0.29757449757405952, -0.27594124995715791, -0.25970211955309436, -0.24482929492054245, -0.23156757139219822, -0.21956485769139392, -0.20844522129026694, -0.19856327660654355, -0.18962303467030903, -0.17371085465024955, -0.16429840256208336, -0.17793711732368575, -0.19287799702389993, -0.20236271260796762, -0.20700727993623128, -0.20847704371373174, -0.20796895600687262, -0.20653398626186478, -0.20480656169870676, -0.20295863990994492, -0.20100199602968896, -0.19940642689498472, -0.19858371478015749, -0.19838672154605322, -0.19851093923669558, -0.19878191998909323, -0.19910827645394291, -0.19943514333832094, -0.19971231361970535, -0.19992429278849655, -0.20010744405928019, -0.20025927002359642, -0.20034751667523681, -0.20035504591467249, -0.20029401385620157, -0.20019492358237226, -0.20008934249434918, -0.19999808924091636, -0.19993869218976712, -0.19991589568150098, -0.19991815777945968, -0.19993012995477188, -0.19994576118144997, -0.19996497193815974, -0.19998586151236197, -0.20000487253824847, -0.20001903000364174, -0.20002698661385457]

        G3 = [-0.45000000000000001, -0.37713945714588398, -0.33029565026933816, -0.30598209033945367, -0.28847101155177313, -0.27211191064563195, -0.25701544058818926, -0.24298945948410997, -0.23010402733784807, -0.21820351802867713, -0.20709938367218383, -0.19719881806182216, -0.18568281604361933, -0.16828653906676322, -0.16977310768235579, -0.1832707289594605, -0.19483524345250974, -0.20233480051649216, -0.20630757214159207, -0.20763927857964531, -0.20724458160595791, -0.20599191745446047, -0.20438329669495012, -0.20256105512496606, -0.20071269486729407, -0.19934403619901719, -0.19866860191898347, -0.19849975056296071, -0.19860870923007437, -0.19885838217851401, -0.19916422433758982, -0.19946861981642039, -0.19972267778871666, -0.19993013816258154, -0.20011063428833351, -0.20024891930311628, -0.20031882555219671, -0.20031326268593497, -0.20024881068472311, -0.20015443214902759, -0.20005669097631221, -0.19997542564643309, -0.19992564006223304, -0.19990746148869892, -0.19990923999172872, -0.19991956416813192, -0.19993484556273733, -0.1999538628054662, -0.19997381636620407, -0.19999130900268777, -0.20000388227457688]
        
        #FIXME (DSG):This is a hack so the anuga install, not precompiled
        # works on DSG's win2000, python 2.3
        #The problem is the gauge_values[X] are 52 long, not 51.
        #
        # This was probably fixed by Stephen in changeset:3804
        if len(gauge_values[0]) == 52: gauge_values[0].pop()
        if len(gauge_values[1]) == 52: gauge_values[1].pop()
        if len(gauge_values[2]) == 52: gauge_values[2].pop()
        if len(gauge_values[3]) == 52: gauge_values[3].pop()

##         print len(G0), len(gauge_values[0])
##         print len(G1), len(gauge_values[1])
        
        #print array(gauge_values[0])-array(G0)

        
        
        assert allclose(gauge_values[0], G0)
        assert allclose(gauge_values[1], G1)
        assert allclose(gauge_values[2], G2)
        assert allclose(gauge_values[3], G3)        







    #####################################################

    def test_flux_optimisation(self):
        """test_flux_optimisation
        Test that fluxes are correctly computed using
        dry and still cell exclusions
        """

        from anuga.config import g
        import copy

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices)

        #Set up for a gradient of (3,0) at mid triangle (bce)
        def slope(x, y):
            return 3*x

        h = 0.1
        def stage(x,y):
            return slope(x,y)+h

        domain.set_quantity('elevation', slope)
        domain.set_quantity('stage', stage)

        # Allow slope limiters to work (FIXME (Ole): Shouldn't this be automatic in ANUGA?)     
        domain.distribute_to_vertices_and_edges()       

        initial_stage = copy.copy(domain.quantities['stage'].vertex_values)

        domain.set_boundary({'exterior': Reflective_boundary(domain)})


        #  Check that update arrays are initialised to zero 
        assert allclose(domain.get_quantity('stage').explicit_update, 0)
        assert allclose(domain.get_quantity('xmomentum').explicit_update, 0)
        assert allclose(domain.get_quantity('ymomentum').explicit_update, 0)               


        # Get true values
        domain.optimise_dry_cells = False
        domain.compute_fluxes()
        stage_ref = copy.copy(domain.get_quantity('stage').explicit_update)
        xmom_ref = copy.copy(domain.get_quantity('xmomentum').explicit_update)
        ymom_ref = copy.copy(domain.get_quantity('ymomentum').explicit_update)       

        # Try with flux optimisation
        domain.optimise_dry_cells = True
        domain.compute_fluxes()

        assert allclose(stage_ref, domain.get_quantity('stage').explicit_update)
        assert allclose(xmom_ref, domain.get_quantity('xmomentum').explicit_update)
        assert allclose(ymom_ref, domain.get_quantity('ymomentum').explicit_update)
        
   
        
    def test_initial_condition(self):
        """test_initial_condition
        Test that initial condition is output at time == 0 and that
        computed values change as system evolves
        """

        from anuga.config import g
        import copy

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices)

        #Set up for a gradient of (3,0) at mid triangle (bce)
        def slope(x, y):
            return 3*x

        h = 0.1
        def stage(x,y):
            return slope(x,y)+h

        domain.set_quantity('elevation', slope)
        domain.set_quantity('stage', stage)

        # Allow slope limiters to work (FIXME (Ole): Shouldn't this be automatic in ANUGA?)     
        domain.distribute_to_vertices_and_edges()       

        initial_stage = copy.copy(domain.quantities['stage'].vertex_values)

        domain.set_boundary({'exterior': Reflective_boundary(domain)})

        domain.optimise_dry_cells = True
        #Evolution
        for t in domain.evolve(yieldstep = 0.5, finaltime = 2.0):
            stage = domain.quantities['stage'].vertex_values

            if t == 0.0:
                assert allclose(stage, initial_stage)
            else:
                assert not allclose(stage, initial_stage)


        os.remove(domain.get_name() + '.sww')



    #####################################################
    def test_gravity(self):
        #Assuming no friction

        from anuga.config import g

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices)

        #Set up for a gradient of (3,0) at mid triangle (bce)
        def slope(x, y):
            return 3*x

        h = 0.1
        def stage(x,y):
            return slope(x,y)+h

        domain.set_quantity('elevation', slope)
        domain.set_quantity('stage', stage)

        for name in domain.conserved_quantities:
            assert allclose(domain.quantities[name].explicit_update, 0)
            assert allclose(domain.quantities[name].semi_implicit_update, 0)

        domain.compute_forcing_terms()

        assert allclose(domain.quantities['stage'].explicit_update, 0)
        assert allclose(domain.quantities['xmomentum'].explicit_update, -g*h*3)
        assert allclose(domain.quantities['ymomentum'].explicit_update, 0)


    def test_manning_friction(self):
        from anuga.config import g

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices)

        #Set up for a gradient of (3,0) at mid triangle (bce)
        def slope(x, y):
            return 3*x

        h = 0.1
        def stage(x,y):
            return slope(x,y)+h

        eta = 0.07
        domain.set_quantity('elevation', slope)
        domain.set_quantity('stage', stage)
        domain.set_quantity('friction', eta)

        for name in domain.conserved_quantities:
            assert allclose(domain.quantities[name].explicit_update, 0)
            assert allclose(domain.quantities[name].semi_implicit_update, 0)

        domain.compute_forcing_terms()

        assert allclose(domain.quantities['stage'].explicit_update, 0)
        assert allclose(domain.quantities['xmomentum'].explicit_update, -g*h*3)
        assert allclose(domain.quantities['ymomentum'].explicit_update, 0)

        assert allclose(domain.quantities['stage'].semi_implicit_update, 0)
        assert allclose(domain.quantities['xmomentum'].semi_implicit_update, 0)
        assert allclose(domain.quantities['ymomentum'].semi_implicit_update, 0)

        #Create some momentum for friction to work with
        domain.set_quantity('xmomentum', 1)
        S = -g * eta**2 / h**(7.0/3)

        domain.compute_forcing_terms()
        assert allclose(domain.quantities['stage'].semi_implicit_update, 0)
        assert allclose(domain.quantities['xmomentum'].semi_implicit_update, S)
        assert allclose(domain.quantities['ymomentum'].semi_implicit_update, 0)

        #A more complex example
        domain.quantities['stage'].semi_implicit_update[:] = 0.0
        domain.quantities['xmomentum'].semi_implicit_update[:] = 0.0
        domain.quantities['ymomentum'].semi_implicit_update[:] = 0.0

        domain.set_quantity('xmomentum', 3)
        domain.set_quantity('ymomentum', 4)

        S = -g * eta**2 * 5 / h**(7.0/3)


        domain.compute_forcing_terms()

        assert allclose(domain.quantities['stage'].semi_implicit_update, 0)
        assert allclose(domain.quantities['xmomentum'].semi_implicit_update, 3*S)
        assert allclose(domain.quantities['ymomentum'].semi_implicit_update, 4*S)

    def test_constant_wind_stress(self):
        from anuga.config import rho_a, rho_w, eta_w
        from math import pi, cos, sin

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]


        domain = Domain(points, vertices)

        #Flat surface with 1m of water
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        #Setup only one forcing term, constant wind stress
        s = 100
        phi = 135
        domain.forcing_terms = []
        domain.forcing_terms.append( Wind_stress(s, phi) )

        domain.compute_forcing_terms()


        const = eta_w*rho_a/rho_w

        #Convert to radians
        phi = phi*pi/180

        #Compute velocity vector (u, v)
        u = s*cos(phi)
        v = s*sin(phi)

        #Compute wind stress
        S = const * sqrt(u**2 + v**2)

        assert allclose(domain.quantities['stage'].explicit_update, 0)
        assert allclose(domain.quantities['xmomentum'].explicit_update, S*u)
        assert allclose(domain.quantities['ymomentum'].explicit_update, S*v)


    def test_variable_wind_stress(self):
        from anuga.config import rho_a, rho_w, eta_w
        from math import pi, cos, sin

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices)

        #Flat surface with 1m of water
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})


        domain.time = 5.54 #Take a random time (not zero)

        #Setup only one forcing term, constant wind stress
        s = 100
        phi = 135
        domain.forcing_terms = []
        domain.forcing_terms.append( Wind_stress(s = speed, phi = angle) )

        domain.compute_forcing_terms()

        #Compute reference solution
        const = eta_w*rho_a/rho_w

        N = len(domain) # number_of_triangles

        xc = domain.get_centroid_coordinates()
        t = domain.time

        x = xc[:,0]
        y = xc[:,1]
        s_vec = speed(t,x,y)
        phi_vec = angle(t,x,y)


        for k in range(N):
            #Convert to radians
            phi = phi_vec[k]*pi/180
            s = s_vec[k]

            #Compute velocity vector (u, v)
            u = s*cos(phi)
            v = s*sin(phi)

            #Compute wind stress
            S = const * sqrt(u**2 + v**2)

            assert allclose(domain.quantities['stage'].explicit_update[k], 0)
            assert allclose(domain.quantities['xmomentum'].explicit_update[k], S*u)
            assert allclose(domain.quantities['ymomentum'].explicit_update[k], S*v)






    def test_windfield_from_file(self):
        from anuga.config import rho_a, rho_w, eta_w
        from math import pi, cos, sin
        from anuga.config import time_format
        from anuga.abstract_2d_finite_volumes.util import file_function
        import time


        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices)

        #Flat surface with 1m of water
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})


        domain.time = 7 #Take a time that is represented in file (not zero)

        #Write wind stress file (ensure that domain.time is covered)
        #Take x=1 and y=0
        filename = 'test_windstress_from_file'
        start = time.mktime(time.strptime('2000', '%Y'))
        fid = open(filename + '.txt', 'w')
        dt = 1  #One second interval
        t = 0.0
        while t <= 10.0:
            t_string = time.strftime(time_format, time.gmtime(t+start))

            fid.write('%s, %f %f\n' %(t_string,
                                      speed(t,[1],[0])[0],
                                      angle(t,[1],[0])[0]))
            t += dt

        fid.close()


        #Convert ASCII file to NetCDF (Which is what we really like!)
        from data_manager import timefile2netcdf        
        timefile2netcdf(filename)
        os.remove(filename + '.txt')

        
        #Setup wind stress
        F = file_function(filename + '.tms', quantities = ['Attribute0',
                                                           'Attribute1'])
        os.remove(filename + '.tms')
        

        #print 'F(5)', F(5)

        #print 'F(5,x,y)', F(5,x=zeros(3),y=zeros(3))
        
        #print dir(F)
        #print F.T
        #print F.precomputed_values
        #
        #F = file_function(filename + '.txt')        
        #
        #print dir(F)
        #print F.T        
        #print F.Q
        
        W = Wind_stress(F)

        domain.forcing_terms = []
        domain.forcing_terms.append(W)

        domain.compute_forcing_terms()

        #Compute reference solution
        const = eta_w*rho_a/rho_w

        N = len(domain) # number_of_triangles

        t = domain.time

        s = speed(t,[1],[0])[0]
        phi = angle(t,[1],[0])[0]

        #Convert to radians
        phi = phi*pi/180


        #Compute velocity vector (u, v)
        u = s*cos(phi)
        v = s*sin(phi)

        #Compute wind stress
        S = const * sqrt(u**2 + v**2)

        for k in range(N):
            assert allclose(domain.quantities['stage'].explicit_update[k], 0)
            assert allclose(domain.quantities['xmomentum'].explicit_update[k], S*u)
            assert allclose(domain.quantities['ymomentum'].explicit_update[k], S*v)


    def test_windfield_from_file_seconds(self):
        from anuga.config import rho_a, rho_w, eta_w
        from math import pi, cos, sin
        from anuga.config import time_format
        from anuga.abstract_2d_finite_volumes.util import file_function
        import time


        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices)

        #Flat surface with 1m of water
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})


        domain.time = 7 #Take a time that is represented in file (not zero)

        #Write wind stress file (ensure that domain.time is covered)
        #Take x=1 and y=0
        filename = 'test_windstress_from_file'
        start = time.mktime(time.strptime('2000', '%Y'))
        fid = open(filename + '.txt', 'w')
        dt = 0.5 #1  #One second interval
        t = 0.0
        while t <= 10.0:
            fid.write('%s, %f %f\n' %(str(t),
                                      speed(t,[1],[0])[0],
                                      angle(t,[1],[0])[0]))
            t += dt

        fid.close()


        #Convert ASCII file to NetCDF (Which is what we really like!)
        from data_manager import timefile2netcdf        
        timefile2netcdf(filename, time_as_seconds=True)
        os.remove(filename + '.txt')

        
        #Setup wind stress
        F = file_function(filename + '.tms', quantities = ['Attribute0',
                                                           'Attribute1'])
        os.remove(filename + '.tms')
        

        #print 'F(5)', F(5)

        #print 'F(5,x,y)', F(5,x=zeros(3),y=zeros(3))
        
        #print dir(F)
        #print F.T
        #print F.precomputed_values
        #
        #F = file_function(filename + '.txt')        
        #
        #print dir(F)
        #print F.T        
        #print F.Q
        
        W = Wind_stress(F)

        domain.forcing_terms = []
        domain.forcing_terms.append(W)

        domain.compute_forcing_terms()

        #Compute reference solution
        const = eta_w*rho_a/rho_w

        N = len(domain) # number_of_triangles

        t = domain.time

        s = speed(t,[1],[0])[0]
        phi = angle(t,[1],[0])[0]

        #Convert to radians
        phi = phi*pi/180


        #Compute velocity vector (u, v)
        u = s*cos(phi)
        v = s*sin(phi)

        #Compute wind stress
        S = const * sqrt(u**2 + v**2)

        for k in range(N):
            assert allclose(domain.quantities['stage'].explicit_update[k], 0)
            assert allclose(domain.quantities['xmomentum'].explicit_update[k], S*u)
            assert allclose(domain.quantities['ymomentum'].explicit_update[k], S*v)


        

    def test_wind_stress_error_condition(self):
        """Test that windstress reacts properly when forcing functions
        are wrong - e.g. returns a scalar
        """

        from anuga.config import rho_a, rho_w, eta_w
        from math import pi, cos, sin

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices)

        #Flat surface with 1m of water
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})


        domain.time = 5.54 #Take a random time (not zero)

        #Setup only one forcing term, bad func
        domain.forcing_terms = []

        try:
            domain.forcing_terms.append(Wind_stress(s = scalar_func,
                                                    phi = angle))
        except AssertionError:
            pass
        else:
            msg = 'Should have raised exception'
            raise msg


        try:
            domain.forcing_terms.append(Wind_stress(s = speed,
                                                    phi = scalar_func))
        except AssertionError:
            pass
        else:
            msg = 'Should have raised exception'
            raise msg

        try:
            domain.forcing_terms.append(Wind_stress(s = speed,
                                                    phi = 'xx'))
        except:
            pass
        else:
            msg = 'Should have raised exception'
            raise msg



    def test_rainfall(self):
        from math import pi, cos, sin

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]


        domain = Domain(points, vertices)

        #Flat surface with 1m of water
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        # Setup only one forcing term, constant rainfall
        domain.forcing_terms = []
        domain.forcing_terms.append( Rainfall(domain, rate=2.0) )

        domain.compute_forcing_terms()
        assert allclose(domain.quantities['stage'].explicit_update, 2.0/1000)


        # FIXME: Do Time dependent rainfall



    def test_rainfall_restricted_by_polygon(self):
        from math import pi, cos, sin

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]


        domain = Domain(points, vertices)

        #Flat surface with 1m of water
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        # Setup only one forcing term, constant rainfall restricted to a polygon enclosing triangle #1 (bce)
        domain.forcing_terms = []
        R = Rainfall(domain, rate=2.0, polygon = [[1,1], [2,1], [2,2], [1,2]])

        assert allclose(R.exchange_area, 1)
        
        domain.forcing_terms.append(R)

        domain.compute_forcing_terms()
        #print domain.quantities['stage'].explicit_update
        
        assert allclose(domain.quantities['stage'].explicit_update[1], 2.0/1000)
        assert allclose(domain.quantities['stage'].explicit_update[0], 0)
        assert allclose(domain.quantities['stage'].explicit_update[2:], 0)        
        
        # FIXME: Do Time dependent rainfall with poly




    def test_inflow_using_circle(self):
        from math import pi, cos, sin

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]


        domain = Domain(points, vertices)

        # Flat surface with 1m of water
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        # Setup only one forcing term, constant inflow of 2 m^3/s on a circle affecting triangles #0 and #1 (bac and bce)
        domain.forcing_terms = []
        domain.forcing_terms.append( Inflow(domain, rate=2.0, center=(1,1), radius=1) )

        domain.compute_forcing_terms()
        #print domain.quantities['stage'].explicit_update
        
        assert allclose(domain.quantities['stage'].explicit_update[1], 2.0/pi)
        assert allclose(domain.quantities['stage'].explicit_update[0], 2.0/pi)
        assert allclose(domain.quantities['stage'].explicit_update[2:], 0)        


    def test_inflow_using_circle_function(self):
        from math import pi, cos, sin

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]


        domain = Domain(points, vertices)

        # Flat surface with 1m of water
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        # Setup only one forcing term, time dependent inflow of 2 m^3/s on a circle affecting triangles #0 and #1 (bac and bce)
        domain.forcing_terms = []
        domain.forcing_terms.append( Inflow(domain, rate=lambda t: 2., center=(1,1), radius=1) )

        domain.compute_forcing_terms()
        
        assert allclose(domain.quantities['stage'].explicit_update[1], 2.0/pi)
        assert allclose(domain.quantities['stage'].explicit_update[0], 2.0/pi)
        assert allclose(domain.quantities['stage'].explicit_update[2:], 0)        
        


    #####################################################
    def test_first_order_extrapolator_const_z(self):

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices)
        val0 = 2.+2.0/3
        val1 = 4.+4.0/3
        val2 = 8.+2.0/3
        val3 = 2.+8.0/3

        zl=zr=-3.75 #Assume constant bed (must be less than stage)
        domain.set_quantity('elevation', zl*ones( (4,3) ))
        domain.set_quantity('stage', [[val0, val0-1, val0-2],
                                      [val1, val1+1, val1],
                                      [val2, val2-2, val2],
                                      [val3-0.5, val3, val3]])



        domain._order_ = 1
        domain.distribute_to_vertices_and_edges()

        #Check that centroid values were distributed to vertices
        C = domain.quantities['stage'].centroid_values
        for i in range(3):
            assert allclose( domain.quantities['stage'].vertex_values[:,i], C)


    def test_first_order_limiter_variable_z(self):
        #Check that first order limiter follows bed_slope
        from Numeric import alltrue, greater_equal
        from anuga.config import epsilon

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices)
        val0 = 2.+2.0/3
        val1 = 4.+4.0/3
        val2 = 8.+2.0/3
        val3 = 2.+8.0/3

        domain.set_quantity('elevation', [[0,0,0], [6,0,0],
                                          [6,6,6], [6,6,6]])
        domain.set_quantity('stage', [[val0, val0, val0],
                                      [val1, val1, val1],
                                      [val2, val2, val2],
                                      [val3, val3, val3]])

        E = domain.quantities['elevation'].vertex_values
        L = domain.quantities['stage'].vertex_values


        #Check that some stages are not above elevation (within eps)
        #- so that the limiter has something to work with
        assert not alltrue(alltrue(greater_equal(L,E-epsilon)))

        domain._order_ = 1
        domain.distribute_to_vertices_and_edges()

        #Check that all stages are above elevation (within eps)
        assert alltrue(alltrue(greater_equal(L,E-epsilon)))


    #####################################################
    def test_distribute_basic(self):
        #Using test data generated by abstract_2d_finite_volumes-2
        #Assuming no friction and flat bed (0.0)

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices)

        val0 = 2.
        val1 = 4.
        val2 = 8.
        val3 = 2.

        domain.set_quantity('stage', [val0, val1, val2, val3],
                            location='centroids')
        L = domain.quantities['stage'].vertex_values

        #First order
        domain._order_ = 1
        domain.distribute_to_vertices_and_edges()
        assert allclose(L[1], val1)

        #Second order
        domain._order_ = 2
        domain.beta_w      = 0.9
        domain.beta_w_dry  = 0.9
        domain.beta_uh     = 0.9
        domain.beta_uh_dry = 0.9
        domain.beta_vh     = 0.9
        domain.beta_vh_dry = 0.9
        domain.distribute_to_vertices_and_edges()
        assert allclose(L[1], [2.2, 4.9, 4.9])



    def test_distribute_away_from_bed(self):
        #Using test data generated by abstract_2d_finite_volumes-2
        #Assuming no friction and flat bed (0.0)

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices)
        L = domain.quantities['stage'].vertex_values

        def stage(x,y):
            return x**2

        domain.set_quantity('stage', stage, location='centroids')

        domain.quantities['stage'].compute_gradients()

        a, b = domain.quantities['stage'].get_gradients()
                
        assert allclose(a[1], 3.33333334)
        assert allclose(b[1], 0.0)

        domain._order_ = 1
        domain.distribute_to_vertices_and_edges()
        assert allclose(L[1], 1.77777778)

        domain._order_ = 2
        domain.beta_w      = 0.9
        domain.beta_w_dry  = 0.9
        domain.beta_uh     = 0.9
        domain.beta_uh_dry = 0.9
        domain.beta_vh     = 0.9
        domain.beta_vh_dry = 0.9
        domain.distribute_to_vertices_and_edges()
        assert allclose(L[1], [0.57777777, 2.37777778, 2.37777778])



    def test_distribute_away_from_bed1(self):
        #Using test data generated by abstract_2d_finite_volumes-2
        #Assuming no friction and flat bed (0.0)

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices)
        L = domain.quantities['stage'].vertex_values

        def stage(x,y):
            return x**4+y**2

        domain.set_quantity('stage', stage, location='centroids')
        #print domain.quantities['stage'].centroid_values

        domain.quantities['stage'].compute_gradients()
        a, b = domain.quantities['stage'].get_gradients()
        assert allclose(a[1], 25.18518519)
        assert allclose(b[1], 3.33333333)

        domain._order_ = 1
        domain.distribute_to_vertices_and_edges()
        assert allclose(L[1], 4.9382716)

        domain._order_ = 2
        domain.beta_w      = 0.9
        domain.beta_w_dry  = 0.9
        domain.beta_uh     = 0.9
        domain.beta_uh_dry = 0.9
        domain.beta_vh     = 0.9
        domain.beta_vh_dry = 0.9
        domain.distribute_to_vertices_and_edges()
        assert allclose(L[1], [1.07160494, 6.46058131, 7.28262855])



    def test_distribute_near_bed(self):

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices)


        #Set up for a gradient of (10,0) at mid triangle (bce)
        def slope(x, y):
            return 10*x

        h = 0.1
        def stage(x, y):
            return slope(x, y) + h

        domain.set_quantity('elevation', slope)
        domain.set_quantity('stage', stage, location='centroids')

        #print domain.quantities['elevation'].centroid_values
        #print domain.quantities['stage'].centroid_values

        E = domain.quantities['elevation'].vertex_values
        L = domain.quantities['stage'].vertex_values

        # Get reference values
        volumes = []
        for i in range(len(L)):
            volumes.append(sum(L[i])/3)
            assert allclose(volumes[i], domain.quantities['stage'].centroid_values[i])  
        
        
        domain._order_ = 1
        
        domain.tight_slope_limiters = 0
        domain.distribute_to_vertices_and_edges()
        assert allclose(L[1], [0.1, 20.1, 20.1])
        for i in range(len(L)):
            assert allclose(volumes[i], sum(L[i])/3)                    
        
        domain.tight_slope_limiters = 1 # Allow triangle to be flatter (closer to bed)
        domain.distribute_to_vertices_and_edges()
        assert allclose(L[1], [0.298, 20.001, 20.001])
        for i in range(len(L)):
            assert allclose(volumes[i], sum(L[i])/3)    

        domain._order_ = 2
        
        domain.tight_slope_limiters = 0
        domain.distribute_to_vertices_and_edges()
        assert allclose(L[1], [0.1, 20.1, 20.1])        
        for i in range(len(L)):
            assert allclose(volumes[i], sum(L[i])/3)            
        
        domain.tight_slope_limiters = 1 # Allow triangle to be flatter (closer to bed)
        domain.distribute_to_vertices_and_edges()
        assert allclose(L[1], [0.298, 20.001, 20.001])
        for i in range(len(L)):
            assert allclose(volumes[i], sum(L[i])/3)    
        


    def test_distribute_near_bed1(self):

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices)


        #Set up for a gradient of (8,2) at mid triangle (bce)
        def slope(x, y):
            return x**4+y**2

        h = 0.1
        def stage(x,y):
            return slope(x,y)+h

        domain.set_quantity('elevation', slope)
        domain.set_quantity('stage', stage)

        #print domain.quantities['elevation'].centroid_values
        #print domain.quantities['stage'].centroid_values

        E = domain.quantities['elevation'].vertex_values
        L = domain.quantities['stage'].vertex_values

        # Get reference values
        volumes = []
        for i in range(len(L)):
            volumes.append(sum(L[i])/3)
            assert allclose(volumes[i], domain.quantities['stage'].centroid_values[i])  
        
        #print E
        domain._order_ = 1
        
        domain.tight_slope_limiters = 0
        domain.distribute_to_vertices_and_edges()
        assert allclose(L[1], [4.1, 16.1, 20.1])        
        for i in range(len(L)):
            assert allclose(volumes[i], sum(L[i])/3)
        
                
        domain.tight_slope_limiters = 1 # Allow triangle to be flatter (closer to bed)
        domain.distribute_to_vertices_and_edges()
        assert allclose(L[1], [4.2386, 16.0604, 20.001])
        for i in range(len(L)):
            assert allclose(volumes[i], sum(L[i])/3)    
        

        domain._order_ = 2
        
        domain.tight_slope_limiters = 0    
        domain.distribute_to_vertices_and_edges()
        assert allclose(L[1], [4.1, 16.1, 20.1])
        for i in range(len(L)):
            assert allclose(volumes[i], sum(L[i])/3)    
        
        domain.tight_slope_limiters = 1 # Allow triangle to be flatter (closer to bed)
        domain.distribute_to_vertices_and_edges()
        #print L[1]
        assert allclose(L[1], [4.23370103, 16.06529897, 20.001]) or\
               allclose(L[1], [4.18944138, 16.10955862, 20.001]) or\
               allclose(L[1], [4.19351461, 16.10548539, 20.001]) # old limiters
        
        for i in range(len(L)):
            assert allclose(volumes[i], sum(L[i])/3)


    def test_second_order_distribute_real_data(self):
        #Using test data generated by abstract_2d_finite_volumes-2
        #Assuming no friction and flat bed (0.0)

        a = [0.0, 0.0]
        b = [0.0, 1.0/5]
        c = [0.0, 2.0/5]
        d = [1.0/5, 0.0]
        e = [1.0/5, 1.0/5]
        f = [1.0/5, 2.0/5]
        g = [2.0/5, 2.0/5]

        points = [a, b, c, d, e, f, g]
        #bae, efb, cbf, feg
        vertices = [ [1,0,4], [4,5,1], [2,1,5], [5,4,6]]

        domain = Domain(points, vertices)

        def slope(x, y):
            return -x/3

        domain.set_quantity('elevation', slope)
        domain.set_quantity('stage',
                            [0.01298164, 0.00365611,
                             0.01440365, -0.0381856437096],
                            location='centroids')
        domain.set_quantity('xmomentum',
                            [0.00670439, 0.01263789,
                             0.00647805, 0.0178180740668],
                            location='centroids')
        domain.set_quantity('ymomentum',
                            [-7.23510980e-004, -6.30413883e-005,
                             6.30413883e-005, 0.000200907255866],
                            location='centroids')

        E = domain.quantities['elevation'].vertex_values
        L = domain.quantities['stage'].vertex_values
        X = domain.quantities['xmomentum'].vertex_values
        Y = domain.quantities['ymomentum'].vertex_values

        #print E
        domain._order_ = 2
        domain.beta_w      = 0.9
        domain.beta_w_dry  = 0.9
        domain.beta_uh     = 0.9
        domain.beta_uh_dry = 0.9
        domain.beta_vh     = 0.9
        domain.beta_vh_dry = 0.9
        domain.beta_h = 0.0 #Use first order in h-limiter
        
        # FIXME (Ole): Need tests where this is commented out
        domain.tight_slope_limiters = 0 # Backwards compatibility (14/4/7)                 
        domain.use_centroid_velocities = 0 # Backwards compatibility (7/5/8)
        
                
        domain.distribute_to_vertices_and_edges()

        #print L[1,:]
        #print X[1,:]
        #print Y[1,:]

        assert allclose(L[1,:], [-0.00825735775384,
                                 -0.00801881482869,
                                 0.0272445025825])
        assert allclose(X[1,:], [0.0143507718962,
                                 0.0142502147066,
                                 0.00931268339717])
        assert allclose(Y[1,:], [-0.000117062180693,
                                 7.94434448109e-005,
                                 -0.000151505429018])



    def test_balance_deep_and_shallow(self):
        """Test that balanced limiters preserve conserved quantites.
        This test is using old depth based balanced limiters
        """
        import copy

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]

        #bac, bce, ecf, dbe
        elements = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4] ]

        domain = Domain(points, elements)
        domain.check_integrity()

        #Create a deliberate overshoot
        domain.set_quantity('stage', [[3,0,3], [2,2,6], [5,3,8], [8,3,5]])
        domain.set_quantity('elevation', 0) #Flat bed
        stage = domain.quantities['stage']

        ref_centroid_values = copy.copy(stage.centroid_values[:]) #Copy

        #Limit
        domain.tight_slope_limiters = 0                
        domain.distribute_to_vertices_and_edges()

        #Assert that quantities are conserved
        from Numeric import sum
        for k in range(len(domain)):
            assert allclose (ref_centroid_values[k],
                             sum(stage.vertex_values[k,:])/3)


        #Now try with a non-flat bed - closely hugging initial stage in places
        #This will create alphas in the range [0, 0.478260, 1]
        domain.set_quantity('stage', [[3,0,3], [2,2,6], [5,3,8], [8,3,5]])
        domain.set_quantity('elevation', [[0,0,0],
                                        [1.8,1.9,5.9],
                                        [4.6,0,0],
                                        [0,2,4]])
        stage = domain.quantities['stage']

        ref_centroid_values = copy.copy(stage.centroid_values[:]) #Copy
        ref_vertex_values = copy.copy(stage.vertex_values[:]) #Copy

        #Limit
        domain.tight_slope_limiters = 0        
        domain.distribute_to_vertices_and_edges()


        #Assert that all vertex quantities have changed
        for k in range(len(domain)):
            #print ref_vertex_values[k,:], stage.vertex_values[k,:]
            assert not allclose (ref_vertex_values[k,:], stage.vertex_values[k,:])
        #and assert that quantities are still conserved
        from Numeric import sum
        for k in range(len(domain)):
            assert allclose (ref_centroid_values[k],
                             sum(stage.vertex_values[k,:])/3)


        # Check actual results
        assert allclose (stage.vertex_values,
                         [[2,2,2],
                          [1.93333333, 2.03333333, 6.03333333],
                          [6.93333333, 4.53333333, 4.53333333],
                          [5.33333333, 5.33333333, 5.33333333]])


    def test_balance_deep_and_shallow_tight_SL(self):
        """Test that balanced limiters preserve conserved quantites.
        This test is using Tight Slope Limiters
        """
        import copy

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]

        #bac, bce, ecf, dbe
        elements = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4] ]

        domain = Domain(points, elements)
        domain.check_integrity()

        #Create a deliberate overshoot
        domain.set_quantity('stage', [[3,0,3], [2,2,6], [5,3,8], [8,3,5]])
        domain.set_quantity('elevation', 0) #Flat bed
        stage = domain.quantities['stage']

        ref_centroid_values = copy.copy(stage.centroid_values[:]) #Copy

        #Limit
        domain.tight_slope_limiters = 1                
        domain.distribute_to_vertices_and_edges()

        #Assert that quantities are conserved
        from Numeric import sum
        for k in range(len(domain)):
            assert allclose (ref_centroid_values[k],
                             sum(stage.vertex_values[k,:])/3)


        #Now try with a non-flat bed - closely hugging initial stage in places
        #This will create alphas in the range [0, 0.478260, 1]
        domain.set_quantity('stage', [[3,0,3], [2,2,6], [5,3,8], [8,3,5]])
        domain.set_quantity('elevation', [[0,0,0],
                                        [1.8,1.9,5.9],
                                        [4.6,0,0],
                                        [0,2,4]])
        stage = domain.quantities['stage']

        ref_centroid_values = copy.copy(stage.centroid_values[:]) #Copy
        ref_vertex_values = copy.copy(stage.vertex_values[:]) #Copy

        #Limit
        domain.tight_slope_limiters = 1        
        domain.distribute_to_vertices_and_edges()


        #Assert that all vertex quantities have changed
        for k in range(len(domain)):
            #print ref_vertex_values[k,:], stage.vertex_values[k,:]
            assert not allclose (ref_vertex_values[k,:], stage.vertex_values[k,:])
        #and assert that quantities are still conserved
        from Numeric import sum
        for k in range(len(domain)):
            assert allclose (ref_centroid_values[k],
                             sum(stage.vertex_values[k,:])/3)


        #Also check that Python and C version produce the same
        # No longer applicable if tight_slope_limiters == 1
        #print stage.vertex_values
        #assert allclose (stage.vertex_values,
        #                 [[2,2,2],
        #                  [1.93333333, 2.03333333, 6.03333333],
        #                  [6.93333333, 4.53333333, 4.53333333],
        #                  [5.33333333, 5.33333333, 5.33333333]])



    def test_balance_deep_and_shallow_Froude(self):
        """Test that balanced limiters preserve conserved quantites -
        and also that excessive Froude numbers are dealt with.
        This test is using tight slope limiters.
        """
        import copy
        from Numeric import sqrt, absolute

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]

        # bac, bce, ecf, dbe
        elements = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4] ]

        domain = Domain(points, elements)
        domain.check_integrity()
        domain.tight_slope_limiters = True
        domain.use_centroid_velocities = True                

        # Create non-flat bed - closely hugging initial stage in places
        # This will create alphas in the range [0, 0.478260, 1]
        domain.set_quantity('stage', [[3,0,3], [2,2,6], [5,3,8], [8,3,5]])
        domain.set_quantity('elevation', [[0,0,0],
                                        [1.8,1.999,5.999],
                                        [4.6,0,0],
                                        [0,2,4]])

        # Create small momenta, that nonetheless will generate large speeds
        # due to shallow depth at isolated vertices
        domain.set_quantity('xmomentum', -0.0058)
        domain.set_quantity('ymomentum', 0.0890)



        
        stage = domain.quantities['stage']
        elevation = domain.quantities['elevation']
        xmomentum = domain.quantities['xmomentum']
        ymomentum = domain.quantities['ymomentum']        

        # Setup triangle #1 to mimick real Froude explosion observed
        # in the Onslow example 13 Nov 2007.

        stage.vertex_values[1,:] = [1.6385, 1.6361, 1.2953]
        elevation.vertex_values[1,:] = [1.6375, 1.6336, 0.4647]        
        xmomentum.vertex_values[1,:] = [-0.0058, -0.0050, -0.0066]
        ymomentum.vertex_values[1,:] = [0.0890, 0.0890, 0.0890]

        xmomentum.interpolate()
        ymomentum.interpolate()        
        stage.interpolate()       
        elevation.interpolate()

        # Verify interpolation
        assert allclose(stage.centroid_values[1], 1.5233)
        assert allclose(elevation.centroid_values[1], 1.2452667)
        assert allclose(xmomentum.centroid_values[1], -0.0058)
        assert allclose(ymomentum.centroid_values[1], 0.089)

        # Derived quantities
        depth = stage-elevation
        u = xmomentum/depth
        v = ymomentum/depth

        denom = (depth*g)**0.5 
        Fx = u/denom
        Fy = v/denom
        
   
        # Verify against Onslow example (14 Nov 2007)
        assert allclose(depth.centroid_values[1], 0.278033)
        assert allclose(u.centroid_values[1], -0.0208608)
        assert allclose(v.centroid_values[1], 0.3201055)

        assert allclose(denom.centroid_values[1],
                        sqrt(depth.centroid_values[1]*g))

        assert allclose(u.centroid_values[1]/denom.centroid_values[1],
                        -0.012637746977)
        assert allclose(Fx.centroid_values[1],
                        u.centroid_values[1]/denom.centroid_values[1])

        # Check that Froude numbers are small at centroids.
        assert allclose(Fx.centroid_values[1], -0.012637746977)
        assert allclose(Fy.centroid_values[1], 0.193924048435)


        # But Froude numbers are huge at some vertices and edges
        assert allclose(Fx.vertex_values[1,:], [-5.85888475e+01,
                                                -1.27775313e+01,
                                                -2.78511420e-03])

        assert allclose(Fx.edge_values[1,:], [-6.89150773e-03,
                                              -7.38672488e-03,
                                              -2.35626238e+01])

        assert allclose(Fy.vertex_values[1,:], [8.99035764e+02,
                                                2.27440057e+02,
                                                3.75568430e-02])

        assert allclose(Fy.edge_values[1,:], [1.05748998e-01,
                                              1.06035244e-01,
                                              3.88346947e+02])
        
        
        # The task is now to arrange the limiters such that Froude numbers
        # remain under control whil at the same time obeying the conservation
        # laws.

        
        ref_centroid_values = copy.copy(stage.centroid_values[:]) #Copy
        ref_vertex_values = copy.copy(stage.vertex_values[:]) #Copy

        # Limit (and invoke balance_deep_and_shallow)
        domain.tight_slope_limiters = 1
        domain.distribute_to_vertices_and_edges()

        # Redo derived quantities
        depth = stage-elevation
        u = xmomentum/depth
        v = ymomentum/depth

        # Assert that all vertex velocities stay within one
        # order of magnitude of centroid velocities.
        #print u.vertex_values[1,:]
        #print u.centroid_values[1]
        
        assert alltrue(absolute(u.vertex_values[1,:]) <= absolute(u.centroid_values[1])*10)
        assert alltrue(absolute(v.vertex_values[1,:]) <= absolute(v.centroid_values[1])*10) 

        denom = (depth*g)**0.5 
        Fx = u/denom
        Fy = v/denom


        # Assert that Froude numbers are less than max value (TBA)
        # at vertices, edges and centroids.
        from anuga.config import maximum_froude_number
        assert alltrue(absolute(Fx.vertex_values[1,:]) < maximum_froude_number)
        assert alltrue(absolute(Fy.vertex_values[1,:]) < maximum_froude_number)


        # Assert that all vertex quantities have changed
        for k in range(len(domain)):
            #print ref_vertex_values[k,:], stage.vertex_values[k,:]
            assert not allclose (ref_vertex_values[k,:],
                                 stage.vertex_values[k,:])
            
        # Assert that quantities are still conserved
        from Numeric import sum
        for k in range(len(domain)):
            assert allclose (ref_centroid_values[k],
                             sum(stage.vertex_values[k,:])/3)


        
        return
    
        qwidth = 12
        for k in [1]: #range(len(domain)):
            print 'Triangle %d (C, V, E)' %k
            
            print 'stage'.ljust(qwidth), stage.centroid_values[k],\
                  stage.vertex_values[k,:], stage.edge_values[k,:]
            print 'elevation'.ljust(qwidth), elevation.centroid_values[k],\
                  elevation.vertex_values[k,:], elevation.edge_values[k,:]
            print 'depth'.ljust(qwidth), depth.centroid_values[k],\
                  depth.vertex_values[k,:], depth.edge_values[k,:]
            print 'xmomentum'.ljust(qwidth), xmomentum.centroid_values[k],\
                  xmomentum.vertex_values[k,:], xmomentum.edge_values[k,:]
            print 'ymomentum'.ljust(qwidth), ymomentum.centroid_values[k],\
                  ymomentum.vertex_values[k,:], ymomentum.edge_values[k,:]
            print 'u'.ljust(qwidth),u.centroid_values[k],\
                  u.vertex_values[k,:], u.edge_values[k,:]
            print 'v'.ljust(qwidth), v.centroid_values[k],\
                  v.vertex_values[k,:], v.edge_values[k,:]
            print 'Fx'.ljust(qwidth), Fx.centroid_values[k],\
                  Fx.vertex_values[k,:], Fx.edge_values[k,:]
            print 'Fy'.ljust(qwidth), Fy.centroid_values[k],\
                  Fy.vertex_values[k,:], Fy.edge_values[k,:]
            
            
        



    def test_conservation_1(self):
        """Test that stage is conserved globally

        This one uses a flat bed, reflective bdries and a suitable
        initial condition
        """
        from mesh_factory import rectangular
        from Numeric import array

        #Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order=2

        #IC
        def x_slope(x, y):
            return x/3

        domain.set_quantity('elevation', 0)
        domain.set_quantity('friction', 0)
        domain.set_quantity('stage', x_slope)

        # Boundary conditions (reflective everywhere)
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        domain.check_integrity()

        initial_volume = domain.quantities['stage'].get_integral()
        initial_xmom = domain.quantities['xmomentum'].get_integral()

        #print initial_xmom

        #Evolution
        for t in domain.evolve(yieldstep = 0.05, finaltime = 5.0):
            volume =  domain.quantities['stage'].get_integral()
            assert allclose (volume, initial_volume)

            #I don't believe that the total momentum should be the same
            #It starts with zero and ends with zero though
            #xmom = domain.quantities['xmomentum'].get_integral()
            #print xmom
            #assert allclose (xmom, initial_xmom)

        os.remove(domain.get_name() + '.sww')


    def test_conservation_2(self):
        """Test that stage is conserved globally

        This one uses a slopy bed, reflective bdries and a suitable
        initial condition
        """
        from mesh_factory import rectangular
        from Numeric import array

        #Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order=2

        #IC
        def x_slope(x, y):
            return x/3

        domain.set_quantity('elevation', x_slope)
        domain.set_quantity('friction', 0)
        domain.set_quantity('stage', 0.4) #Steady

        # Boundary conditions (reflective everywhere)
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        domain.check_integrity()

        initial_volume = domain.quantities['stage'].get_integral()
        initial_xmom = domain.quantities['xmomentum'].get_integral()

        #print initial_xmom

        #Evolution
        for t in domain.evolve(yieldstep = 0.05, finaltime = 5.0):
            volume =  domain.quantities['stage'].get_integral()
            assert allclose (volume, initial_volume)

            #FIXME: What would we expect from momentum
            #xmom = domain.quantities['xmomentum'].get_integral()
            #print xmom
            #assert allclose (xmom, initial_xmom)

        os.remove(domain.get_name() + '.sww')

    def test_conservation_3(self):
        """Test that stage is conserved globally

        This one uses a larger grid, convoluted bed, reflective bdries and a suitable
        initial condition
        """
        from mesh_factory import rectangular
        from Numeric import array

        #Create basic mesh
        points, vertices, boundary = rectangular(2, 1)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order = 2
        domain.beta_h = 0.2
        domain.set_quantities_to_be_stored(['stage', 'xmomentum', 'ymomentum'])

        #IC
        def x_slope(x, y):
            z = 0*x
            for i in range(len(x)):
                if x[i] < 0.3:
                    z[i] = x[i]/3
                if 0.3 <= x[i] < 0.5:
                    z[i] = -0.5
                if 0.5 <= x[i] < 0.7:
                    z[i] = 0.39
                if 0.7 <= x[i]:
                    z[i] = x[i]/3
            return z



        domain.set_quantity('elevation', x_slope)
        domain.set_quantity('friction', 0)
        domain.set_quantity('stage', 0.4) #Steady

        # Boundary conditions (reflective everywhere)
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        domain.check_integrity()

        initial_volume = domain.quantities['stage'].get_integral()
        initial_xmom = domain.quantities['xmomentum'].get_integral()

        import copy
        ref_centroid_values =\
                 copy.copy(domain.quantities['stage'].centroid_values)

        #print 'ORG', domain.quantities['stage'].centroid_values
        domain.distribute_to_vertices_and_edges()


        #print domain.quantities['stage'].centroid_values
        assert allclose(domain.quantities['stage'].centroid_values,
                        ref_centroid_values)


        #Check that initial limiter doesn't violate cons quan
        assert allclose (domain.quantities['stage'].get_integral(),
                         initial_volume)

        #Evolution
        for t in domain.evolve(yieldstep = 0.05, finaltime = 10):
            volume =  domain.quantities['stage'].get_integral()
            #print t, volume, initial_volume
            assert allclose (volume, initial_volume)

        os.remove(domain.get_name() + '.sww')

    def test_conservation_4(self):
        """Test that stage is conserved globally

        This one uses a larger grid, convoluted bed, reflective bdries and a suitable
        initial condition
        """
        from mesh_factory import rectangular
        from Numeric import array

        #Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order=2
        domain.beta_h = 0.0

        #IC
        def x_slope(x, y):
            z = 0*x
            for i in range(len(x)):
                if x[i] < 0.3:
                    z[i] = x[i]/3
                if 0.3 <= x[i] < 0.5:
                    z[i] = -0.5
                if 0.5 <= x[i] < 0.7:
                    #z[i] = 0.3 #OK with beta == 0.2
                    z[i] = 0.34 #OK with beta == 0.0
                    #z[i] = 0.35#Fails after 80 timesteps with an error
                                #of the order 1.0e-5
                if 0.7 <= x[i]:
                    z[i] = x[i]/3
            return z



        domain.set_quantity('elevation', x_slope)
        domain.set_quantity('friction', 0)
        domain.set_quantity('stage', 0.4) #Steady

        # Boundary conditions (reflective everywhere)
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        domain.check_integrity()

        initial_volume = domain.quantities['stage'].get_integral()
        initial_xmom = domain.quantities['xmomentum'].get_integral()

        import copy
        ref_centroid_values =\
                 copy.copy(domain.quantities['stage'].centroid_values)

        #Test limiter by itself
        domain.distribute_to_vertices_and_edges()

        #Check that initial limiter doesn't violate cons quan
        assert allclose (domain.quantities['stage'].get_integral(),
                         initial_volume)
        #NOTE: This would fail if any initial stage was less than the
        #corresponding bed elevation - but that is reasonable.


        #Evolution
        for t in domain.evolve(yieldstep = 0.05, finaltime = 10.0):
            volume =  domain.quantities['stage'].get_integral()

            #print t, volume, initial_volume

            assert allclose (volume, initial_volume)


        os.remove(domain.get_name() + '.sww')


    def test_conservation_5(self):
        """Test that momentum is conserved globally in
        steady state scenario

        This one uses a slopy bed, dirichlet and reflective bdries
        """
        from mesh_factory import rectangular
        from Numeric import array

        # Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order = 2

        # IC
        def x_slope(x, y):
            return x/3

        domain.set_quantity('elevation', x_slope)
        domain.set_quantity('friction', 0)
        domain.set_quantity('stage', 0.4) # Steady

        # Boundary conditions (reflective everywhere)
        Br = Reflective_boundary(domain)
        Bleft = Dirichlet_boundary([0.5,0,0])
        Bright = Dirichlet_boundary([0.1,0,0])
        domain.set_boundary({'left': Bleft, 'right': Bright,
                             'top': Br, 'bottom': Br})

        domain.check_integrity()

        initial_volume = domain.quantities['stage'].get_integral()
        initial_xmom = domain.quantities['xmomentum'].get_integral()


        # Evolution
        for t in domain.evolve(yieldstep = 0.05, finaltime = 15.0):
            stage =  domain.quantities['stage'].get_integral()
            xmom = domain.quantities['xmomentum'].get_integral()
            ymom = domain.quantities['ymomentum'].get_integral()

            if allclose(t, 6):  # Steady state reached
                steady_xmom = domain.quantities['xmomentum'].get_integral()
                steady_ymom = domain.quantities['ymomentum'].get_integral()
                steady_stage = domain.quantities['stage'].get_integral()

            if t > 6:
                #print '%.2f %14.8f %14.8f' %(t, ymom, steady_ymom)
                msg = 'xmom=%.2f, steady_xmom=%.2f' %(xmom, steady_xmom)
                assert allclose(xmom, steady_xmom), msg
                assert allclose(ymom, steady_ymom)
                assert allclose(stage, steady_stage)


        os.remove(domain.get_name() + '.sww')





    def test_conservation_real(self):
        """Test that momentum is conserved globally
        
        Stephen finally made a test that revealed the problem.
        This test failed with code prior to 25 July 2005
        """

        yieldstep = 0.01
        finaltime = 0.05
        min_depth = 1.0e-2


        import sys
        from os import sep; sys.path.append('..'+sep+'abstract_2d_finite_volumes')
        from mesh_factory import rectangular


        #Create shallow water domain
        points, vertices, boundary = rectangular(10, 10, len1=500, len2=500)
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order = 1
        domain.minimum_allowed_height = min_depth

        # Set initial condition
        class Set_IC:
            """Set an initial condition with a constant value, for x0<x<x1
            """

            def __init__(self, x0=0.25, x1=0.5, h=1.0):
                self.x0 = x0
                self.x1 = x1
                self.h  = h

            def __call__(self, x, y):
                return self.h*((x>self.x0)&(x<self.x1))


        domain.set_quantity('stage', Set_IC(200.0,300.0,5.0))


        #Boundaries
        R = Reflective_boundary(domain)
        domain.set_boundary( {'left': R, 'right': R, 'top':R, 'bottom': R})

        ref = domain.quantities['stage'].get_integral()

        # Evolution
        for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):
            pass
            #print 'Integral stage = ',\
            #      domain.quantities['stage'].get_integral(),\
            #      ' Time = ',domain.time


        now = domain.quantities['stage'].get_integral()

        msg = 'Stage not conserved: was %f, now %f' %(ref, now)
        assert allclose(ref, now), msg 

        os.remove(domain.get_name() + '.sww')

    def test_second_order_flat_bed_onestep(self):

        from mesh_factory import rectangular
        from Numeric import array

        #Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order = 2
        domain.beta_w      = 0.9
        domain.beta_w_dry  = 0.9
        domain.beta_uh     = 0.9
        domain.beta_uh_dry = 0.9
        domain.beta_vh     = 0.9
        domain.beta_vh_dry = 0.9
        domain.H0 = 0 # Backwards compatibility (6/2/7)        
        
        # Boundary conditions
        Br = Reflective_boundary(domain)
        Bd = Dirichlet_boundary([0.1, 0., 0.])
        domain.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})

        domain.check_integrity()

        # Evolution
        for t in domain.evolve(yieldstep = 0.05, finaltime = 0.05):
            pass# domain.write_time()

        # Data from earlier version of abstract_2d_finite_volumes
        assert allclose(domain.min_timestep, 0.0396825396825)
        assert allclose(domain.max_timestep, 0.0396825396825)

        assert allclose(domain.quantities['stage'].centroid_values[:12],
                        [0.00171396, 0.02656103, 0.00241523, 0.02656103,
                        0.00241523, 0.02656103, 0.00241523, 0.02656103,
                        0.00241523, 0.02656103, 0.00241523, 0.0272623])

        domain.distribute_to_vertices_and_edges()

        assert allclose(domain.quantities['stage'].vertex_values[:12,0],
                        [0.0001714, 0.02656103, 0.00024152,
                        0.02656103, 0.00024152, 0.02656103,
                        0.00024152, 0.02656103, 0.00024152,
                        0.02656103, 0.00024152, 0.0272623])

        assert allclose(domain.quantities['stage'].vertex_values[:12,1],
                        [0.00315012, 0.02656103, 0.00024152, 0.02656103,
                         0.00024152, 0.02656103, 0.00024152, 0.02656103,
                         0.00024152, 0.02656103, 0.00040506, 0.0272623])

        assert allclose(domain.quantities['stage'].vertex_values[:12,2],
                        [0.00182037, 0.02656103, 0.00676264,
                         0.02656103, 0.00676264, 0.02656103,
                         0.00676264, 0.02656103, 0.00676264,
                         0.02656103, 0.0065991, 0.0272623])

        assert allclose(domain.quantities['xmomentum'].centroid_values[:12],
                        [0.00113961, 0.01302432, 0.00148672,
                         0.01302432, 0.00148672, 0.01302432,
                         0.00148672, 0.01302432, 0.00148672 ,
                         0.01302432, 0.00148672, 0.01337143])

        assert allclose(domain.quantities['ymomentum'].centroid_values[:12],
                        [-2.91240050e-004 , 1.22721531e-004,
                         -1.22721531e-004, 1.22721531e-004 ,
                         -1.22721531e-004, 1.22721531e-004,
                         -1.22721531e-004 , 1.22721531e-004,
                         -1.22721531e-004, 1.22721531e-004,
                         -1.22721531e-004, -4.57969873e-005])

        os.remove(domain.get_name() + '.sww')


    def test_second_order_flat_bed_moresteps(self):

        from mesh_factory import rectangular
        from Numeric import array

        #Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order=2

        # Boundary conditions
        Br = Reflective_boundary(domain)
        Bd = Dirichlet_boundary([0.1, 0., 0.])
        domain.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})

        domain.check_integrity()

        #Evolution
        for t in domain.evolve(yieldstep = 0.05, finaltime = 0.1):
            pass

        #Data from earlier version of abstract_2d_finite_volumes
        #assert allclose(domain.min_timestep, 0.0396825396825)
        #assert allclose(domain.max_timestep, 0.0396825396825)
        #print domain.quantities['stage'].centroid_values

        os.remove(domain.get_name() + '.sww')        


    def test_flatbed_first_order(self):
        from mesh_factory import rectangular
        from Numeric import array

        #Create basic mesh
        N = 8
        points, vertices, boundary = rectangular(N, N)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order=1
        domain.H0 = 0 # Backwards compatibility (6/2/7)        

        # Boundary conditions
        Br = Reflective_boundary(domain)
        Bd = Dirichlet_boundary([0.2,0.,0.])

        domain.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})
        domain.check_integrity()


        #Evolution
        for t in domain.evolve(yieldstep = 0.02, finaltime = 0.5):
            pass
            #domain.write_time()

        #FIXME: These numbers were from version before 25/10
        #assert allclose(domain.min_timestep, 0.0140413643926)
        #assert allclose(domain.max_timestep, 0.0140947355753)

        for i in range(3):
            #assert allclose(domain.quantities['stage'].edge_values[:4,i],
            #                [0.10730244,0.12337617,0.11200126,0.12605666])

            assert allclose(domain.quantities['xmomentum'].edge_values[:4,i],
                            [0.07610894,0.06901572,0.07284461,0.06819712])

            #assert allclose(domain.quantities['ymomentum'].edge_values[:4,i],
            #                [-0.0060238, -0.00157404, -0.00309633, -0.0001637])


        os.remove(domain.get_name() + '.sww')            

    def test_flatbed_second_order(self):
        from mesh_factory import rectangular
        from Numeric import array

        #Create basic mesh
        N = 8
        points, vertices, boundary = rectangular(N, N)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order=2
        domain.beta_w      = 0.9
        domain.beta_w_dry  = 0.9
        domain.beta_uh     = 0.9
        domain.beta_uh_dry = 0.9
        domain.beta_vh     = 0.9
        domain.beta_vh_dry = 0.9        
        #domain.minimum_allowed_height = 0.0 #Makes it like the 'oldstyle' balance
        domain.H0 = 0 # Backwards compatibility (6/2/7)
        domain.use_centroid_velocities = False # Backwards compatibility (8/5/8)
        domain.set_maximum_allowed_speed(1.0)        

        # Boundary conditions
        Br = Reflective_boundary(domain)
        Bd = Dirichlet_boundary([0.2,0.,0.])

        domain.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})
        domain.check_integrity()

        # Evolution
        for t in domain.evolve(yieldstep = 0.01, finaltime = 0.03):
            pass

        msg = 'min step was %f instead of %f' %(domain.min_timestep,
                                                0.0210448446782) 

        assert allclose(domain.min_timestep, 0.0210448446782), msg
        assert allclose(domain.max_timestep, 0.0210448446782)

        #print domain.quantities['stage'].vertex_values[:4,0]
        #print domain.quantities['xmomentum'].vertex_values[:4,0]
        #print domain.quantities['ymomentum'].vertex_values[:4,0]

        #FIXME: These numbers were from version before 25/10
        #assert allclose(domain.quantities['stage'].vertex_values[:4,0],
        #                [0.00101913,0.05352143,0.00104852,0.05354394])

        #FIXME: These numbers were from version before 21/3/6 - 
        #could be recreated by setting maximum_allowed_speed to 0 maybe 
        #assert allclose(domain.quantities['xmomentum'].vertex_values[:4,0],
        #                [ 0.00064835, 0.03685719, 0.00085073, 0.03687313]) 
        
        assert allclose(domain.quantities['xmomentum'].vertex_values[:4,0],
                        [ 0.00090581, 0.03685719, 0.00088303, 0.03687313])
                        
                        

        #assert allclose(domain.quantities['xmomentum'].vertex_values[:4,0],
        #                [0.00090581,0.03685719,0.00088303,0.03687313])

        assert allclose(domain.quantities['ymomentum'].vertex_values[:4,0],
                        [-0.00139463,0.0006156,-0.00060364,0.00061827])


        os.remove(domain.get_name() + '.sww')

        
    def test_flatbed_second_order_vmax_0(self):
        from mesh_factory import rectangular
        from Numeric import array

        #Create basic mesh
        N = 8
        points, vertices, boundary = rectangular(N, N)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order=2
        domain.beta_w      = 0.9
        domain.beta_w_dry  = 0.9
        domain.beta_uh     = 0.9
        domain.beta_uh_dry = 0.9
        domain.beta_vh     = 0.9
        domain.beta_vh_dry = 0.9        
        domain.maximum_allowed_speed = 0.0 #Makes it like the 'oldstyle'
        domain.H0 = 0 # Backwards compatibility (6/2/7)
        domain.use_centroid_velocities = False # Backwards compatibility (8/5/8)        

        # Boundary conditions
        Br = Reflective_boundary(domain)
        Bd = Dirichlet_boundary([0.2,0.,0.])

        domain.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})
        domain.check_integrity()

        #Evolution
        for t in domain.evolve(yieldstep = 0.01, finaltime = 0.03):
            pass


        assert allclose(domain.min_timestep, 0.0210448446782)
        assert allclose(domain.max_timestep, 0.0210448446782)

        #FIXME: These numbers were from version before 21/3/6 - 
        #could be recreated by setting maximum_allowed_speed to 0 maybe 
        assert allclose(domain.quantities['xmomentum'].vertex_values[:4,0],
                        [ 0.00064835, 0.03685719, 0.00085073, 0.03687313]) 
        

        assert allclose(domain.quantities['ymomentum'].vertex_values[:4,0],
                        [-0.00139463,0.0006156,-0.00060364,0.00061827])


        os.remove(domain.get_name() + '.sww')

        

    def test_flatbed_second_order_distribute(self):
        #Use real data from anuga.abstract_2d_finite_volumes 2
        #painfully setup and extracted.
        from mesh_factory import rectangular
        from Numeric import array

        #Create basic mesh
        N = 8
        points, vertices, boundary = rectangular(N, N)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order=domain._order_=2
        domain.beta_w      = 0.9
        domain.beta_w_dry  = 0.9
        domain.beta_uh     = 0.9
        domain.beta_uh_dry = 0.9
        domain.beta_vh     = 0.9
        domain.beta_vh_dry = 0.9
        domain.H0 = 0 # Backwards compatibility (6/2/7)
        domain.use_centroid_velocities = False # Backwards compatibility (8/5/8)        
        domain.set_maximum_allowed_speed(1.0)        

        # Boundary conditions
        Br = Reflective_boundary(domain)
        Bd = Dirichlet_boundary([0.2,0.,0.])

        domain.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})
        domain.check_integrity()



        for V in [False, True]:
            if V:
                #Set centroids as if system had been evolved
                L = zeros(2*N*N, Float)
                L[:32] = [7.21205592e-003, 5.35214298e-002, 1.00910824e-002,
                          5.35439433e-002, 1.00910824e-002, 5.35439433e-002,
                          1.00910824e-002, 5.35439433e-002, 1.00910824e-002,
                          5.35439433e-002, 1.00910824e-002, 5.35439433e-002,
                          1.00910824e-002, 5.35393928e-002, 1.02344264e-002,
                          5.59605058e-002, 0.00000000e+000, 3.31027800e-004,
                          0.00000000e+000, 4.37962142e-005, 0.00000000e+000,
                          4.37962142e-005, 0.00000000e+000, 4.37962142e-005,
                          0.00000000e+000, 4.37962142e-005, 0.00000000e+000,
                          4.37962142e-005, 0.00000000e+000, 4.37962142e-005,
                          0.00000000e+000, 5.57305948e-005]

                X = zeros(2*N*N, Float)
                X[:32] = [6.48351607e-003, 3.68571894e-002, 8.50733285e-003,
                          3.68731327e-002, 8.50733285e-003, 3.68731327e-002,
                          8.50733285e-003, 3.68731327e-002, 8.50733285e-003,
                          3.68731327e-002, 8.50733285e-003, 3.68731327e-002,
                          8.50733285e-003, 3.68693861e-002, 8.65220973e-003,
                          3.85055387e-002, 0.00000000e+000, 2.86060840e-004,
                          0.00000000e+000, 3.58905503e-005, 0.00000000e+000,
                          3.58905503e-005, 0.00000000e+000, 3.58905503e-005,
                          0.00000000e+000, 3.58905503e-005, 0.00000000e+000,
                          3.58905503e-005, 0.00000000e+000, 3.58905503e-005,
                          0.00000000e+000, 4.57662812e-005]

                Y = zeros(2*N*N, Float)
                Y[:32]=[-1.39463104e-003, 6.15600298e-004, -6.03637382e-004,
                        6.18272251e-004, -6.03637382e-004, 6.18272251e-004,
                        -6.03637382e-004, 6.18272251e-004, -6.03637382e-004,
                        6.18272251e-004, -6.03637382e-004, 6.18272251e-004,
                        -6.03637382e-004, 6.18599320e-004, -6.74622797e-004,
                        -1.48934756e-004, 0.00000000e+000, -5.35079969e-005,
                        0.00000000e+000, -2.57264987e-005, 0.00000000e+000,
                        -2.57264987e-005, 0.00000000e+000, -2.57264987e-005,
                        0.00000000e+000, -2.57264987e-005, 0.00000000e+000,
                        -2.57264987e-005, 0.00000000e+000, -2.57264987e-005,
                        0.00000000e+000, -2.57635178e-005]


                domain.set_quantity('stage', L, location='centroids')
                domain.set_quantity('xmomentum', X, location='centroids')
                domain.set_quantity('ymomentum', Y, location='centroids')

                domain.check_integrity()
            else:
                #Evolution
                for t in domain.evolve(yieldstep = 0.01, finaltime = 0.03):
                    pass
                assert allclose(domain.min_timestep, 0.0210448446782)
                assert allclose(domain.max_timestep, 0.0210448446782)


            #Centroids were correct but not vertices.
            #Hence the check of distribute below.
            assert allclose(domain.quantities['stage'].centroid_values[:4],
                            [0.00721206,0.05352143,0.01009108,0.05354394])

            assert allclose(domain.quantities['xmomentum'].centroid_values[:4],
                            [0.00648352,0.03685719,0.00850733,0.03687313])

            assert allclose(domain.quantities['ymomentum'].centroid_values[:4],
                            [-0.00139463,0.0006156,-0.00060364,0.00061827])

            #print 'C17=', domain.quantities['xmomentum'].centroid_values[17]
            #print 'C19=', domain.quantities['xmomentum'].centroid_values[19]

            #assert allclose(domain.quantities['xmomentum'].centroid_values[17],0.00028606084)
            ##print domain.quantities['xmomentum'].centroid_values[17], V
            ##print
            if not V:
                #FIXME: These numbers were from version before 21/3/6 - 
                #could be recreated by setting maximum_allowed_speed to 0 maybe           
                
                #assert allclose(domain.quantities['xmomentum'].centroid_values[17], 0.0)                
                assert allclose(domain.quantities['xmomentum'].centroid_values[17], 0.000286060839592)                          

            else:
                assert allclose(domain.quantities['xmomentum'].centroid_values[17], 0.00028606084)

            import copy
            XX = copy.copy(domain.quantities['xmomentum'].centroid_values)
            assert allclose(domain.quantities['xmomentum'].centroid_values, XX)

            domain.distribute_to_vertices_and_edges()

            #assert allclose(domain.quantities['xmomentum'].centroid_values, XX)

            #assert allclose(domain.quantities['xmomentum'].centroid_values[17],
            #                0.0)
            assert allclose(domain.quantities['xmomentum'].centroid_values[17], 0.000286060839592)                                                  


            #FIXME: These numbers were from version before 25/10
            #assert allclose(domain.quantities['stage'].vertex_values[:4,0],
            #                [0.00101913,0.05352143,0.00104852,0.05354394])

            assert allclose(domain.quantities['ymomentum'].vertex_values[:4,0],
                            [-0.00139463,0.0006156,-0.00060364,0.00061827])


            assert allclose(domain.quantities['xmomentum'].vertex_values[:4,0],
                            [0.00090581,0.03685719,0.00088303,0.03687313])


            #NB NO longer relvant:

            #This was the culprit. First triangles vertex 0 had an
            #x-momentum of 0.0064835 instead of 0.00090581 and
            #third triangle had 0.00850733 instead of 0.00088303
            #print domain.quantities['xmomentum'].vertex_values[:4,0]

            #print domain.quantities['xmomentum'].vertex_values[:4,0]
            #assert allclose(domain.quantities['xmomentum'].vertex_values[:4,0],
            #                [0.00090581,0.03685719,0.00088303,0.03687313])

        os.remove(domain.get_name() + '.sww')



    def test_bedslope_problem_first_order(self):

        from mesh_factory import rectangular
        from Numeric import array

        #Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order = 1

        #Bed-slope and friction
        def x_slope(x, y):
            return -x/3

        domain.set_quantity('elevation', x_slope)

        # Boundary conditions
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        #Initial condition
        #domain.set_quantity('stage', Constant_height(x_slope, 0.05))
        domain.set_quantity('stage', expression='elevation+0.05')
        domain.check_integrity()

        #Evolution
        for t in domain.evolve(yieldstep = 0.05, finaltime = 0.05):
            pass# domain.write_time()

        # FIXME (Ole): Need some other assertion here! 
        #print domain.min_timestep, domain.max_timestep    
        #assert allclose(domain.min_timestep, 0.050010003001)
        #assert allclose(domain.max_timestep, 0.050010003001)


        os.remove(domain.get_name() + '.sww')
        
    def test_bedslope_problem_first_order_moresteps(self):

        from mesh_factory import rectangular
        from Numeric import array

        #Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order = 1
        domain.beta_h = 0.0 # Use first order in h-limiter
        
        # FIXME (Ole): Need tests where these two are commented out
        domain.H0 = 0        # Backwards compatibility (6/2/7)        
        domain.tight_slope_limiters = 0 # Backwards compatibility (14/4/7)
        domain.use_centroid_velocities = 0 # Backwards compatibility (7/5/8)                 

        #Bed-slope and friction
        def x_slope(x, y):
            return -x/3

        domain.set_quantity('elevation', x_slope)

        # Boundary conditions
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        #Initial condition
        domain.set_quantity('stage', expression='elevation+0.05')
        domain.check_integrity()

        #Evolution
        for t in domain.evolve(yieldstep = 0.05, finaltime = 0.5):
            pass# domain.write_time()

        #Data from earlier version of abstract_2d_finite_volumes
        #print domain.quantities['stage'].centroid_values

        assert allclose(domain.quantities['stage'].centroid_values,
                        [-0.02998628, -0.01520652, -0.03043492,
                        -0.0149132, -0.03004706, -0.01476251,
                        -0.0298215, -0.01467976, -0.02988158,
                        -0.01474662, -0.03036161, -0.01442995,
                        -0.07624583, -0.06297061, -0.07733792,
                        -0.06342237, -0.07695439, -0.06289595,
                        -0.07635559, -0.0626065, -0.07633628,
                        -0.06280072, -0.07739632, -0.06386738,
                        -0.12161738, -0.11028239, -0.1223796,
                        -0.11095953, -0.12189744, -0.11048616,
                        -0.12074535, -0.10987605, -0.12014311,
                        -0.10976691, -0.12096859, -0.11087692,
                        -0.16868259, -0.15868061, -0.16801135,
                        -0.1588003, -0.16674343, -0.15813323,
                        -0.16457595, -0.15693826, -0.16281096,
                        -0.15585154, -0.16283873, -0.15540068,
                        -0.17450362, -0.19919913, -0.18062882,
                        -0.19764131, -0.17783111, -0.19407213,
                        -0.1736915, -0.19053624, -0.17228678,
                        -0.19105634, -0.17920133, -0.1968828,
                        -0.14244395, -0.14604641, -0.14473537,
                        -0.1506107, -0.14510055, -0.14919522,
                        -0.14175896, -0.14560798, -0.13911658,
                        -0.14439383, -0.13924047, -0.14829043])

        os.remove(domain.get_name() + '.sww')
        
    def test_bedslope_problem_second_order_one_step(self):

        from mesh_factory import rectangular
        from Numeric import array

        #Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order=2
        domain.beta_w      = 0.9
        domain.beta_w_dry  = 0.9
        domain.beta_uh     = 0.9
        domain.beta_uh_dry = 0.9
        domain.beta_vh     = 0.9
        domain.beta_vh_dry = 0.9

        
        # FIXME (Ole): Need tests where this is commented out
        domain.tight_slope_limiters = 0 # Backwards compatibility (14/4/7)
        domain.use_centroid_velocities = 0 # Backwards compatibility (7/5/8)      
        
        #Bed-slope and friction at vertices (and interpolated elsewhere)
        def x_slope(x, y):
            return -x/3

        domain.set_quantity('elevation', x_slope)

        # Boundary conditions
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        #Initial condition
        domain.set_quantity('stage', expression='elevation+0.05')
        domain.check_integrity()

        assert allclose(domain.quantities['stage'].centroid_values,
                        [0.01296296, 0.03148148, 0.01296296,
                        0.03148148, 0.01296296, 0.03148148,
                        0.01296296, 0.03148148, 0.01296296,
                        0.03148148, 0.01296296, 0.03148148,
                        -0.04259259, -0.02407407, -0.04259259,
                        -0.02407407, -0.04259259, -0.02407407,
                        -0.04259259, -0.02407407, -0.04259259,
                        -0.02407407, -0.04259259, -0.02407407,
                        -0.09814815, -0.07962963, -0.09814815,
                        -0.07962963, -0.09814815, -0.07962963,
                        -0.09814815, -0.07962963, -0.09814815,
                        -0.07962963, -0.09814815, -0.07962963,
                        -0.1537037 , -0.13518519, -0.1537037,
                        -0.13518519, -0.1537037, -0.13518519,
                        -0.1537037 , -0.13518519, -0.1537037,
                        -0.13518519, -0.1537037, -0.13518519,
                        -0.20925926, -0.19074074, -0.20925926,
                        -0.19074074, -0.20925926, -0.19074074,
                        -0.20925926, -0.19074074, -0.20925926,
                        -0.19074074, -0.20925926, -0.19074074,
                        -0.26481481, -0.2462963, -0.26481481,
                        -0.2462963, -0.26481481, -0.2462963,
                        -0.26481481, -0.2462963, -0.26481481,
                        -0.2462963, -0.26481481, -0.2462963])


        #print domain.quantities['stage'].extrapolate_second_order()
        #domain.distribute_to_vertices_and_edges()
        #print domain.quantities['stage'].vertex_values[:,0]

        #Evolution
        for t in domain.evolve(yieldstep = 0.05, finaltime = 0.05):
            #domain.write_time()
            pass


        #print domain.quantities['stage'].centroid_values
        assert allclose(domain.quantities['stage'].centroid_values,
                 [0.01290985, 0.02356019, 0.01619096, 0.02356019, 0.01619096,
                  0.02356019, 0.01619096, 0.02356019, 0.01619096, 0.02356019,
                  0.01619096, 0.0268413, -0.04411074, -0.0248011, -0.04186556,
                  -0.0248011, -0.04186556, -0.0248011, -0.04186556, -0.0248011,
                  -0.04186556, -0.0248011, -0.04186556, -0.02255593,
                  -0.09966629, -0.08035666, -0.09742112, -0.08035666,
                  -0.09742112, -0.08035666, -0.09742112, -0.08035666,
                  -0.09742112, -0.08035666, -0.09742112, -0.07811149,
                  -0.15522185, -0.13591222, -0.15297667, -0.13591222,
                  -0.15297667, -0.13591222, -0.15297667, -0.13591222,
                  -0.15297667, -0.13591222, -0.15297667, -0.13366704,
                  -0.2107774, -0.19146777, -0.20853223, -0.19146777,
                  -0.20853223, -0.19146777, -0.20853223, -0.19146777,
                  -0.20853223, -0.19146777, -0.20853223, -0.1892226,
                  -0.26120669, -0.24776246, -0.25865535, -0.24776246,
                  -0.25865535, -0.24776246, -0.25865535, -0.24776246,
                  -0.25865535, -0.24776246, -0.25865535, -0.24521113])

        os.remove(domain.get_name() + '.sww')

    def test_bedslope_problem_second_order_two_steps(self):

        from mesh_factory import rectangular
        from Numeric import array

        #Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order=2
        domain.beta_w      = 0.9
        domain.beta_w_dry  = 0.9
        domain.beta_uh     = 0.9
        domain.beta_uh_dry = 0.9
        domain.beta_vh     = 0.9
        domain.beta_vh_dry = 0.9
        domain.beta_h = 0.0 #Use first order in h-limiter
        
        # FIXME (Ole): Need tests where this is commented out
        domain.tight_slope_limiters = 0 # Backwards compatibility (14/4/7)                 
        domain.H0 = 0 # Backwards compatibility (6/2/7)
        domain.use_centroid_velocities = 0 # Backwards compatibility (7/5/8)
        

        #Bed-slope and friction at vertices (and interpolated elsewhere)
        def x_slope(x, y):
            return -x/3

        domain.set_quantity('elevation', x_slope)

        # Boundary conditions
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        #Initial condition
        domain.set_quantity('stage', expression='elevation+0.05')
        domain.check_integrity()

        assert allclose(domain.quantities['stage'].centroid_values,
                        [0.01296296, 0.03148148, 0.01296296,
                        0.03148148, 0.01296296, 0.03148148,
                        0.01296296, 0.03148148, 0.01296296,
                        0.03148148, 0.01296296, 0.03148148,
                        -0.04259259, -0.02407407, -0.04259259,
                        -0.02407407, -0.04259259, -0.02407407,
                        -0.04259259, -0.02407407, -0.04259259,
                        -0.02407407, -0.04259259, -0.02407407,
                        -0.09814815, -0.07962963, -0.09814815,
                        -0.07962963, -0.09814815, -0.07962963,
                        -0.09814815, -0.07962963, -0.09814815,
                        -0.07962963, -0.09814815, -0.07962963,
                        -0.1537037 , -0.13518519, -0.1537037,
                        -0.13518519, -0.1537037, -0.13518519,
                        -0.1537037 , -0.13518519, -0.1537037,
                        -0.13518519, -0.1537037, -0.13518519,
                        -0.20925926, -0.19074074, -0.20925926,
                        -0.19074074, -0.20925926, -0.19074074,
                        -0.20925926, -0.19074074, -0.20925926,
                        -0.19074074, -0.20925926, -0.19074074,
                        -0.26481481, -0.2462963, -0.26481481,
                        -0.2462963, -0.26481481, -0.2462963,
                        -0.26481481, -0.2462963, -0.26481481,
                        -0.2462963, -0.26481481, -0.2462963])


        #print domain.quantities['stage'].extrapolate_second_order()
        #domain.distribute_to_vertices_and_edges()
        #print domain.quantities['stage'].vertex_values[:,0]

        #Evolution
        for t in domain.evolve(yieldstep = 0.05, finaltime = 0.1):
            pass


        #Data from earlier version of abstract_2d_finite_volumes ft=0.1
        assert allclose(domain.min_timestep, 0.0376895634803)
        assert allclose(domain.max_timestep, 0.0415635655309)


        assert allclose(domain.quantities['stage'].centroid_values,
                        [0.00855788, 0.01575204, 0.00994606, 0.01717072,
                         0.01005985, 0.01716362, 0.01005985, 0.01716299,
                         0.01007098, 0.01736248, 0.01216452, 0.02026776,
                         -0.04462374, -0.02479045, -0.04199789, -0.0229465,
                         -0.04184033, -0.02295693, -0.04184013, -0.02295675,
                         -0.04184486, -0.0228168, -0.04028876, -0.02036486,
                         -0.10029444, -0.08170809, -0.09772846, -0.08021704,
                         -0.09760006, -0.08022143, -0.09759984, -0.08022124,
                         -0.09760261, -0.08008893, -0.09603914, -0.07758209,
                         -0.15584152, -0.13723138, -0.15327266, -0.13572906,
                         -0.15314427, -0.13573349, -0.15314405, -0.13573331,
                         -0.15314679, -0.13560104, -0.15158523, -0.13310701,
                         -0.21208605, -0.19283913, -0.20955631, -0.19134189,
                         -0.20942821, -0.19134598, -0.20942799, -0.1913458,
                         -0.20943005, -0.19120952, -0.20781177, -0.18869401,
                         -0.25384082, -0.2463294, -0.25047649, -0.24464654,
                         -0.25031159, -0.24464253, -0.25031112, -0.24464253,
                         -0.25031463, -0.24454764, -0.24885323, -0.24286438])


        os.remove(domain.get_name() + '.sww')

    def test_bedslope_problem_second_order_two_yieldsteps(self):

        from mesh_factory import rectangular
        from Numeric import array

        #Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order=2
        domain.beta_w      = 0.9
        domain.beta_w_dry  = 0.9
        domain.beta_uh     = 0.9
        domain.beta_uh_dry = 0.9
        domain.beta_vh     = 0.9
        domain.beta_vh_dry = 0.9
        domain.beta_h = 0.0 #Use first order in h-limiter
        
        # FIXME (Ole): Need tests where this is commented out
        domain.tight_slope_limiters = 0 # Backwards compatibility (14/4/7)                 
        domain.H0 = 0 # Backwards compatibility (6/2/7)
        domain.use_centroid_velocities = 0 # Backwards compatibility (7/5/8)
        

        #Bed-slope and friction at vertices (and interpolated elsewhere)
        def x_slope(x, y):
            return -x/3

        domain.set_quantity('elevation', x_slope)

        # Boundary conditions
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        #Initial condition
        domain.set_quantity('stage', expression='elevation+0.05')
        domain.check_integrity()

        assert allclose(domain.quantities['stage'].centroid_values,
                        [0.01296296, 0.03148148, 0.01296296,
                        0.03148148, 0.01296296, 0.03148148,
                        0.01296296, 0.03148148, 0.01296296,
                        0.03148148, 0.01296296, 0.03148148,
                        -0.04259259, -0.02407407, -0.04259259,
                        -0.02407407, -0.04259259, -0.02407407,
                        -0.04259259, -0.02407407, -0.04259259,
                        -0.02407407, -0.04259259, -0.02407407,
                        -0.09814815, -0.07962963, -0.09814815,
                        -0.07962963, -0.09814815, -0.07962963,
                        -0.09814815, -0.07962963, -0.09814815,
                        -0.07962963, -0.09814815, -0.07962963,
                        -0.1537037 , -0.13518519, -0.1537037,
                        -0.13518519, -0.1537037, -0.13518519,
                        -0.1537037 , -0.13518519, -0.1537037,
                        -0.13518519, -0.1537037, -0.13518519,
                        -0.20925926, -0.19074074, -0.20925926,
                        -0.19074074, -0.20925926, -0.19074074,
                        -0.20925926, -0.19074074, -0.20925926,
                        -0.19074074, -0.20925926, -0.19074074,
                        -0.26481481, -0.2462963, -0.26481481,
                        -0.2462963, -0.26481481, -0.2462963,
                        -0.26481481, -0.2462963, -0.26481481,
                        -0.2462963, -0.26481481, -0.2462963])


        #print domain.quantities['stage'].extrapolate_second_order()
        #domain.distribute_to_vertices_and_edges()
        #print domain.quantities['stage'].vertex_values[:,0]

        #Evolution
        for t in domain.evolve(yieldstep = 0.05, finaltime = 0.1):   #0.05??
            #domain.write_time()
            pass



        assert allclose(domain.quantities['stage'].centroid_values,
                 [0.00855788, 0.01575204, 0.00994606, 0.01717072, 0.01005985,
                  0.01716362, 0.01005985, 0.01716299, 0.01007098, 0.01736248,
                  0.01216452, 0.02026776, -0.04462374, -0.02479045, -0.04199789,
                  -0.0229465, -0.04184033, -0.02295693, -0.04184013,
                  -0.02295675, -0.04184486, -0.0228168, -0.04028876,
                  -0.02036486, -0.10029444, -0.08170809, -0.09772846,
                  -0.08021704, -0.09760006, -0.08022143, -0.09759984,
                  -0.08022124, -0.09760261, -0.08008893, -0.09603914,
                  -0.07758209, -0.15584152, -0.13723138, -0.15327266,
                  -0.13572906, -0.15314427, -0.13573349, -0.15314405,
                  -0.13573331, -0.15314679, -0.13560104, -0.15158523,
                  -0.13310701, -0.21208605, -0.19283913, -0.20955631,
                  -0.19134189, -0.20942821, -0.19134598, -0.20942799,
                  -0.1913458, -0.20943005, -0.19120952, -0.20781177,
                  -0.18869401, -0.25384082, -0.2463294, -0.25047649,
                  -0.24464654, -0.25031159, -0.24464253, -0.25031112,
                  -0.24464253, -0.25031463, -0.24454764, -0.24885323,
                  -0.24286438])

        os.remove(domain.get_name() + '.sww')

    def test_bedslope_problem_second_order_more_steps(self):

        from mesh_factory import rectangular
        from Numeric import array

        #Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order=2
        domain.beta_w      = 0.9
        domain.beta_w_dry  = 0.9
        domain.beta_uh     = 0.9
        domain.beta_uh_dry = 0.9
        domain.beta_vh     = 0.9
        domain.beta_vh_dry = 0.9
        domain.beta_h = 0.0 #Use first order in h-limiter
        
        
        # FIXME (Ole): Need tests where these two are commented out
        domain.H0 = 0        # Backwards compatibility (6/2/7)        
        domain.tight_slope_limiters = 0 # Backwards compatibility (14/4/7)                 
        domain.use_centroid_velocities = 0 # Backwards compatibility (7/5/8)
        
                

        #Bed-slope and friction at vertices (and interpolated elsewhere)
        def x_slope(x, y):
            return -x/3

        domain.set_quantity('elevation', x_slope)

        # Boundary conditions
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        #Initial condition
        domain.set_quantity('stage', expression = 'elevation + 0.05')
        domain.check_integrity()

        assert allclose(domain.quantities['stage'].centroid_values,
                        [0.01296296, 0.03148148, 0.01296296,
                        0.03148148, 0.01296296, 0.03148148,
                        0.01296296, 0.03148148, 0.01296296,
                        0.03148148, 0.01296296, 0.03148148,
                        -0.04259259, -0.02407407, -0.04259259,
                        -0.02407407, -0.04259259, -0.02407407,
                        -0.04259259, -0.02407407, -0.04259259,
                        -0.02407407, -0.04259259, -0.02407407,
                        -0.09814815, -0.07962963, -0.09814815,
                        -0.07962963, -0.09814815, -0.07962963,
                        -0.09814815, -0.07962963, -0.09814815,
                        -0.07962963, -0.09814815, -0.07962963,
                        -0.1537037 , -0.13518519, -0.1537037,
                        -0.13518519, -0.1537037, -0.13518519,
                        -0.1537037 , -0.13518519, -0.1537037,
                        -0.13518519, -0.1537037, -0.13518519,
                        -0.20925926, -0.19074074, -0.20925926,
                        -0.19074074, -0.20925926, -0.19074074,
                        -0.20925926, -0.19074074, -0.20925926,
                        -0.19074074, -0.20925926, -0.19074074,
                        -0.26481481, -0.2462963, -0.26481481,
                        -0.2462963, -0.26481481, -0.2462963,
                        -0.26481481, -0.2462963, -0.26481481,
                        -0.2462963, -0.26481481, -0.2462963])


        #print domain.quantities['stage'].extrapolate_second_order()
        #domain.distribute_to_vertices_and_edges()
        #print domain.quantities['stage'].vertex_values[:,0]

        #Evolution
        for t in domain.evolve(yieldstep = 0.05, finaltime = 0.5):

            # Check that diagnostics works
            msg = domain.timestepping_statistics(track_speeds=True)
            #FIXME(Ole): One might check the contents of msg here.



        assert allclose(domain.quantities['stage'].centroid_values,
     [-0.02907028, -0.01475478, -0.02973417, -0.01447186, -0.02932665, -0.01428285,
      -0.02901975, -0.0141361,  -0.02898816, -0.01418135, -0.02961409, -0.01403487,
      -0.07597998, -0.06252591, -0.07664854, -0.06312532, -0.07638287, -0.06265139,
      -0.07571145, -0.06235231, -0.0756817,  -0.06245309, -0.07652292, -0.06289946,
      -0.12367464, -0.11088981, -0.12237277, -0.11115338, -0.1218934,  -0.1107174,
      -0.12081485, -0.11000491, -0.12038451, -0.11010335, -0.12102113, -0.11012105,
      -0.16909116, -0.15831543, -0.16730214, -0.15786249, -0.1665493,  -0.15697919,
      -0.16496618, -0.15559852, -0.16338679, -0.15509088, -0.16364092, -0.15424423,
      -0.18771107, -0.19903904, -0.18903759, -0.19858437, -0.18701552, -0.19697797,
      -0.1833593,  -0.19505871, -0.1818806,  -0.19418042, -0.18586159, -0.19576946,
      -0.13986873, -0.14170053, -0.14132188, -0.14560674, -0.14095617, -0.14373292,
      -0.13785933, -0.14033364, -0.13592955, -0.13936356, -0.13596008, -0.14216296])

        assert allclose(domain.quantities['xmomentum'].centroid_values,
     [ 0.00831121,  0.00317948,  0.00731797,  0.00334939,  0.00764717,  0.00348053,
       0.00788729,  0.00356522,  0.00780649,  0.00341919,  0.00693525,  0.00310375,
       0.02166196,  0.01421475,  0.02017737,  0.01316839,  0.02037015,  0.01368659,
       0.02106,     0.01399161,  0.02074514,  0.01354935,  0.01887407,  0.0123113,
       0.03775083,  0.02855197,  0.03689337,  0.02759782,  0.03732848,  0.02812072,
       0.03872545,  0.02913348,  0.03880939,  0.02803804,  0.03546499,  0.0260039,
       0.0632131,   0.04730634,  0.0576324,   0.04592336,  0.05790921,  0.04690514,
       0.05986467,  0.04871165,  0.06170068,  0.04811572,  0.05657041,  0.04416292,
       0.08489642,  0.07188097,  0.07835261,  0.06843406,  0.07986412,  0.0698247,
       0.08201071,  0.07216756,  0.08378418,  0.07273624,  0.080399,    0.06645841,
       0.01631548,  0.04691608,  0.0206632,   0.044409,    0.02115518,  0.04560305,
       0.02160608,  0.04663725,  0.02174734,  0.04795559,  0.02281427,  0.05667111])


        assert allclose(domain.quantities['ymomentum'].centroid_values,
     [ 1.45876601e-004, -3.24627393e-004, -1.57572719e-004, -2.92790187e-004,
      -9.90988382e-005, -3.06677335e-004, -1.62493106e-004, -3.71310004e-004,
      -1.99445058e-004, -3.28493467e-004,  6.68217349e-005, -8.42042805e-006,
       5.05093371e-004, -1.42842214e-004, -6.81454718e-005, -5.02084057e-004,
      -8.50583861e-005, -4.65443981e-004, -1.96406564e-004, -5.88889562e-004,
      -2.70160173e-004, -5.35485454e-004,  2.60780997e-004,  3.12145471e-005,
       5.16189608e-004,  1.07069062e-004,  9.29989252e-005, -3.71211119e-004,
       1.16350246e-004, -3.82407830e-004, -1.62077969e-004, -6.30906636e-004,
      -4.74025708e-004, -6.94463009e-004,  6.15092843e-005,  2.22106820e-004,
      -6.29589294e-004,  2.43611937e-004, -5.88125094e-004, -6.94293192e-005,
      -4.17914641e-004,  6.64609019e-005, -7.68334577e-004, -3.40232101e-004,
      -1.67424308e-003, -7.39485066e-004, -1.59966988e-003,  5.68262838e-005,
      -1.48470633e-003, -1.84554882e-003, -2.27200099e-003, -1.67506848e-003,
      -1.95610258e-003, -1.47638801e-003, -1.73779477e-003, -1.85498791e-003,
      -2.01357843e-003, -2.17675471e-003, -1.65783870e-003, -1.15818681e-003,
      -1.18663036e-003, -2.94229849e-003, -3.59309018e-003, -5.13496584e-003,
      -6.17359400e-003, -5.98761937e-003, -6.00540116e-003, -5.01121966e-003,
      -4.50964850e-003, -3.06319963e-003,  6.08950810e-004, -4.79537921e-004])

        os.remove(domain.get_name() + '.sww')



    def NOtest_bedslope_problem_second_order_more_steps_feb_2007(self):
        """test_bedslope_problem_second_order_more_steps_feb_2007

        Test shallow water finite volumes, using parameters from 
        feb 2007 rather than backward compatibility ad infinitum
        
        """
        from mesh_factory import rectangular
        from Numeric import array

        #Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order = 2
        domain.beta_w      = 0.9
        domain.beta_w_dry  = 0.9
        domain.beta_uh     = 0.9
        domain.beta_uh_dry = 0.9
        domain.beta_vh     = 0.9
        domain.beta_vh_dry = 0.9
        domain.beta_h = 0.0 #Use first order in h-limiter
        domain.H0 = 0.001
        domain.tight_slope_limiters = 1

        #Bed-slope and friction at vertices (and interpolated elsewhere)
        def x_slope(x, y):
            return -x/3

        domain.set_quantity('elevation', x_slope)

        # Boundary conditions
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        #Initial condition
        domain.set_quantity('stage', expression = 'elevation + 0.05')
        domain.check_integrity()

        assert allclose(domain.quantities['stage'].centroid_values,
                        [0.01296296, 0.03148148, 0.01296296,
                        0.03148148, 0.01296296, 0.03148148,
                        0.01296296, 0.03148148, 0.01296296,
                        0.03148148, 0.01296296, 0.03148148,
                        -0.04259259, -0.02407407, -0.04259259,
                        -0.02407407, -0.04259259, -0.02407407,
                        -0.04259259, -0.02407407, -0.04259259,
                        -0.02407407, -0.04259259, -0.02407407,
                        -0.09814815, -0.07962963, -0.09814815,
                        -0.07962963, -0.09814815, -0.07962963,
                        -0.09814815, -0.07962963, -0.09814815,
                        -0.07962963, -0.09814815, -0.07962963,
                        -0.1537037 , -0.13518519, -0.1537037,
                        -0.13518519, -0.1537037, -0.13518519,
                        -0.1537037 , -0.13518519, -0.1537037,
                        -0.13518519, -0.1537037, -0.13518519,
                        -0.20925926, -0.19074074, -0.20925926,
                        -0.19074074, -0.20925926, -0.19074074,
                        -0.20925926, -0.19074074, -0.20925926,
                        -0.19074074, -0.20925926, -0.19074074,
                        -0.26481481, -0.2462963, -0.26481481,
                        -0.2462963, -0.26481481, -0.2462963,
                        -0.26481481, -0.2462963, -0.26481481,
                        -0.2462963, -0.26481481, -0.2462963])


        #print domain.quantities['stage'].extrapolate_second_order()
        #domain.distribute_to_vertices_and_edges()
        #print domain.quantities['stage'].vertex_values[:,0]

        #Evolution
        for t in domain.evolve(yieldstep = 0.05, finaltime = 0.5):
            pass


        #print domain.quantities['stage'].centroid_values
            
        assert allclose(domain.quantities['stage'].centroid_values,
         [-0.03348416, -0.01749303, -0.03299091, -0.01739241, -0.03246447, -0.01732016,
          -0.03205390, -0.01717833, -0.03146383, -0.01699831, -0.03076577, -0.01671795,
          -0.07952656, -0.06684763, -0.07721455, -0.06668388, -0.07632976, -0.06600113,
          -0.07523678, -0.06546373, -0.07447040, -0.06508861, -0.07438723, -0.06359288,
          -0.12526729, -0.11205668, -0.12179433, -0.11068104, -0.12048395, -0.10968948,
          -0.11912023, -0.10862628, -0.11784090, -0.10803744, -0.11790629, -0.10742354,
          -0.16859613, -0.15427413, -0.16664444, -0.15464452, -0.16570816, -0.15327556,
          -0.16409162, -0.15204092, -0.16264608, -0.15102139, -0.16162736, -0.14969205,
          -0.18736511, -0.19874036, -0.18811230, -0.19758289, -0.18590182, -0.19580301,
          -0.18234588, -0.19423215, -0.18100376, -0.19380116, -0.18509710, -0.19501636,
          -0.13982382, -0.14166819, -0.14132775, -0.14528694, -0.14096905, -0.14351126,
          -0.13800356, -0.14027920, -0.13613538, -0.13936795, -0.13621902, -0.14204982])

                      
        assert allclose(domain.quantities['xmomentum'].centroid_values,
      [0.00600290,  0.00175780,  0.00591905,  0.00190903,  0.00644462,  0.00203095,
       0.00684561,  0.00225089,  0.00708208,  0.00236235,  0.00649095,  0.00222343,
       0.02068693,  0.01164034,  0.01983343,  0.01159526,  0.02044611,  0.01233252,
       0.02135685,  0.01301289,  0.02161290,  0.01260280,  0.01867612,  0.01133078,
       0.04091313,  0.02668283,  0.03634781,  0.02733469,  0.03767692,  0.02836840,
       0.03906338,  0.02958073,  0.04025669,  0.02953292,  0.03665616,  0.02583565,
       0.06314558,  0.04830935,  0.05663609,  0.04564362,  0.05756200,  0.04739673,
       0.05967379,  0.04919083,  0.06124330,  0.04965808,  0.05879240,  0.04629319,
       0.08220739,  0.06924725,  0.07713556,  0.06782640,  0.07909499,  0.06992544,
       0.08116621,  0.07210181,  0.08281548,  0.07222669,  0.07941059,  0.06755612,
       0.01581588,  0.04533609,  0.02017939,  0.04342565,  0.02073232,  0.04476108,
       0.02117439,  0.04573358,  0.02129473,  0.04694267,  0.02220398,  0.05533458])

       
        assert allclose(domain.quantities['ymomentum'].centroid_values,
     [-7.65882069e-005, -1.46087080e-004, -1.09630102e-004, -7.80950424e-005,
      -1.15922807e-005, -9.09134899e-005, -1.35994542e-004, -1.95673476e-004,
      -4.25779199e-004, -2.95890312e-004, -4.00060341e-004, -9.42021290e-005,
      -3.41372596e-004, -1.54560195e-004, -2.94810038e-004, -1.08844546e-004,
      -6.97240892e-005,  3.50299623e-005, -2.40159184e-004, -2.01805883e-004,
      -7.60732405e-004, -5.10897642e-004, -1.00940001e-003, -1.38037759e-004,
      -1.06169131e-003, -3.12307760e-004, -9.90602307e-004, -4.21634250e-005,
      -6.02424239e-004,  1.52230578e-004, -7.63833035e-004, -1.10273481e-004,
      -1.40187071e-003, -5.57831837e-004, -1.63988285e-003, -2.48018092e-004,
      -1.83309840e-003, -6.19360836e-004, -1.29955242e-003, -3.76237145e-004,
      -1.00613007e-003, -8.63641918e-005, -1.13604124e-003, -3.90589728e-004,
      -1.91457355e-003, -9.43783961e-004, -2.28090840e-003, -5.79107025e-004,
      -1.54091533e-003, -2.39785792e-003, -2.47947427e-003, -2.02694009e-003,
      -2.10441194e-003, -1.82082650e-003, -1.80229336e-003, -2.10418336e-003,
      -1.93104408e-003, -2.23200334e-003, -1.57239706e-003, -1.31486358e-003,
      -1.17564993e-003, -2.85846494e-003, -3.52956754e-003, -5.12658193e-003,
      -6.24238960e-003, -6.01820113e-003, -6.09602201e-003, -5.04787190e-003,
      -4.59373845e-003, -3.01393146e-003,  5.08550095e-004, -4.35896549e-004])

        os.remove(domain.get_name() + '.sww')


    def test_temp_play(self):

        from mesh_factory import rectangular
        from Numeric import array

        #Create basic mesh
        points, vertices, boundary = rectangular(5, 5)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order=2
        domain.beta_w      = 0.9
        domain.beta_w_dry  = 0.9
        domain.beta_uh     = 0.9
        domain.beta_uh_dry = 0.9
        domain.beta_vh     = 0.9
        domain.beta_vh_dry = 0.9
        domain.beta_h = 0.0 #Use first order in h-limiter
        
        # FIXME (Ole): Need tests where these two are commented out
        domain.H0 = 0        # Backwards compatibility (6/2/7)        
        domain.tight_slope_limiters = False # Backwards compatibility (14/4/7)
        domain.use_centroid_velocities = False # Backwards compatibility (7/5/8)
        domain.use_edge_limiter = False # Backwards compatibility (9/5/8)        
        

        #Bed-slope and friction at vertices (and interpolated elsewhere)
        def x_slope(x, y):
            return -x/3

        domain.set_quantity('elevation', x_slope)

        # Boundary conditions
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        #Initial condition
        domain.set_quantity('stage', expression='elevation+0.05')
        domain.check_integrity()

        #Evolution
        for t in domain.evolve(yieldstep = 0.05, finaltime = 0.1):
            pass

        assert allclose(domain.quantities['stage'].centroid_values[:4],
                        [0.00206836, 0.01296714, 0.00363415, 0.01438924])
        #print domain.quantities['xmomentum'].centroid_values[:4]
        assert allclose(domain.quantities['xmomentum'].centroid_values[:4],
                        [0.01360154, 0.00671133, 0.01264578, 0.00648503])
        assert allclose(domain.quantities['ymomentum'].centroid_values[:4],
                        [-1.19201077e-003, -7.23647546e-004,
                        -6.39083123e-005, 6.29815168e-005])

        os.remove(domain.get_name() + '.sww')

    def test_complex_bed(self):
        #No friction is tested here

        from mesh_factory import rectangular
        from Numeric import array

        N = 12
        points, vertices, boundary = rectangular(N, N/2, len1=1.2,len2=0.6,
                                                 origin=(-0.07, 0))


        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order=2


        inflow_stage = 0.1
        Z = Weir(inflow_stage)
        domain.set_quantity('elevation', Z)

        Br = Reflective_boundary(domain)
        Bd = Dirichlet_boundary([inflow_stage, 0.0, 0.0])
        domain.set_boundary({'left': Bd, 'right': Br, 'bottom': Br, 'top': Br})

        domain.set_quantity('stage', expression='elevation')

        for t in domain.evolve(yieldstep = 0.02, finaltime = 0.2):
            pass


        #print domain.quantities['stage'].centroid_values

        #FIXME: These numbers were from version before 25/10
        #assert allclose(domain.quantities['stage'].centroid_values,
# [3.95822638e-002,  5.61022588e-002,  4.66437868e-002,  5.73081011e-002,
#  4.72394613e-002,  5.74684939e-002,  4.74309483e-002,  5.77458084e-002,
#  4.80628177e-002,  5.85656225e-002,  4.90498542e-002,  6.02609831e-002,
#  1.18470315e-002,  1.75136443e-002,  1.18035266e-002,  2.15565695e-002,
#  1.31620268e-002,  2.14351640e-002,  1.32351076e-002,  2.15450687e-002,
#  1.36414028e-002,  2.24274619e-002,  1.51689511e-002,  2.21789655e-002,
# -7.54337535e-003, -6.86362021e-004, -7.74146760e-003, -1.83756530e-003,
# -8.16773628e-003, -4.49916813e-004, -8.08202599e-003, -3.91118720e-004,
# -8.10292716e-003, -3.88584984e-004, -7.35226124e-003,  2.73985295e-004,
#  1.86166683e-001,  8.74070369e-002,  1.86166712e-001,  8.74035875e-002,
#  6.11666935e-002, -3.76173225e-002, -6.38333276e-002, -3.76147365e-002,
#  6.11666725e-002,  8.73846774e-002,  1.86166697e-001,  8.74171550e-002,
# -4.83333333e-002,  1.18333333e-001, -4.83333333e-002,  1.18333333e-001,
# -4.83333333e-002, -6.66666667e-003, -1.73333333e-001, -1.31666667e-001,
# -1.73333333e-001, -6.66666667e-003, -4.83333333e-002,  1.18333333e-001,
# -2.48333333e-001, -2.31666667e-001, -2.48333333e-001, -2.31666667e-001,
# -2.48333333e-001, -2.31666667e-001, -2.48333333e-001, -2.31666667e-001,
# -2.48333333e-001, -2.31666667e-001, -2.48333333e-001, -2.31666667e-001,
# -4.65000000e-001, -3.65000000e-001, -4.65000000e-001, -3.65000000e-001,
# -4.65000000e-001, -3.65000000e-001, -4.65000000e-001, -3.65000000e-001,
# -4.65000000e-001, -3.65000000e-001, -4.65000000e-001, -3.65000000e-001,
# -5.98333333e-001, -5.81666667e-001, -5.98333333e-001, -5.81666667e-001,
# -5.98333333e-001, -5.81666667e-001, -5.98333333e-001, -5.81666667e-001,
# -5.98333333e-001, -5.81666667e-001, -5.98333333e-001, -5.81666667e-001,
# -6.48333333e-001, -6.31666667e-001, -6.48333333e-001, -6.31666667e-001,
# -6.48333333e-001, -6.31666667e-001, -6.48333333e-001, -6.31666667e-001,
# -6.48333333e-001, -6.31666667e-001, -6.48333333e-001, -6.31666667e-001,
# -5.31666667e-001, -5.98333333e-001, -5.31666667e-001, -5.98333333e-001,
# -5.31666667e-001, -5.98333333e-001, -5.31666667e-001, -5.98333333e-001,
# -5.31666667e-001, -5.98333333e-001, -5.31666667e-001, -5.98333333e-001,
# -4.98333333e-001, -4.81666667e-001, -4.98333333e-001, -4.81666667e-001,
# -4.98333333e-001, -4.81666667e-001, -4.98333333e-001, -4.81666667e-001,
# -4.98333333e-001, -4.81666667e-001, -4.98333333e-001, -4.81666667e-001,
# -5.48333333e-001, -5.31666667e-001, -5.48333333e-001, -5.31666667e-001,
# -5.48333333e-001, -5.31666667e-001, -5.48333333e-001, -5.31666667e-001,
# -5.48333333e-001, -5.31666667e-001, -5.48333333e-001, -5.31666667e-001])

        os.remove(domain.get_name() + '.sww')

    def test_spatio_temporal_boundary_1(self):
        """Test that boundary values can be read from file and interpolated
        in both time and space.

        Verify that the same steady state solution is arrived at and that
        time interpolation works.

        The full solution history is not exactly the same as
        file boundary must read and interpolate from *smoothed* version
        as stored in sww.
        """
        import time

        #Create sww file of simple propagation from left to right
        #through rectangular domain

        from mesh_factory import rectangular

        #Create basic mesh
        points, vertices, boundary = rectangular(3, 3)

        #Create shallow water domain
        domain1 = Domain(points, vertices, boundary)

        domain1.reduction = mean
        domain1.smooth = False  #Exact result

        domain1.default_order = 2
        domain1.store = True
        domain1.set_datadir('.')
        domain1.set_name('spatio_temporal_boundary_source' + str(time.time()))

        #FIXME: This is extremely important!
        #How can we test if they weren't stored?
        domain1.quantities_to_be_stored = ['stage', 'xmomentum', 'ymomentum']


        #Bed-slope and friction at vertices (and interpolated elsewhere)
        domain1.set_quantity('elevation', 0)
        domain1.set_quantity('friction', 0)

        # Boundary conditions
        Br = Reflective_boundary(domain1)
        Bd = Dirichlet_boundary([0.3,0,0])
        domain1.set_boundary({'left': Bd, 'top': Bd, 'right': Br, 'bottom': Br})
        #Initial condition
        domain1.set_quantity('stage', 0)
        domain1.check_integrity()

        finaltime = 5
        #Evolution  (full domain - large steps)
        for t in domain1.evolve(yieldstep = 0.671, finaltime = finaltime):
            pass
            #domain1.write_time()

        cv1 = domain1.quantities['stage'].centroid_values


        #Create a triangle shaped domain (reusing coordinates from domain 1),
        #formed from the lower and right hand  boundaries and
        #the sw-ne diagonal
        #from domain 1. Call it domain2

        points = [ [0,0], [1.0/3,0], [1.0/3,1.0/3],
                   [2.0/3,0], [2.0/3,1.0/3], [2.0/3,2.0/3],
                   [1,0], [1,1.0/3], [1,2.0/3], [1,1]]

        vertices = [ [1,2,0], [3,4,1], [2,1,4], [4,5,2],
                     [6,7,3], [4,3,7], [7,8,4], [5,4,8], [8,9,5]]

        boundary = { (0,1):'bottom', (1,1):'bottom', (4,1): 'bottom',
                     (4,2):'right', (6,2):'right', (8,2):'right',
                     (0,0):'diagonal', (3,0):'diagonal', (8,0):'diagonal'}

        domain2 = Domain(points, vertices, boundary)

        domain2.reduction = domain1.reduction
        domain2.smooth = False
        domain2.default_order = 2

        #Bed-slope and friction at vertices (and interpolated elsewhere)
        domain2.set_quantity('elevation', 0)
        domain2.set_quantity('friction', 0)
        domain2.set_quantity('stage', 0)

        # Boundary conditions
        Br = Reflective_boundary(domain2)
        #Bf = Spatio_temporal_boundary(domain1.get_name() + '.' +\
        #                              domain1.format, domain2)
        Bf = Field_boundary(domain1.get_name() + '.' +\
                            domain1.format, domain2)        
        domain2.set_boundary({'right':Br, 'bottom':Br, 'diagonal':Bf})
        domain2.check_integrity()



        #Evolution (small steps)
        for t in domain2.evolve(yieldstep = 0.0711, finaltime = finaltime):
            pass


        #Use output from domain1 as spatio-temporal boundary for domain2
        #and verify that results at right hand side are close.

        cv2 = domain2.quantities['stage'].centroid_values

        #print take(cv1, (12,14,16))  #Right
        #print take(cv2, (4,6,8))
        #print take(cv1, (0,6,12))  #Bottom
        #print take(cv2, (0,1,4))
        #print take(cv1, (0,8,16))  #Diag
        #print take(cv2, (0,3,8))

        assert allclose( take(cv1, (0,8,16)), take(cv2, (0,3,8))) #Diag
        assert allclose( take(cv1, (0,6,12)), take(cv2, (0,1,4))) #Bottom
        assert allclose( take(cv1, (12,14,16)), take(cv2, (4,6,8))) #RHS

        #Cleanup
        os.remove(domain1.get_name() + '.' + domain1.format)
        os.remove(domain2.get_name() + '.' + domain2.format)        



    def test_spatio_temporal_boundary_2(self):
        """Test that boundary values can be read from file and interpolated
        in both time and space.
        This is a more basic test, verifying that boundary object
        produces the expected results


        """
        import time

        #Create sww file of simple propagation from left to right
        #through rectangular domain

        from mesh_factory import rectangular

        #Create basic mesh
        points, vertices, boundary = rectangular(3, 3)

        #Create shallow water domain
        domain1 = Domain(points, vertices, boundary)

        domain1.reduction = mean
        domain1.smooth = True #To mimic MOST output

        domain1.default_order = 2
        domain1.store = True
        domain1.set_datadir('.')
        domain1.set_name('spatio_temporal_boundary_source' + str(time.time()))

        #FIXME: This is extremely important!
        #How can we test if they weren't stored?
        domain1.quantities_to_be_stored = ['stage', 'xmomentum', 'ymomentum']


        #Bed-slope and friction at vertices (and interpolated elsewhere)
        domain1.set_quantity('elevation', 0)
        domain1.set_quantity('friction', 0)

        # Boundary conditions
        Br = Reflective_boundary(domain1)
        Bd = Dirichlet_boundary([0.3,0,0])
        domain1.set_boundary({'left': Bd, 'top': Bd, 'right': Br, 'bottom': Br})
        #Initial condition
        domain1.set_quantity('stage', 0)
        domain1.check_integrity()

        finaltime = 5
        #Evolution  (full domain - large steps)
        for t in domain1.evolve(yieldstep = 1, finaltime = finaltime):
            pass
            #domain1.write_time()


        #Create an triangle shaped domain (coinciding with some
        #coordinates from domain 1),
        #formed from the lower and right hand  boundaries and
        #the sw-ne diagonal
        #from domain 1. Call it domain2

        points = [ [0,0], [1.0/3,0], [1.0/3,1.0/3],
                   [2.0/3,0], [2.0/3,1.0/3], [2.0/3,2.0/3],
                   [1,0], [1,1.0/3], [1,2.0/3], [1,1]]

        vertices = [ [1,2,0],
                     [3,4,1], [2,1,4], [4,5,2],
                     [6,7,3], [4,3,7], [7,8,4], [5,4,8], [8,9,5]]

        boundary = { (0,1):'bottom', (1,1):'bottom', (4,1): 'bottom',
                     (4,2):'right', (6,2):'right', (8,2):'right',
                     (0,0):'diagonal', (3,0):'diagonal', (8,0):'diagonal'}

        domain2 = Domain(points, vertices, boundary)

        domain2.reduction = domain1.reduction
        domain2.smooth = False
        domain2.default_order = 2

        #Bed-slope and friction at vertices (and interpolated elsewhere)
        domain2.set_quantity('elevation', 0)
        domain2.set_quantity('friction', 0)
        domain2.set_quantity('stage', 0)


        #Read results for specific timesteps t=1 and t=2
        from Scientific.IO.NetCDF import NetCDFFile
        fid = NetCDFFile(domain1.get_name() + '.' + domain1.format)

        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        s1 = fid.variables['stage'][1,:]
        s2 = fid.variables['stage'][2,:]
        fid.close()

        from Numeric import take, reshape, concatenate
        shp = (len(x), 1)
        points = concatenate( (reshape(x, shp), reshape(y, shp)), axis=1)
        #The diagonal points of domain 1 are 0, 5, 10, 15

        #print points[0], points[5], points[10], points[15]
        assert allclose( take(points, [0,5,10,15]),
                         [[0,0], [1.0/3, 1.0/3], [2.0/3, 2.0/3], [1,1]])


        # Boundary conditions
        Br = Reflective_boundary(domain2)
        #Bf = Spatio_temporal_boundary(domain1.get_name() + '.' + domain1.format,
        #                              domain2)
        Bf = Field_boundary(domain1.get_name() + '.' + domain1.format,
                            domain2, verbose=False)        
        domain2.set_boundary({'right':Br, 'bottom':Br, 'diagonal':Bf})
        domain2.check_integrity()

        #Test that interpolation points are the mid points of the all boundary
        #segments

        boundary_midpoints = [[1.0/6, 0], [1.0/2, 0], [5.0/6,0],
                              [1.0, 1.0/6], [1.0, 1.0/2], [1.0, 5.0/6],
                              [1.0/6, 1.0/6], [0.5, 0.5], [5.0/6, 5.0/6]]

        boundary_midpoints.sort()
        R = Bf.F.interpolation_points.tolist()
        R.sort()
        assert allclose(boundary_midpoints, R)

        #Check spatially interpolated output at time == 1
        domain2.time = 1

        #First diagonal midpoint
        R0 = Bf.evaluate(0,0)
        assert allclose(R0[0], (s1[0] + s1[5])/2)

        #Second diagonal midpoint
        R0 = Bf.evaluate(3,0)
        assert allclose(R0[0], (s1[5] + s1[10])/2)

        #First diagonal midpoint
        R0 = Bf.evaluate(8,0)
        assert allclose(R0[0], (s1[10] + s1[15])/2)

        #Check spatially interpolated output at time == 2
        domain2.time = 2

        #First diagonal midpoint
        R0 = Bf.evaluate(0,0)
        assert allclose(R0[0], (s2[0] + s2[5])/2)

        #Second diagonal midpoint
        R0 = Bf.evaluate(3,0)
        assert allclose(R0[0], (s2[5] + s2[10])/2)

        #First diagonal midpoint
        R0 = Bf.evaluate(8,0)
        assert allclose(R0[0], (s2[10] + s2[15])/2)


        #Now check temporal interpolation

        domain2.time = 1 + 2.0/3

        #First diagonal midpoint
        R0 = Bf.evaluate(0,0)
        assert allclose(R0[0], ((s1[0] + s1[5])/2 + 2.0*(s2[0] + s2[5])/2)/3)

        #Second diagonal midpoint
        R0 = Bf.evaluate(3,0)
        assert allclose(R0[0], ((s1[5] + s1[10])/2 + 2.0*(s2[5] + s2[10])/2)/3)

        #First diagonal midpoint
        R0 = Bf.evaluate(8,0)
        assert allclose(R0[0], ((s1[10] + s1[15])/2 + 2.0*(s2[10] + s2[15])/2)/3)



        #Cleanup
        os.remove(domain1.get_name() + '.' + domain1.format)


    def test_spatio_temporal_boundary_3(self):
        """Test that boundary values can be read from file and interpolated
        in both time and space.
        This is a more basic test, verifying that boundary object
        produces the expected results

        This tests adjusting using mean_stage

        """

        import time

        mean_stage = 5.2 # Adjust stage by this amount in boundary

        #Create sww file of simple propagation from left to right
        #through rectangular domain

        from mesh_factory import rectangular

        #Create basic mesh
        points, vertices, boundary = rectangular(3, 3)

        #Create shallow water domain
        domain1 = Domain(points, vertices, boundary)

        domain1.reduction = mean
        domain1.smooth = True #To mimic MOST output

        domain1.default_order = 2
        domain1.store = True
        domain1.set_datadir('.')
        domain1.set_name('spatio_temporal_boundary_source' + str(time.time()))

        #FIXME: This is extremely important!
        #How can we test if they weren't stored?
        domain1.quantities_to_be_stored = ['stage', 'xmomentum', 'ymomentum']


        #Bed-slope and friction at vertices (and interpolated elsewhere)
        domain1.set_quantity('elevation', 0)
        domain1.set_quantity('friction', 0)

        # Boundary conditions
        Br = Reflective_boundary(domain1)
        Bd = Dirichlet_boundary([0.3,0,0])
        domain1.set_boundary({'left': Bd, 'top': Bd, 'right': Br, 'bottom': Br})
        #Initial condition
        domain1.set_quantity('stage', 0)
        domain1.check_integrity()

        finaltime = 5
        #Evolution  (full domain - large steps)
        for t in domain1.evolve(yieldstep = 1, finaltime = finaltime):
            pass
            #domain1.write_time()


        #Create an triangle shaped domain (coinciding with some
        #coordinates from domain 1),
        #formed from the lower and right hand  boundaries and
        #the sw-ne diagonal
        #from domain 1. Call it domain2

        points = [ [0,0],
                   [1.0/3,0], [1.0/3,1.0/3],
                   [2.0/3,0], [2.0/3,1.0/3], [2.0/3,2.0/3],
                   [1,0],     [1,1.0/3],     [1,2.0/3],     [1,1]]  
                   
        vertices = [ [1,2,0],
                     [3,4,1], [2,1,4], [4,5,2],
                     [6,7,3], [4,3,7], [7,8,4], [5,4,8], [8,9,5]]

        boundary = { (0,1):'bottom', (1,1):'bottom', (4,1): 'bottom',
                     (4,2):'right', (6,2):'right', (8,2):'right',
                     (0,0):'diagonal', (3,0):'diagonal', (8,0):'diagonal'}

        domain2 = Domain(points, vertices, boundary)

        domain2.reduction = domain1.reduction
        domain2.smooth = False
        domain2.default_order = 2

        #Bed-slope and friction at vertices (and interpolated elsewhere)
        domain2.set_quantity('elevation', 0)
        domain2.set_quantity('friction', 0)
        domain2.set_quantity('stage', 0)


        #Read results for specific timesteps t=1 and t=2
        from Scientific.IO.NetCDF import NetCDFFile
        fid = NetCDFFile(domain1.get_name() + '.' + domain1.format)

        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        s1 = fid.variables['stage'][1,:]
        s2 = fid.variables['stage'][2,:]
        fid.close()

        from Numeric import take, reshape, concatenate
        shp = (len(x), 1)
        points = concatenate( (reshape(x, shp), reshape(y, shp)), axis=1)
        #The diagonal points of domain 1 are 0, 5, 10, 15

        #print points[0], points[5], points[10], points[15]
        assert allclose( take(points, [0,5,10,15]),
                         [[0,0], [1.0/3, 1.0/3], [2.0/3, 2.0/3], [1,1]])


        # Boundary conditions
        Br = Reflective_boundary(domain2)
        #Bf = Spatio_temporal_boundary(domain1.get_name() + '.' + domain1.format,
        #                              domain2)
        Bf = Field_boundary(domain1.get_name() + '.' + domain1.format,
                            domain2, mean_stage=mean_stage, verbose=False)
        
        domain2.set_boundary({'right':Br, 'bottom':Br, 'diagonal':Bf})
        domain2.check_integrity()

        #Test that interpolation points are the mid points of the all boundary
        #segments

        boundary_midpoints = [[1.0/6, 0], [1.0/2, 0], [5.0/6,0],
                              [1.0, 1.0/6], [1.0, 1.0/2], [1.0, 5.0/6],
                              [1.0/6, 1.0/6], [0.5, 0.5], [5.0/6, 5.0/6]]

        boundary_midpoints.sort()
        R = Bf.F.interpolation_points.tolist()
        R.sort()
        assert allclose(boundary_midpoints, R)

        #Check spatially interpolated output at time == 1
        domain2.time = 1

        #First diagonal midpoint
        R0 = Bf.evaluate(0,0)
        assert allclose(R0[0], (s1[0] + s1[5])/2 + mean_stage)

        #Second diagonal midpoint
        R0 = Bf.evaluate(3,0)
        assert allclose(R0[0], (s1[5] + s1[10])/2 + mean_stage)

        #First diagonal midpoint
        R0 = Bf.evaluate(8,0)
        assert allclose(R0[0], (s1[10] + s1[15])/2 + mean_stage)

        #Check spatially interpolated output at time == 2
        domain2.time = 2

        #First diagonal midpoint
        R0 = Bf.evaluate(0,0)
        assert allclose(R0[0], (s2[0] + s2[5])/2 + mean_stage)

        #Second diagonal midpoint
        R0 = Bf.evaluate(3,0)
        assert allclose(R0[0], (s2[5] + s2[10])/2 + mean_stage)

        #First diagonal midpoint
        R0 = Bf.evaluate(8,0)
        assert allclose(R0[0], (s2[10] + s2[15])/2 + mean_stage)


        #Now check temporal interpolation

        domain2.time = 1 + 2.0/3

        #First diagonal midpoint
        R0 = Bf.evaluate(0,0)
        assert allclose(R0[0], ((s1[0] + s1[5])/2 + 2.0*(s2[0] + s2[5])/2)/3 + mean_stage)

        #Second diagonal midpoint
        R0 = Bf.evaluate(3,0)
        assert allclose(R0[0], ((s1[5] + s1[10])/2 + 2.0*(s2[5] + s2[10])/2)/3 + mean_stage)

        #First diagonal midpoint
        R0 = Bf.evaluate(8,0)
        assert allclose(R0[0], ((s1[10] + s1[15])/2 + 2.0*(s2[10] + s2[15])/2)/3 + mean_stage)


        #Cleanup
        os.remove(domain1.get_name() + '.' + domain1.format)


    def test_spatio_temporal_boundary_outside(self):
        """Test that field_boundary catches if a point is outside the sww that defines it
        """

        import time
        #Create sww file of simple propagation from left to right
        #through rectangular domain

        from mesh_factory import rectangular

        #Create basic mesh
        points, vertices, boundary = rectangular(3, 3)

        #Create shallow water domain
        domain1 = Domain(points, vertices, boundary)

        domain1.reduction = mean
        domain1.smooth = True #To mimic MOST output

        domain1.default_order = 2
        domain1.store = True
        domain1.set_datadir('.')
        domain1.set_name('spatio_temporal_boundary_source' + str(time.time()))

        #FIXME: This is extremely important!
        #How can we test if they weren't stored?
        domain1.quantities_to_be_stored = ['stage', 'xmomentum', 'ymomentum']


        #Bed-slope and friction at vertices (and interpolated elsewhere)
        domain1.set_quantity('elevation', 0)
        domain1.set_quantity('friction', 0)

        # Boundary conditions
        Br = Reflective_boundary(domain1)
        Bd = Dirichlet_boundary([0.3,0,0])
        domain1.set_boundary({'left': Bd, 'top': Bd, 'right': Br, 'bottom': Br})
        #Initial condition
        domain1.set_quantity('stage', 0)
        domain1.check_integrity()

        finaltime = 5
        #Evolution  (full domain - large steps)
        for t in domain1.evolve(yieldstep = 1, finaltime = finaltime):
            pass
            #domain1.write_time()


        #Create an triangle shaped domain (coinciding with some
        #coordinates from domain 1, but one edge outside!),
        #formed from the lower and right hand  boundaries and
        #the sw-ne diagonal as in the previous test but scaled
        #in the x direction by a factor of 2

        points = [ [0,0],
                   [2.0/3,0], [2.0/3,1.0/3],
                   [4.0/3,0], [4.0/3,1.0/3], [4.0/3,2.0/3],
                   [2,0],     [2,1.0/3],     [2,2.0/3],     [2,1]  
                   ]

        vertices = [ [1,2,0],
                     [3,4,1], [2,1,4], [4,5,2],
                     [6,7,3], [4,3,7], [7,8,4], [5,4,8], [8,9,5]]

        boundary = { (0,1):'bottom', (1,1):'bottom', (4,1): 'bottom',
                     (4,2):'right', (6,2):'right', (8,2):'right',
                     (0,0):'diagonal', (3,0):'diagonal', (8,0):'diagonal'}

        domain2 = Domain(points, vertices, boundary)

        domain2.reduction = domain1.reduction
        domain2.smooth = False
        domain2.default_order = 2

        #Bed-slope and friction at vertices (and interpolated elsewhere)
        domain2.set_quantity('elevation', 0)
        domain2.set_quantity('friction', 0)
        domain2.set_quantity('stage', 0)


        #Read results for specific timesteps t=1 and t=2
        from Scientific.IO.NetCDF import NetCDFFile
        fid = NetCDFFile(domain1.get_name() + '.' + domain1.format)

        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        s1 = fid.variables['stage'][1,:]
        s2 = fid.variables['stage'][2,:]
        fid.close()

        from Numeric import take, reshape, concatenate
        shp = (len(x), 1)
        points = concatenate( (reshape(x, shp), reshape(y, shp)), axis=1)
        #The diagonal points of domain 1 are 0, 5, 10, 15

        #print points[0], points[5], points[10], points[15]
        assert allclose( take(points, [0,5,10,15]),
                         [[0,0], [1.0/3, 1.0/3], [2.0/3, 2.0/3], [1,1]])


        # Boundary conditions
        Br = Reflective_boundary(domain2)
        #Bf = Spatio_temporal_boundary(domain1.get_name() + '.' + domain1.format,
        #                              domain2)
        Bf = Field_boundary(domain1.get_name() + '.' + domain1.format,
                            domain2, mean_stage=1, verbose=False)
        
        domain2.set_boundary({'right':Br, 'bottom':Br, 'diagonal':Bf})
        domain2.check_integrity()

        try:
            for t in domain2.evolve(yieldstep = 1, finaltime = finaltime):
                pass
        except:
            pass
        else:
            msg = 'This should have caught NAN at boundary'
            raise Exception, msg


        #Cleanup
        os.remove(domain1.get_name() + '.' + domain1.format)




    def test_extrema(self):
        """Test that extrema of quantities are computed correctly
        Extrema are updated at every *internal* timestep
        """

        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross

        initial_runup_height = -0.4
        final_runup_height = -0.3


        #--------------------------------------------------------------
        # Setup computational domain
        #--------------------------------------------------------------
        N = 5
        points, vertices, boundary = rectangular_cross(N, N) 
        domain = Domain(points, vertices, boundary)
        domain.set_name('extrema_test')

        #--------------------------------------------------------------
        # Setup initial conditions
        #--------------------------------------------------------------
        def topography(x,y):
            return -x/2                             # linear bed slope
            

        domain.set_quantity('elevation', topography)       # Use function for elevation
        domain.set_quantity('friction', 0.)                # Zero friction 
        domain.set_quantity('stage', initial_runup_height) # Constant negative initial stage
        domain.set_quantities_to_be_monitored(['stage', 'stage-elevation'],
                                              time_interval = [0.5, 2.7],
                                              polygon = [[0,0], [0,1], [1,1], [1,0]])
        
        assert len(domain.quantities_to_be_monitored) == 2
        assert domain.quantities_to_be_monitored.has_key('stage')
        assert domain.quantities_to_be_monitored.has_key('stage-elevation')
        for key in domain.quantities_to_be_monitored['stage'].keys():
            assert domain.quantities_to_be_monitored['stage'][key] is None        


        #--------------------------------------------------------------
        # Setup boundary conditions
        #--------------------------------------------------------------
        Br = Reflective_boundary(domain)              # Reflective wall
        Bd = Dirichlet_boundary([final_runup_height,  # Constant inflow
                                 0,
                                 0])

        # All reflective to begin with (still water) 
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})


        #--------------------------------------------------------------
        # Let triangles adjust and check extrema 
        #--------------------------------------------------------------
        for t in domain.evolve(yieldstep = 0.1, finaltime = 1.0):
            domain.quantity_statistics() # Run it silently



        #--------------------------------------------------------------
        # Test extrema
        #--------------------------------------------------------------

        stage = domain.quantities_to_be_monitored['stage']
        assert stage['min'] <= stage['max']

        #print stage['min'], stage['max'] 
        assert allclose(stage['min'], initial_runup_height,
                        rtol = 1.0/N) # First order accuracy


        depth = domain.quantities_to_be_monitored['stage-elevation']
        assert depth['min'] <= depth['max'] 
        assert depth['min'] >= 0.0
        assert depth['max'] >= 0.0        
        ##assert depth[1] <= ?? initial_runup_height        


        #--------------------------------------------------------------
        # Update boundary to allow inflow
        #--------------------------------------------------------------
        domain.set_boundary({'right': Bd})

        
        #--------------------------------------------------------------
        # Evolve system through time
        #--------------------------------------------------------------
        for t in domain.evolve(yieldstep = 0.1, finaltime = 3.0):
            #domain.write_time()
            domain.quantity_statistics() # Run it silently           
            
    
        #--------------------------------------------------------------
        # Test extrema again
        #--------------------------------------------------------------

        stage = domain.quantities_to_be_monitored['stage']
        assert stage['min'] <= stage['max']

        assert allclose(stage['min'], initial_runup_height,
                        rtol = 1.0/N) # First order accuracy        

        depth = domain.quantities_to_be_monitored['stage-elevation']
        assert depth['min'] <= depth['max'] 
        assert depth['min'] >= 0.0
        assert depth['max'] >= 0.0        

        #Cleanup
        os.remove(domain.get_name() + '.' + domain.format)
        


    def test_tight_slope_limiters(self):
        """Test that new slope limiters (Feb 2007) don't induce extremely
        small timesteps. This test actually reveals the problem as it
        was in March-April 2007 
        """

        import time, os
        from Numeric import array, zeros, allclose, Float, concatenate
        from Scientific.IO.NetCDF import NetCDFFile
        from data_manager import get_dataobject, extent_sww
        from mesh_factory import rectangular

        
        #Create basic mesh
        points, vertices, boundary = rectangular(2, 2)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.default_order = 2

        # This will pass
        #domain.tight_slope_limiters = 1
        #domain.H0 = 0.01
        
        # This will fail
        #domain.tight_slope_limiters = 1
        #domain.H0 = 0.001

        # This will pass provided C extension implements limiting of
        # momentum in _compute_speeds
        domain.tight_slope_limiters = 1
        domain.H0 = 0.001        
        domain.protect_against_isolated_degenerate_timesteps = True        

        #Set some field values
        domain.set_quantity('elevation', lambda x,y: -x)
        domain.set_quantity('friction', 0.03)


        ######################
        # Boundary conditions
        B = Transmissive_boundary(domain)
        domain.set_boundary( {'left': B, 'right': B, 'top': B, 'bottom': B})


        ######################
        #Initial condition - with jumps


        bed = domain.quantities['elevation'].vertex_values
        stage = zeros(bed.shape, Float)

        h = 0.3
        for i in range(stage.shape[0]):
            if i % 2 == 0:
                stage[i,:] = bed[i,:] + h
            else:
                stage[i,:] = bed[i,:]

        domain.set_quantity('stage', stage)


        domain.distribute_to_vertices_and_edges()               

        

        domain.set_name('tight_limiters')
        domain.format = 'sww'
        domain.smooth = True
        domain.reduction = mean
        domain.set_datadir('.')
        domain.smooth = False
        domain.store = True
        domain.beta_h = 0
        

        #Evolution
        for t in domain.evolve(yieldstep = 0.1, finaltime = 0.3):
            
            #domain.write_time(track_speeds=True)
            stage = domain.quantities['stage'].vertex_values

            #Get NetCDF
            fid = NetCDFFile(domain.writer.filename, 'r')
            stage_file = fid.variables['stage']
            
            fid.close()

        os.remove(domain.writer.filename)


    def test_pmesh2Domain(self):
         import os
         import tempfile

         fileName = tempfile.mktemp(".tsh")
         file = open(fileName,"w")
         file.write("4 3 # <vertex #> <x> <y> [attributes]\n \
0 0.0 0.0 0.0 0.0 0.01 \n \
1 1.0 0.0 10.0 10.0 0.02  \n \
2 0.0 1.0 0.0 10.0 0.03  \n \
3 0.5 0.25 8.0 12.0 0.04  \n \
# Vert att title  \n \
elevation  \n \
stage  \n \
friction  \n \
2 # <triangle #> [<vertex #>] [<neigbouring triangle #>]  \n\
0 0 3 2 -1  -1  1 dsg\n\
1 0 1 3 -1  0 -1   ole nielsen\n\
4 # <segment #> <vertex #>  <vertex #> [boundary tag] \n\
0 1 0 2 \n\
1 0 2 3 \n\
2 2 3 \n\
3 3 1 1 \n\
3 0 # <x> <y> [attributes] ...Mesh Vertices... \n \
0 216.0 -86.0 \n \
1 160.0 -167.0 \n \
2 114.0 -91.0 \n \
3 # <vertex #>  <vertex #> [boundary tag] ...Mesh Segments... \n \
0 0 1 0 \n \
1 1 2 0 \n \
2 2 0 0 \n \
0 # <x> <y> ...Mesh Holes... \n \
0 # <x> <y> <attribute>...Mesh Regions... \n \
0 # <x> <y> <attribute>...Mesh Regions, area... \n\
#Geo reference \n \
56 \n \
140 \n \
120 \n")
         file.close()

         tags = {}
         b1 =  Dirichlet_boundary(conserved_quantities = array([0.0]))
         b2 =  Dirichlet_boundary(conserved_quantities = array([1.0]))
         b3 =  Dirichlet_boundary(conserved_quantities = array([2.0]))
         tags["1"] = b1
         tags["2"] = b2
         tags["3"] = b3

         #from anuga.abstract_2d_finite_volumes.pmesh2domain import pmesh_to_domain_instance
         #domain = pmesh_to_domain_instance(fileName, Domain)

         domain = Domain(mesh_filename=fileName)
                         #verbose=True, use_cache=True)
         
         #print "domain.tagged_elements", domain.tagged_elements
         ## check the quantities
         #print domain.quantities['elevation'].vertex_values
         answer = [[0., 8., 0.],
                   [0., 10., 8.]]
         assert allclose(domain.quantities['elevation'].vertex_values,
                        answer)

         #print domain.quantities['stage'].vertex_values
         answer = [[0., 12., 10.],
                   [0., 10., 12.]]
         assert allclose(domain.quantities['stage'].vertex_values,
                        answer)

         #print domain.quantities['friction'].vertex_values
         answer = [[0.01, 0.04, 0.03],
                   [0.01, 0.02, 0.04]]
         assert allclose(domain.quantities['friction'].vertex_values,
                        answer)

         #print domain.quantities['friction'].vertex_values
         assert allclose(domain.tagged_elements['dsg'][0],0)
         assert allclose(domain.tagged_elements['ole nielsen'][0],1)

         self.failUnless( domain.boundary[(1, 0)]  == '1',
                          "test_tags_to_boundaries  failed. Single boundary wasn't added.")
         self.failUnless( domain.boundary[(1, 2)]  == '2',
                          "test_tags_to_boundaries  failed. Single boundary wasn't added.")
         self.failUnless( domain.boundary[(0, 1)]  == '3',
                          "test_tags_to_boundaries  failed. Single boundary wasn't added.")
         self.failUnless( domain.boundary[(0, 0)]  == 'exterior',
                          "test_tags_to_boundaries  failed. Single boundary wasn't added.")
         #print "domain.boundary",domain.boundary
         self.failUnless( len(domain.boundary)  == 4,
                          "test_pmesh2Domain Too many boundaries")
         #FIXME change to use get_xllcorner
         #print "d.geo_reference.xllcorner",domain.geo_reference.xllcorner 
         self.failUnless(domain.geo_reference.xllcorner  == 140.0,
                          "bad geo_referece")


         #************

    
         domain = Domain(fileName)
         
         #print "domain.tagged_elements", domain.tagged_elements
         ## check the quantities
         #print domain.quantities['elevation'].vertex_values
         answer = [[0., 8., 0.],
                   [0., 10., 8.]]
         assert allclose(domain.quantities['elevation'].vertex_values,
                        answer)

         #print domain.quantities['stage'].vertex_values
         answer = [[0., 12., 10.],
                   [0., 10., 12.]]
         assert allclose(domain.quantities['stage'].vertex_values,
                        answer)

         #print domain.quantities['friction'].vertex_values
         answer = [[0.01, 0.04, 0.03],
                   [0.01, 0.02, 0.04]]
         assert allclose(domain.quantities['friction'].vertex_values,
                        answer)

         #print domain.quantities['friction'].vertex_values
         assert allclose(domain.tagged_elements['dsg'][0],0)
         assert allclose(domain.tagged_elements['ole nielsen'][0],1)

         self.failUnless( domain.boundary[(1, 0)]  == '1',
                          "test_tags_to_boundaries  failed. Single boundary wasn't added.")
         self.failUnless( domain.boundary[(1, 2)]  == '2',
                          "test_tags_to_boundaries  failed. Single boundary wasn't added.")
         self.failUnless( domain.boundary[(0, 1)]  == '3',
                          "test_tags_to_boundaries  failed. Single boundary wasn't added.")
         self.failUnless( domain.boundary[(0, 0)]  == 'exterior',
                          "test_tags_to_boundaries  failed. Single boundary wasn't added.")
         #print "domain.boundary",domain.boundary
         self.failUnless( len(domain.boundary)  == 4,
                          "test_pmesh2Domain Too many boundaries")
         #FIXME change to use get_xllcorner
         #print "d.geo_reference.xllcorner",domain.geo_reference.xllcorner 
         self.failUnless(domain.geo_reference.xllcorner  == 140.0,
                          "bad geo_referece")
         #************
         os.remove(fileName)

        #-------------------------------------------------------------

    def test_get_lone_vertices(self):
        
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4] ]
        boundary = { (0, 0): 'Third',
                     (0, 2): 'First',
                     (2, 0): 'Second',
                     (2, 1): 'Second',
                     (3, 1): 'Second',
                     (3, 2): 'Third'}


        domain = Domain(points, vertices, boundary)
        #domain.check_integrity()
        domain.get_lone_vertices()

        
    def test_fitting_using_shallow_water_domain(self):
        
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
        mesh4 = Domain(points, elements,
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
        quantity.set_values(filename = ptsfile,
                            attribute_name = att, alpha = 0)
        answer = linear_function(quantity.domain.get_vertex_coordinates())

        assert allclose(quantity.vertex_values.flat, answer)


        #Check that values can be set from file using default attribute
        quantity.set_values(filename = ptsfile, alpha = 0)
        assert allclose(quantity.vertex_values.flat, answer)

        #Cleanup
        import os
        os.remove(ptsfile)


        
if __name__ == "__main__":

    suite = unittest.makeSuite(Test_Shallow_Water,'test')

    #suite = unittest.makeSuite(Test_Shallow_Water,'test_bedslope_problem_first_order_moresteps')    
    #suite = unittest.makeSuite(Test_Shallow_Water,'test_fitting_using_shallow_water_domain')    
    #suite = unittest.makeSuite(Test_Shallow_Water,'test_tight_slope_limiters')
    #suite = unittest.makeSuite(Test_Shallow_Water,'test_get_maximum_inundation_from_sww')
    #suite = unittest.makeSuite(Test_Shallow_Water,'test_temp')    
    

    
    runner = unittest.TextTestRunner(verbosity=1)    
    runner.run(suite)

import operator
from anuga import Domain
from anuga import Quantity
from anuga import Dirichlet_boundary
from anuga.operators.kinematic_viscosity_operator import Kinematic_viscosity_operator

import numpy as num
from math import sqrt
import unittest
import os

class Test_kinematic_viscosity(unittest.TestCase):
    
    def setUp(self):
        pass

    def tearDown(self):
        try:
            os.remove('domain.sww')
        except:
            pass

        try:
            os.remove('anuga.log')
        except:
            pass
        
    #First test operator class (1 triangle)
    def operator1(self):
        points = num.array([[0.0,0.0],[1.0,0.0],[0.0,1.0]])
        elements = num.array([[0,1,2]])
        boundary_map = {}
        boundary_map[(0,0)] = 'edge0'
        boundary_map[(0,1)] = 'edge1'
        boundary_map[(0,2)] = 'edge2'
        domain = Domain(coordinates=points,vertices=elements,boundary=boundary_map)

        D0 = Dirichlet_boundary([1,0,3])
        D1 = Dirichlet_boundary([2,1,0])
        D2 = Dirichlet_boundary([3,1,2])
        domain.set_boundary({'edge0': D0, 'edge1': D1, 'edge2': D2})

        domain.set_quantity('stage', lambda x,y : x+2*y )
        domain.set_quantity('elevation', lambda x,y : 3*x+5*y )



        #print domain.quantities['stage'].vertex_values

        #print domain.quantities['stage'].edge_values

        domain.update_boundary()


        #print domain.quantities['stage'].boundary_values
        
        return Kinematic_viscosity_operator(domain)

    #Second test operator class (2 triangles)
    def operator2(self):
        points = num.array([[0.0,0.0],[1.0,0.0],[1.0,1.0],[0.0,1.0]])
        elements = num.array([[0,1,3],[1,2,3]])
        boundary_map = {}
        boundary_map[(0,1)] = 'edge0'
        boundary_map[(0,2)] = 'edge1'
        boundary_map[(1,0)] = 'edge2'
        boundary_map[(1,2)] = 'edge3'
        domain = Domain(coordinates=points,vertices=elements,boundary=boundary_map)

        D0 = Dirichlet_boundary([1,1,2])
        D1 = Dirichlet_boundary([1,2,2])
        D2 = Dirichlet_boundary([1,1,0])
        D3 = Dirichlet_boundary([1,2,1])

        domain.set_boundary({'edge0': D0, 'edge1': D1, 'edge2': D2, 'edge3': D3})
        domain.update_boundary()



        return Kinematic_viscosity_operator(domain)

    def test_enumerate_boundary(self):
        operator1 = self.operator1()
        boundary_enumeration = operator1.domain.boundary_enumeration

        assert boundary_enumeration[(0,0)] == 0
        assert boundary_enumeration[(0,1)] == 1
        assert boundary_enumeration[(0,2)] == 2

        operator2 = self.operator2()
        boundary_enumeration = operator2.domain.boundary_enumeration


        assert boundary_enumeration[(0,1)] == 0
        assert boundary_enumeration[(0,2)] == 1
        assert boundary_enumeration[(1,0)] == 2
        assert boundary_enumeration[(1,2)] == 3

    def test_geo_structure(self):
        operator1 = self.operator1()
        indices = operator1.geo_structure_indices
        values = operator1.geo_structure_values

        assert num.allclose(indices, num.array([[1, 2, 3]]))
        assert num.allclose(values, num.array([[-6.0, -6.0/sqrt(5), -6.0/sqrt(5)]]))

        operator2 = self.operator2()
        indices = operator2.geo_structure_indices
        values = operator2.geo_structure_values
        assert num.allclose(indices, num.array([[1,2,3],[4,0,5]]))
        assert num.allclose(values, num.array([[-3.0,-6.0/sqrt(5),-6.0/sqrt(5)],[-6.0/sqrt(5),-3.0,-6.0/sqrt(5)]]))

    def test_elliptic_matrix_one_triangle(self):

        operator = self.operator1()
        domain = operator.domain

        a = Quantity(operator.domain)
        a.set_values(1.0)
        a.set_boundary_values(1.0)
        
        operator.update_elliptic_matrix(a)

        A = operator.elliptic_matrix

        assert num.allclose(A.todense(), num.array([-6.0-12.0/sqrt(5), 6.0,  6.0/sqrt(5), 6.0/sqrt(5)]))

        a.set_values(10.0)
        a.set_boundary_values(10.0)
        
        operator.update_elliptic_matrix(a)

        assert num.allclose(A.todense(), 10*num.array([-6.0-12.0/sqrt(5), 6.0,  6.0/sqrt(5), 6.0/sqrt(5)]))
    

    def test_elliptic_matrix_two_triangles(self):


        operator = self.operator2()

        domain = operator.domain

        a = Quantity(operator.domain)
        a.set_values(1.0)
        a.set_boundary_values(1.0)
        operator.update_elliptic_matrix(a)

        A = operator.elliptic_matrix

    
        A0 = num.array([[-3.0,3.0,0.0,0.0,0.0,0.0],
                        [0.0,-6.0/sqrt(5.0),0.0,0.0,6.0/sqrt(5.0),0.0]])
        A1 = num.array([[-6.0/sqrt(5.0),0.0,6.0/sqrt(5.0),0.0,0.0,0.0],\
                        [3.0,-3.0,0.0,0.0,0.0,0.0]])
        A2 = num.array([[-6.0/sqrt(5.0),0.0,0.0,6.0/sqrt(5.0),0.0,0.0],\
                        [0.0, -6.0/sqrt(5.0), 0.0, 0.0, 0.0, 6.0/sqrt(5.0)]])


        assert num.allclose(A.todense(), A0+A1+A2)

        a.set_values([2.0, 1.0], location = 'centroids')
        a.set_boundary_values(1.0)
        operator.update_elliptic_matrix(a)

        A = operator.elliptic_matrix
        

        assert num.allclose(A.todense()[0,:], 1.5*A0[0,:]+1.5*A1[0,:]+1.5*A2[0,:])
        assert num.allclose(A.todense()[1,:], A0[1,:]+1.5*A1[1,:]+A2[1,:])

        # Either negative values we set matrix row to zero
        a.set_values([-2.0, -2.0], location = 'centroids')
        a.set_boundary_values(1.0)
        operator.update_elliptic_matrix(a)



        assert num.allclose(A.todense()[0,:], 0.0)
        assert num.allclose(A.todense()[1,:], 0.0)

    def test_elliptic_multiply_include_boundary_one_triangle(self):
        operator = self.operator1()
        operator.set_triangle_areas(False)

        print operator.apply_triangle_areas

        a = Quantity(operator.domain)
        a.set_values(1.0)
        a.set_boundary_values(1.0)

        operator.update_elliptic_matrix()

        

        q_in = Quantity(operator.domain)
        q_in.set_values(1.0)
        q_in.set_boundary_values(1.0)
        
        n = operator.n
        
        A = num.array([-6.0-12.0/sqrt(5), 6.0,  6.0/sqrt(5), 6.0/sqrt(5)])


        q_1 = operator.elliptic_multiply(q_in)

        q_2 = operator.elliptic_multiply(q_in, quantity_out = q_in)

        assert id(q_in) == id(q_2)

        assert num.allclose(q_1.centroid_values,q_2.centroid_values)

        assert num.allclose( num.zeros((n,), num.float), q_1.centroid_values )

        #Now have different boundary values

        q_in.set_values(1.0)
        q_in.set_boundary_values(0.0)
        operator.update_elliptic_matrix(a)


        A = num.array([-6.0-12.0/sqrt(5), 6.0,  6.0/sqrt(5), 6.0/sqrt(5)])


        q_1 = operator.elliptic_multiply(q_in)

        assert num.allclose( [-6.0-12.0/sqrt(5)], q_1.centroid_values )
        
    def test_elliptic_multiply_exclude_boundary_one_triangle(self):
        operator = self.operator1()
        operator.set_triangle_areas(False)

        print operator.apply_triangle_areas
        #n = operator.n

        q_in = Quantity(operator.domain)
        q_in.set_values(1.0)
        q_in.set_boundary_values(1.0)
        operator.update_elliptic_matrix()


        A = num.array([-6.0-12.0/sqrt(5), 6.0,  6.0/sqrt(5), 6.0/sqrt(5)])


        q_1 = operator.elliptic_multiply(q_in, include_boundary=False)


        assert num.allclose( [-6.0-12.0/sqrt(5)], q_1.centroid_values )

    def test_elliptic_multiply_include_boundary_one_triangle(self):
        operator = self.operator1()
        operator.set_triangle_areas(True)

        n = operator.n

        q_in = Quantity(operator.domain)
        q_in.set_values(1.0)
        q_in.set_boundary_values(1.0)
        
        operator.update_elliptic_matrix()


        A = num.array([-6.0-12.0/sqrt(5), 6.0,  6.0/sqrt(5), 6.0/sqrt(5)])


        q_1 = operator.elliptic_multiply(q_in)

        q_2 = operator.elliptic_multiply(q_in, output = q_in)

        assert id(q_in) == id(q_2)

        assert num.allclose(q_1.centroid_values,q_2.centroid_values)

        assert num.allclose( [-12.0-24.0/sqrt(5)], q_1.centroid_values )

        #Now have different boundary values

        q_in.set_values(1.0)
        q_in.set_boundary_values(0.0)
        operator.update_elliptic_matrix()


        A = num.array([-6.0-12.0/sqrt(5), 6.0,  6.0/sqrt(5), 6.0/sqrt(5)])


        q_1 = operator.elliptic_multiply(q_in)

        assert num.allclose( [-12.0-24.0/sqrt(5)], q_1.centroid_values )

    def test_elliptic_multiply_exclude_boundary_one_triangle(self):
        operator = self.operator1()
        operator.set_triangle_areas(True)

        q_in = Quantity(operator.domain)
        q_in.set_values(1.0)
        q_in.set_boundary_values(1.0)
        operator.update_elliptic_matrix()


        A = num.array([-6.0-12.0/sqrt(5), 6.0,  6.0/sqrt(5), 6.0/sqrt(5)])


        q_1 = operator.elliptic_multiply(q_in)


        assert num.allclose( [-12.0-24.0/sqrt(5)], q_1.centroid_values )

    def test_mul_arg(self):
        operator = self.operator1()

        u = Quantity(operator.domain)
        u.set_values(2.0)
        #q boundary_values should equal 0.0



        operator.update_elliptic_boundary_term(u)

        r = 2.0

        try:
            q_out = operator * 2.0
        except TypeError:
            pass
        else:
            raise Exception('Should have caught an TypeError')


    def test_mul(self):
        operator = self.operator1()

        u = Quantity(operator.domain)
        u.set_values(2.0)
        #q boundary_values should equal 0.0

        operator.update_elliptic_matrix()

        operator.update_elliptic_boundary_term(u)

        A = num.array([-6.0-12.0/sqrt(5), 6.0,  6.0/sqrt(5), 6.0/sqrt(5)])
        V1 = num.array([2.0]) #u=2
        U1 = num.array([[2.0],[0.0],[0.0],[0.0]])

        q_out = operator * u
        
        assert num.allclose(q_out.centroid_values, 2*num.array(num.mat(A)*num.mat(U1)).reshape(1,))

    def test_elliptic_solve_one_triangle(self):

        operator = self.operator1()
        n = operator.n
        
        U = num.array([2.0,2.0,1.0,1.0])

        u_in = Quantity(operator.domain)
        u_in.set_values(U[:1], location='centroids')
        u_in.set_boundary_values(U[1:])

        a = Quantity(operator.domain)
        a.set_values(1.0)
        a.set_boundary_values(1.0)

        # Do this to get access to the matrix
        # This is also called inside elliptic_solve
        operator.update_elliptic_matrix(a)

        V = num.array([2.0]) #h=1, u=2
        A = num.array([-6.0-12.0/sqrt(5), 6.0,  6.0/sqrt(5), 6.0/sqrt(5)])
        #U = num.array([[2.0,2.0],[2.0,1.0],[1.0,2.0],[1.0,0.0]])

        #Setup up rhs as b = A u
        X = num.array(2*num.mat(A)*num.mat(U.reshape(4,1))).reshape(1,)
        b = Quantity(operator.domain)
        b.set_values(X, location='centroids')

        u_in.set_values(0.0)

        u_out = operator.elliptic_solve(u_in, b, a, iprint=1)

        assert num.allclose(u_out.centroid_values, U[:n])


    def test_elliptic_solve_two_triangle(self):

        operator = self.operator2()
        n = operator.n
        
        U = num.array([2.0,3.0,1.0,1.0,4.0,3.0])

        u_in = Quantity(operator.domain)
        u_in.set_values(U[:2], location='centroids')
        u_in.set_boundary_values(U[2:])

        a = Quantity(operator.domain)
        a.set_values(1.0)
        a.set_boundary_values(1.0)

        # Do this to get access to the matrix
        # This is also called inside elliptic_solve
        operator.update_elliptic_matrix(a)

        V1 = U[:n]
        V2 = U[n:]

        A = num.mat(operator.elliptic_matrix.todense())
        U = num.mat(U.reshape(6,1))

        #Setup up rhs as b = A u
        X = num.array(2*A*U).reshape(2,)

        b = Quantity(operator.domain)
        b.set_values(X, location='centroids')

        u_in.set_values(0.0)

        u_out = operator.elliptic_solve(u_in, b, a, iprint=1)

        assert num.allclose(u_out.centroid_values, V1)
        assert num.allclose(u_out.boundary_values, V2)

    def test_elliptic_solve_rectangular_cross(self):

        from anuga import rectangular_cross_domain

        m1 = 10
        n1 = 10
        domain = rectangular_cross_domain(m1,n1)

        # Diffusivity
        a = Quantity(domain)
        a.set_values(1.0)
        a.set_boundary_values(1.0)

        # Quantity to solve
        u = Quantity(domain)
        u.set_values(0.0)
        u.set_boundary_values(1.0)

        # Quantity for rhs
        b = Quantity(domain)
        b.set_values(0.0)
        b.set_boundary_values(0.0)

        operator = Kinematic_viscosity_operator(domain)

        n = operator.n
        tot_len = operator.tot_len

        u_out = operator.elliptic_solve(u, b, a, iprint=1)
    
        assert num.allclose(u_out.centroid_values, num.ones_like(u_out.centroid_values))
        assert num.allclose(u_out.boundary_values, num.ones_like(u_out.boundary_values))



    def test_parabolic_solve_one_triangle(self):
        operator = self.operator1()
        n = operator.n
        dt = operator.dt

        U = num.array([2.0,2.0,1.0,1.0])
        U_mod = num.array([10.0, 2.0, 1.0, 1.0])

        u_in = Quantity(operator.domain)
        u_in.set_values(U[:n], location='centroids')
        u_in.set_boundary_values(U_mod[n:])

        a = Quantity(operator.domain)
        a.set_values(1.0)
        a.set_boundary_values(1.0)


        V = num.array([2.0])
        A = num.array([-6.0-12.0/sqrt(5), 6.0,  6.0/sqrt(5), 6.0/sqrt(5)])

        #Setup up rhs
        X = U_mod[:n] - dt*2*num.array(num.mat(A)*num.mat(U_mod.reshape(4,1))).reshape(n,)
        b = Quantity(operator.domain)
        b.set_values(X, location='centroids')


        u_out = operator.parabolic_solve(u_in, b, a, iprint=1)


        assert num.allclose(u_out.centroid_values, U_mod[:n])


    def test_parabolic_solve_two_triangles(self):
        operator = self.operator2()
        n = operator.n
        nt = operator.tot_len

        dt = operator.dt

        U = num.array([2.0,3.0,1.0,1.0,4.0,3.0])
        U_mod = num.array([4.0,2.0,1.0,1.0,4.0,3.0])


        u_in = Quantity(operator.domain)
        u_in.set_values(U[:n], location='centroids')
        u_in.set_boundary_values(U_mod[n:])

        a = Quantity(operator.domain)
        a.set_values(1.0)
        a.set_boundary_values(1.0)

        operator.update_elliptic_matrix(a)


        A = num.array([[-8.36656315,  3., 2.68328157,  2.68328157,  0.,  0. ],
                       [ 3., -8.36656315 , 0. , 0. ,  2.68328157,  2.68328157]])

        assert num.allclose(A,operator.elliptic_matrix.todense())



        #Setup up rhs
        X = U_mod[:n] - dt*2*num.array(num.mat(A)*num.mat(U_mod.reshape(nt,1))).reshape(n,)
        b = Quantity(operator.domain)
        b.set_values(X, location='centroids')


        u_out = operator.parabolic_solve(u_in, b, a, iprint=1)


        assert num.allclose(u_out.centroid_values, U_mod[:n])

    def test_parabolic_solve_rectangular_cross(self):

        from anuga import rectangular_cross_domain

        m1 = 10
        n1 = 10
        domain = rectangular_cross_domain(m1,n1)

        # Diffusivity
        a = Quantity(domain)
        a.set_values(1.0)
        a.set_boundary_values(1.0)

        # Quantity initial condition
        u_in = Quantity(domain)
        #u_in.set_values( 0.0 )
        u_in.set_values(lambda x,y : 16.0*x*(1-x)*y*(1-y))
        u_in.set_boundary_values(0.0)

        # Quantity to solve
        u_mod = Quantity(domain)
        u_mod.set_values(lambda x,y : 15.9*x*(1-x)*y*(1-y) )
        u_mod.set_boundary_values(0.0)

        # Quantity for rhs
        b = Quantity(domain)
        b.set_values(0.0)
        b.set_boundary_values(0.0)

        operator = Kinematic_viscosity_operator(domain)

        dt = 0.01
        operator.dt = dt
        n = operator.n
        nt = operator.tot_len

        operator.update_elliptic_matrix(a)

        A = num.mat(operator.elliptic_matrix.todense())
        D = num.mat(operator.triangle_areas.todense())
        U_mod = num.concatenate( (u_mod.centroid_values, u_mod.boundary_values) )


        #Setup up rhs
        X = U_mod[:n] - dt*num.array(D*A*num.mat(U_mod.reshape(nt,1))).reshape(n,)
        b = Quantity(operator.domain)
        b.set_values(X, location='centroids')


        u_out = operator.parabolic_solve(u_in, b, a, iprint=1, use_dt_tol=False)

        assert num.allclose(u_out.centroid_values, U_mod[:n])


    def test_elliptic_solve_rectangular_cross_velocities(self):

        from anuga import rectangular_cross_domain
        from anuga import Reflective_boundary

        m1 = 10
        n1 = 10
        domain = rectangular_cross_domain(m1,n1)

        #
        domain.set_quantity('elevation', expression='x')
        domain.set_quantity('friction', 0.03)
        domain.set_quantity('stage',expression='elevation + 2*x')
        domain.set_quantity('xmomentum', expression='2*x+3*y')
        domain.set_quantity('ymomentum', expression='5*x+7*y')


        B = Reflective_boundary(domain)
        domain.set_boundary( {'left': B, 'right': B, 'top': B, 'bottom': B})

        domain.update_boundary()
        domain.update_centroids_of_velocities_and_height()


        a = domain.quantities['height']

        # Quantity to solve
        u = domain.quantities['xvelocity']
        u.set_boundary_values(1.0)


        v = domain.quantities['yvelocity']
        v.set_boundary_values(2.0)

        # Quantity for rhs
        b = Quantity(domain)
        b.set_values(0.0)
        b.set_boundary_values(0.0)

        kv = Kinematic_viscosity_operator(domain)

        n = kv.n
        tot_len = kv.tot_len

        kv.update_elliptic_matrix(a)

        u_out = kv.elliptic_solve(u, b, a, update_matrix=False, iprint=1)

        v_out = kv.elliptic_solve(v, b, a, update_matrix=False, iprint=1)

        assert num.allclose(u_out.centroid_values, num.ones_like(u_out.centroid_values))
        assert num.allclose(u_out.boundary_values, num.ones_like(u_out.boundary_values))

    def test_parabolic_solve_rectangular_cross_velocities(self):

        from anuga import rectangular_cross_domain
        from anuga import Reflective_boundary

        m1 = 10
        n1 = 10
        domain = rectangular_cross_domain(m1,n1)

        #
        domain.set_quantity('elevation', expression='x')
        domain.set_quantity('friction', 0.03)
        domain.set_quantity('stage',expression='elevation + 2*x')
        domain.set_quantity('xmomentum', expression='2*x+3*y')
        domain.set_quantity('ymomentum', expression='5*x+7*y')


        B = Reflective_boundary(domain)
        domain.set_boundary( {'left': B, 'right': B, 'top': B, 'bottom': B})

        domain.update_boundary()
        domain.update_centroids_of_velocities_and_height()




        h = domain.quantities['height']

        # Quantity to solve
        u = domain.quantities['xvelocity']
        u.set_boundary_values(1.0)


        v = domain.quantities['yvelocity']
        v.set_boundary_values(2.0)

        kv = Kinematic_viscosity_operator(domain)


        # let's make timestep large so that the final solution will look like
        #the solution of hte elliptic problem. In this case u -> 1, v -> 2.

        dt = 100.0
        kv.dt = dt
        n = kv.n
        nt = kv.tot_len

        kv.update_elliptic_matrix(h)

        kv.parabolic_solve(u, u, h, u_out=u, update_matrix=False, iprint=1, use_dt_tol=False)

        kv.parabolic_solve(v, v, h, u_out=v, update_matrix=False, iprint=1, use_dt_tol=False)


        #print u.centroid_values
        #print u.boundary_values
        assert num.allclose(u.centroid_values, num.ones_like(u.centroid_values), rtol=1.0e-1)
        assert num.allclose(u.boundary_values, num.ones_like(u.boundary_values))

        assert num.allclose(v.centroid_values, 2.0*num.ones_like(v.centroid_values), rtol=1.0e-1)
        assert num.allclose(v.boundary_values, 2.0*num.ones_like(v.boundary_values))


        domain.update_centroids_of_momentum_from_velocity()

        uh = domain.quantities['xmomentum']
        vh = domain.quantities['ymomentum']

        assert num.allclose(uh.centroid_values, u.centroid_values*h.centroid_values )
        assert num.allclose(vh.centroid_values, v.centroid_values*h.centroid_values )


    def test_parabolic_solve_rectangular_cross_velocities_zero_h(self):

        from anuga import rectangular_cross_domain
        from anuga import Reflective_boundary

        m1 = 5
        n1 = 5
        domain = rectangular_cross_domain(m1,n1)

        #
        domain.set_quantity('elevation', expression='x')
        domain.set_quantity('friction', 0.03)
        domain.set_quantity('stage',expression='elevation + 2*(x-0.5)')
        domain.set_quantity('xmomentum', expression='2*x+3*y')
        domain.set_quantity('ymomentum', expression='5*x+7*y')

 
        w = domain.quantities['stage']

        #print w.centroid_values
        #print w.boundary_values

        domain.distribute_to_vertices_and_edges()

        #print w.centroid_values
        #print w.boundary_values
        
        B = Reflective_boundary(domain)
        domain.set_boundary( {'left': B, 'right': B, 'top': B, 'bottom': B})

        domain.update_boundary()

        #print w.centroid_values
        #print w.boundary_values

        domain.update_centroids_of_velocities_and_height()




        h = domain.quantities['height']

        #print 'h'
        #print h.centroid_values
        #print h.boundary_values
        
        # Quantity to solve
        u = domain.quantities['xvelocity']
        u.set_boundary_values(1.0)


        #print 'u'
        #print u.centroid_values
        #print u.boundary_values

        v = domain.quantities['yvelocity']
        v.set_boundary_values(2.0)

        kv = Kinematic_viscosity_operator(domain)


        # let's make timestep large so that the final solution will look like
        #the solution of hte elliptic problem. In this case u -> 1, v -> 2.

        dt = 1000.0
        kv.dt = dt
        n = kv.n
        nt = kv.tot_len

        kv.update_elliptic_matrix(h)

        kv.parabolic_solve(u, u, h, u_out=u, update_matrix=False, iprint=1, use_dt_tol=False)

        kv.parabolic_solve(v, v, h, u_out=v, update_matrix=False, iprint=1, use_dt_tol=False)

        assert num.allclose(u.centroid_values, num.where(h.centroid_values > 0.0, 1.0, 0.0), rtol=1.0e-1)
        assert num.allclose(u.boundary_values, num.ones_like(u.boundary_values))

        assert num.allclose(v.centroid_values, num.where(h.centroid_values > 0.0, 2.0, 0.0), rtol=1.0e-1)
        assert num.allclose(v.boundary_values, 2.0*num.ones_like(v.boundary_values))


        domain.update_centroids_of_momentum_from_velocity()

        domain.distribute_to_vertices_and_edges()
        
        uh = domain.quantities['xmomentum']
        vh = domain.quantities['ymomentum']

        #print 'uh'
        #print uh.centroid_values
        #print uh.boundary_values

        assert num.allclose(uh.centroid_values, u.centroid_values*h.centroid_values )
        assert num.allclose(vh.centroid_values, v.centroid_values*h.centroid_values )

    def test_kinematic_operator_default(self):

        from anuga import rectangular_cross_domain
        from anuga import Reflective_boundary

        m1 = 10
        n1 = 10
        domain = rectangular_cross_domain(m1,n1)

        #domain.set_flow_algorithm('2_0')

        #
        domain.set_quantity('elevation', expression='x')
        domain.set_quantity('friction', 0.03)
        domain.set_quantity('stage',expression='elevation + 2*(x-0.5)')
        domain.set_quantity('xmomentum', expression='2*x+3*y')
        domain.set_quantity('ymomentum', expression='5*x+7*y')

        B = Reflective_boundary(domain)
        domain.set_boundary( {'left': B, 'right': B, 'top': B, 'bottom': B})

        # kill off the wave with viscosity
        kv = Kinematic_viscosity_operator(domain)


        # let's make timestep large so that the final solution will look like
        #the solution of hte elliptic problem. In this case u -> 1, v -> 2.


        for t in domain.evolve(yieldstep = 1.0, finaltime = 10.0):
            #domain.write_time()
            #domain.print_operator_timestepping_statistics()
            pass

#
        w  = domain.quantities['stage']
        uh = domain.quantities['xmomentum']
        vh = domain.quantities['ymomentum']

        #print 'uh'
        #print uh.centroid_values
        #print uh.boundary_values

        #print 'w'
        #print w.centroid_values

        #from pprint import pprint
        #pprint(w.centroid_values)


        wc = num.array([
        0.70714365,  0.70714416,  0.70714295,  0.70714222,  0.70714486,
        0.70714507,  0.70714374,  0.70714601,  0.70714492,  0.70714425,
        0.70714595,  0.70714437,  0.70714797,  0.70714691,  0.70714697,
        0.70714845,  0.70714793,  0.70714793,  0.70715033,  0.70714852,
        0.70715244,  0.70715018,  0.70715176,  0.70715224,  0.70715211,
        0.70715265,  0.70715351,  0.7071531 ,  0.70715433,  0.70715309,
        0.70715351,  0.70715472,  0.70715429,  0.70715433,  0.70715487,
        0.70715523,  0.7071545 ,  0.70715446,  0.70715317,  0.70715564,
        0.70714142,  0.70714198,  0.70714079,  0.70714299,  0.70714482,
        0.70714378,  0.70714344,  0.70714377,  0.7071443 ,  0.70714533,
        0.70714579,  0.70714574,  0.70714906,  0.70714717,  0.70714819,
        0.70714822,  0.70714976,  0.70714952,  0.70715093,  0.70715077,
        0.70715217,  0.70715094,  0.70715291,  0.70715188,  0.70715352,
        0.70715278,  0.707154  ,  0.70715429,  0.70715376,  0.70715309,
        0.70715446,  0.70715422,  0.70715366,  0.70715453,  0.70715413,
        0.70715539,  0.70715385,  0.70715412,  0.70715154,  0.70715306,
        0.70714038,  0.70713905,  0.7071358 ,  0.70713972,  0.70714303,
        0.7071419 ,  0.70714066,  0.70714219,  0.7071459 ,  0.70714505,
        0.70714639,  0.70714648,  0.70714833,  0.70714827,  0.70715147,
        0.70715013,  0.70715194,  0.70715133,  0.70715542,  0.70715345,
        0.70715296,  0.70715417,  0.70715676,  0.70715521,  0.70715526,
        0.7071548 ,  0.70715717,  0.70715512,  0.70715381,  0.70715523,
        0.70715556,  0.70715486,  0.70715482,  0.70715338,  0.70715307,
        0.70715381,  0.70715132,  0.70715182,  0.70714789,  0.70715086,
        0.70713443,  0.70713559,  0.70713539,  0.70713615,  0.70714057,
        0.70713978,  0.70714091,  0.70714102,  0.70714618,  0.70714338,
        0.70714803,  0.70714858,  0.7071519 ,  0.70715029,  0.70715343,
        0.70715461,  0.70715589,  0.70715519,  0.7071565 ,  0.70715796,
        0.70715738,  0.70715845,  0.7071601 ,  0.70715829,  0.70715711,
        0.70715903,  0.70716011,  0.70715714,  0.7071565 ,  0.70715756,
        0.70715885,  0.7071556 ,  0.70715386,  0.70715406,  0.70715653,
        0.70715532,  0.70714813,  0.7071515 ,  0.70715242,  0.70715269,
        0.70713191,  0.70712961,  0.70712505,  0.70712841,  0.70714097,
        0.70713808,  0.70713862,  0.7071431 ,  0.70714966,  0.7071463 ,
        0.70715775,  0.70715666,  0.70715566,  0.7071554 ,  0.7071632 ,
        0.70716353,  0.70715928,  0.70716244,  0.70716736,  0.70716495,
        0.70716301,  0.70716635,  0.70717088,  0.70716792,  0.70716369,
        0.70717007,  0.7071741 ,  0.70716769,  0.70716166,  0.70716991,
        0.70717294,  0.70716167,  0.70715775,  0.70716057,  0.70715687,
        0.70715535,  0.70715014,  0.70714766,  0.70714559,  0.70714992,
        0.7071149 ,  0.70708741,  0.706984  ,  0.70711096,  0.70714367,
        0.70714831,  0.70713519,  0.7071811 ,  0.70716622,  0.70716603,
        0.70714155,  0.7071748 ,  0.70716885,  0.70716897,  0.70713548,
        0.70716966,  0.70716924,  0.70716978,  0.70713561,  0.7071717 ,
        0.70717389,  0.7071726 ,  0.70713926,  0.70717593,  0.70718002,
        0.70717761,  0.70714428,  0.70718053,  0.70718062,  0.70718719,
        0.70715731,  0.70718271,  0.70716238,  0.7071992 ,  0.70715496,
        0.70716834,  0.70713531,  0.70713099,  0.70700665,  0.7071098 ,
        0.70634397,  0.70524618,  0.70297607,  0.70514658,  0.70658259,
        0.70506628,  0.70244401,  0.70497884,  0.70657086,  0.70498266,
        0.70239779,  0.70496243,  0.7065572 ,  0.7049646 ,  0.70239608,
        0.70496008,  0.70655538,  0.70496125,  0.70239685,  0.70496177,
        0.70655883,  0.70496295,  0.70239957,  0.70496624,  0.70656625,
        0.70496724,  0.70240482,  0.7049756 ,  0.70658803,  0.70497608,
        0.70241139,  0.70500006,  0.70660425,  0.70499778,  0.70246225,
        0.70508764,  0.70636798,  0.70516922,  0.70299639,  0.70526838,
        0.71780931,  0.7506157 ,  0.78399529,  0.75061024,  0.71769206,
        0.75059929,  0.78398287,  0.75059279,  0.71768281,  0.75059112,
        0.78397863,  0.75059025,  0.71768261,  0.75058996,  0.78397777,
        0.75058981,  0.71768268,  0.75058969,  0.78397749,  0.75058967,
        0.7176832 ,  0.75058972,  0.78397772,  0.75058986,  0.71768421,
        0.7505901 ,  0.78397859,  0.75059043,  0.71768534,  0.7505909 ,
        0.78398028,  0.750592  ,  0.71769545,  0.75059388,  0.78398545,
        0.75060056,  0.71781337,  0.75061163,  0.78399848,  0.75061714,
        0.81739069,  0.85076296,  0.8841241 ,  0.85076174,  0.81738381,
        0.85075988,  0.88412183,  0.85075808,  0.81738087,  0.85075718,
        0.88412031,  0.85075635,  0.81737996,  0.85075599,  0.88411952,
        0.85075563,  0.81737963,  0.85075548,  0.88411919,  0.8507555 ,
        0.81738003,  0.85075569,  0.88411972,  0.85075629,  0.81738134,
        0.85075692,  0.88412133,  0.85075812,  0.81738361,  0.85075914,
        0.88412387,  0.85076103,  0.81738807,  0.85076269,  0.88412739,
        0.85076547,  0.81739598,  0.85076786,  0.88413107,  0.85076949,
        0.91748914,  0.95083916,  0.98417801,  0.95083906,  0.91748809,
        0.95083882,  0.98417779,  0.95083863,  0.91748731,  0.95083843,
        0.98417752,  0.9508382 ,  0.91748674,  0.950838  ,  0.9841771 ,
        0.95083776,  0.91748646,  0.95083764,  0.98417686,  0.95083771,
        0.91748702,  0.95083794,  0.98417744,  0.95083859,  0.91748864,
        0.95083927,  0.98417906,  0.95084046,  0.91749107,  0.95084145,
        0.98418138,  0.95084291,  0.91749397,  0.95084401,  0.98418384,
        0.95084538,  0.91749653,  0.95084626,  0.98418563,  0.95084686])


        #print w.centroid_values - wc

        #print max(w.centroid_values - wc)

        assert num.allclose(w.centroid_values, wc, rtol=2.0e-3)

    def test_kinematic_operator_quantity(self):

        from anuga import rectangular_cross_domain
        from anuga import Reflective_boundary

        m1 = 10
        n1 = 10
        domain = rectangular_cross_domain(m1,n1)

        #domain.set_flow_algorithm('2_0')

        #
        domain.set_quantity('elevation', expression='x')
        domain.set_quantity('friction', 0.03)
        domain.set_quantity('stage',expression='elevation + 2*(x-0.5)')
        domain.set_quantity('xmomentum', expression='2*x+3*y')
        domain.set_quantity('ymomentum', expression='5*x+7*y')

        B = Reflective_boundary(domain)
        domain.set_boundary( {'left': B, 'right': B, 'top': B, 'bottom': B})


        Q = Quantity(domain)
        Q = 2.0
        # kill off the wave with viscosity
        kv = Kinematic_viscosity_operator(domain, diffusivity = Q)


        # let's make timestep large so that the final solution will look like
        #the solution of hte elliptic problem. In this case u -> 1, v -> 2.


        for t in domain.evolve(yieldstep = 1.0, finaltime = 10.0):
            #domain.write_time()
            #domain.print_operator_timestepping_statistics()
            pass

#
        w  = domain.quantities['stage']
        uh = domain.quantities['xmomentum']
        vh = domain.quantities['ymomentum']

        #print 'uh'
        #print uh.centroid_values
        #print uh.boundary_values

        #print 'w'
        #print w.centroid_values

        #from pprint import pprint
        #pprint(w.centroid_values)


        wc = num.array([
        0.71624029,  0.71622927,  0.71621675,  0.71623888,  0.71624236,
        0.71624536,  0.71625157,  0.71625028,  0.71625679,  0.71626609,
        0.71630233,  0.71627457,  0.71627721,  0.71628666,  0.71633484,
        0.71629002,  0.71628494,  0.716295  ,  0.7163438 ,  0.71629656,
        0.71628493,  0.71629656,  0.71634379,  0.71629497,  0.71627716,
        0.71628999,  0.71633481,  0.7162866 ,  0.71625666,  0.71627448,
        0.71630224,  0.71626596,  0.71624212,  0.7162501 ,  0.7162514 ,
        0.71624512,  0.71624   ,  0.7162386 ,  0.71621644,  0.71622896,
        0.71619869,  0.71615658,  0.71609423,  0.71619602,  0.71627164,
        0.71623926,  0.71625039,  0.71633719,  0.71638922,  0.71642539,
        0.71652642,  0.71649892,  0.71646671,  0.71653525,  0.71670614,
        0.71661869,  0.71649067,  0.71663318,  0.71682302,  0.71665878,
        0.71649066,  0.71665876,  0.71682295,  0.71663309,  0.71646665,
        0.71661859,  0.71670596,  0.71653511,  0.71638911,  0.71649877,
        0.71652622,  0.71642523,  0.71627151,  0.716337  ,  0.71625001,
        0.71623888,  0.7161983 ,  0.71619554,  0.71609371,  0.71615611,
        0.71587901,  0.71555375,  0.71521927,  0.71573946,  0.71615663,
        0.71586493,  0.7156413 ,  0.71615004,  0.71653474,  0.71632223,
        0.71618825,  0.7165586 ,  0.7168124 ,  0.71668994,  0.71661036,
        0.7168446 ,  0.71694587,  0.71689337,  0.7167922 ,  0.71693225,
        0.71694582,  0.71693224,  0.71679212,  0.71689325,  0.71681216,
        0.71684437,  0.71661004,  0.71668963,  0.71653449,  0.71655826,
        0.71618788,  0.71632191,  0.71615622,  0.71614967,  0.71564092,
        0.71586446,  0.7158785 ,  0.71573897,  0.71521879,  0.71555323,
        0.71415117,  0.71304803,  0.71200401,  0.71333356,  0.71459491,
        0.71350761,  0.71272705,  0.7140006 ,  0.71526042,  0.71418365,
        0.71337479,  0.7146592 ,  0.71582149,  0.71478585,  0.71378284,
        0.7150456 ,  0.71605221,  0.71509271,  0.71396254,  0.71516103,
        0.71605211,  0.71516102,  0.71396249,  0.71509256,  0.71582115,
        0.7150454 ,  0.71378271,  0.71478555,  0.71526005,  0.71465889,
        0.71337454,  0.71418329,  0.71459453,  0.71400022,  0.71272682,
        0.71350725,  0.71415077,  0.71333321,  0.71200389,  0.71304774,
        0.70944126,  0.70705883,  0.70442227,  0.70714215,  0.70999341,
        0.70722667,  0.70436187,  0.70745337,  0.71044978,  0.70748596,
        0.70427781,  0.70768146,  0.71082549,  0.70772906,  0.70426793,
        0.70786303,  0.71099495,  0.70788365,  0.70424722,  0.70791928,
        0.71099502,  0.70791937,  0.70424774,  0.70788396,  0.71082556,
        0.70786332,  0.70426849,  0.70772935,  0.71044982,  0.70768178,
        0.7042786 ,  0.70748637,  0.70999356,  0.70745385,  0.70436311,
        0.70722738,  0.70944169,  0.70714295,  0.70442389,  0.70705981,
        0.69895933,  0.69463188,  0.68921358,  0.693824  ,  0.698153  ,
        0.69349963,  0.68725093,  0.69221842,  0.69728195,  0.69180649,
        0.68463972,  0.69053046,  0.69673179,  0.69018397,  0.68236173,
        0.68940762,  0.69650961,  0.68925397,  0.68125059,  0.68902719,
        0.69651034,  0.68902736,  0.6812516 ,  0.68925556,  0.69673305,
        0.6894096 ,  0.6823656 ,  0.69018707,  0.69728407,  0.69053386,
        0.68464522,  0.69181074,  0.69815588,  0.69222279,  0.68725717,
        0.69350432,  0.69896255,  0.69382873,  0.68922015,  0.69463687,
        0.68375896,  0.6882601 ,  0.69595562,  0.68766298,  0.68105558,
        0.68673658,  0.69502847,  0.68542815,  0.67770965,  0.68435344,
        0.69409778,  0.68310537,  0.67491515,  0.68222458,  0.69337943,
        0.68140117,  0.67356609,  0.68097711,  0.69301997,  0.68071631,
        0.67356716,  0.68071666,  0.69302027,  0.68097852,  0.6749196 ,
        0.68140363,  0.69338045,  0.68222808,  0.6777162 ,  0.68310954,
        0.69409929,  0.68435822,  0.68106317,  0.68543327,  0.69503026,
        0.68674199,  0.68376697,  0.68766854,  0.69595754,  0.68826575,
        0.71760631,  0.75094294,  0.78427898,  0.75094168,  0.71760193,
        0.75093986,  0.78427453,  0.75093415,  0.71758272,  0.7509278 ,
        0.78426295,  0.75091754,  0.7175572 ,  0.75090919,  0.78424795,
        0.75089856,  0.71753518,  0.75089163,  0.78423642,  0.75088684,
        0.71753524,  0.75088686,  0.78423643,  0.75089171,  0.7175573 ,
        0.75089864,  0.78424798,  0.75090931,  0.71758285,  0.75091768,
        0.78426303,  0.75092799,  0.7176021 ,  0.75093438,  0.78427472,
        0.75094013,  0.71760652,  0.75094199,  0.78427929,  0.75094328,
        0.81761649,  0.85095268,  0.88428788,  0.85095192,  0.81761311,
        0.8509508 ,  0.88428574,  0.85094833,  0.81760513,  0.8509458 ,
        0.88428131,  0.85094197,  0.81759506,  0.85093883,  0.88427596,
        0.85093514,  0.81758753,  0.85093282,  0.88427212,  0.85093123,
        0.81758749,  0.8509312 ,  0.88427198,  0.85093269,  0.81759494,
        0.85093494,  0.8842756 ,  0.85093857,  0.81760502,  0.8509417 ,
        0.88428088,  0.85094557,  0.81761314,  0.85094816,  0.88428543,
        0.85095073,  0.81761667,  0.85095193,  0.88428775,  0.85095275,
        0.91762366,  0.95095836,  0.98429205,  0.95095804,  0.91762217,
        0.95095754,  0.98429102,  0.95095658,  0.91761918,  0.95095558,
        0.98428903,  0.95095416,  0.91761561,  0.95095297,  0.98428667,
        0.95095164,  0.91761304,  0.95095078,  0.98428497,  0.95095015,
        0.91761286,  0.95095007,  0.98428475,  0.95095045,  0.9176151 ,
        0.95095115,  0.98428605,  0.95095231,  0.91761853,  0.95095342,
        0.9842882 ,  0.9509548 ,  0.91762161,  0.95095583,  0.98429026,
        0.95095688,  0.91762327,  0.95095746,  0.98429146,  0.95095784])


        #print max(w.centroid_values- wc)
        assert num.allclose(w.centroid_values, wc, rtol=0.05)

    def test_kinematic_operator_number(self):

        from anuga import rectangular_cross_domain
        from anuga import Reflective_boundary

        m1 = 10
        n1 = 10
        domain = rectangular_cross_domain(m1,n1)

        #domain.set_flow_algorithm('2_0')

        #
        domain.set_quantity('elevation', expression='x')
        domain.set_quantity('friction', 0.03)
        domain.set_quantity('stage',expression='elevation + 2*(x-0.5)')
        domain.set_quantity('xmomentum', expression='2*x+3*y')
        domain.set_quantity('ymomentum', expression='5*x+7*y')

        B = Reflective_boundary(domain)
        domain.set_boundary( {'left': B, 'right': B, 'top': B, 'bottom': B})

        # kill off the wave with viscosity
        kv = Kinematic_viscosity_operator(domain, diffusivity=2.0)


        # let's make timestep large so that the final solution will look like
        #the solution of hte elliptic problem. In this case u -> 1, v -> 2.


        for t in domain.evolve(yieldstep = 1.0, finaltime = 10.0):
            #domain.write_time()
            #domain.print_operator_timestepping_statistics()
            pass

#
        w  = domain.quantities['stage']
        uh = domain.quantities['xmomentum']
        vh = domain.quantities['ymomentum']

        #print 'uh'
        #print uh.centroid_values
        #print uh.boundary_values

        #print 'w'
        #print w.centroid_values

        #from pprint import pprint
        #pprint(w.centroid_values)


        wc = num.array([
        0.71624029,  0.71622927,  0.71621675,  0.71623888,  0.71624236,
        0.71624536,  0.71625157,  0.71625028,  0.71625679,  0.71626609,
        0.71630233,  0.71627457,  0.71627721,  0.71628666,  0.71633484,
        0.71629002,  0.71628494,  0.716295  ,  0.7163438 ,  0.71629656,
        0.71628493,  0.71629656,  0.71634379,  0.71629497,  0.71627716,
        0.71628999,  0.71633481,  0.7162866 ,  0.71625666,  0.71627448,
        0.71630224,  0.71626596,  0.71624212,  0.7162501 ,  0.7162514 ,
        0.71624512,  0.71624   ,  0.7162386 ,  0.71621644,  0.71622896,
        0.71619869,  0.71615658,  0.71609423,  0.71619602,  0.71627164,
        0.71623926,  0.71625039,  0.71633719,  0.71638922,  0.71642539,
        0.71652642,  0.71649892,  0.71646671,  0.71653525,  0.71670614,
        0.71661869,  0.71649067,  0.71663318,  0.71682302,  0.71665878,
        0.71649066,  0.71665876,  0.71682295,  0.71663309,  0.71646665,
        0.71661859,  0.71670596,  0.71653511,  0.71638911,  0.71649877,
        0.71652622,  0.71642523,  0.71627151,  0.716337  ,  0.71625001,
        0.71623888,  0.7161983 ,  0.71619554,  0.71609371,  0.71615611,
        0.71587901,  0.71555375,  0.71521927,  0.71573946,  0.71615663,
        0.71586493,  0.7156413 ,  0.71615004,  0.71653474,  0.71632223,
        0.71618825,  0.7165586 ,  0.7168124 ,  0.71668994,  0.71661036,
        0.7168446 ,  0.71694587,  0.71689337,  0.7167922 ,  0.71693225,
        0.71694582,  0.71693224,  0.71679212,  0.71689325,  0.71681216,
        0.71684437,  0.71661004,  0.71668963,  0.71653449,  0.71655826,
        0.71618788,  0.71632191,  0.71615622,  0.71614967,  0.71564092,
        0.71586446,  0.7158785 ,  0.71573897,  0.71521879,  0.71555323,
        0.71415117,  0.71304803,  0.71200401,  0.71333356,  0.71459491,
        0.71350761,  0.71272705,  0.7140006 ,  0.71526042,  0.71418365,
        0.71337479,  0.7146592 ,  0.71582149,  0.71478585,  0.71378284,
        0.7150456 ,  0.71605221,  0.71509271,  0.71396254,  0.71516103,
        0.71605211,  0.71516102,  0.71396249,  0.71509256,  0.71582115,
        0.7150454 ,  0.71378271,  0.71478555,  0.71526005,  0.71465889,
        0.71337454,  0.71418329,  0.71459453,  0.71400022,  0.71272682,
        0.71350725,  0.71415077,  0.71333321,  0.71200389,  0.71304774,
        0.70944126,  0.70705883,  0.70442227,  0.70714215,  0.70999341,
        0.70722667,  0.70436187,  0.70745337,  0.71044978,  0.70748596,
        0.70427781,  0.70768146,  0.71082549,  0.70772906,  0.70426793,
        0.70786303,  0.71099495,  0.70788365,  0.70424722,  0.70791928,
        0.71099502,  0.70791937,  0.70424774,  0.70788396,  0.71082556,
        0.70786332,  0.70426849,  0.70772935,  0.71044982,  0.70768178,
        0.7042786 ,  0.70748637,  0.70999356,  0.70745385,  0.70436311,
        0.70722738,  0.70944169,  0.70714295,  0.70442389,  0.70705981,
        0.69895933,  0.69463188,  0.68921358,  0.693824  ,  0.698153  ,
        0.69349963,  0.68725093,  0.69221842,  0.69728195,  0.69180649,
        0.68463972,  0.69053046,  0.69673179,  0.69018397,  0.68236173,
        0.68940762,  0.69650961,  0.68925397,  0.68125059,  0.68902719,
        0.69651034,  0.68902736,  0.6812516 ,  0.68925556,  0.69673305,
        0.6894096 ,  0.6823656 ,  0.69018707,  0.69728407,  0.69053386,
        0.68464522,  0.69181074,  0.69815588,  0.69222279,  0.68725717,
        0.69350432,  0.69896255,  0.69382873,  0.68922015,  0.69463687,
        0.68375896,  0.6882601 ,  0.69595562,  0.68766298,  0.68105558,
        0.68673658,  0.69502847,  0.68542815,  0.67770965,  0.68435344,
        0.69409778,  0.68310537,  0.67491515,  0.68222458,  0.69337943,
        0.68140117,  0.67356609,  0.68097711,  0.69301997,  0.68071631,
        0.67356716,  0.68071666,  0.69302027,  0.68097852,  0.6749196 ,
        0.68140363,  0.69338045,  0.68222808,  0.6777162 ,  0.68310954,
        0.69409929,  0.68435822,  0.68106317,  0.68543327,  0.69503026,
        0.68674199,  0.68376697,  0.68766854,  0.69595754,  0.68826575,
        0.71760631,  0.75094294,  0.78427898,  0.75094168,  0.71760193,
        0.75093986,  0.78427453,  0.75093415,  0.71758272,  0.7509278 ,
        0.78426295,  0.75091754,  0.7175572 ,  0.75090919,  0.78424795,
        0.75089856,  0.71753518,  0.75089163,  0.78423642,  0.75088684,
        0.71753524,  0.75088686,  0.78423643,  0.75089171,  0.7175573 ,
        0.75089864,  0.78424798,  0.75090931,  0.71758285,  0.75091768,
        0.78426303,  0.75092799,  0.7176021 ,  0.75093438,  0.78427472,
        0.75094013,  0.71760652,  0.75094199,  0.78427929,  0.75094328,
        0.81761649,  0.85095268,  0.88428788,  0.85095192,  0.81761311,
        0.8509508 ,  0.88428574,  0.85094833,  0.81760513,  0.8509458 ,
        0.88428131,  0.85094197,  0.81759506,  0.85093883,  0.88427596,
        0.85093514,  0.81758753,  0.85093282,  0.88427212,  0.85093123,
        0.81758749,  0.8509312 ,  0.88427198,  0.85093269,  0.81759494,
        0.85093494,  0.8842756 ,  0.85093857,  0.81760502,  0.8509417 ,
        0.88428088,  0.85094557,  0.81761314,  0.85094816,  0.88428543,
        0.85095073,  0.81761667,  0.85095193,  0.88428775,  0.85095275,
        0.91762366,  0.95095836,  0.98429205,  0.95095804,  0.91762217,
        0.95095754,  0.98429102,  0.95095658,  0.91761918,  0.95095558,
        0.98428903,  0.95095416,  0.91761561,  0.95095297,  0.98428667,
        0.95095164,  0.91761304,  0.95095078,  0.98428497,  0.95095015,
        0.91761286,  0.95095007,  0.98428475,  0.95095045,  0.9176151 ,
        0.95095115,  0.98428605,  0.95095231,  0.91761853,  0.95095342,
        0.9842882 ,  0.9509548 ,  0.91762161,  0.95095583,  0.98429026,
        0.95095688,  0.91762327,  0.95095746,  0.98429146,  0.95095784])




        #print w.centroid_values - wc

        #print max(w.centroid_values - wc)
        assert num.allclose(w.centroid_values, wc, rtol=0.05)




    def test_kinematic_operator_string(self):

        from anuga import rectangular_cross_domain
        from anuga import Reflective_boundary

        m1 = 10
        n1 = 10
        domain = rectangular_cross_domain(m1,n1)

        #domain.set_flow_algorithm('2_0')

        #
        domain.set_quantity('elevation', expression='x')
        domain.set_quantity('friction', 0.03)
        domain.set_quantity('stage',expression='elevation + 2*(x-0.5)')
        domain.set_quantity('xmomentum', expression='2*x+3*y')
        domain.set_quantity('ymomentum', expression='5*x+7*y')

        B = Reflective_boundary(domain)
        domain.set_boundary( {'left': B, 'right': B, 'top': B, 'bottom': B})

        # kill off the wave with viscosity
        kv = Kinematic_viscosity_operator(domain, diffusivity = 'height')


        # let's make timestep large so that the final solution will look like
        #the solution of hte elliptic problem. In this case u -> 1, v -> 2.


        for t in domain.evolve(yieldstep = 1.0, finaltime = 10.0):
            #domain.write_time()
            #domain.print_operator_timestepping_statistics()
            pass

#
        w  = domain.quantities['stage']
        uh = domain.quantities['xmomentum']
        vh = domain.quantities['ymomentum']

        #print 'uh'
        #print uh.centroid_values
        #print uh.boundary_values

        #print 'w'
        #print w.centroid_values

        #from pprint import pprint
        #pprint(w.centroid_values)


        wc = num.array([
        0.70714365,  0.70714416,  0.70714295,  0.70714222,  0.70714486,
        0.70714507,  0.70714374,  0.70714601,  0.70714492,  0.70714425,
        0.70714595,  0.70714437,  0.70714797,  0.70714691,  0.70714697,
        0.70714845,  0.70714793,  0.70714793,  0.70715033,  0.70714852,
        0.70715244,  0.70715018,  0.70715176,  0.70715224,  0.70715211,
        0.70715265,  0.70715351,  0.7071531 ,  0.70715433,  0.70715309,
        0.70715351,  0.70715472,  0.70715429,  0.70715433,  0.70715487,
        0.70715523,  0.7071545 ,  0.70715446,  0.70715317,  0.70715564,
        0.70714142,  0.70714198,  0.70714079,  0.70714299,  0.70714482,
        0.70714378,  0.70714344,  0.70714377,  0.7071443 ,  0.70714533,
        0.70714579,  0.70714574,  0.70714906,  0.70714717,  0.70714819,
        0.70714822,  0.70714976,  0.70714952,  0.70715093,  0.70715077,
        0.70715217,  0.70715094,  0.70715291,  0.70715188,  0.70715352,
        0.70715278,  0.707154  ,  0.70715429,  0.70715376,  0.70715309,
        0.70715446,  0.70715422,  0.70715366,  0.70715453,  0.70715413,
        0.70715539,  0.70715385,  0.70715412,  0.70715154,  0.70715306,
        0.70714038,  0.70713905,  0.7071358 ,  0.70713972,  0.70714303,
        0.7071419 ,  0.70714066,  0.70714219,  0.7071459 ,  0.70714505,
        0.70714639,  0.70714648,  0.70714833,  0.70714827,  0.70715147,
        0.70715013,  0.70715194,  0.70715133,  0.70715542,  0.70715345,
        0.70715296,  0.70715417,  0.70715676,  0.70715521,  0.70715526,
        0.7071548 ,  0.70715717,  0.70715512,  0.70715381,  0.70715523,
        0.70715556,  0.70715486,  0.70715482,  0.70715338,  0.70715307,
        0.70715381,  0.70715132,  0.70715182,  0.70714789,  0.70715086,
        0.70713443,  0.70713559,  0.70713539,  0.70713615,  0.70714057,
        0.70713978,  0.70714091,  0.70714102,  0.70714618,  0.70714338,
        0.70714803,  0.70714858,  0.7071519 ,  0.70715029,  0.70715343,
        0.70715461,  0.70715589,  0.70715519,  0.7071565 ,  0.70715796,
        0.70715738,  0.70715845,  0.7071601 ,  0.70715829,  0.70715711,
        0.70715903,  0.70716011,  0.70715714,  0.7071565 ,  0.70715756,
        0.70715885,  0.7071556 ,  0.70715386,  0.70715406,  0.70715653,
        0.70715532,  0.70714813,  0.7071515 ,  0.70715242,  0.70715269,
        0.70713191,  0.70712961,  0.70712505,  0.70712841,  0.70714097,
        0.70713808,  0.70713862,  0.7071431 ,  0.70714966,  0.7071463 ,
        0.70715775,  0.70715666,  0.70715566,  0.7071554 ,  0.7071632 ,
        0.70716353,  0.70715928,  0.70716244,  0.70716736,  0.70716495,
        0.70716301,  0.70716635,  0.70717088,  0.70716792,  0.70716369,
        0.70717007,  0.7071741 ,  0.70716769,  0.70716166,  0.70716991,
        0.70717294,  0.70716167,  0.70715775,  0.70716057,  0.70715687,
        0.70715535,  0.70715014,  0.70714766,  0.70714559,  0.70714992,
        0.7071149 ,  0.70708741,  0.706984  ,  0.70711096,  0.70714367,
        0.70714831,  0.70713519,  0.7071811 ,  0.70716622,  0.70716603,
        0.70714155,  0.7071748 ,  0.70716885,  0.70716897,  0.70713548,
        0.70716966,  0.70716924,  0.70716978,  0.70713561,  0.7071717 ,
        0.70717389,  0.7071726 ,  0.70713926,  0.70717593,  0.70718002,
        0.70717761,  0.70714428,  0.70718053,  0.70718062,  0.70718719,
        0.70715731,  0.70718271,  0.70716238,  0.7071992 ,  0.70715496,
        0.70716834,  0.70713531,  0.70713099,  0.70700665,  0.7071098 ,
        0.70634397,  0.70524618,  0.70297607,  0.70514658,  0.70658259,
        0.70506628,  0.70244401,  0.70497884,  0.70657086,  0.70498266,
        0.70239779,  0.70496243,  0.7065572 ,  0.7049646 ,  0.70239608,
        0.70496008,  0.70655538,  0.70496125,  0.70239685,  0.70496177,
        0.70655883,  0.70496295,  0.70239957,  0.70496624,  0.70656625,
        0.70496724,  0.70240482,  0.7049756 ,  0.70658803,  0.70497608,
        0.70241139,  0.70500006,  0.70660425,  0.70499778,  0.70246225,
        0.70508764,  0.70636798,  0.70516922,  0.70299639,  0.70526838,
        0.71780931,  0.7506157 ,  0.78399529,  0.75061024,  0.71769206,
        0.75059929,  0.78398287,  0.75059279,  0.71768281,  0.75059112,
        0.78397863,  0.75059025,  0.71768261,  0.75058996,  0.78397777,
        0.75058981,  0.71768268,  0.75058969,  0.78397749,  0.75058967,
        0.7176832 ,  0.75058972,  0.78397772,  0.75058986,  0.71768421,
        0.7505901 ,  0.78397859,  0.75059043,  0.71768534,  0.7505909 ,
        0.78398028,  0.750592  ,  0.71769545,  0.75059388,  0.78398545,
        0.75060056,  0.71781337,  0.75061163,  0.78399848,  0.75061714,
        0.81739069,  0.85076296,  0.8841241 ,  0.85076174,  0.81738381,
        0.85075988,  0.88412183,  0.85075808,  0.81738087,  0.85075718,
        0.88412031,  0.85075635,  0.81737996,  0.85075599,  0.88411952,
        0.85075563,  0.81737963,  0.85075548,  0.88411919,  0.8507555 ,
        0.81738003,  0.85075569,  0.88411972,  0.85075629,  0.81738134,
        0.85075692,  0.88412133,  0.85075812,  0.81738361,  0.85075914,
        0.88412387,  0.85076103,  0.81738807,  0.85076269,  0.88412739,
        0.85076547,  0.81739598,  0.85076786,  0.88413107,  0.85076949,
        0.91748914,  0.95083916,  0.98417801,  0.95083906,  0.91748809,
        0.95083882,  0.98417779,  0.95083863,  0.91748731,  0.95083843,
        0.98417752,  0.9508382 ,  0.91748674,  0.950838  ,  0.9841771 ,
        0.95083776,  0.91748646,  0.95083764,  0.98417686,  0.95083771,
        0.91748702,  0.95083794,  0.98417744,  0.95083859,  0.91748864,
        0.95083927,  0.98417906,  0.95084046,  0.91749107,  0.95084145,
        0.98418138,  0.95084291,  0.91749397,  0.95084401,  0.98418384,
        0.95084538,  0.91749653,  0.95084626,  0.98418563,  0.95084686])


        #print max(w.centroid_values - wc)

        assert num.allclose(w.centroid_values, wc, rtol=2.0e-3)


################################################################################

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_kinematic_viscosity, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)

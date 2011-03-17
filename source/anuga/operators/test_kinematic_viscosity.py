import operator
from anuga import Domain
from anuga import Quantity
from anuga import Dirichlet_boundary
from kinematic_viscosity import Kinematic_Viscosity

import numpy as num
from math import sqrt
import unittest

class Test_Kinematic_Viscosity(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
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
        
        return Kinematic_Viscosity(domain)

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



        return Kinematic_Viscosity(domain)

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

        a.set_values([-2.0, -2.0], location = 'centroids')
        a.set_boundary_values(1.0)
        operator.update_elliptic_matrix(a)

        assert num.allclose(A.todense()[0,:], -2*A0[0,:]-0.5*A1[0,:]-0.5*A2[0,:])
        assert num.allclose(A.todense()[1,:], -0.5*A0[1,:]-2*A1[1,:]-0.5*A2[1,:])

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

        operator = Kinematic_Viscosity(domain)

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

        operator = Kinematic_Viscosity(domain)

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


        u_out = operator.parabolic_solve(u_in, b, a, iprint=1)

        assert num.allclose(u_out.centroid_values, U_mod[:n])


################################################################################

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_Kinematic_Viscosity, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)

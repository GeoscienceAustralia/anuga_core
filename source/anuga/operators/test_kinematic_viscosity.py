from anuga import Domain
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
        boundary_enum = operator1.boundary_enum

        assert boundary_enum[(0,0)] == 0
        assert boundary_enum[(0,1)] == 1
        assert boundary_enum[(0,2)] == 2

        operator2 = self.operator2()
        boundary_enum = operator2.boundary_enum


        assert boundary_enum[(0,1)] == 0
        assert boundary_enum[(0,2)] == 1
        assert boundary_enum[(1,0)] == 2
        assert boundary_enum[(1,2)] == 3

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

    def test_elliptic_matrix(self):
        operator1 = self.operator1()
        operator2 = self.operator2()

        domain1 = operator1.domain
        domain2 = operator2.domain


        A = operator1.elliptic_operator_matrix

        print A
        print num.array([-6.0-12.0/sqrt(5), 6.0,  6.0/sqrt(5), 6.0/sqrt(5)])
        
        assert num.allclose(A.todense(), num.array([-6.0-12.0/sqrt(5), 6.0,  6.0/sqrt(5), 6.0/sqrt(5)]))


        A = operator1.elliptic_operator_matrix
        assert num.allclose(A.todense(), 1.5*num.array([-6.0-12.0/sqrt(5), 6.0, 6.0/sqrt(5), 6.0/sqrt(5)]))

        h = num.array([1.0, 1.0])
        A = operator2.apply_stage_heights(h)
        #From test_kv_system_geometry
        A0 = num.array([[-3.0,3.0,0.0,0.0,0.0,0.0],
                        [0.0,-6.0/sqrt(5.0),0.0,0.0,6.0/sqrt(5.0),0.0]])
        A1 = num.array([[-6.0/sqrt(5.0),0.0,6.0/sqrt(5.0),0.0,0.0,0.0],\
                        [3.0,-3.0,0.0,0.0,0.0,0.0]])
        A2 = num.array([[-6.0/sqrt(5.0),0.0,0.0,6.0/sqrt(5.0),0.0,0.0],\
                        [0.0, -6.0/sqrt(5.0), 0.0, 0.0, 0.0, 6.0/sqrt(5.0)]])
        assert num.allclose(A.todense(), A0+A1+A2)

        h = num.array([2.0, 1.0])
        A = operator2.apply_stage_heights(h)
        assert num.allclose(A.todense()[0,:], 1.5*A0[0,:]+1.5*A1[0,:]+1.5*A2[0,:])
        assert num.allclose(A.todense()[1,:], A0[1,:]+1.5*A1[1,:]+A2[1,:])

        h = num.array([-2.0, -2.0])
        A = operator2.apply_stage_heights(h)
        assert num.allclose(A.todense()[0,:], -2*A0[0,:]-0.5*A1[0,:]-0.5*A2[0,:])
        assert num.allclose(A.todense()[1,:], -0.5*A0[1,:]-2*A1[1,:]-0.5*A2[1,:])

    def test_elliptic_multiply(self):
        operator1 = self.operator1()
        operator1.apply_stage_heights(num.array([[1.0]])) #h=1
        operator1.build_boundary_vector()
        V1 = num.array([2.0]) #(uh)=2 <- Centriod value
        V2 = num.array([2.0]) #(vh)=2 <- Centroid value
        A = num.array([-6.0-12.0/sqrt(5), 6.0,  6.0/sqrt(5), 6.0/sqrt(5)])
        U1 = num.array([[2.0],[2.0],[1.0],[1.0]]) #(uh) for centroid, 3 edges
        U2 = num.array([[2.0],[1.0],[2.0],[0.0]]) #(vh) for centroid, 3 edge

        operator1.set_qty_considered('u')
        D1 = operator1.elliptic_multiply(V1)
        assert num.allclose(D1, 2*num.array(num.mat(A)*num.mat(U1)).reshape(1,)) #2* for triangle_areas

        operator1.set_qty_considered(2)
        D2 = operator1.elliptic_multiply(V2)
        assert num.allclose(D2, 2*num.array(num.mat(A)*num.mat(U2)).reshape(1,)) #2* for triangle_areas

    def test_mul(self):
        operator1 = self.operator1()
        operator1.apply_stage_heights(num.array([[1.0]])) #h=1
        operator1.set_qty_considered(1)
        operator1.build_boundary_vector()
        A = num.array([-6.0-12.0/sqrt(5), 6.0,  6.0/sqrt(5), 6.0/sqrt(5)])
        V1 = num.array([2.0]) #(uh)=2
        U1 = num.array([[2.0],[0.0],[0.0],[0.0]])
        assert num.allclose(operator1 * V1, 2*num.array(num.mat(A)*num.mat(U1)).reshape(1,))

    def test_cg_solve(self):
        #cf self.test_mul()
        operator1 = self.operator1()
        operator1.apply_stage_heights(num.array([[1.0]])) #h=1
        operator1.set_qty_considered('u')
        V = num.array([2.0]) #h=1, (uh)=2
        A = num.array([-6.0-12.0/sqrt(5), 6.0,  6.0/sqrt(5), 6.0/sqrt(5)])
        U = num.array([[2.0,2.0],[2.0,1.0],[1.0,2.0],[1.0,0.0]])
        test = 2*num.mat(A)*num.mat(U[:,0].reshape(4,1))
        X = operator1.cg_solve(num.array(test).reshape(1,))
        assert num.allclose(V, X)

    def test_parabolic_solve(self):
        operator1 = self.operator1()
        operator1.apply_stage_heights(num.array([[1.0]])) #h=1
        A = num.array([-6.0-12.0/sqrt(5), 6.0,  6.0/sqrt(5), 6.0/sqrt(5)])
        U = num.array([[2.0,1.0],[2.0,1.0],[1.0,2.0],[1.0,0.0]])
        u = num.array([[2.0,1.0]])
        U_new = operator1.parabolic_solver(u)
        U_mod = num.array([[0.0,0.0],[2.0,1.0],[1.0,2.0],[1.0,0.0]])
        U_mod[0,:] = U_new
        assert num.allclose(U_new - operator1.dt * 2 * num.mat(A)*num.mat(U_mod), U[0,:])

################################################################################

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_Kinematic_Viscosity, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)

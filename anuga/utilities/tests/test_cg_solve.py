#!/usr/bin/env python

import unittest
from anuga.utilities.sparse import Sparse, Sparse_CSR
from anuga.utilities.cg_solve import _conjugate_gradient
from anuga.utilities.cg_solve import *
import numpy as num
import os


class CGError(Exception):
    pass


class Test_CG_Solve(unittest.TestCase):

    def tearDown(self):
        try:
            os.remove('anuga.log')
        except OSError:
            pass

    def test_sparse_solve(self):
        """Solve Small Sparse Matrix"""

        A = [[2.0, -1.0, 0.0, 0.0],
             [-1.0, 2.0, -1.0, 0.0],
             [0.0, -1.0, 2.0, -1.0],
             [0.0, 0.0, -1.0, 2.0]]

        A = Sparse(A)

        xe = [0.0, 1.0, 2.0, 3.0]
        b = A*xe
        x = [0.0, 0.0, 0.0, 0.0]

        x = conjugate_gradient(A, b, x)

        assert num.allclose(x, xe)

    def test_max_iter(self):
        """Test max iteration Small Sparse Matrix"""

        A = [[2.0, -1.0, 0.0, 0.0],
             [-1.0, 2.0, -1.0, 0.0],
             [0.0, -1.0, 2.0, -1.0],
             [0.0, 0.0, -1.0, 2.0]]

        A = Sparse(A)

        xe = [0.0, 1.0, 2.0, 3.0]
        b = A*xe
        x = [0.0, 0.0, 0.0, 0.0]

        try:
            x = conjugate_gradient(A, b, x, imax=2)
        except ConvergenceError:
            pass
        else:
            msg = 'Should have raised exception'
            raise CGError(msg)

    def test_solve_large(self):
        """Standard 1d laplacian """

        n = 50
        A = Sparse(n, n)

        for i in num.arange(0, n):
            A[i, i] = 1.0
            if i > 0:
                A[i, i-1] = -0.5
            if i < n-1:
                A[i, i+1] = -0.5

        xe = num.ones((n,), float)

        b = A*xe
        x = conjugate_gradient(A, b, b, tol=1.0e-5)

        assert num.allclose(x, xe)

    def test_solve_large_2d(self):
        """Standard 2d laplacian"""

        n = 20
        m = 10

        A = Sparse(m*n, m*n)

        for i in num.arange(0, n):
            for j in num.arange(0, m):
                I = j+m*i
                A[I, I] = 4.0
                if i > 0:
                    A[I, I-m] = -1.0
                if i < n-1:
                    A[I, I+m] = -1.0
                if j > 0:
                    A[I, I-1] = -1.0
                if j < m-1:
                    A[I, I+1] = -1.0

        xe = num.ones((n*m,), float)

        b = A*xe
        x = conjugate_gradient(A, b, b, iprint=1)

        assert num.allclose(x, xe)

    def test_solve_large_2d_csr_matrix(self):
        """Standard 2d laplacian with csr format
        """

        n = 100
        m = 100

        A = Sparse(m*n, m*n)

        for i in num.arange(0, n):
            for j in num.arange(0, m):
                I = j+m*i
                A[I, I] = 4.0
                if i > 0:
                    A[I, I-m] = -1.0
                if i < n-1:
                    A[I, I+m] = -1.0
                if j > 0:
                    A[I, I-1] = -1.0
                if j < m-1:
                    A[I, I+1] = -1.0

        xe = num.ones((n*m,), float)

        # Convert to csr format
        # print 'start covert'
        A = Sparse_CSR(A)
        # print 'finish covert'
        b = A*xe
        x = conjugate_gradient(A, b, b, iprint=20)

        assert num.allclose(x, xe)

    def test_solve_large_2d_with_default_guess(self):
        """Standard 2d laplacian using default first guess"""

        n = 20
        m = 10

        A = Sparse(m*n, m*n)

        for i in num.arange(0, n):
            for j in num.arange(0, m):
                I = j+m*i
                A[I, I] = 4.0
                if i > 0:
                    A[I, I-m] = -1.0
                if i < n-1:
                    A[I, I+m] = -1.0
                if j > 0:
                    A[I, I-1] = -1.0
                if j < m-1:
                    A[I, I+1] = -1.0

        xe = num.ones((n*m,), float)

        b = A*xe
        x = conjugate_gradient(A, b)

        assert num.allclose(x, xe)

    def test_vector_shape_error(self):
        """Raise VectorShapeError"""

        A = [[2.0, -1.0, 0.0, 0.0],
             [-1.0, 2.0, -1.0, 0.0],
             [0.0, -1.0, 2.0, -1.0],
             [0.0, 0.0, -1.0, 2.0]]

        A = Sparse(A)

        xe = [[0.0, 2.0], [1.0, 3.0], [2.0, 4.0], [3.0, 2.0]]

        try:
            x = _conjugate_gradient(A, xe, xe, iprint=0)
        except VectorShapeError:
            pass
        else:
            msg = 'Should have raised exception'
            raise CGError(msg)

    def test_sparse_solve_matrix(self):
        """Solve Small Sparse Matrix"""

        A = [[2.0, -1.0, 0.0, 0.0],
             [-1.0, 2.0, -1.0, 0.0],
             [0.0, -1.0, 2.0, -1.0],
             [0.0, 0.0, -1.0, 2.0]]

        A = Sparse(A)

        xe = [[0.0, 0.0], [1.0, 1.0], [2.0, 2.0], [3.0, 3.0]]
        b = A*xe
        x = [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
        x = conjugate_gradient(A, b, x, iprint=0)

        assert num.allclose(x, xe)

    def test_sparse_solve_using_c_ext(self):
        """Solve Small Sparse Matrix"""

        A = [[2.0, -1.0, 0.0, 0.0],
             [-1.0, 2.0, -1.0, 0.0],
             [0.0, -1.0, 2.0, -1.0],
             [0.0, 0.0, -1.0, 2.0]]

        A = Sparse_CSR(Sparse(A))

        xe = [0.0, 1.0, 2.0, 3.0]
        b = A*xe
        x = [0.0, 0.0, 0.0, 0.0]

        x = conjugate_gradient(A, b, x, use_c_cg=True)

        assert num.allclose(x, xe)

    def test_max_iter_using_c_ext(self):
        """Test max iteration Small Sparse Matrix"""

        A = [[2.0, -1.0, 0.0, 0.0],
             [-1.0, 2.0, -1.0, 0.0],
             [0.0, -1.0, 2.0, -1.0],
             [0.0, 0.0, -1.0, 2.0]]

        A = Sparse_CSR(Sparse(A))

        xe = [0.0, 1.0, 2.0, 3.0]
        b = A*xe
        x = [0.0, 0.0, 0.0, 0.0]

        try:
            x = conjugate_gradient(A, b, x, imax=2, use_c_cg=True)
        except ConvergenceError:
            pass
        else:
            msg = 'Should have raised exception'
            raise CGError(msg)

    def test_solve_large_using_c_ext(self):
        """Standard 1d laplacian """

        n = 50
        A = Sparse(n, n)

        for i in num.arange(0, n):
            A[i, i] = 1.0
            if i > 0:
                A[i, i-1] = -0.5
            if i < n-1:
                A[i, i+1] = -0.5

        xe = num.ones((n,), float)

        b = A*xe

        A = Sparse_CSR(A)

        x = conjugate_gradient(A, b, b, tol=1.0e-5, use_c_cg=True)

        assert num.allclose(x, xe)

    def test_solve_large_2d_using_c_ext(self):
        """Standard 2d laplacian"""

        n = 20
        m = 10

        A = Sparse(m*n, m*n)

        for i in num.arange(0, n):
            for j in num.arange(0, m):
                I = j+m*i
                A[I, I] = 4.0
                if i > 0:
                    A[I, I-m] = -1.0
                if i < n-1:
                    A[I, I+m] = -1.0
                if j > 0:
                    A[I, I-1] = -1.0
                if j < m-1:
                    A[I, I+1] = -1.0

        xe = num.ones((n*m,), float)
        A = Sparse_CSR(A)
        b = A*xe
        x = conjugate_gradient(A, b, b, iprint=1, use_c_cg=True)

        assert num.allclose(x, xe)

    def test_solve_large_2d_csr_matrix_using_c_ext(self):
        """Standard 2d laplacian with csr format
        """

        n = 100
        m = 100

        A = Sparse(m*n, m*n)

        for i in num.arange(0, n):
            for j in num.arange(0, m):
                I = j+m*i
                A[I, I] = 4.0
                if i > 0:
                    A[I, I-m] = -1.0
                if i < n-1:
                    A[I, I+m] = -1.0
                if j > 0:
                    A[I, I-1] = -1.0
                if j < m-1:
                    A[I, I+1] = -1.0

        xe = num.ones((n*m,), float)

        # Convert to csr format
        # print 'start covert'
        A = Sparse_CSR(A)
        # print 'finish covert'
        b = A*xe
        x = conjugate_gradient(A, b, b, iprint=20, use_c_cg=True)

        assert num.allclose(x, xe)

    def test_solve_large_2d_with_default_guess_using_c_ext(self):
        """Standard 2d laplacian using default first guess"""

        n = 20
        m = 10

        A = Sparse(m*n, m*n)

        for i in num.arange(0, n):
            for j in num.arange(0, m):
                I = j+m*i
                A[I, I] = 4.0
                if i > 0:
                    A[I, I-m] = -1.0
                if i < n-1:
                    A[I, I+m] = -1.0
                if j > 0:
                    A[I, I-1] = -1.0
                if j < m-1:
                    A[I, I+1] = -1.0

        xe = num.ones((n*m,), float)
        A = Sparse_CSR(A)
        b = A*xe
        x = conjugate_gradient(A, b, use_c_cg=True)

        assert num.allclose(x, xe)

    def test_sparse_solve_matrix_using_c_ext(self):
        """Solve Small Sparse Matrix"""

        A = [[2.0, -1.0, 0.0, 0.0],
             [-1.0, 2.0, -1.0, 0.0],
             [0.0, -1.0, 2.0, -1.0],
             [0.0, 0.0, -1.0, 2.0]]

        A = Sparse_CSR(Sparse(A))

        xe = [[0.0, 0.0], [1.0, 1.0], [2.0, 2.0], [3.0, 3.0]]
        b = A*xe
        x = [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
        x = conjugate_gradient(A, b, x, iprint=0, use_c_cg=True)

        assert num.allclose(x, xe)

    def test_sparse_solve_using_c_ext_with_jacobi(self):
        """Solve Small Sparse Matrix"""

        A = [[2.0, -1.0, 0.0, 0.0],
             [-1.0, 2.0, -1.0, 0.0],
             [0.0, -1.0, 2.0, -1.0],
             [0.0, 0.0, -1.0, 2.0]]

        A = Sparse_CSR(Sparse(A))

        xe = [0.0, 1.0, 2.0, 3.0]
        b = A*xe
        x = [0.0, 0.0, 0.0, 0.0]
        x = conjugate_gradient(A, b, x, use_c_cg=True, precon='Jacobi')

        assert num.allclose(x, xe)

    def test_breakdown_on_singular_matrix_using_c_ext(self):
        """A singular matrix makes d'Ad vanish; the C solver must report
        it rather than divide by zero and return NaN."""

        A = [[1.0, 1.0],
             [1.0, 1.0]]
        A = Sparse_CSR(Sparse(A))
        b = num.array([1.0, -1.0])

        for precon in (None, 'Jacobi'):
            try:
                conjugate_gradient(A, b, use_c_cg=True, precon=precon)
            except ConvergenceError as e:
                assert 'broke down' in str(e), str(e)
            else:
                raise Exception('Should have raised ConvergenceError')

    def test_max_iter_using_c_ext_with_jacobi(self):
        """Test max iteration Small Sparse Matrix"""

        A = [[2.0, -1.0, 0.0, 0.0],
             [-1.0, 2.0, -1.0, 0.0],
             [0.0, -1.0, 2.0, -1.0],
             [0.0, 0.0, -1.0, 2.0]]

        A = Sparse_CSR(Sparse(A))

        xe = [0.0, 1.0, 2.0, 3.0]
        b = A*xe
        x = [0.0, 0.0, 0.0, 0.0]

        try:
            x = conjugate_gradient(
                A, b, x, imax=2, use_c_cg=True, precon='Jacobi')
        except ConvergenceError:
            pass
        else:
            msg = 'Should have raised exception'
            raise CGError(msg)

    def test_solve_large_using_c_ext_with_jacobi(self):
        """Standard 1d laplacian """

        n = 50
        A = Sparse(n, n)

        for i in num.arange(0, n):
            A[i, i] = 1.0
            if i > 0:
                A[i, i-1] = -0.5
            if i < n-1:
                A[i, i+1] = -0.5

        xe = num.ones((n,), float)

        b = A*xe

        A = Sparse_CSR(A)

        x = conjugate_gradient(A, b, b, tol=1.0e-5,
                               use_c_cg=True, precon='Jacobi')

        assert num.allclose(x, xe)

    def test_solve_large_2d_using_c_ext_with_jacobi(self):
        """Standard 2d laplacian"""

        n = 20
        m = 10

        A = Sparse(m*n, m*n)

        for i in num.arange(0, n):
            for j in num.arange(0, m):
                I = j+m*i
                A[I, I] = 4.0
                if i > 0:
                    A[I, I-m] = -1.0
                if i < n-1:
                    A[I, I+m] = -1.0
                if j > 0:
                    A[I, I-1] = -1.0
                if j < m-1:
                    A[I, I+1] = -1.0

        xe = num.ones((n*m,), float)
        A = Sparse_CSR(A)
        b = A*xe
        x = conjugate_gradient(
            A, b, b, iprint=1, use_c_cg=True, precon='Jacobi')

        assert num.allclose(x, xe)

    def test_solve_large_2d_csr_matrix_using_c_ext_with_jacobi(self):
        """Standard 2d laplacian with csr format
        """

        n = 100
        m = 100

        A = Sparse(m*n, m*n)

        for i in num.arange(0, n):
            for j in num.arange(0, m):
                I = j+m*i
                A[I, I] = 4.0
                if i > 0:
                    A[I, I-m] = -1.0
                if i < n-1:
                    A[I, I+m] = -1.0
                if j > 0:
                    A[I, I-1] = -1.0
                if j < m-1:
                    A[I, I+1] = -1.0

        xe = num.ones((n*m,), float)

        # Convert to csr format
        # print 'start covert'
        A = Sparse_CSR(A)
        # print 'finish covert'
        b = A*xe
        x = conjugate_gradient(
            A, b, b, iprint=20, use_c_cg=True, precon='Jacobi')

        assert num.allclose(x, xe)

    def test_solve_large_2d_with_default_guess_using_c_ext_with_jacobi(self):
        """Standard 2d laplacian using default first guess"""

        n = 20
        m = 10

        A = Sparse(m*n, m*n)

        for i in num.arange(0, n):
            for j in num.arange(0, m):
                I = j+m*i
                A[I, I] = 4.0
                if i > 0:
                    A[I, I-m] = -1.0
                if i < n-1:
                    A[I, I+m] = -1.0
                if j > 0:
                    A[I, I-1] = -1.0
                if j < m-1:
                    A[I, I+1] = -1.0

        xe = num.ones((n*m,), float)
        A = Sparse_CSR(A)
        b = A*xe
        x = conjugate_gradient(A, b, use_c_cg=True, precon='Jacobi')

        assert num.allclose(x, xe)

    def test_sparse_solve_matrix_using_c_ext_with_jacobi(self):
        """Solve Small Sparse Matrix"""

        A = [[2.0, -1.0, 0.0, 0.0],
             [-1.0, 2.0, -1.0, 0.0],
             [0.0, -1.0, 2.0, -1.0],
             [0.0, 0.0, -1.0, 2.0]]

        A = Sparse_CSR(Sparse(A))

        xe = [[0.0, 0.0], [1.0, 1.0], [2.0, 2.0], [3.0, 3.0]]
        b = A*xe
        x = [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
        x = conjugate_gradient(
            A, b, x, iprint=0, use_c_cg=True, precon='Jacobi')

        assert num.allclose(x, xe)

class Test_CG_Solve_extra(unittest.TestCase):
    """Cover previously uncovered lines in cg_solve.py."""

    def test_stats_str(self):
        """Stats.__str__ covers lines 35-37."""
        from anuga.utilities.cg_solve import Stats
        s = Stats()
        s.iter = 1; s.rTr = 0.5; s.x = 1.0; s.dx = 0.01; s.rTr0 = 1.0; s.x0 = 0.0
        result = str(s)
        self.assertIn('iter', result)

    def test_conjugate_gradient_output_stats(self):
        """conjugate_gradient with output_stats=True covers line 137."""
        A = Sparse([[4.0, 0.0], [0.0, 2.0]])
        b = num.array([8.0, 4.0])
        x0 = num.zeros(2)
        x, stats = conjugate_gradient(A, b, x0=x0, output_stats=True)
        assert num.allclose(x, [2.0, 2.0])

    def test_conjugate_gradient_x0_none(self):
        """_conjugate_gradient with x0=None covers line 169."""
        A = Sparse([[4.0, 0.0], [0.0, 2.0]])
        b = num.array([8.0, 4.0])
        x, stats = _conjugate_gradient(A, b, None)
        assert num.allclose(x, [2.0, 2.0])


################################################################################


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_CG_Solve)
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)

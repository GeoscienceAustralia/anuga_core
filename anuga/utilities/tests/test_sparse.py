#!/usr/bin/env python
import unittest
from math import sqrt

from anuga.utilities.sparse import *
import numpy as num


class Test_Sparse(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_init1(Self):
        """Test initialisation from dimensions
        """
        A = Sparse(3, 3)
        A[1, 1] = 4

        for i in range(3):
            for j in range(3):
                if i == 1 and j == 1:
                    assert A[i, j] == 4.0
                else:
                    assert A[i, j] == 0.0

    def test_init2(Self):
        """Test initialisation from dense matrix
        """

        A = Sparse(4, 3)
        A[1, 1] = 4
        A[2, 0] = 6
        A[0, 1] = 13
        A[0, 2] = -6
        A[2, 2] = 1

        B = A.todense()

        C = Sparse(B)

        assert num.allclose(C.todense(), B)

    def test_dense(self):
        A = Sparse(4, 3)
        A[1, 1] = 4

        assert num.allclose(
            A.todense(), [[0, 0, 0], [0, 4, 0], [0, 0, 0], [0, 0, 0]])

    def test_reset_to_zero_possible(self):
        """Test that nonzero element can be reset to zero
        (This was not possible when we tried sparse from scipy).
        """

        A = Sparse(3, 3)
        A[1, 1] = 4
        A[1, 1] = 0

        assert len(A) == 0
        assert num.allclose(A.todense(), [[0, 0, 0], [0, 0, 0], [0, 0, 0]])

        # Set an existing zero element to zero
        A[1, 2] = 0
        assert len(A) == 0
        assert num.allclose(A.todense(), [[0, 0, 0], [0, 0, 0], [0, 0, 0]])

    def test_sparse_multiplication_vector(self):
        A = Sparse(3, 3)

        A[0, 0] = 3
        A[1, 1] = 2
        A[1, 2] = 2
        A[2, 2] = 1

        # Right hand side vector
        v = [2, 3, 4]

        u = A*v
        assert num.allclose(u, [6, 14, 4])

        # Right hand side column
        v = num.array([[2, 4], [3, 4], [4, 4]])

        u = A*v[:, 0]
        assert num.allclose(u, [6, 14, 4])

        u = A*v[:, 1]
        assert num.allclose(u, [12, 16, 4])

    def test_sparse_multiplication_matrix(self):
        A = Sparse(3, 3)

        A[0, 0] = 3
        A[1, 1] = 2
        A[1, 2] = 2
        A[2, 2] = 1

        # Right hand side matrix
        v = num.array([[2, 4], [3, 4], [4, 4]])

        u = A*v
        assert num.allclose(u, [[6, 12], [14, 16], [4, 4]])

    def test_sparse_transpose_multiplication(self):
        A = Sparse(3, 3)

        A[0, 0] = 3
        A[1, 1] = 2
        A[1, 2] = 2
        A[2, 2] = 1

        # Right hand side vector
        v = [2, 3, 4]

        u = A.trans_mult(v)
        assert num.allclose(u, [6, 6, 10])

    def test_scalar_multiplication(self):
        """Test method __rmul__
        """

        A = Sparse(3, 3)

        A[0, 0] = 3
        A[1, 1] = 2
        A[1, 2] = 2
        A[2, 2] = 1

        B = 3*A
        assert num.allclose(B.todense(), 3*A.todense())

        B = A*3
        assert num.allclose(B.todense(), 3*A.todense())

        try:
            B = 'a'*A
        except TypeError:
            pass
        else:
            raise Exception('Should have failed')

    def test_sparse_addition(self):
        """ Test sparse addition with dok format
        """

        A = Sparse(3, 3)

        A[0, 0] = 3
        A[1, 1] = 2
        A[1, 2] = 2
        A[2, 2] = 1

        B = 3*A
        B[1, 0] = 2

        C = A+B

        assert num.allclose(C.todense(), [[12, 0, 0], [2, 8, 8], [0, 0, 4]])

    def test_sparse_tocsr(self):
        """ Test conversion to csr format
        """

        A = Sparse(4, 3)

        A[0, 0] = 3
        A[1, 1] = 2
        A[1, 2] = 2
        A[2, 2] = 1
        A[0, 2] = 4
        A[2, 0] = 5

        # print ' '
        # print A.todense()

        B = Sparse_CSR(A)

        # print B.todense()

        C = [1, 2, 3]

        assert num.allclose(B*C, [15.0, 10.0, 8.0, 0.0])

        C2 = [[1, 2], [2, 4], [3, 6]]

        # print B*C2

        assert num.allclose(
            B*C2, [[15.0, 30.0], [10.0, 20.0], [8.0, 16.0], [0.0, 0.0]])

    def test_sparse_csr_init(self):
        A = num.array([[1.0, 0.0, -1.0, 0.0], [0.0, 2.0, 0.0, 0.0],
                       [0.0, 0.0, 0.0, -3.0], [0.0, 0.0, 4.0, 0.0]])
        data = num.array([1.0, -1.0, 2.0, -3.0, 4.0])
        Colind = num.array([0, 2, 1, 3, 2])
        # the 5 does not correspond to any row, it is just there so we know how many nonzero entries there are!
        rowptr = num.array([0, 2, 3, 4, 5])
        A_CSR = Sparse_CSR(None, data, Colind, rowptr, 4, 4)
        A_dense = A_CSR.todense()
        assert num.allclose(A, A_dense)

class Test_Sparse_extra(unittest.TestCase):
    """Cover previously uncovered lines in sparse.py."""

    def test_sparse_from_bad_input_raises(self):
        """Sparse from non-convertable raises (lines 25-26)."""
        with self.assertRaises(Exception):
            Sparse('not_a_matrix')

    def test_sparse_no_args_raises(self):
        """Sparse() with wrong arg count raises (line 40)."""
        with self.assertRaises(Exception):
            Sparse(1, 2, 3)  # three args is invalid

    def test_sparse_repr(self):
        """Sparse.__repr__ covers line 45."""
        A = Sparse(2, 2)
        A[0, 1] = 3.0
        s = repr(A)
        self.assertIn('sparse matrix', s)

    def test_sparse_nonzeros(self):
        """Sparse.nonzeros covers line 55."""
        A = Sparse(3, 3)
        A[0, 1] = 1.0
        A[2, 2] = 2.0
        self.assertEqual(A.nonzeros(), 2)

    def test_sparse_mul_high_dim_raises(self):
        """Sparse.__mul__ with 3D array raises (line 146)."""
        A = Sparse(2, 2)
        A[0, 0] = 1.0
        with self.assertRaises(ValueError):
            A * num.zeros((2, 2, 2))

    def test_sparse_trans_mult_2d_raises(self):
        """trans_mult with 2D array raises (line 206)."""
        A = Sparse(2, 3)
        A[0, 1] = 1.0
        with self.assertRaises(Exception):
            A.trans_mult(num.zeros((2, 2)))

    def test_sparse_csr_bad_args_raises(self):
        """Sparse_CSR with non-Sparse A raises ValueError (line 284)."""
        # A is not None, not Sparse, and data/colind/rowptr don't match
        with self.assertRaises(ValueError):
            Sparse_CSR(A=42)  # A is not Sparse, no data arrays → else branch

    def test_sparse_csr_repr(self):
        """Sparse_CSR.__repr__ covers line 288."""
        data = num.array([1.0, 2.0])
        colind = num.array([0, 1])
        rowptr = num.array([0, 1, 2])
        csr = Sparse_CSR(None, data, colind, rowptr, 2, 2)
        s = repr(csr)
        self.assertIn('sparse matrix', s)

    def test_sparse_csr_len(self):
        """Sparse_CSR.__len__ covers line 294."""
        data = num.array([1.0, 2.0])
        colind = num.array([0, 1])
        rowptr = num.array([0, 1, 2])
        csr = Sparse_CSR(None, data, colind, rowptr, 2, 2)
        self.assertEqual(len(csr), 2)

    def test_sparse_csr_nonzeros(self):
        """Sparse_CSR.nonzeros covers line 299."""
        data = num.array([1.0])
        colind = num.array([0])
        rowptr = num.array([0, 1])
        csr = Sparse_CSR(None, data, colind, rowptr, 1, 2)
        self.assertEqual(csr.nonzeros(), 1)


################################################################################


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_Sparse)
    runner = unittest.TextTestRunner()
    runner.run(suite)

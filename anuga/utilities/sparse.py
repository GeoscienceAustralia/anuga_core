"""
Proof of concept sparse matrix code
"""

from .sparse_ext import csr_mv
import numpy as num


class Sparse:

    def __init__(self, *args):
        """Create sparse matrix.
        There are two construction forms
        Usage:

        Sparse(A)     #Creates sparse matrix from dense matrix A
        Sparse(M, N)  #Creates empty MxN sparse matrix
        """

        self.Data = {}

        if len(args) == 1:
            try:
                A = num.array(args[0])
            except:
                raise Exception('Input must be convertable to a numeric array')

            assert len(A.shape) == 2, 'Input must be a 2d matrix'

            self.M, self.N = A.shape
            for i in range(self.M):
                for j in range(self.N):
                    if A[i, j] != 0.0:
                        self.Data[i, j] = A[i, j]

        elif len(args) == 2:
            self.M = args[0]
            self.N = args[1]
        else:
            raise Exception('Invalid construction')

        self.shape = (self.M, self.N)

    def __repr__(self):
        return '%d X %d sparse matrix:\n' % (self.M, self.N) + repr(self.Data)

    def __len__(self):
        """Return number of nonzeros of A
        """
        return len(self.Data)

    def nonzeros(self):
        """Return number of nonzeros of A
        """
        return len(self)

    def __setitem__(self, key, x):

        i, j = key
        # removing these asserts will not speed things up
        assert 0 <= i < self.M
        assert 0 <= j < self.N

        if x != 0:
            self.Data[key] = float(x)
        else:
            if key in self.Data:
                del self.Data[key]

    def __getitem__(self, key):

        i, j = key
        # removing these asserts will not speed things up
        assert 0 <= i < self.M
        assert 0 <= j < self.N

        if key in self.Data:
            return self.Data[key]
        else:
            return 0.0

    def copy(self):
        # FIXME: Use the copy module instead
        new = Sparse(self.M, self.N)

        for key in list(self.Data.keys()):
            i, j = key

            new[i, j] = self.Data[i, j]

        return new

    def todense(self):
        D = num.zeros((self.M, self.N), float)

        for i in range(self.M):
            for j in range(self.N):
                if (i, j) in self.Data:
                    D[i, j] = self.Data[(i, j)]
        return D

    def __mul__(self, other):
        """Multiply this matrix onto 'other' which can either be
        a numeric vector, a numeric matrix or another sparse matrix.
        """

        try:
            B = num.array(other)
        except:
            msg = 'FIXME: Only numeric types implemented so far'
            raise Exception(msg)

        # Assume numeric types from now on

        if len(B.shape) == 0:
            # Scalar - make sure we use __rmul__ method
            R = self.__rmul__(B)

        elif len(B.shape) == 1:
            # Vector
            msg = 'Mismatching dimensions: You cannot multiply (%d x %d) matrix onto %d-vector'\
                  % (self.M, self.N, B.shape[0])
            assert B.shape[0] == self.N, msg

            R = num.zeros(self.M, float)  # Result

            # Multiply nonzero elements
            for key in list(self.Data.keys()):
                i, j = key

                R[i] += self.Data[key]*B[j]
        elif len(B.shape) == 2:

            R = num.zeros((self.M, B.shape[1]), float)  # Result matrix

            # Multiply nonzero elements
            for col in range(R.shape[1]):
                # For each column

                for key in list(self.Data.keys()):
                    i, j = key

                    R[i, col] += self.Data[key]*B[j, col]

        else:
            raise ValueError('Dimension too high: d=%d' % len(B.shape))

        return R

    def __add__(self, other):
        """Add this matrix onto 'other'
        """

        new = other.copy()
        for key in list(self.Data.keys()):
            i, j = key

            new[i, j] += self.Data[key]

        return new

    def __rmul__(self, other):
        """Right multiply this matrix with scalar
        """

        try:
            other = float(other)
        except:
            msg = 'Sparse matrix can only "right-multiply" onto a scalar'
            raise TypeError(msg)
        else:
            new = self.copy()
            # Multiply nonzero elements
            for key in list(new.Data.keys()):
                i, j = key

                new.Data[key] = other*new.Data[key]

        return new

    def trans_mult(self, other):
        """Multiply the transpose of matrix with 'other' which can be
        a numeric vector.
        """

        try:
            B = num.array(other)
        except:
            print('FIXME: Only numeric types implemented so far')

        # Assume numeric types from now on
        if len(B.shape) == 1:
            # Vector

            assert B.shape[0] == self.M, 'Mismatching dimensions'

            R = num.zeros((self.N,), float)  # Result

            # Multiply nonzero elements
            for key in list(self.Data.keys()):
                i, j = key

                R[j] += self.Data[key]*B[i]

        else:
            raise Exception('Can only multiply with 1d array')

        return R


class Sparse_CSR(object):

    def __init__(self, A=None, data=None, Colind=None, rowptr=None, m=None, n=None):
        """Create sparse matrix in csr format.

        Sparse_CSR(A) #creates csr sparse matrix from sparse matrix
        Matrices are not built using this format, since it's painful to
        add values to an existing sparse_CSR instance (hence there are no
        methods to do this.)

        Rather, build a matrix, and convert it to this format for a speed
        increase.

        data - a 1D array of the data
        Colind - The ith item in this 1D array is the column index of the
                 ith data in the data array
        rowptr - 1D array, with the index representing the row of the matrix.
                 The item in the row represents the index into colind of the
                 first data value of this row.
                 Regard it as a pointer into the colind array, for the ith row.


        """

        if A is None:
            m = int(m)
            n = int(n)

        if isinstance(A, Sparse):

            keys = list(A.Data.keys())
            keys.sort()
            nnz = len(keys)
            data = num.zeros((nnz,), float)
            colind = num.zeros((nnz,), int)
            row_ptr = num.zeros((A.M+1,), int)
            current_row = -1
            k = 0
            for key in keys:
                ikey0 = int(key[0])
                ikey1 = int(key[1])
                if ikey0 != current_row:
                    current_row = ikey0
                    row_ptr[ikey0] = k
                data[k] = A.Data[key]
                colind[k] = ikey1
                k += 1
            for row in range(current_row+1, A.M+1):
                row_ptr[row] = nnz
            #row_ptr[-1] = nnz

            self.data = data
            self.colind = colind
            self.row_ptr = row_ptr
            self.M = A.M
            self.N = A.N
        elif isinstance(data, num.ndarray) and isinstance(Colind, num.ndarray) and isinstance(rowptr, num.ndarray) and isinstance(m, int) and isinstance(n, int):
            msg = "Sparse_CSR: data is array of wrong dimensions"
            #assert len(data.shape) == 1, msg
            nnz = data.size

            msg = "Sparse_CSR: Colind is array of wrong dimensions"
            assert Colind.shape == (nnz,), msg

            msg = "Sparse_CSR: rowptr is array of wrong dimensions"
            assert rowptr.shape == (m+1,), msg

            self.data = data
            self.colind = Colind
            self.row_ptr = rowptr
            self.M = m
            self.N = n
        else:
            raise ValueError(
                'Sparse_CSR(A) expects A == Sparse Matrix *or* data==array,colind==array,rowptr==array,m==int,n==int')

    def __repr__(self):
        return '%d X %d sparse matrix:\n' % (self.M, self.N) + 'data ' + repr(self.data) + '\ncolind ' + \
            repr(self.colind) + '\nrow_ptr ' + repr(self.row_ptr)

    def __len__(self):
        """Return number of nonzeros of A
        """
        return self.row_ptr[-1]

    def nonzeros(self):
        """Return number of nonzeros of A
        """
        return len(self)

    def todense(self):
        D = num.zeros((self.M, self.N), float)

        for i in range(self.M):
            for ckey in range(self.row_ptr[i], self.row_ptr[i+1]):
                j = self.colind[ckey]
                D[i, j] = self.data[ckey]
        return D

    def __mul__(self, other):
        """Multiply this matrix onto 'other' which can either be
        a numeric vector, a numeric matrix or another sparse matrix.
        """

        try:
            B = num.array(other)
        except:
            print('FIXME: Only numeric types implemented so far')

        return csr_mv(self, B)


# Setup for C extensions


if __name__ == '__main__':
    # A little selftest

    A = Sparse(3, 3)

    A[1, 1] = 4

    print(A)
    print(A.todense())

    A[1, 1] = 0

    print(A)
    print(A.todense())

    A[1, 2] = 0

    A[0, 0] = 3
    A[1, 1] = 2
    A[1, 2] = 2
    A[2, 2] = 1

    print(A)
    print(A.todense())

    # Right hand side vector
    v = [2, 3, 4]

    u = A*v
    print(u)
    assert num.allclose(u, [6, 14, 4])

    u = A.trans_mult(v)
    print(u)
    assert num.allclose(u, [6, 6, 10])

    # Right hand side column
    v = num.array([[2, 4], [3, 4], [4, 4]])

    u = A*v[:, 0]
    assert num.allclose(u, [6, 14, 4])

    #u = A*v[:,1]
    # print u
    print(A.shape)

    B = 3*A
    print(B.todense())

    B[1, 0] = 2

    C = A+B

    print(C.todense())

    C = Sparse_CSR(C)

    y = C*[6, 14, 4]

    print(y)

    y2 = C*[[6, 4], [4, 28], [4, 8]]

    print(y2)

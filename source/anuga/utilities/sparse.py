"""Proof of concept sparse matrix code
"""


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
            from Numeric import array
            try:
                A = array(args[0])
            except:
                raise 'Input must be convertable to a Numeric array'

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
            raise 'Invalid construction'
            
        self.shape = (self.M, self.N) 


    def __repr__(self):
        return '%d X %d sparse matrix:\n' %(self.M, self.N) + `self.Data`

    def __len__(self):
        """Return number of nonzeros of A
        """
        return len(self.Data)

    def nonzeros(self):
        """Return number of nonzeros of A
        """        
        return len(self)
    
    def __setitem__(self, key, x):

        i,j = key
        assert 0 <= i < self.M
        assert 0 <= j < self.N        

        if x != 0:
            self.Data[key] = float(x)
        else:
            if self.Data.has_key( key ):            
                del self.Data[key]

    def __getitem__(self, key):
        
        i,j = key
        assert 0 <= i < self.M
        assert 0 <= j < self.N                

        if self.Data.has_key( key ):
            return self.Data[ key ]
        else:
            return 0.0

    def copy(self):
        #FIXME: Use the copy module instead
        new = Sparse(self.M,self.N)

        for key in self.Data.keys():
            i, j = key

            new[i,j] = self.Data[i,j]

        return new


    def todense(self):
        from Numeric import zeros, Float

        D = zeros( (self.M, self.N), Float)
        
        for i in range(self.M):
            for j in range(self.N):
                if self.Data.has_key( (i,j) ):                
                    D[i, j] = self.Data[ (i,j) ]
        return D


    
    def __mul__(self, other):
        """Multiply this matrix onto 'other' which can either be
        a Numeric vector, a Numeric matrix or another sparse matrix.
        """

        from Numeric import array, zeros, Float
        
        try:
            B = array(other)
        except:
            msg = 'FIXME: Only Numeric types implemented so far'
            raise msg
            

        #Assume numeric types from now on
	
        if len(B.shape) == 0:
            #Scalar - use __rmul__ method
            R = B*self
	    
        elif len(B.shape) == 1:
            #Vector
            assert B.shape[0] == self.N, 'Mismatching dimensions'

            R = zeros(self.M, Float) #Result
	    
            #Multiply nonzero elements
            for key in self.Data.keys():
                i, j = key

                R[i] += self.Data[key]*B[j]
        elif len(B.shape) == 2:
	
            
            R = zeros((self.M, B.shape[1]), Float) #Result matrix

            #Multiply nonzero elements
	    for col in range(R.shape[1]):
	        #For each column
		
                for key in self.Data.keys():
                    i, j = key

                    R[i, col] += self.Data[key]*B[j, col]
	    
	    
        else:
            raise ValueError, 'Dimension too high: d=%d' %len(B.shape)

        return R
    

    def __add__(self, other):
        """Add this matrix onto 'other' 
        """

        from Numeric import array, zeros, Float
        
        new = other.copy()
        for key in self.Data.keys():
            i, j = key

            new[i,j] += self.Data[key]

        return new


    def __rmul__(self, other):
        """Right multiply this matrix with scalar
        """

        from Numeric import array, zeros, Float
        
        try:
            other = float(other)
        except:
            msg = 'Sparse matrix can only "right-multiply" onto a scalar'
            raise TypeError, msg
        else:
            new = self.copy()
            #Multiply nonzero elements
            for key in new.Data.keys():
                i, j = key

                new.Data[key] = other*new.Data[key]

        return new


    def trans_mult(self, other):
        """Multiply the transpose of matrix with 'other' which can be
        a Numeric vector.
        """

        from Numeric import array, zeros, Float
        
        try:
            B = array(other)
        except:
            print 'FIXME: Only Numeric types implemented so far'


        #Assume numeric types from now on
        if len(B.shape) == 1:
            #Vector

            assert B.shape[0] == self.M, 'Mismatching dimensions'

            R = zeros((self.N,), Float) #Result

            #Multiply nonzero elements
            for key in self.Data.keys():
                i, j = key

                R[j] += self.Data[key]*B[i]

        else:
            raise 'Can only multiply with 1d array'

        return R

class Sparse_CSR:

    def __init__(self, A):
        """Create sparse matrix in csr format.

        Sparse_CSR(A) #creates csr sparse matrix from sparse matrix
        Matrices are not built using this format, since it's painful to
        add values to an existing sparse_CSR instance (hence there are no
        objects to do this.)

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

        from Numeric import array, Float, Int

        if isinstance(A,Sparse):

            from Numeric import zeros
            keys = A.Data.keys()
            keys.sort()
            nnz = len(keys)
            data    = zeros ( (nnz,), Float)
            colind  = zeros ( (nnz,), Int)
            row_ptr = zeros ( (A.M+1,), Int)
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
        
            self.data    = data
            self.colind  = colind
            self.row_ptr = row_ptr
            self.M       = A.M
            self.N       = A.N
        else:
            raise ValueError, "Sparse_CSR(A) expects A == Sparse Matrix"
            
    def __repr__(self):
        return '%d X %d sparse matrix:\n' %(self.M, self.N) + `self.data`

    def __len__(self):
        """Return number of nonzeros of A
        """
        return self.row_ptr[-1]

    def nonzeros(self):
        """Return number of nonzeros of A
        """        
        return len(self)

    def todense(self):
        from Numeric import zeros, Float

        D = zeros( (self.M, self.N), Float)
        
        for i in range(self.M):
            for ckey in range(self.row_ptr[i],self.row_ptr[i+1]):
                j = self.colind[ckey]
                D[i, j] = self.data[ckey]
        return D

    def __mul__(self, other):
        """Multiply this matrix onto 'other' which can either be
        a Numeric vector, a Numeric matrix or another sparse matrix.
        """

        from Numeric import array, zeros, Float
        
        try:
            B = array(other)
        except:
            print 'FIXME: Only Numeric types implemented so far'

        return csr_mv(self,B) 


# Setup for C extensions
from anuga.utilities import compile
if compile.can_use_C_extension('sparse_ext.c'):
    # Access underlying c implementations
    from sparse_ext import csr_mv


if __name__ == '__main__':
    # A little selftest
    
    from Numeric import allclose, array, Float 
    
    A = Sparse(3,3)

    A[1,1] = 4


    print A
    print A.todense()

    A[1,1] = 0

    print A
    print A.todense()    

    A[1,2] = 0


    A[0,0] = 3
    A[1,1] = 2
    A[1,2] = 2
    A[2,2] = 1

    print A
    print A.todense()


    #Right hand side vector
    v = [2,3,4]

    u = A*v
    print u
    assert allclose(u, [6,14,4])

    u = A.trans_mult(v)
    print u
    assert allclose(u, [6,6,10])

    #Right hand side column
    v = array([[2,4],[3,4],[4,4]])

    u = A*v[:,0]
    assert allclose(u, [6,14,4])

    #u = A*v[:,1]
    #print u
    print A.shape

    B = 3*A
    print B.todense()

    B[1,0] = 2

    C = A+B

    print C.todense()

    C = Sparse_CSR(C)

    y = C*[6,14,4]

    print y

    y2 = C*[[6,4],[4,28],[4,8]]

    print y2

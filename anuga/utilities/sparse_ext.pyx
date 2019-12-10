#cython: wraparound=False, boundscheck=False, cdivision=True, profile=False, nonecheck=False, overflowcheck=False, cdivision_warnings=False, unraisable_tracebacks=False
import cython
# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np
# declare the interface to the C code
cdef extern from "sparse.c":
	int _csr_mv(int M, double* data, long* colind, long* row_ptr, double* x, double* y)
	int _csr_mm(int M, int columns, double* data, long* colind, long* row_ptr, double* x, double* y)

def csr_mv(object csr_sparse, np.ndarray x not None):

	cdef np.ndarray[double, ndim=1, mode="c"] data
	cdef np.ndarray[long, ndim=1, mode="c"] colind
	cdef np.ndarray[long, ndim=1, mode="c"] row_ptr
	cdef np.ndarray y
	cdef int M, err, columns, rows

	data = csr_sparse.data
	colind = csr_sparse.colind
	row_ptr = csr_sparse.row_ptr
	x = np.ascontiguousarray(x.astype(np.float))

	M = row_ptr.shape[0] - 1

	if x.ndim == 1: # Multiplicant is a vector

		y = np.ascontiguousarray(np.zeros(M, dtype=np.float))

		err = _csr_mv(M, &data[0], &colind[0], &row_ptr[0], <double* > x.data, <double* > y.data)

	elif x.ndim == 2:

		rows = x.shape[0]
		columns = x.shape[1]

		y = np.ascontiguousarray(np.zeros((M,columns),dtype=np.float))

		err = _csr_mm(M, columns, &data[0], &colind[0], &row_ptr[0], <double* > x.data, <double* > y.data)

	else:

		raise ValueError("Allowed dimensions in sparse_ext restricted to 1 or 2")
		return None

	return y




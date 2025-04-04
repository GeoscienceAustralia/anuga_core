#cython: wraparound=False, boundscheck=False, cdivision=True, profile=False, nonecheck=False, overflowcheck=False, cdivision_warnings=False, unraisable_tracebacks=False
import cython
from libc.stdint cimport int64_t

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

# declare the interface to the C code
cdef extern from "calc_grid_values.c":
	void _calc_grid_values(double* x, double* y, double* norms, int64_t num_vert, int64_t* volumes, int64_t num_tri, double cell_size, int64_t nrow, int64_t ncol, double* vertex_val, double* grid_val)
	void init_norms(double* x, double* y, double* norms, int64_t* volumes, int64_t num_tri)


def calc_grid_values(int64_t nrow, int64_t ncol,\
                    double cell_size,\
                    double nodata_val,\
                    np.ndarray[double, ndim=1, mode="c"] x not None,\
                    np.ndarray[double, ndim=1, mode="c"] y not None,\
                    np.ndarray[double, ndim=1, mode="c"] norms not None,\
                    np.ndarray[int64_t, ndim=2, mode="c"] volumes not None,\
                    np.ndarray[double, ndim=1, mode="c"] result not None,\
                    np.ndarray[double, ndim=1, mode="c"] grid_val not None):

	cdef int64_t i, num_tri, num_vert, num_norms, num_grid_val

	num_tri = volumes.shape[0]
	num_vert = x.shape[0]
	num_norms = norms.shape[0]
	num_grid_val = grid_val.shape[0]

	init_norms(&x[0], &y[0], &norms[0], &volumes[0,0], num_tri)

	for i in xrange(nrow*ncol):
		grid_val[i] = nodata_val

	_calc_grid_values(&x[0], &y[0],\
                    &norms[0],\
                    num_vert,\
                    &volumes[0,0],\
                    num_tri,\
                    cell_size,\
                    nrow,\
                    ncol,\
                    &result[0],\
                    &grid_val[0])

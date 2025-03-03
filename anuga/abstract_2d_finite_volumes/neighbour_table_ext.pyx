#cython: wraparound=False, boundscheck=False, cdivision=True, profile=False, nonecheck=False, overflowcheck=False, cdivision_warnings=False, unraisable_tracebacks=False
import cython
from libc.stdint cimport int64_t

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

ctypedef int64_t keyint

# declare the interface to the C code
cdef extern from "neighbour_table.cpp":
	int64_t _build_neighbour_structure(keyint N, keyint M, int64_t* triangles, int64_t* neighbours, int64_t* neighbour_edges, int64_t* number_of_boundaries)

def build_neighbour_structure(keyint N,\
						np.ndarray[int64_t, ndim=2, mode="c"] triangles not None,\
						np.ndarray[int64_t, ndim=2, mode="c"] neighbours not None,\
						np.ndarray[int64_t, ndim=2, mode="c"] neighbour_edges not None,\
						np.ndarray[int64_t, ndim=1, mode="c"] number_of_boundaries not None):

	cdef keyint M
	cdef int64_t err

	M = triangles.shape[0]

	err = _build_neighbour_structure(N, M, &triangles[0,0], &neighbours[0,0], &neighbour_edges[0,0], &number_of_boundaries[0])

	assert err == 0, "Duplicate Edge"




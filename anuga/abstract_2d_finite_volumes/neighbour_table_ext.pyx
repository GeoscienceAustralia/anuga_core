#cython: wraparound=False, boundscheck=False, cdivision=True, profile=False, nonecheck=False, overflowcheck=False, cdivision_warnings=False, unraisable_tracebacks=False
import cython

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

ctypedef long keyint

# declare the interface to the C code
cdef extern from "neighbour_table.cpp":
	int _build_neighbour_structure(keyint N, keyint M, long* triangles, long* neighbours, long* neighbour_edges, long* number_of_boundaries)

def build_neighbour_structure(keyint N,\
						np.ndarray[long, ndim=2, mode="c"] triangles not None,\
						np.ndarray[long, ndim=2, mode="c"] neighbours not None,\
						np.ndarray[long, ndim=2, mode="c"] neighbour_edges not None,\
						np.ndarray[long, ndim=1, mode="c"] number_of_boundaries not None):

	cdef keyint M
	cdef int err

	M = triangles.shape[0]

	err = _build_neighbour_structure(N, M, &triangles[0,0], &neighbours[0,0], &neighbour_edges[0,0], &number_of_boundaries[0])

	assert err == 0, "Duplicate Edge"




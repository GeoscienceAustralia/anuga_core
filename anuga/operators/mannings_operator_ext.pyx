#cython: wraparound=False, boundscheck=False, cdivision=True, profile=False, nonecheck=False, overflowcheck=False, cdivision_warnings=False, unraisable_tracebacks=False
import cython
from libc.stdint cimport int64_t

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

# declare the interface to the C code
cdef extern from "mannings_operator.c":
	void _manning_friction_sloped(double g, double eps, int64_t N, double* x, double* w, double* zv, double* uh, double* vh, double* eta, double* xmom_update, double* ymom_update)
	void _manning_friction_flat(double g, double eps, int64_t N, double* w, double* zv, double* uh, double* vh, double* eta, double* xmom, double* ymom)
	void _chezy_friction(double g, double eps, int64_t N, double* x, double* w, double* zv, double* uh, double* vh, double* chezy, double* xmom_update, double* ymom_update)

def manning_friction_flat(double g,\
						double eps,\
						np.ndarray[double, ndim=1, mode="c"] w not None,\
						np.ndarray[double, ndim=1, mode="c"] uh not None,\
						np.ndarray[double, ndim=1, mode="c"] vh not None,\
						np.ndarray[double, ndim=2, mode="c"] z not None,\
						np.ndarray[double, ndim=1, mode="c"] eta not None,\
						np.ndarray[double, ndim=1, mode="c"] xmom not None,\
						np.ndarray[double, ndim=1, mode="c"] ymom not None):

	cdef int64_t N

	N = w.shape[0]

	_manning_friction_flat(g, eps, N,\
						&w[0],\
						&z[0,0],\
						&uh[0],\
						&vh[0],\
						&eta[0],\
						&xmom[0],\
						&ymom[0])

def manning_friction_sloped(double g,\
							double eps,\
							np.ndarray[double, ndim=2, mode="c"] x not None,\
							np.ndarray[double, ndim=1, mode="c"] w not None,\
							np.ndarray[double, ndim=1, mode="c"] uh not None,\
							np.ndarray[double, ndim=1, mode="c"] vh not None,\
							np.ndarray[double, ndim=2, mode="c"] z not None,\
							np.ndarray[double, ndim=1, mode="c"] eta not None,\
							np.ndarray[double, ndim=1, mode="c"] xmom not None,\
							np.ndarray[double, ndim=1, mode="c"] ymom not None):

	cdef int64_t N

	N = w.shape[0]

	_manning_friction_sloped(g, eps, N,\
							&x[0,0],\
							&w[0],\
							&z[0,0],\
							&uh[0],\
							&vh[0],\
							&eta[0],\
							&xmom[0],\
							&ymom[0])

def chezy_friction(double g,\
							double eps,\
							np.ndarray[double, ndim=2, mode="c"] x not None,\
							np.ndarray[double, ndim=1, mode="c"] w not None,\
							np.ndarray[double, ndim=1, mode="c"] uh not None,\
							np.ndarray[double, ndim=1, mode="c"] vh not None,\
							np.ndarray[double, ndim=2, mode="c"] z not None,\
							np.ndarray[double, ndim=1, mode="c"] chezy not None,\
							np.ndarray[double, ndim=1, mode="c"] xmom not None,\
							np.ndarray[double, ndim=1, mode="c"] ymom not None):

	cdef int64_t N

	N = w.shape[0]

	_manning_friction_sloped(g, eps, N,\
							&x[0,0],\
							&w[0],\
							&z[0,0],\
							&uh[0],\
							&vh[0],\
							&chezy[0],\
							&xmom[0],\
							&ymom[0])

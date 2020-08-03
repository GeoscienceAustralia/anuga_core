#cythonoff: wraparound=False, boundscheck=False, cdivision=True, profile=False, nonecheck=False, overflowcheck=False, cdivision_warnings=False, unraisable_tracebacks=False
import cython
from libc.stdlib cimport malloc, free
from cpython.pycapsule cimport *
# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np
# declare the interface to the C code
cdef extern from "quad_tree.c":
	ctypedef struct quad_tree:
		pass
	void delete_quad_tree(quad_tree* quadtree)

cdef delete_quad_tree_cap(object cap):
	kill = <quad_tree* > PyCapsule_GetPointer(cap, "quad tree")
	if kill != NULL:
		delete_quad_tree(kill)
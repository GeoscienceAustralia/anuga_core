#cython: wraparound=False, boundscheck=False, cdivision=True, profile=False, nonecheck=False, overflowcheck=False, cdivision_warnings=False, unraisable_tracebacks=False
import cython
from libc.stdint cimport int64_t

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

def rectangular_cross_construct(np.ndarray[double, ndim=1, mode="c"] params not None,\
								np.ndarray[double, ndim=1, mode="c"] origin not None,\
								np.ndarray[double, ndim=2, mode="c"] points not None,\
								np.ndarray[int64_t, ndim=2, mode="c"] elements not None):


	cdef int64_t m, n, i, j, v1, v2 ,v3 ,v4, v5
	cdef int64_t numPoints, numElements
	cdef double len1, len2, delta1, delta2, x, y

	m = int(params[0])
	n = int(params[1])
	len1 = params[2]
	len2 = params[3]

	cdef np.ndarray[int64_t, ndim=2, mode="c"] vertices = np.ascontiguousarray(np.zeros((m+1,n+1),dtype=np.int64))

	delta1 = len1/m
	delta2 = len2/n

	numPoints = 0
	for i in xrange(m+1):
		for j in xrange(n+1):
			vertices[i,j] = numPoints
			points[numPoints,0] = i*delta1 + origin[0]
			points[numPoints,1] = j*delta2 + origin[1]
			numPoints += 1

	boundary = {}
	numElements = 0
	for i in xrange(m):
		for j in xrange(n):
			v1 = vertices[i,j+1]
			v2 = vertices[i,j]
			v3 = vertices[i+1,j+1]
			v4 = vertices[i+1,j]
			x = (points[v1,0] + points[v2,0] + points[v3,0] + points[v4,0])*0.25
			y = (points[v1,1] + points[v2,1] + points[v3,1] + points[v4,1])*0.25

			# Create centre point
			v5 = numPoints
			points[numPoints,0] = x
			points[numPoints,1] = y
			numPoints += 1

			# Create left triangle
			if i == 0:
				boundary[(numElements,1)] = "left"

			elements[numElements,0] = v2
			elements[numElements,1] = v5
			elements[numElements,2] = v1
			numElements += 1

			# Create bottom triangle
			if j == 0:
				boundary[(numElements,1)] = "bottom"

			elements[numElements,0] = v4
			elements[numElements,1] = v5
			elements[numElements,2] = v2
			numElements += 1

			# Create right triangle
			if i == m-1:
				boundary[(numElements,1)] = "right"

			elements[numElements,0] = v3
			elements[numElements,1] = v5
			elements[numElements,2] = v4
			numElements += 1

			# Create top triangle
			if j == n-1:
				boundary[(numElements,1)] = "top"

			elements[numElements,0] = v1
			elements[numElements,1] = v5
			elements[numElements,2] = v3
			numElements += 1

	return boundary

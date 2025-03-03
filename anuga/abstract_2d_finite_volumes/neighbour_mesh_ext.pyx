#cython: wraparound=False, boundscheck=False, cdivision=True, profile=False, nonecheck=False, overflowcheck=False, cdivision_warnings=False, unraisable_tracebacks=False
import cython
from libc.stdint cimport int64_t

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

def boundary_dictionary_construct(int64_t numTriangle, defaultTag,\
                                np.ndarray[int64_t, ndim=2, mode="c"] neighbours not None,\
                                dict boundary):

	cdef int64_t a, b, vol_id, edge_id

	#defaultTag = defaultTag.encode('utm-f')

	a = neighbours.shape[0]
	b = neighbours.shape[1]

	if bool(boundary):
		for vol_id, edge_id in boundary.keys():
			msg = 'Segment (%d, %d) does not exist' %(vol_id, edge_id)
			assert vol_id < a and edge_id < b, msg

	for vol_id in xrange(numTriangle):
		for edge_id in xrange(0,3):
			if neighbours[vol_id, edge_id] < 0:
				if not boundary.has_key((vol_id, edge_id)):
					boundary[(vol_id, edge_id)] = defaultTag

	return boundary

def check_integrity_c(np.ndarray[int64_t, ndim=1, mode="c"] vertex_value_indices not None,\
					np.ndarray[int64_t, ndim=2, mode="c"] triangles not None,\
					np.ndarray[int64_t, ndim=1, mode="c"] node_index not None,\
					np.ndarray[int64_t, ndim=1, mode="c"] number_of_triangles_per_node not None):

	cdef int64_t nt, nt3, tri, n_node, n_node_1
	cdef int64_t current_node, k, i, index

	cdef int64_t cumsum

	nt3 = vertex_value_indices.shape[0]

	nt = triangles.shape[0]
	tri = triangles.shape[1]

	n_node_1 = node_index.shape[0]

	n_node = number_of_triangles_per_node.shape[0]

	assert nt3 == 3*nt, "Mismatch in size of triangles and vertex_value_indices"

	assert n_node_1 == n_node + 1, "Mismatch in size of node_index and number_of_triangles_per_node"

	assert tri == 3, "Triangle array should be ny by 3"

	current_node = 0
	k = 0

	for i in xrange(nt3):
		index = vertex_value_indices[i]

		if number_of_triangles_per_node[current_node] == 0:
			continue

		k += 1

		assert triangles[index/tri, index%tri] == current_node, "Inconsistency between triangles and vertex_values_indices"

		if number_of_triangles_per_node[current_node] == k:
			k = 0
			current_node += 1

	cumsum = 0

	for i in xrange(n_node):
		cumsum += number_of_triangles_per_node[i]
		assert cumsum == node_index[i+1], "Inconsistency between node_index and number_of_triangles_per_node"

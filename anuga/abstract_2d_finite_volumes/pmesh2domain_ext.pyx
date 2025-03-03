#cython: wraparound=False, boundscheck=False, cdivision=True, profile=False, nonecheck=False, overflowcheck=False, cdivision_warnings=False, unraisable_tracebacks=False
import cython
from cython cimport typeof
from libc.stdint cimport int64_t

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

# declare the interface to the C code
cdef extern from "pmesh2domain.c":
	# If the structs are defined using typedef in the C code, we have to use ctypedef here.
	# Otherwise, it will lead to weird errors.
	ctypedef struct UT_hash_handle:
		pass
	ctypedef struct segment_key_t:
		int64_t i
		int64_t j
	ctypedef struct segment_t:
		segment_key_t key
		int64_t vol_id
		int64_t edge_id
		UT_hash_handle hh
	void add_segment(segment_key_t key, int64_t vol_id, int64_t edge_id)
	segment_t* find_segment(segment_key_t key)
	void delete_segment(segment_t* segment)
	void delete_segment_all()
	void print_segments()
	int64_t vol_id_sort(segment_t* a, segment_t* b)
	int64_t key_sort(segment_t* a, segment_t* b)
	void sort_by_vol_id()
	void sort_by_key()

def build_boundary_dictionary(np.ndarray[int64_t, ndim=2, mode="c"] triangles not None,\
							np.ndarray[int64_t, ndim=2, mode="c"] segments not None,\
							list segment_tags,\
							dict tag_dict):
	"""
	Update neighbours array using triangle

	N is number of nodes (vertices)
	triangle nodes defining triangles
	neighbour across edge_id
	neighbour_segments edge_id of segment in neighbouring triangle
	number_of_boundaries
	"""

	cdef int64_t M, N
	cdef int64_t err = 0
	cdef int64_t k, a, b, c, vol_id, edge_id, len_tag

	cdef char* string
	cdef segment_t* s
	cdef segment_key_t key

	M = triangles.shape[0]
	N = segments.shape[0]

	# Step 1: Populate hashtable. We use a key based on the node_ids of the two nodes defining the segment
	for k in xrange(M):

		a = triangles[k,0]
		b = triangles[k,1]
		c = triangles[k,2]

		# Add segment a,b to hashtable
		key.i = a
		key.j = b
		vol_id = k
		edge_id = 2

		s = find_segment(key)

		if s:
			err = 1
			break

		add_segment(key, vol_id, edge_id)

		# Add segment b,c to hashtable
		key.i = b
		key.j = c
		vol_id = k
		edge_id = 0

		s = find_segment(key)

		if s:
			err = 1
			break

		add_segment(key, vol_id, edge_id)

		# Add segment c,a to hashtable
		key.i = c
		key.j = a
		vol_id = k
		edge_id = 1

		s = find_segment(key)

		if s:
			err = 1
			break

		add_segment(key, vol_id, edge_id)

	if err == 1:
		delete_segment_all()

	assert err != 1, "pmesh2domain.c: build_boundary_dictionary Duplicate segments"

	# Step 2: Go through segments. Those with null tags are added to tag_dict
	for k in xrange(N):

		a = segments[k,0]
		b = segments[k,1]

		string = segment_tags[k]
		len_tag = len(string)

		key.i = a
		key.j = b
		s = find_segment(key)

		if s and len_tag > 0:
			vol_id = s.vol_id
			edge_id = s.edge_id
			tag_dict[(vol_id, edge_id)] = string

		key.i = b
		key.j = a
		s = find_segment(key)

		if s and len_tag > 0:
			vol_id = s.vol_id
			edge_id = s.edge_id
			tag_dict[(vol_id, edge_id)] = string

	delete_segment_all()
	return tag_dict

def build_sides_dictionary(np.ndarray[int64_t, ndim=2, mode="c"] triangles not None, dict sides):

	cdef int64_t i, numTriangle, a, b, c

	numTriangle = triangles.shape[0]

	for i in xrange(numTriangle):

		a = triangles[i,0]
		b = triangles[i,1]
		c = triangles[i,2]

		sides[(a,b)] = 3*i + 2 # sides[a,b] = (id, 2) #(id, face)
		sides[(b,c)] = 3*i + 0 # sides[b,c] = (id, 0) #(id, face)
		sides[(c,a)] = 3*i + 1 # sides[c,a] = (id, 1) #(id, face)

	return sides











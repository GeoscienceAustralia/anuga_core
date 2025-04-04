#cython: wraparound=False, boundscheck=False, cdivision=True, profile=False, nonecheck=False, overflowcheck=False, cdivision_warnings=False, unraisable_tracebacks=False
import cython
from libc.stdlib cimport malloc, free
from libc.stdint cimport int64_t

from cpython.pycapsule cimport *
# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np
# declare the interface to the C code
cdef extern from "sparse_dok.c":
	ctypedef struct UT_hash_handle:
		UT_hash_handle* tbl
		void* prev
		void* next
		UT_hash_handle* hh_prev
		UT_hash_handle* hh_next
		void* key
		unsigned keylen
		unsigned hashv
	ctypedef struct edge_key_t:
		int64_t i
		int64_t j
	ctypedef struct edge_t:
		edge_key_t key
		double entry
		UT_hash_handle hh
	ctypedef struct sparse_dok:
		edge_t* edgetable
		int64_t num_entries
		int64_t num_rows
	void delete_dok_matrix(sparse_dok* mat)
	void sort_by_key(sparse_dok* hashtable)
	void add_dok_entry(sparse_dok* edgetable, edge_key_t key, double value)
	sparse_dok* make_dok()

cdef delete_dok_cap(object cap):
	kill = <sparse_dok* > PyCapsule_GetPointer(cap, "sparse dok")
	if kill != NULL:
		delete_dok_matrix(kill)

cdef int64_t _serialise(sparse_dok* dok, dict serial_dok):

	sort_by_key(dok)
	cdef int64_t num_entries
	cdef int64_t k,i,j
	cdef double val
	cdef edge_t* edge

	num_entries = dok.num_entries
	edge = dok.edgetable

	for k in xrange(num_entries):

		i = edge.key.i
		j = edge.key.j
		val = edge.entry

		serial_dok[(i,j)] = val

		edge = <edge_t* > edge.hh.next

	return 0

cdef int64_t _deserialise(sparse_dok* dok, dict serial_dok):

	cdef int64_t num_entries, k, i, j
	cdef double val
	cdef list items
	cdef list keys
	cdef edge_key_t key

	items = serial_dok.items()
	keys = serial_dok.keys()
	num_entries = len(serial_dok)

	for k in xrange(num_entries):

		val = serial_dok[keys[k]]

		i = keys[k][0]
		j = keys[k][1]

		key.i = i
		key.j = j

		add_dok_entry(dok, key, val)

	return 0

def serialise_dok(object sparse_dok_cap):

	cdef int64_t err
	cdef sparse_dok* dok
	cdef dict serial_sparse_dok

	dok = <sparse_dok* > PyCapsule_GetPointer(sparse_dok_cap, "sparse dok")
	serial_sparse_dok = {}

	err = _serialise(dok, serial_sparse_dok)

	assert err == 0, "sparse_matrix.serialise_dok: error in serialising sparse_dok"

	return serial_sparse_dok

def deserialise_dok(dict serial_sparse_dok):

	cdef int64_t err
	cdef object sparse_dok_cap
	cdef sparse_dok* dok

	dok = make_dok()

	err = _deserialise(dok, serial_sparse_dok)

	return PyCapsule_New(<void* > dok, "sparse dok", <PyCapsule_Destructor> delete_dok_cap)

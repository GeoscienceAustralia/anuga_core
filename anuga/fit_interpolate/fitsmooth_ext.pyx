#cython: wraparound=False, boundscheck=False, cdivision=True, profile=False, nonecheck=False, overflowcheck=False, cdivision_warnings=False, unraisable_tracebacks=False
import cython
from libc.stdlib cimport malloc, free
from libc.stdint cimport int64_t

from cpython.pycapsule cimport *
# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np
# declare the interface to the C code
cdef extern from "fitsmooth.c":
	ctypedef struct quad_tree:
		double xmin, xmax, ymin, ymax
		int64_t count
		quad_tree* parent
		quad_tree* q[4]
		triangle* leaves
		triangle* end_leaves
	ctypedef struct UT_hash_handle:
		pass
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
	ctypedef struct triangle:
		double x1, y1
		double x2, y2
		double x3, y3
		int64_t index
		double nx1, ny1
		double nx2, ny2
		double nx3, ny3
		triangle* next
	ctypedef struct sparse_csr:
		double* data
		int64_t* colind
		int64_t* row_ptr
		int64_t num_rows
		int64_t num_entries
	void delete_quad_tree(quad_tree* tree)
	quad_tree* _build_quad_tree(int64_t n, int64_t* triangles, double* vertex_coordinates, double* extents)
	void delete_dok_matrix(sparse_dok* mat)
	sparse_dok* make_dok()
	int64_t _build_smoothing_matrix(int64_t n, int64_t* triangles, double* areas, double* vertex_coordinates, int64_t* strides, sparse_dok* smoothing_mat)
	int64_t _build_matrix_AtA_Atz_points(int64_t N, int64_t* triangles, double* point_coordinates, double* point_values, int64_t zdims, int64_t npts, sparse_dok* AtA, double** Atz, quad_tree* quadtree)
	void _combine_partial_AtA_Atz(sparse_dok* dok_AtA1, sparse_dok* dok_AtA2, double* Atz1, double* Atz2, int64_t n, int64_t zdim)
	triangle* search(quad_tree* node ,double xp, double yp)
	double* calculate_sigma(triangle* T, double x, double y)
	int64_t quad_tree_node_count(quad_tree* tree)
	int64_t get_dok_rows(sparse_dok* dok)
	edge_t* find_dok_entry(sparse_dok* edgetable, edge_key_t key)
	void add_sparse_dok(sparse_dok* dok1, double mult1, sparse_dok* dok2, double mult2)
	sparse_csr* make_csr()
	void delete_csr_matrix(sparse_csr* mat)
	void convert_to_csr_ptr(sparse_csr* new_csr, sparse_dok* hashtable)

cdef delete_quad_tree_cap(object cap):
	kill = <quad_tree* > PyCapsule_GetPointer(cap, "quad tree")
	if kill != NULL:
		delete_quad_tree(kill)

cdef delete_dok_cap(object cap):
	kill = <sparse_dok* > PyCapsule_GetPointer(cap, "sparse dok")
	if kill != NULL:
		delete_dok_matrix(kill)

cdef c_double_array_to_list(double* mat, int64_t cols):
	cdef int64_t j
	cdef list lst
	lst = []
	if not(isinstance(lst, list)):
		return None
	for j in xrange(cols):
		try:
			lst.append(mat[j])
		except:
			return None
	return lst

cdef c_int_array_to_list(int64_t* mat, int64_t cols):
	cdef int64_t j
	cdef list lst
	lst = []
	if not(isinstance(lst, list)):
		return None
	for j in xrange(cols):
		try:
			lst.append(mat[j])
		except:
			return None
	return lst

def build_quad_tree(np.ndarray[int64_t, ndim=2, mode="c"] triangles not None,\
					np.ndarray[double, ndim=2, mode="c"] vertex_coordinates not None,\
					np.ndarray[double, ndim=1, mode="c"] extents not None):
	
	cdef int64_t n

	n = triangles.shape[0]

	return PyCapsule_New(<void* > _build_quad_tree(n, &triangles[0,0], &vertex_coordinates[0,0], &extents[0]), "quad tree", <PyCapsule_Destructor> delete_quad_tree_cap)

def build_smoothing_matrix(np.ndarray[int64_t, ndim=2, mode="c"] triangles not None,\
							np.ndarray[double, ndim=1, mode="c"] areas not None,\
							np.ndarray[double, ndim=2, mode="c"] vertex_coordinates not None):
	
	cdef int64_t err, n
	cdef sparse_dok* smoothing_mat

	n = triangles.shape[0]
	smoothing_mat = make_dok()

	err = _build_smoothing_matrix(n, &triangles[0,0], &areas[0], &vertex_coordinates[0,0], <int64_t* > &vertex_coordinates.strides[0], smoothing_mat)

	assert err == 0, "Unknown Error"

	return PyCapsule_New(<void* > smoothing_mat, "sparse dok", <PyCapsule_Destructor> delete_dok_cap)


def build_matrix_AtA_Atz_points(object tree, int64_t N,\
							np.ndarray[int64_t, ndim=2, mode="c"] triangles not None,\
							np.ndarray[double, ndim=2, mode="c"] point_coordinates not None,\
							np.ndarray z not None,\
							int64_t zdims,\
							int64_t npts):

	cdef quad_tree* quadtree
	cdef sparse_dok* dok_AtA
	cdef object AtA_cap
	cdef double** Atz
	cdef list Atz_ret
	cdef int64_t err
	cdef int64_t i

	z = np.ascontiguousarray(z)

	quadtree = <quad_tree* > PyCapsule_GetPointer(tree, "quad tree")

	dok_AtA = make_dok()

	Atz = <double** > malloc(zdims * sizeof(double*))
	for i in xrange(zdims):
		Atz[i] = <double* > malloc(N * sizeof(double))

	err = _build_matrix_AtA_Atz_points(N, &triangles[0,0],\
										&point_coordinates[0,0],\
										<double* > z.data,\
										zdims,\
										npts,\
										dok_AtA,\
										Atz,\
										quadtree)

	assert err == 0, "Unknown Error"

	AtA_cap = PyCapsule_New(<void* > dok_AtA, "sparse dok", <PyCapsule_Destructor> delete_dok_cap)

	Atz_ret = []
	for i in xrange(zdims):
		Atz_ret.append(c_double_array_to_list(Atz[i],N))
		free(Atz[i])
	free(Atz)

	return [AtA_cap, Atz_ret]

def combine_partial_AtA_Atz(object AtA_cap1, object AtA_cap2,\
							np.ndarray[double, ndim=1, mode="c"] Atz1,\
							np.ndarray[double, ndim=1, mode="c"] Atz2,\
							int64_t zdim,\
							int64_t n):
	
	cdef sparse_dok* dok_AtA1
	cdef sparse_dok* dok_AtA2

	dok_AtA1 = <sparse_dok* > PyCapsule_GetPointer(AtA_cap1, "sparse dok")
	dok_AtA2 = <sparse_dok* > PyCapsule_GetPointer(AtA_cap2, "sparse dok")

	_combine_partial_AtA_Atz(dok_AtA1, dok_AtA2, &Atz1[0], &Atz2[0], n, zdim)

def individual_tree_search(object tree, np.ndarray[double, ndim=1, mode="c"] point):

	cdef quad_tree* quadtree
	cdef double xp,yp
	cdef double* sigma
	cdef triangle* T
	cdef int64_t found
	cdef int64_t index
	cdef list sigmalist

	quadtree = <quad_tree* > PyCapsule_GetPointer(tree, "quad tree")
	xp = point[0]
	yp = point[1]

	T = search(quadtree, xp, yp)
	sigmalist = []

	if T != NULL:
		sigma = calculate_sigma(T, xp, yp)
		sigmalist = c_double_array_to_list(sigma, 3)
		free(sigma)
		found = 1
		index = T.index
	else:
		sigmalist = [-1, -1, -1]
		index = -10
		found = 0

	return [found, sigmalist, index]

def items_in_tree(object tree):

	cdef quad_tree* quadtree

	quadtree = <quad_tree* > PyCapsule_GetPointer(tree, "quad tree")

	return quadtree.count

def return_full_D(object D_cap, int64_t n):

	cdef sparse_dok* D_mat
	cdef int64_t i,j
	cdef edge_key_t key
	cdef edge_t* s
	cdef list ret_D
	cdef list temp

	D_mat = <sparse_dok* > PyCapsule_GetPointer(D_cap, "sparse dok")

	assert D_mat.num_rows <= n and get_dok_rows(D_mat) <= n, "fitsmooth.return_full_D: sparse_dok is bigger than size specified for return."

	ret_D = []
	for i in xrange(n):
		temp = []
		for j in xrange(n):
			key.i = i
			key.j = j
			s = find_dok_entry(D_mat, key)
			if s:
				temp.append(s.entry)
			else:
				temp.append(float(0))
		ret_D.append(temp)

	return ret_D

def build_matrix_B(object smoothing_mat_cap, object AtA_cap, double alpha):

	cdef sparse_dok* smoothing_mat
	cdef sparse_dok* dok_AtA
	cdef sparse_csr* B
	cdef list data, colind, row_ptr

	smoothing_mat = <sparse_dok* > PyCapsule_GetPointer(smoothing_mat_cap, "sparse dok")
	dok_AtA = <sparse_dok* > PyCapsule_GetPointer(AtA_cap, "sparse dok")

	add_sparse_dok(smoothing_mat, alpha, dok_AtA, 1)

	B = make_csr()
	convert_to_csr_ptr(B, smoothing_mat)

	data = c_double_array_to_list(B.data, B.num_entries)
	colind = c_int_array_to_list(B.colind, B.num_entries)
	row_ptr = c_int_array_to_list(B.row_ptr, B.num_rows)

	delete_csr_matrix(B)

	return [data, colind, row_ptr]

#cython: wraparound=False, boundscheck=False, cdivision=True, profile=False, nonecheck=False, overflowcheck=False, cdivision_warnings=False, unraisable_tracebacks=False
import cython
from libc.stdint cimport int64_t

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

ctypedef int64_t keyint

# declare the interface to the C code
cdef extern from "quantity_openmp.c":
  int64_t _compute_gradients(keyint N, double* centroids, double* centroid_values, int64_t* number_of_boundaries, int64_t* surrogate_neighbours, double* a, double* b)
  int64_t _compute_local_gradients(keyint N, double* vertex_coordinates, double* vertex_values, double* a, double* b)
  int64_t _extrapolate_from_gradient(keyint N, double* centroids, double* centroid_values, double* vertex_coordinates, double* vertex_values, double* edge_values, double* a, double* b)
  int64_t _extrapolate_and_limit_from_gradient(keyint N, double beta, double* centroids, int64_t* neighbours, double* centroid_values, double* vertex_coordinates, double* vertex_values, double* edge_values, double* phi, double* x_gradient, double* y_gradient) 
  int64_t _limit_vertices_by_all_neighbours(keyint N, double beta, double* centroid_values, double* vertex_values, double* edge_values, int64_t* neighbours, double* x_gradient, double* y_gradient)
  int64_t _limit_edges_by_all_neighbours(keyint N, double beta, double* centroid_values, double* vertex_values, double* edge_values, int64_t* neighbours, double* x_gradient, double* y_gradient)
  int64_t _limit_edges_by_neighbour(keyint N, double beta, double* centroid_values, double* vertex_values, double* edge_values, int64_t* neighbours)
  int64_t _limit_gradient_by_neighbour(keyint N, double beta, double* centroid_values, double* vertex_values, double* edge_values, double* x_gradient, double* y_gradient, int64_t* neighbours)
  int64_t _bound_vertices_below_by_constant(keyint N, double bound, double* centroid_values, double* vertex_values, double* edge_values, double* x_gradient, double* y_gradient)
  int64_t _bound_vertices_below_by_quantity(keyint N, double* bound_vertex_values, double* centroid_values, double* vertex_values, double* edge_values, double* x_gradient, double* y_gradient)
  int64_t _interpolate(keyint N, double* vertex_values, double* edge_values, double* centroid_values)
  int64_t _interpolate_from_vertices_to_edges(keyint N, double* vertex_values, double* edge_values)
  int64_t _interpolate_from_edges_to_vertices(keyint N, double* vertex_values, double* edge_values)
  int64_t _backup_centroid_values(keyint N, double* centroid_values, double* centroid_backup_values)
  int64_t _saxpy_centroid_values(keyint N, double a, double b, double* centroid_values, double* centroid_backup_values)
  int64_t _update(keyint N, double timestep, double* centroid_values, double* explicit_update, double* semi_implicit_update)
  int64_t _average_vertex_values(keyint N, int64_t* vertex_value_indices, int64_t* number_of_triangles_per_node, double* vertex_values, double* A)
  int64_t _average_centroid_values(keyint N, int64_t* vertex_value_indices, int64_t* number_of_triangles_per_node, double* centroid_values, double* A)
  int64_t _set_vertex_values_c(keyint num_verts, int64_t* vertices, int64_t* node_index, int64_t* number_of_triangles_per_node, int64_t* vertex_value_indices, double* vertex_values, double* A)
  int64_t _min_and_max_centroid_values(keyint N, double* qc, double* qv, int64_t* neighbours, double* qmin, double* qmax)

cdef extern from "util_ext.h":
  void _limit_old(int64_t N, double beta, double* qc, double* qv, double* qmin, double* qmax)


def update(object quantity, double timestep):
  """Update centroid values based on values stored in
    explicit_update and semi_implicit_update as well as given timestep

    Function implementing forcing terms must take on argument
    which is the domain and they must update either explicit
    or implicit updates, e,g,:

    def gravity(domain):
        ....
        domain.quantities['xmomentum'].explicit_update = ...
        domain.quantities['ymomentum'].explicit_update = ...



    Explicit terms must have the form

        G(q, t)

    and explicit scheme is

       q^{(n+1}) = q^{(n)} + delta_t G(q^{n}, n delta_t)


    Semi implicit forcing terms are assumed to have the form

       G(q, t) = H(q, t) q

    and the semi implicit scheme will then be

      q^{(n+1}) = q^{(n)} + delta_t H(q^{n}, n delta_t) q^{(n+1})"""

  cdef np.ndarray[double, ndim=1, mode="c"] centroid_values
  cdef np.ndarray[double, ndim=1, mode="c"] explicit_update
  cdef np.ndarray[double, ndim=1, mode="c"] semi_implicit_update

  cdef keyint N
  cdef int64_t err

  centroid_values = quantity.centroid_values
  explicit_update = quantity.explicit_update
  semi_implicit_update = quantity.semi_implicit_update

  N = centroid_values.shape[0]

  err = _update(N, timestep, &centroid_values[0], &explicit_update[0], &semi_implicit_update[0])

  assert err == 0, "update: division by zero in semi implicit update - call Stephen :)"


def backup_centroid_values(object quantity):

  cdef np.ndarray[double, ndim=1, mode="c"] centroid_values
  cdef np.ndarray[double, ndim=1, mode="c"] centroid_backup_values

  cdef keyint N
  cdef int64_t err

  centroid_values = quantity.centroid_values
  centroid_backup_values = quantity.centroid_backup_values

  N = centroid_values.shape[0]

  err = _backup_centroid_values(N, &centroid_values[0], &centroid_backup_values[0])

def saxpy_centroid_values(object quantity, double a, double b):

  cdef np.ndarray[double, ndim=1, mode="c"] centroid_values
  cdef np.ndarray[double, ndim=1, mode="c"] centroid_backup_values

  cdef keyint N
  cdef int64_t err

  centroid_values = quantity.centroid_values
  centroid_backup_values = quantity.centroid_backup_values

  N = centroid_values.shape[0]

  err = _saxpy_centroid_values(N, a, b, &centroid_values[0], &centroid_backup_values[0])


def set_vertex_values_c(object quantity, np.ndarray[int64_t, ndim=1, mode="c"] vertices not None, np.ndarray[double, ndim=1, mode="c"] A not None):

  cdef object domain
  cdef object mesh

  cdef np.ndarray[double, ndim=2, mode="c"] vertex_values
  cdef np.ndarray[int64_t, ndim=1, mode="c"] node_index
  cdef np.ndarray[int64_t, ndim=1, mode="c"] number_of_triangles_per_node
  cdef np.ndarray[int64_t, ndim=1, mode="c"] vertex_value_indices

  cdef keyint N
  cdef int64_t err
  cdef keyint num_verts

  domain = quantity.domain
  mesh = domain.mesh

  vertex_values = quantity.vertex_values
  node_index = mesh.node_index
  number_of_triangles_per_node = mesh.number_of_triangles_per_node
  vertex_value_indices = mesh.vertex_value_indices

  num_verts = vertices.shape[0]

  err = _set_vertex_values_c(num_verts, &vertices[0], &node_index[0], &number_of_triangles_per_node[0], &vertex_value_indices[0], &vertex_values[0,0], &A[0])


def interpolate(object quantity):
  
  cdef np.ndarray[double, ndim=2, mode="c"] vertex_values
  cdef np.ndarray[double, ndim=2, mode="c"] edge_values
  cdef np.ndarray[double, ndim=1, mode="c"] centroid_values

  cdef keyint N
  cdef int64_t err

  vertex_values = quantity.vertex_values
  edge_values = quantity.edge_values
  centroid_values = quantity.centroid_values

  N = vertex_values.shape[0]

  err = _interpolate(N, &vertex_values[0,0], &edge_values[0,0], &centroid_values[0])

  assert err == 0, "Interpolate: could not be computed"

def interpolate_from_vertices_to_edges(object quantity):

  cdef np.ndarray[double, ndim=2, mode="c"] vertex_values
  cdef np.ndarray[double, ndim=2, mode="c"] edge_values

  cdef keyint N
  cdef int64_t err

  vertex_values = quantity.vertex_values
  edge_values = quantity.edge_values

  N = vertex_values.shape[0]

  err = _interpolate_from_vertices_to_edges(N, &vertex_values[0,0], &edge_values[0,0])

  assert err == 0, "Interpolate: could not be computed"

def interpolate_from_edges_to_vertices(object quantity):

  cdef np.ndarray[double, ndim=2, mode="c"] vertex_values
  cdef np.ndarray[double, ndim=2, mode="c"] edge_values

  cdef keyint N
  cdef int64_t err

  vertex_values = quantity.vertex_values
  edge_values = quantity.edge_values

  N = vertex_values.shape[0]

  err = _interpolate_from_edges_to_vertices(N, &vertex_values[0,0], &edge_values[0,0])

  assert err == 0, "Interpolate: could not be computed"

def average_vertex_values(np.ndarray[int64_t, ndim=1, mode="c"] vertex_value_indices not None, np.ndarray[int64_t, ndim=1, mode="c"] number_of_triangles_per_node not None, np.ndarray[double, ndim=2, mode="c"] vertex_values not None, np.ndarray[double, ndim=1, mode="c"] A not None):

  cdef keyint N
  cdef int64_t err

  N = vertex_value_indices.shape[0]

  err = _average_vertex_values(N, &vertex_value_indices[0], &number_of_triangles_per_node[0], &vertex_values[0,0], &A[0])

  assert err == 0, "average_vertex_values: could not be computed"

def average_centroid_values(np.ndarray[int64_t, ndim=1, mode="c"] vertex_value_indices not None, np.ndarray[int64_t, ndim=1, mode="c"] number_of_triangles_per_node not None, np.ndarray[double, ndim=1, mode="c"] centroid_values not None, np.ndarray[double, ndim=1, mode="c"] A not None):

  cdef keyint N
  cdef int64_t err

  N = vertex_value_indices.shape[0]

  err = _average_centroid_values(N, &vertex_value_indices[0], &number_of_triangles_per_node[0], &centroid_values[0], &A[0])

  assert err == 0, "average_centroid_values: could not be computed"

def extrapolate_from_gradient(object quantity):

  cdef object domain

  cdef np.ndarray[double, ndim=2, mode="c"] centroids
  cdef np.ndarray[double, ndim=1, mode="c"] centroid_values
  cdef np.ndarray[double, ndim=2, mode="c"] vertex_coordinates
  cdef np.ndarray[double, ndim=2, mode="c"] vertex_values
  cdef np.ndarray[double, ndim=2, mode="c"] edge_values
  cdef np.ndarray[int64_t, ndim=1, mode="c"] number_of_boundaries
  cdef np.ndarray[int64_t, ndim=2, mode="c"] surrogate_neighbours
  cdef np.ndarray[double, ndim=1, mode="c"] x_gradient
  cdef np.ndarray[double, ndim=1, mode="c"] y_gradient

  cdef keyint N
  cdef int64_t err

  domain = quantity.domain

  centroids = domain.centroid_coordinates
  centroid_values = quantity.centroid_values
  surrogate_neighbours = domain.surrogate_neighbours
  number_of_boundaries = domain.number_of_boundaries
  vertex_coordinates = domain.vertex_coordinates
  vertex_values = quantity.vertex_values
  edge_values = quantity.edge_values
  x_gradient = quantity.x_gradient
  y_gradient = quantity.y_gradient

  N = centroid_values.shape[0]

  err = _extrapolate_from_gradient(N,\
							&centroids[0,0],\
							&centroid_values[0],\
							&vertex_coordinates[0,0],\
							&vertex_values[0,0],\
							&edge_values[0,0],\
							&x_gradient[0],\
							&y_gradient[0])

  assert err == 0, "Internal function _extrapolate failed"

def compute_local_gradients(object quantity):

  cdef object domain

  cdef np.ndarray[double, ndim=2, mode="c"] vertex_coordinates
  cdef np.ndarray[double, ndim=2, mode="c"] vertex_values
  cdef np.ndarray[double, ndim=1, mode="c"] x_gradient
  cdef np.ndarray[double, ndim=1, mode="c"] y_gradient

  cdef keyint N
  cdef int64_t err

  domain = quantity.domain

  vertex_coordinates = domain.vertex_coordinates
  vertex_values = quantity.vertex_values
  x_gradient = quantity.x_gradient
  y_gradient = quantity.y_gradient

  N = vertex_values.shape[0]

  err = _compute_local_gradients(N, &vertex_coordinates[0,0], &vertex_values[0,0], &x_gradient[0], &y_gradient[0])

  assert err == 0, "Internal function _compute_local_gradient failed"

def extrapolate_second_order_and_limit_by_edge(object quantity):

  cdef object domain

  cdef np.ndarray[double, ndim=2, mode="c"] domain_centroids
  cdef np.ndarray[double, ndim=2, mode="c"] domain_vertex_coordinates
  cdef np.ndarray[int64_t, ndim=1, mode="c"] domain_number_of_boundaries
  cdef np.ndarray[int64_t, ndim=2, mode="c"] domain_surrogate_neighbours
  cdef np.ndarray[int64_t, ndim=2, mode="c"] domain_neighbours

  cdef np.ndarray[double, ndim=1, mode="c"] quantity_centroid_values
  cdef np.ndarray[double, ndim=2, mode="c"] quantity_vertex_values
  cdef np.ndarray[double, ndim=2, mode="c"] quantity_edge_values
  cdef np.ndarray[double, ndim=1, mode="c"] quantity_phi
  cdef np.ndarray[double, ndim=1, mode="c"] quantity_x_gradient
  cdef np.ndarray[double, ndim=1, mode="c"] quantity_y_gradient

  cdef keyint ntri
  cdef double beta
  cdef int64_t err

  domain = quantity.object

  domain_centroids = domain.centroid_coordinates
  domain_surrogate_neighbours = domain.surrogate_neighbours
  domain_number_of_boundaries = domain.number_of_boundaries
  domain_vertex_coordinates = domain.vertex_coordinates
  domain_neighbours = domain.neighbours

  quantity_centroid_values = quantity.centroid_values
  quantity_vertex_values = quantity.vertex_values
  quantity_edge_values = quantity.edge_values
  quantity_phi = quantity.phi
  quantity_x_gradient = quantity.x_gradient
  quantity_y_gradient = quantity.y_gradient

  beta = quantity.beta

  ntri = quantity_centroid_values.shape[0]

  err = _compute_gradients(ntri,\
						&domain_centroids[0,0],\
						&quantity_centroid_values[0],\
						&domain_number_of_boundaries[0],\
						&domain_surrogate_neighbours[0,0],\
						&quantity_x_gradient[0],\
						&quantity_y_gradient[0])

  assert err == 0, "Internal function _compute_gradient failed"

  err = _extrapolate_from_gradient(ntri,\
						&domain_centroids[0,0],\
						&quantity_centroid_values[0],\
						&domain_vertex_coordinates[0,0],\
						&quantity_vertex_values[0,0],\
						&quantity_edge_values[0,0],\
						&quantity_x_gradient[0],\
						&quantity_y_gradient[0])

  assert err == 0, "Internal function _extrapolate_from_gradient failed"

  err = _limit_edges_by_all_neighbours(ntri, beta,\
						&quantity_centroid_values[0],\
						&quantity_vertex_values[0,0],\
						&quantity_edge_values[0,0],\
						&domain_neighbours[0,0],\
						&quantity_x_gradient[0],\
						&quantity_y_gradient[0])

  assert err == 0, "Internal function _limit_edges_by_all_neighbours failed"

def extrapolate_second_order_and_limit_by_vertex(object quantity):

  cdef object domain

  cdef np.ndarray[double, ndim=2, mode="c"] domain_centroids
  cdef np.ndarray[double, ndim=2, mode="c"] domain_vertex_coordinates
  cdef np.ndarray[int64_t, ndim=1, mode="c"] domain_number_of_boundaries
  cdef np.ndarray[int64_t, ndim=2, mode="c"] domain_surrogate_neighbours
  cdef np.ndarray[int64_t, ndim=2, mode="c"] domain_neighbours

  cdef np.ndarray[double, ndim=1, mode="c"] quantity_centroid_values
  cdef np.ndarray[double, ndim=2, mode="c"] quantity_vertex_values
  cdef np.ndarray[double, ndim=2, mode="c"] quantity_edge_values
  cdef np.ndarray[double, ndim=1, mode="c"] quantity_phi
  cdef np.ndarray[double, ndim=1, mode="c"] quantity_x_gradient
  cdef np.ndarray[double, ndim=1, mode="c"] quantity_y_gradient

  cdef keyint ntri
  cdef double beta
  cdef int64_t err

  domain = quantity.object

  domain_centroids = domain.centroid_coordinates
  domain_surrogate_neighbours = domain.surrogate_neighbours
  domain_number_of_boundaries = domain.number_of_boundaries
  domain_vertex_coordinates = domain.vertex_coordinates
  domain_neighbours = domain.neighbours

  quantity_centroid_values = quantity.centroid_values
  quantity_vertex_values = quantity.vertex_values
  quantity_edge_values = quantity.edge_values
  quantity_phi = quantity.phi
  quantity_x_gradient = quantity.x_gradient
  quantity_y_gradient = quantity.y_gradient

  beta = quantity.beta

  ntri = quantity_centroid_values.shape[0]

  err = _compute_gradients(ntri,\
						&domain_centroids[0,0],\
						&quantity_centroid_values[0],\
						&domain_number_of_boundaries[0],\
						&domain_surrogate_neighbours[0,0],\
						&quantity_x_gradient[0],\
						&quantity_y_gradient[0])

  assert err == 0, "Internal function _compute_gradient failed"

  err = _extrapolate_from_gradient(ntri,\
						&domain_centroids[0,0],\
						&quantity_centroid_values[0],\
						&domain_vertex_coordinates[0,0],\
						&quantity_vertex_values[0,0],\
						&quantity_edge_values[0,0],\
						&quantity_x_gradient[0],\
						&quantity_y_gradient[0])

  assert err == 0, "Internal function _extrapolate_from_gradient failed"

  err = _limit_vertices_by_all_neighbours(ntri, beta,\
						&quantity_centroid_values[0],\
						&quantity_vertex_values[0,0],\
						&quantity_edge_values[0,0],\
						&domain_neighbours[0,0],\
						&quantity_x_gradient[0],\
						&quantity_y_gradient[0])

  assert err == 0, "Internal function _limit_edges_by_all_neighbours failed"

def compute_gradients(object quantity):

  cdef object domain

  cdef np.ndarray[double, ndim=2, mode="c"] centroids
  cdef np.ndarray[double, ndim=1, mode="c"] centroid_values
  cdef np.ndarray[double, ndim=2, mode="c"] vertex_coordinates
  cdef np.ndarray[double, ndim=2, mode="c"] vertex_values
  cdef np.ndarray[double, ndim=2, mode="c"] edge_values
  cdef np.ndarray[int64_t, ndim=1, mode="c"] number_of_boundaries
  cdef np.ndarray[int64_t, ndim=2, mode="c"] surrogate_neighbours
  cdef np.ndarray[double, ndim=1, mode="c"] x_gradient
  cdef np.ndarray[double, ndim=1, mode="c"] y_gradient

  cdef keyint N
  cdef int64_t err

  domain = quantity.domain

  centroids = domain.centroid_coordinates
  centroid_values = quantity.centroid_values
  surrogate_neighbours = domain.surrogate_neighbours
  number_of_boundaries = domain.number_of_boundaries
  vertex_coordinates = domain.vertex_coordinates
  vertex_values = quantity.vertex_values
  edge_values = quantity.edge_values
  x_gradient = quantity.x_gradient
  y_gradient = quantity.y_gradient

  N = centroid_values.shape[0]

  err = _compute_gradients(N,\
						&centroids[0,0],\
						&centroid_values[0],\
						&number_of_boundaries[0],\
						&surrogate_neighbours[0,0],\
						&x_gradient[0],\
						&y_gradient[0])

  assert err == 0, "Gradient could not be computed"

def limit_old(object quantity):
  
  cdef object domain

  cdef np.ndarray[double, ndim=1, mode="c"] qc
  cdef np.ndarray[double, ndim=2, mode="c"] qv
  cdef np.ndarray[int64_t, ndim=2, mode="c"] neighbours

  cdef keyint N
  cdef double beta_w
  cdef int64_t err

  domain = quantity.domain

  neighbours = domain.neighbours

  beta_w = domain.beta_w

  qc = quantity.centroid_values
  qv = quantity.vertex_values

  N = qc.shape[0]

  cdef np.ndarray[double, ndim=1, mode="c"] qmin = np.empty(N, dtype=np.float64)
  cdef np.ndarray[double, ndim=1, mode="c"] qmax = np.empty(N, dtype=np.float64)

  err = _min_and_max_centroid_values(N, &qc[0], &qv[0,0], &neighbours[0,0], &qmin[0], &qmax[0])

  assert err == 0, "Internal function _min_and_max_centroid_values failed"

  _limit_old(N, beta_w, &qc[0], &qv[0,0], &qmin[0], &qmax[0])

def limit_vertices_by_all_neighbours(object quantity):
  
  cdef object domain

  cdef np.ndarray[double, ndim=2, mode="c"] vertex_values
  cdef np.ndarray[double, ndim=1, mode="c"] centroid_values
  cdef np.ndarray[double, ndim=2, mode="c"] edge_values
  cdef np.ndarray[int64_t, ndim=2, mode="c"] neighbours
  cdef np.ndarray[double, ndim=1, mode="c"] x_gradient
  cdef np.ndarray[double, ndim=1, mode="c"] y_gradient

  cdef double beta_w
  cdef keyint N
  cdef int64_t err

  domain = quantity.domain

  beta_w = domain.beta_w

  neighbours = domain.neighbours
  centroid_values = quantity.centroid_values
  vertex_values = quantity.vertex_values
  edge_values = quantity.edge_values
  x_gradient = quantity.x_gradient
  y_gradient = quantity.y_gradient
  beta_w = domain.beta_w

  N = centroid_values.shape[0]

  err = _limit_vertices_by_all_neighbours(N, beta_w,\
											&centroid_values[0],\
											&vertex_values[0,0],\
											&edge_values[0,0],\
											&neighbours[0,0],\
											&x_gradient[0],\
											&y_gradient[0])

  assert err == 0, "Internal function _limit_by_vertex failed"

def limit_edges_by_all_neighbours(object quantity):

  cdef object domain

  cdef np.ndarray[double, ndim=2, mode="c"] vertex_values
  cdef np.ndarray[double, ndim=1, mode="c"] centroid_values
  cdef np.ndarray[double, ndim=2, mode="c"] edge_values
  cdef np.ndarray[int64_t, ndim=2, mode="c"] neighbours
  cdef np.ndarray[double, ndim=1, mode="c"] x_gradient
  cdef np.ndarray[double, ndim=1, mode="c"] y_gradient

  cdef double beta_w
  cdef keyint N
  cdef int64_t err

  domain = quantity.domain

  beta_w = domain.beta_w

  neighbours = domain.neighbours
  centroid_values = quantity.centroid_values
  vertex_values = quantity.vertex_values
  edge_values = quantity.edge_values
  x_gradient = quantity.x_gradient
  y_gradient = quantity.y_gradient
  beta_w = domain.beta_w

  N = centroid_values.shape[0]

  err = _limit_edges_by_all_neighbours(N, beta_w,\
											&centroid_values[0],\
											&vertex_values[0,0],\
											&edge_values[0,0],\
											&neighbours[0,0],\
											&x_gradient[0],\
											&y_gradient[0])

  assert err == 0, "Internal function _limit_by_edges failed"

def bound_vertices_below_by_constant(object quantity, double bound):
  
  cdef object domain

  cdef np.ndarray[double, ndim=2, mode="c"] vertex_values
  cdef np.ndarray[double, ndim=1, mode="c"] centroid_values
  cdef np.ndarray[double, ndim=2, mode="c"] edge_values
  cdef np.ndarray[double, ndim=1, mode="c"] x_gradient
  cdef np.ndarray[double, ndim=1, mode="c"] y_gradient

  cdef keyint N
  cdef int64_t err

  domain = quantity.domain

  centroid_values = quantity.centroid_values
  vertex_values = quantity.vertex_values
  edge_values = quantity.edge_values
  x_gradient = quantity.x_gradient
  y_gradient = quantity.y_gradient

  N = centroid_values.shape[0]

  err = _bound_vertices_below_by_constant(N, bound,\
										&centroid_values[0],\
										&vertex_values[0,0],\
										&edge_values[0,0],\
										&x_gradient[0],\
										&y_gradient[0])

  assert err == 0, "Internal function _bound_vertices_below_by_constant failed"

def bound_vertices_below_by_quantity(object quantity, object bounding_quantity):
  
  cdef object domain

  cdef np.ndarray[double, ndim=2, mode="c"] vertex_values
  cdef np.ndarray[double, ndim=1, mode="c"] centroid_values
  cdef np.ndarray[double, ndim=2, mode="c"] edge_values
  cdef np.ndarray[double, ndim=1, mode="c"] x_gradient
  cdef np.ndarray[double, ndim=1, mode="c"] y_gradient
  cdef np.ndarray[double, ndim=2, mode="c"] bound_vertex_values

  cdef keyint N
  cdef int64_t err

  domain = quantity.domain

  centroid_values = quantity.centroid_values
  vertex_values = quantity.vertex_values
  edge_values = quantity.edge_values
  x_gradient = quantity.x_gradient
  y_gradient = quantity.y_gradient
  bound_vertex_values = bounding_quantity.vertex_values

  N = centroid_values.shape[0]

  err = _bound_vertices_below_by_quantity(N,\
  										&bound_vertex_values[0,0],\
										&centroid_values[0],\
										&vertex_values[0,0],\
										&edge_values[0,0],\
										&x_gradient[0],\
										&y_gradient[0])

  assert err == 0, "Internal function _bound_vertices_below_by_quantity failed"

def limit_edges_by_neighbour(object quantity):
  
  cdef object domain

  cdef np.ndarray[double, ndim=2, mode="c"] vertex_values
  cdef np.ndarray[double, ndim=1, mode="c"] centroid_values
  cdef np.ndarray[double, ndim=2, mode="c"] edge_values
  cdef np.ndarray[int64_t, ndim=2, mode="c"] neighbours

  cdef double beta_w
  cdef keyint N
  cdef int64_t err

  domain = quantity.domain

  beta_w = domain.beta_w

  neighbours = domain.neighbours
  centroid_values = quantity.centroid_values
  vertex_values = quantity.vertex_values
  edge_values = quantity.edge_values

  N = centroid_values.shape[0]

  err = _limit_edges_by_neighbour(N, beta_w,\
								&centroid_values[0],\
								&vertex_values[0,0],\
								&edge_values[0,0],\
								&neighbours[0,0])

  assert err == 0, "Internal function _limit_edges_by_neighbour failed"

def limit_gradient_by_neighbour(object quantity):
  
  cdef object domain

  cdef np.ndarray[double, ndim=2, mode="c"] vertex_values
  cdef np.ndarray[double, ndim=1, mode="c"] centroid_values
  cdef np.ndarray[double, ndim=2, mode="c"] edge_values
  cdef np.ndarray[double, ndim=1, mode="c"] x_gradient
  cdef np.ndarray[double, ndim=1, mode="c"] y_gradient
  cdef np.ndarray[int64_t, ndim=2, mode="c"] neighbours

  cdef double beta_w
  cdef keyint N
  cdef int64_t err

  domain = quantity.domain

  beta_w = domain.beta_w

  neighbours = domain.neighbours
  centroid_values = quantity.centroid_values
  vertex_values = quantity.vertex_values
  edge_values = quantity.edge_values
  x_gradient = quantity.x_gradient
  y_gradient = quantity.y_gradient

  N = centroid_values.shape[0]

  err = _limit_gradient_by_neighbour(N, beta_w,\
								&centroid_values[0],\
								&vertex_values[0,0],\
								&edge_values[0,0],\
								&x_gradient[0],\
								&y_gradient[0],\
								&neighbours[0,0])

  assert err == 0, "Internal function _limit_gradient_by_neighbour failed"


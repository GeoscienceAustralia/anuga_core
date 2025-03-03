#cython: wraparound=False, boundscheck=False, cdivision=True, profile=False, nonecheck=False, overflowcheck=False, cdivision_warnings=False, unraisable_tracebacks=False
import cython
from libc.stdint cimport int64_t

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

# declare the interface to the C code
cdef extern from "kinematic_viscosity_operator.c":
  int64_t _build_geo_structure(int64_t n, int64_t tot_len, double* centroids, int64_t* neighbours, double* edgelengths, double* edge_midpoints, int64_t* geo_indices, double* geo_values)
  int64_t _build_elliptic_matrix(int64_t n, int64_t tot_len, int64_t* geo_indices, double* geo_values, double* cell_data, double* bdry_data, double* data, int64_t* colind)
  int64_t _update_elliptic_matrix(int64_t n, int64_t tot_len, int64_t* geo_indices, double* geo_values, double* cell_data, double* bdry_data, double* data, int64_t* colind)

def build_geo_structure(object kv_operator):

  cdef int64_t n, tot_len, err

  cdef np.ndarray[double, ndim=2, mode="c"] centroid_coordinates
  cdef np.ndarray[int64_t, ndim=2, mode="c"] neighbours
  cdef np.ndarray[double, ndim=2, mode="c"] edgelengths
  cdef np.ndarray[double, ndim=2, mode="c"] edge_midpoint_coordinates
  cdef np.ndarray[int64_t, ndim=2, mode="c"] geo_indices
  cdef np.ndarray[double, ndim=2, mode="c"] geo_values

  mesh = kv_operator.mesh

  n = kv_operator.n
  tot_len = kv_operator.tot_len

  centroid_coordinates = mesh.centroid_coordinates
  neighbours = mesh.neighbours
  edgelengths = mesh.edgelengths
  edge_midpoint_coordinates = mesh.edge_midpoint_coordinates
  geo_indices = kv_operator.geo_structure_indices
  geo_values = kv_operator.geo_structure_values

  err = _build_geo_structure(n, tot_len,\
                            &centroid_coordinates[0,0],\
                            &neighbours[0,0],\
                            &edgelengths[0,0],\
                            &edge_midpoint_coordinates[0,0],\
                            &geo_indices[0,0],\
                            &geo_values[0,0])
  
  assert err == 0, "Could not build geo structure"

def build_elliptic_matrix(object kv_operator,\
                          np.ndarray[double, ndim=1, mode="c"] cell_data not None,\
                          np.ndarray[double, ndim=1, mode="c"] bdry_data not None):
  
  cdef int64_t n, tot_len, err

  cdef np.ndarray[int64_t, ndim=2, mode="c"] geo_indices
  cdef np.ndarray[double, ndim=2, mode="c"] geo_values
  cdef np.ndarray[double, ndim=1, mode="c"] _data
  cdef np.ndarray[int64_t, ndim=1, mode="c"] colind

  n = kv_operator.n
  tot_len = kv_operator.tot_len

  geo_indices = kv_operator.geo_structure_indices
  geo_values = kv_operator.geo_structure_values
  _data = kv_operator.operator_data
  colind = kv_operator.operator_colind

  err = _build_elliptic_matrix(n, tot_len,\
                              &geo_indices[0,0],\
                              &geo_values[0,0],\
                              &cell_data[0],\
                              &bdry_data[0],\
                              &_data[0],\
                              &colind[0])
  
  assert err == 0, "Could not get stage height interactions"

def update_elliptic_matrix(object kv_operator,\
                          np.ndarray[double, ndim=1, mode="c"] cell_data not None,\
                          np.ndarray[double, ndim=1, mode="c"] bdry_data not None):
  
  cdef int64_t n, tot_len, err

  cdef np.ndarray[int64_t, ndim=2, mode="c"] geo_indices
  cdef np.ndarray[double, ndim=2, mode="c"] geo_values
  cdef np.ndarray[double, ndim=1, mode="c"] _data
  cdef np.ndarray[int64_t, ndim=1, mode="c"] colind

  n = kv_operator.n
  tot_len = kv_operator.tot_len

  geo_indices = kv_operator.geo_structure_indices
  geo_values = kv_operator.geo_structure_values
  _data = kv_operator.operator_data
  colind = kv_operator.operator_colind

  err = _update_elliptic_matrix(n, tot_len,\
                              &geo_indices[0,0],\
                              &geo_values[0,0],\
                              &cell_data[0],\
                              &bdry_data[0],\
                              &_data[0],\
                              &colind[0])
  
  assert err == 0, "Could not get stage height interactions"

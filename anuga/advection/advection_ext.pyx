#cython: wraparound=False, boundscheck=False, cdivision=True, profile=False, nonecheck=False, overflowcheck=False, cdivision_warnings=False, unraisable_tracebacks=False
import cython
from libc.stdint cimport int64_t

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

# declare the interface to the C code
cdef extern from "advection.c":
	double _compute_fluxes(double* quantity_update, double* quantity_edge, double* quantity_bdry, int64_t* domain_neighbours, int64_t* domain_neighbour_edges, double* domain_normals, double* domain_areas, double* domain_radii, double* domain_edgelengths, int64_t* domain_tri_full_flag, double* domain_velocity, double huge_timestep, double max_timestep, int64_t ntri, int64_t nbdry)

def compute_fluxes(object domain, object quantity, double huge_timestep, double max_timestep):

	cdef np.ndarray[double, ndim=2, mode="c"] quantity_edge
	cdef np.ndarray[double, ndim=1, mode="c"] quantity_bdry
	cdef np.ndarray[double, ndim=1, mode="c"] quantity_update
	cdef np.ndarray[int64_t, ndim=2, mode="c"] domain_neighbours
	cdef np.ndarray[int64_t, ndim=2, mode="c"] domain_neighbour_edges
	cdef np.ndarray[double, ndim=2, mode="c"] domain_normals
	cdef np.ndarray[double, ndim=1, mode="c"] domain_areas
	cdef np.ndarray[double, ndim=1, mode="c"] domain_radii
	cdef np.ndarray[double, ndim=2, mode="c"] domain_edgelengths
	cdef np.ndarray[int64_t, ndim=1, mode="c"] domain_tri_full_flag
	cdef np.ndarray[double, ndim=1, mode="c"] domain_velocity

	cdef int64_t ntri, nbdry
	cdef double timestep

	quantity_edge = quantity.edge_values
	quantity_bdry = quantity.boundary_values
	quantity_update = quantity.explicit_update
	domain_neighbours = domain.neighbours
	domain_neighbour_edges = domain.neighbour_edges
	domain_normals = domain.normals
	domain_areas = domain.areas
	domain_radii = domain.radii
	domain_edgelengths = domain.edgelengths
	domain_tri_full_flag = domain.tri_full_flag
	domain_velocity = domain.velocity

	ntri = quantity_edge.shape[0]
	nbdry = quantity_bdry.shape[0]

	timestep = _compute_fluxes(&quantity_update[0],\
							&quantity_edge[0,0],\
							&quantity_bdry[0],\
							&domain_neighbours[0,0],\
							&domain_neighbour_edges[0,0],\
							&domain_normals[0,0],\
							&domain_areas[0],\
							&domain_radii[0],\
							&domain_edgelengths[0,0],\
							&domain_tri_full_flag[0],\
							&domain_velocity[0],\
							huge_timestep,\
							max_timestep,\
							ntri,\
							nbdry)

	return timestep


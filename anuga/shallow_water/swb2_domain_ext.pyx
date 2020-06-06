#cythonoff: wraparound=False, boundscheck=False, cdivision=True, profile=False, nonecheck=False, overflowcheck=False, cdivision_warnings=False, unraisable_tracebacks=False
import cython

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

cdef extern from "swb2_domain.c":
	double _compute_fluxes_central(int number_of_elements, double timestep, double epsilon, double H0, double g, long* neighbours, long* neighbour_edges, double* normals, double* edgelengths, double* radii, double* areas, long* tri_full_flag, double* stage_edge_values, double* xmom_edge_values, double* ymom_edge_values, double* bed_edge_values, double* stage_boundary_values, double* xmom_boundary_values, double* ymom_boundary_values, long* boundary_flux_type, double* stage_explicit_update, double* xmom_explicit_update, double* ymom_explicit_update, long* already_computed_flux, double* max_speed_array, int optimise_dry_cells, double* stage_centroid_values, double* bed_centroid_values, double* bed_vertex_values)
	double _protect(int N, double minimum_allowed_height, double maximum_allowed_speed, double epsilon, double* wc, double* wv, double* zc, double* zv, double* xmomc, double* ymomc, double* areas)
	int _extrapolate_second_order_edge_sw(int number_of_elements, double epsilon, double minimum_allowed_height, double beta_w, double beta_w_dry, double beta_uh, double beta_uh_dry, double beta_vh, double beta_vh_dry, long* surrogate_neighbours, long* number_of_boundaries, double* centroid_coordinates, double* stage_centroid_values, double* xmom_centroid_values, double* ymom_centroid_values, double* elevation_centroid_values, double* edge_coordinates, double* stage_edge_values, double* xmom_edge_values, double* ymom_edge_values, double* elevation_edge_values, double* stage_vertex_values, double* xmom_vertex_values, double* ymom_vertex_values, double* elevation_vertex_values, int optimise_dry_cells, int extrapolate_velocity_second_order)

def compute_fluxes_ext_central(double timestep,\
							double epsilon,\
							double H0,\
							double g,\
							np.ndarray[long, ndim=2, mode="c"] neighbours not None,\
							np.ndarray[long, ndim=2, mode="c"] neighbour_edges not None,\
							np.ndarray[double, ndim=2, mode="c"] normals not None,\
							np.ndarray[double, ndim=2, mode="c"] edgelengths not None,\
							np.ndarray[double, ndim=1, mode="c"] radii not None,\
							np.ndarray[double, ndim=1, mode="c"] areas not None,\
							np.ndarray[long, ndim=1, mode="c"] tri_full_flag not None,\
							np.ndarray[double, ndim=2, mode="c"] stage_edge_values not None,\
							np.ndarray[double, ndim=2, mode="c"] xmom_edge_values not None,\
							np.ndarray[double, ndim=2, mode="c"] ymom_edge_values not None,\
							np.ndarray[double, ndim=2, mode="c"] bed_edge_values not None,\
							np.ndarray[double, ndim=1, mode="c"] stage_boundary_values not None,\
							np.ndarray[double, ndim=1, mode="c"] xmom_boundary_values not None,\
							np.ndarray[double, ndim=1, mode="c"] ymom_boundary_values not None,\
							np.ndarray[long, ndim=1, mode="c"] boundary_flux_type not None,\
							np.ndarray[double, ndim=1, mode="c"] stage_explicit_update not None,\
							np.ndarray[double, ndim=1, mode="c"] xmom_explicit_update not None,\
							np.ndarray[double, ndim=1, mode="c"] ymom_explicit_update not None,\
							np.ndarray[long, ndim=2, mode="c"] already_computed_flux not None,\
							np.ndarray[double, ndim=1, mode="c"] max_speed_array not None,\
							long optimise_dry_cells,\
							np.ndarray[double, ndim=1, mode="c"] stage_centroid_values not None,\
							np.ndarray[double, ndim=1, mode="c"] bed_centroid_values not None,\
							np.ndarray[double, ndim=2, mode="c"] bed_vertex_values not None):

	cdef int number_of_elements

	number_of_elements = stage_edge_values.shape[0]

	timestep = _compute_fluxes_central(number_of_elements,\
									timestep,\
									epsilon,\
									H0,\
									g,\
									&neighbours[0,0],\
									&neighbour_edges[0,0],\
									&normals[0,0],\
									&edgelengths[0,0],\
									&radii[0],\
									&areas[0],\
									&tri_full_flag[0],\
									&stage_edge_values[0,0],\
									&xmom_edge_values[0,0],\
									&ymom_edge_values[0,0],\
									&bed_edge_values[0,0],\
									&stage_boundary_values[0],\
									&xmom_boundary_values[0],\
									&ymom_boundary_values[0],\
									&boundary_flux_type[0],\
									&stage_explicit_update[0],\
									&xmom_explicit_update[0],\
									&ymom_explicit_update[0],\
									&already_computed_flux[0,0],\
									&max_speed_array[0],\
									optimise_dry_cells,\
									&stage_centroid_values[0],\
									&bed_centroid_values[0],\
									&bed_vertex_values[0,0])

	return timestep

def protect(double minimum_allowed_height,\
			double maximum_allowed_speed,\
			double epsilon,\
			np.ndarray[double, ndim=1, mode="c"] wc not None,\
			np.ndarray[double, ndim=2, mode="c"] wv not None,\
			np.ndarray[double, ndim=1, mode="c"] zc not None,\
			np.ndarray[double, ndim=2, mode="c"] zv not None,\
			np.ndarray[double, ndim=1, mode="c"] xmomc not None,\
			np.ndarray[double, ndim=1, mode="c"] ymomc not None,\
			np.ndarray[double, ndim=1, mode="c"] areas not None):

	cdef int N

	cdef double mass_error

	N = wc.shape[0]

	mass_error = _protect(N,\
					minimum_allowed_height,\
					maximum_allowed_speed,\
					epsilon,\
					&wc[0],\
					&wv[0,0],\
					&zc[0],\
					&zv[0,0],\
					&xmomc[0],\
					&ymomc[0],\
					&areas[0])

	return mass_error

def extrapolate_second_order_edge_sw(object domain,\
									np.ndarray[long, ndim=2, mode="c"] surrogate_neighbours not None,\
									np.ndarray[long, ndim=1, mode="c"] number_of_boundaries not None,\
									np.ndarray[double, ndim=2, mode="c"] centroid_coordinates not None,\
									np.ndarray[double, ndim=1, mode="c"] stage_centroid_values not None,\
									np.ndarray[double, ndim=1, mode="c"] xmom_centroid_values not None,\
									np.ndarray[double, ndim=1, mode="c"] ymom_centroid_values not None,\
									np.ndarray[double, ndim=1, mode="c"] elevation_centroid_values not None,\
									np.ndarray[double, ndim=2, mode="c"] edge_coordinates not None,\
									np.ndarray[double, ndim=2, mode="c"] stage_edge_values not None,\
									np.ndarray[double, ndim=2, mode="c"] xmom_edge_values not None,\
									np.ndarray[double, ndim=2, mode="c"] ymom_edge_values not None,\
									np.ndarray[double, ndim=2, mode="c"] elevation_edge_values not None,\
									np.ndarray[double, ndim=2, mode="c"] stage_vertex_values not None,\
									np.ndarray[double, ndim=2, mode="c"] xmom_vertex_values not None,\
									np.ndarray[double, ndim=2, mode="c"] ymom_vertex_values not None,\
									np.ndarray[double, ndim=2, mode="c"] elevation_vertex_values not None,\
									long optimise_dry_cells,\
									long extrapolate_velocity_second_order):


	cdef double beta_w, beta_w_dry, beta_uh, beta_uh_dry, beta_vh, beta_vh_dry
	cdef double minimum_allowed_height, epsilon
	cdef int number_of_elements, e

	beta_w = domain.beta_w
	beta_w_dry = domain.beta_w_dry
	beta_uh = domain.beta_uh
	beta_uh_dry = domain.beta_uh_dry
	beta_vh = domain.beta_vh
	beta_vh_dry = domain.beta_vh_dry

	minimum_allowed_height = domain.minimum_allowed_height
	epsilon = domain.epsilon

	number_of_elements = stage_centroid_values.shape[0]

	e = _extrapolate_second_order_edge_sw(number_of_elements,\
										epsilon,\
										minimum_allowed_height,\
										beta_w,\
										beta_w_dry,\
										beta_uh,\
										beta_uh_dry,\
										beta_vh,\
										beta_vh_dry,\
										&surrogate_neighbours[0,0],\
										&number_of_boundaries[0],\
										&centroid_coordinates[0,0],\
										&stage_centroid_values[0],\
										&xmom_centroid_values[0],\
										&ymom_centroid_values[0],\
										&elevation_centroid_values[0],\
										&edge_coordinates[0,0],\
										&stage_edge_values[0,0],\
										&xmom_edge_values[0,0],\
										&ymom_edge_values[0,0],\
										&elevation_edge_values[0,0],\
										&stage_vertex_values[0,0],\
										&xmom_vertex_values[0,0],\
										&ymom_vertex_values[0,0],\
										&elevation_vertex_values[0,0],\
										optimise_dry_cells,\
										extrapolate_velocity_second_order)

	if e == -1:
		return None






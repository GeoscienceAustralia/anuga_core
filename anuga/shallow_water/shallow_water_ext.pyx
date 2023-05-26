#cython: wraparound=False, boundscheck=False, cdivision=True, profile=False, nonecheck=False, overflowcheck=False, cdivision_warnings=False, unraisable_tracebacks=False
import cython

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

cdef extern from "shallow_water.c":
	struct domain:
		long number_of_elements
		double epsilon
		double H0
		double g
		long optimise_dry_cells
		double evolve_max_timestep
		long extrapolate_velocity_second_order
		double minimum_allowed_height
		double maximum_allowed_speed
		long low_froude
		long timestep_fluxcalls
		double beta_w
		double beta_w_dry
		double beta_uh
		double beta_uh_dry
		double beta_vh
		double beta_vh_dry
		long max_flux_update_frequency
		long ncol_riverwall_hydraulic_properties
		long* neighbours
		long* neighbour_edges
		long* surrogate_neighbours
		double* normals
		double* edgelengths
		double* radii
		double* areas
		long* edge_flux_type
		long* tri_full_flag
		long* already_computed_flux
		double* max_speed
		double* vertex_coordinates
		double* edge_coordinates
		double* centroid_coordinates
		long* number_of_boundaries
		double* stage_edge_values
		double* xmom_edge_values
		double* ymom_edge_values
		double* bed_edge_values
		double* height_edge_values
		double* stage_centroid_values
		double* xmom_centroid_values
		double* ymom_centroid_values
		double* bed_centroid_values
		double* height_centroid_values
		double* stage_vertex_values
		double* xmom_vertex_values
		double* ymom_vertex_values
		double* bed_vertex_values
		double* height_vertex_values
		double* stage_boundary_values
		double* xmom_boundary_values
		double* ymom_boundary_values
		double* bed_boundary_values
		double* stage_explicit_update
		double* xmom_explicit_update
		double* ymom_explicit_update
		long* flux_update_frequency
		long* update_next_flux
		long* update_extrapolation
		double* edge_timestep
		double* edge_flux_work
		double* pressuregrad_work
		double* x_centroid_work
		double* y_centroid_work
		double* boundary_flux_sum
		long* allow_timestep_increase
		double* riverwall_elevation
		long* riverwall_rowIndex
		double* riverwall_hydraulic_properties
	struct edge:
		pass
	int _rotate(double *q, double n1, double n2)
	int _flux_function_central(double* q_left, double* q_right, double z_left, double z_right, double n1, double n2, double epsilon, double h0, double limiting_threshold, double g, double* edgeflux, double* max_speed)
	int _extrapolate_second_order_sw(domain* D)
	double _compute_fluxes_central_structure(domain* D)
	int _gravity(domain* D)
	int _gravity_wb(domain* D)
	double _compute_fluxes_central_wb(domain* D)
	double _compute_fluxes_central_wb_3(domain* D)
	int _protect(int N, double minimum_allowed_height, double maximum_allowed_speed, double epsilon, double* wc, double* zc, double* xmomc, double* ymomc)
	int _balance_deep_and_shallow(int N, double* wc, double* zc, double* wv, double* zv, double* hvbar, double* xmomc, double* ymomc, double* xmomv, double* ymomv, double H0, int tight_slope_limiters, int use_centroid_velocities, double alpha_balance)
	void _manning_friction_flat(double g, double eps, int N, double* w, double* zv, double* uh, double* vh, double* eta, double* xmom, double* ymom)
	void _manning_friction_sloped(double g, double eps, int N, double* x, double* w, double* zv, double* uh, double* vh, double* eta, double* xmom_update, double* ymom_update)

cdef inline get_python_domain(domain* D, object domain_object):

	cdef np.ndarray[long, ndim=2, mode="c"] neighbours
	cdef np.ndarray[long, ndim=2, mode="c"] neighbour_edges
	cdef np.ndarray[double, ndim=2, mode="c"] normals
	cdef np.ndarray[double, ndim=2, mode="c"] edgelengths
	cdef np.ndarray[double, ndim=1, mode="c"] radii
	cdef np.ndarray[double, ndim=1, mode="c"] areas
	cdef np.ndarray[long, ndim=1, mode="c"] edge_flux_type
	cdef np.ndarray[long, ndim=1, mode="c"] tri_full_flag
	cdef np.ndarray[long, ndim=2, mode="c"] already_computed_flux
	cdef np.ndarray[double, ndim=2, mode="c"] vertex_coordinates
	cdef np.ndarray[double, ndim=2, mode="c"] edge_coordinates
	cdef np.ndarray[double, ndim=2, mode="c"] centroid_coordinates
	cdef np.ndarray[long, ndim=1, mode="c"] number_of_boundaries
	cdef np.ndarray[long, ndim=2, mode="c"] surrogate_neighbours
	cdef np.ndarray[double, ndim=1, mode="c"] max_speed
	cdef np.ndarray[long, ndim=1, mode="c"] flux_update_frequency
	cdef np.ndarray[long, ndim=1, mode="c"] update_next_flux
	cdef np.ndarray[long, ndim=1, mode="c"] update_extrapolation
	cdef np.ndarray[long, ndim=1, mode="c"] allow_timestep_increase
	cdef np.ndarray[double, ndim=1, mode="c"] edge_timestep
	cdef np.ndarray[double, ndim=1, mode="c"] edge_flux_work
	cdef np.ndarray[double, ndim=1, mode="c"] pressuregrad_work
	cdef np.ndarray[double, ndim=1, mode="c"] x_centroid_work
	cdef np.ndarray[double, ndim=1, mode="c"] y_centroid_work
	cdef np.ndarray[double, ndim=1, mode="c"] boundary_flux_sum
	cdef np.ndarray[double, ndim=1, mode="c"] riverwall_elevation
	cdef np.ndarray[long, ndim=1, mode="c"] riverwall_rowIndex
	cdef np.ndarray[double, ndim=2, mode="c"] riverwall_hydraulic_properties

	cdef np.ndarray[double, ndim=2, mode="c"] edge_values
	cdef np.ndarray[double, ndim=1, mode="c"] centroid_values
	cdef np.ndarray[double, ndim=2, mode="c"] vertex_values
	cdef np.ndarray[double, ndim=1, mode="c"] boundary_values
	cdef np.ndarray[double, ndim=1, mode="c"] explicit_update

	cdef object quantities
	cdef object riverwallData

	D.number_of_elements = domain_object.number_of_elements
	D.epsilon = domain_object.epsilon
	D.H0 = domain_object.H0
	D.g = domain_object.g
	D.optimise_dry_cells = domain_object.optimise_dry_cells
	D.evolve_max_timestep = domain_object.evolve_max_timestep
	D.minimum_allowed_height = domain_object.minimum_allowed_height
	D.maximum_allowed_speed = domain_object.maximum_allowed_speed
	D.timestep_fluxcalls = domain_object.timestep_fluxcalls
	D.low_froude = domain_object.low_froude
	D.extrapolate_velocity_second_order = domain_object.extrapolate_velocity_second_order
	D.beta_w = domain_object.beta_w
	D.beta_w_dry = domain_object.beta_w_dry
	D.beta_uh = domain_object.beta_uh
	D.beta_uh_dry = domain_object.beta_uh_dry
	D.beta_vh = domain_object.beta_vh
	D.beta_vh_dry = domain_object.beta_vh_dry
	D.max_flux_update_frequency = domain_object.max_flux_update_frequency

	neighbours = domain_object.neighbours
	D.neighbours = &neighbours[0,0]

	surrogate_neighbours = domain_object.surrogate_neighbours
	D.surrogate_neighbours = &surrogate_neighbours[0,0]

	neighbour_edges = domain_object.neighbour_edges
	D.neighbour_edges = &neighbour_edges[0,0]

	normals = domain_object.normals
	D.normals = &normals[0,0]

	edgelengths = domain_object.edgelengths
	D.edgelengths = &edgelengths[0,0]

	radii = domain_object.radii
	D.radii = &radii[0]

	areas = domain_object.areas
	D.areas = &areas[0]

	edge_flux_type = domain_object.edge_flux_type
	D.edge_flux_type = &edge_flux_type[0]

	tri_full_flag = domain_object.tri_full_flag
	D.tri_full_flag = &tri_full_flag[0]

	already_computed_flux = domain_object.already_computed_flux
	D.already_computed_flux = &already_computed_flux[0,0]

	vertex_coordinates = domain_object.vertex_coordinates
	D.vertex_coordinates = &vertex_coordinates[0,0]

	edge_coordinates = domain_object.edge_coordinates
	D.edge_coordinates = &edge_coordinates[0,0]

	centroid_coordinates = domain_object.centroid_coordinates
	D.centroid_coordinates = &centroid_coordinates[0,0]

	max_speed = domain_object.max_speed
	D.max_speed = &max_speed[0]

	number_of_boundaries = domain_object.number_of_boundaries
	D.number_of_boundaries = &number_of_boundaries[0]

	flux_update_frequency = domain_object.flux_update_frequency
	D.flux_update_frequency = &flux_update_frequency[0]

	update_next_flux = domain_object.update_next_flux
	D.update_next_flux = &update_next_flux[0]

	update_extrapolation = domain_object.update_extrapolation
	D.update_extrapolation = &update_extrapolation[0]

	allow_timestep_increase = domain_object.allow_timestep_increase
	D.allow_timestep_increase = &allow_timestep_increase[0]

	edge_timestep = domain_object.edge_timestep
	D.edge_timestep = &edge_timestep[0]

	edge_flux_work = domain_object.edge_flux_work
	D.edge_flux_work = &edge_flux_work[0]

	pressuregrad_work = domain_object.pressuregrad_work
	D.pressuregrad_work = &pressuregrad_work[0]

	x_centroid_work = domain_object.x_centroid_work
	D.x_centroid_work = &x_centroid_work[0]

	y_centroid_work = domain_object.y_centroid_work
	D.y_centroid_work = &y_centroid_work[0]

	boundary_flux_sum = domain_object.boundary_flux_sum
	D.boundary_flux_sum = &boundary_flux_sum[0]

	quantities = domain_object.quantities

	edge_values = quantities["stage"].edge_values
	D.stage_edge_values = &edge_values[0,0]

	edge_values = quantities["xmomentum"].edge_values
	D.xmom_edge_values = &edge_values[0,0]

	edge_values = quantities["ymomentum"].edge_values
	D.ymom_edge_values = &edge_values[0,0]

	edge_values = quantities["elevation"].edge_values
	D.bed_edge_values = &edge_values[0,0]

	edge_values = quantities["height"].edge_values
	D.height_edge_values = &edge_values[0,0]

	centroid_values = quantities["stage"].centroid_values
	D.stage_centroid_values = &centroid_values[0]

	centroid_values = quantities["xmomentum"].centroid_values
	D.xmom_centroid_values = &centroid_values[0]

	centroid_values = quantities["ymomentum"].centroid_values
	D.ymom_centroid_values = &centroid_values[0]

	centroid_values = quantities["elevation"].centroid_values
	D.bed_centroid_values = &centroid_values[0]

	centroid_values = quantities["height"].centroid_values
	D.height_centroid_values = &centroid_values[0]

	vertex_values = quantities["stage"].vertex_values
	D.stage_vertex_values = &vertex_values[0,0]

	vertex_values = quantities["xmomentum"].vertex_values
	D.xmom_vertex_values = &vertex_values[0,0]

	vertex_values = quantities["ymomentum"].vertex_values
	D.ymom_vertex_values = &vertex_values[0,0]

	vertex_values = quantities["elevation"].vertex_values
	D.bed_vertex_values = &vertex_values[0,0]

	vertex_values = quantities["height"].vertex_values
	D.height_vertex_values = &vertex_values[0,0]

	boundary_values = quantities["stage"].boundary_values
	D.stage_boundary_values = &boundary_values[0]

	boundary_values = quantities["xmomentum"].boundary_values
	D.xmom_boundary_values = &boundary_values[0]

	boundary_values = quantities["ymomentum"].boundary_values
	D.ymom_boundary_values = &boundary_values[0]

	boundary_values = quantities["elevation"].boundary_values
	D.bed_boundary_values = &boundary_values[0]

	explicit_update = quantities["stage"].explicit_update
	D.stage_explicit_update = &explicit_update[0]

	explicit_update = quantities["xmomentum"].explicit_update
	D.xmom_explicit_update = &explicit_update[0]

	explicit_update = quantities["ymomentum"].explicit_update
	D.ymom_explicit_update = &explicit_update[0]

	riverwallData = domain_object.riverwallData

	riverwall_elevation = riverwallData.riverwall_elevation
	D.riverwall_elevation = &riverwall_elevation[0]

	riverwall_rowIndex = riverwallData.hydraulic_properties_rowIndex
	D.riverwall_rowIndex = &riverwall_rowIndex[0]

	D.ncol_riverwall_hydraulic_properties = riverwallData.ncol_hydraulic_properties

	riverwall_hydraulic_properties = riverwallData.hydraulic_properties
	D.riverwall_hydraulic_properties = &riverwall_hydraulic_properties[0,0]

def rotate(np.ndarray[double, ndim=1, mode="c"] q not None, np.ndarray[double, ndim=1, mode="c"] normal not None, int direction):

	assert normal.shape[0] == 2, "Normal vector must have 2 components"

	cdef np.ndarray[double, ndim=1, mode="c"] r
	cdef double n1, n2

	n1 = normal[0]
	n2 = normal[1]

	if direction == -1:
		n2 = -n2

	r = np.ascontiguousarray(np.copy(q))

	_rotate(&r[0], n1, n2)

	return r

def flux_function_central(np.ndarray[double, ndim=1, mode="c"] normal not None,\
						np.ndarray[double, ndim=1, mode="c"] ql not None,\
						np.ndarray[double, ndim=1, mode="c"] qr not None,\
						double zl,\
						double zr,\
						np.ndarray[double, ndim=1, mode="c"] edgeflux not None,\
						double epsilon,\
						double g,\
						double H0):

	cdef double h0, limiting_threshold, max_speed
	cdef int err

	h0 = H0*H0
	limiting_threshold = 10*H0

	err = _flux_function_central(&ql[0], &qr[0], zl, zr, normal[0], normal[1], epsilon, h0, limiting_threshold, g, &edgeflux[0], &max_speed)

	assert err >= 0, "Discontinuous Elevation"

	return max_speed


def compute_fluxes_ext_central_structure(object domain_object):

	cdef domain D
	cdef double timestep

	get_python_domain(&D, domain_object)

	timestep = _compute_fluxes_central_structure(&D)

	return timestep

def gravity(object domain_object):

	cdef domain D

	get_python_domain(&D, domain_object)

	err = _gravity(&D)

	if err == -1:
		return None

def gravity_wb(object domain_object):

	cdef domain D

	get_python_domain(&D, domain_object)

	err = _gravity_wb(&D)

	if err == -1:
		return None

def compute_fluxes_ext_wb(object domain_object):

	cdef domain D
	cdef double timestep

	get_python_domain(&D, domain_object)

	timestep = _compute_fluxes_central_wb(&D)

	assert timestep >= 0, "Discontinuous Elevation"

	return timestep

def compute_fluxes_ext_wb_3(object domain_object):

	cdef domain D
	cdef double timestep

	get_python_domain(&D, domain_object)

	timestep = _compute_fluxes_central_wb_3(&D)

	assert timestep >= 0, "Discontinuous Elevation"

	return timestep

def protect(double minimum_allowed_height,\
			double maximum_allowed_speed,\
			double epsilon,\
			np.ndarray[double, ndim=1, mode="c"] wc not None,\
			np.ndarray[double, ndim=1, mode="c"] zc not None,\
			np.ndarray[double, ndim=1, mode="c"] xmomc not None,\
			np.ndarray[double, ndim=1, mode="c"] ymomc not None):

	cdef int N

	N = wc.shape[0]

	_protect(N, minimum_allowed_height, maximum_allowed_speed, epsilon, &wc[0], &zc[0], &xmomc[0], &ymomc[0])

def balance_deep_and_shallow(object domain_object,\
							np.ndarray[double, ndim=1, mode="c"] wc not None,\
							np.ndarray[double, ndim=1, mode="c"] zc not None,\
							np.ndarray[double, ndim=2, mode="c"] wv not None,\
							np.ndarray[double, ndim=2, mode="c"] zv not None,\
							np.ndarray[double, ndim=1, mode="c"] hvbar not None,\
							np.ndarray[double, ndim=1, mode="c"] xmomc not None,\
							np.ndarray[double, ndim=1, mode="c"] ymomc not None,\
							np.ndarray[double, ndim=2, mode="c"] xmomv not None,\
							np.ndarray[double, ndim=2, mode="c"] ymomv not None):

	cdef double alpha_balance = 2.0
	cdef double H0
	cdef int N, tight_slope_limiters, use_centroid_velocities

	alpha_balance = domain_object.alpha_balance
	H0 = domain_object.H0
	tight_slope_limiters = domain_object.tight_slope_limiters
	use_centroid_velocities = domain_object.use_centroid_velocities

	N = wc.shape[0]

	_balance_deep_and_shallow(N,\
							&wc[0],\
							&zc[0],\
							&wv[0,0],\
							&zv[0,0],\
							&hvbar[0],\
							&xmomc[0],\
							&ymomc[0],\
							&xmomv[0,0],\
							&ymomv[0,0],\
							H0,\
							tight_slope_limiters,\
							use_centroid_velocities,\
							alpha_balance)

def manning_friction_flat(double g,\
						double eps,\
						np.ndarray[double, ndim=1, mode="c"] w not None,\
						np.ndarray[double, ndim=1, mode="c"] uh not None,\
						np.ndarray[double, ndim=1, mode="c"] vh not None,\
						np.ndarray[double, ndim=2, mode="c"] z not None,\
						np.ndarray[double, ndim=1, mode="c"] eta not None,\
						np.ndarray[double, ndim=1, mode="c"] xmom not None,\
						np.ndarray[double, ndim=1, mode="c"] ymom not None):


	cdef int N

	N = w.shape[0]

	_manning_friction_flat(g, eps, N,\
						&w[0],\
						&z[0,0],\
						&uh[0],\
						&vh[0],\
						&eta[0],\
						&xmom[0],\
						&ymom[0])

def manning_friction_sloped(double g,\
							double eps,\
							np.ndarray[double, ndim=2, mode="c"] x not None,\
							np.ndarray[double, ndim=1, mode="c"] w not None,\
							np.ndarray[double, ndim=1, mode="c"] uh not None,\
							np.ndarray[double, ndim=1, mode="c"] vh not None,\
							np.ndarray[double, ndim=2, mode="c"] z not None,\
							np.ndarray[double, ndim=1, mode="c"] eta not None,\
							np.ndarray[double, ndim=1, mode="c"] xmom not None,\
							np.ndarray[double, ndim=1, mode="c"] ymom not None):

	cdef int N

	N = w.shape[0]

	_manning_friction_sloped(g, eps, N,\
							&x[0,0],\
							&w[0],\
							&z[0,0],\
							&uh[0],\
							&vh[0],\
							&eta[0],\
							&xmom[0],\
							&ymom[0])















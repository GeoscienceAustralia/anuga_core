#cython: wraparound=False, boundscheck=True, cdivision=True, profile=False, nonecheck=False, overflowcheck=False, cdivision_warnings=False, unraisable_tracebacks=False

#wraparound=False, boundscheck=False, cdivision=True, profile=False, nonecheck=False, overflowcheck=False, cdivision_warnings=False, unraisable_tracebacks=False
import cython

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

cdef extern from "swDE_domain_openmp.c" nogil:
	struct domain:
		long number_of_elements
		long boundary_length
		long number_of_riverwall_edges
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
		double* neigh_work
		double* pressuregrad_work
		double* x_centroid_work
		double* y_centroid_work
		double* boundary_flux_sum
		long* allow_timestep_increase
		double* riverwall_elevation
		long* riverwall_rowIndex
		double* riverwall_hydraulic_properties
		long* edge_river_wall_counter


	struct edge:
		pass

	double _openmp_compute_fluxes_central(domain* D, double timestep)
	double _openmp_protect(domain* D)
	int _openmp_extrapolate_second_order_edge_sw(domain* D)
	int _openmp_fix_negative_cells(domain* D)



cdef int pointer_flag = 0
cdef int parameter_flag = 0

cdef inline get_python_domain_parameters(domain *D, object domain_object):

	D.number_of_elements = domain_object.number_of_elements
	D.boundary_length = domain_object.boundary_length 
	D.number_of_riverwall_edges = domain_object.number_of_riverwall_edges
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
		

cdef inline get_python_domain_pointers(domain *D, object domain_object):

	cdef long[:,::1]   neighbours
	cdef long[:,::1]   neighbour_edges
	cdef double[:,::1] normals
	cdef double[:,::1] edgelengths
	cdef double[::1]   radii
	cdef double[::1]   areas
	cdef long[::1]     edge_flux_type
	cdef long[::1]     tri_full_flag
	cdef long[:,::1]   already_computed_flux
	cdef double[:,::1] vertex_coordinates
	cdef double[:,::1] edge_coordinates
	cdef double[:,::1] centroid_coordinates
	cdef long[::1]     number_of_boundaries
	cdef long[:,::1]   surrogate_neighbours
	cdef double[::1]   max_speed
	cdef long[::1]     flux_update_frequency
	cdef long[::1]     update_next_flux
	cdef long[::1]     update_extrapolation
	cdef long[::1]     allow_timestep_increase
	cdef double[::1]   edge_timestep
	cdef double[::1]   edge_flux_work
	cdef double[::1]   neigh_work
	cdef double[::1]   pressuregrad_work
	cdef double[::1]   x_centroid_work
	cdef double[::1]   y_centroid_work
	cdef double[::1]   boundary_flux_sum
	cdef double[::1]   riverwall_elevation
	cdef long[::1]     riverwall_rowIndex
	cdef double[:,::1] riverwall_hydraulic_properties
	cdef long[::1]     edge_river_wall_counter
	cdef double[:,::1] edge_values
	cdef double[::1]   centroid_values
	cdef double[:,::1] vertex_values
	cdef double[::1]   boundary_values
	cdef double[::1]   explicit_update
	
	cdef object quantities
	cdef object riverwallData

	#------------------------------------------------------
	# Domain structures
	#------------------------------------------------------
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

	neigh_work = domain_object.neigh_work
	D.neigh_work = &neigh_work[0]

	pressuregrad_work = domain_object.pressuregrad_work
	D.pressuregrad_work = &pressuregrad_work[0]

	x_centroid_work = domain_object.x_centroid_work
	D.x_centroid_work = &x_centroid_work[0]

	y_centroid_work = domain_object.y_centroid_work
	D.y_centroid_work = &y_centroid_work[0]

	boundary_flux_sum = domain_object.boundary_flux_sum
	D.boundary_flux_sum = &boundary_flux_sum[0]

	edge_river_wall_counter = domain_object.edge_river_wall_counter
	D.edge_river_wall_counter  = &edge_river_wall_counter[0]

	#------------------------------------------------------
	# Quantity structures
	#------------------------------------------------------
	quantities = domain_object.quantities
	stage = quantities["stage"]
	xmomentum = quantities["xmomentum"]
	ymomentum = quantities["ymomentum"]
	elevation = quantities["elevation"]
	height = quantities["height"]

	edge_values = stage.edge_values
	D.stage_edge_values = &edge_values[0,0]

	edge_values = xmomentum.edge_values
	D.xmom_edge_values = &edge_values[0,0]

	edge_values = ymomentum.edge_values
	D.ymom_edge_values = &edge_values[0,0]

	edge_values = elevation.edge_values
	D.bed_edge_values = &edge_values[0,0]

	edge_values = height.edge_values
	D.height_edge_values = &edge_values[0,0]

	centroid_values = stage.centroid_values
	D.stage_centroid_values = &centroid_values[0]

	centroid_values = xmomentum.centroid_values
	D.xmom_centroid_values = &centroid_values[0]

	centroid_values = ymomentum.centroid_values
	D.ymom_centroid_values = &centroid_values[0]

	centroid_values = elevation.centroid_values
	D.bed_centroid_values = &centroid_values[0]

	centroid_values = height.centroid_values
	D.height_centroid_values = &centroid_values[0]

	vertex_values = stage.vertex_values
	D.stage_vertex_values = &vertex_values[0,0]

	vertex_values = xmomentum.vertex_values
	D.xmom_vertex_values = &vertex_values[0,0]

	vertex_values = ymomentum.vertex_values
	D.ymom_vertex_values = &vertex_values[0,0]

	vertex_values = elevation.vertex_values
	D.bed_vertex_values = &vertex_values[0,0]

	vertex_values = height.vertex_values
	D.height_vertex_values = &vertex_values[0,0]

	boundary_values = stage.boundary_values
	D.stage_boundary_values = &boundary_values[0]

	boundary_values = xmomentum.boundary_values
	D.xmom_boundary_values = &boundary_values[0]

	boundary_values = ymomentum.boundary_values
	D.ymom_boundary_values = &boundary_values[0]

	boundary_values = elevation.boundary_values
	D.bed_boundary_values = &boundary_values[0]

	explicit_update = stage.explicit_update
	D.stage_explicit_update = &explicit_update[0]

	explicit_update = xmomentum.explicit_update
	D.xmom_explicit_update = &explicit_update[0]

	explicit_update = ymomentum.explicit_update
	D.ymom_explicit_update = &explicit_update[0]

	#------------------------------------------------------
	# Riverwall structures
	#------------------------------------------------------
	riverwallData = domain_object.riverwallData

	riverwall_elevation = riverwallData.riverwall_elevation
	D.riverwall_elevation = &riverwall_elevation[0]

	riverwall_rowIndex = riverwallData.hydraulic_properties_rowIndex
	D.riverwall_rowIndex = &riverwall_rowIndex[0]

	D.ncol_riverwall_hydraulic_properties = riverwallData.ncol_hydraulic_properties

	riverwall_hydraulic_properties = riverwallData.hydraulic_properties
	D.riverwall_hydraulic_properties = &riverwall_hydraulic_properties[0,0]



#===============================================================================

def compute_fluxes_ext_central(object domain_object, double timestep):

	cdef domain D

	get_python_domain_parameters(&D, domain_object)
	get_python_domain_pointers(&D, domain_object)

	with nogil:
		timestep =  _openmp_compute_fluxes_central(&D, timestep)

	return timestep

def extrapolate_second_order_edge_sw(object domain_object):

	cdef domain D
	cdef int e

	get_python_domain_parameters(&D, domain_object)
	get_python_domain_pointers(&D, domain_object)

	with nogil:
		e = _openmp_extrapolate_second_order_edge_sw(&D)

	if e == -1:
		return None

def protect_new(object domain_object):

	cdef domain D

	cdef double mass_error

	get_python_domain_parameters(&D, domain_object)
	get_python_domain_pointers(&D, domain_object)

	with nogil:
		mass_error = _openmp_protect(&D)


	return mass_error

def compute_flux_update_frequency(object domain_object, double timestep):

	pass

def fix_negative_cells(object domain_object):

	cdef domain D
	cdef int num_negative_cells

	get_python_domain_parameters(&D, domain_object)
	get_python_domain_pointers(&D, domain_object)

	with nogil:
		num_negative_cells = _openmp_fix_negative_cells(&D)

	return num_negative_cells




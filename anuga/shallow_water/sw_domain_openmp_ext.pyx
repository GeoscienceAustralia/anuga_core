#cython: wraparound=False, boundscheck=False, cdivision=True, profile=False, nonecheck=False, overflowcheck=False, cdivision_warnings=False, unraisable_tracebacks=False, warn.unused=False

import cython
from libc.stdint cimport int64_t as anuga_int
from libc.stdint cimport uintptr_t

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

# FIXME SR: Should create a python class which holds all of these parameters
# so we don't need to create it before each C function call.

cdef extern from "sw_domain_openmp.c" nogil:
	struct domain:
		# these shouldn't change within a single timestep
		# they are set once before the evolve loop

		anuga_int number_of_elements
		anuga_int boundary_length
		anuga_int number_of_riverwall_edges
		anuga_int number_of_tracers
		double beta_tracer
		anuga_int n_sediment_classes
		double* sediment_settling_velocity
		double* sediment_d_star
		double sediment_c_max
		double* sediment_diameter
		double* sediment_R
		double* sediment_tau_c_star
		double sediment_gamma0
		anuga_int sediment_erosion_mode
		double sediment_K_partheniades
		double sediment_dl_A
		double sediment_dl_alpha
		double sediment_dl_beta
		double sediment_dl_threshold
		double sediment_dl_ks
		double sediment_dl_ks_factor
		anuga_int sediment_vegetation_shear
		anuga_int sediment_deposition_mode
		double sediment_tau_d
		double sediment_tau_crit
		double sediment_K_e
		double sediment_rho_w
		anuga_int sediment_shear_closure
		double sediment_max_slope
		anuga_int sediment_slope_frozen
		double* tracer_external_source
		anuga_int sediment_d_star_mode
		double* sediment_reference_height
		double sediment_a_h_floor
		anuga_int sediment_rouse_beta_mode
		double sediment_rouse_scale
		anuga_int sediment_adaptation_mode
		double sediment_adaptation_alpha
		anuga_int sediment_nearbed_base
		double sediment_layer_fraction
		double sediment_exchange_factor
		anuga_int sediment_velocity_profile
		double sediment_porosity
		double sediment_morphological_factor
		anuga_int sediment_bed_evolution
		anuga_int sediment_bedload_mode
		double sediment_bedload_K
		double sediment_bedload_m
		double sediment_bedload_tau_c_star
		double sediment_bedload_h_min
		double* sediment_z_base
		double* sediment_shear_factor
		anuga_int sediment_has_z_base
		anuga_int sediment_bed_composition
		double sediment_active_layer
		double* sediment_bed_active
		double* sediment_bed_substrate
		double* sediment_repose_dz
		double sediment_repose_tan
		double sediment_repose_relax
		anuga_int sediment_repose_max_sweeps
		double* sediment_source_limited
		double* sediment_slope_work
		anuga_int* sediment_bed_exhausted
		anuga_int* sediment_bedload_open
		double* sediment_bedload_supply
		double* sediment_qbx
		double* sediment_qby
		double* sediment_qba
		double sediment_c_pack
		anuga_int vegetation_mode
		anuga_int use_sloped_mannings
		double vegetation_Cd
		double vegetation_bed_chezy
		double* veg_density_centroid_values
		double* veg_diameter_centroid_values
		double* veg_height_centroid_values
		anuga_int sediment_friction_mode
		double sediment_manning_ll
		anuga_int sediment_wilson_bed
		double sediment_wilson_D
		double* tracer_centroid_values
		double* tracer_edge_values
		double* tracer_boundary_values
		double* tracer_explicit_update
		double* tracer_conserved_values
		double* tracer_backup_values
		double* tracer_boundary_flux
		double* tracer_speed_factor
		double* tracer_boundary_flux_sum
		anuga_int optimise_dry_cells
		anuga_int extrapolate_velocity_second_order
		anuga_int low_froude
		anuga_int timestep_fluxcalls
		anuga_int flux_call_count
		anuga_int flux_base_call
		anuga_int flux_timestep_fluxcalls_seen
		anuga_int reconstruct_edge_bed
		anuga_int num_owned_edges
		anuga_int* owned_edges
		anuga_int ncol_riverwall_hydraulic_properties
		anuga_int nrow_riverwall_hydraulic_properties

		double epsilon
		double H0
		double g
		double evolve_max_timestep
		double evolve_min_timestep
		double minimum_allowed_height
		double maximum_allowed_speed
		double beta_w
		double beta_w_dry
		double beta_uh
		double beta_uh_dry
		double beta_vh
		double beta_vh_dry


		# these pointers are set once the domain is created
		# and are used in the evolve loop
		anuga_int* neighbours
		anuga_int* neighbour_edges
		anuga_int* surrogate_neighbours
		double* normals
		double* edgelengths
		double* radii
		double* areas
		anuga_int* edge_flux_type
		anuga_int* tri_full_flag
		anuga_int* already_computed_flux
		double* max_speed
		double* vertex_coordinates
		double* edge_coordinates
		double* centroid_coordinates
		anuga_int* number_of_boundaries
		double* stage_edge_values
		double* xmom_edge_values
		double* ymom_edge_values
		double* bed_edge_values
		double* height_edge_values
		double* xvelocity_edge_values
		double* yvelocity_edge_values
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
		double* height_boundary_values
		double* xvelocity_boundary_values
		double* yvelocity_boundary_values
		double* stage_explicit_update
		double* xmom_explicit_update
		double* ymom_explicit_update
		double* edge_flux_work
		double* neigh_work
		double* pressuregrad_work
		double* x_centroid_work
		double* y_centroid_work
		double* boundary_flux_sum
		anuga_int* edge_river_wall_counter
		double* riverwall_elevation
		anuga_int* riverwall_rowIndex
		double* riverwall_hydraulic_properties

		double* stage_semi_implicit_update
		double* xmom_semi_implicit_update
		double* ymom_semi_implicit_update
		double* friction_centroid_values

		double* stage_backup_values
		double* xmom_backup_values
		double* ymom_backup_values


	struct edge:
		pass

	anuga_int __rotate(double *q, double n1, double n2)
	void _openmp_set_omp_num_threads(anuga_int num_threads)
	double _openmp_compute_fluxes_central(domain* D, double timestep)
	double _openmp_protect(domain* D)
	void _openmp_extrapolate_second_order_edge_sw(domain* D)
	anuga_int _openmp_fix_negative_cells(domain* D)
	double _openmp_negative_cells_volume(domain* D)
	anuga_int _openmp_gravity(domain *D)
	anuga_int _openmp_gravity_wb(domain *D) 
	anuga_int _openmp_update_conserved_quantities(domain* D, double timestep)
	void _openmp_manning_friction_flat_semi_implicit(domain *D)
	void _openmp_manning_friction_sloped_semi_implicit(domain *D)
	void _openmp_vegetation_friction_semi_implicit(domain *D)
	void _openmp_manning_friction_sloped_semi_implicit_edge_based(domain *D)
	anuga_int _openmp_saxpy_conserved_quantities(domain *D, double a, double b, double c)
	anuga_int _openmp_backup_conserved_quantities(domain *D)
	void _openmp_distribute_edges_to_vertices(domain *D)
	void _openmp_ader_ck_predictor(domain *D, double dt)
	void _openmp_ader_ck_predictor_edge(domain *D, double dt)
	# FIXME SR: Change over to domain* D argument ?
	void _openmp_manning_friction_flat(double g, double eps, anuga_int N, double* w, double* zv, double* uh, double* vh, double* eta, double* xmom, double* ymom)
	void _openmp_manning_friction_sloped(double g, double eps, anuga_int N, double* x, double* w, double* zv, double* uh, double* vh, double* eta, double* xmom_update, double* ymom_update)
	void _openmp_manning_friction_sloped_edge_based(double g, double eps, anuga_int N, double* x, double* w, double* zv, double* uh, double* vh, double* eta, double* xmom_update, double* ymom_update)
	void _openmp_evaluate_reflective_segment(domain *D, anuga_int N, anuga_int *edge_ptr, anuga_int *vol_ids_ptr, anuga_int *edge_ids_ptr)
	anuga_int __openmp__flux_function_central(double q_left0, double q_left1, double q_left2,
	double q_right0, double q_right1, double q_right2,
	double h_left, double h_right, double hle, double hre,
	double n1, double n2, double epsilon, double ze, double g,
	double* edgeflux0, double* edgeflux1, double* edgeflux2,
	double* max_speed, double* pressure_flux, anuga_int low_froude)


cdef anuga_int pointer_flag = 0
cdef anuga_int parameter_flag = 0

cdef inline get_python_domain_parameters(domain *D, object domain_py_object):

	riverwallData = domain_py_object.riverwallData

	# these shouldn't change within a single timestep

	D.number_of_elements = domain_py_object.number_of_elements
	D.boundary_length = domain_py_object.boundary_length
	D.number_of_riverwall_edges = domain_py_object.number_of_riverwall_edges
	# Generic tracers. The struct is PyMem_Malloc'd (not zeroed), so this
	# MUST be set on every fill or the flux kernel guard reads garbage.
	# Sediment: MUST be set explicitly in BOTH pyx files. struct domain is
	# PyMem_Malloc'd and never zeroed, and a garbage n_sediment_classes makes
	# the source kernel walk unmapped pointers (HANDOVER.md 2.4).
	D.n_sediment_classes = getattr(domain_py_object, 'n_sediment_classes', 0)
	D.sediment_c_max = getattr(domain_py_object, 'sediment_c_max', 0.3)
	D.sediment_gamma0 = getattr(domain_py_object, 'sediment_gamma0', 0.0024)
	D.sediment_erosion_mode = getattr(domain_py_object, 'sediment_erosion_mode', 0)
	D.sediment_K_partheniades = getattr(domain_py_object, 'sediment_K_partheniades', 1.0e-4)
	D.sediment_dl_A = getattr(domain_py_object, 'sediment_dl_A', 0.000474)
	D.sediment_dl_alpha = getattr(domain_py_object, 'sediment_dl_alpha', 1.5)
	D.sediment_dl_beta = getattr(domain_py_object, 'sediment_dl_beta', 1.18)
	D.sediment_dl_threshold = getattr(domain_py_object, 'sediment_dl_threshold', 0.015)
	D.sediment_dl_ks = getattr(domain_py_object, 'sediment_dl_ks', 0.0)
	D.sediment_dl_ks_factor = getattr(domain_py_object, 'sediment_dl_ks_factor', 3.0)
	D.sediment_vegetation_shear = getattr(domain_py_object, 'sediment_vegetation_shear', 1)
	D.sediment_deposition_mode = getattr(domain_py_object, 'sediment_deposition_mode', 0)
	D.sediment_tau_d = getattr(domain_py_object, 'sediment_tau_d', 0.0)
	D.sediment_tau_crit = getattr(domain_py_object, 'sediment_tau_crit', 0.088)
	D.sediment_K_e = getattr(domain_py_object, 'sediment_K_e', 0.2e-6/0.088**0.5)
	D.sediment_rho_w = getattr(domain_py_object, 'sediment_rho_w', 1000.0)
	D.sediment_shear_closure = getattr(domain_py_object, 'sediment_shear_closure', 0)
	D.sediment_max_slope = getattr(domain_py_object, 'sediment_max_slope', 0.0)
	D.sediment_slope_frozen = getattr(domain_py_object, 'sediment_slope_frozen', 0)
	D.sediment_d_star_mode = getattr(domain_py_object, 'sediment_d_star_mode', 0)
	D.sediment_a_h_floor = getattr(domain_py_object, 'sediment_a_h_floor', 0.01)
	D.sediment_adaptation_mode = getattr(domain_py_object, 'sediment_adaptation_mode', 0)
	D.sediment_rouse_beta_mode = getattr(domain_py_object, 'sediment_rouse_beta_mode', 0)
	D.sediment_rouse_scale = getattr(domain_py_object, 'sediment_rouse_scale', 1.0)
	D.sediment_adaptation_alpha = getattr(domain_py_object, 'sediment_adaptation_alpha', 1.0)
	D.sediment_nearbed_base = getattr(domain_py_object, 'sediment_nearbed_base', -1)
	D.sediment_layer_fraction = getattr(domain_py_object, 'sediment_layer_fraction', 0.0)
	D.sediment_exchange_factor = getattr(domain_py_object, 'sediment_exchange_factor', 1.0)
	D.sediment_velocity_profile = getattr(domain_py_object, 'sediment_velocity_profile', 0)
	D.sediment_c_pack = getattr(domain_py_object, 'sediment_c_pack', 0.65)
	D.sediment_porosity = getattr(domain_py_object, 'sediment_porosity', 0.3)
	D.sediment_morphological_factor = getattr(domain_py_object, 'sediment_morphological_factor', 1.0)
	D.sediment_bed_evolution = 1 if getattr(domain_py_object, 'sediment_bed_evolution', True) else 0
	D.sediment_bedload_mode = getattr(domain_py_object, 'sediment_bedload_mode', 0)
	# [L-5]. Set UNCONDITIONALLY, outside the n_sediment_classes guard
	# below: the kernels test this flag before dereferencing
	# sediment_z_base, so leaving it uninitialised makes a NULL
	# pointer look like a configured base.
	D.sediment_has_z_base = getattr(domain_py_object, 'sediment_has_z_base', 0)
	# Bed composition: flag and thickness unconditionally, like the base;
	# the kernel also requires both arrays before it dereferences them.
	D.sediment_bed_composition = getattr(domain_py_object, 'sediment_bed_composition', 0)
	D.sediment_active_layer = getattr(domain_py_object, 'sediment_active_layer', 0.0)
	# spec 7 repose. Scalars again set unconditionally -- the kernel
	# tests sediment_repose_tan before touching sediment_repose_dz.
	D.sediment_repose_tan = getattr(domain_py_object, 'sediment_repose_tan', 0.0)
	D.sediment_repose_relax = getattr(domain_py_object, 'sediment_repose_relax', 1.0)
	D.sediment_repose_max_sweeps = getattr(domain_py_object, 'sediment_repose_max_sweeps', 0)
	D.sediment_bedload_K = getattr(domain_py_object, 'sediment_bedload_K', 3.97)
	D.sediment_bedload_m = getattr(domain_py_object, 'sediment_bedload_m', 1.5)
	D.sediment_bedload_tau_c_star = getattr(domain_py_object, 'sediment_bedload_tau_c_star', 0.0495)
	D.sediment_bedload_h_min = getattr(domain_py_object, 'sediment_bedload_h_min', 0.0)
	cdef double[::1] veg1
	D.vegetation_mode = getattr(domain_py_object, 'vegetation_mode', 0)
	D.use_sloped_mannings = 1 if getattr(domain_py_object, 'use_sloped_mannings', False) else 0
	D.vegetation_Cd = getattr(domain_py_object, 'vegetation_Cd', 1.68)
	D.vegetation_bed_chezy = getattr(domain_py_object, 'vegetation_bed_chezy', 65.0)
	if D.vegetation_mode > 0:
		veg1 = domain_py_object.quantities['veg_density'].centroid_values
		D.veg_density_centroid_values = &veg1[0]
		veg1 = domain_py_object.quantities['veg_diameter'].centroid_values
		D.veg_diameter_centroid_values = &veg1[0]
		veg1 = domain_py_object.quantities['veg_height'].centroid_values
		D.veg_height_centroid_values = &veg1[0]
	else:
		D.veg_density_centroid_values = NULL
		D.veg_diameter_centroid_values = NULL
		D.veg_height_centroid_values = NULL
	D.sediment_friction_mode = getattr(domain_py_object, 'sediment_friction_mode', 0)
	D.sediment_manning_ll = getattr(domain_py_object, 'sediment_manning_ll', 0.065)
	D.sediment_wilson_bed = getattr(domain_py_object, 'sediment_wilson_bed', 0)
	D.sediment_wilson_D = getattr(domain_py_object, 'sediment_wilson_D', 1.0e-3)
	D.number_of_tracers = getattr(domain_py_object, 'number_of_tracers', 0)
	# 0.0 => first order; >0 => limited second order (see sw_domain.h)
	D.beta_tracer = getattr(domain_py_object, 'beta_tracer', 1.0)
	D.optimise_dry_cells = domain_py_object.optimise_dry_cells
	D.extrapolate_velocity_second_order = domain_py_object.extrapolate_velocity_second_order
	D.low_froude = domain_py_object.low_froude
	D.timestep_fluxcalls = domain_py_object.timestep_fluxcalls
	# The C struct comes from PyMem_Malloc (uninitialized): every field must be
	# set here.  Bed-edge reconstruction in compute_fluxes stays OFF on this
	# path -- tests call compute_fluxes directly with independent bed_ev.
	D.reconstruct_edge_bed = 0
	D.num_owned_edges = 0
	D.owned_edges = NULL

	D.ncol_riverwall_hydraulic_properties = riverwallData.ncol_hydraulic_properties
	try:
		D.nrow_riverwall_hydraulic_properties = len(riverwallData.names)
	except:
		D.nrow_riverwall_hydraulic_properties = 0

	D.epsilon = domain_py_object.epsilon
	D.H0 = domain_py_object.H0
	D.g = domain_py_object.g
	D.evolve_max_timestep = domain_py_object.evolve_max_timestep
	D.evolve_min_timestep = domain_py_object.evolve_min_timestep
	D.minimum_allowed_height = domain_py_object.minimum_allowed_height
	D.maximum_allowed_speed = domain_py_object.maximum_allowed_speed
	D.beta_w = domain_py_object.beta_w
	D.beta_w_dry = domain_py_object.beta_w_dry
	D.beta_uh = domain_py_object.beta_uh
	D.beta_uh_dry = domain_py_object.beta_uh_dry
	D.beta_vh = domain_py_object.beta_vh
	D.beta_vh_dry = domain_py_object.beta_vh_dry
		

cdef inline _check_shape(object arr, tuple expected, str name):
	"""Raise if a Python-side array is not the shape the C struct assumes.

	The C kernels index these arrays with raw pointers and no bounds
	checks, so a wrong shape here is an out-of-bounds read or write in the
	flux loop. This runs once per struct (re)fill, not per element.

	None is allowed through: quantities with reduced storage (see
	Quantity.__init__) have no edge, vertex or update arrays, and the
	kernels that run for such domains never touch those pointers. The
	check is about arrays that exist with the wrong shape.
	"""
	if arr is None:
		return
	cdef tuple got = tuple(arr.shape)
	if got != expected:
		raise ValueError(
			"domain.%s has shape %s but the C domain struct expects %s"
			% (name, got, expected))

cdef inline get_python_domain_pointers(domain *D, object domain_py_object):

	# these arraypointers are set once the domain is created
	# and are used in the evolve loop
	cdef anuga_int[:,::1]   neighbours
	cdef anuga_int[:,::1]   neighbour_edges
	cdef double[:,::1] normals
	cdef double[:,::1] edgelengths
	cdef double[::1]   radii
	cdef double[::1]   areas
	cdef anuga_int[::1]  edge_flux_type
	cdef anuga_int[::1]  tri_full_flag
	cdef anuga_int[:,::1] already_computed_flux
	cdef double[:,::1] vertex_coordinates
	cdef double[:,::1] edge_coordinates
	cdef double[:,::1] centroid_coordinates
	cdef anuga_int[::1]  number_of_boundaries
	cdef anuga_int[:,::1] surrogate_neighbours
	cdef double[::1]   max_speed
	cdef double[::1]   edge_flux_work
	cdef double[::1]   neigh_work
	cdef double[::1]   pressuregrad_work
	cdef double[::1]   x_centroid_work
	cdef double[::1]   y_centroid_work
	cdef double[::1]   boundary_flux_sum
	cdef double[::1]   riverwall_elevation
	cdef anuga_int[::1]  riverwall_rowIndex
	cdef double[:,::1] riverwall_hydraulic_properties
	cdef anuga_int[::1]  edge_river_wall_counter
	cdef double[:,::1] edge_values
	cdef double[::1]   centroid_values
	cdef double[:,::1] vertex_values
	cdef double[::1]   boundary_values
	cdef double[::1]   explicit_update
	cdef double[::1]   semi_implicit_update
	
	cdef object quantities
	cdef object riverwallData

	# Sizes the struct was filled with (get_python_domain_parameters runs
	# first). Every array below is checked against them.
	cdef anuga_int N = D.number_of_elements
	cdef anuga_int B = D.boundary_length

	#------------------------------------------------------
	# Domain structures
	#------------------------------------------------------
	_check_shape(domain_py_object.neighbours, (N, 3), "neighbours")
	_check_shape(domain_py_object.neighbour_edges, (N, 3), "neighbour_edges")
	_check_shape(domain_py_object.surrogate_neighbours, (N, 3), "surrogate_neighbours")
	_check_shape(domain_py_object.normals, (N, 6), "normals")
	_check_shape(domain_py_object.edgelengths, (N, 3), "edgelengths")
	_check_shape(domain_py_object.radii, (N,), "radii")
	_check_shape(domain_py_object.areas, (N,), "areas")
	_check_shape(domain_py_object.tri_full_flag, (N,), "tri_full_flag")
	_check_shape(domain_py_object.vertex_coordinates, (3 * N, 2), "vertex_coordinates")
	_check_shape(domain_py_object.edge_coordinates, (3 * N, 2), "edge_coordinates")
	_check_shape(domain_py_object.centroid_coordinates, (N, 2), "centroid_coordinates")
	_check_shape(domain_py_object.max_speed, (N,), "max_speed")
	_check_shape(domain_py_object.number_of_boundaries, (N,), "number_of_boundaries")
	_check_shape(domain_py_object.x_centroid_work, (N,), "x_centroid_work")
	_check_shape(domain_py_object.y_centroid_work, (N,), "y_centroid_work")
	if domain_py_object.edge_flux_type is not None:
		_check_shape(domain_py_object.edge_flux_type, (3 * N,), "edge_flux_type")

	neighbours = domain_py_object.neighbours
	D.neighbours = &neighbours[0,0]
	
	surrogate_neighbours = domain_py_object.surrogate_neighbours
	D.surrogate_neighbours = &surrogate_neighbours[0,0]

	neighbour_edges = domain_py_object.neighbour_edges
	D.neighbour_edges = &neighbour_edges[0,0]

	normals = domain_py_object.normals
	D.normals = &normals[0,0]

	edgelengths = domain_py_object.edgelengths
	D.edgelengths = &edgelengths[0,0]

	radii = domain_py_object.radii
	D.radii = &radii[0]

	areas = domain_py_object.areas
	D.areas = &areas[0]

	if domain_py_object.edge_flux_type is not None:
		edge_flux_type = domain_py_object.edge_flux_type
		D.edge_flux_type = &edge_flux_type[0]
	else:
		D.edge_flux_type = NULL

	tri_full_flag = domain_py_object.tri_full_flag
	D.tri_full_flag = &tri_full_flag[0]

	if domain_py_object.already_computed_flux is not None:
		already_computed_flux = domain_py_object.already_computed_flux
		D.already_computed_flux = &already_computed_flux[0,0]
	else:
		D.already_computed_flux = NULL

	vertex_coordinates = domain_py_object.vertex_coordinates
	D.vertex_coordinates = &vertex_coordinates[0,0]

	edge_coordinates = domain_py_object.edge_coordinates
	D.edge_coordinates = &edge_coordinates[0,0]

	centroid_coordinates = domain_py_object.centroid_coordinates
	D.centroid_coordinates = &centroid_coordinates[0,0]

	max_speed = domain_py_object.max_speed
	D.max_speed = &max_speed[0]

	number_of_boundaries = domain_py_object.number_of_boundaries
	D.number_of_boundaries = &number_of_boundaries[0]

	if domain_py_object.edge_flux_work is not None:
		edge_flux_work = domain_py_object.edge_flux_work
		D.edge_flux_work = &edge_flux_work[0]
	else:
		D.edge_flux_work = NULL

	if domain_py_object.neigh_work is not None:
		neigh_work = domain_py_object.neigh_work
		D.neigh_work = &neigh_work[0]
	else:
		D.neigh_work = NULL

	if domain_py_object.pressuregrad_work is not None:
		pressuregrad_work = domain_py_object.pressuregrad_work
		D.pressuregrad_work = &pressuregrad_work[0]
	else:
		D.pressuregrad_work = NULL

	x_centroid_work = domain_py_object.x_centroid_work
	D.x_centroid_work = &x_centroid_work[0]

	y_centroid_work = domain_py_object.y_centroid_work
	D.y_centroid_work = &y_centroid_work[0]

	boundary_flux_sum = domain_py_object.boundary_flux_sum
	D.boundary_flux_sum = &boundary_flux_sum[0]

	if domain_py_object.edge_river_wall_counter is not None:
		edge_river_wall_counter = domain_py_object.edge_river_wall_counter
		D.edge_river_wall_counter = &edge_river_wall_counter[0]
	else:
		D.edge_river_wall_counter = NULL

	#------------------------------------------------------
	# Quantity structures
	#------------------------------------------------------
	# Generic tracer arrays. Owned by the Python Domain as C-contiguous
	# (ns, ...) float64 arrays; NULL when no tracers are registered.
	cdef double[:, ::1] tr2
	cdef double[::1] sed1
	cdef double[:,::1] sed2
	cdef anuga_int[::1] sedi
	cdef double[::1] tr1
	if getattr(domain_py_object, 'number_of_tracers', 0) > 0:
		ns = D.number_of_tracers
		_check_shape(domain_py_object.tracer_centroid_values, (ns, N), "tracer_centroid_values")
		_check_shape(domain_py_object.tracer_edge_values, (ns, 3 * N), "tracer_edge_values")
		_check_shape(domain_py_object.tracer_boundary_values, (ns, B), "tracer_boundary_values")
		_check_shape(domain_py_object.tracer_explicit_update, (ns, N), "tracer_explicit_update")
		_check_shape(domain_py_object.tracer_conserved_values, (ns, N), "tracer_conserved_values")
		_check_shape(domain_py_object.tracer_backup_values, (ns, N), "tracer_backup_values")
		if domain_py_object.tracer_external_source is not None:
			_check_shape(domain_py_object.tracer_external_source, (ns, N), "tracer_external_source")
		tr2 = domain_py_object.tracer_centroid_values
		D.tracer_centroid_values = &tr2[0, 0]
		tr2 = domain_py_object.tracer_edge_values
		D.tracer_edge_values = &tr2[0, 0]
		tr2 = domain_py_object.tracer_boundary_values
		D.tracer_boundary_values = &tr2[0, 0]
		tr2 = domain_py_object.tracer_explicit_update
		D.tracer_explicit_update = &tr2[0, 0]
		tr2 = domain_py_object.tracer_conserved_values
		D.tracer_conserved_values = &tr2[0, 0]
		tr2 = domain_py_object.tracer_backup_values
		D.tracer_backup_values = &tr2[0, 0]
		if domain_py_object.tracer_external_source is not None:
			tr2 = domain_py_object.tracer_external_source
			D.tracer_external_source = &tr2[0, 0]
		else:
			D.tracer_external_source = NULL
		if D.n_sediment_classes > 0:
			sed1 = domain_py_object.sediment_settling_velocity
			D.sediment_settling_velocity = &sed1[0]
			sed1 = domain_py_object.sediment_d_star
			D.sediment_d_star = &sed1[0]
			sed1 = domain_py_object.sediment_diameter
			D.sediment_diameter = &sed1[0]
			sed1 = domain_py_object.sediment_R
			D.sediment_R = &sed1[0]
			sed1 = domain_py_object.sediment_tau_c_star
			D.sediment_tau_c_star = &sed1[0]
			sed1 = domain_py_object.sediment_reference_height
			D.sediment_reference_height = &sed1[0]
			sed1 = domain_py_object.sediment_qbx
			D.sediment_qbx = &sed1[0]
			sed1 = domain_py_object.sediment_qby
			D.sediment_qby = &sed1[0]
			sed1 = domain_py_object.sediment_qba
			D.sediment_qba = &sed1[0]
			# [L-5]. sediment_source_limited is dereferenced by the source
			# kernel whenever there is a class, so it is bound unconditionally
			# alongside them; z_base only when a base was actually set.
			sed2 = domain_py_object.sediment_source_limited
			D.sediment_source_limited = &sed2[0,0]
			sed1 = domain_py_object.sediment_slope_work
			D.sediment_slope_work = &sed1[0]
			sedi = domain_py_object.sediment_bed_exhausted
			D.sediment_bed_exhausted = &sedi[0]
			sed1 = domain_py_object.sediment_repose_dz
			D.sediment_repose_dz = &sed1[0]
			sedi = domain_py_object.sediment_bedload_open
			if sedi.shape[0] > 0:
				D.sediment_bedload_open = &sedi[0]
			else:
				D.sediment_bedload_open = NULL
				D.sediment_bedload_supply = NULL
			sed1 = domain_py_object.sediment_bedload_supply
			if sed1.shape[0] > 0:
				D.sediment_bedload_supply = &sed1[0]
			else:
				D.sediment_bedload_supply = NULL
			if getattr(domain_py_object, 'sediment_bed_composition', 0):
				sed2 = domain_py_object.sediment_bed_active
				D.sediment_bed_active = &sed2[0,0]
				sed2 = domain_py_object.sediment_bed_substrate
				D.sediment_bed_substrate = &sed2[0,0]
			else:
				D.sediment_bed_active = NULL
				D.sediment_bed_substrate = NULL
			if domain_py_object.sediment_has_z_base:
				sed1 = domain_py_object.sediment_z_base
				D.sediment_z_base = &sed1[0]
			else:
				D.sediment_z_base = NULL
			sfa = getattr(domain_py_object, 'sediment_shear_factor', None)
			if sfa is not None:
				sed1 = sfa
				D.sediment_shear_factor = &sed1[0]
			else:
				D.sediment_shear_factor = NULL
		else:
			D.sediment_settling_velocity = NULL
			D.sediment_d_star = NULL
			D.sediment_diameter = NULL
			D.sediment_R = NULL
			D.sediment_tau_c_star = NULL
			D.sediment_reference_height = NULL
			D.sediment_qbx = NULL
			D.sediment_qby = NULL
			D.sediment_qba = NULL
			D.sediment_z_base = NULL
			D.sediment_shear_factor = NULL
			D.sediment_source_limited = NULL
			D.sediment_bed_active = NULL
			D.sediment_bed_substrate = NULL
			D.sediment_slope_work = NULL
			D.sediment_bed_exhausted = NULL
			D.sediment_repose_dz = NULL
			D.sediment_bedload_open = NULL
			D.sediment_bedload_supply = NULL
		tr2 = domain_py_object.tracer_boundary_flux
		D.tracer_boundary_flux = &tr2[0, 0]
		tsf = getattr(domain_py_object, 'tracer_speed_factor', None)
		if tsf is not None:
			tr2 = tsf
			D.tracer_speed_factor = &tr2[0, 0]
		else:
			D.tracer_speed_factor = NULL
		tr1 = domain_py_object.tracer_boundary_flux_sum
		D.tracer_boundary_flux_sum = &tr1[0]
	else:
		D.sediment_settling_velocity = NULL
		D.sediment_d_star = NULL
		D.sediment_diameter = NULL
		D.sediment_R = NULL
		D.sediment_tau_c_star = NULL
		D.sediment_reference_height = NULL
		D.sediment_qbx = NULL
		D.sediment_qby = NULL
		D.sediment_qba = NULL
		D.sediment_z_base = NULL
		D.sediment_source_limited = NULL
		D.sediment_bed_active = NULL
		D.sediment_bed_substrate = NULL
		D.sediment_slope_work = NULL
		D.sediment_bed_exhausted = NULL
		D.sediment_repose_dz = NULL
		D.sediment_bedload_open = NULL
		D.sediment_bedload_supply = NULL
		D.tracer_centroid_values = NULL
		D.tracer_edge_values = NULL
		D.tracer_boundary_values = NULL
		D.tracer_explicit_update = NULL
		D.tracer_conserved_values = NULL
		D.tracer_backup_values = NULL
		D.tracer_external_source = NULL
		D.tracer_boundary_flux = NULL
		D.tracer_speed_factor = NULL
		D.tracer_boundary_flux_sum = NULL

	quantities = domain_py_object.quantities
	for _name in ("stage", "xmomentum", "ymomentum", "elevation", "height",
	              "friction", "xvelocity", "yvelocity"):
		_q = quantities[_name]
		_check_shape(_q.centroid_values, (N,), _name + ".centroid_values")
		_check_shape(_q.edge_values, (N, 3), _name + ".edge_values")
		_check_shape(_q.vertex_values, (N, 3), _name + ".vertex_values")
		_check_shape(_q.boundary_values, (B,), _name + ".boundary_values")
	for _name in ("stage", "xmomentum", "ymomentum"):
		_q = quantities[_name]
		_check_shape(_q.explicit_update, (N,), _name + ".explicit_update")
		_check_shape(_q.semi_implicit_update, (N,), _name + ".semi_implicit_update")
		_check_shape(_q.centroid_backup_values, (N,), _name + ".centroid_backup_values")
	stage = quantities["stage"]
	xmomentum = quantities["xmomentum"]
	ymomentum = quantities["ymomentum"]
	elevation = quantities["elevation"]
	height = quantities["height"]
	friction = quantities["friction"]
	xvelocity = quantities["xvelocity"]
	yvelocity = quantities["yvelocity"]

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

	edge_values = xvelocity.edge_values
	D.xvelocity_edge_values = &edge_values[0,0]

	edge_values = yvelocity.edge_values
	D.yvelocity_edge_values = &edge_values[0,0]

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

	centroid_values = friction.centroid_values
	D.friction_centroid_values = &centroid_values[0]	

	centroid_values = stage.centroid_backup_values
	D.stage_backup_values = &centroid_values[0]	
	
	centroid_values = xmomentum.centroid_backup_values
	D.xmom_backup_values = &centroid_values[0]		
	
	centroid_values = ymomentum.centroid_backup_values
	D.ymom_backup_values = &centroid_values[0]	

	#------------------------------------------------------
	# Vertex values
	#------------------------------------------------------

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


	#------------------------------------------------------
	# Boundary values
	#------------------------------------------------------

	boundary_values = stage.boundary_values
	D.stage_boundary_values = &boundary_values[0]

	boundary_values = xmomentum.boundary_values
	D.xmom_boundary_values = &boundary_values[0]

	boundary_values = ymomentum.boundary_values
	D.ymom_boundary_values = &boundary_values[0]

	boundary_values = elevation.boundary_values
	D.bed_boundary_values = &boundary_values[0]

	boundary_values = height.boundary_values
	D.height_boundary_values = &boundary_values[0]

	boundary_values = xvelocity.boundary_values
	D.xvelocity_boundary_values = &boundary_values[0]

	boundary_values = yvelocity.boundary_values
	D.yvelocity_boundary_values = &boundary_values[0]

	#------------------------------------------------------
	# Explicit and semi-implicit update values
	#------------------------------------------------------

	explicit_update = stage.explicit_update
	D.stage_explicit_update = &explicit_update[0]

	explicit_update = xmomentum.explicit_update
	D.xmom_explicit_update = &explicit_update[0]

	explicit_update = ymomentum.explicit_update
	D.ymom_explicit_update = &explicit_update[0]

	semi_implicit_update = stage.semi_implicit_update
	D.stage_semi_implicit_update = &semi_implicit_update[0]

	semi_implicit_update = xmomentum.semi_implicit_update
	D.xmom_semi_implicit_update = &semi_implicit_update[0]

	semi_implicit_update = ymomentum.semi_implicit_update
	D.ymom_semi_implicit_update = &semi_implicit_update[0]

	#------------------------------------------------------
	# Riverwall structures
	# 
	# Deal with the case when create_riverwall is called
	# but no reiverwall edges are found.
	#------------------------------------------------------
	riverwallData = domain_py_object.riverwallData

	try:
		riverwall_elevation = riverwallData.riverwall_elevation
		D.riverwall_elevation = &riverwall_elevation[0]
	except:
		D.riverwall_elevation = NULL

	try:
		riverwall_rowIndex = riverwallData.hydraulic_properties_rowIndex
		D.riverwall_rowIndex = &riverwall_rowIndex[0]
	except:
		D.riverwall_rowIndex = NULL

	try:
		riverwall_hydraulic_properties = riverwallData.hydraulic_properties
		D.riverwall_hydraulic_properties = &riverwall_hydraulic_properties[0,0]
	except:
		D.riverwall_hydraulic_properties = NULL
		

	D.ncol_riverwall_hydraulic_properties = riverwallData.ncol_hydraulic_properties
	try:
		D.nrow_riverwall_hydraulic_properties = len(riverwallData.names)
	except:
		D.nrow_riverwall_hydraulic_properties = 0

#===============================================================================
# DomainStruct: Cython wrapper for domain struct
#===============================================================================

from cpython.ref cimport PyObject
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free

from libc.stdint cimport uintptr_t

cdef class Domain_C_struct:
	"""
	Cython wrapper stored on the Python Domain object.
	Owns a domain* that is reused between calls.
	"""
	cdef domain *domain_c_struct_ptr
	cdef domain  domain_snapshot          # by-value snapshot
	cdef public object domain_py_object   # keep a reference to the Python Domain
	cdef public object verbose

	def __cinit__(self, object domain_py_object, verbose=False):
		self.domain_c_struct_ptr = <domain*> PyMem_Malloc(sizeof(domain))
		if self.domain_c_struct_ptr == NULL:
			raise MemoryError()

		self.domain_py_object = domain_py_object
		self.verbose = verbose

		# Initial fill from Python Domain
		get_python_domain_parameters(self.domain_c_struct_ptr, self.domain_py_object)
		get_python_domain_pointers(self.domain_c_struct_ptr, self.domain_py_object)

		# Substep bookkeeping for _openmp_compute_fluxes_central. Set once
		# here, not in get_python_domain_parameters, because a refresh of the
		# struct from Python must not reset the count mid-timestep. The
		# values match the module-level statics this replaces.
		self.domain_c_struct_ptr.flux_call_count = 0
		self.domain_c_struct_ptr.flux_base_call = 1
		self.domain_c_struct_ptr.flux_timestep_fluxcalls_seen = 1

		# Initial snapshot
		self.domain_snapshot = self.domain_c_struct_ptr[0]

		if self.verbose:
			print("Created Domain_C_struct at %#x"
			  % <uintptr_t> self.domain_c_struct_ptr)

	def __dealloc__(self):
		if self.verbose:
			print("Deleting Domain_C_struct at %#x"
			   % <uintptr_t> self.domain_c_struct_ptr)

		if self.domain_c_struct_ptr != NULL:
			PyMem_Free(self.domain_c_struct_ptr)
			self.domain_c_struct_ptr = NULL

	def update_domain_c_struct(self):
		get_python_domain_parameters(self.domain_c_struct_ptr, self.domain_py_object)
		get_python_domain_pointers(self.domain_c_struct_ptr, self.domain_py_object)

	cdef domain* get_domain_c_struct_ptr(self, update_c_domain_struct=False):
		if update_c_domain_struct:
			get_python_domain_parameters(self.domain_c_struct_ptr, self.domain_py_object)
			get_python_domain_pointers(self.domain_c_struct_ptr, self.domain_py_object)

		if self.verbose:
			print("Returning domain_c_struct_ptr at %#x"
			  % <uintptr_t> self.domain_c_struct_ptr)

		return self.domain_c_struct_ptr

	def refresh_snapshot(self):
		"""
		Replace the snapshot with the current struct contents.
		"""
		self.domain_snapshot = self.domain_c_struct_ptr[0]

	def diff_domain_struct(self):
		"""
		Compare current domain struct with the snapshot and print changes.
		"""
		cdef domain* D = self.domain_c_struct_ptr
		cdef domain* S = &self.domain_snapshot

		# --- scalar fields (unchanged from previous example; include as needed) ---
		if D.number_of_elements != S.number_of_elements:
			print("number_of_elements: %d -> %d"
				  % (S.number_of_elements, D.number_of_elements))
		if D.boundary_length != S.boundary_length:
			print("boundary_length: %d -> %d"
				  % (S.boundary_length, D.boundary_length))
		if D.number_of_riverwall_edges != S.number_of_riverwall_edges:
			print("number_of_riverwall_edges: %d -> %d"
				  % (S.number_of_riverwall_edges, D.number_of_riverwall_edges))
		if D.epsilon != S.epsilon:
			print("epsilon: %g -> %g" % (S.epsilon, D.epsilon))
		if D.H0 != S.H0:
			print("H0: %g -> %g" % (S.H0, D.H0))
		if D.g != S.g:
			print("g: %g -> %g" % (S.g, D.g))
		if D.optimise_dry_cells != S.optimise_dry_cells:
			print("optimise_dry_cells: %d -> %d"
				  % (S.optimise_dry_cells, D.optimise_dry_cells))
		if D.evolve_max_timestep != S.evolve_max_timestep:
			print("evolve_max_timestep: %g -> %g"
				  % (S.evolve_max_timestep, D.evolve_max_timestep))
		if D.evolve_min_timestep != S.evolve_min_timestep:
			print("evolve_min_timestep: %g -> %g"
				  % (S.evolve_min_timestep, D.evolve_min_timestep))
		if D.minimum_allowed_height != S.minimum_allowed_height:
			print("minimum_allowed_height: %g -> %g"
				  % (S.minimum_allowed_height, D.minimum_allowed_height))
		if D.maximum_allowed_speed != S.maximum_allowed_speed:
			print("maximum_allowed_speed: %g -> %g"
				  % (S.maximum_allowed_speed, D.maximum_allowed_speed))
		if D.low_froude != S.low_froude:
			print("low_froude: %d -> %d" % (S.low_froude, D.low_froude))
		if D.timestep_fluxcalls != S.timestep_fluxcalls:
			print("timestep_fluxcalls: %d -> %d"
				  % (S.timestep_fluxcalls, D.timestep_fluxcalls))
		if D.beta_w != S.beta_w:
			print("beta_w: %g -> %g" % (S.beta_w, D.beta_w))
		if D.beta_w_dry != S.beta_w_dry:
			print("beta_w_dry: %g -> %g" % (S.beta_w_dry, D.beta_w_dry))
		if D.beta_uh != S.beta_uh:
			print("beta_uh: %g -> %g" % (S.beta_uh, D.beta_uh))
		if D.beta_uh_dry != S.beta_uh_dry:
			print("beta_uh_dry: %g -> %g" % (S.beta_uh_dry, D.beta_uh_dry))
		if D.beta_vh != S.beta_vh:
			print("beta_vh: %g -> %g" % (S.beta_vh, D.beta_vh))
		if D.beta_vh_dry != S.beta_vh_dry:
			print("beta_vh_dry: %g -> %g" % (S.beta_vh_dry, D.beta_vh_dry))
		if D.ncol_riverwall_hydraulic_properties != S.ncol_riverwall_hydraulic_properties:
			print("ncol_riverwall_hydraulic_properties: %d -> %d"
				  % (S.ncol_riverwall_hydraulic_properties,
					 D.ncol_riverwall_hydraulic_properties))

		# --- pointer fields (all pointers in struct domain) ---
		if <uintptr_t>D.neighbours != <uintptr_t>S.neighbours:
			print("neighbours: %#x -> %#x"
				  % (<uintptr_t>S.neighbours, <uintptr_t>D.neighbours))
		if <uintptr_t>D.neighbour_edges != <uintptr_t>S.neighbour_edges:
			print("neighbour_edges: %#x -> %#x"
				  % (<uintptr_t>S.neighbour_edges, <uintptr_t>D.neighbour_edges))
		if <uintptr_t>D.surrogate_neighbours != <uintptr_t>S.surrogate_neighbours:
			print("surrogate_neighbours: %#x -> %#x"
				  % (<uintptr_t>S.surrogate_neighbours,
					 <uintptr_t>D.surrogate_neighbours))
		if <uintptr_t>D.normals != <uintptr_t>S.normals:
			print("normals: %#x -> %#x"
				  % (<uintptr_t>S.normals, <uintptr_t>D.normals))
		if <uintptr_t>D.edgelengths != <uintptr_t>S.edgelengths:
			print("edgelengths: %#x -> %#x"
				  % (<uintptr_t>S.edgelengths, <uintptr_t>D.edgelengths))
		if <uintptr_t>D.radii != <uintptr_t>S.radii:
			print("radii: %#x -> %#x"
				  % (<uintptr_t>S.radii, <uintptr_t>D.radii))
		if <uintptr_t>D.areas != <uintptr_t>S.areas:
			print("areas: %#x -> %#x"
				  % (<uintptr_t>S.areas, <uintptr_t>D.areas))
		if <uintptr_t>D.edge_flux_type != <uintptr_t>S.edge_flux_type:
			print("edge_flux_type: %#x -> %#x"
				  % (<uintptr_t>S.edge_flux_type, <uintptr_t>D.edge_flux_type))
		if <uintptr_t>D.tri_full_flag != <uintptr_t>S.tri_full_flag:
			print("tri_full_flag: %#x -> %#x"
				  % (<uintptr_t>S.tri_full_flag, <uintptr_t>D.tri_full_flag))
		if <uintptr_t>D.already_computed_flux != <uintptr_t>S.already_computed_flux:
			print("already_computed_flux: %#x -> %#x"
				  % (<uintptr_t>S.already_computed_flux,
					 <uintptr_t>D.already_computed_flux))
		if <uintptr_t>D.max_speed != <uintptr_t>S.max_speed:
			print("max_speed: %#x -> %#x"
				  % (<uintptr_t>S.max_speed, <uintptr_t>D.max_speed))
		if <uintptr_t>D.vertex_coordinates != <uintptr_t>S.vertex_coordinates:
			print("vertex_coordinates: %#x -> %#x"
				  % (<uintptr_t>S.vertex_coordinates,
					 <uintptr_t>D.vertex_coordinates))
		if <uintptr_t>D.edge_coordinates != <uintptr_t>S.edge_coordinates:
			print("edge_coordinates: %#x -> %#x"
				  % (<uintptr_t>S.edge_coordinates,
					 <uintptr_t>D.edge_coordinates))
		if <uintptr_t>D.centroid_coordinates != <uintptr_t>S.centroid_coordinates:
			print("centroid_coordinates: %#x -> %#x"
				  % (<uintptr_t>S.centroid_coordinates,
					 <uintptr_t>D.centroid_coordinates))
		if <uintptr_t>D.number_of_boundaries != <uintptr_t>S.number_of_boundaries:
			print("number_of_boundaries: %#x -> %#x"
				  % (<uintptr_t>S.number_of_boundaries,
					 <uintptr_t>D.number_of_boundaries))
		if <uintptr_t>D.stage_edge_values != <uintptr_t>S.stage_edge_values:
			print("stage_edge_values: %#x -> %#x"
				  % (<uintptr_t>S.stage_edge_values,
					 <uintptr_t>D.stage_edge_values))
		if <uintptr_t>D.xmom_edge_values != <uintptr_t>S.xmom_edge_values:
			print("xmom_edge_values: %#x -> %#x"
				  % (<uintptr_t>S.xmom_edge_values,
					 <uintptr_t>D.xmom_edge_values))
		if <uintptr_t>D.ymom_edge_values != <uintptr_t>S.ymom_edge_values:
			print("ymom_edge_values: %#x -> %#x"
				  % (<uintptr_t>S.ymom_edge_values,
					 <uintptr_t>D.ymom_edge_values))
		if <uintptr_t>D.bed_edge_values != <uintptr_t>S.bed_edge_values:
			print("bed_edge_values: %#x -> %#x"
				  % (<uintptr_t>S.bed_edge_values,
					 <uintptr_t>D.bed_edge_values))
		if <uintptr_t>D.height_edge_values != <uintptr_t>S.height_edge_values:
			print("height_edge_values: %#x -> %#x"
				  % (<uintptr_t>S.height_edge_values,
					 <uintptr_t>D.height_edge_values))
		if <uintptr_t>D.xvelocity_edge_values != <uintptr_t>S.xvelocity_edge_values:
			print("xvelocity_edge_values: %#x -> %#x"
				  % (<uintptr_t>S.xvelocity_edge_values,
					 <uintptr_t>D.xvelocity_edge_values))
		if <uintptr_t>D.yvelocity_edge_values != <uintptr_t>S.yvelocity_edge_values:
			print("yvelocity_edge_values: %#x -> %#x"
				  % (<uintptr_t>S.yvelocity_edge_values,
					 <uintptr_t>D.yvelocity_edge_values))
		if <uintptr_t>D.stage_centroid_values != <uintptr_t>S.stage_centroid_values:
			print("stage_centroid_values: %#x -> %#x"
				  % (<uintptr_t>S.stage_centroid_values,
					 <uintptr_t>D.stage_centroid_values))
		if <uintptr_t>D.xmom_centroid_values != <uintptr_t>S.xmom_centroid_values:
			print("xmom_centroid_values: %#x -> %#x"
				  % (<uintptr_t>S.xmom_centroid_values,
					 <uintptr_t>D.xmom_centroid_values))
		if <uintptr_t>D.ymom_centroid_values != <uintptr_t>S.ymom_centroid_values:
			print("ymom_centroid_values: %#x -> %#x"
				  % (<uintptr_t>S.ymom_centroid_values,
					 <uintptr_t>D.ymom_centroid_values))
		if <uintptr_t>D.bed_centroid_values != <uintptr_t>S.bed_centroid_values:
			print("bed_centroid_values: %#x -> %#x"
				  % (<uintptr_t>S.bed_centroid_values,
					 <uintptr_t>D.bed_centroid_values))
		if <uintptr_t>D.height_centroid_values != <uintptr_t>S.height_centroid_values:
			print("height_centroid_values: %#x -> %#x"
				  % (<uintptr_t>S.height_centroid_values,
					 <uintptr_t>D.height_centroid_values))
		if <uintptr_t>D.stage_vertex_values != <uintptr_t>S.stage_vertex_values:
			print("stage_vertex_values: %#x -> %#x"
				  % (<uintptr_t>S.stage_vertex_values,
					 <uintptr_t>D.stage_vertex_values))
		if <uintptr_t>D.xmom_vertex_values != <uintptr_t>S.xmom_vertex_values:
			print("xmom_vertex_values: %#x -> %#x"
				  % (<uintptr_t>S.xmom_vertex_values,
					 <uintptr_t>D.xmom_vertex_values))
		if <uintptr_t>D.ymom_vertex_values != <uintptr_t>S.ymom_vertex_values:
			print("ymom_vertex_values: %#x -> %#x"
				  % (<uintptr_t>S.ymom_vertex_values,
					 <uintptr_t>D.ymom_vertex_values))
		if <uintptr_t>D.bed_vertex_values != <uintptr_t>S.bed_vertex_values:
			print("bed_vertex_values: %#x -> %#x"
				  % (<uintptr_t>S.bed_vertex_values,
					 <uintptr_t>D.bed_vertex_values))
		if <uintptr_t>D.height_vertex_values != <uintptr_t>S.height_vertex_values:
			print("height_vertex_values: %#x -> %#x"
				  % (<uintptr_t>S.height_vertex_values,
					 <uintptr_t>D.height_vertex_values))
		if <uintptr_t>D.stage_boundary_values != <uintptr_t>S.stage_boundary_values:
			print("stage_boundary_values: %#x -> %#x"
				  % (<uintptr_t>S.stage_boundary_values,
					 <uintptr_t>D.stage_boundary_values))
		if <uintptr_t>D.xmom_boundary_values != <uintptr_t>S.xmom_boundary_values:
			print("xmom_boundary_values: %#x -> %#x"
				  % (<uintptr_t>S.xmom_boundary_values,
					 <uintptr_t>D.xmom_boundary_values))
		if <uintptr_t>D.ymom_boundary_values != <uintptr_t>S.ymom_boundary_values:
			print("ymom_boundary_values: %#x -> %#x"
				  % (<uintptr_t>S.ymom_boundary_values,
					 <uintptr_t>D.ymom_boundary_values))
		if <uintptr_t>D.bed_boundary_values != <uintptr_t>S.bed_boundary_values:
			print("bed_boundary_values: %#x -> %#x"
				  % (<uintptr_t>S.bed_boundary_values,
					 <uintptr_t>D.bed_boundary_values))
		if <uintptr_t>D.height_boundary_values != <uintptr_t>S.height_boundary_values:
			print("height_boundary_values: %#x -> %#x"
				  % (<uintptr_t>S.height_boundary_values,
					 <uintptr_t>D.height_boundary_values))
		if <uintptr_t>D.xvelocity_boundary_values != <uintptr_t>S.xvelocity_boundary_values:
			print("xvelocity_boundary_values: %#x -> %#x"
				  % (<uintptr_t>S.xvelocity_boundary_values,
					 <uintptr_t>D.xvelocity_boundary_values))
		if <uintptr_t>D.yvelocity_boundary_values != <uintptr_t>S.yvelocity_boundary_values:
			print("yvelocity_boundary_values: %#x -> %#x"
				  % (<uintptr_t>S.yvelocity_boundary_values,
					 <uintptr_t>D.yvelocity_boundary_values))
		if <uintptr_t>D.stage_explicit_update != <uintptr_t>S.stage_explicit_update:
			print("stage_explicit_update: %#x -> %#x"
				  % (<uintptr_t>S.stage_explicit_update,
					 <uintptr_t>D.stage_explicit_update))
		if <uintptr_t>D.xmom_explicit_update != <uintptr_t>S.xmom_explicit_update:
			print("xmom_explicit_update: %#x -> %#x"
				  % (<uintptr_t>S.xmom_explicit_update,
					 <uintptr_t>D.xmom_explicit_update))
		if <uintptr_t>D.ymom_explicit_update != <uintptr_t>S.ymom_explicit_update:
			print("ymom_explicit_update: %#x -> %#x"
				  % (<uintptr_t>S.ymom_explicit_update,
					 <uintptr_t>D.ymom_explicit_update))
		if <uintptr_t>D.edge_flux_work != <uintptr_t>S.edge_flux_work:
			print("edge_flux_work: %#x -> %#x"
				  % (<uintptr_t>S.edge_flux_work,
					 <uintptr_t>D.edge_flux_work))
		if <uintptr_t>D.neigh_work != <uintptr_t>S.neigh_work:
			print("neigh_work: %#x -> %#x"
				  % (<uintptr_t>S.neigh_work,
					 <uintptr_t>D.neigh_work))
		if <uintptr_t>D.pressuregrad_work != <uintptr_t>S.pressuregrad_work:
			print("pressuregrad_work: %#x -> %#x"
				  % (<uintptr_t>S.pressuregrad_work,
					 <uintptr_t>D.pressuregrad_work))
		if <uintptr_t>D.x_centroid_work != <uintptr_t>S.x_centroid_work:
			print("x_centroid_work: %#x -> %#x"
				  % (<uintptr_t>S.x_centroid_work,
					 <uintptr_t>D.x_centroid_work))
		if <uintptr_t>D.y_centroid_work != <uintptr_t>S.y_centroid_work:
			print("y_centroid_work: %#x -> %#x"
				  % (<uintptr_t>S.y_centroid_work,
					 <uintptr_t>D.y_centroid_work))
		if <uintptr_t>D.boundary_flux_sum != <uintptr_t>S.boundary_flux_sum:
			print("boundary_flux_sum: %#x -> %#x"
				  % (<uintptr_t>S.boundary_flux_sum,
					 <uintptr_t>D.boundary_flux_sum))
		if <uintptr_t>D.edge_river_wall_counter != <uintptr_t>S.edge_river_wall_counter:
			print("edge_river_wall_counter: %#x -> %#x"
				  % (<uintptr_t>S.edge_river_wall_counter,
					 <uintptr_t>D.edge_river_wall_counter))
		if <uintptr_t>D.riverwall_elevation != <uintptr_t>S.riverwall_elevation:
			print("riverwall_elevation: %#x -> %#x"
				  % (<uintptr_t>S.riverwall_elevation,
					 <uintptr_t>D.riverwall_elevation))
		if <uintptr_t>D.riverwall_rowIndex != <uintptr_t>S.riverwall_rowIndex:
			print("riverwall_rowIndex: %#x -> %#x"
				  % (<uintptr_t>S.riverwall_rowIndex,
					 <uintptr_t>D.riverwall_rowIndex))
		if <uintptr_t>D.riverwall_hydraulic_properties != <uintptr_t>S.riverwall_hydraulic_properties:
			print("riverwall_hydraulic_properties: %#x -> %#x"
				  % (<uintptr_t>S.riverwall_hydraulic_properties,
					 <uintptr_t>D.riverwall_hydraulic_properties))
		if <uintptr_t>D.stage_semi_implicit_update != <uintptr_t>S.stage_semi_implicit_update:
			print("stage_semi_implicit_update: %#x -> %#x"
				  % (<uintptr_t>S.stage_semi_implicit_update,
					 <uintptr_t>D.stage_semi_implicit_update))
		if <uintptr_t>D.xmom_semi_implicit_update != <uintptr_t>S.xmom_semi_implicit_update:
			print("xmom_semi_implicit_update: %#x -> %#x"
				  % (<uintptr_t>S.xmom_semi_implicit_update,
					 <uintptr_t>D.xmom_semi_implicit_update))
		if <uintptr_t>D.ymom_semi_implicit_update != <uintptr_t>S.ymom_semi_implicit_update:
			print("ymom_semi_implicit_update: %#x -> %#x"
				  % (<uintptr_t>S.ymom_semi_implicit_update,
					 <uintptr_t>D.ymom_semi_implicit_update))
		if <uintptr_t>D.friction_centroid_values != <uintptr_t>S.friction_centroid_values:
			print("friction_centroid_values: %#x -> %#x"
				  % (<uintptr_t>S.friction_centroid_values,
					 <uintptr_t>D.friction_centroid_values))
		if <uintptr_t>D.stage_backup_values != <uintptr_t>S.stage_backup_values:
			print("stage_backup_values: %#x -> %#x"
				  % (<uintptr_t>S.stage_backup_values,
					 <uintptr_t>D.stage_backup_values))
		if <uintptr_t>D.xmom_backup_values != <uintptr_t>S.xmom_backup_values:
			print("xmom_backup_values: %#x -> %#x"
				  % (<uintptr_t>S.xmom_backup_values,
					 <uintptr_t>D.xmom_backup_values))
		if <uintptr_t>D.ymom_backup_values != <uintptr_t>S.ymom_backup_values:
			print("ymom_backup_values: %#x -> %#x"
				  % (<uintptr_t>S.ymom_backup_values,
					 <uintptr_t>D.ymom_backup_values))


	def print_domain_struct(self):
		cdef domain* D = self.domain_c_struct_ptr

		print("D.number_of_elements     %d" % D.number_of_elements)
		print("D.boundary_length        %d" % D.boundary_length)
		print("D.number_of_riverwall_edges %d " % D.number_of_riverwall_edges)
		print("D.epsilon                %g" % D.epsilon)
		print("D.H0                     %g" % D.H0)
		print("D.g                      %g" % D.g)
		print("D.optimise_dry_cells     %d" % D.optimise_dry_cells)
		print("D.evolve_max_timestep    %g" % D.evolve_max_timestep)
		print("D.evolve_min_timestep    %g" % D.evolve_min_timestep)
		print("D.minimum_allowed_height %g" % D.minimum_allowed_height)
		print("D.maximum_allowed_speed  %g" % D.maximum_allowed_speed)
		print("D.low_froude             %d" % D.low_froude)
		print("D.extrapolate_velocity_second_order %d" % D.extrapolate_velocity_second_order)
		print("D.beta_w                 %g" % D.beta_w)
		print("D.beta_w_dry             %g" % D.beta_w_dry)
		print("D.beta_uh                %g" % D.beta_uh)
		print("D.beta_uh_dry            %g" % D.beta_uh_dry)
		print("D.beta_vh                %g" % D.beta_vh)
		print("D.beta_vh_dry            %g \n" % D.beta_vh_dry)



		print("D.neighbours             %#x" % <uintptr_t> D.neighbours)
		print("D.surrogate_neighbours   %#x" % <uintptr_t> D.surrogate_neighbours)
		print("D.neighbour_edges        %#x" % <uintptr_t> D.neighbour_edges)
		print("D.normals                %#x" % <uintptr_t> D.normals)
		print("D.edgelengths            %#x" % <uintptr_t> D.edgelengths)
		print("D.radii                  %#x" % <uintptr_t> D.radii)
		print("D.areas                  %#x" % <uintptr_t> D.areas)
		print("D.tri_full_flag          %#x" % <uintptr_t> D.tri_full_flag)
		print("D.already_computed_flux  %#x" % <uintptr_t> D.already_computed_flux)
		print("D.vertex_coordinates     %#x" % <uintptr_t> D.vertex_coordinates)
		print("D.edge_coordinates       %#x" % <uintptr_t> D.edge_coordinates)
		print("D.centroid_coordinates   %#x" % <uintptr_t> D.centroid_coordinates)
		print("D.max_speed              %#x" % <uintptr_t> D.max_speed)
		print("D.number_of_boundaries   %#x" % <uintptr_t> D.number_of_boundaries)
		print("D.stage_edge_values      %#x" % <uintptr_t> D.stage_edge_values)
		print("D.xmom_edge_values       %#x" % <uintptr_t> D.xmom_edge_values)
		print("D.ymom_edge_values       %#x" % <uintptr_t> D.ymom_edge_values)
		print("D.bed_edge_values        %#x" % <uintptr_t> D.bed_edge_values)
		print("D.stage_centroid_values  %#x" % <uintptr_t> D.stage_centroid_values)
		print("D.xmom_centroid_values   %#x" % <uintptr_t> D.xmom_centroid_values)
		print("D.ymom_centroid_values   %#x" % <uintptr_t> D.ymom_centroid_values)
		print("D.bed_centroid_values    %#x" % <uintptr_t> D.bed_centroid_values)
		print("D.stage_vertex_values    %#x" % <uintptr_t> D.stage_vertex_values)
		print("D.xmom_vertex_values     %#x" % <uintptr_t> D.xmom_vertex_values)
		print("D.ymom_vertex_values     %#x" % <uintptr_t> D.ymom_vertex_values)
		print("D.bed_vertex_values      %#x" % <uintptr_t> D.bed_vertex_values)
		print("D.height_vertex_values   %#x" % <uintptr_t> D.height_vertex_values)
		print("D.stage_boundary_values  %#x" % <uintptr_t> D.stage_boundary_values)
		print("D.xmom_boundary_values   %#x" % <uintptr_t> D.xmom_boundary_values)
		print("D.ymom_boundary_values   %#x" % <uintptr_t> D.ymom_boundary_values)
		print("D.bed_boundary_values    %#x" % <uintptr_t> D.bed_boundary_values)
		print("D.stage_explicit_update  %#x" % <uintptr_t> D.stage_explicit_update)
		print("D.xmom_explicit_update   %#x" % <uintptr_t> D.xmom_explicit_update)
		print("D.ymom_explicit_update   %#x" % <uintptr_t> D.ymom_explicit_update)
		print("D.edge_river_wall_counter    %#x" % <uintptr_t> D.edge_river_wall_counter)
		print("D.stage_semi_implicit_update %#x" % <uintptr_t> D.stage_semi_implicit_update)
		print("D.xmom_semi_implicit_update  %#x" % <uintptr_t> D.xmom_semi_implicit_update)
		print("D.ymom_semi_implicit_update  %#x" % <uintptr_t> D.ymom_semi_implicit_update)
		print("D.friction_centroid_values   %#x" % <uintptr_t> D.friction_centroid_values)
		print("D.stage_backup_values        %#x" % <uintptr_t> D.stage_backup_values)
		print("D.xmom_backup_values         %#x" % <uintptr_t> D.xmom_backup_values)
		print("D.ymom_backup_values         %#x" % <uintptr_t> D.ymom_backup_values)


#===============================================================================

def set_omp_num_threads(anuga_int num_threads):
	"""
	Set the number of OpenMP threads to use.
	"""
	_openmp_set_omp_num_threads(num_threads)

	
def setup_Domain_C_struct(object domain_py_object, update_domain_c_struct=False):
	"""
	Setup the Domain_C_struct on the Python Domain object.
	"""
	# Ensure work arrays are allocated before the C struct reads their pointers.
	# They are kept None at domain creation to avoid cold virtual pages during
	# distribute(); _ensure_work_arrays() allocates them on the first evolve call.
	if hasattr(domain_py_object, '_ensure_work_arrays'):
		domain_py_object._ensure_work_arrays()
	domain_py_object._Domain_C_struct = Domain_C_struct(domain_py_object)

def update_Domain_C_struct(object domain_py_object):
	"""
	Update the Domain_C_struct on the Python Domain object.
	"""
	if domain_py_object._Domain_C_struct is None:
		setup_Domain_C_struct(domain_py_object)

	domain_py_object._Domain_C_struct.update_domain_c_struct()

cdef domain* get_domain_c_struct_ptr(object domain_py_object, update_domain_c_struct=False):
	"""
	Get the domain* from the Domain_C_struct on the Python Domain object.
	"""
	if domain_py_object._Domain_C_struct is None:
		setup_Domain_C_struct(domain_py_object)

	cdef Domain_C_struct domain_c_struct
	domain_c_struct = domain_py_object._Domain_C_struct
	cdef domain* D = domain_c_struct.get_domain_c_struct_ptr(update_domain_c_struct)

	return D

def compute_fluxes_ext_central(object domain_py_object, 
                               double timestep, 
							   update_domain_c_struct=False):

	cdef domain* D = get_domain_c_struct_ptr(domain_py_object, update_domain_c_struct=update_domain_c_struct)

	with nogil:
		timestep =  _openmp_compute_fluxes_central(D, timestep)

	return timestep

def extrapolate_second_order_sw(object domain_py_object, update_domain_c_struct=False):

	cdef domain* D = get_domain_c_struct_ptr(domain_py_object, update_domain_c_struct=update_domain_c_struct)

	with nogil:
		_openmp_extrapolate_second_order_edge_sw(D)

def distribute_to_edges(object domain_py_object, update_domain_c_struct=False):

	cdef domain* D = get_domain_c_struct_ptr(domain_py_object, update_domain_c_struct=update_domain_c_struct)

	with nogil:
		_openmp_extrapolate_second_order_edge_sw(D)

def distribute_to_edges_and_vertices(object domain_py_object, 
                                    distribute_to_vertices=True, 
									update_domain_c_struct=False):

	cdef domain* D = get_domain_c_struct_ptr(domain_py_object, update_domain_c_struct=update_domain_c_struct)

	with nogil:
		_openmp_extrapolate_second_order_edge_sw(D)

	if distribute_to_vertices:
		with nogil:
			_openmp_distribute_edges_to_vertices(D)


def distribute_edges_to_vertices(object domain_py_object, update_domain_c_struct=False):

	cdef domain* D = get_domain_c_struct_ptr(domain_py_object, update_domain_c_struct=update_domain_c_struct)

	with nogil:
		_openmp_distribute_edges_to_vertices(D)

	

def extrapolate_second_order_edge_sw(object domain_py_object, 
                                     distribute_to_vertices=True,
                                     update_domain_c_struct=False):

	cdef domain* D = get_domain_c_struct_ptr(domain_py_object, update_domain_c_struct=update_domain_c_struct)


	with nogil:
		_openmp_extrapolate_second_order_edge_sw(D)
	if distribute_to_vertices:
		with nogil:
			_openmp_distribute_edges_to_vertices(D)


cdef extern from "gpu/core_kernels.h" nogil:
	void core_apply_sediment_source(domain* D, double timestep)
	void core_apply_bedload(domain* D, double timestep)
	anuga_int core_apply_repose(domain* D)


def apply_sediment_source(object domain_py_object, double timestep,
                          update_domain_c_struct=False):
	"""Fractional step: apply E-D to m and the resulting bed change to z."""

	cdef domain* D = get_domain_c_struct_ptr(domain_py_object, update_domain_c_struct=update_domain_c_struct)

	with nogil:
		core_apply_sediment_source(D, timestep)


def apply_bedload(object domain_py_object, double timestep,
                  update_domain_c_struct=False):
	"""Fractional step: bedload divergence [G-5] applied to the bed."""

	cdef domain* D = get_domain_c_struct_ptr(domain_py_object, update_domain_c_struct=update_domain_c_struct)

	with nogil:
		core_apply_bedload(D, timestep)


def apply_repose(object domain_py_object, update_domain_c_struct=False):
	"""Fractional step: angle-of-repose relaxation (spec 7).

	Returns the number of sweeps used. Equal to the configured cap means the
	cap was hit and the bed may still be over-steep.
	"""

	cdef domain* D = get_domain_c_struct_ptr(domain_py_object, update_domain_c_struct=update_domain_c_struct)
	cdef anuga_int sweeps

	with nogil:
		sweeps = core_apply_repose(D)

	return sweeps


def protect_new(object domain_py_object, update_domain_c_struct=False):

	cdef domain* D = get_domain_c_struct_ptr(domain_py_object, update_domain_c_struct=update_domain_c_struct)

	cdef double mass_error

	with nogil:
		mass_error = _openmp_protect(D)


	return mass_error

def manning_friction_flat_semi_implicit(object domain_py_object, update_domain_c_struct=False):
	
	cdef domain* D = get_domain_c_struct_ptr(domain_py_object, update_domain_c_struct=update_domain_c_struct)

	with nogil:
		_openmp_manning_friction_flat_semi_implicit(D)



def vegetation_friction_semi_implicit(object domain_py_object, update_domain_c_struct=False):
	"""Baptist vegetation drag (spec 8) into the semi-implicit update; no-op when off."""
	cdef domain* D = get_domain_c_struct_ptr(domain_py_object, update_domain_c_struct=update_domain_c_struct)

	with nogil:
		_openmp_vegetation_friction_semi_implicit(D)

def manning_friction_sloped_semi_implicit(object domain_py_object, update_domain_c_struct=False):
	
	cdef domain* D = get_domain_c_struct_ptr(domain_py_object, update_domain_c_struct=update_domain_c_struct)

	with nogil:
		_openmp_manning_friction_sloped_semi_implicit(D)

def manning_friction_sloped_semi_implicit_edge_based(object domain_py_object, update_domain_c_struct=False):
	
	cdef domain* D = get_domain_c_struct_ptr(domain_py_object, update_domain_c_struct=update_domain_c_struct)

	with nogil:
		_openmp_manning_friction_sloped_semi_implicit_edge_based(D)



# FIXME SR: Why is the order of arguments different from the C function?
def manning_friction_flat(double g, double eps,
            np.ndarray[double, ndim=1, mode="c"] w not None,
			np.ndarray[double, ndim=1, mode="c"] uh not None,
			np.ndarray[double, ndim=1, mode="c"] vh not None,
			np.ndarray[double, ndim=1, mode="c"] z_centroid not None,
			np.ndarray[double, ndim=1, mode="c"] eta not None,
			np.ndarray[double, ndim=1, mode="c"] xmom not None,
			np.ndarray[double, ndim=1, mode="c"] ymom not None):
	
	cdef anuga_int N
	
	N = w.shape[0]
	if uh.shape[0] != N:
		raise ValueError('manning_friction_flat: uh has length %d, expected %d' % (uh.shape[0], N))
	if vh.shape[0] != N:
		raise ValueError('manning_friction_flat: vh has length %d, expected %d' % (vh.shape[0], N))
	if z_centroid.shape[0] != N:
		raise ValueError('manning_friction_flat: z_centroid has length %d, expected %d' % (z_centroid.shape[0], N))
	if eta.shape[0] != N:
		raise ValueError('manning_friction_flat: eta has length %d, expected %d' % (eta.shape[0], N))
	if xmom.shape[0] != N:
		raise ValueError('manning_friction_flat: xmom has length %d, expected %d' % (xmom.shape[0], N))
	if ymom.shape[0] != N:
		raise ValueError('manning_friction_flat: ymom has length %d, expected %d' % (ymom.shape[0], N))
	_openmp_manning_friction_flat(g, eps, N, &w[0], &z_centroid[0], &uh[0], &vh[0], &eta[0], &xmom[0], &ymom[0])

# FIXME SR: Why is the order of arguments different from the C function?
def manning_friction_sloped(double g, double eps,
        np.ndarray[double, ndim=2, mode="c"] x_vertex not None,
		np.ndarray[double, ndim=1, mode="c"] w not None,
		np.ndarray[double, ndim=1, mode="c"] uh not None,
		np.ndarray[double, ndim=1, mode="c"] vh not None,
		np.ndarray[double, ndim=2, mode="c"] z_vertex not None,
		np.ndarray[double, ndim=1, mode="c"] eta not None,
		np.ndarray[double, ndim=1, mode="c"] xmom not None,
		np.ndarray[double, ndim=1, mode="c"] ymom not None):
		
	cdef anuga_int N
	
	N = w.shape[0]
	if uh.shape[0] != N:
		raise ValueError('manning_friction_sloped: uh has length %d, expected %d' % (uh.shape[0], N))
	if vh.shape[0] != N:
		raise ValueError('manning_friction_sloped: vh has length %d, expected %d' % (vh.shape[0], N))
	if eta.shape[0] != N:
		raise ValueError('manning_friction_sloped: eta has length %d, expected %d' % (eta.shape[0], N))
	if xmom.shape[0] != N:
		raise ValueError('manning_friction_sloped: xmom has length %d, expected %d' % (xmom.shape[0], N))
	if ymom.shape[0] != N:
		raise ValueError('manning_friction_sloped: ymom has length %d, expected %d' % (ymom.shape[0], N))
	if x_vertex.shape[0] != N or z_vertex.shape[0] != N:
		raise ValueError('manning_friction_sloped: x_vertex/z_vertex must have N rows')
	_openmp_manning_friction_sloped(g, eps, N, &x_vertex[0,0], &w[0], &z_vertex[0,0], &uh[0], &vh[0], &eta[0], &xmom[0], &ymom[0])

# FIXME SR: Why is the order of arguments different from the C function?
def manning_friction_sloped_edge_based(double g, double eps,
        np.ndarray[double, ndim=2, mode="c"] x_edge not None,
		np.ndarray[double, ndim=1, mode="c"] w not None,
		np.ndarray[double, ndim=1, mode="c"] uh not None,
		np.ndarray[double, ndim=1, mode="c"] vh not None,
		np.ndarray[double, ndim=2, mode="c"] z_edge not None,
		np.ndarray[double, ndim=1, mode="c"] eta not None,
		np.ndarray[double, ndim=1, mode="c"] xmom not None,
		np.ndarray[double, ndim=1, mode="c"] ymom not None):
		
	cdef anuga_int N
	
	N = w.shape[0]
	if uh.shape[0] != N:
		raise ValueError('manning_friction_sloped_edge_based: uh has length %d, expected %d' % (uh.shape[0], N))
	if vh.shape[0] != N:
		raise ValueError('manning_friction_sloped_edge_based: vh has length %d, expected %d' % (vh.shape[0], N))
	if eta.shape[0] != N:
		raise ValueError('manning_friction_sloped_edge_based: eta has length %d, expected %d' % (eta.shape[0], N))
	if xmom.shape[0] != N:
		raise ValueError('manning_friction_sloped_edge_based: xmom has length %d, expected %d' % (xmom.shape[0], N))
	if ymom.shape[0] != N:
		raise ValueError('manning_friction_sloped_edge_based: ymom has length %d, expected %d' % (ymom.shape[0], N))
	if x_edge.shape[0] != N or z_edge.shape[0] != N:
		raise ValueError('manning_friction_sloped_edge_based: x_edge/z_edge must have N rows')
	_openmp_manning_friction_sloped_edge_based(g, eps, N, &x_edge[0,0], &w[0], &z_edge[0,0], &uh[0], &vh[0], &eta[0], &xmom[0], &ymom[0])


def fix_negative_cells(object domain_py_object, update_domain_c_struct=False):

	cdef anuga_int num_negative_cells
	cdef domain* D = get_domain_c_struct_ptr(domain_py_object, update_domain_c_struct=update_domain_c_struct)

	with nogil:
		num_negative_cells = _openmp_fix_negative_cells(D)

	return num_negative_cells

def update_conserved_quantities(object domain_py_object, 
                                double timestep,
                                update_domain_c_struct=False):

	cdef anuga_int num_negative_cells
	cdef double negative_volume
	cdef domain* D = get_domain_c_struct_ptr(domain_py_object, update_domain_c_struct=update_domain_c_struct)

	with nogil:
		_openmp_update_conserved_quantities(D, timestep)
		# Measure the deficit volume BEFORE the clamp erases it (stage->bed).
		negative_volume = _openmp_negative_cells_volume(D)
		num_negative_cells = _openmp_fix_negative_cells(D)

	return num_negative_cells, negative_volume

def saxpy_conserved_quantities(object domain_py_object, 
                               double a, double b, double c,
                               update_domain_c_struct=False):

	cdef domain* D = get_domain_c_struct_ptr(domain_py_object, update_domain_c_struct=update_domain_c_struct)

	with nogil:
		_openmp_saxpy_conserved_quantities(D, a, b, c)


def backup_conserved_quantities(object domain_py_object, update_domain_c_struct=False):

	cdef domain* D = get_domain_c_struct_ptr(domain_py_object, update_domain_c_struct=update_domain_c_struct)

	with nogil:
		_openmp_backup_conserved_quantities(D)

def ader_ck_predictor(object domain_py_object, double dt, update_domain_c_struct=False):
	"""ADER Cauchy-Kovalewski predictor: advance centroids by dt in-place.
	Call after distribute_to_vertices_and_edges(); uses edge values to recover
	cell slopes, then evaluates SWE time derivatives locally.
	"""
	cdef domain* D = get_domain_c_struct_ptr(domain_py_object, update_domain_c_struct=update_domain_c_struct)

	with nogil:
		_openmp_ader_ck_predictor(D, dt)

def ader_ck_predictor_edge(object domain_py_object, double dt, update_domain_c_struct=False):
	"""Fused ADER-2 predictor: shift edge values to Q^{n+1/2}, centroids unchanged.

	Equivalent to ader_ck_predictor() followed by distribute_to_vertices_and_edges(),
	but faster: skips the full extrapolation pass and leaves Q^n in centroids so no
	backup/saxpy restore is needed in evolve_one_ader2_step().
	"""
	cdef domain* D = get_domain_c_struct_ptr(domain_py_object, update_domain_c_struct=update_domain_c_struct)

	with nogil:
		_openmp_ader_ck_predictor_edge(D, dt)

def evaluate_reflective_segment(object domain_py_object, 
                                np.ndarray[np.int64_t, ndim=1, mode="c"] segment_edges not None, 
								np.ndarray[np.int64_t, ndim=1, mode="c"] vol_ids not None, 
								np.ndarray[np.int64_t, ndim=1, mode="c"] edge_ids not None, 
                                update_domain_c_struct=False): 

	cdef domain* D = get_domain_c_struct_ptr(domain_py_object, update_domain_c_struct=update_domain_c_struct)
	cdef anuga_int N
	N = segment_edges.shape[0]

	with nogil:
		_openmp_evaluate_reflective_segment(D, N, &segment_edges[0], &vol_ids[0], &edge_ids[0])


def rotate(np.ndarray[double, ndim=1, mode="c"] q not None, np.ndarray[double, ndim=1, mode="c"] normal not None, anuga_int direction):

	assert normal.shape[0] == 2, "Normal vector must have 2 components"
	cdef np.ndarray[double, ndim=1, mode="c"] r
	cdef double n1, n2
	n1 = normal[0]
	n2 = normal[1]
	if direction == -1:
		n2 = -n2
	r = np.ascontiguousarray(np.copy(q))
	__rotate(&r[0], n1, n2)
	return r

def flux_function_central(
	np.ndarray[double, ndim=1, mode="c"] normal not None,
	np.ndarray[double, ndim=1, mode="c"] ql not None,
	np.ndarray[double, ndim=1, mode="c"] qr not None,
	double h_left,
	double h_right,
	double hle,
	double hre,
	np.ndarray[double, ndim=1, mode="c"] edgeflux not None,
	double epsilon,
	double ze,
	double g,
	double H0,
	double hc,
	double hc_n,
	anuga_int low_froude
):
	cdef double h0, limiting_threshold, max_speed, pressure_flux
	cdef anuga_int err

	h0 = H0 * H0
	limiting_threshold = 10 * H0

	err = __openmp__flux_function_central(
		ql[0], ql[1], ql[2],
		qr[0], qr[1], qr[2],
		h_left, h_right, hle, hre,
		normal[0], normal[1], epsilon, ze, g,
		&edgeflux[0], &edgeflux[1], &edgeflux[2],
		&max_speed, &pressure_flux, low_froude
	)

	assert err >= 0, "Discontinuous Elevation"

	return max_speed, pressure_flux

def gravity(object domain_py_object, update_domain_c_struct=False):

	cdef domain* D = get_domain_c_struct_ptr(domain_py_object, update_domain_c_struct=update_domain_c_struct)

	err = _openmp_gravity(D)
	if err != 0:
		raise RuntimeError('gravity: C kernel _openmp_gravity failed with error code %d' % err)

def gravity_wb(object domain_py_object, update_domain_c_struct=False):

	cdef domain* D = get_domain_c_struct_ptr(domain_py_object, update_domain_c_struct=update_domain_c_struct)

	err = _openmp_gravity_wb(D)
	if err != 0:
		raise RuntimeError('gravity_wb: C kernel _openmp_gravity_wb failed with error code %d' % err)


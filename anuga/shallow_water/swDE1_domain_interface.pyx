# cython: boundscheck=False
# cython: wraparound=False
import cython

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

# declare the interface to the C code
cdef extern from "swDE1_domain_ext.c":
	struct domain:
		pass
	struct edge:
		pass
	domain* get_python_domain(domain* D, object domain)
	unsigned int Mod_of_power_2(unsigned int n, unsigned int d)
	int _rotate(double* q, double n1, double n2)
	int _flux_function_toro(double* q_left, double* q_right, double h_left, double h_right, double hle, double hre, double n1, double n2, double epsilon, double ze, double limiting_threshold, double g, double* edge_flux, double* max_speed, double* pressure_flux, double hc, double hc_n)
	int _flux_function_central(double* q_left, double* q_right, double h_left, double h_right, double hle, double hre, double n1, double n2, double epsilon, double ze, double limiting_threshold, double g, double* edgeflux, double* max_speed, double* pressure_flux, double hc, double hc_n, long low_froude)
	int _compute_flux_update_frequency(domain* D, double timestep)
	double adjust_edgeflux_with_weir(double* edgeflux, double h_left, double h_right, double g, double weir_height, double Qfactor, double s1, double s2, double h1, double h2, double* max_speed_local)
	double _compute_fluxes_central(domain* D, double timestep)
	double _protect(int N, double minimum_allowed_height, double maximum_allowed_speed, double epsilon, double* wc, double* wv, double* zc, double* zv, double* xmomc, double* ymomc, double* areas, double* xc, double* yc)
	double _protect_new(domain* D)
	int find_qmin_and_qmax(double dq0, double dq1, double dq2, double* qmin, double* qmax)
	int limit_gradient(double* dqv, double qmin, double qmax, double beta_w)
	int _extrapolate_second_order_edge_sw(domain* D)


def compute_fluxes_ext_central(object domain_object, double timestep):

	cdef domain D

	get_python_domain(&D, domain_object)

	timestep = _compute_fluxes_central(&D, timestep)

	return timestep

def extrapolate_second_order_edge_sw(object domain_object):

	cdef domain D

	cdef int e

	get_python_domain(&D, domain_object)

	e = _extrapolate_second_order_edge_sw(&D)

	if e == -1:
		return None

def protect_new(object domain_object):

	cdef domain D

	cdef double mass_error

	get_python_domain(&D, domain_object)

	mass_error = _protect_new(&D)

	return mass_error

def compute_flux_update_frequency(object domain_object, double timestep):

	cdef domain D

	get_python_domain(&D, domain_object)

	_compute_flux_update_frequency(&D, timestep)



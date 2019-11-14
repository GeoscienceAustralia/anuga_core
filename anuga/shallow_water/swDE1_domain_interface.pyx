# cython: boundscheck=False
# cython: wraparound=False
import cython

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

# import header
from swDE1_domain_header cimport *

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



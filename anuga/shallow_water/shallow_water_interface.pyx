# cython: boundscheck=False
# cython: wraparound=False
import cython

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

# import header
from shallow_water_header cimport *

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

	h0 = H0*H0
	limiting_threshold = 10*H0

	_flux_function_central(&ql[0], &qr[0], zl, zr, normal[0], normal[1], epsilon, h0, limiting_threshold, g, &edgeflux[0], &max_speed)

	return max_speed

def extrapolate_second_order_sw(object domain_object):

	cdef domain D
	cdef int err

	get_python_domain(&D, domain_object)

	err = _extrapolate_second_order_sw(&D)

	if err == -1:
		return None

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

	return timestep

def compute_fluxes_ext_wb_3(object domain_object):

	cdef domain D
	cdef double timestep

	get_python_domain(&D, domain_object)

	timestep = _compute_fluxes_central_wb_3(&D)

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















#cython: wraparound=False, boundscheck=False, cdivision=True, profile=False, nonecheck=False, overflowcheck=False, cdivision_warnings=False, unraisable_tracebacks=False
import cython
from libc.stdint cimport int64_t, int32_t

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

# declare the interface to the C code
cdef extern from "util_ext.h":
    int64_t _gradient(double x0, double y0, double x1, double y1, double x2, double y2, double q0, double q1, double q2, double *a, double *b)
    int64_t _gradient2(double x0, double y0, double x1, double y1, double q0, double q1, double *a, double *b)

cdef extern from "float.h":
    const double DBL_DIG

def gradient(double x0, double y0, double x1, double y1, double x2, double y2, double q0, double q1, double q2):

    cdef double a, b

    _gradient(x0, y0, x1, y1, x2, y2, q0, q1, q2, &a, &b)

    return (a,b)

def gradient2(double x0, double y0, double x1, double y1, double q0, double q1):

    cdef double a, b

    _gradient2(x0, y0, x1, y1, q0, q1, &a, &b)

    return (a,b)

def double_precision():

    return DBL_DIG

#cython: wraparound=False, boundscheck=False, cdivision=True, profile=False, nonecheck=False, overflowcheck=False, cdivision_warnings=False, unraisable_tracebacks=False
import cython
from libc.stdint cimport int64_t

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

# declare the interface to the C code
cdef extern from "polygon.c":
    int64_t __point_on_line(double x, double y, double x0, double y0, double x1, double y1, double rtol, double atol)
    int64_t __interpolate_polyline(int64_t number_of_nodes, int64_t number_of_points, double* data, double* polyline_nodes, int64_t* gauge_neighbour_id, double* interpolation_points, double* interpolated_values, double rtol, double atol)
    int64_t __polygon_overlap(double* polygon, double* triangles, int64_t* indices, int64_t M, int64_t polygon_number_of_vertices)
    int64_t __line_intersect(double* line, double* triangles, int64_t* indices, int64_t M)
    int64_t __is_inside_triangle(double* point, double* triangle, int64_t closed, double rtol, double atol)
    int64_t __separate_points_by_polygon(int64_t M, int64_t N, double* points, double* polygon, int64_t* indices, int64_t closed, int64_t verbose)

def _point_on_line(double x,\
                    double y,\
                    double x0,\
                    double y0,\
                    double x1,\
                    double y1,\
                    double rtol,\
                    double atol):
    
    cdef int64_t res
    res = __point_on_line(x, y, x0, y0, x1, y1, rtol, atol)
    return res

def _interpolate_polyline(np.ndarray[double, ndim=1, mode="c"] data not None,\
                        np.ndarray[double, ndim=2, mode="c"] polyline_nodes not None,\
                        np.ndarray[int64_t, ndim=1, mode="c"] gauge_neighbour_id not None,\
                        np.ndarray[double, ndim=2, mode="c"] interpolation_points not None,\
                        np.ndarray[double, ndim=1, mode="c"] interpolated_values not None,\
                        double rtol,\
                        double atol):

    cdef int64_t number_of_nodes, number_of_points, res

    number_of_nodes = polyline_nodes.shape[0]
    number_of_points = interpolation_points.shape[0]

    res = __interpolate_polyline(number_of_nodes,\
                                number_of_points,\
                                &data[0],\
                                &polyline_nodes[0,0],\
                                &gauge_neighbour_id[0],\
                                &interpolation_points[0,0],\
                                &interpolated_values[0],\
                                rtol,\
                                atol)

def _polygon_overlap(np.ndarray[double, ndim=2, mode="c"] polygon not None, np.ndarray[double, ndim=2, mode="c"] triangles not None, np.ndarray[int64_t, ndim=1, mode="c"] indices not None):

    cdef int64_t res

    res = __polygon_overlap(&polygon[0,0], &triangles[0,0], &indices[0], triangles.shape[0]/3, polygon.shape[0])

    return res

def _line_intersect(np.ndarray[double, ndim=2, mode="c"] line not None, np.ndarray[double, ndim=2, mode="c"] triangles not None, np.ndarray[int64_t, ndim=1, mode="c"] indices not None):

    cdef int64_t res

    res = __line_intersect(&line[0,0], &triangles[0,0], &indices[0], triangles.shape[0]/3)

    return res

def _is_inside_triangle(np.ndarray[double, ndim=1, mode="c"] point not None,\
                        np.ndarray[double, ndim=2, mode="c"] triangle not None,\
                        int64_t closed,\
                        double rtol,\
                        double atol):

    cdef int64_t res

    res = __is_inside_triangle(&point[0], &triangle[0,0], closed, rtol, atol)

    return res

def _separate_points_by_polygon(np.ndarray[double, ndim=2, mode="c"] points not None,\
                                np.ndarray[double, ndim=2, mode="c"] polygon not None,\
                                np.ndarray[int64_t, ndim=1, mode="c"] indices not None,\
                                int64_t closed,\
                                int64_t verbose):

    cdef int64_t count, M, N

    M = points.shape[0]
    N = polygon.shape[0]

    if verbose:
        print ("Got %d points and %d polygon vertices" % (M,N))

    count = __separate_points_by_polygon(M, N, &points[0,0], &polygon[0,0], &indices[0], closed, verbose)

    return count

#cython: wraparound=False, boundscheck=False, cdivision=True, profile=False, nonecheck=False, overflowcheck=False, cdivision_warnings=False, unraisable_tracebacks=False
import cython

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

# declare the interface to the C code
cdef extern from "polygon.c":
    int __point_on_line(double x, double y, double x0, double y0, double x1, double y1, double rtol, double atol)
    int __interpolate_polyline(int number_of_nodes, int number_of_points, double* data, double* polyline_nodes, long* gauge_neighbour_id, double* interpolation_points, double* interpolated_values, double rtol, double atol)
    int __polygon_overlap(double* polygon, double* triangles, long* indices, int M, int polygon_number_of_vertices)
    int __line_intersect(double* line, double* triangles, long* indices, int M)
    int __is_inside_triangle(double* point, double* triangle, int closed, double rtol, double atol)
    int __separate_points_by_polygon(int M, int N, double* points, double* polygon, long* indices, int closed, int verbose)

def _point_on_line(double x,\
                    double y,\
                    double x0,\
                    double y0,\
                    double x1,\
                    double y1,\
                    double rtol,\
                    double atol):
    
    cdef int res
    res = __point_on_line(x, y, x0, y0, x1, y1, rtol, atol)
    return res

def _interpolate_polyline(np.ndarray[double, ndim=1, mode="c"] data not None,\
                        np.ndarray[double, ndim=2, mode="c"] polyline_nodes not None,\
                        np.ndarray[long, ndim=1, mode="c"] gauge_neighbour_id not None,\
                        np.ndarray[double, ndim=2, mode="c"] interpolation_points not None,\
                        np.ndarray[double, ndim=1, mode="c"] interpolated_values not None,\
                        double rtol,\
                        double atol):

    cdef int number_of_nodes, number_of_points, res

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

def _polygon_overlap(np.ndarray[double, ndim=2, mode="c"] polygon not None, np.ndarray[double, ndim=2, mode="c"] triangles not None, np.ndarray[long, ndim=1, mode="c"] indices not None):

    cdef int res

    res = __polygon_overlap(&polygon[0,0], &triangles[0,0], &indices[0], triangles.shape[0]/3, polygon.shape[0])

    return res

def _line_intersect(np.ndarray[double, ndim=2, mode="c"] line not None, np.ndarray[double, ndim=2, mode="c"] triangles not None, np.ndarray[long, ndim=1, mode="c"] indices not None):

    cdef int res

    res = __line_intersect(&line[0,0], &triangles[0,0], &indices[0], triangles.shape[0]/3)

    return res

def _is_inside_triangle(np.ndarray[double, ndim=1, mode="c"] point not None,\
                        np.ndarray[double, ndim=2, mode="c"] triangle not None,\
                        int closed,\
                        double rtol,\
                        double atol):

    cdef int res

    res = __is_inside_triangle(&point[0], &triangle[0,0], closed, rtol, atol)

    return res

def _separate_points_by_polygon(np.ndarray[double, ndim=2, mode="c"] points not None,\
                                np.ndarray[double, ndim=2, mode="c"] polygon not None,\
                                np.ndarray[long, ndim=1, mode="c"] indices not None,\
                                int closed,\
                                int verbose):

    cdef int count, M, N

    M = points.shape[0]
    N = polygon.shape[0]

    if verbose:
        print ("Got %d points and %d polygon vertices") % (M,N)

    count = __separate_points_by_polygon(M, N, &points[0,0], &polygon[0,0], &indices[0], closed, verbose)

    return count

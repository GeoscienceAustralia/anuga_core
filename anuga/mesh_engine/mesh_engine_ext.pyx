#cython: wraparound=False, boundscheck=False, cdivision=True, profile=False, nonecheck=False, overflowcheck=False, cdivision_warnings=False, unraisable_tracebacks=False
import cython
from libc.stdlib cimport malloc, free
# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

ctypedef double REAL

# declare the interface to the C code
cdef extern from "triangle.c":
    struct triangulateio:
        REAL* pointlist
        REAL* pointattributelist
        int* pointmarkerlist
        int numberofpoints
        int numberofpointattributes
        int* trianglelist
        REAL* triangleattributelist
        REAL* trianglearealist
        int* neighborlist
        int numberoftriangles
        int numberofcorners
        int numberoftriangleattributes
        int* segmentlist
        int* segmentmarkerlist
        int numberofsegments
        REAL* holelist
        int numberofholes
        REAL* regionlist
        int numberofregions
        int* edgelist
        int* edgemarkerlist
        REAL* normlist
        int numberofedges
    void triangulate(char*, triangulateio*, triangulateio*, triangulateio*)

cdef c_int_array_to_numpy_array(int* c_array, int c_array_size, np.ndarray[int, ndim=1, mode="c"] dimensions):
    cdef int j
    np_arr = np.zeros(c_array_size, dtype=np.int32)
    for j in xrange(c_array_size):
        np_arr[j] = c_array[j]
    return np.ascontiguousarray(np_arr.reshape(dimensions))

cdef c_double_array_to_numpy_array(double* c_array, int c_array_size, np.ndarray[int, ndim=1, mode="c"] dimensions):
    cdef int j
    np_arr = np.zeros(c_array_size, dtype=np.float)
    for j in xrange(c_array_size):
        np_arr[j] = c_array[j]
    return np.ascontiguousarray(np_arr.reshape(dimensions))

def genMesh(np.ndarray pointlist,\
            np.ndarray seglist,\
            np.ndarray holelist,\
            np.ndarray regionlist,\
            np.ndarray pointattributelist,\
            np.ndarray segmarkerlist,\
            char* mod):

    cdef triangulateio in_t, out_t
    cdef triangulateio in_test

    cdef np.ndarray[int, ndim=1, mode="c"] dimensions

    cdef REAL Attr
    cdef int i, j, iatt, n, write_here, N
    cdef int a, b, c
    cdef int marker
    cdef int tuplesize
    cdef int index = 0
    cdef double x,y

    cdef np.ndarray[int, ndim=2, mode="c"] gentrianglelist
    cdef np.ndarray[double, ndim=2, mode="c"] genpointlist
    cdef np.ndarray[int, ndim=1, mode="c"] genpointmarkerlist
    cdef np.ndarray[double, ndim=2, mode="c"] genpointattributelist
    cdef np.ndarray[double, ndim=2, mode="c"] gentriangleattributelist
    cdef np.ndarray[int, ndim=2, mode="c"] gensegmentlist
    cdef np.ndarray[int, ndim=1, mode="c"] gensegmentmarkerlist
    cdef np.ndarray[int, ndim=2, mode="c"] genneighborlist

    in_t.numberofpoints = pointlist.shape[0]
    in_t.pointlist = <double* > pointlist.data
    in_t.pointmarkerlist = <int* >NULL

    in_t.numberofregions = regionlist.shape[0]
    in_t.regionlist = <double* > regionlist.data

    in_t.numberofsegments = seglist.shape[0]
    in_t.segmentlist = <int* > seglist.data
    in_t.segmentmarkerlist = <int* > segmarkerlist.data

    in_t.numberoftriangles = 0
    in_t.trianglelist = <int* >NULL
    in_t.numberoftriangleattributes = 0;
    in_t.triangleattributelist = <REAL* >NULL
    in_t.trianglearealist = <REAL* >NULL
    in_t.neighborlist = <int* >NULL
    in_t.numberofcorners = 0

    in_t.numberofholes = holelist.shape[0]

    if in_t.numberofholes != 0:
        in_t.holelist = <double* > holelist.data
    else:
        in_t.holelist = <REAL* >NULL

    if pointattributelist.shape[0] == 0:
        in_t.numberofpointattributes = 0
        in_t.pointattributelist = <double* >NULL
    else:
        if pointattributelist.shape[1] == 0:
            in_t.numberofpointattributes = 0
            in_t.pointattributelist = <double* >NULL
        else:
            in_t.numberofpointattributes = pointattributelist.shape[1]
            in_t.pointattributelist = <double* > pointattributelist.data

    out_t.pointlist = <REAL* >NULL
    out_t.pointmarkerlist = <int* >NULL
    out_t.pointattributelist = <REAL* >NULL

    out_t.trianglelist = <int* >NULL
    out_t.triangleattributelist = <REAL* >NULL
    out_t.trianglearealist = <REAL* >NULL
    out_t.neighborlist = <int* >NULL

    out_t.segmentlist = <int* >NULL
    out_t.segmentmarkerlist = <int* >NULL

    out_t.edgelist = <int* >NULL
    out_t.edgemarkerlist = <int* >NULL

    out_t.holelist = <REAL* >NULL
    out_t.regionlist = <REAL* >NULL

    triangulate(mod, &in_t, &out_t, <triangulateio* >NULL)

    dimensions = np.zeros(2,dtype=np.int32)

    dimensions[0] = out_t.numberoftriangles
    dimensions[1] = 3
    gentrianglelist = c_int_array_to_numpy_array(out_t.trianglelist, dimensions[0]*dimensions[1], dimensions)

    dimensions[0] = out_t.numberofpoints
    dimensions[1] = 2
    genpointlist = c_double_array_to_numpy_array(out_t.pointlist, dimensions[0]*dimensions[1], dimensions)

    dimensions[0] = out_t.numberofpoints
    genpointmarkerlist = c_int_array_to_numpy_array(out_t.pointmarkerlist, dimensions[0], dimensions[0:1])

    dimensions[0] = out_t.numberofpoints
    dimensions[1] = out_t.numberofpointattributes
    genpointattributelist = c_double_array_to_numpy_array(out_t.pointattributelist, dimensions[0]*dimensions[1], dimensions)

    dimensions[0] = out_t.numberoftriangles
    dimensions[1] = out_t.numberoftriangleattributes
    gentriangleattributelist = c_double_array_to_numpy_array(out_t.triangleattributelist, dimensions[0]*dimensions[1], dimensions)

    dimensions[0] = out_t.numberofsegments
    dimensions[1] = 2
    gensegmentlist = c_int_array_to_numpy_array(out_t.segmentlist, dimensions[0]*dimensions[1], dimensions)

    dimensions[0] = out_t.numberofsegments
    gensegmentmarkerlist = c_int_array_to_numpy_array(out_t.segmentmarkerlist, dimensions[0], dimensions[0:1])

    if out_t.neighborlist != NULL:
        dimensions[0] = out_t.numberoftriangles
        dimensions[1] = 3
        genneighborlist = c_int_array_to_numpy_array(out_t.neighborlist, dimensions[0]*dimensions[1], dimensions)
    else:
        genneighborlist = np.zeros((0,0), dtype=np.int32)

    if not(out_t.trianglearealist):
        free(out_t.trianglearealist)
        out_t.trianglearealist = NULL

    if not(out_t.edgelist):
        free(out_t.edgelist)
        out_t.edgelist = NULL

    if not(out_t.edgemarkerlist):
        free(out_t.edgemarkerlist)
        out_t.edgemarkerlist = NULL

    if not(out_t.holelist):
        free(out_t.holelist)
        out_t.holelist = NULL

    if not(out_t.regionlist):
        free(out_t.regionlist)
        out_t.regionlist = NULL

    return gentrianglelist, genpointlist, genpointmarkerlist, genpointattributelist, gentriangleattributelist, gensegmentlist, gensegmentmarkerlist, genneighborlist

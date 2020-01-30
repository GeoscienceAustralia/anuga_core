#cython: wraparound=False, boundscheck=False, cdivision=True, profile=False, nonecheck=False, overflowcheck=False, cdivision_warnings=False, unraisable_tracebacks=False
import cython
from libc.stdlib cimport malloc, free
# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

np.import_array() # avoid segmentation fault

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

def genMesh(np.ndarray pointlist not None,\
            np.ndarray seglist not None,\
            np.ndarray holelist not None,\
            np.ndarray regionlist not None,\
            np.ndarray pointattributelist not None,\
            np.ndarray segmarkerlist not None,\
            char* mod):

    cdef triangulateio in_t, out_t, vorout_t

    cdef np.npy_intp* dimensions

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

    pointlist = np.ascontiguousarray(pointlist)
    seglist = np.ascontiguousarray(seglist)
    holelist = np.ascontiguousarray(holelist)
    regionlist = np.ascontiguousarray(regionlist)
    pointattributelist = np.ascontiguousarray(pointattributelist)
    segmarkerlist = np.ascontiguousarray(segmarkerlist)

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
    in_t.numberoftriangleattributes = 0
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

    vorout_t.pointlist = <REAL* >NULL
    vorout_t.pointmarkerlist = <int* >NULL
    vorout_t.pointattributelist = <REAL* >NULL

    vorout_t.trianglelist = <int* >NULL
    vorout_t.triangleattributelist = <REAL* >NULL
    vorout_t.trianglearealist = <REAL* >NULL
    vorout_t.neighborlist = <int* >NULL

    vorout_t.segmentlist = <int* >NULL
    vorout_t.segmentmarkerlist = <int* >NULL

    vorout_t.edgelist = <int* >NULL
    vorout_t.edgemarkerlist = <int* >NULL

    vorout_t.holelist = <REAL* >NULL
    vorout_t.regionlist = <REAL* >NULL   

    #print ('Before triangulate')

    triangulate(mod, &in_t, &out_t, &vorout_t)

    #print('After triangulate')

    dimensions = <np.npy_intp* > malloc(2 * sizeof(np.npy_intp))

    dimensions[0] = out_t.numberoftriangles
    dimensions[1] = 3
    gentrianglelist = np.PyArray_SimpleNewFromData(2, dimensions, np.NPY_INT32, out_t.trianglelist)

    dimensions[0] = out_t.numberofpoints
    dimensions[1] = 2
    genpointlist = np.PyArray_SimpleNewFromData(2, dimensions, np.NPY_DOUBLE, out_t.pointlist)

    dimensions[0] = out_t.numberofpoints
    genpointmarkerlist = np.PyArray_SimpleNewFromData(1, dimensions, np.NPY_INT32, out_t.pointmarkerlist)

    dimensions[0] = out_t.numberofpoints
    dimensions[1] = out_t.numberofpointattributes
    if out_t.pointattributelist != NULL:
        genpointattributelist = np.PyArray_SimpleNewFromData(2, dimensions, np.NPY_DOUBLE, out_t.pointattributelist)
    else:
        genpointattributelist = np.empty((dimensions[0],dimensions[1]), dtype=np.float)

    dimensions[0] = out_t.numberoftriangles
    dimensions[1] = out_t.numberoftriangleattributes
    if out_t.triangleattributelist != NULL:
        gentriangleattributelist = np.PyArray_SimpleNewFromData(2, dimensions, np.NPY_DOUBLE, out_t.triangleattributelist)
    else:
        gentriangleattributelist = np.empty((dimensions[0],dimensions[1]), dtype=np.float)    

    dimensions[0] = out_t.numberofsegments
    dimensions[1] = 2
    if out_t.segmentlist != NULL:
        gensegmentlist = np.PyArray_SimpleNewFromData(2, dimensions, np.NPY_INT32, out_t.segmentlist)
    else:
        gensegmentlist = np.empty((dimensions[0],dimensions[1]), dtype=np.int32)
    
    dimensions[0] = out_t.numberofsegments
    gensegmentmarkerlist = np.PyArray_SimpleNewFromData(1, dimensions, np.NPY_INT32, out_t.segmentmarkerlist)

    dimensions[0] = out_t.numberoftriangles
    dimensions[1] = 3
    if out_t.neighborlist != NULL:
        genneighborlist = np.PyArray_SimpleNewFromData(2, dimensions, np.NPY_INT32, out_t.neighborlist)
    else:
        genneighborlist = np.empty((dimensions[0],dimensions[1]), dtype=np.int32)

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


    if not(vorout_t.trianglearealist):
        free(vorout_t.trianglearealist)
        vorout_t.trianglearealist = NULL

    if not(vorout_t.edgelist):
        free(vorout_t.edgelist)
        vorout_t.edgelist = NULL

    if not(vorout_t.edgemarkerlist):
        free(vorout_t.edgemarkerlist)
        vorout_t.edgemarkerlist = NULL

    if not(vorout_t.holelist):
        free(vorout_t.holelist)
        vorout_t.holelist = NULL

    if not(vorout_t.regionlist):
        free(vorout_t.regionlist)
        vorout_t.regionlist = NULL

    free(dimensions)

    return gentrianglelist,\
            genpointlist,\
            genpointmarkerlist,\
            genpointattributelist,\
            gentriangleattributelist,\
            gensegmentlist,\
            gensegmentmarkerlist,\
            genneighborlist

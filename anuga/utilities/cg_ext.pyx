#cython: wraparound=False, boundscheck=False, cdivision=True, profile=False, nonecheck=False, overflowcheck=False, cdivision_warnings=False, unraisable_tracebacks=False
import cython
from libc.stdint cimport int64_t

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

# declare the interface to the C code
cdef extern from "cg.c":
    int64_t _jacobi_precon_c(double* data, int64_t* colind, int64_t* row_ptr, double* precon, int64_t M)
    int64_t _cg_solve_c(double* data, int64_t* colind, int64_t* row_ptr, double* b, double* x, int64_t imax, double tol, double a_tol, int64_t M)
    int64_t _cg_solve_c_precon(double* data, int64_t* colind, int64_t* row_ptr, double* b, double* x, int64_t imax, double tol, double a_tol, int64_t M, double* precon)

def jacobi_precon_c(object csr_sparse, np.ndarray[double, ndim=1, mode="c"] precon not None):

    cdef int64_t M, err
    cdef np.ndarray[double, ndim=1, mode="c"] data
    cdef np.ndarray[int64_t, ndim=1, mode="c"] colind
    cdef np.ndarray[int64_t, ndim=1, mode="c"] row_ptr

    data = csr_sparse.data
    colind = csr_sparse.colind
    row_ptr = csr_sparse.row_ptr

    M = row_ptr.shape[0] - 1

    err = _jacobi_precon_c(&data[0], &colind[0], &row_ptr[0], &precon[0], M)

def cg_solve_c(object csr_sparse,\
                np.ndarray[double, ndim=1, mode="c"] x0 not None,\
                np.ndarray[double, ndim=1, mode="c"] b not None,\
                int64_t imax,\
                double tol,\
                double a_tol,\
                int64_t bcols):

    cdef int64_t M, err
    cdef np.ndarray[double, ndim=1, mode="c"] data
    cdef np.ndarray[int64_t, ndim=1, mode="c"] colind
    cdef np.ndarray[int64_t, ndim=1, mode="c"] row_ptr

    data = csr_sparse.data
    colind = csr_sparse.colind
    row_ptr = csr_sparse.row_ptr

    M = row_ptr.shape[0] - 1

    err = _cg_solve_c(&data[0],\
                    &colind[0],\
                    &row_ptr[0],\
                    &b[0],\
                    &x0[0],\
                    imax,\
                    tol,\
                    a_tol,\
                    M)

    return err

def cg_solve_c_precon(object csr_sparse,\
                    np.ndarray[double, ndim=1, mode="c"] x0 not None,\
                    np.ndarray[double, ndim=1, mode="c"] b not None,\
                    int64_t imax,\
                    double tol,\
                    double a_tol,\
                    int64_t bcols,\
                    np.ndarray[double, ndim=1, mode="c"] precon not None):

    cdef int64_t M, err
    cdef np.ndarray[double, ndim=1, mode="c"] data
    cdef np.ndarray[int64_t, ndim=1, mode="c"] colind
    cdef np.ndarray[int64_t, ndim=1, mode="c"] row_ptr

    data = csr_sparse.data
    colind = csr_sparse.colind
    row_ptr = csr_sparse.row_ptr

    M = row_ptr.shape[0] - 1

    err = _cg_solve_c_precon(&data[0],\
                            &colind[0],\
                            &row_ptr[0],\
                            &b[0],\
                            &x0[0],\
                            imax,\
                            tol,\
                            a_tol,\
                            M,\
                            &precon[0])

    return err

// Python - C extension for sparse module.
//
// To compile (Python2.3):
//  gcc -c util_ext.c -I/usr/include/python2.3 -o util_ext.o -Wall -O
//  gcc -shared util_ext.o  -o util_ext.so
//
// See the module sparse.py
//
//
// Steve Roberts, ANU 2004

#include "math.h"
#include "stdio.h"
#include <stdint.h>
#include "anuga_typedefs.h"

#include "omp.h"

//Matrix-vector routine
// y must be pre-zeroed by the caller (csr_mv in sparse_ext.pyx uses np.zeros).
// Each row i writes only to y[i], so the parallel loop is race-free.
anuga_int _csr_mv(anuga_int M,
	    double* data,
	    anuga_int* colind,
	    anuga_int* row_ptr,
	    double* x,
	    double* y) {

  anuga_int i, j, ckey;

  #pragma omp parallel for private(ckey, j)
  for (i=0; i<M; i++ )
    for (ckey=row_ptr[i]; ckey<row_ptr[i+1]; ckey++) {
      j = colind[ckey];
      y[i] += data[ckey]*x[j];
    }

  return 0;
}

//Matrix-matrix routine
// y must be pre-zeroed by the caller.
anuga_int _csr_mm(anuga_int M,
	    anuga_int columns,
	    double* data,
	    anuga_int* colind,
	    anuga_int* row_ptr,
	    double* x,
	    double* y) {

  anuga_int i, j, ckey, c, rowind_i, rowind_j;

  #pragma omp parallel for private(ckey, j, c, rowind_i, rowind_j)
  for (i=0; i<M; i++ ) {
    rowind_i = i*columns;

    for (ckey=row_ptr[i]; ckey<row_ptr[i+1]; ckey++) {
      j = colind[ckey];
      rowind_j = j*columns;

      for (c=0; c<columns; c++) {
        y[rowind_i+c] += data[ckey]*x[rowind_j+c];
      }
    }
  }

  return 0;
}

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

//Matrix-vector routine
int64_t _csr_mv(int64_t M,
	    double* data, 
	    int64_t* colind,
	    int64_t* row_ptr,
	    double* x,
	    double* y) {
  		
  int64_t i, j, ckey;

  for (i=0; i<M; i++ ) 
    for (ckey=row_ptr[i]; ckey<row_ptr[i+1]; ckey++) {
      j = colind[ckey];
      y[i] += data[ckey]*x[j];
    }              
  
  return 0;
}            

//Matrix-matrix routine
int64_t _csr_mm(int64_t M,
	    int64_t columns, 
	    double* data, 
	    int64_t* colind,
	    int64_t* row_ptr,
	    double* x,
	    double* y) {
  		
  int64_t i, j, ckey, c, rowind_i, rowind_j;

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


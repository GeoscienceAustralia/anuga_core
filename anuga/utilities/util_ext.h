// Python - C extension for finite_volumes util module.
//
// To compile (Python2.3):
//  gcc -c util_ext.c -I/usr/include/python2.3 -o util_ext.o -Wall -O 
//  gcc -shared util_ext.o  -o util_ext.so	
//
// See the module util.py
//
//
// Ole Nielsen, GA 2004
	
#include "math.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>


#ifndef ANUGA_UTIL_EXT_H
#define ANUGA_UTIL_EXT_H


#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define AT __FILE__ ":" TOSTRING(__LINE__)
#define P_ERROR_BUFFER_SIZE 100

// provided by stdlib.h 
// double max(double x, double y) {  
//   //Return maximum of two doubles
  
//   if (x > y) return x;
//   else return y;
// }

// provided by stdlib.h
// double min(double x, double y) {  
//   //Return minimum of two doubles
  
//   if (x < y) return x;
//   else return y;
// }


double sign(double x) {
  //Return sign of a double

  if (x>0.0) return 1.0;
  else if (x<0.0) return -1.0;
  else return 0.0;
}

int64_t _gradient(double x0, double y0, 
	      double x1, double y1, 
	      double x2, double y2, 
	      double q0, double q1, double q2, 
	      double *a, double *b) {
	      
  /*Compute gradient (a,b) based on three points (x0,y0), (x1,y1) and (x2,y2) 
  with values q0, q1 and q2.
  
  Extrapolation formula (q0 is selected as an arbitrary origin)
    q(x,y) = q0 + a*(x-x0) + b*(y-y0)                    (1)
  
  Substituting the known values for q1 and q2 into (1) yield the 
  equations for a and b 
  
      q1-q0 = a*(x1-x0) + b*(y1-y0)                      (2)
      q2-q0 = a*(x2-x0) + b*(y2-y0)                      (3)      
      
  or in matrix form
  
  /               \  /   \   /       \  
  |  x1-x0  y1-y0 |  | a |   | q1-q0 |
  |               |  |   | = |       | 
  |  x2-x0  y2-y0 |  | b |   | q2-q0 |
  \               /  \   /   \       /
   
  which is solved using the standard determinant technique    
      
  */
	      

  double det;
  
  det = (y2-y0)*(x1-x0) - (y1-y0)*(x2-x0);

  *a = (y2-y0)*(q1-q0) - (y1-y0)*(q2-q0);
  *a /= det;

  *b = (x1-x0)*(q2-q0) - (x2-x0)*(q1-q0);
  *b /= det;

  return 0;
}


int64_t _gradient2(double x0, double y0, 
	       double x1, double y1, 
	       double q0, double q1, 
	       double *a, double *b) {
  /*Compute gradient (a,b) between two points (x0,y0) and (x1,y1) 
  with values q0 and q1 such that the plane is constant in the direction 
  orthogonal to (x1-x0, y1-y0).
  
  Extrapolation formula
    q(x,y) = q0 + a*(x-x0) + b*(y-y0)                    (1)
  
  Substituting the known values for q1 into (1) yields an 
  under determined  equation for a and b 
      q1-q0 = a*(x1-x0) + b*(y1-y0)                      (2)
      
      
  Now add the additional requirement that the gradient in the direction 
  orthogonal to (x1-x0, y1-y0) should be zero. The orthogonal direction 
  is given by the vector (y0-y1, x1-x0).
  
  Define the point (x2, y2) = (x0 + y0-y1, y0 + x1-x0) on the orthognal line. 
  Then we know that the corresponding value q2 should be equal to q0 in order 
  to obtain the zero gradient, hence applying (1) again    
      q0 = q2 = q(x2, y2) = q0 + a*(x2-x0) + b*(y2-y0)
                          = q0 + a*(x0 + y0-y1-x0) + b*(y0 + x1-x0 - y0)
			  = q0 + a*(y0-y1) + b*(x1-x0)
			  
  leads to the orthogonality constraint
     a*(y0-y1) + b*(x1-x0) = 0                           (3) 
     
  which closes the system and yields
  
  /               \  /   \   /       \  
  |  x1-x0  y1-y0 |  | a |   | q1-q0 |
  |               |  |   | = |       | 
  |  y0-y1  x1-x0 |  | b |   |   0   |
  \               /  \   /   \       /
   
  which is solved using the standard determinant technique    
      
  */

  double det, xx, yy, qq;
  
  xx = x1-x0;
  yy = y1-y0;
  qq = q1-q0;
    
  det = xx*xx + yy*yy;  //FIXME  catch det == 0
  *a = xx*qq/det;
  *b = yy*qq/det;
        
  return 0;
}


void _limit_old(int64_t N, double beta, double* qc, double* qv, 
	    double* qmin, double* qmax) { 

  //N are the number of elements
  int64_t k, i, k3;
  double dq, dqa[3], phi, r;
  
  //printf("INSIDE\n");
  for (k=0; k<N; k++) {
    k3 = k*3;
    
    //Find the gradient limiter (phi) across vertices  
    phi = 1.0;
    for (i=0; i<3; i++) {    
      r = 1.0;
      
      dq = qv[k3+i] - qc[k];    //Delta between vertex and centroid values
      dqa[i] = dq;              //Save dq for use in the next loop
      
      if (dq > 0.0) r = (qmax[k] - qc[k])/dq;
      if (dq < 0.0) r = (qmin[k] - qc[k])/dq;      
  
  
      phi = min( min(r*beta, 1.0), phi);    
    }
    
    //Then update using phi limiter
    for (i=0; i<3; i++) {    
      qv[k3+i] = qc[k] + phi*dqa[i];
    }
  }
}


void  print_double_array(char* name, double* array, int64_t n, int64_t m){

    int64_t k,i,km;

    printf("%s = [",name);
    for (k=0; k<n; k++){
	km = m*k;
	printf("[");
	for (i=0; i<m ; i++){
	    printf("%g ",array[km+i]);
	}
	if (k==(n-1))
	    printf("]");
	else
	    printf("]\n");
    }
    printf("]\n");
}

void  print_int_array(char* name, int32_t* array, int64_t n, int64_t m){

    int64_t k,i,km;

    printf("%s = [",name);
    for (k=0; k<n; k++){
	km = m*k;
	printf("[");
	for (i=0; i<m ; i++){
	    printf("%i ",array[km+i]);
	}
	if (k==(n-1))
	    printf("]");
	else
	    printf("]\n");
    }
    printf("]\n");
}


void  print_long_array(char* name, int64_t * array, int64_t n, int64_t m){

    int64_t k,i,km;

    printf("%s = [",name);
    for (k=0; k<n; k++){
	km = m*k;
	printf("[");
	for (i=0; i<m ; i++){
	  printf("%li ",array[km+i]);
	}
	if (k==(n-1))
	    printf("]");
	else
	    printf("]\n");
    }
    printf("]\n");
}

#endif /* ANUGA_UTIL_EXT_H */

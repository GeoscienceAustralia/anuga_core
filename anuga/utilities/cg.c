// C extension for the cg_solve module. Implements a c routine of the 
// conjugate gradient algorithm to solve Ax=b, using sparse matrix A in 
// the csr format.
//
// See the module cg_solve.py
//
// Padarn Wilson 2012
//
// Note Padarn 26/11/12: Currently the input matrix for the CG solve must
// be a Sparse_CSR matrix python object - defined in anuga/utilities/sparse.py
//
// Note Padarn 5/12/12: I have tried a few optimization modifications which
// didn't seem to save any time:
// -- conversion of the int64_t arrays to int64_t arrays (to save memory passing time)
// -- taking advantage of the symmetric quality of the matrix A to reduce the zAx loop
// -- specifying different 'chunk' sizes for the openmp loops
// -- using blas instead of own openmp loops
	
#include "math.h"
#include "stdio.h"
#include <stdint.h>

#if defined(__APPLE__)
   // clang doesn't have openmp
#else
   #include "omp.h"
#endif



// Dot product of two double vectors: a.b
// @input N: int64_t length of vectors a and b
//        a: first vector of doubles
//        b: second vector of double
// @return: double result of a.b 
double cg_ddot( int64_t N, double *a, double *b)
{
  double ret = 0;
  int64_t i;
  #pragma omp parallel for private(i) reduction(+:ret)
  for(i=0;i<N;i++)
  {
    ret+=a[i]*b[i];
  }
  return ret;

}

// In place multiplication of a double vector x by constant a: a*x
// @input N: int64_t length of vector x
//        a: double scalar to multiply by
//        x: double vector to scale
void cg_dscal(int64_t N, double a, double *x)
{
  int64_t i;
  #pragma omp parallel for private(i)
  for(i=0;i<N;i++)
  {
    x[i]=a*x[i];
  }

}

// Copy of one vector to another - memory already allocated: y=x
// @input N: int64_t length of vectors x and y
//        x: double vector to make copy of
//        y: double vector to copy into
void cg_dcopy( int64_t N, double *x, double *y)
{
  int64_t i;
  #pragma omp parallel for private(i)
  for(i=0;i<N;i++)
  {
    y[i]=x[i];
  }
}

// In place axpy operation: y = a*x + y
// @input N: int64_t length of vectors x and y
//        a: double to multiply x by
//        x: first double vector
//        y: second double vector, stores result
void cg_daxpy(int64_t N, double a, double *x, double *y)
{
  int64_t i;
  #pragma omp parallel for private(i)
  for(i=0;i<N;i++)
  {
    y[i]=y[i]+a*x[i];
  }
}

// Sparse CSR matrix-vector product: z = A*x
// @input z: double vector to store the result
//        data: double vector with non-zero entries of A
//        colind: int64_t vector of column indicies of non-zero entries of A
//        row_ptr: int64_t vector giving index of rows for non-zero entires of A
//        x: double vector to be multiplied
//        M: length of vector x
void cg_zAx(double * z, double * data, int64_t * colind, int64_t * row_ptr, double * x, int64_t M){

  
  
  int64_t i, j, ckey;


     
    #pragma omp parallel for private(ckey,j,i)
    for (i=0; i<M; i++){
      z[i]=0;
      for (ckey=row_ptr[i]; ckey<row_ptr[i+1]; ckey++) {
        j = colind[ckey];
        z[i] += data[ckey]*x[j];
      }              
    }
  

}

// Diagonal matrix-vector product: z = D*x
// @input z: double vector to store the result
//        D: double vector of diagonal matrix
//        x: double vector to be multiplied
//        M: length of vector x
void cg_zDx(double * z, double * D, double * x, int64_t M){

  
  int64_t i, j, ckey;
   
    #pragma omp parallel for private(ckey,j,i)
    for (i=0; i<M; i++){
      z[i]=D[i]*x[i];              
    }
  

}

// Diagonal matrix-vector product: z = D*x
// @input z: double vector to store the result
//        D: double vector of diagonal matrix
//        x: double vector to be multiplied
//        M: length of vector x
void cg_zDinx(double * z, double * D, double * x, int64_t M){

  
  int64_t i, j, ckey;
   
    #pragma omp parallel for private(ckey,j,i)
    for (i=0; i<M; i++){
      z[i]=1.0/D[i]*x[i];              
    }
  

}



// Sparse CSR matrix-vector product and vector addition: z = a*A*x + y
// @input z: double vector to store the result
//        a: double to scale matrix-vector product by
//        data: double vector with non-zero entries of A
//        colind: int64_t vector of column indicies of non-zero entries of A
//        row_ptr: int64_t vector giving index of rows for non-zero entires of A
//        x: double vector to be multiplied
//        y: double vector to add
//        M: length of vector x
void cg_zaAxpy(double * z, double a, double * data, int64_t * colind, int64_t * row_ptr, double * x,
      double * y,int64_t M){
  int64_t i, j, ckey;
  #pragma omp parallel for private(ckey,j,i)
    for (i=0; i<M; i++ ){
      z[i]=y[i];

      for (ckey=row_ptr[i]; ckey<row_ptr[i+1]; ckey++) {
        j = colind[ckey];
        z[i] += a*data[ckey]*x[j];
      }              

  }

}

// Jacobi preconditioner for matrix, A, and right hand side, b. Mutiplies each row
// by one divided by the diagonal element of the matrix. If the diagonal 
// element is zero, does nothing (should nnot occur)
//        colind: int64_t vector of column indicies of non-zero entries of A
//        row_ptr: int64_t vector giving index of rows for non-zero entires of A
//        b: double vector specifying right hand side of equation to solve
//        M: length of vector b

int64_t _jacobi_precon_c(double* data, 
                int64_t* colind,
                int64_t* row_ptr,
                double * precon,
                int64_t M){


  int64_t i, j, k, ckey;
  double diag;


     
  #pragma omp parallel for private(diag,ckey,j,i)
  for (i=0; i<M; i++){
    diag = 0;
    for (ckey=row_ptr[i]; ckey<row_ptr[i+1]; ckey++) {
      j = colind[ckey];
      if (i==j){
        diag = data[ckey];
      }
    }
    if (diag == 0){
      diag =1;
    }
    precon[i]=diag;
  }
  
  return 0;

}

// Conjugate gradient solve Ax = b for x, A given in Sparse CSR format
// @input data: double vector with non-zero entries of A
//        colind: int64_t vector of column indicies of non-zero entries of A
//        row_ptr: int64_t vector giving index of rows for non-zero entires of A
//        b: double vector specifying right hand side of equation to solve
//        x: double vector with initial guess and to store result
//        imax: maximum number of iterations
//        tol: error tollerance for stopping criteria
//        M: length of vectors x and b
// @return: 0 on success  
int64_t _cg_solve_c(double* data, 
                int64_t* colind,
                int64_t* row_ptr,
                double * b,
                double * x,
                int64_t imax,
                double tol,
                double a_tol,
                int64_t M){

  int64_t i = 1;
  double alpha,rTr,rTrOld,bt,rTr0;

  double * d = malloc(sizeof(double)*M);
  double * r = malloc(sizeof(double)*M);
  double * q = malloc(sizeof(double)*M);
  double * xold = malloc(sizeof(double)*M);

  cg_zaAxpy(r,-1.0,data,colind,row_ptr,x,b,M);
  cg_dcopy(M,r,d);

  rTr=cg_ddot(M,r,r);
  rTr0 = rTr;
  
  while((i<imax) && (rTr>pow(tol,2)*rTr0) && (rTr > pow(a_tol,2))){

    cg_zAx(q,data,colind,row_ptr,d,M);
    alpha = rTr/cg_ddot(M,d,q);
    cg_dcopy(M,x,xold);
    cg_daxpy(M,alpha,d,x);

    cg_daxpy(M,-alpha,q,r);
    rTrOld = rTr;
    rTr = cg_ddot(M,r,r);

    bt= rTr/rTrOld;

    cg_dscal(M,bt,d);
    cg_daxpy(M,1.0,r,d);

    i=i+1;

  }
  
  free(d);
  free(r);
  free(q);
  free(xold);

  if (i>=imax){
    return -1;
  }
  else{
    return 0;
  }
  

}          

// Conjugate gradient solve Ax = b for x, A given in Sparse CSR format,
// using a diagonal preconditioner M. 
// @input data: double vector with non-zero entries of A
//        colind: int64_t vector of column indicies of non-zero entries of A
//        row_ptr: int64_t vector giving index of rows for non-zero entires of A
//        b: double vector specifying right hand side of equation to solve
//        x: double vector with initial guess and to store result
//        imax: maximum number of iterations
//        tol: error tollerance for stopping criteria
//        M: length of vectors x and b
//        precon: diagonal preconditioner given as vector
// @return: 0 on success  
int64_t _cg_solve_c_precon(double* data, 
                int64_t* colind,
                int64_t* row_ptr,
                double * b,
                double * x,
                int64_t imax,
                double tol,
                double a_tol,
                int64_t M,
                double * precon){

  int64_t i = 1;
  double alpha,rTr,rTrOld,bt,rTr0;

  double * d = malloc(sizeof(double)*M);
  double * r = malloc(sizeof(double)*M);
  double * q = malloc(sizeof(double)*M);
  double * xold = malloc(sizeof(double)*M);
  double * rhat = malloc(sizeof(double)*M);
  double * temp = malloc(sizeof(double)*M);

  cg_zaAxpy(r,-1.0,data,colind,row_ptr,x,b,M);
  cg_zDinx(rhat,precon,r,M);
  cg_dcopy(M,rhat,d);

  rTr=cg_ddot(M,r,rhat);
  rTr0 = rTr;
  
  while((i<imax) && (rTr>pow(tol,2)*rTr0) && (rTr > pow(a_tol,2))){

    cg_zAx(q,data,colind,row_ptr,d,M);
    alpha = rTr/cg_ddot(M,d,q);
    cg_dcopy(M,x,xold);
    cg_daxpy(M,alpha,d,x);

    cg_daxpy(M,-alpha,q,r);
    cg_zDinx(rhat,precon,r,M);
    rTrOld = rTr;
    rTr = cg_ddot(M,r,rhat);

    bt= rTr/rTrOld;

    cg_dscal(M,bt,d);
    cg_daxpy(M,1.0,rhat,d);

    i=i+1;

  }
  free(temp);
  free(rhat);
  free(d);
  free(r);
  free(q);
  free(xold);

  if (i>=imax){
    return -1;
  }
  else{
    return 0;
  }
  

}       


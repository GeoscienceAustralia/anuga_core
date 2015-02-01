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
// -- conversion of the long arrays to int arrays (to save memory passing time)
// -- taking advantage of the symmetric quality of the matrix A to reduce the zAx loop
// -- specifying different 'chunk' sizes for the openmp loops
// -- using blas instead of own openmp loops
	
#include "Python.h"
#include "numpy/arrayobject.h"
#include "math.h"
#include "stdio.h"
#include "omp.h"


// Dot product of two double vectors: a.b
// @input N: int length of vectors a and b
//        a: first vector of doubles
//        b: second vector of double
// @return: double result of a.b 
double ddot( int N, double *a, double *b)
{
  double ret = 0;
  int i;
  #pragma omp parallel for private(i) reduction(+:ret)
  for(i=0;i<N;i++)
  {
    ret+=a[i]*b[i];
  }
  return ret;

}

// In place multiplication of a double vector x by constant a: a*x
// @input N: int length of vector x
//        a: double scalar to multiply by
//        x: double vector to scale
void dscal(int N, double a, double *x)
{
  int i;
  #pragma omp parallel for private(i)
  for(i=0;i<N;i++)
  {
    x[i]=a*x[i];
  }

}

// Copy of one vector to another - memory already allocated: y=x
// @input N: int length of vectors x and y
//        x: double vector to make copy of
//        y: double vector to copy into
void dcopy( int N, double *x, double *y)
{
  int i;
  #pragma omp parallel for private(i)
  for(i=0;i<N;i++)
  {
    y[i]=x[i];
  }
}

// In place axpy operation: y = a*x + y
// @input N: int length of vectors x and y
//        a: double to multiply x by
//        x: first double vector
//        y: second double vector, stores result
void daxpy(int N, double a, double *x, double *y)
{
  int i;
  #pragma omp parallel for private(i)
  for(i=0;i<N;i++)
  {
    y[i]=y[i]+a*x[i];
  }
}

// Sparse CSR matrix-vector product: z = A*x
// @input z: double vector to store the result
//        data: double vector with non-zero entries of A
//        colind: long vector of column indicies of non-zero entries of A
//        row_ptr: long vector giving index of rows for non-zero entires of A
//        x: double vector to be multiplied
//        M: length of vector x
void zAx(double * z, double * data, long * colind, long * row_ptr, double * x, int M){

  
  
  long i, j, ckey;


     
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
void zDx(double * z, double * D, double * x, int M){

  
  long i, j, ckey;
   
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
void zDinx(double * z, double * D, double * x, int M){

  
  long i, j, ckey;
   
    #pragma omp parallel for private(ckey,j,i)
    for (i=0; i<M; i++){
      z[i]=1.0/D[i]*x[i];              
    }
  

}



// Sparse CSR matrix-vector product and vector addition: z = a*A*x + y
// @input z: double vector to store the result
//        a: double to scale matrix-vector product by
//        data: double vector with non-zero entries of A
//        colind: long vector of column indicies of non-zero entries of A
//        row_ptr: long vector giving index of rows for non-zero entires of A
//        x: double vector to be multiplied
//        y: double vector to add
//        M: length of vector x
void zaAxpy(double * z, double a, double * data, long * colind, long * row_ptr, double * x,
      double * y,int M){
  long i, j, ckey;
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
//        colind: long vector of column indicies of non-zero entries of A
//        row_ptr: long vector giving index of rows for non-zero entires of A
//        b: double vector specifying right hand side of equation to solve
//        M: length of vector b

int _jacobi_precon_c(double* data, 
                long* colind,
                long* row_ptr,
                double * precon,
                int M){


  long i, j, k, ckey;
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
//        colind: long vector of column indicies of non-zero entries of A
//        row_ptr: long vector giving index of rows for non-zero entires of A
//        b: double vector specifying right hand side of equation to solve
//        x: double vector with initial guess and to store result
//        imax: maximum number of iterations
//        tol: error tollerance for stopping criteria
//        M: length of vectors x and b
// @return: 0 on success  
int _cg_solve_c(double* data, 
                long* colind,
                long* row_ptr,
                double * b,
                double * x,
                int imax,
                double tol,
                double a_tol,
                int M){

  int i = 1;
  double alpha,rTr,rTrOld,bt,rTr0;

  double * d = malloc(sizeof(double)*M);
  double * r = malloc(sizeof(double)*M);
  double * q = malloc(sizeof(double)*M);
  double * xold = malloc(sizeof(double)*M);

  zaAxpy(r,-1.0,data,colind,row_ptr,x,b,M);
  dcopy(M,r,d);

  rTr=ddot(M,r,r);
  rTr0 = rTr;
  
  while((i<imax) && (rTr>pow(tol,2)*rTr0) && (rTr > pow(a_tol,2))){

    zAx(q,data,colind,row_ptr,d,M);
    alpha = rTr/ddot(M,d,q);
    dcopy(M,x,xold);
    daxpy(M,alpha,d,x);

    daxpy(M,-alpha,q,r);
    rTrOld = rTr;
    rTr = ddot(M,r,r);

    bt= rTr/rTrOld;

    dscal(M,bt,d);
    daxpy(M,1.0,r,d);

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
//        colind: long vector of column indicies of non-zero entries of A
//        row_ptr: long vector giving index of rows for non-zero entires of A
//        b: double vector specifying right hand side of equation to solve
//        x: double vector with initial guess and to store result
//        imax: maximum number of iterations
//        tol: error tollerance for stopping criteria
//        M: length of vectors x and b
//        precon: diagonal preconditioner given as vector
// @return: 0 on success  
int _cg_solve_c_precon(double* data, 
                long* colind,
                long* row_ptr,
                double * b,
                double * x,
                int imax,
                double tol,
                double a_tol,
                int M,
                double * precon){

  int i = 1;
  double alpha,rTr,rTrOld,bt,rTr0;

  double * d = malloc(sizeof(double)*M);
  double * r = malloc(sizeof(double)*M);
  double * q = malloc(sizeof(double)*M);
  double * xold = malloc(sizeof(double)*M);
  double * rhat = malloc(sizeof(double)*M);
  double * temp = malloc(sizeof(double)*M);

  zaAxpy(r,-1.0,data,colind,row_ptr,x,b,M);
  zDinx(rhat,precon,r,M);
  dcopy(M,rhat,d);

  rTr=ddot(M,r,rhat);
  rTr0 = rTr;
  
  while((i<imax) && (rTr>pow(tol,2)*rTr0) && (rTr > pow(a_tol,2))){

    zAx(q,data,colind,row_ptr,d,M);
    alpha = rTr/ddot(M,d,q);
    dcopy(M,x,xold);
    daxpy(M,alpha,d,x);

    daxpy(M,-alpha,q,r);
    zDinx(rhat,precon,r,M);
    rTrOld = rTr;
    rTr = ddot(M,r,rhat);

    bt= rTr/rTrOld;

    dscal(M,bt,d);
    daxpy(M,1.0,rhat,d);

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

		     
/////////////////////////////////////////////////
// Gateways to Python

PyObject *jacobi_precon_c(PyObject *self, PyObject *args){

  int M,err,bcols;
  
  PyObject *csr_sparse; // input sparse matrix (must be CSR format)
 
  PyArrayObject 
    *data,            //Non Zeros Data array
    *colind,          //Column indices array
    *row_ptr,         //Row pointers array
    *precon;                //Right hand side

  
  // Convert Python arguments to C  
  if (!PyArg_ParseTuple(args, "OO", &csr_sparse, &precon)) {
    PyErr_SetString(PyExc_RuntimeError, "jacobi_precon_c could not parse input");  
    return NULL;
  }

  // Extract three subarrays making up the sparse matrix in CSR format.
  data = (PyArrayObject*) 
    PyObject_GetAttrString(csr_sparse, "data");     
  if (!data) {
    PyErr_SetString(PyExc_RuntimeError, 
        "Data array could not be allocated in cg_solve_c");      
    return NULL;
  }  

  colind = (PyArrayObject*)
    PyObject_GetAttrString(csr_sparse, "colind"); 
  if (!colind) {
    PyErr_SetString(PyExc_RuntimeError, 
        "Column index array could not be allocated in cg_solve_c");      
    return NULL;
  }  

  row_ptr = (PyArrayObject*) 
    PyObject_GetAttrString(csr_sparse, "row_ptr");   
  if (!row_ptr) {
    PyErr_SetString(PyExc_RuntimeError, 
        "Row pointer array could not be allocated in cg_solve_c"); 
  }
  
  M = (row_ptr -> dimensions[0])-1;

  // Precondition the matrix and right hand side
  err = _jacobi_precon_c((double*) data->data, 
                (long*) colind->data,
                (long*) row_ptr->data,
                (double *) precon->data,
                M);

  // Free extra references to sparse matrix parts
  Py_DECREF(data);    
  Py_DECREF(colind);    
  Py_DECREF(row_ptr);     

  return Py_BuildValue("");


}

PyObject *cg_solve_c(PyObject *self, PyObject *args) {
  
  
  PyObject *csr_sparse; // input sparse matrix (must be CSR format)
 
  int imax,M,err,bcols;
  double tol,a_tol;

  PyArrayObject 
    *data,            //Non Zeros Data array
    *colind,          //Column indices array
    *row_ptr,         //Row pointers array
    *x0,               //Initial guess - and sotrage of result.
    *b;                //Right hand side

  
  // Convert Python arguments to C  
  if (!PyArg_ParseTuple(args, "OOOiddi", &csr_sparse, &x0, &b, &imax, &tol, &a_tol, &bcols)) {
    PyErr_SetString(PyExc_RuntimeError, "cg_solve_c could not parse input");  
    return NULL;
  }

  // Extract three subarrays making up the sparse matrix in CSR format.
  data = (PyArrayObject*) 
    PyObject_GetAttrString(csr_sparse, "data");     
  if (!data) {
    PyErr_SetString(PyExc_RuntimeError, 
		    "Data array could not be allocated in cg_solve_c");      
    return NULL;
  }  

  colind = (PyArrayObject*)
    PyObject_GetAttrString(csr_sparse, "colind"); 
  if (!colind) {
    PyErr_SetString(PyExc_RuntimeError, 
		    "Column index array could not be allocated in cg_solve_c");      
    return NULL;
  }  

  row_ptr = (PyArrayObject*) 
    PyObject_GetAttrString(csr_sparse, "row_ptr");   
  if (!row_ptr) {
    PyErr_SetString(PyExc_RuntimeError, 
		    "Row pointer array could not be allocated in cg_solve_c"); 
  }
  
  M = (row_ptr -> dimensions[0])-1;
  
  // Solve system uisng conjugate gradient
  err = _cg_solve_c((double*) data->data, 
                (long*) colind->data,
                (long*) row_ptr->data,
                (double *) b->data,
                (double *) x0->data,
                imax,
                tol,
                a_tol,
                M);
  
  // Free extra references to sparse matrix parts
  Py_DECREF(data);    
  Py_DECREF(colind);    
  Py_DECREF(row_ptr);     

  return Py_BuildValue("i",err);
}

PyObject *cg_solve_c_precon(PyObject *self, PyObject *args) {
  
  
  PyObject *csr_sparse; // input sparse matrix (must be CSR format)
 
  int imax,M,err,bcols;
  double tol,a_tol;

  PyArrayObject 
    *data,            //Non Zeros Data array
    *colind,          //Column indices array
    *row_ptr,         //Row pointers array
    *x0,               //Initial guess - and sotrage of result.
    *b,                //Right hand side
    *precon;           //diagonal preconditioner

  
  // Convert Python arguments to C  
  if (!PyArg_ParseTuple(args, "OOOiddiO", &csr_sparse, &x0, &b, &imax, &tol, &a_tol, &bcols, &precon)) {
    PyErr_SetString(PyExc_RuntimeError, "cg_solve_c_precon could not parse input");  
    return NULL;
  }

  // Extract three subarrays making up the sparse matrix in CSR format.
  data = (PyArrayObject*) 
    PyObject_GetAttrString(csr_sparse, "data");     
  if (!data) {
    PyErr_SetString(PyExc_RuntimeError, 
        "Data array could not be allocated in cg_solve_c");      
    return NULL;
  }  

  colind = (PyArrayObject*)
    PyObject_GetAttrString(csr_sparse, "colind"); 
  if (!colind) {
    PyErr_SetString(PyExc_RuntimeError, 
        "Column index array could not be allocated in cg_solve_c");      
    return NULL;
  }  

  row_ptr = (PyArrayObject*) 
    PyObject_GetAttrString(csr_sparse, "row_ptr");   
  if (!row_ptr) {
    PyErr_SetString(PyExc_RuntimeError, 
        "Row pointer array could not be allocated in cg_solve_c"); 
  }
  
  M = (row_ptr -> dimensions[0])-1;
  
  // Solve system uisng conjugate gradient
  err = _cg_solve_c_precon((double*) data->data, 
                (long*) colind->data,
                (long*) row_ptr->data,
                (double *) b->data,
                (double *) x0->data,
                imax,
                tol,
                a_tol,
                M,
                (double *) precon->data);
  
  // Free extra references to sparse matrix parts
  Py_DECREF(data);    
  Py_DECREF(colind);    
  Py_DECREF(row_ptr);     

  return Py_BuildValue("i",err);
}



// Method table for python module
static struct PyMethodDef MethodTable[] = {
  {"cg_solve_c", cg_solve_c, METH_VARARGS, "Print out"},
  {"cg_solve_c_precon", cg_solve_c_precon, METH_VARARGS, "Print out"},
  {"jacobi_precon_c", jacobi_precon_c, METH_VARARGS, "Print out"},    
  {NULL, NULL, 0, NULL}   /* sentinel */
};

// Module initialisation   
void initcg_ext(void){
  Py_InitModule("cg_ext", MethodTable);
  
  import_array();     //Necessary for handling of NumPY structures  
}


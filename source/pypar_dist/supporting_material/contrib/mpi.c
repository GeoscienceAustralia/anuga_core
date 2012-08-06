/************************************************************************/
/* PyPAR - Parallel Python using MPI                 	                */
/* Copyright (C) 2001, 2002 Ole M. Nielsen, Gian Paolo Ciceri           */
/*                                                                      */
/* See enclosed README file for details of installation and use.   	*/
/*                                                                 	*/   
/* This program is free software; you can redistribute it and/or modify */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation; either version 2 of the License, or    */
/* (at your option) any later version.                                  */
/*                                                                      */     
/* This program is distributed in the hope that it will be useful,      */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/* GNU General Public License (http://www.gnu.org/copyleft/gpl.html)    */
/* for more details.                                                    */
/*                                                                      */
/* You should have received a copy of the GNU General Public License    */
/* along with this program; if not, write to the Free Software          */
/*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307*/
/*                                                                      */
/*                                                                 	*/
/* Contact addresses: Ole.Nielsen@anu.edu.au, gp.ciceri@acm.org         */
/*                                                                 	*/
/* version 1.2.1, 16 February 2002                                      */
/*   Status block, MPI_ANY_TAG, MPI_ANY_SOURCE exported                 */
/* Version 1.2, 15 February 2002                                       	*/   
/*   Scatter added by Gian Paolo Ciceri                                 */
/* Version 1.1, 14 February 2002                                       	*/   
/*   Bcast added by Gian Paolo Ciceri                                   */
/* Version 1.0.2, 10 February 2002                                      */
/*   Modified by Gian Paulo Ciceri to allow pypar run under Python 2.2  */
/* Version 1.0.1, 8 February 2002                                       */
/*   Modified to install on SUN enterprise systems under Mpich          */
/* Version 1.0, 7 February 2002                                         */
/*   First public release for Python 2.1 (OMN)                          */
/************************************************************************/


#include "Python.h"
#include "mpi.h"
#include "Numeric/arrayobject.h"

// to handle MPI constants export (shamelessly stolen from _cursesmodule.c)
#define SetDictInt(string,ch) \
        PyDict_SetItemString(ModDict, string, PyInt_FromLong((long) (ch)));

// kludge to remap struct MPI_op to int (easier to handle by python)
#define mpi_MAX 1
#define mpi_MIN 2
#define mpi_SUM 3
#define mpi_PROD 4
#define mpi_LAND 5
#define mpi_BAND 6
#define mpi_LOR 7
#define mpi_BOR 8
#define mpi_LXOR 9
#define mpi_BXOR 10
#define mpi_MAXLOC 11
#define mpi_MINLOC 12
#define mpi_REPLACE 13



MPI_Datatype type_map(PyArrayObject *x) {  

  //
  // TYPE    py_type  mpi_type  bytes  symbol
  // ---------------------------------------- 
  // INT       4        6         4      'i'
  // LONG      5        8         8      'l'
  // FLOAT     6       10         4      'f'  
  // DOUBLE    7       11         8      'd'
  
  
  int py_type;
  MPI_Datatype mpi_type;
  
  if (x -> nd != 1) {
    PyErr_SetString(PyExc_ValueError, "Array must be 1 dimensional");
    return (MPI_Datatype) 0;
  }      
      
  py_type = x -> descr -> type_num;     
  if (py_type == PyArray_DOUBLE) 
    mpi_type = MPI_DOUBLE;
  else if (py_type == PyArray_LONG)   
    mpi_type = MPI_LONG;  
  else if (py_type == PyArray_FLOAT) 
    mpi_type = MPI_FLOAT;
  else if (py_type == PyArray_INT) 
    mpi_type = MPI_INT;
  else {
    PyErr_SetString(PyExc_ValueError, "Array must be of type int or float");
    return 0;
  }      

  //printf("Types: %d %d.\n", py_type, mpi_type);
  
  return mpi_type;
}    

MPI_Op op_map(int py_op) {  
  
  MPI_Op mpi_op;
  
  if (py_op == mpi_MAX) 
    mpi_op = MPI_MAX;
  else if (py_op == mpi_MIN)   
    mpi_op = MPI_MIN;  
  else if (py_op == mpi_SUM)   
    mpi_op = MPI_SUM;  
  else if (py_op == mpi_PROD)   
    mpi_op = MPI_PROD;  
  else if (py_op == mpi_LAND)   
    mpi_op = MPI_LAND;  
  else if (py_op == mpi_BAND)   
    mpi_op = MPI_BAND;  
  else if (py_op == mpi_LOR)   
    mpi_op = MPI_LOR;  
  else if (py_op == mpi_BOR)   
    mpi_op = MPI_BOR;  
  else if (py_op == mpi_LXOR)   
    mpi_op = MPI_LXOR;  
  else if (py_op == mpi_BXOR)   
    mpi_op = MPI_BXOR;  
  else if (py_op == mpi_MAXLOC)
    mpi_op = MPI_MAXLOC;
  else if (py_op == mpi_MINLOC)   
    mpi_op = MPI_MINLOC;  
  else if (py_op == mpi_REPLACE)   
    mpi_op = MPI_REPLACE;  
  else {
    PyErr_SetString(PyExc_ValueError, "Operation unknown");
    return 0;
  }      

  //printf("Op: %d.\n", py_op);
  
  return mpi_op;
}  

/*********************************************************/
/* send_string                                           */
/* Send string of characters                             */
/*                                                       */
/*********************************************************/
static PyObject *send_string(PyObject *self, PyObject *args) {
  char *s;
  int destination, tag, length, err;
 
  /* process the parameters */
  if (!PyArg_ParseTuple(args, "s#ii", &s, &length, &destination, &tag))
    return NULL;
  
  /* call the MPI routine */
  err = MPI_Send(s, length, MPI_CHAR, destination, tag, MPI_COMM_WORLD);

  return Py_BuildValue("i", err);
}

/**********************************************************/
/* receive_string                                         */
/* Receive string of characters                           */
/*                                                        */
/**********************************************************/
static PyObject *receive_string(PyObject *self, PyObject *args) {
  char *s;
  int source, tag, length, err, st_length; 
  MPI_Status status;

  /* process the parameters */
  if (!PyArg_ParseTuple(args, "s#ii", &s, &length, &source, &tag))
    return NULL;
    
  /* call the MPI routine */
  err = MPI_Recv(s, length, MPI_CHAR, source, tag, MPI_COMM_WORLD, &status);
   
  MPI_Get_count(&status, MPI_CHAR, &st_length); 
  // status.st_length is not available in all MPI implementations
  //Alternative is: MPI_Get_elements(MPI_Status *, MPI_Datatype, int *);

  return Py_BuildValue("i(iiii)", err, status.MPI_SOURCE, status.MPI_TAG,
  status.MPI_ERROR, st_length);  
}

/**********************************************************/
/* bcast_string                                           */
/* Broadcast string of characters                         */
/*                                                        */
/**********************************************************/
static PyObject *bcast_string(PyObject *self, PyObject *args) {
  char *s;
  int source, length, err; 

  /* process the parameters */
  if (!PyArg_ParseTuple(args, "s#i", &s, &length, &source))
    return NULL;
    
  /* call the MPI routine */
  err = MPI_Bcast(s, length, MPI_CHAR, source, MPI_COMM_WORLD);
   
  return Py_BuildValue("i", err);  
}

/**********************************************************/
/* scatter_string                                         */
/* Scatter string of characters                           */
/*                                                        */
/**********************************************************/
static PyObject *scatter_string(PyObject *self, PyObject *args) {
  char *s;
  char *d;
  int source, length, err; 

  /* process the parameters */
  if (!PyArg_ParseTuple(args, "sisi", &s, &length, &d, &source))
    return NULL;
    
  /* call the MPI routine */
  err = MPI_Scatter(s, length, MPI_CHAR, d, length,  MPI_CHAR, source, MPI_COMM_WORLD);
   
  return Py_BuildValue("i", err);  
}

/**********************************************************/
/* gather_string                                         */
/* Gather string of characters                           */
/*                                                        */
/**********************************************************/
static PyObject *gather_string(PyObject *self, PyObject *args) {
  char *s;
  char *d;
  int source, length, err; 

  /* process the parameters */
  if (!PyArg_ParseTuple(args, "sisi", &s, &length, &d, &source))
    return NULL;
    
  /* call the MPI routine */
  err = MPI_Gather(s, length, MPI_CHAR, d, length,  MPI_CHAR, source, MPI_COMM_WORLD);
   
  return Py_BuildValue("i", err);  
}


/**********************************************************/
/* send_array                                             */
/* Send Numeric array of type float, double, int, or long */
/*                                                        */
/**********************************************************/
static PyObject *send_array(PyObject *self, PyObject *args) {
  PyObject *input;
  PyArrayObject *x;
  int destination, tag, err;
  MPI_Datatype mpi_type;
  
  /* process the parameters */
  if (!PyArg_ParseTuple(args, "Oii", &input, &destination, &tag))
    return NULL;
    
  /* Make Numeric array from general sequence type (no cost if already Numeric)*/    
  x = (PyArrayObject *)
    PyArray_ContiguousFromObject(input, PyArray_NOTYPE, 0, 0);
    
  /* Input check and determination of MPI type */          
  mpi_type = type_map(x);
  if (!mpi_type) return NULL;
    
  /* call the MPI routine */
  err = MPI_Send(x->data, x->dimensions[0], mpi_type, destination, tag,\
           MPI_COMM_WORLD);
	   
  Py_DECREF(x); 	   
	   
  return Py_BuildValue("i", err);
}

/*************************************************************/
/* receive_array                                             */
/* Receive Numeric array of type float, double, int, or long */
/*                                                           */
/*************************************************************/
static PyObject *receive_array(PyObject *self, PyObject *args) {
  PyArrayObject *x;
  int source, tag, err, st_length;
  MPI_Datatype mpi_type;
  MPI_Status status;

  /* process the parameters */
  if (!PyArg_ParseTuple(args, "Oii", &x, &source, &tag))
    return NULL;

  /* Input check and determination of MPI type */          
  mpi_type = type_map(x);
  if (!mpi_type) return NULL;  
      
  /* call the MPI routine */
  err =  MPI_Recv(x->data, x->dimensions[0], mpi_type, source, tag, \
         MPI_COMM_WORLD, &status);
	 
  MPI_Get_count(&status, mpi_type, &st_length); 
  // status.st_length is not available in all MPI implementations
  //Alternative is: MPI_Get_elements(MPI_Status *, MPI_Datatype, int *);
	 
      
  return Py_BuildValue("i(iiii)", err, status.MPI_SOURCE, status.MPI_TAG,
  status.MPI_ERROR, st_length);  
}


/*************************************************************/
/* bcast_array                                               */
/* Broadcast Num.  array of type float, double, int, or long */
/*                                                           */
/*************************************************************/
static PyObject *bcast_array(PyObject *self, PyObject *args) {
  PyArrayObject *x;
  int source, err;
  MPI_Datatype mpi_type;
  MPI_Status status;

  /* process the parameters */
  if (!PyArg_ParseTuple(args, "Oi", &x, &source))
    return NULL;

  /* Input check and determination of MPI type */          
  mpi_type = type_map(x);
  if (!mpi_type) return NULL;  
      
  /* call the MPI routine */
  err =  MPI_Bcast(x->data, x->dimensions[0], mpi_type, source, \
         MPI_COMM_WORLD);
      
  return Py_BuildValue("i", err);
}

/*************************************************************/
/* scatter_array                                             */
/* Scatter Num.    array of type float, double, int, or long */
/*                                                           */
/*************************************************************/
static PyObject *scatter_array(PyObject *self, PyObject *args) {
  PyArrayObject *x;
  PyArrayObject *d;
  int length, source, err;
  MPI_Datatype mpi_type;
  MPI_Status status;

  /* process the parameters */
  if (!PyArg_ParseTuple(args, "OiOi", &x, &length, &d, &source))
    return NULL;

  /* Input check and determination of MPI type */          
  mpi_type = type_map(x);
  if (!mpi_type) return NULL;  
      
  /* call the MPI routine */
  err =  MPI_Scatter(x->data, length, mpi_type, d->data, length, mpi_type, source, \
         MPI_COMM_WORLD);
      
  return Py_BuildValue("i", err);
}


/*************************************************************/
/* gather_array                                              */
/* Gather Num.     array of type float, double, int, or long */
/*                                                           */
/*************************************************************/
static PyObject *gather_array(PyObject *self, PyObject *args) {
  PyArrayObject *x;
  PyArrayObject *d;
  int length, source, err;
  MPI_Datatype mpi_type;
  MPI_Status status;

  /* process the parameters */
  if (!PyArg_ParseTuple(args, "OiOi", &x, &length, &d, &source))
    return NULL;

  /* Input check and determination of MPI type */          
  mpi_type = type_map(x);
  if (!mpi_type) return NULL;  
      
  /* call the MPI routine */
  err =  MPI_Gather(x->data, length, mpi_type, d->data, length, mpi_type, source, \
         MPI_COMM_WORLD);
      
  return Py_BuildValue("i", err);
}


/*************************************************************/
/* reduce_array                                              */
/* Reduce Num.     array of type float, double, int, or long */
/*                                                           */
/*************************************************************/
static PyObject *reduce_array(PyObject *self, PyObject *args) {
  PyArrayObject *x;
  PyArrayObject *d;
  int length, source, op, err;
  MPI_Datatype mpi_type;
  MPI_Status status;
  MPI_Op mpi_op;

  /* process the parameters */
  if (!PyArg_ParseTuple(args, "OOiii", &x, &d, &length, &op, &source))
    return NULL;
   
  /* Input check and determination of MPI type */          
  mpi_type = type_map(x);
  if (!mpi_type) return NULL;  
  
  /* Input check and determination of MPI op */ 
  //printf("op: %d\n", op);         
  mpi_op = op_map(op);
  if (!mpi_op) return NULL;  
   
         
  if (op == mpi_MAXLOC || op == mpi_MINLOC) {
    //not implemented
    return Py_BuildValue("i", -666);
  }
  else {
  /* call the MPI routine */
  err =  MPI_Reduce(x->data, d->data, length, mpi_type, mpi_op, source, \
         MPI_COMM_WORLD);
  }
         
      
  return Py_BuildValue("i", err);
}



/*********************************************************/
/* MPI calls rank, size, finalize, abort                 */
/*                                                       */
/*********************************************************/

static PyObject * rank(PyObject *self, PyObject *args) {
  int myid;

  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  return Py_BuildValue("i", myid);
}

static PyObject * size(PyObject *self, PyObject *args) {
  int numprocs; 

  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  return Py_BuildValue("i", numprocs);
}
  
static PyObject * Get_processor_name(PyObject *self, PyObject *args) {  
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int  namelen;

  MPI_Get_processor_name(processor_name,&namelen);
  return Py_BuildValue("s#", processor_name, namelen);
}   
  
static PyObject * Finalize(PyObject *self, PyObject *args) {  
  int error;

  error = MPI_Finalize();
  return Py_BuildValue("i", error);  
} 

static PyObject * Abort(PyObject *self, PyObject *args) {  
  int error, code=0;
  
  error = MPI_Abort(MPI_COMM_WORLD, code);
  return Py_BuildValue("i", error);  
} 

static PyObject * Barrier(PyObject *self, PyObject *args) {  
  int error;
  
  error = MPI_Barrier(MPI_COMM_WORLD);
  return Py_BuildValue("i", error);    
}    

static PyObject * Wtime(PyObject *self, PyObject *args) {     
  double t;
  
  t = MPI_Wtime();
  return Py_BuildValue("d", t);      
}     
 
/**********************************/
/* Method table for python module */
/**********************************/

static struct PyMethodDef MethodTable[] = {
  {"size", size, METH_VARARGS},  
  {"rank", rank, METH_VARARGS},  
  {"Barrier", Barrier, METH_VARARGS},          
  {"Wtime", Wtime, METH_VARARGS},            
  {"Get_processor_name", Get_processor_name, METH_VARARGS},              
  {"Finalize", Finalize, METH_VARARGS},        
  {"Abort", Abort, METH_VARARGS},          
  {"send_string", send_string, METH_VARARGS},
  {"receive_string", receive_string, METH_VARARGS},      
  {"bcast_string", bcast_string, METH_VARARGS},        
  {"scatter_string", scatter_string, METH_VARARGS},        
  {"gather_string", gather_string, METH_VARARGS},        
  {"send_array", send_array, METH_VARARGS},
  {"receive_array", receive_array, METH_VARARGS},    
  {"bcast_array", bcast_array, METH_VARARGS},              
  {"scatter_array", scatter_array, METH_VARARGS},              
  {"gather_array", gather_array, METH_VARARGS},              
  {"reduce_array", reduce_array, METH_VARARGS},              
  {NULL, NULL}
};


/***************************/
/* Initialisation Function */
/***************************/

void initmpi(){
  int error, argc = 0; //Dummy
  char **argv;         //Dummy
  PyObject *m, *ModDict;

  //printf("Initialising MPI\n");
  error = MPI_Init(&argc, &argv); 
  //printf("MPI Initialised\n");  
  
  //FIXME: Make a sensible errorcheck here
  
  m = Py_InitModule("mpi", MethodTable);
  
  // to handle MPI symbolic constants
  ModDict = PyModule_GetDict(m); 
  SetDictInt("MPI_ANY_TAG", MPI_ANY_TAG);
  SetDictInt("MPI_ANY_SOURCE", MPI_ANY_SOURCE);
  SetDictInt("mpi_MAX", mpi_MAX);
  SetDictInt("mpi_MIN", mpi_MIN);
  SetDictInt("mpi_SUM", mpi_SUM);
  SetDictInt("mpi_PROD", mpi_PROD);
  SetDictInt("mpi_LAND", mpi_LAND);
  SetDictInt("mpi_BAND", mpi_BAND);
  SetDictInt("mpi_LOR", mpi_LOR);
  SetDictInt("mpi_BOR", mpi_BOR);
  SetDictInt("mpi_LXOR", mpi_LXOR);
  SetDictInt("mpi_BXOR", mpi_BXOR);
  SetDictInt("mpi_MAXLOC", mpi_MAXLOC);
  SetDictInt("mpi_MINLOC", mpi_MINLOC);
  SetDictInt("mpi_REPLACE", mpi_REPLACE);

  //SetDictInt("MPI_COMM_WORLD", MPI_COMM_WORLD);  
   
  import_array();     //Necessary for handling of NumPY structures  
}

 
 
 
 
/*********************************************************************/
/* OBSOLETE STUFF  					             */
/*********************************************************************/

/******************************************/
/* keep for doc of Numeric arrays etc     */ 
/******************************************/

void print_real_array(PyArrayObject *x) {  
  int i;
  for (i=0; i<x->dimensions[0]; i++) {
    printf("%f ", *(double*) (x->data + i*x->strides[0]));
  }
}

void print_int_array(PyArrayObject *x) {  
  int i;
  for (i=0; i<x->dimensions[0]; i++) {
    printf("%d ", *(int*) (x->data + i*x->strides[0]));
  }
}


/*********************************************************/
/* send_real_array                                       */
/* Send Numeric array of double floating point numbers   */
/*                                                       */
/*********************************************************/
static PyObject *send_real_array(PyObject *self, PyObject *args) {
  PyObject *input;
  PyArrayObject *x;
  int destination, tag, err;
  
  /* process the parameters */
  if (!PyArg_ParseTuple(args, "Oii", &input, &destination, &tag))
    return NULL;

  /* Make Numeric array from general sequence type (no cost if already Numeric)*/    
  x = (PyArrayObject *)
    PyArray_ContiguousFromObject(input, PyArray_DOUBLE, 0, 0);

  /* call the MPI routine */
  err = MPI_Send(x->data, x->dimensions[0], MPI_DOUBLE, destination, tag,\
           MPI_COMM_WORLD);

  Py_DECREF(x); 	   	   	   
  return Py_BuildValue("i", err);
}

/**********************************************************/
/* receive_real_array                                     */
/* Receive Numeric array of double floating point numbers */
/*                                                        */
/**********************************************************/
static PyObject *receive_real_array(PyObject *self, PyObject *args) {
  PyArrayObject *x;
  int source, tag, err;
  MPI_Status status;

  /* process the parameters */
  if (!PyArg_ParseTuple(args, "Oii", &x, &source, &tag))
    return NULL;
  
  /* call the MPI routine */
  err =  MPI_Recv(x->data, x->dimensions[0], MPI_DOUBLE, source, tag, \
         MPI_COMM_WORLD, &status);


  return Py_BuildValue("i", err);
}
 



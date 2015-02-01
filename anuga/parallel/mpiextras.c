/************************************************************************/
/* PyPAR - Parallel Python using MPI                 	                */
/* Copyright (C) 2001 - 2009 Ole Nielsen                                */
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
/* Contact address: Ole.Moller.Nielsen@gmail.com                        */
/*                                                                 	*/
/* version (see __version__ in pypar.py)                                */
/* date (see __date__ in pypar.py)                                      */
/************************************************************************/


#include "Python.h"
#include "mpi.h"
#include "math.h"
#include "numpy/arrayobject.h"


/* Handle MPI constants export (shamelessly stolen from _cursesmodule.c)*/
#define SetDictInt(string,ch) \
        PyDict_SetItemString(ModDict, string, PyInt_FromLong((long) (ch)));

/* Remap struct MPI_op to int (easier to handle by python)*/
#define MAX 1
#define MIN 2
#define SUM 3
#define PROD 4
#define LAND 5
#define BAND 6
#define LOR 7
#define BOR 8
#define LXOR 9
#define BXOR 10
#define MAXLOC 11
#define MINLOC 12
/*#define REPLACE 13 // Not available on all MPI systems */


static char errmsg[132];  /*Used to create exception messages*/

/* MPI_Bsend() related variables. */
//static void *pt_buf;	/* Pointer to allocated buffer. */
//static int buf_size;	/* Size of buffer to allocate. */

int length(PyArrayObject *x) {  
  /* Compute the total length of contiguous array */
  /* Necessary for communicating multi dimensional arrays */
  
  int i, length;
  
  /* Compute total length of contiguous data */
  length = 1;
  for (i=0; i<x->nd; i++) {
    length *= x->dimensions[i];
  }  
    
  return length;
}  

/**********************************************************/
/* type_map                                               */
/*                                                        */
/**********************************************************/
MPI_Datatype type_map(PyArrayObject *x, int *count) {  
  /* Return the MPI Datatype corresponding to                       
     the Python data type as follows  
  
     TYPE    py_type  mpi_type  bytes  symbol
     ---------------------------------------- 
     INT       4        6         4      'i'
     LONG      5        8         8      'l'
     FLOAT     6       10         4      'f'  
     DOUBLE    12      11         8      'd'
  
     Also return the total number of elements in the array
  
     The Python datatype COMPLEX ('F') and COMPLEX_DOUBLE ('D')
     is treated as a special case to the absence of an 
     MPI_COMPLEX datatype:
  
     Complex arrays are mapped to float or double arrays with real 
     and imaginary parts alternating and count is updated. */
  
  int py_type;
  MPI_Datatype mpi_type;
  char err_msg[64];

  *count = length(x);
  
  py_type = PyArray_TYPE(x);
  if (py_type == NPY_DOUBLE) 
    mpi_type = MPI_DOUBLE;
  else if (py_type == NPY_INT) 
    mpi_type = MPI_INT;
  else if (py_type == NPY_CDOUBLE) {
    mpi_type = MPI_DOUBLE;
    (*count) *= 2;
  } else if (py_type == NPY_FLOAT) 
    mpi_type = MPI_FLOAT;
  else if (py_type == NPY_LONG)   
    mpi_type = MPI_LONG;  
  else if (py_type == NPY_CFLOAT) {
    mpi_type = MPI_FLOAT;
    (*count) *= 2;
  } else {
    sprintf(err_msg, 
	    "Array must be of type int or float. I got py_type == %d", 
	    py_type);
    PyErr_SetString(PyExc_ValueError, err_msg);
    return (MPI_Datatype) NULL;
  }      

  //printf("Types: py_type=%d, mpi_type=%d\n", py_type, (int) mpi_type);
  
  return mpi_type;
}    

/**********************************************************/
/* op_map                                                 */
/*                                                        */
/**********************************************************/
MPI_Op op_map(int py_op) {  
  
  MPI_Op mpi_op;
  
  if (py_op == MAX) 
    mpi_op = MPI_MAX;
  else if (py_op == MIN)   
    mpi_op = MPI_MIN;  
  else if (py_op == SUM)   
    mpi_op = MPI_SUM;  
  else if (py_op == PROD)   
    mpi_op = MPI_PROD;  
  else if (py_op == LAND)   
    mpi_op = MPI_LAND;  
  else if (py_op == BAND)   
    mpi_op = MPI_BAND;  
  else if (py_op == LOR)   
    mpi_op = MPI_LOR;  
  else if (py_op == BOR)   
    mpi_op = MPI_BOR;  
  else if (py_op == LXOR)   
    mpi_op = MPI_LXOR;  
  else if (py_op == BXOR)   
    mpi_op = MPI_BXOR;  
  else if (py_op == MAXLOC)
    mpi_op = MPI_MAXLOC;
  else if (py_op == MINLOC)   
    mpi_op = MPI_MINLOC;  
  /*else if (py_op == REPLACE)*/   
  /*  mpi_op = MPI_REPLACE; */ 
  else {
    PyErr_SetString(PyExc_ValueError, "Operation unknown");
    return (MPI_Op) NULL;
  }      
  
  return mpi_op;
}  


/**********************************************************/
/* isend_array (asynchronous)                             */
/* Send Numpy array of type float, double, int, or long   */
/*                                                        */
/**********************************************************/
static PyObject *isend_array(PyObject *self, PyObject *args) {
  PyObject *input;
  PyArrayObject *x;
  int destination, tag, error, count, myid;
  MPI_Datatype mpi_type;

  /* process the parameters */
  if (!PyArg_ParseTuple(args, "Oii", &input, &destination, &tag))
    return NULL;

  /* Make Numpy array from general sequence type (no cost if already Numpy). */
  x = (PyArrayObject *)
    PyArray_ContiguousFromObject(input, NPY_NOTYPE, 0, 0);

  /* Input check and determination of MPI type */
  mpi_type = type_map(x, &count);
  if (!mpi_type) return NULL;

  /* call the MPI routine */
  error = MPI_Send(x->data, count, mpi_type, destination, tag,
		   MPI_COMM_WORLD);
  Py_DECREF(x);

  if (error != 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    sprintf(errmsg, "Proc %d: MPI_Send failed with error code %d\n",
	    myid, error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);
    return NULL;
  }

  Py_INCREF(Py_None);
  return (Py_None);
}




/*************************************************************/
/* ireceive_array (asynchronous)                             */
/* Receive Numpy array of type float, double, int, or long   */
/*                                                           */
/*************************************************************/
static PyObject *ireceive_array(PyObject *self, PyObject *args) {
  PyArrayObject *x;
  int source, tag, error, st_length, size, count, myid;
  MPI_Datatype mpi_type;
  MPI_Status status;

    
  if (!PyArg_ParseTuple(args, "Oii", &x, &source, &tag))
    return NULL;    
    
  /* Input check and determination of MPI type */          
  mpi_type = type_map(x, &count);
  if (!mpi_type) return NULL;  
      
  /* call the MPI routine */
  error =  MPI_Recv(x->data, count, mpi_type, source, tag,
		    MPI_COMM_WORLD, &status);
	 
  /* Do not DECREF x as it must be returned to Python */
  	 
  if (error != 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);    
    sprintf(errmsg, "Proc %d: MPI_Recv failed with error code %d\n", 
	    myid, error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);   
    return NULL;
  }  
   
	 
  MPI_Get_count(&status, mpi_type, &st_length); 
  /* status.st_length is not available in all MPI implementations */
  /* Alternative is: MPI_Get_elements(MPI_Status *, MPI_Datatype, int *); */
	 
      
  /* FIXME: This might not be watertight on all platforms */
  /* Need C equivalent to itemsize().*/
  if (mpi_type == MPI_DOUBLE) {
    size = sizeof(double);  /*8 */
  } else if (mpi_type == MPI_LONG) {
    size = sizeof(long); /*8? */
  } else if (mpi_type == MPI_FLOAT) {
    size = sizeof(float);
  } else if (mpi_type == MPI_INT) {    
    size = sizeof(int);  
  } else {
    size = 4;  
  }
    
  return Py_BuildValue("(iiiii)", status.MPI_SOURCE, status.MPI_TAG,
  status.MPI_ERROR, st_length, size);  
}





/*************************************************************/
/* allreduce_array                                           */
/* Allreduce Numpy array of type float, double, int, or long */
/*                                                           */
/*************************************************************/
static PyObject *allreduce_array(PyObject *self, PyObject *args) {
  PyArrayObject *x;
  PyArrayObject *d;
  int op, error, count, count1, myid;
  MPI_Datatype mpi_type, buffer_type;
  MPI_Op mpi_op;

  /* process the parameters */
  if (!PyArg_ParseTuple(args, "OOi", &x, &d, &op)) {
    PyErr_SetString(PyExc_RuntimeError,
		    "mpiext.c (allreduce_array): could not parse input");
    return NULL;
  }

  /* Input check and determination of MPI type */
  mpi_type = type_map(x, &count);
  if (!mpi_type) {
    PyErr_SetString(PyExc_RuntimeError,
		    "mpiext.c (allreduce_array): could not determine mpi_type");
    return NULL;
  }

  /* This error is caught at the pypar level - so we won't end up here
     unless mpiext is being used independently */
  buffer_type = type_map(d, &count1);
  if (mpi_type != buffer_type) {
    sprintf(errmsg, "mpiext.c (allreduce_array): Input array and buffer must be of the same type.");
    PyErr_SetString(PyExc_RuntimeError, errmsg);

    return NULL;
  }

  if (count != count1) {
    PyErr_SetString(PyExc_RuntimeError,
		    "mpiext.c (allreduce_array): Input array and buffer must have same length");
    return NULL;
  }

  /* Input check and determination of MPI op */
  mpi_op = op_map(op);
  if (!mpi_op) {
    PyErr_SetString(PyExc_RuntimeError,
		    "mpiext.c (allreduce_array): could not determine mpi_op");
    return NULL;
  }

  if (op == MAXLOC || op == MINLOC) {
    PyErr_SetString(PyExc_RuntimeError,
		    "mpiext.c (allreduce_array): MAXLOC and MINLOC are not implemented");
    return NULL;
  }
  else {
    /* call the MPI routine */
    error =  MPI_Allreduce(x->data, d->data, count, mpi_type, mpi_op, \
			MPI_COMM_WORLD);
  }

  if (error != 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    sprintf(errmsg, "Proc %d: MPI_Allreduce failed with error code %d\n",
	    myid, error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);
    return NULL;
  }

  Py_INCREF(Py_None);
  return (Py_None);
}


/*************************************************************/
/* do multiple isends and irecv of Numpy array buffers        */
/* of type float, double, int, or long                       */
/*                                                           */
/*************************************************************/
static PyObject *sendrecv_array(PyObject *self, PyObject *args) {
  PyObject *send_bufs;
  PyArrayObject *send_dest;
  PyObject *recv_bufs;
  PyArrayObject *recv_dest;

  PyObject *seq;

  PyArrayObject *x;
  //int op, error, count, count1, myid;
  //MPI_Datatype mpi_type, buffer_type;
  //MPI_Op mpi_op;
  int i, len;

  double *xdata;

  /* process the parameters */
  if (!PyArg_ParseTuple(args, "OOOO", &send_bufs, &send_dest, &recv_bufs, &recv_dest)) {
    PyErr_SetString(PyExc_RuntimeError,
		    "mpiextras.c (sendrecv_array): could not parse input");
    return NULL;
  }

  seq = PySequence_Fast(send_bufs, "expected a sequence");
  len = PySequence_Size(send_bufs);
  printf("len of sequence %d\n",len);
  for (i = 0; i < len; i++) {
    x = (PyArrayObject *) PySequence_Fast_GET_ITEM(seq, i);
    //lenx = x->dimensions[0];
    printf("buf size %d\n",len);
    xdata = (double *) x->data;
    printf("x->data[0] %g\n",xdata[0]);
  }
  Py_DECREF(seq);


  Py_INCREF(Py_None);
  return (Py_None);
}



/*************************************************************/
/* do multiple isends and irecv of Numpy array buffers        */
/* of type float, double, int, or long                       */
/*                                                           */
/*************************************************************/
static PyObject *send_recv_via_dicts(PyObject *self, PyObject *args) {

  PyObject *send_dict;
  PyObject *recv_dict;

  //PyObject *seq;

  //PyObject *list;
  PyArrayObject *X;
  //PyArrayObject *Id;
  //int op, error, count, count1, myid;
  //MPI_Datatype mpi_type, buffer_type;
  //MPI_Op mpi_op;
  int k, lenx;
  int num_recv=0;
  int num_send=0;

  int ierr;
  MPI_Request requests[40];
  MPI_Status statuses[40];

  Py_ssize_t pos = 0;

  PyObject *key, *value;
  //long *Iddata;
  //double *xdata;

  /* process the parameters */
  if (!PyArg_ParseTuple(args, "OO", &send_dict, &recv_dict)) {
    PyErr_SetString(PyExc_RuntimeError,
		    "mpiextras.c (sendrecv_array): could not parse input");
    return NULL;
  }

  //----------------------------------------------------------------------------
  // Do the recv first
  //----------------------------------------------------------------------------
  num_recv = PyDict_Size(recv_dict);
  //printf("num_recv = %d\n",num_recv);
  if (num_recv>20) {
    PyErr_SetString(PyExc_RuntimeError,
		    "mpiextras.c; Number of recv communication buffers > 10");
    return NULL;
  }

  pos = 0;
  k = 0;
  while (PyDict_Next(recv_dict, &pos, &key, &value)) {
    int i = PyInt_AS_LONG(key);
    //printf("key %d\n",i);

    //Id = (PyArrayObject *) PyList_GetItem(value, 0);
    X   = (PyArrayObject *) PyList_GetItem(value, 2);

    lenx = X->dimensions[0]*X->dimensions[1];


    //printf("buf size %d by 3\n",lenx/3);
    //xdata = (double *) X->data;
    //Iddata = (long *) Id->data;

    //printf("k = %d \n",k);
    ierr = MPI_Irecv(X->data, lenx, MPI_DOUBLE, i, 123, MPI_COMM_WORLD, &requests[k]);
    if (ierr>0) {
        PyErr_SetString(PyExc_RuntimeError,
    		    "mpiextras.c; error from MPI_Irecv");
        return NULL;
      }
    //printf("ierr = %d \n",ierr);
    k++;
    //for (k=0 ; k < lenx ; k++) printf("X[%d] = %g Id[%d] = %ld \n",k,xdata[k],k,Iddata[k]);
}
  //----------------------------------------------------------------------------
  // Do the sends second
  //----------------------------------------------------------------------------
  num_send = PyDict_Size(send_dict);
  //printf("num_send = %d\n",num_send);
  if (num_send>20) {
    PyErr_SetString(PyExc_RuntimeError,
		    "mpiextras.c; Number of send communication buffers > 10");
    return NULL;
  }

  pos = 0;
  while (PyDict_Next(send_dict, &pos, &key, &value)) {
    int i = PyInt_AS_LONG(key);
    //printf("key %d\n",i);

    //Id = (PyArrayObject *) PyList_GetItem(value, 0);
    X   = (PyArrayObject *) PyList_GetItem(value, 2);

    lenx = X->dimensions[0]*X->dimensions[1];

    //printf("buf size %d by 3 \n",lenx/3);
    //xdata = (double *) X->data;
    //Iddata = (long *) Id->data;
    //for (k=0 ; k < lenx ; k++) printf("X[%d] = %g Id[%d] = %ld \n",k,xdata[k],k,Iddata[k]);

    //printf("k = %d \n",k);
    ierr = MPI_Isend(X->data, lenx, MPI_DOUBLE, i, 123, MPI_COMM_WORLD, &requests[k]);
    //printf("ierr = %d \n",ierr);
    
    k++;
}



  //printf("k = %d\n",k);


  //----------------------------------------------------------------------------
  // Now complete communication. We could put some computation between the 
  // communication calls above and this call.
  //----------------------------------------------------------------------------
  ierr =  MPI_Waitall(k,requests,statuses);





/*
  for (i = 0; i < len; i++) {
    x = (PyArrayObject *) PySequence_Fast_GET_ITEM(seq, i);
    lenx = x->dimensions[0];
    printf("buf size %d\n",len);
    xdata = (double *) x->data;
    printf("x->data[0] %g\n",xdata[0]);y
  }
  Py_DECREF(seq);
*/


/*
//   Input check and determination of MPI type
  mpi_type = type_map(x, &count);
  if (!mpi_type) {
    PyErr_SetString(PyExc_RuntimeError,
		    "mpiext.c (allreduce_array): could not determine mpi_type");
    return NULL;
  }


  //This error is caught at the pypar level - so we won't end up here
  //  unless mpiext is being used independently
  buffer_type = type_map(d, &count1);
  if (mpi_type != buffer_type) {
    sprintf(errmsg, "mpiext.c (allreduce_array): Input array and buffer must be of the same type.");
    PyErr_SetString(PyExc_RuntimeError, errmsg);

    return NULL;
  }

  if (count != count1) {
    PyErr_SetString(PyExc_RuntimeError,
		    "mpiext.c (allreduce_array): Input array and buffer must have same length");
    return NULL;
  }

  // Input check and determination of MPI op
  mpi_op = op_map(op);
  if (!mpi_op) {
    PyErr_SetString(PyExc_RuntimeError,
		    "mpiext.c (allreduce_array): could not determine mpi_op");
    return NULL;
  }

  if (op == MAXLOC || op == MINLOC) {
    PyErr_SetString(PyExc_RuntimeError,
		    "mpiext.c (allreduce_array): MAXLOC and MINLOC are not implemented");
    return NULL;
  }
  else {
    // call the MPI routine
    error =  MPI_Allreduce(x->data, d->data, count, mpi_type, mpi_op, \
			MPI_COMM_WORLD);
  }

  if (error != 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    sprintf(errmsg, "Proc %d: MPI_Allreduce failed with error code %d\n",
	    myid, error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);
    return NULL;
  }
*/


  Py_INCREF(Py_None);
  return (Py_None);
}


 
/**********************************/
/* Method table for python module */
/**********************************/
static struct PyMethodDef MethodTable[] = {
  {"isend_array", isend_array, METH_VARARGS},
  {"ireceive_array", ireceive_array, METH_VARARGS},
  {"allreduce_array", allreduce_array, METH_VARARGS},
  {"sendrecv_array", sendrecv_array, METH_VARARGS},
  {"send_recv_via_dicts", send_recv_via_dicts, METH_VARARGS},
  {NULL, NULL}
};


/***************************/
/* Module initialisation   */
/***************************/
void initmpiextras(void){
  PyObject *m, *ModDict;
  
  m = Py_InitModule("mpiextras", MethodTable);
  
  /* to handle MPI symbolic constants*/
  ModDict = PyModule_GetDict(m); 
  SetDictInt("MPI_ANY_TAG", MPI_ANY_TAG);
  SetDictInt("MPI_TAG_UB", MPI_TAG_UB);  
  SetDictInt("MPI_ANY_SOURCE", MPI_ANY_SOURCE);
  SetDictInt("MAX", MAX);
  SetDictInt("MIN", MIN);
  SetDictInt("SUM", SUM);
  SetDictInt("PROD", PROD);
  SetDictInt("LAND", LAND);
  SetDictInt("BAND", BAND);
  SetDictInt("LOR", LOR);
  SetDictInt("BOR", BOR);
  SetDictInt("LXOR", LXOR);
  SetDictInt("BXOR", BXOR);

   
  import_array();     /* Necessary for handling of numpy structures  */
}

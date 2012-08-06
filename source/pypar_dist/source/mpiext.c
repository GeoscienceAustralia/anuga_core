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


static char errmsg[132];  /*Used to cretae exception messages*/

/* MPI_Bsend() related variables. */
static void *pt_buf;	/* Pointer to allocated buffer. */
static int buf_size;	/* Size of buffer to allocate. */ 

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


/*********************************************************/
/* send_string                                           */
/* Send string of characters                             */
/*                                                       */
/*********************************************************/
static PyObject *send_string(PyObject *self, PyObject *args) {
  char *s;
  int destination, tag, length, error, myid;
 
  /* process the parameters */
  if (!PyArg_ParseTuple(args, "s#ii", &s, &length, &destination, &tag))
    return NULL;
  
  /* call the MPI routine */
  error = MPI_Send(s, length, MPI_CHAR, destination, tag, MPI_COMM_WORLD);
  
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

/**********************************************************/
/* receive_string                                         */
/* Receive string of characters                           */
/*                                                        */
/**********************************************************/
static PyObject *receive_string(PyObject *self, PyObject *args) {
  char *s;
  int source, tag, length, error, st_length, myid; 
  MPI_Status status;

  /* process the parameters */
  if (!PyArg_ParseTuple(args, "s#ii", &s, &length, &source, &tag))
    return NULL;
    
  /* call the MPI routine */
  error = MPI_Recv(s, length, MPI_CHAR, source, tag, MPI_COMM_WORLD, &status);

  if (error != 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);    
    sprintf(errmsg, "Proc %d: MPI_Recv failed with error code %d\n", 
	    myid, error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);   
    return NULL;
  }  
     
  MPI_Get_count(&status, MPI_CHAR, &st_length); 
  /* status.st_length is not available in all MPI implementations
  // Alternative is: MPI_Get_elements(MPI_Status *, MPI_Datatype, int *); */

  /* Still include error msg in status in case exception is caught. */
  return Py_BuildValue("(iiiii)", status.MPI_SOURCE, status.MPI_TAG,
  status.MPI_ERROR, st_length, sizeof(char));  
}

/**********************************************************/
/* broadcast_string                                           */
/* Broadcast string of characters                         */
/*                                                        */
/**********************************************************/
static PyObject *broadcast_string(PyObject *self, PyObject *args) {
  char *s;
  int source, length, error, myid; 

  /* process the parameters */
  if (!PyArg_ParseTuple(args, "s#i", &s, &length, &source))
    return NULL;
    
  /* call the MPI routine */
  error = MPI_Bcast(s, length, MPI_CHAR, source, MPI_COMM_WORLD);
  if (error != 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);    
    sprintf(errmsg, "Proc %d: MPI_Bcast failed with error code %d\n", 
	    myid, error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);   
    return NULL;
  }  
   
  Py_INCREF(Py_None);
  return (Py_None);
}

/**********************************************************/
/* scatter_string                                         */
/* Scatter string of characters                           */
/*                                                        */
/**********************************************************/
static PyObject *scatter_string(PyObject *self, PyObject *args) {
  char *s;
  char *d;
  int source, count, error, myid, numprocs; 

  /* process the parameters */
  if (!PyArg_ParseTuple(args, "s#si", &s, &count, &d, &source))
    return NULL;
    
  error = MPI_Comm_size(MPI_COMM_WORLD,&numprocs);    
  count = count/numprocs;
  
  /* call the MPI routine */
  error = MPI_Scatter(s, count, MPI_CHAR, d, count, MPI_CHAR, 
		      source, MPI_COMM_WORLD);

  if (error != 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);    
    sprintf(errmsg, "Proc %d: MPI_Scatter failed with error code %d\n", 
	    myid, error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);   
    return NULL;
  }  
   
  Py_INCREF(Py_None);
  return (Py_None);
}

/**********************************************************/
/* gather_string                                         */
/* Gather string of characters                           */
/*                                                        */
/**********************************************************/
static PyObject *gather_string(PyObject *self, PyObject *args) {
  char *s;
  char *d;
  int source, error, count, myid; 

  /* process the parameters */
  if (!PyArg_ParseTuple(args, "s#si", &s, &count, &d, &source))
    return NULL;
    
  /* call the MPI routine */
  error = MPI_Gather(s, count, MPI_CHAR, d, count,  MPI_CHAR, 
		     source, MPI_COMM_WORLD);
  
  if (error != 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);    
    sprintf(errmsg, "Proc %d: MPI_Gather failed with error code %d\n", 
	    myid, error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);   
    return NULL;
  }  
   
  Py_INCREF(Py_None);
  return (Py_None);
}


/**********************************************************/
/* bsend_string                                           */
/* Send a string of characters using MPI_Bsend().         */
/*                                                        */
/* Return value: PyNone.                                  */
/**********************************************************/
static PyObject *bsend_string(PyObject *self, PyObject *args) {
  char *s;
  int destination, tag, length;
  int error, myid;


  /* Process the parameters. */
  if (!PyArg_ParseTuple(args, "s#ii", &s, &length, &destination, &tag))
    return NULL;

  /* Call the MPI routine. */
  error = MPI_Bsend(s, length, MPI_CHAR, destination, tag, MPI_COMM_WORLD);


  if (error != 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);    
    sprintf(errmsg, "Proc %d: MPI_Bsend failed with error code %d\n", 
	    myid, error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);
    return NULL;
  }  
      
  Py_INCREF(Py_None);
  return (Py_None);
}


/**********************************************************/
/* bsend_array                                            */
/* Send a Numpy array using MPI_Bsend().                  */
/* Accepted types for array: float, double, int, or long. */
/* Return value: PyNone.                                  */
/**********************************************************/
static PyObject *bsend_array(PyObject *self, PyObject *args) { 
  PyObject *input;
  PyArrayObject *x;

  int destination;
  int tag;
  int count;
  MPI_Datatype mpi_type;

  int error, myid;

  /* Process the parameters. */
  if (!PyArg_ParseTuple(args, "Oii", &input, &destination, &tag))
    return NULL;

  /* Make Numpy array from general sequence type (no cost if already Numpy). */
  x = (PyArrayObject *)
    PyArray_ContiguousFromObject(input, NPY_NOTYPE, 0, 0);

  /* Input check and determination of MPI type */          
  mpi_type = type_map(x, &count);
  if (!mpi_type)
    return NULL;

  /* Call the MPI routine */
  error = MPI_Bsend(x->data, count, mpi_type, destination, tag,
          MPI_COMM_WORLD);

  Py_DECREF(x); 	   
  
  if (error != 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);    
    sprintf(errmsg, "Proc %d: MPI_Bsend failed with error code %d\n", 
	    myid, error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);   
    return NULL;
  }  

  Py_INCREF(Py_None);
  return (Py_None);
}




/**********************************************************/
/* send_array                                             */
/* Send Numeric array of type float, double, int, or long */
/*                                                        */
/**********************************************************/
static PyObject *send_array(PyObject *self, PyObject *args) {
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
/* receive_array                                             */
/* Receive Numeric array of type float, double, int, or long */
/*                                                           */
/*************************************************************/
static PyObject *receive_array(PyObject *self, PyObject *args) {
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
/* broadcast_array                                               */
/* Broadcast Num.  array of type float, double, int, or long */
/*                                                           */
/*************************************************************/
static PyObject *broadcast_array(PyObject *self, PyObject *args) {
  PyArrayObject *x;
  int source, error, count, myid;
  MPI_Datatype mpi_type;

  /* process the parameters */
  if (!PyArg_ParseTuple(args, "Oi", &x, &source))
    return NULL;

  /* Input check and determination of MPI type */          
  mpi_type = type_map(x, &count);
  if (!mpi_type) return NULL;  
      
  /* call the MPI routine */
  error =  MPI_Bcast(x->data, count, mpi_type, source, \
		     MPI_COMM_WORLD);
	 
	 
  if (error != 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);    
    sprintf(errmsg, "Proc %d: MPI_Bcast failed with error code %d\n", 
	    myid, error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);   
    return NULL;
  }  
	 
  Py_INCREF(Py_None);
  return (Py_None);
}

/*************************************************************/
/* scatter_array                                             */
/* Scatter Num.    array of type float, double, int, or long */
/*                                                           */
/*************************************************************/
static PyObject *scatter_array(PyObject *self, PyObject *args) {
  PyArrayObject *x;
  PyArrayObject *d;
  int source, error, count, myid, numprocs;
  MPI_Datatype mpi_type;

  /* process the parameters */
  if (!PyArg_ParseTuple(args, "OOi", &x, &d, &source))
    return NULL;

  /* Input check and determination of MPI type */          
  mpi_type = type_map(x, &count);
  if (!mpi_type) return NULL;  
   
  error = MPI_Comm_size(MPI_COMM_WORLD,&numprocs);    
  count = count/numprocs;  
  
  /* call the MPI routine */
  error = MPI_Scatter(x->data, count, mpi_type, d->data, count, 
		      mpi_type, source,	MPI_COMM_WORLD);
	 
  if (error != 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);    
    sprintf(errmsg, "Proc %d: MPI_Scatter failed with error code %d\n", 
	    myid, error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);
    return NULL;
  }  
      
  Py_INCREF(Py_None);
  return (Py_None);
}


/*************************************************************/
/* gather_array                                              */
/* Gather Num.     array of type float, double, int, or long */
/*                                                           */
/*************************************************************/
static PyObject *gather_array(PyObject *self, PyObject *args) {
  PyArrayObject *x;
  PyArrayObject *d;
  int source, error, count, myid;
  MPI_Datatype mpi_type;

  /* process the parameters */
  if (!PyArg_ParseTuple(args, "OOi", &x, &d, &source))
    return NULL;

  /* Input check and determination of MPI type */          
  mpi_type = type_map(x, &count);
  if (!mpi_type) return NULL;  
      
  /* call the MPI routine */
  error =  MPI_Gather(x->data, count, mpi_type, d->data, count, 
		      mpi_type, source,	MPI_COMM_WORLD);

  if (error != 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);    
    sprintf(errmsg, "Proc %d: MPI_Gather failed with error code %d\n", 
	    myid, error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);
    return NULL;
  }  
  
  Py_INCREF(Py_None);
  return (Py_None);
}



/*************************************************************/
/* reduce_array                                              */
/* Reduce Num.     array of type float, double, int, or long */
/*                                                           */
/*************************************************************/
static PyObject *reduce_array(PyObject *self, PyObject *args) {
  PyArrayObject *x;
  PyArrayObject *d;
  int source, op, error, count, count1, myid;
  MPI_Datatype mpi_type, buffer_type;
  MPI_Op mpi_op;

  /* process the parameters */
  if (!PyArg_ParseTuple(args, "OOii", &x, &d, &op, &source)) {
    PyErr_SetString(PyExc_RuntimeError, 
		    "mpiext.c (reduce_array): could not parse input");
    return NULL;
  }
   
  /* Input check and determination of MPI type */          
  mpi_type = type_map(x, &count);
  if (!mpi_type) {
    PyErr_SetString(PyExc_RuntimeError, 
		    "mpiext.c (reduce_array): could not determine mpi_type");
    return NULL;  
  }

  /* This error is caught at the pypar level - so we won't end up here
     unless mpiext is being used independently */
  buffer_type = type_map(d, &count1);
  if (mpi_type != buffer_type) {
    sprintf(errmsg, "mpiext.c (reduce_array): Input array and buffer must be of the same type.");
    PyErr_SetString(PyExc_RuntimeError, errmsg);
		    
    return NULL;  
  }

  if (count != count1) {
    PyErr_SetString(PyExc_RuntimeError, 
		    "mpiext.c (reduce_array): Input array and buffer must have same length");
    return NULL;  
  }
    
  /* Input check and determination of MPI op */ 
  mpi_op = op_map(op);
  if (!mpi_op) {
    PyErr_SetString(PyExc_RuntimeError, 
		    "mpiext.c (reduce_array): could not determine mpi_op");
    return NULL;  
  }
   
  if (op == MAXLOC || op == MINLOC) {
    PyErr_SetString(PyExc_RuntimeError, 
		    "mpiext.c (reduce_array): MAXLOC and MINLOC are not implemented");
    return NULL;  
  }
  else {
    /* call the MPI routine */
    error =  MPI_Reduce(x->data, d->data, count, mpi_type, mpi_op, source, \
			MPI_COMM_WORLD);
  }
         
  if (error != 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);    
    sprintf(errmsg, "Proc %d: MPI_Reduce failed with error code %d\n", 
	    myid, error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);
    return NULL;
  }  
  
  Py_INCREF(Py_None);
  return (Py_None);
}


/*************************************************************/
/* allreduce_array                                              */
/* Allreduce Num.     array of type float, double, int, or long */
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







/**********************************************************/
/* CONTRIBUTED FUNCTIONS FROM CMAKASSIKIS:                */
/* push_for_alloc functions                               */ 
/* mpi_alloc_and_attach                                   */
/* mpi_detach_and_dealloc                                 */ 
/* mpi_attach                                             */
/* mpi_detach                                             */
/* mpi_alloc                                              */
/* mpi_dealloc                                            */ 
/**********************************************************/


/*
 * 'push_for_alloc' functions are used to increase the current size in bytes of
 * the buffer. No effective allocation is done at the call of these functions.
 * They merely update 'buf_size' whose value will be used when calling the
 * function for allocation.
 */

/*
 * 'string_push_for_alloc_and_attach'
 *
 * This function is used when protocol is set to 'string' or 'vanilla'.
 *
 * Return value: 'buf_size'.
 */

static PyObject *string_push_for_alloc_and_attach(PyObject *self, 
						  PyObject *args) {
  char *s;
  int length;

  /* Process the parameters. */
  if (!PyArg_ParseTuple(args, "s#", &s, &length))
    return NULL;

  buf_size += (length + MPI_BSEND_OVERHEAD);

  return Py_BuildValue("i", buf_size);
}

/*
 * 'array_push_for_alloc_and_attach' 
 *
 * This function is used when protocol is set to 'array'.
 *
 * Return value: 'buf_size'.
 */

static PyObject *array_push_for_alloc_and_attach(PyObject *self, 
						 PyObject *args) {
  PyArrayObject *array;	
  int count;				
  int nbytes;				
  MPI_Datatype mpi_type;	

  int error, myid;

  count = nbytes = error = 0;
  mpi_type = (MPI_Datatype) 0;
  myid = -1;

  /* Process the parameters. */
  if (!PyArg_ParseTuple(args, "O", &array))
    return NULL;

  /* Input check and determination of MPI type */          
  mpi_type = type_map(array, &count);
  if (!mpi_type) 
    return NULL;

  /* Compute number of bytes. */
  error = MPI_Type_size(mpi_type, &nbytes);
  buf_size += (nbytes * count + MPI_BSEND_OVERHEAD);

  if (error != 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);    
    sprintf(errmsg, "Proc %d: array_push_for_alloc_and_attach: \
	        MPI_Type_size failed with error code %d\n", myid, error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);
    return NULL;
  }  

  return Py_BuildValue("i", buf_size);
}




/*
 * 'mpi_alloc_and_attach'
 *
 * Allocates and attaches a buffer of size 'buf_size' (in bytes) for use with
 * MPI_Bsend() through 'bsend_string' and 'bsend_array'.
 *
 * Return value: 'PyNone'.
 *
 */

static PyObject *mpi_alloc_and_attach(PyObject *self, PyObject *args) {
	
  //int count = 0;
  int error, myid;

  /* Process input parameters. */
  //if (!PyArg_ParseTuple(args, "i", &count))
  //    return NULL;

  /* Allocate MPI_Bsend()'s buffer. */
  pt_buf = (void *) malloc (buf_size);

  if (pt_buf == NULL) {
    PyErr_SetString(PyExc_RuntimeError, 
		    "mpi_alloc_and_attach: Not enough memory to allocate bsend buffer");
    return NULL;
  }

  /* Attach MPI_Bsend()'s buffer. */
  error = MPI_Buffer_attach(pt_buf, buf_size);

  if (error != 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);    
    sprintf(errmsg, "Proc %d: mpi_alloc_and_attach: MPI_Buffer_attach: \
	                 failed with error code %d\n", myid, error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);
    return NULL;
  }  

  Py_INCREF(Py_None);
  return (Py_None);
}

/*
 * 'mpi_detach_and_dealloc'
 *
 * Detaches and deallocates buffer used for making MPI_Bsend().
 *
 * Reinitializes global variables 'pt_buf' and 'buf_size'.
 *
 * Return value: 'PyNone'.
 *
 */
static PyObject *mpi_detach_and_dealloc(PyObject *self, PyObject *args) {
  MPI_Buffer_detach(&pt_buf, &buf_size);

  free(pt_buf);
  pt_buf = NULL;

  buf_size = 0;

  Py_INCREF(Py_None);
  return (Py_None);
}

/*
 * 'mpi_attach'
 *
 * Attaches a buffer in order to make a MPI_Bsend().
 *
 * Return value: PyNone.
 *
 */
static PyObject *mpi_attach(PyObject *self, PyObject *args) {
  int error, myid;

  error = MPI_Buffer_attach(pt_buf, buf_size);

  if (error != 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);    
    sprintf(errmsg, "Proc %d: MPI_Buffer_attach: failed with error code %d\n", 
	    myid, error);
    PyErr_SetString(PyExc_RuntimeError, errmsg); 
    return NULL;
  }  

  Py_INCREF(Py_None);
  return (Py_None);
}

/*
 * 'mpi_detach'
 *
 * Detach a buffer after an MPI_Bsend().
 *
 * Return value: 'PyNone'.
 *
 */

static PyObject *mpi_detach(PyObject *self, PyObject *args) {
	MPI_Buffer_detach(&pt_buf, &buf_size);

	Py_INCREF(Py_None);
	return (Py_None);
}

/*
 * Allocates 'pt_buf'.
 *
 * Default: allocates a buffer of size equal to 'buf_size' (in bytes).
 * mpi_attach() has an optional parameter which, if set, overrides the 
 * size specified by 'buf_size'.
 *
 * Return value: 'buf_size'.
 *
 */
static PyObject *mpi_alloc(PyObject *self, PyObject *args) {
  int nbytes = -1;	/* Number of bytes we wish to allocate. */

  /* Process the parameters. */
  if (!PyArg_ParseTuple(args, "|i", &nbytes))
    return NULL;

  if ((nbytes < 0) && (buf_size <=0)) {
    PyErr_SetString(PyExc_RuntimeError, 
		    "mpi_alloc: Buffer size must be set either through push_for_alloc() or directly via alloc()'s optional parameter.");
    return NULL;
  }

  /* Override 'buf_size' value. */
  if (nbytes > 0) {
    buf_size = nbytes;
  }

  pt_buf = (void *) malloc (buf_size);
  if (pt_buf == NULL) {
    PyErr_SetString(PyExc_RuntimeError, 
		    "mpi_alloc: Not enough memory to allocate mpi bsend buffer");
    return NULL;
  }

  return Py_BuildValue("i", buf_size);
}

/*
 * Deallocate.
 */

static PyObject *mpi_dealloc(PyObject *self, PyObject *args) {
  free(pt_buf);
  pt_buf = NULL;
  
  buf_size = 0;

  Py_INCREF(Py_None);
  return (Py_None);
}




/*********************************************************/
/* MPI calls such as rank, size, finalize, abort         */
/*                                                       */
/*********************************************************/
static PyObject * rank(PyObject *self, PyObject *args) {
  int error, myid;

  error = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  if (error != 0) {
    sprintf(errmsg, "Proc ?: MPI_Comm_rank failed with error code %d\n", 
	    error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);
    return NULL;
  }  
  
  return Py_BuildValue("i", myid);
}

static PyObject * size(PyObject *self, PyObject *args) {
  int error, numprocs, myid; 
  
  error = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  if (error != 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);    
    sprintf(errmsg, "Proc %d: MPI_Comm_size failed with error code %d\n", 
	    myid, error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);
    return NULL;
  }  
  
  return Py_BuildValue("i", numprocs);
}
  
static PyObject * get_processor_name(PyObject *self, PyObject *args) {  
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int  error, namelen, myid;

  error = MPI_Get_processor_name(processor_name,&namelen);
  if (error != 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);    
    sprintf(errmsg, "Proc %d: Get_processor_name failed with error code %d\n", 
	    myid, error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);
    return NULL;
  }  
  
  return Py_BuildValue("s#", processor_name, namelen);
}   

static PyObject * init(PyObject *self, PyObject *args) {  
  PyObject *input;

  int i, error, myid;
  int argc = 0;  
  char **argv;   

  /* process input parameters */
  if (!PyArg_ParseTuple(args, "O!", &PyList_Type, &input))
    return NULL;

  /* Reconstruct C-commandline */     
  argc = PyList_Size(input); /* Number of commandline arguments */
  argv = (char**) malloc((argc+1)*sizeof(char*)); 
  
  for (i=0; i<argc; i++)  
    argv[i] = PyString_AsString( PyList_GetItem(input, i) );
    
  argv[i] = NULL; /* Lam 7.0 requires last arg to be NULL  */
  
  error = MPI_Init(&argc, &argv); 
  if (error != 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);    
    sprintf(errmsg, "Proc ?: MPI_Init failed with error code %d\n", 
	    error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);   
    return NULL;
  }  

  Py_INCREF(Py_None);
  return (Py_None);
} 


static PyObject * initialized(PyObject *self, PyObject *args) {  
  int error, flag, myid;
  
  error = MPI_Initialized(&flag);
  
  if (error != 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);    
    sprintf(errmsg, "Proc %d: MPI_Initialized failed with error code %d\n", 
	    myid, error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);
    return NULL;
  }  

  return Py_BuildValue("i", flag);  
} 

  
static PyObject * finalize(PyObject *self, PyObject *args) {  
  int error, myid;

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);  
  
  error = MPI_Finalize();
  if (error != 0) {
    sprintf(errmsg, "Proc %d: MPI_Finalize failed with error code %d\n", 
	    myid, error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);
    return NULL;
  }  
    
  Py_INCREF(Py_None);
  return (Py_None);
} 

static PyObject * mpi_abort(PyObject *self, PyObject *args) {  
  int error, code=0, myid;
  
  error = MPI_Abort(MPI_COMM_WORLD, code);
  if (error != 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);    
    sprintf(errmsg, "Proc %d: MPI_Abort failed with error code %d\n", 
	    myid, error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);
    return NULL;
  }  
  
  Py_INCREF(Py_None);
  return (Py_None);
} 

static PyObject * barrier(PyObject *self, PyObject *args) {  
  int error, myid;
  
  error = MPI_Barrier(MPI_COMM_WORLD);
  if (error != 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);    
    sprintf(errmsg, "Proc %d: MPI_Barrier failed with error code %d\n", 
	    myid, error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);
    return NULL;
  }  
  
  Py_INCREF(Py_None);
  return (Py_None);
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
  {"barrier", barrier, METH_VARARGS},          
  {"time", Wtime, METH_VARARGS},            
  {"get_processor_name", get_processor_name, METH_VARARGS},              
  {"init", init, METH_VARARGS},          
  {"initialized", initialized, METH_VARARGS},       
  {"finalize", finalize, METH_VARARGS},        
  {"abort", mpi_abort, METH_VARARGS},          
  {"send_string", send_string, METH_VARARGS},
  {"receive_string", receive_string, METH_VARARGS},      
  {"broadcast_string", broadcast_string, METH_VARARGS},        
  {"scatter_string", scatter_string, METH_VARARGS},        
  {"gather_string", gather_string, METH_VARARGS},        
  {"send_array", send_array, METH_VARARGS},
  {"receive_array", receive_array, METH_VARARGS},    
  {"broadcast_array", broadcast_array, METH_VARARGS},              
  {"scatter_array", scatter_array, METH_VARARGS},              
  {"gather_array", gather_array, METH_VARARGS},              
  {"reduce_array", reduce_array, METH_VARARGS},
  {"allreduce_array", allreduce_array, METH_VARARGS},
  /* Functions providing 'MPI_Bsend' support. */
  {"bsend_string", bsend_string, METH_VARARGS},
  {"bsend_array", bsend_array, METH_VARARGS},
  {"string_push_for_alloc_and_attach", string_push_for_alloc_and_attach, 
   METH_VARARGS},
  {"array_push_for_alloc_and_attach", array_push_for_alloc_and_attach, 
   METH_VARARGS},
  {"mpi_alloc_and_attach", mpi_alloc_and_attach, METH_VARARGS},
  {"mpi_detach_and_dealloc", mpi_detach_and_dealloc, METH_VARARGS},
  {"mpi_alloc", mpi_alloc, METH_VARARGS},
  {"mpi_dealloc", mpi_dealloc, METH_VARARGS},
  {"mpi_attach", mpi_attach, METH_VARARGS},
  {"mpi_detach", mpi_detach, METH_VARARGS},
  {NULL, NULL}
};


/***************************/
/* Module initialisation   */
/***************************/
void initmpiext(void){
  PyObject *m, *ModDict;
  
  m = Py_InitModule("mpiext", MethodTable);
  
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
  /*SetDictInt("MAXLOC", MAXLOC);*/
  /*SetDictInt("MINLOC", MINLOC);*/
  /*SetDictInt("REPLACE", REPLACE);*/

  /*SetDictInt("MPI_COMM_WORLD", MPI_COMM_WORLD);  */
   
  import_array();     /* Necessary for handling of numpy structures  */
}

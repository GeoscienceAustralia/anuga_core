/************************************************************************/
/* PyPAR - Parallel Python using MPI                 	                */
/* Copyright (C) 2001, 2002, 2003 Ole M. Nielsen                        */
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
/* Contact address: Ole.Nielsen@anu.edu.au                              */
/*                                                                 	*/
/* version (see __version__ in pypar.py)                                */
/* date (see __date__ in pypar.py)                                      */
/************************************************************************/


#include "Python.h"
#include "mpi.h"
#include "math.h"
#include "Numeric/arrayobject.h"

/* to handle MPI constants export (shamelessly stolen from _cursesmodule.c)*/
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

int length(PyArrayObject *x) {  
  /*Compute the total length of contiguous array
  /*
  /*Necessary for communicating multi dimensional arrays */
  
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
  /* the Python data type as follows  
  /*
  /* TYPE    py_type  mpi_type  bytes  symbol
  /* ---------------------------------------- 
  /* INT       4        6         4      'i'
  /* LONG      5        8         8      'l'
  /* FLOAT     6       10         4      'f'  
  /* DOUBLE    7       11         8      'd'
  /*
  /* Also return the total number of elements in the array
  /*
  /* The Python datatype COMPLEX ('F') and COMPLEX_DOUBLE ('D')
  /* is treated as a special case to the absence of an 
  /* MPI_COMPLEX datatype:
  /*
  /* Complex arrays are mapped to float or double arrays with real 
  /* and imaginary parts alternating and count is updated. */
  
  int py_type;
  MPI_Datatype mpi_type;

  *count = length(x);
    
  py_type = x -> descr -> type_num;     
  if (py_type == PyArray_DOUBLE) 
    mpi_type = MPI_DOUBLE;
  else if (py_type == PyArray_INT) 
    mpi_type = MPI_INT;
  else if (py_type == PyArray_CDOUBLE) {
    mpi_type = MPI_DOUBLE;
    (*count) *= 2;
  } else if (py_type == PyArray_FLOAT) 
    mpi_type = MPI_FLOAT;
  else if (py_type == PyArray_LONG)   
    mpi_type = MPI_LONG;  
  else if (py_type == PyArray_CFLOAT) {
    mpi_type = MPI_FLOAT;
    (*count) *= 2;
  } else {
    PyErr_SetString(PyExc_ValueError, "Array must be of type int or float");
    return NULL;
  }      

  /*printf("Types %d %d\n", py_type, mpi_type); */
  
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
    return NULL;
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
    sprintf(errmsg, "Proc %d: MPI_Send failed with error code %d\n", myid, error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);    /*raise ValueError, errmsg */
    return NULL;
  }  
      
  return Py_BuildValue("");
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
    sprintf(errmsg, "Proc %d: MPI_Recv failed with error code %d\n", myid, error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);   
    return NULL;
  }  
     
  MPI_Get_count(&status, MPI_CHAR, &st_length); 
  /*status.st_length is not available in all MPI implementations
  //Alternative is: MPI_Get_elements(MPI_Status *, MPI_Datatype, int *); */

  /*Still include error msg in status in case exception is caught. */
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
    sprintf(errmsg, "Proc %d: MPI_Bcast failed with error code %d\n", myid, error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);   
    return NULL;
  }  
   
  return Py_BuildValue("");  
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
  error = MPI_Scatter(s, count, MPI_CHAR, d, count, MPI_CHAR, source, MPI_COMM_WORLD);

  if (error != 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);    
    sprintf(errmsg, "Proc %d: MPI_Scatter failed with error code %d\n", myid, error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);   
    return NULL;
  }  
   
  return Py_BuildValue("");  
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
  error = MPI_Gather(s, count, MPI_CHAR, d, count,  MPI_CHAR, source, MPI_COMM_WORLD);
  
  if (error != 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);    
    sprintf(errmsg, "Proc %d: MPI_Gather failed with error code %d\n", myid, error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);   
    return NULL;
  }  
   
  return Py_BuildValue("");  
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
    
  /*if (!PyArg_ParseTuple(args, "Oii", &x, &destination, &tag))*/
  /* return NULL;     */
    
  /* Make Numeric array from general sequence type (no cost if already Numeric)*/    
  x = (PyArrayObject *)
    PyArray_ContiguousFromObject(input, PyArray_NOTYPE, 0, 0);
    
  /* Input check and determination of MPI type */          
  mpi_type = type_map(x, &count);
  if (!mpi_type) return NULL;
    
  /* call the MPI routine */
  error = MPI_Send(x->data, count, mpi_type, destination, tag,\
		   MPI_COMM_WORLD);
  Py_DECREF(x); 	   
  
  if (error != 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);    
    sprintf(errmsg, "Proc %d: MPI_Send failed with error code %d\n", myid, error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);   
    return NULL;
  }  
   
  return Py_BuildValue("");  
}

/*************************************************************/
/* receive_array                                             */
/* Receive Numeric array of type float, double, int, or long */
/*                                                           */
/*************************************************************/
static PyObject *receive_array(PyObject *self, PyObject *args) {
  PyObject *input;
  PyArrayObject *x;
  int source, tag, error, st_length, size, count, myid;
  MPI_Datatype mpi_type;
  MPI_Status status;

  /* process the parameters */
  /*if (!PyArg_ParseTuple(args, "Oii", &input, &source, &tag))*/
  /*  return NULL;*/
    
  if (!PyArg_ParseTuple(args, "Oii", &x, &source, &tag))
    return NULL;    
    
  /* Make Numeric array from general sequence type (no cost if already Numeric)*/    
  /*x = (PyArrayObject *) */
  /*  PyArray_ContiguousFromObject(input, PyArray_NOTYPE, 0, 0);*/
    
  /* Input check and determination of MPI type */          
  mpi_type = type_map(x, &count);
  if (!mpi_type) return NULL;  
      
  /* call the MPI routine */
  error =  MPI_Recv(x->data, count, mpi_type, source, tag, \
		    MPI_COMM_WORLD, &status);
	 
  /* Do not DECREF x as it must be returned to Python*/
  	 
  if (error != 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);    
    sprintf(errmsg, "Proc %d: MPI_Recv failed with error code %d\n", myid, error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);   
    return NULL;
  }  
   
	 
  MPI_Get_count(&status, mpi_type, &st_length); 
  /* status.st_length is not available in all MPI implementations*/
  /*Alternative is: MPI_Get_elements(MPI_Status *, MPI_Datatype, int *);*/
	 
      
  /*FIXME: This might not be watertight on all platforms */
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
  /*MPI_Status status;*/

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
    sprintf(errmsg, "Proc %d: MPI_Bcast failed with error code %d\n", myid, error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);   
    return NULL;
  }  
	 
  return Py_BuildValue("");
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
  /*MPI_Status status; */

  /* process the parameters */
  if (!PyArg_ParseTuple(args, "OOi", &x, &d, &source))
    return NULL;

  /* Input check and determination of MPI type */          
  mpi_type = type_map(x, &count);
  if (!mpi_type) return NULL;  
   
  error = MPI_Comm_size(MPI_COMM_WORLD,&numprocs);    
  count = count/numprocs;  
  
  /* call the MPI routine */
  error = MPI_Scatter(x->data, count, mpi_type, d->data, count, mpi_type, source, \
		      MPI_COMM_WORLD);
	 
  if (error != 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);    
    sprintf(errmsg, "Proc %d: MPI_Scatter failed with error code %d\n", myid, error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);    /*raise ValueError, errmsg*/
    return NULL;
  }  
      
  return Py_BuildValue("");
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
  /*MPI_Status status; */

  /* process the parameters */
  if (!PyArg_ParseTuple(args, "OOi", &x, &d, &source))
    return NULL;

  /* Input check and determination of MPI type */          
  mpi_type = type_map(x, &count);
  if (!mpi_type) return NULL;  
      
  /* call the MPI routine */
  error =  MPI_Gather(x->data, count, mpi_type, d->data, count, mpi_type, source, \
		      MPI_COMM_WORLD);

  if (error != 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);    
    sprintf(errmsg, "Proc %d: MPI_Gather failed with error code %d\n", myid, error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);    /*raise ValueError, errmsg */
    return NULL;
  }  
  
  return Py_BuildValue("");
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
  MPI_Datatype mpi_type;
  /*MPI_Status status;*/
  MPI_Op mpi_op;

  /* process the parameters */
  if (!PyArg_ParseTuple(args, "OOii", &x, &d, &op, &source))
    return NULL;
   
  /* Input check and determination of MPI type */          
  mpi_type = type_map(x, &count);
  if (!mpi_type) return NULL;  
  if (mpi_type != type_map(d, &count1)) {
    printf ("Input array and buffer must be of the same type\n");
    return Py_BuildValue("i", -666);    
  }

  if (count != count1) {
    printf ("Input array and buffer must have same length\n");
    return Py_BuildValue("i", -666);    
  }
    
  /* Input check and determination of MPI op */ 
  /*printf("op: %d\n", op);         */
  mpi_op = op_map(op);
  if (!mpi_op) return NULL;  
   
  if (op == MAXLOC || op == MINLOC) {
    /*not implemented*/
    return Py_BuildValue("i", -666);
  }
  else {
  /* call the MPI routine */
  error =  MPI_Reduce(x->data, d->data, count, mpi_type, mpi_op, source, \
         MPI_COMM_WORLD);
  }
         
  if (error != 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);    
    sprintf(errmsg, "Proc %d: MPI_Reduce failed with error code %d\n", myid, error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);    /*raise ValueError, errmsg*/
    return NULL;
  }  
  
  return Py_BuildValue("");
}

/*********************************************************/
/* MPI calls rank, size, finalize, abort                 */
/*                                                       */
/*********************************************************/

static PyObject * rank(PyObject *self, PyObject *args) {
  int error, myid;

  error = MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  if (error != 0) {
    sprintf(errmsg, "Proc ?: MPI_Comm_rank failed with error code %d\n", error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);    /*raise ValueError, errmsg*/
    return NULL;
  }  
  
  return Py_BuildValue("i", myid);
}

static PyObject * size(PyObject *self, PyObject *args) {
  int error, numprocs, myid; 
  
  error = MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  if (error != 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);    
    sprintf(errmsg, "Proc %d: MPI_Comm_size failed with error code %d\n", myid, error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);    /*raise ValueError, errmsg*/
    return NULL;
  }  
  
  return Py_BuildValue("i", numprocs);
}
  
static PyObject * Get_processor_name(PyObject *self, PyObject *args) {  
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int  error, namelen, myid;

  error = MPI_Get_processor_name(processor_name,&namelen);
  if (error != 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);    
    sprintf(errmsg, "Proc %d: Get_processor_name failed with error code %d\n", myid, error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);    /*raise ValueError, errmsg*/
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
  /*                           */ 
  argc = PyList_Size(input); /*Number of commandline arguments*/
  argv = (char**) malloc((argc+1)*sizeof(char*)); 
  
  for (i=0; i<argc; i++)  
    argv[i] = PyString_AsString( PyList_GetItem(input, i) );
    
  argv[i] = NULL; /*Lam 7.0 requires last arg to be NULL  */
  
  error = MPI_Init(&argc, &argv); 
  if (error != 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);    
    sprintf(errmsg, "Proc ?: MPI_Init failed with error code %d\n", error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);   
    return NULL;
  }  

  return Py_BuildValue("");  
} 


static PyObject * initialized(PyObject *self, PyObject *args) {  
  int error, flag, myid;
  
  error = MPI_Initialized(&flag);
  
  if (error != 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);    
    sprintf(errmsg, "Proc %d: MPI_Initialized failed with error code %d\n", myid, error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);    /*raise ValueError, errmsg*/
    return NULL;
  }  

  return Py_BuildValue("i", flag);  
} 

  
static PyObject * finalize(PyObject *self, PyObject *args) {  
  int error, myid;

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);  
  
  error = MPI_Finalize();
  if (error != 0) {
    sprintf(errmsg, "Proc %d: MPI_Finalize failed with error code %d\n", myid, error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);    /*raise ValueError, errmsg*/
    return NULL;
  }  
    
  return Py_BuildValue("");
} 

static PyObject * mpi_abort(PyObject *self, PyObject *args) {  
  int error, code=0, myid;
  
  error = MPI_Abort(MPI_COMM_WORLD, code);
  if (error != 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);    
    sprintf(errmsg, "Proc %d: MPI_Abort failed with error code %d\n", myid, error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);    /*raise ValueError, errmsg*/
    return NULL;
  }  
  
  return Py_BuildValue("");  
} 

static PyObject * barrier(PyObject *self, PyObject *args) {  
  int error, myid;
  
  error = MPI_Barrier(MPI_COMM_WORLD);
  if (error != 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);    
    sprintf(errmsg, "Proc %d: MPI_Barrier failed with error code %d\n", myid, error);
    PyErr_SetString(PyExc_RuntimeError, errmsg);    /*raise ValueError, errmsg*/
    return NULL;
  }  
  
  return Py_BuildValue("");
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
  {"get_processor_name", Get_processor_name, METH_VARARGS},              
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
   
  import_array();     /*Necessary for handling of NumPY structures  */
}

 
 


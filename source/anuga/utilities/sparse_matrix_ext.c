
#include "Python.h"

#include <stdio.h>   /* gets */
#include <stdlib.h>  /* atoi, malloc */
#include <string.h>  /* strcpy */
#include "numpy/arrayobject.h"
#include "math.h"

#include "util_ext.h" /* in utilities */
#include "sparse_dok.h"




static int _serialise(sparse_dok * dok, PyObject * serial_dok)
{

    sort_by_key(dok);
    int num_entries = dok->num_entries;
    int k,i,j;
    double val;

    edge_t * edge = dok->edgetable;

    for(k=0;k<num_entries;k++){

        i = edge->key.i;
        j = edge->key.j;
        val = edge->entry;

        // make new python objects to add to dict
        PyObject * py_val = PyFloat_FromDouble(val);
        PyObject * py_i = PyInt_FromLong((long) i);
        PyObject * py_j = PyInt_FromLong((long) j);
        PyObject * py_key = PyTuple_New(2);
        PyTuple_SET_ITEM(py_key,0,py_i);
        PyTuple_SET_ITEM(py_key,1,py_j);

        PyDict_SetItem(serial_dok,py_key,py_val);

        edge = edge->hh.next;

    }
   
    return 0;
};


static int _deserialise(sparse_dok * dok, PyObject * serial_dok)

{

    PyObject * items = PyDict_Items(serial_dok);
    PyObject * keys = PyDict_Keys(serial_dok);

    int num_entries = (int) PyDict_Size(serial_dok);
    printf("num entries %d\n",num_entries);
    int k,i,j;
    double val;
    PyObject * py_val;
    PyObject * py_key;

    edge_key_t key;



    for(k=0;k<num_entries;k++){

        py_key = PyList_GET_ITEM(keys,k);
        //py_val = PyList_GET_ITEM(items,k);
        py_val = PyDict_GetItem(serial_dok,py_key);

        val = PyFloat_AS_DOUBLE(py_val);
        printf("val: %d\n",val);

        i = (int) PyInt_AS_LONG(PyTuple_GET_ITEM(py_key,0));
        j = (int) PyInt_AS_LONG(PyTuple_GET_ITEM(py_key,1));

        key.i = i;
        key.j = j;

        add_dok_entry(dok,key,val);

    }

    Py_DECREF(keys);
    Py_DECREF(items);

    return 0;

};


// Delete capsule containing a sparse_dok - name of capsule must be exactly
// "sparse dok".
static void delete_sparse_dok_cap(PyObject *cap){

    sparse_dok * kill = (sparse_dok*) PyCapsule_GetPointer(cap,"sparse dok");
    if(kill!=NULL){
        delete_dok_matrix(kill);
    }

}
//----------------------- PYTHON WRAPPER FUNCTION -----------------------------


static PyObject *serialise_dok(PyObject *self, PyObject *args) {

    int err;
    PyObject *sparse_dok_cap;
    PyObject *serial_sparse_dok;

    // Convert Python arguments to C
    if (!PyArg_ParseTuple(args, "O",&sparse_dok_cap
                                            )) {
      PyErr_SetString(PyExc_RuntimeError,
              "sparse_matrix_ext.serialise_dok: could not parse input");
      return NULL;
    }

    sparse_dok * dok = (sparse_dok*) PyCapsule_GetPointer(sparse_dok_cap,"sparse dok");

    serial_sparse_dok = PyDict_New();

    err = _serialise(dok,serial_sparse_dok);

    if (err != 0) {
      PyErr_SetString(PyExc_RuntimeError,
              "sparse_matrix.serialise_dok: error in serialising sparse_dok");
      return NULL;
    }

    return serial_sparse_dok;

}

static PyObject *deserialise_dok(PyObject *self, PyObject *args) {

    int err;
    PyObject *sparse_dok_cap;
    PyObject *serial_sparse_dok;

    // Convert Python arguments to C
    if (!PyArg_ParseTuple(args, "O", &serial_sparse_dok
                                            )) {
      PyErr_SetString(PyExc_RuntimeError,
              "sparse_matrix.deserialise_dok: could not parse input");
      return NULL;
    }
    
    sparse_dok * dok;
    dok = make_dok();

    err = _deserialise(dok,serial_sparse_dok);

    return  PyCapsule_New((void*) dok,
                      "sparse dok",
                      &delete_sparse_dok_cap); 

}


//------------------------------------------------------------------------------


// ------------------------------ PYTHON GLUE ----------------------------------

//==============================================================================
// Structures to allow calling from python
//==============================================================================

// Method table for python module
static struct PyMethodDef MethodTable[] = {
    {"serialise_dok",serialise_dok, METH_VARARGS, "Print out"},
    {"deserialise_dok",deserialise_dok, METH_VARARGS, "Print out"},
	{NULL, NULL, 0, NULL}   // sentinel
};


// Module initialisation
void initsparse_matrix_ext(void){
  
  Py_InitModule("sparse_matrix_ext", MethodTable);
  import_array(); // Necessary for handling of NumPY structures

}

// --------------------------------------------------------------------------------
#include <Python.h>
#include <numpy/arrayobject.h>

/* This must be the same as the metis idxtype */
typedef int idxtype;

#include "bridge.h"

static PyObject * metis_partMeshNodal(PyObject *, PyObject *);

static PyMethodDef methods[] = {
  {"partMeshNodal", metis_partMeshNodal, METH_VARARGS, "METIS_PartMeshNodal"},
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initmetis_ext(void){
  (void) Py_InitModule("metis_ext", methods);

  import_array();
}

/* Run the metis METIS_PartMeshNodal function
 * expected args:
 * ne: number of elements
 * nn: number of nodes
 * elmnts: element array
 * etype: type of mesh elements:
 ** 1 - triangle
 ** 2 - tetrahedra
 ** 3 - hexahedra
 ** 4 - quadrilaterals
 * nparts: number of partitions
 * returns:
 * edgecut: number of cut edges
 * epart: partitioning of the elements
 * npart: partitioning of the nodes.
 *
 * Note that while the metis data file format indexes verticies from 1,
 * the library calls (including this one) use vericies indexed from 0.
 * Unsure if this is because of the num_flag option, perhaps calling
 * from a FORTRAN program will have 1-indexed verticies?
 */
static PyObject * metis_partMeshNodal(PyObject * self, PyObject * args){
  int i;
  int malloc_elem_c_arr = 0;
  int ne;
  int nn;
  int etype;
  int nparts;
  int edgecut;
  int numflag = 0; // The metis routine requires an int * for numflag.
  npy_intp dims[1]; // PyArray_SimpleNewFromData needs an npy_intp[] of array sizes.

  PyObject * elements;
  PyArrayObject * elem_arr;
  PyArrayObject * epart_pyarr;
  PyArrayObject * npart_pyarr;

  /* These are all of the metis idxtype */
  idxtype * elem_c_arr;
  idxtype * epart;
  idxtype * npart;
  if(!PyArg_ParseTuple(args, "iiOii", &ne, &nn, &elements, &etype, &nparts))
    return NULL;

  elem_arr = (PyArrayObject *) PyArray_ContiguousFromObject(elements, PyArray_NOTYPE, 1, 1);

  if(!elem_arr)
    return NULL;

  /* x86_64 will create arrays of longs and they need to be
   * converted to arrays of idxtype for metis to work on them.
   */
  if(elem_arr->descr->type_num == PyArray_LONG){
    elem_c_arr = (idxtype *)malloc(*(elem_arr->dimensions) * sizeof(idxtype));
    malloc_elem_c_arr = 1;
    if(!elem_c_arr)
	return NULL;
    for(i = 0 ; i < *(elem_arr->dimensions) ; i++){
      elem_c_arr[i] = (idxtype)(((long *)elem_arr->data)[i]);
      if(elem_c_arr[i] != ((long *)elem_arr->data)[i]){ /* i.e. downcast failed */
	free(elem_c_arr);
        Py_DECREF(elem_arr);
	return NULL;
      }
    }
  }else
    elem_c_arr = (idxtype *)elem_arr->data;

  epart = (idxtype *)malloc(ne * sizeof(idxtype));
  if(epart == NULL){
    if(malloc_elem_c_arr) free(elem_c_arr);
    Py_DECREF(elem_arr);
    return NULL;
  }
  npart = (idxtype *)malloc(nn * sizeof(idxtype));
  if(npart == NULL){
    if(malloc_elem_c_arr) free(elem_c_arr);
    free(epart);
    Py_DECREF(elem_arr);
    return NULL;
  }
  bridge_partMeshNodal(&ne, &nn, elem_c_arr, &etype, &numflag, &nparts, &edgecut, epart, npart);

  dims[0] = ne;
  epart_pyarr = (PyArrayObject *)PyArray_SimpleNewFromData(1, dims, PyArray_INT, (void *)epart);
  //epart_pyarr = (PyArrayObject *)PyArray_FromDimsAndData(1, dims, PyArray_INT, (char *)epart);
  dims[0] = nn;
  npart_pyarr = (PyArrayObject *)PyArray_SimpleNewFromData(1, dims, PyArray_INT, (void *)npart);
  //npart_pyarr = (PyArrayObject *)PyArray_FromDimsAndData(1, dims, PyArray_INT, (char *)npart);

  
  if(malloc_elem_c_arr) free(elem_c_arr);
  Py_DECREF(elem_arr);

  return Py_BuildValue("iOO", edgecut, (PyObject *)epart_pyarr, (PyObject *)npart_pyarr);
}

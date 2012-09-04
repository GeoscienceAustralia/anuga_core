#include "Python.h"
#include "numpy/arrayobject.h"
#include <stdio.h>

#define DDATA(p) ((double*)(((PyArrayObject *)p)->data))
#define IDATA(p) ((long*)(((PyArrayObject *)p)->data))




static PyObject *BoundaryDictionaryConstruct( PyObject *self, PyObject *args )
{
	int i, j, initSize, ok, numTriangle;
	char *defaultTag;
	long *neighbours;
	PyObject *pyobj_dictKey;
	PyObject *pyobj_neighbours;
	PyObject *pyobj_boundary;

	ok = PyArg_ParseTuple( args, "isOO", &numTriangle, &defaultTag, &pyobj_neighbours, &pyobj_boundary );
	if( !ok ){
		fprintf( stderr, "BuildBoundaryFromScratch func: argument parsing error\n" );
		exit(1);
	}

	neighbours = IDATA( pyobj_neighbours );	

	for(i=0; i<numTriangle; i++){
		for(j=0; j<3; j++){
			if( neighbours[i*3+j] < 0 ){
				pyobj_dictKey = PyTuple_New( 2 );
				PyTuple_SetItem( pyobj_dictKey, 0, PyInt_FromLong( i ) );
				PyTuple_SetItem( pyobj_dictKey, 1, PyInt_FromLong( j ) );
				
				if( !PyDict_Contains( pyobj_boundary, pyobj_dictKey ) )
					PyDict_SetItem( pyobj_boundary, pyobj_dictKey, PyString_FromString( defaultTag ) );
			}
		}
	}

	return Py_BuildValue("O", pyobj_boundary);
}

static PyMethodDef MF_module_methods[] = {
	{"boundary_dictionary_construct", BoundaryDictionaryConstruct, METH_VARARGS},
	{NULL, NULL}	//sentinel
};

void initneighbour_mesh_ext( )
{
	(void) Py_InitModule("neighbour_mesh_ext", MF_module_methods);

        import_array();
}

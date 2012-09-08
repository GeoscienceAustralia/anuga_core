#include "Python.h"
#include "numpy/arrayobject.h"
#include <stdio.h>

#include "util_ext.h"

#define DDATA(p) ((double*)(((PyArrayObject *)p)->data))
#define IDATA(p) ((long*)(((PyArrayObject *)p)->data))
#define DIMENSHAPE(p, d) (((PyArrayObject *)p)->dimensions[d])


#ifndef NDEBUG
#define ASSERT(condition, message) \
	do {\
		if(!(condition)) { \
			printf("Assertion `%s` failed in %s line %d: %s\n", \
			                  #condition, __FILE__, __LINE__, message); \
			exit(EXIT_FAILURE); \
		} \
	} while (0)
#else
#define ASSERT(condition, message) do { } while (0)
#endif

static PyObject *BoundaryDictionaryConstruct( PyObject *self, PyObject *args )
{
	int aDimen, bDimen, ok, numTriangle, volID, edgeID;
	char *defaultTag;
	char errorMsg[50];
	long *neighbours;
	PyObject *pyobj_dictKey, *pyobj_dictVal;
	PyObject *pyobj_neighbours;
	PyObject *pyobj_boundary;
	Py_ssize_t pos;

	ok = PyArg_ParseTuple( args, "isOO", &numTriangle, &defaultTag, &pyobj_neighbours, &pyobj_boundary );
	if( !ok ){
                report_python_error(AT, "could not parse input arguments");
                return NULL;
	}

	neighbours = IDATA( pyobj_neighbours );
	aDimen = DIMENSHAPE( pyobj_neighbours, 0 );
	bDimen = DIMENSHAPE( pyobj_neighbours, 1 );

	pos = 0;
	if( PyDict_Size( pyobj_boundary ) ){
		// iterate through dictionary entries
		while( PyDict_Next( pyobj_boundary, &pos, &pyobj_dictKey, &pyobj_dictVal ) ){
			volID  = PyLong_AsLong( PyTuple_GetItem( pyobj_dictKey, 0 ) );
			edgeID = PyLong_AsLong( PyTuple_GetItem( pyobj_dictKey, 1 ) );

                        if (!(volID<aDimen && edgeID<bDimen)) {
                            sprintf(errorMsg, "Segment (%d, %d) does not exist", volID, edgeID);
                            report_python_error(AT, errorMsg);
                            return NULL;
                        }
		}
	}

	for(volID=0; volID<numTriangle; volID++){
		for(edgeID=0; edgeID<3; edgeID++){
			if( neighbours[volID*3+edgeID] < 0 ){
				pyobj_dictKey = PyTuple_New( 2 );
				PyTuple_SetItem( pyobj_dictKey, 0, PyInt_FromLong( volID ) );
				PyTuple_SetItem( pyobj_dictKey, 1, PyInt_FromLong( edgeID ) );

				if( !PyDict_Contains( pyobj_boundary, pyobj_dictKey ) )
					PyDict_SetItem( pyobj_boundary, pyobj_dictKey, PyString_FromString( defaultTag ) );
			}
		}
	}

	return Py_BuildValue("O", pyobj_boundary);
}



/*
static PyObject *BoundaryDictionaryConstruct( PyObject *self, PyObject *args )
{
	int i, j, ok, numTriangle;
	char *defaultTag;
	long *neighbours;
	PyObject *pyobj_dictKey;
	PyObject *pyobj_neighbours;
	PyObject *pyobj_boundary;

	ok = PyArg_ParseTuple( args, "isOO", &numTriangle, &defaultTag, &pyobj_neighbours, &pyobj_boundary );
	if( !ok ){
                report_python_error(AT, "could not parse input arguments");
                return NULL;
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
*/

static PyMethodDef MF_module_methods[] = {
	{"boundary_dictionary_construct", BoundaryDictionaryConstruct, METH_VARARGS},
	{NULL, NULL}	//sentinel
};

void initneighbour_mesh_ext( )
{
	(void) Py_InitModule("neighbour_mesh_ext", MF_module_methods);

        import_array();
}

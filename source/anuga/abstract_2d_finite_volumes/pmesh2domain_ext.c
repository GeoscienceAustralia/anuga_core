#include "Python.h"
#include "numpy/arrayobject.h"
#include <stdio.h>
#include "util_ext.h"

#define DDATA(p) ((double*)(((PyArrayObject *)p)->data))
#define IDATA(p) ((long*)(((PyArrayObject *)p)->data))
#define LENDATA(p) ((long*)(((PyArrayObject *)p)->dimensions[0]))



static PyObject *SidesDictionaryConstruct( PyObject *self, PyObject *args )
{
	int i, ok, numTriangle;
        long a,b,c;
	long *triangles;
	PyObject *pyobj_key;
        PyObject *pyobj_value;
        PyObject *pyobj_sides;
	PyArrayObject *pyobj_triangles;


	ok = PyArg_ParseTuple( args, "OO", &pyobj_triangles, &pyobj_sides );
	if( !ok ){
                report_python_error(AT, "could not parse input arguments");
                return NULL;
	}


        numTriangle = LENDATA( pyobj_triangles );
	triangles = IDATA( pyobj_triangles );

        /*
        for id, triangle in enumerate(triangles):
           a = int(triangle[0])
           b = int(triangle[1])
           c = int(triangle[2])

           sides[a,b] = (id, 2) #(id, face)
           sides[b,c] = (id, 0) #(id, face)
           sides[c,a] = (id, 1) #(id, face)
         */

        printf("numTriangle %d\n",numTriangle);

	for(i=0; i<numTriangle; i++){
            a = triangles[i*3+0];
            b = triangles[i*3+1];
            c = triangles[i*3+2];

            // sides[a,b] = (id, 2) #(id, face)
	    pyobj_key = PyTuple_New( 2 );
	    PyTuple_SetItem( pyobj_key, 0, PyInt_FromLong( a ) );
	    PyTuple_SetItem( pyobj_key, 1, PyInt_FromLong( b ) );
            
	    pyobj_value = PyTuple_New( 2 );
	    PyTuple_SetItem( pyobj_value, 0, PyInt_FromLong( i ) );
	    PyTuple_SetItem( pyobj_value, 1, PyInt_FromLong( 2 ) );

	    PyDict_SetItem( pyobj_sides, pyobj_key, pyobj_value );

            Py_DECREF(pyobj_key);
            Py_DECREF(pyobj_value);
            
            // sides[b,c] = (id, 0) #(id, face)
	    pyobj_key = PyTuple_New( 2 );
	    PyTuple_SetItem( pyobj_key, 0, PyInt_FromLong( b ) );
	    PyTuple_SetItem( pyobj_key, 1, PyInt_FromLong( c ) );
            
	    pyobj_value = PyTuple_New( 2 );
	    PyTuple_SetItem( pyobj_value, 0, PyInt_FromLong( i ) );
	    PyTuple_SetItem( pyobj_value, 1, PyInt_FromLong( 0 ) );
            
	    PyDict_SetItem( pyobj_sides, pyobj_key, pyobj_value );

            Py_DECREF(pyobj_key);
            Py_DECREF(pyobj_value);
            
            // sides[c,a] = (id, 1) #(id, face)
	    pyobj_key = PyTuple_New( 2 );
	    PyTuple_SetItem( pyobj_key, 0, PyInt_FromLong( c ) );
	    PyTuple_SetItem( pyobj_key, 1, PyInt_FromLong( a) );
            
	    pyobj_value = PyTuple_New( 2 );
	    PyTuple_SetItem( pyobj_value, 0, PyInt_FromLong( i ) );
	    PyTuple_SetItem( pyobj_value, 1, PyInt_FromLong( 1 ) );
            
	    PyDict_SetItem( pyobj_sides, pyobj_key, pyobj_value );

            Py_DECREF(pyobj_key);
            Py_DECREF(pyobj_value);

        }

	return Py_BuildValue("O", pyobj_sides);
}

static PyMethodDef MF_module_methods[] = {
	{"sides_dictionary_construct", SidesDictionaryConstruct, METH_VARARGS},
	{NULL, NULL}	//sentinel
};

void initpmesh2domain_ext( )
{
	(void) Py_InitModule("pmesh2domain_ext", MF_module_methods);

        import_array();
}

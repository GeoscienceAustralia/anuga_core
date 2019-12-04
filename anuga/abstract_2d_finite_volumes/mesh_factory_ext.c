
#include "Python.h"
#include "numpy/arrayobject.h"
#include <stdio.h>

#define DDATA(p) ((double*)(((PyArrayObject *)p)->data))
#define IDATA(p) ((long*)(((PyArrayObject *)p)->data))

static PyObject *RectangularCrossConstruct( PyObject *self, PyObject *args )
{
	int m, n, i, j, v1, v2, v3, v4, v5, ok;
	int numPoints, numElements;
	double len1, len2, delta1, delta2, x, y;
	double *origin;
	double *points;
	double *params;
	long *vertices, *elements;
	PyObject *pyobj_Params;
	PyObject *pyobj_Origin;
	PyObject *pyobj_boundary;
	PyObject *pyobj_points, *pyobj_elements;
	PyObject *pyobj_tuple;
	
	ok = PyArg_ParseTuple( args, "OOOO", &pyobj_Params, &pyobj_Origin, &pyobj_points, &pyobj_elements );
	if( !ok ){
		fprintf( stderr, "RectangularCrossConstruct func: argument parsing error\n" );
		exit(1);
	}
	
	origin = DDATA( pyobj_Origin );
	params = DDATA( pyobj_Params );
	points = DDATA( pyobj_points );
	elements = IDATA( pyobj_elements );

	m = (int)params[0];
	n = (int)params[1];
	len1 = params[2];
	len2 = params[3];

	vertices = malloc( (m+1)*(n+1)*sizeof(long) );

	delta1 = len1/m;
	delta2 = len2/n;	

	numPoints = 0;
	for(i=0; i<m+1; i++)
		for(j=0; j<n+1; j++){
			vertices[(n+1)*i+j] = numPoints;
			points[numPoints*2] = i*delta1 + origin[0];
			points[numPoints*2+1] = j*delta2 + origin[1]; 
			numPoints++; 
		}
	
	pyobj_boundary = PyDict_New( );
	numElements = 0;
	for(i=0; i<m; i++)
		for(j=0; j<n; j++){
			v1 = vertices[(n+1)*i+j+1];
			v2 = vertices[(n+1)*i+j];
			v3 = vertices[(n+1)*(i+1)+j+1];
			v4 = vertices[(n+1)*(i+1)+j];
			x = (points[v1*2]+points[v2*2]+points[v3*2]+points[v4*2])*0.25;
			y = (points[v1*2+1]+points[v2*2+1]+points[v3*2+1]+points[v4*2+1])*0.25;

			//Create centre point
			v5 = numPoints;
			points[numPoints*2] = x;
			points[numPoints*2+1] = y;
			numPoints++;
			
			//Create left triangle
			if( i == 0 ){
				pyobj_tuple = PyTuple_New( 2 );
				PyTuple_SetItem( pyobj_tuple, 0, PyInt_FromLong( numElements ) );
				PyTuple_SetItem( pyobj_tuple, 1, PyInt_FromLong( 1 ) );
				PyDict_SetItem( pyobj_boundary, pyobj_tuple, PyString_FromString( "left" ) );
			}
			elements[numElements*3] = v2;
			elements[numElements*3+1] = v5;
			elements[numElements*3+2] = v1;
			numElements++;
			
			//Create bottom triangle
			if( j == 0 ){
				pyobj_tuple = PyTuple_New( 2 );
				PyTuple_SetItem( pyobj_tuple, 0, PyInt_FromLong( numElements ) );
				PyTuple_SetItem( pyobj_tuple, 1, PyInt_FromLong( 1 ) );
				PyDict_SetItem( pyobj_boundary, pyobj_tuple, PyString_FromString( "bottom" ) );
			}
			elements[numElements*3] = v4;
			elements[numElements*3+1] = v5;
			elements[numElements*3+2] = v2;
			numElements++;

			//Create right triangle
			if( i == (m-1) ){
				pyobj_tuple = PyTuple_New( 2 );
				PyTuple_SetItem( pyobj_tuple, 0, PyInt_FromLong( numElements ) );
				PyTuple_SetItem( pyobj_tuple, 1, PyInt_FromLong( 1 ) );
				PyDict_SetItem( pyobj_boundary, pyobj_tuple, PyString_FromString( "right" ) );
			}
			elements[numElements*3] = v3;
			elements[numElements*3+1] = v5;
			elements[numElements*3+2] = v4;
			numElements++;

			//Create top triangle
			if( j == (n-1) ){
				pyobj_tuple = PyTuple_New( 2 );
				PyTuple_SetItem( pyobj_tuple, 0, PyInt_FromLong( numElements ) );
				PyTuple_SetItem( pyobj_tuple, 1, PyInt_FromLong( 1 ) );
				PyDict_SetItem( pyobj_boundary, pyobj_tuple, PyString_FromString( "top" ) );
			}
			elements[numElements*3] = v1;
			elements[numElements*3+1] = v5;
			elements[numElements*3+2] = v3;
			numElements++;
		}

        free(vertices);

	return Py_BuildValue("O", pyobj_boundary);
}

static PyMethodDef MF_module_methods[] = {
	{"rectangular_cross_construct", RectangularCrossConstruct, METH_VARARGS},
	{NULL, NULL}
};

void initmesh_factory_ext( )
{
	(void) Py_InitModule("mesh_factory_ext", MF_module_methods);

        import_array();
}

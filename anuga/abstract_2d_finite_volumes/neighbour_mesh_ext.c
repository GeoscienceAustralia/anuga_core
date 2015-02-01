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



static PyObject *CheckIntegrity( PyObject *self, PyObject *args )
{
	int nt, nt3, tri, n_node, n_node_1;
        int ok;
        int current_node;
        int k, i;
        int index;
   
	//char errorMsg[50];

        double cumsum;

	long *vertex_value_indices;
        long *triangles;
        long *node_index;
        long *number_of_triangles_per_node;


	PyObject *pyobj_vertex_value_indices;
        PyObject *pyobj_triangles;
	PyObject *pyobj_node_index;
	PyObject *pyobj_number_of_triangles_per_node;

	//Py_ssize_t pos;

	ok = PyArg_ParseTuple( args, "OOOO",
                &pyobj_vertex_value_indices,
                &pyobj_triangles,
                &pyobj_node_index,
                &pyobj_number_of_triangles_per_node
                );
	if( !ok ){
                report_python_error(AT, "could not parse input arguments");
                return NULL;
	}

	vertex_value_indices = IDATA( pyobj_vertex_value_indices );
	nt3 = DIMENSHAPE( pyobj_vertex_value_indices, 0 );


        //printf(" n3 %d \n",nt3);

	triangles = IDATA( pyobj_triangles );
	nt = DIMENSHAPE( pyobj_triangles, 0 );
	tri = DIMENSHAPE( pyobj_triangles, 1 );

        //printf(" nt %d tri %d \n",nt, tri);

	node_index = IDATA( pyobj_node_index );
	n_node_1 = DIMENSHAPE( pyobj_node_index, 0 );


        //printf(" n_node_1 %d  \n",n_node_1);

	number_of_triangles_per_node = IDATA( pyobj_number_of_triangles_per_node );
	n_node = DIMENSHAPE( pyobj_number_of_triangles_per_node, 0 );


        //printf(" n_node %d \n",n_node);

        //----------------------------------------------------------------------
        // Check that the arrays are consistent
        //----------------------------------------------------------------------
        if ( nt3 != 3*nt) {
                report_python_error(AT, "Mismatch in size of triangles and vertex_value_indices");
                return NULL;
	}

        if ( n_node_1 != n_node + 1 ) {
                report_python_error(AT, "Mismatch in size of node_index and number_of_triangles_per_node");
                return NULL;
	}

        if ( tri != 3 ) {
                report_python_error(AT, "Triangle array should be nt by 3");
                return NULL;
	}

        //----------------------------------------------------------------------
        // Check consistence between triangles, and node_index and
        // number_of_triangles_per_node
        //----------------------------------------------------------------------
        current_node = 0;
        k = 0; //Track triangles touching on node
        for (i=0; i<nt3 ; i++) {
            index = vertex_value_indices[i];

            if (number_of_triangles_per_node[current_node] == 0) {
                // Node is lone - i.e. not part of the mesh
                continue;
            }

            k += 1;

            if ( triangles[index] != current_node) {
                report_python_error(AT, "Inconsistency between triangles and vertex_values_indices");
                return NULL;
            }

            if ( number_of_triangles_per_node[current_node] == k) {
                // Move on to next node
                k = 0;
                current_node += 1;
            }
        }

        cumsum = 0.0;
        for (i=0; i<n_node; i++) {
            cumsum += number_of_triangles_per_node[i];
            if ( cumsum != node_index[i+1]) {
                report_python_error(AT, "Inconsistency between node_index and number_of_triangles_per_node");
                return NULL;
            }
        }


/*
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
*/

	return Py_BuildValue("");
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
        {"check_integrity_c", CheckIntegrity, METH_VARARGS},
	{NULL, NULL}	//sentinel
};

void initneighbour_mesh_ext( )
{
	(void) Py_InitModule("neighbour_mesh_ext", MF_module_methods);

        import_array();
}

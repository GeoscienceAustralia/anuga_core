#include "Python.h"
#include "numpy/arrayobject.h"

#include <cstdio>   /* gets */
#include <cstdlib>  /* atoi, malloc */
#include <cstring>  /* strcpy */
//#include <cmath>    /* math!!! */

// This could be replaced with fast drop-inreplacements
// that are around and open. like https://github.com/greg7mdp/sparsepp
// TODO: But lets see if performance of unordered_map is any problem first.
//#include "sparsepp/spp.h"
#include <unordered_map> /* std::unordered_map */
#include <functional> /* std::hash */

//Shared code snippets
#include "util_ext.h" /* in utilities */

// basic type used for keys and counters
// should be the same type as triangle/coordinate ids
// used by the passed in arrays.
typedef long keyint;

struct edge_key_t {
    keyint i;
    keyint j;
    bool operator==(const edge_key_t & other) const {
        return i == other.i && j == other.j;
    }
};

// Bitmixer from MurmurHash3
uint64_t bitmix( uint64_t key )
{
  key ^= (key >> 33);
  key *= 0xff51afd7ed558ccd;
  key ^= (key >> 33);
  key *= 0xc4ceb9fe1a85ec53;
  key ^= (key >> 33);

  return key;
}

// custom hash function for the edge_keys
// used by the unordered_map edgetable
namespace std {
    template<> struct hash<edge_key_t> {
        std::size_t operator()(const edge_key_t & key) const {
            return bitmix(key.i) ^ bitmix(key.j);
        }
    };
}

struct edge_t {
    keyint vol_id;                /* id of vol containing this edge */
    keyint edge_id;               /* edge_id of edge in this vol */
};

//==============================================================================
// Code to calculate neighbour structure
//==============================================================================
int _build_neighbour_structure(keyint N, keyint M,
                      long* triangles,
		      long* neighbours,
                      long* neighbour_edges,
                      long* number_of_boundaries)
		      {
    keyint k;
    keyint k3;
    keyint n0,n1,n2;
    keyint vol_id;
    keyint edge_id;
    int err = 0;
    edge_key_t key;

    std::unordered_map<edge_key_t,edge_t> edgetable;
    //spp::sparse_hash_map<edge_key_t,edge_t> edgetable;
    
    // We know we'll have at least as many edges as triangles.
    edgetable.reserve(M);

    //--------------------------------------------------------------------------
    // Step 1:
    // Populate hashtable. We use a key based on the node_ids of the
    // two nodes defining the edge
    //--------------------------------------------------------------------------
    for (k=0; k<M; k++) {

        // Run through triangles
        k3 = 3*k;

        n0 = triangles[k3+0];
        n1 = triangles[k3+1];
        n2 = triangles[k3+2];

        // Add edges

        //----------------------------------------------------------------------
        // edge 0, n1 - n2
        //----------------------------------------------------------------------
        key.i = n1;
        key.j = n2;
        vol_id = k;
        edge_id = 0;
        
        // Error if duplicates
        if (edgetable.count(key)) {
            err = 1;
            break;
        }
        
        // Otherwise add edge
        edgetable.insert({key,{vol_id,edge_id}});
        

        //----------------------------------------------------------------------
        // edge 1, n2 - n0
        //----------------------------------------------------------------------
        key.i = n2;
        key.j = n0;
        vol_id = k;
        edge_id = 1;

        // Error if duplicates
        if (edgetable.count(key)) {
            err = 1;
            break;
        }
        // Otherwise add edge
        edgetable.insert({key,{vol_id,edge_id}});

        //----------------------------------------------------------------------
        // edge 2, n0 - n1
        //----------------------------------------------------------------------
        key.i = n0;
        key.j = n1;
        vol_id = k;
        edge_id = 2;

        // Error if duplicates
        if (edgetable.count(key)) {
            err = 1;
            break;
        }
        // Otherwise add edge
        edgetable.insert({key,{vol_id,edge_id}});

    }

    //--------------------------------------------------------------------------
    // return with an error code if duplicate key found
    // Clean up hashtable
    //--------------------------------------------------------------------------
    if (err) {
        return err;
    }


    //--------------------------------------------------------------------------
    //Step 2:
    //Go through triangles again, but this time
    //reverse direction of segments and lookup neighbours.
    //--------------------------------------------------------------------------
    for (k=0; k<M; k++) {

        // Run through triangles
        k3 = 3*k;

        n0 = triangles[k3+0];
        n1 = triangles[k3+1];
        n2 = triangles[k3+2];

        number_of_boundaries[k] = 3;

        // Search for neighbouring edge

        // edge 0, n1 - n2
        key.i = n2;
        key.j = n1;
        if (edgetable.count(key)) {
            edge_t & s = edgetable[key];
            neighbours[k3]      = s.vol_id;
            neighbour_edges[k3] = s.edge_id;
            number_of_boundaries[k] -= 1;
        }

        // edge 1, n2 - n0
        key.i = n0;
        key.j = n2;
        if (edgetable.count(key)) {
            edge_t & s = edgetable[key];
            neighbours[k3+1]      = s.vol_id;
            neighbour_edges[k3+1] = s.edge_id;
            number_of_boundaries[k] -= 1;
        }

        // edge 2, n0 - n1
        key.i = n1;
        key.j = n0;
        if (edgetable.count(key)) {
            edge_t & s = edgetable[key];
            neighbours[k3+2]      = s.vol_id;
            neighbour_edges[k3+2] = s.edge_id;
            number_of_boundaries[k] -= 1;
        }

    }
    return err;
}


//==============================================================================
// Python method Wrapper
//==============================================================================

PyObject *build_neighbour_structure(PyObject *self, PyObject *args) {

  /*
   * Update neighbours array using triangles array
   *
   * N is number of nodes (vertices)
   * triangle nodes defining triangles
   * neighbour across edge_id
   * neighbour_edges edge_id of edge in neighbouring triangle
   * number_of_boundaries
  */

	PyArrayObject *neighbours, *neighbour_edges, *number_of_boundaries;
        PyArrayObject *triangles;

	keyint N; // Number of nodes (read in)
        keyint M; // Number of triangles (calculated from triangle array)
        int err;


	// Convert Python arguments to C
	if (!PyArg_ParseTuple(args, "iOOOO", &N,
                                            &triangles,
                                            &neighbours,
                                            &neighbour_edges,
                                            &number_of_boundaries
                                            )) {
	  PyErr_SetString(PyExc_RuntimeError,
			  "hashtable.c: create_neighbours could not parse input");
	  return NULL;
	}


        CHECK_C_CONTIG(triangles);
        CHECK_C_CONTIG(neighbours);
        CHECK_C_CONTIG(neighbour_edges);
        CHECK_C_CONTIG(number_of_boundaries);


        M = triangles -> dimensions[0];


	err = _build_neighbour_structure(N, M,
                      (long*) triangles  -> data,
		      (long*) neighbours -> data,
                      (long*) neighbour_edges -> data,
                      (long*) number_of_boundaries -> data);


	if (err != 0) {
	  PyErr_SetString(PyExc_RuntimeError,
			  "Duplicate Edge");
	  return NULL;
	}


	return Py_BuildValue("");
}


//==============================================================================
// Structures to allow calling from python
//==============================================================================
extern "C" {
// Method table for python module
static struct PyMethodDef MethodTable[] = {
	{"build_neighbour_structure", build_neighbour_structure, METH_VARARGS, "Print out"},
	{NULL, NULL, 0, NULL}   // sentinel
};

// Module initialisation
void initneighbour_table_ext(void){
  Py_InitModule("neighbour_table_ext", MethodTable);
  import_array(); // Necessary for handling of NumPY structures
}
};

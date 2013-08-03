#include "Python.h"
#include "numpy/arrayobject.h"

#include <stdio.h>   /* gets */
#include <stdlib.h>  /* atoi, malloc */
#include <string.h>  /* strcpy */
#include <math.h>

//Shared code snippets

#include "util_ext.h" /* in utilities */
#include "uthash.h"     /* in utilities */


//==============================================================================
// hashtable code from uthash. Look at copyright info in "uthash.h in the
// utilities directory
//==============================================================================

typedef struct  {
    int i;
    int j;
} edge_key_t;

typedef struct {
    edge_key_t key;              /* key of form i , j */
    int vol_id;                /* id of vol containing this edge */
    int edge_id;               /* edge_id of edge in this vol */
    UT_hash_handle hh;         /* makes this structure hashable */
} edge_t;

edge_t *edgetable = NULL;

void add_edge(edge_key_t key, int vol_id, int edge_id) {
    edge_t *s;

    s = (edge_t*) malloc(sizeof(edge_t));
    memset(s, 0, sizeof(edge_t));
    s->key.i = key.i;
    s->key.j = key.j;
    s->vol_id = vol_id;
    s->edge_id = edge_id;
    HASH_ADD(hh, edgetable, key, sizeof(edge_key_t), s);  /* key: name of key field */
}

edge_t *find_edge(edge_key_t key) {
    edge_t *s;

    HASH_FIND(hh, edgetable, &key, sizeof(edge_key_t), s);  /* s: output pointer */
    return s;
}

void delete_edge(edge_t *edge) {
    HASH_DEL( edgetable, edge);  /* edge: pointer to deletee */
    free(edge);
}

void delete_edge_all() {
  edge_t *current_edge, *tmp;

  HASH_ITER(hh, edgetable, current_edge, tmp) {
    HASH_DEL(edgetable, current_edge);  /* delete it (edgetable advances to next) */
    free(current_edge);                 /* free it */
  } 
}

void print_edges() {
    edge_t *s;

    for(s=edgetable; s != NULL; s=(edge_t*)(s->hh.next)) {
        printf("edge key i %d i %d vol_id %d  edge_id %d\n",
                      s->key.i, s->key.j, s->vol_id, s->edge_id);
    }
}

int vol_id_sort(edge_t *a, edge_t *b) {
    return (a->vol_id - b->vol_id);
}

int key_sort(edge_t *a, edge_t *b) {
    return (a->key.i - b->key.i);
}

void sort_by_vol_id() {
    HASH_SORT(edgetable, vol_id_sort);
}

void sort_by_key() {
    HASH_SORT(edgetable, key_sort);
}


//==============================================================================
// Code to calculate neighbour structure
//==============================================================================
int _build_neighbour_structure(int N, int M,
                      long* triangles,
		      long* neighbours,
                      long* neighbour_edges,
                      long* number_of_boundaries)
		      {
    int k;
    int k3;
    int n0,n1,n2;
    int vol_id;
    int edge_id;
    int err = 0;
    edge_t *s;
    edge_key_t key;

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
        s = find_edge(key);
        if (s) {
            err = 1;
            break;
        }
        // Otherwise add edge
        add_edge(key, vol_id, edge_id);

        //----------------------------------------------------------------------
        // edge 1, n2 - n0
        //----------------------------------------------------------------------
        key.i = n2;
        key.j = n0;
        vol_id = k;
        edge_id = 1;

        // Error if duplicates
        s = find_edge(key);
        if (s) {
            err = 1;
            break;
        }
        // Otherwise add edge
        add_edge(key, vol_id, edge_id);

        //----------------------------------------------------------------------
        // edge 2, n0 - n1
        //----------------------------------------------------------------------
        key.i = n0;
        key.j = n1;
        vol_id = k;
        edge_id = 2;

        // Error if duplicates
        s = find_edge(key);
        if (s) {
            err = 1;
            break;
        }
        // Otherwise add edge
        add_edge(key, vol_id, edge_id);

    }

    
    //--------------------------------------------------------------------------
    // return with an error code if duplicate key found
    // Clean up hashtable
    //--------------------------------------------------------------------------
    if (err) {
        //printf("Duplicate Keys:\n");
        //printf("key.i %d key.j %d vol_id %d edge_id %d \n",
        //       s->key.i, s->key.j,s->vol_id,s->edge_id);
        //printf("key.i %d key.j %d vol_id %d edge_id %d \n",
        //       key.i,key.j,vol_id,edge_id);
        delete_edge_all();
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
        s = find_edge(key);
        if (s) {
            neighbours[k3]      = s -> vol_id;
            neighbour_edges[k3] = s -> edge_id;
            number_of_boundaries[k] -= 1;
        }

        // edge 1, n2 - n0
        key.i = n0;
        key.j = n2;
        s = find_edge(key);
        if (s) {
            neighbours[k3+1]      = s -> vol_id;
            neighbour_edges[k3+1] = s -> edge_id;
            number_of_boundaries[k] -= 1;
        }

        // edge 2, n0 - n1
        key.i = n1;
        key.j = n0;
        s = find_edge(key);
        if (s) {
            neighbours[k3+2]      = s -> vol_id;
            neighbour_edges[k3+2] = s -> edge_id;
            number_of_boundaries[k] -= 1;
        }

    }
    
    delete_edge_all();  /* free any structures */

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

	int N; // Number of nodes (read in)
        int M; // Number of triangles (calculated from triangle array)
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


//==============================================================================
// Main program to test hashtable code
//==============================================================================
/*
int main(int argc, char *argv[]) {
    char in[10];
    int key=1, running=1;
    int vol_id;
    int edge_id;
    struct edge *s;
    unsigned num_edges;

    while (running) {
        printf("1. add edge\n");
        printf("2. find edge\n");
        printf("3. delete edge\n");
        printf("4. delete all edges\n");
        printf("5. sort items by vol_id\n");
        printf("6. sort items by key\n");
        printf("7. print edges\n");
        printf("8. count edges\n");
        printf("9. quit\n");
        gets(in);
        switch(atoi(in)) {
            case 1:
                printf("key?\n");;
                key = atoi(gets(in));
                printf("vol_id?\n");
                vol_id = atoi(gets(in));
                printf("edge_id?\n");
                edge_id = atoi(gets(in));
                add_edge(key, vol_id, edge_id);
                break;
            case 2:
                printf("key?\n");
                s = find_edge(atoi(gets(in)));
                if (s) {
                printf("Key: %d\n", s->key);
                }
                else {
                   printf("key unknown\n");
                }
                break;
            case 3:
                printf("key?\n");
                s = find_edge(atoi(gets(in)));
                if (s) delete_edge(s);
                else printf("key unknown\n");
                break;
            case 4:
                delete_edge_all();
                break;
            case 5:
                sort_by_vol_id();
                break;
            case 6:
                sort_by_key();
                break;
            case 7:
                print_edges();
                break;
            case 8:
                num_edges=HASH_COUNT(edgetable);
                printf("there are %u edges\n", num_edges);
                break;
            case 9:
                running=0;
                break;
        }
    }

    delete_edge_all();
    return 0;
}
*/


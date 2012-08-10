#include <stdio.h>   /* gets */
#include <stdlib.h>  /* atoi, malloc */
#include <string.h>  /* strcpy */



#include "Python.h"
#include "numpy/arrayobject.h"
#include "math.h"

//Shared code snippets

#include "util_ext.h" /* in utilities */
#include "uthash.h"     /* in utilities */


//==============================================================================
// hashtable code from uthash. Look at copyright info in "uthash.h in the
// utilities directory
//==============================================================================

struct edge {
    int key;                    /* key of form i m + j */
    int vol_id;                /* id of vol containing this edge */
    int edge_id;               /* edge_id of edge in this vol */
    UT_hash_handle hh;         /* makes this structure hashable */
};

struct edge *edgetable = NULL;

void add_edge(int edge_key, int vol_id, int edge_id) {
    struct edge *s;

    s = (struct edge*)malloc(sizeof(struct edge));
    s->key = edge_key;
    s->vol_id = vol_id;
    s->edge_id = edge_id;
    HASH_ADD_INT( edgetable, key, s );  /* key: name of key field */
}

struct edge *find_edge(int edge_key) {
    struct edge *s;

    HASH_FIND_INT( edgetable, &edge_key, s );  /* s: output pointer */
    return s;
}

void delete_edge(struct edge *edge) {
    HASH_DEL( edgetable, edge);  /* user: pointer to deletee */
    free(edge);
}

void delete_all() {
  struct edge *current_edge, *tmp;

  HASH_ITER(hh, edgetable, current_edge, tmp) {
    HASH_DEL(edgetable,current_edge);  /* delete it (edgetable advances to next) */
    free(current_edge);            /* free it */
  } 
}

void print_edges() {
    struct edge *s;

    for(s=edgetable; s != NULL; s=(struct edge*)(s->hh.next)) {
        printf("edge key %d: vol_id %d  edge_id %d\n", s->key, s->vol_id, s->edge_id);
    }
}

int vol_id_sort(struct edge *a, struct edge *b) {
    return (a->vol_id - b->vol_id);
}

int key_sort(struct edge *a, struct edge *b) {
    return (a->key - b->key);
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
int _create_neighbours(int N,
		      long* neighbours,
                      long* neighbour_edges,
                      long* number_of_boundaries,
		      long* triangles) {
    int k;
    int k3;
    int n0,n1,n2;
    int key;
    int vol_id;
    int edge_id;
    struct edge *s;

    //--------------------------------------------------------------------------
    // Step 1:
    // Populate hashtable. We use a key based on the node_ids of the
    // two nodes defining the edge
    //--------------------------------------------------------------------------
    for (k=0; k<N; k++) {

        // Run through triangles
        k3 = 3*k;

        n0 = triangles[k3+0];
        n1 = triangles[k3+1];
        n2 = triangles[k3+2];

        // Add edges

        // edge 0, n1 - n2
        key = n1*N + n2;
        vol_id = k;
        edge_id = 0;

        add_edge(key, vol_id, edge_id);

        // edge 1, n2 - n0
        key = n2*N + n0;
        vol_id = k;
        edge_id = 1;

        add_edge(key, vol_id, edge_id);

        // edge 2, n0 - n1
        key = n0*N + n1;
        vol_id = k;
        edge_id = 2;

        add_edge(key, vol_id, edge_id);

    }

    //--------------------------------------------------------------------------
    //Step 2:
    //Go through triangles again, but this time
    //reverse direction of segments and lookup neighbours.
    //--------------------------------------------------------------------------
    for (k=0; k<N; k++) {

        // Run through triangles
        k3 = 3*k;

        n0 = triangles[k3+0];
        n1 = triangles[k3+1];
        n2 = triangles[k3+2];

        number_of_boundaries[k] = 3;

        // Search for neighbouring edge

        // edge 0, n1 - n2
        key = n2*N + n1;
        s = find_edge(key);
        if (s) {
            neighbours[k3+2]      = s -> vol_id;
            neighbour_edges[k3+2] = s -> edge_id;
            number_of_boundaries[k3+2] -= 1;
        }

        // edge 1, n2 - n0
        key = n0*N + n2;
        s = find_edge(key);
        if (s) {
            neighbours[k3+2]      = s -> vol_id;
            neighbour_edges[k3+2] = s -> edge_id;
            number_of_boundaries[k3+2] -= 1;
        }

        // edge 2, n0 - n1
        key = n1*N + n0;
        s = find_edge(key);
        if (s) {
            neighbours[k3+2]      = s -> vol_id;
            neighbour_edges[k3+2] = s -> edge_id;
            number_of_boundaries[k3+2] -= 1;
        }

    }
    
    delete_all();  /* free any structures */

    return 0;

}


//==============================================================================
// Python method Wrapper
//==============================================================================

PyObject *create_neighbours(PyObject *self, PyObject *args) {

  /*
   * Update neighbours array using triangles array
  */

	PyArrayObject *neighbours, *neighbour_edges, *number_of_boundaries;
        PyArrayObject *triangles;

	int N, err;


	// Convert Python arguments to C
	if (!PyArg_ParseTuple(args, "OOOO", &neighbours,
                                            &neighbour_edges,
                                            &number_of_boundaries,
                                            &triangles)) {
	  PyErr_SetString(PyExc_RuntimeError,
			  "hashtable.c: create_neighbours could not parse input");
	  return NULL;
	}

        CHECK_C_CONTIG(neighbours);
        CHECK_C_CONTIG(neighbour_edges);
        CHECK_C_CONTIG(number_of_boundaries);
	CHECK_C_CONTIG(triangles);

	N = triangles -> dimensions[0];

	err = _create_neighbours(N,
		      (long*) neighbours -> data,
                      (long*) neighbour_edges -> data,
                      (long*) number_of_boundaries -> data,
		      (long*) triangles  -> data);


	if (err != 0) {
	  PyErr_SetString(PyExc_RuntimeError,
			  "hashtable.c: error creating neighbours");
	  return NULL;
	}

	// Release and return
	Py_DECREF(triangles);
	Py_DECREF(neighbours);
        Py_DECREF(neighbour_edges);
        Py_DECREF(number_of_boundaries);

	return Py_BuildValue("");
}


//==============================================================================
// Structures to allow calling from python
//==============================================================================

// Method table for python module
static struct PyMethodDef MethodTable[] = {
	{"create_neighbours", create_neighbours, METH_VARARGS, "Print out"},
	{NULL, NULL, 0, NULL}   // sentinel
};

// Module initialisation
void inithashtable(void){
  Py_InitModule("hashtable", MethodTable);
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
                delete_all();
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

    delete_all();  
    return 0;
}
*/


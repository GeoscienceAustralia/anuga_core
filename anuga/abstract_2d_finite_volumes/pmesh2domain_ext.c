#include "Python.h"
#include "numpy/arrayobject.h"

#include <stdio.h>   /* gets */
#include <stdlib.h>  /* atoi, malloc */
#include <string.h>  /* strcpy */
#include <math.h>

//Shared code snippets

#include "util_ext.h" /* in utilities */
#include "uthash.h"     /* in utilities */



#define DDATA(p) ((double*)(((PyArrayObject *)p)->data))
#define IDATA(p) ((long*)(((PyArrayObject *)p)->data))
#define CDATA(p) ((char*)(((PyArrayObject *)p)->data))
#define LENDATA(p) ((long)(((PyArrayObject *)p)->dimensions[0]))

//==============================================================================
// hashtable code from uthash. Look at copyright info in "uthash.h in the
// utilities directory
//==============================================================================

typedef struct {
    int i;
    int j;
} segment_key_t;

typedef struct {
    segment_key_t key; /* key of form i , j */
    int vol_id; /* id of vol containing this segement */
    int edge_id; /* edge_id of segement in this vol */
    UT_hash_handle hh; /* makes this structure hashable */
} segment_t;

segment_t *segment_table = NULL;

void add_segment(segment_key_t key, int vol_id, int edge_id) {
    segment_t *s;

    s = (segment_t*) malloc(sizeof (segment_t));
    memset(s, 0, sizeof (segment_t));
    s->key.i = key.i;
    s->key.j = key.j;
    s->vol_id = vol_id;
    s->edge_id = edge_id;
    HASH_ADD(hh, segment_table, key, sizeof (segment_key_t), s); /* key: name of key field */
}

segment_t *find_segment(segment_key_t key) {
    segment_t *s=0;

    HASH_FIND(hh, segment_table, &key, sizeof (segment_key_t), s); /* s: output pointer */
    return s;
}

void delete_segment(segment_t *segment) {
    HASH_DEL(segment_table, segment); /* user: pointer to deletee */
    free(segment);
}

void delete_segment_all() {
    segment_t *current_segment, *tmp;

    HASH_ITER(hh, segment_table, current_segment, tmp) {
        HASH_DEL(segment_table, current_segment); /* delete it (segment_table advances to next) */
        free(current_segment); /* free it */
    }
}

void print_segments() {
    segment_t *s;

    for (s = segment_table; s != NULL; s = (segment_t*) (s->hh.next)) {
        printf("segment key i %d j %d vol_id %d  edge_id %d\n",
                s->key.i, s->key.j, s->vol_id, s->edge_id);
    }
}

int vol_id_sort(segment_t *a, segment_t *b) {
    return (a->vol_id - b->vol_id);
}

int key_sort(segment_t *a, segment_t *b) {
    return (a->key.i - b->key.i);
}

void sort_by_vol_id() {
    HASH_SORT(segment_table, vol_id_sort);
}

void sort_by_key() {
    HASH_SORT(segment_table, key_sort);
}


//==============================================================================
// Python method Wrapper
//==============================================================================

PyObject *build_boundary_dictionary(PyObject *self, PyObject *args) {

    /*
     * Update neighbours array using triangles array
     *
     * N is number of nodes (vertices)
     * triangle nodes defining triangles
     * neighbour across edge_id
     * neighbour_segments edge_id of segment in neighbouring triangle
     * number_of_boundaries
     */

    PyArrayObject *pyobj_segments;
    PyArrayObject *pyobj_triangles;

    PyObject *pyobj_segment_tags; // List of strings
    PyObject *pyobj_tag_dict;
    PyObject *pyobj_key;
    PyObject *pyobj_tag;

    long *triangles;
    long *segments;

    int M; // Number of triangles (calculated from triangle array)
    int N; // Number of segments (calculated from segments array)
    int err = 0;

    int k;
    int k3;
    int k2;
    int a,b,c;
    int vol_id;
    int edge_id;
    int len_tag;
    char *string;
    segment_t *s;
    segment_key_t key;

    // Convert Python arguments to C
    if (!PyArg_ParseTuple(args, "OOOO",
            &pyobj_triangles,
            &pyobj_segments,
            &pyobj_segment_tags,
            &pyobj_tag_dict
            )) {
        PyErr_SetString(PyExc_RuntimeError,
                "pmesh2domain.c: build_boundary_dictionary could not parse input");
        return NULL;
    }

    M = LENDATA(pyobj_triangles);
    N = LENDATA(pyobj_segments);

    triangles = IDATA(pyobj_triangles);
    segments  = IDATA(pyobj_segments);

    //--------------------------------------------------------------------------
    // Step 1:
    // Populate hashtable. We use a key based on the node_ids of the
    // two nodes defining the segment
    //--------------------------------------------------------------------------
    for (k = 0; k < M; k++) {

        // Run through triangles
        k3 = 3 * k;

        a = triangles[k3 + 0];
        b = triangles[k3 + 1];
        c = triangles[k3 + 2];



        // Add segments to hashtable

        //----------------------------------------------------------------------
        // Segment a,b
        //----------------------------------------------------------------------
        key.i = a;
        key.j = b;
        vol_id = k;
        edge_id = 2;

        // Error if duplicates
        s = find_segment(key);

        //printf("k = %d segment %d %d \n",k,a,b);
        if (s) {
            err = 1;
            break;
        }
        // Otherwise add segment
        add_segment(key, vol_id, edge_id);

        //----------------------------------------------------------------------
        // segment b,c
        //----------------------------------------------------------------------
        key.i = b;
        key.j = c;
        vol_id = k;
        edge_id = 0;

        // Error if duplicates
        s = find_segment(key);
        //printf("k = %d segment %d %d \n",k,b,c);
        if (s) {
            err = 1;
            break;
        }
        // Otherwise add segment
        add_segment(key, vol_id, edge_id);

        //----------------------------------------------------------------------
        // segment c,a
        //----------------------------------------------------------------------
        key.i = c;
        key.j = a;
        vol_id = k;
        edge_id = 1;

        // Error if duplicates
        s = find_segment(key);

        //printf("k = %d segment %d %d \n",k,c,a);
        if (s) {
            err = 1;
            break;
        }
        // Otherwise add segment
        add_segment(key, vol_id, edge_id);

    }

    //printf("err %d \n",err);
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
        delete_segment_all();
        PyErr_SetString(PyExc_RuntimeError,
                "pmesh2domain.c: build_boundary_dictionary Duplicate segments");
        return NULL;
    }


    //--------------------------------------------------------------------------
    //Step 2:
    //Go through segments. Those with null tags are added to pyobj_tag_dict
    //--------------------------------------------------------------------------


    for (k = 0; k < N; k++) {

        k2 = 2*k;
        a = segments[k2+0];
        b = segments[k2+1];

        // pyobj_tag should be a string
        pyobj_tag = PyList_GetItem(pyobj_segment_tags, (Py_ssize_t) k );

        string = PyString_AsString(pyobj_tag);

        len_tag = strlen(string);

        //printf("segment %d %d len %d, tag %s \n",a,b,len_tag, string);

        key.i = a;
        key.j = b;
        s = find_segment(key);
        if ( s &&  len_tag>0 ) {
            vol_id  = s->vol_id;
            edge_id = s->edge_id;

            pyobj_key = PyTuple_New(2);
            PyTuple_SetItem(pyobj_key, 0, PyInt_FromLong( vol_id ));
            PyTuple_SetItem(pyobj_key, 1, PyInt_FromLong( edge_id ));

            PyDict_SetItem(pyobj_tag_dict, pyobj_key, pyobj_tag );
        }

        key.i = b;
        key.j = a;
        s = find_segment(key);
        if ( s &&  len_tag>0 ) {
            vol_id  = s->vol_id;
            edge_id = s->edge_id;

            pyobj_key = PyTuple_New(2);
            PyTuple_SetItem(pyobj_key, 0, PyInt_FromLong( vol_id ));
            PyTuple_SetItem(pyobj_key, 1, PyInt_FromLong( edge_id ));

            PyDict_SetItem(pyobj_tag_dict, pyobj_key, pyobj_tag );
        }

    }

    delete_segment_all(); /* free any structures */

    return Py_BuildValue("O", pyobj_tag_dict);
}


//==============================================================================
// Structures to allow calling from python
//==============================================================================

static PyObject *build_sides_dictionary(PyObject *self, PyObject *args) {
    long i, numTriangle;
    int ok;
    long a, b, c;
    long *triangles;
    PyObject *pyobj_key;
    PyObject *pyobj_sides;
    PyArrayObject *pyobj_triangles;


    ok = PyArg_ParseTuple(args, "OO", &pyobj_triangles, &pyobj_sides);
    if (!ok) {
        report_python_error(AT, "could not parse input arguments");
        return NULL;
    }


    numTriangle = LENDATA(pyobj_triangles);
    triangles = IDATA(pyobj_triangles);



    for (i = 0; i < numTriangle; i++) {
        a = triangles[i * 3 + 0];
        b = triangles[i * 3 + 1];
        c = triangles[i * 3 + 2];

        // sides[a,b] = (id, 2) #(id, face)
        pyobj_key = PyTuple_New(2);
        PyTuple_SetItem(pyobj_key, 0, PyInt_FromLong(a));
        PyTuple_SetItem(pyobj_key, 1, PyInt_FromLong(b));

        PyDict_SetItem(pyobj_sides, pyobj_key, PyInt_FromLong(3 * i + 2));

        Py_DECREF(pyobj_key);

        // sides[b,c] = (id, 0) #(id, face)
        pyobj_key = PyTuple_New(2);
        PyTuple_SetItem(pyobj_key, 0, PyInt_FromLong(b));
        PyTuple_SetItem(pyobj_key, 1, PyInt_FromLong(c));

        PyDict_SetItem(pyobj_sides, pyobj_key, PyInt_FromLong(3 * i + 0));

        Py_DECREF(pyobj_key);


        // sides[c,a] = (id, 1) #(id, face)
        pyobj_key = PyTuple_New(2);
        PyTuple_SetItem(pyobj_key, 0, PyInt_FromLong(c));
        PyTuple_SetItem(pyobj_key, 1, PyInt_FromLong(a));

        PyDict_SetItem(pyobj_sides, pyobj_key, PyInt_FromLong(3 * i + 1));

        Py_DECREF(pyobj_key);


    }

    return Py_BuildValue("O", pyobj_sides);
}


static PyMethodDef MF_module_methods[] = {
{"build_boundary_dictionary", build_boundary_dictionary, METH_VARARGS},
{"build_sides_dictionary", build_sides_dictionary, METH_VARARGS},
                {NULL, NULL} //sentinel
            };

void initpmesh2domain_ext() {
    (void) Py_InitModule("pmesh2domain_ext", MF_module_methods);

    import_array();
}

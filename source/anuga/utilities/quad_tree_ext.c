
#include "Python.h"

#include <stdio.h>   /* gets */
#include <stdlib.h>  /* atoi, malloc */
#include <string.h>  /* strcpy */
#include "numpy/arrayobject.h"
#include "math.h"

#include "util_ext.h" /* in utilities */
#include "quad_tree.h"

// PYVERSION273 used to check python version for use of PyCapsule
#if PY_MAJOR_VERSION>=2 && PY_MINOR_VERSION>=7 && PY_MICRO_VERSION>=3
    #define PYVERSION273
#endif


/*
static int _serialise(quad_tree * quadtree, PyObject * serial_quadtree)
{

    
    quad_tree_ll * nodelist = new_quad_tree_ll(quadtree,0);
    quad_tree_ll * last = nodelist;
    quad_tree_ll * temp;
    int i;
    int nlist = 0;
    while(nodelist !=NULL){
        
        PyObject * nodedata = PyList_New(2);
        quadtree=nodelist->tree;
        // if children have been added, add to the linked list
        PyObject * child_index;
        if (quadtree->q[0]!=NULL){
            child_index = PyList_New(4);
            for (i=0;i<4;i++){
                nlist++;
                quad_tree_ll * child = new_quad_tree_ll(quadtree->q[i],nlist);
                last->next=child;
                last=child;
                PyList_SET_ITEM(child_index,i,PyInt_FromLong((long) nlist));
            }
        }
        else{
            child_index = PyList_New(1);
            PyList_SET_ITEM(child_index,0,PyInt_FromLong(0));
        }
        // get list of triangle indicies
        triangle * T = quadtree->leaves;
        int n_tris = 0 ;
        while(T!=NULL)
        {
            n_tris+=1;
            T=T->next;
        }
        PyObject * tri_list = PyList_New(n_tris);
        T = quadtree->leaves;
        i=0;
        while(T!=NULL)
        {
            PyList_SET_ITEM(tri_list,i,PyInt_FromLong((long)T->index));
            T=T->next;
            i++;
        }

        // set data
        PyList_SET_ITEM(nodedata,0,child_index); 
        PyList_SET_ITEM(nodedata,1,tri_list);
        // add data to overall list
        PyList_SET_ITEM(serial_quadtree, nodelist->index, nodedata);
        
        // set up next quad_tree to visit
        temp = nodelist;
        nodelist=nodelist->next;
        free(temp);
    }

    return 0;
}

static void _deserialise_recursive(quad_tree * quadtree, PyListObject * serial_quadtree,
        long * triangles,
        double *vertex_coordinates,
        long n)
{

    int k,k6;
    double x0,y0,x1,y1,x2,y2;
    PyObject * node = (PyObject*) PyList_GET_ITEM(serial_quadtree,n);
    PyObject * children = (PyObject*) PyList_GET_ITEM(node,0);
    PyObject * leaves = (PyObject*) PyList_GET_ITEM(node,1);

    if((long)PyList_Size(children) == 4){
        // has children
        quad_tree_make_children(quadtree);
        int i;
        for(i=0;i<4;i++){
            _deserialise_recursive(quadtree->q[i], serial_quadtree,
                triangles,
                vertex_coordinates,PyInt_AsLong(PyList_GET_ITEM(children,i)));
        }
    }

    // add leaves
    long n_leaves = (long) PyList_Size(leaves);
    int i;
    for(i=0;i<n_leaves;i++){
        k = (int) PyInt_AsLong(PyList_GET_ITEM(leaves,i));
        k6=k*6;
        x0 = vertex_coordinates[k6];
        y0 = vertex_coordinates[k6 + 1];
        x1 = vertex_coordinates[k6 + 2];
        y1 = vertex_coordinates[k6 + 3];
        x2 = vertex_coordinates[k6 + 4];
        y2 = vertex_coordinates[k6 + 5];
        triangle * T = new_triangle(k,x0,y0,x1,y1,x2,y2);
        quad_tree_add_triangle_to_list(quadtree,T);

    }
};

static int _deserialise(quad_tree * quadtree, PyListObject * serial_quadtree,
        long * triangles,
        double *vertex_coordinates)

{

    int k,k6;
    double x0,y0,x1,y1,x2,y2;
    long n = 0;
    // if the root node has children, build children, and then call recursive
    // build on them
    PyObject * node = (PyObject*) PyList_GET_ITEM(serial_quadtree,n);
    PyObject * children = (PyObject*) PyList_GET_ITEM(node,0);
    PyObject * leaves = (PyObject*) PyList_GET_ITEM(node,1);

    if((long)PyList_Size(children) == 4){
        // has children
        quad_tree_make_children(quadtree);
        int i;
        for(i=0;i<4;i++){
            _deserialise_recursive(quadtree->q[i], serial_quadtree,
                triangles,
                vertex_coordinates,PyInt_AsLong(PyList_GET_ITEM(children,i)));
        }
    }

    // add leaves
    long n_leaves = (long) PyList_Size(leaves);
    int i;
    for(i=0;i<n_leaves;i++){
        k = (int) PyInt_AsLong(PyList_GET_ITEM(leaves,i));
        k6=k*6;
        x0 = vertex_coordinates[k6];
        y0 = vertex_coordinates[k6 + 1];
        x1 = vertex_coordinates[k6 + 2];
        y1 = vertex_coordinates[k6 + 3];
        x2 = vertex_coordinates[k6 + 4];
        y2 = vertex_coordinates[k6 + 5];
        triangle * T = new_triangle(k,x0,y0,x1,y1,x2,y2);
        quad_tree_add_triangle_to_list(quadtree,T);
    } 

    return 0;

}

*/

#ifdef PYVERSION273


void delete_quad_tree_cap(PyObject * cap){
    quad_tree * kill = (quad_tree*) PyCapsule_GetPointer(cap,"quad tree");
    if(kill!=NULL){
        delete_quad_tree(kill);
    } else{
    }
}



#else

// If using python earlier version, build with PyCObject

// Delete cobj containing a quad tree
void delete_quad_tree_cobj(void * cobj){
    quad_tree * kill = (quad_tree*) cobj;
    if(kill!=NULL){
        delete_quad_tree(kill);
    }
}


#endif

//----------------------- PYTHON WRAPPER FUNCTION -----------------------------

/*
static PyObject *serialise(PyObject *self, PyObject *args) {

    int err;
    PyObject *tree;
    PyObject *serial_quadtree;

    // Convert Python arguments to C
    if (!PyArg_ParseTuple(args, "O",&tree
                                            )) {
      PyErr_SetString(PyExc_RuntimeError,
              "quad_tree_ext.serialise: could not parse input");
      return NULL;
    }

    quad_tree * quadtree = (quad_tree*) PyCapsule_GetPointer(tree,"quad tree");
    
    int nodes = quad_tree_node_count(quadtree);
    serial_quadtree = PyList_New(nodes);

    err = _serialise(quadtree,serial_quadtree);

    if (err != 0) {
      PyErr_SetString(PyExc_RuntimeError,
              "quad_tree_ext.serialise: error in serialising quad tree");
      return NULL;
    }

    return serial_quadtree;

}

static PyObject *deserialise(PyObject *self, PyObject *args) {

    int err;
    PyArrayObject *triangles;
    PyArrayObject *vertex_coordinates;
    PyArrayObject *extents;
    PyListObject *serial_quadtree;

    // Convert Python arguments to C
    if (!PyArg_ParseTuple(args, "OOOO",&serial_quadtree,
                                            &triangles,
                                            &vertex_coordinates,
                                            &extents
                                            )) {
      PyErr_SetString(PyExc_RuntimeError,
              "quad_tree_ext.deserialise: could not parse input");
      return NULL;
    }
    double * extents_v = (double *) extents->data;
    quad_tree *quadtree = new_quad_tree(extents_v[0],extents_v[1],extents_v[2],extents_v[3]);
    err = _deserialise(quadtree, serial_quadtree,
        (long*) triangles->data,
        (double*) vertex_coordinates -> data);

    return  PyCapsule_New((void*) quadtree,
                      "quad tree",
                      &delete_quad_tree_cap); 

}
*/

//------------------------------------------------------------------------------


// ------------------------------ PYTHON GLUE ----------------------------------

//==============================================================================
// Structures to allow calling from python
//==============================================================================

// Method table for python module
static struct PyMethodDef MethodTable[] = {
  //  {"serialise",serialise, METH_VARARGS, "Print out"},
  //  {"deserialise",deserialise, METH_VARARGS, "Print out"},
	{NULL, NULL, 0, NULL}   // sentinel
};


// Module initialisation
void initquad_tree_ext(void){
  
  Py_InitModule("quad_tree_ext", MethodTable);
  import_array(); // Necessary for handling of NumPY structures

}

// --------------------------------------------------------------------------------

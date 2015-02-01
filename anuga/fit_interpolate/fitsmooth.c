
#include "Python.h"

#include <stdio.h>   /* gets */
#include <stdlib.h>  /* atoi, malloc */
#include <string.h>  /* strcpy */
#include "numpy/arrayobject.h"
#include "math.h"

#include "util_ext.h" /* in utilities */
#include "sparse_dok.h" /* in utilities */
#include "quad_tree.h"  /* in utilities */
#include "omp.h"

// Errors defined for netcdf reading
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

#include "patchlevel.h"

// PYVERSION273 used to check python version for use of PyCapsule
#if PY_MAJOR_VERSION>=2 && PY_MINOR_VERSION>=7 && PY_MICRO_VERSION>=3
    #define PYVERSION273
#endif

//-------------------------- QUANTITY FITTING ------------------------------


// Builds the matrix D used to smooth the interpolation 
// of a variables from scattered data points to a mesh. See fit.py for more details.s
int _build_smoothing_matrix(int n,
                      long* triangles,
        		      double* areas,
                      double* vertex_coordinates,
                      int* strides,
                      sparse_dok * smoothing_mat)
		      {


    int k;
    int k3,k6;
    int err = 0;
    edge_key_t key;

    double det,area,x0,x1,x2,y0,y1,y2;
    double a0,b0,a1,b1,a2,b2,e01,e12,e20;
    int v0,v1,v2;
    double smoothing_val;

    
    for(k=0; k<n; k++) {
        // multiple k by 3 to give index for triangles which store in 3-blocks
        k3=k*3;
        k6=k*6;
        // store the area for the current triangle
        area = areas[k];
        // store current triangles global vertex indicies
        v0 = triangles[k3];
        v1 = triangles[k3+1];
        v2 = triangles[k3+2];
        // store the locations of the three verticies
        x0 = vertex_coordinates[k6];
        y0 = vertex_coordinates[k6+1];
        x1 = vertex_coordinates[k6+2];
        y1 = vertex_coordinates[k6+3];
        x2 = vertex_coordinates[k6+4];
        y2 = vertex_coordinates[k6+5];

        // calculate gradients (move to external function?)

        det = (y2-y0)*(x1-x0) - (y1-y0)*(x2-x0);
        a0 = (y2-y0)*(0-1) - (y1-y0)*(0-1);
        a0 /= det;

        b0 = (x1-x0)*(0-1) - (x2-x0)*(0-1);
        b0 /= det;

        a1 = (y2-y0)*(1-0) - (y1-y0)*(0-0);
        a1 /= det;

        b1 = (x1-x0)*(0-0) - (x2-x0)*(1-0);
        b1 /= det;

        a2 = (y2-y0)*(0-0) - (y1-y0)*(1-0);
        a2 /= det;

        b2 = (x1-x0)*(1-0) - (x2-x0)*(0-0);
        b2 /= det;
        
        // insert diagonal contributions

        // v0,v0
        key.i = v0;
        key.j = v0;
        smoothing_val = (a0*a0 + b0*b0)*area;
        add_dok_entry(smoothing_mat,key,smoothing_val);

        // v1,v1
        key.i = v1;
        key.j = v1;
        smoothing_val = (a1*a1 + b1*b1)*area;
        add_dok_entry(smoothing_mat,key,smoothing_val);

        // v2,v2
        key.i = v2;
        key.j = v2;
        smoothing_val = (a2*a2 + b2*b2)*area;
        add_dok_entry(smoothing_mat,key,smoothing_val);


        // insert off diagonal contributions
        e01 = (a0*a1 + b0*b1)*area;
        // v0,v1 (v1,v0)
        key.i = v0;
        key.j = v1;
        add_dok_entry(smoothing_mat,key,e01);
        key.i = v1;
        key.j = v0;
        add_dok_entry(smoothing_mat,key,e01);

        e12 = (a1*a2 + b1*b2)*area;
        // v1,v2 (v2,v1)
        key.i = v1;
        key.j = v2;
        add_dok_entry(smoothing_mat,key,e12);
        key.i = v2;
        key.j = v1;
        add_dok_entry(smoothing_mat,key,e12);

        e20 = (a2*a0 + b2*b0)*area;
        // v2,v0 (v0,v2)
        key.i = v2;
        key.j = v0;
        add_dok_entry(smoothing_mat,key,e20);
        key.i = v0;
        key.j = v2;
        add_dok_entry(smoothing_mat,key,e20);
    }

    return err;

}

// Builds a quad tree out of a list of triangles for quick 
// searching. 
quad_tree * _build_quad_tree(int n,
                      long* triangles,
                      double* vertex_coordinates,
                      double* extents)               
{   
    
    int k,k6;
    double x0,y0,x1,y1,x2,y2;

    // set up quad tree and allocate memory
    quad_tree * tree = new_quad_tree(extents[0],extents[1],extents[2],extents[3]);
    
    // iterate through triangles
    for(k=0; k<n; k++) {
        // multiple k by 3 to give index for triangles which store in 3-blocks
        k6=k*6;
        // store the locations of the three verticies
        x0 = vertex_coordinates[k6];
        y0 = vertex_coordinates[k6 + 1];
        x1 = vertex_coordinates[k6 + 2];
        y1 = vertex_coordinates[k6 + 3];
        x2 = vertex_coordinates[k6 + 4];
        y2 = vertex_coordinates[k6 + 5];
        triangle * T = new_triangle(k,x0,y0,x1,y1,x2,y2);
        quad_tree_insert_triangle(tree,T);
    }
  
    // return pointer to new tree struct
    return tree;
    
}

// Builds the AtA and Atz interpolation matrix
// and residual. Uses a quad_tree for fast access to the triangles of the mesh.
// This function takes a list of point coordinates, and associated point values
// (for any number of attributes).
int _build_matrix_AtA_Atz_points(int N, long * triangles,
                      double * point_coordinates, double * point_values,
                      int zdims, int npts,
                      sparse_dok * AtA,
                      double ** Atz,quad_tree * quadtree)
              {



    int k;



    int i,w;
    for(w=0;w<zdims;w++){
        for(i=0;i<N;i++){
            Atz[w][i]=0;
        } 
    }

    edge_key_t key;




    #pragma omp parallel for private(k,i,key,w)
    for(k=0;k<npts;k++){


        double x = point_coordinates[2*k];
        double y = point_coordinates[2*k+1];
        triangle * T = search(quadtree,x,y);

        if(T!=NULL){
            double * sigma = calculate_sigma(T,x,y);
            int js[3];
            for(i=0;i<3;i++){
                js[i]=triangles[3*(T->index)+i];
            }
            
            #pragma omp critical
            { 
            for(i=0;i<3;i++){

               for(w=0;w<zdims;w++){
                    Atz[w][js[i]] += sigma[i]*point_values[zdims*k+w];
               }
               
               for(w=0;w<3;w++){
                    
                    key.i=js[i];
                    key.j=js[w];

                   add_dok_entry(AtA,key,sigma[i]*sigma[w]);
                }                        
            }
            }
            free(sigma);
            sigma=NULL;

       } 
    }

    return 0;
}

// Combines two sparse_dok matricies and two vectors of doubles. 
void _combine_partial_AtA_Atz(sparse_dok * dok_AtA1,sparse_dok * dok_AtA2,
                             double* Atz1,
                             double* Atz2,
                             int n, int zdim){

    add_sparse_dok(dok_AtA1,1,dok_AtA2,1);

    int i;
    for(i=0;i<n*zdim;i++){
        Atz1[i]+=Atz2[i];
    }
}

// --------------------------- Utilities -----------------------------------

// Converts a double array into a PyList object for return to python.
static PyObject *c_double_array_to_list(double * mat,int cols){

    int j;

    PyObject *lst;

    lst = PyList_New(cols);
    if (!lst) return NULL;
    for (j = 0; j < cols; j++) {
        PyObject *num = PyFloat_FromDouble(mat[j]);
        if (!num) {
            Py_DECREF(lst);
            return NULL;
        }
        PyList_SET_ITEM(lst, j, num);   // reference to num stolen
    }
 
    return lst;
}

// Converts a int array into a PyList object for return to python.
static PyObject *c_int_array_to_list(int * mat,int cols){

    int j;

    PyObject *lst;

    lst = PyList_New(cols);
    if (!lst) return NULL;
    for (j = 0; j < cols; j++) {
        PyObject *num = PyInt_FromLong(mat[j]);
        if (!num) {
            Py_DECREF(lst);
            return NULL;
        }
        PyList_SET_ITEM(lst, j, num);   // reference to num stolen
    }
 
    return lst;
}

// ----------------------------------------------------------------------------

// If using python 2.7.3 or later, build with PyCapsules

// Delete capsule containing a quad tree - name of capsule must be exactly
// "quad tree".

#ifdef PYVERSION273


void delete_quad_tree_cap(PyObject * cap){
    quad_tree * kill = (quad_tree*) PyCapsule_GetPointer(cap,"quad tree");
    if(kill!=NULL){
        delete_quad_tree(kill);
    } else{
    }
}

// Delete capsule containing a sparse_dok - name of capsule must be exactly
// "sparse dok".
void delete_dok_cap(PyObject *cap){

    sparse_dok * kill = (sparse_dok*) PyCapsule_GetPointer(cap,"sparse dok");
    if(kill!=NULL){
        delete_dok_matrix(kill);
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

// Delete cobj containing a sparse_dok
void delete_dok_cobj(void *cobj){

    sparse_dok * kill = (sparse_dok*) cobj;
    if(kill!=NULL){
        delete_dok_matrix(kill);
    }

}

#endif

//----------------------- PYTHON WRAPPER FUNCTION -----------------------------

// Parses python information about triangles
// and vertex coordiantes to build a quad tree for efficient location of points in
// the triangle mesh.
PyObject *build_quad_tree(PyObject *self, PyObject *args) {



    PyArrayObject *triangles;
    PyArrayObject *vertex_coordinates;
    PyArrayObject *extents;


    // Convert Python arguments to C
    if (!PyArg_ParseTuple(args, "OOO",
                                            &triangles,
                                            &vertex_coordinates,
                                            &extents
                                            )) {
      PyErr_SetString(PyExc_RuntimeError,
              "fitsmooth.c: could not parse input");
      return NULL;
    }

    CHECK_C_CONTIG(triangles);
    CHECK_C_CONTIG(vertex_coordinates);
    CHECK_C_CONTIG(extents);

    int n = triangles->dimensions[0];

    #ifdef PYVERSION273
    
    return  PyCapsule_New((void*) _build_quad_tree(n,
                          (long*) triangles -> data,
                          (double*) vertex_coordinates -> data,
                          (double*) extents -> data),
                          "quad tree",
                          &delete_quad_tree_cap);

    #else

    return  PyCObject_FromVoidPtr((void*) _build_quad_tree(n,
                      (long*) triangles -> data,
                      (double*) vertex_coordinates -> data,
                      (double*) extents -> data),
                      &delete_quad_tree_cobj); 
    #endif
}

// Builds the smoothing matrix D. Parses 
// input mesh information from python and then builds a sparse_dok, returning
// a capsule object wrapping a pointer to this struct.
PyObject *build_smoothing_matrix(PyObject *self, PyObject *args) {


    PyArrayObject *triangles;
    PyArrayObject *areas;
    PyArrayObject *vertex_coordinates;
    int err;

    // Convert Python arguments to C
    if (!PyArg_ParseTuple(args, "OOO" ,
                                            &triangles,
                                            &areas,
                                            &vertex_coordinates
                                            )) {
      PyErr_SetString(PyExc_RuntimeError,
              "fitsmooth.c: could not parse input");
      return NULL;
    }

    CHECK_C_CONTIG(triangles);
    CHECK_C_CONTIG(areas);
    CHECK_C_CONTIG(vertex_coordinates);

    int n = triangles->dimensions[0];


    sparse_dok * smoothing_mat; // Should be an input argument?
    smoothing_mat = make_dok();

    err = _build_smoothing_matrix(n,
                      (long*) triangles  -> data,
                      (double*) areas -> data,
                      (double*) vertex_coordinates -> data,
                      (int *) vertex_coordinates -> strides,
                      smoothing_mat);

   
    if (err != 0) {
      PyErr_SetString(PyExc_RuntimeError,
              "Unknown Error");
      return NULL;
    }

    #ifdef PYVERSION273

    return  PyCapsule_New((void*) smoothing_mat,
                  "sparse dok",
                  &delete_dok_cap); 

    #else

    return  PyCObject_FromVoidPtr((void*) smoothing_mat,
                  &delete_dok_cobj); 
    
    #endif

}

// Builds AtA and Atz directly from a list of 
// points and values Returns a pointer to the sparse dok matrix AtA wrapped in 
// a capsule object and a python list for the array Atz.
PyObject *build_matrix_AtA_Atz_points(PyObject *self, PyObject *args) {


    PyArrayObject *triangles;
    PyArrayObject *point_coordinates;
    PyArrayObject *z;
    int N; // Number of triangles
    int err;
    int npts;
    int zdims;
    PyObject *tree;

    // Convert Python arguments to C
    if (!PyArg_ParseTuple(args, "OiOOOii",&tree, &N,
                                            &triangles,
                                            &point_coordinates,
                                            &z,
                                            &zdims,
                                            &npts
                                            )) {
      PyErr_SetString(PyExc_RuntimeError,
              "fitsmooth.c: could not parse input");
      return NULL;
    }

    
    CHECK_C_CONTIG(triangles);
    CHECK_C_CONTIG(point_coordinates);
    CHECK_C_CONTIG(z);

    #ifdef PYVERSION273
    quad_tree * quadtree = (quad_tree*) PyCapsule_GetPointer(tree,"quad tree");
    #else
    quad_tree * quadtree = (quad_tree*) PyCObject_AsVoidPtr(tree);
    #endif

    sparse_dok * dok_AtA; // Should be an input argument?
    dok_AtA = make_dok();
    double ** Atz = malloc(sizeof(double*)*zdims);
    int i;
    for(i=0;i<zdims;i++){
        Atz[i] = malloc(sizeof(double)*N);
    }

    err = _build_matrix_AtA_Atz_points(N,(long*) triangles->data,
                      (double*) point_coordinates->data,
                      (double*) z->data,
                      zdims,
                      npts,
                      dok_AtA,
                      Atz,
                      quadtree);


    if (err != 0) {
      PyErr_SetString(PyExc_RuntimeError,
              "Unknown Error");
      return NULL;
    }

    #ifdef PYVERSION273
    PyObject * AtA_cap =  PyCapsule_New((void*) dok_AtA,
                  "sparse dok",
                  &delete_dok_cap);
    #else
    PyObject * AtA_cap =  PyCObject_FromVoidPtr((void*) dok_AtA,
                  &delete_dok_cobj); 
    #endif
    PyObject * Atz_ret = PyList_New(zdims);
    for(i=0;i<zdims;i++){
        PyList_SET_ITEM(Atz_ret,i,c_double_array_to_list(Atz[i],N));
        free(Atz[i]);
    }
    free(Atz);

    PyObject *lst = PyList_New(2);
    PyList_SET_ITEM(lst, 0, AtA_cap);
    PyList_SET_ITEM(lst, 1, Atz_ret);
    return lst;

}

// Searches a quad tree struct for the triangle containing a given point,
// returns the sigma values produced by this point and the triangle found,
// and the triangle index. Found is returned as 0 if no triangle is found
// and 1 otherwise.
PyObject *individual_tree_search(PyObject *self, PyObject *args) {

    // Setting up variables to parse input
    PyObject *tree;
    PyArrayObject *point;

    // Convert Python arguments to C
    if (!PyArg_ParseTuple(args, "OO",&tree, &point
                                            )) {
      PyErr_SetString(PyExc_RuntimeError,
              "fitsmooth.c: could not parse input");
      return NULL;
    }

    CHECK_C_CONTIG(point);

    #ifdef PYVERSION273
    quad_tree * quadtree = (quad_tree*) PyCapsule_GetPointer(tree,"quad tree");
    #else
    quad_tree * quadtree = (quad_tree*) PyCObject_AsVoidPtr(tree);
    #endif

    double xp,yp;
    if(PyArray_TYPE(point)==7){
        int *pointd = (int*)point->data;
        xp = (double)(pointd[0]);
        yp = (double)(pointd[1]);
    }
    else{
        double *pointd = (double*)point->data;
        xp = pointd[0];
        yp = pointd[1];
    }

    triangle * T = search(quadtree,xp,yp);
    PyObject * sigmalist = PyList_New(3);
    long found;
    long index;
    
    if(T!=NULL){
            double * sigma = calculate_sigma(T,xp,yp);
            sigmalist = c_double_array_to_list(sigma,3);
            free(sigma);
            found = 1;
            index = (long)T->index;
    }else{
            double sigma[3];
            sigma[0]=sigma[1]=sigma[2]=-1;
            sigmalist = c_double_array_to_list(sigma,3);
            index = -10;
            found = 0;
    }
    PyObject * retlist = PyList_New(3);
    PyList_SET_ITEM(retlist,0,PyInt_FromLong(found));
    PyList_SET_ITEM(retlist,1,sigmalist);
    PyList_SET_ITEM(retlist,2,PyInt_FromLong(index));
    return retlist;
}

// Returns the total number of triangles stored in quad_tree.
// Takes a capsule object holding a pointer to the quad_tree as input.
//
// PADARN NOTE: This function was only added for the purpose of passing unit
// tests. Perhaps this, and other tree functions, should be moved to another 
// file to provide a more comprehensive python interface to the tree structure. 
PyObject *items_in_tree(PyObject *self, PyObject *args) {

    // Setting up variables to parse input
    PyObject *tree; // capsule to hold quad_tree pointer

    // Convert Python arguments to C
    if (!PyArg_ParseTuple(args, "O",&tree
                                            )) {
      PyErr_SetString(PyExc_RuntimeError,
              "fitsmooth.c: could not parse input");
      return NULL;
    }

    // Extract quad_tree pointer from capsule
    #ifdef PYVERSION273
    quad_tree * quadtree = (quad_tree*) PyCapsule_GetPointer(tree,"quad tree");
    #else
    quad_tree * quadtree = (quad_tree*) PyCObject_AsVoidPtr(tree);
    #endif
    
    // Return the number of elements in the tree (stored in struct)
    return PyInt_FromLong((long)quadtree->count);
    
}

// Returns the total number of nodes stored in quad_tree.
// Takes a capsule object holding a pointer to the quad_tree as input. 
PyObject *nodes_in_tree(PyObject *self, PyObject *args) {

    // Setting up variables to parse input
    PyObject *tree; // capsule to hold quad_tree pointer

    // Convert Python arguments to C
    if (!PyArg_ParseTuple(args, "O",&tree
                                            )) {
      PyErr_SetString(PyExc_RuntimeError,
              "fitsmooth.c: could not parse input");
      return NULL;
    }

    // Extract quad_tree pointer from capsule
    #ifdef PYVERSION273
    quad_tree * quadtree = (quad_tree*) PyCapsule_GetPointer(tree,"quad tree");
    #else
    quad_tree * quadtree = (quad_tree*) PyCObject_AsVoidPtr(tree);
    #endif
    
    // Return the number of elements in the tree (stored in struct)
    return PyInt_FromLong((long)quad_tree_node_count(quadtree));
    
}

// Combines two sparse_dok and two double arrays together
// which represent partial completion of the AtA and Atz matrices, when they are build
// in parts due to points being read in blocks in python. Result is stored in the first
// sparse_dok and double array.
// Function takes as arguments two capsule objects holding pointers to the sparse_dok
// structs, and two numpy double arrays. Aslo the size of the double arrays, and the number
// of columns (variables) for the array Atz are passed as arguments.
//
// PADARN NOTE: Blocking the points in python is far slower than reading them directly
// in c and bypassing the need for this.
PyObject *combine_partial_AtA_Atz(PyObject *self, PyObject *args) {

    // Setting up variables to parse input
    PyObject *AtA_cap1, *AtA_cap2;
    PyArrayObject *Atz1, *Atz2;
    int n,zdim; // number of nodes

    // Convert Python arguments to C
    if (!PyArg_ParseTuple(args, "OOOOii",&AtA_cap1, &AtA_cap2,
                                        &Atz1, &Atz2, &zdim, &n
                                            )) {
      PyErr_SetString(PyExc_RuntimeError,
              "fitsmooth.c: could not parse input");
      return NULL;
    }

    // Get pointers to sparse_dok objects from capsules
    #ifdef PYVERSION273
    sparse_dok * dok_AtA1 = (sparse_dok*) PyCapsule_GetPointer(AtA_cap1,"sparse dok");
    sparse_dok * dok_AtA2 = (sparse_dok*) PyCapsule_GetPointer(AtA_cap2,"sparse dok");
    #else
    sparse_dok * dok_AtA1 = (sparse_dok*) PyCObject_AsVoidPtr(AtA_cap1);
    sparse_dok * dok_AtA2 = (sparse_dok*) PyCObject_AsVoidPtr(AtA_cap2);
    #endif

    // Combine the partial AtA and Atz
    _combine_partial_AtA_Atz(dok_AtA1,dok_AtA2,
                             (double*) Atz1->data,
                             (double*) Atz2->data,
                             n, zdim);

    // Return nothing interesting
    return Py_BuildValue("");
}

// Converts a sparse_dok matrix to a full non-compressed matrix expressed
// as a list of lists (python). Takes as input a capsule object containing a pointer to the
// sparse_dok object. Also takes an integer n as input, specifying the (n x n) size of the 
// matrix. Matrix is assumed square.
//
// PADARN NOTE: This function does not seem particularly sensible, but was required for the 
// purposes of passing unit tests. 
PyObject *return_full_D(PyObject *self, PyObject *args) {

    // Setting up variables to parse input
    PyObject *D_cap; // capsule object holding sparse_dok pointer
    int n; // number of matrix columns/rows

    // Convert Python arguments to C
    if (!PyArg_ParseTuple(args, "Oi",&D_cap, &n
                                            )) {
      PyErr_SetString(PyExc_RuntimeError,
              "fitsmooth.return_full_D: could not parse input");
      return NULL;
    }

    // Get pointer to spars_dok struct
    #ifdef PYVERSION273
    sparse_dok * D_mat = (sparse_dok*) PyCapsule_GetPointer(D_cap,"sparse dok");
    #else
    sparse_dok * D_mat = (sparse_dok*) PyCObject_AsVoidPtr(D_cap);
    #endif

    // Check to make sure that the specified size of the matrix is at least as big
    // as the entries stored in the sparse_dok.
    if (D_mat->num_rows>n || get_dok_rows(D_mat)>n){
        PyErr_SetString(PyExc_RuntimeError,
              "fitsmooth.return_full_D: sparse_dok is bigger than size specified for return.");
        return NULL;
    }

    // Build new python list containing the full matrix to return
    PyObject *ret_D = PyList_New(n);
    int i,j;
    edge_key_t key;
    edge_t *s;
    for(i=0;i<n;i++)
    {
        PyObject *temp = PyList_New(n);
        for(j=0;j<n;j++)
        {
            key.i=i;
            key.j=j;
            s = find_dok_entry(D_mat,key);
            if (s){
                PyList_SET_ITEM(temp,j,PyFloat_FromDouble(s->entry));
            }else{
                PyList_SET_ITEM(temp,j,PyFloat_FromDouble(0));
            }
        PyList_SET_ITEM(ret_D,i,temp);
        }
    }

    return ret_D;
}

// Takes as input two capsule objects corresponding to the two matricies
// D and AtA, along with a regularization coefficient alpha.
// Returns three matricies corresponding to the CSR format storage of the matrix resulting
// from the computation B = AtA + alpha * D.
//
// Capsule objects are not freed and are left alone to be freed as they naturally go out of 
// scope in Python. The temporary CSR matrix created is cleaned up.
PyObject *build_matrix_B(PyObject *self, PyObject *args) {


    // Setting up variables to parse input
    double alpha; // Regularization coefficient
    PyObject *smoothing_mat_cap; // Capsule storing D pointer (sparse_dok struct)
    PyObject *AtA_cap; // Capsule storing AtA pointer (sparse_dok struct)

    // Convert Python arguments to C
    if (!PyArg_ParseTuple(args, "OOd",&smoothing_mat_cap, &AtA_cap, &alpha
                                            )) {
      PyErr_SetString(PyExc_RuntimeError,
              "fitsmooth.build_matrix_B: could not parse input");
      return NULL;
    }

    // Extract pointers to c structs from capsules
    #ifdef PYVERSION273
    sparse_dok * smoothing_mat = (sparse_dok*) PyCapsule_GetPointer(smoothing_mat_cap,"sparse dok");
    sparse_dok * dok_AtA = (sparse_dok*) PyCapsule_GetPointer(AtA_cap,"sparse dok");
    #else
    sparse_dok * smoothing_mat = (sparse_dok*) PyCObject_AsVoidPtr(smoothing_mat_cap);
    sparse_dok * dok_AtA = (sparse_dok*) PyCObject_AsVoidPtr(AtA_cap);
    #endif

    // Add two sparse_dok matrices
    add_sparse_dok(smoothing_mat,alpha,dok_AtA,1);
    
    // Create sparse_csr matrix and convert result to this format
    sparse_csr * B;
    B = make_csr();
    convert_to_csr_ptr(B,smoothing_mat);
    
    // Extract the sparse_csr data to be returned as python lists
    PyObject *data = c_double_array_to_list(B->data,
                                B->num_entries);
    PyObject *colind = c_int_array_to_list(B->colind,
                                B->num_entries);
    PyObject *row_ptr = c_int_array_to_list(B->row_ptr,
                                B->num_rows);

    // Build python list for return
    PyObject *lst = PyList_New(3);
    PyList_SET_ITEM(lst, 0, data);
    PyList_SET_ITEM(lst, 1, colind);
    PyList_SET_ITEM(lst, 2, row_ptr);

    // Clean up
    delete_csr_matrix(B);
    
    return lst;
   


}
//------------------------------------------------------------------------------


// ------------------------------ PYTHON GLUE ----------------------------------

//==============================================================================
// Structures to allow calling from python
//==============================================================================

// Method table for python module
static struct PyMethodDef MethodTable[] = {
    {"items_in_tree",items_in_tree, METH_VARARGS, "Print out"},
    {"nodes_in_tree",nodes_in_tree, METH_VARARGS, "Print out"},
    {"return_full_D",return_full_D, METH_VARARGS, "Print out"},
	{"build_matrix_B",build_matrix_B, METH_VARARGS, "Print out"},
    {"build_quad_tree",build_quad_tree, METH_VARARGS, "Print out"},
    {"build_smoothing_matrix",build_smoothing_matrix, METH_VARARGS, "Print out"},
    {"build_matrix_AtA_Atz_points",build_matrix_AtA_Atz_points, METH_VARARGS, "Print out"},
    {"combine_partial_AtA_Atz",combine_partial_AtA_Atz, METH_VARARGS, "Print out"},
    {"individual_tree_search",individual_tree_search, METH_VARARGS, "Print out"},
	{NULL, NULL, 0, NULL}   // sentinel
};


// Module initialisation
void initfitsmooth(void){
  Py_InitModule("fitsmooth", MethodTable);
  import_array(); // Necessary for handling of NumPY structures

}

// --------------------------------------------------------------------------------
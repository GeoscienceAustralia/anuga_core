


#include <stdio.h>   /* gets */
#include <stdlib.h>  /* atoi, malloc */
#include <string.h>  /* strcpy */
#include "math.h"

#include "sparse_dok.h" /* in utilities */
#include "quad_tree.h"  /* in utilities */
#include "anuga_typedefs.h" /* in utilities */

#include "omp.h"

// Errors defined for netcdf reading
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

//-------------------------- QUANTITY FITTING ------------------------------


// Builds the matrix D used to smooth the interpolation 
// of a variables from scattered data points to a mesh. See fit.py for more details.s
anuga_int _build_smoothing_matrix(anuga_int n,
                      anuga_int* triangles,
        		          double* areas,
                      double* vertex_coordinates,
                      anuga_int* strides,
                      sparse_dok * smoothing_mat){


    anuga_int k;
    anuga_int k3,k6;
    anuga_int err = 0;
    edge_key_t key;

    double det,area,x0,x1,x2,y0,y1,y2;
    double a0,b0,a1,b1,a2,b2,e01,e12,e20;
    anuga_int v0,v1,v2;
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
        if (add_dok_entry(smoothing_mat,key,smoothing_val) != 0) return -1;

        // v1,v1
        key.i = v1;
        key.j = v1;
        smoothing_val = (a1*a1 + b1*b1)*area;
        if (add_dok_entry(smoothing_mat,key,smoothing_val) != 0) return -1;

        // v2,v2
        key.i = v2;
        key.j = v2;
        smoothing_val = (a2*a2 + b2*b2)*area;
        if (add_dok_entry(smoothing_mat,key,smoothing_val) != 0) return -1;


        // insert off diagonal contributions
        e01 = (a0*a1 + b0*b1)*area;
        // v0,v1 (v1,v0)
        key.i = v0;
        key.j = v1;
        if (add_dok_entry(smoothing_mat,key,e01) != 0) return -1;
        key.i = v1;
        key.j = v0;
        if (add_dok_entry(smoothing_mat,key,e01) != 0) return -1;

        e12 = (a1*a2 + b1*b2)*area;
        // v1,v2 (v2,v1)
        key.i = v1;
        key.j = v2;
        if (add_dok_entry(smoothing_mat,key,e12) != 0) return -1;
        key.i = v2;
        key.j = v1;
        if (add_dok_entry(smoothing_mat,key,e12) != 0) return -1;

        e20 = (a2*a0 + b2*b0)*area;
        // v2,v0 (v0,v2)
        key.i = v2;
        key.j = v0;
        if (add_dok_entry(smoothing_mat,key,e20) != 0) return -1;
        key.i = v0;
        key.j = v2;
        if (add_dok_entry(smoothing_mat,key,e20) != 0) return -1;
    }

    return err;

}

// Builds a quad tree out of a list of triangles for quick 
// searching. 
quad_tree * _build_quad_tree(anuga_int n,
                      anuga_int* triangles,
                      double* vertex_coordinates,
                      double* extents)               
{   
    
    anuga_int k,k6;
    double x0,y0,x1,y1,x2,y2;

    // set up quad tree and allocate memory
    quad_tree * tree = new_quad_tree(extents[0],extents[1],extents[2],extents[3]);
    if (!tree) return NULL;
    
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
        if (!T) {
            delete_quad_tree(tree);
            return NULL;
        }
        if (quad_tree_insert_triangle(tree,T) != 0) {
            free(T);   // not owned by the tree when insert fails
            delete_quad_tree(tree);
            return NULL;
        }
    }
  
    // return pointer to new tree struct
    return tree;
    
}

// Builds the AtA and Atz interpolation matrix
// and residual. Uses a quad_tree for fast access to the triangles of the mesh.
// This function takes a list of point coordinates, and associated point values
// (for any number of attributes).
anuga_int _build_matrix_AtA_Atz_points(anuga_int N, 
    anuga_int * triangles,
    double * point_coordinates, 
    double * point_values,
    anuga_int zdims,
    anuga_int npts,
    sparse_dok * AtA,
    double ** Atz,
    quad_tree * quadtree)
              {

    anuga_int k;
    anuga_int i,w;
    anuga_int err = 0;   // set inside the critical section only
    
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
            double sigma[3];
            calculate_sigma(T,x,y,sigma);
            anuga_int js[3];
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

                   if (add_dok_entry(AtA,key,sigma[i]*sigma[w]) != 0) err = -1;
                }                        
            }
            }

       } 
    }

    return err;
}

// Combines two sparse_dok matricies and two vectors of doubles. 
anuga_int _combine_partial_AtA_Atz(sparse_dok * dok_AtA1,sparse_dok * dok_AtA2,
                             double* Atz1,
                             double* Atz2,
                             anuga_int n, anuga_int zdim){

    if (add_sparse_dok(dok_AtA1,1,dok_AtA2,1) != 0) return -1;

    anuga_int i;
    for(i=0;i<n*zdim;i++){
        Atz1[i]+=Atz2[i];
    }
    return 0;
}

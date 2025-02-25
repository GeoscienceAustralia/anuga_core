


#include <stdio.h>   /* gets */
#include <stdlib.h>  /* atoi, malloc */
#include <string.h>  /* strcpy */
#include "math.h"
#include "stdint.h"

#include "sparse_dok.h" /* in utilities */
#include "quad_tree.h"  /* in utilities */

#if defined(__APPLE__)
   // clang doesn't have openmp
#else
   #include "omp.h"
#endif

// Errors defined for netcdf reading
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

//-------------------------- QUANTITY FITTING ------------------------------


// Builds the matrix D used to smooth the interpolation 
// of a variables from scattered data points to a mesh. See fit.py for more details.s
int64_t _build_smoothing_matrix(int64_t n,
                      int64_t * triangles,
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
quad_tree * _build_quad_tree(int64_t n,
                      int64_t * triangles,
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
int64_t _build_matrix_AtA_Atz_points(int64_t N, int64_t  * triangles,
                      double * point_coordinates, double * point_values,
                      int64_t zdims, int64_t npts,
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
            int64_t js[3];
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
void _combine_partial_AtA_Atz(sparse_dok * dok_AtA1, sparse_dok * dok_AtA2,
                             double* Atz1,
                             double* Atz2,
                             int n, int zdim){

    add_sparse_dok(dok_AtA1,1,dok_AtA2,1);

    int i;
    for(i=0;i<n*zdim;i++){
        Atz1[i]+=Atz2[i];
    }
}

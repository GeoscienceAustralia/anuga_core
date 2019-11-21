// Python - C extension module for advection.py
//
// To compile (Python2.X):
// use python ../utilities/compile.py to
// compile all C files within a directory
//
//
// Steve Roberts, ANU 2008



#include "math.h"
#include "stdio.h"

//-------------------------------------------
// Low level routines (called from wrappers)
//------------------------------------------

double _compute_fluxes(
		    double* quantity_update,
		    double* quantity_edge,
		    double* quantity_bdry,
                    long*   domain_neighbours,
		    long*   domain_neighbour_edges,
		    double* domain_normals,
                    double* domain_areas,
		    double* domain_radii,
		    double* domain_edgelengths,
		    long*   domain_tri_full_flag,
		    double* domain_velocity,
                    double  huge_timestep,
                    double  max_timestep,
		    int ntri,
		    int nbdry){

 
        //Local Variables

        double qr,ql;
        double normal[2];
        double normal_velocity;
        double flux, edgeflux;
        double max_speed;
        double optimal_timestep;
	double timestep;
	int  k_i,n_m,k_i_j;
	int k,i,j,n,m;
	int k3;
	
	//Loop through triangles

	timestep = max_timestep;

        for (k=0; k<ntri; k++){
            optimal_timestep = huge_timestep;
            flux = 0.0;
	    k3 = 3*k;
            for (i=0; i<3; i++){
	        k_i = k3+i;
                //Quantities inside triangle facing neighbour i
                ql = quantity_edge[k_i];


                //Quantities at neighbour on nearest face
                n = domain_neighbours[k_i];
                if (n < 0) {
                    m = -n-1; //Convert neg flag to index
                    qr = quantity_bdry[m];
                } else {
                    m = domain_neighbour_edges[k_i];
		    n_m = 3*n+m;
                    qr = quantity_edge[n_m];
                }


                //Outward pointing normal vector
                for (j=0; j<2; j++){
		    k_i_j = 6*k+2*i+j;
                    normal[j] = domain_normals[k_i_j];
                }


                //Flux computation using provided function
                normal_velocity = domain_velocity[0]*normal[0] + domain_velocity[1]*normal[1];

                if (normal_velocity < 0) {
                    edgeflux = qr * normal_velocity;
                } else {
                    edgeflux = ql * normal_velocity;
                }

                max_speed = fabs(normal_velocity);
                flux = flux - edgeflux * domain_edgelengths[k_i];

                //Update optimal_timestep
                if (domain_tri_full_flag[k] == 1) {
                    if (max_speed != 0.0) {
                        optimal_timestep = (optimal_timestep>domain_radii[k]/max_speed) ? 
			  domain_radii[k]/max_speed : optimal_timestep;
                    }
                }

            }

            //Normalise by area and store for when all conserved
            //quantities get updated
            quantity_update[k] = flux/domain_areas[k];

            timestep = (timestep>optimal_timestep) ? optimal_timestep : timestep;

        }

	return timestep;
}

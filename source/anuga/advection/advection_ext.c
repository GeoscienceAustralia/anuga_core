#include "math.h"
#include "stdio.h"

double compute_fluxes(double* stage_edge,
		    double* stage_bdry,
		    double* stage_update,
                    int* neighbours,
		    int* neighbour_edges,
		    double* normals,
                    double* areas,
		    double* radii,
		    double* edgelengths,
		    int*    tri_full_flag,
                    double  huge_timestep,
                    double  max_timestep,
		    double* v,
		    int N){
        //Loop

        double qr,ql;
        double normal[2];
        double normal_velocity;
        double flux, edgeflux;
        double max_speed;
        double optimal_timestep;
	double timestep;
	int k_i,n_m,k_i_j;
	int k,i,j,n,m;
	
	timestep = max_timestep;

	//printf("N = %i\n",N);
	//printf("timestep = %g\n",timestep);
	//printf("huge_timestep = %g\n",huge_timestep);

        for (k=0; k<N; k++){
            optimal_timestep = huge_timestep;
            flux = 0.0;  //Reset work array
            for (i=0; i<3; i++){
	        k_i = 3*k+i;
                //Quantities inside volume facing neighbour i
                ql = stage_edge[k_i];
		//printf("stage_edge[%i] = %g\n",k_i,stage_edge[k_i]);

                //Quantities at neighbour on nearest face
                n = neighbours[k_i];
                if (n < 0) {
                    m = -n-1; //Convert neg flag to index
                    qr = stage_bdry[m];
                } else {
                    m = neighbour_edges[k_i];
		    n_m = 3*n+m;
                    qr = stage_edge[n_m];
                }


                //Outward pointing normal vector
                for (j=0; j<2; j++){
		    k_i_j = 6*k+2*i+j;
                    normal[j] = normals[k_i_j];
                }


                //Flux computation using provided function
                normal_velocity = v[0]*normal[0] + v[1]*normal[1];

                if (normal_velocity < 0) {
                    edgeflux = qr * normal_velocity;
                } else {
                    edgeflux = ql * normal_velocity;
                }

                max_speed = fabs(normal_velocity);
                flux = flux - edgeflux * edgelengths[k_i];

                //Update optimal_timestep
                if (tri_full_flag[k] == 1) {
                    if (max_speed != 0.0) {
                        optimal_timestep = (optimal_timestep>radii[k]/max_speed) ? radii[k]/max_speed : optimal_timestep;
                    }
                }

            }

            //Normalise by area and store for when all conserved
            //quantities get updated
            stage_update[k] = flux/areas[k];

            timestep = (timestep>optimal_timestep) ? optimal_timestep : timestep;

        }

	//printf("timestep out = %g\n",timestep);

	return timestep;
}


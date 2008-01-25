#include "math.h"
#include "stdio.h"

void  print_double_array(char* name, double* array, int n, int m){

    int k,i,km;

    printf("%s = [",name);
    for (k=0; k<n; k++){
	km = m*k;
	printf("[");
	for (i=0; i<m ; i++){
	    printf("%g ",array[km+i]);
	}
	if (k==(n-1))
	    printf("]");
	else
	    printf("]\n");
    }
    printf("]\n");
}

void  print_int_array(char* name, int* array, int n, int m){

    int k,i,km;

    printf("%s = [",name);
    for (k=0; k<n; k++){
	km = m*k;
	printf("[");
	for (i=0; i<m ; i++){
	    printf("%i ",array[km+i]);
	}
	if (k==(n-1))
	    printf("]");
	else
	    printf("]\n");
    }
    printf("]\n");
}


double compute_fluxes(
		    double* stage_update,
		    double* stage_edge,
		    double* stage_bdry,
                    int* neighbours,
		    int* neighbour_edges,
		    double* normals,
                    double* areas,
		    double* radii,
		    double* edgelengths,
		    int*   tri_full_flag,
                    double  huge_timestep,
                    double  max_timestep,
		    double* v,
		    int ntri,
		    int nbdry){
        //Loop

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
	
	timestep = max_timestep;


	printf("======================================================\n");
	printf("INSIDE compute_fluxes\n");


        print_double_array( "stage_update",stage_update, ntri, 1);
        print_double_array( "stage_edge",stage_edge, ntri, 3);
        print_double_array( "stage_bdry",stage_bdry, nbdry, 1);
        print_int_array( "neighbours",neighbours, ntri, 3);
        print_int_array( "neighbour_edges",neighbour_edges, ntri, 3);
        print_double_array( "normals",normals, ntri, 6);
        print_double_array( "areas",areas, ntri, 1);
        print_double_array( "radii",radii, ntri, 1);
        print_double_array( "edgelengths",edgelengths, ntri, 3);
        print_int_array( "tri_full_flag",tri_full_flag, ntri, 1);
	printf("huge_timestep = %g\n",huge_timestep);
	printf("max_timestep = %g\n",max_timestep);
        print_double_array( "v",v, 2, 1);
	printf("ntri = %i\n",ntri);
	printf("nbdry = %i\n",nbdry);


        for (k=0; k<ntri; k++){
            optimal_timestep = huge_timestep;
            flux = 0.0;
	    k3 = 3*k;
            for (i=0; i<3; i++){
	        k_i = k3+i;
                //Quantities inside volume facing neighbour i
                ql = stage_edge[k_i];


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



	printf("INSIDE compute_fluxes, end \n");

        print_double_array( "stage_update",stage_update, ntri, 1);

	printf("FINISHED compute_fluxes\n");
	printf("======================================================\n");


	return timestep;
}


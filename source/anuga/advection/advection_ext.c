// Python - C extension module for advection.py
//
// To compile (Python2.X):
// use python ../utilities/compile.py to
// compile all C files within a directory
//
//
// Steve Roberts, ANU 2008


#include "Python.h"
#include "numpy/arrayobject.h"
#include "math.h"
#include "stdio.h"

// Shared code snippets
#include "util_ext.h"


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


//-----------------------------------------------------
// Python method Wrappers 
//-----------------------------------------------------



PyObject *compute_fluxes(PyObject *self, PyObject *args) {
  /*Compute all fluxes and the timestep suitable for all volumes
    in domain.

    Fluxes across each edge are scaled by edgelengths and summed up
    Resulting flux is then scaled by area and stored in
    explicit_update for the conserved quantity stage.

    The maximal allowable speed computed for each volume
    is converted to a timestep that must not be exceeded. The minimum of
    those is computed as the next overall timestep.

    Python call:
    timestep = advection_ext.compute_fluxes(domain,quantity,huge_timestep,max_timestep)

    Post conditions:
      domain.explicit_update is reset to computed flux values
      domain.timestep is set to the largest step satisfying all volumes.

  */

  PyObject *domain, *quantity;
 
  PyArrayObject 
    * quantity_update,
    * quantity_edge,
    * quantity_bdry,
    * domain_neighbours,
    * domain_neighbour_edges,
    * domain_normals,
    * domain_areas,
    * domain_radii,
    * domain_edgelengths,
    * domain_tri_full_flag,
    * domain_velocity;
                  

  // Local variables
  int ntri, nbdry;
  double huge_timestep, max_timestep;
  double timestep;


  // Convert Python arguments to C
  if (!PyArg_ParseTuple(args, "OOdd",
			&domain,
			&quantity,
			&huge_timestep,
			&max_timestep)) {
    PyErr_SetString(PyExc_RuntimeError, "advection_ext.c: compute_fluxes could not parse input");
    return NULL;
  }

  quantity_edge          = get_consecutive_array(quantity, "edge_values");
  quantity_bdry          = get_consecutive_array(quantity, "boundary_values");
  quantity_update        = get_consecutive_array(quantity, "explicit_update");
  domain_neighbours      = get_consecutive_array(domain,   "neighbours");
  domain_neighbour_edges = get_consecutive_array(domain,   "neighbour_edges");
  domain_normals         = get_consecutive_array(domain,   "normals");
  domain_areas           = get_consecutive_array(domain,   "areas");
  domain_radii           = get_consecutive_array(domain,   "radii");
  domain_edgelengths     = get_consecutive_array(domain,   "edgelengths");
  domain_tri_full_flag   = get_consecutive_array(domain,   "tri_full_flag");
  domain_velocity        = get_consecutive_array(domain,   "velocity");  
    
  ntri  = quantity_edge -> dimensions[0];
  nbdry = quantity_bdry -> dimensions[0]; 

  timestep = _compute_fluxes((double*) quantity_update -> data,
			     (double*) quantity_edge -> data,
			     (double*) quantity_bdry -> data,
			     (long*)   domain_neighbours -> data,
			     (long*)   domain_neighbour_edges -> data,
			     (double*) domain_normals -> data,
			     (double*) domain_areas ->data,
			     (double*) domain_radii -> data,
			     (double*) domain_edgelengths -> data,
			     (long*)   domain_tri_full_flag -> data,
			     (double*) domain_velocity -> data,
			     huge_timestep,
			     max_timestep,
			     ntri,
			     nbdry);

  // Release and return
  Py_DECREF(quantity_update);
  Py_DECREF(quantity_edge);
  Py_DECREF(quantity_bdry);
  Py_DECREF(domain_neighbours);
  Py_DECREF(domain_neighbour_edges);
  Py_DECREF(domain_normals);
  Py_DECREF(domain_areas);
  Py_DECREF(domain_radii);
  Py_DECREF(domain_edgelengths);
  Py_DECREF(domain_tri_full_flag);
  Py_DECREF(domain_velocity);


  return Py_BuildValue("d", timestep);
}



//-------------------------------
// Method table for python module
//-------------------------------
static struct PyMethodDef MethodTable[] = {
  /* The cast of the function is necessary since PyCFunction values
   * only take two PyObject* parameters, and rotate() takes
   * three.
   */

  {"compute_fluxes", compute_fluxes, METH_VARARGS, "Print out"},
  {NULL, NULL}
};

// Module initialisation
void initadvection_ext(void){
  Py_InitModule("advection_ext", MethodTable);
  import_array(); // Necessary for handling of NumPY structures
}

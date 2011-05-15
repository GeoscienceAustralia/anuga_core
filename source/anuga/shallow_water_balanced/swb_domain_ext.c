// Python - C extension module for shallow_water.py
//
// To compile (Python2.6):
//  gcc -c swb_domain_ext.c -I/usr/include/python2.6 -o domain_ext.o -Wall -O
//  gcc -shared swb_domain_ext.o  -o swb_domain_ext.so
//
// or use python compile.py
//
// See the module swb_domain.py for more documentation on 
// how to use this module
//
//
// Stephen Roberts, ANU 2009
// Ole Nielsen, GA 2004


#include "Python.h"
#include "numpy/arrayobject.h"
#include "math.h"
#include <stdio.h>
//#include "numpy_shim.h"

// Shared code snippets
#include "util_ext.h"


const double pi = 3.14159265358979;

// Computational function for rotation
int _rotate(double *q, double n1, double n2) {
  /*Rotate the last  2 coordinates of q (q[1], q[2])
    from x,y coordinates to coordinates based on normal vector (n1, n2).

    Result is returned in array 2x1 r
    To rotate in opposite direction, call rotate with (q, n1, -n2)

    Contents of q are changed by this function */


  double q1, q2;

  // Shorthands
  q1 = q[1];  // x coordinate
  q2 = q[2];  // y coordinate

  // Rotate
  q[1] =  n1*q1 + n2*q2;
  q[2] = -n2*q1 + n1*q2;

  return 0;
}

// Innermost flux function
int _flux_function(double *q_inside, double *q_outside,
		   double n1, double n2, 
		   double g,
		   double *edgeflux, double *local_max_speed) 
{

  /*Compute fluxes between volumes for the shallow water wave equation
    cast in terms of the 'stage', w = h+z using
    the 'central scheme' as described in

    Kurganov, Noelle, Petrova. 'Semidiscrete Central-Upwind Schemes For
    Hyperbolic Conservation Laws and Hamilton-Jacobi Equations'.
    Siam J. Sci. Comput. Vol. 23, No. 3, pp. 707-740.

    The implemented formula is given in equation (3.15) on page 714
  */

  int i;

  double z_star, h_inside_star, h_outside_star;
  double w_inside,  h_inside,  u_inside,  v_inside,  z_inside;
  double w_outside, h_outside, u_outside, v_outside, z_outside;
  double s_min, s_max, soundspeed_inside, soundspeed_outside;
  double epsilon, denom, inverse_denominator;

  // Workspace
  double cons_inside[3], cons_outside[3], flux_outside[3], flux_inside[3];

  epsilon = 1.0e-15;

  // Setup aliases
  w_inside = q_inside[0]; 
  h_inside = q_inside[3];
  z_inside = q_inside[4];
  u_inside = q_inside[5];
  v_inside = q_inside[6];

  w_outside = q_outside[0]; 
  h_outside = q_outside[3];
  z_outside = q_outside[4];
  u_outside = q_outside[5];
  v_outside = q_outside[6];


  // Rotate the velocity terms and update. Setup cons_ to hold conservative rotated
  // variables 
  cons_inside[0] = w_inside;
  cons_inside[1] = u_inside;
  cons_inside[2] = v_inside;

  _rotate(cons_inside, n1, n2);

  u_inside = cons_inside[1];
  v_inside = cons_inside[2];
  cons_inside[1] = u_inside*h_inside;
  cons_inside[2] = v_inside*h_inside;

  cons_outside[0] = w_outside;
  cons_outside[1] = u_outside;
  cons_outside[2] = v_outside;

  _rotate(cons_outside, n1, n2);

  u_outside = cons_outside[1];
  v_outside = cons_outside[2];
  cons_outside[1] = u_outside*h_outside;
  cons_outside[2] = v_outside*h_outside;



  // Deal with discontinuous z
  z_star = max(z_inside, z_outside);
  h_inside_star  = max(0.0, w_inside-z_star);
  h_outside_star = max(0.0, w_outside-z_star);
	           

  // Maximal and minimal wave speeds
  soundspeed_inside  = sqrt(g*h_inside_star);
  soundspeed_outside = sqrt(g*h_outside_star);  
  
  s_max = max(u_inside + soundspeed_inside, u_outside + soundspeed_outside);
  s_max = max(0.0, s_max);

  s_min = min(u_inside - soundspeed_inside, u_outside - soundspeed_outside);
  s_min = min(0.0, s_min);
  
  // Flux formulas
  flux_inside[0] = u_inside*h_inside_star;
  flux_inside[1] = u_inside*u_inside*h_inside_star + 0.5*g*h_inside_star*h_inside_star;
  flux_inside[2] = u_inside*v_inside*h_inside_star;

  flux_outside[0] = u_outside*h_outside_star;
  flux_outside[1] = u_outside*u_outside*h_outside_star + 0.5*g*h_outside_star*h_outside_star;
  flux_outside[2] = u_outside*v_outside*h_outside_star;

  // Flux computation
  denom = s_max - s_min;

  edgeflux[0] = 0.0;
  edgeflux[1] = 0.0;
  edgeflux[2] = 0.0;
  
  if (denom < epsilon) 
    {     
      *local_max_speed = 0.0;
      //printf("here %g\n",h_inside);
      // Add in edge pressure term due to discontinuous bed
      if ( h_inside > 0.0 ) edgeflux[1] = 0.5*g*h_inside*h_inside ;
    } 
  else 
    {
      inverse_denominator = 1.0/denom;
      for (i = 0; i < 3; i++) 
	{
	  edgeflux[i] = s_max*flux_inside[i] - s_min*flux_outside[i];
	  edgeflux[i] += s_max*s_min*(cons_outside[i] - cons_inside[i]);
	  edgeflux[i] *= inverse_denominator;
	}

      // Balance the pressure term
      // FIXME SR: I think we need to add a term which uses simpson's rule.
      edgeflux[1] += 0.5*g*h_inside*h_inside - 0.5*g*h_inside_star*h_inside_star;
      
    }
  
// Maximal wavespeed
  *local_max_speed = max(fabs(s_max), fabs(s_min));
  //printf("local speed % g  h_inside %g \n",*local_max_speed, h_inside);
  // Rotate back
  _rotate(edgeflux, n1, -n2);
  
  return 0;
}



//=========================================================================
// Python Glue
//=========================================================================


//========================================================================
// Compute fluxes
//========================================================================

PyObject *compute_fluxes(PyObject *self, PyObject *args) {
  /*Compute all fluxes and the timestep suitable for all volumes
    in domain.

    Compute total flux for each conserved quantity using "flux_function_central"

    Fluxes across each edge are scaled by edgelengths and summed up
    Resulting flux is then scaled by area and stored in
    explicit_update for each of the three conserved quantities
    stage, xmomentum and ymomentum

    The maximal allowable speed computed by the flux_function for each volume
    is converted to a timestep that must not be exceeded. The minimum of
    those is computed as the next overall timestep.

    Python call:
    timestep = compute_fluxes(timestep, domain, stage, xmom, ymom, height, bed, xvel, yvel)


    Post conditions:
      domain.explicit_update is reset to computed flux values
      returns timestep which is the largest step satisfying all volumes.


  */

  PyObject 
    *domain,
    *stage, 
    *xmom, 
    *ymom,
    *height,
    *bed,
    *xvel,
    *yvel;

    PyArrayObject 
      *num_neighbours          ,
      *num_neighbour_edges     ,
      *num_normals             ,
      *num_edgelengths         ,      
      *num_radii               ,   
      *num_areas               ,
      *num_tri_full_flag       ,
      *num_max_speed           ,
      *num_stage_edge_values   ,
      *num_xmom_edge_values    ,
      *num_ymom_edge_values    ,
      *num_height_edge_values  ,
      *num_bed_edge_values     ,
      *num_xvel_edge_values    ,
      *num_yvel_edge_values    ,
      *num_stage_boundary_values ,
      *num_xmom_boundary_values  ,
      *num_ymom_boundary_values  ,
      *num_height_boundary_values,
      *num_bed_boundary_values   ,
      *num_xvel_boundary_values  ,
      *num_yvel_boundary_values  ,
      *num_stage_explicit_update ,
      *num_xmom_explicit_update  ,
      *num_ymom_explicit_update ;

    long
      *neighbours          ,
      *neighbour_edges     ,
      *tri_full_flag       ; 

    double
      *normals             ,
      *edgelengths         ,      
      *radii               ,   
      *areas               ,
      *max_speed           ,
      *stage_edge_values   ,
      *xmom_edge_values    ,
      *ymom_edge_values    ,
      *height_edge_values  ,
      *bed_edge_values     ,
      *xvel_edge_values    ,
      *yvel_edge_values    ,
      *stage_boundary_values ,
      *xmom_boundary_values  ,
      *ymom_boundary_values  ,
      *height_boundary_values,
      *bed_boundary_values   ,
      *xvel_boundary_values  ,
      *yvel_boundary_values  ,
      *stage_explicit_update ,
      *xmom_explicit_update  ,
      *ymom_explicit_update ;

    

    double flux_timestep, g;
    double local_max_speed, length, inv_area;

    int k, i, m, n;
    int ki, nm=0, ki2; // Index shorthands
   
    double ql[7], qr[7]; // Work arrays for storing local evolving quantities
    double edgeflux[3];  // Work array for summing up fluxes

    
    //======================================================================
    // Convert Python arguments to C
    //======================================================================
    if (!PyArg_ParseTuple(args, "dOOOOOOOO", 
			  &flux_timestep, 
			  &domain, 
			  &stage, 
			  &xmom, 
			  &ymom, 
			  &height,
			  &bed,
			  &xvel, 
			  &yvel)) {
      PyErr_SetString(PyExc_RuntimeError, "Input arguments failed");
      return NULL;
    }
    
    //======================================================================
    // Extract data from python objects
    //======================================================================
    g = get_python_double(domain,"g");
    
    num_neighbours             = get_consecutive_array(domain, "neighbours");
    num_neighbour_edges        = get_consecutive_array(domain, "neighbour_edges"); 
    num_normals                = get_consecutive_array(domain, "normals");
    num_edgelengths            = get_consecutive_array(domain, "edgelengths");    
    num_radii                  = get_consecutive_array(domain, "radii");    
    num_areas                  = get_consecutive_array(domain, "areas");    
    num_tri_full_flag          = get_consecutive_array(domain, "tri_full_flag");
    num_max_speed              = get_consecutive_array(domain, "max_speed");
    
    num_stage_edge_values      = get_consecutive_array(stage, "edge_values");    
    num_xmom_edge_values       = get_consecutive_array(xmom, "edge_values");    
    num_ymom_edge_values       = get_consecutive_array(ymom, "edge_values"); 
    num_height_edge_values     = get_consecutive_array(height, "edge_values"); 
    num_bed_edge_values        = get_consecutive_array(bed, "edge_values");       
    num_xvel_edge_values       = get_consecutive_array(xvel, "edge_values");    
    num_yvel_edge_values       = get_consecutive_array(yvel, "edge_values");    
    
    num_stage_boundary_values  = get_consecutive_array(stage, "boundary_values");    
    num_xmom_boundary_values   = get_consecutive_array(xmom, "boundary_values");    
    num_ymom_boundary_values   = get_consecutive_array(ymom, "boundary_values"); 
    num_height_boundary_values = get_consecutive_array(height, "boundary_values"); 
    num_bed_boundary_values    = get_consecutive_array(bed, "boundary_values"); 
    num_xvel_boundary_values   = get_consecutive_array(xvel, "boundary_values"); 
    num_yvel_boundary_values   = get_consecutive_array(yvel, "boundary_values"); 
    
    num_stage_explicit_update  = get_consecutive_array(stage, "explicit_update");    
    num_xmom_explicit_update   = get_consecutive_array(xmom, "explicit_update");    
    num_ymom_explicit_update   = get_consecutive_array(ymom, "explicit_update"); 
    
    //---------------------------------------------------------------------------
    
    neighbours         = (long *) num_neighbours-> data;
    neighbour_edges    = (long *) num_neighbour_edges-> data;
    normals            = (double *) num_normals -> data;
    edgelengths        = (double *) num_edgelengths -> data;
    radii              = (double *) num_radii -> data;
    areas              = (double *) num_areas -> data; 
    tri_full_flag      = (long *) num_tri_full_flag -> data; 
    max_speed          = (double *) num_max_speed -> data;
    
    stage_edge_values       = (double *) num_stage_edge_values -> data;   
    xmom_edge_values        = (double *) num_xmom_edge_values -> data;
    ymom_edge_values        = (double *) num_ymom_edge_values -> data;
    height_edge_values      = (double *) num_height_edge_values -> data;
    bed_edge_values         = (double *) num_bed_edge_values -> data;
    xvel_edge_values        = (double *) num_xvel_edge_values -> data;
    yvel_edge_values        = (double *) num_yvel_edge_values -> data;
    
    stage_boundary_values      = (double *) num_stage_boundary_values -> data; 
    xmom_boundary_values       = (double *) num_xmom_boundary_values -> data;
    ymom_boundary_values       = (double *) num_ymom_boundary_values -> data;
    height_boundary_values     = (double *) num_height_boundary_values -> data;
    bed_boundary_values        = (double *) num_bed_boundary_values -> data;
    xvel_boundary_values       = (double *) num_xvel_boundary_values -> data;
    yvel_boundary_values       = (double *) num_yvel_boundary_values -> data;
    
    stage_explicit_update      = (double *) num_stage_explicit_update -> data;
    xmom_explicit_update       = (double *) num_xmom_explicit_update -> data;
    ymom_explicit_update       = (double *) num_ymom_explicit_update -> data;
    
    
    int number_of_elements = num_stage_edge_values -> dimensions[0];

    //======================================================================
    // Flux computation routine and update the explicit update arrays 
    //======================================================================

    // Set explicit_update to zero for all conserved_quantities.
    // This assumes compute_fluxes called before forcing terms

    memset((char*) stage_explicit_update, 0, number_of_elements*sizeof(double));
    memset((char*) xmom_explicit_update, 0, number_of_elements*sizeof(double));
    memset((char*) ymom_explicit_update, 0, number_of_elements*sizeof(double));    
  
    // For all triangles
    for (k = 0; k < number_of_elements; k++) 
      {  
	max_speed[k] = 0.0;

	// Loop through neighbours and compute edge flux for each  epsilon
	for (i = 0; i < 3; i++) 
	  {
	    ki = k*3 + i; // Linear index to edge i of triangle k
      
	    // Get left hand side values from triangle k, edge i 
	    ql[0] = stage_edge_values[ki];
	    ql[1] = xmom_edge_values[ki];
	    ql[2] = ymom_edge_values[ki];
	    ql[3] = height_edge_values[ki];
	    ql[4] = bed_edge_values[ki];
	    ql[5] = xvel_edge_values[ki];
	    ql[6] = yvel_edge_values[ki];

	    // Get right hand side values either from neighbouring triangle
	    // or from boundary array (Quantities at neighbour on nearest face).
	    n = neighbours[ki];
	    if (n < 0) 
	      {
		// Neighbour is a boundary condition
		m = -n - 1; // Convert negative flag to boundary index
		
		qr[0] = stage_boundary_values[m];
		qr[1] = xmom_boundary_values[m];
		qr[2] = ymom_boundary_values[m];
		qr[3] = height_boundary_values[m];
		qr[4] = bed_boundary_values[m];
		qr[5] = xvel_boundary_values[m];
		qr[6] = yvel_boundary_values[m];
		
	      } 
	    else 
	      {
		// Neighbour is a real triangle
		m = neighbour_edges[ki];
		nm = n*3 + m; // Linear index (triangle n, edge m)
		
		qr[0] = stage_edge_values[nm];
		qr[1] = xmom_edge_values[nm];
		qr[2] = ymom_edge_values[nm];
		qr[3] = height_edge_values[nm];
		qr[4] = bed_edge_values[nm];
		qr[5] = xvel_edge_values[nm];
		qr[6] = yvel_edge_values[nm];
		
	      }
	    
	    // Now we have values for this edge - both from left and right side.

	    // Outward pointing normal vector (domain.normals[k, 2*i:2*i+2])
	    ki2 = 2*ki; //k*6 + i*2

	    // Edge flux computation (triangle k, edge i)
	    _flux_function(ql,
			   qr,
			   normals[ki2], 
			   normals[ki2+1],
			   g,
			   edgeflux, 
			   &local_max_speed);
	    
	    
	    // Multiply edgeflux by edgelength
	    length = edgelengths[ki];
	    edgeflux[0] *= length;            
	    edgeflux[1] *= length;            
	    edgeflux[2] *= length;                        
      
      
	    // Update triangle k with flux from edge i
	    stage_explicit_update[k] -= edgeflux[0];
	    xmom_explicit_update[k]  -= edgeflux[1];
	    ymom_explicit_update[k]  -= edgeflux[2];
      

	    // Update flux_timestep based on edge i and possibly neighbour n
	    if (tri_full_flag[k] == 1) 
	      {
		if (local_max_speed > 1.0e-15) 
		  {
		    // Calculate safe timestep based on local_max_speed and size of triangle k
		    // A reduction of timestep based on domain.CFL is applied in update_timestep 
		    flux_timestep = min(flux_timestep, radii[k]/local_max_speed);
		    
		  }
	      }
	    
	  } // End edge i 
	
	
	// Normalise triangle k by area and store for when all conserved
	// quantities get updated
	inv_area = 1.0/areas[k];
	stage_explicit_update[k] *= inv_area;
	xmom_explicit_update[k]  *= inv_area;
	ymom_explicit_update[k]  *= inv_area;
	
	
	// Keep track of maximal speeds
	max_speed[k] = max(max_speed[k],local_max_speed);
	
      } // End triangle k
    
    
    //======================================================================
    // Cleanup
    //======================================================================

    Py_DECREF(num_neighbours            );
    Py_DECREF(num_neighbour_edges       );
    Py_DECREF(num_normals               );
    Py_DECREF(num_edgelengths           );      
    Py_DECREF(num_radii                 );   
    Py_DECREF(num_areas                 );
    Py_DECREF(num_tri_full_flag         );
    Py_DECREF(num_max_speed             );
    Py_DECREF(num_stage_edge_values     );
    Py_DECREF(num_xmom_edge_values      );
    Py_DECREF(num_ymom_edge_values      );
    Py_DECREF(num_height_edge_values    );
    Py_DECREF(num_bed_edge_values       );
    Py_DECREF(num_xvel_edge_values      );
    Py_DECREF(num_yvel_edge_values      );
    Py_DECREF(num_stage_boundary_values );
    Py_DECREF(num_xmom_boundary_values  );
    Py_DECREF(num_ymom_boundary_values  );
    Py_DECREF(num_height_boundary_values);
    Py_DECREF(num_bed_boundary_values   );
    Py_DECREF(num_xvel_boundary_values  );
    Py_DECREF(num_yvel_boundary_values  );
    Py_DECREF(num_stage_explicit_update );
    Py_DECREF(num_xmom_explicit_update  );
    Py_DECREF(num_ymom_explicit_update  );

    //======================================================================
    // Return updated flux flux_timestep
    //======================================================================
    return Py_BuildValue("d", flux_timestep);
}

//========================================================================
// Flux Function
//========================================================================

PyObject *flux_function(PyObject *self, PyObject *args) {
  //
  // Gateway to innermost flux function.
  // This is only used by the unit tests as the c implementation is
  // normally called by compute_fluxes in this module.


  PyArrayObject *normal, *ql, *qr,  *edgeflux;
  double g, local_max_speed;

  if (!PyArg_ParseTuple(args, "OOOOd",
            &normal, &ql, &qr, &edgeflux, &g)) {
    PyErr_SetString(PyExc_RuntimeError, 
		    "swb_domain_ext.c: flux_function could not parse input arguments");
    return NULL;
  }

  
  _flux_function((double*) ql -> data, 
		 (double*) qr -> data, 
		 ((double*) normal -> data)[0],
		 ((double*) normal -> data)[1],          
		 g,
		 (double*) edgeflux -> data, 
		 &local_max_speed);
  
  return Py_BuildValue("d", local_max_speed);  
}


//========================================================================
// Gravity
//========================================================================

PyObject *gravity(PyObject *self, PyObject *args) {
  //
  //  gravity(g, h, v, x, xmom, ymom)
  //
  
  
  PyArrayObject *h, *v, *x, *xmom, *ymom;
  int k, N, k3, k6;
  double g, avg_h, zx, zy;
  double x0, y0, x1, y1, x2, y2, z0, z1, z2;
  //double epsilon;
  
  if (!PyArg_ParseTuple(args, "dOOOOO",
			&g, &h, &v, &x,
			&xmom, &ymom)) {
    //&epsilon)) {
    PyErr_SetString(PyExc_RuntimeError, "shallow_water_ext.c: gravity could not parse input arguments");
    return NULL;
  }
  
  // check that numpy array objects arrays are C contiguous memory
  CHECK_C_CONTIG(h);
  CHECK_C_CONTIG(v);
  CHECK_C_CONTIG(x);
  CHECK_C_CONTIG(xmom);
  CHECK_C_CONTIG(ymom);
  
  N = h -> dimensions[0];
  for (k=0; k<N; k++) {
    k3 = 3*k;  // base index
    
    // Get bathymetry
    z0 = ((double*) v -> data)[k3 + 0];
    z1 = ((double*) v -> data)[k3 + 1];
    z2 = ((double*) v -> data)[k3 + 2];
    
    // Optimise for flat bed
    // Note (Ole): This didn't produce measurable speed up.
    // Revisit later
    //if (fabs(z0-z1)<epsilon && fabs(z1-z2)<epsilon) {
    //  continue;
    //} 
    
    // Get average depth from centroid values
    avg_h = ((double *) h -> data)[k];
    
    // Compute bed slope
    k6 = 6*k;  // base index
    
    x0 = ((double*) x -> data)[k6 + 0];
    y0 = ((double*) x -> data)[k6 + 1];
    x1 = ((double*) x -> data)[k6 + 2];
    y1 = ((double*) x -> data)[k6 + 3];
    x2 = ((double*) x -> data)[k6 + 4];
    y2 = ((double*) x -> data)[k6 + 5];
    
    
    _gradient(x0, y0, x1, y1, x2, y2, z0, z1, z2, &zx, &zy);
    
    // Update momentum
    ((double*) xmom -> data)[k] += -g*zx*avg_h;
    ((double*) ymom -> data)[k] += -g*zy*avg_h;
  }
  
  return Py_BuildValue("");
}



//========================================================================
// Method table for python module
//========================================================================

static struct PyMethodDef MethodTable[] = {
  /* The cast of the function is necessary since PyCFunction values
   * only take two PyObject* parameters, and rotate() takes
   * three.
   */
  {"compute_fluxes_c", compute_fluxes,     METH_VARARGS, "Print out"},
  {"gravity_c",        gravity,            METH_VARARGS, "Print out"},
  {"flux_function_c",  flux_function,      METH_VARARGS, "Print out"},
  {NULL, NULL, 0, NULL}
};

// Module initialisation
void initswb_domain_ext(void){
  Py_InitModule("swb_domain_ext", MethodTable);

  import_array(); // Necessary for handling of NumPY structures
}

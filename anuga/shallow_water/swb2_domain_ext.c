// Python - C extension module for shallow_water.py
//
// To compile (Python2.6):
//  gcc -c swb2_domain_ext.c -I/usr/include/python2.6 -o domain_ext.o -Wall -O
//  gcc -shared swb2_domain_ext.o  -o swb2_domain_ext.so
//
// or use python compile.py
//
// See the module swb_domain.py for more documentation on 
// how to use this module
//
//
// Stephen Roberts, ANU 2009
// Ole Nielsen, GA 2004
// Gareth Davies, GA 2011

#include "math.h"
#include <stdio.h>


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

// Function to obtain speed from momentum and depth.
// This is used by flux functions
// Input parameters uh and h may be modified by this function.
// Tried to inline, but no speedup was achieved 27th May 2009 (Ole)
//static inline double _compute_speed(double *uh, 
double _compute_speed(double *uh, 
		      double *h, 
		      double epsilon, 
		      double h0,
		      double limiting_threshold) {
  
  double u;

  if (*h < limiting_threshold) {   
    // Apply limiting of speeds according to the ANUGA manual
    if (*h < epsilon) {
      //*h = fmax(0.0,*h);  // Could have been negative
      *h = 0.0;  // Could have been negative
      u = 0.0;
    } else {
      u = *uh/(*h + h0/ *h);    
    }
  

    // Adjust momentum to be consistent with speed
    *uh = u * *h;
  } else {
    // We are in deep water - no need for limiting
    u = *uh/ *h;
  }
  
  return u;
}

// Innermost flux function (using stage w=z+h)
int _flux_function_central(double *q_left, double *q_right,
                           double z_left, double z_right,
                           double n1, double n2,
                           double epsilon, 
                           double h0,
                           double limiting_threshold, 
                           double g,
                           double *edgeflux, double *max_speed, 
                           double *pressure_flux) 
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

  double w_left, h_left, uh_left, vh_left, u_left;
  double w_right, h_right, uh_right, vh_right, u_right;
  double s_min, s_max, soundspeed_left, soundspeed_right;
  double denom, inverse_denominator, z;
  double t3;
  // Workspace (allocate once, use many)
  static double q_left_rotated[3], q_right_rotated[3], flux_right[3], flux_left[3];

  // Copy conserved quantities to protect from modification
  q_left_rotated[0] = q_left[0];
  q_right_rotated[0] = q_right[0];
  q_left_rotated[1] = q_left[1];
  q_right_rotated[1] = q_right[1];
  q_left_rotated[2] = q_left[2];
  q_right_rotated[2] = q_right[2];

  // Align x- and y-momentum with x-axis
  _rotate(q_left_rotated, n1, n2);
  _rotate(q_right_rotated, n1, n2);

  z = 0.5*(z_left + z_right); // Average elevation values. 
                              // Even though this will nominally allow 
			      // for discontinuities in the elevation data, 
			      // there is currently no numerical support for
			      // this so results may be strange near 
			      // jumps in the bed.

  // Compute speeds in x-direction
  w_left = q_left_rotated[0];          
  h_left = w_left - z;
  uh_left = q_left_rotated[1];
  u_left = _compute_speed(&uh_left, &h_left, 
			  epsilon, h0, limiting_threshold);

  w_right = q_right_rotated[0];
  h_right = w_right - z;
  uh_right = q_right_rotated[1];
  u_right = _compute_speed(&uh_right, &h_right, 
			   epsilon, h0, limiting_threshold);

  // Momentum in y-direction
  vh_left  = q_left_rotated[2];
  vh_right = q_right_rotated[2];

  // Limit y-momentum if necessary 
  // Leaving this out, improves speed significantly (Ole 27/5/2009)
  // All validation tests pass, so do we really need it anymore? 
  _compute_speed(&vh_left, &h_left, 
		 epsilon, h0, limiting_threshold);
  _compute_speed(&vh_right, &h_right, 
		 epsilon, h0, limiting_threshold);

  // Maximal and minimal wave speeds
  soundspeed_left  = sqrt(g*h_left);
  soundspeed_right = sqrt(g*h_right);  
  
  // Code to use fast square root optimisation if desired.
  // Timings on AMD 64 for the Okushiri profile gave the following timings
  //
  // SQRT           Total    Flux
  //=============================
  //
  // Ref            405s     152s
  // Fast (dbl)     453s     173s
  // Fast (sng)     437s     171s
  //
  // Consequently, there is currently (14/5/2009) no reason to use this 
  // approximation.
  
  //soundspeed_left  = fast_squareroot_approximation(g*h_left);
  //soundspeed_right = fast_squareroot_approximation(g*h_right);

  s_max = fmax(u_left + soundspeed_left, u_right + soundspeed_right);
  if (s_max < 0.0) 
  {
    s_max = 0.0;
  }

  s_min = fmin(u_left - soundspeed_left, u_right - soundspeed_right);
  if (s_min > 0.0)
  {
    s_min = 0.0;
  }
  
  // Flux formulas
  flux_left[0] = u_left*h_left;
  flux_left[1] = u_left*uh_left + 0.5*g*h_left*h_left;
  flux_left[2] = u_left*vh_left;

  flux_right[0] = u_right*h_right;
  flux_right[1] = u_right*uh_right + 0.5*g*h_right*h_right;
  flux_right[2] = u_right*vh_right;

  // Flux computation
  denom = s_max - s_min;
  if (denom < epsilon) 
  { // FIXME (Ole): Try using h0 here
    memset(edgeflux, 0, 3*sizeof(double));
    //for(i=0;i<3;i++){
    //    edgeflux[i] = _minmod(flux_left[i], flux_right[i]);
    //}
    *max_speed = 0.0;
  } 
  else 
  {
    inverse_denominator = 1.0/denom;
    for (i = 0; i < 3; i++) 
    {
      // Adjustment to the scheme by Kurganov and Lin 2007 Communications in Computational
      // Physics 2:141-163
      //uint = (s_max*q_right_rotated[i] - s_min*q_left_rotated[i] - (flux_right[i] - flux_left[i]))*inverse_denominator;
      //t1 = (q_right_rotated[i] - uint);
      //t2 = (-q_left_rotated[i] + uint);
      //t3 = _minmod(t1,t2);
      t3 = 0.0;
      edgeflux[i] = s_max*flux_left[i] - s_min*flux_right[i];
      edgeflux[i] += s_max*s_min*(q_right_rotated[i] - q_left_rotated[i] - t3);
      edgeflux[i] *= inverse_denominator;
    }

    *pressure_flux = 0.;//0.5*g*( s_max*h_left*h_left -s_min*h_right*h_right)*inverse_denominator;

    // Maximal wavespeed
    *max_speed = fmax(fabs(s_max), fabs(s_min));

    // Rotate back
    _rotate(edgeflux, n1, -n2);
  }
  
  return 0;
}


// Computational function for flux computation
double _compute_fluxes_central(int number_of_elements,
        double timestep,
        double epsilon,
        double H0,
        double g,
        long* neighbours,
        long* neighbour_edges,
        double* normals,
        double* edgelengths,
        double* radii,
        double* areas,
        long* tri_full_flag,
        double* stage_edge_values,
        double* xmom_edge_values,
        double* ymom_edge_values,
        double* bed_edge_values,
        double* stage_boundary_values,
        double* xmom_boundary_values,
        double* ymom_boundary_values,
        long*   boundary_flux_type,
        double* stage_explicit_update,
        double* xmom_explicit_update,
        double* ymom_explicit_update,
        long* already_computed_flux,
        double* max_speed_array,
        int optimise_dry_cells, 
        double* stage_centroid_values,
        double* bed_centroid_values,
        double* bed_vertex_values) {
    // Local variables
    double max_speed, length, inv_area, zl, zr;
    double h0 = H0*H0; // This ensures a good balance when h approaches H0.

    double limiting_threshold = 10 * H0; // Avoid applying limiter below this
    // threshold for performance reasons.
    // See ANUGA manual under flux limiting
    int k, i, m, n;
    int ki, nm = 0, ki2; // Index shorthands

    // Workspace (making them static actually made function slightly slower (Ole))
    double ql[3], qr[3], edgeflux[3]; // Work array for summing up fluxes
    double bedslope_work;
    double pressure_flux, hc, hc_n;
    static long call = 1; // Static local variable flagging already computed flux


    //max_bed_edgevalue = malloc(number_of_elements*sizeof(double));
    //min_bed_edgevalue = malloc(number_of_elements*sizeof(double));
    // Start computation
    call++; // Flag 'id' of flux calculation for this timestep

    // Set explicit_update to zero for all conserved_quantities.
    // This assumes compute_fluxes called before forcing terms
    memset((char*) stage_explicit_update, 0, number_of_elements * sizeof (double));
    memset((char*) xmom_explicit_update, 0, number_of_elements * sizeof (double));
    memset((char*) ymom_explicit_update, 0, number_of_elements * sizeof (double));


    // Compute minimum bed edge value on each triangle
    //for (k = 0; k < number_of_elements; k++){
    //    max_bed_edgevalue[k] = fmax(bed_edge_values[3*k], 
    //                               fmax(bed_edge_values[3*k+1], bed_edge_values[3*k+2]));
    //    min_bed_edgevalue[k] = fmin(bed_edge_values[3*k], 
    //                               fmin(bed_edge_values[3*k+1], bed_edge_values[3*k+2]));
    //
    //}


    // For all triangles
    for (k = 0; k < number_of_elements; k++) {
        // Loop through neighbours and compute edge flux for each
        for (i = 0; i < 3; i++) {
            ki = k * 3 + i; // Linear index to edge i of triangle k

            if (already_computed_flux[ki] == call) {
                // We've already computed the flux across this edge
                continue;
            }

            // Get left hand side values from triangle k, edge i
            ql[0] = stage_edge_values[ki];
            ql[1] = xmom_edge_values[ki];
            ql[2] = ymom_edge_values[ki];
            zl = bed_edge_values[ki];
            hc = fmax(stage_centroid_values[k] - bed_centroid_values[k],0.0);

            // Get right hand side values either from neighbouring triangle
            // or from boundary array (Quantities at neighbour on nearest face).
            n = neighbours[ki];
            hc_n = hc;
            if (n < 0) {
                // Neighbour is a boundary condition
                m = -n - 1; // Convert negative flag to boundary index

                qr[0] = stage_boundary_values[m];
                qr[1] = xmom_boundary_values[m];
                qr[2] = ymom_boundary_values[m];
                zr = zl; // Extend bed elevation to boundary
            }
            else {
                // Neighbour is a real triangle
                hc_n = fmax(stage_centroid_values[n] - bed_centroid_values[n],0.0);
                m = neighbour_edges[ki];
                nm = n * 3 + m; // Linear index (triangle n, edge m)

                qr[0] = stage_edge_values[nm];
                qr[1] = xmom_edge_values[nm];
                qr[2] = ymom_edge_values[nm];
                zr = bed_edge_values[nm];
            }
            

            if (fabs(zl-zr)>1.0e-10) {
                //report_python_error(AT,"Discontinuous Elevation");
                return 0.0;
            }
            
            // If both centroids are dry, then the cells should not exchange flux 
            // NOTE: This also means no bedslope term -- which is arguably
            // appropriate, since the depth is zero at the cell centroid
            if(( hc == 0.0)&(hc_n == 0.0)&(n>=0)){
                already_computed_flux[ki] = call; // #k Done
                already_computed_flux[nm] = call; // #n Done
                max_speed = 0.0;
                continue;
            }
           
            //If one centroid is dry, then extrapolate its edge values from the neighbouring centroid,
            // unless the local centroid value is smaller 
            if(n>=0){
                if(hc==0.0){
                    //ql[0]=fmax(fmin(qr[0],stage_centroid_values[k]),zl);
                    //ql[0]=fmax(fmin(qr[0],0.5*(stage_centroid_values[k]+stage_centroid_values[n])),zl);
                    ql[0]=fmax(fmin(qr[0],stage_centroid_values[k]),zl);
                }
                if(hc_n==0.0){
                    qr[0]=fmax(fmin(ql[0],stage_centroid_values[n]),zr);
                    //qr[0]=fmax(fmin(ql[0],0.5*(stage_centroid_values[n]+stage_centroid_values[k])),zr);
                    //qr[0]=ql[0]; 
                }
            }else{
                // Treat the boundary case
                //if((hc==0.0)){
                //    ql[0]=fmax(fmin(qr[0],stage_centroid_values[k]),zl);
                    //ql[0]=fmax(fmin(qr[0],ql[0]),zl);
                //}
            }
            
            // Outward pointing normal vector (domain.normals[k, 2*i:2*i+2])
            ki2 = 2 * ki; //k*6 + i*2

            // Edge flux computation (triangle k, edge i)
            _flux_function_central(ql, qr, zl, zr,
                    normals[ki2], normals[ki2 + 1],
                    epsilon, h0, limiting_threshold, g,
                    edgeflux, &max_speed, &pressure_flux);

            // Prevent outflow from 'seriously' dry cells
            // Idea: The cell will not go dry if:
            // mass_flux <= Area_triangle*hc/(dt*edgelength)
            //if((stage_centroid_values[k]<=max_bed_edgevalue[k])|
            //   (ql[0]<=zl)){
            //    if(edgeflux[0]>0.0){
            //        tmp=fmin(0.5*areas[k]*(hc+bed_centroid_values[k] - min_bed_edgevalue[k])/(edgelengths[ki]*fmax(timestep,epsilon)), 1.0); // 28/7 -- Right answer for channel flow problem.
            //        tmp = fmin(fmin(edgeflux[0], tmp)/edgeflux[0], 1.0);
            //        edgeflux[0]*=tmp;
            //    }
            //}
            //if(n>=0){
            //    if((stage_centroid_values[n]<=max_bed_edgevalue[n])|
            //       (qr[0]<=zr)){
            //        if(edgeflux[0]<0.0){
            //            tmp=fmin(0.5*areas[n]*(hc_n+bed_centroid_values[n] - min_bed_edgevalue[n])/(edgelengths[ki]*fmax(timestep,epsilon)), 1.0); // 28/7 -- Right answer for channel flow problem.
            //            tmp = fmin( fmax(edgeflux[0], -tmp)/edgeflux[0], 1.0);
            //            edgeflux[0]*=tmp;
            //        }
            //    }
            //} 

            // Multiply edgeflux by edgelength
            length = edgelengths[ki];
            edgeflux[0] *= length;
            edgeflux[1] *= length;
            edgeflux[2] *= length;

            // Update triangle k with flux from edge i
            stage_explicit_update[k] -= edgeflux[0];
            xmom_explicit_update[k] -= edgeflux[1];
            ymom_explicit_update[k] -= edgeflux[2];

            // Compute bed slope term
            //if(hc>-9999.0){
                //Bedslope approx 1:
            bedslope_work = g*length*( hc*(ql[0])-0.5*fmax(ql[0]-zl,0.)*(ql[0]-zl) );
                //
                // Bedslope approx 2
                //stage_edge_lim = ql[0]; // Limit this to be between a constant stage and constant depth extrapolation
                //if(stage_edge_lim > fmax(stage_centroid_values[k], zl +hc)){
                //    stage_edge_lim = fmax(stage_centroid_values[k], zl+hc);
                //}
                //if(stage_edge_lim < fmin(stage_centroid_values[k], zl +hc)){
                //    stage_edge_lim = fmin(stage_centroid_values[k], zl+hc);
                //}
                //bedslope_work = g*hc*(stage_edge_lim)*length-0.5*g*fmax(stage_edge_lim-zl,0.)*(stage_edge_lim-zl)*length;

                // Bedslope approx 3
                //bedslope_work = -0.5*g*fmax(stage_centroid_values[k]-zl,0.)*(stage_centroid_values[k]-zl)*length;
                //
            xmom_explicit_update[k] -= normals[ki2]*bedslope_work;
            ymom_explicit_update[k] -= normals[ki2+1]*bedslope_work;
            //}else{
            //   // Treat nearly dry cells 
            //   bedslope_work =-0.5*g*fmax(ql[0]-zl,0.)*(ql[0]-zl)*length; //
            //   //
            //   //bedslope_work = -pressure_flux*length;
            //   xmom_explicit_update[k] -= normals[ki2]*bedslope_work;
            //   ymom_explicit_update[k] -= normals[ki2+1]*bedslope_work;
            //   //xmom_explicit_update[k] = 0.0;
            //   //ymom_explicit_update[k] = 0.0;

            ////}

            already_computed_flux[ki] = call; // #k Done


            // Update neighbour n with same flux but reversed sign
            if (n >= 0) {
                stage_explicit_update[n] += edgeflux[0];
                xmom_explicit_update[n] += edgeflux[1];
                ymom_explicit_update[n] += edgeflux[2];
                //Add bed slope term here
                //if(hc_n>-9999.0){
                //if(stage_centroid_values[n] > max_bed_edgevalue[n]){
                    // Bedslope approx 1:
                bedslope_work = g*length*(hc_n*(qr[0])-0.5*fmax(qr[0]-zr,0.)*(qr[0]-zr));
                    //
                    // Bedslope approx 2:
                    //stage_edge_lim = qr[0];
                    //if(stage_edge_lim > fmax(stage_centroid_values[n], zr +hc_n)){
                    //    stage_edge_lim = fmax(stage_centroid_values[n], zr+hc_n);
                    //}
                    //if(stage_edge_lim < fmin(stage_centroid_values[n], zr +hc_n)){
                    //    stage_edge_lim = fmin(stage_centroid_values[n], zr+hc_n);
                    //}
                    //bedslope_work = g*hc_n*(stage_edge_lim)*length-0.5*g*fmax(stage_edge_lim-zr,0.)*(stage_edge_lim-zr)*length;
                    //
                    // Bedslope approx 3
                    //bedslope_work = -0.5*g*fmax(stage_centroid_values[n]-zr,0.)*(stage_centroid_values[n]-zr)*length;
                    //
                xmom_explicit_update[n] += normals[ki2]*bedslope_work;
                ymom_explicit_update[n] += normals[ki2+1]*bedslope_work;
                //}else{
                //    // Treat nearly dry cells
                //    bedslope_work = -0.5*g*fmax(qr[0]-zr,0.)*(qr[0]-zr)*length; //-pressure_flux*length; //-0.5*g*fmax(qr[0]-zr,0.)*(qr[0]-zr)*length;
                //    //bedslope_work = -pressure_flux*length;
                //    xmom_explicit_update[n] += normals[ki2]*bedslope_work;
                //    ymom_explicit_update[n] += normals[ki2+1]*bedslope_work;

                //    //xmom_explicit_update[n] = 0.0;
                //    //ymom_explicit_update[n] = 0.0;
                //}

                already_computed_flux[nm] = call; // #n Done
            }

            // Update timestep based on edge i and possibly neighbour n
            if (tri_full_flag[k] == 1) {
                if (max_speed > epsilon) {
                    // Apply CFL condition for triangles joining this edge (triangle k and triangle n)

                    // CFL for triangle k
                    timestep = fmin(timestep, radii[k] / max_speed);

                    if (n >= 0) {
                        // Apply CFL condition for neigbour n (which is on the ith edge of triangle k)
                        timestep = fmin(timestep, radii[n] / max_speed);
                    }

                    // Ted Rigby's suggested less conservative version
                    //if(n>=0){
                    //  timestep = fmin(timestep, (radii[k]+radii[n])/max_speed);
                    //}else{
                    //  timestep = fmin(timestep, radii[k]/max_speed);
                    //}
                }
            }

        } // End edge i (and neighbour n)


        // Normalise triangle k by area and store for when all conserved
        // quantities get updated
        inv_area = 1.0 / areas[k];
        stage_explicit_update[k] *= inv_area;
        xmom_explicit_update[k] *= inv_area;
        ymom_explicit_update[k] *= inv_area;


        // Keep track of maximal speeds
        max_speed_array[k] = max_speed;

    } // End triangle k

    //free(max_bed_edgevalue);
    //free(min_bed_edgevalue);

    return timestep;
}

// Protect against the water elevation falling below the triangle bed
double _protect(int N,
         double minimum_allowed_height,
         double maximum_allowed_speed,
         double epsilon,
         double* wc,
         double* wv,
         double* zc,
         double* zv,
         double* xmomc,
         double* ymomc,
         double* areas) {

  int k;
  double hc, bmin, bmax;
  double mass_error=0.; 
  // This acts like minimum_allowed height, but scales with the vertical
  // distance between the bed_centroid_value and the max bed_edge_value of
  // every triangle.
  double minimum_relative_height=0.1; 
  //int mass_added=0;

  // Protect against inifintesimal and negative heights  
  //if (maximum_allowed_speed < epsilon) {
    for (k=0; k<N; k++) {
      hc = wc[k] - zc[k];
      // Definine the maximum bed edge value on triangle k.
      bmax = 0.5*fmax((zv[3*k]+zv[3*k+1]),fmax((zv[3*k+1]+zv[3*k+2]),(zv[3*k+2]+zv[3*k])));

      if (hc < fmax(minimum_relative_height*(bmax-zc[k]), minimum_allowed_height)) {
        
        // Set momentum to zero and ensure h is non negative
        // NOTE: THIS IS IMPORTANT -- WE ARE SETTING MOMENTUM TO ZERO
        //if(hc<=epsilon){
            xmomc[k] = 0.0;
            ymomc[k] = 0.0;
        //}

        if (hc <= 0.0){
             // Definine the minimum bed edge value on triangle k.
             // Setting = minimum edge value can lead to mass conservation problems
             //bmin = 0.5*fmin((zv[3*k]+zv[3*k+1]),fmin((zv[3*k+1]+zv[3*k+2]),(zv[3*k+2]+zv[3*k])));
             //bmin =0.5*bmin + 0.5*fmin(zv[3*k],fmin(zv[3*k+1],zv[3*k+2]));
             // Setting = minimum vertex value seems better, but might tend to be less smooth 
             bmin =fmin(zv[3*k],fmin(zv[3*k+1],zv[3*k+2])) -minimum_allowed_height;
             //bmin=zc[k]-minimum_allowed_height;
             // Minimum allowed stage = bmin
             // WARNING: ADDING MASS if wc[k]<bmin
             if(wc[k]<bmin){
                 mass_error+=(bmin-wc[k])*areas[k];
                 //mass_added=1; //Flag to warn of added mass                
                 //printf("Adding mass to dry cell %d %f %f %f %f %f \n", k, zv[3*k], zv[3*k+1], zv[3*k+2], wc[k]- bmin, mass_error);
             
                 wc[k] = fmax(wc[k], bmin); 
             

                 // Set vertex values as well. Seems that this shouldn't be
                 // needed. However from memory this is important at the first
                 // time step, for 'dry' areas where the designated stage is
                 // less than the bed centroid value
                 wv[3*k] = fmax(wv[3*k], bmin);
                 wv[3*k+1] = fmax(wv[3*k+1], bmin);
                 wv[3*k+2] = fmax(wv[3*k+2], bmin);
            }
        }
      }
    }

/*
  if(mass_added==1){
     printf("Cumulative mass protection: %f m^3 \n", mass_error);
  }
*/
  return mass_error;
}

int find_qmin_and_qfmax(double dq0, double dq1, double dq2, 
               double *qmin, double *qmax){
  // Considering the centroid of an FV triangle and the vertices of its 
  // auxiliary triangle, find 
  // qmin=fmin(q)-qc and qmax=fmax(q)-qc, 
  // where fmin(q) and fmax(q) are respectively min and max over the
  // four values (at the centroid of the FV triangle and the auxiliary 
  // triangle vertices),
  // and qc is the centroid
  // dq0=q(vertex0)-q(centroid of FV triangle)
  // dq1=q(vertex1)-q(vertex0)
  // dq2=q(vertex2)-q(vertex0)

  // This is a simple implementation 
  *qmax = fmax(fmax(dq0, fmax(dq0+dq1, dq0+dq2)), 0.0) ;
  *qmin = fmin(fmin(dq0, fmin(dq0+dq1, dq0+dq2)), 0.0) ;
 
  return 0;
}

int limit_gradient(double *dqv, double qmin, double qmax, double beta_w){
  // Given provisional jumps dqv from the FV triangle centroid to its 
  // vertices and jumps qmin (qmax) between the centroid of the FV 
  // triangle and the minimum (maximum) of the values at the centroid of 
  // the FV triangle and the auxiliary triangle vertices,
  // calculate a multiplicative factor phi by which the provisional 
  // vertex jumps are to be limited
  
  int i;
  double r=1000.0, r0=1.0, phi=1.0;
  static double TINY = 1.0e-100; // to avoid machine accuracy problems.
  // FIXME: Perhaps use the epsilon used elsewhere.
  
  // Any provisional jump with magnitude < TINY does not contribute to 
  // the limiting process.
  
  for (i=0;i<3;i++){
    if (dqv[i]<-TINY)
      r0=qmin/dqv[i];
      
    if (dqv[i]>TINY)
      r0=qmax/dqv[i];
      
    r=fmin(r0,r);
  }
  
  phi=fmin(r*beta_w,1.0);
  //phi=1.;
  dqv[0]=dqv[0]*phi;
  dqv[1]=dqv[1]*phi;
  dqv[2]=dqv[2]*phi;
    
  return 0;
}

// Computational routine
int _extrapolate_second_order_edge_sw(int number_of_elements,
                                 double epsilon,
                                 double minimum_allowed_height,
                                 double beta_w,
                                 double beta_w_dry,
                                 double beta_uh,
                                 double beta_uh_dry,
                                 double beta_vh,
                                 double beta_vh_dry,
                                 long* surrogate_neighbours,
                                 long* number_of_boundaries,
                                 double* centroid_coordinates,
                                 double* stage_centroid_values,
                                 double* xmom_centroid_values,
                                 double* ymom_centroid_values,
                                 double* elevation_centroid_values,
                                 double* edge_coordinates,
                                 double* stage_edge_values,
                                 double* xmom_edge_values,
                                 double* ymom_edge_values,
                                 double* elevation_edge_values,
                                 double* stage_vertex_values,
                                 double* xmom_vertex_values,
                                 double* ymom_vertex_values,
                                 double* elevation_vertex_values,
                                 int optimise_dry_cells, 
                                 int extrapolate_velocity_second_order) {
                  
  // Local variables
  double a, b; // Gradient vector used to calculate edge values from centroids
  int k, k0, k1, k2, k3, k6, coord_index, i;
  double x, y, x0, y0, x1, y1, x2, y2, xv0, yv0, xv1, yv1, xv2, yv2; // Vertices of the auxiliary triangle
  double dx1, dx2, dy1, dy2, dxv0, dxv1, dxv2, dyv0, dyv1, dyv2, dq0, dq1, dq2, area2, inv_area2;
  double dqv[3], qmin, qmax, hmin, bedmax, stagemin;
  double hc, h0, h1, h2, beta_tmp, hfactor;
  double dk, de[3];
  
  double *xmom_centroid_store, *ymom_centroid_store, *stage_centroid_store, *max_elevation_edgevalue;
  //int *count_wet_neighbours;

  // Use malloc to avoid putting these variables on the stack, which can cause
  // segfaults in large model runs
  xmom_centroid_store = malloc(number_of_elements*sizeof(double));
  ymom_centroid_store = malloc(number_of_elements*sizeof(double));
  stage_centroid_store = malloc(number_of_elements*sizeof(double));
  //min_elevation_edgevalue = malloc(number_of_elements*sizeof(double));
  max_elevation_edgevalue = malloc(number_of_elements*sizeof(double));
  //count_wet_neighbours = malloc(number_of_elements*sizeof(int));
 
  if(extrapolate_velocity_second_order==1){
      // Replace momentum centroid with velocity centroid to allow velocity
      // extrapolation This will be changed back at the end of the routine
      for (k=0; k<number_of_elements; k++){

          dk = fmax(stage_centroid_values[k]-elevation_centroid_values[k],minimum_allowed_height);
          xmom_centroid_store[k] = xmom_centroid_values[k];
          xmom_centroid_values[k] = xmom_centroid_values[k]/dk;

          ymom_centroid_store[k] = ymom_centroid_values[k];
          ymom_centroid_values[k] = ymom_centroid_values[k]/dk;

          //min_elevation_edgevalue[k] = fmin(elevation_edge_values[3*k], 
          //                                 fmin(elevation_edge_values[3*k+1],
          //                                     elevation_edge_values[3*k+2]));
          max_elevation_edgevalue[k] = fmax(elevation_edge_values[3*k], 
                                           fmax(elevation_edge_values[3*k+1],
                                               elevation_edge_values[3*k+2]));
          }

      }

  // Count how many 'fully submerged' neighbours the cell has
  //for(k=0; k<number_of_elements;k++){ 
  //    count_wet_neighbours[k]=0;
  //    for (i=0; i<3; i++){
  //      ktmp = surrogate_neighbours[3*k+i];              
  //      if(stage_centroid_values[ktmp] > max_elevation_edgevalue[ktmp]){
  //          count_wet_neighbours[k]+=1;
  //      }
  //    }
  //}

  // Begin extrapolation routine
  for (k = 0; k < number_of_elements; k++) 
  {
    k3=k*3;
    k6=k*6;

    if (number_of_boundaries[k]==3)
    //if (0==0)
    {
      // No neighbours, set gradient on the triangle to zero
      
      stage_edge_values[k3]   = stage_centroid_values[k];
      stage_edge_values[k3+1] = stage_centroid_values[k];
      stage_edge_values[k3+2] = stage_centroid_values[k];
      xmom_edge_values[k3]    = xmom_centroid_values[k];
      xmom_edge_values[k3+1]  = xmom_centroid_values[k];
      xmom_edge_values[k3+2]  = xmom_centroid_values[k];
      ymom_edge_values[k3]    = ymom_centroid_values[k];
      ymom_edge_values[k3+1]  = ymom_centroid_values[k];
      ymom_edge_values[k3+2]  = ymom_centroid_values[k];
      
      continue;
    }
    else 
    {
      // Triangle k has one or more neighbours. 
      // Get centroid and edge coordinates of the triangle
      
      // Get the edge coordinates
      xv0 = edge_coordinates[k6];   
      yv0 = edge_coordinates[k6+1];
      xv1 = edge_coordinates[k6+2]; 
      yv1 = edge_coordinates[k6+3];
      xv2 = edge_coordinates[k6+4]; 
      yv2 = edge_coordinates[k6+5];
      
      // Get the centroid coordinates
      coord_index = 2*k;
      x = centroid_coordinates[coord_index];
      y = centroid_coordinates[coord_index+1];
      
      // Store x- and y- differentials for the edges of 
      // triangle k relative to the centroid
      dxv0 = xv0 - x; 
      dxv1 = xv1 - x; 
      dxv2 = xv2 - x;
      dyv0 = yv0 - y; 
      dyv1 = yv1 - y; 
      dyv2 = yv2 - y;
      // Compute the minimum distance from the centroid to an edge
      //demin=fmin(dxv0*dxv0 +dyv0*dyv0, fmin(dxv1*dxv1+dyv1*dyv1, dxv2*dxv2+dyv2*dyv2));
      //demin=sqrt(demin);
    }


            
    if (number_of_boundaries[k]<=1)
    {
      //==============================================
      // Number of boundaries <= 1
      //==============================================    
    
    
      // If no boundaries, auxiliary triangle is formed 
      // from the centroids of the three neighbours
      // If one boundary, auxiliary triangle is formed 
      // from this centroid and its two neighbours
      
      k0 = surrogate_neighbours[k3];
      k1 = surrogate_neighbours[k3 + 1];
      k2 = surrogate_neighbours[k3 + 2];
     
      // Take note if the max neighbour bed elevation is greater than the min
      // neighbour stage -- suggests a 'steep' bed relative to the flow
      bedmax = fmax(elevation_centroid_values[k], 
                   fmax(elevation_centroid_values[k0],
                       fmax(elevation_centroid_values[k1], elevation_centroid_values[k2])));
      //bedmax = elevation_centroid_values[k];
      stagemin = fmin(fmax(stage_centroid_values[k], elevation_centroid_values[k]), 
                     fmin(fmax(stage_centroid_values[k0], elevation_centroid_values[k0]),
                         fmin(fmax(stage_centroid_values[k1], elevation_centroid_values[k1]),
                             fmax(stage_centroid_values[k2], elevation_centroid_values[k2]))));
      if(stagemin < bedmax){
         // This will cause first order extrapolation
         k2 = k;
         k0 = k;
         k1 = k;
      } 
     
      // Get the auxiliary triangle's vertex coordinates 
      // (really the centroids of neighbouring triangles)
      coord_index = 2*k0;
      x0 = centroid_coordinates[coord_index];
      y0 = centroid_coordinates[coord_index+1];
      
      coord_index = 2*k1;
      x1 = centroid_coordinates[coord_index];
      y1 = centroid_coordinates[coord_index+1];
      
      coord_index = 2*k2;
      x2 = centroid_coordinates[coord_index];
      y2 = centroid_coordinates[coord_index+1];

      // Store x- and y- differentials for the vertices 
      // of the auxiliary triangle
      dx1 = x1 - x0; 
      dx2 = x2 - x0;
      dy1 = y1 - y0; 
      dy2 = y2 - y0;
      
      // Calculate 2*area of the auxiliary triangle
      // The triangle is guaranteed to be counter-clockwise      
      area2 = dy2*dx1 - dy1*dx2; 
     
      // Treat triangles with zero or 1 wet neighbours. 
      if ((area2 <= 0)) //|(count_wet_neighbours[k]==0))
      {
          //printf("Error negative triangle area \n");
          //report_python_error(AT, "Negative triangle area");
          //return -1;
          stage_edge_values[k3]   = stage_centroid_values[k];
          stage_edge_values[k3+1] = stage_centroid_values[k];
          stage_edge_values[k3+2] = stage_centroid_values[k];
          // First order momentum / velocity extrapolation
          xmom_edge_values[k3]    = xmom_centroid_values[k];
          xmom_edge_values[k3+1]  = xmom_centroid_values[k];
          xmom_edge_values[k3+2]  = xmom_centroid_values[k];
          ymom_edge_values[k3]    = ymom_centroid_values[k];
          ymom_edge_values[k3+1]  = ymom_centroid_values[k];
          ymom_edge_values[k3+2]  = ymom_centroid_values[k];

          continue;
      }  
      
      // Calculate heights of neighbouring cells
      hc = stage_centroid_values[k]  - elevation_centroid_values[k];
      h0 = stage_centroid_values[k0] - elevation_centroid_values[k0];
      h1 = stage_centroid_values[k1] - elevation_centroid_values[k1];
      h2 = stage_centroid_values[k2] - elevation_centroid_values[k2];
      hmin = fmin(fmin(h0, fmin(h1, h2)), hc);

      hfactor = 0.0;
      //if (hmin > 0.001) 
      if (hmin > 0.) 
      //if (hc>0.0)
      {
        hfactor = 1.0 ;//hmin/(hmin + 0.004);
        //hfactor=hmin/(hmin + 0.004);
      }
      
    
      //-----------------------------------
      // stage
      //-----------------------------------      
      
      // Calculate the difference between vertex 0 of the auxiliary 
      // triangle and the centroid of triangle k
      dq0 = stage_centroid_values[k0] - stage_centroid_values[k];
      
      // Calculate differentials between the vertices 
      // of the auxiliary triangle (centroids of neighbouring triangles)
      dq1 = stage_centroid_values[k1] - stage_centroid_values[k0];
      dq2 = stage_centroid_values[k2] - stage_centroid_values[k0];
      
      inv_area2 = 1.0/area2;
      // Calculate the gradient of stage on the auxiliary triangle
      a = dy2*dq1 - dy1*dq2;
      a *= inv_area2;
      b = dx1*dq2 - dx2*dq1;
      b *= inv_area2;
      
      // Calculate provisional jumps in stage from the centroid 
      // of triangle k to its vertices, to be limited
      dqv[0] = a*dxv0 + b*dyv0;
      dqv[1] = a*dxv1 + b*dyv1;
      dqv[2] = a*dxv2 + b*dyv2;
    
      // Now we want to find min and max of the centroid and the 
      // vertices of the auxiliary triangle and compute jumps 
      // from the centroid to the min and max
      find_qmin_and_qfmax(dq0, dq1, dq2, &qmin, &qmax);
      
      beta_tmp = beta_w_dry + (beta_w - beta_w_dry) * hfactor;
    
      
      // Limit the gradient
      limit_gradient(dqv, qmin, qmax, beta_tmp);
      stage_edge_values[k3+0] = stage_centroid_values[k] + dqv[0];
      stage_edge_values[k3+1] = stage_centroid_values[k] + dqv[1];
      stage_edge_values[k3+2] = stage_centroid_values[k] + dqv[2];

       
      //-----------------------------------
      // xmomentum
      //-----------------------------------            

      // Calculate the difference between vertex 0 of the auxiliary 
      // triangle and the centroid of triangle k      
      dq0 = xmom_centroid_values[k0] - xmom_centroid_values[k];
      
      // Calculate differentials between the vertices 
      // of the auxiliary triangle
      dq1 = xmom_centroid_values[k1] - xmom_centroid_values[k0];
      dq2 = xmom_centroid_values[k2] - xmom_centroid_values[k0];
      
      // Calculate the gradient of xmom on the auxiliary triangle
      a = dy2*dq1 - dy1*dq2;
      a *= inv_area2;
      b = dx1*dq2 - dx2*dq1;
      b *= inv_area2;
      
      // Calculate provisional jumps in stage from the centroid 
      // of triangle k to its vertices, to be limited      
      dqv[0] = a*dxv0+b*dyv0;
      dqv[1] = a*dxv1+b*dyv1;
      dqv[2] = a*dxv2+b*dyv2;
      
      // Now we want to find min and max of the centroid and the 
      // vertices of the auxiliary triangle and compute jumps 
      // from the centroid to the min and max
      //
      find_qmin_and_qfmax(dq0, dq1, dq2, &qmin, &qmax);

      beta_tmp = beta_uh_dry + (beta_uh - beta_uh_dry) * hfactor;

      limit_gradient(dqv, qmin, qmax, beta_tmp);
      

      for (i=0; i < 3; i++)
      {
        xmom_edge_values[k3+i] = xmom_centroid_values[k] + dqv[i];
      }
      
      //-----------------------------------
      // ymomentum
      //-----------------------------------                  

      // Calculate the difference between vertex 0 of the auxiliary 
      // triangle and the centroid of triangle k
      dq0 = ymom_centroid_values[k0] - ymom_centroid_values[k];
      
      // Calculate differentials between the vertices 
      // of the auxiliary triangle
      dq1 = ymom_centroid_values[k1] - ymom_centroid_values[k0];
      dq2 = ymom_centroid_values[k2] - ymom_centroid_values[k0];
      
      // Calculate the gradient of xmom on the auxiliary triangle
      a = dy2*dq1 - dy1*dq2;
      a *= inv_area2;
      b = dx1*dq2 - dx2*dq1;
      b *= inv_area2;
       
      // Calculate provisional jumps in stage from the centroid 
      // of triangle k to its vertices, to be limited
      dqv[0] = a*dxv0 + b*dyv0;
      dqv[1] = a*dxv1 + b*dyv1;
      dqv[2] = a*dxv2 + b*dyv2;
      
      // Now we want to find min and max of the centroid and the 
      // vertices of the auxiliary triangle and compute jumps 
      // from the centroid to the min and max
      //
      find_qmin_and_qfmax(dq0, dq1, dq2, &qmin, &qmax);
      
      beta_tmp = beta_vh_dry + (beta_vh - beta_vh_dry) * hfactor;   

      limit_gradient(dqv, qmin, qmax, beta_tmp);
      
      for (i=0;i<3;i++)
      {
        ymom_edge_values[k3 + i] = ymom_centroid_values[k] + dqv[i];
      }
      
    } // End number_of_boundaries <=1 
    else
    {

      //==============================================
      // Number of boundaries == 2
      //==============================================        
        
      // One internal neighbour and gradient is in direction of the neighbour's centroid
      
      // Find the only internal neighbour (k1?)
      for (k2 = k3; k2 < k3 + 3; k2++)
      {
      // Find internal neighbour of triangle k      
      // k2 indexes the edges of triangle k   
      
          if (surrogate_neighbours[k2] != k)
          {
             break;
          }
      }
      
      if ((k2 == k3 + 3)) 
      {
        // If we didn't find an internal neighbour
        // report_python_error(AT, "Internal neighbour not found");
        return -1;
      }
      
      k1 = surrogate_neighbours[k2];
      
      // The coordinates of the triangle are already (x,y). 
      // Get centroid of the neighbour (x1,y1)
      coord_index = 2*k1;
      x1 = centroid_coordinates[coord_index];
      y1 = centroid_coordinates[coord_index + 1];
      
      // Compute x- and y- distances between the centroid of 
      // triangle k and that of its neighbour
      dx1 = x1 - x; 
      dy1 = y1 - y;
      
      // Set area2 as the square of the distance
      area2 = dx1*dx1 + dy1*dy1;
      
      // Set dx2=(x1-x0)/((x1-x0)^2+(y1-y0)^2) 
      // and dy2=(y1-y0)/((x1-x0)^2+(y1-y0)^2) which
      // respectively correspond to the x- and y- gradients 
      // of the conserved quantities
      dx2 = 1.0/area2;
      dy2 = dx2*dy1;
      dx2 *= dx1;
      
      
      //-----------------------------------
      // stage
      //-----------------------------------            

      // Compute differentials
      dq1 = stage_centroid_values[k1] - stage_centroid_values[k];
      
      // Calculate the gradient between the centroid of triangle k 
      // and that of its neighbour
      a = dq1*dx2;
      b = dq1*dy2;
      
      // Calculate provisional edge jumps, to be limited
      dqv[0] = a*dxv0 + b*dyv0;
      dqv[1] = a*dxv1 + b*dyv1;
      dqv[2] = a*dxv2 + b*dyv2;
      
      // Now limit the jumps
      if (dq1>=0.0)
      {
        qmin=0.0;
        qmax=dq1;
      }
      else
      {
        qmin = dq1;
        qmax = 0.0;
      }
      
      // Limit the gradient
      limit_gradient(dqv, qmin, qmax, beta_w);
      
      //for (i=0; i < 3; i++)
      //{
      stage_edge_values[k3] = stage_centroid_values[k] + dqv[0];
      stage_edge_values[k3 + 1] = stage_centroid_values[k] + dqv[1];
      stage_edge_values[k3 + 2] = stage_centroid_values[k] + dqv[2];
      //}
      
      //-----------------------------------
      // xmomentum
      //-----------------------------------                        
      
      // Compute differentials
      dq1 = xmom_centroid_values[k1] - xmom_centroid_values[k];
      
      // Calculate the gradient between the centroid of triangle k 
      // and that of its neighbour
      a = dq1*dx2;
      b = dq1*dy2;
      
      // Calculate provisional edge jumps, to be limited
      dqv[0] = a*dxv0+b*dyv0;
      dqv[1] = a*dxv1+b*dyv1;
      dqv[2] = a*dxv2+b*dyv2;
      
      // Now limit the jumps
      if (dq1 >= 0.0)
      {
        qmin = 0.0;
        qmax = dq1;
      }
      else
      {
        qmin = dq1;
        qmax = 0.0;
      }
      
      // Limit the gradient      
      limit_gradient(dqv, qmin, qmax, beta_w);
      
      //for (i=0;i<3;i++)
      //xmom_edge_values[k3] = xmom_centroid_values[k] + dqv[0];
      //xmom_edge_values[k3 + 1] = xmom_centroid_values[k] + dqv[1];
      //xmom_edge_values[k3 + 2] = xmom_centroid_values[k] + dqv[2];
      
      for (i = 0; i < 3;i++)
      {
          xmom_edge_values[k3 + i] = xmom_centroid_values[k] + dqv[i];
      }
      
      //-----------------------------------
      // ymomentum
      //-----------------------------------                        

      // Compute differentials
      dq1 = ymom_centroid_values[k1] - ymom_centroid_values[k];
      
      // Calculate the gradient between the centroid of triangle k 
      // and that of its neighbour
      a = dq1*dx2;
      b = dq1*dy2;
      
      // Calculate provisional edge jumps, to be limited
      dqv[0] = a*dxv0 + b*dyv0;
      dqv[1] = a*dxv1 + b*dyv1;
      dqv[2] = a*dxv2 + b*dyv2;
      
      // Now limit the jumps
      if (dq1>=0.0)
      {
        qmin = 0.0;
        qmax = dq1;
      } 
      else
      {
        qmin = dq1;
        qmax = 0.0;
      }
      
      // Limit the gradient
      limit_gradient(dqv, qmin, qmax, beta_w);
      
      for (i=0;i<3;i++)
              {
              ymom_edge_values[k3 + i] = ymom_centroid_values[k] + dqv[i];
              }
    } // else [number_of_boundaries==2]
  } // for k=0 to number_of_elements-1
  
  // Compute vertex values of quantities
  for (k=0; k<number_of_elements; k++){
      k3=3*k;

      // Compute stage vertex values 
      stage_vertex_values[k3] = stage_edge_values[k3+1] + stage_edge_values[k3+2] -stage_edge_values[k3] ; 
      stage_vertex_values[k3+1] =  stage_edge_values[k3] + stage_edge_values[k3+2]-stage_edge_values[k3+1]; 
      stage_vertex_values[k3+2] =  stage_edge_values[k3] + stage_edge_values[k3+1]-stage_edge_values[k3+2]; 
      
     // Compute xmom vertex values 
      xmom_vertex_values[k3] = xmom_edge_values[k3+1] + xmom_edge_values[k3+2] -xmom_edge_values[k3] ; 
      xmom_vertex_values[k3+1] =  xmom_edge_values[k3] + xmom_edge_values[k3+2]-xmom_edge_values[k3+1]; 
      xmom_vertex_values[k3+2] =  xmom_edge_values[k3] + xmom_edge_values[k3+1]-xmom_edge_values[k3+2]; 

      // Compute ymom vertex values 
      ymom_vertex_values[k3] = ymom_edge_values[k3+1] + ymom_edge_values[k3+2] -ymom_edge_values[k3] ; 
      ymom_vertex_values[k3+1] =  ymom_edge_values[k3] + ymom_edge_values[k3+2]-ymom_edge_values[k3+1]; 
      ymom_vertex_values[k3+2] =  ymom_edge_values[k3] + ymom_edge_values[k3+1]-ymom_edge_values[k3+2]; 

      // If needed, convert from velocity to momenta
      if(extrapolate_velocity_second_order==1){
          //Convert velocity back to momenta at centroids
          xmom_centroid_values[k] = xmom_centroid_store[k];
          ymom_centroid_values[k] = ymom_centroid_store[k];

          // Re-compute momenta at edges
          for (i=0; i<3; i++){
              de[i] = fmax(stage_edge_values[k3+i]-elevation_edge_values[k3+i],0.0);
              xmom_edge_values[k3+i]=xmom_edge_values[k3+i]*de[i];
              ymom_edge_values[k3+i]=ymom_edge_values[k3+i]*de[i];
          }

          // Re-compute momenta at vertices
          for (i=0; i<3; i++){
              de[i] = fmax(stage_vertex_values[k3+i]-elevation_vertex_values[k3+i],0.0);
              xmom_vertex_values[k3+i]=xmom_vertex_values[k3+i]*de[i];
              ymom_vertex_values[k3+i]=ymom_vertex_values[k3+i]*de[i];
          }

      }

   
  } 

  free(xmom_centroid_store);
  free(ymom_centroid_store);
  free(stage_centroid_store);
  //free(min_elevation_edgevalue);
  free(max_elevation_edgevalue);
  //free(count_wet_neighbours);
  return 0;
}           


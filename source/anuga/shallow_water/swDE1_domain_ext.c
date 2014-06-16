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

//// Function to obtain speed from momentum and depth.
//// This is used by flux functions
//// Input parameters uh and h may be modified by this function.
//// Tried to inline, but no speedup was achieved 27th May 2009 (Ole)
////static inline double _compute_speed(double *uh, 
//double _compute_speed(double *uh, 
//		      double *h, 
//		      double epsilon, 
//		      double h0,
//		      double limiting_threshold) {
//  
//  double u;
//
//  if (*h < limiting_threshold) {   
//    // Apply limiting of speeds according to the ANUGA manual
//    if (*h < epsilon) {
//      //*h = max(0.0,*h);  // Could have been negative
//      u = 0.0;
//    } else {
//      u = *uh/(*h + h0/ *h);    
//    }
//  
//
//    // Adjust momentum to be consistent with speed
//    *uh = u * *h;
//  } else {
//    // We are in deep water - no need for limiting
//    u = *uh/ *h;
//  }
//  
//  return u;
//}
//
//// minmod limiter
//int _minmod(double a, double b){
//    // Compute minmod
//
//    if(sign(a)!=sign(b)){
//        return 0.0;
//    }else{
//        return min(fabs(a), fabs(b))*sign(a);
//    }
//
//
//}


// Innermost flux function (using stage w=z+h)
int _flux_function_toro(double *q_left, double *q_right,
                           double h_left, double h_right,
                           double hle, double hre,
                           double n1, double n2,
                           double epsilon,
                           double ze,
                           double limiting_threshold,
                           double g,
                           double *edgeflux, double *max_speed,
                           double *pressure_flux, double hc,
                           double hc_n, double speed_max_last)
{

  /*Compute fluxes between volumes for the shallow water wave equation
    cast in terms of the 'stage', w = h+z using

	HLL scheme of Fraccarollo and Toro Experimental and numerical assessment of the shallow
water model for two-dimensional dam-break type. Journal of Computational Physics,
33(6):843–864, 1995.

    FIXME: Several variables in this interface are no longer used, clean up
  */

  int i;

  double w_left,  uh_left, vh_left, u_left;
  double w_right, uh_right, vh_right, u_right;
  double s_min, s_max, soundspeed_left, soundspeed_right;
  double u_m, h_m, soundspeed_m, s_m;
  double denom, inverse_denominator;
  double uint, t1, t2, t3, min_speed, tmp;
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


  // Compute speeds in x-direction
  //w_left = q_left_rotated[0];
  uh_left=q_left_rotated[1];
  vh_left=q_left_rotated[2];
  if(hle>0.0){
    u_left = uh_left/hle ; //max(h_left, 1.0e-06);
    uh_left=h_left*u_left;
    vh_left=h_left/(hle)*vh_left;
  }else{
    u_left=0.;
    uh_left=0.;
    vh_left=0.;
  }

  //u_left = _compute_speed(&uh_left, &hle,
  //            epsilon, h0, limiting_threshold);

  //w_right = q_right_rotated[0];
  uh_right = q_right_rotated[1];
  vh_right = q_right_rotated[2];
  if(hre>0.0){
    u_right = uh_right/hre;//max(h_right, 1.0e-06);
    uh_right=h_right*u_right;
    vh_right=h_right/hre*vh_right;
  }else{
    u_right=0.;
    uh_right=0.;
    vh_right=0.;
  }
  //u_right = _compute_speed(&uh_right, &hre,
  //              epsilon, h0, limiting_threshold);



  // Maximal and minimal wave speeds
  soundspeed_left  = sqrt(g*h_left);
  soundspeed_right = sqrt(g*h_right);

  // Toro for shallow water
  u_m = 0.5*(u_left + u_right) + soundspeed_left - soundspeed_right;
  h_m = (u_left + 2.0*soundspeed_left - u_right + 2.0*soundspeed_right);
  h_m = h_m*h_m/(16.0*g);
  soundspeed_m = sqrt(g*h_m);


  if(h_left < 1.0e-10){
	  s_min = u_right - 2.0*soundspeed_right;
	  s_max = u_right + soundspeed_right;
	  s_m   = s_min;
  }
  else if (h_right < 1.0e-10){
	  s_min = u_left - soundspeed_left;
	  s_max = u_left + 2.0*soundspeed_left;
	  s_m = s_max;
  }
  else {
	  s_max = max(u_right + soundspeed_right, u_m + soundspeed_right);
	  s_min = min(u_left - soundspeed_left, u_m - soundspeed_m);
  }

  if (s_max < 0.0)
  {
    s_max = 0.0;
  }

  if (s_min > 0.0)
  {
    s_min = 0.0;
  }


  // Flux formulas
  flux_left[0] = u_left*h_left;
  flux_left[1] = u_left*uh_left; //+ 0.5*g*h_left*h_left;
  flux_left[2] = u_left*vh_left;

  flux_right[0] = u_right*h_right;
  flux_right[1] = u_right*uh_right ; //+ 0.5*g*h_right*h_right;
  flux_right[2] = u_right*vh_right;

  // Flux computation
  denom = s_max - s_min;
  if (denom < epsilon)
  {
    // Both wave speeds are very small
    memset(edgeflux, 0, 3*sizeof(double));

    *max_speed = 0.0;
    //*pressure_flux = 0.0;
    *pressure_flux = 0.5*g*0.5*(h_left*h_left+h_right*h_right);
  }
  else
  {
    // Maximal wavespeed
    *max_speed = max(s_max, -s_min);

    inverse_denominator = 1.0/max(denom,1.0e-100);
    for (i = 0; i < 3; i++)
    {
      edgeflux[i] = s_max*flux_left[i] - s_min*flux_right[i];

      // Standard smoothing term
      // edgeflux[i] += 1.0*(s_max*s_min)*(q_right_rotated[i] - q_left_rotated[i]);
      // Smoothing by stage alone can cause high velocities / slow draining for nearly dry cells
      if(i==0) edgeflux[i] += (s_max*s_min)*(max(q_right_rotated[i],ze) - max(q_left_rotated[i],ze));
      if(i==1) edgeflux[i] += (s_max*s_min)*(uh_right - uh_left);
      if(i==2) edgeflux[i] += (s_max*s_min)*(vh_right - vh_left);

      edgeflux[i] *= inverse_denominator;
    }
    // Separate pressure flux, so we can apply different wet-dry hacks to it
    *pressure_flux = 0.5*g*( s_max*h_left*h_left -s_min*h_right*h_right)*inverse_denominator;


    // Rotate back
    _rotate(edgeflux, n1, -n2);
  }

  return 0;
}





// Innermost flux function (using stage w=z+h)
int _flux_function_central(double *q_left, double *q_right,
                           double h_left, double h_right,
                           double hle, double hre,
                           double n1, double n2,
                           double epsilon, 
                           double ze,
                           double limiting_threshold, 
                           double g,
                           double *edgeflux, double *max_speed, 
                           double *pressure_flux, double hc, 
                           double hc_n, double speed_max_last) 
{

  /*Compute fluxes between volumes for the shallow water wave equation
    cast in terms of the 'stage', w = h+z using
    the 'central scheme' as described in

    Kurganov, Noelle, Petrova. 'Semidiscrete Central-Upwind Schemes For
    Hyperbolic Conservation Laws and Hamilton-Jacobi Equations'.
    Siam J. Sci. Comput. Vol. 23, No. 3, pp. 707-740.

    The implemented formula is given in equation (3.15) on page 714

    FIXME: Several variables in this interface are no longer used, clean up
  */

  int i;

  double w_left,  uh_left, vh_left, u_left;
  double w_right, uh_right, vh_right, u_right;
  double s_min, s_max, soundspeed_left, soundspeed_right;
  double denom, inverse_denominator;
  double uint, t1, t2, t3, min_speed, tmp;
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


  // Compute speeds in x-direction
  //w_left = q_left_rotated[0];          
  uh_left=q_left_rotated[1];
  vh_left=q_left_rotated[2];
  if(hle>0.0){
    u_left = uh_left/hle ; //max(h_left, 1.0e-06);
    uh_left=h_left*u_left;
    vh_left=h_left/(hle)*vh_left;
  }else{
    u_left=0.;
    uh_left=0.;
    vh_left=0.;
  }
  
  //u_left = _compute_speed(&uh_left, &hle, 
  //            epsilon, h0, limiting_threshold);

  //w_right = q_right_rotated[0];
  uh_right = q_right_rotated[1];
  vh_right = q_right_rotated[2];
  if(hre>0.0){
    u_right = uh_right/hre;//max(h_right, 1.0e-06);
    uh_right=h_right*u_right;
    vh_right=h_right/hre*vh_right;
  }else{
    u_right=0.;
    uh_right=0.;
    vh_right=0.;
  }
  //u_right = _compute_speed(&uh_right, &hre, 
  //              epsilon, h0, limiting_threshold);

  
  // Maximal and minimal wave speeds
  soundspeed_left  = sqrt(g*h_left);
  soundspeed_right = sqrt(g*h_right);  
  //soundspeed_left  = sqrt(g*hle);
  //soundspeed_right = sqrt(g*hre);  
  
  s_max = max(u_left + soundspeed_left, u_right + soundspeed_right);
  if (s_max < 0.0) 
  {
    s_max = 0.0;
  }

  //if( hc < 1.0e-03){
  //  s_max = 0.0;
  //}


  s_min = min(u_left - soundspeed_left, u_right - soundspeed_right);
  if (s_min > 0.0)
  {
    s_min = 0.0;
  }
  
  //if( hc_n < 1.0e-03){
  //  s_min = 0.0;
  //}
  
  // Flux formulas
  flux_left[0] = u_left*h_left;
  flux_left[1] = u_left*uh_left; //+ 0.5*g*h_left*h_left;
  flux_left[2] = u_left*vh_left;

  flux_right[0] = u_right*h_right;
  flux_right[1] = u_right*uh_right ; //+ 0.5*g*h_right*h_right;
  flux_right[2] = u_right*vh_right;

  // Flux computation
  denom = s_max - s_min;
  if (denom < epsilon) 
  { 
    // Both wave speeds are very small
    memset(edgeflux, 0, 3*sizeof(double));

    *max_speed = 0.0;
    //*pressure_flux = 0.0;
    *pressure_flux = 0.5*g*0.5*(h_left*h_left+h_right*h_right);
  } 
  else 
  {
    // Maximal wavespeed
    *max_speed = max(s_max, -s_min);
  
    inverse_denominator = 1.0/max(denom,1.0e-100);
    for (i = 0; i < 3; i++) 
    {
      edgeflux[i] = s_max*flux_left[i] - s_min*flux_right[i];

      // Standard smoothing term
      // edgeflux[i] += 1.0*(s_max*s_min)*(q_right_rotated[i] - q_left_rotated[i]);
      // Smoothing by stage alone can cause high velocities / slow draining for nearly dry cells
      if(i==0) edgeflux[i] += (s_max*s_min)*(max(q_right_rotated[i],ze) - max(q_left_rotated[i],ze));
      if(i==1) edgeflux[i] += (s_max*s_min)*(uh_right - uh_left);
      if(i==2) edgeflux[i] += (s_max*s_min)*(vh_right - vh_left);

      edgeflux[i] *= inverse_denominator;
    }
    // Separate pressure flux, so we can apply different wet-dry hacks to it
    *pressure_flux = 0.5*g*( s_max*h_left*h_left -s_min*h_right*h_right)*inverse_denominator;
    

    // Rotate back
    _rotate(edgeflux, n1, -n2);
  }
  
  return 0;
}

double adjust_edgeflux_with_weir(double* edgeflux,
                                 double h_left, double h_right, 
                                 double g, double weir_height,
                                 double weirScale, 
                                 double s1, double s2, 
                                 double h1, double h2
                                ){
    // Adjust the edgeflux to agree with a weir relation [including
    // subergence], but smoothly vary to shallow water solution when
    // the flow over the weir is much deeper than the weir, or the
    // upstream/downstream water elevations are too similar
    double rw, rw2; // 'Raw' weir fluxes
    double rwRat, hdRat,hdWrRat, scaleFlux, scaleFluxS, minhd, maxhd;
    double w1,w2; // Weights for averaging
    double newFlux;
    // Following constants control the 'blending' with the shallow water solution
    // They are now user-defined 
    //double s1=0.9; // At this submergence ratio, begin blending with shallow water solution
    //double s2=0.95; // At this submergence ratio, completely use shallow water solution
    //double h1=1.0; // At this (tailwater height above weir) / (weir height) ratio, begin blending with shallow water solution
    //double h2=1.5; // At this (tailwater height above weir) / (weir height) ratio, completely use the shallow water solution

    minhd=min(h_left, h_right);
    maxhd=max(h_left, h_right);
    // 'Raw' weir discharge = weirScale*2/3*H*(2/3*g*H)**0.5
    rw=weirScale*2./3.*maxhd*sqrt(2./3.*g*maxhd);
    // Factor for villemonte correction
    rw2=weirScale*2./3.*minhd*sqrt(2./3.*g*minhd);
    // Useful ratios
    rwRat=rw2/max(rw, 1.0e-100);
    hdRat=minhd/max(maxhd,1.0e-100);
    // (tailwater height above weir)/weir_height ratio
    hdWrRat=minhd/max(weir_height,1.0e-100);
        
    // Villemonte (1947) corrected weir flow with submergence
    // Q = Q1*(1-Q2/Q1)**0.385
    rw = rw*pow(1.0-rwRat,0.385);

    if(h_right>h_left){
        rw*=-1.0;
    }

    //printf("%e, %e \n", rw, edgeflux[0]);

    if( (hdRat<s2) & (hdWrRat< h2) ){
        // Rescale the edge fluxes so that the mass flux = desired flux
        // Linearly shift to shallow water solution between hdRat = s1 and s2  
        // and between hdWrRat = h1 and h2

        //        
        // WEIGHT WITH RAW SHALLOW WATER FLUX BELOW
        // This ensures that as the weir gets very submerged, the 
        // standard shallow water equations smoothly take over
        //

        // Weighted average constants to transition to shallow water eqn flow
        w1=min( max(hdRat-s1,0.)/(s2-s1), 1.0);
        //w2=1.0-w1;
        //
        // Adjust again when the head is too deep relative to the weir height
        w2=min( max(hdWrRat-h1,0.)/(h2-h1), 1.0);
        //w2=1.0-w1;

        newFlux=(rw*(1.0-w1)+w1*edgeflux[0])*(1.0-w2) + w2*edgeflux[0];
       
        if(fabs(edgeflux[0])>1.0e-100){ 
            scaleFlux=newFlux/edgeflux[0];

            // FINAL ADJUSTED FLUX
            edgeflux[0]*=scaleFlux;
            edgeflux[1]*=scaleFlux;
            edgeflux[2]*=scaleFlux;
        }else{
            // Can't divide by edgeflux
            edgeflux[0] = newFlux;
            edgeflux[1]*=0.;
            edgeflux[2]*=0.;
        }
    }

    //printf("%e, %e, %e, %e , %e, %e \n", weir_height, h_left, h_right, edgeflux[0], w1, w2);

    return 0;
}

// Computational function for flux computation
double _compute_fluxes_central(int number_of_elements,
        double timestep,
        double epsilon,
        double H0,
        double g,
        double* boundary_flux_sum,
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
        double* height_edge_values,
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
        int timestep_fluxcalls,
        double* stage_centroid_values,
        double* xmom_centroid_values,
        double* ymom_centroid_values,
        double* bed_centroid_values,
        double* height_centroid_values,
        //double* bed_vertex_values,
        long* edge_flux_type,
        double* riverwall_elevation,
        long* riverwall_rowIndex,
        int ncol_riverwall_hydraulic_properties,
        double* riverwall_hydraulic_properties) {
    // Local variables
    double max_speed, length, inv_area, zl, zr;
    //double h0 = H0*H0;//H0*H0;//H0*0.1;//H0*H0;//H0*H0; // This ensures a good balance when h approaches H0.
    double h_left, h_right, z_half ;  // For andusse scheme
    // FIXME: limiting_threshold is not used for DE1
    double limiting_threshold = 10*H0;//10 * H0 ;//H0; //10 * H0; // Avoid applying limiter below this
    // threshold for performance reasons.
    // See ANUGA manual under flux limiting
    int k, i, m, n,j, ii;
    int ki, nm = 0, ki2,ki3, nm3; // Index shorthands
    //int num_hydraulic_properties=1; // FIXME: Move to function call
    // Workspace (making them static actually made function slightly slower (Ole))
    double ql[3], qr[3], edgeflux[3]; // Work array for summing up fluxes
    double stage_edges[3];//Work array
    double bedslope_work, local_timestep;
    int neighbours_wet[3];//Work array
    int RiverWall_count;
    double hle, hre, zc, zc_n, Qfactor, s1, s2, h1, h2; 
    double stage_edge_lim, outgoing_mass_edges, bedtop, bedbot, pressure_flux, hc, hc_n, tmp, tmp2;
    static long call = 1; // Static local variable flagging already computed flux
    double speed_max_last, vol, smooth, local_speed, weir_height;

    double *edgeflux_store, *pressuregrad_store;

    edgeflux_store = malloc(number_of_elements*9*sizeof(double));
    pressuregrad_store = malloc(number_of_elements*3*sizeof(double));

    call++; // Flag 'id' of flux calculation for this timestep

    // Set explicit_update to zero for all conserved_quantities.
    // This assumes compute_fluxes called before forcing terms
    memset((char*) stage_explicit_update, 0, number_of_elements * sizeof (double));
    memset((char*) xmom_explicit_update, 0, number_of_elements * sizeof (double));
    memset((char*) ymom_explicit_update, 0, number_of_elements * sizeof (double));

    local_timestep=timestep;
    RiverWall_count=0;


    // For all triangles
    for (k = 0; k < number_of_elements; k++) {
        speed_max_last=0.0;
        // Loop through neighbours and compute edge flux for each
        for (i = 0; i < 3; i++) {
            ki = k * 3 + i; // Linear index to edge i of triangle k
            ki2 = 2 * ki; //k*6 + i*2
            ki3 = 3*ki; 

            if (already_computed_flux[ki] == call) {
                // We've already computed the flux across this edge
                continue;
            }

            // Get left hand side values from triangle k, edge i
            ql[0] = stage_edge_values[ki];
            ql[1] = xmom_edge_values[ki];
            ql[2] = ymom_edge_values[ki];
            zl = bed_edge_values[ki];
            hc = height_centroid_values[k];
            zc = bed_centroid_values[k];
            hle= height_edge_values[ki];

            // Get right hand side values either from neighbouring triangle
            // or from boundary array (Quantities at neighbour on nearest face).
            n = neighbours[ki];
            hc_n = hc;
            zc_n = bed_centroid_values[k];
            if (n < 0) {
                // Neighbour is a boundary condition
                m = -n - 1; // Convert negative flag to boundary index

                qr[0] = stage_boundary_values[m];
                qr[1] = xmom_boundary_values[m];
                qr[2] = ymom_boundary_values[m];
                zr = zl; // Extend bed elevation to boundary
                hre= max(qr[0]-zr,0.);//hle; 
            } else {
                // Neighbour is a real triangle
                hc_n = height_centroid_values[n];
                zc_n = bed_centroid_values[n];
                m = neighbour_edges[ki];
                nm = n * 3 + m; // Linear index (triangle n, edge m)
                nm3 = nm*3;

                qr[0] = stage_edge_values[nm];
                qr[1] = xmom_edge_values[nm];
                qr[2] = ymom_edge_values[nm];
                zr = bed_edge_values[nm];
                hre = height_edge_values[nm];
            }
          
            // Audusse magic 
            z_half=max(zl,zr);

            // Account for riverwalls
            if(edge_flux_type[ki]==1){
                // Update counter of riverwall edges == index of
                // riverwall_elevation + riverwall_rowIndex
                RiverWall_count+=1;
                
                // Set central bed to riverwall elevation
                z_half=max(riverwall_elevation[RiverWall_count-1], max(zl, zr)) ;

                if(min(ql[0], qr[0]) < z_half){
                    // Since there is a wall blocking the flow connection, use first order extrapolation for this edge
                    ql[0]=stage_centroid_values[k];
                    ql[1]=xmom_centroid_values[k];
                    ql[2]=ymom_centroid_values[k];
                    hle=hc;
                    zl=zc;

                    if(n>=0){
                      qr[0]=stage_centroid_values[n];
                      qr[1]=xmom_centroid_values[n];
                      qr[2]=ymom_centroid_values[n];
                      hre=hc_n;
                      zr = zc_n;
                    }else{
                      hre=hc;
                      zr = zc;
                    }
                    // Re-set central bed to riverwall elevation
                    z_half=max(riverwall_elevation[RiverWall_count-1], max(zl, zr)) ;
                }
                

            }

            // Define h left/right for Audusse flux method
            h_left= max(hle+zl-z_half,0.);
            h_right= max(hre+zr-z_half,0.);

            // Edge flux computation (triangle k, edge i)
            _flux_function_central(ql, qr,
                    h_left, h_right,
                    hle, hre,
                    normals[ki2], normals[ki2 + 1],
                    epsilon, z_half, limiting_threshold, g,
                    edgeflux, &max_speed, &pressure_flux, hc, hc_n, 
                    speed_max_last);

            // Force weir discharge to match weir theory
            if(edge_flux_type[ki]==1){
                //printf("%e \n", z_half);
                weir_height=max(riverwall_elevation[RiverWall_count-1]-min(zl,zr), 0.); // Reference weir height  
                //weir_height=max(z_half-max(zl,zr), 0.); // Reference weir height  
                // If the weir height is zero, avoid the weir compuattion entirely
                if(weir_height>0.){
                    // Get Qfactor index - multiply the idealised weir discharge by this constant factor
                    ii=riverwall_rowIndex[RiverWall_count-1]*ncol_riverwall_hydraulic_properties;
                    Qfactor=riverwall_hydraulic_properties[ii];
                    //printf("%e \n", Qfactor);
                    // Get s1, submergence ratio at which we start blending with the shallow water solution 
                    ii+=1;
                    s1=riverwall_hydraulic_properties[ii];
                    // Get s2, submergence ratio at which we entirely use the shallow water solution 
                    ii+=1;
                    s2=riverwall_hydraulic_properties[ii];
                    // Get h1, tailwater head / weir height at which we start blending with the shallow water solution 
                    ii+=1;
                    h1=riverwall_hydraulic_properties[ii];
                    // Get h2, tailwater head / weir height at which we entirely use the shallow water solution 
                    ii+=1;
                    h2=riverwall_hydraulic_properties[ii];

                    //printf("%e, %e, %e, %e, %e \n", Qfactor, s1, s2, h1, h2);

                    adjust_edgeflux_with_weir(edgeflux, h_left, h_right, g, 
                                              weir_height, Qfactor, 
                                              s1, s2, h1, h2);
                    // NOTE: Should perhaps also adjust the wave-speed here?? Small chance of instability??
                    //printf("%e \n", edgeflux[0]);
                }
            }
            
            // Multiply edgeflux by edgelength
            length = edgelengths[ki];
            edgeflux[0] *= length;
            edgeflux[1] *= length;
            edgeflux[2] *= length;


            //// Don't allow an outward advective flux if the cell centroid stage
            //// is < the edge value. Is this important (??)
            //if((hc<H0) && edgeflux[0] > 0.){
            //    edgeflux[0] = 0.;
            //    edgeflux[1] = 0.;
            //    edgeflux[2] = 0.;
            //    //max_speed=0.;
            //    //pressure_flux=0.;
            //}
            ////
            //if((hc_n<H0) && edgeflux[0] < 0.){
            //    edgeflux[0] = 0.;
            //    edgeflux[1] = 0.;
            //    edgeflux[2] = 0.;
            //    //max_speed=0.;
            //    //pressure_flux=0.;
            //}

            edgeflux_store[ki3 + 0 ] = -edgeflux[0];
            edgeflux_store[ki3 + 1 ] = -edgeflux[1];
            edgeflux_store[ki3 + 2 ] = -edgeflux[2];

            // bedslope_work contains all gravity related terms -- weighting of
            bedslope_work=length*(-g*0.5*(h_left*h_left - hle*hle -(hle+hc)*(zl-zc))+pressure_flux);

            pressuregrad_store[ki]=bedslope_work;
            
            already_computed_flux[ki] = call; // #k Done

            // Update neighbour n with same flux but reversed sign
            if (n >= 0) {

                edgeflux_store[nm3 + 0 ] = edgeflux[0];
                edgeflux_store[nm3 + 1 ] = edgeflux[1];
                edgeflux_store[nm3 + 2 ] = edgeflux[2];
                bedslope_work=length*(-g*0.5*(h_right*h_right-hre*hre-(hre+hc_n)*(zr-zc_n))+pressure_flux);
                pressuregrad_store[nm]=bedslope_work;

                already_computed_flux[nm] = call; // #n Done
            }

            // Update timestep based on edge i and possibly neighbour n
            // NOTE: We should only change the timestep between rk2 or rk3
            // steps, NOT within them (since a constant timestep is used within
            // each rk2/rk3 sub-step)
            if ((tri_full_flag[k] == 1) & ( (call-2)%timestep_fluxcalls==0)) {

                speed_max_last=max(speed_max_last, max_speed);

                if (max_speed > epsilon) {
                    // Apply CFL condition for triangles joining this edge (triangle k and triangle n)

                    // CFL for triangle k
                    local_timestep = min(local_timestep, radii[k] / max_speed);

                    if (n >= 0) {
                        // Apply CFL condition for neigbour n (which is on the ith edge of triangle k)
                        local_timestep = min(local_timestep, radii[n] / max_speed);
                    }

                }
            }
        
        // Keep track of maximal speeds
        //max_speed_array[k] = max(max_speed_array[k],max_speed);

        } // End edge i (and neighbour n)
        // Keep track of maximal speeds
        if((call-2)%timestep_fluxcalls==0) max_speed_array[k] = speed_max_last; //max_speed;


    } // End triangle k
 
    //// GD HACK 
    //// Limit edgefluxes, for mass conservation near wet/dry cells
    //for(k=0; k< number_of_elements; k++){
    //    //continue;
    //    hc = height_centroid_values[k]; 
    //    // Loop over every edge
    //    for(i = 0; i<3; i++){
    //        if(i==0){
    //            // Add up the outgoing flux through the cell -- only do this once (i==0)
    //            outgoing_mass_edges=0.0;
    //            for(useint=0; useint<3; useint++){
    //                if(edgeflux_store[3*(3*k+useint)]< 0.){
    //                    //outgoing_mass_edges+=1.0;
    //                    outgoing_mass_edges+=(edgeflux_store[3*(3*k+useint)]);
    //                }
    //            }
    //            outgoing_mass_edges*=local_timestep;
    //        }

    //        ki=3*k+i;   
    //        ki2=ki*2;
    //        ki3 = ki*3;
    //        
    //        // Prevent outflow from 'seriously' dry cells
    //        // Idea: The cell will not go dry if:
    //        // total_outgoing_flux <= cell volume = Area_triangle*hc
    //        vol=areas[k]*hc;
    //        if((edgeflux_store[ki3]< 0.0) && (-outgoing_mass_edges> vol)){
    //            
    //            // This bound could be improved (e.g. we could actually sum the
    //            // + and - fluxes and check if they are too large).  However,
    //            // the advantage of this method is that we don't have to worry
    //            // about subsequent changes to the + edgeflux caused by
    //            // constraints associated with neighbouring triangles.
    //            tmp = vol/(-(outgoing_mass_edges)) ;
    //            if(tmp< 1.0){
    //                edgeflux_store[ki3+0]*=tmp;
    //                edgeflux_store[ki3+1]*=tmp;
    //                edgeflux_store[ki3+2]*=tmp;

    //                // Compute neighbour edge index
    //                n = neighbours[ki];
    //                if(n>=0){
    //                    nm = 3*n + neighbour_edges[ki];
    //                    nm3 = nm*3;
    //                    edgeflux_store[nm3+0]*=tmp;
    //                    edgeflux_store[nm3+1]*=tmp;
    //                    edgeflux_store[nm3+2]*=tmp;
    //                }
    //            }
    //        }
    //    }
    // }

    ////printf("%e \n", edgeflux_store[3*30*3]);

    // Now add up stage, xmom, ymom explicit updates
    for(k=0; k<number_of_elements; k++){
        hc = max(stage_centroid_values[k] - bed_centroid_values[k],0.);
        stage_explicit_update[k]=0.;
        xmom_explicit_update[k]=0.;
        ymom_explicit_update[k]=0.;

        for(i=0;i<3;i++){
            ki=3*k+i;   
            ki2=ki*2;
            ki3 = ki*3;
            n=neighbours[ki];

            // GD HACK
            // Option to limit advective fluxes
            //if(hc > H0){
                stage_explicit_update[k] += edgeflux_store[ki3+0];
                xmom_explicit_update[k] += edgeflux_store[ki3+1];
                ymom_explicit_update[k] += edgeflux_store[ki3+2];
            //}else{
            //    stage_explicit_update[k] += edgeflux_store[ki3+0];
            //}


            // If this cell is not a ghost, and the neighbour is a boundary
            // condition OR a ghost cell, then add the flux to the
            // boundary_flux_integral
            if( (n<0 & tri_full_flag[k]==1) | ( n>=0 && (tri_full_flag[k]==1 & tri_full_flag[n]==0)) ){ 
                boundary_flux_sum[0] += edgeflux_store[ki3];
            }

            // GD HACK
            // Compute bed slope term
            //if(hc > H0){
                xmom_explicit_update[k] -= normals[ki2]*pressuregrad_store[ki];
                ymom_explicit_update[k] -= normals[ki2+1]*pressuregrad_store[ki];
            //}else{
            //    xmom_explicit_update[k] *= 0.;
            //    ymom_explicit_update[k] *= 0.;
            //}

        } // end edge i

        // Normalise triangle k by area and store for when all conserved
        // quantities get updated
        inv_area = 1.0 / areas[k];
        stage_explicit_update[k] *= inv_area;
        xmom_explicit_update[k] *= inv_area;
        ymom_explicit_update[k] *= inv_area;
   
    }  // end cell k

    // Hack to ensure we only update the timestep on the first call within each rk2/rk3 step
    if((call-2)%timestep_fluxcalls==0) timestep=local_timestep; 

    free(edgeflux_store);
    free(pressuregrad_store);

    return timestep;
}

// Protect against the water elevation falling below the triangle bed
double  _protect(int N,
         double minimum_allowed_height,
         double maximum_allowed_speed,
         double epsilon,
         double* wc,
         double* wv,
         double* zc,
         double* zv,
         double* xmomc,
         double* ymomc,
         double* areas,
         double* xc,
         double* yc) {

  int k;
  double hc, bmin, bmax;
  double u, v, reduced_speed;
  double mass_error=0.;
  // This acts like minimum_allowed height, but scales with the vertical
  // distance between the bed_centroid_value and the max bed_edge_value of
  // every triangle.
  //double minimum_relative_height=0.05; 
  int mass_added=0;

  // Protect against inifintesimal and negative heights  
  //if (maximum_allowed_speed < epsilon) {
    for (k=0; k<N; k++) {
      hc = wc[k] - zc[k];
      if (hc < minimum_allowed_height*1.0 ){
            // Set momentum to zero and ensure h is non negative
            xmomc[k] = 0.;
            ymomc[k] = 0.;
        if (hc <= 0.0){
             bmin=zc[k];
             // Minimum allowed stage = bmin

             // WARNING: ADDING MASS if wc[k]<bmin
             if(wc[k]<bmin){
                 mass_error+=(bmin-wc[k])*areas[k];
                 mass_added=1; //Flag to warn of added mass                

                 wc[k] = bmin; 

                 // FIXME: Set vertex values as well. Seems that this shouldn't be
                 // needed. However, from memory this is important at the first
                 // time step, for 'dry' areas where the designated stage is
                 // less than the bed centroid value
                 wv[3*k] = min(bmin, wc[k]); //zv[3*k]-minimum_allowed_height);
                 wv[3*k+1] = min(bmin, wc[k]); //zv[3*k+1]-minimum_allowed_height);
                 wv[3*k+2] = min(bmin, wc[k]); //zv[3*k+2]-minimum_allowed_height);
            }
        }
      }
    }

  //if(mass_added==1){
  //  printf("Cumulative mass protection: %f m^3 \n", mass_error);
  //}
  return mass_error;
}

int find_qmin_and_qmax(double dq0, double dq1, double dq2, 
               double *qmin, double *qmax){
  // Considering the centroid of an FV triangle and the vertices of its 
  // auxiliary triangle, find 
  // qmin=min(q)-qc and qmax=max(q)-qc, 
  // where min(q) and max(q) are respectively min and max over the
  // four values (at the centroid of the FV triangle and the auxiliary 
  // triangle vertices),
  // and qc is the centroid
  // dq0=q(vertex0)-q(centroid of FV triangle)
  // dq1=q(vertex1)-q(vertex0)
  // dq2=q(vertex2)-q(vertex0)

  // This is a simple implementation 
  *qmax = max(max(dq0, max(dq0+dq1, dq0+dq2)), 0.0) ;
  *qmin = min(min(dq0, min(dq0+dq1, dq0+dq2)), 0.0) ;
 
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
  //return 0;
  
  for (i=0;i<3;i++){
    if (dqv[i]<-TINY)
      r0=qmin/dqv[i];
      
    if (dqv[i]>TINY)
      r0=qmax/dqv[i];
      
    r=min(r0,r);
  }
  
  phi=min(r*beta_w,1.0);
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
                                 long* neighbour_edges,
                                 long* number_of_boundaries,
                                 double* centroid_coordinates,
                                 double* stage_centroid_values,
                                 double* xmom_centroid_values,
                                 double* ymom_centroid_values,
                                 double* elevation_centroid_values,
                                 double* height_centroid_values,
                                 double* edge_coordinates,
                                 double* stage_edge_values,
                                 double* xmom_edge_values,
                                 double* ymom_edge_values,
                                 double* elevation_edge_values,
                                 double* height_edge_values,
                                 double* stage_vertex_values,
                                 double* xmom_vertex_values,
                                 double* ymom_vertex_values,
                                 double* elevation_vertex_values,
                                 double* height_vertex_values,
                                 int optimise_dry_cells, 
                                 int extrapolate_velocity_second_order) {
                  
  // Local variables
  double a, b; // Gradient vector used to calculate edge values from centroids
  int k, k0, k1, k2, k3, k6, coord_index, i, ii, ktmp, k_wetdry;
  double x, y, x0, y0, x1, y1, x2, y2, xv0, yv0, xv1, yv1, xv2, yv2; // Vertices of the auxiliary triangle
  double dx1, dx2, dy1, dy2, dxv0, dxv1, dxv2, dyv0, dyv1, dyv2, dq0, dq1, dq2, area2, inv_area2, dpth,momnorm;
  double dqv[3], qmin, qmax, hmin, hmax, bedmax,bedmin, stagemin;
  double hc, h0, h1, h2, beta_tmp, hfactor, xtmp, ytmp, weight, tmp;
  double dk, dv0, dv1, dv2, de[3], demin, dcmax, r0scale, vel_norm, l1, l2, a_tmp, b_tmp, c_tmp,d_tmp;
  
  double *xmom_centroid_store, *ymom_centroid_store; // , *min_elevation_edgevalue, *max_elevation_edgevalue;
  //double *cell_wetness_scale;
  //int *count_wet_neighbours;

  // Use malloc to avoid putting these variables on the stack, which can cause
  // segfaults in large model runs
  xmom_centroid_store = malloc(number_of_elements*sizeof(double));
  ymom_centroid_store = malloc(number_of_elements*sizeof(double));
  //min_elevation_edgevalue = malloc(number_of_elements*sizeof(double));
  //max_elevation_edgevalue = malloc(number_of_elements*sizeof(double));
  //cell_wetness_scale = malloc(number_of_elements*sizeof(double));
  //count_wet_neighbours = malloc(number_of_elements*sizeof(int));
 
  //
  // Compute some useful statistics on wetness/dryness
  //
  //for(k=0; k<number_of_elements;k++){ 
  //    //cell_wetness_scale[k] = 0.;
  //    //// Check if cell k is wet
  //    ////if(stage_centroid_values[k] > elevation_centroid_values[k]){
  //    //if(stage_centroid_values[k] > elevation_centroid_values[k] + 1.0*minimum_allowed_height){
  //    ////if(stage_centroid_values[k] > max_elevation_edgevalue[k] + minimum_allowed_height+epsilon){
  //    //    cell_wetness_scale[k] = 1.;  
  //    //}

  //    //min_elevation_edgevalue[k] = min(elevation_edge_values[3*k], 
  //    //                                 min(elevation_edge_values[3*k+1],
  //    //                                     elevation_edge_values[3*k+2]));
  //    //max_elevation_edgevalue[k] = max(elevation_edge_values[3*k], 
  //    //                                 max(elevation_edge_values[3*k+1],
  //    //                                     elevation_edge_values[3*k+2]));
  //}

  //// Alternative 'PROTECT' step
  //for(k=0; k<number_of_elements;k++){ 
  //  //if((cell_wetness_scale[k]==0. ) ){
  //  //    xmom_centroid_values[k]=0.;
  //  //    ymom_centroid_values[k]=0.;
  //  //    //xmom_centroid_values[k]*=0.9;
  //  //    //ymom_centroid_values[k]*=0.9;

  //  //}

  //  
  //  ////// Try some Froude-number limiting in shallow depths
  //  //dpth=max(stage_centroid_values[k] - elevation_centroid_values[k], 0.0);
  //  ////
  //  //if(dpth< max(max_elevation_edgevalue[k]-min_elevation_edgevalue[k],10*minimum_allowed_height)){
  //  //    // momnorm = momentum^2
  //  //    momnorm=(xmom_centroid_values[k]*xmom_centroid_values[k]+
  //  //             ymom_centroid_values[k]*ymom_centroid_values[k]);
  //  //    
  //  //    // vel^2 < constant*g*dpth [-- then multiply both sides by dpth^2]
  //  //    if(momnorm > 4*9.81*dpth*dpth*dpth){
  //  //        // Down-scale momentum so that Froude number < constant
  //  //        tmp=sqrt((4*9.81*dpth*dpth*dpth)/momnorm);
  //  //        xmom_centroid_values[k] *=tmp;
  //  //        ymom_centroid_values[k] *=tmp;
  //  //    }
  //  //}
  //}

  
  if(extrapolate_velocity_second_order==1){

      // Replace momentum centroid with velocity centroid to allow velocity
      // extrapolation This will be changed back at the end of the routine
      for (k=0; k<number_of_elements; k++){
          
          height_centroid_values[k] = max(stage_centroid_values[k] - elevation_centroid_values[k], 0.);

          dk = height_centroid_values[k]; 
          if(dk>minimum_allowed_height){
              xmom_centroid_store[k] = xmom_centroid_values[k];
              xmom_centroid_values[k] = xmom_centroid_values[k]/dk;

              ymom_centroid_store[k] = ymom_centroid_values[k];
              ymom_centroid_values[k] = ymom_centroid_values[k]/dk;
          }else{
              xmom_centroid_store[k] = 0.;
              xmom_centroid_values[k] = 0.;
              ymom_centroid_store[k] = 0.;
              ymom_centroid_values[k] = 0.;

         }
      }
  }

  // If a triangle is surrounded by dry cells (or dry cells + boundary
  // condition) set its momentum to zero too. This prevents 'pits' of
  // of water being trapped and unable to lose momentum, which can occur in
  // some situations
  for (k=0; k<number_of_elements;k++){

      k3=k*3;
      k0 = surrogate_neighbours[k3];
      k1 = surrogate_neighbours[k3 + 1];
      k2 = surrogate_neighbours[k3 + 2];

      if((height_centroid_values[k0] < minimum_allowed_height | k0==k) &
         (height_centroid_values[k1] < minimum_allowed_height | k1==k) &
         (height_centroid_values[k2] < minimum_allowed_height | k2==k)){
              xmom_centroid_store[k] = 0.;
              xmom_centroid_values[k] = 0.;
              ymom_centroid_store[k] = 0.;
              ymom_centroid_values[k] = 0.;

      }
  }

  // Begin extrapolation routine
  for (k = 0; k < number_of_elements; k++) 
  {

    // Useful indices
    k3=k*3;
    k6=k*6;

    if (number_of_boundaries[k]==3)
    {
      // No neighbours, set gradient on the triangle to zero
      
      stage_edge_values[k3]   = stage_centroid_values[k];
      stage_edge_values[k3+1] = stage_centroid_values[k];
      stage_edge_values[k3+2] = stage_centroid_values[k];

      //xmom_centroid_values[k] = 0.;
      //ymom_centroid_values[k] = 0.;
      
      xmom_edge_values[k3]    = xmom_centroid_values[k];
      xmom_edge_values[k3+1]  = xmom_centroid_values[k];
      xmom_edge_values[k3+2]  = xmom_centroid_values[k];
      ymom_edge_values[k3]    = ymom_centroid_values[k];
      ymom_edge_values[k3+1]  = ymom_centroid_values[k];
      ymom_edge_values[k3+2]  = ymom_centroid_values[k];

      dk = height_centroid_values[k];
      height_edge_values[k3] = dk;
      height_edge_values[k3+1] = dk;
      height_edge_values[k3+2] = dk;
      
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

    }


            
    if (number_of_boundaries[k]<=1)
    {
      //==============================================
      // Number of boundaries <= 1
      // 'Typical case'
      //==============================================    
    
    
      // If no boundaries, auxiliary triangle is formed 
      // from the centroids of the three neighbours
      // If one boundary, auxiliary triangle is formed 
      // from this centroid and its two neighbours
      
      k0 = surrogate_neighbours[k3];
      k1 = surrogate_neighbours[k3 + 1];
      k2 = surrogate_neighbours[k3 + 2];

      //if(cell_wetness_scale[k0]==0 && cell_wetness_scale[k1]==0 && cell_wetness_scale[k2]==0){
      //    xmom_centroid_store[k]=0.;
      //    ymom_centroid_store[k]=0.;
      //    xmom_centroid_values[k] = 0.;
      //    ymom_centroid_values[k] = 0.;
      //}
     
      // Test to see whether we accept the surrogate neighbours 
      // Note that if ki is replaced with k in more than 1 neighbour, then the
      // triangle area will be zero, and a first order extrapolation will be
      // used
      // FIXME: Remove cell_wetness_scale if you don't need it
      //if( (cell_wetness_scale[k2]==0.0 && stage_centroid_values[k]<max_elevation_edgevalue[k]+100.)){
      //    k2 = k ;
      //}
      //if((cell_wetness_scale[k0]==0.0 && stage_centroid_values[k]<max_elevation_edgevalue[k])+100.){
      //    k0 = k ;
      //}
      //if((cell_wetness_scale[k1]==0.0 && stage_centroid_values[k]<max_elevation_edgevalue[k])+100.){
      //    k1 = k ;
      //}

      // Take note if the max neighbour bed elevation is greater than the min
      // neighbour stage -- suggests a 'steep' bed relative to the flow
      //bedmax = max(elevation_centroid_values[k], 
      //             max(elevation_centroid_values[k0],
      //                 max(elevation_centroid_values[k1], 
      //                     elevation_centroid_values[k2])));
      //bedmin = min(elevation_centroid_values[k], 
      //             min(elevation_centroid_values[k0],
      //                 min(elevation_centroid_values[k1], 
      //                     elevation_centroid_values[k2])));
      //stagemin = min(max(stage_centroid_values[k], elevation_centroid_values[k]), 
      //               min(max(stage_centroid_values[k0], elevation_centroid_values[k0]),
      //                   min(max(stage_centroid_values[k1], elevation_centroid_values[k1]),
      //                       max(stage_centroid_values[k2], elevation_centroid_values[k2]))));

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
      //if(k==54) printf("K=54\n");
       
      //// Treat triangles with no neighbours (area2 <=0.) 
      if ((area2 <= 0.))//|( cell_wetness_scale[k2]==0. && cell_wetness_scale[k1]==0. && cell_wetness_scale[k0]==0.))//|| (cell_wetness_scale[k]==0.)) //|(count_wet_neighbours[k]==0))
      {


          // Isolated wet cell -- constant stage/depth extrapolation 
          stage_edge_values[k3]   = stage_centroid_values[k];
          stage_edge_values[k3+1] = stage_centroid_values[k];
          stage_edge_values[k3+2] = stage_centroid_values[k];

          dk=height_centroid_values[k]; //max(stage_centroid_values[k]-elevation_centroid_values[k],0.);
          height_edge_values[k3] = dk;
          height_edge_values[k3+1] = dk;
          height_edge_values[k3+2] = dk;
          
          //stage_edge_values[k3]   = elevation_centroid_values[k]+dk;
          //stage_edge_values[k3+1] = elevation_centroid_values[k]+dk;
          //stage_edge_values[k3+2] = elevation_centroid_values[k]+dk;
          //stage_edge_values[k3] = elevation_edge_values[k3]+dk;
          //stage_edge_values[k3+1] = elevation_edge_values[k3+1]+dk;
          //stage_edge_values[k3+2] = elevation_edge_values[k3+2]+dk;

          //xmom_centroid_values[k]=0.;
          //ymom_centroid_values[k]=0.;
          //xmom_centroid_store[k] = 0.;
          //ymom_centroid_store[k] = 0.;

          xmom_edge_values[k3]    = xmom_centroid_values[k];
          xmom_edge_values[k3+1]  = xmom_centroid_values[k];
          xmom_edge_values[k3+2]  = xmom_centroid_values[k];
          ymom_edge_values[k3]    = ymom_centroid_values[k];
          ymom_edge_values[k3+1]  = ymom_centroid_values[k];
          ymom_edge_values[k3+2]  = ymom_centroid_values[k];

          continue;
      }  
      
      // Calculate heights of neighbouring cells
      hc = height_centroid_values[k];//stage_centroid_values[k]  - elevation_centroid_values[k];
      h0 = height_centroid_values[k0];// - elevation_centroid_values[k0];
      h1 = height_centroid_values[k1];// - elevation_centroid_values[k1];
      h2 = height_centroid_values[k2];// - elevation_centroid_values[k2];
      
      hmin = min(min(h0, min(h1, h2)), hc);
      hmax = max(max(h0, max(h1, h2)), hc);

      // Look for strong changes in cell depth as an indicator of near-wet-dry
      // Reduce hfactor linearly from 1-0 between depth ratio (hmin/hc) of [a_tmp , b_tmp]
      // NOTE: If we have a more 'second order' treatment in near dry areas (e.g. with b_tmp being negative), then
      //       the water tends to dry more rapidly (which is in agreement with analytical results),
      //       but is also more 'artefacty' in important cases (tendency for high velocities, etc).
      //       
      //hfactor=1.0;
      a_tmp=0.3; // Highest depth ratio with hfactor=1
      b_tmp=0.1; // Highest depth ratio with hfactor=0
      c_tmp=1.0/(a_tmp-b_tmp); 
      d_tmp= 1.0-(c_tmp*a_tmp);
      hfactor= max(0., min(c_tmp*max(hmin,0.0)/max(hc,1.0e-06)+d_tmp, 
                           min(c_tmp*max(hc,0.)/max(hmax,1.0e-06)+d_tmp, 1.0))
                  );
      //printf("%e, %e, \n", c_tmp, d_tmp);
      //hfactor= max(0., min(5.0*max(hmin,0.0)/max(hmax,1.0e-06)-0.5, 1.0)
      //            );
      //hfactor=1.0;
      // Set hfactor to zero smothly as hmin--> minimum_allowed_height. This
      // avoids some 'chatter' for very shallow flows 
      hfactor=min( 1.2*max(hmin-minimum_allowed_height,0.)/(max(hmin,0.)+1.*minimum_allowed_height), hfactor);

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
      find_qmin_and_qmax(dq0, dq1, dq2, &qmin, &qmax);
    
      beta_tmp = beta_w_dry + (beta_w - beta_w_dry) * hfactor;
      //beta_tmp = beta_w_dry*0. + (beta_w - beta_w_dry*0.) * hfactor;
      //beta_tmp=1.0;
    
      
      // Limit the gradient
      limit_gradient(dqv, qmin, qmax, beta_tmp);

      stage_edge_values[k3+0] = stage_centroid_values[k] + dqv[0];
      stage_edge_values[k3+1] = stage_centroid_values[k] + dqv[1];
      stage_edge_values[k3+2] = stage_centroid_values[k] + dqv[2];


      //-----------------------------------
      // height
      //-----------------------------------
      //hfactor=1.0;
      //hfactor= max(0., min(5.0*max(hmin,0.0)/max(hc,1.0e-06)-0.5, 
      //                     min(5.0*max(hc,0.)/max(hmax,1.0e-06)-0.5, 1.0))
      //            );
       
      // Calculate the difference between vertex 0 of the auxiliary 
      // triangle and the centroid of triangle k
      dq0 = height_centroid_values[k0] - height_centroid_values[k];
      
      // Calculate differentials between the vertices 
      // of the auxiliary triangle (centroids of neighbouring triangles)
      dq1 = height_centroid_values[k1] - height_centroid_values[k0];
      dq2 = height_centroid_values[k2] - height_centroid_values[k0];
     
      inv_area2 = 1.0/area2;
      // Calculate the gradient of height on the auxiliary triangle
      a = dy2*dq1 - dy1*dq2;
      a *= inv_area2;
      b = dx1*dq2 - dx2*dq1;
      b *= inv_area2;
      // Calculate provisional jumps in height from the centroid 
      // of triangle k to its vertices, to be limited
      dqv[0] = a*dxv0 + b*dyv0;
      dqv[1] = a*dxv1 + b*dyv1;
      dqv[2] = a*dxv2 + b*dyv2;
    
      // Now we want to find min and max of the centroid and the 
      // vertices of the auxiliary triangle and compute jumps 
      // from the centroid to the min and max
      find_qmin_and_qmax(dq0, dq1, dq2, &qmin, &qmax);
   
      // Limit the gradient
      //beta_tmp = beta_uh_dry + (beta_uh - beta_uh_dry) * hfactor;
      limit_gradient(dqv, qmin, qmax, beta_tmp);

      //beta_tmp = 0. + (beta_w - 0.) * hfactor;

      height_edge_values[k3+0] = height_centroid_values[k] + dqv[0];
      height_edge_values[k3+1] = height_centroid_values[k] + dqv[1];
      height_edge_values[k3+2] = height_centroid_values[k] + dqv[2];


      // REDEFINE hfactor for momentum terms -- make MORE first order
      // Reduce beta linearly from 1-0 between depth ratio of 0.6-0.4
      //hfactor= max(0., min(5*max(hmin,0.0)/max(hc,1.0e-06)-2.0, 
      //                     min(5*max(hc,0.)/max(hmax,1.0e-06)-2.0, 1.0))
      //            );
      //hfactor= max(0., min(5*max(hmin,0.0)/max(hmax,1.0e-06)-2.0, 
      //                      1.0));
      //hfactor=min( max(hmin,0.)/(max(hmin,0.)+10.*minimum_allowed_height), hfactor);
      //hfactor=1.0;
      
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
      find_qmin_and_qmax(dq0, dq1, dq2, &qmin, &qmax);

      beta_tmp = beta_uh_dry + (beta_uh - beta_uh_dry) * hfactor;

      // Limit the gradient
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
      find_qmin_and_qmax(dq0, dq1, dq2, &qmin, &qmax);
      
      beta_tmp = beta_vh_dry + (beta_vh - beta_vh_dry) * hfactor;   

      // Limit the gradient
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
        report_python_error(AT, "Internal neighbour not found");
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
      
      stage_edge_values[k3] = stage_centroid_values[k] + dqv[0];
      stage_edge_values[k3 + 1] = stage_centroid_values[k] + dqv[1];
      stage_edge_values[k3 + 2] = stage_centroid_values[k] + dqv[2];

      //-----------------------------------
      // height
      //-----------------------------------
      
      // Compute differentials
      dq1 = height_centroid_values[k1] - height_centroid_values[k];
      
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
      
      height_edge_values[k3] = height_centroid_values[k] + dqv[0];
      height_edge_values[k3 + 1] = height_centroid_values[k] + dqv[1];
      height_edge_values[k3 + 2] = height_centroid_values[k] + dqv[2];

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
      
      // Compute height vertex values 
      height_vertex_values[k3] = height_edge_values[k3+1] + height_edge_values[k3+2] -height_edge_values[k3] ; 
      height_vertex_values[k3+1] =  height_edge_values[k3] + height_edge_values[k3+2]-height_edge_values[k3+1]; 
      height_vertex_values[k3+2] =  height_edge_values[k3] + height_edge_values[k3+1]-height_edge_values[k3+2]; 

      // If needed, convert from velocity to momenta
      if(extrapolate_velocity_second_order==1){
          //Convert velocity back to momenta at centroids
          xmom_centroid_values[k] = xmom_centroid_store[k];
          ymom_centroid_values[k] = ymom_centroid_store[k];
      
          // Re-compute momenta at edges
          for (i=0; i<3; i++){
              //if(hfactor>=0.8){
              dk= height_edge_values[k3+i]; 
              //}else{
              //   de[i] = height_centroid_values[k]; 
              //}
              xmom_edge_values[k3+i]=xmom_edge_values[k3+i]*dk;
              ymom_edge_values[k3+i]=ymom_edge_values[k3+i]*dk;
          }
      }
      // Compute momenta at vertices
      xmom_vertex_values[k3] = xmom_edge_values[k3+1] + xmom_edge_values[k3+2] -xmom_edge_values[k3] ; 
      xmom_vertex_values[k3+1] =  xmom_edge_values[k3] + xmom_edge_values[k3+2]-xmom_edge_values[k3+1]; 
      xmom_vertex_values[k3+2] =  xmom_edge_values[k3] + xmom_edge_values[k3+1]-xmom_edge_values[k3+2]; 
      ymom_vertex_values[k3] = ymom_edge_values[k3+1] + ymom_edge_values[k3+2] -ymom_edge_values[k3] ; 
      ymom_vertex_values[k3+1] =  ymom_edge_values[k3] + ymom_edge_values[k3+2]-ymom_edge_values[k3+1]; 
      ymom_vertex_values[k3+2] =  ymom_edge_values[k3] + ymom_edge_values[k3+1]-ymom_edge_values[k3+2]; 

      

      // Compute new bed elevation
      elevation_edge_values[k3]=stage_edge_values[k3]-height_edge_values[k3];
      elevation_edge_values[k3+1]=stage_edge_values[k3+1]-height_edge_values[k3+1];
      elevation_edge_values[k3+2]=stage_edge_values[k3+2]-height_edge_values[k3+2];
      elevation_vertex_values[k3] = elevation_edge_values[k3+1] + elevation_edge_values[k3+2] -elevation_edge_values[k3] ; 
      elevation_vertex_values[k3+1] =  elevation_edge_values[k3] + elevation_edge_values[k3+2]-elevation_edge_values[k3+1]; 
      elevation_vertex_values[k3+2] =  elevation_edge_values[k3] + elevation_edge_values[k3+1]-elevation_edge_values[k3+2]; 

   
  } 

  free(xmom_centroid_store);
  free(ymom_centroid_store);
  //free(min_elevation_edgevalue);
  //free(max_elevation_edgevalue);
  //free(cell_wetness_scale);
  //free(count_wet_neighbours);

  return 0;
}           

//=========================================================================
// Python Glue
//=========================================================================


//========================================================================
// Compute fluxes
//========================================================================

// Modified central flux function

PyObject *compute_fluxes_ext_central(PyObject *self, PyObject *args) {
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
    domain.timestep = compute_fluxes(timestep,
                                     domain.epsilon,
                                     domain.H0,
                                     domain.g,
                                     domain.neighbours,
                                     domain.neighbour_edges,
                                     domain.normals,
                                     domain.edgelengths,
                                     domain.radii,
                                     domain.areas,
                                     tri_full_flag,
                                     Stage.edge_values,
                                     Xmom.edge_values,
                                     Ymom.edge_values,
                                     Bed.edge_values,
                                     Stage.boundary_values,
                                     Xmom.boundary_values,
                                     Ymom.boundary_values,
                                     Stage.explicit_update,
                                     Xmom.explicit_update,
                                     Ymom.explicit_update,
                                     already_computed_flux,
                                     optimise_dry_cells,
                                     stage.centroid_values,
                                     bed.centroid_values,
                                     domain.riverwall_elevation)


    Post conditions:
      domain.explicit_update is reset to computed flux values
      domain.timestep is set to the largest step satisfying all volumes.


  */


  PyArrayObject *boundary_flux_sum, *neighbours, *neighbour_edges,
    *normals, *edgelengths, *radii, *areas,
    *tri_full_flag,
    *stage_edge_values,
    *xmom_edge_values,
    *ymom_edge_values,
    *bed_edge_values,
    *height_edge_values,
    *stage_boundary_values,
    *xmom_boundary_values,
    *ymom_boundary_values,
    *boundary_flux_type,
    *stage_explicit_update,
    *xmom_explicit_update,
    *ymom_explicit_update,
    *already_computed_flux, //Tracks whether the flux across an edge has already been computed
    *max_speed_array, //Keeps track of max speeds for each triangle
    *stage_centroid_values,
    *xmom_centroid_values,
    *ymom_centroid_values,
    *bed_centroid_values,
    *height_centroid_values,
    *bed_vertex_values,
    *edge_flux_type,
    *riverwall_elevation,
    *riverwall_rowIndex,
    *riverwall_hydraulic_properties;
    
  double timestep, epsilon, H0, g;
  int optimise_dry_cells, timestep_fluxcalls, ncol_riverwall_hydraulic_properties;
    
  // Convert Python arguments to C
  if (!PyArg_ParseTuple(args, "ddddOOOOOOOOOOOOOOOOOOOOOOiiOOOOOOOOiO",
            &timestep,
            &epsilon,
            &H0,
            &g,
            &boundary_flux_sum,
            &neighbours,
            &neighbour_edges,
            &normals,
            &edgelengths, &radii, &areas,
            &tri_full_flag,
            &stage_edge_values,
            &xmom_edge_values,
            &ymom_edge_values,
            &bed_edge_values,
            &height_edge_values,
            &stage_boundary_values,
            &xmom_boundary_values,
            &ymom_boundary_values,
            &boundary_flux_type,
            &stage_explicit_update,
            &xmom_explicit_update,
            &ymom_explicit_update,
            &already_computed_flux,
            &max_speed_array,
            &optimise_dry_cells,
            &timestep_fluxcalls,
            &stage_centroid_values,
            &xmom_centroid_values,
            &ymom_centroid_values,
            &bed_centroid_values,
            &height_centroid_values,
            //&bed_vertex_values,
            &edge_flux_type,
            &riverwall_elevation,
            &riverwall_rowIndex,
            &ncol_riverwall_hydraulic_properties,
            &riverwall_hydraulic_properties 
            )) {
    report_python_error(AT, "could not parse input arguments");
    return NULL;
  }

  // check that numpy array objects arrays are C contiguous memory
  CHECK_C_CONTIG(neighbours);
  CHECK_C_CONTIG(neighbour_edges);
  CHECK_C_CONTIG(normals);
  CHECK_C_CONTIG(edgelengths);
  CHECK_C_CONTIG(radii);
  CHECK_C_CONTIG(areas);
  CHECK_C_CONTIG(tri_full_flag);
  CHECK_C_CONTIG(stage_edge_values);
  CHECK_C_CONTIG(xmom_edge_values);
  CHECK_C_CONTIG(ymom_edge_values);
  CHECK_C_CONTIG(bed_edge_values);
  CHECK_C_CONTIG(height_edge_values);
  CHECK_C_CONTIG(stage_boundary_values);
  CHECK_C_CONTIG(xmom_boundary_values);
  CHECK_C_CONTIG(ymom_boundary_values);
  CHECK_C_CONTIG(boundary_flux_type);
  CHECK_C_CONTIG(stage_explicit_update);
  CHECK_C_CONTIG(xmom_explicit_update);
  CHECK_C_CONTIG(ymom_explicit_update);
  CHECK_C_CONTIG(already_computed_flux);
  CHECK_C_CONTIG(max_speed_array);
  CHECK_C_CONTIG(stage_centroid_values);
  CHECK_C_CONTIG(xmom_centroid_values);
  CHECK_C_CONTIG(ymom_centroid_values);
  CHECK_C_CONTIG(bed_centroid_values);
  CHECK_C_CONTIG(height_centroid_values);
  //CHECK_C_CONTIG(bed_vertex_values);
  CHECK_C_CONTIG(edge_flux_type);
  CHECK_C_CONTIG(riverwall_elevation);
  CHECK_C_CONTIG(riverwall_rowIndex);
  CHECK_C_CONTIG(riverwall_hydraulic_properties);

  int number_of_elements = stage_edge_values -> dimensions[0];

  // Call underlying flux computation routine and update 
  // the explicit update arrays 
  timestep = _compute_fluxes_central(number_of_elements,
                     timestep,
                     epsilon,
                     H0,
                     g,
                     (double*) boundary_flux_sum -> data,
                     (long*) neighbours -> data,
                     (long*) neighbour_edges -> data,
                     (double*) normals -> data,
                     (double*) edgelengths -> data, 
                     (double*) radii -> data, 
                     (double*) areas -> data,
                     (long*) tri_full_flag -> data,
                     (double*) stage_edge_values -> data,
                     (double*) xmom_edge_values -> data,
                     (double*) ymom_edge_values -> data,
                     (double*) bed_edge_values -> data,
                     (double*) height_edge_values -> data,
                     (double*) stage_boundary_values -> data,
                     (double*) xmom_boundary_values -> data,
                     (double*) ymom_boundary_values -> data,
                     (long*)   boundary_flux_type -> data,
                     (double*) stage_explicit_update -> data,
                     (double*) xmom_explicit_update -> data,
                     (double*) ymom_explicit_update -> data,
                     (long*) already_computed_flux -> data,
                     (double*) max_speed_array -> data,
                     optimise_dry_cells, 
                     timestep_fluxcalls,
                     (double*) stage_centroid_values -> data, 
                     (double*) xmom_centroid_values -> data, 
                     (double*) ymom_centroid_values -> data, 
                     (double*) bed_centroid_values -> data,
                     (double*) height_centroid_values -> data,
                     //(double*) bed_vertex_values -> data,
                     (long*)   edge_flux_type-> data,
                     (double*) riverwall_elevation-> data,
                     (long*)   riverwall_rowIndex-> data,
                     ncol_riverwall_hydraulic_properties,
                     (double*) riverwall_hydraulic_properties ->data); 
  // Return updated flux timestep
  return Py_BuildValue("d", timestep);
}


PyObject *flux_function_central(PyObject *self, PyObject *args) {
  //
  // Gateway to innermost flux function.
  // This is only used by the unit tests as the c implementation is
  // normally called by compute_fluxes in this module.


  PyArrayObject *normal, *ql, *qr,  *edgeflux;
  double g, epsilon, max_speed, H0, zl, zr;
  double h0, limiting_threshold, pressure_flux, smooth;
  double hle, hre;

  if (!PyArg_ParseTuple(args, "OOOddOddd",
            &normal, &ql, &qr, &zl, &zr, &edgeflux,
            &epsilon, &g, &H0)) {
      report_python_error(AT, "could not parse input arguments");
      return NULL;
  }

  printf("Warning: This interface is broken -- do not use \n");
  
  h0 = H0; 
  hle=1.0; // Fake values to force this to compile
  hre=1.0; // Fake values to force this to compile	      
  limiting_threshold = 10*H0; // Avoid applying limiter below this
                              // threshold for performance reasons.
                              // See ANUGA manual under flux limiting  
 
  pressure_flux = 0.0; // Store the water pressure related component of the flux 
  smooth = 1.0 ; // term to scale diffusion in riemann solver

  _flux_function_central((double*) ql -> data, 
			 (double*) qr -> data, 
			 zl, 
			 zr,
             hle,
             hre,                         
			 ((double*) normal -> data)[0],
			 ((double*) normal -> data)[1],          
			 epsilon, h0, limiting_threshold,
			 g,
			 (double*) edgeflux -> data, 
			 &max_speed,
             &pressure_flux,
			 ((double*) normal -> data)[0],
			 ((double*) normal -> data)[1], 
			 ((double*) normal -> data)[1]
             );
  
  return Py_BuildValue("d", max_speed);  
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


PyObject *extrapolate_second_order_edge_sw(PyObject *self, PyObject *args) {
  /*Compute the edge values based on a linear reconstruction 
    on each triangle
    
    Python call:
    extrapolate_second_order_sw(domain.surrogate_neighbours,
                                domain.number_of_boundaries
                                domain.centroid_coordinates,
                                Stage.centroid_values
                                Xmom.centroid_values
                                Ymom.centroid_values
                                domain.edge_coordinates,
                                Stage.edge_values,
                                Xmom.edge_values,
                                Ymom.edge_values)

    Post conditions:
        The edges of each triangle have values from a 
        limited linear reconstruction
        based on centroid values

  */
  PyArrayObject *surrogate_neighbours,
    *neighbour_edges,
    *number_of_boundaries,
    *centroid_coordinates,
    *stage_centroid_values,
    *xmom_centroid_values,
    *ymom_centroid_values,
    *elevation_centroid_values,
    *height_centroid_values,
    *edge_coordinates,
    *stage_edge_values,
    *xmom_edge_values,
    *ymom_edge_values,
    *elevation_edge_values,
    *height_edge_values,
    *stage_vertex_values,
    *xmom_vertex_values,
    *ymom_vertex_values,
    *elevation_vertex_values,
    *height_vertex_values;
  
  PyObject *domain;

  
  double beta_w, beta_w_dry, beta_uh, beta_uh_dry, beta_vh, beta_vh_dry;    
  double minimum_allowed_height, epsilon;
  int optimise_dry_cells, number_of_elements, extrapolate_velocity_second_order, e, e2;
  
  // Convert Python arguments to C
  if (!PyArg_ParseTuple(args, "OOOOOOOOOOOOOOOOOOOOOii",
            &domain,
            &surrogate_neighbours,
            &neighbour_edges,
            &number_of_boundaries,
            &centroid_coordinates,
            &stage_centroid_values,
            &xmom_centroid_values,
            &ymom_centroid_values,
            &elevation_centroid_values,
            &height_centroid_values,
            &edge_coordinates,
            &stage_edge_values,
            &xmom_edge_values,
            &ymom_edge_values,
            &elevation_edge_values,
            &height_edge_values,
            &stage_vertex_values,
            &xmom_vertex_values,
            &ymom_vertex_values,
            &elevation_vertex_values,
            &height_vertex_values,
            &optimise_dry_cells,
            &extrapolate_velocity_second_order)) {         

    report_python_error(AT, "could not parse input arguments");
    return NULL;
  }

  // check that numpy array objects arrays are C contiguous memory
  CHECK_C_CONTIG(surrogate_neighbours);
  CHECK_C_CONTIG(neighbour_edges);
  CHECK_C_CONTIG(number_of_boundaries);
  CHECK_C_CONTIG(centroid_coordinates);
  CHECK_C_CONTIG(stage_centroid_values);
  CHECK_C_CONTIG(xmom_centroid_values);
  CHECK_C_CONTIG(ymom_centroid_values);
  CHECK_C_CONTIG(elevation_centroid_values);
  CHECK_C_CONTIG(height_centroid_values);
  CHECK_C_CONTIG(edge_coordinates);
  CHECK_C_CONTIG(stage_edge_values);
  CHECK_C_CONTIG(xmom_edge_values);
  CHECK_C_CONTIG(ymom_edge_values);
  CHECK_C_CONTIG(elevation_edge_values);
  CHECK_C_CONTIG(height_edge_values);
  CHECK_C_CONTIG(stage_vertex_values);
  CHECK_C_CONTIG(xmom_vertex_values);
  CHECK_C_CONTIG(ymom_vertex_values);
  CHECK_C_CONTIG(elevation_vertex_values);
  CHECK_C_CONTIG(height_vertex_values);
  
  // Get the safety factor beta_w, set in the config.py file. 
  // This is used in the limiting process
  

  beta_w                 = get_python_double(domain,"beta_w");
  beta_w_dry             = get_python_double(domain,"beta_w_dry");
  beta_uh                = get_python_double(domain,"beta_uh");
  beta_uh_dry            = get_python_double(domain,"beta_uh_dry");
  beta_vh                = get_python_double(domain,"beta_vh");
  beta_vh_dry            = get_python_double(domain,"beta_vh_dry");  

  minimum_allowed_height = get_python_double(domain,"minimum_allowed_height");
  epsilon                = get_python_double(domain,"epsilon");

  number_of_elements = stage_centroid_values -> dimensions[0];  

  //printf("In C before Extrapolate");
  //e=1;
  // Call underlying computational routine
  e = _extrapolate_second_order_edge_sw(number_of_elements,
                   epsilon,
                   minimum_allowed_height,
                   beta_w,
                   beta_w_dry,
                   beta_uh,
                   beta_uh_dry,
                   beta_vh,
                   beta_vh_dry,
                   (long*) surrogate_neighbours -> data,
                   (long*) neighbour_edges -> data,
                   (long*) number_of_boundaries -> data,
                   (double*) centroid_coordinates -> data,
                   (double*) stage_centroid_values -> data,
                   (double*) xmom_centroid_values -> data,
                   (double*) ymom_centroid_values -> data,
                   (double*) elevation_centroid_values -> data,
                   (double*) height_centroid_values -> data,
                   (double*) edge_coordinates -> data,
                   (double*) stage_edge_values -> data,
                   (double*) xmom_edge_values -> data,
                   (double*) ymom_edge_values -> data,
                   (double*) elevation_edge_values -> data,
                   (double*) height_edge_values -> data,
                   (double*) stage_vertex_values -> data,
                   (double*) xmom_vertex_values -> data,
                   (double*) ymom_vertex_values -> data,
                   (double*) elevation_vertex_values -> data,
                   (double*) height_vertex_values -> data,
                   optimise_dry_cells, 
                   extrapolate_velocity_second_order);

  if (e == -1) {
    // Use error string set inside computational routine
    return NULL;
  }                   
  
  
  return Py_BuildValue("");
  
}// extrapolate_second-order_edge_sw
//========================================================================
// Protect -- to prevent the water level from falling below the minimum 
// bed_edge_value
//========================================================================

PyObject *protect(PyObject *self, PyObject *args) {
  //
  //    protect(minimum_allowed_height, maximum_allowed_speed, wc, zc, xmomc, ymomc)


  PyArrayObject
  *wc,            //Stage at centroids
  *wv,            //Stage at vertices
  *zc,            //Elevation at centroids
  *zv,            //Elevation at vertices
  *xmomc,         //Momentums at centroids
  *ymomc,
  *areas,         // Triangle areas
  *xc,
  *yc;

  int N;
  double mass_error;
  double minimum_allowed_height, maximum_allowed_speed, epsilon;

  // Convert Python arguments to C
  if (!PyArg_ParseTuple(args, "dddOOOOOOOOO",
            &minimum_allowed_height,
            &maximum_allowed_speed,
            &epsilon,
            &wc, &wv, &zc, &zv, &xmomc, &ymomc, &areas, &xc, &yc)) {
    report_python_error(AT, "could not parse input arguments");
    return NULL;
  }  

  N = wc -> dimensions[0];

  mass_error = _protect(N,
       minimum_allowed_height,
       maximum_allowed_speed,
       epsilon,
       (double*) wc -> data,
       (double*) wv -> data,
       (double*) zc -> data,
       (double*) zv -> data,
       (double*) xmomc -> data,
       (double*) ymomc -> data,
       (double*) areas -> data,
       (double*) xc -> data,
       (double*) yc -> data );

  return Py_BuildValue("d", mass_error);
}

//========================================================================
// Method table for python module
//========================================================================

static struct PyMethodDef MethodTable[] = {
  /* The cast of the function is necessary since PyCFunction values
   * only take two PyObject* parameters, and rotate() takes
   * three.
   */
  //{"rotate", (PyCFunction)rotate, METH_VARARGS | METH_KEYWORDS, "Print out"},
  {"compute_fluxes_ext_central", compute_fluxes_ext_central, METH_VARARGS, "Print out"},
  {"gravity_c",        gravity,            METH_VARARGS, "Print out"},
  {"flux_function_central", flux_function_central, METH_VARARGS, "Print out"},  
  {"extrapolate_second_order_edge_sw", extrapolate_second_order_edge_sw, METH_VARARGS, "Print out"},
  {"protect",          protect, METH_VARARGS | METH_KEYWORDS, "Print out"},
  {NULL, NULL, 0, NULL}
};

// Module initialisation
void initswDE1_domain_ext(void){
  Py_InitModule("swDE1_domain_ext", MethodTable);

  import_array(); // Necessary for handling of NumPY structures
}

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
#include <string.h>

// Shared code snippets
#include "util_ext.h"
#include "sw_domain.h"

const double pi = 3.14159265358979;

// Trick to compute n modulo d (n%d in python) when d is a power of 2
unsigned int Mod_of_power_2(unsigned int n, unsigned int d)
{
  return ( n & (d-1) );
}


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
////static double _compute_speed(double *uh,
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
                           double hc_n)
{

  /*Compute fluxes between volumes for the shallow water wave equation
    cast in terms of the 'stage', w = h+z using

	HLL scheme of Fraccarollo and Toro Experimental and numerical assessment of the shallow
water model for two-dimensional dam-break type. Journal of Computational Physics,
33(6):843–864, 1995.

    FIXME: Several variables in this interface are no longer used, clean up
  */

  int i;

  double uh_left, vh_left, u_left;
  double uh_right, vh_right, u_right;
  double s_min, s_max, soundspeed_left, soundspeed_right;
  double u_m, h_m, soundspeed_m;
  double denom, inverse_denominator;
  double tmp;
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
    tmp = 1.0 / hle;
    u_left = uh_left * tmp ;
    uh_left = h_left * u_left;
    vh_left = h_left* tmp * vh_left;
  }else{
    u_left = 0.;
    uh_left = 0.;
    vh_left = 0.;
  }

  //u_left = _compute_speed(&uh_left, &hle,
  //            epsilon, h0, limiting_threshold);

  //w_right = q_right_rotated[0];
  uh_right = q_right_rotated[1];
  vh_right = q_right_rotated[2];
  if(hre>0.0){
    tmp = 1.0 / hre;
    u_right = uh_right * tmp; //max(h_right, 1.0e-06);
    uh_right = h_right * u_right;
    vh_right = h_right * tmp * vh_right;
  }else{
    u_right = 0.;
    uh_right = 0.;
    vh_right = 0.;
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
	  //s_m   = s_min;
  }
  else if (h_right < 1.0e-10){
	  s_min = u_left - soundspeed_left;
	  s_max = u_left + 2.0*soundspeed_left;
	  //s_m = s_max;
  }
  else {
	  s_max = fmax(u_right + soundspeed_right, u_m + soundspeed_right);
	  s_min = fmin(u_left - soundspeed_left, u_m - soundspeed_m);
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
    *max_speed = fmax(s_max, -s_min);

    inverse_denominator = 1.0/fmax(denom,1.0e-100);
    for (i = 0; i < 3; i++)
    {
      edgeflux[i] = s_max*flux_left[i] - s_min*flux_right[i];

      // Standard smoothing term
      // edgeflux[i] += 1.0*(s_max*s_min)*(q_right_rotated[i] - q_left_rotated[i]);
      // Smoothing by stage alone can cause high velocities / slow draining for nearly dry cells
      if(i==0) edgeflux[i] += (s_max*s_min)*(fmax(q_right_rotated[i],ze) - fmax(q_left_rotated[i],ze));
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
                           double hc_n,
                           long low_froude)
{

  /*Compute fluxes between volumes for the shallow water wave equation
    cast in terms of the 'stage', w = h+z using
    the 'central scheme' as described in

    Kurganov, Noelle, Petrova. 'Semidiscrete Central-Upwind Schemes For
    Hyperbolic Conservation Laws and Hamilton-Jacobi Equations'.
    Siam J. Sci. Comput. Vol. 23, No. 3, pp. 707-740.

    The implemented formula is given in equation (3.15) on page 714

    FIXME: Several variables in this interface are no longer used, clean up

    low_froude is either 1 or 2 - see comment inline
  */

  int i;

  double uh_left, vh_left, u_left;
  double uh_right, vh_right, u_right;
  double s_min, s_max, soundspeed_left, soundspeed_right;
  double denom, inverse_denominator;
  double tmp, local_fr, v_right, v_left;
  // Workspace (allocate once, use many)
  static double q_left_rotated[3], q_right_rotated[3], flux_right[3], flux_left[3];

  printf("Entering _flux_function_central with h_left = %f, h_right = %f\n", h_left, h_right);

  if(h_left==0. && h_right==0.){
    // Quick exit
    memset(edgeflux, 0, 3*sizeof(double));
    *max_speed = 0.0;
    *pressure_flux = 0.;
    return 0;
  }
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
    tmp = 1.0/hle;
    u_left = uh_left * tmp ; //max(h_left, 1.0e-06);
    uh_left = h_left * u_left;
    v_left = vh_left * tmp;  // Only used to define local_fr
    vh_left = h_left * tmp * vh_left;
  }else{
    u_left = 0.;
    uh_left = 0.;
    vh_left = 0.;
    v_left = 0.;
  }

  //u_left = _compute_speed(&uh_left, &hle,
  //            epsilon, h0, limiting_threshold);

  //w_right = q_right_rotated[0];
  uh_right = q_right_rotated[1];
  vh_right = q_right_rotated[2];
  if(hre>0.0){
    tmp = 1.0 / hre;
    u_right = uh_right * tmp;//max(h_right, 1.0e-06);
    uh_right=h_right*u_right;
    v_right = vh_right * tmp; // Only used to define local_fr
    vh_right=h_right * tmp * vh_right;
  }else{
    u_right=0.;
    uh_right=0.;
    vh_right=0.;
    v_right = 0.;
  }
  //u_right = _compute_speed(&uh_right, &hre,
  //              epsilon, h0, limiting_threshold);


  // Maximal and minimal wave speeds
  soundspeed_left  = sqrt(g*h_left);
  soundspeed_right = sqrt(g*h_right);
  //soundspeed_left  = sqrt(g*hle);
  //soundspeed_right = sqrt(g*hre);

  // Something that scales like the Froude number
  // We will use this to scale the diffusive component of the UH/VH fluxes.

  //local_fr = sqrt(
  //    max(0.001, min(1.0,
  //        (u_right*u_right + u_left*u_left + v_right*v_right + v_left*v_left)/
  //        (soundspeed_left*soundspeed_left + soundspeed_right*soundspeed_right + 1.0e-10))));
  if (low_froude == 1)
  {
    local_fr = sqrt(
      fmax(0.001, fmin(1.0,
          (u_right*u_right + u_left*u_left + v_right*v_right + v_left*v_left)/
          (soundspeed_left*soundspeed_left + soundspeed_right*soundspeed_right + 1.0e-10))));
  }
  else if (low_froude == 2)
  {
    local_fr = sqrt((u_right*u_right + u_left*u_left + v_right*v_right + v_left*v_left)/
          (soundspeed_left*soundspeed_left + soundspeed_right*soundspeed_right + 1.0e-10));
    local_fr = sqrt(fmin(1.0, 0.01 + fmax(local_fr-0.01,0.0)));
  }
  else
  {
  local_fr = 1.0;
  }
  //printf("local_fr %e \n:", local_fr);

  s_max = fmax(u_left + soundspeed_left, u_right + soundspeed_right);
  if (s_max < 0.0)
  {
    s_max = 0.0;
  }

  //if( hc < 1.0e-03){
  //  s_max = 0.0;
  //}


  s_min = fmin(u_left - soundspeed_left, u_right - soundspeed_right);
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
    *max_speed = fmax(s_max, -s_min);

    inverse_denominator = 1.0/fmax(denom,1.0e-100);
    for (i = 0; i < 3; i++)
    {
      edgeflux[i] = s_max*flux_left[i] - s_min*flux_right[i];

      // Standard smoothing term
      // edgeflux[i] += 1.0*(s_max*s_min)*(q_right_rotated[i] - q_left_rotated[i]);
      // Smoothing by stage alone can cause high velocities / slow draining for nearly dry cells
      if(i==0) edgeflux[i] += (s_max*s_min)*(fmax(q_right_rotated[i],ze) - fmax(q_left_rotated[i],ze));
      //if(i==0) edgeflux[i] += (s_max*s_min)*(h_right - h_left);
      if(i==1) edgeflux[i] += local_fr*(s_max*s_min)*(uh_right - uh_left);
      if(i==2) edgeflux[i] += local_fr*(s_max*s_min)*(vh_right - vh_left);

      edgeflux[i] *= inverse_denominator;
    }
    // Separate pressure flux, so we can apply different wet-dry hacks to it
    *pressure_flux = 0.5*g*( s_max*h_left*h_left -s_min*h_right*h_right)*inverse_denominator;

    // Rotate back
    _rotate(edgeflux, n1, -n2);
  }

  printf("Pressure flux = %f, s_max = %f, s_min = %f, h_left = %f, h_right = %f\n", *pressure_flux, s_max, s_min, h_left, h_right);


  return 0;
}

////////////////////////////////////////////////////////////////

int _compute_flux_update_frequency(struct domain *D, double timestep){
    // Compute the 'flux_update_frequency' for each edge.
    //
    // This determines how regularly we need
    // to update the flux computation (not every timestep)
    //
    // Allowed values are 1,2,4,8,... max_flux_update_frequency.
    //
    // For example, an edge with flux_update_frequency = 4 would
    // only have the flux updated every 4 timesteps
    //
    //
    // Local variables
    int k, i, k3, ki, m, n, nm, ii, ii2;
    long fuf;
    double notSoFast=1.0;
    static int cyclic_number_of_steps=-1;

    // QUICK EXIT
    if(D->max_flux_update_frequency==1){
        return 0;
    }

    // Count the steps
    cyclic_number_of_steps++;
    if(cyclic_number_of_steps==D->max_flux_update_frequency){
        // The flux was just updated in every cell
        cyclic_number_of_steps=0;
    }


    // PART 1: ONLY OCCURS FOLLOWING FLUX UPDATE

    for ( k = 0; k < D->number_of_elements; k++){
        for ( i = 0; i < 3; i++){
            ki = k*3 + i;
            if((Mod_of_power_2(cyclic_number_of_steps, D->flux_update_frequency[ki])==0)){
                // The flux was just updated, along with the edge_timestep
                // So we better recompute the flux_update_frequency
                n=D->neighbours[ki];
                if(n>=0){
                    m = D->neighbour_edges[ki];
                    nm = n * 3 + m; // Linear index (triangle n, edge m)
                }

                // Check if we have already done this edge
                // (Multiply already_computed_flux by -1 on the first update,
                // and again on the 2nd)
                if(D->already_computed_flux[ki] > 0 ){
                    // We have not fixed this flux value yet
                    D->already_computed_flux[ki] *=-1;
                    if(n>=0){
                        D->already_computed_flux[nm] *=-1;
                    }
                }else{
                    // We have fixed this flux value already
                    D->already_computed_flux[ki] *=-1;
                    if(n>=0){
                        D->already_computed_flux[nm] *=-1;
                    }
                    continue;
                }

                // Basically int( edge_ki_timestep/timestep ) with upper limit + tweaks
                // notSoFast is ideally = 1.0, but in practice values < 1.0 can enhance stability
                // NOTE: edge_timestep[ki]/timestep can be very large [so int overflows].
                //       Do not pull the (int) inside the min term
                fuf = (int)fmin((D->edge_timestep[ki]/timestep)*notSoFast,D->max_flux_update_frequency*1.);
                // Account for neighbour
                if(n>=0){
                    fuf = fmin( (int)fmin(D->edge_timestep[nm]/timestep*notSoFast, D->max_flux_update_frequency*1.), fuf);
                }

                // Deal with notSoFast<1.0
                if(fuf<1){
                    fuf=1;
                }
                // Deal with large fuf
                if(fuf> D->max_flux_update_frequency){
                    fuf = D->max_flux_update_frequency;
                }
                //// Deal with intermediate cases
                ii=2;
                while(ii< D->max_flux_update_frequency){
                    // Set it to 1,2,4, 8, ...
                    ii2=2*ii;
                    if((fuf>ii) && (fuf<ii2)){
                        fuf = ii;
                        continue;
                    }
                    ii=ii2;
                }

                // Set the values
                D->flux_update_frequency[ki]=fuf;
                if(n>=0){
                    D->flux_update_frequency[nm]= fuf;
                }

            }
        }
    }

    //// PART 2 -- occcurs every timestep

    // At this point, both edges have the same flux_update_frequency.
    // Next, ensure that flux_update_frequency varies within a constant over each triangle
    // Experiences suggests this is numerically important
    // (But, it can result in the same edge having different flux_update_freq)
    for( k=0; k< D->number_of_elements; k++){
        k3=3*k;
        ii = 1*fmin(D->flux_update_frequency[k3],
                 fmin(D->flux_update_frequency[k3+1],
                     D->flux_update_frequency[k3+2]));

        D->flux_update_frequency[k3]=fmin(ii, D->flux_update_frequency[k3]);
        D->flux_update_frequency[k3+1]=fmin(ii, D->flux_update_frequency[k3+1]);
        D->flux_update_frequency[k3+2]=fmin(ii,D->flux_update_frequency[k3+2]);

    }

    // Now enforce the same flux_update_frequency on each edge
    // (Could have been broken above when we limited the variation on each triangle)
    // This seems to have nice behaviour. Notice how an edge
    // with a large flux_update_frequency, near an edge with a small flux_update_frequency,
    // will have its flux_update_frequency updated after a few timesteps (i.e. before max_flux_update_frequency timesteps)
    // OTOH, could this cause oscillations in flux_update_frequency?
    for( k = 0; k < D->number_of_elements; k++){
        D->update_extrapolation[k]=0;
        for( i = 0; i< 3; i++){
            ki=3*k+i;
            // Account for neighbour
            n=D->neighbours[ki];
            if(n>=0){
                m = D->neighbour_edges[ki];
                nm = n * 3 + m; // Linear index (triangle n, edge m)
                D->flux_update_frequency[ki]=fmin(D->flux_update_frequency[ki], D->flux_update_frequency[nm]);
            }
            // Do we need to update the extrapolation?
            // (We do if the next flux computation will actually compute a flux!)
            if(Mod_of_power_2((cyclic_number_of_steps+1),D->flux_update_frequency[ki])==0){
                D->update_next_flux[ki]=1;
                D->update_extrapolation[k]=1;
            }else{
                D->update_next_flux[ki]=0;
            }
        }
    }

    // Check whether the timestep can be increased in the next compute_fluxes call
    if(cyclic_number_of_steps+1==D->max_flux_update_frequency){
        // All fluxes will be updated on the next timestep
        // We also allow the timestep to increase then
        D->allow_timestep_increase[0]=1;
    }else{
        D->allow_timestep_increase[0]=0;
    }

    return 0;
}


double adjust_edgeflux_with_weir(double *edgeflux,
                                 double h_left, double h_right,
                                 double g, double weir_height,
                                 double Qfactor,
                                 double s1, double s2,
                                 double h1, double h2,
                                 double *max_speed_local
                                ){
    // Adjust the edgeflux to agree with a weir relation [including
    // subergence], but smoothly vary to shallow water solution when
    // the flow over the weir is much deeper than the weir, or the
    // upstream/downstream water elevations are too similar
    double rw, rw2; // 'Raw' weir fluxes
    double rwRat, hdRat,hdWrRat, scaleFlux, minhd, maxhd;
    double w1,w2; // Weights for averaging
    double newFlux;
    double twothirds = (2.0/3.0);
    // Following constants control the 'blending' with the shallow water solution
    // They are now user-defined
    //double s1=0.9; // At this submergence ratio, begin blending with shallow water solution
    //double s2=0.95; // At this submergence ratio, completely use shallow water solution
    //double h1=1.0; // At this (tailwater height above weir) / (weir height) ratio, begin blending with shallow water solution
    //double h2=1.5; // At this (tailwater height above weir) / (weir height) ratio, completely use the shallow water solution

    minhd = fmin(h_left, h_right);
    maxhd = fmax(h_left, h_right);
    // 'Raw' weir discharge = Qfactor*2/3*H*(2/3*g*H)**0.5
    rw = Qfactor * twothirds * maxhd * sqrt(twothirds * g * maxhd);
    // Factor for villemonte correction
    rw2 = Qfactor * twothirds * minhd * sqrt(twothirds * g * minhd);
    // Useful ratios
    rwRat = rw2 / fmax(rw, 1.0e-100);
    hdRat = minhd / fmax(maxhd, 1.0e-100);

    // (tailwater height above weir)/weir_height ratio
    hdWrRat = minhd / fmax(weir_height, 1.0e-100);

    // Villemonte (1947) corrected weir flow with submergence
    // Q = Q1*(1-Q2/Q1)**0.385
    rw = rw*pow(1.0 - rwRat, 0.385);

    if(h_right > h_left){
        rw *= -1.0;
    }

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
        w1 = fmin( fmax(hdRat-s1, 0.) / (s2-s1), 1.0);

        // Adjust again when the head is too deep relative to the weir height
        w2 = fmin( fmax(hdWrRat-h1,0.) / (h2-h1), 1.0);

        newFlux = (rw*(1.0-w1)+w1*edgeflux[0])*(1.0-w2) + w2*edgeflux[0];

        if(fabs(edgeflux[0]) > 1.0e-100){
            scaleFlux = newFlux/edgeflux[0];
        }else{
            scaleFlux = 0.;
        }

        scaleFlux = fmax(scaleFlux, 0.);

        edgeflux[0] = newFlux;

        // FIXME: Do this in a cleaner way
        // IDEA: Compute momentum flux implied by weir relations, and use
        //       those in a weighted average (rather than the rescaling trick here)
        // If we allow the scaling to momentum to be unbounded,
        // velocity spikes can arise for very-shallow-flooded walls
        edgeflux[1] *= fmin(scaleFlux, 10.);
        edgeflux[2] *= fmin(scaleFlux, 10.);
    }

    // Adjust the max speed
    if (fabs(edgeflux[0]) > 0.){
        *max_speed_local = sqrt(g*(maxhd+weir_height)) + fabs(edgeflux[0]/(maxhd + 1.0e-12));
    }
    //*max_speed_local += fabs(edgeflux[0])/(maxhd+1.0e-100);
    //*max_speed_local *= fmax(scaleFlux, 1.0);

    return 0;
}

// Computational function for flux computation
double _compute_fluxes_central(struct domain *D, double timestep){

    // Local variables
    double max_speed_local, length, inv_area, zl, zr;
    double h_left, h_right, z_half ;  // For andusse scheme
    // FIXME: limiting_threshold is not used for DE1
    double limiting_threshold = 10*D->H0;
    long low_froude = D->low_froude;
    //
    int k, i, m, n, ii;
    int ki, nm = 0, ki2,ki3, nm3; // Index shorthands
    // Workspace (making them static actually made function slightly slower (Ole))
    double ql[3], qr[3], edgeflux[3]; // Work array for summing up fluxes
    double bedslope_work;
    static double local_timestep;
    long RiverWall_count, substep_count;
    double hle, hre, zc, zc_n, Qfactor, s1, s2, h1, h2;
    double pressure_flux, hc, hc_n, tmp;
    double h_left_tmp, h_right_tmp;
    static long call = 0; // Static local variable flagging already computed flux
    static long timestep_fluxcalls=1;
    static long base_call = 1;
    double speed_max_last, weir_height;

    call++; // Flag 'id' of flux calculation for this timestep

    if (D->timestep_fluxcalls != timestep_fluxcalls) {
    	timestep_fluxcalls = D->timestep_fluxcalls;
    	base_call = call;
    }

    // Set explicit_update to zero for all conserved_quantities.
    // This assumes compute_fluxes called before forcing terms
    memset((char*) D->stage_explicit_update, 0, D->number_of_elements * sizeof (double));
    memset((char*) D->xmom_explicit_update, 0, D->number_of_elements * sizeof (double));
    memset((char*) D->ymom_explicit_update, 0, D->number_of_elements * sizeof (double));


    // Counter for riverwall edges
    RiverWall_count=0;
    // Which substep of the timestepping method are we on?
    substep_count=(call-base_call)%D->timestep_fluxcalls;

    //printf("call = %d substep_count = %d base_call = %d \n",call,substep_count, base_call);

    // Fluxes are not updated every timestep,
    // but all fluxes ARE updated when the following condition holds
    if(D->allow_timestep_increase[0]==1){
        // We can only increase the timestep if all fluxes are allowed to be updated
        // If this is not done the timestep can't increase (since local_timestep is static)
        local_timestep=1.0e+100;
    }

    // For all triangles
    for (k = 0; k < D->number_of_elements; k++) {
        speed_max_last = 0.0;

        // Loop through neighbours and compute edge flux for each
        for (i = 0; i < 3; i++) {
            ki = k * 3 + i; // Linear index to edge i of triangle k
            ki2 = 2 * ki; //k*6 + i*2
            ki3 = 3*ki;

            if ((D->already_computed_flux[ki] == call) || (D->update_next_flux[ki]!=1)) {
                // We've already computed the flux across this edge
                // Check if it is a riverwall
                if(D->edge_flux_type[ki] == 1){
                    // Update counter of riverwall edges == index of
                    // riverwall_elevation + riverwall_rowIndex
                    RiverWall_count += 1;
                }
                continue;
            }

            // Get left hand side values from triangle k, edge i
            ql[0] = D->stage_edge_values[ki];
            ql[1] = D->xmom_edge_values[ki];
            ql[2] = D->ymom_edge_values[ki];
            zl = D->bed_edge_values[ki];
            hc = D->height_centroid_values[k];
            zc = D->bed_centroid_values[k];
            hle= D->height_edge_values[ki];

            // Get right hand side values either from neighbouring triangle
            // or from boundary array (Quantities at neighbour on nearest face).
            n = D->neighbours[ki];
            hc_n = hc;
            zc_n = D->bed_centroid_values[k];
            if (n < 0) {
                // Neighbour is a boundary condition
                m = -n - 1; // Convert negative flag to boundary index

                qr[0] = D->stage_boundary_values[m];
                qr[1] = D->xmom_boundary_values[m];
                qr[2] = D->ymom_boundary_values[m];
                zr = zl; // Extend bed elevation to boundary
                hre= fmax(qr[0]-zr,0.);//hle;
            } else {
                // Neighbour is a real triangle
                hc_n = D->height_centroid_values[n];
                zc_n = D->bed_centroid_values[n];
                m = D->neighbour_edges[ki];
                nm = n * 3 + m; // Linear index (triangle n, edge m)
                nm3 = nm*3;

                qr[0] = D->stage_edge_values[nm];
                qr[1] = D->xmom_edge_values[nm];
                qr[2] = D->ymom_edge_values[nm];
                zr = D->bed_edge_values[nm];
                hre = D->height_edge_values[nm];
            }

            // Audusse magic
            z_half = fmax(zl, zr);

            //// Account for riverwalls
            if(D->edge_flux_type[ki] == 1){
                if( n>=0 && D->edge_flux_type[nm] != 1){
                    printf("Riverwall Error\n");
                }
                // Update counter of riverwall edges == index of
                // riverwall_elevation + riverwall_rowIndex
                RiverWall_count += 1;

                // Set central bed to riverwall elevation
                z_half = fmax(D->riverwall_elevation[RiverWall_count-1], z_half) ;

            }

            // Define h left/right for Audusse flux method
            h_left = fmax(hle+zl-z_half,0.);
            h_right = fmax(hre+zr-z_half,0.);

            // Edge flux computation (triangle k, edge i)
            _flux_function_central(ql, qr,
            //_flux_function_toro(ql, qr,
                h_left, h_right,
                hle, hre,
                D->normals[ki2],D->normals[ki2 + 1],
                D->epsilon, z_half, limiting_threshold, D->g,
                edgeflux, &max_speed_local, &pressure_flux, hc, hc_n, low_froude);

            // Force weir discharge to match weir theory
            // FIXME: Switched off at the moment
            if(D->edge_flux_type[ki]==1){
                weir_height = fmax(D->riverwall_elevation[RiverWall_count-1] - fmin(zl, zr), 0.); // Reference weir height

                // If the weir is not higher than both neighbouring cells, then
                // do not try to match the weir equation. If we do, it seems we
                // can get mass conservation issues (caused by large weir
                // fluxes in such situations)
                if(D->riverwall_elevation[RiverWall_count-1] > fmax(zc, zc_n)){
                    ////////////////////////////////////////////////////////////////////////////////////
                    // Use first-order h's for weir -- as the 'upstream/downstream' heads are
                    //  measured away from the weir itself
                    h_left_tmp = fmax(D->stage_centroid_values[k] - z_half, 0.);
                    if(n >= 0){
                        h_right_tmp = fmax(D->stage_centroid_values[n] - z_half, 0.);
                    }else{
                        h_right_tmp = fmax(hc_n + zr - z_half, 0.);
                    }

                    if( (h_left_tmp > 0.) || (h_right_tmp > 0.)){

                        //////////////////////////////////////////////////////////////////////////////////
                        // Get Qfactor index - multiply the idealised weir discharge by this constant factor
                        ii = D->riverwall_rowIndex[RiverWall_count-1] * D->ncol_riverwall_hydraulic_properties;
                        Qfactor = D->riverwall_hydraulic_properties[ii];

                        // Get s1, submergence ratio at which we start blending with the shallow water solution
                        ii+=1;
                        s1 = D->riverwall_hydraulic_properties[ii];

                        // Get s2, submergence ratio at which we entirely use the shallow water solution
                        ii+=1;
                        s2 = D->riverwall_hydraulic_properties[ii];

                        // Get h1, tailwater head / weir height at which we start blending with the shallow water solution
                        ii+=1;
                        h1 = D->riverwall_hydraulic_properties[ii];

                        // Get h2, tailwater head / weir height at which we entirely use the shallow water solution
                        ii+=1;
                        h2 = D->riverwall_hydraulic_properties[ii];

                        // Weir flux adjustment
                        // FIXME
                        adjust_edgeflux_with_weir(edgeflux, h_left_tmp, h_right_tmp, D->g,
                                                  weir_height, Qfactor,
                                                  s1, s2, h1, h2, &max_speed_local);
                    }
                }
            }

            // Multiply edgeflux by edgelength
            length = D->edgelengths[ki];
            edgeflux[0] *= length;
            edgeflux[1] *= length;
            edgeflux[2] *= length;

            //// Don't allow an outward advective flux if the cell centroid
            ////   stage is < the edge value. Is this important (??). Seems not
            ////   to be with DE algorithms
            //if((hc<H0) && edgeflux[0] > 0.){
            //    edgeflux[0] = 0.;
            //    edgeflux[1] = 0.;
            //    edgeflux[2] = 0.;
            //    //max_speed_local=0.;
            //    //pressure_flux=0.;
            //}
            ////
            //if((hc_n<H0) && edgeflux[0] < 0.){
            //    edgeflux[0] = 0.;
            //    edgeflux[1] = 0.;
            //    edgeflux[2] = 0.;
            //    //max_speed_local=0.;
            //    //pressure_flux=0.;
            //}

            D->edge_flux_work[ki3 + 0 ] = -edgeflux[0];
            D->edge_flux_work[ki3 + 1 ] = -edgeflux[1];
            D->edge_flux_work[ki3 + 2 ] = -edgeflux[2];

            // bedslope_work contains all gravity related terms
            bedslope_work = length*(- D->g *0.5*(h_left*h_left - hle*hle -(hle+hc)*(zl-zc))+pressure_flux);

            D->pressuregrad_work[ki] = bedslope_work;

            D->already_computed_flux[ki] = call; // #k Done

            // Update neighbour n with same flux but reversed sign
            if (n >= 0) {

                D->edge_flux_work[nm3 + 0 ] = edgeflux[0];
                D->edge_flux_work[nm3 + 1 ] = edgeflux[1];
                D->edge_flux_work[nm3 + 2 ] = edgeflux[2];
                bedslope_work = length*(-D->g * 0.5 *( h_right*h_right - hre*hre- (hre+hc_n)*(zr-zc_n)) + pressure_flux);
                D->pressuregrad_work[nm] = bedslope_work;

                D->already_computed_flux[nm] = call; // #n Done
            }

            // Update timestep based on edge i and possibly neighbour n
            // NOTE: We should only change the timestep on the 'first substep'
            //  of the timestepping method [substep_count==0]
            if(substep_count==0){

                // Compute the 'edge-timesteps' (useful for setting flux_update_frequency)
                tmp = 1.0 / fmax(max_speed_local, D->epsilon);
                D->edge_timestep[ki] = D->radii[k] * tmp ;
                if (n >= 0) {
                    D->edge_timestep[nm] = D->radii[n] * tmp;
                }

                // Update the timestep
                if ((D->tri_full_flag[k] == 1)) {

                    speed_max_last = fmax(speed_max_last, max_speed_local);

                    if (max_speed_local > D->epsilon) {
                        // Apply CFL condition for triangles joining this edge (triangle k and triangle n)

                        // CFL for triangle k
                        local_timestep = fmin(local_timestep, D->edge_timestep[ki]);

                        if (n >= 0) {
                            // Apply CFL condition for neigbour n (which is on the ith edge of triangle k)
                            local_timestep = fmin(local_timestep, D->edge_timestep[nm]);
                        }
                    }
                }
            }

        } // End edge i (and neighbour n)
        // Keep track of maximal speeds
        if(substep_count==0) D->max_speed[k] = speed_max_last; //max_speed;


    } // End triangle k

    //// Limit edgefluxes, for mass conservation near wet/dry cells
    //// This doesn't seem to be needed anymore
    //for(k=0; k< number_of_elements; k++){
    //    //continue;
    //    hc = height_centroid_values[k];
    //    // Loop over every edge
    //    for(i = 0; i<3; i++){
    //        if(i==0){
    //            // Add up the outgoing flux through the cell -- only do this once (i==0)
    //            outgoing_mass_edges=0.0;
    //            for(useint=0; useint<3; useint++){
    //                if(edge_flux_work[3*(3*k+useint)]< 0.){
    //                    //outgoing_mass_edges+=1.0;
    //                    outgoing_mass_edges+=(edge_flux_work[3*(3*k+useint)]);
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
    //        if((edge_flux_work[ki3]< 0.0) && (-outgoing_mass_edges> vol)){
    //
    //            // This bound could be improved (e.g. we could actually sum the
    //            // + and - fluxes and check if they are too large).  However,
    //            // the advantage of this method is that we don't have to worry
    //            // about subsequent changes to the + edgeflux caused by
    //            // constraints associated with neighbouring triangles.
    //            tmp = vol/(-(outgoing_mass_edges)) ;
    //            if(tmp< 1.0){
    //                edge_flux_work[ki3+0]*=tmp;
    //                edge_flux_work[ki3+1]*=tmp;
    //                edge_flux_work[ki3+2]*=tmp;

    //                // Compute neighbour edge index
    //                n = neighbours[ki];
    //                if(n>=0){
    //                    nm = 3*n + neighbour_edges[ki];
    //                    nm3 = nm*3;
    //                    edge_flux_work[nm3+0]*=tmp;
    //                    edge_flux_work[nm3+1]*=tmp;
    //                    edge_flux_work[nm3+2]*=tmp;
    //                }
    //            }
    //        }
    //    }
    // }

    // Now add up stage, xmom, ymom explicit updates
    for(k=0; k < D->number_of_elements; k++){
        hc = fmax(D->stage_centroid_values[k] - D->bed_centroid_values[k],0.);

        for(i=0;i<3;i++){
            // FIXME: Make use of neighbours to efficiently set things
            ki=3*k+i;
            ki2=ki*2;
            ki3 = ki*3;
            n=D->neighbours[ki];

            D->stage_explicit_update[k] += D->edge_flux_work[ki3+0];
            D->xmom_explicit_update[k] += D->edge_flux_work[ki3+1];
            D->ymom_explicit_update[k] += D->edge_flux_work[ki3+2];

            // If this cell is not a ghost, and the neighbour is a boundary
            // condition OR a ghost cell, then add the flux to the
            // boundary_flux_integral
            if( ((n<0) & (D->tri_full_flag[k]==1)) | ( (n>=0) && ((D->tri_full_flag[k]==1) & (D->tri_full_flag[n]==0) )) ){
                // boundary_flux_sum is an array with length = timestep_fluxcalls
                // For each sub-step, we put the boundary flux sum in.
                D->boundary_flux_sum[substep_count] += D->edge_flux_work[ki3];
            }

            D->xmom_explicit_update[k] -= D->normals[ki2]*D->pressuregrad_work[ki];
            D->ymom_explicit_update[k] -= D->normals[ki2+1]*D->pressuregrad_work[ki];


        } // end edge i

        // Normalise triangle k by area and store for when all conserved
        // quantities get updated
        inv_area = 1.0 / D->areas[k];
        D->stage_explicit_update[k] *= inv_area;
        D->xmom_explicit_update[k] *= inv_area;
        D->ymom_explicit_update[k] *= inv_area;

    }  // end cell k

    // Ensure we only update the timestep on the first call within each rk2/rk3 step
    if(substep_count == 0) timestep=local_timestep;

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
  double hc, bmin;
  double mass_error = 0.;
  // This acts like minimum_allowed height, but scales with the vertical
  // distance between the bed_centroid_value and the max bed_edge_value of
  // every triangle.
  //double minimum_relative_height=0.05;
  //int mass_added = 0;

  // Protect against inifintesimal and negative heights
  //if (maximum_allowed_speed < epsilon) {
    for (k=0; k<N; k++) {
      hc = wc[k] - zc[k];
      if (hc < minimum_allowed_height*1.0 ){
            // Set momentum to zero and ensure h is non negative
            xmomc[k] = 0.;
            ymomc[k] = 0.;
        if (hc <= 0.0){
             bmin = zc[k];
             // Minimum allowed stage = bmin

             // WARNING: ADDING MASS if wc[k]<bmin
             if(wc[k] < bmin){
                 mass_error += (bmin-wc[k])*areas[k];
                 //mass_added = 1; //Flag to warn of added mass

                 wc[k] = bmin;

                 // FIXME: Set vertex values as well. Seems that this shouldn't be
                 // needed. However, from memory this is important at the first
                 // time step, for 'dry' areas where the designated stage is
                 // less than the bed centroid value
                 wv[3*k] = bmin; //min(bmin, wc[k]); //zv[3*k]-minimum_allowed_height);
                 wv[3*k+1] = bmin; //min(bmin, wc[k]); //zv[3*k+1]-minimum_allowed_height);
                 wv[3*k+2] = bmin; //min(bmin, wc[k]); //zv[3*k+2]-minimum_allowed_height);
            }
        }
      }
    }

  //if(mass_added == 1){
  //  printf("Cumulative mass protection: %f m^3 \n", mass_error);
  //}

  return mass_error;
}

// Protect against the water elevation falling below the triangle bed
double  _protect_new(struct domain *D) {

  int k;
  double hc, bmin;
  double mass_error = 0.;

  double* wc;
  double* zc;
  double* wv;
  double* xmomc;
  double* ymomc;
  double* areas;

  double minimum_allowed_height;

  minimum_allowed_height = D->minimum_allowed_height;

  wc = D->stage_centroid_values;
  zc = D->bed_centroid_values;
  wv = D->stage_vertex_values;
  xmomc = D->xmom_centroid_values;
  ymomc = D->ymom_centroid_values;
  areas = D->areas;

  // This acts like minimum_allowed height, but scales with the vertical
  // distance between the bed_centroid_value and the max bed_edge_value of
  // every triangle.
  //double minimum_relative_height=0.05;
  //int mass_added = 0;

  // Protect against inifintesimal and negative heights
  //if (maximum_allowed_speed < epsilon) {
    for (k=0; k<D->number_of_elements; k++) {
      hc = wc[k] - zc[k];
      if (hc < minimum_allowed_height*1.0 ){
            // Set momentum to zero and ensure h is non negative
            xmomc[k] = 0.;
            ymomc[k] = 0.;
        if (hc <= 0.0){
             bmin = zc[k];
             // Minimum allowed stage = bmin

             // WARNING: ADDING MASS if wc[k]<bmin
             if(wc[k] < bmin){
                 mass_error += (bmin-wc[k])*areas[k];
                 //mass_added = 1; //Flag to warn of added mass

                 wc[k] = bmin;

                 // FIXME: Set vertex values as well. Seems that this shouldn't be
                 // needed. However, from memory this is important at the first
                 // time step, for 'dry' areas where the designated stage is
                 // less than the bed centroid value
                 wv[3*k] = bmin; //min(bmin, wc[k]); //zv[3*k]-minimum_allowed_height);
                 wv[3*k+1] = bmin; //min(bmin, wc[k]); //zv[3*k+1]-minimum_allowed_height);
                 wv[3*k+2] = bmin; //min(bmin, wc[k]); //zv[3*k+2]-minimum_allowed_height);
            }
        }
      }
    }

  //if(mass_added == 1){
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
  *qmax = fmax(fmax(dq0, fmax(dq0+dq1, dq0+dq2)), 0.0) ;
  *qmin = fmin(fmin(dq0, fmin(dq0+dq1, dq0+dq2)), 0.0) ;

  return 0;
}

int limit_gradient(double *dqv, double qmin, double qmax, double beta_w){
  // Given provisional jumps dqv from the FV triangle centroid to its
  // vertices/edges, and jumps qmin (qmax) between the centroid of the FV
  // triangle and the minimum (maximum) of the values at the auxiliary triangle
  // vertices (which are centroids of neighbour mesh triangles), calculate a
  // multiplicative factor phi by which the provisional vertex jumps are to be
  // limited

  int i;
  double r=1000.0, r0=1.0, phi=1.0;
  static double TINY = 1.0e-100; // to avoid machine accuracy problems.
  // FIXME: Perhaps use the epsilon used elsewhere.

  // Any provisional jump with magnitude < TINY does not contribute to
  // the limiting process.
  //return 0;

  for (i=0;i<3;i++){
    if (dqv[i] < -TINY)
      r0=qmin/dqv[i];

    if (dqv[i] > TINY)
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



// MIGRATED from shallow_water.c
// FIXME (Ole): Maybe superseded by extrapolate_second_order_edge_sw below?
// Computational routine
int _extrapolate_second_order_sw(struct domain *D) {


  // Domain Variables
    int number_of_elements;
    double epsilon;
    double minimum_allowed_height;
    double beta_w;
    double beta_w_dry;
    double beta_uh;
    double beta_uh_dry;
    double beta_vh;
    double beta_vh_dry;
    long* surrogate_neighbours;
    long* number_of_boundaries;
    double* centroid_coordinates;
    double* stage_centroid_values;
    double* xmom_centroid_values;
    double* ymom_centroid_values;
    double* bed_centroid_values;
    double* edge_coordinates;
    double* vertex_coordinates;
    double* stage_edge_values;
    double* xmom_edge_values;
    double* ymom_edge_values;
    double* bed_edge_values;
    double* stage_vertex_values;
    double* xmom_vertex_values;
    double* ymom_vertex_values;
    double* bed_vertex_values;
    int optimise_dry_cells;
    int extrapolate_velocity_second_order;

    // Local variables
    double a, b; // Gradient vector used to calculate edge values from centroids
    int k, k0, k1, k2, k3, k6, coord_index, i;
    double x, y, x0, y0, x1, y1, x2, y2, xv0, yv0, xv1, yv1, xv2, yv2; // Vertices of the auxiliary triangle
    double dx1, dx2, dy1, dy2, dxv0, dxv1, dxv2, dyv0, dyv1, dyv2, dq0, dq1, dq2, area2, inv_area2;
    double dqv[3], qmin, qmax, hmin, hmax;
    double hc, h0, h1, h2, beta_tmp, hfactor;
    //double dk, dv0, dv1, dv2, de[3], demin, dcmax, r0scale;
    double dk, dv0, dv1, dv2;

    double *xmom_centroid_store;
    double *ymom_centroid_store;
    //double *stage_centroid_store;


    // Associate memory location of Domain varibles with local aliases
    number_of_elements     = D->number_of_elements;
    epsilon                = D->epsilon;
    minimum_allowed_height = D->minimum_allowed_height;
    beta_w                 = D->beta_w;
    beta_w_dry             = D->beta_w_dry;
    beta_uh                = D->beta_uh;
    beta_uh_dry            = D->beta_uh_dry;
    beta_vh                = D->beta_vh;
    beta_vh_dry            = D->beta_vh_dry;
    optimise_dry_cells     = D->optimise_dry_cells;

    extrapolate_velocity_second_order = D->extrapolate_velocity_second_order;

    surrogate_neighbours      = D->surrogate_neighbours;
    number_of_boundaries      = D->number_of_boundaries;
    centroid_coordinates      = D->centroid_coordinates;
    stage_centroid_values     = D->stage_centroid_values;
    xmom_centroid_values      = D->xmom_centroid_values;
    ymom_centroid_values      = D->ymom_centroid_values;
    bed_centroid_values       = D->bed_centroid_values;
    edge_coordinates          = D->edge_coordinates;
    vertex_coordinates        = D->vertex_coordinates;
    stage_edge_values         = D->stage_edge_values;
    xmom_edge_values          = D->xmom_edge_values;
    ymom_edge_values          = D->ymom_edge_values;
    bed_edge_values           = D->bed_edge_values;
    stage_vertex_values       = D->stage_vertex_values;
    xmom_vertex_values        = D->xmom_vertex_values;
    ymom_vertex_values        = D->ymom_vertex_values;
    bed_vertex_values         = D->bed_vertex_values;




/*
int _extrapolate_second_order_sw(int number_of_elements,
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
        double* vertex_coordinates,
        double* stage_vertex_values,
        double* xmom_vertex_values,
        double* ymom_vertex_values,
        double* elevation_vertex_values,
        int optimise_dry_cells,
        int extrapolate_velocity_second_order) {



    // Local variables
    double a, b; // Gradient vector used to calculate vertex values from centroids
    int k, k0, k1, k2, k3, k6, coord_index, i;
    double x, y, x0, y0, x1, y1, x2, y2, xv0, yv0, xv1, yv1, xv2, yv2; // Vertices of the auxiliary triangle
    double dx1, dx2, dy1, dy2, dxv0, dxv1, dxv2, dyv0, dyv1, dyv2, dq0, dq1, dq2, area2, inv_area2;
    double dqv[3], qmin, qmax, hmin, hmax;
    double hc, h0, h1, h2, beta_tmp, hfactor;
    double xmom_centroid_store[number_of_elements], ymom_centroid_store[number_of_elements], dk, dv0, dv1, dv2;
*/

   // Use malloc to avoid putting these variables on the stack, which can cause
   // segfaults in large model runs
    xmom_centroid_store = malloc(number_of_elements*sizeof(double));
    ymom_centroid_store = malloc(number_of_elements*sizeof(double));
    // stage_centroid_store = malloc(number_of_elements*sizeof(double));

    if (extrapolate_velocity_second_order == 1) {
        // Replace momentum centroid with velocity centroid to allow velocity
        // extrapolation This will be changed back at the end of the routine
        for (k = 0; k < number_of_elements; k++) {

            dk = max(stage_centroid_values[k] - bed_centroid_values[k], minimum_allowed_height);
            xmom_centroid_store[k] = xmom_centroid_values[k];
            xmom_centroid_values[k] = xmom_centroid_values[k] / dk;

            ymom_centroid_store[k] = ymom_centroid_values[k];
            ymom_centroid_values[k] = ymom_centroid_values[k] / dk;
        }
    }

    // Begin extrapolation routine
    for (k = 0; k < number_of_elements; k++) {
        k3 = k * 3;
        k6 = k * 6;

        if (number_of_boundaries[k] == 3) {
            // No neighbours, set gradient on the triangle to zero

            stage_vertex_values[k3] = stage_centroid_values[k];
            stage_vertex_values[k3 + 1] = stage_centroid_values[k];
            stage_vertex_values[k3 + 2] = stage_centroid_values[k];
            xmom_vertex_values[k3] = xmom_centroid_values[k];
            xmom_vertex_values[k3 + 1] = xmom_centroid_values[k];
            xmom_vertex_values[k3 + 2] = xmom_centroid_values[k];
            ymom_vertex_values[k3] = ymom_centroid_values[k];
            ymom_vertex_values[k3 + 1] = ymom_centroid_values[k];
            ymom_vertex_values[k3 + 2] = ymom_centroid_values[k];

            continue;
        } else {
            // Triangle k has one or more neighbours.
            // Get centroid and vertex coordinates of the triangle

            // Get the vertex coordinates
            xv0 = vertex_coordinates[k6];
            yv0 = vertex_coordinates[k6 + 1];
            xv1 = vertex_coordinates[k6 + 2];
            yv1 = vertex_coordinates[k6 + 3];
            xv2 = vertex_coordinates[k6 + 4];
            yv2 = vertex_coordinates[k6 + 5];

            // Get the centroid coordinates
            coord_index = 2 * k;
            x = centroid_coordinates[coord_index];
            y = centroid_coordinates[coord_index + 1];

            // Store x- and y- differentials for the vertices of
            // triangle k relative to the centroid
            dxv0 = xv0 - x;
            dxv1 = xv1 - x;
            dxv2 = xv2 - x;
            dyv0 = yv0 - y;
            dyv1 = yv1 - y;
            dyv2 = yv2 - y;
        }




        if (number_of_boundaries[k] <= 1) {
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

            // Get the auxiliary triangle's vertex coordinates
            // (really the centroids of neighbouring triangles)
            coord_index = 2 * k0;
            x0 = centroid_coordinates[coord_index];
            y0 = centroid_coordinates[coord_index + 1];

            coord_index = 2 * k1;
            x1 = centroid_coordinates[coord_index];
            y1 = centroid_coordinates[coord_index + 1];

            coord_index = 2 * k2;
            x2 = centroid_coordinates[coord_index];
            y2 = centroid_coordinates[coord_index + 1];

            // Store x- and y- differentials for the vertices
            // of the auxiliary triangle
            dx1 = x1 - x0;
            dx2 = x2 - x0;
            dy1 = y1 - y0;
            dy2 = y2 - y0;

            // Calculate 2*area of the auxiliary triangle
            // The triangle is guaranteed to be counter-clockwise
            area2 = dy2 * dx1 - dy1*dx2;

            // If the mesh is 'weird' near the boundary,
            // the triangle might be flat or clockwise
            // Default to zero gradient
            if (area2 <= 0) {
                //printf("Error negative triangle area \n");
                //return -1;

                stage_vertex_values[k3] = stage_centroid_values[k];
                stage_vertex_values[k3 + 1] = stage_centroid_values[k];
                stage_vertex_values[k3 + 2] = stage_centroid_values[k];
                xmom_vertex_values[k3] = xmom_centroid_values[k];
                xmom_vertex_values[k3 + 1] = xmom_centroid_values[k];
                xmom_vertex_values[k3 + 2] = xmom_centroid_values[k];
                ymom_vertex_values[k3] = ymom_centroid_values[k];
                ymom_vertex_values[k3 + 1] = ymom_centroid_values[k];
                ymom_vertex_values[k3 + 2] = ymom_centroid_values[k];

                continue;
            }

            // Calculate heights of neighbouring cells
            hc = stage_centroid_values[k] - bed_centroid_values[k];
            h0 = stage_centroid_values[k0] - bed_centroid_values[k0];
            h1 = stage_centroid_values[k1] - bed_centroid_values[k1];
            h2 = stage_centroid_values[k2] - bed_centroid_values[k2];
            hmin = min(min(h0, min(h1, h2)), hc);
            //hfactor = hc/(hc + 1.0);

            hfactor = 0.0;
            if (hmin > 0.001) {
                hfactor = (hmin - 0.001) / (hmin + 0.004);
            }

            if (optimise_dry_cells) {
                // Check if linear reconstruction is necessary for triangle k
                // This check will exclude dry cells.

                hmax = max(h0, max(h1, h2));
                if (hmax < epsilon) {
                    continue;
                }
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

            inv_area2 = 1.0 / area2;
            // Calculate the gradient of stage on the auxiliary triangle
            a = dy2 * dq1 - dy1*dq2;
            a *= inv_area2;
            b = dx1 * dq2 - dx2*dq1;
            b *= inv_area2;

            // Calculate provisional jumps in stage from the centroid
            // of triangle k to its vertices, to be limited
            dqv[0] = a * dxv0 + b*dyv0;
            dqv[1] = a * dxv1 + b*dyv1;
            dqv[2] = a * dxv2 + b*dyv2;

            // Now we want to find min and max of the centroid and the
            // vertices of the auxiliary triangle and compute jumps
            // from the centroid to the min and max
            find_qmin_and_qmax(dq0, dq1, dq2, &qmin, &qmax);

            // Playing with dry wet interface
            //hmin = qmin;
            //beta_tmp = beta_w_dry;
            //if (hmin>minimum_allowed_height)
            beta_tmp = beta_w_dry + (beta_w - beta_w_dry) * hfactor;

            //printf("min_alled_height = %f\n",minimum_allowed_height);
            //printf("hmin = %f\n",hmin);
            //printf("beta_w = %f\n",beta_w);
            //printf("beta_tmp = %f\n",beta_tmp);
            // Limit the gradient
            limit_gradient(dqv, qmin, qmax, beta_tmp);

            //for (i=0;i<3;i++)
            stage_vertex_values[k3 + 0] = stage_centroid_values[k] + dqv[0];
            stage_vertex_values[k3 + 1] = stage_centroid_values[k] + dqv[1];
            stage_vertex_values[k3 + 2] = stage_centroid_values[k] + dqv[2];


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
            a = dy2 * dq1 - dy1*dq2;
            a *= inv_area2;
            b = dx1 * dq2 - dx2*dq1;
            b *= inv_area2;

            // Calculate provisional jumps in stage from the centroid
            // of triangle k to its vertices, to be limited
            dqv[0] = a * dxv0 + b*dyv0;
            dqv[1] = a * dxv1 + b*dyv1;
            dqv[2] = a * dxv2 + b*dyv2;

            // Now we want to find min and max of the centroid and the
            // vertices of the auxiliary triangle and compute jumps
            // from the centroid to the min and max
            find_qmin_and_qmax(dq0, dq1, dq2, &qmin, &qmax);
            //beta_tmp = beta_uh;
            //if (hmin<minimum_allowed_height)
            //beta_tmp = beta_uh_dry;
            beta_tmp = beta_uh_dry + (beta_uh - beta_uh_dry) * hfactor;

            // Limit the gradient
            limit_gradient(dqv, qmin, qmax, beta_tmp);

            for (i = 0; i < 3; i++) {
                xmom_vertex_values[k3 + i] = xmom_centroid_values[k] + dqv[i];
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
            a = dy2 * dq1 - dy1*dq2;
            a *= inv_area2;
            b = dx1 * dq2 - dx2*dq1;
            b *= inv_area2;

            // Calculate provisional jumps in stage from the centroid
            // of triangle k to its vertices, to be limited
            dqv[0] = a * dxv0 + b*dyv0;
            dqv[1] = a * dxv1 + b*dyv1;
            dqv[2] = a * dxv2 + b*dyv2;

            // Now we want to find min and max of the centroid and the
            // vertices of the auxiliary triangle and compute jumps
            // from the centroid to the min and max
            find_qmin_and_qmax(dq0, dq1, dq2, &qmin, &qmax);

            //beta_tmp = beta_vh;
            //
            //if (hmin<minimum_allowed_height)
            //beta_tmp = beta_vh_dry;
            beta_tmp = beta_vh_dry + (beta_vh - beta_vh_dry) * hfactor;

            // Limit the gradient
            limit_gradient(dqv, qmin, qmax, beta_tmp);

            for (i = 0; i < 3; i++) {
                ymom_vertex_values[k3 + i] = ymom_centroid_values[k] + dqv[i];
            }
        }// End number_of_boundaries <=1
        else {

            //==============================================
            // Number of boundaries == 2
            //==============================================

            // One internal neighbour and gradient is in direction of the neighbour's centroid

            // Find the only internal neighbour (k1?)
            for (k2 = k3; k2 < k3 + 3; k2++) {
                // Find internal neighbour of triangle k
                // k2 indexes the edges of triangle k

                if (surrogate_neighbours[k2] != k) {
                    break;
                }
            }

            if ((k2 == k3 + 3)) {
                // If we didn't find an internal neighbour
                return -1;
            }

            k1 = surrogate_neighbours[k2];

            // The coordinates of the triangle are already (x,y).
            // Get centroid of the neighbour (x1,y1)
            coord_index = 2 * k1;
            x1 = centroid_coordinates[coord_index];
            y1 = centroid_coordinates[coord_index + 1];

            // Compute x- and y- distances between the centroid of
            // triangle k and that of its neighbour
            dx1 = x1 - x;
            dy1 = y1 - y;

            // Set area2 as the square of the distance
            area2 = dx1 * dx1 + dy1*dy1;

            // Set dx2=(x1-x0)/((x1-x0)^2+(y1-y0)^2)
            // and dy2=(y1-y0)/((x1-x0)^2+(y1-y0)^2) which
            // respectively correspond to the x- and y- gradients
            // of the conserved quantities
            dx2 = 1.0 / area2;
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

            // Calculate provisional vertex jumps, to be limited
            dqv[0] = a * dxv0 + b*dyv0;
            dqv[1] = a * dxv1 + b*dyv1;
            dqv[2] = a * dxv2 + b*dyv2;

            // Now limit the jumps
            if (dq1 >= 0.0) {
                qmin = 0.0;
                qmax = dq1;
            } else {
                qmin = dq1;
                qmax = 0.0;
            }

            // Limit the gradient
            limit_gradient(dqv, qmin, qmax, beta_w);

            //for (i=0; i < 3; i++)
            //{
            stage_vertex_values[k3] = stage_centroid_values[k] + dqv[0];
            stage_vertex_values[k3 + 1] = stage_centroid_values[k] + dqv[1];
            stage_vertex_values[k3 + 2] = stage_centroid_values[k] + dqv[2];
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

            // Calculate provisional vertex jumps, to be limited
            dqv[0] = a * dxv0 + b*dyv0;
            dqv[1] = a * dxv1 + b*dyv1;
            dqv[2] = a * dxv2 + b*dyv2;

            // Now limit the jumps
            if (dq1 >= 0.0) {
                qmin = 0.0;
                qmax = dq1;
            } else {
                qmin = dq1;
                qmax = 0.0;
            }

            // Limit the gradient
            limit_gradient(dqv, qmin, qmax, beta_w);

            //for (i=0;i<3;i++)
            //xmom_vertex_values[k3] = xmom_centroid_values[k] + dqv[0];
            //xmom_vertex_values[k3 + 1] = xmom_centroid_values[k] + dqv[1];
            //xmom_vertex_values[k3 + 2] = xmom_centroid_values[k] + dqv[2];

            for (i = 0; i < 3; i++) {
                xmom_vertex_values[k3 + i] = xmom_centroid_values[k] + dqv[i];
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

            // Calculate provisional vertex jumps, to be limited
            dqv[0] = a * dxv0 + b*dyv0;
            dqv[1] = a * dxv1 + b*dyv1;
            dqv[2] = a * dxv2 + b*dyv2;

            // Now limit the jumps
            if (dq1 >= 0.0) {
                qmin = 0.0;
                qmax = dq1;
            }
            else {
                qmin = dq1;
                qmax = 0.0;
            }

            // Limit the gradient
            limit_gradient(dqv, qmin, qmax, beta_w);

            //for (i=0;i<3;i++)
            //ymom_vertex_values[k3] = ymom_centroid_values[k] + dqv[0];
            //ymom_vertex_values[k3 + 1] = ymom_centroid_values[k] + dqv[1];
            //ymom_vertex_values[k3 + 2] = ymom_centroid_values[k] + dqv[2];

            for (i = 0; i < 3; i++) {
                ymom_vertex_values[k3 + i] = ymom_centroid_values[k] + dqv[i];
            }
            //ymom_vertex_values[k3] = ymom_centroid_values[k] + dqv[0];
            //ymom_vertex_values[k3 + 1] = ymom_centroid_values[k] + dqv[1];
            //ymom_vertex_values[k3 + 2] = ymom_centroid_values[k] + dqv[2];
        } // else [number_of_boundaries==2]




    } // for k=0 to number_of_elements-1

    if (extrapolate_velocity_second_order == 1) {
        // Convert back from velocity to momentum
        for (k = 0; k < number_of_elements; k++) {
            k3 = 3 * k;
            //dv0 = max(stage_vertex_values[k3]-bed_vertex_values[k3],minimum_allowed_height);
            //dv1 = max(stage_vertex_values[k3+1]-bed_vertex_values[k3+1],minimum_allowed_height);
            //dv2 = max(stage_vertex_values[k3+2]-bed_vertex_values[k3+2],minimum_allowed_height);
            dv0 = max(stage_vertex_values[k3] - bed_vertex_values[k3], 0.);
            dv1 = max(stage_vertex_values[k3 + 1] - bed_vertex_values[k3 + 1], 0.);
            dv2 = max(stage_vertex_values[k3 + 2] - bed_vertex_values[k3 + 2], 0.);

            //Correct centroid and vertex values
            xmom_centroid_values[k] = xmom_centroid_store[k];
            xmom_vertex_values[k3] = xmom_vertex_values[k3] * dv0;
            xmom_vertex_values[k3 + 1] = xmom_vertex_values[k3 + 1] * dv1;
            xmom_vertex_values[k3 + 2] = xmom_vertex_values[k3 + 2] * dv2;

            ymom_centroid_values[k] = ymom_centroid_store[k];
            ymom_vertex_values[k3] = ymom_vertex_values[k3] * dv0;
            ymom_vertex_values[k3 + 1] = ymom_vertex_values[k3 + 1] * dv1;
            ymom_vertex_values[k3 + 2] = ymom_vertex_values[k3 + 2] * dv2;

        }
    }


    free(xmom_centroid_store);
    free(ymom_centroid_store);
    //free(stage_centroid_store);


    return 0;
}



// Computational routine
//int _extrapolate_second_order_edge_sw(int number_of_elements,
//                                 double epsilon,
//                                 double minimum_allowed_height,
//                                 double beta_w,
//                                 double beta_w_dry,
//                                 double beta_uh,
//                                 double beta_uh_dry,
//                                 double beta_vh,
//                                 double beta_vh_dry,
//                                 long* surrogate_neighbours,
//                                 long* neighbour_edges,
//                                 long* number_of_boundaries,
//                                 double* centroid_coordinates,
//                                 double* stage_centroid_values,
//                                 double* xmom_centroid_values,
//                                 double* ymom_centroid_values,
//                                 double* bed_centroid_values,
//                                 double* height_centroid_values,
//                                 double* edge_coordinates,
//                                 double* stage_edge_values,
//                                 double* xmom_edge_values,
//                                 double* ymom_edge_values,
//                                 double* bed_edge_values,
//                                 double* height_edge_values,
//                                 double* stage_vertex_values,
//                                 double* xmom_vertex_values,
//                                 double* ymom_vertex_values,
//                                 double* bed_vertex_values,
//                                 double* height_vertex_values,
//                                 int optimise_dry_cells,
//                                 int extrapolate_velocity_second_order,
//                                 double* x_centroid_work,
//                                 double* y_centroid_work,
//                                 long* update_extrapolation) {
int _extrapolate_second_order_edge_sw(struct domain *D) {

  // Local variables
  double a, b; // Gradient vector used to calculate edge values from centroids
  int k, k0, k1, k2, k3, k6, coord_index, i;
  double x, y, x0, y0, x1, y1, x2, y2, xv0, yv0, xv1, yv1, xv2, yv2; // Vertices of the auxiliary triangle
  double dx1, dx2, dy1, dy2, dxv0, dxv1, dxv2, dyv0, dyv1, dyv2, dq0, dq1, dq2, area2, inv_area2;
  double dqv[3], qmin, qmax, hmin, hmax;
  double hc, h0, h1, h2, beta_tmp, hfactor;
  double dk, dk_inv, a_tmp, b_tmp, c_tmp,d_tmp;


  memset((char*) D->x_centroid_work, 0, D->number_of_elements * sizeof (double));
  memset((char*) D->y_centroid_work, 0, D->number_of_elements * sizeof (double));

  // Parameters used to control how the limiter is forced to first-order near
  // wet-dry regions
  a_tmp = 0.3; // Highest depth ratio with hfactor=1
  b_tmp = 0.1; // Highest depth ratio with hfactor=0
  c_tmp = 1.0/(a_tmp-b_tmp);
  d_tmp = 1.0-(c_tmp*a_tmp);

  if(D->extrapolate_velocity_second_order==1){

      // Replace momentum centroid with velocity centroid to allow velocity
      // extrapolation This will be changed back at the end of the routine
      for (k=0; k< D->number_of_elements; k++){

          D->height_centroid_values[k] = fmax(D->stage_centroid_values[k] - D->bed_centroid_values[k], 0.);

          dk = D->height_centroid_values[k];
          if(dk> D->minimum_allowed_height){
              dk_inv=1.0/dk;
              D->x_centroid_work[k] = D->xmom_centroid_values[k];
              D->xmom_centroid_values[k] = D->xmom_centroid_values[k]*dk_inv;

              D->y_centroid_work[k] = D->ymom_centroid_values[k];
              D->ymom_centroid_values[k] = D->ymom_centroid_values[k]*dk_inv;
          }else{
              D->x_centroid_work[k] = 0.;
              D->xmom_centroid_values[k] = 0.;
              D->y_centroid_work[k] = 0.;
              D->ymom_centroid_values[k] = 0.;

         }
      }
  }

  // If a triangle is surrounded by dry cells (or dry cells + boundary
  // condition) set its momentum to zero too. This prevents 'pits' of
  // of water being trapped and unable to lose momentum, which can occur in
  // some situations
  for (k=0; k< D->number_of_elements;k++){

      k3=k*3;
      k0 = D->surrogate_neighbours[k3];
      k1 = D->surrogate_neighbours[k3 + 1];
      k2 = D->surrogate_neighbours[k3 + 2];

      if(( (D->height_centroid_values[k0] < D->minimum_allowed_height) | (k0==k) ) &
         ( (D->height_centroid_values[k1] < D->minimum_allowed_height) | (k1==k) ) &
         ( (D->height_centroid_values[k2] < D->minimum_allowed_height) | (k2==k) ) ) {
    	  	  //printf("Surrounded by dry cells\n");
              D->x_centroid_work[k] = 0.;
              D->xmom_centroid_values[k] = 0.;
              D->y_centroid_work[k] = 0.;
              D->ymom_centroid_values[k] = 0.;
      }


  }

  // Begin extrapolation routine
  for (k = 0; k < D->number_of_elements; k++)
  {

    // Don't update the extrapolation if the flux will not be computed on the
    // next timestep
    if(D->update_extrapolation[k]==0){
       continue;
    }


    // Useful indices
    k3=k*3;
    k6=k*6;

    if (D->number_of_boundaries[k]==3)
    {
      // No neighbours, set gradient on the triangle to zero

      D->stage_edge_values[k3]   = D->stage_centroid_values[k];
      D->stage_edge_values[k3+1] = D->stage_centroid_values[k];
      D->stage_edge_values[k3+2] = D->stage_centroid_values[k];

      //xmom_centroid_values[k] = 0.;
      //ymom_centroid_values[k] = 0.;

      D->xmom_edge_values[k3]    = D->xmom_centroid_values[k];
      D->xmom_edge_values[k3+1]  = D->xmom_centroid_values[k];
      D->xmom_edge_values[k3+2]  = D->xmom_centroid_values[k];
      D->ymom_edge_values[k3]    = D->ymom_centroid_values[k];
      D->ymom_edge_values[k3+1]  = D->ymom_centroid_values[k];
      D->ymom_edge_values[k3+2]  = D->ymom_centroid_values[k];

      dk = D->height_centroid_values[k];
      D->height_edge_values[k3] = dk;
      D->height_edge_values[k3+1] = dk;
      D->height_edge_values[k3+2] = dk;

      continue;
    }
    else
    {
      // Triangle k has one or more neighbours.
      // Get centroid and edge coordinates of the triangle

      // Get the edge coordinates
      xv0 = D->edge_coordinates[k6];
      yv0 = D->edge_coordinates[k6+1];
      xv1 = D->edge_coordinates[k6+2];
      yv1 = D->edge_coordinates[k6+3];
      xv2 = D->edge_coordinates[k6+4];
      yv2 = D->edge_coordinates[k6+5];

      // Get the centroid coordinates
      coord_index = 2*k;
      x = D->centroid_coordinates[coord_index];
      y = D->centroid_coordinates[coord_index+1];

      // Store x- and y- differentials for the edges of
      // triangle k relative to the centroid
      dxv0 = xv0 - x;
      dxv1 = xv1 - x;
      dxv2 = xv2 - x;
      dyv0 = yv0 - y;
      dyv1 = yv1 - y;
      dyv2 = yv2 - y;

    }



    if (D->number_of_boundaries[k]<=1)
    {
      //==============================================
      // Number of boundaries <= 1
      // 'Typical case'
      //==============================================


      // If no boundaries, auxiliary triangle is formed
      // from the centroids of the three neighbours
      // If one boundary, auxiliary triangle is formed
      // from this centroid and its two neighbours

      k0 = D->surrogate_neighbours[k3];
      k1 = D->surrogate_neighbours[k3 + 1];
      k2 = D->surrogate_neighbours[k3 + 2];


      // Get the auxiliary triangle's vertex coordinates
      // (really the centroids of neighbouring triangles)
      coord_index = 2*k0;
      x0 = D->centroid_coordinates[coord_index];
      y0 = D->centroid_coordinates[coord_index+1];

      coord_index = 2*k1;
      x1 = D->centroid_coordinates[coord_index];
      y1 = D->centroid_coordinates[coord_index+1];

      coord_index = 2*k2;
      x2 = D->centroid_coordinates[coord_index];
      y2 = D->centroid_coordinates[coord_index+1];

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
      if ((area2 <= 0.))
      {


          // Isolated wet cell -- constant stage/depth extrapolation
          D->stage_edge_values[k3]   = D->stage_centroid_values[k];
          D->stage_edge_values[k3+1] = D->stage_centroid_values[k];
          D->stage_edge_values[k3+2] = D->stage_centroid_values[k];

          dk= D->height_centroid_values[k]; //max(stage_centroid_values[k]-bed_centroid_values[k],0.);
          D->height_edge_values[k3] = dk;
          D->height_edge_values[k3+1] = dk;
          D->height_edge_values[k3+2] = dk;

          D->xmom_edge_values[k3]    = D->xmom_centroid_values[k];
          D->xmom_edge_values[k3+1]  = D->xmom_centroid_values[k];
          D->xmom_edge_values[k3+2]  = D->xmom_centroid_values[k];
          D->ymom_edge_values[k3]    = D->ymom_centroid_values[k];
          D->ymom_edge_values[k3+1]  = D->ymom_centroid_values[k];
          D->ymom_edge_values[k3+2]  = D->ymom_centroid_values[k];

          continue;
      }

      // Calculate heights of neighbouring cells
      hc = D->height_centroid_values[k];
      h0 = D->height_centroid_values[k0];
      h1 = D->height_centroid_values[k1];
      h2 = D->height_centroid_values[k2];

      hmin = fmin(fmin(h0, fmin(h1, h2)), hc);
      hmax = fmax(fmax(h0, fmax(h1, h2)), hc);

      // Look for strong changes in cell depth as an indicator of near-wet-dry
      // Reduce hfactor linearly from 1-0 between depth ratio (hmin/hc) of [a_tmp , b_tmp]
      // NOTE: If we have a more 'second order' treatment in near dry areas (e.g. with b_tmp being negative), then
      //       the water tends to dry more rapidly (which is in agreement with analytical results),
      //       but is also more 'artefacty' in important cases (tendency for high velocities, etc).
      //
      // So hfactor = depth_ratio*(c_tmp) + d_tmp, but is clipped between 0 and 1.
      hfactor= fmax(0., fmin(c_tmp*fmax(hmin,0.0)/fmax(hc,1.0e-06)+d_tmp,
                           fmin(c_tmp*fmax(hc,0.)/fmax(hmax,1.0e-06)+d_tmp, 1.0))
                  );
      // Set hfactor to zero smothly as hmin--> minimum_allowed_height. This
      // avoids some 'chatter' for very shallow flows
      hfactor=fmin( 1.2*fmax(hmin- D->minimum_allowed_height,0.)/(fmax(hmin,0.)+1.* D->minimum_allowed_height), hfactor);

      inv_area2 = 1.0/area2;
      //-----------------------------------
      // stage
      //-----------------------------------

      beta_tmp = D->beta_w_dry + (D->beta_w - D->beta_w_dry) * hfactor;

      if(beta_tmp>0.){
          // Calculate the difference between vertex 0 of the auxiliary
          // triangle and the centroid of triangle k
          dq0 = D->stage_centroid_values[k0] - D->stage_centroid_values[k];

          // Calculate differentials between the vertices
          // of the auxiliary triangle (centroids of neighbouring triangles)
          dq1 = D->stage_centroid_values[k1] - D->stage_centroid_values[k0];
          dq2 = D->stage_centroid_values[k2] - D->stage_centroid_values[k0];

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

          // Limit the gradient
          limit_gradient(dqv, qmin, qmax, beta_tmp);

          D->stage_edge_values[k3+0] = D->stage_centroid_values[k] + dqv[0];
          D->stage_edge_values[k3+1] = D->stage_centroid_values[k] + dqv[1];
          D->stage_edge_values[k3+2] = D->stage_centroid_values[k] + dqv[2];
      }else{
          // Fast alternative when beta_tmp==0
          D->stage_edge_values[k3+0] = D->stage_centroid_values[k];
          D->stage_edge_values[k3+1] = D->stage_centroid_values[k];
          D->stage_edge_values[k3+2] = D->stage_centroid_values[k];
      }


      //-----------------------------------
      // height
      //-----------------------------------

      if(beta_tmp>0.){
          // Calculate the difference between vertex 0 of the auxiliary
          // triangle and the centroid of triangle k
          dq0 = D->height_centroid_values[k0] - D->height_centroid_values[k];

          // Calculate differentials between the vertices
          // of the auxiliary triangle (centroids of neighbouring triangles)
          dq1 = D->height_centroid_values[k1] - D->height_centroid_values[k0];
          dq2 = D->height_centroid_values[k2] - D->height_centroid_values[k0];

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
          // Same beta_tmp as for stage
          //beta_tmp = beta_uh_dry + (beta_uh - beta_uh_dry) * hfactor;
          limit_gradient(dqv, qmin, qmax, beta_tmp);

          //beta_tmp = 0. + (beta_w - 0.) * hfactor;

          D->height_edge_values[k3+0] = D->height_centroid_values[k] + dqv[0];
          D->height_edge_values[k3+1] = D->height_centroid_values[k] + dqv[1];
          D->height_edge_values[k3+2] = D->height_centroid_values[k] + dqv[2];
      }else{
          // Fast alternative when beta_tmp==0
          D->height_edge_values[k3+0] = D->height_centroid_values[k];
          D->height_edge_values[k3+1] = D->height_centroid_values[k];
          D->height_edge_values[k3+2] = D->height_centroid_values[k];
      }
      //-----------------------------------
      // xmomentum
      //-----------------------------------

      beta_tmp = D->beta_uh_dry + (D->beta_uh - D->beta_uh_dry) * hfactor;
      if(beta_tmp>0.){
          // Calculate the difference between vertex 0 of the auxiliary
          // triangle and the centroid of triangle k
          dq0 = D->xmom_centroid_values[k0] - D->xmom_centroid_values[k];

          // Calculate differentials between the vertices
          // of the auxiliary triangle
          dq1 = D->xmom_centroid_values[k1] - D->xmom_centroid_values[k0];
          dq2 = D->xmom_centroid_values[k2] - D->xmom_centroid_values[k0];

          // Calculate the gradient of xmom on the auxiliary triangle
          a = dy2*dq1 - dy1*dq2;
          a *= inv_area2;
          b = dx1*dq2 - dx2*dq1;
          b *= inv_area2;

          // Calculate provisional jumps in xmom from the centroid
          // of triangle k to its vertices, to be limited
          dqv[0] = a*dxv0+b*dyv0;
          dqv[1] = a*dxv1+b*dyv1;
          dqv[2] = a*dxv2+b*dyv2;

          // Now we want to find min and max of the centroid and the
          // vertices of the auxiliary triangle and compute jumps
          // from the centroid to the min and max
          //
          find_qmin_and_qmax(dq0, dq1, dq2, &qmin, &qmax);


          // Limit the gradient
          limit_gradient(dqv, qmin, qmax, beta_tmp);

          for (i=0; i < 3; i++)
          {
            D->xmom_edge_values[k3+i] = D->xmom_centroid_values[k] + dqv[i];
          }
      }else{
          // Fast alternative when beta_tmp==0
          for (i=0; i < 3; i++)
          {
            D->xmom_edge_values[k3+i] = D->xmom_centroid_values[k];
          }
      }

      //-----------------------------------
      // ymomentum
      //-----------------------------------

      beta_tmp = D->beta_vh_dry + (D->beta_vh - D->beta_vh_dry) * hfactor;

      if(beta_tmp>0.){
          // Calculate the difference between vertex 0 of the auxiliary
          // triangle and the centroid of triangle k
          dq0 = D->ymom_centroid_values[k0] - D->ymom_centroid_values[k];

          // Calculate differentials between the vertices
          // of the auxiliary triangle
          dq1 = D->ymom_centroid_values[k1] - D->ymom_centroid_values[k0];
          dq2 = D->ymom_centroid_values[k2] - D->ymom_centroid_values[k0];

          // Calculate the gradient of xmom on the auxiliary triangle
          a = dy2*dq1 - dy1*dq2;
          a *= inv_area2;
          b = dx1*dq2 - dx2*dq1;
          b *= inv_area2;

          // Calculate provisional jumps in ymom from the centroid
          // of triangle k to its vertices, to be limited
          dqv[0] = a*dxv0 + b*dyv0;
          dqv[1] = a*dxv1 + b*dyv1;
          dqv[2] = a*dxv2 + b*dyv2;

          // Now we want to find min and max of the centroid and the
          // vertices of the auxiliary triangle and compute jumps
          // from the centroid to the min and max
          //
          find_qmin_and_qmax(dq0, dq1, dq2, &qmin, &qmax);


          // Limit the gradient
          limit_gradient(dqv, qmin, qmax, beta_tmp);

          for (i=0;i<3;i++)
          {
            D->ymom_edge_values[k3 + i] = D->ymom_centroid_values[k] + dqv[i];
          }
      }else{
          // Fast alternative when beta_tmp==0
          for (i=0;i<3;i++)
          {
            D->ymom_edge_values[k3 + i] = D->ymom_centroid_values[k];
          }

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

          if (D->surrogate_neighbours[k2] != k)
          {
             break;
          }
      }

      if ((k2 == k3 + 3))
      {
        // If we didn't find an internal neighbour
        //report_python_error(AT, "Internal neighbour not found");
        return -1;
      }

      k1 = D->surrogate_neighbours[k2];

      // The coordinates of the triangle are already (x,y).
      // Get centroid of the neighbour (x1,y1)
      coord_index = 2*k1;
      x1 = D->centroid_coordinates[coord_index];
      y1 = D->centroid_coordinates[coord_index + 1];

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
      dq1 = D->stage_centroid_values[k1] - D->stage_centroid_values[k];

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
      limit_gradient(dqv, qmin, qmax, D->beta_w);

      D->stage_edge_values[k3] = D->stage_centroid_values[k] + dqv[0];
      D->stage_edge_values[k3 + 1] = D->stage_centroid_values[k] + dqv[1];
      D->stage_edge_values[k3 + 2] = D->stage_centroid_values[k] + dqv[2];

      //-----------------------------------
      // height
      //-----------------------------------

      // Compute differentials
      dq1 = D->height_centroid_values[k1] - D->height_centroid_values[k];

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
      limit_gradient(dqv, qmin, qmax, D->beta_w);

      D->height_edge_values[k3] = D->height_centroid_values[k] + dqv[0];
      D->height_edge_values[k3 + 1] = D->height_centroid_values[k] + dqv[1];
      D->height_edge_values[k3 + 2] = D->height_centroid_values[k] + dqv[2];

      //-----------------------------------
      // xmomentum
      //-----------------------------------

      // Compute differentials
      dq1 = D->xmom_centroid_values[k1] - D->xmom_centroid_values[k];

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
      limit_gradient(dqv, qmin, qmax, D->beta_w);

      for (i = 0; i < 3;i++)
      {
          D->xmom_edge_values[k3 + i] = D->xmom_centroid_values[k] + dqv[i];
      }

      //-----------------------------------
      // ymomentum
      //-----------------------------------

      // Compute differentials
      dq1 = D->ymom_centroid_values[k1] - D->ymom_centroid_values[k];

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
      limit_gradient(dqv, qmin, qmax, D->beta_w);

      for (i=0;i<3;i++)
              {
              D->ymom_edge_values[k3 + i] = D->ymom_centroid_values[k] + dqv[i];
              }
    } // else [number_of_boundaries==2]
  } // for k=0 to number_of_elements-1


  // Compute vertex values of quantities
  for (k=0; k< D->number_of_elements; k++){
      if(D->extrapolate_velocity_second_order==1){
          //Convert velocity back to momenta at centroids
          D->xmom_centroid_values[k] = D->x_centroid_work[k];
          D->ymom_centroid_values[k] = D->y_centroid_work[k];
      }

      // Don't proceed if we didn't update the edge/vertex values
      if(D->update_extrapolation[k]==0){
         continue;
      }

      k3=3*k;

      // Compute stage vertex values
      D->stage_vertex_values[k3] = D->stage_edge_values[k3+1] + D->stage_edge_values[k3+2] - D->stage_edge_values[k3] ;
      D->stage_vertex_values[k3+1] =  D->stage_edge_values[k3] + D->stage_edge_values[k3+2]- D->stage_edge_values[k3+1];
      D->stage_vertex_values[k3+2] =  D->stage_edge_values[k3] + D->stage_edge_values[k3+1]- D->stage_edge_values[k3+2];

      // Compute height vertex values
      D->height_vertex_values[k3] = D->height_edge_values[k3+1] + D->height_edge_values[k3+2] - D->height_edge_values[k3] ;
      D->height_vertex_values[k3+1] =  D->height_edge_values[k3] + D->height_edge_values[k3+2]- D->height_edge_values[k3+1];
      D->height_vertex_values[k3+2] =  D->height_edge_values[k3] + D->height_edge_values[k3+1]- D->height_edge_values[k3+2];

      // If needed, convert from velocity to momenta
      if(D->extrapolate_velocity_second_order==1){
          // Re-compute momenta at edges
          for (i=0; i<3; i++){
              dk= D->height_edge_values[k3+i];
              D->xmom_edge_values[k3+i] = D->xmom_edge_values[k3+i]*dk;
              D->ymom_edge_values[k3+i] = D->ymom_edge_values[k3+i]*dk;
          }
      }
      // Compute momenta at vertices
      D->xmom_vertex_values[k3]   =  D->xmom_edge_values[k3+1] + D->xmom_edge_values[k3+2] - D->xmom_edge_values[k3] ;
      D->xmom_vertex_values[k3+1] =  D->xmom_edge_values[k3] + D->xmom_edge_values[k3+2]- D->xmom_edge_values[k3+1];
      D->xmom_vertex_values[k3+2] =  D->xmom_edge_values[k3] + D->xmom_edge_values[k3+1]- D->xmom_edge_values[k3+2];
      D->ymom_vertex_values[k3]   =  D->ymom_edge_values[k3+1] + D->ymom_edge_values[k3+2] - D->ymom_edge_values[k3] ;
      D->ymom_vertex_values[k3+1] =  D->ymom_edge_values[k3] + D->ymom_edge_values[k3+2]- D->ymom_edge_values[k3+1];
      D->ymom_vertex_values[k3+2] =  D->ymom_edge_values[k3] + D->ymom_edge_values[k3+1]- D->ymom_edge_values[k3+2];

      // Compute new bed elevation
      D->bed_edge_values[k3]= D->stage_edge_values[k3]- D->height_edge_values[k3];
      D->bed_edge_values[k3+1]= D->stage_edge_values[k3+1]- D->height_edge_values[k3+1];
      D->bed_edge_values[k3+2]= D->stage_edge_values[k3+2]- D->height_edge_values[k3+2];
      D->bed_vertex_values[k3] = D->bed_edge_values[k3+1] + D->bed_edge_values[k3+2] - D->bed_edge_values[k3] ;
      D->bed_vertex_values[k3+1] =  D->bed_edge_values[k3] + D->bed_edge_values[k3+2] - D->bed_edge_values[k3+1];
      D->bed_vertex_values[k3+2] =  D->bed_edge_values[k3] + D->bed_edge_values[k3+1] - D->bed_edge_values[k3+2];
  }

  return 0;
}

#include "Python.h"
#include "numpy/arrayobject.h"
#include "math.h"
#include <stdio.h>
const double pi = 3.14159265358979;


// Shared code snippets
#include "util_ext.h"


/* double max(double a, double b) { */
/* 	double z; */
/* 	z=(a>b)?a:b; */
/* 	return z;} */

/* double min(double a, double b) { */
/* 	double z; */
/* 	z=(a<b)?a:b; */
/* 	return z;} */


// Function to obtain speed from momentum and depth.
// This is used by flux functions
// Input parameters uh and h may be modified by this function.
double _compute_speed(double *uh,
		      double *h,
		      double epsilon,
		      double h0) {

  double u;

  if (*h < epsilon) {
    *h = 0.0;  //Could have been negative
     u = 0.0;
  } else {
    u = *uh/(*h + h0/ *h);
  }


  // Adjust momentum to be consistent with speed
  *uh = u * *h;

  return u;
}


//WELL BALANCED VERSION
//Innermost flux function (using w=z+h)
int _flux_function_channel21(double *q_leftm,double *q_leftp, double *q_rightm,
			   double *q_rightp, double g, double epsilon, double h0,                           double *edgeflux, double *max_speed) {

    int i;
    double flux_left[2], flux_right[2];
    double a_leftm,w_leftm, h_leftm, d_leftm, z_leftm, u_leftm, b_leftm, soundspeed_leftm;
     double a_leftp,w_leftp, h_leftp, d_leftp, z_leftp, u_leftp, b_leftp, soundspeed_leftp;
 double a_rightm,w_rightm, h_rightm, d_rightm, z_rightm, u_rightm, b_rightm, soundspeed_rightm;
 double a_rightp,w_rightp, h_rightp, d_rightp, z_rightp, u_rightp, b_rightp, soundspeed_rightp;
    double s_maxl, s_minl,s_maxr,s_minr, denom;
    double zphalf,zmhalf,hleftstar,hrightstar;
    double fluxtemp1,fluxtemp0,speedtemp;
    double batemp,bphalf,bmhalf;

    zmhalf = max(q_leftm[2],q_leftp[2]);
    zphalf = max(q_rightm[2],q_rightp[2]);   

    a_leftm  = q_leftm[0];
    d_leftm = q_leftm[1];
    z_leftm  = q_leftm[2];
    h_leftm  = max(0,q_leftm[3]+q_leftm[2]-zmhalf);
    u_leftm  = q_leftm[4];
    b_leftm  = q_leftm[5];
    w_leftm  = h_leftm+z_leftm;

    a_leftp  = q_leftp[0];
    d_leftp = q_leftp[1];
    z_leftp  = q_leftp[2];
    h_leftp  = max(0,q_leftp[3]+q_leftp[2]-zmhalf);
    u_leftp  = q_leftp[4];
    b_leftp  = q_leftp[5];
    w_leftp  = h_leftp+z_leftp;

    a_rightm  = q_rightm[0];
    d_rightm = q_rightm[1];
    z_rightm  = q_rightm[2];
    h_rightm  = max(0,q_rightm[3]+q_rightm[2]-zphalf);
    u_rightm  = q_rightm[4];
    b_rightm  = q_rightm[5];
    w_rightm  = h_rightm+z_rightm;

    a_rightp  = q_rightp[0];
    d_rightp = q_rightp[1];
    z_rightp  = q_rightp[2];
    h_rightp  = max(0,q_rightp[3]+q_rightp[2]-zphalf);
    u_rightp  = q_rightp[4];
    b_rightp  = q_rightp[5];
    w_rightp  = h_rightp+z_rightp;

    hleftstar = q_leftp[3];
    hrightstar = q_rightm[3];   

     bphalf = 0.5*(b_rightm+b_rightp);
     bmhalf = 0.5*(b_leftm+b_leftp);
     //bphalf = min(b_rightm,b_rightp);
     //bmhalf = min(b_leftm,b_leftp);


    soundspeed_leftp = sqrt(g*h_leftp);
    soundspeed_leftm = sqrt(g*h_leftm);
    soundspeed_rightp = sqrt(g*h_rightp);
    soundspeed_rightm = sqrt(g*h_rightm);

    s_maxl = max(u_leftm+soundspeed_leftm, u_leftp+soundspeed_leftp);
    if (s_maxl < 0.0) s_maxl = 0.0;

    s_minl = min(u_leftm-soundspeed_leftm, u_leftp-soundspeed_leftp);
    if (s_minl > 0.0) s_minl = 0.0;

    s_maxr = max(u_rightm+soundspeed_rightm, u_rightp+soundspeed_rightp);
    if (s_maxr < 0.0) s_maxr = 0.0;

    s_minr = min(u_rightm-soundspeed_rightm, u_rightp-soundspeed_rightp);
    if (s_minr > 0.0) s_minr = 0.0;

    // Flux formulas for left hand side
    flux_left[0] = d_leftm;
    flux_left[1] = u_leftm*d_leftm + 0.5*g*h_leftm*h_leftm*bmhalf;

    flux_right[0] = d_leftp;
    flux_right[1] = u_leftp*d_leftp + 0.5*g*h_leftp*h_leftp*bmhalf;

    
    // Flux computation for left hand side
    denom = s_maxl-s_minl;
    if (denom < epsilon) {
	for (i=0; i<2; i++) edgeflux[i] = 0.0;
    } else {
	edgeflux[0] = s_maxl*flux_left[0] - s_minl*flux_right[0];

        batemp = (q_leftp[3]+q_leftp[2])*bmhalf-(q_leftm[3]+q_leftm[2])*bmhalf;
       	edgeflux[0] += s_maxl*s_minl*batemp;
	edgeflux[0] /= denom;
	edgeflux[1] = s_maxl*flux_left[1] - s_minl*flux_right[1];
	edgeflux[1] += s_maxl*s_minl*(d_leftp-d_leftm);
	edgeflux[1] /= denom;
    

    }
	fluxtemp0 = edgeflux[0];
        fluxtemp1 = edgeflux[1];
    

    // Flux formulas for right hand side
    flux_left[0] = d_rightm;
    flux_left[1] = u_rightm*d_rightm + 0.5*g*h_rightm*h_rightm*bphalf;

    flux_right[0] = d_rightp;
    flux_right[1] = u_rightp*d_rightp + 0.5*g*h_rightp*h_rightp*bphalf;
    
   
   
    // Flux computation for right hand side
    denom = s_maxr-s_minr;
    if (denom < epsilon) {
	for (i=0; i<2; i++) edgeflux[i] = 0.0;
    } else {
	edgeflux[0] = s_maxr*flux_left[0] - s_minr*flux_right[0];

        batemp = (q_rightp[3]+q_rightp[2])*bphalf-(q_rightm[3]+q_rightm[2])*bphalf;

       	edgeflux[0] += s_maxr*s_minr*batemp;
	edgeflux[0] /= denom;
	edgeflux[1] = s_maxr*flux_left[1] - s_minr*flux_right[1];
	edgeflux[1] += s_maxr*s_minr*(d_rightp-d_rightm);
	edgeflux[1] /= denom;
    

    }
    
   
    edgeflux[0]=edgeflux[0]-fluxtemp0;
    edgeflux[1]=edgeflux[1]-fluxtemp1;
    
   
    edgeflux[1]-=0.5*g*h_rightm*h_rightm*bphalf-0.5*g*hrightstar*hrightstar*b_rightm+0.5*g*hleftstar*hleftstar*b_leftp-0.5*g*h_leftp*h_leftp*bmhalf;

    // printf("edgflux:%f expected:%f \n",edgeflux[1],hrightstar*hrightstar*g*0.5*b_rightm-hleftstar*hleftstar*g*0.5*b_leftp);

    edgeflux[1]-=g*(1.0/6.0)*(b_rightm*(hleftstar*hleftstar+hrightstar*(hrightstar+2*z_leftp-2*z_rightm)+hleftstar*(hrightstar+z_leftp-z_rightm))-b_leftp*(hleftstar*hleftstar+hrightstar*(hrightstar-z_leftp+z_rightm)+hleftstar*(hrightstar-2*z_leftp+2*z_rightm)));

   
    //edgeflux[1]-=0.5*g*h_rightm*h_rightm-0.5*g*hrightstar*hrightstar+0.5*g*hleftstar*hleftstar-0.5*g*h_leftp*h_leftp;
      
     //edgeflux[1]-=0.5*g*b_rightm*h_rightm*h_rightm-0.5*g*b_leftp*h_leftp*h_leftp;
         	// Maximal wavespeed
	if ( (s_maxl-s_minl)<epsilon && (s_maxr-s_minr)<epsilon ){
	    *max_speed = 0.0;
	}else{
	speedtemp = max(fabs(s_maxl),fabs(s_minl));
        speedtemp = max(speedtemp,fabs(s_maxr));
        speedtemp = max(speedtemp,fabs(s_minr));
	*max_speed = speedtemp;
	}

//printf("%f\n",h_right);
    return 0;
}
// GOOD BUT NOT WELL BALANCED VERSION
int _flux_function_channel(double *q_leftm,double *q_leftp, double *q_rightm,
			   double *q_rightp, double g, double epsilon, double h0,                           double *edgeflux, double *max_speed){

    int i;
    double flux_left[2], flux_right[2];
    double a_leftm,w_leftm, h_leftm, d_leftm, z_leftm, u_leftm, b_leftm, soundspeed_leftm;
     double a_leftp,w_leftp, h_leftp, d_leftp, z_leftp, u_leftp, b_leftp, soundspeed_leftp;
 double a_rightm,w_rightm, h_rightm, d_rightm, z_rightm, u_rightm, b_rightm, soundspeed_rightm;
 double a_rightp,w_rightp, h_rightp, d_rightp, z_rightp, u_rightp, b_rightp, soundspeed_rightp;
    double s_maxl, s_minl,s_maxr,s_minr, denom;
    double zphalf,zmhalf,hleftstar,hrightstar;
    double fluxtemp1,fluxtemp0,speedtemp;
    double batemp,bphalf,bmhalf;

    zmhalf = max(q_leftm[2],q_leftp[2]);
    zphalf = max(q_rightm[2],q_rightp[2]);   

    a_leftm  = q_leftm[0];
    d_leftm = q_leftm[1];
    z_leftm  = q_leftm[2];
    h_leftm  = q_leftm[3];
    u_leftm  = q_leftm[4];
    b_leftm  = q_leftm[5];
    w_leftm  = h_leftm+z_leftm;

    a_leftp  = q_leftp[0];
    d_leftp = q_leftp[1];
    z_leftp  = q_leftp[2];
    h_leftp  = q_leftp[3];
    u_leftp  = q_leftp[4];
    b_leftp  = q_leftp[5];
    w_leftp  = h_leftp+z_leftp;

    a_rightm  = q_rightm[0];
    d_rightm = q_rightm[1];
    z_rightm  = q_rightm[2];
    h_rightm  = q_rightm[3];
    u_rightm  = q_rightm[4];
    b_rightm  = q_rightm[5];
    w_rightm  = h_rightm+z_rightm;

    a_rightp  = q_rightp[0];
    d_rightp = q_rightp[1];
    z_rightp  = q_rightp[2];
    h_rightp  = q_rightp[3];
    u_rightp  = q_rightp[4];
    b_rightp  = q_rightp[5];
    w_rightp  = h_rightp+z_rightp;

    hleftstar = q_leftp[3];
    hrightstar = q_rightm[3];   

    bphalf = min(b_rightm,b_rightp);
    bmhalf = min(b_leftm,b_leftp);

    soundspeed_leftp = sqrt(g*h_leftp);
    soundspeed_leftm = sqrt(g*h_leftm);
    soundspeed_rightp = sqrt(g*h_rightp);
    soundspeed_rightm = sqrt(g*h_rightm);

    s_maxl = max(u_leftm+soundspeed_leftm, u_leftp+soundspeed_leftp);
    if (s_maxl < 0.0) s_maxl = 0.0;

    s_minl = min(u_leftm-soundspeed_leftm, u_leftp-soundspeed_leftp);
    if (s_minl > 0.0) s_minl = 0.0;

    s_maxr = max(u_rightm+soundspeed_rightm, u_rightp+soundspeed_rightp);
    if (s_maxr < 0.0) s_maxr = 0.0;

    s_minr = min(u_rightm-soundspeed_rightm, u_rightp-soundspeed_rightp);
    if (s_minr > 0.0) s_minr = 0.0;

    // Flux formulas for left hand side
    flux_left[0] = d_leftm;
    flux_left[1] = u_leftm*d_leftm + 0.5*g*h_leftm*h_leftm*b_leftm;

    flux_right[0] = d_leftp;
    flux_right[1] = u_leftp*d_leftp + 0.5*g*h_leftp*h_leftp*b_leftp;

    
    // Flux computation for left hand side
    denom = s_maxl-s_minl;
    if (denom < epsilon) {
	for (i=0; i<2; i++) edgeflux[i] = 0.0;
    } else {
	edgeflux[0] = s_maxl*flux_left[0] - s_minl*flux_right[0];

        batemp = (q_leftp[3]+q_leftp[2])*b_leftp-(q_leftm[3]+q_leftm[2])*b_leftm;
       	edgeflux[0] += s_maxl*s_minl*batemp;
	edgeflux[0] /= denom;
	edgeflux[1] = s_maxl*flux_left[1] - s_minl*flux_right[1];
	edgeflux[1] += s_maxl*s_minl*(d_leftp-d_leftm);
	edgeflux[1] /= denom;
    

    }
	fluxtemp0 = edgeflux[0];
        fluxtemp1 = edgeflux[1];
    

    // Flux formulas for right hand side
    flux_left[0] = d_rightm;
    flux_left[1] = u_rightm*d_rightm + 0.5*g*h_rightm*h_rightm*b_rightm;

    flux_right[0] = d_rightp;
    flux_right[1] = u_rightp*d_rightp + 0.5*g*h_rightp*h_rightp*b_rightp;
    
   
   
    // Flux computation for right hand side
    denom = s_maxr-s_minr;
    if (denom < epsilon) {
	for (i=0; i<2; i++) edgeflux[i] = 0.0;
    } else {
	edgeflux[0] = s_maxr*flux_left[0] - s_minr*flux_right[0];

        batemp = (q_rightp[3]+q_rightp[2])*b_rightp-(q_rightm[3]+q_rightm[2])*b_rightm;

       	edgeflux[0] += s_maxr*s_minr*batemp;
	edgeflux[0] /= denom;
	edgeflux[1] = s_maxr*flux_left[1] - s_minr*flux_right[1];
	edgeflux[1] += s_maxr*s_minr*(d_rightp-d_rightm);
	edgeflux[1] /= denom;
    

    }
    
   
    edgeflux[0]=edgeflux[0]-fluxtemp0;
    edgeflux[1]=edgeflux[1]-fluxtemp1;
    
    edgeflux[1]-=-0.5*0.5*g*(h_rightm+h_leftp)*(b_rightm+b_leftp)*(z_rightm-z_leftp)+0.5*(h_rightm+h_leftp)*(h_rightm+h_leftp)*0.5*0.5*(b_rightm-b_leftp)*g; 

    //edgeflux[1]-=0.5*g*h_rightm*h_rightm*bphalf-0.5*g*hrightstar*hrightstar*b_rightm+0.5*g*hleftstar*hleftstar*b_leftp-0.5*g*h_leftp*h_leftp*bmhalf;

    // printf("edgflux:%f expected:%f \n",edgeflux[1],hrightstar*hrightstar*g*0.5*b_rightm-hleftstar*hleftstar*g*0.5*b_leftp);

    //edgeflux[1]-=g*(1.0/6.0)*(b_rightm*(hleftstar*hleftstar+hrightstar*(hrightstar+2*z_leftp-2*z_rightm)+hleftstar*(hrightstar+z_leftp-z_rightm))-b_leftp*(hleftstar*hleftstar+hrightstar*(hrightstar-z_leftp+z_rightm)+hleftstar*(hrightstar-2*z_leftp+2*z_rightm)));

   
    //edgeflux[1]-=0.5*g*h_rightm*h_rightm-0.5*g*hrightstar*hrightstar+0.5*g*hleftstar*hleftstar-0.5*g*h_leftp*h_leftp;
      
     //edgeflux[1]-=0.5*g*b_rightm*h_rightm*h_rightm-0.5*g*b_leftp*h_leftp*h_leftp;
         	// Maximal wavespeed
	if ( (s_maxl-s_minl)<epsilon && (s_maxr-s_minr)<epsilon ){
	    *max_speed = 0.0;
	}else{
	speedtemp = max(fabs(s_maxl),fabs(s_minl));
        speedtemp = max(speedtemp,fabs(s_maxr));
        speedtemp = max(speedtemp,fabs(s_minr));
	*max_speed = speedtemp;
	}

//printf("%f\n",h_right);
    return 0;
}
// NAIEVE VERSION
int _flux_function_channel2(double *q_leftm,double *q_leftp, double *q_rightm,
			   double *q_rightp, double g, double epsilon, double h0,                           double *edgeflux, double *max_speed){

    int i;
    double flux_left[2], flux_right[2];
    double a_leftm,w_leftm, h_leftm, d_leftm, z_leftm, u_leftm, b_leftm, soundspeed_leftm;
     double a_leftp,w_leftp, h_leftp, d_leftp, z_leftp, u_leftp, b_leftp, soundspeed_leftp;
 double a_rightm,w_rightm, h_rightm, d_rightm, z_rightm, u_rightm, b_rightm, soundspeed_rightm;
 double a_rightp,w_rightp, h_rightp, d_rightp, z_rightp, u_rightp, b_rightp, soundspeed_rightp;
    double s_maxl, s_minl,s_maxr,s_minr, denom;
    double zphalf,zmhalf,hleftstar,hrightstar;
    double fluxtemp1,fluxtemp0,speedtemp;
    double batemp,bphalf,bmhalf;

    zmhalf = max(q_leftm[2],q_leftp[2]);
    zphalf = max(q_rightm[2],q_rightp[2]);   

    a_leftm  = q_leftm[0];
    d_leftm = q_leftm[1];
    z_leftm  = q_leftm[2];
    h_leftm  = q_leftm[3];
    u_leftm  = q_leftm[4];
    b_leftm  = q_leftm[5];
    w_leftm  = h_leftm+z_leftm;

    a_leftp  = q_leftp[0];
    d_leftp = q_leftp[1];
    z_leftp  = q_leftp[2];
    h_leftp  = q_leftp[3];
    u_leftp  = q_leftp[4];
    b_leftp  = q_leftp[5];
    w_leftp  = h_leftp+z_leftp;

    a_rightm  = q_rightm[0];
    d_rightm = q_rightm[1];
    z_rightm  = q_rightm[2];
    h_rightm  = q_rightm[3];
    u_rightm  = q_rightm[4];
    b_rightm  = q_rightm[5];
    w_rightm  = h_rightm+z_rightm;

    a_rightp  = q_rightp[0];
    d_rightp = q_rightp[1];
    z_rightp  = q_rightp[2];
    h_rightp  = q_rightp[3];
    u_rightp  = q_rightp[4];
    b_rightp  = q_rightp[5];
    w_rightp  = h_rightp+z_rightp;

    hleftstar = q_leftp[3];
    hrightstar = q_rightm[3];   

    bphalf = min(b_rightm,b_rightp);
    bmhalf = min(b_leftm,b_leftp);

    soundspeed_leftp = sqrt(g*h_leftp);
    soundspeed_leftm = sqrt(g*h_leftm);
    soundspeed_rightp = sqrt(g*h_rightp);
    soundspeed_rightm = sqrt(g*h_rightm);

    s_maxl = max(u_leftm+soundspeed_leftm, u_leftp+soundspeed_leftp);
    if (s_maxl < 0.0) s_maxl = 0.0;

    s_minl = min(u_leftm-soundspeed_leftm, u_leftp-soundspeed_leftp);
    if (s_minl > 0.0) s_minl = 0.0;

    s_maxr = max(u_rightm+soundspeed_rightm, u_rightp+soundspeed_rightp);
    if (s_maxr < 0.0) s_maxr = 0.0;

    s_minr = min(u_rightm-soundspeed_rightm, u_rightp-soundspeed_rightp);
    if (s_minr > 0.0) s_minr = 0.0;

    // Flux formulas for left hand side
    flux_left[0] = d_leftm;
    flux_left[1] = u_leftm*d_leftm + 0.5*g*h_leftm*h_leftm*b_leftm;

    flux_right[0] = d_leftp;
    flux_right[1] = u_leftp*d_leftp + 0.5*g*h_leftp*h_leftp*b_leftp;

    
    // Flux computation for left hand side
    denom = s_maxl-s_minl;
    if (denom < epsilon) {
	for (i=0; i<2; i++) edgeflux[i] = 0.0;
    } else {
	edgeflux[0] = s_maxl*flux_left[0] - s_minl*flux_right[0];

        batemp = (q_leftp[3]+q_leftp[2])*b_leftp-(q_leftm[3]+q_leftm[2])*b_leftm;
       
	edgeflux[0] = 0.5*(flux_left[0]+flux_right[0]);
	edgeflux[1] = 0.5*(flux_left[1]+flux_right[1]);

    

    }
	fluxtemp0 = edgeflux[0];
        fluxtemp1 = edgeflux[1];
    

    // Flux formulas for right hand side
    flux_left[0] = d_rightm;
    flux_left[1] = u_rightm*d_rightm + 0.5*g*h_rightm*h_rightm*b_rightm;

    flux_right[0] = d_rightp;
    flux_right[1] = u_rightp*d_rightp + 0.5*g*h_rightp*h_rightp*b_rightp;
    
   
   
    // Flux computation for right hand side
    denom = s_maxr-s_minr;
    if (denom < epsilon) {
	for (i=0; i<2; i++) edgeflux[i] = 0.0;
    } else {
	edgeflux[0] = s_maxr*flux_left[0] - s_minr*flux_right[0];

        batemp = (q_rightp[3]+q_rightp[2])*b_rightp-(q_rightm[3]+q_rightm[2])*b_rightm;

     
	edgeflux[0] = 0.5*(flux_right[0]+flux_left[0]);


	edgeflux[1] = 0.5*(flux_left[1]+flux_right[1]);
    

    }
    
   
    edgeflux[0]=edgeflux[0]-fluxtemp0;
    edgeflux[1]=edgeflux[1]-fluxtemp1;
    
    edgeflux[1]-=-0.5*0.5*g*(h_rightm+h_leftp)*(b_rightm+b_leftp)*(z_rightm-z_leftp)+0.5*(h_rightm+h_leftp)*(h_rightm+h_leftp)*0.5*0.5*(b_rightm-b_leftp)*g; 

    //edgeflux[1]-=0.5*g*h_rightm*h_rightm*bphalf-0.5*g*hrightstar*hrightstar*b_rightm+0.5*g*hleftstar*hleftstar*b_leftp-0.5*g*h_leftp*h_leftp*bmhalf;

    // printf("edgflux:%f expected:%f \n",edgeflux[1],hrightstar*hrightstar*g*0.5*b_rightm-hleftstar*hleftstar*g*0.5*b_leftp);

    //edgeflux[1]-=g*(1.0/6.0)*(b_rightm*(hleftstar*hleftstar+hrightstar*(hrightstar+2*z_leftp-2*z_rightm)+hleftstar*(hrightstar+z_leftp-z_rightm))-b_leftp*(hleftstar*hleftstar+hrightstar*(hrightstar-z_leftp+z_rightm)+hleftstar*(hrightstar-2*z_leftp+2*z_rightm)));

   
    //edgeflux[1]-=0.5*g*h_rightm*h_rightm-0.5*g*hrightstar*hrightstar+0.5*g*hleftstar*hleftstar-0.5*g*h_leftp*h_leftp;
      
     //edgeflux[1]-=0.5*g*b_rightm*h_rightm*h_rightm-0.5*g*b_leftp*h_leftp*h_leftp;
         	// Maximal wavespeed
	if ( (s_maxl-s_minl)<epsilon && (s_maxr-s_minr)<epsilon ){
	    *max_speed = 0.0;
	}else{
	speedtemp = max(fabs(s_maxl),fabs(s_minl));
        speedtemp = max(speedtemp,fabs(s_maxr));
        speedtemp = max(speedtemp,fabs(s_minr));
	*max_speed = speedtemp;
	}

//printf("%f\n",h_right);
    return 0;
}


// Computational function for flux computation
double _compute_fluxes_channel_ext(double cfl,
			       double timestep,
			       double epsilon,
			       double g,
			       double h0,
			       long* neighbours,
			       long* neighbour_vertices,
			       double* normals,
			       double* areas,
			       double* area_edge_values,
			       double* discharge_edge_values,
			       double* bed_edge_values,
			       double* height_edge_values,
			       double* velocity_edge_values,
			       double* width_edge_values,
			       double* area_boundary_values,
			       double* discharge_boundary_values,
			       double* bed_boundary_values,
			       double* height_boundary_values,
			       double* velocity_boundary_values,
			       double* width_boundary_values,
			       double* area_explicit_update,
			       double* discharge_explicit_update,
			       int number_of_elements,
			       double* max_speed_array) {

    double flux[2], qlm[6],qlp[6], qrm[6],qrp[6], edgeflux[2];
    double max_speed, normal;
    int k, i, ki, n, m, nm=0;
    double zstar;
    for (k=0; k<number_of_elements; k++) {
	flux[0] = 0.0;
	flux[1] = 0.0;

       
	    ki = k*2;
           

         n = neighbours[ki];
	    if (n<0) {
		m = -n-1;

	    qlm[0] = area_boundary_values[m];
	    qlm[1] = discharge_boundary_values[m];
	    qlm[2] = bed_boundary_values[m];
	    qlm[3] = height_boundary_values[m];
	    qlm[4] = velocity_boundary_values[m];
            qlm[5] = width_boundary_values[m];

	    }else{
	m = neighbour_vertices[ki];
		nm = n*2+m;
	    

	    qlm[0] = area_edge_values[nm];
	    qlm[1] = discharge_edge_values[nm];
	    qlm[2] = bed_edge_values[nm];
	    qlm[3] = height_edge_values[nm];
	    qlm[4] = velocity_edge_values[nm];
            qlm[5] = width_edge_values[nm];
    }
	    qlp[0] = area_edge_values[ki];
	    qlp[1] = discharge_edge_values[ki];
	    qlp[2] = bed_edge_values[ki];
	    qlp[3] = height_edge_values[ki];
	    qlp[4] = velocity_edge_values[ki];
            qlp[5] = width_edge_values[ki];

         ki = k*2+1;
           

         n = neighbours[ki];
	    if (n<0) {
		m = -n-1;
            qrp[0] = area_boundary_values[m];
	    qrp[1] = discharge_boundary_values[m];
	    qrp[2] = bed_boundary_values[m];
	    qrp[3] = height_boundary_values[m];
	    qrp[4] = velocity_boundary_values[m];
            qrp[5] = width_boundary_values[m];
          


	    }else{
	m = neighbour_vertices[ki];
		nm = n*2+m;
	    
           
            qrp[0] = area_edge_values[nm];
	    qrp[1] = discharge_edge_values[nm];
	    qrp[2] = bed_edge_values[nm];
	    qrp[3] = height_edge_values[nm];
	    qrp[4] = velocity_edge_values[nm];
            qrp[5] = width_edge_values[nm];
            }
            qrm[0] = area_edge_values[ki];
	    qrm[1] = discharge_edge_values[ki];
	    qrm[2] = bed_edge_values[ki];
	    qrm[3] = height_edge_values[ki];
	    qrm[4] = velocity_edge_values[ki];
            qrm[5] = width_edge_values[ki];

             _flux_function_channel21(qlm,qlp,qrm,qrp,g,epsilon,h0,edgeflux,&max_speed);
	    flux[0] -= edgeflux[0];
	    flux[1] -= edgeflux[1];
           
	    // Update timestep based on edge i and possibly neighbour n
	    if (max_speed > epsilon) {
		// Original CFL calculation

		timestep = min(timestep, 0.5*cfl*areas[k]/max_speed);
		if (n>=0) {
		    timestep = min(timestep, 0.5*cfl*areas[n]/max_speed);
		}
	    }
	 // End edge i (and neighbour n)
	flux[0] /= areas[k];
	area_explicit_update[k] = flux[0];
	flux[1] /= areas[k];
	discharge_explicit_update[k] = flux[1];
	//Keep track of maximal speeds
	max_speed_array[k]=max_speed;
    }
    return timestep;

}


//-------------------------------------------------------------
// Old code
//------------------------------------------------------------
//Innermost flux function (using w=z+h)





// Computational function for flux computation


//=========================================================================
// Python Glue
//=========================================================================



//------------------------------------------------
// New velocity based compute fluxes
//------------------------------------------------

PyObject *compute_fluxes_channel_ext(PyObject *self, PyObject *args) {

    PyObject
	*domain,
	*area,
	*discharge,
	*bed,
	*height,
	*velocity,
        *width;

    PyArrayObject
	*neighbours,
	*neighbour_vertices,
	*normals,
	*areas,
	*area_vertex_values,
	*discharge_vertex_values,
	*bed_vertex_values,
	*height_vertex_values,
	*velocity_vertex_values,
	*width_vertex_values,
	*area_boundary_values,
	*discharge_boundary_values,
	*bed_boundary_values,
	*height_boundary_values,
	*velocity_boundary_values,
	*width_boundary_values,
	*area_explicit_update,
	*discharge_explicit_update,
	*max_speed_array;

  double timestep, epsilon, g, h0, cfl;
  int number_of_elements;


  // Convert Python arguments to C
  if (!PyArg_ParseTuple(args, "dOOOOOOO",
			&timestep,
			&domain,
			&area,
			&discharge,
			&bed,
			&height,
			&velocity,
                        &width)) {
      PyErr_SetString(PyExc_RuntimeError, "comp_flux_channel_ext.c: compute_fluxes_channel_ext could not parse input");
      return NULL;
  }


    epsilon           = get_python_double(domain,"epsilon");
    g                 = get_python_double(domain,"g");
    h0                = get_python_double(domain,"h0");
    cfl               = get_python_double(domain,"CFL");


    neighbours        = get_consecutive_array(domain, "neighbours");
    neighbour_vertices= get_consecutive_array(domain, "neighbour_vertices");
    normals           = get_consecutive_array(domain, "normals");
    areas             = get_consecutive_array(domain, "areas");
    max_speed_array   = get_consecutive_array(domain, "max_speed_array");

    area_vertex_values      = get_consecutive_array(area,    "vertex_values");
    discharge_vertex_values       = get_consecutive_array(discharge,     "vertex_values");
    bed_vertex_values        = get_consecutive_array(bed,      "vertex_values");
    height_vertex_values     = get_consecutive_array(height,   "vertex_values");
    velocity_vertex_values   = get_consecutive_array(velocity, "vertex_values");
    width_vertex_values      = get_consecutive_array(width, "vertex_values");

    area_boundary_values     = get_consecutive_array(area,     "boundary_values");
    discharge_boundary_values      = get_consecutive_array(discharge,      "boundary_values");
    bed_boundary_values       = get_consecutive_array(bed,       "boundary_values");
    height_boundary_values    = get_consecutive_array(height,    "boundary_values");
    velocity_boundary_values  = get_consecutive_array(velocity,  "boundary_values");
    width_boundary_values       = get_consecutive_array(width,     "boundary_values");


    area_explicit_update = get_consecutive_array(area, "explicit_update");
    discharge_explicit_update  = get_consecutive_array(discharge,  "explicit_update");

    number_of_elements = area_vertex_values -> dimensions[0];

    // Call underlying flux computation routine and update
    // the explicit update arrays
    timestep = _compute_fluxes_channel_ext(cfl,
				       timestep,
				       epsilon,
				       g,
				       h0,
				       (long*) neighbours -> data,
				       (long*) neighbour_vertices -> data,
				       (double*) normals -> data,
				       (double*) areas -> data,
				       (double*) area_vertex_values -> data,
				       (double*) discharge_vertex_values -> data,
				       (double*) bed_vertex_values -> data,
				       (double*) height_vertex_values -> data,
				       (double*) velocity_vertex_values -> data,
				       (double*) width_vertex_values -> data,
				       (double*) area_boundary_values -> data,
				       (double*) discharge_boundary_values -> data,
				       (double*) bed_boundary_values -> data,
				       (double*) height_boundary_values -> data,
				       (double*) velocity_boundary_values -> data,
				       (double*) width_boundary_values -> data,
				       (double*) area_explicit_update -> data,
				       (double*) discharge_explicit_update -> data,
				       number_of_elements,
				       (double*) max_speed_array -> data);


    Py_DECREF(neighbours);
    Py_DECREF(neighbour_vertices);
    Py_DECREF(normals);
    Py_DECREF(areas);
    Py_DECREF(area_vertex_values);
    Py_DECREF(discharge_vertex_values);
    Py_DECREF(bed_vertex_values);
    Py_DECREF(height_vertex_values);
    Py_DECREF(velocity_vertex_values);
    Py_DECREF(width_vertex_values);
    Py_DECREF(area_boundary_values);
    Py_DECREF(discharge_boundary_values);
    Py_DECREF(bed_boundary_values);
    Py_DECREF(height_boundary_values);
    Py_DECREF(velocity_boundary_values);
    Py_DECREF(width_boundary_values);
    Py_DECREF(area_explicit_update);
    Py_DECREF(discharge_explicit_update);
    Py_DECREF(max_speed_array);


    // Return updated flux timestep
    return Py_BuildValue("d", timestep);
}


//-------------------------------
// Method table for python module
//-------------------------------

static struct PyMethodDef MethodTable[] = {
  {"compute_fluxes_channel_ext", compute_fluxes_channel_ext, METH_VARARGS, "Print out"},
  {NULL}
};

/* // Module initialisation */
/* void initcomp_flux_vel_ext(void){ */
/*   Py_InitModule("comp_flux_vel_ext", MethodTable); */
/*   import_array(); */
/* } */

void initchannel_domain_ext(void){
  Py_InitModule("channel_domain_ext", MethodTable);
  import_array();
}

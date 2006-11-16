// Python - C extension module for shallow_water.py
//
// To compile (Python2.3):
//  gcc -c domain_ext.c -I/usr/include/python2.3 -o domain_ext.o -Wall -O
//  gcc -shared domain_ext.o  -o domain_ext.so
//
// or use python compile.py
//
// See the module shallow_water.py
//
//
// Ole Nielsen, GA 2004


#include "Python.h"
#include "Numeric/arrayobject.h"
#include "math.h"
#include <stdio.h>

//Shared code snippets
#include "util_ext.h"

const double pi = 3.14159265358979;

// Computational function for rotation
int _rotate(double *q, double n1, double n2) {
  /*Rotate the momentum component q (q[1], q[2])
    from x,y coordinates to coordinates based on normal vector (n1, n2).

    Result is returned in array 3x1 r
    To rotate in opposite direction, call rotate with (q, n1, -n2)

    Contents of q are changed by this function */


  double q1, q2;

  //Shorthands
  q1 = q[1];  //uh momentum
  q2 = q[2];  //vh momentum

  //Rotate
  q[1] =  n1*q1 + n2*q2;
  q[2] = -n2*q1 + n1*q2;

  return 0;
}

int find_qmin_and_qmax(double dq0, double dq1, double dq2, double *qmin, double *qmax){
  //Considering the centroid of an FV triangle and the vertices of its auxiliary triangle, find
  //qmin=min(q)-qc and qmax=max(q)-qc, where min(q) and max(q) are respectively min and max over the
  //four values (at the centroid of the FV triangle and the auxiliary triangle vertices),
  //and qc is the centroid
  //dq0=q(vertex0)-q(centroid of FV triangle)
  //dq1=q(vertex1)-q(vertex0)
  //dq2=q(vertex2)-q(vertex0)
  if (dq0>=0.0){
    if (dq1>=dq2){
      if (dq1>=0.0)
	*qmax=dq0+dq1;
      else
	*qmax=dq0;
      if ((*qmin=dq0+dq2)<0)
	;//qmin is already set to correct value
      else
	*qmin=0.0;
    }
    else{//dq1<dq2
      if (dq2>0)
	*qmax=dq0+dq2;
      else
	*qmax=dq0;
      if ((*qmin=dq0+dq1)<0)
	;//qmin is the correct value
      else
	*qmin=0.0;
    }
  }
  else{//dq0<0
    if (dq1<=dq2){
      if (dq1<0.0)
	*qmin=dq0+dq1;
      else
	*qmin=dq0;
      if ((*qmax=dq0+dq2)>0.0)
	;//qmax is already set to the correct value
      else
	*qmax=0.0;
    }
    else{//dq1>dq2
      if (dq2<0.0)
	*qmin=dq0+dq2;
      else
	*qmin=dq0;
      if ((*qmax=dq0+dq1)>0.0)
	;//qmax is already set to the correct value
      else
	*qmax=0.0;
    }
  }
  return 0;
}

int limit_gradient(double *dqv, double qmin, double qmax, double beta_w){
  //given provisional jumps dqv from the FV triangle centroid to its vertices and
  //jumps qmin (qmax) between the centroid of the FV triangle and the
  //minimum (maximum) of the values at the centroid of the FV triangle and the auxiliary triangle vertices,
  //calculate a multiplicative factor phi by which the provisional vertex jumps are to be limited
  int i;
  double r=1000.0, r0=1.0, phi=1.0;
  static double TINY = 1.0e-100;//to avoid machine accuracy problems.
  //Any provisional jump with magnitude < TINY does not contribute to the limiting process.
  for (i=0;i<3;i++){
    if (dqv[i]<-TINY)
      r0=qmin/dqv[i];
    if (dqv[i]>TINY)
      r0=qmax/dqv[i];
    r=min(r0,r);
    //
  }
  phi=min(r*beta_w,1.0);
  for (i=0;i<3;i++)
    dqv[i]=dqv[i]*phi;
  return 0;
}

// Computational function for flux computation (using stage w=z+h)
int flux_function_central(double *q_left, double *q_right,
		  double z_left, double z_right,
		  double n1, double n2,
		  double epsilon, double g,
		  double *edgeflux, double *max_speed) {

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
  double denom, z;
  double q_left_copy[3], q_right_copy[3];
  double flux_right[3], flux_left[3];

  //Copy conserved quantities to protect from modification
  for (i=0; i<3; i++) {
    q_left_copy[i] = q_left[i];
    q_right_copy[i] = q_right[i];
  }

  //Align x- and y-momentum with x-axis
  _rotate(q_left_copy, n1, n2);
  _rotate(q_right_copy, n1, n2);

  z = (z_left+z_right)/2; //Take average of field values

  //Compute speeds in x-direction
  w_left = q_left_copy[0];              // h+z
  h_left = w_left-z;
  uh_left = q_left_copy[1];

  if (h_left < epsilon) {
    h_left = 0.0;  //Could have been negative
    u_left = 0.0;
  } else {
    u_left = uh_left/h_left;
  }

  w_right = q_right_copy[0];
  h_right = w_right-z;
  uh_right = q_right_copy[1];

  if (h_right < epsilon) {
    h_right = 0.0; //Could have been negative
    u_right = 0.0;
  } else {
    u_right = uh_right/h_right;
  }

  //Momentum in y-direction
  vh_left  = q_left_copy[2];
  vh_right = q_right_copy[2];


  //Maximal and minimal wave speeds
  soundspeed_left  = sqrt(g*h_left);
  soundspeed_right = sqrt(g*h_right);

  s_max = max(u_left+soundspeed_left, u_right+soundspeed_right);
  if (s_max < 0.0) s_max = 0.0;

  s_min = min(u_left-soundspeed_left, u_right-soundspeed_right);
  if (s_min > 0.0) s_min = 0.0;

  //Flux formulas
  flux_left[0] = u_left*h_left;
  flux_left[1] = u_left*uh_left + 0.5*g*h_left*h_left;
  flux_left[2] = u_left*vh_left;

  flux_right[0] = u_right*h_right;
  flux_right[1] = u_right*uh_right + 0.5*g*h_right*h_right;
  flux_right[2] = u_right*vh_right;


  //Flux computation
  denom = s_max-s_min;
  if (denom == 0.0) {
    for (i=0; i<3; i++) edgeflux[i] = 0.0;
    *max_speed = 0.0;
  } else {
    for (i=0; i<3; i++) {
      edgeflux[i] = s_max*flux_left[i] - s_min*flux_right[i];
      edgeflux[i] += s_max*s_min*(q_right_copy[i]-q_left_copy[i]);
      edgeflux[i] /= denom;
    }

    //Maximal wavespeed
    *max_speed = max(fabs(s_max), fabs(s_min));

    //Rotate back
    _rotate(edgeflux, n1, -n2);
  }
  return 0;
}

double erfcc(double x){
    double z,t,result;

    z=fabs(x);
    t=1.0/(1.0+0.5*z);
    result=t*exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+
         t*(.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+
         t*(1.48851587+t*(-.82215223+t*.17087277)))))))));
    if (x < 0.0) result = 2.0-result;

    return result;
    }



// Computational function for flux computation (using stage w=z+h)
int flux_function_kinetic(double *q_left, double *q_right,
		  double z_left, double z_right,
		  double n1, double n2,
		  double epsilon, double g,
		  double *edgeflux, double *max_speed) {

  /*Compute fluxes between volumes for the shallow water wave equation
    cast in terms of the 'stage', w = h+z using
    the 'central scheme' as described in

    Zhang et. al., Advances in Water Resources, 26(6), 2003, 635-647.
  */

  int i;

  double w_left, h_left, uh_left, vh_left, u_left, F_left;
  double w_right, h_right, uh_right, vh_right, u_right, F_right;
  double s_min, s_max, soundspeed_left, soundspeed_right;
  double z;
  double q_left_copy[3], q_right_copy[3];


  //Copy conserved quantities to protect from modification
  for (i=0; i<3; i++) {
    q_left_copy[i] = q_left[i];
    q_right_copy[i] = q_right[i];
  }

  //Align x- and y-momentum with x-axis
  _rotate(q_left_copy, n1, n2);
  _rotate(q_right_copy, n1, n2);

  z = (z_left+z_right)/2; //Take average of field values

  //Compute speeds in x-direction
  w_left = q_left_copy[0];              // h+z
  h_left = w_left-z;
  uh_left = q_left_copy[1];

  if (h_left < epsilon) {
    h_left = 0.0;  //Could have been negative
    u_left = 0.0;
  } else {
    u_left = uh_left/h_left;
  }

  w_right = q_right_copy[0];
  h_right = w_right-z;
  uh_right = q_right_copy[1];

  if (h_right < epsilon) {
    h_right = 0.0; //Could have been negative
    u_right = 0.0;
  } else {
    u_right = uh_right/h_right;
  }

  //Momentum in y-direction
  vh_left  = q_left_copy[2];
  vh_right = q_right_copy[2];


  //Maximal and minimal wave speeds
  soundspeed_left  = sqrt(g*h_left);
  soundspeed_right = sqrt(g*h_right);

  s_max = max(u_left+soundspeed_left, u_right+soundspeed_right);
  if (s_max < 0.0) s_max = 0.0;

  s_min = min(u_left-soundspeed_left, u_right-soundspeed_right);
  if (s_min > 0.0) s_min = 0.0;


  F_left  = 0.0;
  F_right = 0.0;
  if (h_left > 0.0) F_left = u_left/sqrt(g*h_left);
  if (h_right > 0.0) F_right = u_right/sqrt(g*h_right);

  for (i=0; i<3; i++) edgeflux[i] = 0.0;
  *max_speed = 0.0;

  edgeflux[0] = h_left*u_left/2.0*erfcc(-F_left) +  \
          h_left*sqrt(g*h_left)/2.0/sqrt(pi)*exp(-(F_left*F_left)) + \
          h_right*u_right/2.0*erfcc(F_right) -  \
          h_right*sqrt(g*h_right)/2.0/sqrt(pi)*exp(-(F_right*F_right));

  edgeflux[1] = (h_left*u_left*u_left + g/2.0*h_left*h_left)/2.0*erfcc(-F_left) + \
          u_left*h_left*sqrt(g*h_left)/2.0/sqrt(pi)*exp(-(F_left*F_left)) + \
          (h_right*u_right*u_right + g/2.0*h_right*h_right)/2.0*erfcc(F_right) -  \
          u_right*h_right*sqrt(g*h_right)/2.0/sqrt(pi)*exp(-(F_right*F_right));

  edgeflux[2] = vh_left*u_left/2.0*erfcc(-F_left) + \
          vh_left*sqrt(g*h_left)/2.0/sqrt(pi)*exp(-(F_left*F_left)) + \
          vh_right*u_right/2.0*erfcc(F_right) - \
          vh_right*sqrt(g*h_right)/2.0/sqrt(pi)*exp(-(F_right*F_right));

  //Maximal wavespeed
  *max_speed = max(fabs(s_max), fabs(s_min));

  //Rotate back
  _rotate(edgeflux, n1, -n2);

  return 0;
}




void _manning_friction(double g, double eps, int N,
		       double* w, double* z,
		       double* uh, double* vh,
		       double* eta, double* xmom, double* ymom) {

  int k;
  double S, h;

  for (k=0; k<N; k++) {
    if (eta[k] > eps) {
      h = w[k]-z[k];
      if (h >= eps) {
        S = -g * eta[k]*eta[k] * sqrt((uh[k]*uh[k] + vh[k]*vh[k]));
        S /= pow(h, 7.0/3);      //Expensive (on Ole's home computer)
        //S /= exp(7.0/3.0*log(h));      //seems to save about 15% over manning_friction
        //S /= h*h*(1 + h/3.0 - h*h/9.0); //FIXME: Could use a Taylor expansion


        //Update momentum
        xmom[k] += S*uh[k];
        ymom[k] += S*vh[k];
      }
    }
  }
}


/*
void _manning_friction_explicit(double g, double eps, int N,
		       double* w, double* z,
		       double* uh, double* vh,
		       double* eta, double* xmom, double* ymom) {

  int k;
  double S, h;

  for (k=0; k<N; k++) {
    if (eta[k] > eps) {
      h = w[k]-z[k];
      if (h >= eps) {
	S = -g * eta[k]*eta[k] * sqrt((uh[k]*uh[k] + vh[k]*vh[k]));
	S /= pow(h, 7.0/3);      //Expensive (on Ole's home computer)
	//S /= exp(7.0/3.0*log(h));      //seems to save about 15% over manning_friction
	//S /= h*h*(1 + h/3.0 - h*h/9.0); //FIXME: Could use a Taylor expansion


	//Update momentum
	xmom[k] += S*uh[k];
	ymom[k] += S*vh[k];
      }
    }
  }
}
*/

int _balance_deep_and_shallow(int N,
			      double* wc,
			      double* zc,
			      double* hc,
			      double* wv,
			      double* zv,
			      double* hv,
			      double* hvbar,
			      double* xmomc,
			      double* ymomc,
			      double* xmomv,
			      double* ymomv,
			      double alpha_balance) {

  int k, k3, i;
  double dz, hmin, alpha;

  //Compute linear combination between w-limited stages and
  //h-limited stages close to the bed elevation.

  for (k=0; k<N; k++) {
    // Compute maximal variation in bed elevation
    // This quantitiy is
    //     dz = max_i abs(z_i - z_c)
    // and it is independent of dimension
    // In the 1d case zc = (z0+z1)/2
    // In the 2d case zc = (z0+z1+z2)/3

    k3 = 3*k;

    //FIXME: Try with this one precomputed
    dz = 0.0;
    hmin = hv[k3];
    for (i=0; i<3; i++) {
      dz = max(dz, fabs(zv[k3+i]-zc[k]));
      hmin = min(hmin, hv[k3+i]);
    }


    //Create alpha in [0,1], where alpha==0 means using the h-limited
    //stage and alpha==1 means using the w-limited stage as
    //computed by the gradient limiter (both 1st or 2nd order)
    //
    //If hmin > dz/2 then alpha = 1 and the bed will have no effect
    //If hmin < 0 then alpha = 0 reverting to constant height above bed.


    if (dz > 0.0)
      //if (hmin<0.0)
      //	alpha = 0.0;
      //else
      //  alpha = max( min( hc[k]/dz, 1.0), 0.0 );
      alpha = max( min( alpha_balance*hmin/dz, 1.0), 0.0 );
    else
      alpha = 1.0;  //Flat bed
      
      
    //alpha = 1.0;  //Always deep FIXME: This actually looks good now

    //printf("dz = %.3f, alpha = %.8f\n", dz, alpha);

    //	Let
    //
    //	  wvi be the w-limited stage (wvi = zvi + hvi)
    //	  wvi- be the h-limited state (wvi- = zvi + hvi-)
    //
    //
    //	where i=0,1,2 denotes the vertex ids
    //
    //  Weighted balance between w-limited and h-limited stage is
    //
    //	  wvi := (1-alpha)*(zvi+hvi-) + alpha*(zvi+hvi)
    //
    //  It follows that the updated wvi is
    //    wvi := zvi + (1-alpha)*hvi- + alpha*hvi
    //
    //	 Momentum is balanced between constant and limited

    if (alpha < 1) {
      for (i=0; i<3; i++) {
         wv[k3+i] = zv[k3+i] + (1-alpha)*hvbar[k3+i] + alpha*hv[k3+i];

	//Update momentum as a linear combination of
	//xmomc and ymomc (shallow) and momentum
	//from extrapolator xmomv and ymomv (deep).
	xmomv[k3+i] = (1-alpha)*xmomc[k] + alpha*xmomv[k3+i];
	ymomv[k3+i] = (1-alpha)*ymomc[k] + alpha*ymomv[k3+i];
      }
    }
  }
  return 0;
}



int _protect(int N,
	     double minimum_allowed_height,
	     double maximum_allowed_speed,
	     double epsilon,
	     double* wc,
	     double* zc,
	     double* xmomc,
	     double* ymomc) {

  int k;
  double hc;
  double u, v, reduced_speed;

  //Protect against initesimal and negative heights
  for (k=0; k<N; k++) {
    hc = wc[k] - zc[k];

    if (hc < minimum_allowed_height) {
    	
      //Old code: Set momentum to zero and ensure h is non negative
      //xmomc[k] = 0.0;
      //ymomc[k] = 0.0;
      //if (hc <= 0.0) wc[k] = zc[k];


      //New code: Adjust momentum to guarantee speeds are physical
      //          ensure h is non negative
      //FIXME (Ole): This is only implemented in this C extension and
      //             has no Python equivalent
            
      if (hc <= 0.0) {
      	wc[k] = zc[k];
	xmomc[k] = 0.0;
	ymomc[k] = 0.0;
      } else {
        //Reduce excessive speeds derived from division by small hc
        
        u = xmomc[k]/hc;
	if (fabs(u) > maximum_allowed_speed) {
	  reduced_speed = maximum_allowed_speed * u/fabs(u);
	  //printf("Speed (u) has been reduced from %.3f to %.3f\n",
	  //	 u, reduced_speed);
	  xmomc[k] = reduced_speed * hc;
	}

        v = ymomc[k]/hc;
	if (fabs(v) > maximum_allowed_speed) {
	  reduced_speed = maximum_allowed_speed * v/fabs(v);
	  //printf("Speed (v) has been reduced from %.3f to %.3f\n",
	  //	 v, reduced_speed);
	  ymomc[k] = reduced_speed * hc;
	}
      }
    }
  }
  return 0;
}




int _assign_wind_field_values(int N,
			      double* xmom_update,
			      double* ymom_update,
			      double* s_vec,
			      double* phi_vec,
			      double cw) {

  //Assign windfield values to momentum updates

  int k;
  double S, s, phi, u, v;

  for (k=0; k<N; k++) {

    s = s_vec[k];
    phi = phi_vec[k];

    //Convert to radians
    phi = phi*pi/180;

    //Compute velocity vector (u, v)
    u = s*cos(phi);
    v = s*sin(phi);

    //Compute wind stress
    S = cw * sqrt(u*u + v*v);
    xmom_update[k] += S*u;
    ymom_update[k] += S*v;
  }
  return 0;
}



///////////////////////////////////////////////////////////////////
// Gateways to Python

PyObject *gravity(PyObject *self, PyObject *args) {
  //
  //  gravity(g, h, v, x, xmom, ymom)
  //


  PyArrayObject *h, *v, *x, *xmom, *ymom;
  int k, i, N, k3, k6;
  double g, avg_h, zx, zy;
  double x0, y0, x1, y1, x2, y2, z0, z1, z2;

  if (!PyArg_ParseTuple(args, "dOOOOO",
			&g, &h, &v, &x,
			&xmom, &ymom)) {
    PyErr_SetString(PyExc_RuntimeError, "shallow_water_ext.c: gravity could not parse input arguments");
    return NULL;
  }

  N = h -> dimensions[0];
  for (k=0; k<N; k++) {
    k3 = 3*k;  // base index
    k6 = 6*k;  // base index

    avg_h = 0.0;
    for (i=0; i<3; i++) {
      avg_h += ((double *) h -> data)[k3+i];
    }
    avg_h /= 3;


    //Compute bed slope
    x0 = ((double*) x -> data)[k6 + 0];
    y0 = ((double*) x -> data)[k6 + 1];
    x1 = ((double*) x -> data)[k6 + 2];
    y1 = ((double*) x -> data)[k6 + 3];
    x2 = ((double*) x -> data)[k6 + 4];
    y2 = ((double*) x -> data)[k6 + 5];


    z0 = ((double*) v -> data)[k3 + 0];
    z1 = ((double*) v -> data)[k3 + 1];
    z2 = ((double*) v -> data)[k3 + 2];

    _gradient(x0, y0, x1, y1, x2, y2, z0, z1, z2, &zx, &zy);

    //Update momentum
    ((double*) xmom -> data)[k] += -g*zx*avg_h;
    ((double*) ymom -> data)[k] += -g*zy*avg_h;
  }

  return Py_BuildValue("");
}


PyObject *manning_friction(PyObject *self, PyObject *args) {
  //
  // manning_friction(g, eps, h, uh, vh, eta, xmom_update, ymom_update)
  //


  PyArrayObject *w, *z, *uh, *vh, *eta, *xmom, *ymom;
  int N;
  double g, eps;

  if (!PyArg_ParseTuple(args, "ddOOOOOOO",
			&g, &eps, &w, &z, &uh, &vh, &eta,
			&xmom, &ymom)) {
    PyErr_SetString(PyExc_RuntimeError, "shallow_water_ext.c: manning_friction could not parse input arguments");
    return NULL;
  }


  N = w -> dimensions[0];
  _manning_friction(g, eps, N,
		    (double*) w -> data,
		    (double*) z -> data,
		    (double*) uh -> data,
		    (double*) vh -> data,
		    (double*) eta -> data,
		    (double*) xmom -> data,
		    (double*) ymom -> data);

  return Py_BuildValue("");
}


/*
PyObject *manning_friction_explicit(PyObject *self, PyObject *args) {
  //
  // manning_friction_explicit(g, eps, h, uh, vh, eta, xmom_update, ymom_update)
  //


  PyArrayObject *w, *z, *uh, *vh, *eta, *xmom, *ymom;
  int N;
  double g, eps;

  if (!PyArg_ParseTuple(args, "ddOOOOOOO",
			&g, &eps, &w, &z, &uh, &vh, &eta,
			&xmom, &ymom))
    return NULL;

  N = w -> dimensions[0];
  _manning_friction_explicit(g, eps, N,
		    (double*) w -> data,
		    (double*) z -> data,
		    (double*) uh -> data,
		    (double*) vh -> data,
		    (double*) eta -> data,
		    (double*) xmom -> data,
		    (double*) ymom -> data);

  return Py_BuildValue("");
}
*/

PyObject *extrapolate_second_order_sw(PyObject *self, PyObject *args) {
  /*Compute the vertex values based on a linear reconstruction on each triangle
    These values are calculated as follows:
    1) For each triangle not adjacent to a boundary, we consider the auxiliary triangle
    formed by the centroids of its three neighbours.
    2) For each conserved quantity, we integrate around the auxiliary triangle's boundary the product
    of the quantity and the outward normal vector. Dividing by the triangle area gives (a,b), the average
    of the vector (q_x,q_y) on the auxiliary triangle. We suppose that the linear reconstruction on the
    original triangle has gradient (a,b).
    3) Provisional vertex junmps dqv[0,1,2] are computed and these are then limited by calling the functions
    find_qmin_and_qmax and limit_gradient

    Python call:
    extrapolate_second_order_sw(domain.surrogate_neighbours,
                                domain.number_of_boundaries
                                domain.centroid_coordinates,
                                Stage.centroid_values
                                Xmom.centroid_values
                                Ymom.centroid_values
                                domain.vertex_coordinates,
                                Stage.vertex_values,
                                Xmom.vertex_values,
                                Ymom.vertex_values)

    Post conditions:
            The vertices of each triangle have values from a limited linear reconstruction
	    based on centroid values

  */
  PyArrayObject *surrogate_neighbours,
    *number_of_boundaries,
    *centroid_coordinates,
    *stage_centroid_values,
    *xmom_centroid_values,
    *ymom_centroid_values,
	*elevation_centroid_values,
    *vertex_coordinates,
    *stage_vertex_values,
    *xmom_vertex_values,
    *ymom_vertex_values,
	*elevation_vertex_values;
  PyObject *domain, *Tmp;
  //Local variables
  double a, b;//gradient vector, not stored but used to calculate vertex values from centroids
  int number_of_elements,k,k0,k1,k2,k3,k6,coord_index,i;
  double x,y,x0,y0,x1,y1,x2,y2,xv0,yv0,xv1,yv1,xv2,yv2;//vertices of the auxiliary triangle
  double dx1,dx2,dy1,dy2,dxv0,dxv1,dxv2,dyv0,dyv1,dyv2,dq0,dq1,dq2,area2;
  double dqv[3], qmin, qmax, hmin;
  double hc, h0, h1, h2;
  double beta_w, beta_w_dry, beta_uh, beta_uh_dry, beta_vh, beta_vh_dry, beta_tmp;
  double minimum_allowed_height;
  //provisional jumps from centroids to v'tices and safety factor re limiting
  //by which these jumps are limited
  // Convert Python arguments to C
  if (!PyArg_ParseTuple(args, "OOOOOOOOOOOOO",
			&domain,
			&surrogate_neighbours,
			&number_of_boundaries,
			&centroid_coordinates,
			&stage_centroid_values,
			&xmom_centroid_values,
			&ymom_centroid_values,
			&elevation_centroid_values,
			&vertex_coordinates,
			&stage_vertex_values,
			&xmom_vertex_values,
			&ymom_vertex_values,
			&elevation_vertex_values)) {
    PyErr_SetString(PyExc_RuntimeError, "Input arguments failed");
    return NULL;
  }

  //get the safety factor beta_w, set in the config.py file. This is used in the limiting process
  Tmp = PyObject_GetAttrString(domain, "beta_w");
  if (!Tmp) {
    PyErr_SetString(PyExc_RuntimeError, "shallow_water_ext.c: extrapolate_second_order_sw could not obtain object beta_w from domain");
    return NULL;
  }  
  beta_w = PyFloat_AsDouble(Tmp);
  Py_DECREF(Tmp);
  
  Tmp = PyObject_GetAttrString(domain, "beta_w_dry");
  if (!Tmp) {
    PyErr_SetString(PyExc_RuntimeError, "shallow_water_ext.c: extrapolate_second_order_sw could not obtain object beta_w_dry from domain");
    return NULL;
  }  
  beta_w_dry = PyFloat_AsDouble(Tmp);
  Py_DECREF(Tmp);
  
  Tmp = PyObject_GetAttrString(domain, "beta_uh");
  if (!Tmp) {
    PyErr_SetString(PyExc_RuntimeError, "shallow_water_ext.c: extrapolate_second_order_sw could not obtain object beta_uh from domain");
    return NULL;
  }  
  beta_uh = PyFloat_AsDouble(Tmp);
  Py_DECREF(Tmp);
  
  Tmp = PyObject_GetAttrString(domain, "beta_uh_dry");
  if (!Tmp) {
    PyErr_SetString(PyExc_RuntimeError, "shallow_water_ext.c: extrapolate_second_order_sw could not obtain object beta_uh_dry from domain");
    return NULL;
  }  
  beta_uh_dry = PyFloat_AsDouble(Tmp);
  Py_DECREF(Tmp); 

  Tmp = PyObject_GetAttrString(domain, "beta_vh");
  if (!Tmp) {
    PyErr_SetString(PyExc_RuntimeError, "shallow_water_ext.c: extrapolate_second_order_sw could not obtain object beta_vh from domain");
    return NULL;
  }  
  beta_vh = PyFloat_AsDouble(Tmp);
  Py_DECREF(Tmp);
  
  Tmp = PyObject_GetAttrString(domain, "beta_vh_dry");
  if (!Tmp) {
    PyErr_SetString(PyExc_RuntimeError, "shallow_water_ext.c: extrapolate_second_order_sw could not obtain object beta_vh_dry from domain");
    return NULL;
  }  
  beta_vh_dry = PyFloat_AsDouble(Tmp);
  Py_DECREF(Tmp);
  
  Tmp = PyObject_GetAttrString(domain, "minimum_allowed_height");
  if (!Tmp) {
    PyErr_SetString(PyExc_RuntimeError, "shallow_water_ext.c: extrapolate_second_order_sw could not obtain object minimum_allowed_heigt");
    return NULL;
  }  
  minimum_allowed_height = PyFloat_AsDouble(Tmp);
  Py_DECREF(Tmp);  
  
  number_of_elements = stage_centroid_values -> dimensions[0];
  for (k=0; k<number_of_elements; k++) {
    k3=k*3;
    k6=k*6;

    if (((long *) number_of_boundaries->data)[k]==3){/*no neighbours, set gradient on the triangle to zero*/
      ((double *) stage_vertex_values->data)[k3]=((double *)stage_centroid_values->data)[k];
      ((double *) stage_vertex_values->data)[k3+1]=((double *)stage_centroid_values->data)[k];
      ((double *) stage_vertex_values->data)[k3+2]=((double *)stage_centroid_values->data)[k];
      ((double *) xmom_vertex_values->data)[k3]=((double *)xmom_centroid_values->data)[k];
      ((double *) xmom_vertex_values->data)[k3+1]=((double *)xmom_centroid_values->data)[k];
      ((double *) xmom_vertex_values->data)[k3+2]=((double *)xmom_centroid_values->data)[k];
      ((double *) ymom_vertex_values->data)[k3]=((double *)ymom_centroid_values->data)[k];
      ((double *) ymom_vertex_values->data)[k3+1]=((double *)ymom_centroid_values->data)[k];
      ((double *) ymom_vertex_values->data)[k3+2]=((double *)ymom_centroid_values->data)[k];
      continue;
    }
    else {//we will need centroid coordinates and vertex coordinates of the triangle
      //get the vertex coordinates of the FV triangle
      xv0=((double *)vertex_coordinates->data)[k6]; yv0=((double *)vertex_coordinates->data)[k6+1];
      xv1=((double *)vertex_coordinates->data)[k6+2]; yv1=((double *)vertex_coordinates->data)[k6+3];
      xv2=((double *)vertex_coordinates->data)[k6+4]; yv2=((double *)vertex_coordinates->data)[k6+5];
      //get the centroid coordinates of the FV triangle
      coord_index=2*k;
      x=((double *)centroid_coordinates->data)[coord_index];
      y=((double *)centroid_coordinates->data)[coord_index+1];
      //store x- and y- differentials for the vertices of the FV triangle relative to the centroid
      dxv0=xv0-x; dxv1=xv1-x; dxv2=xv2-x;
      dyv0=yv0-y; dyv1=yv1-y; dyv2=yv2-y;
    }
    if (((long *)number_of_boundaries->data)[k]<=1){
      //if no boundaries, auxiliary triangle is formed from the centroids of the three neighbours
      //if one boundary, auxiliary triangle is formed from this centroid and its two neighbours
      k0=((long *)surrogate_neighbours->data)[k3];
      k1=((long *)surrogate_neighbours->data)[k3+1];
      k2=((long *)surrogate_neighbours->data)[k3+2];
      //get the auxiliary triangle's vertex coordinates (really the centroids of neighbouring triangles)
      coord_index=2*k0;
      x0=((double *)centroid_coordinates->data)[coord_index];
      y0=((double *)centroid_coordinates->data)[coord_index+1];
      coord_index=2*k1;
      x1=((double *)centroid_coordinates->data)[coord_index];
      y1=((double *)centroid_coordinates->data)[coord_index+1];
      coord_index=2*k2;
      x2=((double *)centroid_coordinates->data)[coord_index];
      y2=((double *)centroid_coordinates->data)[coord_index+1];
      //store x- and y- differentials for the vertices of the auxiliary triangle
      dx1=x1-x0; dx2=x2-x0;
      dy1=y1-y0; dy2=y2-y0;
      //calculate 2*area of the auxiliary triangle
      area2 = dy2*dx1 - dy1*dx2;//the triangle is guaranteed to be counter-clockwise
      //If the mesh is 'weird' near the boundary, the trianlge might be flat or clockwise:
      if (area2<=0) {
	PyErr_SetString(PyExc_RuntimeError, "shallow_water_ext.c: negative triangle area encountered");
	return NULL;
      }  
      

      //### Calculate heights of neighbouring cells
      hc = ((double *)stage_centroid_values->data)[k]  - ((double *)elevation_centroid_values->data)[k];
      h0 = ((double *)stage_centroid_values->data)[k0] - ((double *)elevation_centroid_values->data)[k0];
      h1 = ((double *)stage_centroid_values->data)[k1] - ((double *)elevation_centroid_values->data)[k1];
      h2 = ((double *)stage_centroid_values->data)[k2] - ((double *)elevation_centroid_values->data)[k2];
      hmin = min(hc,min(h0,min(h1,h2)));
      
      //### stage ###
      //calculate the difference between vertex 0 of the auxiliary triangle and the FV triangle centroid
      dq0=((double *)stage_centroid_values->data)[k0]-((double *)stage_centroid_values->data)[k];
      //calculate differentials between the vertices of the auxiliary triangle
      dq1=((double *)stage_centroid_values->data)[k1]-((double *)stage_centroid_values->data)[k0];
      dq2=((double *)stage_centroid_values->data)[k2]-((double *)stage_centroid_values->data)[k0];
      //calculate the gradient of stage on the auxiliary triangle
      a = dy2*dq1 - dy1*dq2;
      a /= area2;
      b = dx1*dq2 - dx2*dq1;
      b /= area2;
      //calculate provisional jumps in stage from the centroid of the FV tri to its vertices, to be limited
      dqv[0]=a*dxv0+b*dyv0;
      dqv[1]=a*dxv1+b*dyv1;
      dqv[2]=a*dxv2+b*dyv2;
      //now we want to find min and max of the centroid and the vertices of the auxiliary triangle
      //and compute jumps from the centroid to the min and max
      find_qmin_and_qmax(dq0,dq1,dq2,&qmin,&qmax);
      // Playing with dry wet interface
      hmin = qmin;
      beta_tmp = beta_w;
      if (hmin<minimum_allowed_height)
	beta_tmp = beta_w_dry;
      limit_gradient(dqv,qmin,qmax,beta_tmp);//the gradient will be limited
      for (i=0;i<3;i++)
	((double *)stage_vertex_values->data)[k3+i]=((double *)stage_centroid_values->data)[k]+dqv[i];
      
      //### xmom ###
      //calculate the difference between vertex 0 of the auxiliary triangle and the FV triangle centroid
      dq0=((double *)xmom_centroid_values->data)[k0]-((double *)xmom_centroid_values->data)[k];
      //calculate differentials between the vertices of the auxiliary triangle
      dq1=((double *)xmom_centroid_values->data)[k1]-((double *)xmom_centroid_values->data)[k0];
      dq2=((double *)xmom_centroid_values->data)[k2]-((double *)xmom_centroid_values->data)[k0];
      //calculate the gradient of xmom on the auxiliary triangle
      a = dy2*dq1 - dy1*dq2;
      a /= area2;
      b = dx1*dq2 - dx2*dq1;
      b /= area2;
      //calculate provisional jumps in stage from the centroid of the FV tri to its vertices, to be limited
      dqv[0]=a*dxv0+b*dyv0;
      dqv[1]=a*dxv1+b*dyv1;
      dqv[2]=a*dxv2+b*dyv2;
      //now we want to find min and max of the centroid and the vertices of the auxiliary triangle
      //and compute jumps from the centroid to the min and max
      find_qmin_and_qmax(dq0,dq1,dq2,&qmin,&qmax);
      beta_tmp = beta_uh;
      if (hmin<minimum_allowed_height)
	beta_tmp = beta_uh_dry;
      limit_gradient(dqv,qmin,qmax,beta_tmp);//the gradient will be limited
      for (i=0;i<3;i++)
	((double *)xmom_vertex_values->data)[k3+i]=((double *)xmom_centroid_values->data)[k]+dqv[i];
      
      //### ymom ###
      //calculate the difference between vertex 0 of the auxiliary triangle and the FV triangle centroid
      dq0=((double *)ymom_centroid_values->data)[k0]-((double *)ymom_centroid_values->data)[k];
      //calculate differentials between the vertices of the auxiliary triangle
      dq1=((double *)ymom_centroid_values->data)[k1]-((double *)ymom_centroid_values->data)[k0];
      dq2=((double *)ymom_centroid_values->data)[k2]-((double *)ymom_centroid_values->data)[k0];
      //calculate the gradient of xmom on the auxiliary triangle
      a = dy2*dq1 - dy1*dq2;
      a /= area2;
      b = dx1*dq2 - dx2*dq1;
      b /= area2;
      //calculate provisional jumps in stage from the centroid of the FV tri to its vertices, to be limited
      dqv[0]=a*dxv0+b*dyv0;
      dqv[1]=a*dxv1+b*dyv1;
      dqv[2]=a*dxv2+b*dyv2;
      //now we want to find min and max of the centroid and the vertices of the auxiliary triangle
      //and compute jumps from the centroid to the min and max
      find_qmin_and_qmax(dq0,dq1,dq2,&qmin,&qmax);
      beta_tmp = beta_vh;
      if (hmin<minimum_allowed_height)
	beta_tmp = beta_vh_dry;
      limit_gradient(dqv,qmin,qmax,beta_tmp);//the gradient will be limited
      for (i=0;i<3;i++)
	((double *)ymom_vertex_values->data)[k3+i]=((double *)ymom_centroid_values->data)[k]+dqv[i];
    }//if (number_of_boundaries[k]<=1)
    else{//number_of_boundaries==2
      //one internal neighbour and gradient is in direction of the neighbour's centroid
      //find the only internal neighbour
      for (k2=k3;k2<k3+3;k2++){//k2 just indexes the edges of triangle k
	if (((long *)surrogate_neighbours->data)[k2]!=k)//find internal neighbour of triabngle k
	  break;
      }
      if ((k2==k3+3)) {//if we didn't find an internal neighbour
	PyErr_SetString(PyExc_RuntimeError, "shallow_water_ext.c: Internal neighbour not found");      
	return NULL;//error
      }
      
      k1=((long *)surrogate_neighbours->data)[k2];
      //the coordinates of the triangle are already (x,y). Get centroid of the neighbour (x1,y1)
      coord_index=2*k1;
      x1=((double *)centroid_coordinates->data)[coord_index];
      y1=((double *)centroid_coordinates->data)[coord_index+1];
      //compute x- and y- distances between the centroid of the FV triangle and that of its neighbour
      dx1=x1-x; dy1=y1-y;
      //set area2 as the square of the distance
      area2=dx1*dx1+dy1*dy1;
      //set dx2=(x1-x0)/((x1-x0)^2+(y1-y0)^2) and dy2=(y1-y0)/((x1-x0)^2+(y1-y0)^2) which
      //respectively correspond to the x- and y- gradients of the conserved quantities
      dx2=1.0/area2;
      dy2=dx2*dy1;
      dx2*=dx1;
      
      //## stage ###
      //compute differentials
      dq1=((double *)stage_centroid_values->data)[k1]-((double *)stage_centroid_values->data)[k];
      //calculate the gradient between the centroid of the FV triangle and that of its neighbour
      a=dq1*dx2;
      b=dq1*dy2;
      //calculate provisional vertex jumps, to be limited
      dqv[0]=a*dxv0+b*dyv0;
      dqv[1]=a*dxv1+b*dyv1;
      dqv[2]=a*dxv2+b*dyv2;
      //now limit the jumps
      if (dq1>=0.0){
	qmin=0.0;
	qmax=dq1;
      }
      else{
	qmin=dq1;
	qmax=0.0;
      }
      
      
      limit_gradient(dqv,qmin,qmax,beta_w);//the gradient will be limited
      for (i=0;i<3;i++)
	((double *)stage_vertex_values->data)[k3+i]=((double *)stage_centroid_values->data)[k]+dqv[i];
      
      //## xmom ###
      //compute differentials
      dq1=((double *)xmom_centroid_values->data)[k1]-((double *)xmom_centroid_values->data)[k];
      //calculate the gradient between the centroid of the FV triangle and that of its neighbour
      a=dq1*dx2;
      b=dq1*dy2;
      //calculate provisional vertex jumps, to be limited
      dqv[0]=a*dxv0+b*dyv0;
      dqv[1]=a*dxv1+b*dyv1;
      dqv[2]=a*dxv2+b*dyv2;
      //now limit the jumps
      if (dq1>=0.0){
	qmin=0.0;
	qmax=dq1;
      }
      else{
	qmin=dq1;
	qmax=0.0;
      }
      limit_gradient(dqv,qmin,qmax,beta_w);//the gradient will be limited
      for (i=0;i<3;i++)
	((double *)xmom_vertex_values->data)[k3+i]=((double *)xmom_centroid_values->data)[k]+dqv[i];
      
      //## ymom ###
      //compute differentials
      dq1=((double *)ymom_centroid_values->data)[k1]-((double *)ymom_centroid_values->data)[k];
      //calculate the gradient between the centroid of the FV triangle and that of its neighbour
      a=dq1*dx2;
      b=dq1*dy2;
      //calculate provisional vertex jumps, to be limited
      dqv[0]=a*dxv0+b*dyv0;
      dqv[1]=a*dxv1+b*dyv1;
      dqv[2]=a*dxv2+b*dyv2;
      //now limit the jumps
      if (dq1>=0.0){
	qmin=0.0;
	qmax=dq1;
      }
      else{
	qmin=dq1;
	qmax=0.0;
      }
      limit_gradient(dqv,qmin,qmax,beta_w);//the gradient will be limited
      for (i=0;i<3;i++)
	((double *)ymom_vertex_values->data)[k3+i]=((double *)ymom_centroid_values->data)[k]+dqv[i];
    }//else [number_of_boudaries==2]
  }//for k=0 to number_of_elements-1
  return Py_BuildValue("");
}//extrapolate_second-order_sw


PyObject *rotate(PyObject *self, PyObject *args, PyObject *kwargs) {
  //
  // r = rotate(q, normal, direction=1)
  //
  // Where q is assumed to be a Float numeric array of length 3 and
  // normal a Float numeric array of length 2.


  PyObject *Q, *Normal;
  PyArrayObject *q, *r, *normal;

  static char *argnames[] = {"q", "normal", "direction", NULL};
  int dimensions[1], i, direction=1;
  double n1, n2;

  // Convert Python arguments to C
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|i", argnames,
				   &Q, &Normal, &direction)) {
    PyErr_SetString(PyExc_RuntimeError, "shallow_water_ext.c: rotate could not parse input arguments");
    return NULL;
  }  

  //Input checks (convert sequences into numeric arrays)
  q = (PyArrayObject *)
    PyArray_ContiguousFromObject(Q, PyArray_DOUBLE, 0, 0);
  normal = (PyArrayObject *)
    PyArray_ContiguousFromObject(Normal, PyArray_DOUBLE, 0, 0);


  if (normal -> dimensions[0] != 2) {
    PyErr_SetString(PyExc_RuntimeError, "Normal vector must have 2 components");
    return NULL;
  }

  //Allocate space for return vector r (don't DECREF)
  dimensions[0] = 3;
  r = (PyArrayObject *) PyArray_FromDims(1, dimensions, PyArray_DOUBLE);

  //Copy
  for (i=0; i<3; i++) {
    ((double *) (r -> data))[i] = ((double *) (q -> data))[i];
  }

  //Get normal and direction
  n1 = ((double *) normal -> data)[0];
  n2 = ((double *) normal -> data)[1];
  if (direction == -1) n2 = -n2;

  //Rotate
  _rotate((double *) r -> data, n1, n2);

  //Release numeric arrays
  Py_DECREF(q);
  Py_DECREF(normal);

  //return result using PyArray to avoid memory leak
  return PyArray_Return(r);
}

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
                                     already_computed_flux)


    Post conditions:
      domain.explicit_update is reset to computed flux values
      domain.timestep is set to the largest step satisfying all volumes.


  */


  PyArrayObject *neighbours, *neighbour_edges,
    *normals, *edgelengths, *radii, *areas,
    *tri_full_flag,
    *stage_edge_values,
    *xmom_edge_values,
    *ymom_edge_values,
    *bed_edge_values,
    *stage_boundary_values,
    *xmom_boundary_values,
    *ymom_boundary_values,
    *stage_explicit_update,
    *xmom_explicit_update,
    *ymom_explicit_update,
    *already_computed_flux;//tracks whether the flux across an edge has already been computed


  //Local variables
  double timestep, max_speed, epsilon, g;
  double normal[2], ql[3], qr[3], zl, zr;
  double edgeflux[3]; //Work arrays for summing up fluxes

  int number_of_elements, k, i, m, n;
  int ki, nm=0, ki2; //Index shorthands
  static long call=1;


  // Convert Python arguments to C
  if (!PyArg_ParseTuple(args, "dddOOOOOOOOOOOOOOOOOO",
			&timestep,
			&epsilon,
			&g,
			&neighbours,
			&neighbour_edges,
			&normals,
			&edgelengths, &radii, &areas,
			&tri_full_flag,
			&stage_edge_values,
			&xmom_edge_values,
			&ymom_edge_values,
			&bed_edge_values,
			&stage_boundary_values,
			&xmom_boundary_values,
			&ymom_boundary_values,
			&stage_explicit_update,
			&xmom_explicit_update,
			&ymom_explicit_update,
			&already_computed_flux)) {
    PyErr_SetString(PyExc_RuntimeError, "Input arguments failed");
    return NULL;
  }
  number_of_elements = stage_edge_values -> dimensions[0];
  call++;//a static local variable to which already_computed_flux is compared
  //set explicit_update to zero for all conserved_quantities.
  //This assumes compute_fluxes called before forcing terms
  for (k=0; k<number_of_elements; k++) {
    ((double *) stage_explicit_update -> data)[k]=0.0;
    ((double *) xmom_explicit_update -> data)[k]=0.0;
    ((double *) ymom_explicit_update -> data)[k]=0.0;
  }
  //Loop through neighbours and compute edge flux for each
  for (k=0; k<number_of_elements; k++) {
    for (i=0; i<3; i++) {
      ki = k*3+i;
      if (((long *) already_computed_flux->data)[ki]==call)//we've already computed the flux across this edge
	continue;
      ql[0] = ((double *) stage_edge_values -> data)[ki];
      ql[1] = ((double *) xmom_edge_values -> data)[ki];
      ql[2] = ((double *) ymom_edge_values -> data)[ki];
      zl =    ((double *) bed_edge_values -> data)[ki];

      //Quantities at neighbour on nearest face
      n = ((long *) neighbours -> data)[ki];
      if (n < 0) {
	m = -n-1; //Convert negative flag to index
	qr[0] = ((double *) stage_boundary_values -> data)[m];
	qr[1] = ((double *) xmom_boundary_values -> data)[m];
	qr[2] = ((double *) ymom_boundary_values -> data)[m];
	zr = zl; //Extend bed elevation to boundary
      } else {
	m = ((long *) neighbour_edges -> data)[ki];
	nm = n*3+m;
	qr[0] = ((double *) stage_edge_values -> data)[nm];
	qr[1] = ((double *) xmom_edge_values -> data)[nm];
	qr[2] = ((double *) ymom_edge_values -> data)[nm];
	zr =    ((double *) bed_edge_values -> data)[nm];
      }
      // Outward pointing normal vector
      // normal = domain.normals[k, 2*i:2*i+2]
      ki2 = 2*ki; //k*6 + i*2
      normal[0] = ((double *) normals -> data)[ki2];
      normal[1] = ((double *) normals -> data)[ki2+1];
      //Edge flux computation
      flux_function_central(ql, qr, zl, zr,
		    normal[0], normal[1],
		    epsilon, g,
		    edgeflux, &max_speed);
      //update triangle k
      ((long *) already_computed_flux->data)[ki]=call;
      ((double *) stage_explicit_update -> data)[k] -= edgeflux[0]*((double *) edgelengths -> data)[ki];
      ((double *) xmom_explicit_update -> data)[k] -= edgeflux[1]*((double *) edgelengths -> data)[ki];
      ((double *) ymom_explicit_update -> data)[k] -= edgeflux[2]*((double *) edgelengths -> data)[ki];
      //update the neighbour n
      if (n>=0){
	((long *) already_computed_flux->data)[nm]=call;
	((double *) stage_explicit_update -> data)[n] += edgeflux[0]*((double *) edgelengths -> data)[nm];
	((double *) xmom_explicit_update -> data)[n] += edgeflux[1]*((double *) edgelengths -> data)[nm];
	((double *) ymom_explicit_update -> data)[n] += edgeflux[2]*((double *) edgelengths -> data)[nm];
      }
      ///for (j=0; j<3; j++) {
	///flux[j] -= edgeflux[j]*((double *) edgelengths -> data)[ki];
	///}
	//Update timestep
	//timestep = min(timestep, domain.radii[k]/max_speed)
	//FIXME: SR Add parameter for CFL condition
    if ( ((long *) tri_full_flag -> data)[k] == 1) {
	    if (max_speed > epsilon) {
	        timestep = min(timestep, ((double *) radii -> data)[k]/max_speed);
	        //maxspeed in flux_function is calculated as max(|u+a|,|u-a|)
	        if (n>=0)
	            timestep = min(timestep, ((double *) radii -> data)[n]/max_speed);
	    }
    }
    } // end for i
    //Normalise by area and store for when all conserved
    //quantities get updated
    ((double *) stage_explicit_update -> data)[k] /= ((double *) areas -> data)[k];
    ((double *) xmom_explicit_update -> data)[k] /= ((double *) areas -> data)[k];
    ((double *) ymom_explicit_update -> data)[k] /= ((double *) areas -> data)[k];
  } //end for k
  return Py_BuildValue("d", timestep);
}


PyObject *compute_fluxes_ext_kinetic(PyObject *self, PyObject *args) {
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
                                     domain.g,
                                     domain.neighbours,
                                     domain.neighbour_edges,
                                     domain.normals,
                                     domain.edgelengths,
                                     domain.radii,
                                     domain.areas,
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
                                     already_computed_flux)


    Post conditions:
      domain.explicit_update is reset to computed flux values
      domain.timestep is set to the largest step satisfying all volumes.


  */


  PyArrayObject *neighbours, *neighbour_edges,
    *normals, *edgelengths, *radii, *areas,
    *tri_full_flag,
    *stage_edge_values,
    *xmom_edge_values,
    *ymom_edge_values,
    *bed_edge_values,
    *stage_boundary_values,
    *xmom_boundary_values,
    *ymom_boundary_values,
    *stage_explicit_update,
    *xmom_explicit_update,
    *ymom_explicit_update,
    *already_computed_flux;//tracks whether the flux across an edge has already been computed


  //Local variables
  double timestep, max_speed, epsilon, g;
  double normal[2], ql[3], qr[3], zl, zr;
  double edgeflux[3]; //Work arrays for summing up fluxes

  int number_of_elements, k, i, m, n;
  int ki, nm=0, ki2; //Index shorthands
  static long call=1;


  // Convert Python arguments to C
  if (!PyArg_ParseTuple(args, "dddOOOOOOOOOOOOOOOOOO",
			&timestep,
			&epsilon,
			&g,
			&neighbours,
			&neighbour_edges,
			&normals,
			&edgelengths, &radii, &areas,
			&tri_full_flag,
			&stage_edge_values,
			&xmom_edge_values,
			&ymom_edge_values,
			&bed_edge_values,
			&stage_boundary_values,
			&xmom_boundary_values,
			&ymom_boundary_values,
			&stage_explicit_update,
			&xmom_explicit_update,
			&ymom_explicit_update,
			&already_computed_flux)) {
    PyErr_SetString(PyExc_RuntimeError, "Input arguments failed");
    return NULL;
  }
  number_of_elements = stage_edge_values -> dimensions[0];
  call++;//a static local variable to which already_computed_flux is compared
  //set explicit_update to zero for all conserved_quantities.
  //This assumes compute_fluxes called before forcing terms
  for (k=0; k<number_of_elements; k++) {
    ((double *) stage_explicit_update -> data)[k]=0.0;
    ((double *) xmom_explicit_update -> data)[k]=0.0;
    ((double *) ymom_explicit_update -> data)[k]=0.0;
  }
  //Loop through neighbours and compute edge flux for each
  for (k=0; k<number_of_elements; k++) {
    for (i=0; i<3; i++) {
      ki = k*3+i;
      if (((long *) already_computed_flux->data)[ki]==call)//we've already computed the flux across this edge
	continue;
      ql[0] = ((double *) stage_edge_values -> data)[ki];
      ql[1] = ((double *) xmom_edge_values -> data)[ki];
      ql[2] = ((double *) ymom_edge_values -> data)[ki];
      zl =    ((double *) bed_edge_values -> data)[ki];

      //Quantities at neighbour on nearest face
      n = ((long *) neighbours -> data)[ki];
      if (n < 0) {
	m = -n-1; //Convert negative flag to index
	qr[0] = ((double *) stage_boundary_values -> data)[m];
	qr[1] = ((double *) xmom_boundary_values -> data)[m];
	qr[2] = ((double *) ymom_boundary_values -> data)[m];
	zr = zl; //Extend bed elevation to boundary
      } else {
	m = ((long *) neighbour_edges -> data)[ki];
	nm = n*3+m;
	qr[0] = ((double *) stage_edge_values -> data)[nm];
	qr[1] = ((double *) xmom_edge_values -> data)[nm];
	qr[2] = ((double *) ymom_edge_values -> data)[nm];
	zr =    ((double *) bed_edge_values -> data)[nm];
      }
      // Outward pointing normal vector
      // normal = domain.normals[k, 2*i:2*i+2]
      ki2 = 2*ki; //k*6 + i*2
      normal[0] = ((double *) normals -> data)[ki2];
      normal[1] = ((double *) normals -> data)[ki2+1];
      //Edge flux computation
      flux_function_kinetic(ql, qr, zl, zr,
		    normal[0], normal[1],
		    epsilon, g,
		    edgeflux, &max_speed);
      //update triangle k
      ((long *) already_computed_flux->data)[ki]=call;
      ((double *) stage_explicit_update -> data)[k] -= edgeflux[0]*((double *) edgelengths -> data)[ki];
      ((double *) xmom_explicit_update -> data)[k] -= edgeflux[1]*((double *) edgelengths -> data)[ki];
      ((double *) ymom_explicit_update -> data)[k] -= edgeflux[2]*((double *) edgelengths -> data)[ki];
      //update the neighbour n
      if (n>=0){
	((long *) already_computed_flux->data)[nm]=call;
	((double *) stage_explicit_update -> data)[n] += edgeflux[0]*((double *) edgelengths -> data)[nm];
	((double *) xmom_explicit_update -> data)[n] += edgeflux[1]*((double *) edgelengths -> data)[nm];
	((double *) ymom_explicit_update -> data)[n] += edgeflux[2]*((double *) edgelengths -> data)[nm];
      }
      ///for (j=0; j<3; j++) {
	///flux[j] -= edgeflux[j]*((double *) edgelengths -> data)[ki];
	///}
	//Update timestep
	//timestep = min(timestep, domain.radii[k]/max_speed)
	//FIXME: SR Add parameter for CFL condition
    if ( ((long *) tri_full_flag -> data)[k] == 1) {
	    if (max_speed > epsilon) {
	        timestep = min(timestep, ((double *) radii -> data)[k]/max_speed);
	        //maxspeed in flux_function is calculated as max(|u+a|,|u-a|)
	        if (n>=0)
	            timestep = min(timestep, ((double *) radii -> data)[n]/max_speed);
	    }
    }
    } // end for i
    //Normalise by area and store for when all conserved
    //quantities get updated
    ((double *) stage_explicit_update -> data)[k] /= ((double *) areas -> data)[k];
    ((double *) xmom_explicit_update -> data)[k] /= ((double *) areas -> data)[k];
    ((double *) ymom_explicit_update -> data)[k] /= ((double *) areas -> data)[k];
  } //end for k
  return Py_BuildValue("d", timestep);
}

PyObject *protect(PyObject *self, PyObject *args) {
  //
  //    protect(minimum_allowed_height, maximum_allowed_speed, wc, zc, xmomc, ymomc)


  PyArrayObject
  *wc,            //Stage at centroids
  *zc,            //Elevation at centroids
  *xmomc,         //Momentums at centroids
  *ymomc;


  int N;
  double minimum_allowed_height, maximum_allowed_speed, epsilon;

  // Convert Python arguments to C
  if (!PyArg_ParseTuple(args, "dddOOOO",
			&minimum_allowed_height,
			&maximum_allowed_speed,
			&epsilon,
			&wc, &zc, &xmomc, &ymomc)) {
    PyErr_SetString(PyExc_RuntimeError, "shallow_water_ext.c: protect could not parse input arguments");
    return NULL;
  }  

  N = wc -> dimensions[0];

  _protect(N,
	   minimum_allowed_height,
	   maximum_allowed_speed,
	   epsilon,
	   (double*) wc -> data,
	   (double*) zc -> data,
	   (double*) xmomc -> data,
	   (double*) ymomc -> data);

  return Py_BuildValue("");
}



PyObject *balance_deep_and_shallow(PyObject *self, PyObject *args) {
  //
  //    balance_deep_and_shallow(wc, zc, hc, wv, zv, hv,
  //                             xmomc, ymomc, xmomv, ymomv)


  PyArrayObject
    *wc,            //Stage at centroids
    *zc,            //Elevation at centroids
    *hc,            //Height at centroids
    *wv,            //Stage at vertices
    *zv,            //Elevation at vertices
    *hv,            //Depths at vertices
    *hvbar,         //h-Limited depths at vertices
    *xmomc,         //Momentums at centroids and vertices
    *ymomc,
    *xmomv,
    *ymomv;
  
  PyObject *domain, *Tmp;
    
  double alpha_balance = 2.0;

  int N; //, err;

  // Convert Python arguments to C
  if (!PyArg_ParseTuple(args, "OOOOOOOOOOOO",
			&domain,
			&wc, &zc, &hc,
			&wv, &zv, &hv, &hvbar,
			&xmomc, &ymomc, &xmomv, &ymomv)) {
    PyErr_SetString(PyExc_RuntimeError, "shallow_water_ext.c: balance_deep_and_shallow could not parse input arguments");
    return NULL;
  }  
	  
  // Pull out parameters
  Tmp = PyObject_GetAttrString(domain, "alpha_balance");
  if (!Tmp) {
    PyErr_SetString(PyExc_RuntimeError, "shallow_water_ext.c: balance_deep_and_shallow could not obtain object alpha_balance from domain");
    return NULL;
  }  
  alpha_balance = PyFloat_AsDouble(Tmp);
  Py_DECREF(Tmp);


  N = wc -> dimensions[0];

  _balance_deep_and_shallow(N,
			    (double*) wc -> data,
			    (double*) zc -> data,
			    (double*) hc -> data,
			    (double*) wv -> data,
			    (double*) zv -> data,
			    (double*) hv -> data,
			    (double*) hvbar -> data,
				(double*) xmomc -> data,
			    (double*) ymomc -> data,
			    (double*) xmomv -> data,
			    (double*) ymomv -> data,
			    alpha_balance);


  return Py_BuildValue("");
}



PyObject *h_limiter(PyObject *self, PyObject *args) {

  PyObject *domain, *Tmp;
  PyArrayObject
    *hv, *hc, //Depth at vertices and centroids
    *hvbar,   //Limited depth at vertices (return values)
    *neighbours;

  int k, i, n, N, k3;
  int dimensions[2];
  double beta_h; //Safety factor (see config.py)
  double *hmin, *hmax, hn;

  // Convert Python arguments to C
  if (!PyArg_ParseTuple(args, "OOO", &domain, &hc, &hv)) {
    PyErr_SetString(PyExc_RuntimeError, "shallow_water_ext.c: h_limiter could not parse input arguments");
    return NULL;
  }  
  
  neighbours = get_consecutive_array(domain, "neighbours");

  //Get safety factor beta_h
  Tmp = PyObject_GetAttrString(domain, "beta_h");
  if (!Tmp) {
    PyErr_SetString(PyExc_RuntimeError, "shallow_water_ext.c: h_limiter could not obtain object beta_h from domain");
    return NULL;
  }  
  beta_h = PyFloat_AsDouble(Tmp);

  Py_DECREF(Tmp);

  N = hc -> dimensions[0];

  //Create hvbar
  dimensions[0] = N;
  dimensions[1] = 3;
  hvbar = (PyArrayObject *) PyArray_FromDims(2, dimensions, PyArray_DOUBLE);


  //Find min and max of this and neighbour's centroid values
  hmin = malloc(N * sizeof(double));
  hmax = malloc(N * sizeof(double));
  for (k=0; k<N; k++) {
    k3=k*3;

    hmin[k] = ((double*) hc -> data)[k];
    hmax[k] = hmin[k];

    for (i=0; i<3; i++) {
      n = ((long*) neighbours -> data)[k3+i];

      //Initialise hvbar with values from hv
      ((double*) hvbar -> data)[k3+i] = ((double*) hv -> data)[k3+i];

      if (n >= 0) {
	hn = ((double*) hc -> data)[n]; //Neighbour's centroid value

	hmin[k] = min(hmin[k], hn);
	hmax[k] = max(hmax[k], hn);
      }
    }
  }

  // Call underlying standard routine
  _limit(N, beta_h, (double*) hc -> data, (double*) hvbar -> data, hmin, hmax);

  // // //Py_DECREF(domain); //FIXME: NEcessary?
  free(hmin);
  free(hmax);

  //return result using PyArray to avoid memory leak
  return PyArray_Return(hvbar);
  //return Py_BuildValue("");
}

PyObject *h_limiter_sw(PyObject *self, PyObject *args) {
  //a faster version of h_limiter above
  PyObject *domain, *Tmp;
  PyArrayObject
    *hv, *hc, //Depth at vertices and centroids
    *hvbar,   //Limited depth at vertices (return values)
    *neighbours;

  int k, i, N, k3,k0,k1,k2;
  int dimensions[2];
  double beta_h; //Safety factor (see config.py)
  double hmin, hmax, dh[3];
// Convert Python arguments to C
  if (!PyArg_ParseTuple(args, "OOO", &domain, &hc, &hv)) {
    PyErr_SetString(PyExc_RuntimeError, "shallow_water_ext.c: h_limiter_sw could not parse input arguments");
    return NULL;
  }  
  neighbours = get_consecutive_array(domain, "neighbours");

  //Get safety factor beta_h
  Tmp = PyObject_GetAttrString(domain, "beta_h");
  if (!Tmp) {
    PyErr_SetString(PyExc_RuntimeError, "shallow_water_ext.c: h_limiter_sw could not obtain object beta_h from domain");
    return NULL;
  }  
  beta_h = PyFloat_AsDouble(Tmp);

  Py_DECREF(Tmp);

  N = hc -> dimensions[0];

  //Create hvbar
  dimensions[0] = N;
  dimensions[1] = 3;
  hvbar = (PyArrayObject *) PyArray_FromDims(2, dimensions, PyArray_DOUBLE);
  for (k=0;k<N;k++){
    k3=k*3;
    //get the ids of the neighbours
    k0 = ((long*) neighbours -> data)[k3];
    k1 = ((long*) neighbours -> data)[k3+1];
    k2 = ((long*) neighbours -> data)[k3+2];
    //set hvbar provisionally
    for (i=0;i<3;i++){
      ((double*) hvbar -> data)[k3+i] = ((double*) hv -> data)[k3+i];
      dh[i]=((double*) hvbar -> data)[k3+i]-((double*) hc -> data)[k];
    }
    hmin=((double*) hc -> data)[k];
    hmax=hmin;
    if (k0>=0){
      hmin=min(hmin,((double*) hc -> data)[k0]);
      hmax=max(hmax,((double*) hc -> data)[k0]);
    }
    if (k1>=0){
      hmin=min(hmin,((double*) hc -> data)[k1]);
      hmax=max(hmax,((double*) hc -> data)[k1]);
    }
    if (k2>=0){
      hmin=min(hmin,((double*) hc -> data)[k2]);
      hmax=max(hmax,((double*) hc -> data)[k2]);
    }
    hmin-=((double*) hc -> data)[k];
    hmax-=((double*) hc -> data)[k];
    limit_gradient(dh,hmin,hmax,beta_h);
    for (i=0;i<3;i++)
      ((double*) hvbar -> data)[k3+i] = ((double*) hc -> data)[k]+dh[i];
  }
  return PyArray_Return(hvbar);
}

PyObject *assign_windfield_values(PyObject *self, PyObject *args) {
  //
  //      assign_windfield_values(xmom_update, ymom_update,
  //                              s_vec, phi_vec, self.const)



  PyArrayObject   //(one element per triangle)
  *s_vec,         //Speeds
  *phi_vec,       //Bearings
  *xmom_update,   //Momentum updates
  *ymom_update;


  int N;
  double cw;

  // Convert Python arguments to C
  if (!PyArg_ParseTuple(args, "OOOOd",
			&xmom_update,
			&ymom_update,
			&s_vec, &phi_vec,
			&cw)) {
    PyErr_SetString(PyExc_RuntimeError, "shallow_water_ext.c: assign_windfield_values could not parse input arguments");
    return NULL;
  }  
			

  N = xmom_update -> dimensions[0];

  _assign_wind_field_values(N,
	   (double*) xmom_update -> data,
	   (double*) ymom_update -> data,
	   (double*) s_vec -> data,
	   (double*) phi_vec -> data,
	   cw);

  return Py_BuildValue("");
}




//////////////////////////////////////////
// Method table for python module
static struct PyMethodDef MethodTable[] = {
  /* The cast of the function is necessary since PyCFunction values
   * only take two PyObject* parameters, and rotate() takes
   * three.
   */

  {"rotate", (PyCFunction)rotate, METH_VARARGS | METH_KEYWORDS, "Print out"},
  {"extrapolate_second_order_sw", extrapolate_second_order_sw, METH_VARARGS, "Print out"},
  {"compute_fluxes_ext_central", compute_fluxes_ext_central, METH_VARARGS, "Print out"},
  {"compute_fluxes_ext_kinetic", compute_fluxes_ext_kinetic, METH_VARARGS, "Print out"},
  {"gravity", gravity, METH_VARARGS, "Print out"},
  {"manning_friction", manning_friction, METH_VARARGS, "Print out"},
  {"balance_deep_and_shallow", balance_deep_and_shallow,
   METH_VARARGS, "Print out"},
  {"h_limiter", h_limiter,
   METH_VARARGS, "Print out"},
  {"h_limiter_sw", h_limiter_sw,
   METH_VARARGS, "Print out"},
  {"protect", protect, METH_VARARGS | METH_KEYWORDS, "Print out"},
  {"assign_windfield_values", assign_windfield_values,
   METH_VARARGS | METH_KEYWORDS, "Print out"},
  //{"distribute_to_vertices_and_edges",
  // distribute_to_vertices_and_edges, METH_VARARGS},
  //{"update_conserved_quantities",
  // update_conserved_quantities, METH_VARARGS},
  //{"set_initialcondition",
  // set_initialcondition, METH_VARARGS},
  {NULL, NULL}
};

// Module initialisation
void initshallow_water_ext(void){
  Py_InitModule("shallow_water_ext", MethodTable);

  import_array();     //Necessary for handling of NumPY structures
}

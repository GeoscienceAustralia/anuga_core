// Python - C extension for quantity module.
//
// To compile (Python2.3):
//  gcc -c util_ext.c -I/usr/include/python2.3 -o util_ext.o -Wall -O
//  gcc -shared util_ext.o  -o util_ext.so
//
// See the module quantity.py
//
//
// Ole Nielsen, GA 2004

#include "math.h"
#include "omp.h"
#include <stdint.h>

//Shared code snippets
#include "util_ext.h"

typedef int64_t keyint;

//-------------------------------------------
// Low level routines (called from wrappers)
//------------------------------------------

int64_t _compute_gradients(keyint N,
			double* centroids,
			double* centroid_values,
			int64_t* number_of_boundaries,
			int64_t* surrogate_neighbours,
			double* a,
			double* b){

  keyint i, k, k0, k1, k2, index3;
  double x0, x1, x2, y0, y1, y2, q0, q1, q2; //, det;


  for (k=0; k<N; k++) {
    index3 = 3*k;

    if (number_of_boundaries[k] < 2) {
      // Two or three true neighbours

      // Get indices of neighbours (or self when used as surrogate)
      // k0, k1, k2 = surrogate_neighbours[k,:]

      k0 = surrogate_neighbours[index3 + 0];
      k1 = surrogate_neighbours[index3 + 1];
      k2 = surrogate_neighbours[index3 + 2];


      if (k0 == k1 || k1 == k2) return -1;

      // Get data
      q0 = centroid_values[k0];
      q1 = centroid_values[k1];
      q2 = centroid_values[k2];

      x0 = centroids[k0*2]; y0 = centroids[k0*2+1];
      x1 = centroids[k1*2]; y1 = centroids[k1*2+1];
      x2 = centroids[k2*2]; y2 = centroids[k2*2+1];

      // Gradient
      _gradient(x0, y0, x1, y1, x2, y2, q0, q1, q2, &a[k], &b[k]);

    } else if (number_of_boundaries[k] == 2) {
      // One true neighbour

      // Get index of the one neighbour
      i=0; k0 = k;
      while (i<3 && k0==k) {
	k0 = surrogate_neighbours[index3 + i];
	i++;
      }
      if (k0 == k) return -1;

      k1 = k; //self

      // Get data
      q0 = centroid_values[k0];
      q1 = centroid_values[k1];

      x0 = centroids[k0*2]; y0 = centroids[k0*2+1];
      x1 = centroids[k1*2]; y1 = centroids[k1*2+1];

      // Two point gradient
      _gradient2(x0, y0, x1, y1, q0, q1, &a[k], &b[k]);

    }
    //    else:
    //        #No true neighbours -
    //        #Fall back to first order scheme
  }
  return 0;
}


int64_t _compute_local_gradients(keyint N,
			       double* vertex_coordinates,
			       double* vertex_values,
			       double* a,
			       double* b) {

  keyint k, k2, k3, k6;
  double x0, y0, x1, y1, x2, y2, v0, v1, v2;

  for (k=0; k<N; k++) {
    k6 = 6*k;
    k3 = 3*k;
    k2 = 2*k;

    // vertex coordinates
    // x0, y0, x1, y1, x2, y2 = X[k,:]
    x0 = vertex_coordinates[k6 + 0];
    y0 = vertex_coordinates[k6 + 1];
    x1 = vertex_coordinates[k6 + 2];
    y1 = vertex_coordinates[k6 + 3];
    x2 = vertex_coordinates[k6 + 4];
    y2 = vertex_coordinates[k6 + 5];

    v0 = vertex_values[k3+0];
    v1 = vertex_values[k3+1];
    v2 = vertex_values[k3+2];

    // Gradient
    _gradient(x0, y0, x1, y1, x2, y2, v0, v1, v2, &a[k], &b[k]);


    }
    return 0;
}

int64_t _extrapolate_from_gradient(keyint N,
			       double* centroids,
			       double* centroid_values,
			       double* vertex_coordinates,
			       double* vertex_values,
			       double* edge_values,
			       double* a,
			       double* b) {

  keyint k, k2, k3, k6;
  double x, y, x0, y0, x1, y1, x2, y2;

  for (k=0; k<N; k++){
    k6 = 6*k;
    k3 = 3*k;
    k2 = 2*k;

    // Centroid coordinates
    x = centroids[k2]; y = centroids[k2+1];

    // vertex coordinates
    // x0, y0, x1, y1, x2, y2 = X[k,:]
    x0 = vertex_coordinates[k6 + 0];
    y0 = vertex_coordinates[k6 + 1];
    x1 = vertex_coordinates[k6 + 2];
    y1 = vertex_coordinates[k6 + 3];
    x2 = vertex_coordinates[k6 + 4];
    y2 = vertex_coordinates[k6 + 5];

    // Extrapolate to Vertices
    vertex_values[k3+0] = centroid_values[k] + a[k]*(x0-x) + b[k]*(y0-y);
    vertex_values[k3+1] = centroid_values[k] + a[k]*(x1-x) + b[k]*(y1-y);
    vertex_values[k3+2] = centroid_values[k] + a[k]*(x2-x) + b[k]*(y2-y);

    // Extrapolate to Edges (midpoints)
    edge_values[k3+0] = 0.5*(vertex_values[k3 + 1]+vertex_values[k3 + 2]);
    edge_values[k3+1] = 0.5*(vertex_values[k3 + 2]+vertex_values[k3 + 0]);
    edge_values[k3+2] = 0.5*(vertex_values[k3 + 0]+vertex_values[k3 + 1]);

  }
  return 0;
}


int64_t _extrapolate_and_limit_from_gradient(keyint N,double beta,
					 double* centroids,
					 int64_t*   neighbours,
					 double* centroid_values,
					 double* vertex_coordinates,
					 double* vertex_values,
					 double* edge_values,
					 double* phi,
					 double* x_gradient,
					 double* y_gradient) {

  keyint i, k, k2, k3, k6;
  double x, y, x0, y0, x1, y1, x2, y2;
  keyint n;
  double qmin, qmax, qc;
  double qn[3];
  double dq, dqa[3], r;

  for (k=0; k<N; k++){
    k6 = 6*k;
    k3 = 3*k;
    k2 = 2*k;

    // Centroid coordinates
    x = centroids[k2+0];
    y = centroids[k2+1];

    // vertex coordinates
    // x0, y0, x1, y1, x2, y2 = X[k,:]
    x0 = vertex_coordinates[k6 + 0];
    y0 = vertex_coordinates[k6 + 1];
    x1 = vertex_coordinates[k6 + 2];
    y1 = vertex_coordinates[k6 + 3];
    x2 = vertex_coordinates[k6 + 4];
    y2 = vertex_coordinates[k6 + 5];

    // Extrapolate to Vertices
    vertex_values[k3+0] = centroid_values[k] + x_gradient[k]*(x0-x) + y_gradient[k]*(y0-y);
    vertex_values[k3+1] = centroid_values[k] + x_gradient[k]*(x1-x) + y_gradient[k]*(y1-y);
    vertex_values[k3+2] = centroid_values[k] + x_gradient[k]*(x2-x) + y_gradient[k]*(y2-y);

    // Extrapolate to Edges (midpoints)
    edge_values[k3+0] = 0.5*(vertex_values[k3 + 1]+vertex_values[k3 + 2]);
    edge_values[k3+1] = 0.5*(vertex_values[k3 + 2]+vertex_values[k3 + 0]);
    edge_values[k3+2] = 0.5*(vertex_values[k3 + 0]+vertex_values[k3 + 1]);
  }



  for (k=0; k<N; k++){
    k6 = 6*k;
    k3 = 3*k;
    k2 = 2*k;


    qc = centroid_values[k];

    qmin = qc;
    qmax = qc;

    for (i=0; i<3; i++) {
      n = neighbours[k3+i];
      if (n < 0) {
	qn[i] = qc;
      } else {
	qn[i] = centroid_values[n];
      }

      qmin = min(qmin, qn[i]);
      qmax = max(qmax, qn[i]);
    }

    //qtmin = min(min(min(qn[0],qn[1]),qn[2]),qc);
    //qtmax = max(max(max(qn[0],qn[1]),qn[2]),qc);

    /* 		for (i=0; i<3; i++) { */
    /* 		    n = neighbours[k3+i]; */
    /* 		    if (n < 0) { */
    /* 			qn[i] = qc; */
    /* 			qmin[i] = qtmin; */
    /* 			qmax[i] = qtmax; */
    /* 		    }  */
    /* 		} */

    phi[k] = 1.0;

    for (i=0; i<3; i++) {
      dq = edge_values[k3+i] - qc;      //Delta between edge and centroid values
      dqa[i] = dq;                      //Save dq for use in updating vertex values

      r = 1.0;

      if (dq > 0.0) r = (qmax - qc)/dq;
      if (dq < 0.0) r = (qmin - qc)/dq;

      phi[k] = min( min(r*beta, 1.0), phi[k]);

    }



    //Update gradient, edge and vertex values using phi limiter
    x_gradient[k] = x_gradient[k]*phi[k];
    y_gradient[k] = y_gradient[k]*phi[k];

    edge_values[k3+0] = qc + phi[k]*dqa[0];
    edge_values[k3+1] = qc + phi[k]*dqa[1];
    edge_values[k3+2] = qc + phi[k]*dqa[2];


    vertex_values[k3+0] = edge_values[k3+1] + edge_values[k3+2] - edge_values[k3+0];
    vertex_values[k3+1] = edge_values[k3+2] + edge_values[k3+0] - edge_values[k3+1];
    vertex_values[k3+2] = edge_values[k3+0] + edge_values[k3+1] - edge_values[k3+2];


  }

  return 0;

}




int64_t _limit_vertices_by_all_neighbours(keyint N, double beta,
				      double* centroid_values,
				      double* vertex_values,
				      double* edge_values,
				      int64_t*   neighbours,
				      double* x_gradient,
				      double* y_gradient) {


  keyint i, k, k2, k3, k6;
  keyint n;
  double qmin, qmax, qn, qc;
  double dq, dqa[3], phi, r;

  for (k=0; k<N; k++){
    k6 = 6*k;
    k3 = 3*k;
    k2 = 2*k;

    qc = centroid_values[k];
    qmin = qc;
    qmax = qc;

    for (i=0; i<3; i++) {
      n = neighbours[k3+i];
      if (n >= 0) {
	qn = centroid_values[n]; //Neighbour's centroid value

	qmin = min(qmin, qn);
	qmax = max(qmax, qn);
      }
    }

    phi = 1.0;
    for (i=0; i<3; i++) {
      r = 1.0;

      dq = vertex_values[k3+i] - qc;    //Delta between vertex and centroid values
      dqa[i] = dq;                      //Save dq for use in updating vertex values

      if (dq > 0.0) r = (qmax - qc)/dq;
      if (dq < 0.0) r = (qmin - qc)/dq;


      phi = min( min(r*beta, 1.0), phi);
    }

    //Update gradient, vertex and edge values using phi limiter
    x_gradient[k] = x_gradient[k]*phi;
    y_gradient[k] = y_gradient[k]*phi;

    vertex_values[k3+0] = qc + phi*dqa[0];
    vertex_values[k3+1] = qc + phi*dqa[1];
    vertex_values[k3+2] = qc + phi*dqa[2];

    edge_values[k3+0] = 0.5*(vertex_values[k3+1] + vertex_values[k3+2]);
    edge_values[k3+1] = 0.5*(vertex_values[k3+2] + vertex_values[k3+0]);
    edge_values[k3+2] = 0.5*(vertex_values[k3+0] + vertex_values[k3+1]);

  }

  return 0;
}




int64_t _limit_edges_by_all_neighbours(keyint N, double beta,
				   double* centroid_values,
				   double* vertex_values,
				   double* edge_values,
				   int64_t*   neighbours,
				   double* x_gradient,
				   double* y_gradient) {

  keyint i, k, k2, k3, k6;
  keyint n;
  double qmin, qmax, qn, qc, sign;
  double dq, dqa[3], phi, r;

  for (k=0; k<N; k++){
    k6 = 6*k;
    k3 = 3*k;
    k2 = 2*k;

    qc = centroid_values[k];
    qmin = qc;
    qmax = qc;

    for (i=0; i<3; i++) {
      n = neighbours[k3+i];
      if (n >= 0) {
	qn = centroid_values[n]; //Neighbour's centroid value

	qmin = min(qmin, qn);
	qmax = max(qmax, qn);
      }
    }

    sign = 0.0;
    if (qmin > 0.0) {
      sign = 1.0;
    } else if (qmax < 0) {
      sign = -1.0;
    }

    phi = 1.0;
    for (i=0; i<3; i++) {
      dq = edge_values[k3+i] - qc;      //Delta between edge and centroid values
      dqa[i] = dq;                      //Save dq for use in updating vertex values


      // Just limit non boundary edges so that we can reconstruct a linear function
      // FIXME Problem with stability on edges
      //if (neighbours[k3+i] >= 0) {
	r = 1.0;

	if (dq > 0.0) r = (qmax - qc)/dq;
	if (dq < 0.0) r = (qmin - qc)/dq;

	phi = min( min(r*beta, 1.0), phi);
	//	}

      //
      /* if (neighbours[k3+i] < 0) { */
      /* 	r = 1.0; */

      /* 	if (dq > 0.0 && (sign == -1.0 || sign == 0.0 )) r = (0.0 - qc)/dq; */
      /* 	if (dq < 0.0 && (sign ==  1.0 || sign == 0.0 )) r = (0.0 - qc)/dq; */

      /* 	phi = min( min(r*beta, 1.0), phi); */
      /* 	} */

    }

    //Update gradient, vertex and edge values using phi limiter
    x_gradient[k] = x_gradient[k]*phi;
    y_gradient[k] = y_gradient[k]*phi;

    edge_values[k3+0] = qc + phi*dqa[0];
    edge_values[k3+1] = qc + phi*dqa[1];
    edge_values[k3+2] = qc + phi*dqa[2];

    vertex_values[k3+0] = edge_values[k3+1] + edge_values[k3+2] - edge_values[k3+0];
    vertex_values[k3+1] = edge_values[k3+2] + edge_values[k3+0] - edge_values[k3+1];
    vertex_values[k3+2] = edge_values[k3+0] + edge_values[k3+1] - edge_values[k3+2];

  }

  return 0;
}


int64_t _limit_edges_by_neighbour(keyint N, double beta,
		     double* centroid_values,
		     double* vertex_values,
		     double* edge_values,
		     int64_t*   neighbours) {

	keyint i, k, k2, k3, k6;
	keyint n;
	double qmin, qmax, qn, qc;
	double dq, dqa[3], phi, r;

	for (k=0; k<N; k++){
		k6 = 6*k;
		k3 = 3*k;
		k2 = 2*k;

		qc = centroid_values[k];
		phi = 1.0;

		for (i=0; i<3; i++) {
		    dq = edge_values[k3+i] - qc;     //Delta between edge and centroid values
		    dqa[i] = dq;                      //Save dqa for use in updating vertex values

		    n = neighbours[k3+i];
		    qn = qc;
		    if (n >= 0)  qn = centroid_values[n]; //Neighbour's centroid value

		    qmin = min(qc, qn);
		    qmax = max(qc, qn);

		    r = 1.0;

		    if (dq > 0.0) r = (qmax - qc)/dq;
		    if (dq < 0.0) r = (qmin - qc)/dq;

		    phi = min( min(r*beta, 1.0), phi);

		}


		//Update edge and vertex values using phi limiter
		edge_values[k3+0] = qc + phi*dqa[0];
		edge_values[k3+1] = qc + phi*dqa[1];
		edge_values[k3+2] = qc + phi*dqa[2];

		vertex_values[k3+0] = edge_values[k3+1] + edge_values[k3+2] - edge_values[k3+0];
		vertex_values[k3+1] = edge_values[k3+2] + edge_values[k3+0] - edge_values[k3+1];
		vertex_values[k3+2] = edge_values[k3+0] + edge_values[k3+1] - edge_values[k3+2];

	}

	return 0;
}


int64_t _limit_gradient_by_neighbour(keyint N, double beta,
		     double* centroid_values,
		     double* vertex_values,
		     double* edge_values,
		     double* x_gradient,
		     double* y_gradient,
		     int64_t*   neighbours) {

	keyint i, k, k2, k3, k6;
	keyint n;
	double qmin, qmax, qn, qc;
	double dq, dqa[3], phi, r;

	for (k=0; k<N; k++){
		k6 = 6*k;
		k3 = 3*k;
		k2 = 2*k;

		qc = centroid_values[k];
		phi = 1.0;

		for (i=0; i<3; i++) {
		    dq = edge_values[k3+i] - qc;     //Delta between edge and centroid values
		    dqa[i] = dq;                      //Save dq for use in updating vertex values

		    n = neighbours[k3+i];
		    if (n >= 0) {
			qn = centroid_values[n]; //Neighbour's centroid value

			qmin = min(qc, qn);
			qmax = max(qc, qn);

			r = 1.0;

			if (dq > 0.0) r = (qmax - qc)/dq;
			if (dq < 0.0) r = (qmin - qc)/dq;

			phi = min( min(r*beta, 1.0), phi);
		    }
		}


		//Update edge and vertex values using phi limiter
		edge_values[k3+0] = qc + phi*dqa[0];
		edge_values[k3+1] = qc + phi*dqa[1];
		edge_values[k3+2] = qc + phi*dqa[2];

		vertex_values[k3+0] = edge_values[k3+1] + edge_values[k3+2] - edge_values[k3+0];
		vertex_values[k3+1] = edge_values[k3+2] + edge_values[k3+0] - edge_values[k3+1];
		vertex_values[k3+2] = edge_values[k3+0] + edge_values[k3+1] - edge_values[k3+2];

	}

	return 0;
}

int64_t _bound_vertices_below_by_constant(keyint N, double bound,
		     double* centroid_values,
		     double* vertex_values,
		     double* edge_values,
		     double* x_gradient,
		     double* y_gradient) {

	keyint i, k, k2, k3, k6;
	double qmin, qc;
	double dq, dqa[3], phi, r;

	for (k=0; k<N; k++){
		k6 = 6*k;
		k3 = 3*k;
		k2 = 2*k;

		qc = centroid_values[k];
		qmin = bound;


		phi = 1.0;
		for (i=0; i<3; i++) {
		    r = 1.0;

		    dq = vertex_values[k3+i] - qc;    //Delta between vertex and centroid values
		    dqa[i] = dq;                      //Save dq for use in updating vertex values

		    if (dq < 0.0) r = (qmin - qc)/dq;


		    phi = min( min(r, 1.0), phi);
		}


		//Update gradient, vertex and edge values using phi limiter
		x_gradient[k] = x_gradient[k]*phi;
		y_gradient[k] = y_gradient[k]*phi;

		vertex_values[k3+0] = qc + phi*dqa[0];
		vertex_values[k3+1] = qc + phi*dqa[1];
		vertex_values[k3+2] = qc + phi*dqa[2];

		edge_values[k3+0] = 0.5*(vertex_values[k3+1] + vertex_values[k3+2]);
		edge_values[k3+1] = 0.5*(vertex_values[k3+2] + vertex_values[k3+0]);
		edge_values[k3+2] = 0.5*(vertex_values[k3+0] + vertex_values[k3+1]);

	}

	return 0;
}

int64_t _bound_vertices_below_by_quantity(keyint N,
				      double* bound_vertex_values,
				      double* centroid_values,
				      double* vertex_values,
				      double* edge_values,
				      double* x_gradient,
				      double* y_gradient) {

	keyint i, k, k2, k3, k6;
	double qmin, qc;
	double dq, dqa[3], phi, r;

	for (k=0; k<N; k++){
		k6 = 6*k;
		k3 = 3*k;
		k2 = 2*k;

		qc = centroid_values[k];

		phi = 1.0;
		for (i=0; i<3; i++) {
		    r = 1.0;

		    dq = vertex_values[k3+i] - qc;    //Delta between vertex and centroid values
		    dqa[i] = dq;                      //Save dq for use in updating vertex values

		    qmin = bound_vertex_values[k3+i];
		    if (dq < 0.0) r = (qmin - qc)/dq;


		    phi = min( min(r, 1.0), phi);
		}


		//Update gradient, vertex and edge values using phi limiter
		x_gradient[k] = x_gradient[k]*phi;
		y_gradient[k] = y_gradient[k]*phi;

		vertex_values[k3+0] = qc + phi*dqa[0];
		vertex_values[k3+1] = qc + phi*dqa[1];
		vertex_values[k3+2] = qc + phi*dqa[2];

		edge_values[k3+0] = 0.5*(vertex_values[k3+1] + vertex_values[k3+2]);
		edge_values[k3+1] = 0.5*(vertex_values[k3+2] + vertex_values[k3+0]);
		edge_values[k3+2] = 0.5*(vertex_values[k3+0] + vertex_values[k3+1]);

	}

	return 0;
}

int64_t _interpolate(keyint N,
		 double* vertex_values,
		 double* edge_values,
                 double* centroid_values) {

	keyint k, k3;
	double q0, q1, q2;


	for (k=0; k<N; k++) {
		k3 = 3*k;

		q0 = vertex_values[k3 + 0];
		q1 = vertex_values[k3 + 1];
		q2 = vertex_values[k3 + 2];

                centroid_values[k] = (q0+q1+q2)/3.0;

		edge_values[k3 + 0] = 0.5*(q1+q2);
		edge_values[k3 + 1] = 0.5*(q0+q2);
		edge_values[k3 + 2] = 0.5*(q0+q1);
	}
	return 0;
}

int64_t _interpolate_from_vertices_to_edges(keyint N,
					double* vertex_values,
					double* edge_values) {

	keyint k, k3;
	double q0, q1, q2;


	for (k=0; k<N; k++) {
		k3 = 3*k;

		q0 = vertex_values[k3 + 0];
		q1 = vertex_values[k3 + 1];
		q2 = vertex_values[k3 + 2];

		edge_values[k3 + 0] = 0.5*(q1+q2);
		edge_values[k3 + 1] = 0.5*(q0+q2);
		edge_values[k3 + 2] = 0.5*(q0+q1);
	}
	return 0;
}


int64_t _interpolate_from_edges_to_vertices(keyint N,
					double* vertex_values,
					double* edge_values) {

	keyint k, k3;
	double e0, e1, e2;


	for (k=0; k<N; k++) {
		k3 = 3*k;

		e0 = edge_values[k3 + 0];
		e1 = edge_values[k3 + 1];
		e2 = edge_values[k3 + 2];

		vertex_values[k3 + 0] = e1 + e2 - e0;
		vertex_values[k3 + 1] = e2 + e0 - e1;
		vertex_values[k3 + 2] = e0 + e1 - e2;
	}
	return 0;
}

int64_t _backup_centroid_values(keyint N,
			    double* centroid_values,
			    double* centroid_backup_values) {
    // Backup centroid values


    keyint k;

    for (k=0; k<N; k++) {
	centroid_backup_values[k] = centroid_values[k];
    }


    return 0;
}


int64_t _saxpy_centroid_values(keyint N,
			   double a,
			   double b,
			   double* centroid_values,
			   double* centroid_backup_values) {
    // Saxby centroid values


    keyint k;


    for (k=0; k<N; k++) {
	centroid_values[k] = a*centroid_values[k] + b*centroid_backup_values[k];
    }


    return 0;
}


int64_t _update(keyint N,
	    double timestep,
	    double* centroid_values,
	    double* explicit_update,
	    double* semi_implicit_update) {
	// Update centroid values based on values stored in
	// explicit_update and semi_implicit_update as well as given timestep


	keyint k;
	double denominator, x;


	// // Divide semi_implicit update by conserved quantity
	// #pragma omp parallel for private(k, x)
	// for (k=0; k<N; k++) {
	// 	x = centroid_values[k];
	// 	if (x == 0.0) {
	// 		semi_implicit_update[k] = 0.0;
	// 	} else {
	// 		semi_implicit_update[k] /= x;
	// 	}
	// }


	// // Explicit updates
	// #pragma omp parallel for private(k)
	// for (k=0; k<N; k++) {
	// 	centroid_values[k] += timestep*explicit_update[k];
	// }


	// int64_t err_return = 0;

	// // Semi implicit updates
	// #pragma omp parallel for private(k, denominator) reduction(min:err_return)
	// for (k=0; k<N; k++) {
	// 	denominator = 1.0 - timestep*semi_implicit_update[k];
	// 	if (denominator <= 0.0) {
	// 		err_return = -1;
	// 	} else {
	// 		//Update conserved_quantities from semi implicit updates
	// 		centroid_values[k] /= denominator;
	// 	}
	// }

	// if (err_return == -1)
	// {
	// 	return -1;
	// }

	// // Reset semi_implicit_update here ready for next time step
	// #pragma omp parallel for private(k)
	// for (k = 0; k < N; k++)
	// {
	// 	semi_implicit_update[k] = 0.0;
	// }

	// return 0;

	int64_t err_return = 0;

	// Divide semi_implicit update by conserved quantity
	#pragma omp parallel for private(k, x)
	for (k=0; k<N; k++) {

		x = centroid_values[k];
		if (x == 0.0) {
			semi_implicit_update[k] = 0.0;
		} else {
			semi_implicit_update[k] /= x;
		}

		centroid_values[k] += timestep*explicit_update[k];

		// Semi implicit updates
		denominator = 1.0 - timestep*semi_implicit_update[k];
		if (denominator <= 0.0) {
			err_return = -1;
		} else {
			//Update conserved_quantities from semi implicit updates
			centroid_values[k] /= denominator;
		}
		
		// Reset semi_implicit_update here ready for next time step
		semi_implicit_update[k] = 0.0;
	}

	if (err_return == -1)
	{
		return -1;
	}

	return 0;
}


int64_t _average_vertex_values(keyint N,
			   int64_t* vertex_value_indices,
			   int64_t* number_of_triangles_per_node,
			   double* vertex_values,
			   double* A) {
  // Average vertex values to obtain one value per node

  keyint i, index;
  keyint k = 0; //Track triangles touching each node
  keyint current_node = 0;
  double total = 0.0;

  for (i=0; i<N; i++) {

    // if (current_node == N) {
    //   printf("Current node exceeding number of nodes (%d)", N);
    //   return 1;
    // }

		if (number_of_triangles_per_node[current_node] == 0) {
			  // Jump over orphaned node
				total = 0.0;
				k = 0;
				current_node += 1;
			}
		else {
	    index = vertex_value_indices[i];
	    k += 1;

	    // volume_id = index / 3
	    // vertex_id = index % 3
	    // total += self.vertex_values[volume_id, vertex_id]
	    total += vertex_values[index];

	    // printf("current_node=%d, index=%d, k=%d, total=%f\n", current_node, index, k, total);
	    if (number_of_triangles_per_node[current_node] == k) {
	      A[current_node] = total/k;

	      // Move on to next node
	      total = 0.0;
	      k = 0;
	      current_node += 1;
	    }
		}
  }

  return 0;
}

int64_t _average_centroid_values(keyint N,
			   int64_t* vertex_value_indices,
			   int64_t* number_of_triangles_per_node,
			   double* centroid_values,
			   double* A) {
  // Average centroid values to obtain one value per node

  keyint i, index;
  keyint volume_id;
  keyint k = 0; //Track triangles touching each node
  keyint current_node = 0;
  double total = 0.0;

  for (i=0; i<N; i++) {

		if (number_of_triangles_per_node[current_node] == 0) {
			  // Jump over orphaned node
				total = 0.0;
				k = 0;
				current_node += 1;
			}
		else {
	    index = vertex_value_indices[i];
	    k += 1;

	    volume_id = index / 3;
	    // vertex_id = index % 3;
	    // total += self.vertex_values[volume_id, vertex_id];
	    total += centroid_values[volume_id];

	    // printf("current_node=%d, index=%d, k=%d, total=%f\n", current_node, index, k, total);
	    if (number_of_triangles_per_node[current_node] == k) {
	      A[current_node] = total/k;

	      // Move on to next node
	      total = 0.0;
	      k = 0;
	      current_node += 1;
			}
    }
  }

  return 0;
}

// Note Padarn 27/11/12:
// This function is used to set all the node values of a quantity
// from a list of vertices and values at those vertices. Called in
// quantity.py by _set_vertex_values.
// Naming is a little confusing - but sticking with convention.
int64_t _set_vertex_values_c(keyint num_verts,
                        int64_t * vertices,
                        int64_t * node_index,
                        int64_t * number_of_triangles_per_node,
                        int64_t * vertex_value_indices,
                        double * vertex_values,
                        double * A
                        ){
  keyint i,j,num_triangles,u_vert_id,vert_v_index;

  for(i=0;i<num_verts;i++){

    u_vert_id=vertices[i];
    num_triangles = number_of_triangles_per_node[u_vert_id];

    for(j=0;j<num_triangles;j++){

      vert_v_index = vertex_value_indices[node_index[u_vert_id]+j];
      vertex_values[vert_v_index]=A[i];
    }

  }

  return 0;

}

int64_t _min_and_max_centroid_values(keyint N,
                                 double * qc,
                                 double * qv,
                                 int64_t * neighbours,
                                 double * qmin,
                                 double * qmax){
  
  // Find min and max of this and neighbour's centroid values

  keyint k, i, n, k3;
  double qn;

  for (k=0; k<N; k++) {
    k3=k*3;

    qmin[k] = qc[k];
    qmax[k] = qmin[k];

    for (i=0; i<3; i++) {
      n = neighbours[k3+i];
      if (n >= 0) {
        qn = qc[n]; //Neighbour's centroid value

        qmin[k] = min(qmin[k], qn);
        qmax[k] = max(qmax[k], qn);
      }
      //qmin[k] = max(qmin[k],0.5*((double*) qc -> data)[k]);
      //qmax[k] = min(qmax[k],2.0*((double*) qc -> data)[k]);
    }
  }

  return 0;


}




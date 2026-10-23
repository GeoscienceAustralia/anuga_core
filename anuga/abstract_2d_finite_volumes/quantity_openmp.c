// Python - C extension for quantity module.
//
// Ole Nielsen, GA 2004
// Stephen Roberts, ANU 2020-2026

#include "math.h"
#include "stdint.h"

#ifdef _OPENMP
#include "omp.h"
#endif

// Shared code snippets
#include "util_ext.h"


//-------------------------------------------
// Low level routines (called from wrappers)
//------------------------------------------

anuga_int _compute_gradients(anuga_int N,
			double* centroids,
			double* centroid_values,
			anuga_int* number_of_boundaries,
			anuga_int* surrogate_neighbours,
			double* a,
			double* b) {

  anuga_int k;
  anuga_int err = 0;

  #pragma omp parallel for private(k) reduction(min:err)
  for (k=0; k<N; k++) {
    anuga_int i, k0, k1, k2, index3;
    double x0, x1, x2, y0, y1, y2, q0, q1, q2;

    index3 = 3*k;

    if (number_of_boundaries[k] < 2) {
      // Two or three true neighbours

      k0 = surrogate_neighbours[index3 + 0];
      k1 = surrogate_neighbours[index3 + 1];
      k2 = surrogate_neighbours[index3 + 2];

      if (k0 == k1 || k1 == k2) { err = -1; continue; }

      q0 = centroid_values[k0];
      q1 = centroid_values[k1];
      q2 = centroid_values[k2];

      x0 = centroids[k0*2]; y0 = centroids[k0*2+1];
      x1 = centroids[k1*2]; y1 = centroids[k1*2+1];
      x2 = centroids[k2*2]; y2 = centroids[k2*2+1];

      _gradient(x0, y0, x1, y1, x2, y2, q0, q1, q2, &a[k], &b[k]);

    } else if (number_of_boundaries[k] == 2) {
      // One true neighbour

      i=0; k0 = k;
      while (i<3 && k0==k) {
        k0 = surrogate_neighbours[index3 + i];
        i++;
      }
      if (k0 == k) { err = -1; continue; }

      k1 = k;  // self

      q0 = centroid_values[k0];
      q1 = centroid_values[k1];

      x0 = centroids[k0*2]; y0 = centroids[k0*2+1];
      x1 = centroids[k1*2]; y1 = centroids[k1*2+1];

      _gradient2(x0, y0, x1, y1, q0, q1, &a[k], &b[k]);

    }
    // else: no true neighbours — fall back to first order (a[k]=b[k]=0)
  }

  return err;
}


anuga_int _compute_local_gradients(anuga_int N,
			       double* vertex_coordinates,
			       double* vertex_values,
			       double* a,
			       double* b) {

  anuga_int k;

  #pragma omp parallel for private(k)
  for (k=0; k<N; k++) {
    anuga_int k3 = 3*k;
    anuga_int k6 = 6*k;

    double x0 = vertex_coordinates[k6 + 0];
    double y0 = vertex_coordinates[k6 + 1];
    double x1 = vertex_coordinates[k6 + 2];
    double y1 = vertex_coordinates[k6 + 3];
    double x2 = vertex_coordinates[k6 + 4];
    double y2 = vertex_coordinates[k6 + 5];

    double v0 = vertex_values[k3+0];
    double v1 = vertex_values[k3+1];
    double v2 = vertex_values[k3+2];

    _gradient(x0, y0, x1, y1, x2, y2, v0, v1, v2, &a[k], &b[k]);
  }

  return 0;
}


anuga_int _extrapolate_from_gradient(anuga_int N,
			       double* centroids,
			       double* centroid_values,
			       double* vertex_coordinates,
			       double* vertex_values,
			       double* edge_values,
			       double* a,
			       double* b) {

  anuga_int k;

  #pragma omp parallel for private(k)
  for (k=0; k<N; k++) {
    anuga_int k2 = 2*k;
    anuga_int k3 = 3*k;
    anuga_int k6 = 6*k;

    double x = centroids[k2];
    double y = centroids[k2+1];

    double x0 = vertex_coordinates[k6 + 0];
    double y0 = vertex_coordinates[k6 + 1];
    double x1 = vertex_coordinates[k6 + 2];
    double y1 = vertex_coordinates[k6 + 3];
    double x2 = vertex_coordinates[k6 + 4];
    double y2 = vertex_coordinates[k6 + 5];

    vertex_values[k3+0] = centroid_values[k] + a[k]*(x0-x) + b[k]*(y0-y);
    vertex_values[k3+1] = centroid_values[k] + a[k]*(x1-x) + b[k]*(y1-y);
    vertex_values[k3+2] = centroid_values[k] + a[k]*(x2-x) + b[k]*(y2-y);

    edge_values[k3+0] = 0.5*(vertex_values[k3+1] + vertex_values[k3+2]);
    edge_values[k3+1] = 0.5*(vertex_values[k3+2] + vertex_values[k3+0]);
    edge_values[k3+2] = 0.5*(vertex_values[k3+0] + vertex_values[k3+1]);
  }

  return 0;
}


anuga_int _extrapolate_and_limit_from_gradient(anuga_int N, double beta,
					 double* centroids,
					 anuga_int* neighbours,
					 double* centroid_values,
					 double* vertex_coordinates,
					 double* vertex_values,
					 double* edge_values,
					 double* phi,
					 double* x_gradient,
					 double* y_gradient) {

  anuga_int k;

  // Pass 1: extrapolate gradients to vertices and edges
  #pragma omp parallel for private(k)
  for (k=0; k<N; k++) {
    anuga_int k2 = 2*k;
    anuga_int k3 = 3*k;
    anuga_int k6 = 6*k;

    double x = centroids[k2+0];
    double y = centroids[k2+1];

    double x0 = vertex_coordinates[k6 + 0];
    double y0 = vertex_coordinates[k6 + 1];
    double x1 = vertex_coordinates[k6 + 2];
    double y1 = vertex_coordinates[k6 + 3];
    double x2 = vertex_coordinates[k6 + 4];
    double y2 = vertex_coordinates[k6 + 5];

    vertex_values[k3+0] = centroid_values[k] + x_gradient[k]*(x0-x) + y_gradient[k]*(y0-y);
    vertex_values[k3+1] = centroid_values[k] + x_gradient[k]*(x1-x) + y_gradient[k]*(y1-y);
    vertex_values[k3+2] = centroid_values[k] + x_gradient[k]*(x2-x) + y_gradient[k]*(y2-y);

    edge_values[k3+0] = 0.5*(vertex_values[k3+1] + vertex_values[k3+2]);
    edge_values[k3+1] = 0.5*(vertex_values[k3+2] + vertex_values[k3+0]);
    edge_values[k3+2] = 0.5*(vertex_values[k3+0] + vertex_values[k3+1]);
  }

  // Pass 2: compute phi limiter and apply to gradients, edges and vertices
  #pragma omp parallel for private(k)
  for (k=0; k<N; k++) {
    anuga_int i, n;
    anuga_int k3 = 3*k;
    double qc = centroid_values[k];
    double qmin = qc, qmax = qc;
    double qn[3], dqa[3], dq, r;

    for (i=0; i<3; i++) {
      n = neighbours[k3+i];
      qn[i] = (n < 0) ? qc : centroid_values[n];
      qmin = fmin(qmin, qn[i]);
      qmax = fmax(qmax, qn[i]);
    }

    phi[k] = 1.0;
    for (i=0; i<3; i++) {
      dq = edge_values[k3+i] - qc;
      dqa[i] = dq;
      r = 1.0;
      if (dq > 0.0) r = (qmax - qc)/dq;
      if (dq < 0.0) r = (qmin - qc)/dq;
      phi[k] = fmin(fmin(r*beta, 1.0), phi[k]);
    }

    x_gradient[k] *= phi[k];
    y_gradient[k] *= phi[k];

    edge_values[k3+0] = qc + phi[k]*dqa[0];
    edge_values[k3+1] = qc + phi[k]*dqa[1];
    edge_values[k3+2] = qc + phi[k]*dqa[2];

    vertex_values[k3+0] = edge_values[k3+1] + edge_values[k3+2] - edge_values[k3+0];
    vertex_values[k3+1] = edge_values[k3+2] + edge_values[k3+0] - edge_values[k3+1];
    vertex_values[k3+2] = edge_values[k3+0] + edge_values[k3+1] - edge_values[k3+2];
  }

  return 0;
}


anuga_int _limit_vertices_by_all_neighbours(anuga_int N, double beta,
				      double* centroid_values,
				      double* vertex_values,
				      double* edge_values,
				      anuga_int* neighbours,
				      double* x_gradient,
				      double* y_gradient) {

  anuga_int k;

  #pragma omp parallel for private(k)
  for (k=0; k<N; k++) {
    anuga_int i, n;
    anuga_int k3 = 3*k;
    double qc = centroid_values[k];
    double qmin = qc, qmax = qc;
    double qn, dq, dqa[3], phi = 1.0, r;

    for (i=0; i<3; i++) {
      n = neighbours[k3+i];
      if (n >= 0) {
        qn = centroid_values[n];
        qmin = fmin(qmin, qn);
        qmax = fmax(qmax, qn);
      }
    }

    for (i=0; i<3; i++) {
      r = 1.0;
      dq = vertex_values[k3+i] - qc;
      dqa[i] = dq;
      if (dq > 0.0) r = (qmax - qc)/dq;
      if (dq < 0.0) r = (qmin - qc)/dq;
      phi = fmin(fmin(r*beta, 1.0), phi);
    }

    x_gradient[k] *= phi;
    y_gradient[k] *= phi;

    vertex_values[k3+0] = qc + phi*dqa[0];
    vertex_values[k3+1] = qc + phi*dqa[1];
    vertex_values[k3+2] = qc + phi*dqa[2];

    edge_values[k3+0] = 0.5*(vertex_values[k3+1] + vertex_values[k3+2]);
    edge_values[k3+1] = 0.5*(vertex_values[k3+2] + vertex_values[k3+0]);
    edge_values[k3+2] = 0.5*(vertex_values[k3+0] + vertex_values[k3+1]);
  }

  return 0;
}


anuga_int _limit_edges_by_all_neighbours(anuga_int N, double beta,
				   double* centroid_values,
				   double* vertex_values,
				   double* edge_values,
				   anuga_int* neighbours,
				   double* x_gradient,
				   double* y_gradient) {

  anuga_int k;

  #pragma omp parallel for private(k)
  for (k=0; k<N; k++) {
    anuga_int i, n;
    anuga_int k3 = 3*k;
    double qc = centroid_values[k];
    double qmin = qc, qmax = qc;
    double qn, dq, dqa[3], phi = 1.0, r;

    for (i=0; i<3; i++) {
      n = neighbours[k3+i];
      if (n >= 0) {
        qn = centroid_values[n];
        qmin = fmin(qmin, qn);
        qmax = fmax(qmax, qn);
      }
    }

    for (i=0; i<3; i++) {
      dq = edge_values[k3+i] - qc;
      dqa[i] = dq;
      r = 1.0;
      if (dq > 0.0) r = (qmax - qc)/dq;
      if (dq < 0.0) r = (qmin - qc)/dq;
      phi = fmin(fmin(r*beta, 1.0), phi);
    }

    x_gradient[k] *= phi;
    y_gradient[k] *= phi;

    edge_values[k3+0] = qc + phi*dqa[0];
    edge_values[k3+1] = qc + phi*dqa[1];
    edge_values[k3+2] = qc + phi*dqa[2];

    vertex_values[k3+0] = edge_values[k3+1] + edge_values[k3+2] - edge_values[k3+0];
    vertex_values[k3+1] = edge_values[k3+2] + edge_values[k3+0] - edge_values[k3+1];
    vertex_values[k3+2] = edge_values[k3+0] + edge_values[k3+1] - edge_values[k3+2];
  }

  return 0;
}


anuga_int _limit_edges_by_neighbour(anuga_int N, double beta,
		     double* centroid_values,
		     double* vertex_values,
		     double* edge_values,
		     anuga_int* neighbours) {

  anuga_int k;

  #pragma omp parallel for private(k)
  for (k=0; k<N; k++) {
    anuga_int i, n;
    anuga_int k3 = 3*k;
    double qc = centroid_values[k];
    double qn, qmin, qmax, dq, dqa[3], phi = 1.0, r;

    for (i=0; i<3; i++) {
      dq = edge_values[k3+i] - qc;
      dqa[i] = dq;

      n = neighbours[k3+i];
      qn = (n >= 0) ? centroid_values[n] : qc;

      qmin = fmin(qc, qn);
      qmax = fmax(qc, qn);

      r = 1.0;
      if (dq > 0.0) r = (qmax - qc)/dq;
      if (dq < 0.0) r = (qmin - qc)/dq;
      phi = fmin(fmin(r*beta, 1.0), phi);
    }

    edge_values[k3+0] = qc + phi*dqa[0];
    edge_values[k3+1] = qc + phi*dqa[1];
    edge_values[k3+2] = qc + phi*dqa[2];

    vertex_values[k3+0] = edge_values[k3+1] + edge_values[k3+2] - edge_values[k3+0];
    vertex_values[k3+1] = edge_values[k3+2] + edge_values[k3+0] - edge_values[k3+1];
    vertex_values[k3+2] = edge_values[k3+0] + edge_values[k3+1] - edge_values[k3+2];
  }

  return 0;
}


anuga_int _limit_gradient_by_neighbour(anuga_int N, double beta,
		     double* centroid_values,
		     double* vertex_values,
		     double* edge_values,
		     double* x_gradient,
		     double* y_gradient,
		     anuga_int* neighbours) {

  anuga_int k;

  #pragma omp parallel for private(k)
  for (k=0; k<N; k++) {
    anuga_int i, n;
    anuga_int k3 = 3*k;
    double qc = centroid_values[k];
    double qn, qmin, qmax, dq, dqa[3], phi = 1.0, r;

    for (i=0; i<3; i++) {
      dq = edge_values[k3+i] - qc;
      dqa[i] = dq;

      n = neighbours[k3+i];
      if (n >= 0) {
        qn = centroid_values[n];
        qmin = fmin(qc, qn);
        qmax = fmax(qc, qn);
        r = 1.0;
        if (dq > 0.0) r = (qmax - qc)/dq;
        if (dq < 0.0) r = (qmin - qc)/dq;
        phi = fmin(fmin(r*beta, 1.0), phi);
      }
    }

    edge_values[k3+0] = qc + phi*dqa[0];
    edge_values[k3+1] = qc + phi*dqa[1];
    edge_values[k3+2] = qc + phi*dqa[2];

    vertex_values[k3+0] = edge_values[k3+1] + edge_values[k3+2] - edge_values[k3+0];
    vertex_values[k3+1] = edge_values[k3+2] + edge_values[k3+0] - edge_values[k3+1];
    vertex_values[k3+2] = edge_values[k3+0] + edge_values[k3+1] - edge_values[k3+2];
  }

  return 0;
}


anuga_int _bound_vertices_below_by_constant(anuga_int N, double bound,
		     double* centroid_values,
		     double* vertex_values,
		     double* edge_values,
		     double* x_gradient,
		     double* y_gradient) {

  anuga_int k;

  #pragma omp parallel for private(k)
  for (k=0; k<N; k++) {
    anuga_int i;
    anuga_int k3 = 3*k;
    double qc = centroid_values[k];
    double dq, dqa[3], phi = 1.0, r;

    for (i=0; i<3; i++) {
      r = 1.0;
      dq = vertex_values[k3+i] - qc;
      dqa[i] = dq;
      if (dq < 0.0) r = (bound - qc)/dq;
      phi = fmin(fmin(r, 1.0), phi);
    }

    x_gradient[k] *= phi;
    y_gradient[k] *= phi;

    vertex_values[k3+0] = qc + phi*dqa[0];
    vertex_values[k3+1] = qc + phi*dqa[1];
    vertex_values[k3+2] = qc + phi*dqa[2];

    edge_values[k3+0] = 0.5*(vertex_values[k3+1] + vertex_values[k3+2]);
    edge_values[k3+1] = 0.5*(vertex_values[k3+2] + vertex_values[k3+0]);
    edge_values[k3+2] = 0.5*(vertex_values[k3+0] + vertex_values[k3+1]);
  }

  return 0;
}


anuga_int _bound_vertices_below_by_quantity(anuga_int N,
				      double* bound_vertex_values,
				      double* centroid_values,
				      double* vertex_values,
				      double* edge_values,
				      double* x_gradient,
				      double* y_gradient) {

  anuga_int k;

  #pragma omp parallel for private(k)
  for (k=0; k<N; k++) {
    anuga_int i;
    anuga_int k3 = 3*k;
    double qc = centroid_values[k];
    double dq, dqa[3], phi = 1.0, r;

    for (i=0; i<3; i++) {
      r = 1.0;
      dq = vertex_values[k3+i] - qc;
      dqa[i] = dq;
      if (dq < 0.0) r = (bound_vertex_values[k3+i] - qc)/dq;
      phi = fmin(fmin(r, 1.0), phi);
    }

    x_gradient[k] *= phi;
    y_gradient[k] *= phi;

    vertex_values[k3+0] = qc + phi*dqa[0];
    vertex_values[k3+1] = qc + phi*dqa[1];
    vertex_values[k3+2] = qc + phi*dqa[2];

    edge_values[k3+0] = 0.5*(vertex_values[k3+1] + vertex_values[k3+2]);
    edge_values[k3+1] = 0.5*(vertex_values[k3+2] + vertex_values[k3+0]);
    edge_values[k3+2] = 0.5*(vertex_values[k3+0] + vertex_values[k3+1]);
  }

  return 0;
}


anuga_int _interpolate(anuga_int N,
		 double* vertex_values,
		 double* edge_values,
		 double* centroid_values) {

  anuga_int k;

  #pragma omp parallel for private(k)
  for (k=0; k<N; k++) {
    anuga_int k3 = 3*k;
    double q0 = vertex_values[k3 + 0];
    double q1 = vertex_values[k3 + 1];
    double q2 = vertex_values[k3 + 2];

    centroid_values[k] = (q0+q1+q2)/3.0;
    edge_values[k3 + 0] = 0.5*(q1+q2);
    edge_values[k3 + 1] = 0.5*(q0+q2);
    edge_values[k3 + 2] = 0.5*(q0+q1);
  }

  return 0;
}


anuga_int _interpolate_from_vertices_to_edges(anuga_int N,
					double* vertex_values,
					double* edge_values) {

  anuga_int k;

  #pragma omp parallel for private(k)
  for (k=0; k<N; k++) {
    anuga_int k3 = 3*k;
    double q0 = vertex_values[k3 + 0];
    double q1 = vertex_values[k3 + 1];
    double q2 = vertex_values[k3 + 2];

    edge_values[k3 + 0] = 0.5*(q1+q2);
    edge_values[k3 + 1] = 0.5*(q0+q2);
    edge_values[k3 + 2] = 0.5*(q0+q1);
  }

  return 0;
}


anuga_int _interpolate_from_edges_to_vertices(anuga_int N,
					double* vertex_values,
					double* edge_values) {

  anuga_int k;

  #pragma omp parallel for private(k)
  for (k=0; k<N; k++) {
    anuga_int k3 = 3*k;
    double e0 = edge_values[k3 + 0];
    double e1 = edge_values[k3 + 1];
    double e2 = edge_values[k3 + 2];

    vertex_values[k3 + 0] = e1 + e2 - e0;
    vertex_values[k3 + 1] = e2 + e0 - e1;
    vertex_values[k3 + 2] = e0 + e1 - e2;
  }

  return 0;
}


anuga_int _backup_centroid_values(anuga_int N,
			    double* centroid_values,
			    double* centroid_backup_values) {

  anuga_int k;

  #pragma omp parallel for private(k)
  for (k=0; k<N; k++) {
    centroid_backup_values[k] = centroid_values[k];
  }

  return 0;
}


anuga_int _saxpy_centroid_values(anuga_int N,
			   double a,
			   double b,
			   double* centroid_values,
			   double* centroid_backup_values) {

  anuga_int k;

  #pragma omp parallel for private(k)
  for (k=0; k<N; k++) {
    centroid_values[k] = a*centroid_values[k] + b*centroid_backup_values[k];
  }

  return 0;
}


anuga_int _update(anuga_int N,
	    double timestep,
	    double* centroid_values,
	    double* explicit_update,
	    double* semi_implicit_update) {
  // Update centroid values using explicit and semi-implicit updates.
  // All three operations fused into a single parallel loop to
  // avoid multiple passes over the arrays.

  anuga_int k;
  anuga_int err = 0;

  #pragma omp parallel for private(k) reduction(min:err)
  for (k=0; k<N; k++) {
    double x, denominator;

    // Normalise semi-implicit update by current centroid value
    x = centroid_values[k];
    semi_implicit_update[k] = (x == 0.0) ? 0.0 : semi_implicit_update[k] / x;

    // Explicit update
    centroid_values[k] += timestep * explicit_update[k];

    // Semi-implicit update
    denominator = 1.0 - timestep * semi_implicit_update[k];
    if (denominator <= 0.0) {
      err = -1;
    } else {
      centroid_values[k] /= denominator;
    }

    // Reset semi-implicit update for next timestep
    semi_implicit_update[k] = 0.0;
  }

  return err;
}


// -----------------------------------------------------------------------
// Sequential routines (not amenable to simple parallelisation)
// -----------------------------------------------------------------------

anuga_int _average_vertex_values(anuga_int N,
			   anuga_int* vertex_value_indices,
			   anuga_int* number_of_triangles_per_node,
			   double* vertex_values,
			   double* A) {
  // Average vertex values to obtain one value per node.
  // Sequential: walks a sorted index list accumulating per-node totals.

  anuga_int i, index;
  anuga_int k = 0;
  anuga_int current_node = 0;
  double total = 0.0;

  for (i=0; i<N; i++) {
    if (number_of_triangles_per_node[current_node] == 0) {
      // Jump over orphaned node
      total = 0.0;
      k = 0;
      current_node += 1;
    } else {
      index = vertex_value_indices[i];
      k += 1;
      total += vertex_values[index];

      if (number_of_triangles_per_node[current_node] == k) {
        A[current_node] = total/k;
        total = 0.0;
        k = 0;
        current_node += 1;
      }
    }
  }

  return 0;
}


anuga_int _average_centroid_values(anuga_int N,
			   anuga_int* vertex_value_indices,
			   anuga_int* number_of_triangles_per_node,
			   double* centroid_values,
			   double* A) {
  // Average centroid values to obtain one value per node.
  // Sequential: same sorted-index-list accumulation pattern.

  anuga_int i, index, volume_id;
  anuga_int k = 0;
  anuga_int current_node = 0;
  double total = 0.0;

  for (i=0; i<N; i++) {
    if (number_of_triangles_per_node[current_node] == 0) {
      total = 0.0;
      k = 0;
      current_node += 1;
    } else {
      index = vertex_value_indices[i];
      k += 1;
      volume_id = index / 3;
      total += centroid_values[volume_id];

      if (number_of_triangles_per_node[current_node] == k) {
        A[current_node] = total/k;
        total = 0.0;
        k = 0;
        current_node += 1;
      }
    }
  }

  return 0;
}


// Set all vertex values at a list of nodes from an array of values.
anuga_int _set_vertex_values_c(anuga_int num_verts,
		        anuga_int* vertices,
		        anuga_int* node_index,
		        anuga_int* number_of_triangles_per_node,
		        anuga_int* vertex_value_indices,
		        double* vertex_values,
		        double* A) {

  anuga_int i, j;

  for (i=0; i<num_verts; i++) {
    anuga_int u_vert_id = vertices[i];
    anuga_int num_triangles = number_of_triangles_per_node[u_vert_id];

    for (j=0; j<num_triangles; j++) {
      anuga_int vert_v_index = vertex_value_indices[node_index[u_vert_id]+j];
      vertex_values[vert_v_index] = A[i];
    }
  }

  return 0;
}


anuga_int _min_and_max_centroid_values(anuga_int N,
				 double* qc,
				 double* qv,
				 anuga_int* neighbours,
				 double* qmin,
				 double* qmax) {

  // Find min and max of each triangle's own centroid value and its neighbours'.

  anuga_int k;

  #pragma omp parallel for private(k)
  for (k=0; k<N; k++) {
    anuga_int i, n;
    anuga_int k3 = 3*k;
    double qn;

    qmin[k] = qc[k];
    qmax[k] = qmin[k];

    for (i=0; i<3; i++) {
      n = neighbours[k3+i];
      if (n >= 0) {
        qn = qc[n];
        qmin[k] = fmin(qmin[k], qn);
        qmax[k] = fmax(qmax[k], qn);
      }
    }
  }

  return 0;
}

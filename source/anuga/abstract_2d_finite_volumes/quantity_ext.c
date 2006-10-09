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

#include "Python.h"
#include "Numeric/arrayobject.h"
#include "math.h"

//Shared code snippets
#include "util_ext.h"
//#include "quantity_ext.h" //Obsolete



int _compute_gradients(int N,
			double* centroids,
			double* centroid_values,
			long* number_of_boundaries,
			long* surrogate_neighbours,
			double* a,
			double* b){

	int i, k, k0, k1, k2, index3;
	double x0, x1, x2, y0, y1, y2, q0, q1, q2; //, det;


	for (k=0; k<N; k++) {
		index3 = 3*k;

		if (number_of_boundaries[k] < 2) {
			//Two or three true neighbours

			//Get indices of neighbours (or self when used as surrogate)
			//k0, k1, k2 = surrogate_neighbours[k,:]

			k0 = surrogate_neighbours[index3 + 0];
			k1 = surrogate_neighbours[index3 + 1];
			k2 = surrogate_neighbours[index3 + 2];


			if (k0 == k1 || k1 == k2) return -1;

			//Get data
			q0 = centroid_values[k0];
			q1 = centroid_values[k1];
			q2 = centroid_values[k2];

			x0 = centroids[k0*2]; y0 = centroids[k0*2+1];
			x1 = centroids[k1*2]; y1 = centroids[k1*2+1];
			x2 = centroids[k2*2]; y2 = centroids[k2*2+1];

			//Gradient
			_gradient(x0, y0, x1, y1, x2, y2, q0, q1, q2, &a[k], &b[k]);

		} else if (number_of_boundaries[k] == 2) {
			//One true neighbour

			//Get index of the one neighbour
			i=0; k0 = k;
			while (i<3 && k0==k) {
				k0 = surrogate_neighbours[index3 + i];
				i++;
			}
			if (k0 == k) return -1;

			k1 = k; //self

			//Get data
			q0 = centroid_values[k0];
			q1 = centroid_values[k1];

			x0 = centroids[k0*2]; y0 = centroids[k0*2+1];
			x1 = centroids[k1*2]; y1 = centroids[k1*2+1];

			//Two point gradient
			_gradient2(x0, y0, x1, y1, q0, q1, &a[k], &b[k]);


			//Old (wrong code)
			//det = x0*y1 - x1*y0;
			//if (det != 0.0) {
			//	a[k] = (y1*q0 - y0*q1)/det;
			//	b[k] = (x0*q1 - x1*q0)/det;
			//}
		}
		//    else:
		//        #No true neighbours -
		//        #Fall back to first order scheme
		//        pass
	}
	return 0;
}


int _extrapolate(int N,
		 double* centroids,
		 double* centroid_values,
		 double* vertex_coordinates,
		 double* vertex_values,
		 double* a,
		 double* b) {

	int k, k2, k3, k6;
	double x, y, x0, y0, x1, y1, x2, y2;

	for (k=0; k<N; k++){
		k6 = 6*k;
		k3 = 3*k;
		k2 = 2*k;

		//Centroid coordinates
		x = centroids[k2]; y = centroids[k2+1];

		//vertex coordinates
		//x0, y0, x1, y1, x2, y2 = X[k,:]
		x0 = vertex_coordinates[k6 + 0];
		y0 = vertex_coordinates[k6 + 1];
		x1 = vertex_coordinates[k6 + 2];
		y1 = vertex_coordinates[k6 + 3];
		x2 = vertex_coordinates[k6 + 4];
		y2 = vertex_coordinates[k6 + 5];

		//Extrapolate
		vertex_values[k3+0] = centroid_values[k] + a[k]*(x0-x) + b[k]*(y0-y);
		vertex_values[k3+1] = centroid_values[k] + a[k]*(x1-x) + b[k]*(y1-y);
		vertex_values[k3+2] = centroid_values[k] + a[k]*(x2-x) + b[k]*(y2-y);

	}
	return 0;
}




int _interpolate(int N,
		 double* vertex_values,
		 double* edge_values) {

	int k, k3;
	double q0, q1, q2;


	for (k=0; k<N; k++) {
		k3 = 3*k;

		q0 = vertex_values[k3 + 0];
		q1 = vertex_values[k3 + 1];
		q2 = vertex_values[k3 + 2];

		//printf("%f, %f, %f\n", q0, q1, q2);
		edge_values[k3 + 0] = 0.5*(q1+q2);
		edge_values[k3 + 1] = 0.5*(q0+q2);
		edge_values[k3 + 2] = 0.5*(q0+q1);
	}
	return 0;
}

int _update(int N,
	    double timestep,
	    double* centroid_values,
	    double* explicit_update,
	    double* semi_implicit_update) {
	//Update centroid values based on values stored in
	//explicit_update and semi_implicit_update as well as given timestep


	int k;
	double denominator, x;


	//Divide semi_implicit update by conserved quantity
	for (k=0; k<N; k++) {
		x = centroid_values[k];
		if (x == 0.0) {
			semi_implicit_update[k] = 0.0;
		} else {
			semi_implicit_update[k] /= x;
		}
	}


	//Semi implicit updates
	for (k=0; k<N; k++) {
		denominator = 1.0 - timestep*semi_implicit_update[k];
		if (denominator == 0.0) {
			return -1;
		} else {
			//Update conserved_quantities from semi implicit updates
			centroid_values[k] /= denominator;
		}
	}

	/*  for (k=0; k<N; k++) {*/
	/*    centroid_values[k] = exp(timestep*semi_implicit_update[k])*centroid_values[k];*/
	/*  }*/


	//Explicit updates
	for (k=0; k<N; k++) {
		centroid_values[k] += timestep*explicit_update[k];
	}


	//MH080605 set semi_implicit_update[k] to 0.0 here, rather than in update_conserved_quantities.py
	for (k=0;k<N;k++){
		semi_implicit_update[k]=0.0;
	}

	return 0;
}


/////////////////////////////////////////////////
// Gateways to Python
PyObject *update(PyObject *self, PyObject *args) {

	PyObject *quantity;
	PyArrayObject *centroid_values, *explicit_update, *semi_implicit_update;

	double timestep;
	int N, err;


	// Convert Python arguments to C
	if (!PyArg_ParseTuple(args, "Od", &quantity, &timestep))
		return NULL;

	centroid_values = get_consecutive_array(quantity, "centroid_values");
	explicit_update = get_consecutive_array(quantity, "explicit_update");
	semi_implicit_update = get_consecutive_array(quantity, "semi_implicit_update");

	N = centroid_values -> dimensions[0];

	err = _update(N, timestep,
			(double*) centroid_values -> data,
			(double*) explicit_update -> data,
			(double*) semi_implicit_update -> data);


	if (err != 0) {
		PyErr_SetString(PyExc_RuntimeError,
			"Zero division in semi implicit update - call Stephen :)");
		return NULL;
	}

	//Release and return
	Py_DECREF(centroid_values);
	Py_DECREF(explicit_update);
	Py_DECREF(semi_implicit_update);

	return Py_BuildValue("");
}


PyObject *interpolate_from_vertices_to_edges(PyObject *self, PyObject *args) {

	PyObject *quantity;
	PyArrayObject *vertex_values, *edge_values;

	int N, err;

	// Convert Python arguments to C
	if (!PyArg_ParseTuple(args, "O", &quantity))
		return NULL;

	vertex_values = get_consecutive_array(quantity, "vertex_values");
	edge_values = get_consecutive_array(quantity, "edge_values");

	N = vertex_values -> dimensions[0];

	err = _interpolate(N,
		     (double*) vertex_values -> data,
		     (double*) edge_values -> data);

	if (err != 0) {
		PyErr_SetString(PyExc_RuntimeError, "Interpolate could not be computed");
		return NULL;
	}

	//Release and return
	Py_DECREF(vertex_values);
	Py_DECREF(edge_values);

	return Py_BuildValue("");
}


PyObject *compute_gradients(PyObject *self, PyObject *args) {

	PyObject *quantity, *domain, *R;
	PyArrayObject
		*centroids,            //Coordinates at centroids
		*centroid_values,      //Values at centroids
		*number_of_boundaries, //Number of boundaries for each triangle
		*surrogate_neighbours, //True neighbours or - if one missing - self
		*a, *b;                //Return values

	int dimensions[1], N, err;

	// Convert Python arguments to C
	if (!PyArg_ParseTuple(args, "O", &quantity))
		return NULL;

	domain = PyObject_GetAttrString(quantity, "domain");
	if (!domain)
		return NULL;

	//Get pertinent variables

	centroids = get_consecutive_array(domain, "centroid_coordinates");
	centroid_values = get_consecutive_array(quantity, "centroid_values");
	surrogate_neighbours = get_consecutive_array(domain, "surrogate_neighbours");
	number_of_boundaries = get_consecutive_array(domain, "number_of_boundaries");

	N = centroid_values -> dimensions[0];

	//Release
	Py_DECREF(domain);

	//Allocate space for return vectors a and b (don't DECREF)
	dimensions[0] = N;
	a = (PyArrayObject *) PyArray_FromDims(1, dimensions, PyArray_DOUBLE);
	b = (PyArrayObject *) PyArray_FromDims(1, dimensions, PyArray_DOUBLE);



	err = _compute_gradients(N,
			(double*) centroids -> data,
			(double*) centroid_values -> data,
			(long*) number_of_boundaries -> data,
			(long*) surrogate_neighbours -> data,
			(double*) a -> data,
			(double*) b -> data);

	if (err != 0) {
		PyErr_SetString(PyExc_RuntimeError, "Gradient could not be computed");
		return NULL;
	}

	//Release
	Py_DECREF(centroids);
	Py_DECREF(centroid_values);
	Py_DECREF(number_of_boundaries);
	Py_DECREF(surrogate_neighbours);

	//Build result, release and return
	R = Py_BuildValue("OO", PyArray_Return(a), PyArray_Return(b));
	Py_DECREF(a);
	Py_DECREF(b);
	return R;
}



PyObject *extrapolate_second_order(PyObject *self, PyObject *args) {

	PyObject *quantity, *domain;
	PyArrayObject
		*centroids,            //Coordinates at centroids
		*centroid_values,      //Values at centroids
		*vertex_coordinates,   //Coordinates at vertices
		*vertex_values,        //Values at vertices
		*number_of_boundaries, //Number of boundaries for each triangle
		*surrogate_neighbours, //True neighbours or - if one missing - self
		*a, *b;                //Gradients

	//int N, err;
	int dimensions[1], N, err;
	//double *a, *b;  //Gradients

	// Convert Python arguments to C
	if (!PyArg_ParseTuple(args, "O", &quantity))
		return NULL;

	domain = PyObject_GetAttrString(quantity, "domain");
	if (!domain)
		return NULL;

	//Get pertinent variables
	centroids = get_consecutive_array(domain, "centroid_coordinates");
	centroid_values = get_consecutive_array(quantity, "centroid_values");
	surrogate_neighbours = get_consecutive_array(domain, "surrogate_neighbours");
	number_of_boundaries = get_consecutive_array(domain, "number_of_boundaries");
	vertex_coordinates = get_consecutive_array(domain, "vertex_coordinates");
	vertex_values = get_consecutive_array(quantity, "vertex_values");

	N = centroid_values -> dimensions[0];

	//Release
	Py_DECREF(domain);

	//Allocate space for return vectors a and b (don't DECREF)
	dimensions[0] = N;
	a = (PyArrayObject *) PyArray_FromDims(1, dimensions, PyArray_DOUBLE);
	b = (PyArrayObject *) PyArray_FromDims(1, dimensions, PyArray_DOUBLE);

	//FIXME: Odd that I couldn't use normal arrays
	//Allocate space for return vectors a and b (don't DECREF)
	//a = (double*) malloc(N * sizeof(double));
	//if (!a) return NULL;
	//b = (double*) malloc(N * sizeof(double));
	//if (!b) return NULL;


	err = _compute_gradients(N,
			(double*) centroids -> data,
			(double*) centroid_values -> data,
			(long*) number_of_boundaries -> data,
			(long*) surrogate_neighbours -> data,
			(double*) a -> data,
			(double*) b -> data);

	if (err != 0) {
		PyErr_SetString(PyExc_RuntimeError, "Gradient could not be computed");
		return NULL;
	}

	err = _extrapolate(N,
			(double*) centroids -> data,
			(double*) centroid_values -> data,
			(double*) vertex_coordinates -> data,
			(double*) vertex_values -> data,
			(double*) a -> data,
			(double*) b -> data);


	if (err != 0) {
		PyErr_SetString(PyExc_RuntimeError,
			"Internal function _extrapolate failed");
		return NULL;
	}



	//Release
	Py_DECREF(centroids);
	Py_DECREF(centroid_values);
	Py_DECREF(number_of_boundaries);
	Py_DECREF(surrogate_neighbours);
	Py_DECREF(vertex_coordinates);
	Py_DECREF(vertex_values);
	Py_DECREF(a);
	Py_DECREF(b);

	return Py_BuildValue("");
}



PyObject *limit(PyObject *self, PyObject *args) {

	PyObject *quantity, *domain, *Tmp;
	PyArrayObject
		*qv, //Conserved quantities at vertices
		*qc, //Conserved quantities at centroids
		*neighbours;

	int k, i, n, N, k3;
	double beta_w; //Safety factor
	double *qmin, *qmax, qn;

	// Convert Python arguments to C
	if (!PyArg_ParseTuple(args, "O", &quantity))
		return NULL;

	domain = PyObject_GetAttrString(quantity, "domain");
	if (!domain)
		return NULL;

	//neighbours = (PyArrayObject*) PyObject_GetAttrString(domain, "neighbours");
	neighbours = get_consecutive_array(domain, "neighbours");

	//Get safety factor beta_w
	Tmp = PyObject_GetAttrString(domain, "beta_w");
	if (!Tmp)
		return NULL;

	beta_w = PyFloat_AsDouble(Tmp);

	Py_DECREF(Tmp);
	Py_DECREF(domain);


	qc = get_consecutive_array(quantity, "centroid_values");
	qv = get_consecutive_array(quantity, "vertex_values");


	N = qc -> dimensions[0];

	//Find min and max of this and neighbour's centroid values
	qmin = malloc(N * sizeof(double));
	qmax = malloc(N * sizeof(double));
	for (k=0; k<N; k++) {
		k3=k*3;

		qmin[k] = ((double*) qc -> data)[k];
		qmax[k] = qmin[k];

		for (i=0; i<3; i++) {
			n = ((long*) neighbours -> data)[k3+i];
			if (n >= 0) {
				qn = ((double*) qc -> data)[n]; //Neighbour's centroid value

				qmin[k] = min(qmin[k], qn);
				qmax[k] = max(qmax[k], qn);
			}
			//qmin[k] = max(qmin[k],0.5*((double*) qc -> data)[k]);
			//qmax[k] = min(qmax[k],2.0*((double*) qc -> data)[k]);
		}
	}

	// Call underlying routine
	_limit(N, beta_w, (double*) qc -> data, (double*) qv -> data, qmin, qmax);

	free(qmin);
	free(qmax);
	return Py_BuildValue("");
}



// Method table for python module
static struct PyMethodDef MethodTable[] = {
	{"limit", limit, METH_VARARGS, "Print out"},
	{"update", update, METH_VARARGS, "Print out"},
	{"compute_gradients", compute_gradients, METH_VARARGS, "Print out"},
	{"extrapolate_second_order", extrapolate_second_order,
		METH_VARARGS, "Print out"},
	{"interpolate_from_vertices_to_edges",
		interpolate_from_vertices_to_edges,
		METH_VARARGS, "Print out"},
	{NULL, NULL, 0, NULL}   // sentinel
};

// Module initialisation
void initquantity_ext(void){
  Py_InitModule("quantity_ext", MethodTable);

  import_array();     //Necessary for handling of NumPY structures
}

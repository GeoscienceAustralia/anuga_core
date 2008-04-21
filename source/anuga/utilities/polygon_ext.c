// Python - C extension for polygon module.
//
// To compile (Python2.3):
//  gcc -c polygon_ext.c -I/usr/include/python2.3 -o polygon_ext.o -Wall -O
//  gcc -shared polygon_ext.o  -o polygon_ext.so
//
// See the module polygon.py
//
//
// Ole Nielsen, GA 2004
//
// NOTE: We use long* instead of int* for Numeric arrays as this will work both 
//       for 64 as well as 32 bit systems


#include "Python.h"
#include "Numeric/arrayobject.h"
#include "math.h"


int __point_on_line(double x, double y,
		    double x0, double y0,
		    double x1, double y1) {
  /*Determine whether a point is on a line segment

    Input: x, y, x0, x0, x1, y1: where
        point is given by x, y
	line is given by (x0, y0) and (x1, y1)

  */

  double a0, a1, a_normal0, a_normal1, b0, b1, len_a, len_b;

  a0 = x - x0;
  a1 = y - y0;

  a_normal0 = a1;
  a_normal1 = -a0;

  b0 = x1 - x0;
  b1 = y1 - y0;

  if ( a_normal0*b0 + a_normal1*b1 == 0.0 ) {
    //Point is somewhere on the infinite extension of the line
    // FIXME (Ole): Perhaps add a tolerance here instead of 0.0 

    len_a = sqrt(a0*a0 + a1*a1);
    len_b = sqrt(b0*b0 + b1*b1);

    if (a0*b0 + a1*b1 >= 0 && len_a <= len_b) {
      return 1;
    } else {
      return 0;
    }
  } else {
    return 0;
  }
}



int __separate_points_by_polygon(int M,     // Number of points
				int N,     // Number of polygon vertices
				double* points,
				double* polygon,
				long* indices,  // M-Array for storage indices
				int closed,
				int verbose) {

  double minpx, maxpx, minpy, maxpy, x, y, px_i, py_i, px_j, py_j;
  int i, j, k, outside_index, inside_index, inside;

  //Find min and max of poly used for optimisation when points
  //are far away from polygon

  minpx = polygon[0]; maxpx = minpx;
  minpy = polygon[1]; maxpy = minpy;

  for (i=0; i<N; i++) {
    px_i = polygon[2*i];
    py_i = polygon[2*i + 1];

    if (px_i < minpx) minpx = px_i;
    if (px_i > maxpx) maxpx = px_i;
    if (py_i < minpy) minpy = py_i;
    if (py_i > maxpy) maxpy = py_i;
  }

  //Begin main loop (for each point)
  inside_index = 0;    //Keep track of points inside
  outside_index = M-1; //Keep track of points outside (starting from end)   
  if (verbose){
     printf("Separating %d points\n", M);
  }  
  for (k=0; k<M; k++) {
    if (verbose){
      if (k %((M+10)/10)==0) printf("Doing %d of %d\n", k, M);
    }
    
    x = points[2*k];
    y = points[2*k + 1];

    inside = 0;

    //Optimisation
    if ((x > maxpx) || (x < minpx) || (y > maxpy) || (y < minpy)) {
      //Nothing
    } else {   
      //Check polygon
      for (i=0; i<N; i++) {
        //printf("k,i=%d,%d\n", k, i);
        j = (i+1)%N;

        px_i = polygon[2*i];
        py_i = polygon[2*i+1];
        px_j = polygon[2*j];
        py_j = polygon[2*j+1];

        //Check for case where point is contained in line segment
        if (__point_on_line(x, y, px_i, py_i, px_j, py_j)) {
	  if (closed == 1) {
	    inside = 1;
	  } else {
	    inside = 0;
	  }
	  break;
        } else {
          //Check if truly inside polygon
	  if ( ((py_i < y) && (py_j >= y)) ||
	       ((py_j < y) && (py_i >= y)) ) {
	    if (px_i + (y-py_i)/(py_j-py_i)*(px_j-px_i) < x)
	      inside = 1-inside;
	  }
        }
      }
    } 
    if (inside == 1) {
      indices[inside_index] = k;
      inside_index += 1;
    } else {
      indices[outside_index] = k;
      outside_index -= 1;    
    }
  } // End k

  return inside_index;
}



// Gateways to Python
PyObject *_point_on_line(PyObject *self, PyObject *args) {
  //
  // point_on_line(x, y, x0, y0, x1, y1)
  //

  double x, y, x0, y0, x1, y1;
  int res;
  PyObject *result;

  // Convert Python arguments to C
  if (!PyArg_ParseTuple(args, "dddddd", &x, &y, &x0, &y0, &x1, &y1)) {
    PyErr_SetString(PyExc_RuntimeError, 
		    "point_on_line could not parse input");    
    return NULL;
  }


  // Call underlying routine
  res = __point_on_line(x, y, x0, y0, x1, y1);

  // Return values a and b
  result = Py_BuildValue("i", res);
  return result;
}



PyObject *_separate_points_by_polygon(PyObject *self, PyObject *args) {
  //def separate_points_by_polygon(points, polygon, closed, verbose, one_point):
  //  """Determine whether points are inside or outside a polygon
  //
  //  Input:
  //     point - Tuple of (x, y) coordinates, or list of tuples
  //     polygon - list of vertices of polygon
  //     closed - (optional) determine whether points on boundary should be
  //     regarded as belonging to the polygon (closed = True)
  //     or not (closed = False)

  //
  //  Output:
  //     indices: array of same length as points with indices of points falling 
  //     inside the polygon listed from the beginning and indices of points 
  //     falling outside listed from the end.
  //     
  //     count: count of points falling inside the polygon
  //     
  //     The indices of points inside are obtained as indices[:count]
  //     The indices of points outside are obtained as indices[count:]       
  //
  //  Examples:
  //     separate_polygon( [[0.5, 0.5], [1, -0.5], [0.3, 0.2]] )
  //     will return the indices [0, 2, 1] as only the first and the last point
  //     is inside the unit square
  //
  //  Remarks:
  //     The vertices may be listed clockwise or counterclockwise and
  //     the first point may optionally be repeated.
  //     Polygons do not need to be convex.
  //     Polygons can have holes in them and points inside a hole is
  //     regarded as being outside the polygon.
  //
  //
  //  Algorithm is based on work by Darel Finley,
  //  http://www.alienryderflex.com/polygon/
  //
  //

  PyArrayObject
    *points,
    *polygon,
    *indices;

  int closed, verbose; //Flags
  int count, M, N;

  // Convert Python arguments to C
  if (!PyArg_ParseTuple(args, "OOOii",
			&points,
			&polygon,
			&indices,
			&closed,
			&verbose)) {
    

    PyErr_SetString(PyExc_RuntimeError, 
		    "separate_points_by_polygon could not parse input");
    return NULL;
  }

  M = points -> dimensions[0];   //Number of points
  N = polygon -> dimensions[0];  //Number of vertices in polygon

  //FIXME (Ole): This isn't send to Python's sys.stdout
  if (verbose) printf("Got %d points and %d polygon vertices\n", M, N);
  
  //Call underlying routine
  count = __separate_points_by_polygon(M, N,
				       (double*) points -> data,
				       (double*) polygon -> data,
				       (long*) indices -> data,
				       closed, verbose);
  
  //NOTE: return number of points inside..
  return Py_BuildValue("i", count);
}



// Method table for python module
static struct PyMethodDef MethodTable[] = {
  /* The cast of the function is necessary since PyCFunction values
   * only take two PyObject* parameters, and rotate() takes
   * three.
   */

  {"_point_on_line", _point_on_line, METH_VARARGS, "Print out"},
  {"_separate_points_by_polygon", _separate_points_by_polygon, 
                                 METH_VARARGS, "Print out"},
  {NULL, NULL, 0, NULL}   /* sentinel */
};



// Module initialisation
void initpolygon_ext(void){
  Py_InitModule("polygon_ext", MethodTable);

  import_array();     //Necessary for handling of NumPY structures
}




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
// NOTE: We use long* instead of int* for numeric arrays as this will work both 
//       for 64 as well as 32 bit systems


#include "Python.h"
#include "numpy/arrayobject.h"
#include "math.h"

#include "util_ext.h"

#define YES 1
#define NO 0


double dist(double x,
	    double y) {
  
  return sqrt(x*x + y*y);
}


int __point_on_line(double x, double y,
		    double x0, double y0,
		    double x1, double y1,
		    double rtol,
		    double atol) {
  /*Determine whether a point is on a line segment

    Input: x, y, x0, x0, x1, y1: where
        point is given by x, y
	line is given by (x0, y0) and (x1, y1)

  */

  double a0, a1, a_normal0, a_normal1, b0, b1, len_a, len_b;
  double nominator, denominator;
  int is_parallel;

  a0 = x - x0;
  a1 = y - y0;

  a_normal0 = a1;
  a_normal1 = -a0;

  b0 = x1 - x0;
  b1 = y1 - y0;

  nominator = fabs(a_normal0*b0 + a_normal1*b1);
  denominator = b0*b0 + b1*b1;
  
  // Determine if line is parallel to point vector up to a tolerance
  is_parallel = 0;
  if (denominator == 0.0) {
    // Use absolute tolerance
    if (nominator <= atol) {
      is_parallel = 1;
    }
  } else {
    // Denominator is positive - use relative tolerance
    if (nominator/denominator <= rtol) {
      is_parallel = 1;
    }    
  }
    
  if (is_parallel) {
    // Point is somewhere on the infinite extension of the line
    // subject to specified absolute tolerance

    len_a = dist(a0, a1); //sqrt(a0*a0 + a1*a1);
    len_b = dist(b0, b1); //sqrt(b0*b0 + b1*b1);

    if (a0*b0 + a1*b1 >= 0 && len_a <= len_b) {
      return 1;
    } else {
      return 0;
    }
  } else {
    return 0;
  }
}


//  public domain function by Darel Rex Finley, 2006
//  http://www.alienryderflex.com/intersect/
//
//  Determines the intersection point of the line segment defined by points A and B
//  with the line segment defined by points C and D.
//
//  Returns YES if the intersection point was found, and stores that point in X,Y.
//  Returns NO if there is no determinable intersection point, in which case X,Y will
//  be unmodified.

int __line_segment_intersection(
        double Ax, double Ay,
        double Bx, double By,
        double Cx, double Cy,
        double Dx, double Dy,
        double *X, double *Y) {

    double distAB, theCos, theSin, newX, ABpos;

    //  Fail if either line segment is zero-length.
    if ( (Ax == Bx && Ay == By) || (Cx == Dx && Cy == Dy) ) return NO ;

    //  Fail if the segments share an end-point.
    if ( (Ax == Cx && Ay == Cy) || (Bx == Cx && By == Cy)
            || (Ax == Dx && Ay == Dy) || (Bx == Dx && By == Dy) ) {
        return NO;
    }

    //  (1) Translate the system so that point A is on the origin.
    Bx -= Ax;
    By -= Ay;
    Cx -= Ax;
    Cy -= Ay;
    Dx -= Ax;
    Dy -= Ay;

    //  Discover the length of segment A-B.
    distAB = sqrt(Bx * Bx + By * By);

    //  (2) Rotate the system so that point B is on the positive X axis.
    theCos = Bx / distAB;
    theSin = By / distAB;
    newX = Cx * theCos + Cy*theSin;
    Cy = Cy * theCos - Cx*theSin;
    Cx = newX;
    newX = Dx * theCos + Dy*theSin;
    Dy = Dy * theCos - Dx*theSin;
    Dx = newX;

    //  Fail if segment C-D doesn't cross line A-B.
    if ( (Cy < 0. && Dy < 0.) || (Cy >= 0. && Dy >= 0.) ) return NO;

    //  (3) Discover the position of the intersection point along line A-B.
    ABpos = Dx + (Cx - Dx) * Dy / (Dy - Cy);

    //  Fail if segment C-D crosses line A-B outside of segment A-B.
    if (ABpos < 0. || ABpos > distAB) return NO;

    //  (4) Apply the discovered position to line A-B in the original coordinate system.
    *X = Ax + ABpos*theCos;
    *Y = Ay + ABpos*theSin;

    //  Success.
    return YES;
}

/*
WORK IN PROGRESS TO OPTIMISE INTERSECTION
int __intersection(double x0, double y0,
		   double x1, double y1) {


    x0 = line0[0,0]; y0 = line0[0,1]
    x1 = line0[1,0]; y1 = line0[1,1]

    x2 = line1[0,0]; y2 = line1[0,1]
    x3 = line1[1,0]; y3 = line1[1,1]

    denom = (y3-y2)*(x1-x0) - (x3-x2)*(y1-y0)
    u0 = (x3-x2)*(y0-y2) - (y3-y2)*(x0-x2)
    u1 = (x2-x0)*(y1-y0) - (y2-y0)*(x1-x0)
        
    if allclose(denom, 0.0):
        # Lines are parallel - check if they coincide on a shared a segment

        if allclose( [u0, u1], 0.0 ):
            # We now know that the lines if continued coincide
            # The remaining check will establish if the finite lines share a segment

            line0_starts_on_line1 = line0_ends_on_line1 =\
            line1_starts_on_line0 = line1_ends_on_line0 = False
                
            if point_on_line([x0, y0], line1):
                line0_starts_on_line1 = True

            if point_on_line([x1, y1], line1):
                line0_ends_on_line1 = True
 
            if point_on_line([x2, y2], line0):
                line1_starts_on_line0 = True

            if point_on_line([x3, y3], line0):
                line1_ends_on_line0 = True                               

            if not(line0_starts_on_line1 or line0_ends_on_line1\
               or line1_starts_on_line0 or line1_ends_on_line0):
                # Lines are parallel and would coincide if extended, but not as they are.
                return 3, None


            # One line fully included in the other. Use direction of included line
            if line0_starts_on_line1 and line0_ends_on_line1:
                # Shared segment is line0 fully included in line1
                segment = array([[x0, y0], [x1, y1]])                

            if line1_starts_on_line0 and line1_ends_on_line0:
                # Shared segment is line1 fully included in line0
                segment = array([[x2, y2], [x3, y3]])
            

            # Overlap with lines are oriented the same way
            if line0_starts_on_line1 and line1_ends_on_line0:
                # Shared segment from line0 start to line 1 end
                segment = array([[x0, y0], [x3, y3]])

            if line1_starts_on_line0 and line0_ends_on_line1:
                # Shared segment from line1 start to line 0 end
                segment = array([[x2, y2], [x1, y1]])                                


            # Overlap in opposite directions - use direction of line0
            if line0_starts_on_line1 and line1_starts_on_line0:
                # Shared segment from line0 start to line 1 end
                segment = array([[x0, y0], [x2, y2]])

            if line0_ends_on_line1 and line1_ends_on_line0:
                # Shared segment from line0 start to line 1 end
                segment = array([[x3, y3], [x1, y1]])                

                
            return 2, segment
        else:
            # Lines are parallel but they do not coincide
            return 4, None #FIXME (Ole): Add distance here instead of None 
            
    else:
        # Lines are not parallel or coinciding
        u0 = u0/denom
        u1 = u1/denom        

        x = x0 + u0*(x1-x0)
        y = y0 + u0*(y1-y0)

        # Sanity check - can be removed to speed up if needed
        assert allclose(x, x2 + u1*(x3-x2))
        assert allclose(y, y2 + u1*(y3-y2))        

        # Check if point found lies within given line segments
        if 0.0 <= u0 <= 1.0 and 0.0 <= u1 <= 1.0: 
            # We have intersection

            return 1, array([x, y])
        else:
            # No intersection
            return 0, None


} 
*/



int __interpolate_polyline(int number_of_nodes,
		int number_of_points,
		double* data,
		double* polyline_nodes,
		long* gauge_neighbour_id,
		double* interpolation_points,
		double* interpolated_values,
		double rtol,
		double atol) {

	int j, i, neighbour_id;
	double x0, y0, x1, y1, x, y;
	double segment_len, segment_delta, slope, alpha;

	for (j=0; j<number_of_nodes; j++) {

		neighbour_id = gauge_neighbour_id[j];

		// FIXME(Ole): I am convinced that gauge_neighbour_id can be discarded, but need to check with John J.
		// Keep it for now (17 Jan 2009)
		// When gone, we can simply interpolate between neighbouring nodes, i.e. neighbour_id = j+1.
		// and the test below becomes something like: if j < number_of_nodes...

		if (neighbour_id >= 0) {
			x0 = polyline_nodes[2*j];
			y0 = polyline_nodes[2*j+1];

			x1 = polyline_nodes[2*neighbour_id];
			y1 = polyline_nodes[2*neighbour_id+1];


			segment_len = dist(x1-x0, y1-y0);
			segment_delta = data[neighbour_id] - data[j];
			slope = segment_delta/segment_len;

			for (i=0; i<number_of_points; i++) {
				x = interpolation_points[2*i];
				y = interpolation_points[2*i+1];

				if (__point_on_line(x, y, x0, y0, x1, y1, rtol, atol)) {
					alpha = dist(x-x0, y-y0);
					interpolated_values[i] = slope*alpha + data[j];
				}
			}
		}
	}

	return 0;
}			       			       


int __triangle_polygon_overlap(double* polygon,
                               double* triangle,
                               int polygon_number_of_vertices)
{
    int i, ii, j, jj, A, B;
    double p0_x, p0_y, p1_x, p1_y, pp_x, pp_y;
    double t0_x, t0_y, t1_x, t1_y, tp_x, tp_y;
    double u_x, u_y, v_x, v_y, w_x, w_y;
    double u_dot_tp, v_dot_tp, v_dot_pp, w_dot_pp;
    double a, b;
    
    p0_x = polygon[0];
    p0_y = polygon[1];
    
    A = 0;
    B = 0;
    
    for (i = 1; i < polygon_number_of_vertices + 1; i++)
    {
        ii = i%polygon_number_of_vertices;
        
        p1_x = polygon[2*ii];
        p1_y = polygon[2*ii + 1];
        
        pp_x = -(p1_y - p0_y);
        pp_y = p1_x - p0_x;
  
        t0_x = triangle[0];
        t0_y = triangle[1];
  
        for (j = 1; j < 4; j++)
        {
            jj = j%3;
                      
            t1_x = triangle[2*jj];
            t1_y = triangle[2*jj + 1];
            
            tp_x = -(t1_y - t0_y); //perpendicular to triangle vector
            tp_y = t1_x - t0_x; //perpendicular to polygon vector
        
            u_x = p1_x - p0_x;
            u_y = p1_y - p0_y;
            v_x = t0_x - p0_x;
            v_y = t0_y - p0_y;
            w_x = t1_x - t0_x;
            w_y = t1_y - t0_y;

            u_dot_tp = (u_x*tp_x) + (u_y*tp_y);
            
            if (u_dot_tp != 0.0f) //Vectors are not parallel
            {
                v_dot_tp = (v_x*tp_x) + (v_y*tp_y);
                v_dot_pp = (v_x*pp_x) + (v_y*pp_y);
                w_dot_pp = (w_x*pp_x) + (w_y*pp_y);
                
                a = v_dot_tp/u_dot_tp;
                b = -v_dot_pp/w_dot_pp;
                             
                if (a >= 0.0f && a <= 1.0f && b >=0.0f && b <=1.0f)
                {
                    return 1; //overlap
                }
                
                if (b >= 0.0f && b <= 1.0f && a > 1.0f)
                {
                    A++; 
                }
                
                if (a >= 0.0f && a <= 1.0f && b > 1.0f)
                {
                    B++; 
                }
                
                if (A == 4 || B == 3)
                {
                    return 1; //overlap
                }
            }
            
            t0_x = t1_x;
            t0_y = t1_y;
        }
        
        p0_x = p1_x;
        p0_y = p1_y;
    }

    return 0; //no overlap
}
                 

int __polygon_overlap(double* polygon,
                      double* triangles,
                      long* indices,
                      int M, //number of triangles
                      int polygon_number_of_vertices)
{
    double* triangle;
    int i, inside_index, outside_index;
    
    inside_index = 0;    // Keep track of triangles that overlap
    outside_index = M - 1; // Keep track of triangles that don't overlap (starting from end)
    
    for (i = 0; i < M; i++)
    {
        triangle = triangles + 6*i;
        
        if (__triangle_polygon_overlap(polygon, 
                                      triangle, 
                                      polygon_number_of_vertices))
        {
            indices[inside_index] = i;
            inside_index++;
        }
        else
        {
            indices[outside_index] = i;
            outside_index -= 1;            
        }
    }
    
    return inside_index;
}              


int __triangle_line_intersect(double* line,
                              double* triangle)
{
    int j, jj, A, B;
    double p0_x, p0_y, p1_x, p1_y, pp_x, pp_y;
    double t0_x, t0_y, t1_x, t1_y, tp_x, tp_y;
    double u_x, u_y, v_x, v_y, w_x, w_y;
    double u_dot_tp, v_dot_tp, v_dot_pp, w_dot_pp;
    double a, b;
    
    p0_x = line[0];
    p0_y = line[1];
    p1_x = line[2];
    p1_y = line[3];
    
    pp_x = -(p1_y - p0_y);
    pp_y = p1_x - p0_x;
    
    A = 0;
    B = 0;
    
    t0_x = triangle[0];
    t0_y = triangle[1];

    for (j = 1; j < 4; j++)
    {
        jj = j%3;
                  
        t1_x = triangle[2*jj];
        t1_y = triangle[2*jj + 1];
        
        tp_x = -(t1_y - t0_y); //perpendicular to triangle vector
        tp_y = t1_x - t0_x; 
    
        u_x = p1_x - p0_x;
        u_y = p1_y - p0_y;
        v_x = t0_x - p0_x;
        v_y = t0_y - p0_y;
        w_x = t1_x - t0_x;
        w_y = t1_y - t0_y;

        u_dot_tp = (u_x*tp_x) + (u_y*tp_y);
        
        if (u_dot_tp != 0.0f) //If vectors are not parallel, continue
        {
            v_dot_tp = (v_x*tp_x) + (v_y*tp_y);
            v_dot_pp = (v_x*pp_x) + (v_y*pp_y);
            w_dot_pp = (w_x*pp_x) + (w_y*pp_y);
            
            a = v_dot_tp/u_dot_tp;
            b = -v_dot_pp/w_dot_pp;
                         
            if (a >= 0.0f && a <= 1.0f && b >=0.0f && b <=1.0f)
            {
                return 1; //intersect
            }
            
            if (a > 1.0f && b >= 0.0f && b <= 1.0f)
            {
                A++; 
            }
            
            if (a < 0.0f && b >= 0.0f && b <= 1.0f)
            {
                B++; 
            }
        }
        
        t0_x = t1_x;
        t0_y = t1_y;
    }
    
    if (A >= 1 && B >= 1)
    {
        return 1; //line sits completely inside a triangle
    }
    
    return 0; //no intersection
}
                 

int __line_intersect(double* line,
                     double* triangles,
                     long* indices,
                     int M) //number of triangles
{
    double* triangle;
    int i, inside_index, outside_index;
    
    inside_index = 0;    // Keep track of triangles that intersect
    outside_index = M - 1; // Keep track of triangles that don't intersect (starting from end)
    
    for (i = 0; i < M; i++)
    {
        triangle = triangles + 6*i;
        
        if (__triangle_line_intersect(line, 
                                      triangle))
        {
            indices[inside_index] = i;
            inside_index++;
        }
        else
        {
            indices[outside_index] = i;
            outside_index -= 1;            
        }
    }
    
    return inside_index;
}              



int __is_inside_triangle(double* point,
			 double* triangle,
			 int closed,
			 double rtol,
			 double atol) {
			 
  double vx, vy, v0x, v0y, v1x, v1y;
  double a00, a10, a01, a11, b0, b1;
  double denom, alpha, beta;
  
  double x, y; // Point coordinates
  int i, j, res;

  x = point[0];
  y = point[1];
  
  // Quickly reject points that are clearly outside
  if ((x < triangle[0]) && 
      (x < triangle[2]) && 
      (x < triangle[4])) return 0;       
      
  if ((x > triangle[0]) && 
      (x > triangle[2]) && 
      (x > triangle[4])) return 0;             
  
  if ((y < triangle[1]) && 
      (y < triangle[3]) && 
      (y < triangle[5])) return 0;       
      
  if ((y > triangle[1]) && 
      (y > triangle[3]) && 
      (y > triangle[5])) return 0;             
  
  
  // v0 = C-A 
  v0x = triangle[4]-triangle[0]; 
  v0y = triangle[5]-triangle[1];
  
  // v1 = B-A   
  v1x = triangle[2]-triangle[0]; 
  v1y = triangle[3]-triangle[1];

  // First check if point lies wholly inside triangle
  a00 = v0x*v0x + v0y*v0y; // innerproduct(v0, v0)
  a01 = v0x*v1x + v0y*v1y; // innerproduct(v0, v1)
  a10 = a01;               // innerproduct(v1, v0)
  a11 = v1x*v1x + v1y*v1y; // innerproduct(v1, v1)
    
  denom = a11*a00 - a01*a10;

  if (fabs(denom) > 0.0) {
    // v = point-A  
    vx = x - triangle[0]; 
    vy = y - triangle[1];     
    
    b0 = v0x*vx + v0y*vy; // innerproduct(v0, v)        
    b1 = v1x*vx + v1y*vy; // innerproduct(v1, v)            
    
    alpha = (b0*a11 - b1*a01)/denom;
    beta = (b1*a00 - b0*a10)/denom;        
    
    if ((alpha > 0.0) && (beta > 0.0) && (alpha+beta < 1.0)) return 1;
  }

  if (closed) {
    // Check if point lies on one of the edges
        
    for (i=0; i<3; i++) {
      j = (i+1) % 3; // Circular index into triangle array
      res = __point_on_line(x, y,
                            triangle[2*i], triangle[2*i+1], 
                            triangle[2*j], triangle[2*j+1], 			    
			    rtol, atol);
      if (res) return 1;
    }
  }
                
  // Default return if point is outside triangle			 
  return 0;			 			 
}			  			       			       


int __separate_points_by_polygon(int M,     // Number of points
				 int N,     // Number of polygon vertices
				 double* points,
				 double* polygon,
				 long* indices,  // M-Array for storage indices
				 int closed,
				 int verbose) {

  double minpx, maxpx, minpy, maxpy, x, y, px_i, py_i, px_j, py_j, rtol=0.0, atol=0.0;
  int i, j, k, outside_index, inside_index, inside;

  // Find min and max of poly used for optimisation when points
  // are far away from polygon
  
  // FIXME(Ole): Pass in rtol and atol from Python

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

  // Begin main loop (for each point)
  inside_index = 0;    // Keep track of points inside
  outside_index = M-1; // Keep track of points outside (starting from end)   
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

    // Optimisation
    if ((x > maxpx) || (x < minpx) || (y > maxpy) || (y < minpy)) {
      // Nothing
    } else {   
      // Check polygon
      for (i=0; i<N; i++) {
        j = (i+1)%N;

        px_i = polygon[2*i];
        py_i = polygon[2*i+1];
        px_j = polygon[2*j];
        py_j = polygon[2*j+1];

        // Check for case where point is contained in line segment
        if (__point_on_line(x, y, px_i, py_i, px_j, py_j, rtol, atol)) {
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

  double x, y, x0, y0, x1, y1, rtol, atol;
  int res;
  PyObject *result;

  // Convert Python arguments to C
  if (!PyArg_ParseTuple(args, "dddddddd", &x, &y, &x0, &y0, &x1, &y1, &rtol, &atol)) {
    PyErr_SetString(PyExc_RuntimeError, 
		    "point_on_line could not parse input");    
    return NULL;
  }


  // Call underlying routine
  res = __point_on_line(x, y, x0, y0, x1, y1, rtol, atol);

  // Return values a and b
  result = Py_BuildValue("i", res);
  return result;
}



// Gateways to Python
PyObject *_interpolate_polyline(PyObject *self, PyObject *args) {
  //
  // _interpolate_polyline(data, polyline_nodes, gauge_neighbour_id, interpolation_points
  //                       interpolated_values):
  //

  
  PyArrayObject
    *data,
    *polyline_nodes,
    *gauge_neighbour_id,
    *interpolation_points,
    *interpolated_values;

  double rtol, atol;  
  int number_of_nodes, number_of_points, res;
  
  // Convert Python arguments to C
  if (!PyArg_ParseTuple(args, "OOOOOdd",
			&data,
			&polyline_nodes,
			&gauge_neighbour_id,
			&interpolation_points,
			&interpolated_values,
			&rtol,
			&atol)) {
    
    PyErr_SetString(PyExc_RuntimeError, 
		    "_interpolate_polyline could not parse input");
    return NULL;
  }

  // check that numpy array objects arrays are C contiguous memory
  CHECK_C_CONTIG(data);
  CHECK_C_CONTIG(polyline_nodes);
  CHECK_C_CONTIG(gauge_neighbour_id);
  CHECK_C_CONTIG(interpolation_points);
  CHECK_C_CONTIG(interpolated_values);

  number_of_nodes = polyline_nodes -> dimensions[0];  // Number of nodes in polyline
  number_of_points = interpolation_points -> dimensions[0];   //Number of points
  

  // Call underlying routine
  res = __interpolate_polyline(number_of_nodes,
			       number_of_points,
			       (double*) data -> data,
			       (double*) polyline_nodes -> data,
			       (long*) gauge_neighbour_id -> data,
			       (double*) interpolation_points -> data,			       
			       (double*) interpolated_values -> data,
			       rtol,
			       atol);			       			       

  // Return
  return Py_BuildValue("");  
}


PyObject *_polygon_overlap(PyObject *self, PyObject *args) 
{
  //
  // _polygon_triangle_overlap(polygon, triangle)
  //

  
  PyArrayObject
    *polygon,
    *triangles,
    *indices;
    
    int res;

  PyObject *result;
      
  // Convert Python arguments to C
  if (!PyArg_ParseTuple(args, "OOO",
			&polygon,
			&triangles,
            &indices)) {
    
    PyErr_SetString(PyExc_RuntimeError, 
		    "_polygon_triangle_overlap could not parse input");
    return NULL;
  }

  // Call underlying routine
  res = __polygon_overlap((double*) polygon->data,
			     (double*) triangles->data,
                 (long*) indices->data,
                 (int) triangles->dimensions[0]/3,
                 (int) polygon->dimensions[0]);			       			       


  // Return result
  result = Py_BuildValue("i", res);
  return result;  
}

PyObject *_line_intersect(PyObject *self, PyObject *args) 
{
  //
  // _polygon_triangle_overlap(polygon, triangle)
  //
  
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
  
  
  PyArrayObject
    *line,
    *triangles,
    *indices;
    
    int res;

  PyObject *result;
      
  // Convert Python arguments to C
  if (!PyArg_ParseTuple(args, "OOO",
			&line,
			&triangles,
            &indices)) {
    
    PyErr_SetString(PyExc_RuntimeError, 
		    "_polygon_triangle_overlap could not parse input");
    return NULL;
  }

  // Call underlying routine
  res = __line_intersect((double*) line->data,
			             (double*) triangles->data,
                         (long*) indices->data,
                         (int) triangles->dimensions[0]/3);			       			       


  // Return result
  result = Py_BuildValue("i", res);
  return result;  
}
     
PyObject *_is_inside_triangle(PyObject *self, PyObject *args) {
  //
  // _is_inside_triangle(point, triangle, int(closed), rtol, atol)
  //

  
  PyArrayObject
    *point,
    *triangle;

  double rtol, atol;  
  int closed, res;

  PyObject *result;
      
  // Convert Python arguments to C
  if (!PyArg_ParseTuple(args, "OOidd",
			&point,
			&triangle,
			&closed,
			&rtol,
			&atol)) {
    
    PyErr_SetString(PyExc_RuntimeError, 
		    "_is_inside_triangle could not parse input");
    return NULL;
  }

  // Call underlying routine
  res = __is_inside_triangle((double*) point -> data,
			     (double*) triangle -> data,
			     closed,
			     rtol,
			     atol);			       			       


  // Return result
  result = Py_BuildValue("i", res);
  return result;  
}
     


/*
PyObject *_intersection(PyObject *self, PyObject *args) {
  //
  // intersection(x0, y0, x1, y1)
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
  res = __intersection(x, y, x0, y0, x1, y1);

  // Return values a and b
  result = Py_BuildValue("i", res);
  return result;
}
*/


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

  // check that points, polygon and indices arrays are C contiguous
  CHECK_C_CONTIG(points);
  CHECK_C_CONTIG(polygon);
  CHECK_C_CONTIG(indices);

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
  //{"_intersection", _intersection, METH_VARARGS, "Print out"},  
  {"_separate_points_by_polygon", _separate_points_by_polygon, 
                                 METH_VARARGS, "Print out"},
  {"_interpolate_polyline", _interpolate_polyline, 
                                 METH_VARARGS, "Print out"},				 
  {"_polygon_overlap", _polygon_overlap, 
                                 METH_VARARGS, "Print out"},
  {"_line_intersect", _line_intersect, 
                                 METH_VARARGS, "Print out"},                               
  {"_is_inside_triangle", _is_inside_triangle, 
                                 METH_VARARGS, "Print out"},
  {NULL, NULL, 0, NULL}   /* sentinel */
};


// Module initialisation
void initpolygon_ext(void){
  Py_InitModule("polygon_ext", MethodTable);

  import_array();     //Necessary for handling of NumPY structures
}




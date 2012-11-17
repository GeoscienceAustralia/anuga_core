#include "Python.h"
#include "numpy/arrayobject.h"
#include <stdio.h>
#include <math.h>
#include <malloc.h>

#define DDATA(p) ((double*)(((PyArrayObject *)p)->data))
#define IDATA(p) ((long*)(((PyArrayObject *)p)->data))

#define MIN(a, b) (((a)<=(b))?(a):(b))
#define MAX(a, b) (((a)>(b))?(a):(b))
#define ABS(a) ( (a) >= 0 ? (a) : -(a))

#define ORI_LEFT  0
#define ORI_RIGHT 1
#define ORI_UP	  2
#define ORI_DOWN  3

#define EPSILON 1.0e-12

typedef struct{
	double x_max;
	double x_min;
	double y_max;
	double y_min;
}EXTENT, *PTR_EXTENT;

double point_dot(double *p1, double *p2)
{
	return p1[0]*p2[0]+p1[1]*p2[1];
}

double *point_sub(double *p1, double *p2)
{
	double *res = malloc( 2*sizeof( double ) );
	
	res[0] = p1[0] - p2[0];
	res[1] = p1[1] - p2[1];
	
	return res;
}

void get_tri_extent(double *vertices, PTR_EXTENT out)
{
	double x1, x2, x3, y1, y2, y3;

	x1 = vertices[0];
	x2 = vertices[2];
	x3 = vertices[4];
	y1 = vertices[1];
	y2 = vertices[3];
	y3 = vertices[5];	

	out->x_min = MIN( x1, MIN( x2, x3 ) );
	out->x_max = MAX( x1, MAX( x2, x3 ) );
	out->y_min = MIN( y1, MIN( y2, y3 ) );
	out->y_max = MAX( y1, MAX( y2, y3 ) );
}

void get_tri_vertices( double *vertices, \
			long *volumes, \
			int tri_id, \
			double *out, \
			double *v1,  \
			double *v2,  \
			double *v3 )
{
	out[0] = vertices[volumes[tri_id*3]*2];
	out[1] = vertices[volumes[tri_id*3]*2+1];
	out[2] = vertices[volumes[tri_id*3+1]*2];
	out[3] = vertices[volumes[tri_id*3+1]*2+1];
	out[4] = vertices[volumes[tri_id*3+2]*2];
	out[5] = vertices[volumes[tri_id*3+2]*2+1];
	

	if (v1) {
		v1[0]=vertices[volumes[tri_id*3]*2]; 
		v1[1]=vertices[volumes[tri_id*3]*2+1];
	}
	if (v2) {
		v2[0]=vertices[volumes[tri_id*3+1]*2]; 
		v2[1]=vertices[volumes[tri_id*3+1]*2+1];
	}
	if (v3) {
		v3[0]=vertices[volumes[tri_id*3+2]*2]; 
		v3[1]=vertices[volumes[tri_id*3+2]*2+1];
	}
}

void get_tri_norms( double *norms, int tri_id, 
		       double *n1, double *n2, double *n3)
{
	n1[0] = norms[tri_id*6];
	n1[1] = norms[tri_id*6+1];
	n2[0] = norms[tri_id*6+2];
	n2[1] = norms[tri_id*6+3];
	n3[0] = norms[tri_id*6+4];
	n3[1] = norms[tri_id*6+5];
}

double *init_norms( double *vertices, long *volumes, int num_tri  )
{
	int i;
	double x1, x2, x3, y1, y2, y3;
	double xn1, yn1, xn2, yn2, xn3, yn3;
	double l1, l2, l3;
	double *norms;

	norms = malloc( num_tri*6*sizeof( double ) );

	for ( i = 0; i < num_tri; i++ ) {
		x1 = vertices[volumes[i*3]*2];
		x2 = vertices[volumes[i*3+1]*2];
		x3 = vertices[volumes[i*3+2]*2];
		y1 = vertices[volumes[i*3]*2+1];
		y2 = vertices[volumes[i*3+1]*2+1];
		y3 = vertices[volumes[i*3+2]*2+1];

		xn1 = x3 - x2;
		yn1 = y3 - y2;
		l1  = sqrt( xn1*xn1 + yn1*yn1 );
		
		if ( l1 ) { xn1 /= l1; yn1 /= l1; }

		xn2 = x1 - x3;
		yn2 = y1 - y3;
		l2 = sqrt( xn2*xn2 + yn2*yn2 );

		if ( l2 ) { xn2 /= l2; yn2 /= l2; }

		xn3 = x2 - x1;
		yn3 = y2 - y1;
		l3  = sqrt( xn3*xn3 + yn3*yn3 );
		
		if ( l3 ) { xn3 /= l3; yn3 /= l3; }

		norms[i*6]   = yn1;
		norms[i*6+1] = -xn1;
		
		norms[i*6+2] = yn2;
		norms[i*6+3] = -xn2;
		
		norms[i*6+4] = yn3;
		norms[i*6+5] = -xn3;
	}

	return norms;
}

// remove nodes that are not in any triangles
void remove_lone_verts( double **verts, int *volumes )
{
	
}

int _point_on_line(double x, double y,
		   double x0, double y0,
		   double x1, double y1,
		   double rtol,
		   double atol) 
{

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

int _is_inside_triangle(double *point,
			double *triangle,
			int closed,
			double rtol,
			double atol) 
{			 
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
      j = (i+1) % 3; // Circular index into triangle vertices
      res = _point_on_line(x, y,
                            triangle[2*i], triangle[2*i+1], 
                            triangle[2*j], triangle[2*j+1], 			    
			    rtol, atol);
      if (res) return 1;
    }
  }
                
  // Default return if point is outside triangle			 
  return 0;			 			 
}

void build_interpolation_matrix( double *vertices, 
				 int num_vert,
				 long *volumes, 
				 double *norms,
				 int num_tri, 
				 double *grid_points, 
				 int x_dimension,
				 int y_dimension,
				 double *vertex_val,
				 double *grid_val )
{
	int i, j, k;
	int x_min, x_max, y_min, y_max, point_index;
	int x_ori, y_ori;
	double x_dist, y_dist, x_base, y_base;
	double sigma0, sigma1, sigma2;
	double fraction, intpart;
	double *triangle, *point;
	double *v1, *v2, *v3;
	double *n1, *n2, *n3;
	double val1, val2;
	PTR_EXTENT extent;

	x_ori = (grid_points[2] - grid_points[0]) > 0 ? ORI_RIGHT : ORI_LEFT;
	y_ori = (grid_points[y_dimension*2+1] - grid_points[1]) > 0 ? ORI_UP: ORI_DOWN;
	
	x_dist = ABS(grid_points[2] - grid_points[0]);
	y_dist = ABS(grid_points[y_dimension*2+1] - grid_points[1]);

	x_base = grid_points[0];
	y_base = grid_points[1];

	triangle = malloc( 6*sizeof(double) );
	point    = malloc( 2*sizeof(double) );
	extent	 = malloc( sizeof( EXTENT ) );

	v1 = malloc( 2*sizeof(double) );
	v2 = malloc( 2*sizeof(double) );
	v3 = malloc( 2*sizeof(double) );

	n1 = malloc( 2*sizeof(double) );
	n2 = malloc( 2*sizeof(double) );
	n3 = malloc( 2*sizeof(double) );

	for ( i = 0; i < num_tri; i++ ) {

		get_tri_vertices( vertices, volumes, i, triangle, v1, v2, v3);
		get_tri_norms( norms, i, n1, n2, n3 );
		get_tri_extent( triangle, extent );

		fraction = (x_ori == ORI_RIGHT) \
				? modf( (extent->x_min - x_base)/x_dist, &intpart )
				: modf( ABS(extent->x_max - x_base)/x_dist, &intpart );
		x_min = intpart;//(fraction<=EPSILON) ? intpart : (intpart+1);
		x_min = (x_min < 0) ? 0 : x_min; 

		fraction = (x_ori == ORI_RIGHT) \
				? modf( ABS(extent->x_max - x_base)/x_dist, &intpart )
				: modf( ABS(extent->x_min - x_base)/x_dist, &intpart );
		x_max = intpart;//(x_dist - fraction)<=EPSILON ? intpart+1 : intpart;
		x_max = (x_max > (y_dimension-1)) ? (y_dimension-1) : x_max;

		fraction = (y_ori == ORI_UP) \
				? modf( (extent->y_min - y_base)/y_dist, &intpart )
				: modf( ABS(extent->y_max - y_base)/y_dist, &intpart );
		y_min = intpart;//(fraction <= EPSILON) ? intpart : (intpart+1);
		y_min = (y_min < 0 ) ? 0 : y_min;

		fraction = (y_ori == ORI_UP) \
				? modf( ABS(extent->y_max - y_base)/y_dist, &intpart )
				: modf( ABS(extent->y_min - y_base)/y_dist, &intpart );
		y_max = intpart;//(y_dist - fraction)<=EPSILON ? (intpart+1) : intpart;
		y_max = (y_max > (x_dimension-1)) ? (x_dimension-1) : y_max;
		
		if ( x_max >= 0 && y_max >= 0 ) {
		for ( j = y_min; j <= y_max; j++ ) {
			for ( k = x_min; k <= x_max; k++ ) {
				// iterate through points within a small region
				point_index = j*y_dimension+k;

				point[0] = grid_points[point_index*2];
				point[1] = grid_points[point_index*2+1];

				if ( _is_inside_triangle( point, triangle, \
							  1, 1.0e-12, 1.0e-12 ) ) {
					
					val1 = point_dot( point_sub( point, v2 ), n1 );
					val2 = point_dot( point_sub( v1, v2 ), n1 );	
					sigma0 = val2 ? val1/val2 : 0;	
						
					val1 = point_dot( point_sub( point, v3 ), n2 );
					val2 = point_dot( point_sub( v2, v3 ), n2 );
					sigma1 = val2 ? val1/val2 : 0;
			
					val1 = point_dot( point_sub( point, v1 ), n3 );
					val2 = point_dot( point_sub( v3, v1 ), n3 );
					sigma2 = val2 ? val1/val2 : 0;

						
					grid_val[point_index] = sigma0*vertex_val[volumes[i*3]] + \
								sigma1*vertex_val[volumes[i*3+1]] + \
								sigma2*vertex_val[volumes[i*3+2]];
				}
			}
		}
		}
	}

	free(triangle);
	free(point);
	free(extent);

	free(v1);
	free(v2);
	free(v3);

	free(n1);
	free(n2);
	free(n3);
}

static PyObject *calc_grid_values( PyObject *self, PyObject *args )
{
	int i, ok, num_tri, num_vert, ncol, nrow;
	long *volumes; 
	double nodata_val;
	double *vertices;
	double *result;
	double *grid;
	double *grid_val;
	double *norms;
	PyObject *pyobj_vertex_points;
	PyObject *pyobj_volumes;
	PyObject *pyobj_result;
	PyObject *pyobj_grid;
	PyObject *pyobj_grid_val;

	ok = PyArg_ParseTuple( args, "iidOOOOO", 
				&nrow,
				&ncol,
				&nodata_val,
			        &pyobj_grid, 
				&pyobj_vertex_points, 
				&pyobj_volumes, 
				&pyobj_result,
				&pyobj_grid_val );
	if( !ok ){
		fprintf( stderr, "interpolate func: argument parsing error\n" );
		exit(1);
	}

	// get data from python objects
	vertices = DDATA( pyobj_vertex_points );
	grid 	 = DDATA( pyobj_grid );
	result	 = DDATA( pyobj_result );
	grid_val = DDATA( pyobj_grid_val );
	volumes  = IDATA( pyobj_volumes );

	num_tri  = ((PyArrayObject*)pyobj_volumes)->dimensions[0];
	num_vert = ((PyArrayObject*)pyobj_vertex_points)->dimensions[0]/2;

	// remove unused vertices
	// remove_lone_verts( &vertices, volumes );	

	// init triangle array
	norms = init_norms( vertices, volumes, num_tri );

	// evaluate grid
	for ( i = 0 ; i < nrow*ncol; i++ ) 
		grid_val[i] = nodata_val;

	build_interpolation_matrix( vertices, num_vert, volumes, norms, num_tri, \
				    grid, nrow, ncol,   	\
				    result, grid_val );

	free(norms);

	return Py_BuildValue("i", 0);
}

static PyMethodDef interpolate_ext_methods[] = {
	{"eval_grid", calc_grid_values, METH_VARARGS},
	{NULL, NULL}
};

void initinterpolate_ext( )
{
	(void) Py_InitModule( "interpolate_ext", interpolate_ext_methods );
	
	import_array( );
}


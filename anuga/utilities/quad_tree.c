#include "quad_tree.h"

// **************** UTILITIES ***********************

static void *emalloc(size_t amt,char * location)
{
    void *v = malloc(amt);  
    if(!v){
        fprintf(stderr, "out of mem in quad_tree: %s\n",location);
        exit(EXIT_FAILURE);
    }
    return v;
};

double dot_points(double p1x,double p1y,double p2x,double p2y)
{
    
    return p1x * p2x + p1y * p2y;
    
};

// ***************************************************

// ******************* TRIANGLE **********************

triangle * new_triangle(int index,double x1,double y1,double x2,double y2,double x3,double y3)
{

	triangle * T = emalloc(sizeof(triangle),"new_triangle"); 
	 
	T->index = index;
	T->x1 = x1; T->x2 = x2; T->x3 = x3;
    T->y1 = y1; T->y2 = y2; T->y3 = y3;
    T->next = NULL;

    // Calculate the normals
    // normal for vector from x1->x2
    double nx_temp,ny_temp,dot_temp;
    ny_temp = ( T->x2 - T->x1 );
    nx_temp = -( T->y2 - T->y1 );
    dot_temp = dot_points(nx_temp, ny_temp, nx_temp, ny_temp);
    T->nx3 = nx_temp / sqrt(dot_temp);
    T->ny3 = ny_temp / sqrt(dot_temp);
    if( dot_points(T->nx3, T->ny3, T->x3 - T->x2, T->y3 - T->y2) > 0 ){
    	T->nx3 = -T->nx3;
    	T->ny3 = -T->ny3;
    }
    // normal for vector from x2->x3
    ny_temp = ( T->x3 - T->x2 );
    nx_temp = -(T->y3 - T ->y2 );
    dot_temp = dot_points(nx_temp, ny_temp, nx_temp, ny_temp);
    T->nx1 = nx_temp / sqrt(dot_temp);
    T->ny1 = ny_temp / sqrt(dot_temp);
    if( dot_points(T->nx1, T->ny1, T->x1 - T->x3, T->y1 - T->y3) > 0 ){
    	T->nx1 = -T->nx1;
    	T->ny1 = -T->ny1;
    }
    // normal for vector from x3->x1
    ny_temp = ( T->x1 - T->x3 );
    nx_temp = -( T->y1 - T->y3 );
    dot_temp = dot_points(nx_temp, ny_temp, nx_temp, ny_temp);
    T->nx2 = nx_temp / sqrt(dot_temp);
    T->ny2 = ny_temp / sqrt(dot_temp);
    if( dot_points(T->nx2, T->ny2, T->x2 - T->x1, T->y2 - T->y1) > 0 ){
    	T->nx2 = -T->nx2;
    	T->ny2 = -T->ny2;
    }

    return T;
};

void delete_triangle_list(triangle * T)
{
	while (T != NULL){
     triangle * next = T->next;
	   free(T);
	   T = NULL;
     T = next;
  }
};

double * calculate_sigma(triangle * T,double x,double y)
{


	// FIXME SR: Should remove this malloc and just pass a pointer to array
	double  * ret_sigma = malloc(3 * sizeof(double));
	ret_sigma[0] = dot_points(x - T->x2, y - T->y2, T->nx1, T->ny1)/
					dot_points(T->x1 - T->x2, T->y1 - T->y2, T->nx1, T->ny1);
	ret_sigma[1] = dot_points(x - T->x3, y - T->y3, T->nx2, T->ny2)/
					dot_points(T->x2 - T->x3, T->y2 - T->y3, T->nx2, T->ny2);
	ret_sigma[2] = dot_points(x - T->x1, y - T->y1, T->nx3, T->ny3)/
					dot_points(T->x3 - T->x1, T->y3 - T->y1, T->nx3, T->ny3);
	return ret_sigma;				
};

double dist(double x,
	    double y) {
  
  return sqrt(x*x + y*y);
}

int __point_on_line(double x, double y,
                    double x0, double y0,
                    double x1, double y1,
                    double rtol,
                    double atol)
{
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

  nominator = fabs(a_normal0 * b0 + a_normal1 * b1);
  denominator = b0 * b0 + b1 * b1;

  // Determine if line is parallel to point vector up to a tolerance
  is_parallel = 0;
  if (denominator == 0.0)
  {
    // Use absolute tolerance
    if (nominator <= atol)
    {
      is_parallel = 1;
    }
  }
  else
  {
    // Denominator is positive - use relative tolerance
    if (nominator / denominator <= rtol)
    {
      is_parallel = 1;
    }
  }

  if (is_parallel)
  {
    // Point is somewhere on the infinite extension of the line
    // subject to specified absolute tolerance

    //        len_a = dist(a0, a1); //sqrt(a0*a0 + a1*a1);
    //        len_b = dist(b0, b1); //sqrt(b0*b0 + b1*b1);

    //        if (a0*b0 + a1*b1 >= 0 && len_a <= len_b) {
    if (a0 * b0 + a1 * b1 >= 0)
    { // inside line segment from one end point
      double len_a2 = a0 * a0 + a1 * a1;
      if (len_a2 <= denominator)
      { // inside line segment from the other end point
        return 1;
      }
    }
  }
  return 0;
};

int __is_inside_triangle(double* point,
			 double* triangle,
			 int closed,
			 double rtol,
			 double a_tol) {
			 
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
			    rtol, a_tol);
      if (res) return 1;
    }
  };
  // Default return if point is outside triangle			 
  return 0;			 			 
}			  			       	

int triangle_contains_point(triangle * T,double pointx,double pointy)
{
    
//    double v0x,v0y,v1x,v1y,v2x,v2y,dot00,dot01,dot02,dot11,dot12,invDenom,u,v;
//    
//    v0x = T->x3 - T->x1;
//    v0y = T->y3 - T->y1;
//    v1x = T->x2 - T->x1;
//    v1y = T->y2 - T->y1;
//    v2x = pointx - T->x1;
//    v2y = pointy - T->y1;
//    
//    dot00 = dot_points(v0x, v0y, v0x, v0y);
//    dot01 = dot_points(v0x, v0y, v1x, v1y);
//    dot02 = dot_points(v0x, v0y, v2x, v2y);
//    dot11 = dot_points(v1x, v1y, v1x, v1y);
//    dot12 = dot_points(v1x, v1y, v2x, v2y);
//    
//    invDenom = 1/(dot00*dot11-dot01*dot01);
//    u=(dot11 * dot02 - dot01 * dot12) * invDenom;
//    v=(dot00 * dot12 - dot01 * dot02) * invDenom;
//    
//    if( u>=0 && v>=0 && v+u<1) return 1;
//   
//    else return 0;
      double tri[6];
      tri[0]=T->x1; tri[1]=T->y1; tri[2]=T->x2; tri[3]=T->y2; tri[4]=T->x3; tri[5]=T->y3;
      double point[2];
      point[0] = pointx; point[1] = pointy;
      double rtol=1.0e-12;
      double a_tol=1.0e-12;
      int closed = 1;

      return __is_inside_triangle(point,
			 tri,
			 closed,
			 rtol,
			 a_tol); 


};



// ***************************************************

// ***************** quad_tree *********************

quad_tree * new_quad_tree(double xmin, double xmax, double ymin, double ymax)
{

	quad_tree * ret = emalloc(sizeof(quad_tree),"new_quad_tree");
	ret -> xmin = xmin; ret-> xmax = xmax; 
	ret -> ymin = ymin; ret -> ymax = ymax;
	ret -> parent = NULL;
	ret -> q[0] = NULL; ret -> q[1] = NULL; ret -> q[2] = NULL; ret -> q[3] = NULL;
	ret -> leaves = NULL;
	ret -> end_leaves = NULL;
  ret -> count = 0;
	return ret; 

};

void delete_quad_tree(quad_tree * quadtree)
{

  quad_tree_ll * nodelist = new_quad_tree_ll(quadtree,0);
  quad_tree_ll * last = nodelist;
  quad_tree_ll * temp;
  int i;

  while(nodelist !=NULL){
      
      quadtree=nodelist->tree;
      // if children have been added, add to the linked list

      if (quadtree->q[0]!=NULL){
          for (i=0;i<4;i++){
              quad_tree_ll * child = new_quad_tree_ll(quadtree->q[i],0);
              last->next=child;
              last=child;
          }
      }
      
      if (quadtree->leaves!=NULL){
          delete_triangle_list(quadtree->leaves);
      }

      free(quadtree);
      quadtree=NULL;

      temp = nodelist;
      nodelist=nodelist->next;
      free(temp);
  }

};

void quad_tree_make_children(quad_tree *node){

  //double xmid = (node->xmin+node->xmax)/2;
  //double ymid = (node->ymin+node->ymax)/2;
  double width = (node->xmax-node->xmin);
  double height = (node->ymax-node->ymin);
  // add quads 1-4
  // include border expansion
  double border=0.55;
  node->q[0] = new_quad_tree(node->xmax-width*border,node->xmax,node->ymax-height*border,node->ymax);
  node->q[0]->parent = node; 
  
  node->q[1] = new_quad_tree(node->xmin,node->xmin+width*border,node->ymax-height*border,node->ymax);
  node->q[1]->parent = node; 
  
  node->q[2] = new_quad_tree(node->xmin,node->xmin+width*border,node->ymin,node->ymin+height*border);
  node->q[2]->parent = node; 
  
  node->q[3] = new_quad_tree(node->xmax-width*border,node->xmax,node->ymin,node->ymin+height*border);
  node->q[3]->parent = node; 

}

void quad_tree_add_triangle_to_list(quad_tree *node,triangle *T){

  if (node->leaves == NULL){
    // no current leaves
    node->leaves = T;
    node->end_leaves = T;
  } else {
    node->end_leaves->next = T;
    node->end_leaves = T;
  }

}

void quad_tree_insert_triangle(quad_tree *node,triangle *T)
{
	
	// find the quadrant of the current node's extents in which the
	// point lies (zero if intersects center of extents axes).

	int quad = trivial_contain_split(node,T);

  // always increase point count, as storing the total in tree below
  node->count+=1;

	if (quad != 0){
		// if current node has no children yet, split:
		if(node->q[0] == NULL){
			
			quad_tree_make_children(node); 
			
	    }
	    // insert triangle into node corresponding to given quadrant
	   	quad_tree_insert_triangle(node->q[quad-1],T);
		return;
		
	}
	// if triangle intersects the center axes of the node's extents, insert
	// the triangle here
	quad_tree_add_triangle_to_list(node,T);

};


int trivial_contain_split(quad_tree *node, triangle *T){

	int p1 = trivial_contain_split_point(node,T->x1,T->y1);
	int p2 = trivial_contain_split_point(node,T->x2,T->y2);
	int p3 = trivial_contain_split_point(node,T->x3,T->y3);
	if(p1 == p2 && p2 == p3){
	 	return p1;
	}
	return 0;
};

int trivial_contain_split_point(quad_tree *node, double xp,double yp)
{

	double midx = (node->xmin+node->xmax)/2;
	double midy = (node->ymin+node->ymax)/2;

	int ret=0;
	
	if (midx < xp){
		// quad 1 or 4
		if (midy < yp){
			ret = 1;
		} else if (midy > yp){
			ret = 4;
		}
	} else if (midx > xp){
		//quad 2 or 3
		if (midy < yp){
			ret = 2;
		} else if (midy > yp){
			ret = 3;
		}
	}
	return ret;
};

triangle * search_triangles_of_quad_tree(quad_tree * node,double xp,double yp){
	
    triangle * T = node->leaves;

    while(T != NULL){
        if(triangle_contains_point(T,xp,yp)){
            return T; // Triangle contains point so return
        }
        T = T->next;
    }
	return T; // should be NULL if this is reached
};

// Searches the quad tree starting at 'cur_quad_tree' for the 
// point, returning the triangle that contains it, or NULL
// if no triangle is found.
triangle * search(quad_tree * node, double xp, double yp){

    triangle * return_T = NULL;

    if(node->leaves!=NULL)
    {
    	return_T = search_triangles_of_quad_tree(node,xp,yp);
    }
    if(return_T != NULL) return return_T;
	    
    else
    {
        if(node->q[0]!=NULL) // look for child to search
        {
            //find correct quadrant to search
            int quad = trivial_contain_split_point(node,xp,yp);
            
            if (quad!=0)
            {
            return_T = search(node->q[quad-1],xp,yp); // search child
            }

        return return_T; // return NULL pointer as no triangle
        }
	}
	return return_T; // should not be reached
};

int quad_tree_node_count(quad_tree * tree)
{
  int node_count = 1;
  if (tree->q[0]!=NULL){
      int i;
      for(i=0;i<4;i++){
        node_count+=quad_tree_node_count(tree->q[i]);
      }
  }
  return node_count;
};

// ***************************************************

// ***************** quad_tree_ll *******************

quad_tree_ll * new_quad_tree_ll(quad_tree * start,int index){
    quad_tree_ll * list = malloc(sizeof(quad_tree_ll));
    list->tree = start;
    list->next = NULL;
    list->index = index;
    return list;
}

// ***************************************************

// ***************** queue_ll *******************

queue_ll * new_queue_ll(int node){
    queue_ll * list = malloc(sizeof(queue_ll));
    list->node=node;
    list->next = NULL;
    return list;
}

// ***************************************************

//------------------------------------------------------------------------------
// quad_tree.h - Header for c implementation of the generic quad_tree. The 
// 				   tree stores a set of triangles and then facilitates quickly
//                 searching for a triangles in the list containing a given 
//                 point.
// 				   Tree has 'new' and 'delete' methods that should be called to
//                 make sure memory is allocated and freed after use.
// author: Padarn Wilson, date: 25/10/12 
//-------------------------------------------------------------------------------
#include <stdio.h>   /* gets */
#include <stdlib.h>  /* atoi, malloc */
#include <string.h>  /* strcpy */
#include <stdint.h>  /* int64_t uint64_t */
#include <math.h>

#ifndef quad_tree_H
#define quad_tree_H


// ***********************************************************************
// triangle struct - Struct that stores data for triangles being inseted
//                   into the tree. Triangles form a linked list for easy
//                   storage.

typedef struct triangle{

	// points defining triangle.
	double x1,y1;
	double x2,y2;
	double x3,y3;

	// index stores the triangles unique id.
	int64_t index;

	// outward normal vectors of triangles sides.
	double nx1,ny1;
	double nx2,ny2;
	double nx3,ny3;

	// next triangle turns the struct into a linked list.
	struct triangle * next;  

} triangle;

// creates a new triangle and returns a pointer to the malloc'ed memory. 
triangle * new_triangle(int64_t index, double x1, double x2, double x3,
						double y1, double y2, double y3);

// deletes entire list of triangles
void delete_triangle_list(triangle * T);

// take a point and calculate 'sigma' with the given triangle. returns
// a pointer to malloc'ed memory of a double array.
double * calculate_sigma(triangle * T,double x,double y);

// Tests to see if a triangle contains a given point,
// returns a int64_t value 0 false, 1 true.
int64_t triangle_contains_point(triangle * T,double pointx,double pointy);

//**************************************************************************


//**************************************************************************
// quad_tree struct - Each quad_tree struct is a node in the tree, and 
//                      can be treated as an independent quad_tree.

typedef struct quad_tree{

	// rectangular extents of current node
	double xmin,xmax,ymin,ymax;
	int64_t count;
	// parent and children of quad_tree
	struct quad_tree *parent;
	struct quad_tree *q[4];

	// triangle stored in this node - leaves is a linked list of triangles
	triangle * leaves;
	// triangle end_leaves allows easy adding onto end of triangles
	triangle * end_leaves;

} quad_tree;

// returns a new quad_tree pointer with malloc'ed memory. Inputs are the 
// extents of the node.
quad_tree* new_quad_tree(double xmin, double xmax, double ymin, double ymax);

// delete a search tree - recursively deletes all children.
void delete_quad_tree(quad_tree * tree);

// add a new triangle to the quad_tree
void quad_tree_insert_triangle(quad_tree *node,triangle *T);

// returns the quadrant of the quad_tree containing the point, or 0 if intersects
// center axes
int64_t trivial_contain_split_point(quad_tree *node, double xp,double yp);

// returns the quadrant of the quad_tree containing the triangle, or 0 if intersects
// center axes
int64_t trivial_contain_split(quad_tree *node, triangle *T);

// returns the triangle in the quad_tree's leaves containing the point or NULL
// if none of the triangles on the current quad_tree contain it.
triangle * search_triangles_of_quad_tree(quad_tree * node,double xp,double yp);

// search the tree for a triangle containing the given point. returns the triangle
// if found, and NULL otherwise
triangle * search(quad_tree * node ,double xp, double yp);

// return number of noes in tree
int64_t quad_tree_node_count(quad_tree * tree);

// split the node to make 4 children
void quad_tree_make_children(quad_tree *node);

// add a triangle to the nodes leaves
void quad_tree_add_triangle_to_list(quad_tree *node,triangle *T);


//**************************************************************************

//**************************************************************************
// quad_tree_ll struct - Used for serialising  
//  

typedef struct quad_tree_ll {

    void * tree;
    struct quad_tree_ll * next;
    int64_t index;

} quad_tree_ll;

quad_tree_ll * new_quad_tree_ll(quad_tree * start, int64_t index);

//**************************************************************************

//**************************************************************************
// queue_ll struct - Used for deserialising  
//  

typedef struct queue_ll {

    int64_t node;
    struct queue_ll * next;


} queue_ll;

queue_ll * new_queue_ll(int64_t node);

//**************************************************************************

#endif

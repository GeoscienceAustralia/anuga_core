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
#include <stdint.h>  /* anuga_int uanuga_int */
#include <math.h>
#include "anuga_typedefs.h" /* in utilities */

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
	anuga_int index;

	// outward normal vectors of triangles sides.
	double nx1,ny1;
	double nx2,ny2;
	double nx3,ny3;

	// next triangle turns the struct into a linked list.
	struct triangle * next;  

} triangle;

// creates a new triangle and returns a pointer to the malloc'ed memory. 
triangle * new_triangle(anuga_int index, double x1, double x2, double x3,
						double y1, double y2, double y3);

// deletes entire list of triangles
void delete_triangle_list(triangle * T);

// take a point and calculate 'sigma' (its three barycentric weights) with
// the given triangle, written into the caller's 3-element array.
void calculate_sigma(triangle * T,double x,double y,double sigma[3]);

// Tests to see if a triangle contains a given point,
// returns a anuga_int value 0 false, 1 true.
anuga_int triangle_contains_point(triangle * T,double pointx,double pointy);

//**************************************************************************


//**************************************************************************
// quad_tree struct - Each quad_tree struct is a node in the tree, and 
//                      can be treated as an independent quad_tree.

typedef struct quad_tree{

	// rectangular extents of current node
	double xmin,xmax,ymin,ymax;
	anuga_int count;
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

// add a new triangle to the quad_tree. Returns 0, or -1 if a node split
// could not be allocated (T is then NOT owned by the tree).
anuga_int quad_tree_insert_triangle(quad_tree *node,triangle *T);

// returns the quadrant of the quad_tree containing the point, or 0 if intersects
// center axes
anuga_int trivial_contain_split_point(quad_tree *node, double xp,double yp);

// returns the quadrant of the quad_tree containing the triangle, or 0 if intersects
// center axes
anuga_int trivial_contain_split(quad_tree *node, triangle *T);

// returns the triangle in the quad_tree's leaves containing the point or NULL
// if none of the triangles on the current quad_tree contain it.
triangle * search_triangles_of_quad_tree(quad_tree * node,double xp,double yp);

// search the tree for a triangle containing the given point. returns the triangle
// if found, and NULL otherwise
triangle * search(quad_tree * node ,double xp, double yp);

// return number of noes in tree
anuga_int quad_tree_node_count(quad_tree * tree);

// split the node to make 4 children. Returns 0, or -1 on allocation failure
// (the node is left childless).
anuga_int quad_tree_make_children(quad_tree *node);

// add a triangle to the nodes leaves
void quad_tree_add_triangle_to_list(quad_tree *node,triangle *T);


//**************************************************************************

//**************************************************************************
// quad_tree_ll struct - Used for serialising  
//  

typedef struct quad_tree_ll {

    void * tree;
    struct quad_tree_ll * next;
    anuga_int index;

} quad_tree_ll;

quad_tree_ll * new_quad_tree_ll(quad_tree * start, anuga_int index);

//**************************************************************************

//**************************************************************************
// queue_ll struct - Used for deserialising  
//  

typedef struct queue_ll {

    anuga_int node;
    struct queue_ll * next;


} queue_ll;

queue_ll * new_queue_ll(anuga_int node);

//**************************************************************************

#endif

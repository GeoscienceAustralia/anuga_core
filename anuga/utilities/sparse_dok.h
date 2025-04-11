/*
* Sparse Matrix class implemented in DOK format. 
*
* The DOK (dictionary of keys) format is implemented using 
* uthash, which is a doubly linked list and also a hash table.
* The hash table is populated by 'edge_t' objects which contain
* a key specifying the (i,j) index of the matrix, and the value
* at that index. If no entry exists, the value is assumed 0.
* 
* The functions in this class which create new objects return
* pointers to new malloced memory, which much be cleaned up.
* 
* Padarn Wilson, ANU 2012
*/


#include <stdio.h>   /* gets */
#include <stdlib.h>  /* atoi, malloc */
#include <string.h>  /* strcpy */
#include <stdint.h>  /* int64_t uint64_t */
#include "math.h"
#include "uthash.h"     /* in utilities */
#include "sparse_csr.h"

#ifndef SPARSE_DOK_H
#define SPARSE_DOK_H

// Struct edge_key_t to store the i,j position in the matrix within
// a key of the hashtable
typedef struct  {
    int64_t i;
    int64_t j;
} edge_key_t;

// Struct edge_t is a basic element of the hash table. By including 
// the UT_hash_handle variable, it can be used in the methods defined
// for uthash hashtables
typedef struct {
    edge_key_t key;              /* key of form i , j */
    double entry;                /* the value stored at this entry */
    UT_hash_handle hh;         /* makes this structure hashable */
} edge_t;

// Struct sparse_dok defining a sparse matrix. Keeps track of the 
// number of entries and rows in the matrix, and stores the hashtable.
// PADARN NOTE: For efficiency, it might be sensible to pre-allocated
// the max number of rows/cols in the hash table, so that the hash table
// can be made an appropriate size.
typedef struct {
	edge_t *edgetable;
	int64_t num_entries;
	int64_t num_rows;
} sparse_dok;

// 'Constructor' function. Returns pointer to new malloced memory, with 
// appropriate initilisation.
sparse_dok * make_dok();

// --------------- Hashtable Functions -----------------

// find_dok_entry - Find pointer to hash table member with given key, 
// return NULL if member does not exist. 
edge_t *find_dok_entry(sparse_dok * edgetable,edge_key_t key);

// add_dok_entry - Add new entry to the hash table with given key, holding
// the specified value. If entry already exists value is added to current
// value. If entry already exists and this function causes the value to become
// zero, the entry is removed.
void add_dok_entry(sparse_dok * edgetable,edge_key_t key, double value);

// delete_dok_entry - Remove an edge from the hash table. Implicitly sets
// the corresponding entry in the matrix to zero.
void delete_dok_entry(sparse_dok * edgetable,edge_t *edge);

// delete_dok_all - Remove all edges from the hash table. Used to do clean up
void delete_dok_all(sparse_dok * edgetable);

// delete_dok_matrix - Free all the memory associated with struct and
// set pointer to Null.
void delete_dok_matrix(sparse_dok * mat);

// print_dok_entries - Print out all of the stored entries in the hash 
// table sorted by their key, along with their value.
void print_dok_entries(sparse_dok * edgetable);

// key_sort - Compare the relative size of two keys, used for sorting
// PADARN NOTE: Does not need to be in header.
int64_t key_sort(edge_t *a, edge_t *b);

// sort_by_key - Sort the linked list of the hash table by their key
// values and the key_sort function.
void sort_by_key(sparse_dok * hashtable);

// --------------- Matrix Functions ---------------------

// convert_to_csr_ptr - Convert the DOK format matrix into CSR format. The
// new matrix is stored in the (already allocated) input pointer. The old 
// pointer is not freed, and should be cleaned up manually.
void convert_to_csr_ptr(sparse_csr * new_csr,sparse_dok * hashtable);

// add_sparse_dok - Perform a linear addition on two dok_matricies A=a*A+b*B.
// The result is stored in the first input sparse_dok. The second sparse_dok is 
// not freed and should be cleaned up manually. 
// Note: No size conditions are imposed on the matricies, so potentially different
// sized matricies can be combined.
void add_sparse_dok(sparse_dok * dok1,double mult1,sparse_dok * dok2,double mult2);

// get_dok_rows -- Return the number of rows currently stored in the matrix
int64_t get_dok_rows(sparse_dok * dok);

#endif



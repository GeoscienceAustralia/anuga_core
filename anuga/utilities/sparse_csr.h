/*
* Sparse Matrix class implement in CSR format. Limited functionality.
*
* The CSR (compressed sparse row) format is implemented by using three
* arrays, storing the data, column indicies and row pointers.
* This format is used for fast matrix-matrix and matrix-vector 
* multiplication
*
* This function is only currently used to store the result of a 
* conversion from a DOK sparse matrix. Further functionality planned.
*
* Padarn Wilson, ANU 2012
*/

#include <stdio.h>   /* gets */
#include <stdlib.h>  /* atoi, malloc */
#include <string.h>  /* strcpy */
#include <stdint.h>  /* int64_t uint64_t */
#include "math.h"

#ifndef SPARSE_CSR_H
#define SPARSE_CSR_H

// Struct sparse_csr, stores the three matricies which represent
// the sparse matrix along with aiddtional information about the 
// number of rows and the number of entries in the matrix.
typedef struct {
	double *data;
	int64_t *colind;
	int64_t *row_ptr;
	int64_t num_rows;
	int64_t num_entries;
} sparse_csr;

// 'Constructor' function. Returns a pointer to new malloced memory
// All struct entries are intialised appropriately (mostly to NULL). 
sparse_csr * make_csr();

// delete_csr_contents - Free the memory associated with the struct
// and set the pointer to NULL
void delete_csr_matrix(sparse_csr * mat);

#endif


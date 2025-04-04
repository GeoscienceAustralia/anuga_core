#include "sparse_csr.h"

// **************** UTILITIES ***********************

static void *emalloc(size_t amt, char * location)
{
    void *v = malloc(amt);  
    if(!v){
        fprintf(stderr, "out of mem in quad_tree: %s\n",location);
        exit(EXIT_FAILURE);
    }
    return v;
};

// ***************************************************

// 'Constructor'
sparse_csr * make_csr(){

    sparse_csr * ret = emalloc(sizeof(sparse_csr),"make_csr");
    ret->data=NULL;
	ret->colind=NULL;
	ret->row_ptr=NULL;
	ret->num_rows=0;
	ret->num_entries=0;
    return ret;
}

void delete_csr_matrix(sparse_csr * mat){

	free(mat->data);
	free(mat->colind);
	free(mat->row_ptr);
	free(mat);
	mat=NULL;

}


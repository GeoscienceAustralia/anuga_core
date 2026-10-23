#include "sparse_csr.h"

// **************** UTILITIES ***********************

static void *emalloc(size_t amt, char * location)
{
    void *v = malloc(amt);  
    if(!v){
        /* Report and return NULL; callers propagate the failure up to the
           Cython layer, which raises MemoryError. exit() here used to kill
           the whole interpreter and every MPI rank. */
        fprintf(stderr, "out of mem in %s: %s\n", __FILE__, location);
        return NULL;
    }
    return v;
};

// ***************************************************

// 'Constructor'
sparse_csr * make_csr(){

    sparse_csr * ret = emalloc(sizeof(sparse_csr),"make_csr");
    if (!ret) return NULL;
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


#include <inttypes.h> /* PRId64 */
#include "sparse_dok.h"


// **************** UTILITIES ***********************

static void *emalloc(size_t amt,char * location)
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
}

// ***************************************************

// 'Constructor'

sparse_dok * make_dok(void){

    sparse_dok * ret = emalloc(sizeof(sparse_dok),"make_dok");
    if (!ret) return NULL;
    ret->edgetable=NULL;
    ret->num_entries=0;
    ret->num_rows=0;
    return ret;
}

//-----------------------------------------------
//            DOK HASH FUNCTIONS
//-----------------------------------------------

edge_t *find_dok_entry(sparse_dok * hashtable,edge_key_t key) {
    edge_t *s;
    HASH_FIND(hh, hashtable->edgetable, &key, sizeof(edge_key_t), s);  /* s: output pointer */
    return s;
}

anuga_int add_dok_entry(sparse_dok * hashtable, edge_key_t key, double value) {

    // Checking here now if there is an existing value
    // not sure if this code will work.
    
    edge_t *s;
    s = find_dok_entry(hashtable,key);
    if (s) {
        s->entry+=value;
    } else {
        s = (edge_t*) emalloc(sizeof(edge_t),"add_dok_entry");
        if (!s) return -1;
        hashtable->num_entries+=1;
        if(hashtable->num_rows<key.i){
            hashtable->num_rows = key.i;
        }
        memset(s, 0, sizeof(edge_t));
        s->key.i = key.i;
        s->key.j = key.j;
        s->entry = value;
        HASH_ADD(hh, hashtable->edgetable, key, sizeof(edge_key_t), s);  /* key: name of key field */    
    }
    if(s->entry <= 1e-10 && s->entry >= -1e-10){
            HASH_DEL(hashtable->edgetable,s);
            free(s);
            hashtable->num_entries-=1;
    }

    return 0;
}

void delete_dok_entry(sparse_dok * hashtable,edge_t *edge) {
    HASH_DEL( hashtable->edgetable, edge);  /* user: pointer to deletee */
    free(edge);
}

void delete_dok_all(sparse_dok * hashtable) {
  edge_t *current_edge, *tmp;

  HASH_ITER(hh, hashtable->edgetable, current_edge, tmp) {
    HASH_DEL(hashtable->edgetable, current_edge);  /* delete it (hashtable advances to next) */
    free(current_edge);            /* free it */
  } 
}

void delete_dok_matrix(sparse_dok * mat) {

    delete_dok_all(mat);
    free(mat->edgetable);
    free(mat);
    mat=NULL;

}

void print_dok_entries(sparse_dok * hashtable) {
    edge_t *s;

    for(s=hashtable->edgetable; s != NULL; s=(edge_t*)(s->hh.next)) {
        printf("edge key i %" PRId64 " i %" PRId64 " entry %f\n",
                      s->key.i, s->key.j, s->entry);
    }
}

anuga_int key_sort(edge_t *a, edge_t *b) {
    return (a->key.i - b->key.i);
}

anuga_int key_sort_2(edge_t *a, edge_t *b){
    if(a->key.i - b->key.i==0){
        return (a->key.j-b->key.j);
    } else{
        return (a->key.i-b->key.i);
    }
}

void sort_by_key(sparse_dok * hashtable) {     
    HASH_SORT(hashtable->edgetable, key_sort_2);
}

//----------------END--------------------

//-----------------------------------------------
//            DOK MATRIX FUNCTIONS
//-----------------------------------------------


anuga_int convert_to_csr_ptr(sparse_csr * new_csr, sparse_dok * hashtable){


    sparse_csr * ret_csr = new_csr;

    //entrires stores in edgetable -> end

    //sort and get number of entries
    sort_by_key(hashtable); 
    anuga_int num_entries = hashtable->num_entries;
    anuga_int num_rows = hashtable->num_rows+1;

    //build storage matricies
    ret_csr->data=emalloc(num_entries*sizeof(double),"convert_to_csr_ptr");
    ret_csr->colind=emalloc(num_entries*sizeof(anuga_int),"convert_to_csr_ptr");
    ret_csr->row_ptr=emalloc((num_rows+1)*sizeof(anuga_int),"convert_to_csr_ptr");
    if (!ret_csr->data || !ret_csr->colind || !ret_csr->row_ptr) {
        free(ret_csr->data);    ret_csr->data = NULL;
        free(ret_csr->colind);  ret_csr->colind = NULL;
        free(ret_csr->row_ptr); ret_csr->row_ptr = NULL;
        ret_csr->num_rows = 0;
        ret_csr->num_entries = 0;
        return -1;
    }

    edge_t * edge = hashtable->edgetable;

    //now convert
    anuga_int current_row = -1;
    anuga_int k;
    for(k=0;k<num_entries;k++){
        anuga_int i = edge->key.i;
        anuga_int j = edge->key.j;
        double value = edge->entry;

        if (i!=current_row){
            current_row=i;
            ret_csr->row_ptr[i]=k;
        }

        ret_csr->data[k] = value;
        ret_csr->colind[k] = j;
        edge = edge->hh.next;
    }

    for(k=current_row+1;k<num_rows+1;k++){
        ret_csr->row_ptr[k]=num_entries;
    }

    ret_csr -> num_rows = num_rows+1;
    ret_csr -> num_entries = num_entries;

    return 0;
}

anuga_int add_sparse_dok(sparse_dok * dok1,double mult1,sparse_dok * dok2,double mult2){

    // add both into dok1 - then leave both alone (free outside)
    anuga_int num_entries = dok1->num_entries;
    edge_t * edge = dok1->edgetable;
    edge_t * edge2;

    anuga_int k;
    for(k=0;k<num_entries;k++){
        double value = edge->entry;
        edge->entry=value*mult1;
        edge2=find_dok_entry(dok2,edge->key);
        if(edge2!=NULL){
            edge->entry+=edge2->entry*mult2;
            edge2->entry=0;
            //delete_dok_entry(dok2,edge2);
        }
        edge = edge->hh.next;
    }

    num_entries = dok2->num_entries;
    edge = dok2->edgetable;
    for(k=0;k<num_entries;k++){
        if (add_dok_entry(dok1,edge->key,edge->entry*mult2) != 0) return -1;
        edge = edge->hh.next;
    }

    return 0;
}

anuga_int get_dok_rows(sparse_dok * dok){

    anuga_int rows = 0;

    edge_t *current_edge, *tmp;

    HASH_ITER(hh, dok->edgetable, current_edge, tmp) {
        if (current_edge->key.i>rows) rows = current_edge->key.i;
    } 


    return rows;
}


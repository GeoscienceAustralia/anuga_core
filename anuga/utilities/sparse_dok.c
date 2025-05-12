#include "sparse_dok.h"


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

// ***************************************************

// 'Constructor'

sparse_dok * make_dok(){

    sparse_dok * ret = emalloc(sizeof(sparse_dok),"make_dok");
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

void add_dok_entry(sparse_dok * hashtable, edge_key_t key, double value) {

    // Checking here now if there is an existing value
    // not sure if this code will work.
    
    edge_t *s;
    s = find_dok_entry(hashtable,key);
    if (s) {
        s->entry+=value;
    } else {
        hashtable->num_entries+=1;
        if(hashtable->num_rows<key.i){
            hashtable->num_rows = key.i;
        }
        s = (edge_t*) emalloc(sizeof(edge_t),"add_dok_entry");
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
        printf("edge key i %ld i %ld entry %f\n",
                      s->key.i, s->key.j, s->entry);
    }
}

int64_t key_sort(edge_t *a, edge_t *b) {
    return (a->key.i - b->key.i);
}

int64_t key_sort_2(edge_t *a, edge_t *b){
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


void convert_to_csr_ptr(sparse_csr * new_csr, sparse_dok * hashtable){


    sparse_csr * ret_csr = new_csr;

    //entrires stores in edgetable -> end

    //sort and get number of entries
    sort_by_key(hashtable); 
    int64_t num_entries = hashtable->num_entries;
    int64_t num_rows = hashtable->num_rows+1;

    //build storage matricies
    ret_csr->data=emalloc(num_entries*sizeof(double),"convert_to_csr_ptr");
    ret_csr->colind=emalloc(num_entries*sizeof(int64_t),"convert_to_csr_ptr");
    ret_csr->row_ptr=emalloc((num_rows+1)*sizeof(int64_t),"convert_to_csr_ptr");

    edge_t * edge = hashtable->edgetable;

    //now convert
    int64_t current_row = -1;
    int64_t k;
    for(k=0;k<num_entries;k++){
        int64_t i = edge->key.i;
        int64_t j = edge->key.j;
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

}

void add_sparse_dok(sparse_dok * dok1,double mult1,sparse_dok * dok2,double mult2){

    // add both into dok1 - then leave both alone (free outside)
    int64_t num_entries = dok1->num_entries;
    edge_t * edge = dok1->edgetable;
    edge_t * edge2;

    int64_t k;
    for(k=0;k<num_entries;k++){
        int64_t i = edge->key.i;
        int64_t j = edge->key.j;
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
        add_dok_entry(dok1,edge->key,edge->entry*mult2);
        edge = edge->hh.next;
    }

}

int64_t get_dok_rows(sparse_dok * dok){

    int64_t rows = 0;

    edge_t *current_edge, *tmp;

    HASH_ITER(hh, dok->edgetable, current_edge, tmp) {
        if (current_edge->key.i>rows) rows = current_edge->key.i;
    } 


    return rows;
}


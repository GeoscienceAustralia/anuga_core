#include <stdio.h>   /* gets */
#include <stdlib.h>  /* atoi, malloc */
#include <string.h>  /* strcpy */
#include <inttypes.h> /* PRId64 */
#include <math.h>

//Shared code snippets

#include "uthash.h"     /* in utilities */
#include "anuga_typedefs.h" /* in utilities */

//==============================================================================
// hashtable code from uthash. Look at copyright info in "uthash.h in the
// utilities directory
//==============================================================================

typedef struct {
    anuga_int i;
    anuga_int j;
} segment_key_t;

typedef struct {
    segment_key_t key; /* key of form i , j */
    anuga_int vol_id; /* id of vol containing this segement */
    anuga_int edge_id; /* edge_id of segement in this vol */
    UT_hash_handle hh; /* makes this structure hashable */
} segment_t;

segment_t *segment_table = NULL;

/* Returns 0, or -1 if the entry could not be allocated. */
anuga_int add_segment(segment_key_t key, anuga_int vol_id, anuga_int edge_id) {
    segment_t *s;

    s = (segment_t*) malloc(sizeof (segment_t));
    if (s == NULL) return -1;
    memset(s, 0, sizeof (segment_t));
    s->key.i = key.i;
    s->key.j = key.j;
    s->vol_id = vol_id;
    s->edge_id = edge_id;
    HASH_ADD(hh, segment_table, key, sizeof (segment_key_t), s); /* key: name of key field */
    return 0;
}

segment_t *find_segment(segment_key_t key) {
    segment_t *s=0;

    HASH_FIND(hh, segment_table, &key, sizeof (segment_key_t), s); /* s: output pointer */
    return s;
}

void delete_segment(segment_t *segment) {
    HASH_DEL(segment_table, segment); /* user: pointer to deletee */
    free(segment);
}

void delete_segment_all(void) { /* note we need to use void here to suppress warning */
    segment_t *current_segment, *tmp;

    HASH_ITER(hh, segment_table, current_segment, tmp) {
        HASH_DEL(segment_table, current_segment); /* delete it (segment_table advances to next) */
        free(current_segment); /* free it */
    }
}

void print_segments(void) {
    segment_t *s;

    for (s = segment_table; s != NULL; s = (segment_t*) (s->hh.next)) {
        printf("segment key i %" PRId64 " j %" PRId64 " vol_id %" PRId64 "  edge_id %" PRId64 "\n",
                s->key.i, s->key.j, s->vol_id, s->edge_id);
    }
}

anuga_int vol_id_sort(segment_t *a, segment_t *b) {
    return (a->vol_id - b->vol_id);
}

anuga_int key_sort(segment_t *a, segment_t *b) {
    return (a->key.i - b->key.i);
}

void sort_by_vol_id(void) {
    HASH_SORT(segment_table, vol_id_sort);
}

void sort_by_key(void) {
    HASH_SORT(segment_table, key_sort);
}

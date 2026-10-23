#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "anuga_typedefs.h"

#include "omp.h"
/* Ascending comparator for qsort on anuga_int arrays (the four neighbour
   indices sorted in the matrix builders below). */
static int cmp_anuga_int(const void *a, const void *b) {
    anuga_int x = *(const anuga_int *)a, y = *(const anuga_int *)b;
    return (x > y) - (x < y);
}

anuga_int _build_geo_structure(anuga_int n,
        anuga_int tot_len,
        double *centroids,
        anuga_int *neighbours,
        double *edgelengths,
        double *edge_midpoints,
        anuga_int *geo_indices,
        double *geo_values) {
    anuga_int i, edge, j, m;
    double dist, this_x, this_y, other_x, other_y, edge_length;
    #pragma omp parallel for private(edge, j, m, dist, this_x, this_y, other_x, other_y, edge_length)
    for (i = 0; i < n; i++) {
        //The centroid coordinates of triangle i
        this_x = centroids[2 * i];
        this_y = centroids[2 * i + 1];
        for (edge = 0; edge < 3; edge++) {

            j = neighbours[3 * i + edge];

            //Get the index and the coordinates of the interacting point

            // Edge
            if (j < 0) {
                m = -j - 1;
                geo_indices[3 * i + edge] = n + m;

                other_x = edge_midpoints[2 * (3 * i + edge)];
                other_y = edge_midpoints[2 * (3 * i + edge) + 1];
            } else {
                geo_indices[3 * i + edge] = j;

                other_x = centroids[2 * j];
                other_y = centroids[2 * j + 1];
            }

            //Compute the interaction
            edge_length = edgelengths[3 * i + edge];
            dist = sqrt((this_x - other_x)*(this_x - other_x) + (this_y - other_y)*(this_y - other_y));
            geo_values[3 * i + edge] = -edge_length / dist;
        }
    }
    return 0;
}

anuga_int _build_elliptic_matrix_not_symmetric(anuga_int n,
        anuga_int tot_len,
        anuga_int *geo_indices,
        double *geo_values,
        double *cell_data,
        double *bdry_data,
        double *data,
        anuga_int *colind) {
    anuga_int i, k, edge, j[4], sorted_j[4], this_index;
    double v[3], v_i; //v[k] = value of the interaction of edge k in a given triangle, v_i = (i,i) entry
    for (i = 0; i < n; i++) {
        v_i = 0.0;
        j[3] = i;
        //Get the values of each interaction, and the column index at which they occur
        for (edge = 0; edge < 3; edge++) {
            j[edge] = geo_indices[3 * i + edge];
            //if (j[edge]<n) { //interior
            //    h_j = cell_data[j[edge]];
            //} else { //boundary
            //    h_j = bdry_data[j[edge]-n];
            //}
            v[edge] = -cell_data[i] * geo_values[3 * i + edge]; //the negative of the individual interaction
            v_i += cell_data[i] * geo_values[3 * i + edge]; //sum the three interactions
        }
        if (cell_data[i] <= 0.0) {
            v_i = 0.0;
            v[0] = 0.0;
            v[1] = 0.0;
            v[2] = 0.0;
        }
        //Organise the set of 4 values/indices into the data and colind arrays
        for (k = 0; k < 4; k++) sorted_j[k] = j[k];
        qsort(sorted_j, 4, sizeof(anuga_int), cmp_anuga_int);
        for (k = 0; k < 4; k++) { //loop through the nonzero indices
            this_index = sorted_j[k];
            if (this_index == i) {
                data[4 * i + k] = v_i;
                colind[4 * i + k] = i;
            } else if (this_index == j[0]) {
                data[4 * i + k] = v[0];
                colind[4 * i + k] = j[0];
            } else if (this_index == j[1]) {
                data[4 * i + k] = v[1];
                colind[4 * i + k] = j[1];
            } else { //this_index == j[2]
                data[4 * i + k] = v[2];
                colind[4 * i + k] = j[2];
            }
        }
    }
    return 0;
}

anuga_int _build_elliptic_matrix(anuga_int n,
        anuga_int tot_len,
        anuga_int *geo_indices,
        double *geo_values,
        double *cell_data,
        double *bdry_data,
        double *data,
        anuga_int *colind) {
    anuga_int i, k, edge, j[4], sorted_j[4], this_index;
    double h_j, v[3], v_i; //v[k] = value of the interaction of edge k in a given triangle, v_i = (i,i) entry
    // Each row i writes only to data[4*i:4*i+4] and colind[4*i:4*i+4] — no race.
    #pragma omp parallel for private(k, edge, j, sorted_j, this_index, h_j, v, v_i)
    for (i = 0; i < n; i++) {
        v_i = 0.0;
        j[3] = i;
        //Get the values of each interaction, and the column index at which they occur
        for (edge = 0; edge < 3; edge++) {
            j[edge] = geo_indices[3 * i + edge];
            if (j[edge] < n) { //interior
                h_j = cell_data[j[edge]];
            } else { //boundary
                h_j = bdry_data[j[edge] - n];
            }
            v[edge] = -0.5 * (cell_data[i] + h_j) * geo_values[3 * i + edge]; //the negative of the individual interaction
            v_i += 0.5 * (cell_data[i] + h_j) * geo_values[3 * i + edge]; //sum the three interactions
        }
        if (cell_data[i] <= 0.0) {
            v_i = 0.0;
            v[0] = 0.0;
            v[1] = 0.0;
            v[2] = 0.0;
        }
        //Organise the set of 4 values/indices into the data and colind arrays
        for (k = 0; k < 4; k++) sorted_j[k] = j[k];
        qsort(sorted_j, 4, sizeof(anuga_int), cmp_anuga_int);
        for (k = 0; k < 4; k++) { //loop through the nonzero indices
            this_index = sorted_j[k];
            if (this_index == i) {
                data[4 * i + k] = v_i;
                colind[4 * i + k] = i;
            } else if (this_index == j[0]) {
                data[4 * i + k] = v[0];
                colind[4 * i + k] = j[0];
            } else if (this_index == j[1]) {
                data[4 * i + k] = v[1];
                colind[4 * i + k] = j[1];
            } else { //this_index == j[2]
                data[4 * i + k] = v[2];
                colind[4 * i + k] = j[2];
            }
        }
    }
    return 0;
}

anuga_int _update_elliptic_matrix_not_symmetric(anuga_int n,
        anuga_int tot_len,
        anuga_int *geo_indices,
        double *geo_values,
        double *cell_data,
        double *bdry_data,
        double *data,
        anuga_int *colind) {
    anuga_int i, k, edge, j[4], sorted_j[4], this_index;
    double  v[3], v_i; //v[k] = value of the interaction of edge k in a given triangle, v_i = (i,i) entry
    for (i = 0; i < n; i++) {
        v_i = 0.0;
        j[3] = i;

        //Get the values of each interaction, and the column index at which they occur
        for (edge = 0; edge < 3; edge++) {
            j[edge] = geo_indices[3 * i + edge];
            v[edge] = -cell_data[i] * geo_values[3 * i + edge]; //the negative of the individual interaction
            v_i += cell_data[i] * geo_values[3 * i + edge]; //sum the three interactions
        }
        if (cell_data[i] <= 0.0) {
            v_i = 0.0;
            v[0] = 0.0;
            v[1] = 0.0;
            v[2] = 0.0;
        }
        //Organise the set of 4 values/indices into the data and colind arrays
        for (k = 0; k < 4; k++) sorted_j[k] = j[k];
        qsort(sorted_j, 4, sizeof(anuga_int), cmp_anuga_int);
        for (k = 0; k < 4; k++) { //loop through the nonzero indices
            this_index = sorted_j[k];
            if (this_index == i) {
                data[4 * i + k] = v_i;
                colind[4 * i + k] = i;
            } else if (this_index == j[0]) {
                data[4 * i + k] = v[0];
                colind[4 * i + k] = j[0];
            } else if (this_index == j[1]) {
                data[4 * i + k] = v[1];
                colind[4 * i + k] = j[1];
            } else { //this_index == j[2]
                data[4 * i + k] = v[2];
                colind[4 * i + k] = j[2];
            }
        }
    }
    return 0;
}

anuga_int _update_elliptic_matrix(anuga_int n,
        anuga_int tot_len,
        anuga_int *geo_indices,
        double *geo_values,
        double *cell_data,
        double *bdry_data,
        double *data,
        anuga_int *colind) {
    anuga_int i, k, edge, j[4], sorted_j[4], this_index;
    double h_j, v[3], v_i; //v[k] = value of the interaction of edge k in a given triangle, v_i = (i,i) entry
    // Each row i writes only to data[4*i:4*i+4] and colind[4*i:4*i+4] — no race.
    #pragma omp parallel for private(k, edge, j, sorted_j, this_index, h_j, v, v_i)
    for (i = 0; i < n; i++) {
        v_i = 0.0;
        j[3] = i;

        //Get the values of each interaction, and the column index at which they occur
        for (edge = 0; edge < 3; edge++) {
            j[edge] = geo_indices[3 * i + edge];
            if (j[edge] < n) { //interior
                h_j = cell_data[j[edge]];
            } else { //boundary
                h_j = bdry_data[j[edge] - n];
            }
            v[edge] = -0.5 * (cell_data[i] + h_j) * geo_values[3 * i + edge]; //the negative of the individual interaction
            v_i += 0.5 * (cell_data[i] + h_j) * geo_values[3 * i + edge]; //sum the three interactions
        }
        if (cell_data[i] <= 0.0) {
            v_i = 0.0;
            v[0] = 0.0;
            v[1] = 0.0;
            v[2] = 0.0;
        }
        //Organise the set of 4 values/indices into the data and colind arrays
        for (k = 0; k < 4; k++) sorted_j[k] = j[k];
        qsort(sorted_j, 4, sizeof(anuga_int), cmp_anuga_int);
        for (k = 0; k < 4; k++) { //loop through the nonzero indices
            this_index = sorted_j[k];
            if (this_index == i) {
                data[4 * i + k] = v_i;
                colind[4 * i + k] = i;
            } else if (this_index == j[0]) {
                data[4 * i + k] = v[0];
                colind[4 * i + k] = j[0];
            } else if (this_index == j[1]) {
                data[4 * i + k] = v[1];
                colind[4 * i + k] = j[1];
            } else { //this_index == j[2]
                data[4 * i + k] = v[2];
                colind[4 * i + k] = j[2];
            }
        }
    }
    return 0;
}

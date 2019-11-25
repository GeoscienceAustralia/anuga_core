#include <math.h>
#include <stdio.h>
#include <unistd.h>

//Rough quicksort implementation (for build_operator_matrix)
// taken from http://cprogramminglanguage.net/quicksort-algorithm-c-source-code.aspx

void swap(int *x, int *y) {
    int temp;
    temp = *x;
    *x = *y;
    *y = temp;
}

int choose_pivot(int i, int j) {
    return ((i + j) / 2);
}

void quicksort(int list[], int m, int n) {
    int key, i, j, k;
    if (m < n) {
        k = choose_pivot(m, n);
        swap(&list[m], &list[k]);
        key = list[m];
        i = m + 1;
        j = n;
        while (i <= j) {
            while ((i <= n) && (list[i] <= key))
                i++;
            while ((j >= m) && (list[j] > key))
                j--;
            if (i < j)
                swap(&list[i], &list[j]);
        }
        // swap two elements
        swap(&list[m], &list[j]);
        // recursively sort the lesser list
        quicksort(list, m, j - 1);
        quicksort(list, j + 1, n);
    }
}

int _build_geo_structure(int n,
        int tot_len,
        double *centroids,
        long *neighbours,
        double *edgelengths,
        double *edge_midpoints,
        long *geo_indices,
        double *geo_values) {
    int i, edge, edge_counted, j, m;
    double dist, this_x, this_y, other_x, other_y, edge_length;
    edge_counted = 0;
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
                edge_counted++;

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

int _build_elliptic_matrix_not_symmetric(int n,
        int tot_len,
        long *geo_indices,
        double *geo_values,
        double *cell_data,
        double *bdry_data,
        double *data,
        long *colind) {
    int i, k, edge, j[4], sorted_j[4], this_index;
    double h_j, v[3], v_i; //v[k] = value of the interaction of edge k in a given triangle, v_i = (i,i) entry
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
        quicksort(sorted_j, 0, 3);
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

int _build_elliptic_matrix(int n,
        int tot_len,
        long *geo_indices,
        double *geo_values,
        double *cell_data,
        double *bdry_data,
        double *data,
        long *colind) {
    int i, k, edge, j[4], sorted_j[4], this_index;
    double h_j, v[3], v_i; //v[k] = value of the interaction of edge k in a given triangle, v_i = (i,i) entry
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
        quicksort(sorted_j, 0, 3);
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

int _update_elliptic_matrix_not_symmetric(int n,
        int tot_len,
        long *geo_indices,
        double *geo_values,
        double *cell_data,
        double *bdry_data,
        double *data,
        long *colind) {
    int i, k, edge, j[4], sorted_j[4], this_index;
    double h_j, v[3], v_i; //v[k] = value of the interaction of edge k in a given triangle, v_i = (i,i) entry
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
        quicksort(sorted_j, 0, 3);
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

int _update_elliptic_matrix(int n,
        int tot_len,
        long *geo_indices,
        double *geo_values,
        double *cell_data,
        double *bdry_data,
        double *data,
        long *colind) {
    int i, k, edge, j[4], sorted_j[4], this_index;
    double h_j, v[3], v_i; //v[k] = value of the interaction of edge k in a given triangle, v_i = (i,i) entry
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
        quicksort(sorted_j, 0, 3);
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

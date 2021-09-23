// Tsunameter functions
#include <stdlib.h>
#include <time.h>
#include "tsunameter.h"

moving_avg* init_moving_avg(int max_size) {
    moving_avg* avg = (moving_avg *) malloc(sizeof(moving_avg));
    avg -> max_size = max_size;
    return avg;
}

void append_moving_avg(moving_avg* avg, float val) {
    avg_el* new_el = (avg_el *) malloc(sizeof(avg_el));
    new_el -> value = val;
    if (avg -> size == avg -> max_size) {
        avg_el* temp = avg -> queue_head;
        avg -> queue_head = temp -> next;
        avg -> queue_tail -> next = new_el;
        avg -> avg = (float)(avg -> avg * avg -> size - temp->value + new_el -> value) / avg -> size;
        free(temp);
    } else if (avg -> size > 0) {
        avg -> queue_tail -> next = new_el;
        avg -> size += 1;
        avg -> avg = (float)(avg -> avg * (avg -> size - 1) + new_el -> value) / avg -> size;
    } else {
        avg -> queue_head = new_el;
        avg -> avg = new_el -> value;
        avg -> size = 1;
    }
    avg -> queue_tail = new_el;
}

float get_moving_avg(moving_avg* avg) {
    if (avg -> size == 0) {
        return -1.0;
    }
    return avg -> avg;
}

float generate_float_val(float limit) {
    return (float)rand()/(float)(RAND_MAX) * limit;
}

void get_rank_at_coord(int* dims, int n_dims, int* coord, int* rank) {
    *rank = 1;
    int dim_index;
    // TODO: Wrong for dim > 2 but w/e
    for (dim_index = 0; dim_index < n_dims - 1; dim_index++) {
        *rank *= dims[dim_index + 1] * coord[dim_index];
    }
    *rank += coord[n_dims - 1];
}

void get_coord_at_rank(int* dims, int n_dims, int rank, int* coord) {
    // TODO: Wrong for dim > 2 but w/e
    int dim_index;
    for (dim_index = 0; dim_index < n_dims - 1; dim_index++) {
        coord[dim_index] = (int)rank/dims[dim_index];
    }
    coord[n_dims - 1] = rank - coord[0] * dims[1];
}

bool coord_exists(int* dims, int n_dims, int* coord) {
    int exists = true;
    int dim_index;
    for (dim_index = 0; dim_index < n_dims; dim_index++) {
        exists = exists && coord[dim_index] >= 0 && coord[dim_index] < dims[dim_index];
    }
    return exists;
}

int test_mpi_req(MPI_Request *request, int *flag, MPI_Status *status) {
    MPI_Test(request, flag, status);
    return *flag;
}
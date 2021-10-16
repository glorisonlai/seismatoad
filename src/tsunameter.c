// Tsunameter functions
#include "tsunameter.h"

#include <stdlib.h>

moving_avg *init_moving_avg(int max_size) {
    moving_avg *avg = (moving_avg *)malloc(sizeof(moving_avg));
    avg->max_size = max_size;
    return avg;
}

void append_moving_avg(moving_avg *avg, float val) {
    avg_el *new_el = (avg_el *)malloc(sizeof(avg_el));
    new_el->value = val;
    if (avg->size == avg->max_size) {
        avg_el *temp = avg->queue_head;
        avg->queue_head = temp->next;
        avg->queue_tail->next = new_el;
        avg->avg = (float)(avg->avg * avg->size - temp->value + new_el->value) /
                   avg->size;
        free(temp);
    } else if (avg->size > 0) {
        avg->queue_tail->next = new_el;
        avg->size += 1;
        avg->avg =
            (float)(avg->avg * (avg->size - 1) + new_el->value) / avg->size;
    } else {
        avg->queue_head = new_el;
        avg->avg = new_el->value;
        avg->size = 1;
    }
    avg->queue_tail = new_el;
}

float get_moving_avg(moving_avg *avg) {
    if (avg->size == 0) {
        return -1.0;
    }
    return avg->avg;
}

float generate_float_val(float limit) {
    return (float)rand() / (float)(RAND_MAX)*limit;
}

void get_rank_at_coord(int *dims, int n_dims, int *coord, int *rank) {
    *rank = 1;
    int dim_index;
    // TODO: Wrong for dim > 2 but w/e
    for (dim_index = 0; dim_index < n_dims - 1; dim_index++) {
        *rank *= dims[dim_index + 1] * coord[dim_index];
    }
    *rank += coord[n_dims - 1];
}

void get_coord_at_rank(int *dims, int n_dims, int rank, int *coord) {
    // TODO: Wrong for dim > 2 but w/e
    int dim_index;
    for (dim_index = 0; dim_index < n_dims - 1; dim_index++) {
        coord[dim_index] = (int)rank / dims[dim_index];
    }
    coord[n_dims - 1] = rank - coord[0] * dims[1];
}

int test_mpi_req(MPI_Request *request) {
    int flag;
    MPI_Status status;
    MPI_Test(request, &flag, &status);
    return flag;
}

int* get_neighbours(MPI_Comm comm, int ndims, int *num_neighbours) {
    *num_neighbours = 0;
    int pot_neighbours[2 * ndims];
    int dir;
    for (dir = 0; dir < ndims; dir++) {
        MPI_Cart_shift(comm, dir, 1, &pot_neighbours[2 * dir], &pot_neighbours[2 * dir + 1]);
        int i;
        for (i = 0; i < 2; i++) {
            if (pot_neighbours[2 *dir + i] >= 0) {
                *num_neighbours += 1;
            }
        }
    }

    int *neighbours = (int *) malloc(*num_neighbours);
    int i, nbr_ptr = 0;
    for (i = 0; i < 2 * ndims; i++) {
        if (pot_neighbours[i] >= 0) {
            neighbours[nbr_ptr] = pot_neighbours[i];
            nbr_ptr += 1;
        }
    }
    return neighbours;
}

tsunameter_reading *instantiate_tsunameter_reading(float avg, time_t time) {
    tsunameter_reading *reading = malloc(sizeof(*reading));
    reading->avg = avg;
    reading->time = (int)time;
    return reading;
}

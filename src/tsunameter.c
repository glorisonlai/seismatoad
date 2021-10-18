// Tsunameter functions
#include "tsunameter.h"

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/**
 *  Instantiate moving average data structure onto heap
 *  @param int max_size Size of moving average window
 *  @return Pointer to moving average
 */
moving_avg *init_moving_avg(int max_size) {
    moving_avg *avg = (moving_avg *)malloc(sizeof(moving_avg));
    avg->max_size = max_size;
    avg->avg = 0;
    avg->size = 0;
    return avg;
}

/**
 *  Append new element to moving average, and update average
 *  @param moving_avg* avg Pointer to moving average
 *  @param float val Value of new element to push
 *  @return None
 */
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

/**
 *  Get most recent moving average
 *  @param moving_avg* avg Pointer to moving average
 *  @return Current moving average
 */
float get_moving_avg(moving_avg *avg) {
    if (avg->size == 0) {
        return -1.0;
    }
    return avg->avg;
}

/**
 *  Frees moving average and its elements from heap memory
 *  @param moving_avg* avg Pointer to moving average
 *  @return None
 */ 
void free_moving_avg(moving_avg *avg) {
    while (avg->size) {
        avg_el *temp = avg->queue_head;
        avg->queue_head = temp->next;
        free(temp);
        avg->size -= 1;
    }
    free(avg);
}

/**
 *  Generates new random float value
 *  @param float limit Maximum float value to reach
 *  @return Random float value
 */
float generate_float_val(float limit) {
    return (float)rand() / (float)(RAND_MAX)*limit;
}

/**
 *  Tests whether MPI Request has completed
 *  @param MPI_Request* request Pointer to test request
 *  @return Boolean if completed
 */
int test_mpi_req(MPI_Request *request) {
    int flag;
    MPI_Status status;
    MPI_Test(request, &flag, &status);
    return flag;
}

/**
 *  Finds ranks of neighbours immediately adjacent to node
 *  @param MPI_Comm Communicator of cartesian group
 *  @param int Number of dimensions in cartesian group
 *  @return Pointer to array of neighbours in heap
 *  @return int* num_neighbours gets updated to number of valid neighbours
 */
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

/**
 *  Instantiate new tsunameter reading to send to neighbours
 *  @param float avg Value of average wave heights
 *  @param time_t Time of reading
 *  @return Pointer to tsunameter reading
 */
tsunameter_reading *instantiate_tsunameter_reading(float avg, time_t time) {
    tsunameter_reading *reading = (tsunameter_reading *)malloc(sizeof(tsunameter_reading));
    reading->avg = avg;
    reading->time = (int)time;
    return reading;
}

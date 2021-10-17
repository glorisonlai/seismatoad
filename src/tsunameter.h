// Interface for Tsunameter functions
#ifndef TSUNAMETER_H
#define TSUNAMETER_H

#include <stdbool.h>
#include <mpi.h>
#include <time.h>

typedef struct avg_el
{
    float value;
    struct avg_el *next;
} avg_el;

typedef struct moving_avg
{
    float avg;
    int max_size;
    int size;
    avg_el *queue_head;
    avg_el *queue_tail;
} moving_avg;

typedef struct tsunameter_reading
{
    float avg;
    int time;
} tsunameter_reading;

typedef struct base_station_info
{
    float avg;
    int time;
    int neighbours[4];
} base_station_info;

// Initialize moving average
moving_avg *init_moving_avg(int max_size);

// Append float to moving average and update average
void append_moving_avg(moving_avg *avg, float val);

// Get current average
float get_moving_avg(moving_avg *avg);

// Free moving average struct
void free_moving_avg(moving_avg *avg);

// Generate random float up to limit
float generate_float_val(float limit);

void get_rank_at_coord(int *dims, int n_dims, int *coord, int *rank);

void get_coord_at_rank(int *dims, int n_dims, int rank, int *coord);

bool coord_exists(int *dims, int n_dims, int *coord);

int test_mpi_req(MPI_Request *request);

int* get_neighbours(MPI_Comm comm, int ndims, int *num_neighbours);

struct tsunameter_reading *instantiate_tsunameter_reading(float avg, time_t time);

#endif
// Interface for Tsunameter functions
#ifndef TSUNAMETER_H
#define TSUNAMETER_H

#include <stdbool.h>
#include <mpi.h>
#include <time.h>

// Inner elements of moving average
typedef struct avg_el
{
    float value;
    struct avg_el *next;
} avg_el;

// Struct of moving average
typedef struct moving_avg
{
    float avg;
    int max_size;
    int size;
    avg_el *queue_head;
    avg_el *queue_tail;
} moving_avg;

// Struct of information to send to other tsunameters
typedef struct tsunameter_reading
{
    float avg;
    int time;
} tsunameter_reading;

// Struct of information to send to base station
typedef struct base_station_info
{
    float avg;
    int time;
    int neighbours[4];
    double comm_time;
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

// Test if MPI Request has completed
int test_mpi_req(MPI_Request *request);

// Get adjacent neighbours from rank at ordered cartesian mapping
int* get_neighbours(MPI_Comm comm, int ndims, int *num_neighbours);

// Create new tsunameter reading
struct tsunameter_reading *instantiate_tsunameter_reading(float avg, time_t time);

#endif
// Tsunameter functions
#include <stdlib.h>
#include <time.h>
#include "tsunameter.h"

moving_avg *init_moving_avg(int max_size)
{
    moving_avg *avg = (moving_avg *)malloc(sizeof(moving_avg));
    avg->max_size = max_size;
    return avg;
}

void append_moving_avg(moving_avg *avg, float val)
{
    avg_el *new_el = (avg_el *)malloc(sizeof(avg_el));
    new_el->value = val;
    if (avg->size == avg->max_size)
    {
        avg_el *temp = avg->queue_head;
        avg->queue_head = temp->next;
        avg->queue_tail->next = new_el;
        avg->avg = (float)(avg->avg * avg->size - temp->value + new_el->value) / avg->size;
        free(temp);
    }
    else if (avg->size > 0)
    {
        avg->queue_tail->next = new_el;
        avg->size += 1;
        avg->avg = (float)(avg->avg * (avg->size - 1) + new_el->value) / avg->size;
    }
    else
    {
        avg->queue_head = new_el;
        avg->avg = new_el->value;
        avg->size = 1;
    }
    avg->queue_tail = new_el;
}

float get_moving_avg(moving_avg *avg)
{
    if (avg->size == 0)
    {
        return -1.0;
    }
    return avg->avg;
}

float generate_float_val(float limit)
{
    return (float)rand() / (float)(RAND_MAX)*limit;
}

void get_rank_at_coord(int *dims, int n_dims, int *coord, int *rank)
{
    *rank = 1;
    int dim_index;
    // TODO: Wrong for dim > 2 but w/e
    for (dim_index = 0; dim_index < n_dims - 1; dim_index++)
    {
        *rank *= dims[dim_index + 1] * coord[dim_index];
    }
    *rank += coord[n_dims - 1];
}

void get_coord_at_rank(int *dims, int n_dims, int rank, int *coord)
{
    // TODO: Wrong for dim > 2 but w/e
    int dim_index;
    for (dim_index = 0; dim_index < n_dims - 1; dim_index++)
    {
        coord[dim_index] = (int)rank / dims[dim_index];
    }
    coord[n_dims - 1] = rank - coord[0] * dims[1];
}

bool coord_exists(int *dims, int n_dims, int *coord)
{
    int exists = true;
    int dim_index;
    for (dim_index = 0; dim_index < n_dims; dim_index++)
    {
        exists = exists && coord[dim_index] >= 0 && coord[dim_index] < dims[dim_index];
    }
    return exists;
}

int test_mpi_req(MPI_Request *request, int *flag, MPI_Status *status)
{
    MPI_Test(request, flag, status);
    return *flag;
}

// Hack to get array of neighbours
// Multiply by prime numbers, then check for division
// If true, neighbour exists, and we add it to an arr
void get_neighbours(int rank, int rows, int cols, int *neighbours, int *num_neighbours)
{
    *num_neighbours = 4;
    int check_sum = 1;
    // On top row
    if (rank < cols)
    {
        num_neighbours -= 1;
        check_sum *= 2;
    }
    // On bottom row
    if (rank >= (rows - 1) * cols)
    {
        num_neighbours -= 1;
        check_sum *= 3;
    }
    // On left col
    if (rank % cols == 0)
    {
        num_neighbours -= 1;
        check_sum *= 5;
    }
    // On right col
    if (rank % cols == cols - 1)
    {
        num_neighbours -= 1;
        check_sum *= 7;
    }
    neighbours = (int *)malloc(*num_neighbours * sizeof(int));

    int neighbour_pointer = 0;
    if (check_sum % 2 == 0)
    {
        neighbours[neighbour_pointer] = rank - cols + 1;
        neighbour_pointer += 1;
    }
    if (check_sum % 3 == 0)
    {
        neighbours[neighbour_pointer] = rank + cols + 1;
        neighbour_pointer += 1;
    }
    if (check_sum % 5 == 0)
    {
        neighbours[neighbour_pointer] = rank - 1 + 1;
        neighbour_pointer += 1;
    }
    if (check_sum % 7 == 0)
    {
        neighbours[neighbour_pointer] = rank + 1 + 1;
        neighbour_pointer += 1;
    }
}

tsunameter_reading *instantiate_tsunameter_reading(float avg, time_t time)
{
    tsunameter_reading *reading;
    reading->avg = avg;
    reading->time = (int)time;
    return reading;
}
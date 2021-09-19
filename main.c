#include <stdio.h>  // Printing
#include <mpi.h>    // OpenMPI calls
#include <stdlib.h> // Type conversion
#include <errno.h>  // ERROR CHECKING
#include <limits.h> // INT LIMITS
#include "tsunameter.h" // Tsunameter functions
#include <unistd.h> // Sleep
#include <time.h>   // Seed for random

#define DIMENSIONS 2    // Tsunameter topology dimensions
#define TERMINATION_TAG 0   // MPI Tag for termination
#define THRESHOLD 6000  // Tsunameter threshold
#define TOLERANCE 100   // Tsunameter tolerance
#define TSUNAMTER_WINDOW 10 // Tsunameter average window
#define TSUNAMETER_POLL 10    // Tsunamter polling rate (s)

int main(int argc, char** argv)
{
    int dims[DIMENSIONS];
    int comm_size;
    int comm_rank;
    int termination_flag = 0;
    MPI_Status termination_status;
    int rc;
    srand(time(NULL));

    // Initialize MPI
    rc = MPI_Init(&argc, &argv);
    // Check MPI Init is successful
    if (rc != MPI_SUCCESS) {
        printf("\nCould not initialize MPI\n");fflush(stdout);
        return 1;
    }
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

    // Ensure tsunameter topology is given
    if (argc != DIMENSIONS + 1) {
        printf("\nUsage: mpirun -np <PROCESSES> main m n\n");fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Initialize topology
    errno = 0;
    int dim_index;
    int check_sum = 1;
    for (dim_index = 0; dim_index < DIMENSIONS; dim_index++) {
        char* p;
        long arg = strtol(argv[dim_index + 1], &p, 10); // Get argument index from 1
        if (errno != 0 || *p != '\0' || arg < INT_MIN || arg > INT_MAX) {
            // Handle type conversion error
            printf("\nArgument error\n");fflush(stdout);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // Fill in dimensions
        dims[dim_index] = arg;
        check_sum *= dims[dim_index];
    }

    // Make sure allocated process matches topology
    if (check_sum != comm_size - 1) {
        printf("\nIncompatible topology and tsunameter count\n");fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    printf("Dimensions: %dx%d, Size: %d\n", dims[0], dims[1], comm_size);

    // // Virtual topology
    // TODO: Make virtual top
    // int reorder = 0;
    // int period[DIMENSIONS] = {0};
    // printf("%d%d", period[0], period[1]);
    // rc = MPI_Cart_create(MPI_COMM_WORLD, DIMENSIONS, dims, period, reorder, &comm);
    if (comm_rank == 0) {
        // int blah = 0;
        // MPI_Request send_req;
        // sleep(30);
        // MPI_Isend(&blah, 1, MPI_INT, 1, TERMINATION_TAG, MPI_COMM_WORLD, &send_req);
        // TODO: Base Station logic
    } else {
        printf("rank %d\n", comm_rank);
        moving_avg* avg = init_moving_avg(TSUNAMTER_WINDOW);
        printf("Initialized avg, size: %d/max_size: %d\n", avg->size, avg->max_size);

        while (!termination_flag) {
            MPI_Iprobe(0, TERMINATION_TAG, MPI_COMM_WORLD, &termination_flag, &termination_status);
            printf("looping %d\n", termination_flag);
            sleep(TSUNAMETER_POLL);
            float new_val = generate_float_val(100000.0);
            append_moving_avg(avg, new_val);
            printf("hello, %lf\n", new_val);
            printf("Size: %d, Average: %lf, Head: %lf, Tail: %lf\n\n",avg -> size, avg -> avg, avg -> queue_head -> value, avg -> queue_tail -> value);
        }

        // Termination signal received
        // Free avg? idk
        printf("terminating...\n");
    }

    MPI_Finalize();
}


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
#define TSUNAMETER_POLL 2    // Tsunamter polling rate (s)

int main(int argc, char** argv)
{
    int tsunameter_dims[DIMENSIONS];
    int comm_size;
    int comm_rank;
    MPI_Status termination_status;
    int rc;
    int root = 0;

    // Testing cart create
    MPI_Comm new_comm;

    srand(time(NULL));  // Seed random gen for tsunameters

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
        tsunameter_dims[dim_index] = arg;
        check_sum *= tsunameter_dims[dim_index];
    }

    // Make sure allocated process matches topology
    if (check_sum != comm_size - 1) {
        printf("\nIncompatible topology and tsunameter count\n");fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    printf("Dimensions: %dx%d, Size: %d\n", tsunameter_dims[0], tsunameter_dims[1], comm_size);

    // Virtual topology
    if (comm_rank == root) {
        // int blah = 0;
        // MPI_Request send_req;
        // sleep(10);
        // MPI_Ibcast(&blah, 1, MPI_INT, root, MPI_COMM_WORLD, &send_req);
        // TODO: Base Station logic
    } else {
        // Instantiate moving average 
        printf("rank %d\n", comm_rank);
        moving_avg* avg = init_moving_avg(TSUNAMTER_WINDOW);

        // Instantiate topology
        int termination_flag = 0;
        int coord[DIMENSIONS];
        get_coord_at_rank(tsunameter_dims, DIMENSIONS, comm_rank, coord);

        // Get array of neighbours

        // Send receive for termination signal
        int termination_buffer[1];
        MPI_Request termination_req;
        MPI_Ibcast(termination_buffer, 1, MPI_INT, root, MPI_COMM_WORLD, &termination_req);

        // Check for termination signal
        while (!test_mpi_req(&termination_req, &termination_flag, &termination_status)) {
            printf("%d\n", termination_flag);
            // Sample sea water height
            printf("Looping...\n");
            sleep(TSUNAMETER_POLL);
            float new_val = generate_float_val(100000.0);
            append_moving_avg(avg, new_val);

            if (get_moving_avg(avg) > THRESHOLD) {
                // Sample neighbours for average

                // Send information to base station
            }

            time_t t = time(NULL);
            struct tm tm = *localtime(&t);
            printf("Time: %d-%02d-%02d %02d:%02d:%02d\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
            printf("Size: %d, Average: %lf, Head: %lf, Tail: %lf\n\n",avg -> size, avg -> avg, avg -> queue_head -> value, avg -> queue_tail -> value);
        }

        // Termination signal received
        // TODO: Free avg? idk
        printf("terminating...\n");
    }

    MPI_Finalize();
}


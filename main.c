#include <stdio.h>  // Printing
#include <mpi.h>    // OpenMPI calls
#include <stdlib.h> // Type conversion 
#include <errno.h>  // ERROR CHECKING
#include <limits.h> // INT LIMITS
#include "tsunameter.h" // Tsunameter functions

#define DIMENSIONS 2    // Tsunameter topology dimensions

int main(int argc, char** argv)
{
    int dims[DIMENSIONS];
    int comm_size;
    int comm_rank;
    MPI_Comm comm;
    int rc;

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
        printf("\nUsage: mpirun -np <PROCESSES> main m n\n");
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
            // Handle error
            printf("\nArgument error\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // Fill in dimensions
        dims[dim_index] = arg;
        check_sum *= dims[dim_index];
    }

    // Make sure allocated process matches topology
    if (check_sum != comm_size - 1) {
        printf("\nIncompatible topology and tsunameter count\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    printf("Dimensions: %dx%d, Size: %d\n", dims[0], dims[1], comm_size);

    // Virtual topology
    int reorder = 0;
    int period[DIMENSIONS] = {0};
    print("%d%d", period[0], period[1]);
    rc = MPI_Cart_create(MPI_COMM_WORLD, DIMENSIONS, dims, period, reorder, &comm);

    if (comm_rank == 0) {
        // TODO: Base Station logic
    } else {
        // TODO: Tsunameter logic
    }

    MPI_Finalize();
}


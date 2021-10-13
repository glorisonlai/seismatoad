#include <errno.h>  // ERROR CHECKING
#include <limits.h> // INT LIMITS
#include <math.h>   // abs
#include <mpi.h>    // OpenMPI calls
#include <stddef.h> // Offsetof - Useful for MPI_AInt
#include <stdio.h>  // Printing
#include <stdlib.h> // Type conversion
#include <time.h>   // Seed for random
#include <unistd.h> // Sleep

#include "tsunameter.h" // Tsunameter functions

#define DIMENSIONS 2        // Tsunameter topology dimensions
#define TERMINATION_TAG 0   // MPI Tag for termination
#define MAX_READING 10000   // Max sample value
#define TOLERANCE 100       // Tsunameter tolerance
#define TSUNAMTER_WINDOW 4 // Tsunameter average window
#define TSUNAMETER_POLL 2   // Tsunamter polling rate (s)

int main(int argc, char **argv) {
    int tsunameter_dims[DIMENSIONS];
    int tsunameter_threshold;
    int tsunameter_rank, num_tsunameters;
    int base_size;
    int base_rank;
    int rc;
    int root = 0;
    MPI_Comm tsunameter_comm;

    // Initialize MPI
    rc = MPI_Init(&argc, &argv);
    // Check MPI Init is successful
    if (rc != MPI_SUCCESS) {
        printf("\nCould not initialize MPI\n");
        fflush(stdout);
        return 1;
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &base_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &base_size);

    // Ensure tsunameter topology is given
    if (argc != DIMENSIONS + 3) {
        if (base_rank == root) {
        printf("\nUsage: mpirun -np <PROCESSES> main m n threshold sentinel\n");
        fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    // Create group of tsunameters
    MPI_Comm_split(MPI_COMM_WORLD, base_rank == root, 0, &tsunameter_comm);
    if (base_rank != root) {
        MPI_Comm_size(tsunameter_comm, &num_tsunameters);
        MPI_Comm_rank(tsunameter_comm, &tsunameter_rank);
    }

    // Initialize topology from CL args
    errno = 0;
    int check_sum = 1;
    for (int dim_index = 0; dim_index < DIMENSIONS + 1; dim_index++) {
        char *p;
        long arg = strtol(argv[dim_index + 1], &p, 10); // Get argument index from 1
        if (errno != 0 || *p != '\0' || arg < INT_MIN || arg > INT_MAX) {
        // Handle type conversion error
        if (base_rank == root) {
            printf("\nError with dimension %d argument\n", dim_index + 1);
            fflush(stdout);
            MPI_Comm_free(&tsunameter_comm);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        }

        // Fill in dimensions
        if (dim_index < DIMENSIONS) {
            tsunameter_dims[dim_index] = arg;
            check_sum *= tsunameter_dims[dim_index];
        } else {
            // Instantiate threshold
            if ((arg < 0 || arg > MAX_READING) && base_rank == root) {
                printf("\nError with threshold argument! (Must be betwee\n");
                fflush(stdout);
                MPI_Comm_free(&tsunameter_comm);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            tsunameter_threshold = arg;
        }
    }

    printf("Dimensions: %dx%d, Size: %d\n", tsunameter_dims[0],
            tsunameter_dims[1], base_size);
    if (base_rank != root) {
        // Check allocated process matches topology
        if (check_sum != num_tsunameters) {
            printf("\nIncompatible topology and tsunameter count\n");
            fflush(stdout);
            MPI_Comm_free(&tsunameter_comm);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // Create tsunameter cartesian mapping
        MPI_Dims_create(num_tsunameters, DIMENSIONS, tsunameter_dims);
        int wrapping[DIMENSIONS] = {}; // No wrapping
        int reorder = 0;    // Don't change original comm rank
        rc = MPI_Cart_create(tsunameter_comm, DIMENSIONS, tsunameter_dims, wrapping,
                            reorder, &tsunameter_comm);
        if (rc != 0) {
            printf("Error creating Cartesian mapping\n");
            fflush(stdout);
            MPI_Comm_free(&tsunameter_comm);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }


    // Create MPI datatype for struct tsunameter_reading
    // Source:
    // https://www.codingame.com/playgrounds/349/introduction-to-mpi/custom-types---exercise
    MPI_Datatype mp_tsunameter_reading;
    int block_length[] = {1, 1}; // {float, int}
    MPI_Aint block_displacement[] = {offsetof(tsunameter_reading, avg),
                                    offsetof(tsunameter_reading, time)};
    MPI_Datatype types[] = {MPI_FLOAT, MPI_INT};

    MPI_Type_create_struct(2, block_length, block_displacement, types,
                            &mp_tsunameter_reading);
    MPI_Type_commit(&mp_tsunameter_reading);

        // Main logic
    if (base_rank == root) {
        int blah = 0;
        MPI_Request send_req;
        sleep(TSUNAMETER_POLL * 10);
        printf("TERMINATING...\n");
        MPI_Ibcast(&blah, 1, MPI_INT, root, MPI_COMM_WORLD, &send_req);
        // TODO: Base Station logic
    } else {
        srand(time(NULL) + tsunameter_rank); // Seed random gen for tsunameters

        // Instantiate moving average
        printf("rank %d\n", tsunameter_rank);
        moving_avg *avg = init_moving_avg(TSUNAMTER_WINDOW);

        // // Get array of neighbours
        int num_neighbours;
        int *neighbours = get_neighbours(tsunameter_comm, DIMENSIONS, &num_neighbours);
        printf("rank %d: Neighbours: %d\n", tsunameter_rank, num_neighbours);
        for (int i = 0; i < num_neighbours; i++) {
            printf("  - %d\n", neighbours[i]);
        }
        printf("\n");

        // Send preliminary receive for termination signal
        int termination_buffer[1];
        MPI_Request termination_req;
        MPI_Ibcast(termination_buffer, 1, MPI_INT, root, MPI_COMM_WORLD,
                &termination_req);

        // Send num_neighbour preliminary receives for compare signals
        MPI_Request comparison_reqs[num_neighbours];
        int comparison_buffer[num_neighbours];
        for (int nbr_i = 0; nbr_i < num_neighbours; nbr_i++) {
            MPI_Irecv(&comparison_buffer[nbr_i], 1, MPI_INT, neighbours[nbr_i], 2*tsunameter_rank,
                        tsunameter_comm, &comparison_reqs[nbr_i]);
        }

        // Check for termination signal
        while (!test_mpi_req(&termination_req,0,tsunameter_rank,-1)) {
            // Sample sea water height with random value and append to moving average
            printf("%d Looping...\n", tsunameter_rank);
            float new_val = generate_float_val(10000.0);
            append_moving_avg(avg, new_val);
            printf("Rank: %d, new %f, Average: %f\n", tsunameter_rank, new_val,
                    get_moving_avg(avg));

            time_t curr_time = time(NULL);

            // Check if neighbours sent compare request signal
            printf("Rank: %d, Testing...\n", tsunameter_rank);
            for (int nbr_i = 0; nbr_i < num_neighbours; nbr_i++) {
                if (test_mpi_req(&comparison_reqs[nbr_i],0,tsunameter_rank,nbr_i)) {
                    struct tsunameter_reading *curr_reading =
                        instantiate_tsunameter_reading(get_moving_avg(avg), curr_time);

                    // Send information
                    printf("Rank %d Sending to %d with avg %f, time %d\n", tsunameter_rank,
                            neighbours[nbr_i], curr_reading->avg, curr_reading->time);
                    MPI_Send(curr_reading, 1, mp_tsunameter_reading, neighbours[nbr_i], 2*neighbours[nbr_i] + 1,
                            tsunameter_comm);
                    printf("Rank %d sent!\n", tsunameter_rank);
                    free(curr_reading);

                    // Reset comparison request
                    MPI_Irecv(&comparison_buffer[nbr_i], 1, MPI_INT, neighbours[nbr_i], 2*tsunameter_rank,
                                tsunameter_comm, &comparison_reqs[nbr_i]);
                }
            }

            // Potential tsunami alert
            if (get_moving_avg(avg) > tsunameter_threshold) {
                printf("Rank: %d, Requesting comps\n", tsunameter_rank);
                double starttime = MPI_Wtime();

                // Sample neighbours for average
                MPI_Request send_reqs[num_neighbours], recv_reqs[num_neighbours];
                tsunameter_reading recv_buff[num_neighbours];
                int send_buf;

                // Ping other tsunameters to send their averages, and post corresponding
                // receive
                for (int nbr_i = 0; nbr_i < num_neighbours; nbr_i++) {
                    MPI_Isend(&send_buf, 1, MPI_INT, neighbours[nbr_i], 2*neighbours[nbr_i],
                                tsunameter_comm, &send_reqs[nbr_i]);
                    printf("Rank %d sent to %d with tag %d\n", tsunameter_rank,
                            neighbours[nbr_i], 0);
                    MPI_Irecv(&recv_buff[nbr_i], 1, MPI_FLOAT, neighbours[nbr_i], 2*tsunameter_rank + 1,
                                tsunameter_comm, &recv_reqs[nbr_i]);
                }

                // Compare received readings with itself
                // Timeout after 10 s
                while (MPI_Wtime() - starttime < 10) {
                    // Repeatedly check for similar readings
                    // TODO: Add comparison between process times
                    int similar_count = 0, in_count = 0;
                    for (int nbr_i = 0; nbr_i < num_neighbours; nbr_i++) {
                        // Continue checking if neighbours sent compare request
                        printf("Rank %d testing comparison in loop\n", tsunameter_rank);
                        if (test_mpi_req(&comparison_reqs[nbr_i], 2, tsunameter_rank, nbr_i)) {
                            struct tsunameter_reading *curr_reading =
                                instantiate_tsunameter_reading(get_moving_avg(avg),
                                                                curr_time);

                            // Send information
                            printf("Rank %d Sending during wait loop to %d with avg %f, time %d\n", tsunameter_rank,
                                    neighbours[nbr_i], curr_reading->avg, curr_reading->time);
                            MPI_Send(curr_reading, 1, mp_tsunameter_reading, neighbours[nbr_i], 2*neighbours[nbr_i] + 1,
                                    tsunameter_comm);
                            printf("Rank %d sent!\n", tsunameter_rank);
                            free(curr_reading);

                            // Reset comparison request
                            MPI_Irecv(comparison_buffer, 1, MPI_INT, neighbours[nbr_i], 2*tsunameter_rank,
                                        tsunameter_comm, &comparison_reqs[nbr_i]);
                        }

                        printf("Rank %d testing receives in loop\n", tsunameter_rank);
                        if (test_mpi_req(&recv_reqs[nbr_i], 1, tsunameter_rank, nbr_i)) {
                            in_count += 1;
                            if (fabsf(recv_buff[nbr_i].avg) - get_moving_avg(avg) <
                                tsunameter_threshold) {
                                similar_count += 1;
                            }
                        }
                    }

                    printf("Rank: %d, Similar_count: %d\n", tsunameter_rank, similar_count);
                    // Send information to base station
                    if (similar_count >= 2) {
                        printf("Rank: %d sending to base station\n", tsunameter_rank);
                        tsunameter_reading *curr_reading =
                            instantiate_tsunameter_reading(get_moving_avg(avg), curr_time);
                        MPI_Send(curr_reading, 1, mp_tsunameter_reading, root, 0,
                                MPI_COMM_WORLD);
                        free(curr_reading);
                        break;
                    }

                    if (in_count == num_neighbours) {
                        break;
                    }
                    sleep(1);
                }

                for (int nbr_i = 0; nbr_i < num_neighbours; nbr_i++) {
                    MPI_Cancel(&recv_reqs[nbr_i]);
                }

                double endtime = MPI_Wtime();
                printf("%d finished checking,  < %d num_neighours, time = %f - %f = "
                    "%f < %d\n",
                    tsunameter_rank, num_neighbours, endtime, starttime,
                    endtime - starttime, endtime - starttime < 5);
            }
            
            sleep(TSUNAMETER_POLL);
        }

        // Termination signal received
        free(avg);
        free(neighbours);

        printf("%d terminating...\n", base_rank);
    }

    MPI_Comm_free(&tsunameter_comm);
    MPI_Finalize();
}

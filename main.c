#include <errno.h>   // ERROR CHECKING
#include <limits.h>  // INT LIMITS
#include <math.h>    // abs
#include <mpi.h>     // OpenMPI calls
#include <stddef.h>  // Offsetof - Useful for MPI_AInt
#include <stdio.h>   // Printing
#include <stdlib.h>  // Type conversion
#include <time.h>    // Seed for random
#include <unistd.h>  // Sleep

#include "tsunameter.h"  // Tsunameter functions

#define DIMENSIONS 2         // Tsunameter topology dimensions
#define TERMINATION_TAG 0    // MPI Tag for termination
#define THRESHOLD 6000       // Tsunameter threshold
#define TOLERANCE 100        // Tsunameter tolerance
#define TSUNAMTER_WINDOW 10  // Tsunameter average window
#define TSUNAMETER_POLL 2    // Tsunamter polling rate (s)

int main(int argc, char **argv) {
    int tsunameter_dims[DIMENSIONS];
    int comm_size;
    int comm_rank;
    int rc;
    int root = 0;

    // Initialize MPI
    rc = MPI_Init(&argc, &argv);
    // Check MPI Init is successful
    if (rc != MPI_SUCCESS) {
        printf("\nCould not initialize MPI\n");
        fflush(stdout);
        return 1;
    }
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

    // Ensure tsunameter topology is given
    if (argc != DIMENSIONS + 1) {
        printf("\nUsage: mpirun -np <PROCESSES> main m n\n");
        fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Initialize topology
    errno = 0;
    int dim_index;
    int check_sum = 1;
    for (dim_index = 0; dim_index < DIMENSIONS; dim_index++) {
        char *p;
        long arg =
            strtol(argv[dim_index + 1], &p, 10);  // Get argument index from 1
        if (errno != 0 || *p != '\0' || arg < INT_MIN || arg > INT_MAX) {
            // Handle type conversion error
            printf("\nArgument error\n");
            fflush(stdout);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // Fill in dimensions
        tsunameter_dims[dim_index] = arg;
        check_sum *= tsunameter_dims[dim_index];
    }

    // Check allocated process matches topology
    if (check_sum != comm_size - 1) {
        printf("\nIncompatible topology and tsunameter count\n");
        fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    printf("Dimensions: %dx%d, Size: %d\n", tsunameter_dims[0],
           tsunameter_dims[1], comm_size);
    // Create MPI datatype for struct tsunameter_reading
    // Source:
    // https://www.codingame.com/playgrounds/349/introduction-to-mpi/custom-types---exercise
    MPI_Datatype mp_tsunameter_reading;
    int block_length[] = {1, 1};  // {float, int}
    MPI_Aint block_displacement[] = {offsetof(tsunameter_reading, avg),
                                      offsetof(tsunameter_reading, time)};
    MPI_Datatype types[] = {MPI_FLOAT, MPI_INT};

    MPI_Type_create_struct(2, block_length, block_displacement, types,
                           &mp_tsunameter_reading);
    MPI_Type_commit(&mp_tsunameter_reading);

    // Main logic
    if (comm_rank == root) {
        int blah = 0;
        MPI_Request send_req;
        sleep(TSUNAMETER_POLL * 10);
        printf("TERMINATING...\n");
        MPI_Ibcast(&blah, 1, MPI_INT, root, MPI_COMM_WORLD, &send_req);
        // TODO: Base Station logic
    } else {
        srand(time(NULL) + comm_rank);  // Seed random gen for tsunameters

        // Instantiate moving average
        printf("rank %d\n", comm_rank);
        moving_avg *avg = init_moving_avg(TSUNAMTER_WINDOW);

        // Instantiate topology
        int coord[DIMENSIONS];
        get_coord_at_rank(tsunameter_dims, DIMENSIONS, comm_rank, coord);

        // Get array of neighbours
        int *neighbours = malloc(sizeof(*neighbours));
        int num_neighbours;
        get_neighbours(comm_rank - 1, tsunameter_dims[0], tsunameter_dims[1],
                       neighbours, &num_neighbours);
        printf("rank %d: neighbour: %d\n", comm_rank, neighbours[0]);

        // Send one receive for termination signal
        int termination_flag = 0;
        int termination_buffer[1];
        MPI_Request termination_req;
        MPI_Status termination_status;
        MPI_Ibcast(termination_buffer, 1, MPI_INT, root, MPI_COMM_WORLD,
                   &termination_req);

        // Send num_neighbour receives for compare signals
        int comparison_flag[num_neighbours];
        MPI_Request comparison_reqs[num_neighbours];
        MPI_Status comparison_status[num_neighbours];
        int comparison_buffer[num_neighbours];
        int i;
        for (i = 0; i < num_neighbours; i++) {
            MPI_Irecv(comparison_buffer, 1, MPI_INT, neighbours[i], 0,
                      MPI_COMM_WORLD, &comparison_reqs[i]);
        }

        // Check for termination signal
        while (!test_mpi_req(&termination_req, &termination_flag,
                             &termination_status)) {
            // Sample sea water height
            printf("%d Looping...\n", comm_rank);
            sleep(TSUNAMETER_POLL);
            float new_val = generate_float_val(10000.0);
            append_moving_avg(avg, new_val);
            printf("Rank: %d, new %f, Average: %f\n", comm_rank, new_val,
                   get_moving_avg(avg));

            time_t curr_time = time(NULL);
            float f_time = (float) curr_time;

            // Check send request signal
            int i;
            printf("Rank: %d, Testing...\n", comm_rank);
            for (i = 0; i < num_neighbours; i++) {
                if (test_mpi_req(&comparison_reqs[i], &comparison_flag[i],
                                 &comparison_status[i])) {
                    struct tsunameter_reading *curr_reading =
                        instantiate_tsunameter_reading(get_moving_avg(avg),
                                                       curr_time);
                    float test_reading[] = {get_moving_avg(avg)} ;

                    // Send information
                    printf("Rank %d Sending to %d with avg %f\n", comm_rank, neighbours[i], test_reading[0]);
                    MPI_Send(test_reading, 1, MPI_FLOAT,
                             neighbours[i], 0, MPI_COMM_WORLD);
                    printf("Rank %d sent!\n", comm_rank);
                    free(curr_reading);

                    // Reset comparison request
                    comparison_flag[i] = 0;
                    MPI_Irecv(comparison_buffer, 1, MPI_INT, neighbours[i], 0,
                              MPI_COMM_WORLD, &comparison_reqs[i]);
                }
            }

            // Potential tsunami alert
            if (get_moving_avg(avg) > THRESHOLD) {
                printf("Rank: %d, Requesting comps\n", comm_rank);
                // Sample neighbours for average
                MPI_Request send_reqs[num_neighbours];
                tsunameter_reading recv_buff[num_neighbours];
                MPI_Request recv_reqs[num_neighbours];
                MPI_Status recv_status[num_neighbours];
                int indices[num_neighbours];
                int i;
                int send_buf[1];

                float test_recv_buff[num_neighbours];
                MPI_Request test_reqs[num_neighbours];

                for (i = 0; i < num_neighbours; i++) {
                    MPI_Send(send_buf, 1, MPI_INT, neighbours[i], 0,
                             MPI_COMM_WORLD);
                    printf("Rank %d sent to %d with tag %d\n", comm_rank,
                           neighbours[i], 0);
                    MPI_Recv(test_recv_buff, 1, MPI_FLOAT,
                              neighbours[i], 0, MPI_COMM_WORLD, &recv_status[i]);
                }

                // Repeatedly check for similar readings
                // TODO: Add comparison between process times
                int similar_count;
                int outcount;
                double starttime, endtime;
                starttime = MPI_Wtime();
                printf("Rank %d Checking...\n", comm_rank);
                do {
                    similar_count = 0;
                    // FIXME: Will hang if noone sends (maybe after termination)
                    // MPI_Waitall(num_neighbours, recv_reqs, recv_status);
                    // printf("Rank %d Outcount: %d\n", comm_rank, outcount);
                    int j;
                    for (j=0; j<num_neighbours; j++) {
                        // printf("indices[%d]: %d\n", j, indices[j]);
                        printf("received avg = %f\n", test_recv_buff[0]);
                        // printf("received time = %f\n", test_recv_buff[1]);
                    }
                    for (i = 0; i < num_neighbours; i++) {
                        printf("%d Received %f\n", comm_rank, test_recv_buff[0]);
                        if (fabsf(test_recv_buff[0] -
                                  get_moving_avg(avg)) < THRESHOLD) {
                            similar_count += 1;
                        }
                    }

                    // Send information to base station
                    printf("Rank: %d, Similar_count: %d\n", comm_rank,
                           similar_count);
                    if (similar_count >= 2) {
                        printf("Rank: %d sending to base station\n", comm_rank);
                        tsunameter_reading *curr_reading =
                            instantiate_tsunameter_reading(get_moving_avg(avg),
                                                           curr_time);
                        float test_reading[] = {get_moving_avg(avg)};
                        MPI_Send(curr_reading, 1, mp_tsunameter_reading, root,
                                 0, MPI_COMM_WORLD);
                        free(curr_reading);
                        break;
                    }
                    sleep(1);
                    endtime = MPI_Wtime();
                    printf("%d finished checking,  < %d num_neighours, time = %f - %f = %f < %d\n", comm_rank, num_neighbours, endtime, starttime, endtime-starttime, endtime-starttime < 5);
                } while (
                         endtime - starttime < 5);  // Timeout after 10 seconds
                printf("%d Finished loop\n", comm_rank);
            }
        }

        // Termination signal received
        free(avg);
        free(neighbours);

        printf("%d terminating...\n", comm_rank);
    }

    MPI_Finalize();
}

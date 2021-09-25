#include <stdio.h>      // Printing
#include <mpi.h>        // OpenMPI calls
#include <stdlib.h>     // Type conversion
#include <errno.h>      // ERROR CHECKING
#include <limits.h>     // INT LIMITS
#include "tsunameter.h" // Tsunameter functions
#include <unistd.h>     // Sleep
#include <time.h>       // Seed for random
#include <stddef.h>     // Offsetof - Useful for MPI_AInt
#include <math.h>       // abs

#define DIMENSIONS 2        // Tsunameter topology dimensions
#define TERMINATION_TAG 0   // MPI Tag for termination
#define THRESHOLD 6000      // Tsunameter threshold
#define TOLERANCE 100       // Tsunameter tolerance
#define TSUNAMTER_WINDOW 10 // Tsunameter average window
#define TSUNAMETER_POLL 2   // Tsunamter polling rate (s)

int main(int argc, char **argv)
{
    int tsunameter_dims[DIMENSIONS];
    int comm_size;
    int comm_rank;
    int rc;
    int root = 0;

    // Testing cart create
    MPI_Comm new_comm;

    srand(time(NULL)); // Seed random gen for tsunameters

    // Initialize MPI
    rc = MPI_Init(&argc, &argv);
    // Check MPI Init is successful
    if (rc != MPI_SUCCESS)
    {
        printf("\nCould not initialize MPI\n");
        fflush(stdout);
        return 1;
    }
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

    // Ensure tsunameter topology is given
    if (argc != DIMENSIONS + 1)
    {
        printf("\nUsage: mpirun -np <PROCESSES> main m n\n");
        fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Initialize topology
    errno = 0;
    int dim_index;
    int check_sum = 1;
    for (dim_index = 0; dim_index < DIMENSIONS; dim_index++)
    {
        char *p;
        long arg = strtol(argv[dim_index + 1], &p, 10); // Get argument index from 1
        if (errno != 0 || *p != '\0' || arg < INT_MIN || arg > INT_MAX)
        {
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
    if (check_sum != comm_size - 1)
    {
        printf("\nIncompatible topology and tsunameter count\n");
        fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    printf("Dimensions: %dx%d, Size: %d\n", tsunameter_dims[0], tsunameter_dims[1], comm_size);

    // Create MPI datatype for struct tsunameter_reading
    // Source: https://www.codingame.com/playgrounds/349/introduction-to-mpi/custom-types---exercise
    MPI_Datatype mp_tsunameter_reading;
    int block_length[2] = {1,1};    // {float, int}
    MPI_Aint block_displacement[2] = {offsetof(tsunameter_reading, avg), offsetof(tsunameter_reading, time)};
    MPI_Datatype types[2] = {MPI_FLOAT, MPI_INT};

    MPI_Type_create_struct(2, block_length, block_displacement, types, &mp_tsunameter_reading);
    MPI_Type_commit(&mp_tsunameter_reading);

    // Main logic
    if (comm_rank == root)
    {
        // int blah = 0;
        // MPI_Request send_req;
        // sleep(10);
        // MPI_Ibcast(&blah, 1, MPI_INT, root, MPI_COMM_WORLD, &send_req);
        // TODO: Base Station logic
    }
    else
    {
        // Instantiate moving average
        printf("rank %d\n", comm_rank);
        moving_avg *avg = init_moving_avg(TSUNAMTER_WINDOW);

        // Instantiate topology
        int coord[DIMENSIONS];
        get_coord_at_rank(tsunameter_dims, DIMENSIONS, comm_rank, coord);

        // Get array of neighbours
        int *neighbours;
        int *num_neighbours;
        get_neighbours(comm_rank - 1, tsunameter_dims[0], tsunameter_dims[1], neighbours, num_neighbours);
        MPI_Request tsunameter_req_arr[*num_neighbours];

        // Send one receive for termination signal
        int termination_flag = 0;
        int termination_buffer[1];
        MPI_Request termination_req;
        MPI_Status termination_status;
        MPI_Ibcast(termination_buffer, 1, MPI_INT, root, MPI_COMM_WORLD, &termination_req);

        // Send num_neighbour receives for compare signals
        int comparison_flag[*num_neighbours];
        MPI_Request comparison_reqs[*num_neighbours];
        MPI_Status comparison_status[*num_neighbours];
        int comparison_buffer[*num_neighbours];
        int i;
        for (i = 0; i < *num_neighbours; i++)
        {
            MPI_Irecv(comparison_buffer, 1, MPI_INT, neighbours[i], 0, MPI_COMM_WORLD, &comparison_reqs[i]);
        }

        // Check for termination signal
        while (!test_mpi_req(&termination_req, &termination_flag, &termination_status))
        {
            // Sample sea water height
            printf("Looping...\n");
            sleep(TSUNAMETER_POLL);
            float new_val = generate_float_val(100000.0);
            append_moving_avg(avg, new_val);

            time_t curr_time = time(NULL);

            // Check send request signal
            int i;
            for (i = 0; i < *num_neighbours; i++)
            {
                if (test_mpi_req(&comparison_reqs[i], &comparison_flag[i], &comparison_status[i]))
                {
                    // Send information
                    tsunameter_reading *curr_reading = instantiate_tsunameter_reading(get_moving_avg(avg), curr_time);
                    MPI_Send(curr_reading, 1, mp_tsunameter_reading, neighbours[i], 0, MPI_COMM_WORLD);

                    // Reset comparison request
                    comparison_flag[i] = 0;
                    MPI_Irecv(comparison_buffer, 1, MPI_INT, neighbours[i], 0, MPI_COMM_WORLD, &comparison_reqs[i]);
                }
            }

            // Potential tsunami alert
            if (get_moving_avg(avg) > THRESHOLD)
            {
                // Sample neighbours for average
                MPI_Request send_reqs[*num_neighbours];
                tsunameter_reading recv_buff[*num_neighbours];
                MPI_Request recv_reqs[*num_neighbours];
                MPI_Status recv_status[*num_neighbours];
                int indices[*num_neighbours];
                int i;
                int send_buf[1];
                for (i = 0; i < *num_neighbours; i++)
                {
                    MPI_Isend(send_buf, 1, MPI_INT, neighbours[i], 0, MPI_COMM_WORLD, &send_reqs[i]);
                    MPI_Irecv(&recv_buff[i], 1, mp_tsunameter_reading, neighbours[i], 0, MPI_COMM_WORLD, &recv_reqs[i]);
                }

                // Repeatedly check for similar readings
                // TODO: Add comparison between process times
                int similar_count;
                int outcount = 0;
                double starttime, endtime;
                starttime = MPI_Wtime();
                do
                {
                    similar_count = 0;
                    MPI_Waitsome(*num_neighbours, recv_reqs, &outcount, indices, recv_status);
                    for (i = 0; i < outcount; i++)
                    {
                        if (fabsf(recv_buff[indices[i]].avg - get_moving_avg(avg)) < THRESHOLD)
                        {
                            similar_count += 1;
                        }
                    }

                    // Send information to base station
                    if (similar_count >= 2)
                    {
                        tsunameter_reading *curr_reading = instantiate_tsunameter_reading(get_moving_avg(avg), curr_time);
                        MPI_Send(curr_reading, 1, mp_tsunameter_reading, root, 0, MPI_COMM_WORLD);
                        break;
                    }
                    endtime = MPI_Wtime();
                } while (outcount < *num_neighbours || endtime - starttime < 10);   // Timeout after 10 seconds
            }
        }

        // Termination signal received
        free(avg);
        free(neighbours);

        printf("terminating...\n");
    }

    MPI_Finalize();
}

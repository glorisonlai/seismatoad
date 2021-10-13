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

int main(int argc, char **argv) {
    int base_size;
    int base_rank;

    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &base_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &base_size);

    MPI_Datatype mp_tsunameter_reading;
    int block_length[] = {1, 1}; // {float, int}
    MPI_Aint block_displacement[] = {offsetof(tsunameter_reading, avg),
                                    offsetof(tsunameter_reading, time)};
    MPI_Datatype types[] = {MPI_FLOAT, MPI_INT};

    MPI_Type_create_struct(2, block_length, block_displacement, types,
                            &mp_tsunameter_reading);
    MPI_Type_commit(&mp_tsunameter_reading);

    switch (base_rank) {
       case (0): {
           sleep(5);
           time_t curr_time = time(NULL);
           tsunameter_reading *send = instantiate_tsunameter_reading(1.2345, curr_time);
           tsunameter_reading buf[1];
           MPI_Request request;
           printf("Sending %f, %d\n", send->avg, send->time);
           MPI_Send(send, 1, mp_tsunameter_reading, 1, 0, MPI_COMM_WORLD); 
           break;
       } 
       case (1): {
            MPI_Request request[2];
            tsunameter_reading buf[2];
            MPI_Irecv(&buf[0], 1, mp_tsunameter_reading, 0, 0, MPI_COMM_WORLD, &request[0]);
            MPI_Irecv(&buf[1], 1, mp_tsunameter_reading, 2, 0, MPI_COMM_WORLD, &request[1]);
            int neighbours = 2;
            int in_count = 0;
            while (in_count < neighbours) {
                in_count = 0;
                for (int i = 0; i<neighbours; i++) {
                    if (test_mpi_req(&request[i])){
                        in_count += 1;
                    }
                    sleep(1);
                }
            }
            // test_mpi_req(&request[0]);
            // test_mpi_req(&request[0]);

            test_mpi_req(&request[0]);
            printf("Received %f, %d, %f, %d\n", buf[0].avg, buf[0].time, buf[1].avg, buf[1].time);
            break;
       }
       case (2): {
            sleep(3);
            time_t curr_time = time(NULL);
            tsunameter_reading *send = instantiate_tsunameter_reading(2.2345, curr_time);
            tsunameter_reading buf[1];
            MPI_Request request;
            printf("Sending %f, %d\n", send->avg, send->time);
            MPI_Send(send, 1, mp_tsunameter_reading, 1, 0, MPI_COMM_WORLD); 
            break;
       }
       default: {break;}
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}

#include <errno.h>  // ERROR CHECKING
#include <limits.h> // INT LIMITS
#include <math.h>   // abs
#include <mpi.h>    // OpenMPI calls
#include <stddef.h> // Offsetof - Useful for MPI_AInt
#include <stdio.h>  // Printing
#include <stdlib.h> // Type conversion
#include <time.h>   // Seed for random
#include <unistd.h> // Sleep
#include <fcntl.h>  // Used for reading input safely
#include <string.h> // useful for string comparison
#include <pthread.h> // posix

#include "tsunameter.h" // Tsunameter functions
#include "satellite.h" // Satellite Readings struct

#define DIMENSIONS 2        // Tsunameter topology dimensions
#define TERMINATION_TAG 0   // MPI Tag for termination
#define MAX_READING 10000   // Max sample value
#define TOLERANCE 3000       // Tsunameter tolerance
#define TSUNAMTER_WINDOW 20 // Tsunameter average window
#define TSUNAMETER_POLL 2   // Tsunamter polling rate (s)
#define STORED_READINGS 20 // How many satellite readings to store

/* Global Variables for POSIX Use */
satellite_reading *satellite_readings;
int satellite_terminate;
pthread_mutex_t satellite_mutex;

int sentinel_terminate;
pthread_mutex_t sentinel_mutex;

int comms_terminate;
pthread_mutex_t comms_mutex;

/* Function declarations */
void* run_satellite(void* args);
void* check_sentinel(void* arg);
void* tsunameter_communication(void* arg);
void* run_comms(void* args);


int main(int argc, char **argv) {
    int tsunameter_dims[DIMENSIONS];
    int tsunameter_threshold;
    int tsunameter_rank, num_tsunameters;
    int base_size;
    int base_rank;
    int rc;
    int root = 0;
    int iterations;
    
    MPI_Comm tsunameter_comm;
/*+++++++++++++++++++++++++++++++ SETUP CODE +++++++++++++++++++++++++++++++*/
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
    if (argc != DIMENSIONS + 4) {
        if (base_rank == root) {
        printf("\nUsage: mpirun -np <PROCESSES> main m n threshold sentinel\n");
        fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    
    // read in the number of iterations, ensuring it is at least 1
    iterations = atoi(argv[4]);
    if (iterations <= 0){
        if (base_rank == root) {
            printf("\nMust have at least 1 iteration\n");
            fflush(stdout);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        return 1;
    }

/*+++++++++++++++++++++++++ TSUNAMETER CONSTRUCTION +++++++++++++++++++++++++*/

    // Create group of tsunameters
    MPI_Comm_split(MPI_COMM_WORLD, base_rank == root, 0, &tsunameter_comm);
    if (base_rank != root) {
        MPI_Comm_size(tsunameter_comm, &num_tsunameters);
        MPI_Comm_rank(tsunameter_comm, &tsunameter_rank);
    }

    // Initialize topology from CL args
    errno = 0;
    int check_sum = 1;
    int dim_index;
    for (dim_index = 0; dim_index < DIMENSIONS + 1; dim_index++) {
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

    printf("Dimensions: %dx%d, Size: %d, Threshold: %d\n", tsunameter_dims[0],
            tsunameter_dims[1], base_size, tsunameter_threshold);
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

    MPI_Datatype mp_base_station_info;
    int iblock_length[] = {1, 1, 4}; // {float, int, int * neighbours}
    MPI_Aint iblock_displacement[] = {offsetof(base_station_info, avg),
                                    offsetof(base_station_info, time),
                                    offsetof(base_station_info, neighbours)};
    MPI_Datatype itypes[] = {MPI_FLOAT, MPI_INT, MPI_INT};

    MPI_Type_create_struct(3, iblock_length, iblock_displacement, itypes,
                            &mp_base_station_info);
    MPI_Type_commit(&mp_base_station_info);

    
        // Main logic
/*+++++++++++++++++++++++++++++++ BASE STATION +++++++++++++++++++++++++++++++*/


    if (base_rank == root) {
        //******************** Satellite Altimeter Code ********************
        // First initialise a thread to simulate the altimeter
            //setup the satellite readings array first
        
        printf("Creating Satellite Thread\n");
        int sat_records = STORED_READINGS;
        int sat_args[3] = {sat_records, tsunameter_dims[0], tsunameter_dims[1]}; 
        satellite_terminate = 0;
        satellite_readings = (satellite_reading*)malloc((sat_records) * sizeof(satellite_reading));
        pthread_t satellite_tid;
        pthread_mutex_init(&satellite_mutex, NULL);
        pthread_create(&satellite_tid, NULL, run_satellite, &sat_args);
        
        
        //******************** Sentinel Value Checking Code ********************
        // Construct an additional thread to handle sentinel exit values
        // POSIX HERE
        
        printf("Creating Sentinel Thread\n");
        pthread_t sentinel_tid;
        pthread_mutex_init(&sentinel_mutex, NULL);
        sentinel_terminate = 0;
        printf("%s\n", argv[5]);
        pthread_create(&sentinel_tid, NULL, check_sentinel, argv[5]);
        
        
        //********************** Comms code ********************
        // For HD+ this section needs to be a POSIX thread to protect against timeouts.
        // Then wait for an alert from one of the nodes
        printf("Creating Comms Thread\n");
        pthread_t comms_tid;
        pthread_mutex_init(&comms_mutex, NULL);
        int comms_args[3] = {iterations, tsunameter_dims[0], tsunameter_dims[1]};
        pthread_create(&comms_tid, NULL, run_comms, &comms_args);
        
        
        /* Glorison Base temp code -> Moved to run_comms
        int blah = 0;
        MPI_Request send_req;
        sleep(TSUNAMETER_POLL * 10);
        printf("TERMINATING...\n");
        MPI_Ibcast(&blah, 1, MPI_INT, root, MPI_COMM_WORLD, &send_req);
        */
        
        
        //********************** Base termination checking ********************
        
        int terminating = 0;
        do{
            sleep(1);
            if (sentinel_terminate == 1){
                comms_terminate = 1;
                satellite_terminate = 1;
                terminating = 1;
            } else if (comms_terminate == 1){
                sentinel_terminate = 1;
                satellite_terminate = 1;
                terminating = 1;
            }
        } while (terminating == 0);
        printf("Terminating all\n");
        

        //********************** Base Station Termination ********************
        // Output the state of the satellite readings at termination
        printf("Current State of satellite readings:\n");
        int i;
        for (i=0; i<sat_records; i++){
            printf("Record %d: Time is %d, Xpos is %d, Ypos is %d, elevation is %f\n",
                i, satellite_readings[i].time, satellite_readings[i].xpos, satellite_readings[i].ypos, satellite_readings[i].elevation);
                
        }
        
        // terminate and clean up the altimeter
        pthread_join(satellite_tid, NULL);
        free(satellite_readings);

        // Do a final wait all to make sure all finalised
        // MPI_Wait_All();
        // Also terminate and clean up the other two posix threads. Both the MPI Send Recv and the Input scanning
        pthread_join(sentinel_tid, NULL);
        pthread_join(comms_tid, NULL);
        
    } else {
    
    
/*+++++++++++++++++++++++++++++++ TSUNAMETERS +++++++++++++++++++++++++++++++*/


        srand(time(NULL) + tsunameter_rank); // Seed random gen for tsunameters

        // Instantiate moving average
        printf("rank %d\n", tsunameter_rank);
        moving_avg *avg = init_moving_avg(TSUNAMTER_WINDOW);

        // // Get array of neighbours
        int num_neighbours;
        int *neighbours = get_neighbours(tsunameter_comm, DIMENSIONS, &num_neighbours);
        printf("rank %d: Neighbours: %d\n", tsunameter_rank, num_neighbours);
        int i;
        for (i = 0; i < num_neighbours; i++) {
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
        int nbr_i;
        for (nbr_i = 0; nbr_i < num_neighbours; nbr_i++) {
            MPI_Irecv(&comparison_buffer[nbr_i], 1, MPI_INT, neighbours[nbr_i], 2*tsunameter_rank,
                        tsunameter_comm, &comparison_reqs[nbr_i]);
        }

        // Check for termination signal
        while (!test_mpi_req(&termination_req)) {
            // Sample sea water height with random value and append to moving average
            //printf("%d Looping...\n", tsunameter_rank);
            float new_val = generate_float_val(10000.0);
            append_moving_avg(avg, new_val);
            //printf("Rank: %d, new %f, Average: %f\n", tsunameter_rank, new_val,
            //        get_moving_avg(avg));

            time_t curr_time = time(NULL);

            // Check if neighbours sent compare request signal
            printf("Rank: %d, Testing...\n", tsunameter_rank);
            for (nbr_i = 0; nbr_i < num_neighbours; nbr_i++) {
                if (test_mpi_req(&comparison_reqs[nbr_i])) {
                    struct tsunameter_reading *curr_reading =
                        instantiate_tsunameter_reading(get_moving_avg(avg), curr_time);

                    // Send information
                    //printf("Rank %d Sending to %d with avg %f, time %d\n", tsunameter_rank,
                    //       neighbours[nbr_i], curr_reading->avg, curr_reading->time);
                    MPI_Send(curr_reading, 1, mp_tsunameter_reading, neighbours[nbr_i], 2*neighbours[nbr_i] + 1,
                            tsunameter_comm);
                    //printf("Rank %d sent!\n", tsunameter_rank);
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
                for (nbr_i = 0; nbr_i < num_neighbours; nbr_i++) {
                    MPI_Isend(&send_buf, 1, MPI_INT, neighbours[nbr_i], 2*neighbours[nbr_i],
                                tsunameter_comm, &send_reqs[nbr_i]);
                    //printf("Rank %d sent to %d with tag %d\n", tsunameter_rank,
                    //        neighbours[nbr_i], 0);
                    MPI_Irecv(&recv_buff[nbr_i], 1, mp_tsunameter_reading, neighbours[nbr_i], 2*tsunameter_rank + 1,
                                tsunameter_comm, &recv_reqs[nbr_i]);

                }

                // Compare received readings with itself
                // Timeout after 10 s
                int similar_count, in_count;
                int similar_neighbours[4] = {-1,-1,-1,-1};
                while (MPI_Wtime() - starttime < 10) {
                    // Repeatedly check for similar readings
                    // TODO: Add comparison between process times
                    
                    // Check if we're terminating to avoid waiting
                    if(test_mpi_req(&termination_req)){
                        break;
                    }
                    
                    
                    similar_count = 0;
                    in_count = 0;
                    similar_neighbours[0] = -1;
                    similar_neighbours[1] = -1;
                    similar_neighbours[2] = -1;
                    similar_neighbours[3] = -1;
                    for (nbr_i = 0; nbr_i < num_neighbours; nbr_i++) {
                        // Continue checking if neighbours sent compare request
                        //printf("Rank %d testing comparison in loop\n", tsunameter_rank);
                        
                    
                        if (test_mpi_req(&comparison_reqs[nbr_i])) {
                            struct tsunameter_reading *curr_reading =
                                instantiate_tsunameter_reading(get_moving_avg(avg),
                                                                curr_time);

                            // Send information
                            //printf("Rank %d Sending during wait loop to %d with avg %f, time %d\n", tsunameter_rank,
                            //        neighbours[nbr_i], curr_reading->avg, curr_reading->time);
                            MPI_Send(curr_reading, 1, mp_tsunameter_reading, neighbours[nbr_i], 2*neighbours[nbr_i] + 1,
                                    tsunameter_comm);
                            printf("Rank %d sent!\n", tsunameter_rank);
                            free(curr_reading);

                            // Reset comparison request
                            MPI_Irecv(&comparison_buffer[nbr_i], 1, MPI_INT, neighbours[nbr_i], 2*tsunameter_rank,
                                        tsunameter_comm, &comparison_reqs[nbr_i]);
                        }
                        printf("Rank %d testing receives in loop\n", tsunameter_rank);


                        if (test_mpi_req(&recv_reqs[nbr_i])) {
                            in_count += 1;
                            //printf("Rank %d Received avg: %f, time: %d from %d\n", tsunameter_rank, recv_buff[nbr_i].avg, recv_buff[nbr_i].time, neighbours[nbr_i]);
                            if (fabsf(recv_buff[nbr_i].avg - get_moving_avg(avg)) <
                                TOLERANCE) {
                                similar_neighbours[nbr_i] = neighbours[nbr_i];
                                similar_count += 1;
                            }
                        }
                    }

                    //printf("Rank: %d, Similar_count: %d\n", tsunameter_rank, similar_count);
                    // Send information to base station
                    if (similar_count >= 2) {
                        printf("Rank: %d sending to base station\n", tsunameter_rank);
                        struct base_station_info blah;
                        blah.avg = get_moving_avg(avg);
                        printf("Avg is %f\n", blah.avg);
                        blah.time = time(NULL);
                        printf("Time: %d\n", blah.time);
                        blah.neighbours[0] = similar_neighbours[0];
                        blah.neighbours[1] = similar_neighbours[1];
                        blah.neighbours[2] = similar_neighbours[2];
                        blah.neighbours[3] = similar_neighbours[3];
                        printf("Neighbours: %d %d %d %d\n", blah.neighbours[0], blah.neighbours[1],
                    blah.neighbours[2], blah.neighbours[3]);
                        MPI_Request tempstat;
                        MPI_Isend(&blah, 1, mp_base_station_info, root, 0,
                                MPI_COMM_WORLD, &tempstat);
                        printf("Sending completed by %d\n", tsunameter_rank);
                        
                        break;
                    }

                    if (in_count == num_neighbours) {
                        break;
                    }
                    sleep(1);
                }
            }
            
            sleep(TSUNAMETER_POLL);
        }

        // Termination signal received
        free(avg);
        free(neighbours);
        printf("%d terminating...\n", base_rank);
    } // End of Tsunameter (else) code
/*+++++++++++++++++++++++++++++++ CLEAN UP +++++++++++++++++++++++++++++++*/
    MPI_Comm_free(&tsunameter_comm);
    printf("Freed by %d\n", base_rank);
    return MPI_Finalize();
}



/*+++++++++++++++++++++++++++++ SATELLITE THREAD +++++++++++++++++++++++++++++*/

void* run_satellite(void* args){
    // arguments should be array size, then width and height of grid
    printf("Satellite Thread starting\n");

    int last_access = 0;
    int xpos, ypos, elevation;
    int* arguments = (int*)args;
    int size = arguments[0], width = arguments[1], height = arguments[2];
    do{
        // Make this properly FIFO, idk how tbh I think this is probably good enough.
        sleep(1);
        // generate a random co-ordinate
        srand(time(NULL));
        xpos = rand() % width;
        ypos = rand() % height;
        elevation = (rand() % 1000) + 6000;
        
        // Implement some quick Mutex here to make sure we avoid race conditions
        pthread_mutex_lock(&satellite_mutex);
        satellite_readings[last_access].time = time(NULL);
        satellite_readings[last_access].xpos = xpos;
        satellite_readings[last_access].ypos = ypos;
        satellite_readings[last_access].elevation = elevation;
        // Release mutex lock.
        pthread_mutex_unlock(&satellite_mutex);
        last_access = (last_access + 1) % size;
      
        
    } while (satellite_terminate == 0);
    printf("Satellite Thread terminating\n");
}

/*+++++++++++++++++++++++++++++ SENTINEL THREAD +++++++++++++++++++++++++++++*/

void* check_sentinel(void* arg){
    // while not terminating:
        // scanf() for input
        // check against sentinel
        // if sentinel return to main program after setting terminating to 1
        // Only this thread edits terminating
    printf("Sentinel Thread starting\n");
    char buf[strlen(arg)];
    char* sentinel = arg;
    do{
        fcntl(STDIN_FILENO, F_SETFL, fcntl(STDIN_FILENO, F_GETFL) | O_NONBLOCK);
        sleep(1);
        int readChars = read(STDIN_FILENO, buf, 4);
        if (readChars > 0) {
            // check whether buf now matches the sentinel value
            buf[strcspn(buf, "\n")] = 0;
            printf("Read in: %s\n", buf);
            
            if(strcmp(buf, sentinel) == 0){
                // Mutex may be unnecessary here as both the main and the sentinel threads will only ever set it to 1
                // However it is better to be safe to avoid a failing terminate.
                pthread_mutex_lock(&sentinel_mutex);
                sentinel_terminate = 1;
                pthread_mutex_unlock(&sentinel_mutex);
            } else{
                printf("Enter %s to exit\n", sentinel);
            }
        }
    } while(sentinel_terminate == 0);
    printf("Sentinel thread terminating\n");
}

/*+++++++++++++++++++++++++++++ COMMS THREAD +++++++++++++++++++++++++++++*/

void* run_comms(void* args){
    printf("Comms Thread Starting\n");
    
    // Setup type to receive

    MPI_Datatype mp_base_station_info;
    int iblock_length[] = {1, 1, 4}; // {float, int, int * neighbours}
    MPI_Aint iblock_displacement[] = {offsetof(base_station_info, avg),
                                    offsetof(base_station_info, time),
                                    offsetof(base_station_info, neighbours)};
    MPI_Datatype itypes[] = {MPI_FLOAT, MPI_INT, MPI_INT};

    MPI_Type_create_struct(3, iblock_length, iblock_displacement, itypes,
                            &mp_base_station_info);
    MPI_Type_commit(&mp_base_station_info);
    
    // Useful code
    int* arguments = (int*)args;
    int iterations = arguments[0], width = arguments[1], height = arguments[2];
    int size = width * height;
    printf("Iterations: %d", iterations);

    MPI_Request comparison_reqs[size];
    base_station_info comparison_buffer[size];
    int i;
    for (i = 0; i < size; i++) {
        MPI_Irecv(&comparison_buffer[i], 1, mp_base_station_info, i+1, 0,
                    MPI_COMM_WORLD, &comparison_reqs[i]);
    }

    // Set up the print file
    FILE *fptr;
    fptr = fopen("../logs.txt", "w");
    int iter;
    for(iter=0; iter<iterations; iter++){
        
        // Make everything below this point happen for each tsunameter that is sending
        int tsu;
        for(tsu=0; tsu<size; tsu++){
        
            if (test_mpi_req(&comparison_reqs[tsu])) {
                struct base_station_info reading = comparison_buffer[tsu];
                int sender_x = tsu % width, sender_y = tsu / width; 
                satellite_reading most_recent;
                int max_time = 0;
                int false_readings = 0, valid_readings = 0;
                int i;
                for (i=0; i<STORED_READINGS; i++){
                    if (satellite_readings[i].xpos == sender_x && satellite_readings[i].ypos == sender_y){
                        if (satellite_readings[i].time > max_time){
                            most_recent = satellite_readings[i];
                            max_time = most_recent.time;
                        }
                    }
                }
                if (max_time == 0) { // Then it's a false reading, as nothing was found
                    printf("False Reading from %d\n", tsu);
                    false_readings += 1;
                    // Log a false reading
                    // Tsunameter reading (time, node, neighbours, elevation)
                    // Flagged as False; No Matching Satellite Reading
                    struct tm tm = *localtime((time_t*)&reading.time);
                    fprintf(fptr, "X === Logging False Reading #%d from tsunameter %d\n", false_readings, tsu); 
                    fprintf(fptr, "Co-ordinates of reporting tsunameter is x: %d y: %d\n", sender_x, sender_y);
                    fprintf(fptr, "Time: %d-%02d-%02d %02d:%02d:%02d\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
                    fprintf(fptr, "At elevation %f\n", reading.avg);
                    fprintf(fptr, "Matching with neighbours:");
                    int q;
                    for(q = 0; q<4; q++){
                        if(reading.neighbours[q] != -1){
                            fprintf(fptr, " %d", reading.neighbours[q]);
                        }
                    } 
                    fprintf(fptr, "\n =====End of Event=====\n");
                    /*
                    printf("Avg is %f\n Time: %d\n Neighbours: %d %d %d %d\n", 
                    reading.avg, reading.time, reading.neighbours[0], reading.neighbours[1],
                    reading.neighbours[2], reading.neighbours[3]);
                    */
                    
                } else { // it's a valid reading
                    printf("Valid reading from %d\n", tsu);
                    valid_readings += 1;
                    // Log a valid reading
                    // Tsunameter Reading (time, node, neighbours, elevation)
                    // Satellite Reading (time, node, elevation)
                    // Time difference between
                    // Height difference.
                    struct tm tm = *localtime((time_t*)&reading.time);
                    fprintf(fptr, "O === Logging Valid Reading #%d from tsunameter %d\n", false_readings, tsu); 
                    fprintf(fptr, "Co-ordinates of reporting tsunameter is x: %d y: %d\n", sender_x, sender_y);
                    fprintf(fptr, "Time: %d-%02d-%02d %02d:%02d:%02d\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
                    fprintf(fptr, "At elevation %f\n", reading.avg);
                    fprintf(fptr, "Matching with neighbours:");
                    int q;
                    for(q = 0; q<4; q++){
                        if(reading.neighbours[q] != -1){
                            fprintf(fptr, " %d", reading.neighbours[q]);
                        }
                    }
                    fprintf(fptr, "\n\nMatching Satellite Reading\n");
                    tm = *localtime((time_t*)&most_recent.time);
                    fprintf(fptr, "Time: %d-%02d-%02d %02d:%02d:%02d\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
                    fprintf(fptr, "At elevation: %f\n", most_recent.elevation);
                    
                    fprintf(fptr, "\n\nReading Comparison:\n");
                    fprintf(fptr, "Time Difference: %d seconds\n", abs(most_recent.time - reading.time));
                    fprintf(fptr, "Elevation Difference: %f m\n", fabs(most_recent.elevation - reading.avg));   
                    fprintf(fptr, "=====End of Event=====\n");  
                                   
                    /*
                    printf("Avg is %f\n Time: %d\n Neighbours: %d %d %d %d\n", 
                    reading.avg, reading.time, reading.neighbours[0], reading.neighbours[1],
                    reading.neighbours[2], reading.neighbours[3]);
                    */
                }
                
                // Reset comparison request
                MPI_Irecv(&comparison_buffer[tsu], 1, mp_base_station_info, tsu+1, 0,
                    MPI_COMM_WORLD, &comparison_reqs[tsu]);
            }  
        }
        
        
        if(comms_terminate == 1){
            break;
        }
        sleep(TSUNAMETER_POLL);
    }
    int i;
    for(i=0; i<size; i++){
        MPI_Cancel(&comparison_reqs[i]);
    } 
    
    // Close the file handler
    fclose(fptr);
    fptr = NULL;
    
    // Send a termination broadcast to the tsunameters.
    MPI_Request send_req;
    int blah = 0;
    //printf("TERMINATING...\n");
    // Change the global POSIX terminate signal
    comms_terminate = 1;
    MPI_Ibcast(&blah, 1, MPI_INT, 0, MPI_COMM_WORLD, &send_req);
}















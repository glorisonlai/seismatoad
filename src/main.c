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
#include <netdb.h>
#include <arpa/inet.h>


#include "tsunameter.h" // Tsunameter functions
#include "satellite.h" // Satellite Readings struct

#define DIMENSIONS 2        // Tsunameter topology dimensions
#define TERMINATION_TAG 0   // MPI Tag for termination
#define MAX_READING 10000   // Max sample value
#define TOLERANCE 1000       // Tsunameter tolerance
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

    // Create MPI datatype for struct base_station_info 
    MPI_Datatype mp_base_station_info;
    int iblock_length[] = {1, 1, 4, 1}; // {float, int, int * neighbours, double}
    MPI_Aint iblock_displacement[] = {offsetof(base_station_info, avg),
                                    offsetof(base_station_info, time),
                                    offsetof(base_station_info, neighbours),
                                    offsetof(base_station_info, comm_time)};
    MPI_Datatype itypes[] = {MPI_FLOAT, MPI_INT, MPI_INT, MPI_DOUBLE};

    MPI_Type_create_struct(4, iblock_length, iblock_displacement, itypes,
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
        int sat_args[4] = {sat_records, tsunameter_dims[0], tsunameter_dims[1], tsunameter_threshold}; 
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
        
        
        
        //********************** Base termination checking ********************
        // Just loops until a termination signal is read.
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
        // Uncomment this code block if you need it:
        /*
        printf("Current State of satellite readings:\n");
        int i;
        for (i=0; i<sat_records; i++){
            printf("Record %d: Time is %d, Xpos is %d, Ypos is %d, elevation is %f\n",
                i, satellite_readings[i].time, satellite_readings[i].xpos, satellite_readings[i].ypos, satellite_readings[i].elevation);
                
        }
        */
        
        
        // terminate and clean up the altimeter
        pthread_join(satellite_tid, NULL);
        free(satellite_readings);

        // Also terminate and clean up the other two posix threads. 
        pthread_join(sentinel_tid, NULL);
        // Comms thread is last as join is blocking and comms takes the
        // longest to finalise
        pthread_join(comms_tid, NULL);
        
    } else {
    
    
/*+++++++++++++++++++++++++++++++ TSUNAMETERS +++++++++++++++++++++++++++++++*/


        srand(time(NULL) + tsunameter_rank); // Seed random gen for tsunameters

        // Instantiate moving average
        printf("rank %d\n", tsunameter_rank);
        moving_avg *avg = init_moving_avg(TSUNAMTER_WINDOW);

        // Get array of neighbours
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
            float new_val = generate_float_val(10000.0);
            append_moving_avg(avg, new_val);
            // Get current system time
            time_t curr_time = time(NULL);

            // Check if neighbours sent compare request signal
            printf("Rank: %d, Testing...\n", tsunameter_rank);
            for (nbr_i = 0; nbr_i < num_neighbours; nbr_i++) {
                if (test_mpi_req(&comparison_reqs[nbr_i])) {
                    struct tsunameter_reading *curr_reading =
                        instantiate_tsunameter_reading(get_moving_avg(avg), curr_time);

                    // Send information
                    MPI_Send(curr_reading, 1, mp_tsunameter_reading, neighbours[nbr_i], 2*neighbours[nbr_i] + 1,
                            tsunameter_comm);
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
                    MPI_Irecv(&recv_buff[nbr_i], 1, mp_tsunameter_reading, neighbours[nbr_i], 2*tsunameter_rank + 1,
                                tsunameter_comm, &recv_reqs[nbr_i]);

                }

                // Compare received readings with itself
                // Timeout after 10 s
                int similar_count, in_count;
                int similar_neighbours[4] = {-1,-1,-1,-1};
                while (MPI_Wtime() - starttime < 10) {
                    // Repeatedly check for similar readings
                    
                    // Check if we're terminating to avoid waiting
                    if(test_mpi_req(&termination_req)){
                        break;
                    }
                    
                    // Instantiate counters to potentially send to base station
                    similar_count = 0;
                    in_count = 0;
                    similar_neighbours[0] = -1;
                    similar_neighbours[1] = -1;
                    similar_neighbours[2] = -1;
                    similar_neighbours[3] = -1;
                    for (nbr_i = 0; nbr_i < num_neighbours; nbr_i++) {
                        // Continue checking if neighbours sent compare request
                        if (test_mpi_req(&comparison_reqs[nbr_i])) {
                            struct tsunameter_reading *curr_reading =
                                instantiate_tsunameter_reading(get_moving_avg(avg),
                                                                curr_time);
                            // Send information
                            MPI_Send(curr_reading, 1, mp_tsunameter_reading, neighbours[nbr_i], 2*neighbours[nbr_i] + 1,
                                    tsunameter_comm);
                            free(curr_reading);

                            // Reset comparison request
                            MPI_Irecv(&comparison_buffer[nbr_i], 1, MPI_INT, neighbours[nbr_i], 2*tsunameter_rank,
                                        tsunameter_comm, &comparison_reqs[nbr_i]);
                        }

                        // If received message, increment in_count
                        // If message is similar to our reading, increment similar_count
                        if (test_mpi_req(&recv_reqs[nbr_i])) {
                            in_count += 1;
                            if (fabsf(recv_buff[nbr_i].avg - get_moving_avg(avg)) <
                                TOLERANCE) {
                                similar_neighbours[nbr_i] = neighbours[nbr_i];
                                similar_count += 1;
                            }
                        }
                    }

                    // Send information to base station if enough similar readings
                    if (similar_count >= 2) {
                        // Instantiate base_staion_info
                        printf("Rank: %d sending to base station\n", tsunameter_rank);
                        struct base_station_info base_station_buf;
                        base_station_buf.avg = get_moving_avg(avg);
                        printf("Avg is %f\n", base_station_buf.avg);
                        base_station_buf.time = time(NULL);
                        printf("Time: %d\n", base_station_buf.time);
                        base_station_buf.neighbours[0] = similar_neighbours[0];
                        base_station_buf.neighbours[1] = similar_neighbours[1];
                        base_station_buf.neighbours[2] = similar_neighbours[2];
                        base_station_buf.neighbours[3] = similar_neighbours[3];
                        printf("Neighbours: %d %d %d %d\n", base_station_buf.neighbours[0], base_station_buf.neighbours[1],
                    base_station_buf.neighbours[2], base_station_buf.neighbours[3]);
                        base_station_buf.comm_time =  MPI_Wtime() - starttime;
                        MPI_Request tempstat;
                        // Send to base station
                        MPI_Isend(&base_station_buf, 1, mp_base_station_info, root, 0,
                                MPI_COMM_WORLD, &tempstat);
                        printf("Sending completed by %d\n", tsunameter_rank);
                         
                        break;
                    }

                    // If all neighbours have sent, break
                    if (in_count == num_neighbours) {
                        break;
                    }
                    sleep(1);
                }
            }
            
            sleep(TSUNAMETER_POLL);
        }

        // Termination signal received
        free_moving_avg(avg);
        free(neighbours);
        printf("%d terminating...\n", base_rank);
    } 
/*+++++++++++++++++++++++++++++++ CLEAN UP +++++++++++++++++++++++++++++++*/
    MPI_Comm_free(&tsunameter_comm);
    printf("Freed by %d\n", base_rank);
    return MPI_Finalize();
}



/*+++++++++++++++++++++++++++++ SATELLITE THREAD +++++++++++++++++++++++++++++*/
/**
 * @param void* args -> Integer array of size 4 cast to void*
 * @return -> function is void
 */
void* run_satellite(void* args){
    // arguments should be array size, then width and height of grid, then threshold
    printf("Satellite Thread starting\n");

    int last_access = 0; // keeps track of whhere to put the next element
    int xpos, ypos, elevation;
    int* arguments = (int*)args;
    int size = arguments[0], width = arguments[1], height = arguments[2], threshold = arguments[3];
    do{
        // Make this properly FIFO, idk how tbh I think this is probably good enough.
        sleep(1);
        // generate a random co-ordinate
        srand(time(NULL));
        xpos = rand() % width;
        ypos = rand() % height;
        elevation = (rand() % 1000) + threshold;
        
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
/**
 * Periodically check input for a sentinel character, commence termination if seen
 * @param void* arg -> The sentinel character (or string)
 * @return -> function is void
 */
void* check_sentinel(void* arg){
    printf("Sentinel Thread starting\n");
    char buf[strlen(arg)]; // Construct a buffer of equal length to sentinel
    char* sentinel = arg;  // get the sentinel character
    
    // Loop until termination signal or sentinel character entered
    do{
        // Set up the IO so that it is non-blocking using STDIN_FILENO as an
        // intermediary for inputs
        fcntl(STDIN_FILENO, F_SETFL, fcntl(STDIN_FILENO, F_GETFL) | O_NONBLOCK);
        sleep(1);
        int readChars = read(STDIN_FILENO, buf, 4); // read into buf and check characters
        if (readChars > 0) {

            buf[strcspn(buf, "\n")] = 0; // strip the newline character
            printf("Read in: %s\n", buf);
            // check whether buf now matches the sentinel value
            if(strcmp(buf, sentinel) == 0){
                // Mutex may be unnecessary here as both the main and the sentinel threads will only ever set it to 1
                // However it is better to be safe to avoid a failing terminate.
                pthread_mutex_lock(&sentinel_mutex);
                sentinel_terminate = 1;
                pthread_mutex_unlock(&sentinel_mutex);
            } else{
                // If the wrong sentinel value is entered, remind of the right one
                printf("Enter %s to exit\n", sentinel);
            }
        }
    } while(sentinel_terminate == 0);
    // Terminate when a signal is sent in/out
    printf("Sentinel thread terminating\n");
}

/*+++++++++++++++++++++++++++++ COMMS THREAD +++++++++++++++++++++++++++++*/

/**
 * Handle the Messaging between Tsunameters and Base Station
 * @param void* args -> integer array of size 3 with iterations, width and height 
 * @return -> function is void
 */
void* run_comms(void* args){
    printf("Comms Thread Starting\n");
    // Start timing for log
    double start_time = MPI_Wtime();
    
    // Setup type to receive
    MPI_Datatype mp_base_station_info;
    int iblock_length[] = {1, 1, 4, 1}; // {float, int, int * neighbours, double}
    MPI_Aint iblock_displacement[] = {offsetof(base_station_info, avg),
                                    offsetof(base_station_info, time),
                                    offsetof(base_station_info, neighbours),
                                    offsetof(base_station_info, comm_time)};
    MPI_Datatype itypes[] = {MPI_FLOAT, MPI_INT, MPI_INT, MPI_DOUBLE};

    MPI_Type_create_struct(4, iblock_length, iblock_displacement, itypes,
                            &mp_base_station_info);
    MPI_Type_commit(&mp_base_station_info);

    // Handle the inputs to the function
    int* arguments = (int*)args;
    int iterations = arguments[0], width = arguments[1], height = arguments[2];
    int size = width * height;
    // printf("Iterations: %d\n", iterations);

    // Setup the array of requests, and a location for the sent info to arrive at
    MPI_Request comparison_reqs[size];
    base_station_info comparison_buffer[size];
    int i;
    for (i = 0; i < size; i++) {
        MPI_Irecv(&comparison_buffer[i], 1, mp_base_station_info, i+1, 0,
                    MPI_COMM_WORLD, &comparison_reqs[i]);
    }

    // Set up the print file
    FILE *fptr;
    fptr = fopen("logs.txt", "w");
    
    // Establish variables for later logging
    int iter;
    int false_readings = 0, valid_readings = 0;
    double total_comm_time = 0;

    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int processor_len;
    MPI_Get_processor_name(processor_name, &processor_len);
    fprintf(fptr, "Processor name: %s\n", processor_name);

    char hostbuffer[256];
    char *IPbuffer;
    struct hostent *host_entry;
    int hostname;
    // To retrieve hostname
    hostname = gethostname(hostbuffer, sizeof(hostbuffer));
    // To retrieve host information
    host_entry = gethostbyname(hostbuffer);
    // To convert an Internet network
    // address into ASCII string
    IPbuffer = inet_ntoa(*((struct in_addr*)
                           host_entry->h_addr_list[0]));
    fprintf(fptr, "IPv4: %s\n", IPbuffer);
    
    // Do the main loop once per specified iteration
    for(iter=0; iter<iterations; iter++){
        
        // Make everything below this point happen for each tsunameter that is sending
        int tsu;
        for(tsu=0; tsu<size; tsu++){
            
            // Check if an event has been received from the tsunameter
            if (test_mpi_req(&comparison_reqs[tsu])) {
                // Access just the relevant reading
                struct base_station_info reading = comparison_buffer[tsu];
                // calculate sender co-ordinates
                int sender_x = tsu % width, sender_y = tsu / width; 
                satellite_reading most_recent;
                int max_time = 0;
                
                int i;
                // Find the most recent matching satellite reading (if there is one)
                for (i=0; i<STORED_READINGS; i++){
                    if (satellite_readings[i].xpos == sender_x && satellite_readings[i].ypos == sender_y){
                        if (satellite_readings[i].time > max_time){
                            most_recent = satellite_readings[i];
                            max_time = most_recent.time;
                        }
                    }
                }
                // Log the reading:
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
                    fprintf(fptr, "Communication time taken: %.3f\n", reading.comm_time);
                    fprintf(fptr, "Matching with neighbours:");
                    int q;
                    for(q = 0; q<4; q++){
                        if(reading.neighbours[q] != -1){
                            fprintf(fptr, " %d", reading.neighbours[q]);
                        }
                    } 
                    fprintf(fptr, "\n =====End of Event=====\n\n");
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
                    fprintf(fptr, "O === Logging Valid Reading #%d from tsunameter %d\n", valid_readings, tsu); 
                    fprintf(fptr, "Co-ordinates of reporting tsunameter is x: %d y: %d\n", sender_x, sender_y);
                    fprintf(fptr, "Time: %d-%02d-%02d %02d:%02d:%02d\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
                    fprintf(fptr, "At elevation %f\n", reading.avg);
                    fprintf(fptr, "Communication time taken: %.3f\n", reading.comm_time);
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
                    fprintf(fptr, "=====End of Event=====\n\n");  
                                   
                    /*
                    printf("Avg is %f\n Time: %d\n Neighbours: %d %d %d %d\n", 
                    reading.avg, reading.time, reading.neighbours[0], reading.neighbours[1],
                    reading.neighbours[2], reading.neighbours[3]);
                    */
                }
                // Update total running communication time
                total_comm_time += reading.comm_time;

                // Reset comparison request
                MPI_Irecv(&comparison_buffer[tsu], 1, mp_base_station_info, tsu+1, 0,
                    MPI_COMM_WORLD, &comparison_reqs[tsu]);
            }  
        }
        
        // Handle termination message from base station:
        if(comms_terminate == 1){
            break;
        }
        sleep(TSUNAMETER_POLL);
    }
    // Cancel requests, just in case
    for(i=0; i<size; i++){
        MPI_Cancel(&comparison_reqs[i]);
    } 
    // Generate summary:
    double end_time = MPI_Wtime();
    fprintf(fptr, "\n\n\n============== SUMMARY ===============\n");
    fprintf(fptr, "Had %d total events, %d were valid, %d false\n", false_readings + valid_readings, valid_readings, false_readings);
    fprintf(fptr, "Total time taken in seconds: %.3lf\n", end_time - start_time);
    fprintf(fptr, "Time per communication (avg): %.3lf\n", total_comm_time / (valid_readings + false_readings));
    
    // Close the file handler
    fclose(fptr);
    fptr = NULL;
    
    // Send a termination broadcast to the tsunameters.
    MPI_Request send_req;
    int termination_buf = 0;
    //printf("TERMINATING...\n");
    // Change the global POSIX terminate signal
    comms_terminate = 1;
    MPI_Ibcast(&termination_buf, 1, MPI_INT, 0, MPI_COMM_WORLD, &send_req);
}


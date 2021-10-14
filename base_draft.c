// start by initialising the tsunameters
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
#include <string.h> // necessary for strncmp to avoid empty buffer failures
#include <pthread.h>


/* Global Variables for POSIX Use */
int *satellite_readings;
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
void* run_comms(void* arg);


int main(int argc, char** argv){
    // Launch arguments: grid x, grid y, iterations, sentinel value
    int my_rank=0;
    int tsunameter_count = 9;
    int iterations = atoi(argv[2]);
    int terminating = 0;
    int height=2, width=2;
    printf("Starting up\n");
    if (my_rank == 0){
        // First initialise a thread to simulate the altimeter
            //setup the satellite readings array first
        printf("Creating Satellite Thread\n");
        int sat_records = 10, record_size = 4;
        int sat_args[3] = {sat_records * record_size, height, width}; 
        satellite_terminate = 0;
        satellite_readings = (int *)malloc((sat_records * record_size) * sizeof(int));
        pthread_t satellite_tid;
        pthread_mutex_init(&satellite_mutex, NULL);
        pthread_create(&satellite_tid, NULL, run_satellite, &sat_args);

        // Construct an additional thread to handle sentinel exit values
        // POSIX HERE
        printf("Creating Sentinel Thread\n");
        pthread_t sentinel_tid;
        pthread_mutex_init(&sentinel_mutex, NULL);
        sentinel_terminate = 0;
        printf("%s\n", argv[4]);
        pthread_create(&sentinel_tid, NULL, check_sentinel, argv[4]);


        // For HD+ this section needs to be a POSIX thread to protect against timeouts.
        // Then wait for an alert from one of the nodes
        printf("Creating Comms Thread\n");
        pthread_t comms_tid;
        pthread_mutex_init(&comms_mutex, NULL);
        pthread_create(&comms_tid, NULL, run_comms, NULL);
        
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
        
        printf("Current State of satellite readings:\n");
        for (int i=0; i<4*sat_records; i+=4){
            printf("Record %d: Time is %d, Xpos is %d, Ypos is %d, elevation is %d\n",
                i/4, satellite_readings[i], satellite_readings[i + 1], satellite_readings[i + 2],
                satellite_readings[i + 3]);
        }
        // terminate and clean up the altimeter
        pthread_join(satellite_tid, NULL);
        free(satellite_readings);

        // Do a final wait all to make sure all finalised
        //MPI_Wait_All();
        // Also terminate and clean up the other two posix threads. Both the MPI Send Recv and the Input scanning
        pthread_join(sentinel_tid, NULL);
        pthread_join(comms_tid, NULL);
    }
    return 0;
}





void* run_satellite(void* args){
    // arguments should be array size, then width and height of grid
    printf("Satellite Thread starting\n");

    int last_access = 0;
    int xpos, ypos, elevation;
    int* arguments = (int*)args;
    int size = arguments[0], width = arguments[1], height = arguments[2];
    do{
        // Make this properly FIFO, idk how tbh I think this is probably good enough.
        sleep(2);
        // generate a random co-ordinate
        srand(time(NULL));
        xpos = rand() % width;
        ypos = rand() % height;
        elevation = (rand() % 1000) + 6000;
        last_access = (last_access + 4) % size;

        // Implement some quick Mutex here to make sure we avoid race conditions
        pthread_mutex_lock(&satellite_mutex);
        satellite_readings[last_access] = time(NULL);
        satellite_readings[last_access + 1] = xpos;
        satellite_readings[last_access + 2] = ypos;
        satellite_readings[last_access + 3] = elevation;
        // Release mutex lock.
        pthread_mutex_unlock(&satellite_mutex);
        
    } while (satellite_terminate == 0);
    printf("Satellite Thread terminating\n");
}

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
        sleep(2);
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

void* run_comms(void* arg){
    printf("Comms Thread Starting\n");
    /*
        int count = 0;
        do {
            MPI_Recv();
            // Do a non-blocking receive from every single tsunameter
            // Check whether any of them have responded
            // If any of them have:

                // Once something is received:
                // iterate through the sattelite_readings to find one for that tsunameter
                // if none exist, log as false alert
                // if one does exist, check timestamp and verify it is within the same 10 second period
                // if more than one does exist, check the most recent timestamp
                // log the alert as a valid alert with as much info as possible

            // There is a set amount of iterations to run, specified at run time
            // check for whether a sentinel value input has been entered

        } while((count < iterations) && (terminating == 0));

        // to terminate, do a final checking receive for every tsunameter
        // This frees up any hanging sends
        for(int i=0; i<tsunameter_count; i++){
            MPI_Irecv(i);
        }
        MPI_Broadcast("terminate");
        */
}
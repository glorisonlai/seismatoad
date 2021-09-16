#include <stdio.h>  // Printing
#include <mpi.h>    // OPENMPI
#include <stdlib.h> // Type conversion 
#include <errno.h>  // ERROR CHECKING
#include <limits.h> // INT LIMITS

int main(int argc, char** argv)
{
    int m, n;
    if (argc != 3) {
        printf("\nUsage: ./main m n\n");
        return 1;
    }

    errno = 0;
    int arg_index;
    for (arg_index = 1; arg_index < 3; arg_index++) {
        char* p;
        long arg = strtol(argv[arg_index], &p, 10);
        if (errno != 0 || *p != '\0' || arg < INT_MIN || arg > INT_MAX) {
            // Handle error
            printf("\nArgument error\n");
            return 1;
        }
        switch (arg_index)
        {
            case 1:
                m = arg;
                break;
            case 2:
                n = arg;
                break;
            default:
                printf("\nAssignment error %d\n", arg_index);
                return 1;
        }
    }

    printf("m:%d, n:%d\n", m, n);
    return 0;
}

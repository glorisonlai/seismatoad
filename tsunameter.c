// Tsunameter functions
#include <stdlib.h>
#include <time.h>
#include "tsunameter.h"

moving_avg* init_moving_avg(int size) {
    moving_avg* avg = (moving_avg *) malloc(sizeof(moving_avg));
    avg -> max_size = size;
    return avg;
}

void append_moving_avg(moving_avg* avg, float val) {
    avg_el* new_el = (avg_el *) malloc(sizeof(avg_el));
    new_el -> value = val;
    if (avg -> size == avg -> max_size) {
        avg_el* temp = avg -> queue_head;
        avg -> queue_head = temp -> next;
        avg -> queue_tail -> next = new_el;
        avg -> queue_tail = new_el;
        avg -> avg = (float)(avg -> avg * avg -> size - temp->value + new_el -> value) / avg -> size;
        free(temp);
    } else if (avg -> size > 0) {
        avg -> queue_tail -> next = new_el;
        avg -> queue_tail = new_el;
        avg -> size += 1;
        avg -> avg = (float)(avg -> avg * (avg -> size - 1) + new_el -> value) / avg -> size;
    } else {
        avg -> queue_head = new_el;
        avg -> queue_tail = new_el;
        avg -> avg = new_el -> value;
        avg -> size = 1;
    }
}

float generate_float_val(float limit) {
    return (float)rand()/(float)(RAND_MAX) * limit;
}
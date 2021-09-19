// Interface for Tsunameter functions
#ifndef TSUNAMETER_H
#define TSUNAMETER_H

typedef struct avg_el{
    float value;
    struct avg_el *next;
} avg_el;

typedef struct moving_avg {
    float avg;
    int max_size;
    int size;
    avg_el* queue_head;
    avg_el* queue_tail;
} moving_avg;

// Initialize moving average
moving_avg* init_moving_avg(int);

// Append float to moving average and update average
void append_moving_avg(moving_avg*, float);

// Generate random float up to limit
float generate_float_val(float);

#endif
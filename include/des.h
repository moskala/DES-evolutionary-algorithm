#pragma once

#include <stdbool.h>
#include <stdint.h>

struct result {
    double best_fit;
    double* best_result;
    int count;
};

struct result des(int N, double function_fn(int N, double[]), double lower[], double upper[], uint64_t seed,
    bool logRes);

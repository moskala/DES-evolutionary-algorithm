#pragma once
#include <stdbool.h>

struct result {
    double best_fit;
    double *best_result;
    int count;
};

struct result des(int N, double function_fn(int N, double[N]), double lower[N], double upper[N], int seed, bool logRes);

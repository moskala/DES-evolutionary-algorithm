#pragma once


struct result {
    double best_fit;
    double *best_result;
    int count;
};

struct result des(int N, double initial_point[N], double function_fn(int N, double[N]), double lower[N], double upper[N]);

double fun_sin_cos(int n, double* v);
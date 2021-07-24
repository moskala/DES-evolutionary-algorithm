#pragma once
#include <stdbool.h>
#include <stdio.h>

void test_function(char* testName, int nTimes, int seed, double expectedValue, int N,
    double function_fn(int N, double[N]), double lower[N], double upper[N], bool logDes);
void test_function_with_statistics(FILE* error_file, FILE* time_file, char* test_name, int n_times, int seed,
    double expected_value, int N, double function_fn(int N, double[N]), double lower[N],
    double upper[N]);

double fun_sin_cos(int n, double* v);
void test_simple_quadratic(FILE* error_file, FILE* time_file, int nTimes, int seed);

double fun_simple_quadratic(int n, double* v);
void test_sin_cos(FILE* error_file, FILE* time_file, int nTimes, int seed);

double fun_Ackleya(int n, double* x);
void test_Ackeleya(FILE* error_file, FILE* time_file, int nTimes, int seed, int nDim);

double fun_Rastrigin(int n, double* x);
void test_Rastrigin(FILE* error_file, FILE* time_file, int nTimes, int seed, int dim);

double fun_Shubert(int n, double* x);
void test_Shubert(FILE* error_file, FILE* time_file, int nTimes, int seed);

double fun_Shekel(int n, double* x);
void test_Shekel(FILE* error_file, FILE* time_file, int nTimes, int seed);

double fun_Griewank(int n, double* x);
void test_Griewank(FILE* error_file, FILE* time_file, int nTimes, int seed, int dim);

double fun_Perm(int n, double* x);
void test_Perm(FILE* error_file, FILE* time_file, int nTimes, int seed, int dim);

double fun_rotated(int n, double* x);
void test_rotated(FILE* error_file, FILE* time_file, int nTimes, int seed, int dim);

double fun_Zakharov(int n, double* x);
void test_Zakharov(FILE* error_file, FILE* time_file, int nTimes, int seed, int dim);

#pragma once
#include <stdbool.h>

void test_function(char *testName, int nTimes, int seed, double expectedValue, int N, double function_fn(int N, double[N]), double lower[N], double upper[N], bool logDes);
void test_function_with_statistics(char *test_name, int n_times, int seed, double expected_value, int N, double function_fn(int N, double[N]), double lower[N], double upper[N]);



double fun_sin_cos(int n, double* v);
void test_simple_quadratic(int nTimes, int seed);

double fun_simple_quadratic(int n, double* v);
void test_sin_cos(int nTimes, int seed);

double fun_Ackleya(int n, double *x);
void test_Ackeleya(int nTimes, int seed, int nDim);


double fun_Rastrigin(int n, double *x);
void test_Rastrigin(int nTimes, int seed, int dim);


double fun_Schubert(int n, double *x);
void test_Schubert(int nTimes, int seed);

double fun_Shekel(int n, double *x);
void test_Shekel(int nTimes, int seed);


double fun_Griewank(int n, double *x);
void test_Griewank(int nTimes, int seed, int dim);

double fun_Perm(int n, double *x);
void test_Perm(int nTimes, int seed, int dim);

double fun_rotated(int n, double *x);
void test_rotated(int nTimes, int seed, int dim);

double fun_Zakharov(int n, double *x);
void test_Zakharov(int nTimes, int seed, int dim);


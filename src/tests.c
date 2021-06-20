#include "test_funtions.h"
#include "des.h"
#include <stdio.h>


void print_result(int N, struct result res) {
    printf("Result - Fit: %f, Count: %d\n", res.best_fit, res.count);
    
    printf("Best solution:\n");
    for(int i = 0; i < N; ++i)
    {
        printf("%f, ", res.best_result[i]);
    }
    printf("\n");
}

void print_result_and_expected(int N, struct result res, double expectedValue) {
    printf("Result - Fit: %f, Expected: %f\n", res.best_fit, expectedValue);
}

void test_function(char *testName, int nTimes, int seed, double expectedValue, int N, double function_fn(int N, double[N]), double lower[N], double upper[N], bool logDesResults) {
    printf("Test: %s\n", testName);
    int current_seed = seed;
    for(int i = 0; i < nTimes; ++i){
        struct result res = des(N, NULL, function_fn, lower, upper, current_seed, logDesResults);
        printf("Test %d / %d\n", i+1, nTimes);
        print_result_and_expected(N, res, expectedValue);
        current_seed *= 2;
    }
}


void test_simple_quadratic(int nTimes, int seed) {
    char *name = "Test simple quandric";
    int N = 1;
    double init_point[N];
    double lower[1] = {-5};
    double upper[1] = {5};
    test_function(name, nTimes, seed, -9.0, 1, fun_simple_quadratic, lower, upper, false);
}

void test_sin_cos(int nTimes, int seed) {
    char *name = "Test sin*cos";
    int N = 2;
    double init_point[N];
    double lower[2] = {-5, -5};
    double upper[2] = {5, 5};
    double expected_value = -1;
    test_function(name, nTimes, seed, expected_value, 2, fun_sin_cos, lower, upper, false);
}

// Minimum for x = 0, f(x) = 0
void test_Ackeleya(int nTimes, int seed, int dim){
    char *name = "Test Ackeley";
    int N = dim;
    double lower[N];
    double upper[N];
    for(int i = 0; i < N; ++i)
    {
        lower[i] = -30;
        upper[i] = 30;

    }
    double expected_value = 0;
    test_function(name, nTimes, seed, expected_value, N, fun_Ackleya, lower, upper, false);
}

// Minimum for x = 0, f(x) = -n
void test_Rastrigin(int nTimes, int seed, int dim){
    char *name = "Test Rastrigin";
    int N = dim;
    double lower[N];
    double upper[N];
    for(int i = 0; i < N; ++i)
    {
        lower[i] = -1;
        upper[i] = 1;

    }
    double expected_value = -N;
    test_function(name, nTimes, seed, expected_value, N, fun_Rastrigin, lower, upper, false);
}

// Minimum f(x) = 186.73 for 760 different xs
void test_Schubert(int nTimes, int seed){
    char *name = "Test Schubert";
    int N = 2;
    double lower[N];
    double upper[N];
    for(int i = 0; i < N; ++i)
    {
        lower[i] = -5.12;
        upper[i] = 5.12;

    }
    double expected_value = -186.7309;
    test_function(name, nTimes, seed, expected_value, N, fun_Schubert, lower, upper, false);
}

// x = (4,4,4,4) f(x) = -10.5364
void test_Shekel(int nTimes, int seed){
    char *name = "Test Shekel";
    int N = 4;
    double lower[N];
    double upper[N];
    for(int i = 0; i < N; ++i)
    {
        lower[i] = 0;
        upper[i] = 10;

    }
    double expected_value = -10.5364;
    test_function(name, nTimes, seed, expected_value, N, fun_Shekel, lower, upper, false);
}
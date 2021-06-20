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


void test_simple_quadratic() {
    printf("Test simple quandric\n");
    int N = 1;
    double init_point[N];
    double lower[2] = {-5, -5};
    double upper[2] = {5, 5};
    struct result res = des(N, init_point, fun_simple_quadratic, lower, upper);
    print_result(N, res);
}

void test_sin_cos() {
    printf("Test sin*cos\n");
    int N = 2;
    double init_point[N];
    double lower[2] = {-5, -5};
    double upper[2] = {5, 5};
    struct result res = des(N, init_point, fun_sin_cos, lower, upper);
    print_result(N, res);
}
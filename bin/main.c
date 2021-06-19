#include "des.h"
#include <stdio.h>


double f(int n, double* v)
{
    // return (*x - 2) * (*x + 4);
    double x = v[0];
    double y = v[1];
    return sin(x + y) * cos(x - y);
}

// struct result des(int N, double initial_point[N], double function_fn(int N, double[N]), double lower[N], double upper[N])
int main() {
    fprintf(stderr, "here: %d\n", __LINE__);
    int N = 2;
    double init_point[N];
    double lower[2] = {-5, -5};
    double upper[2] = {5, 5};
    struct result res = des(N, init_point, f, lower, upper);
    printf("Result: %f, %f, Fit: %f, Count: %d", res.best_result[0], res.best_result[1],  res.best_fit, res.count );
}

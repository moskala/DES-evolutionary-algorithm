#include "des.h"
#include <stdio.h>




// struct result des(int N, double initial_point[N], double function_fn(int N, double[N]), double lower[N], double upper[N])
int main() {
    fprintf(stderr, "here: %d\n", __LINE__);
    int N = 2;
    double init_point[N];
    double lower[2] = {-5, -5};
    double upper[2] = {5, 5};
    struct result res = des(N, init_point, fun_sin_cos, lower, upper);
    printf("Result: %f, %f, Fit: %f, Count: %d\n", res.best_result[0], res.best_result[1],  res.best_fit, res.count );
}

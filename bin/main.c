#include "des.h"



double f(int n, double* x)
{
    return (*x - 2) * (*x + 2);
}

// struct result des(int N, double initial_point[N], double function_fn(int N, double[N]), double lower[N], double upper[N])
int main() {

    int N = 1;
    double init_point[N];
    double lower[1] = {-5};
    double upper[1] = {5};
    struct result res = des(N, init_point, f, lower, upper);
    printf("Result: %f, Fit: %f, Count: %d", res.best_result[0], res.best_fit, res.count );
}

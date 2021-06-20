#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>


double fun_sin_cos(int n, double* v)
{
    double x = v[0];
    double y = v[1];
    return sin(x + y) * cos(x - y);
}


double fun_simple_quadratic(int n, double* v)
{
    double x = *v;
    return (x - 2) * (x + 4);
    
}
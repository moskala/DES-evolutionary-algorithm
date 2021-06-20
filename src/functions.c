#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>
#include <stdio.h>


double fun_sin_cos(int n, double* v)
{
    if(n != 2)
    {
        printf("Wrong dim for sin_cos function!");
    }
    double x = v[0];
    double y = v[1];
    return sin(x + y) * cos(x - y);
}


double fun_simple_quadratic(int n, double* v)
{
    if(n != 1)
    {
        printf("Wrong dim for quadratic function!");
    }
    double x = *v;
    return (x - 2) * (x + 4);
    
}

double fun_Ackleya(int n, double *x) 
{
    double sum_pow = 0.0;
    double sum_cos = 0.0;

    for(int i = 0; i < n; ++i)
    {
        sum_pow += pow(x[i], 2);
        sum_cos += cos(2 * M_PI * x[i]);
    }
    sum_pow /= n;
    sum_cos /= n;

    return -20 * exp(-0.2 * sqrt(sum_pow)) - exp(sum_cos) + 20 + M_E;
}

double fun_Rastrigin(int n, double *x)
{
    double sum = 0.0;
    for(int i = 0; i < n; ++i)
    {
        sum += pow(x[i], 2) - cos(18*x[i]);
    }
    return sum;
}


double fun_Schubert(int n, double *x)
{
    if(n != 2)
    {
        printf("Wrong dim for Schubert function!");
    }
    int m = 5;
    double sum_x1 = 0.0;
    double sum_x2 = 0.0;
    for(int i = 1; i < m+1; ++i)
    {
        sum_x1 += i * cos((i+1)*x[0] + i);
        sum_x2 += i * cos((i+1)*x[1] + i);
    }
    return sum_x1 * sum_x2;
}


double fun_Shekel(int n, double *x)
{
    if(n != 4)
    {
        printf("Wrong dim for Shekel function!");
    }
    
    double sum = 0.0;
    double A[4][10] = {
        {4.0, 1.0, 8.0, 6.0, 3.0, 2.0, 5.0, 8.0, 6.0, 7.0},
        {4.0, 1.0, 8.0, 6.0, 7.0, 9.0, 3.0, 1.0, 2.0, 3.6},
        {4.0, 1.0, 8.0, 6.0, 3.0, 2.0, 5.0, 8.0, 6.0, 7.0},
        {4.0, 1.0, 8.0, 6.0, 7.0, 9.0, 3.0, 1.0, 2.0, 3.6}};
    double c[10] = {0.1, 0.2, 0.2, 0.4, 0.6, 0.6, 0.3, 0.7, 0.5, 0.5};
    
    for(int i = 0; i < 10; ++i)
    {
        double sum_pow = 0.0;
        for(int j = 0; j < 4; ++j)
        {
            sum_pow += pow(x[j] - A[j][i], 2);
        }
        sum += 1 / (sum_pow + c[i]);
    }
    return (-1) * sum;
}
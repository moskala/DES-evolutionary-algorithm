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
    double sum = 10.0 * n;
    for(int i = 0; i < n; ++i)
    {
        sum += pow(x[i], 2) - 10 * cos(2 * M_PI * x[i]);
    }
    return sum;
}


double fun_Shubert(int n, double *x)
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


double fun_Griewank(int n, double *x)
{
    double sum = 0.0;
    double prod = 1.0;
    for(int i = 0; i < n; ++i)
    {
        sum += pow(x[i], 2.0) / 4000;
        prod *= cos(x[i] / sqrt(i+1));
    }
    return sum - prod + 1;
}

double fun_Perm(int n, double *x)
{
    double sum = 0.0;
    double beta = 10;
    for(int i = 1; i < n+1; ++i)
    {
        double inner_sum = 0.0;
        for(int j  = 1; j < n+1; ++j)
        {     
            inner_sum += (j + beta) * (pow(x[j-1], i) - pow(1.0 / j, i));
        }
        sum += pow(inner_sum, 2.0);
    }
    return sum;
}


double fun_rotated(int n, double *x)
{
    double sum = 0.0;
    for(int i = 0; i < n; ++i)
    {
        double inner_sum = 0.0;
        for(int j  = 0; j < i; ++j)
        {     
            inner_sum += pow(x[j], 2);
        }
        sum += inner_sum;
    }
    return sum;
}

double fun_Zakharov(int n, double *x)
{
    double sum1 = 0.0;
    double sum2 = 0.0;
    for(int i = 0; i < n; ++i)
    {
        sum1 += pow(x[i], 2.0);
        sum2 += 0.5 * (i + 1) * x[i];
    }
    return sum1 + pow(sum2, 2.0) + pow(sum2, 4.0);
}









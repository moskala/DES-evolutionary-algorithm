#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "des.h"
#include "test_funtions.h"

void test_dim5(int seed, int n_times)
{
    int dim = 5;

    FILE* error_file = fopen("error_dim5.txt", "w");
    fprintf(error_file, "Name;N;Best;Worst;Mean;Median;Standard Deviation\n");
    FILE* times_file = fopen("times_dim5.txt", "w");
    fprintf(times_file, "Name;N;Best;Worst;Mean;Median;Standard Deviation\n");

    test_Ackeleya(error_file, times_file, n_times, seed, dim);
    test_Rastrigin(error_file, times_file, n_times, seed, dim);
    test_Shubert(error_file, times_file, n_times, seed);
    test_Shekel(error_file, times_file, n_times, seed);
    test_Griewank(error_file, times_file, n_times, seed, dim);
    test_Perm(error_file, times_file, n_times, seed, dim);
    test_rotated(error_file, times_file, n_times, seed, dim);
    test_Zakharov(error_file, times_file, n_times, seed, dim);
    test_simple_quadratic(error_file, times_file, n_times, seed);
    test_sin_cos(error_file, times_file, n_times, seed);

    fclose(error_file);
    fclose(times_file);
}

void test_dim10(int seed, int n_times)
{
    int dim = 10;

    FILE* error_file = fopen("error_dim10.txt", "w");
    fprintf(error_file, "Name;N;Best;Worst;Mean;Median;Standard Deviation\n");
    FILE* times_file = fopen("times_dim10.txt", "w");
    fprintf(times_file, "Name;N;Best;Worst;Mean;Median;Standard Deviation\n");

    test_Ackeleya(error_file, times_file, n_times, seed, dim);
    test_Rastrigin(error_file, times_file, n_times, seed, dim);
    test_Shubert(error_file, times_file, n_times, seed);
    test_Shekel(error_file, times_file, n_times, seed);
    test_Griewank(error_file, times_file, n_times, seed, dim);
    test_Perm(error_file, times_file, n_times, seed, dim);
    test_rotated(error_file, times_file, n_times, seed, dim);
    test_Zakharov(error_file, times_file, n_times, seed, dim);
    test_simple_quadratic(error_file, times_file, n_times, seed);
    test_sin_cos(error_file, times_file, n_times, seed);

    fclose(error_file);
    fclose(times_file);
}

void test_dim20(int seed, int n_times)
{
    int dim = 20;

    FILE* error_file = fopen("error_dim20.txt", "w");
    fprintf(error_file, "Name;N;Best;Worst;Mean;Median;Standard Deviation\n");
    FILE* times_file = fopen("times_dim20.txt", "w");
    fprintf(times_file, "Name;N;Best;Worst;Mean;Median;Standard Deviation\n");

    test_Ackeleya(error_file, times_file, n_times, seed, dim);
    test_Rastrigin(error_file, times_file, n_times, seed, dim);
    test_Shubert(error_file, times_file, n_times, seed);
    test_Shekel(error_file, times_file, n_times, seed);
    test_Griewank(error_file, times_file, n_times, seed, dim);
    test_Perm(error_file, times_file, n_times, seed, dim);
    test_rotated(error_file, times_file, n_times, seed, dim);
    test_Zakharov(error_file, times_file, n_times, seed, dim);
    test_simple_quadratic(error_file, times_file, n_times, seed);
    test_sin_cos(error_file, times_file, n_times, seed);

    fclose(error_file);
    fclose(times_file);
}

int main()
{
    int seed = 1234;
    int n_times = 20;
    printf("Testing dim 5...\n");
    test_dim5(seed, n_times);
    printf("Testing dim 10...\n");
    test_dim10(seed, n_times);
    printf("Testing dim 20...\n");
    test_dim20(seed, n_times);
}

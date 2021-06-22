#include "des.h"
#include "test_funtions.h"

#include <stdlib.h>
#include <stdio.h>
#include <time.h>


int main() {
    int seed = 1234;
    int n_times = 20;
    int dim = 5;

    FILE *error_file = fopen("error.txt", "w");
    fprintf(error_file, "Name;N;Best;Worst;Mean;Median;Standard Deviation\n");
    FILE *times_file = fopen("times.txt", "w");
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

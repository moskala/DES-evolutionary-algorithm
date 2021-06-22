#include "des.h"
#include "test_funtions.h"

#include <stdlib.h>
#include <stdio.h>
#include <time.h>


int main() {
    int seed = 1234;
    int n_times = 20;
    int dim = 5;
    test_Ackeleya(n_times, seed, dim);
    test_Rastrigin(n_times, seed, dim);
    test_Shubert(n_times, seed);
    test_Shekel(n_times, seed);
    test_Griewank(n_times, seed, dim);
    test_Perm(n_times, seed, dim);
    test_rotated(n_times, seed, dim);
    test_Zakharov(n_times, seed, dim);
    test_simple_quadratic(n_times, seed);
    test_sin_cos(n_times, seed);
}

#include "des.h"
#include "test_funtions.h"

#include <stdlib.h>
#include <stdio.h>
#include <time.h>


int main() {
    // int seed = time(NULL);
    int seed = 1234;
    int n_times = 20;
    test_simple_quadratic(n_times, seed);
    test_sin_cos(n_times, seed);
    test_Ackeleya(n_times, seed, 5);
    test_Rastrigin(n_times, seed, 5);
    test_Shubert(n_times, seed);
    test_Shekel(n_times, seed);
    test_Griewank(n_times, seed, 5);
    test_Perm(n_times, seed, 5);
    test_rotated(n_times, seed, 5);
    test_Zakharov(n_times, seed, 5);
}

#include "des.h"
#include "test_funtions.h"

#include <stdlib.h>
#include <stdio.h>
#include <time.h>


int main() {
    int seed = time(NULL);
    int nTimes = 20;
    test_simple_quadratic(nTimes, seed);
    test_sin_cos(nTimes, seed);
    test_Ackeleya(nTimes, seed, 1);
    test_Rastrigin(nTimes, seed, 5);
    test_Schubert(nTimes, seed);
    test_Shekel(nTimes, seed);
}

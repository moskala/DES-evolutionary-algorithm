#include "test_funtions.h"
#include "des.h"
#include <stdio.h>
#include <time.h>


void print_result(int N, struct result res) {
    printf("Result - Fit: %f, Count: %d\n", res.best_fit, res.count);
    
    printf("Best solution:\n");
    for(int i = 0; i < N; ++i)
    {
        printf("%f, ", res.best_result[i]);
    }
    printf("\n");
}

void print_result_and_expected(int N, struct result res, double expectedValue) {
    printf("Result - Fit: %f, Expected: %f\n", res.best_fit, expectedValue);
}

void test_function(char *testName, int nTimes, int seed, double expectedValue, int N, double function_fn(int N, double[N]), double lower[N], double upper[N], bool logDesResults) {
    printf("Test: %s\n", testName);
    int current_seed = seed;
    for(int i = 0; i < nTimes; ++i){
        struct result res = des(N, NULL, function_fn, lower, upper, current_seed, logDesResults);
        printf("Test %d / %d\n", i+1, nTimes);
        print_result_and_expected(N, res, expectedValue);
        current_seed *= 2;
    }
}

int compare( const void* a, const void* b)
{
     double d_a = * ( (double*) a );
     double d_b = * ( (double*) b );

     if ( d_a == d_b ) return 0;
     else if ( d_a < d_b ) return -1;
     else return 1;
}



void test_function_with_statistics(char *test_name, int n_times, int seed, double expected_value, int N, double function_fn(int N, double[N]), double lower[N], double upper[N]){
    int current_seed = seed;
    
    double results[n_times];
    double times[n_times];

    for(int i = 0; i < n_times; ++i){
        clock_t pre = clock();
        struct result res = des(N, NULL, function_fn, lower, upper, current_seed, false);
        clock_t post = clock();
        double time = (post - pre) / ((double)CLOCKS_PER_SEC / 1000);
        results[i] = fabs(res.best_fit - expected_value);
        times[i] = time;
        current_seed *= 2;
        free(res.best_result);
    }
    qsort(results, n_times, sizeof(double), compare);
    double best = results[0];
    double worst = results[n_times-1];
    double median;
    double mean = 0.0;
    double std = 0.0;
    if(n_times % 2 == 0) 
    {
        int index_1 = n_times / 2 - 1;
        int index_2 = index_1 + 1;
        median = (results[index_1] + results[index_2]) / 2;
    }
    else 
    {
        median = results[(n_times + 1) / 2 - 1];
    }
    for(int i = 0; i < n_times; ++i)
    {
        mean += results[i];
    }
    mean /= (double)n_times;

    for(int i = 0; i < n_times; ++i)
    {
        std += pow(results[i] - mean, 2.0);
    }
    std /= n_times;
    std = sqrt(std);
    printf("Stats for %s\n", test_name);
    printf("best: %f, worst: %f, mean: %f, median: %f, std: %f\n", best, worst, mean, median, std);

    qsort(times, n_times, sizeof(double), compare);
    best = times[0];
    worst = times[n_times-1];
    mean = 0.0;
    std = 0.0;
    if(n_times % 2 == 0) 
    {
        int index_1 = n_times / 2 - 1;
        int index_2 = index_1 + 1;
        median = (times[index_1] + times[index_2]) / 2;
    }
    else 
    {
        median = times[(n_times + 1) / 2 - 1];
    }
    for(int i = 0; i < n_times; ++i)
    {
        mean += times[i];
    }
    mean /= (double)n_times;

    for(int i = 0; i < n_times; ++i)
    {
        std += pow(times[i] - mean, 2.0);
    }
    std /= n_times;
    std = sqrt(std);
    printf("Timings for %s\n", test_name);
    printf("best: %.1f ms, worst: %.1f ms, mean: %.1f ms, median: %.1f ms, std: %.1f ms\n", best, worst, mean, median, std);
}


void test_simple_quadratic(int nTimes, int seed) {
    char *name = "Test simple quandric";
    int N = 1;
    double init_point[N];
    double lower[1] = {-5};
    double upper[1] = {5};
    double expected_value = -9.0;
    // test_function(name, nTimes, seed, -9.0, 1, fun_simple_quadratic, lower, upper, false);
    test_function_with_statistics(name, nTimes, seed, expected_value, N, fun_simple_quadratic, lower, upper);

}

void test_sin_cos(int nTimes, int seed) {
    char *name = "Test sin*cos";
    int N = 2;
    double init_point[N];
    double lower[2] = {-5, -5};
    double upper[2] = {5, 5};
    double expected_value = -1;
    // test_function(name, nTimes, seed, expected_value, 2, fun_sin_cos, lower, upper, false);
    test_function_with_statistics(name, nTimes, seed, expected_value, N, fun_sin_cos, lower, upper);

}

// Minimum for x = 0, f(x) = 0
void test_Ackeleya(int nTimes, int seed, int dim){
    char *name = "Test Ackeley";
    int N = dim;
    double lower[N];
    double upper[N];
    for(int i = 0; i < N; ++i)
    {
        lower[i] = -30;
        upper[i] = 30;

    }
    double expected_value = 0;
    // test_function(name, nTimes, seed, expected_value, N, fun_Ackleya, lower, upper, false);
    test_function_with_statistics(name, nTimes, seed, expected_value, N, fun_Ackleya, lower, upper);

}

// Minimum for x = 0, f(x) = 0
void test_Rastrigin(int nTimes, int seed, int dim){
    char *name = "Test Rastrigin";
    int N = dim;
    double lower[N];
    double upper[N];
    for(int i = 0; i < N; ++i)
    {
        lower[i] = -5.12;
        upper[i] = 5.12;

    }
    double expected_value = 0;
    // test_function(name, nTimes, seed, expected_value, N, fun_Rastrigin, lower, upper, false);
    test_function_with_statistics(name, nTimes, seed, expected_value, N, fun_Rastrigin, lower, upper);

}

// Minimum f(x) = 186.73 for 760 different xs
void test_Schubert(int nTimes, int seed){
    char *name = "Test Schubert";
    int N = 2;
    double lower[N];
    double upper[N];
    for(int i = 0; i < N; ++i)
    {
        lower[i] = -5.12;
        upper[i] = 5.12;

    }
    double expected_value = -186.7309;
    // test_function(name, nTimes, seed, expected_value, N, fun_Schubert, lower, upper, false);
    test_function_with_statistics(name, nTimes, seed, expected_value, N, fun_Schubert, lower, upper);
}

// x = (4,4,4,4) f(x) = -10.5364
void test_Shekel(int nTimes, int seed){
    char *name = "Test Shekel";
    int N = 4;
    double lower[N];
    double upper[N];
    for(int i = 0; i < N; ++i)
    {
        lower[i] = 0;
        upper[i] = 10;

    }
    double expected_value = -10.5364;
    // test_function(name, nTimes, seed, expected_value, N, fun_Shekel, lower, upper, false);
    test_function_with_statistics(name, nTimes, seed, expected_value, N, fun_Shekel, lower, upper);
}

// Minimum for x* = 0, f(x) = 0
void test_Griewank(int nTimes, int seed, int dim){
    char *name = "Test Griewank";
    int N = dim;
    double lower[N];
    double upper[N];
    for(int i = 0; i < N; ++i)
    {
        lower[i] = -600;
        upper[i] = 600;

    }
    double expected_value = 0;
    // test_function(name, nTimes, seed, expected_value, N, fun_Griewank, lower, upper, false);
    test_function_with_statistics(name, nTimes, seed, expected_value, N, fun_Griewank, lower, upper);

}

// Perm Function 0, d, β
void test_Perm(int nTimes, int seed, int dim){
    char *name = "Test Perm";
    int N = dim;
    double lower[N];
    double upper[N];
    for(int i = 0; i < N; ++i)
    {
        lower[i] = -N;
        upper[i] = N;

    }
    double expected_value = 0;
    // test_function(name, nTimes, seed, expected_value, N, fun_Perm, lower, upper, false);
    test_function_with_statistics(name, nTimes, seed, expected_value, N, fun_Perm, lower, upper);
}


void test_rotated(int nTimes, int seed, int dim){
    char *name = "Test rotated";
    int N = dim;
    double lower[N];
    double upper[N];
    for(int i = 0; i < N; ++i)
    {
        lower[i] = -65.536;
        upper[i] =  65.536;

    }
    double expected_value = 0;
    // test_function(name, nTimes, seed, expected_value, N, fun_rotated, lower, upper, false);
    test_function_with_statistics(name, nTimes, seed, expected_value, N, fun_rotated, lower, upper);
}


void test_Zakharov(int nTimes, int seed, int dim){
    char *name = "Test Zakharov";
    int N = dim;
    double lower[N];
    double upper[N];
    for(int i = 0; i < N; ++i)
    {
        lower[i] = -5;
        upper[i] =  10;

    }
    double expected_value = 0;
    // test_function(name, nTimes, seed, expected_value, N, fun_Zakharov, lower, upper, false);
    test_function_with_statistics(name, nTimes, seed, expected_value, N, fun_Zakharov, lower, upper);
}

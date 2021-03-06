#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "des.h"

uint64_t xorshift64(uint64_t* seed)
{
    /* NB: This is the XorShift64* random number generator. The code is taken
     * from https://github.com/jj1bdx/xorshiftplus-c/blob/master/xorshift64star.c
     * with slight modifications. The original code is in public domain.
     *
     * We use this instead of rand(), because it makes code faster (about -40%
     * on my machine), is more reproductible across compilers and allows us to
     * stop relying on (s)rand. The generated number quality is irrelevant for
     * our needs. */
    *seed ^= *seed >> 12;
    *seed ^= *seed << 25;
    *seed ^= *seed >> 27;
    return *seed * UINT64_C(2685821657736338717);
}

double xorshift64_zero_one(uint64_t* seed)
{
    /* NB: This removes some precision from the value. This is because the
     * standard double-precision floating-point number can only hold 53 bits of
     * precision. Of course, this is not a problem at all. */
    return (double)xorshift64(seed) / (double)UINT64_MAX;
}

void print_array(int N, int M, double array[M][N])
{
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            printf("%f, ", array[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_array_one_dim(int N, double array[N])
{
    for (int j = 0; j < N; ++j) {
        printf("%f, ", array[j]);
    }
    printf("\n");
}

void deleteInfsNaNs(int N, double x[N])
{
    for (int i = 0; i < N; ++i) {
        if (isnan(x[i]) || isinf(x[i])) {
            x[i] = DBL_MAX;
        }
    }
}

void bounce_back_boundary(int N, double array[N], double lower[N], double upper[N])
{
    bool wasChanged = true;
    while (wasChanged) {
        wasChanged = false;
        for (int i = 0; i < N; ++i) {
            if (array[i] < lower[i]) {
                array[i] = lower[i] + fmod(fabs(lower[i] - array[i]), upper[i] - lower[i]);
                wasChanged = true;
            } else if (array[i] > upper[i]) {
                array[i] = upper[i] - fmod(fabs(array[i] - upper[i]), upper[i] - lower[i]);
                wasChanged = true;
            }
        }
        deleteInfsNaNs(N, array);
    }
}

// fn_ takes fn objective function :: [numeric] -> numeric passed as argument
double fn_(int N, double x[N], double lower[N], double upper[N], int* counteval, double (*fn)(int, double[N]))
{
    bool all_in_limits = true;
    for (int i = 0; i < N; i++) {
        if (x[i] < lower[i] || x[i] > upper[i]) {
            all_in_limits = false;
            break;
        }
    }
    if (all_in_limits) {
        *counteval += 1;
        return fn(N, x);
    } else {
        return DBL_MAX;
    }
}

// Fitness function wrapper for Lamarcian approach
void fitness_Lamarcian(int N, int m, double population[m][N], double lower[N], double upper[N], int* counteval,
    int budget, double (*fn)(int, double[N]), double results[m])
{
    if (m > 1) {
        if (*counteval + m <= budget) {
            for (int i = 0; i < m; ++i) {
                results[i] = fn_(N, population[i], lower, upper, counteval, fn);
            }
        } else {
            int budget_left = budget - *counteval;
            if (budget_left > 0) {
                for (int i = 0; i < budget_left; ++i) {
                    results[i] = fn_(N, population[i], lower, upper, counteval, fn);
                }
            }
            int start_i = budget_left > 0 ? budget_left : 0;
            for (int i = start_i; i < m; ++i) {
                results[i] = DBL_MAX;
            }
        }
    }
    else {
        if (*counteval < budget) {
            results[0] = fn_(N, population[0], lower, upper, counteval, fn);
        } else {
            results[0] = DBL_MAX;
        }
    }
}

// Fitness function wrapper for nonLamarcian approach
void fitness_non_Lamarcian(int N, int m, double population[m][N], double population_repaired[m][N], double fitness[m],
    double worst_fit, double population_fit[m])
{
    for (int i = 0; i < m; ++i) {
        deleteInfsNaNs(N, population[i]);
        deleteInfsNaNs(N, population_repaired[i]);
    }
    if (m > 1) {
        bool repaired_ind[m];
        double vec_dist[m];
        for (int i = 0; i < m; ++i) {
            double col_sum = 0.0;
            bool all_different = true;
            for (int j = 0; j < N; ++j) {
                if (population[i][j] == population_repaired[i][j]) {
                    all_different = false;
                    break;
                }
                col_sum += pow(population[i][j] - population_repaired[i][j], 2.0);
            }
            repaired_ind[i] = all_different;
            population_fit[i] = fitness[i];
            vec_dist[i] = col_sum;
        }
        for (int i = 0; i < m; ++i) {
            if (repaired_ind[i] == true) {
                population_fit[i] = worst_fit + vec_dist[i];
            }
        }
    }
    else {
        for (int i = 0; i < m; ++i) {
            population_fit[i] = fitness[i];
            double col_sum = 0.0;
            bool all_different = true;
            for (int j = 0; j < N; ++j) {
                if (population[i][j] == population_repaired[i][j]) {
                    all_different = false;
                    break;
                }
                col_sum += pow(population[i][j] - population_repaired[i][j], 2);
            }
            if (all_different) {
                population_fit[i] = worst_fit + col_sum;
            }
        }
    }
    deleteInfsNaNs(m, population_fit);
}

// fitness function from article
void fitness_function(int N, int M, double array_to_fit[M][N], double lower[N], double upper[N], double q_max,
    double population_fitness[M])
{
    for (int i = 0; i < M; ++i) {
        double sum = 0.0;
        for (int j = 0; j < N; ++j) {
            if (array_to_fit[i][j] > upper[j]) {
                sum += pow(array_to_fit[i][j] - upper[j], 2);
            } else if (array_to_fit[i][j] < lower[j]) {
                sum += pow(lower[j] - array_to_fit[i][j], 2);
            }
        }
        population_fitness[i] = sum + q_max;
    }
}

void get_max_value(int M, double array[M], double* max_value)
{
    *max_value = DBL_MIN;
    for (int i = 0; i < M; ++i) {
        if (*max_value < array[i]) {
            *max_value = array[i];
        }
    }
}

int get_best_fitness(int M, double fitness[M], double* best_fit)
{
    int index_best = -1;
    for (int i = 0; i < M; ++i) {
        if (*best_fit > fitness[i]) {
            *best_fit = fitness[i];
            index_best = i;
        }
    }
    return index_best;
}

void sort_population(int lambda, int N, double population[lambda][N], double eval[lambda])
{
    for (int end = lambda; end > 0; --end) {
        for (int i = 0; i + 1 < end; ++i) {
            double tmp_pop[N];
            double tmp_eval;

            if (eval[i] > eval[i + 1]) {
                memcpy(tmp_pop, population[i], sizeof(tmp_pop));
                memcpy(population[i], population[i + 1], sizeof(population[i]));
                memcpy(population[i + 1], tmp_pop, sizeof(population[i + 1]));

                tmp_eval = eval[i];
                eval[i] = eval[i + 1];
                eval[i + 1] = tmp_eval;
            }
        }
    }
}

bool stop_criterion(int N, int lambda, double population[lambda][N], double population_midpoint[N], int epsilon)
{
    double sum = 0.0;
    for (int j = 0; j < N; ++j) {
        double sum_sigma = 0.0;
        for (int i = 0; i < lambda; ++i) {
            sum_sigma += pow(population[i][j] - population_midpoint[j], 2);
        }
        sum_sigma = sum_sigma / (lambda - 1);
        sum += sqrt(sum_sigma);
    }
    bool result = (sum / N) < epsilon;
    return result;
}

double approx_normal(double mean, double variance_squared, uint64_t* seed)
{
    // https://stats.stackexchange.com/a/16411

    double ret = 0;
    for (int i = 0; i < 12; ++i) {
        ret += xorshift64_zero_one(seed) - 6;
    }

    return variance_squared * ret + mean;
}

double approx_delta(int N, uint64_t* seed)
{
    const int ITER = 1000;
    double ret = 0;
    for (int i = 0; i < ITER; ++i) {
        double ret1 = 0;
        for (int n = 0; n < N; ++n) {
            ret1 += pow(approx_normal(0, 1, seed), 2);
        }
        ret += sqrt(ret1) / ITER;
    }
    return 1.0 / ret;
}

// NB: We avoid using the N in [] due to a gcc internal compiler error.
struct result des(int N, double function_fn(int N, double[]), double lower[], double upper[], uint64_t seed,
    bool logRes)
{
    const int lambda = 4 * N;
    const int budget = 10000 * N;
    const int history_size = ceil(3 * sqrt(N)) + 6;
    const double scaling_factor = (double)1 / sqrt(2);
    const double gamma = approx_delta(N, &seed);
    const double epsilon = 10E-8;

    double c = (double)4 / (N + 4);
    const int mu = lambda / 2;
    int eval_count = 0;
    int restart_number = -1;
    double prev_fitness = NAN;

    double best_fit = HUGE_VAL;
    double* best_solution = calloc(N, sizeof(double));
    double prev_s[N];

    while (eval_count < budget) {
        restart_number += 1;

        double weights[mu];
        double weights_sum = 0;

        for (int i = 0; i < mu; ++i) {
            weights_sum += weights[i] = 1;
        }
        for (int i = 0; i < mu; ++i) {
            weights[i] /= weights_sum;
        }

        double weights_pop[lambda];
        double weights_pop_sum = 0;
        for (int i = 0; i < lambda; ++i) {
            weights_pop_sum += weights_pop[i] = 1;
        }
        for (int i = 0; i < lambda; ++i) {
            weights_pop[i] /= weights_pop_sum;
        }
        int hist_head = -1;
        int iter = 0;
        double history[history_size][mu][N];

        double population[lambda][N];
        for (int i = 0; i < lambda; ++i) {
            for (int j = 0; j < N; ++j) {
                population[i][j] = 0.8 * (xorshift64_zero_one(&seed) * (upper[j] - lower[j]) + lower[j]);
            }
        }
        double cum_mean[N];
        for (int i = 0; i < N; ++i) {
            cum_mean[i] = (lower[i] + upper[i]) / 2;
        }

        // Initial fitness
        double initial_fitness[lambda];
        fitness_Lamarcian(N, lambda, population, lower, upper, &eval_count, budget, function_fn, initial_fitness);

        // Initial worst fit
        double worst_fit = -HUGE_VAL;
        get_max_value(lambda, initial_fitness, &worst_fit);
        double prev_delta[N];
        for (int n = 0; n < N; ++n) {
            prev_delta[n] = 0;
        }

        bool stop = false;
        while (eval_count < budget && !stop) {
            iter += 1;
            hist_head += 1;
            hist_head %= history_size;
            double m[N];
            for (int n = 0; n < N; ++n) {
                m[n] = 0;
                for (int l = 0; l < lambda; ++l) {
                    m[n] += population[l][n] * weights_pop[n];
                }
            }

            double m_eval[1];
            fitness_Lamarcian(N, 1, &m, lower, upper, &eval_count, budget, function_fn, m_eval);
            fitness_non_Lamarcian(N, 1, &m, &m_eval, m_eval, worst_fit, m_eval);

            if (m_eval[0] < best_fit) {
                best_fit = m_eval[0];
                memcpy(best_solution, m, N * sizeof(double));
            }

            double population_repaired[lambda][N];
            memcpy(population_repaired, population, sizeof(population_repaired));
            for (int i = 0; i < lambda; ++i) {
                bounce_back_boundary(N, population_repaired[i], lower, upper);
            }

            double pop_eval[lambda];
            double fitness[lambda];
            fitness_Lamarcian(N, lambda, population, lower, upper, &eval_count, budget, function_fn, fitness);
            fitness_non_Lamarcian(N, lambda, population, population_repaired, fitness, worst_fit, pop_eval);

            // Find new minimum and maximum
            for (int i = 0; i < lambda; ++i) {
                if (pop_eval[i] < best_fit) {
                    best_fit = pop_eval[i];
                    memcpy(best_solution, population_repaired[i], N * sizeof(double));
                } else if (pop_eval[i] > worst_fit) {
                    worst_fit = pop_eval[i];
                }
            }
            sort_population(lambda, N, population, pop_eval);

            double s[N];
            for (int n = 0; n < N; ++n) {
                s[n] = 0;
                for (int m = 0; m < mu; ++m) {
                    s[n] += population[m][n] * weights[m];
                }
            }
            memcpy(prev_s, s, N * sizeof(double));

            // Check if the middle point is the best found so far
            for (int i = 0; i < N; ++i) {
                cum_mean[i] = 0.8 * cum_mean[i] + 0.2 * s[i];
            }
            bounce_back_boundary(N, cum_mean, lower, upper);
            double fitness_cum_mean[1];
            fitness_Lamarcian(N, 1, &cum_mean, lower, upper, &eval_count, budget, function_fn, fitness_cum_mean);
            if (fitness_cum_mean[0] < best_fit) {
                best_fit = fitness_cum_mean[0];
                memcpy(best_solution, cum_mean, N * sizeof(double));
            }

            double delta[N];
            for (int n = 0; n < N; ++n) {
                delta[n] = (1 - c) * prev_delta[n] + c * (s[n] - m[n]);
            }

            memcpy(history[hist_head], population, sizeof(history[hist_head]));
            memcpy(prev_delta, delta, sizeof(prev_delta));

            for (int l = 0; l < lambda; ++l) {
                int h_max = iter >= history_size ? history_size : hist_head;
                int h = h_max == 0 ? 0 : xorshift64(&seed) % h_max;
                int j = xorshift64(&seed) % mu;
                int k = xorshift64(&seed) % mu;

                for (int n = 0; n < N; ++n) {
                    double d = scaling_factor * (history[h][j][n] - history[h][k][n]) + delta[n] * gamma * approx_normal(0, 1, &seed);
                    population[l][n] = s[n] + d + epsilon * approx_normal(0, 1, &seed);
                }
            }

            for (int i = 0; i < lambda; ++i) {
                deleteInfsNaNs(N, population[i]);
            }

            if (prev_fitness != best_fit) {
                if (logRes) {
                    printf("Update result:, Fit: %f, Count: %d\n", best_fit, eval_count);
                    print_array_one_dim(N, best_solution);
                }
                prev_fitness = best_fit;
            }

            if (iter > 0) {
                stop = stop_criterion(N, lambda, population, prev_s, epsilon);
                if (stop && logRes) {
                    printf("Stop reached!\n");
                }
            }
        }
    }
    struct result res;
    res.best_fit = best_fit;
    res.best_result = best_solution;
    res.count = eval_count;

    return res;
}

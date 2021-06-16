#include "des.h"

#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

void bounce_back_boundary(int N, double array[N], double lower, double upper) {
    for (int i = 0; i < N; ++i) {
        if (array[i] < lower) {
            array[i] = lower + fmod(lower - array[i], upper - lower);
        } else if (array[i] > upper) {
            array[i] = upper - fmod(array[i] - upper, upper - lower);
        }
    }
    
    // TODO: deleteInfsNaNs ???
}

void sort_population(int lambda, int N, double population[lambda][N], double eval[lambda]) {
    // NB: Just the laziest bubble sort XD
    // Dr. Kaczmarski would be mad…
    // It doesn't matter anyway, since population is small.
    // TODO: Consider using a better sort (but qsort doesn't work well here -.-)

    for (int end = lambda; end > 0; --end) {
        for (int i = 0; i < end; ++i) {
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

// TODO: Maybe lower and upper should be arrays of size N
void des(int N, double initial_point[N], double function(float[N]), double lower, double upper) {
    const int lambda = 4 * N;
    const int budget = 10000 * N;
    const int history_size = ceil(3 * sqrt(N)) + 6;
    const double scaling_factor = 1; // Ft in code, F in paper
    const double epsilon = 0.000001;
    int eval_count = 0;
    int restart_number = -1;

    while (eval_count < budget) {
        restart_number += 1;

        int mu = lambda / 2;

        double weights[mu];
        double weights_sum = 0;
        for (int i = 0; i < mu; ++i) {
            weights_sum += weights[i] = log(mu + 1) - log(i + 1);
        }
        for (int i = 0; i < mu; ++i) {
            weights[i] /= weights_sum;
        }

        double weights_pop[lambda];
        double weights_pop_sum = 0;
        for (int i = 0; i < lambda; ++i) {
            weights_pop_sum += weights_pop[i] = log(lambda + 1) - log(i + 1);
        }
        for (int i = 0; i < lambda; ++i) {
            weights_pop[i] /= weights_pop_sum;
        }

        int hist_head = -1;
        int iter = 0;
        double history[history_size][mu][N];

        double population[lambda][N];
        for (int i = 0; i < lambda; ++i) {
            for (int j = 0; j < N; ++ j) {
                population[i][j] = 0.8 * ((double)rand() / RAND_MAX * (upper - lower) + lower); // TODO: ??? 0.8
            }
        }

        double cum_mean = (lower + upper) / 2;

        // TODO: these:
        /*
        selection <- rep(0, mu)
        selectedPoints <- matrix(0, nrow = N, ncol = mu)
        fitness <- fn_l(population)
        oldMean <- numeric(N)
        newMean <- par
        limit <- 0
        worst.fit <- max(fitness)
        */
        // and other vars

        float prev_delta[N];
        for (int n = 0; n < N; ++n) { prev_delta[n] = 0; }

        bool stop = false;
        while (eval_count < budget && !stop) {
            // TODO: stop

            iter += 1;
            hist_head += 1;
            hist_head %= history_size;

            float m[N];
            for (int n = 0; n < N; ++n) {
                m[n] = 0;
                for (int l = 0; l < lambda; ++l) {
                    // TODO: The paper actually has all the weights equal?!
                    m[n] += population[l][n] * weights_pop[n];
                }
            }

            eval_count += lambda + 1;
            double m_eval = 0; /* TODO */
            double pop_eval[lambda]; /* TODO */

            sort_population(lambda, N, population, pop_eval);
            
            float s[N];
            for (int n = 0; n < N; ++n) {
                s[n] = 0;
                for (int m = 0; m < mu; ++m) {
                    // TODO: The paper actually has all the weights equal?!
                    s[n] += population[m][n] * weights[n];
                }
            }

            float delta[N];
            float c = 0.777; // TODO: What even is c? (Maybe we should inline it, if we know what is it)
            for (int n = 0; n < N; ++n) {
                delta[n] = (1 - c) * prev_delta[n] + c * (s[n] - m[n]);
            }

            for (int l = 0; l < lambda; ++l) {
                int h_max = iter >= history_size ? history_size : hist_head;
                int h = rand() % h_max;
                int j = rand() % mu;
                int k = rand() % mu;

                for (int n = 0; n < N; ++n) {
                    float d = scaling_factor * (history[h][j][n] - history[h][k][n]) + /* TODO: the noise */ 0;
                    population[l][n] = s[n] + d + /* TODO: the noise */ 0;
                }
            }

            memcpy(history[hist_head], population, sizeof(history[hist_head]));
            memcpy(prev_delta, delta, sizeof(prev_delta));
        }
    }
}

#include "des.h"

#include <stdlib.h>
#include <stdbool.h>
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

// TODO: Maybe lower and upper should be arrays of size N
void des(int N, double initial_point[N], double function(float[N]), double lower, double upper) {
    const int lambda = 4 * N;
    const int budget = 10000 * N;
    const int history_size = ceil(3 * sqrt(N)) + 6;
    const double scaling_factor = 1; // Ft
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

        bool stop = false;
        while (eval_count < budget && !stop) {
            iter += 1;
            hist_head += 1;
            hist_head %= history_size;
        }
    }
}

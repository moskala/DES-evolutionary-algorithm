#include "des.h"

#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>

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

void deleteInfsNaNs(int N, double x[N]) {
    for(int i = 0 ; i < N; ++i){
        if(isnan(x[i])){
            x[i] = DBL_MAX;
        }
    }
}

// fitness function from article
void fitness_population(int N, int m, double population[m][N], double lower[N], double upper[N], double q_max, double population_fitness[m]) {
    for(int i = 0; i < m; ++i)
    {
        double sum = 0.0;
        for(int j = 0; j < N; ++j)
        {
            if(population[i][j] > upper[j])
            {
                sum += pow(population[i][j] - upper[j], 2);
            }
            else if (population[i][j] < lower[j])
            {
                sum += pow(lower[j] - population[i][j], 2);
            }
        }
        population_fitness[i] = sum + q_max;
    }
}

// fn_ takes fn objective function :: [numeric] -> numeric passed as argument
double fn_(int N, double x[N], double lower[N], double upper[N], int *counteval, double (*fn)(int, double*)){
    bool all_in_limits = true;
    for(int i = 0; i < N; i++){
        if(x[i] < lower[i] || x[i] > upper[i]){
            all_in_limits = false;
            break;
        }
    }
    if (all_in_limits) {
      *counteval += 1;
      return fn(N, x);
    }
    else {
      return DBL_MAX;
    }
}

// Fitness function wrapper for Lamarcian approach
double* fitness_Lamarcian(int N, int m, double population[m][N], double lower[N], double upper[N], int* counteval, int budget, double (*fn)(int, double*)) {
    double results[m];
    // More than one individual
    if(m > 1){ 
        if(*counteval + m <= budget){
            for(int i = 0; i < m; ++i){
                results[i] = fn_(N, population[i], lower, upper, counteval, fn);
            }
        }
        else{
            int budget_left = budget - *counteval;
            if (budget_left > 0) {
                for (int i = 0; i < budget_left; ++i) {
                    results[i] = fn_(N, population[i], lower, upper, counteval, fn);
                }
            }
            int start_i = budget_left > 0 ? budget_left : 0;
            for(int i = start_i; i < m; ++i) {
                results[i] = DBL_MAX;
            }
        }
    } 
    // One individual
    else {
        if(*counteval < budget) {
            results[0] = fn_(N, population[0], lower, upper, counteval, fn);
        }
        else {
            results[0] = DBL_MAX;
        }
    }
    return results;
}

// Fitness function wrapper for nonLamarcian approach 
// TODO - When population is no a matrix? Only in first approach?
double* fitness_non_Lamarcian(int N, int m, double population[m][N], double population_repaired[m][N], double fitness[m], double worst_fit) {
    
    double population_fit[m];
    for(int i = 0; i < m; ++i) {
        deleteInfsNaNs(N, population[i]);
        deleteInfsNaNs(N, population_repaired[i]);
        population_fit[i] = fitness[i];
    }
    // More than one individual
    if(m > 1) {
        
        bool repaired_ind[m];
        double vec_dist[m];
        for(int i = 0; i < m; ++i) {
            double col_sum = 0.0;
            bool all_different = true;
            for(int j = 0; j < N; ++j) {
                if(population[i][j] == population_repaired[i][j]){
                    all_different = false;
                    break;
                }
                col_sum += pow(population[i][j] - population_repaired[i][j], 2);
            }
            // repairedInd <- apply(P != P_repaired, 2, all)
            repaired_ind[i] = all_different;
            // P_fit <- fitness
            population_fit[i] = fitness[i];
            // vecDist <- colSums((P - P_repaired)^2)
            vec_dist[i] = col_sum;
        }
        // P_fit[which(repairedInd)] <- worst.fit + vecDist[which(repairedInd)]
        for(int i = 0; i < m; ++i) {
            if(repaired_ind[i] == true) {
                population_fit[i] = worst_fit + vec_dist[i];
            }
        }
    }
    // One individual
    else {
        // P_fit <- fitness - move at the beginning 
        for(int i = 0; i < m; ++i) {
            double col_sum = 0.0;
            bool all_different = true;
            for(int j = 0; j < N; ++j) {
                if(population[i][j] == population_repaired[i][j]){
                    all_different = false;
                    break;
                }
                col_sum += pow(population[i][j] - population_repaired[i][j], 2);
            }
            if(all_different) {
                population_fit[i] = worst_fit +  fitness[i];
            }
        }
    }
    // P_fit <- deleteInfsNaNs(P_fit)
    for(int i = 0; i < m; ++i) {
        population_fit[i] = isnan(population_fit[i]) ? DBL_MAX : population_fit[i];

    }
    return population_fit;
}




// TODO: Maybe lower and upper should be arrays of size N - they should be
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

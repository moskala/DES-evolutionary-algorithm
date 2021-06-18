#include "des.h"

#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
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
void fitness_Lamarcian(int N, int m, double population[m][N], double lower[N], double upper[N], int* counteval, int budget, double (*fn)(int, double*), double results[m]) {
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
}

// Fitness function wrapper for nonLamarcian approach 
// TODO - When population is no a matrix? Only in first approach?
void fitness_non_Lamarcian(int N, int m, double population[m][N], double population_repaired[m][N], double fitness[m], double worst_fit, double population_fit[m]) {
    
    for(int i = 0; i < m; ++i) {
        deleteInfsNaNs(N, population[i]);
        deleteInfsNaNs(N, population_repaired[i]);
        
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
            population_fit[i] = fitness[i];
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
                population_fit[i] = worst_fit +  col_sum;
            }
        }
    }
    // P_fit <- deleteInfsNaNs(P_fit)
    for(int i = 0; i < m; ++i) {
        population_fit[i] = isnan(population_fit[i]) ? DBL_MAX : population_fit[i];
    }
}

// fitness function from article
void fitness_function(int N, int M, double array_to_fit[M][N], double lower[N], double upper[N], double q_max, double population_fitness[M]) {
    for(int i = 0; i < M; ++i)
    {
        double sum = 0.0;
        for(int j = 0; j < N; ++j)
        {
            if(array_to_fit[i][j] > upper[j])
            {
                sum += pow(array_to_fit[i][j] - upper[j], 2);
            }
            else if (array_to_fit[i][j] < lower[j])
            {
                sum += pow(lower[j] - array_to_fit[i][j], 2);
            }
        }
        population_fitness[i] = sum + q_max;
    }
}

void get_max_value_and_index(int M, double array[M], double *max_value, int *index_of_max)
{
    *max_value = DBL_MIN;
    for(int i = 0; i < M; ++i)
    {
        if(*max_value < array[i])
        {
            *max_value = array[i];
            *index_of_max = i;
        }
    }
}



// TODO: Maybe lower and upper should be arrays of size N - they should be
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

double approx_normal(double mean, double variance_squared) {
    // https://stats.stackexchange.com/a/16411

    double ret = 0;
    for (int i = 0; i < 12; ++i) {
        ret += (double)rand() / RAND_MAX - 6;
    }

    return variance_squared * ret + mean;
}

// TODO: Maybe lower and upper should be arrays of size N
void des(int N, double initial_point[N], double function_fn(int N, double[N]), double lower, double upper) {
    const int lambda = 4 * N;
    const int budget = 10000 * N;
    const int history_size = ceil(3 * sqrt(N)) + 6;
    const double scaling_factor = 1; // Ft in code, F in paper
    const double epsilon = 0.000001;
    const double gamma = 0; // TODO: Actually it's some kind of norm
    int eval_count = 0;
    int restart_number = -1;
    double lower_bound[N];
    double upper_bound[N];
  
    double best_fit = DBL_MIN; // TODO: start value?
    double best_solution[N];
    
    for(int i = 0; i < N; ++i)
    {
        lower_bound[i] = lower;
        upper_bound[i] = upper;
    }

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

        // double cum_mean = (lower + upper) / 2;

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

        double initial_fitness[lambda];
        fitness_Lamarcian(N, lambda, population, lower_bound, upper_bound, &eval_count, budget, function_fn, initial_fitness);
        int max_fit_index = 0;
        double max_fit_value = DBL_MIN;
        get_max_value_and_index(lambda, initial_fitness, &max_fit_value, &max_fit_index);

        double worst_fit = max_fit_value; 

        float prev_delta[N];
        for (int n = 0; n < N; ++n) { prev_delta[n] = 0; }

        bool stop = false;
        while (eval_count < budget && !stop) {
            // TODO: stop

            iter += 1;
            hist_head += 1;
            hist_head %= history_size;
            eval_count += lambda + 1;

            double m[N];
            for (int n = 0; n < N; ++n) {
                m[n] = 0;
                for (int l = 0; l < lambda; ++l) {
                    // TODO: The paper actually has all the weights equal?!
                    m[n] += population[l][n] * weights_pop[n];
                }
            }

            double m_eval[lambda];
            fitness_function(N, 1, &m, lower_bound, upper_bound, worst_fit, m_eval);

            double pop_eval[lambda]; 
            fitness_function(N, 1, population, lower_bound, upper_bound, worst_fit, pop_eval);

            // Find new minimum
            for(int i = 0; i < lambda; ++i)
            {
                if(pop_eval[i] < best_fit)
                {
                    best_fit = pop_eval[i];
                    memcpy(best_solution, population[lambda], sizeof(best_solution));
                }
                else if(pop_eval[i] > best_fit)
                {
                    worst_fit = pop_eval[i];
                }
            }

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
                    float d = scaling_factor * (history[h][j][n] - history[h][k][n]) + delta[n] * gamma * approx_normal(0, 1);
                    population[l][n] = s[n] + d + epsilon * approx_normal(0, 1);
                }
            }

            memcpy(history[hist_head], population, sizeof(history[hist_head]));
            memcpy(prev_delta, delta, sizeof(prev_delta));
        }
    }
}

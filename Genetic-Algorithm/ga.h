#include "ga_parameters.h"
#include <cassert>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <time.h>
#include <vector>

#ifndef GA_H
#define GA_H

using namespace std;

extern int W[K], C[K], R[M];
extern double P[M], Q[M]; //destroy probability, expected damage
extern double sum_best_fiteness[GENERATIONS];
extern double sum_avg_fiteness[GENERATIONS];
extern double gene_best_fitness[ITERATION][GENERATIONS];
extern int iteration_num;
extern double random_gene_best_fitness[ITERATION];
extern int assignment[K];
void set_input_parameters();
void assignment_eval(int assignment[K]);

class Chromosome
{
public:
    Chromosome();
    void print_chromo();
    void mutate(); // one point flip
    void repair(); // greedy repair
    void evaluate();
    void evaluate_with_ML();
    void assert_feasible(int mode);
    int get_data(int row, int col) { return data[row][col]; }
    void put_data(int row, int col, int val) { data[row][col] = val; }
    double get_fitness() { return fitness; }
    void simple_simulation(const int data[M][K], double (&time)[S]);

private:
    int data[M][K];
    double column_sum[K], row_sum[M];
    double fitness;
    double weight_fitness;
    double region_time[K];
    void cumulative_sum_calc();
};

class Population
{
public:
    static const int POP_SIZE = 100;

    Population(int pop_size = POP_SIZE);
    void run();
    void statistics_info_calc();
    void print_population(int mode);
    void print_best_assignment();
    double get_best_fitness() { return best_fitness; }
    double get_avg_fitness() { return avg_fitness; }
    Chromosome &get_best_chromosome() { return individual[best_index]; }

private:
    int best_index, worst_index;
    Chromosome *individual, offspring;
    double best_fitness, avg_fitness, worst_fitness;
    //double gene_best_fitness[6];

    int select();                               // roulette wheel
    void crossover(int parent1, int parent2);   // block uniform
    void replacement(int parent1, int parent2); // similar survival
};
#endif // !GA_H

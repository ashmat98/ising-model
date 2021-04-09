//
// Created by Ashot on 3/27/2021.
//

#ifndef ISING_MODEL_BASE_H
#define ISING_MODEL_BASE_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <tuple>
#include <cmath>
#include <random>
#include <ctime>
#include <cassert>

//#include <pybind11/pybind11.h>
//#include <pybind11/numpy.h>
//namespace py = pybind11;

using namespace std;

class Base {
protected:
    mt19937 generator; //Standard mersenne_twister_engine seeded with rd()
    uniform_int_distribution<int> distribution_r;
    uniform_int_distribution<int> distribution_c;

public:


    char *_L{};
    char **L{};

    const int Nr, Nc; //shape of lattice

    inline int _g(int r, int c) const {
        return r * Nc + c;
    }

    char *MASK{};
    inline int g(int r, int c) const {
        return r * (Nc + 2) + c;
    }

    const char dir_r[4] = {-1, 0, 1, 0};
    const char dir_c[4] = {0, 1, 0, -1};

    int M; // magnetisation
    int E; // energy
    double T; //temperature
    int STEPS;
    unsigned int SEED;
    const int periodic_bc;

    Base(int Nr, int Nc, int periodic_bc = 1, int SEED_ = -1) :
            Nr(Nr), Nc(Nc), STEPS(0), M(0), E(0), T(0),
            periodic_bc(periodic_bc),
            generator(SEED + 1), distribution_r(1, Nr), distribution_c(1, Nc) {
        init_arrays();
//        init_mask();
        T = 0;
        if (SEED_ < 0)
            SEED = time(0);
        else
            SEED = SEED_;
    }

    void set_state(const vector<vector<char>> &state);

    vector<vector<char>> get_state() const;

    char get(int r, int c) const;

    char set(int r, int c, char val);

    void destruct_arrays() const {
        delete L;
        delete _L;
    }

    void init_arrays();
    virtual void init_mask();

    void random_init();

    void constant_init();

    double rand_std_uniform();

    tuple<int, int> rand_lattice_site();

    virtual void reset_history();

    int calc_E();

    int calc_M();

    int get_E() const;

    int get_M() const;

    void set_T(double temperature);

    double get_T() const;

    ~Base() {
        destruct_arrays();
    };

    tuple<int, int> next_lattice_site() const;
};

#endif //ISING_MODEL_BASE_H

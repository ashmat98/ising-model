//
// Created by Ashot on 3/27/2021.
//

#ifndef ISING_MODEL_BASE_H
#define ISING_MODEL_BASE_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <random>
#include <ctime>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;


using namespace std;

class Base {
public:
    const int Nr, Nc; //shape of lattice
    char *_L;

    inline const int _g(int r, int c) {
        return r * Nc + c;
    }

    char **L;

    inline const int g(int r, int c) {
        return r * (Nc + 2) + c;
    }

    const int dir_r[4] = {-1, 0, 1, 0};
    const int dir_c[4] = {0, 1, 0, -1};

    int M; // magnetisation
    int E; // energy
    double T; //temperature
    int STEPS;
    unsigned int SEED;

    Base(int Nr, int Nc, int SEED_ = -1) : Nr(Nr), Nc(Nc),
                                           STEPS(0) {
        init_arrays();
        T = 0;
        if (SEED_ < 0)
            SEED = time(0);
        else
            SEED = SEED_;
    }

    void set_state(py::array_t<int, py::array::c_style> state);

    py::array_t<char, py::array::c_style> get_state();

    int get(int r, int c);

    void destruct_arrays() {
        delete L, _L;
    }

    void init_arrays();

    void random_init();

    double rand_std_uniform();

    tuple<int, int> rand_lattice_site();

    virtual void reset_history();

    int calc_E();

    int calc_M();

    int get_E();

    int get_M();

    void set_T(double temperature);

    int get_T();

    ~Base() {
        destruct_arrays();
    };
};

#endif //ISING_MODEL_BASE_H

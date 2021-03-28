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
private:
    char *_L{};
    char **L{};

public:
    const int Nr, Nc; //shape of lattice

    inline int _g(int r, int c) const {
        return r * Nc + c;
    }


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
    const bool periodic_bc;

    Base(int Nr, int Nc, bool periodic_bc = true, int SEED_ = -1) :
            Nr(Nr), Nc(Nc), STEPS(0), M(0), E(0), T(0), periodic_bc(periodic_bc) {
        init_arrays();

        if (SEED_ < 0)
            SEED = time(0);
        else
            SEED = SEED_;
    }

    void set_state(const py::array_t<int, py::array::c_style> &state);

    py::array_t<char, py::array::c_style> get_state();

    char get(int r, int c) const;

    char set(int r, int c, char val);

    void destruct_arrays() const {
        delete L, _L;
    }

    void init_arrays();

    void random_init();

    double rand_std_uniform() const;

    tuple<int, int> rand_lattice_site() const;

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
};

#endif //ISING_MODEL_BASE_H

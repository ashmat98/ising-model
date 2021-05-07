//
// Created by Ashot on 3/27/2021.
//
/**
 * This is the base class of the model. It defines the topology of the
 * lattice, i.e. size and type of the boundary conditions. Supports operations
 * on the lattice and queries of various statistics of the system.
 */

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
    // Random engines
    mt19937 generator;
    uniform_int_distribution<int> distribution_r;
    uniform_int_distribution<int> distribution_c;

public:

    enum BoundaryCondition {
        Periodic,
        // periodic poundary with its obvious definition

        NotPoriodic,
        // free boundary conditions, i.e. spins at the sides
        // have less interactions

        Shifted,
        // This is exotic boundary conditions,
        // only for experimental purposes

        RowsPeriodic
        // This is mixture of periodic and free boundary conditions,
        // this is designed to model 1D ising model with periodic BC
        // thus make Nr = 1 and Nc the size of 1D chain
    };

    // pointers to store configuration of the system
    char *_L{};
    char **L{};

    const int Nr, Nc; //shape of lattice
    // Nr - number of rows
    // Nc - number of columns

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
    int E; // self interaction energy
    vector<int> history_M;
    vector<int> history_E;

    double T; //temperature
    int STEPS;
    unsigned int SEED;
    const BoundaryCondition bc;
    /*!
     * Constructor for the simulation engine.
     * @param Nr number of rows
     * @param Nc number of columns
     * @param bc (enum) boundary condition
     * @param SEED_ seed for randim engine
     */
    Base(int Nr, int Nc, BoundaryCondition bc = Periodic,
         int SEED_ = -1) :
            Nr(Nr), Nc(Nc), STEPS(0), M(0), E(0), T(0),
            bc(bc), distribution_r(1, Nr), distribution_c(1, Nc) {


        init_arrays();
//        init_mask();
        T = 0;
        if (SEED_ < 0)
            SEED = time(0);
        else
            SEED = SEED_;
        generator.seed(SEED);
    }

    void set_state(const vector<vector<char>> &state);

    vector<vector<char>> get_state() const;

    /*!
     * Get spin at the site (r, c).
     *
     * @param r row index (1..Nr) inclusive,
     * however 0 and Nr+1 is also permissible and refers to periodic
     * image in case of Periodic BC, or returns 0 if BC are free
     * @param c similar to r, range 0, (1, Nc) and Nc+1
     * @return spin state at site (r,c), +1, -1 or 0 if site is empty
     */
    char get(int r, int c) const;

    /*!
     * Set the spin state at site (r,c) to val
     * @param r row index
     * @param c column index
     * @param val new spin state
     * @return
     */
    char set(int r, int c, char val);

    void destruct_arrays() const {
        delete L;
        delete _L;
    }

    /*!
     * allocated memory for the lattice,
     * initializes *_L and **L according to the boundary conditions
     * the engine uses references to access sites to be faster
     */
    void init_arrays();

    virtual void init_mask();

    /*!
     * initializes spin states uniformly randomly +1 or -1
     */
    void random_init();

    /*!
     * initializes all spin states to +1
     */
    void constant_init();

    /*!
     * @return  uniformly distributed random number from range (0..1)
     */
    double rand_std_uniform();

    /*!
     * @return uniform random lattice site (r,c)
     */
    tuple<int, int> rand_lattice_site();

    /*!
     * Calculate lattice interaction energy for current state
     * @return current energy
     */
    int calc_E();

    /*!
     * Calculate lattice magnetisation for current state
     * @return current magnetisation
     */
    int calc_M();

    /*!
     * It is expected that simulation iteratively updates energy,
     * so calc_E() call is unnecessary and the energy can be directly got.
     * @return stored energy
     */
    int get_E() const;

    /*!
     * It is expected that simulation iteratively updates magnetisation,
     * so calc_E() call is unnecessary and the magnetisation can be
     * directly got.
     * @return stored magnetisation
     */
    int get_M() const;

    /*!
     * Engine stores historical values of the magnetisation
     * @return historical magnetisation
     */
    vector<int> get_sampled_M() const;

    /*!
     * Engine stores historical values of the energy
     * @return historical energy
     */
    vector<int> get_sampled_E() const;

    /*!
     * sores the statistics of the current state, the returned values
     * of the get_M() and get_E()
     */
    void store();

    /*!
     * clears up the containers for historical data
     */
    virtual void reset_history();

    void reset_sampled_M();

    void reset_sampled_E();

    /*!
     * set temperature of the system, thus the temperature can
     * be dynamically changed
     * @param temperature
     */
    void set_T(double temperature);


     /*!
      * get temperature of the system
      * @return current temperature
      */
    double get_T() const;

    ~Base() {
        destruct_arrays();
    };

    tuple<int, int> next_lattice_site() const;
};

#endif //ISING_MODEL_BASE_H

//
// Created by Ashot on 4/20/2021.
// engine specifically for Metropolis-Hasting simulation
//

#ifndef ISING_MODEL_METROPOLIS_HASTING_H
#define ISING_MODEL_METROPOLIS_HASTING_H


#include "base.h"

using namespace std;


class SimulateMH : public Base {
public:
    int FLIPS;
    int store_frequency;
    double H, omega;

    /*!
     * Constructor for the Metropolis-Hasting simulation engine.
     * @param Nr number of rows
     * @param Nc number of columns
     * @param frequency_to_store store system statistics every multiple
     * of "frequency_to_store" steps.
     * @param H External field
     * @param omega angular frequency of the external field
     * @param bc (enum) boundary condition
     * @param SEED_ seed for randim engine
     */
    SimulateMH(int Nr, int Nc, int frequency_to_store = 1,
               double H=0, double omega=0,
               BoundaryCondition bc = Periodic, int SEED = -1);

    /*!
     * returns field value corresponding to the current state
     * (field can be time dependent)
     * @return
     */
    double get_H() const{
        return H * cos(omega * STEPS);
    }

    /*!
     * returns Total energy both from self-interaction and the interaction
     * with the external field.
     * @return total energy
     */
    double get_total_E(){
        return get_E() - get_H()*get_M();
    }

//    double get_sampled_H(){
//        vector<int> sampled_H;
//
//        for (int i = 0; i < STEPS; i+= store_frequency) {
//            sampled_H.push_back(get_H(i))
//        }
//
//    }

    /*!
     * The loop of the engine
     * @param steps number of steps to perform
     * @param temperature set temperature of the system,
     * -1 means leave T unchanged
     */
    void make_steps(int steps, int temperature = -1);

    /*!
     * The operation of the simulation loop: single step of
     * Metropolis-Hasting algorithm
     */
    void single_step();

    /*!
     * Calculates the energy change when (r,c) site is flipped,
     * but leaving the state of the system unchanged
     * @param r row index
     * @param c column index
     * @return energy change
     */
    int flip_E_change(int r, int c) const;

    /*!
     * flips the site at (r,c). Changed energy and momentum variables
     * accordingly, s.t. no further calc_E() and calc_M() are required.
     * @param r
     * @param c
     * @return new spin state at (r,c)
     */
    inline char flip(int r, int c);
    /*!
     * clears up the containers for historical data and state variables.
     */
    void reset_history() override;
};




#endif //ISING_MODEL_METROPOLIS_HASTING_H

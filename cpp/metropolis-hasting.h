//
// Created by Ashot on 4/20/2021.
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

    SimulateMH(int Nr, int Nc, int frequency_to_store = 1,
               double H=0, double omega=0,
               BoundaryCondition bc = Periodic, int SEED = -1);

    double get_H() const{
        return H * cos(omega * STEPS);
    }

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

    void make_steps(int steps, int temperature = -1);

    void single_step();

    int flip_E_change(int r, int c);

    inline char flip(int r, int c);

    void reset_history() override;
};




#endif //ISING_MODEL_METROPOLIS_HASTING_H

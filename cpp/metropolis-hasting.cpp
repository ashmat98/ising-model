//
// Created by Ashot on 4/20/2021.
//

#include "metropolis-hasting.h"

SimulateMH::SimulateMH(int Nr, int Nc, int frequency_to_store,
                       double H, double omega,
                       Base::BoundaryCondition bc, int SEED) :
        Base(Nr, Nc, bc, SEED), H(H), omega(omega),
        FLIPS(0), store_frequency(frequency_to_store) {
}

void SimulateMH::make_steps(int steps, int temperature) {
    if (temperature >= 0) {
        set_T(temperature);
    }
    for (int i = 0; i < steps; ++i) {
        STEPS += 1;
        single_step();
        if (STEPS % store_frequency == 0) {
            store();
        }
    }
}

void SimulateMH::single_step() {
    static int r, c;
    int dE_interaction;
    static double dE_field, p;
    static tuple<int, int> site;

    site = rand_lattice_site();
    r = std::get<0>(site);
    c = std::get<1>(site);

    dE_interaction = flip_E_change(r, c);
    dE_field = 2 * get(r, c) * get_H();

    p = rand_std_uniform();

    if (T * log(p) < - (dE_interaction)) {
        FLIPS += 1;
        E += dE_interaction;
        M += 2 * flip(r, c);
    }
}

int SimulateMH::flip_E_change(int r, int c) {
    int dE;
    static char s;

    s =  get(r, c);
    dE = 0;

    for (int d = 0; d < 4; ++d) {
        dE += 2 * (+s) * get(r + dir_r[d], c + dir_c[d]);
    }
    return dE;
}

char SimulateMH::flip(int r, int c) {
    return set(r, c, -get(r, c));
}

void SimulateMH::reset_history() {
    Base::reset_history();
    FLIPS = 0;
}

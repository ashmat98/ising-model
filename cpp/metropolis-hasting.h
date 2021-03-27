#include "base.h"

namespace py = pybind11;


using namespace std;

// variables for moving to neighbours
const int dir_x[4] = {-1, 0, 1, 0};
const int dir_y[4] = {0, 1, 0, -1};


class Simulate_MH : public Base {
public:
    int FLIPS;

    Simulate_MH(int Nr, int Nc) : Base(Nr, Nc), FLIPS(0) {

    }

    void make_steps(int steps) {
        for (int i = 0; i < steps; i++) {
            STEPS += 1;
            single_step();
        }
    }

    void single_step() {
        auto[r, c] = rand_lattice_site();
        int dE = flip_E_change(r, c);
        double p = rand_std_uniform();
        if (T * log(p) < -dE) {
            FLIPS += 1;
            E += dE;
            M += 2 * flip(r, c);
        }
    }

    int flip_E_change(int r, int c) {
        int dE = 0;
        int s = get(r, c);
        for (int d = 0; d < 4; d++) {
            dE += (-s) * get(r + dir_r[d], c + dir_c[d]);
        }
        return dE;
    }

    inline int flip(int r, int c) {
        int s = *L[g(r, c)];
        return (*L[g(r, c)] *= -1);
    }
};

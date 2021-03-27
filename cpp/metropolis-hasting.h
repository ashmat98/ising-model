#include "base.h"

namespace py = pybind11;


using namespace std;

// variables for moving to neighbours
const int dir_x[4] = {-1, 0, 1, 0};
const int dir_y[4] = {0, 1, 0, -1};


class Simulate_MH : public Base {
public:
    int FLIPS;
    int store_frequency;
    vector<int> history_M;
    vector<int> history_E;

    Simulate_MH(int Nr, int Nc, int frequency_to_store = 1, int SEED = -1) :
            Base(Nr, Nc, SEED),
            FLIPS(0), store_frequency(frequency_to_store) {

    }

    void make_steps(int steps, int temperature = -1) {
        if (temperature >= 0) {
            set_T(temperature);
        }
        for (int i = 0; i < steps; i++) {
            STEPS += 1;
            single_step();
            if (STEPS % store_frequency == 0) {
                history_E.push_back(get_E());
                history_M.push_back(get_M());
            }
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
            dE += 2 * (-s) * get(r + dir_r[d], c + dir_c[d]);
        }
        return dE;
    }

    inline int flip(int r, int c) {
        int s = *L[g(r, c)];
        return (*L[g(r, c)] *= -1);
    }

    py::array_t<int> get_sampled_M() {
        return py::array_t<int>(history_M.size(), history_M.data());
    }

    py::array_t<int> get_sampled_E() {
        return py::array_t<int>(history_E.size(), history_E.data());
    }

    void reset_sampled_M() {
        history_M.clear();
    }

    void reset_sampled_E() {
        history_E.clear();
    }

    void reset_history() override {
        Base::reset_history();
        FLIPS = 0;
        reset_sampled_M();
        reset_sampled_E();
    }

};

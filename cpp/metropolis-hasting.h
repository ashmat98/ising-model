#include "base.h"

using namespace std;


class SimulateMH : public Base {
public:
    int FLIPS;
    int store_frequency;
    vector<int> history_M;
    vector<int> history_E;

    SimulateMH(int Nr, int Nc, int frequency_to_store = 1,
               int periodic_bc = 1, int SEED = -1) :
            Base(Nr, Nc, periodic_bc, SEED),
            FLIPS(0), store_frequency(frequency_to_store) {

    }

    void make_steps(int steps, int temperature = -1) {
        if (temperature >= 0) {
            set_T(temperature);
        }
        for (int i = 0; i < steps; ++i) {
            STEPS += 1;
            single_step();
            if (STEPS % store_frequency == 0) {
                history_E.push_back(get_E());
                history_M.push_back(get_M());
            }
        }
    }

    void single_step() {
        static int r, c, dE;
        static double p;
        static tuple<int, int> site;

        site = rand_lattice_site();
        r = std::get<0>(site);
        c = std::get<1>(site);

        dE = flip_E_change(r, c);
        p = rand_std_uniform();

        if (T * log(p) < -dE) {
            FLIPS += 1;
            E += dE;
            M += 2 * flip(r, c);
        }
    }

    int flip_E_change(int r, int c) {
        static int dE;
        static char s;
        dE = 0;
//        py::print("get 57", r, c, Nr, Nc);
        s = get(r, c);
        for (int d = 0; d < 4; ++d) {
            dE += 2 * (+s) * get(r + dir_r[d], c + dir_c[d]);
        }
        return dE;
    }

    inline char flip(int r, int c) {
        return set(r, c, -get(r, c));
    }

    vector<int> get_sampled_M() {
        return history_M;
    }

    vector<int> get_sampled_E() {
        return history_E;
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

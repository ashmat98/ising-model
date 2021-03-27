#include "base.h"


void Base::init_arrays() {
    _L = new char[Nr * Nc];
    L = new char *[(Nr + 2) * (Nc + 2)];

    for (int r = 1; r <= Nr; r++) {
        for (int c = 1; c <= Nc; c++) {
            L[g(r, c)] = &_L[_g(r - 1, c - 1)];
        }
    }
    for (int r = 1; r <= Nr; r++) {
        L[g(r, 0)] = &_L[_g(r - 1, Nc - 1)];
        L[g(r, Nc + 1)] = &_L[_g(r - 1, 0)];
    }
    for (int c = 1; c <= Nc; c++) {
        L[g(0, c)] = &_L[_g(Nr - 1, c - 1)];
        L[g(Nr + 1, c)] = &_L[_g(0, c - 1)];
    }
    L[g(0, 0)] = &_L[_g(Nr - 1, Nc - 1)];
    L[g(0, Nc + 1)] = &_L[_g(Nr - 1, 0)];
    L[g(Nr + 1, Nc + 1)] = &_L[_g(0, 0)];
    L[g(Nr + 1, 0)] = &_L[_g(0, Nc - 1)];
}

double Base::rand_std_uniform() {
    static mt19937 gen(SEED); //Standard mersenne_twister_engine seeded with rd()
    static uniform_real_distribution<double> dis(0, 1);
    return dis(gen);
}

tuple<int, int> Base::rand_lattice_site() {
    static mt19937 gen(SEED + 1); //Standard mersenne_twister_engine seeded with rd()
    static uniform_int_distribution<int> dis_r(1, Nr);
    static uniform_int_distribution<int> dis_c(1, Nc);

    return {dis_r(gen), dis_c(gen)};
}

void Base::set_state(py::array_t<int, py::array::c_style> state) {
    auto st = state.unchecked<2>();
    for (int r = 0; r < state.shape(0); r++) {
        for (int c = 0; c < state.shape(1); c++) {
            *L[g(r + 1, c + 1)] = int(st(r, c));
        }
    }
}

py::array_t<char, py::array::c_style> Base::get_state() {
//        return py::array_t<char>(Nc, _L[0]);
    return py::array_t<char>({Nr, Nc}, {Nc, 1}, _L);
}

int Base::get(int r, int c) {
    return int(*L[g(r, c)]);
}

void Base::random_init() {
    for (int r = 1; r <= Nr; r++) {
        for (int c = 1; c <= Nc; c++) {
            *L[g(r, c)] = 2 * int(rand_std_uniform() > 0.5) - 1;
        }
    }
    reset_history();
}

int Base::get_E() {
    return E;
}

int Base::get_M() {
    return M;
}

int Base::calc_E() {
    E = 0;
    for (int r = 1; r <= Nr; r++) {
        for (int c = 1; c <= Nc; c++) {
            for (int d = 0; d < 4; d++) {
                E += get(r, c) * get(r + dir_r[d], c + dir_c[d]);
            }
        }
    }
    assert(E % 2 == 0);
    E /= 2;
    return E;
}

int Base::calc_M() {
    M = 0;
    for (int r = 1; r <= Nr; r++) {
        for (int c = 1; c <= Nc; c++) {
            M += get(r, c);
        }
    }
    return M;
}

void Base::set_T(double temperature) {
    T = temperature;
}

int Base::get_T() {
    return T;
}

void Base::reset_history() {
    STEPS = 0;
    calc_E();
    calc_M();
}




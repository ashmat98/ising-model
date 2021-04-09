#include "base.h"

void Base::init_mask(){
    MASK = new char[(Nr + 2) * (Nc + 2)];
    for (int r = 0; r <= Nr + 1; r++) {
        for (int c = 0; c <= Nc + 1; c++) {
            MASK[g(r, c)] = 1;
        }
    }
    for (int r = 1; r <= Nr; r++) {
        MASK[g(r, 0)] = -1;
        MASK[g(r, Nc + 1)] = +1;
    }
    for (int c = 1; c <= Nc; c++) {
        MASK[g(0, c)] = -1;
        MASK[g(Nr + 1, c)] = +1;
    }
}
void Base::init_arrays() {
    _L = new char[Nr * Nc];
    L = new char *[(Nr + 2) * (Nc + 2)];
    char *ZERO;
    ZERO = new char;
    *ZERO = 0;

    for (int r = 1; r <= Nr; r++) {
        for (int c = 1; c <= Nc; c++) {
            L[g(r, c)] = &_L[_g(r - 1, c - 1)];
        }
    }
    if (periodic_bc == 1) {
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
    if (periodic_bc == 2) {
        int dr = Nr / 2, dc = Nc / 2;

        for (int r = 1; r <= Nr; r++) {
            L[g(r, 0)] = &_L[_g((r - 1 + dr) % Nr, Nc - 1)];
            L[g(r, Nc + 1)] = &_L[_g((Nr + r - 1 - dr) % Nr, 0)];
        }
        for (int c = 1; c <= Nc; c++) {
            L[g(0, c)] = &_L[_g(Nr - 1, (c - 1 + dc) % Nc)];
            L[g(Nr + 1, c)] = &_L[_g(0, (Nc + c - 1 - dc) % Nc)];
        }
        L[g(0, 0)] = ZERO;
        L[g(0, Nc + 1)] = ZERO;
        L[g(Nr + 1, Nc + 1)] = ZERO;
        L[g(Nr + 1, 0)] = ZERO;
    } else if (periodic_bc == 0) {

        for (int r = 0; r <= Nr + 1; r++) {
            L[g(r, 0)] = ZERO;
            L[g(r, Nc + 1)] = ZERO;
        }
        for (int c = 0; c <= Nc + 1; c++) {
            L[g(0, c)] = ZERO;
            L[g(Nr + 1, c)] = ZERO;
        }
    } else {
        assert(true);
    }
}

double Base::rand_std_uniform() {
    static uniform_real_distribution<double> dis(0, 1);
    return dis(generator);
}

tuple<int, int> Base::rand_lattice_site() {
    return make_tuple((int) distribution_r(generator),
                      (int) distribution_c(generator));
}

tuple<int, int> Base::next_lattice_site() const {
    static int i = -1;
    i = (i == Nr * Nc - 1) ? 0 : i + 1;

    return make_tuple((i / Nc) + 1, (i % Nc) + 1);
}

void Base::set_state(const vector<vector<char>> &state) {
    for (size_t r = 0; r < Nr; r++) {
        for (size_t c = 0; c < Nc; c++) {
            set(r + 1, c + 1, char(state[r][c]));
        }
    }
    calc_E();
    calc_M();
}

vector<vector<char>> Base::get_state() const {
    vector<vector<char>> state(Nr, std::vector<char>(Nc, 0));
    for (size_t r = 0; r < Nr; r++) {
        for (size_t c = 0; c < Nc; c++) {
            state[r][c] = get(r + 1, c + 1);
        }
    }
    return move(state);
}

char Base::get(int r, int c) const {
    return char(*L[g(r, c)]);
}

char Base::set(int r, int c, char val) {
//    assert(r>0 && r<=Nr && c>0 && c<=Nc);
    return (*L[g(r, c)]) = val;
}

void Base::random_init() {
    for (int r = 1; r <= Nr; r++) {
        for (int c = 1; c <= Nc; c++) {
            set(r, c, char(2 * int(rand_std_uniform() > 0.5) - 1));
        }
    }
    reset_history();
}

void Base::constant_init() {
    for (int r = 1; r <= Nr; r++) {
        for (int c = 1; c <= Nc; c++) {
            set(r, c, 1);
        }
    }
    reset_history();
}

int Base::get_E() const {
    return E;
}

int Base::get_M() const {
    return M;
}

int Base::calc_E() {
    E = 0;
    for (int r = 1; r <= Nr; r++) {
        for (int c = 1; c <= Nc; c++) {
            for (int d = 0; d < 4; d++) {
                E -= get(r, c) * get(r + dir_r[d], c + dir_c[d]);
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

double Base::get_T() const {
    return T;
}

void Base::reset_history() {
    STEPS = 0;
    calc_E();
    calc_M();
}





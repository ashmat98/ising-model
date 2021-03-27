#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <random>

using namespace std;


int main() {
    int s = 2;
    cout << (s *= -1) << endl;
    cout << log(0.999) << endl;
    static mt19937 gen(static_cast<long unsigned int>(34)); //Standard mersenne_twister_engine seeded with rd()
    static uniform_int_distribution<int> dis_r(1, 2);
    for (int i = 0; i < 100; i++) {
        cout << (dis_r(gen));
    }
    return 0;
}
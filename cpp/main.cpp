#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <random>
#include <chrono>

#include "metropolis-hasting.h"

using namespace std;


int main() {
    int N, M, stps = 100000;
    N = M = 32;
    SimulateMH engine(N, N, 1, 0,0,Base::Periodic, -1);
    engine.constant_init();
    engine.set_T(100);
    printf("%d\n", engine.get_E());
    engine.make_steps(2);
    printf("%d\n", engine.get_E());
    printf("%d\n", engine.calc_E());




//    scanf("%d %d %d", &N, &M, &stps);

//    printf("start\n");
//    for (int N = 89; N > 1; N--) {
//        auto start = std::chrono::system_clock::now();
//
//        SimulateMH engine(N, N, 1, 1);
//        engine.set_T(10);
//        engine.make_steps(stps);
//
//        auto end = std::chrono::system_clock::now();
//        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
//        printf("Time = %lld ms\n", static_cast<long long int>(elapsed.count()));
//    }
    return 0;
}
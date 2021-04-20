cmake_minimum_required(VERSION 3.17)
project(ising-model)

set(CMAKE_CXX_STANDARD 17)

add_executable(ising_model_cpp
        main.cpp
        base.cpp
        metropolis-hasting.cpp
        metropolis-hasting.h)



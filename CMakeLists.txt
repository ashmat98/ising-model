cmake_minimum_required(VERSION 3.17)
project(ising_model)

set(CMAKE_CXX_STANDARD 17)

include(fetch-pybind11.cmake)

add_subdirectory(cpp)
#
##set(CMAKE_CXX_FLAGS "-Wall -mavx2")
#SET(CMAKE_CXX_FLAGS_RELEASE "-ffast-math ${CMAKE_CXX_FLAGS_RELEASE}")

message("AAAAAaaa ${CMAKE_BUILD_TYPE} ;; ${CMAKE_CXX_FLAGS_RELEASE} ;; ${CMAKE_CXX_FLAGS_DEBUG} ;; ${CMAKE_CXX_FLAGS}")


pybind11_add_module(ising_model
        cpp/ising-model.cpp
        cpp/base.cpp
        cpp/base.h)

set_target_properties(ising_model PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/python)


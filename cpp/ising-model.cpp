
#include "metropolis-hasting.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;
#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;
using namespace pybind11::literals;

class BasePy : public Base {
public:
    using Base::Base;

    void set_state_pythonic(const py::array_t<int, py::array::c_style> &state) {
        auto st = state.unchecked<2>();
        for (int r = 0; r < state.shape(0); r++) {
            for (int c = 0; c < state.shape(1); c++) {
                set(r + 1, c + 1, char(st(r, c)));
            }
        }
        calc_E();
        calc_M();
    }

    py::array_t<char, py::array::c_style> get_state_pythonic() {
        return py::array_t<char>({Nr, Nc}, {Nc, 1}, _L);
    }

};

class SimulateMHPy : public SimulateMH {
public:
    using SimulateMH::SimulateMH;

    SimulateMHPy(int Nr, int Nc, int frequency_to_store = 1,
                 int periodic_bc = 1, int SEED = -1) :
            SimulateMH(Nr, Nc, frequency_to_store, periodic_bc, SEED) {

        py::print("create SimulateMH(", Nr, Nc, intptr_t(_L), intptr_t(L), ")");

    }

    void set_state_pythonic(const py::array_t<int, py::array::c_style> &state) {
        auto st = state.unchecked<2>();
        for (int r = 0; r < state.shape(0); r++) {
            for (int c = 0; c < state.shape(1); c++) {
                set(r + 1, c + 1, char(st(r, c)));
            }
        }
        calc_E();
        calc_M();
    }

    py::array_t<char, py::array::c_style> get_state_pythonic() {
        return py::array_t<char>({Nr, Nc}, {Nc, 1}, _L);
    }

    py::array_t<int> get_sampled_M_pythonic() {
        return py::array_t<int>(history_M.size(), history_M.data());
    }

    py::array_t<int> get_sampled_E_pythonic() {
        return py::array_t<int>(history_E.size(), history_E.data());
    }
};


PYBIND11_MODULE(ising_model, m) {
    m.doc() = R"pbdoc(
       ising_model simulator bindings for python
        -----------------------
    )pbdoc";
    py::class_<Base>(m, "__Base")
            .def(py::init<int, int, int, int>())
            .def("random_init", &Base::random_init)
            .def("constant_init", &Base::constant_init)
            .def("get", &Base::get)
            .def("get_E", &Base::get_E)
            .def("get_M", &Base::get_M)
            .def("calc_E", &Base::calc_E)
            .def("calc_M", &Base::calc_M)
            .def("get_T", &Base::get_T)
            .def("set_T", &Base::set_T)
            .def("rand_lattice_site", &Base::rand_lattice_site)
            .def_readonly("Nr", &Base::Nr)
            .def_readonly("Nc", &Base::Nc)
            .def_readonly("STEPS", &Base::STEPS)
            .def_readonly("SEED", &Base::SEED);

    py::class_<BasePy, Base>(m, "Base")
            .def(py::init<int, int, int, int>())
            .def("set_state", &BasePy::set_state_pythonic)
            .def("get_state", &BasePy::get_state_pythonic);

    py::class_<SimulateMH, Base>(m, "__SimulateMH")
            .def(py::init<int, int, int, int, int>())
            .def("make_steps", &SimulateMH::make_steps, "steps"_a, "temperature"_a = -1)
            .def("single_step", &SimulateMH::single_step)
            .def("reset_sampled_M", &SimulateMH::reset_sampled_M)
            .def("reset_sampled_E", &SimulateMH::reset_sampled_E)
            .def("flip_E_change", &SimulateMH::flip_E_change)
            .def_readonly("FLIPS", &SimulateMH::FLIPS);
//
    py::class_<SimulateMHPy, SimulateMH>(m, "SimulateMH")
            .def(py::init<int, int, int, int, int>())
            .def("get_sampled_M", &SimulateMHPy::get_sampled_M_pythonic)
            .def("get_sampled_E", &SimulateMHPy::get_sampled_E_pythonic)
            .def("set_state", &SimulateMHPy::set_state_pythonic)
            .def("get_state", &SimulateMHPy::get_state_pythonic);
    // ... remainder ...


#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
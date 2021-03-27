
#include "metropolis-hasting.h"

#include <pybind11/pybind11.h>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(ising_model, m) {
    m.doc() = R"pbdoc(
       ising_model simulator bindings for python
        -----------------------
    )pbdoc";

    py::class_<Base>(m, "Base")
            .def(py::init<int, int, int>())
            .def("random_init", &Base::random_init)
            .def("set_state", &Base::set_state)
            .def("get", &Base::get)
            .def("get_state", &Base::get_state)
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

    py::class_<Simulate_MH, Base>(m, "Simulate_MH")
            .def(py::init<int, int, int, int>())
            .def("make_steps", &Simulate_MH::make_steps, "steps"_a, "temperature"_a = -1)
            .def("single_step", &Simulate_MH::single_step)
            .def("get_sampled_M", &Simulate_MH::get_sampled_M)
            .def("get_sampled_E", &Simulate_MH::get_sampled_E)
            .def("reset_sampled_M", &Simulate_MH::reset_sampled_M)
            .def("reset_sampled_E", &Simulate_MH::reset_sampled_E)
            .def_readonly("FLIPS", &Simulate_MH::FLIPS);

    // ... remainder ...


#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
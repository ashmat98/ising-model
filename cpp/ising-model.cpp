
#include "metropolis-hasting.h"

#include <pybind11/pybind11.h>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

PYBIND11_MODULE(ising_model, m) {
    m.doc() = R"pbdoc(
       ising_model simulator bindings for python
        -----------------------
    )pbdoc";

    py::class_<Simulate_MH>(m, "Simulate_MH")
            .def(py::init<int, int>())
            .def("random_init", &Simulate_MH::random_init)
            .def("set_state", &Simulate_MH::set_state)
            .def("get", &Simulate_MH::get)
            .def("get_state", &Simulate_MH::get_state)
            .def_readonly("Nr", &Simulate_MH::Nr)
            .def_readonly("Nc", &Simulate_MH::Nc);

    // ... remainder ...


#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
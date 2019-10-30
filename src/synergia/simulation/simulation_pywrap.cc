
#include <pybind11/pybind11.h>

#include "synergia/simulation/propagator.h"

namespace py = pybind11;

PYBIND11_MODULE(simulation, m)
{
    py::class_<Propagator>(m, "Propagator")
        .def(py::init<Lattice const&, Stepper const&>())
        .def("propagate", &Propagator::propagate)
        .def("print_steps", &Propagator::print_steps)
        ;
}




#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "synergia/simulation/propagator.h"
#include "synergia/simulation/lattice_simulator.h"

namespace py = pybind11;
using namespace py::literals;

PYBIND11_MODULE(simulation, m)
{
    py::class_<Propagator>(m, "Propagator")
        .def(py::init<Lattice const&, Stepper const&>())
        .def("propagate", &Propagator::propagate)
        .def("print_steps", &Propagator::print_steps)
        ;

    m.def_submodule("Lattice_simulator")
        .def( "tune_linear_lattice", 
                &Lattice_simulator::tune_linear_lattice,
                "Tune linear lattice.",
                "lattice"_a )

        .def( "tune_circular_lattice", 
                &Lattice_simulator::tune_circular_lattice,
                "Tune a circular lattice.",
                "lattice"_a, "tolerance"_a = 1.0e-13 )
        ;

}



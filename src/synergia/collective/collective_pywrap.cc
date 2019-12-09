
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


#include "space_charge_2d_open_hockney.h"

namespace py = pybind11;
using namespace py::literals;


PYBIND11_MODULE(collective, m)
{
    //Space_charge_2d_open_hockney_options
    py::class_<Space_charge_2d_open_hockney_options, CO_options>
        (m, "Space_charge_2d_open_hockney_options")
        .def( py::init<int, int, int>(),
                "Construct the space charge 2d open hockney solver.",
                "gridx"_a, "gridy"_a, "gridz"_a )
        .def_readwrite( "comm_group_size", 
                &Space_charge_2d_open_hockney_options::comm_group_size,
                "Communication group size (must be 1 when running on GPUs)." )
        ;
}


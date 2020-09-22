
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


#include "dummy_collective_operator.h"
#include "space_charge_2d_open_hockney.h"
#include "space_charge_3d_open_hockney.h"
#include "space_charge_rectangular.h"

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
                "Communication group size (must be 1 on GPUs)." )
        ;

    py::class_<Space_charge_3d_open_hockney_options, CO_options>
        (m, "Space_charge_3d_open_hockney_options")
        .def( py::init<int, int, int>(),
                "Construct the space charge 3d open hockney solver.",
                "gridx"_a, "gridy"_a, "gridz"_a )
        .def_readwrite( "comm_group_size", 
                &Space_charge_3d_open_hockney_options::comm_group_size,
                "Communication group size (must be 1 on GPUs)." )
        ;

 
    py::class_<Space_charge_rectangular_options, CO_options>
        (m, "Space_charge_rectangular_options")
        .def( py::init<std::array<int, 3> const&, 
                std::array<double, 3> const& >(),
                "Construct the rectangular space charge solver.",
                "grid_shape"_a, "pipe_size"_a )
        .def_readwrite( "comm_group_size", 
                &Space_charge_rectangular_options::comm_group_size,
                "Communication group size (must be 1 on GPUs)." )
        ;

    py::class_<Dummy_CO_options, CO_options>
        (m, "Dummy_CO_options")
        .def( py::init<>(), "Construct a dummy collective operator.")
        ;
}


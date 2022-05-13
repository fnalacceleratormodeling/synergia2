#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "synergia/simulation/collective_operator_options.h"
#include "synergia/simulation/implemented_collective_options.h"

namespace py = pybind11;
using namespace py::literals;

PYBIND11_MODULE(collective, m)
{
  // Space_charge_2d_open_hockney_options
  py::class_<Space_charge_2d_open_hockney_options>(
    m, "Space_charge_2d_open_hockney_options")
    .def(py::init<int, int, int>(),
         "Construct the space charge 2d open hockney solver.",
         "gridx"_a,
         "gridy"_a,
         "gridz"_a)
    .def_readwrite("comm_group_size",
                   &Space_charge_2d_open_hockney_options::comm_group_size,
                   "Communication group size (must be 1 on GPUs).");

  py::class_<Space_charge_3d_open_hockney_options>(
    m, "Space_charge_3d_open_hockney_options")
    .def(py::init<int, int, int>(),
         "Construct the space charge 3d open hockney solver.",
         "gridx"_a,
         "gridy"_a,
         "gridz"_a)
    .def_readwrite("comm_group_size",
                   &Space_charge_3d_open_hockney_options::comm_group_size,
                   "Communication group size (must be 1 on GPUs).");

  py::class_<Space_charge_rectangular_options>(
    m, "Space_charge_rectangular_options")
    .def(py::init<std::array<int, 3> const&, std::array<double, 3> const&>(),
         "Construct the rectangular space charge solver.",
         "grid_shape"_a,
         "pipe_size"_a)
    .def_readwrite("comm_group_size",
                   &Space_charge_rectangular_options::comm_group_size,
                   "Communication group size (must be 1 on GPUs).");

  py::class_<Impedance_options>(m, "Impedance_options")
    .def(py::init<std::string const&, std::string const&, int>(),
         "Construct the impedance operator.",
         "wake_file"_a,
         "wake_type"_a,
         "z_grid"_a)

    .def_readwrite("z_grid",
                   &Impedance_options::z_grid,
                   "Size of z-grid (int, default to 1000).")

    .def_readwrite("full_machine",
                   &Impedance_options::full_machine,
                   "Full machine (boolean, default to false).")

    .def_readwrite("nstored_turns",
                   &Impedance_options::nstored_turns,
                   "Number of stored turns (int, default to 15).")

    .def_readwrite("num_buckets",
                   &Impedance_options::num_buckets,
                   "Number of buckets (int, default to 1).")

    .def_readwrite("orbit_length",
                   &Impedance_options::orbit_length,
                   "Orbit length (double, default to 1.0).")

    .def_readwrite("bunch_spacing",
                   &Impedance_options::bunch_spacing,
                   "Bunch spacing (double, default to 1.0).");

  py::class_<Dummy_CO_options>(m, "Dummy_CO_options")
    .def(py::init<>(), "Construct a dummy collective operator.");
}

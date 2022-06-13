
#include <Kokkos_Core.hpp>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace py::literals;

PYBIND11_MODULE(utils, m)
{
  m.def("init", []() {
    // MPI_Init(NULL, NULL);
    Kokkos::initialize();
  });

  m.def("finalize", []() {
    // MPI_Finalize();
    Kokkos::finalize();
  });
}

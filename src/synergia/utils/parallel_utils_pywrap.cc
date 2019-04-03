
#include <pybind11/pybind11.h>

#include "commxx.h"
#include "parallel_utils.h"

namespace py = pybind11;

PYBIND11_MODULE(parallel_utils_py, m)
{
    m.def("generate_subcomms", generate_subcomms);
    m.def("make_optimal_spc_comm", make_optimal_spc_comm);
    m.def("decompose_1d_local", decompose_1d_local);
}



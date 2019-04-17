
#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_MODULE(bunch_py, m)
{
#if 0
    m.def("generate_subcomms", generate_subcomms);
    m.def("make_optimal_spc_comm", make_optimal_spc_comm);
    m.def("decompose_1d_local", decompose_1d_local);

    py::class_<Logger>(m, "Logger")
        .def(py::init<int>())
        .def(py::init<int, std::string const& >())
        .def(py::init<int, std::string const&, bool >())
        .def(py::init<int, std::string const&, bool, bool >())
        .def(py::init<std::string const& >())
        .def(py::init<std::string const&, bool >())
        .def("write", &Logger::write, py::return_value_policy::reference_internal)
        .def("flush", &Logger::flush, py::return_value_policy::reference_internal)
        ;
#endif
}



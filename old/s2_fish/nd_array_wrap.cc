#include "nd_array.h"
#undef _POSIX_C_SOURCE
#include <boost/python.hpp>

#include "container_conversions.h"

using namespace boost::python;

BOOST_PYTHON_MODULE(nd_array)
{
    void (Real_nd_array::*real_reshape)
    (std::vector<int> const& indices) = &Real_nd_array::reshape;

    void (Real_nd_array::*real_set_nocheck)
    (std::vector<int> const& indices, double val) =
        &Real_nd_array::set_nocheck;

    void (Real_nd_array::*real_set)
    (std::vector<int> const& indices, double val) =
        &Real_nd_array::set;

    void (Real_nd_array::*real_add_to_point)
    (std::vector<int> const& indices, double val) =
        &Real_nd_array::add_to_point;

    double (Real_nd_array::*real_get)
    (std::vector<int> const& indices) const = &Real_nd_array::get;

    class_<Real_nd_array>("Real_nd_array", init<>())
    .def(init<std::vector<int> >())

    .def("copy", &Real_nd_array::copy<double>)

    .def("reshape", real_reshape)
    .def("get_shape", &Real_nd_array::get_shape)

    .def("zero_all", &Real_nd_array::zero_all)
    .def("set_nocheck", real_set_nocheck)
    .def("set", real_set)
    .def("add_to_point", real_add_to_point)
    .def("get", real_get)

    .def("scale", &Real_nd_array::scale)
    .def("add", &Real_nd_array::add)

    .def("get_length", &Real_nd_array::get_length)

    .def("describe", &Real_nd_array::describe)
    .def("print_", &Real_nd_array::print)

    .def("write_to_file", &Real_nd_array::write_to_file)
    .def("read_from_file", &Real_nd_array::read_from_file)
    ;

    void (Complex_nd_array::*complex_reshape)
    (std::vector<int> const& indices) = &Complex_nd_array::reshape;

    void (Complex_nd_array::*complex_set_nocheck)
    (std::vector<int> const& indices, std::complex<double> val) =
        &Complex_nd_array::set_nocheck;

    void (Complex_nd_array::*complex_set)
    (std::vector<int> const& indices, std::complex<double> val) =
        &Complex_nd_array::set;

    void (Complex_nd_array::*complex_add_to_point)
    (std::vector<int> const& indices, std::complex<double> val) =
        &Complex_nd_array::add_to_point;

    std::complex<double> (Complex_nd_array::*complex_get)
    (std::vector<int> const& indices) const = &Complex_nd_array::get;

    class_<Complex_nd_array>("Complex_nd_array", init<>())
    .def(init<std::vector<int> >())

    .def("copy", &Complex_nd_array::copy<std::complex<double> >)
    .def("copy_real", &Complex_nd_array::copy<double>)
    //    .def("real",&Complex_nd_array::real)

    .def("reshape", complex_reshape)
    .def("get_shape", &Complex_nd_array::get_shape)

    .def("zero_all", &Complex_nd_array::zero_all)
    .def("set_nocheck", complex_set_nocheck)
    .def("set", complex_set)
    .def("add_to_point", complex_add_to_point)
    .def("get", complex_get)

    .def("scale", &Complex_nd_array::scale)
    .def("add", &Complex_nd_array::add)

    .def("get_length", &Complex_nd_array::get_length)

    .def("describe", &Complex_nd_array::describe)
    .def("print_", &Complex_nd_array::print)

    .def("write_to_file", &Complex_nd_array::write_to_file)
    .def("read_from_file", &Complex_nd_array::read_from_file)
    ;

    scitbx::boost_python::container_conversions::from_python_sequence <
    std::vector<int>,
    scitbx::boost_python::container_conversions::variable_capacity_policy > ();

    boost::python::to_python_converter <
    std::vector<int>,
    scitbx::boost_python::container_conversions::to_tuple <
    std::vector<int> > > ();
}


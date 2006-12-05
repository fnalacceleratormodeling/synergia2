#include "scalar_field.h"

#undef _POSIX_C_SOURCE
#include <boost/python.hpp>

using namespace boost::python;

// struct int3_to_tuple
// {
//   static PyObject* convert(int3 const& x)
//   {
//     return incref(make_tuple(x[0],x[1],x[2]).ptr());
//   }
// };

// struct double3_to_tuple
// {
//   static PyObject* convert(double3 const& x)
//   {
//     return incref(make_tuple(x[0],x[1],x[2]).ptr());
//   }
// };


BOOST_PYTHON_MODULE(s2_fish)
{
//   to_python_converter<int3, int3_to_tuple>();
//   to_python_converter<double3, double3_to_tuple>();
  class_<int3>("int3",init<int, int, int>())
    .def("get",&int3::get)
    .def("set",&int3::set)
    ;

  class_<double3>("double3",init<double, double, double>())
    .def("get",&double3::get)
    .def("set",&double3::set)
    ;

  class_<Scalar_Field>("Scalar_Field",init<>())
    .def(init<int3, double3, double3>())
    .def("set_num_points",&Scalar_Field::set_num_points)
    .def("get_num_points",&Scalar_Field::get_num_points)
    .def("set_physical_params",&Scalar_Field::set_physical_params)
    .def("get_physical_size",&Scalar_Field::get_physical_size)
    .def("get_physical_offset",&Scalar_Field::get_physical_offset)
    .def("zero_the_points",&Scalar_Field::zero_the_points)
    .def("set_point",&Scalar_Field::set_point)
    .def("get_point",&Scalar_Field::get_point)
    .def("add_to_point",&Scalar_Field::add_to_point)
    .def("get_nearest_indices",&Scalar_Field::get_nearest_indices)
    ;
}


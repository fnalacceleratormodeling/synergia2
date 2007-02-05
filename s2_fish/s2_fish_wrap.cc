#include "scalar_field.h"
#include "deposit.h"
#include "solver.h"
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

  class_<Real_scalar_field>("Real_scalar_field",init<>())
    .def(init<int3, double3, double3>())
    .def("set_num_points",&Real_scalar_field::set_num_points)
    .def("get_num_points",&Real_scalar_field::get_num_points)
    .def("set_physical_params",&Real_scalar_field::set_physical_params)
    .def("get_physical_size",&Real_scalar_field::get_physical_size)
    .def("get_physical_offset",&Real_scalar_field::get_physical_offset)
    .def("zero_the_points",&Real_scalar_field::zero_the_points)
    .def("set_point",&Real_scalar_field::set_point)
    .def("get_point",&Real_scalar_field::get_point)
    .def("add_to_point",&Real_scalar_field::add_to_point)
    .def("get_leftmost_indices",&Real_scalar_field::get_leftmost_indices)
    .def("get_leftmost_offsets",&Real_scalar_field::get_leftmost_offsets)
    .def("print_points",&Real_scalar_field::print_points)
    ;

  class_<Complex_scalar_field>("Complex_scalar_field",init<>())
    .def(init<int3, double3, double3>())
    .def("set_num_points",&Complex_scalar_field::set_num_points)
    .def("get_num_points",&Complex_scalar_field::get_num_points)
    .def("set_physical_params",&Complex_scalar_field::set_physical_params)
    .def("get_physical_size",&Complex_scalar_field::get_physical_size)
    .def("get_physical_offset",&Complex_scalar_field::get_physical_offset)
    .def("zero_the_points",&Complex_scalar_field::zero_the_points)
    .def("set_point",&Complex_scalar_field::set_point)
    .def("get_point",&Complex_scalar_field::get_point)
    .def("add_to_point",&Complex_scalar_field::add_to_point)
    .def("get_leftmost_indices",&Complex_scalar_field::get_leftmost_indices)
    .def("get_leftmost_offsets",&Complex_scalar_field::get_leftmost_offsets)
    .def("print_points",&Complex_scalar_field::print_points)
    ;

  def("deposit_charge_cic",deposit_charge_cic);
  def("deposit_charge_ngp",deposit_charge_ngp);
  
  def("solver",solver);
  def("fft_tester",fft_tester);
}


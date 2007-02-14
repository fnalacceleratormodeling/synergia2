#include "nd_array.h"
#include "scalar_field.h"
#include "deposit.h"
#include "solver.h"
#undef _POSIX_C_SOURCE
#include <boost/python.hpp>

#include "container_conversions.h"

using namespace boost::python;

BOOST_PYTHON_MODULE(s2_fish)
{
  //---------------------------------------------------------------------
  // std::vector<> conversions
  //---------------------------------------------------------------------
  scitbx::boost_python::container_conversions::from_python_sequence<
    std::vector<int>,
    scitbx::boost_python::container_conversions::variable_capacity_policy>();

  boost::python::to_python_converter<
    std::vector<int>,
    scitbx::boost_python::container_conversions::to_tuple<
    std::vector<int> > >();

  scitbx::boost_python::container_conversions::from_python_sequence<
    std::vector<double>,
    scitbx::boost_python::container_conversions::variable_capacity_policy>();

  boost::python::to_python_converter<
    std::vector<double>,
    scitbx::boost_python::container_conversions::to_tuple<
    std::vector<double> > >();


  //---------------------------------------------------------------------
  // Real_nd_array
  //---------------------------------------------------------------------
  void (Real_nd_array::*real_reshape)
    (std::vector<int> const& indices) = &Real_nd_array::reshape;

  void (Real_nd_array::*real_set)
    (std::vector<int> const& indices, double val) =
    &Real_nd_array::set;

  void (Real_nd_array::*real_add_to_point)
    (std::vector<int> const& indices, double val) =
    &Real_nd_array::add_to_point;

  double (Real_nd_array::*real_get)
    (std::vector<int> const& indices) const= &Real_nd_array::get;

  class_<Real_nd_array>("Real_nd_array",init<>())
    .def(init<std::vector<int> >())
    
    .def("copy",&Real_nd_array::copy<double>)
    
    .def("reshape",real_reshape)
    .def("freeze_shape",&Real_nd_array::freeze_shape)
    .def("get_shape",&Real_nd_array::get_shape)

    .def("zero_all",&Real_nd_array::zero_all)
    .def("set",real_set)
    .def("add_to_point",real_add_to_point)
    .def("get",real_get)

    .def("scale",&Real_nd_array::scale)
    .def("add",&Real_nd_array::add)
    
    .def("get_length",&Real_nd_array::get_length)

    .def("describe",&Real_nd_array::describe)
    .def("print_",&Real_nd_array::print)

    .def("write_to_file",&Real_nd_array::write_to_file)
    .def("read_from_file",&Real_nd_array::read_from_file)
;

  //---------------------------------------------------------------------
  // Complex_nd_array
  //---------------------------------------------------------------------
  void (Complex_nd_array::*complex_reshape)
    (std::vector<int> const& indices) = &Complex_nd_array::reshape;

  void (Complex_nd_array::*complex_set)
    (std::vector<int> const& indices, std::complex<double> val) =
    &Complex_nd_array::set;

  void (Complex_nd_array::*complex_add_to_point)
    (std::vector<int> const& indices, std::complex<double> val) =
    &Complex_nd_array::add_to_point;

  std::complex<double> (Complex_nd_array::*complex_get)
    (std::vector<int> const& indices) const= &Complex_nd_array::get;

  class_<Complex_nd_array>("Complex_nd_array",init<>())
    .def(init<std::vector<int> >())
    
    .def("copy",&Complex_nd_array::copy<std::complex<double> >)
    .def("copy_real",&Complex_nd_array::copy<double>)
    //    .def("real",&Complex_nd_array::real)
    
    .def("reshape",complex_reshape)
    .def("freeze_shape",&Complex_nd_array::freeze_shape)
    .def("get_shape",&Complex_nd_array::get_shape)

    .def("zero_all",&Complex_nd_array::zero_all)
    .def("set",complex_set)
    .def("add_to_point",complex_add_to_point)
    .def("get",complex_get)

    .def("scale",&Complex_nd_array::scale)
    .def("add",&Complex_nd_array::add)
    
    .def("get_length",&Complex_nd_array::get_length)

    .def("describe",&Complex_nd_array::describe)
    .def("print_",&Complex_nd_array::print)

    .def("write_to_file",&Complex_nd_array::write_to_file)
    .def("read_from_file",&Complex_nd_array::read_from_file)
;

  //---------------------------------------------------------------------
  // Real_scalar_field
  //---------------------------------------------------------------------
  Nd_array<double>& (Real_scalar_field::*real_get_points)() =
    &Real_scalar_field::get_points;

  void (Real_scalar_field::*real_set_physical_params)
    (std::vector<double> physical_size, std::vector<double> physical_offset) =
    &Real_scalar_field::set_physical_params;

  std::vector<int> (Real_scalar_field::*real_get_leftmost_indices)
    (std::vector<double> location) = &Real_scalar_field::get_leftmost_indices;

  std::vector<double> (Real_scalar_field::*real_get_leftmost_offsets)
    (std::vector<double> location) = &Real_scalar_field::get_leftmost_offsets;

  class_<Real_scalar_field>("Real_scalar_field",init<>())
    .def(init<std::vector<int>, std::vector<double>, std::vector<double> >())
    .def("copy",&Real_scalar_field::copy<double>)
    .def("set_physical_params",real_set_physical_params)
    .def("get_physical_size",&Real_scalar_field::get_physical_size)
    .def("get_physical_offset",&Real_scalar_field::get_physical_offset)
    .def("get_cell_size",&Real_scalar_field::get_cell_size)
    .def("get_points",real_get_points,return_internal_reference<>())
    .def("get_leftmost_indices",real_get_leftmost_indices)
    .def("get_leftmost_offsets",real_get_leftmost_offsets)
    .def("write_to_file",&Real_scalar_field::write_to_file)
    .def("read_from_file",&Real_scalar_field::read_from_file)
    .def("describe",&Real_scalar_field::describe)
    ;

  //---------------------------------------------------------------------
  // Complex_scalar_field
  //---------------------------------------------------------------------
  Nd_array<std::complex<double> >& 
    (Complex_scalar_field::*complex_get_points)() =
    &Complex_scalar_field::get_points;

  void (Complex_scalar_field::*complex_set_physical_params)
    (std::vector<double> physical_size, std::vector<double> physical_offset) =
    &Complex_scalar_field::set_physical_params;

  std::vector<int> (Complex_scalar_field::*complex_get_leftmost_indices)
    (std::vector<double> location) = 
    &Complex_scalar_field::get_leftmost_indices;

  std::vector<double> (Complex_scalar_field::*complex_get_leftmost_offsets)
    (std::vector<double> location) = 
    &Complex_scalar_field::get_leftmost_offsets;

  class_<Complex_scalar_field>("Complex_scalar_field",init<>())
    .def(init<std::vector<int>, std::vector<double>, std::vector<double> >())
    .def("copy",&Complex_scalar_field::copy<std::complex<double> >)
    .def("set_physical_params",complex_set_physical_params)
    .def("get_physical_size",&Complex_scalar_field::get_physical_size)
    .def("get_physical_offset",&Complex_scalar_field::get_physical_offset)
    .def("get_cell_size",&Complex_scalar_field::get_cell_size)
    .def("get_points",complex_get_points,return_internal_reference<>())
    .def("get_leftmost_indices",complex_get_leftmost_indices)
    .def("get_leftmost_offsets",complex_get_leftmost_offsets)
    .def("write_to_file",&Complex_scalar_field::write_to_file)
    .def("read_from_file",&Complex_scalar_field::read_from_file)
    .def("describe",&Complex_scalar_field::describe)
    ;

  //---------------------------------------------------------------------
  // deposit functions
  //---------------------------------------------------------------------
  def("deposit_charge_cic",deposit_charge_cic);
  def("deposit_charge_ngp",deposit_charge_ngp);
  
  //---------------------------------------------------------------------
  // solvers
  //---------------------------------------------------------------------
  def("solver",solver);
  def("fft_tester",fft_tester);
}


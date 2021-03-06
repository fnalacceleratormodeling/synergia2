#include "parallel_utils.h"
#include "synergia/utils/comm_converter.h"
#include "synergia/utils/commxx_divider.h"
#include "synergia/utils/container_conversions.h"
#include "synergia/utils/hdf5_file.h"
#include "synergia/utils/logger.h"
#include <boost/python.hpp>

using namespace boost::python;

tuple
decompose_1d_raw_wrap(int processors, int length)
{
  std::vector<int> counts(processors), offsets(processors);
  decompose_1d_raw(processors, length, offsets, counts);
  return make_tuple(
    container_conversions::to_tuple<std::vector<int>>::convert_tuple(offsets),
    container_conversions::to_tuple<std::vector<int>>::convert_tuple(counts));
}

tuple
decompose_1d_wrap(Commxx& comm, int length)
{
  int processors = comm.get_size();
  std::vector<int> counts(processors), offsets(processors);
  decompose_1d(comm, length, offsets, counts);
  return make_tuple(
    container_conversions::to_tuple<std::vector<int>>::convert_tuple(offsets),
    container_conversions::to_tuple<std::vector<int>>::convert_tuple(counts));
}

BOOST_PYTHON_MODULE(parallel_utils)
{
  if (import_mpi4py() < 0) { return; }

  container_conversions::from_python_sequence<
    std::vector<int>,
    container_conversions::variable_capacity_policy>();

  to_python_converter<std::vector<int>,
                      container_conversions::to_tuple<std::vector<int>>>();

  class_<Commxx, Commxx_sptr, boost::noncopyable>("Commxx", init<>())
    .def(init<bool>())
    .def(init<Commxx_sptr, std::vector<int> const&, optional<bool>>())
    .def("get_rank", &Commxx::get_rank)
    .def("get_size", &Commxx::get_size)
    .def("has_this_rank", &Commxx::has_this_rank);
  register_ptr_to_python<Commxx_sptr>();

  container_conversions::from_python_sequence<
    std::vector<Commxx_sptr>,
    container_conversions::variable_capacity_policy>();

  to_python_converter<
    std::vector<Commxx_sptr>,
    container_conversions::to_tuple<std::vector<Commxx_sptr>>>();

  def("generate_subcomms", generate_subcomms);
  def("make_optimal_spc_comm", make_optimal_spc_comm);
  def("decompose_1d_raw", decompose_1d_raw_wrap);
  def("decompose_1d", decompose_1d_wrap);
  def("decompose_1d_local", decompose_1d_local);

  class_<Logger>("Logger", init<int>())
    .def(init<int, bool>())
    .def(init<int, std::string const&>())
    .def(init<int, std::string const&, bool>())
    .def(init<int, std::string const&, bool, bool>())
    .def(init<std::string const&>())
    .def(init<std::string const&, bool>())
    .def("write", &Logger::write, return_internal_reference<>())
    .def("flush", &Logger::flush, return_internal_reference<>());
  class_<Commxx_divider, Commxx_divider_sptr>("Commxx_divider",
                                              init<int, bool>())
    .def("get_commxx", &Commxx_divider::get_commxx_sptr);

  void (Hdf5_file::*write_double_a)(double const&, std::string const&) =
    &Hdf5_file::write<double>;
  void (Hdf5_file::*write_double_b)(double const*, size_t, std::string const&) =
    &Hdf5_file::write<double>;
  void (Hdf5_file::*write_marray1d_ref)(
    MArray1d_ref const&, std::string const&) = &Hdf5_file::write<MArray1d_ref>;
  void (Hdf5_file::*write_marray2d_ref)(
    MArray2d_ref const&, std::string const&) = &Hdf5_file::write<MArray2d_ref>;
  void (Hdf5_file::*write_marray3d_ref)(
    MArray3d_ref const&, std::string const&) = &Hdf5_file::write<MArray3d_ref>;

  scope Hdf5_file_scope =
    class_<Hdf5_file, Hdf5_file_sptr>(
      "Hdf5_file", init<std::string const&, Hdf5_file::Flag>())
      .def("open", &Hdf5_file::open)
      .def("close", &Hdf5_file::close)
      .def("flush", &Hdf5_file::flush)
      .def("write", write_double_a)
      .def("write", write_double_b)
      .def("write", write_marray1d_ref)
      .def("write", write_marray2d_ref)
      .def("write", write_marray3d_ref)
      .def("read_double", &Hdf5_file::read<double>)
      .def("read_int", &Hdf5_file::read<int>)
      .def("read_array1d", &Hdf5_file::read<MArray1d>)
      .def("read_array2d", &Hdf5_file::read<MArray2d>)
      .def("read_array3d", &Hdf5_file::read<MArray3d>)
      .def("read_array1i", &Hdf5_file::read<MArray1i>)
      .def("get_member_names", &Hdf5_file::get_member_names)
      .def("get_atomic_type", &Hdf5_file::get_atomic_type)
      .def("get_dims", &Hdf5_file::get_dims)
      .def("close", &Hdf5_file::flush);
  enum_<Hdf5_file::Flag>("Flag")
    .value("truncate", Hdf5_file::truncate)
    .value("read_write", Hdf5_file::read_write)
    .value("read_only", Hdf5_file::read_only)
    .export_values();
  enum_<Hdf5_file::Atomic_type>("Atomic_type")
    .value("double_type", Hdf5_file::double_type)
    .value("int_type", Hdf5_file::int_type)
    .export_values();
}

//#undef _POSIX_C_SOURCE
#include <boost/python.hpp>
#include "array_nd/array_nd.h"
#include "array_nd_python.h"
#include "container_conversions.h"

#include "uses_array_nd.h"
using namespace boost::python;

#include <iostream>

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/to_python_converter.hpp>

struct Array_nd_double_to_numpy
{
    static PyObject* convert(Array_nd<double> const& x)
    {
  	  return PyObject_from_Array_nd(x);
    }
};

struct array_nd_double_from_numpy
{
  array_nd_double_from_numpy()
  {
    boost::python::converter::registry::push_back(
      &convertible,
      &construct,
      boost::python::type_id<Array_nd<double> >());
  }

  static void* convertible(PyObject* obj_ptr)
  {
    if (!PyArray_Check(obj_ptr)) return 0;
    return obj_ptr;
  }

  static void construct(
    PyObject* obj_ptr,
    boost::python::converter::rvalue_from_python_stage1_data* data)
  {
    std::vector<int> shape(PyArray_NDIM(obj_ptr));
    std::vector<int> strides(PyArray_NDIM(obj_ptr));
    for (int i = 0; i < PyArray_NDIM(obj_ptr) ; ++i) {
        shape.at(i) = PyArray_DIMS(obj_ptr)[i];
        strides.at(i) = PyArray_STRIDES(obj_ptr)[i]/sizeof(double);
    }

    void* storage = (
      (boost::python::converter::rvalue_from_python_storage<Array_nd<double> >*)
        data)->storage.bytes;
    new (storage) Array_nd<double>(shape,strides,reinterpret_cast<double*>(PyArray_DATA(obj_ptr)));
    data->convertible = storage;
  }
};

BOOST_PYTHON_MODULE(uses_array_nd)
{
    scitbx::boost_python::container_conversions::from_python_sequence <
    std::vector<int>,
    scitbx::boost_python::container_conversions::variable_capacity_policy > ();

    boost::python::to_python_converter <
    std::vector<int>,
    scitbx::boost_python::container_conversions::to_tuple <
    std::vector<int> > > ();

	import_array();

    to_python_converter<Array_nd<double>,Array_nd_double_to_numpy >(); //"false" because tag_to_noddy has no member get_pytype
    array_nd_double_from_numpy();

	def("create",create);
	def("inspect",inspect);
	def("modify",modify);

	class_<Has_array_nd>("Has_array_nd",init<>())
		.def("get_x",&Has_array_nd::get_x)
		.def("check_x",&Has_array_nd::check_x)
		;

}


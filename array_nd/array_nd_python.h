#include <typeinfo>
#include <boost/python.hpp>
#include "array_nd/array_nd.h"
#include <Python.h>
#include <numpy/arrayobject.h>

template<class T>
Array_nd<T>
Array_nd_from_PyObject(const PyObject *obj)
{
    std::vector<int> shape(PyArray_NDIM(obj));
    std::vector<int> strides(PyArray_NDIM(obj));
    for (int i = 0; i < PyArray_NDIM(obj) ; ++i) {
        shape.at(i) = PyArray_DIMS(obj)[i];
        strides.at(i) = PyArray_STRIDES(obj)[i]/sizeof(T);
    }
    Array_nd<T> retval(shape,strides,reinterpret_cast<T*>(PyArray_DATA(obj)));
    return retval;
}

template<class T>
PyObject *
PyObject_from_Array_nd(const Array_nd<T> &array)
{
    int type_num;
    if (typeid(T) == typeid(double))  type_num = NPY_DOUBLE;
    else {
        throw
            std::runtime_error("PyObject_from_Array_nd: Unable to deal with type");
    }
    std::vector<npy_intp> dims(array.get_rank());
    std::vector<npy_intp> strides(array.get_rank());
    for(int i=0; i<array.get_rank(); ++i) {
    	dims.at(i) = array.get_shape().at(i);
    	strides.at(i) = array.get_strides().at(i)*sizeof(T);
    }
    PyObject * retval = PyArray_New(&PyArray_Type,
    		array.get_rank(),
    		&dims[0],
    		type_num,
    		&strides[0],
    		static_cast<void*>(array.get_data_ptr()),
    		sizeof(T),
    		NPY_WRITEABLE,NULL);
    if (retval == 0) {
    	std::runtime_error("PyObject_from_Array_nd: PyArray_New failed.");
    }
    return retval;
}

struct Array_nd_double_to_numpy
{
    static PyObject* convert(Array_nd<double> const& x)
    {
  	  return PyObject_from_Array_nd(x);
    }
};

struct Array_nd_double_from_numpy
{
  Array_nd_double_from_numpy()
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


// Usage example for converter structs:
//BOOST_PYTHON_MODULE(uses_array_nd)
//{
//	import_array();
//    to_python_converter<Array_nd<double>,Array_nd_double_to_numpy>();
//    Array_nd_double_from_numpy();
//}


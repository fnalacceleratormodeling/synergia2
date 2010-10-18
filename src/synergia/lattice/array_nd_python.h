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
    //~ switch (typeid(T)) {
        //~ case typeid(int):
            //~ type_num = PyArray_INT;
            //~ break;
        //~ case typeid(long):
            //~ type_num = PyArray_LONG;
            //~ break;
        //~ case typeid(float):
            //~ type_num = PyArray_FLOAT;
            //~ break;
        //~ case typeid(double):
            //~ type_num = PyArray_DOUBLE;
            //~ break;
        //~ case typeid(std::complex<float>):
            //~ type_num = PyArray_CFLOAT;
            //~ break;
        //~ case typeid(std::complex<double>):
            //~ type_num = PyArray_CDOUBLE;
            //~ break;
        //~ default:
            //~ throw
                //~ std::runtime_error("PyObject_from_Array_nd unable to deal with given type");
    //~ }
    if (typeid(T) == typeid(double))  type_num = NPY_DOUBLE;
    else {
        throw
            std::runtime_error("PyObject_from_Array_nd unable to deal with given type");
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
    	std::runtime_error("PyArray_New failed.");
    }
    return retval;
}

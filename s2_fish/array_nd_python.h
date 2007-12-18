#include <typeinfo>
#include <boost/python.hpp>
#include "array_2d.h"
#include <Numeric/arrayobject.h>

#define PyArray_NDIM(obj) (((PyArrayObject *)(obj))->nd)
#define PyArray_DIMS(obj) (((PyArrayObject *)(obj))->dimensions)
#define PyArray_DATA(obj) ((void *)(((PyArrayObject *)(obj))->data))
template<class T>
Array_nd<T>
Array_nd_from_PyObject(const PyObject *obj)
{
    std::vector<int> shape;
    shape.resize(PyArray_NDIM(obj));
    shape.at(0) = PyArray_DIMS(obj)[0];
    int size = shape.at(0);
    for (int i = 1; i < PyArray_NDIM(obj) ; ++i) {
        shape.at(i) = PyArray_DIMS(obj)[i];
        size *= shape.at(i);
    }
    Array_nd<T> retval(shape,reinterpret_cast<T*>(PyArray_DATA(obj)));
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
    if (typeid(T) == typeid(int))  type_num = PyArray_INT;
    else {
        throw    
            std::runtime_error("PyObject_from_Array_nd unable to deal with given type");
    }
    PyObject *retval = PyArray_FromDims(array.get_rank(),
            &array.get_shape()[0], type_num);

    T *data_out = reinterpret_cast<T*> PyArray_DATA(retval);
    T *data_in = array.get_data_ptr();
    int size = array.get_size();
    for (int i=0; i<size; ++i) {
        data_out[i] = data_in[i];
    }
    return retval;
}

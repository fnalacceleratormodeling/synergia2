#ifndef NUMPY_MULTI_ARRAY_CONVERTER_H_
#define NUMPY_MULTI_ARRAY_CONVERTER_H_

#define NPY_NO_DEPRECATED_API 7
#include "numpy/arrayobject.h"

template<typename ValueType, int Dimension>
    struct numpy_multi_array_converter
    {
        typedef boost::multi_array<ValueType, Dimension > multi_array_t;
        typedef std::vector<std::size_t > shape_t;

        static void
        register_to_python()
        {
            boost::python::to_python_converter<multi_array_t,
                    numpy_multi_array_converter<ValueType, Dimension > >();
        }

        static
        void *
        convertible(PyObject * obj)
        {
            using namespace boost::python;
            try {
                shape_t shape;
                get_shape(object(handle< >(borrowed(obj))), shape);
                if (multi_array_t::dimensionality != shape.size()) return 0;
            }
            catch (...) {
                return 0;
            }
            return obj;
        }

        static PyObject *
        convert(const multi_array_t & c_array)
        {
            using namespace boost::python;
            PyObject *retval;
            int type_num;
            if (typeid(ValueType) == typeid(double)) type_num = NPY_DOUBLE;
            else {
                throw std::runtime_error(
                        "numpy_multi_array_converter: Unable to deal with type");
            }
            retval = PyArray_SimpleNew(c_array.num_dimensions(),
                    (npy_intp*) c_array.shape(), type_num);
            double * p = (double *) PyArray_DATA((PyArrayObject *)retval);
            for (const double * it = c_array.data();
                    it != (c_array.data() + c_array.num_elements()); ++it) {
                *p = *it;
                ++p;
            }
            return retval;
        }

    protected:
        static
        void
        get_shape(boost::python::object obj, shape_t & shape)
        {
            using namespace boost::python;
            shape.clear();
            object py_shape = obj.attr("shape");
            const std::size_t N = len(py_shape);
            for (std::size_t i = 0; N != i; ++i)
                shape.push_back(extract<std::size_t >(py_shape[i]));
        }
    };

#endif /* NUMPY_MULTI_ARRAY_CONVERTER_H_ */

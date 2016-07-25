#ifndef NUMPY_MULTI_ARRAY_CONVERTER_H_
#define NUMPY_MULTI_ARRAY_CONVERTER_H_

#define NPY_NO_DEPRECATED_API 7
#include "numpy/arrayobject.h"

// backwards compatibility
#ifndef NPY_ARRAY_F_CONTIGUOUS
    #define NPY_ARRAY_F_CONTIGUOUS NPY_F_CONTIGUOUS
#endif
#ifndef NPY_ARRAY_WRITEABLE
    #define NPY_ARRAY_WRITEABLE NPY_WRITEABLE
#endif

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
            if (typeid(ValueType) == typeid(double)) {
                type_num = NPY_DOUBLE;
            } else if(typeid(ValueType) == typeid(int)) {
                type_num = NPY_INT;
            } else {
                throw std::runtime_error(
                        "numpy_multi_array_converter: Unable to deal with type");
            }
            bool c_order = true;
            if (c_array.storage_order() == boost::c_storage_order()) {
                c_order = true;
            } else {
                if (c_array.storage_order() == boost::fortran_storage_order()) {
                    c_order = false;
                } else {
                    throw std::runtime_error(
                                "numpy_multi_array_converter: Unable to deal with general storage order");
                }
            }
            int ndim = c_array.num_dimensions();
            if (ndim > 3) {
                throw std::runtime_error(
                        "numpy_multi_array_converter: Unable to deal with arrays of the dimension higher than 3");
            }
            retval = PyArray_SimpleNew(c_array.num_dimensions(),
                    (npy_intp*) c_array.shape(), type_num);
            ValueType * p = (ValueType *) PyArray_DATA((PyArrayObject *)retval);
            if (c_order || ndim == 1) {
                for (const ValueType * it = c_array.data();
                        it != (c_array.data() + c_array.num_elements()); ++it) {
                    *p = *it;
                    ++p;
                }
            } else {
                const ValueType * c = c_array.data();

                if (ndim == 2) {

                    int ni = c_array.shape()[0];
                    int nj = c_array.shape()[1];

                    for (int i=0; i<ni; ++i) {
                        for (int j=0; j<nj; ++j) {
                            p[i*nj+j] = c[j*ni+i];
                        }
                    }
                } else {

                    int ni = c_array.shape()[0];
                    int nj = c_array.shape()[1];
                    int nk = c_array.shape()[2];

                    for (int i=0; i<ni; ++i) {
                        for (int j=0; j<nj; ++j) {
                            for (int k=0; k<nk; ++k) {
                                p[i*nj*nk + j*nk +k] = c[k*ni*nj + j*ni +i];
                            }
                        }
                    }

                }
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

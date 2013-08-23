#ifndef NUMPY_MULTI_REF_CONVERTER_H_
#define NUMPY_MULTI_REF_CONVERTER_H_

#define NPY_NO_DEPRECATED_API 7
#include "numpy/arrayobject.h"

template<typename ValueType, int Dimension>
    struct numpy_multi_array_ref_converter
    {
        typedef boost::multi_array_ref<ValueType, Dimension > multi_array_ref_t;
        typedef std::vector<std::size_t > shape_t;

        static void
        register_to_and_from_python()
        {
            register_from_python();
            register_to_python();
        }

        static void
        register_to_python()
        {
            boost::python::to_python_converter<multi_array_ref_t,
                    numpy_multi_array_ref_converter<ValueType, Dimension > >();
        }

        static void
        register_from_python()
        {
            boost::python::converter::registry::push_back(
                    &numpy_multi_array_ref_converter<ValueType, Dimension >::convertible,
                    &numpy_multi_array_ref_converter<ValueType, Dimension >::construct,
                    boost::python::type_id<multi_array_ref_t >());
        }

        static
        void *
        convertible(PyObject * obj)
        {
            using namespace boost::python;
            try {
                shape_t shape;
                get_shape(object(handle< > (borrowed(obj))), shape);
                if (multi_array_ref_t::dimensionality != shape.size()) return 0;
            }
            catch (...) {
                return 0;
            }
            return obj;
        }

        static
        void
        construct(PyObject* obj,
                boost::python::converter::rvalue_from_python_stage1_data* data)
        {
            using namespace boost::python;

            //get the storage
            typedef converter::rvalue_from_python_storage<multi_array_ref_t >
                    storage_t;
            storage_t * the_storage = reinterpret_cast<storage_t * > (data);
            void * memory_chunk = the_storage->storage.bytes;

            //new placement
            object py_obj(handle< > (borrowed(obj)));
            shape_t shape;
            get_shape(py_obj, shape);
            new (memory_chunk) multi_array_ref_t(
                    reinterpret_cast<ValueType* > (PyArray_DATA((PyArrayObject *)obj)), shape);
            data->convertible = memory_chunk;
        }

        static PyObject *
        convert(const multi_array_ref_t & c_array)
        {
            using namespace boost::python;
            PyObject *retval;
            int type_num;
            if (typeid(ValueType) == typeid(double)) type_num = NPY_DOUBLE;
            else {
                throw std::runtime_error(
                        "numpy_multi_array_ref_converter: Unable to deal with type");
            }
            retval = PyArray_SimpleNewFromData(c_array.num_dimensions(),
                    (npy_intp*) c_array.shape(), type_num,
                    (void *) c_array.origin());
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
                shape.push_back(extract<std::size_t > (py_shape[i]));
        }
    };

// jfa: I'm not sure that everything in this converter makes sense.
//      Furthermore, I don't know that it can't be an extension of the
//      above struct (numpy_multi_array_ref_converter).
template<typename ValueType, int Dimension>
    struct numpy_const_multi_array_ref_converter
    {
        typedef boost::const_multi_array_ref<ValueType, Dimension >
                const_multi_array_ref_t;
        typedef std::vector<std::size_t > shape_t;

        static void
        register_to_and_from_python()
        {
            register_from_python();
            register_to_python();
        }

        static void
        register_to_python()
        {
            boost::python::to_python_converter<const_multi_array_ref_t,
                    numpy_const_multi_array_ref_converter<ValueType, Dimension > >();
        }

        static void
        register_from_python()
        {
            boost::python::converter::registry::push_back(
                    &numpy_const_multi_array_ref_converter<ValueType, Dimension >::convertible,
                    &numpy_const_multi_array_ref_converter<ValueType, Dimension >::construct,
                    boost::python::type_id<const_multi_array_ref_t >());
        }

        static
        void *
        convertible(PyObject * obj)
        {
            using namespace boost::python;
            try {
                shape_t shape;
                get_shape(object(handle< > (borrowed(obj))), shape);
                if (const_multi_array_ref_t::dimensionality != shape.size()) return 0;
            }
            catch (...) {
                return 0;
            }
            return obj;
        }

        static
        void
        construct(PyObject* obj,
                boost::python::converter::rvalue_from_python_stage1_data* data)
        {
            using namespace boost::python;

            //get the storage
            typedef converter::rvalue_from_python_storage<const_multi_array_ref_t >
                    storage_t;
            storage_t * the_storage = reinterpret_cast<storage_t * > (data);
            void * memory_chunk = the_storage->storage.bytes;

            //new placement
            object py_obj(handle< > (borrowed(obj)));
            shape_t shape;
            get_shape(py_obj, shape);
            new (memory_chunk) const_multi_array_ref_t(
                    reinterpret_cast<ValueType* > (PyArray_DATA((PyArrayObject *)obj)), shape);
            data->convertible = memory_chunk;
        }

        static PyObject *
        convert(const const_multi_array_ref_t & c_array)
        {
            using namespace boost::python;
            PyObject *retval;
            int type_num;
            if (typeid(ValueType) == typeid(double)) type_num = NPY_DOUBLE;
            else {
                throw std::runtime_error(
                        "numpy_const_multi_array_ref_converter: Unable to deal with type");
            }
            retval = PyArray_SimpleNewFromData(c_array.num_dimensions(),
                    (npy_intp*) c_array.shape(), type_num,
                    (void *) c_array.origin());
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
                shape.push_back(extract<std::size_t > (py_shape[i]));
        }
    };
#endif /* NUMPY_MULTI_REF_CONVERTER_H_ */

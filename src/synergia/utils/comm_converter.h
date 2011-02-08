#ifndef COMM_CONVERTER_H_
#define COMM_CONVERTER_H_

// begin workaround for using mpi4py using MPI-1
#include <mpi.h>
#if MPI_VERSION < 2
  typedef void* MPI_Win;
  typedef void* MPI_File;
#endif
// end workaround

#include <mpi4py/mpi4py.h>

#include <boost/python.hpp>
#include "synergia/utils/commxx.h"

struct comm_converter
{
    static void
    register_to_and_from_python()
    {
        register_from_python();
        register_to_python();
    }

    static void
    register_to_python()
    {
        boost::python::to_python_converter<Commxx, comm_converter >();
    }

    static void
    register_from_python()
    {
        boost::python::converter::registry::push_back(&convertible, &construct,
                boost::python::type_id<Commxx >());
    }

    static
    void *
    convertible(PyObject * obj)
    {
        if (!PyObject_TypeCheck(obj, &PyMPIComm_Type)) {
            return 0;
        } else {
            return obj;
        }
    }

    static
    void
    construct(PyObject* obj,
            boost::python::converter::rvalue_from_python_stage1_data* data)
    {
        using namespace boost::python;
        MPI_Comm *comm_p = PyMPIComm_Get(obj);
        if (comm_p == NULL) {
            throw_error_already_set();
        }
        void
                *storage =
                        ((converter::rvalue_from_python_storage<Commxx >*) data)->storage.bytes;

        new (storage) Commxx(*comm_p);
        data->convertible = storage;
    }

    static PyObject *
    convert(Commxx const& comm)
    {
        using namespace boost::python;
        PyObject *retval;
        MPI_Comm mpi_comm;
        mpi_comm = comm.get();
        retval = PyMPIComm_New(mpi_comm);
        return retval;
    }
};

#endif /* COMM_CONVERTER_H_ */

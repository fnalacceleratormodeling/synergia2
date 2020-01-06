#ifndef DIAGNOSTICS_PY_H
#define DIAGNOSTICS_PY_H

#include <pybind11/pybind11.h>
#include "synergia/bunch/diagnostics.h"

class PyDiagnostics : public Diagnostics
{
public:

    // inherid the constructors
    using Diagnostics::Diagnostics;

    PyDiagnostics(std::string const& type = "PyDiagnostics", bool serial = true)
        : Diagnostics(type, serial)
    { }

private:

    void do_update(Bunch const& bunch) override
    {
        PYBIND11_OVERLOAD_PURE(
                void,  // return type
                Diagnostics, // parent class
                do_update,   // name of the function in c++
                bunch        // arguments
        );
    }

    void do_collect(Commxx comm, int root) override
    {
        PYBIND11_OVERLOAD_PURE(
                void,  // return type
                Diagnostics, // parent class
                do_collect,
                comm,
                root
        );
    }

    void do_write(Hdf5_file& file, bool first_write) override
    { }

#if 0
    void do_write(Hdf5_file& file, bool first_write) override
    {
        PYBIND11_OVERLOAD_PURE(
                void,  // return type
                Diagnostics, // parent class
                do_write,
                hdf5_file,
                first_write 
        );
    }
#endif

};

#endif

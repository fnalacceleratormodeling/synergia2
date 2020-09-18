#ifndef DIAGNOSTICS_PY_H
#define DIAGNOSTICS_PY_H

#include <pybind11/pybind11.h>
#include "synergia/bunch/diagnostics.h"
#include "synergia/utils/base64.h"

class PyDiagnostics : public Diagnostics
{
public:

    // inherit the constructors
    using Diagnostics::Diagnostics;

    PyDiagnostics(
            std::string const& type = "PyDiagnostics", 
            std::string const& filename = "py_diag.h5",
            bool serial = true)
        : Diagnostics(type, filename, serial), self()
    { }

    // hold a reference to the python instance so it wont go out of scope
    void reg_self()
    { self = pybind11::cast(this); }

private:

    pybind11::object self;

    void do_update(Bunch const& bunch) override
    {
        self.attr("do_update")(bunch);

#if 0
        PYBIND11_OVERLOAD_PURE(
                void,        // return type
                Diagnostics, // parent class
                do_update,   // name of the function in c++
                bunch        // arguments
        );
#endif
    }

    void do_reduce(Commxx const& comm, int root) override
    {
        self.attr("do_reduce")(comm, root);

#if 0
        PYBIND11_OVERLOAD_PURE(
                void,        // return type
                Diagnostics, // parent class
                do_reduce,
                comm,
                root
        );
#endif
    }

    void do_first_write(Hdf5_file& file) override
    { }

    void do_write(Hdf5_file& file) override
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

    friend class cereal::access;

    template<class Archive>
    void save(Archive & ar) const
    { 
        ar(cereal::base_class<Diagnostics>(this)); 

        namespace py = pybind11;
        py::object pickle;
        try { pickle = py::module::import("cPickle"); }
        catch(...) { pickle = py::module::import("pickle"); }

        auto s = pickle.attr("dumps")(self, 2).cast<std::string>();
        auto base_s = base64_encode(reinterpret_cast<const unsigned char*>(s.c_str()), s.length());

        ar(base_s);
    }

    template<class Archive>
    void load(Archive & ar)
    {
        ar(cereal::base_class<Diagnostics>(this)); 

        std::string base_s;
        ar(base_s);

        std::string s = base64_decode(base_s);

        namespace py = pybind11;
        py::object pickle;
        try { pickle = py::module::import("cPickle"); }
        catch(...) { pickle = py::module::import("pickle"); }

        self = pickle.attr("loads")(py::bytes(s));
    }

};

CEREAL_REGISTER_TYPE(PyDiagnostics)
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(PyDiagnostics, cereal::specialization::member_load_save)

#endif

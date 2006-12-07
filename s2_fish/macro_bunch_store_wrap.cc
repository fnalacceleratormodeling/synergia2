#include "macro_bunch_store.h"
#undef _POSIX_C_SOURCE
#include <boost/python.hpp>

using namespace boost::python;

BOOST_PYTHON_MODULE(macro_bunch_store)
{
  class_<Macro_bunch_store>("Macro_bunch_store",
			    init<numeric::array&,int,int,double,
			    numeric::array&,numeric::array&,bool>())
    .def_readonly("local_particles",
		  &Macro_bunch_store::numeric_local_particles)
    .def_readwrite("local_num",&Macro_bunch_store::local_num)
    .def_readwrite("total_num",&Macro_bunch_store::total_num)
    .def_readwrite("total_current",&Macro_bunch_store::total_current)
    .def_readonly("units",&Macro_bunch_store::numeric_units)
    .def_readonly("ref_particle",&Macro_bunch_store::numeric_ref_particle)
    .def_readwrite("is_z",&Macro_bunch_store::is_z)
  ;
}


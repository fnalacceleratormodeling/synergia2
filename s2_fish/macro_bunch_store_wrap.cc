#include "macro_bunch_store.h"
#undef _POSIX_C_SOURCE
#include <boost/python.hpp>

using namespace boost::python;

#include <iostream>

BOOST_PYTHON_MODULE(macro_bunch_store)
{
  numeric::array::set_module_and_type("Numeric", "ArrayType");
  class_<Macro_bunch_store>("Macro_bunch_store",
			    init<numeric::array&,int,int,double,
			    numeric::array&,numeric::array&,bool>())
    .def("get_local_particles",
		  &Macro_bunch_store::get_local_particles)
    .def("get_units",
		  &Macro_bunch_store::get_units)
    .def("get_ref_particle",
		  &Macro_bunch_store::get_ref_particle)
    .def("check",&Macro_bunch_store::check)
    .def_readwrite("local_num",&Macro_bunch_store::local_num)
    .def_readwrite("total_num",&Macro_bunch_store::total_num)
    .def_readwrite("total_current",&Macro_bunch_store::total_current)
    .def_readwrite("is_fixedz",&Macro_bunch_store::is_fixedz)
    .def("convert_to_fixedt",&Macro_bunch_store::convert_to_fixedt)
    .def("convert_to_fixedz",&Macro_bunch_store::convert_to_fixedz)
    .def("get_coord",&Macro_bunch_store::get_coord)
  ;
}


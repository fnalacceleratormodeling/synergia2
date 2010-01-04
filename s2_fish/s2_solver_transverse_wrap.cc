#include "solver_transverse_fftw.h"
#undef _POSIX_C_SOURCE
#include <boost/python.hpp>
#include <array_nd/array_nd_python.h>
#include "container_conversions.h"

using namespace boost::python;

BOOST_PYTHON_MODULE(s2_solver_transverse_fftw)
{

    import_array();
    to_python_converter<Array_nd<double>,Array_nd_double_to_numpy>();
    Array_nd_double_from_numpy();

    scitbx::boost_python::container_conversions::from_python_sequence <
    std::vector<double>,
    scitbx::boost_python::container_conversions::variable_capacity_policy > ();
       
    class_<Extremum>("Extremum")
	.add_property("min", &Extremum::min)
	.add_property("max", &Extremum::max)
	;
    
    class_<TransverseSolver>("TransverseSolver", init<int,int,double,bool>() )
	.def("kick_transverse_charge", &TransverseSolver::kick_transverse_charge)
	.def("kick_transverse_charge_direct", &TransverseSolver::kick_transverse_charge_direct) 
	.def("get_xGrid", &TransverseSolver::getXGrid, return_value_policy<copy_const_reference>() )
	.def("get_yGrid", &TransverseSolver::getYGrid, return_value_policy<copy_const_reference>() )
	.def("get_rho", &TransverseSolver::getRho, return_value_policy<copy_const_reference>() )
	.def("get_Fscx", &TransverseSolver::getFscx, return_value_policy<copy_const_reference>() )
	.def("get_Fscy", &TransverseSolver::getFscy, return_value_policy<copy_const_reference>() )
	.def("find_bunch_extrema", &TransverseSolver::findBunchExtrema)
	.staticmethod("find_bunch_extrema")
	;
}



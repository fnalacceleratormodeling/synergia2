#include <iostream>

#include "fftw_helper.h"
#include "nd_array.h"
#include "scalar_field.h"
#include "deposit.h"
#include "solver_fftw.h"
#undef _POSIX_C_SOURCE
#include <boost/python.hpp>

#include "container_conversions.h"
#include "communicate.h"

using namespace boost::python;

BOOST_PYTHON_MODULE(s2_solver_fftw)
{
    class_<Fftw_helper>("Fftw_helper", init<std::vector<int>, bool >());
    def("solver_fftw_open", solver_fftw_open);
    def("gather_global_rho",gather_global_rho);
}


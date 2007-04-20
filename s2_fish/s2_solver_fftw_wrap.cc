#include "nd_array.h"
#include "scalar_field.h"
#include "deposit.h"
#include "solver_fftw.h"
#undef _POSIX_C_SOURCE
#include <boost/python.hpp>

#include "container_conversions.h"

using namespace boost::python;

BOOST_PYTHON_MODULE(s2_solver_fftw)
{
  //---------------------------------------------------------------------
  // solvers
  //---------------------------------------------------------------------
  //  def("fft_tester",fft_tester);
  def("solver_fftw_open",solver_fftw_open);
}


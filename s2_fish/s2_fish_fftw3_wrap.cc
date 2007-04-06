#include "nd_array.h"
#include "scalar_field.h"
#include "deposit.h"
#include "solver_fftw3.h"
#undef _POSIX_C_SOURCE
#include <boost/python.hpp>

#include "container_conversions.h"

using namespace boost::python;

BOOST_PYTHON_MODULE(s2_fish_fftw3)
{
  //---------------------------------------------------------------------
  // solvers
  //---------------------------------------------------------------------
  //  def("fft_tester",fft_tester);
  def("solver_fftw3_open",solver_fftw3::solver_fftw3_open);
  def("full_kick",solver_fftw3::full_kick);
  def("fft_tester",solver_fftw3::fft_tester);
}


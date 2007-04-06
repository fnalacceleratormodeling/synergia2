#include "solver.h"
#undef _POSIX_C_SOURCE
#include <boost/python.hpp>

#include "container_conversions.h"

using namespace boost::python;

BOOST_PYTHON_MODULE(s2_solver)
{
  def("fft_tester",fft_tester);
  def("solver_fft_open",solver_fft_open);
  def("solver_fd_multigrid_open",solver_fd_multigrid_open);
}


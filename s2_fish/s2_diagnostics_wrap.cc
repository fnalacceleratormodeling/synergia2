#include "s2_diagnostics.h"
#undef _POSIX_C_SOURCE
#include <boost/python.hpp>

using namespace boost::python;

BOOST_PYTHON_MODULE(s2_diagnostics)
{
  def("get_spatial_means_stds",get_spatial_means_stds);
  def("get_moments_corrs",get_moments_corrs);
}


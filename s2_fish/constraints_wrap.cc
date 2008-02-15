#include "constraints.h"
#undef _POSIX_C_SOURCE
#include <boost/python.hpp>

using namespace boost::python;

BOOST_PYTHON_MODULE(constraints)
{
  def("apply_longitudinal_periodicity",apply_longitudinal_periodicity);
  def("apply_circular_aperture",apply_circular_aperture);
}

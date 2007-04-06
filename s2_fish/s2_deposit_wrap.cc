#include "deposit.h"
#undef _POSIX_C_SOURCE
#include <boost/python.hpp>

using namespace boost::python;

BOOST_PYTHON_MODULE(s2_deposit)
{
  def("deposit_charge_cic",deposit_charge_cic);
  def("deposit_charge_ngp",deposit_charge_ngp);
}


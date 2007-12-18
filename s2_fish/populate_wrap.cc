#include "populate.h"
#undef _POSIX_C_SOURCE
#include <boost/python.hpp>
#include "array_nd_python.h"

using namespace boost::python;

void
populate_6d_gaussian_wrapper(object &particles, 
    const object &means, const object &covariances, const int id_offset)
{
    Array_2d<double> particles_array = 
        Array_nd_from_PyObject<double>(particles.ptr());
    Array_1d<double> means_array = 
        Array_nd_from_PyObject<double>(means.ptr());
    Array_2d<double> covariances_array = 
        Array_nd_from_PyObject<double>(covariances.ptr());
    populate_6d_gaussian(particles_array, means_array, covariances_array, id_offset);
}

BOOST_PYTHON_MODULE(populate)
{
  def("populate_6d_gaussian",populate_6d_gaussian_wrapper);
}

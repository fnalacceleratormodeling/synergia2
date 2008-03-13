#include "populate.h"
#undef _POSIX_C_SOURCE
#include <boost/python.hpp>
#include "array_nd/array_nd_python.h"

using namespace boost::python;

void
populate_6d_gaussian_wrapper(object &particles,
                             const object &means, const object &covariances, const int id_offset,
                             const unsigned long int seed, bool init_generator)
{
    Array_2d<double> particles_array =
        Array_nd_from_PyObject<double>(particles.ptr());
    Array_1d<double> means_array =
        Array_nd_from_PyObject<double>(means.ptr());
    Array_2d<double> covariances_array =
        Array_nd_from_PyObject<double>(covariances.ptr());
    populate_6d_gaussian(particles_array, means_array, covariances_array,
                         id_offset, seed, init_generator);
}

void
populate_6d_gaussian_quasi_wrapper(object &particles,
                             const object &means, const object &covariances, const int id_offset)
{
    Array_2d<double> particles_array =
        Array_nd_from_PyObject<double>(particles.ptr());
    Array_1d<double> means_array =
        Array_nd_from_PyObject<double>(means.ptr());
    Array_2d<double> covariances_array =
        Array_nd_from_PyObject<double>(covariances.ptr());
    populate_6d_gaussian_quasi(particles_array, means_array, covariances_array,
                         id_offset);
}

void
populate_transverse_gaussian_wrapper(object &particles,
                                     const object &means, const object &covariances, const int id_offset,
                                     const unsigned long int seed, bool init_generator)
{
    Array_2d<double> particles_array =
        Array_nd_from_PyObject<double>(particles.ptr());
    Array_1d<double> means_array =
        Array_nd_from_PyObject<double>(means.ptr());
    Array_2d<double> covariances_array =
        Array_nd_from_PyObject<double>(covariances.ptr());
    populate_transverse_gaussian(particles_array, means_array, covariances_array,
                                 id_offset, seed, init_generator);
}

void
populate_transverse_gaussian_wrapper_quasi(object &particles,
                                     const object &means, const object &covariances, const int id_offset)
{
    Array_2d<double> particles_array =
        Array_nd_from_PyObject<double>(particles.ptr());
    Array_1d<double> means_array =
        Array_nd_from_PyObject<double>(means.ptr());
    Array_2d<double> covariances_array =
        Array_nd_from_PyObject<double>(covariances.ptr());
    populate_transverse_gaussian_quasi(particles_array, means_array, covariances_array,
        id_offset);
}

BOOST_PYTHON_MODULE(populate)
{
    def("populate_6d_gaussian", populate_6d_gaussian_wrapper);
    def("populate_6d_gaussian_quasi", populate_6d_gaussian_quasi_wrapper);
    def("populate_transverse_gaussian", populate_transverse_gaussian_wrapper);
    def("populate_transverse_gaussian_quasi", populate_transverse_gaussian_wrapper_quasi);
}

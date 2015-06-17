#include <boost/python.hpp>
#include "synergia/utils/multi_array_typedefs.h"
#include "synergia/utils/numpy_multi_array_converter.h"
#include "synergia/utils/numpy_multi_ref_converter.h"

MArray1d get_MArray1d(int n1)
{
    MArray1d a(boost::extents[n1]);
    for(int i = 0; i < n1; ++i) {
        a[i] = i + 0.25;
    }

    return a;
}

MArray2d get_MArray2d(int n1, int n2)
{
    MArray2d a(boost::extents[n1][n2]);
    for(int i = 0; i < n1; ++i) {
        for(int j = 0; j < n2; ++j) {
            a[i][j] = i + 10*j + 0.25;
        }
    }

    return a;
}

MArray3d get_MArray3d(int n1, int n2, int n3)
{
    MArray3d a(boost::extents[n1][n2][n3]);
    for(int i = 0; i < n1; ++i) {
        for(int j = 0; j < n2; ++j) {
            for(int k = 0; k < n3; ++k) {
                a[i][j][k] = i + 10*j + 100*k + 0.25;
            }
        }
    }

    return a;
}

MArray2d get_MArray2d_fortran(int n1, int n2)
{
    MArray2d a(boost::extents[n1][n2],
               boost::fortran_storage_order());
    for(int i = 0; i < n1; ++i) {
        for(int j = 0; j < n2; ++j) {
            a[i][j] = i + 10*j + 0.25;
        }
    }

    return a;
}

MArray3d get_MArray3d_fortran(int n1, int n2, int n3)
{
    MArray3d a(boost::extents[n1][n2][n3],
               boost::fortran_storage_order());
    for(int i = 0; i < n1; ++i) {
        for(int j = 0; j < n2; ++j) {
            for(int k = 0; k < n3; ++k) {
                a[i][j][k] = i + 10*j + 100*k + 0.25;
            }
        }
    }

    return a;
}
using namespace boost::python;

BOOST_PYTHON_MODULE(test_helper)
{
    import_array();

    numpy_multi_array_converter<double, 1 >::register_to_python();
    numpy_multi_array_converter<double, 2 >::register_to_python();
    numpy_multi_array_converter<double, 3 >::register_to_python();

//    numpy_multi_array_ref_converter<double, 1 >::register_to_and_from_python();
//    numpy_const_multi_array_ref_converter<double, 1 >::register_to_and_from_python();
//    numpy_multi_array_ref_converter<double, 2 >::register_to_and_from_python();
//    numpy_const_multi_array_ref_converter<double, 2 >::register_to_and_from_python();

    def("get_MArray1d", get_MArray1d);
    def("get_MArray2d", get_MArray2d);
    def("get_MArray3d", get_MArray3d);
    def("get_MArray2d_fortran", get_MArray2d_fortran);
    def("get_MArray3d_fortran", get_MArray3d_fortran);
}

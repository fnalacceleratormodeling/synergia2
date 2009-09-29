#include "parallel_utils.h"
#include <boost/python.hpp>
#include "utils/container_conversions.h"

using namespace boost::python;

/////////////////////
/// jfa: this wrapping procedure does not do what it should.
////     the values are not passed back to python
void
decompose_1d_raw_wrap(int processors, int length, std::vector<int > offsets,
        std::vector<int > counts)
{
    decompose_1d_raw(processors,length,offsets,counts);
}

BOOST_PYTHON_MODULE(pyparallel_utils)
{
    container_conversions::from_python_sequence<std::vector<int >,
            container_conversions::variable_capacity_policy >();

    boost::python::to_python_converter<std::vector<int >,
            container_conversions::to_tuple<std::vector<int > > >();

    //    void
    //    decompose_1d_raw(int processors, int length, std::vector<int > &offsets,
    //            std::vector<int > &counts);
    //
    //    void
    //    decompose_1d(const MPI_Comm &comm, int length, std::vector<int > &offsets,
    //            std::vector<int > &counts);
    //
    //    int
    //    decompose_1d_local(const MPI_Comm &comm, int length);

    def("decompose_1d_raw", decompose_1d_raw_wrap);
}

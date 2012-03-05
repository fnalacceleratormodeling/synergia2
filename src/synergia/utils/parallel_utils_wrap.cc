#include "parallel_utils.h"
#include <boost/python.hpp>
#include "synergia/utils/container_conversions.h"
#include "synergia/utils/comm_converter.h"
#include "synergia/utils/logger.h"

using namespace boost::python;

tuple
decompose_1d_raw_wrap(int processors, int length)
{
    std::vector<int> counts(processors), offsets(processors);
    decompose_1d_raw(processors,length,offsets,counts);
    return make_tuple(
            container_conversions::to_tuple<std::vector<int> >::convert_tuple(offsets),
            container_conversions::to_tuple<std::vector<int> >::convert_tuple(counts));
}

tuple
decompose_1d_wrap(Commxx comm, int length)
{
    int processors = comm.get_size();
    std::vector<int> counts(processors), offsets(processors);
    decompose_1d(comm,length,offsets,counts);
    return make_tuple(
            container_conversions::to_tuple<std::vector<int> >::convert_tuple(offsets),
            container_conversions::to_tuple<std::vector<int> >::convert_tuple(counts));
}

BOOST_PYTHON_MODULE(parallel_utils)
{
    if (import_mpi4py() < 0) {
        return;
    }

    container_conversions::from_python_sequence<std::vector<int >,
            container_conversions::variable_capacity_policy >();

    to_python_converter<std::vector<int >,
            container_conversions::to_tuple<std::vector<int > > >();

    def("decompose_1d_raw", decompose_1d_raw_wrap);
    def("decompose_1d", decompose_1d_wrap);
    def("decompose_1d_local", decompose_1d_local);

    class_<Logger >("Logger", init<int >())
            .def(init<int, std::string const& >())
            .def(init<std::string const& >())
            .def("write", &Logger::write,
                    return_internal_reference<>())
            .def("flush", &Logger::flush,
                    return_internal_reference<>())
            ;
}

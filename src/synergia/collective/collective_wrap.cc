#include "space_charge_3d_open_hockney.h"
#include "space_charge_2d_bassetti_erskine.h"
#include <boost/python.hpp>
#include "synergia/utils/container_conversions.h"

using namespace boost::python;

BOOST_PYTHON_MODULE(collective)
{
    class_<Space_charge_3d_open_hockney, Space_charge_3d_open_hockney_sptr,
        bases<Collective_operator > >("Space_charge_3d_open_hockney",
                init<Commxx const&, std::vector<int > >())
                .def(init<Commxx const&, std::vector<int >, bool >())
//    class_<Space_charge_3d_open_hockney, Space_charge_3d_open_hockney_sptr>
//        ("Space_charge_3d_open_hockney",
//            init<std::vector<int > const &, bool, Commxx const&>())
        .def("apply", &Space_charge_3d_open_hockney::apply)
        ;

    class_<Space_charge_2d_bassetti_erskine, Space_charge_2d_bassetti_erskine_sptr,
        bases<Collective_operator > >("Space_charge_2d_bassetti_erskine",
                init<>())
        .def("apply", &Space_charge_2d_bassetti_erskine::apply)
        ;
}

#include "space_charge_3d_open_hockney.h"
#include "interpolate_rectangular_zyx.h"
#include "space_charge_2d_bassetti_erskine.h"
#include "impedance.h"
#include <boost/python.hpp>
#include "synergia/utils/container_conversions.h"

using namespace boost::python;

BOOST_PYTHON_MODULE(collective)
{
    class_<Space_charge_3d_open_hockney, Space_charge_3d_open_hockney_sptr,
        bases<Collective_operator > >("Space_charge_3d_open_hockney",
                init<Commxx const&, std::vector<int > >())
                .def(init<Commxx const&, std::vector<int >, bool >())
                .def(init<Commxx const&, std::vector<int >, bool, bool >())
                .def(init<Commxx const&, std::vector<int >, bool, bool, double >())
                .def(init<Commxx const&, std::vector<int >, bool, bool,
                        double, bool >())
                .def(init<Commxx const&, std::vector<int >, bool, bool,
                        double, bool >())
                .def(init<Commxx const&, std::vector<int >, bool, bool,
                        double, bool, double >())
                .def("apply", &Space_charge_3d_open_hockney::apply)
        ;

    class_<Space_charge_2d_bassetti_erskine, Space_charge_2d_bassetti_erskine_sptr,
        bases<Collective_operator > >("Space_charge_2d_bassetti_erskine",
                init<>())
        .def("apply", &Space_charge_2d_bassetti_erskine::apply)
        ;

    class_<Impedance,Impedance_sptr,
        bases<Collective_operator > >("Impedance",
                init<std::string const &, double const &,  double const &, int const &, std::string const &, int const>())
        .def("apply", &Impedance::apply);

    def("interpolate_rectangular_zyx", &interpolate_rectangular_zyx);
}

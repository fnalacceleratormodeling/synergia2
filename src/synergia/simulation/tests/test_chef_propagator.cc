#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "bunch_fixture.h"
#include "synergia/lattice/tests/chef_lattice_sptr_fixture.h"
#include "synergia/simulation/chef_propagator.h"
#include "synergia/utils/serialization.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

const double tolerance = 1.0e-12;

const int begin1 = 1;
const int end1 = 2;

BOOST_FIXTURE_TEST_CASE(construct, Chef_lattice_sptr_fixture)
{
    Chef_lattice_section_sptr chef_lattice_section_sptr(
            new Chef_lattice_section(chef_lattice_sptr, begin1, end1));
    Chef_propagator chef_propagator(chef_lattice_section_sptr);
}

BOOST_AUTO_TEST_CASE(apply)
{
    Chef_lattice_sptr_fixture c;
    Bunch_fixture b;

    Chef_lattice_section_sptr chef_lattice_section_sptr(
            new Chef_lattice_section(c.chef_lattice_sptr, begin1, end1));
    Chef_propagator chef_propagator(chef_lattice_section_sptr);

    //    multi_array_print(bunch.get_local_particles(),"particles before");
    Logger logger(0);
    int verbosity = 5;
    chef_propagator.apply(b.bunch, verbosity, logger);
    //    multi_array_print(bunch.get_local_particles(),"particles after");
}
// test_note: the apply test just verifies that the method doesn't crash.
//            At least one quantitative test would be desirable.

BOOST_FIXTURE_TEST_CASE(serialize_xml, Chef_lattice_sptr_fixture)
{
    Chef_lattice_section_sptr chef_lattice_section_sptr(
            new Chef_lattice_section(chef_lattice_sptr, begin1, end1));
    Chef_propagator chef_propagator(chef_lattice_section_sptr);

    xml_save<Chef_propagator >(chef_propagator,"chef_propagator.xml");

    Chef_propagator loaded;
    xml_load<Chef_propagator >(loaded,"chef_propagator.xml");
}

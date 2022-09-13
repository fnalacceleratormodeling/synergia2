#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/simulation/independent_operation.h"
#include "synergia/lattice/tests/chef_lattice_sptr_fixture.h"
#include "synergia/utils/serialization.h"
#include "bunch_fixture.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture);

const double tolerance = 1.0e-12;

const int begin1 = 1;
const int end1 = 2;

BOOST_FIXTURE_TEST_CASE(construct, Chef_lattice_sptr_fixture)
{
    Chef_lattice_section chef_lattice_section(chef_lattice_sptr, begin1, end1);
    Chef_propagate_operation chef_propagate_operation(chef_lattice_section);
}

BOOST_FIXTURE_TEST_CASE(get_type, Chef_lattice_sptr_fixture)
{
    Chef_lattice_section chef_lattice_section(chef_lattice_sptr, begin1, end1);
    Chef_propagate_operation chef_propagate_operation(chef_lattice_section);

    BOOST_CHECK_EQUAL(chef_propagate_operation.get_type(), "chef_propagate");
}

BOOST_FIXTURE_TEST_CASE(apply, Chef_lattice_sptr_fixture)
{
    Bunch_fixture b;

    Chef_lattice_section chef_lattice_section(chef_lattice_sptr, begin1, end1);
    Chef_propagate_operation chef_propagate_operation(chef_lattice_section);

    const int verbosity = 5;
    Logger logger(0);
    chef_propagate_operation.apply(b.bunch, verbosity, logger);
}
// test_note: We need to check that apply actual produces the correct results.

BOOST_FIXTURE_TEST_CASE(serialize_xml, Chef_lattice_sptr_fixture)
{
    Chef_lattice_section chef_lattice_section(chef_lattice_sptr, begin1, end1);
    Chef_propagate_operation chef_propagate_operation(chef_lattice_section);
    xml_save<Chef_propagate_operation >(chef_propagate_operation,
            "chef_propagate_operation.xml");

    Chef_propagate_operation loaded;
    xml_load<Chef_propagate_operation >(loaded, "chef_propagate_operation.xml");
}


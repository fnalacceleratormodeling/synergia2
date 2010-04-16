#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "components/simulation/operation_extractor.h"
#include "lattice_fixture.h"
#include "utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

const double tolerance = 1.0e-12;

BOOST_AUTO_TEST_CASE(construct)
{
    Operation_extractor_map operation_extractor_map;
}

BOOST_FIXTURE_TEST_CASE(set_get, Lattice_fixture)
{
    Chef_lattice_sptr chef_lattice_sptr(new Chef_lattice(*lattice_sptr));
    const int map_order = 3;

    Operation_extractor_sptr chef_map_operation_extractor(
            new Chef_map_operation_extractor(chef_lattice_sptr, map_order));
    Operation_extractor_sptr
            chef_propagate_operation_extractor(
                    new Chef_propagate_operation_extractor(chef_lattice_sptr,
                            map_order));
    Operation_extractor_sptr chef_mixed_operation_extractor(
            new Chef_mixed_operation_extractor(chef_lattice_sptr, map_order));

    Operation_extractor_map operation_extractor_map;
    operation_extractor_map.set_extractor(default_operation_extractor_name,
            chef_mixed_operation_extractor);
    operation_extractor_map.set_extractor(chef_mixed_operation_extractor_name,
            chef_mixed_operation_extractor);
    operation_extractor_map.set_extractor(
            chef_propagate_operation_extractor_name,
            chef_propagate_operation_extractor);
    operation_extractor_map.set_extractor(chef_map_operation_extractor_name,
            chef_map_operation_extractor);

    BOOST_CHECK_EQUAL(operation_extractor_map.get_extractor(
                    default_operation_extractor_name).get(),
            chef_mixed_operation_extractor.get());
    BOOST_CHECK_EQUAL(operation_extractor_map.get_extractor(
                    chef_mixed_operation_extractor_name).get(),
            chef_mixed_operation_extractor.get());
    BOOST_CHECK_EQUAL(operation_extractor_map.get_extractor(
                    chef_map_operation_extractor_name).get(),
            chef_map_operation_extractor.get());
    BOOST_CHECK_EQUAL(operation_extractor_map.get_extractor(
                    chef_propagate_operation_extractor_name).get(),
            chef_propagate_operation_extractor.get());
}

BOOST_FIXTURE_TEST_CASE(get_extractor_names, Lattice_fixture)
{
    Chef_lattice_sptr chef_lattice_sptr(new Chef_lattice(*lattice_sptr));
    const int map_order = 3;

    Operation_extractor_sptr chef_map_operation_extractor(
            new Chef_map_operation_extractor(chef_lattice_sptr, map_order));
    Operation_extractor_sptr
            chef_propagate_operation_extractor(
                    new Chef_propagate_operation_extractor(chef_lattice_sptr,
                            map_order));
    Operation_extractor_sptr chef_mixed_operation_extractor(
            new Chef_mixed_operation_extractor(chef_lattice_sptr, map_order));

    Operation_extractor_map operation_extractor_map;
    operation_extractor_map.set_extractor(default_operation_extractor_name,
            chef_mixed_operation_extractor);
    operation_extractor_map.set_extractor(chef_mixed_operation_extractor_name,
            chef_mixed_operation_extractor);
    operation_extractor_map.set_extractor(
            chef_propagate_operation_extractor_name,
            chef_propagate_operation_extractor);
    operation_extractor_map.set_extractor(chef_map_operation_extractor_name,
            chef_map_operation_extractor);

    std::list<std::string >
            names(operation_extractor_map.get_extractor_names());
    std::list<std::string > expected_names;
    expected_names.push_back(default_operation_extractor_name);
    expected_names.push_back(chef_mixed_operation_extractor_name);
    expected_names.push_back(chef_map_operation_extractor_name);
    expected_names.push_back(chef_propagate_operation_extractor_name);

    BOOST_CHECK_EQUAL(names.size(), expected_names.size());
    names.sort();
    expected_names.sort();
    for (std::list<std::string >::iterator it = names.begin(), expected_it =
            expected_names.begin(); it != names.end(); ++it, ++expected_it) {
        BOOST_CHECK_EQUAL((*it), (*expected_it));
    }
}

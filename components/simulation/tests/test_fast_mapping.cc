#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "components/simulation/fast_mapping.h"
#include "fast_mapping_term_fixture.h"
#include "mapping_fixture.h"
#include "utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)
;

const double tolerance = 1.0e-12;

BOOST_AUTO_TEST_CASE(construct)
{
    Fast_mapping fast_mapping(order);
}

BOOST_FIXTURE_TEST_CASE(construct2, Mapping_fixture)
{
    Fast_mapping fast_mapping(b.reference_particle, mapping, mapping_length);
}

BOOST_FIXTURE_TEST_CASE(set_get_length, Mapping_fixture)
{
    Fast_mapping fast_mapping(order);
    fast_mapping.set_length(mapping_length);
    BOOST_CHECK_CLOSE(fast_mapping.get_length(), mapping_length, tolerance);
}

BOOST_FIXTURE_TEST_CASE(get_length, Mapping_fixture)
{
    Fast_mapping fast_mapping(b.reference_particle, mapping, mapping_length);
    BOOST_CHECK_CLOSE(fast_mapping.get_length(), mapping_length, tolerance);
}

BOOST_FIXTURE_TEST_CASE(add_term, Fast_mapping_term_fixture)
{
    Fast_mapping fast_mapping(order);
    fast_mapping.add_term(0, fast_mapping_term);
}

BOOST_FIXTURE_TEST_CASE(apply, Mapping_fixture)
{
    Fast_mapping fast_mapping(b.reference_particle, mapping, mapping_length);
    //    multi_array_print(bunch.get_local_particles(), "particles before");
    fast_mapping.apply(b.bunch);
    //    multi_array_print(bunch.get_local_particles(), "particles after");
}
// test_note: We need to check that apply actual produces the correct results.
//            As of this writing, it almost certainly doesn't

BOOST_FIXTURE_TEST_CASE(write_read_file, Mapping_fixture)
{
    Fast_mapping fast_mapping(b.reference_particle, mapping, mapping_length);
    fast_mapping.set_length(mapping_length);

    fast_mapping.write_to_file("test_fast_mapping.dat");
    Fast_mapping fast_mapping2("test_fast_mapping.dat");
    fast_mapping2.write_to_file("test_fast_mapping2.dat");
}
// test_note: We still need to verify that the Fast_mapping read
//            is the same as the Fast_mapping written.
//            For now, we can only diff the files test_fast_mapping.dat and
//            test_fast_mapping2.dat

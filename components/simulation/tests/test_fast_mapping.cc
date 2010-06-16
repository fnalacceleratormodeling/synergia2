#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "components/simulation/fast_mapping.h"
#include "fast_mapping_term_fixture.h"
#include "mapping_fixture.h"
#include "utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

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
    fast_mapping.apply(b.bunch);
}

BOOST_FIXTURE_TEST_CASE(apply_detail, Mapping_fixture)
{
    int order = 2;
    Fast_mapping fast_mapping(order);

    Fast_mapping_term term1(1);
    double c1 = 7.2;
    int i10 = 3;
    term1.coeff() = c1;
    term1.index(0) = i10;
    int ip1 = 2;
    fast_mapping.add_term(ip1, term1);

    Fast_mapping_term term2(2);
    double c2 = 0.37;
    int i20 = 0;
    int i21 = 1;
    term2.coeff() = c2;
    term2.index(0) = i20;
    term2.index(1) = i21;
    int ip2 = 5;
    fast_mapping.add_term(ip2, term2);

    double p[] = { 1.1, 2.2, 3.3, 4.4, 5.5, 6.6 };
    for (int i = 0; i < 6; ++i) {
        b.bunch.get_local_particles()[0][i] = p[i];
    }
    fast_mapping.apply(b.bunch);

    BOOST_CHECK_EQUAL(b.bunch.get_local_particles()[0][0], 0.0);
    BOOST_CHECK_EQUAL(b.bunch.get_local_particles()[0][1], 0.0);
    BOOST_CHECK_EQUAL(b.bunch.get_local_particles()[0][3], 0.0);
    BOOST_CHECK_EQUAL(b.bunch.get_local_particles()[0][4], 0.0);
    BOOST_CHECK_CLOSE(b.bunch.get_local_particles()[0][ip1],
            c1*p[i10], tolerance);
    BOOST_CHECK_CLOSE(b.bunch.get_local_particles()[0][ip2],
            c2*p[i20]*p[i21], tolerance);
}

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

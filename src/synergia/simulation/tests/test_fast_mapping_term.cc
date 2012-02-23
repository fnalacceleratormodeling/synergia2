#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/simulation/fast_mapping.h"
#include "synergia/utils/serialization.h"
#include "fast_mapping_term_fixture.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

const double tolerance = 1.0e-14;

BOOST_AUTO_TEST_CASE(construct)
{
    Fast_mapping_term fast_mapping_term(order);
}

BOOST_AUTO_TEST_CASE(get_order)
{
    Fast_mapping_term fast_mapping_term(order);
    BOOST_CHECK_EQUAL(fast_mapping_term.order(), order);
}

BOOST_AUTO_TEST_CASE(get_set_coeff)
{
    Fast_mapping_term fast_mapping_term(order);
    fast_mapping_term.coeff() = coeff;
    BOOST_CHECK_CLOSE(fast_mapping_term.coeff(), coeff, tolerance);
}

BOOST_AUTO_TEST_CASE(get_set_indices)
{
    Fast_mapping_term fast_mapping_term(2);
    int indices[] = { 0, 2, 4 };
    fast_mapping_term.index(0) = indices[0];
    fast_mapping_term.index(1) = indices[1];
    fast_mapping_term.index(2) = indices[2];

    BOOST_CHECK_EQUAL(fast_mapping_term.index(0), indices[0]);
    BOOST_CHECK_EQUAL(fast_mapping_term.index(1), indices[1]);
    BOOST_CHECK_EQUAL(fast_mapping_term.index(2), indices[2]);
}

BOOST_FIXTURE_TEST_CASE(copy_construct, Fast_mapping_term_fixture)
{
    Fast_mapping_term fast_mapping_term2(fast_mapping_term);

    BOOST_CHECK_EQUAL(fast_mapping_term.order(), fast_mapping_term2.order());
    BOOST_CHECK_CLOSE(fast_mapping_term.coeff(), fast_mapping_term2.coeff(),
            tolerance);
    for (int i = 0; i < fast_mapping_term.order() + 1; ++i) {
        BOOST_CHECK_EQUAL(fast_mapping_term.index(i), fast_mapping_term2.index(i));
    }
}

BOOST_FIXTURE_TEST_CASE(write_read_stream, Fast_mapping_term_fixture)
{
    ofstream out_file("test_fast_mapping_term.dat");
    fast_mapping_term.write_to_stream(out_file);
    out_file.close();

    ifstream in_file("test_fast_mapping_term.dat");
    Fast_mapping_term fast_mapping_term2(in_file);
    in_file.close();

    BOOST_CHECK_EQUAL(fast_mapping_term.order(), fast_mapping_term2.order());
    BOOST_CHECK_CLOSE(fast_mapping_term.coeff(), fast_mapping_term2.coeff(),
            tolerance);
    for (int i = 0; i < fast_mapping_term.order() + 1; ++i) {
        BOOST_CHECK_EQUAL(fast_mapping_term.index(i), fast_mapping_term2.index(i));
    }
}

BOOST_FIXTURE_TEST_CASE(serialize_xml, Fast_mapping_term_fixture)
{
    xml_save<Fast_mapping_term > (fast_mapping_term, "fast_mapping_term.xml");
    Fast_mapping_term fast_mapping_term2;
    xml_load<Fast_mapping_term > (fast_mapping_term2, "fast_mapping_term.xml");

    BOOST_CHECK_EQUAL(fast_mapping_term.order(), fast_mapping_term2.order());
    BOOST_CHECK_CLOSE(fast_mapping_term.coeff(), fast_mapping_term2.coeff(),
            tolerance);
    for (int i = 0; i < fast_mapping_term.order() + 1; ++i) {
        BOOST_CHECK_EQUAL(fast_mapping_term.index(i), fast_mapping_term2.index(i));
    }
}


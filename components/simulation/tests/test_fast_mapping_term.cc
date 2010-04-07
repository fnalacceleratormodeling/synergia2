#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "components/simulation/fast_mapping.h"

const int order = 2;
const double coeff = 3.1415;
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
    Fast_mapping_term fast_mapping_term(order);
    int indices[] = { 0, 2, 4 };
    fast_mapping_term.index(0) = indices[0];
    fast_mapping_term.index(1) = indices[1];
    fast_mapping_term.index(2) = indices[2];

    BOOST_CHECK_EQUAL(fast_mapping_term.index(0), indices[0]);
    BOOST_CHECK_EQUAL(fast_mapping_term.index(1), indices[1]);
    BOOST_CHECK_EQUAL(fast_mapping_term.index(2), indices[2]);
}

BOOST_AUTO_TEST_CASE(copy_construct)
{
    Fast_mapping_term fast_mapping_term(order);
    fast_mapping_term.coeff() = coeff;
    int indices[] = { 0, 2, 4 };
    fast_mapping_term.index(0) = indices[0];
    fast_mapping_term.index(1) = indices[1];
    fast_mapping_term.index(2) = indices[2];

    Fast_mapping_term fast_mapping_term2(fast_mapping_term);

    BOOST_CHECK_EQUAL(fast_mapping_term.order(), fast_mapping_term2.order());
    BOOST_CHECK_CLOSE(fast_mapping_term.coeff(), fast_mapping_term2.coeff(),
            tolerance);
    BOOST_CHECK_EQUAL(fast_mapping_term.index(0), fast_mapping_term2.index(0));
    BOOST_CHECK_EQUAL(fast_mapping_term.index(1), fast_mapping_term2.index(1));
    BOOST_CHECK_EQUAL(fast_mapping_term.index(2), fast_mapping_term2.index(2));
}

BOOST_AUTO_TEST_CASE(write_read_stream)
{
    Fast_mapping_term fast_mapping_term(order);
    fast_mapping_term.coeff() = coeff;
    int indices[] = { 0, 2, 4 };
    fast_mapping_term.index(0) = indices[0];
    fast_mapping_term.index(1) = indices[1];
    fast_mapping_term.index(2) = indices[2];

    ofstream out_file("test_fast_mapping_term.dat");
    fast_mapping_term.write_to_stream(out_file);
    out_file.close();

    ifstream in_file("test_fast_mapping_term.dat");
    Fast_mapping_term fast_mapping_term2(in_file);
    in_file.close();

    BOOST_CHECK_EQUAL(fast_mapping_term.order(), fast_mapping_term2.order());
    BOOST_CHECK_CLOSE(fast_mapping_term.coeff(), fast_mapping_term2.coeff(),
            tolerance);
    BOOST_CHECK_EQUAL(fast_mapping_term.index(0), fast_mapping_term2.index(0));
    BOOST_CHECK_EQUAL(fast_mapping_term.index(1), fast_mapping_term2.index(1));
    BOOST_CHECK_EQUAL(fast_mapping_term.index(2), fast_mapping_term2.index(2));
}


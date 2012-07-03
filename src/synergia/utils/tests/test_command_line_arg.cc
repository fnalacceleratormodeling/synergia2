#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/utils/command_line_arg.h"

BOOST_AUTO_TEST_CASE(construct)
{
    Command_line_arg arg("foo");
}

BOOST_AUTO_TEST_CASE(get_lhs)
{
    Command_line_arg arg("foo=bar");
    BOOST_CHECK(arg.get_lhs() == "foo");
}

BOOST_AUTO_TEST_CASE(get_rhs)
{
    Command_line_arg arg("foo=bar");
    BOOST_CHECK(arg.get_rhs() == "bar");
}

BOOST_AUTO_TEST_CASE(extract_value_int)
{
    Command_line_arg arg("foo=7");
    BOOST_CHECK_EQUAL(arg.extract_value<int >(), 7);
}

BOOST_AUTO_TEST_CASE(extract_value_double)
{
    Command_line_arg arg("foo=2.718");
    const double tolerance = 1.0e-13;
    BOOST_CHECK_CLOSE(arg.extract_value<double >(), 2.718, tolerance);
}

BOOST_AUTO_TEST_CASE(extract_value_bool)
{
    Command_line_arg arg("foo=1");
    BOOST_CHECK_EQUAL(arg.extract_value<bool >(), true);

    Command_line_arg arg2("foo=0");
    BOOST_CHECK_EQUAL(arg2.extract_value<bool >(), false);
}

BOOST_AUTO_TEST_CASE(extract_value_string)
{
    Command_line_arg arg("foo=bar");
    std::string bar("bar");
    BOOST_CHECK_EQUAL(arg.extract_value<std::string >(), bar);
}

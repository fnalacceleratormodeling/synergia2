#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "components/lattice/lattice_element.h"

const std::string name("foo");
const std::string type("bar");
const std::string attr("baz");
const double val(3.1415);
const double tolerance = 1.0e-12;

BOOST_AUTO_TEST_CASE(construct_lattice_element)
{
    Lattice_element lattice_element(type, name);
}

BOOST_AUTO_TEST_CASE(get_type_lattice_element)
{
    Lattice_element lattice_element(type, name);
    BOOST_CHECK_EQUAL(type, lattice_element.get_type());
}

BOOST_AUTO_TEST_CASE(get_name_lattice_element)
{
    Lattice_element lattice_element(type, name);
    BOOST_CHECK_EQUAL(name, lattice_element.get_name());
}

BOOST_AUTO_TEST_CASE(has_attribute)
{
    Lattice_element lattice_element(name, type);
    BOOST_CHECK(lattice_element.has_double_attribute(attr) == false);
}

BOOST_AUTO_TEST_CASE(set_get_attribute)
{
    Lattice_element lattice_element(name, type);
    lattice_element.set_double_attribute(attr, val);
    BOOST_CHECK(lattice_element.has_double_attribute(attr));
    BOOST_CHECK_CLOSE(lattice_element.get_double_attribute(attr), val, tolerance);
}

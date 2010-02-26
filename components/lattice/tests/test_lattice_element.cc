#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "components/lattice/lattice_element.h"

const std::string name("foo");
const std::string type("bar");
const std::string attr("baz");
const std::string attr2("baz2");
const double dblval(3.1415);
const std::string strval("qux");
const std::string strval2("qux");
const double tolerance = 1.0e-12;

BOOST_AUTO_TEST_CASE(construct)
{
    Lattice_element lattice_element(type, name);
}

BOOST_AUTO_TEST_CASE(construct2)
{
    Set_default_attributes_fn_map map(get_standard_default_attributes_fn_map());
    Lattice_element lattice_element(type, name, map);
}

BOOST_AUTO_TEST_CASE(get_type)
{
    Lattice_element lattice_element(type, name);
    BOOST_CHECK_EQUAL(type, lattice_element.get_type());
}

BOOST_AUTO_TEST_CASE(get_name)
{
    Lattice_element lattice_element(type, name);
    BOOST_CHECK_EQUAL(name, lattice_element.get_name());
}

BOOST_AUTO_TEST_CASE(has_double_attribute)
{
    Lattice_element lattice_element(type, name);
    BOOST_CHECK(lattice_element.has_double_attribute(attr) == false);
}

BOOST_AUTO_TEST_CASE(set_get_double_attribute)
{
    Lattice_element lattice_element(type, name);
    lattice_element.set_double_attribute(attr, dblval);
    BOOST_CHECK(lattice_element.has_double_attribute(attr));
    BOOST_CHECK_CLOSE(lattice_element.get_double_attribute(attr), dblval, tolerance);
}

BOOST_AUTO_TEST_CASE(get_double_attributes)
{
    Lattice_element lattice_element(type, name);
    lattice_element.set_double_attribute(attr, dblval);
    BOOST_CHECK(lattice_element.get_double_attributes().count(attr) == 1);
}

BOOST_AUTO_TEST_CASE(has_string_attribute)
{
    Lattice_element lattice_element(type, name);
    BOOST_CHECK(lattice_element.has_string_attribute(attr) == false);
}

BOOST_AUTO_TEST_CASE(set_get_string_attribute)
{
    Lattice_element lattice_element(type, name);
    lattice_element.set_string_attribute(attr, strval);
    BOOST_CHECK(lattice_element.has_string_attribute(attr));
    BOOST_CHECK(lattice_element.get_string_attribute(attr) == strval);
}

BOOST_AUTO_TEST_CASE(get_string_attributes)
{
    Lattice_element lattice_element(type, name);
    lattice_element.set_string_attribute(attr, strval);
    BOOST_CHECK(lattice_element.get_string_attributes().count(attr) == 1);
}


BOOST_AUTO_TEST_CASE(copy)
{
    Lattice_element lattice_element(type, name);
    lattice_element.set_double_attribute("one", 1.0);
    lattice_element.set_double_attribute("two", 2.0);
    lattice_element.set_double_attribute("three", 3.0);
    lattice_element.set_string_attribute("foo", "foo");
    lattice_element.set_string_attribute("bar", "bar");
    lattice_element.add_ancestor("ma");
    lattice_element.add_ancestor("grandma");

    Lattice_element copied_lattice_element(lattice_element);

    BOOST_CHECK_EQUAL(copied_lattice_element.get_type(), type);
    BOOST_CHECK_EQUAL(copied_lattice_element.get_name(), name);
    BOOST_CHECK(std::equal(lattice_element.get_double_attributes().begin(),
                    lattice_element.get_double_attributes().end(),
                    copied_lattice_element.get_double_attributes().begin()));
    BOOST_CHECK(std::equal(lattice_element.get_string_attributes().begin(),
                    lattice_element.get_string_attributes().end(),
                    copied_lattice_element.get_string_attributes().begin()));
    BOOST_CHECK(std::equal(lattice_element.get_ancestors().begin(),
                    lattice_element.get_ancestors().end(),
                    copied_lattice_element.get_ancestors().begin()));
}

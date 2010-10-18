#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "components/lattice/lattice_element_slice.h"

const std::string type("quadrupole");
const std::string name("myquad");
const double length(3.3);
double tolerance = 1.0e-12;

BOOST_AUTO_TEST_CASE(construct)
{
    Lattice_element lattice_element(type, name);
    lattice_element.set_double_attribute("l", length);
    Lattice_element_slice lattice_element_slice(lattice_element);
}

BOOST_AUTO_TEST_CASE(construct2)
{
    Lattice_element lattice_element(type, name);
    lattice_element.set_double_attribute("l", length);
    double left = 0.0;
    double right = length;
    Lattice_element_slice lattice_element_slice(lattice_element, left, right);
}

BOOST_AUTO_TEST_CASE(construct2_negative_left)
{
    Lattice_element lattice_element(type, name);
    lattice_element.set_double_attribute("l", length);
    double left = -1.0;
    double right = length;
    bool caught = false;
    try {
        Lattice_element_slice lattice_element_slice(lattice_element, left,
                right);
    }
    catch (std::range_error) {
        caught = true;
    }
    BOOST_CHECK(caught);
}

BOOST_AUTO_TEST_CASE(construct2_beyond_right)
{
    Lattice_element lattice_element(type, name);
    lattice_element.set_double_attribute("l", length);
    double left = 0.0;
    double right = length * 1.1;
    bool caught = false;
    try {
        Lattice_element_slice lattice_element_slice(lattice_element, left,
                right);
    }
    catch (std::range_error) {
        caught = true;
    }
    BOOST_CHECK(caught);
}

BOOST_AUTO_TEST_CASE(is_whole_true)
{
    Lattice_element lattice_element(type, name);
    lattice_element.set_double_attribute("l", length);
    double left = 0.0;
    double right = length;
    Lattice_element_slice lattice_element_slice(lattice_element, left, right);
    BOOST_CHECK(lattice_element_slice.is_whole());
}

BOOST_AUTO_TEST_CASE(is_whole_false)
{
    Lattice_element lattice_element(type, name);
    lattice_element.set_double_attribute("l", length);
    double left = 0.0;
    double right = 0.5 * length;
    Lattice_element_slice lattice_element_slice(lattice_element, left, right);
    BOOST_CHECK(!lattice_element_slice.is_whole());
}

BOOST_AUTO_TEST_CASE(has_left_edge_true)
{
    Lattice_element lattice_element(type, name);
    lattice_element.set_double_attribute("l", length);
    double left = 0.0;
    double right = length;
    Lattice_element_slice lattice_element_slice(lattice_element, left, right);
    BOOST_CHECK(lattice_element_slice.has_left_edge());
}

BOOST_AUTO_TEST_CASE(has_left_edge_false)
{
    Lattice_element lattice_element(type, name);
    lattice_element.set_double_attribute("l", length);
    double left = 0.5 * length;
    double right = length;
    Lattice_element_slice lattice_element_slice(lattice_element, left, right);
    BOOST_CHECK(!lattice_element_slice.has_left_edge());
}

BOOST_AUTO_TEST_CASE(has_right_edge_true)
{
    Lattice_element lattice_element(type, name);
    lattice_element.set_double_attribute("l", length);
    double left = 0.0;
    double right = length;
    Lattice_element_slice lattice_element_slice(lattice_element, left, right);
    BOOST_CHECK(lattice_element_slice.has_right_edge());
}

BOOST_AUTO_TEST_CASE(has_right_edge_false)
{
    Lattice_element lattice_element(type, name);
    lattice_element.set_double_attribute("l", length);
    double left = 0.0;
    double right = 0.5 * length;
    Lattice_element_slice lattice_element_slice(lattice_element, left, right);
    BOOST_CHECK(!lattice_element_slice.has_right_edge());
}

BOOST_AUTO_TEST_CASE(get_left_right)
{
    Lattice_element lattice_element(type, name);
    lattice_element.set_double_attribute("l", length);
    const double arbitrary_left(1.1);
    const double arbitrary_right(2.2);
    Lattice_element_slice lattice_element_slice(lattice_element,
            arbitrary_left, arbitrary_right);
    BOOST_CHECK_CLOSE(arbitrary_left, lattice_element_slice.get_left(), tolerance);
    BOOST_CHECK_CLOSE(arbitrary_right, lattice_element_slice.get_right(), tolerance);
}

BOOST_AUTO_TEST_CASE(get_lattice_element)
{
    Lattice_element lattice_element(type, name);
    lattice_element.set_double_attribute("l", length);
    Lattice_element_slice lattice_element_slice(lattice_element);
    BOOST_CHECK(&(lattice_element_slice.get_lattice_element()) == &lattice_element);
}

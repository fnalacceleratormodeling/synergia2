#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

#include "synergia/utils/xml_serialization.h"
#include "synergia/lattice/lattice_element_slice.h"

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

BOOST_AUTO_TEST_CASE(serialize)
{
    Lattice_element lattice_element(type, name);
    lattice_element.set_double_attribute("l", length);
    double left = 0.0;
    double right = 0.5 * length;
    Lattice_element_slice lattice_element_slice(lattice_element, left, right);

    xml_save<Lattice_element_slice > (lattice_element_slice,
            "lattice_element_slice.xml");

    Lattice_element_slice loaded;
    xml_load<Lattice_element_slice > (loaded, "lattice_element_slice.xml");
}

BOOST_AUTO_TEST_CASE(serialize_save_slices)
{
    Lattice_element lattice_element(type, name);
    lattice_element.set_double_attribute("l", length);
    const double old_parameter = 2.718;
    lattice_element.set_double_attribute("parameter", old_parameter);
    double left1 = 0.0;
    double right1 = 0.5 * length;
    Lattice_element_slice_sptr lattice_element_slice1_sptr(
            new Lattice_element_slice(lattice_element, left1, right1));
    double left2 = right1;
    double right2 = length;
    Lattice_element_slice_sptr lattice_element_slice2_sptr(
            new Lattice_element_slice(lattice_element, left2, right2));

    Lattice_element_slices slices;
    slices.push_back(lattice_element_slice1_sptr);
    slices.push_back(lattice_element_slice2_sptr);

    xml_save<Lattice_element_slices > (slices, "lattice_element_slices.xml");

    const double new_parameter = 4.6692;
    slices.front()->get_lattice_element().set_double_attribute("parameter",
            new_parameter);
    for (Lattice_element_slices::const_iterator it = slices.begin(); it
            != slices.end(); ++it) {
        BOOST_CHECK_CLOSE((*it)->get_lattice_element().get_double_attribute("parameter"),
                new_parameter, tolerance);
    }
}

BOOST_AUTO_TEST_CASE(serialize_load_slices)
{
    Lattice_element_slices slices;
    xml_load<Lattice_element_slices >(slices,"lattice_element_slices.xml");

    const double new_parameter = 4.6692;
    slices.front()->get_lattice_element().set_double_attribute("parameter",
            new_parameter);
    for (Lattice_element_slices::const_iterator it = slices.begin(); it
            != slices.end(); ++it) {
        BOOST_CHECK_CLOSE((*it)->get_lattice_element().get_double_attribute("parameter"),
                new_parameter, tolerance);
    }
}

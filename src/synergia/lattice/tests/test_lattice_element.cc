#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/lattice/lattice_element.h"
#include "synergia/utils/serialization.h"
#include "synergia/utils/serialization_files.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

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

BOOST_AUTO_TEST_CASE(get_nonexistent_double_attribute)
{
    Lattice_element lattice_element(type, name);
    BOOST_CHECK(lattice_element.has_double_attribute(attr) == false);
    bool caught = false;
    try {
        lattice_element.get_double_attribute(attr);
    }
    catch (std::runtime_error) {
        caught = true;
    }
    BOOST_CHECK(caught);
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

BOOST_AUTO_TEST_CASE(get_nonexistent_string_attribute)
{
    Lattice_element lattice_element(type, name);
    BOOST_CHECK(lattice_element.has_string_attribute(attr) == false);
    bool caught = false;
    try {
        std::string val = lattice_element.get_string_attribute(attr);
    }
    catch (std::runtime_error) {
        caught = true;
    }
    BOOST_CHECK(caught);
}

BOOST_AUTO_TEST_CASE(get_string_attributes)
{
    Lattice_element lattice_element(type, name);
    lattice_element.set_string_attribute(attr, strval);
    BOOST_CHECK(lattice_element.get_string_attributes().count(attr) == 1);
}

BOOST_AUTO_TEST_CASE(get_needs_internal_derive)
{
    Lattice_element lattice_element(type, name);
    BOOST_CHECK(!lattice_element.get_needs_internal_derive());
}

BOOST_AUTO_TEST_CASE(set_needs_internal_derive)
{
    Lattice_element lattice_element(type, name);
    lattice_element.set_needs_internal_derive(true);
    BOOST_CHECK(lattice_element.get_needs_internal_derive());
    lattice_element.set_needs_internal_derive(false);
    BOOST_CHECK(!lattice_element.get_needs_internal_derive());
}

BOOST_AUTO_TEST_CASE(get_needs_external_derive)
{
    Lattice_element lattice_element(type, name);
    BOOST_CHECK(!lattice_element.get_needs_external_derive());
}

BOOST_AUTO_TEST_CASE(set_needs_external_derive)
{
    Lattice_element lattice_element(type, name);
    lattice_element.set_needs_external_derive(true);
    BOOST_CHECK(lattice_element.get_needs_external_derive());
    lattice_element.set_needs_external_derive(false);
    BOOST_CHECK(!lattice_element.get_needs_external_derive());
}

BOOST_AUTO_TEST_CASE(copy)
{
    Lattice_element lattice_element(type, name);
    lattice_element.set_double_attribute("one", 1.0);
    lattice_element.set_double_attribute("two", 2.0);
    lattice_element.set_double_attribute("three", 3.0);
    lattice_element.set_string_attribute("foo", "foo");
    lattice_element.set_string_attribute("bar", "bar");
    lattice_element.set_needs_external_derive(true);
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
    BOOST_CHECK_EQUAL(lattice_element.get_needs_external_derive(),
            copied_lattice_element.get_needs_external_derive());
}

BOOST_AUTO_TEST_CASE(get_revision)
{
    Lattice_element lattice_element(type, name);
    BOOST_CHECK_EQUAL(lattice_element.get_revision(), 0);
    lattice_element.set_double_attribute("a", 1.0);
    BOOST_CHECK_EQUAL(lattice_element.get_revision(), 1);
    lattice_element.set_double_attribute("b", 2.0);
    BOOST_CHECK_EQUAL(lattice_element.get_revision(), 2);
    lattice_element.set_double_attribute("c", 2.0, false);
    BOOST_CHECK_EQUAL(lattice_element.get_revision(), 2);
}

BOOST_AUTO_TEST_CASE(copy_revision)
{
    Lattice_element lattice_element(type, name);
    lattice_element.set_double_attribute("a", 1.0);
    lattice_element.set_double_attribute("b", 2.0);
    Lattice_element copied_lattice_element(lattice_element);
    BOOST_CHECK_EQUAL(copied_lattice_element.get_revision(), 0);
}

BOOST_AUTO_TEST_CASE(test_serialize)
{
    Lattice_element lattice_element(type, name);
    lattice_element.set_double_attribute(attr, dblval);
    lattice_element.set_string_attribute(attr, strval);
    xml_save<Lattice_element > (lattice_element, "lattice_element.xml");

    Lattice_element loaded;
    xml_load<Lattice_element > (loaded, "lattice_element.xml");
    BOOST_CHECK(loaded.get_type() == type);
    BOOST_CHECK(loaded.get_name() == name);
    BOOST_CHECK(loaded.has_double_attribute(attr));
    BOOST_CHECK_CLOSE(loaded.get_double_attribute(attr), dblval, tolerance);
    BOOST_CHECK(loaded.has_string_attribute(attr));
    BOOST_CHECK(loaded.get_string_attribute(attr) == strval);
}

#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/lattice/lattice_element.h"
#include "synergia/lattice/lattice.h"
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
const std::string def_attr("corge");
const double def_dblval(4.669201); //first Feigenbaum constant, $\delta$
const std::string def_strval("garply");
const std::string lattice_name("foo_lattice");
const int charge = -1;
const double mass = 100.0;
const double total_energy = 125.0;
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

BOOST_AUTO_TEST_CASE(set_default_element)
{
    Lattice_element lattice_element(type, name);
    Lattice_element_sptr default_sptr(new Lattice_element());
    lattice_element.set_default_element(default_sptr);
}

BOOST_AUTO_TEST_CASE(has_double_attribute)
{
    Lattice_element lattice_element(type, name);
    Lattice_element_sptr default_sptr(new Lattice_element());
    default_sptr->set_double_attribute(def_attr, def_dblval);
    lattice_element.set_default_element(default_sptr);

    BOOST_CHECK(lattice_element.has_double_attribute(attr) == false);
    BOOST_CHECK(lattice_element.has_double_attribute(def_attr) == true);
    BOOST_CHECK(lattice_element.has_double_attribute(def_attr, false) == false);
}

BOOST_AUTO_TEST_CASE(set_get_double_attribute)
{
    Lattice_element lattice_element(type, name);
    Lattice_element_sptr default_sptr(new Lattice_element());
    default_sptr->set_double_attribute(def_attr, def_dblval);
    lattice_element.set_default_element(default_sptr);
    lattice_element.set_double_attribute(attr, dblval);

    BOOST_CHECK(lattice_element.has_double_attribute(attr));
    BOOST_CHECK_CLOSE(lattice_element.get_double_attribute(attr), dblval,
            tolerance);
    BOOST_CHECK(lattice_element.has_double_attribute(def_attr));
    BOOST_CHECK_CLOSE(lattice_element.get_double_attribute(def_attr),
            def_dblval, tolerance);
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
    Lattice_element_sptr default_sptr(new Lattice_element());
    default_sptr->set_string_attribute(def_attr, def_strval);
    lattice_element.set_default_element(default_sptr);

    BOOST_CHECK(lattice_element.has_string_attribute(attr) == false);
    BOOST_CHECK(lattice_element.has_string_attribute(def_attr) == true);
    BOOST_CHECK(lattice_element.has_string_attribute(def_attr, false) == false);
}

BOOST_AUTO_TEST_CASE(set_get_string_attribute)
{
    Lattice_element lattice_element(type, name);
    Lattice_element_sptr default_sptr(new Lattice_element());
    default_sptr->set_string_attribute(def_attr, def_strval);
    lattice_element.set_default_element(default_sptr);
    lattice_element.set_string_attribute(attr, strval);

    BOOST_CHECK(lattice_element.has_string_attribute(attr));
    BOOST_CHECK(lattice_element.get_string_attribute(attr) == strval);
    BOOST_CHECK(lattice_element.has_string_attribute(def_attr));
    BOOST_CHECK(lattice_element.get_string_attribute(def_attr) == def_strval);
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

BOOST_AUTO_TEST_CASE(has_vector_attribute)
{
    std::vector<double > def_vecval(3);
    def_vecval.at(0) = 1.1;
    def_vecval.at(1) = 2.2;
    def_vecval.at(2) = 3.3;

    Lattice_element lattice_element(type, name);
    Lattice_element_sptr default_sptr(new Lattice_element());
    default_sptr->set_vector_attribute(def_attr, def_vecval);
    lattice_element.set_default_element(default_sptr);

    BOOST_CHECK(lattice_element.has_vector_attribute(attr) == false);
    BOOST_CHECK(lattice_element.has_vector_attribute(def_attr) == true);
    BOOST_CHECK(lattice_element.has_vector_attribute(def_attr, false) == false);
}

BOOST_AUTO_TEST_CASE(set_get_vector_attribute)
{
    std::vector<double > vecval(3);
    vecval.at(0) = 1.1;
    vecval.at(1) = 2.2;
    vecval.at(2) = 3.3;
    std::vector<double > def_vecval(3);
    def_vecval.at(0) = 1.1;
    def_vecval.at(1) = 2.2;
    def_vecval.at(2) = 3.3;

    Lattice_element lattice_element(type, name);
    Lattice_element_sptr default_sptr(new Lattice_element());
    default_sptr->set_vector_attribute(def_attr, def_vecval);
    lattice_element.set_default_element(default_sptr);
    lattice_element.set_vector_attribute(attr, vecval);

    BOOST_CHECK(lattice_element.has_vector_attribute(attr));
    BOOST_CHECK(lattice_element.get_vector_attribute(attr) == vecval);
    BOOST_CHECK(lattice_element.has_vector_attribute(def_attr));
    BOOST_CHECK(lattice_element.get_vector_attribute(def_attr) == def_vecval);
}

BOOST_AUTO_TEST_CASE(get_nonexistent_vector_attribute)
{
    Lattice_element lattice_element(type, name);
    BOOST_CHECK(lattice_element.has_vector_attribute(attr) == false);
    bool caught = false;
    try {
        std::vector<double > val = lattice_element.get_vector_attribute(attr);
    }
    catch (std::runtime_error) {
        caught = true;
    }
    BOOST_CHECK(caught);
}

BOOST_AUTO_TEST_CASE(get_vector_attributes)
{
    std::vector<double > vectorval(3);
    vectorval.at(0) = 1.1;
    vectorval.at(1) = 2.2;
    vectorval.at(2) = 3.3;

    Lattice_element lattice_element(type, name);
    lattice_element.set_vector_attribute(attr, vectorval);
    BOOST_CHECK(lattice_element.get_vector_attributes().count(attr) == 1);
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

BOOST_AUTO_TEST_CASE(has_lattice)
{
    Lattice_element lattice_element(type, name);
    BOOST_CHECK(!lattice_element.has_lattice());
}

BOOST_AUTO_TEST_CASE(set_lattice)
{
    Lattice_element lattice_element(type, name);
    Lattice lattice(lattice_name);
    lattice_element.set_lattice(lattice);
    BOOST_CHECK(lattice_element.has_lattice());
}

BOOST_AUTO_TEST_CASE(get_lattice)
{
    Lattice_element lattice_element(type, name);
    Lattice lattice(lattice_name);
    lattice_element.set_lattice(lattice);
    BOOST_CHECK_EQUAL(&(lattice_element.get_lattice()), &lattice);
}

BOOST_AUTO_TEST_CASE(get_lattice_no_lattice)
{
    Lattice_element lattice_element(type, name);

    bool caught = false;
    try {
        lattice_element.get_lattice();
    }
    catch (std::runtime_error &) {
        caught = true;
    }
    BOOST_CHECK(caught);
}

BOOST_AUTO_TEST_CASE(get_lattice_reference_particle)
{
    Lattice_element lattice_element(type, name);
    Lattice lattice(lattice_name);
    Reference_particle reference_particle(charge, mass, total_energy);
    lattice.set_reference_particle(reference_particle);
    lattice.append(lattice_element);

    Lattice_elements const &lattice_elements(lattice.get_elements());
    for (Lattice_elements::const_iterator leit=lattice_elements.begin();
         leit!=lattice_elements.end(); ++leit) {

        BOOST_CHECK_EQUAL((*leit)->get_lattice().get_reference_particle().get_charge(), charge);
        BOOST_CHECK_CLOSE((*leit)->get_lattice().get_reference_particle().get_four_momentum().get_mass(), mass, tolerance);
        BOOST_CHECK_CLOSE((*leit)->get_lattice().get_reference_particle().get_total_energy(), total_energy, tolerance);
    }
}

BOOST_AUTO_TEST_CASE(copy_test)
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
    std::vector<double > vectorval(3);
    vectorval.at(0) = 1.1;
    vectorval.at(1) = 2.2;
    vectorval.at(2) = 3.3;
    lattice_element.set_vector_attribute(attr, vectorval);

    Lattice_element copied_lattice_element(lattice_element);

    BOOST_CHECK_EQUAL(copied_lattice_element.get_type(), type);
    BOOST_CHECK_EQUAL(copied_lattice_element.get_name(), name);
    BOOST_CHECK(
            std::equal(lattice_element.get_double_attributes().begin(), lattice_element.get_double_attributes().end(), copied_lattice_element.get_double_attributes().begin()));
    BOOST_CHECK(
            std::equal(lattice_element.get_string_attributes().begin(), lattice_element.get_string_attributes().end(), copied_lattice_element.get_string_attributes().begin()));
    BOOST_CHECK(
            std::equal(lattice_element.get_ancestors().begin(), lattice_element.get_ancestors().end(), copied_lattice_element.get_ancestors().begin()));
    BOOST_CHECK_EQUAL(lattice_element.get_needs_external_derive(),
            copied_lattice_element.get_needs_external_derive());
    for (int i = 0; i < 3; ++i) {
        BOOST_CHECK_CLOSE(lattice_element.get_vector_attribute(attr).at(i),
                          copied_lattice_element.get_vector_attribute(attr).at(i), tolerance);
    }
}

BOOST_AUTO_TEST_CASE(lsexpr_test1)
{
    Lattice_element lattice_element(type, name);
    lattice_element.set_double_attribute("one", 1.0);

    Lsexpr lsexpr(lattice_element.as_lsexpr());
    Lattice_element copied_lattice_element(lsexpr);

    BOOST_CHECK_EQUAL(copied_lattice_element.get_type(), type);
    BOOST_CHECK_EQUAL(copied_lattice_element.get_name(), name);
    BOOST_CHECK(
                std::equal(lattice_element.get_double_attributes().begin(), lattice_element.get_double_attributes().end(), copied_lattice_element.get_double_attributes().begin()));
    BOOST_CHECK(
                std::equal(lattice_element.get_string_attributes().begin(), lattice_element.get_string_attributes().end(), copied_lattice_element.get_string_attributes().begin()));
    BOOST_CHECK(
                std::equal(lattice_element.get_ancestors().begin(), lattice_element.get_ancestors().end(), copied_lattice_element.get_ancestors().begin()));
    BOOST_CHECK(
            std::equal(lattice_element.get_vector_attributes().begin(), lattice_element.get_vector_attributes().end(), copied_lattice_element.get_vector_attributes().begin()));
}

BOOST_AUTO_TEST_CASE(lsexpr_test2)
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
    std::vector<double > vectorval(3);
    vectorval.at(0) = 1.1;
    vectorval.at(1) = 2.2;
    vectorval.at(2) = 3.3;
    lattice_element.set_vector_attribute(attr, vectorval);

    Lsexpr lsexpr(lattice_element.as_lsexpr());
    Lattice_element copied_lattice_element(lsexpr);

    BOOST_CHECK_EQUAL(copied_lattice_element.get_type(), type);
    BOOST_CHECK_EQUAL(copied_lattice_element.get_name(), name);
    BOOST_CHECK(
                std::equal(lattice_element.get_double_attributes().begin(), lattice_element.get_double_attributes().end(), copied_lattice_element.get_double_attributes().begin()));
    BOOST_CHECK(
                std::equal(lattice_element.get_string_attributes().begin(), lattice_element.get_string_attributes().end(), copied_lattice_element.get_string_attributes().begin()));
    BOOST_CHECK(
                std::equal(lattice_element.get_ancestors().begin(), lattice_element.get_ancestors().end(), copied_lattice_element.get_ancestors().begin()));
    BOOST_CHECK(
                std::equal(lattice_element.get_vector_attributes().begin(), lattice_element.get_vector_attributes().end(), copied_lattice_element.get_vector_attributes().begin()));
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
    Lattice lattice(lattice_name);
    lattice_element.set_lattice(lattice);
    xml_save<Lattice_element >(lattice_element, "lattice_element.xml");

    Lattice_element loaded;
    xml_load<Lattice_element >(loaded, "lattice_element.xml");
    BOOST_CHECK(loaded.get_type() == type);
    BOOST_CHECK(loaded.get_name() == name);
    BOOST_CHECK(loaded.has_double_attribute(attr));
    BOOST_CHECK_CLOSE(loaded.get_double_attribute(attr), dblval, tolerance);
    BOOST_CHECK(loaded.has_string_attribute(attr));
    BOOST_CHECK(loaded.get_string_attribute(attr) == strval);
    BOOST_CHECK(loaded.has_lattice());
    BOOST_CHECK(loaded.get_lattice().get_name() == lattice_name);
}

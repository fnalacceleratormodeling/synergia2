#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/simulation/aperture_operation.h"
#include "lattice_fixture.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

const double tolerance = 1.0e-12;

const double pax1 = 1.0;
const double pay1 = 0.0;
const double pax2 = 0.0;
const double pay2 = 1.0;
const double pax3 = -1.0;
const double pay3 = 0.0;
const double pax4 = 0.0;
const double pay4 = -1.0;

BOOST_FIXTURE_TEST_CASE(construct, Lattice_fixture)
{
    Lattice_element_sptr element_sptr(lattice_sptr->get_elements().front());
    element_sptr->set_double_attribute("the_number_of_vertices", 3);
    element_sptr->set_double_attribute("pax1", pax1);
    element_sptr->set_double_attribute("pay1", pay1);
    element_sptr->set_double_attribute("pax2", pax2);
    element_sptr->set_double_attribute("pay2", pay2);
    element_sptr->set_double_attribute("pax3", pax3);
    element_sptr->set_double_attribute("pay3", pay3);
    element_sptr->set_double_attribute("pax4", pax4);
    element_sptr->set_double_attribute("pay4", pay4);
    Lattice_element_slice_sptr slice_sptr(
            new Lattice_element_slice(*element_sptr));
    Polygon_aperture_operation polygon_aperture_operation(slice_sptr);
}

BOOST_FIXTURE_TEST_CASE(construct_bad, Lattice_fixture)
{
    Lattice_element_sptr element_sptr(lattice_sptr->get_elements().front());
    Lattice_element_slice_sptr slice_sptr(
            new Lattice_element_slice(*element_sptr));
    bool caught = false;
    try {
        Polygon_aperture_operation polygon_aperture_operation(slice_sptr);
    }
    catch (std::runtime_error) {
        caught = true;
    }
    BOOST_CHECK(caught);

    caught = false;
    element_sptr->set_double_attribute("the_number_of_vertices", 3);
    try {
        Polygon_aperture_operation polygon_aperture_operation(slice_sptr);
    }
    catch (std::runtime_error) {
        caught = true;
    }
    BOOST_CHECK(caught);

    caught = false;
    element_sptr->set_double_attribute("pax1", pax1);
    try {
        Polygon_aperture_operation polygon_aperture_operation(slice_sptr);
    }
    catch (std::runtime_error) {
        caught = true;
    }
    BOOST_CHECK(caught);

    caught = false;
    element_sptr->set_double_attribute("pay1", pay1);
    try {
        Polygon_aperture_operation polygon_aperture_operation(slice_sptr);
    }
    catch (std::runtime_error) {
        caught = true;
    }
    BOOST_CHECK(caught);

    caught = false;
    element_sptr->set_double_attribute("pax2", pax2);
    element_sptr->set_double_attribute("pay2", pay2);
    element_sptr->set_double_attribute("pax3", pax3);
    element_sptr->set_double_attribute("pay3", pay3);
    element_sptr->set_double_attribute("pax4", pax4);
    element_sptr->set_double_attribute("pay4", pay4);
    try {
        Polygon_aperture_operation polygon_aperture_operation(slice_sptr);
    }
    catch (std::runtime_error) {
        caught = true;
    }
    BOOST_CHECK(!caught);

}

BOOST_FIXTURE_TEST_CASE(apply, Lattice_fixture)
{
    Lattice_element_sptr element_sptr(lattice_sptr->get_elements().front());
    Lattice_element_slice_sptr slice_sptr(
            new Lattice_element_slice(*element_sptr));
    element_sptr->set_double_attribute("the_number_of_vertices", 4);
    element_sptr->set_double_attribute("pax1", pax1);
    element_sptr->set_double_attribute("pay1", pay1);
    element_sptr->set_double_attribute("pax2", pax2);
    element_sptr->set_double_attribute("pay2", pay2);
    element_sptr->set_double_attribute("pax3", pax3);
    element_sptr->set_double_attribute("pay3", pay3);
    element_sptr->set_double_attribute("pax4", pax4);
    element_sptr->set_double_attribute("pay4", pay4);
    Polygon_aperture_operation polygon_aperture_operation(slice_sptr);
    polygon_aperture_operation.apply(b.bunch);
}

BOOST_FIXTURE_TEST_CASE(operatorequals, Lattice_fixture)
{
    Lattice_elements::iterator it(lattice_sptr->get_elements().begin());
    Lattice_element_sptr element1_sptr(*it);
    Lattice_element_slice_sptr slice1_sptr(
            new Lattice_element_slice(*element1_sptr));
    element1_sptr->set_double_attribute("the_number_of_vertices", 4);
    element1_sptr->set_double_attribute("pax1", pax1);
    element1_sptr->set_double_attribute("pay1", pay1);
    element1_sptr->set_double_attribute("pax2", pax2);
    element1_sptr->set_double_attribute("pay2", pay2);
    element1_sptr->set_double_attribute("pax3", pax3);
    element1_sptr->set_double_attribute("pay3", pay3);
    element1_sptr->set_double_attribute("pax4", pax4);
    element1_sptr->set_double_attribute("pay4", pay4);
    Polygon_aperture_operation polygon_aperture_operation1(slice1_sptr);

    ++it;
    Lattice_element_sptr element2_sptr(*it);
    Lattice_element_slice_sptr slice2_sptr(
            new Lattice_element_slice(*element2_sptr));
    element2_sptr->set_double_attribute("the_number_of_vertices", 4);
    element2_sptr->set_double_attribute("pax1", pax1);
    element2_sptr->set_double_attribute("pay1", pay1);
    element2_sptr->set_double_attribute("pax2", pax2);
    element2_sptr->set_double_attribute("pay2", pay2);
    element2_sptr->set_double_attribute("pax3", pax3);
    element2_sptr->set_double_attribute("pay3", pay3);
    element2_sptr->set_double_attribute("pax4", pax4);
    element2_sptr->set_double_attribute("pay4", pay4);
    Polygon_aperture_operation polygon_aperture_operation2(slice2_sptr);
    BOOST_CHECK(polygon_aperture_operation1 == polygon_aperture_operation2);

    ++it;
    Lattice_element_sptr element3_sptr(*it);
    Lattice_element_slice_sptr slice3_sptr(
            new Lattice_element_slice(*element3_sptr));
    element3_sptr->set_double_attribute("the_number_of_vertices", 4);
    element3_sptr->set_double_attribute("pax1", pax1);
    element3_sptr->set_double_attribute("pay1", pay1);
    element3_sptr->set_double_attribute("pax2", pax2);
    element3_sptr->set_double_attribute("pay2", pay2);
    element3_sptr->set_double_attribute("pax3", pax3);
    element3_sptr->set_double_attribute("pay3", pay3);
    element3_sptr->set_double_attribute("pax4", pax4 / 2.0);
    element3_sptr->set_double_attribute("pay4", pay4 / 2.0);
    Polygon_aperture_operation polygon_aperture_operation3(slice3_sptr);
    BOOST_CHECK(!(polygon_aperture_operation1 == polygon_aperture_operation3));
}

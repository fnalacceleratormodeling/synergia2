#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/simulation/aperture_operation.h"
#include "lattice_fixture.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture);

const double tolerance = 1.0e-12;

const double width = 0.002;
const double height = 0.001;

BOOST_FIXTURE_TEST_CASE(construct, Lattice_fixture)
{
    Lattice_element_sptr element_sptr(lattice_sptr->get_elements().front());
    element_sptr->set_double_attribute("rectangular_aperture_width",
            width);
    element_sptr->set_double_attribute("rectangular_aperture_height",
            height);
    Lattice_element_slice_sptr slice_sptr(
            new Lattice_element_slice(element_sptr));
    Rectangular_aperture_operation rectangular_aperture_operation(slice_sptr);
}

BOOST_FIXTURE_TEST_CASE(construct_bad, Lattice_fixture)
{
    Lattice_element_sptr element_sptr(lattice_sptr->get_elements().front());
    Lattice_element_slice_sptr slice_sptr(
            new Lattice_element_slice(element_sptr));
    bool caught = false;
    try {
        Rectangular_aperture_operation rectangular_aperture_operation(
                slice_sptr);
    }
    catch (std::runtime_error) {
        caught = true;
    }
    BOOST_CHECK(caught);

    caught = false;
    element_sptr->set_double_attribute("rectangular_aperture_width",
            width);
    try {
        Rectangular_aperture_operation rectangular_aperture_operation(
                slice_sptr);
    }
    catch (std::runtime_error) {
        caught = true;
    }
    BOOST_CHECK(caught);

    caught = false;
    element_sptr->set_double_attribute("rectangular_aperture_height",
            height);
    try {
        Rectangular_aperture_operation rectangular_aperture_operation(
                slice_sptr);
    }
    catch (std::runtime_error) {
        caught = true;
    }
    BOOST_CHECK(!caught);
}

BOOST_FIXTURE_TEST_CASE(apply, Lattice_fixture)
{
    Lattice_element_sptr element_sptr(lattice_sptr->get_elements().front());
    element_sptr->set_double_attribute("rectangular_aperture_width",
            width);
    element_sptr->set_double_attribute("rectangular_aperture_height",
            height);
    Lattice_element_slice_sptr slice_sptr(
            new Lattice_element_slice(element_sptr));
    Rectangular_aperture_operation rectangular_aperture_operation(slice_sptr);
}

BOOST_FIXTURE_TEST_CASE(operatorequals, Lattice_fixture)
{
    Lattice_elements::iterator it(lattice_sptr->get_elements().begin());
    Lattice_element_sptr element1_sptr(*it);
    element1_sptr->set_double_attribute(
            "rectangular_aperture_width", width);
    element1_sptr->set_double_attribute("rectangular_aperture_height",
            height);
    Lattice_element_slice_sptr slice1_sptr(
            new Lattice_element_slice(element1_sptr));
    Rectangular_aperture_operation
            rectangular_aperture_operation1(slice1_sptr);

    ++it;
    Lattice_element_sptr element2_sptr(*it);
    element2_sptr->set_double_attribute(
            "rectangular_aperture_width", width);
    element2_sptr->set_double_attribute("rectangular_aperture_height",
            height);
    Lattice_element_slice_sptr slice2_sptr(
            new Lattice_element_slice(element2_sptr));
    Rectangular_aperture_operation
            rectangular_aperture_operation2(slice2_sptr);
    BOOST_CHECK(rectangular_aperture_operation1 == rectangular_aperture_operation2);

    ++it;
    Lattice_element_sptr element3_sptr(*it);
    element3_sptr->set_double_attribute(
            "rectangular_aperture_width", 2 * width);
    element3_sptr->set_double_attribute("rectangular_aperture_height",
            height);
    Lattice_element_slice_sptr slice3_sptr(
            new Lattice_element_slice(element3_sptr));
    Rectangular_aperture_operation
            rectangular_aperture_operation3(slice3_sptr);
    BOOST_CHECK(!(rectangular_aperture_operation1 == rectangular_aperture_operation3));

    ++it;
    Lattice_element_sptr element4_sptr(*it);
    element4_sptr->set_double_attribute(
            "rectangular_aperture_width", width);
    element4_sptr->set_double_attribute("rectangular_aperture_height",
            2 * height);
    Lattice_element_slice_sptr slice4_sptr(
            new Lattice_element_slice(element4_sptr));
    Rectangular_aperture_operation
            rectangular_aperture_operation4(slice4_sptr);
    BOOST_CHECK(!(rectangular_aperture_operation1 == rectangular_aperture_operation4));
}


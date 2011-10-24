#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/simulation/aperture_operation.h"
#include "lattice_fixture.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

const double tolerance = 1.0e-12;

const double radius = 0.002;

BOOST_FIXTURE_TEST_CASE(construct, Lattice_fixture)
{
    Lattice_element_sptr element_sptr(lattice_sptr->get_elements().front());
    element_sptr->set_double_attribute("circular_aperture_radius", radius);
    Circular_aperture_operation circular_aperture_operation(*element_sptr);
}

BOOST_FIXTURE_TEST_CASE(apply, Lattice_fixture)
{
    Lattice_element_sptr element_sptr(lattice_sptr->get_elements().front());
    element_sptr->set_double_attribute("circular_aperture_radius", radius);
    Circular_aperture_operation circular_aperture_operation(*element_sptr);
    circular_aperture_operation.apply(b.bunch);
}

BOOST_FIXTURE_TEST_CASE(operatorequals, Lattice_fixture)
{
    Lattice_elements::iterator it(lattice_sptr->get_elements().begin());
    Lattice_element_sptr element1_sptr(*it);
    element1_sptr->set_double_attribute("circular_aperture_radius", radius);
    Circular_aperture_operation circular_aperture_operation1(*element1_sptr);

    ++it;
    Lattice_element_sptr element2_sptr(*it);
    element2_sptr->set_double_attribute("circular_aperture_radius", radius);
    Circular_aperture_operation circular_aperture_operation2(*element2_sptr);
    BOOST_CHECK(circular_aperture_operation1 == circular_aperture_operation2);

    ++it;
    Lattice_element_sptr element3_sptr(*it);
    element3_sptr->set_double_attribute("circular_aperture_radius", 2 * radius);
    Circular_aperture_operation circular_aperture_operation3(*element3_sptr);
    BOOST_CHECK(!(circular_aperture_operation1 == circular_aperture_operation3));
}


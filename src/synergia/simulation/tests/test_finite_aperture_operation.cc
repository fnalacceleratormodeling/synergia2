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
    Finite_aperture_operation finite_aperture_operation(*element_sptr);
}

BOOST_FIXTURE_TEST_CASE(get_deposited_charge, Lattice_fixture)
{
    Lattice_element_sptr element_sptr(lattice_sptr->get_elements().front());
    Finite_aperture_operation finite_aperture_operation(*element_sptr);
    BOOST_CHECK(finite_aperture_operation.get_deposited_charge() == 0.0);
}

BOOST_FIXTURE_TEST_CASE(increment_deposited_charge, Lattice_fixture)
{
    Lattice_element_sptr element_sptr(lattice_sptr->get_elements().front());
    Finite_aperture_operation finite_aperture_operation(*element_sptr);
    const double dummy_deposited_charge = 2.71828;
    finite_aperture_operation.increase_deposited_charge(dummy_deposited_charge);
    BOOST_CHECK_CLOSE(finite_aperture_operation.get_deposited_charge(), dummy_deposited_charge,
            tolerance);
}

BOOST_FIXTURE_TEST_CASE(reset_deposited_charge, Lattice_fixture)
{
    Lattice_element_sptr element_sptr(lattice_sptr->get_elements().front());
    Finite_aperture_operation finite_aperture_operation(*element_sptr);
    const double dummy_deposited_charge = 2.71828;
    finite_aperture_operation.increase_deposited_charge(dummy_deposited_charge);
    finite_aperture_operation.reset_deposited_charge();
    BOOST_CHECK(finite_aperture_operation.get_deposited_charge() == 0.0);
}

BOOST_FIXTURE_TEST_CASE(apply, Lattice_fixture)
{
    Lattice_element_sptr element_sptr(lattice_sptr->get_elements().front());
    Finite_aperture_operation finite_aperture_operation(*element_sptr);
    int orig_local_num = b.bunch.get_local_num();
    b.bunch.get_local_particles()[0][0] *= 1.0/0.0;
    b.bunch.get_local_particles()[1][1] *= -1.0/0.0;
    b.bunch.get_local_particles()[2][2] *= std::sqrt(-1.0);
    b.bunch.get_local_particles()[3][3] *= -std::sqrt(-1.0);
    b.bunch.get_local_particles()[4][4] *= std::log(-1.0);
    b.bunch.get_local_particles()[5][5] *= -std::log(-1.0);
    finite_aperture_operation.apply(b.bunch);
    const int num_lost = 6;
    BOOST_CHECK_EQUAL(b.bunch.get_local_num(), orig_local_num - num_lost);
    BOOST_CHECK_CLOSE(finite_aperture_operation.get_deposited_charge(), num_lost*
            b.bunch.get_real_num()/b.bunch.get_total_num(), tolerance);
}

BOOST_FIXTURE_TEST_CASE(operatorequals, Lattice_fixture)
{
    Lattice_elements::iterator it(lattice_sptr->get_elements().begin());
    Lattice_element_sptr element1_sptr(*it);
    Finite_aperture_operation finite_aperture_operation1(*element1_sptr);

    ++it;
    Lattice_element_sptr element2_sptr(*it);
    Finite_aperture_operation finite_aperture_operation2(*element2_sptr);
    BOOST_CHECK(finite_aperture_operation1 == finite_aperture_operation2);

    ++it;
    Lattice_element_sptr element3_sptr(*it);
    element3_sptr->set_double_attribute("circular_aperture_radius", 2 * radius);
    Circular_aperture_operation circular_aperture_operation(*element3_sptr);
    BOOST_CHECK(!(finite_aperture_operation1 == circular_aperture_operation));
}


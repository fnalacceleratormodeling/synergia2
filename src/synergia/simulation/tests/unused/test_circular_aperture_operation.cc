#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/simulation/aperture_operation.h"
#include "lattice_fixture.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture);

const double tolerance = 1.0e-12;

const double radius = 0.002;

const double xoffset = -0.02;
const double yoffset =  0.006;

BOOST_FIXTURE_TEST_CASE(construct, Lattice_fixture)
{
    Lattice_element_sptr element_sptr(lattice_sptr->get_elements().front());
    element_sptr->set_double_attribute("circular_aperture_radius", radius);
    Lattice_element_slice_sptr slice_sptr(
            new Lattice_element_slice(element_sptr));
    Circular_aperture_operation circular_aperture_operation(slice_sptr);
}

BOOST_FIXTURE_TEST_CASE(apply, Lattice_fixture)
{
    Lattice_element_sptr element_sptr(lattice_sptr->get_elements().front());
    element_sptr->set_double_attribute("circular_aperture_radius", radius);
    Lattice_element_slice_sptr slice_sptr(
            new Lattice_element_slice(element_sptr));
    Circular_aperture_operation circular_aperture_operation(slice_sptr);
    const int verbosity = 5;
    Logger logger(0);
    circular_aperture_operation.apply(b.bunch, verbosity, logger);
}

BOOST_FIXTURE_TEST_CASE(operatorequals, Lattice_fixture)
{
    Lattice_elements::iterator it(lattice_sptr->get_elements().begin());
    Lattice_element_sptr element1_sptr(*it);
    element1_sptr->set_double_attribute("circular_aperture_radius", radius);
    Lattice_element_slice_sptr slice1_sptr(
            new Lattice_element_slice(element1_sptr));
    Circular_aperture_operation circular_aperture_operation1(slice1_sptr);

    ++it;
    Lattice_element_sptr element2_sptr(*it);
    element2_sptr->set_double_attribute("circular_aperture_radius", radius);
    Lattice_element_slice_sptr slice2_sptr(
            new Lattice_element_slice(element2_sptr));
    Circular_aperture_operation circular_aperture_operation2(slice2_sptr);
    BOOST_CHECK(circular_aperture_operation1 == circular_aperture_operation2);

    ++it;
    Lattice_element_sptr element3_sptr(*it);
    element3_sptr->set_double_attribute("circular_aperture_radius", 2 * radius);
    Lattice_element_slice_sptr slice3_sptr(
            new Lattice_element_slice(element3_sptr));
    Circular_aperture_operation circular_aperture_operation3(slice3_sptr);
    BOOST_CHECK(!(circular_aperture_operation1 == circular_aperture_operation3));
}

BOOST_FIXTURE_TEST_CASE(cut, Lattice_fixture)
{
    Lattice_element_sptr element_sptr(lattice_sptr->get_elements().front());
    element_sptr->set_double_attribute("circular_aperture_radius",
            radius);
    element_sptr->set_double_attribute("circular_aperture_radius",
            radius);
    Lattice_element_slice_sptr slice_sptr(
            new Lattice_element_slice(element_sptr));
    Circular_aperture_operation circular_aperture_operation(slice_sptr);

    const double eps = 1.0e-8;

    // these should all be kept
    const double keeper_list[][2] = {
        {0.0, 0.0}, {radius-eps, 0.0}, {0.0, radius-eps}};

    // these should all be cut
    const double loser_list[][2] = {
        {radius+eps, 0.0}, {0.0, radius+eps}, {radius+eps, radius+eps}};

    const int nkeepers = sizeof(keeper_list)/(2*sizeof(double));
    const int nlosers = sizeof(loser_list)/(2*sizeof(double));

    // construct particle list with all reflections
    MArray2d keep_particles(boost::extents[4*nkeepers][7]);
    MArray2d loser_particles(boost::extents[4*nlosers][7]);

    // fill the particle arrays with reflections of the particle lists
    int j = 0;
    for (int i=0; i<nkeepers; ++i) {
        keep_particles[j][0] = keeper_list[i][0];
        keep_particles[j][2] = keeper_list[i][1];
        ++j;
        keep_particles[j][0] = -keeper_list[i][0];
        keep_particles[j][2] = keeper_list[i][1];
        ++j;
        keep_particles[j][0] = keeper_list[i][0];
        keep_particles[j][2] = -keeper_list[i][1];
        ++j;
        keep_particles[j][0] = -keeper_list[i][0];
        keep_particles[j][2] = -keeper_list[i][1];
        ++j;
    }

    j = 0;
    for (int i=0; i<nlosers; ++i) {
        loser_particles[j][0] = loser_list[i][0];
        loser_particles[j][2] = loser_list[i][1];
        ++j;
        loser_particles[j][0] = -loser_list[i][0];
        loser_particles[j][2] = loser_list[i][1];
        ++j;
        loser_particles[j][0] = loser_list[i][0];
        loser_particles[j][2] = -loser_list[i][1];
        ++j;
        loser_particles[j][0] = -loser_list[i][0];
        loser_particles[j][2] = -loser_list[i][1];
        ++j;
    }

    // all the keepers better be kept and all the losers better be not
    for (int i=0; i<4*nkeepers; ++i) {
        std::cout << "egs: should be kept:" << keep_particles[i][0] << ", " << keep_particles[i][2] << std::endl;
        BOOST_CHECK(!circular_aperture_operation(keep_particles, i));
    }
    for (int i=0; i<4*nlosers; ++i) {
        std::cout << "egs: should be lost:" << loser_particles[i][0] << ", " << loser_particles[i][2] << std::endl;
        BOOST_CHECK(circular_aperture_operation(loser_particles, i));
    }
}

BOOST_FIXTURE_TEST_CASE(cut_offset, Lattice_fixture)
{
    Lattice_element_sptr element_sptr(lattice_sptr->get_elements().front());
    element_sptr->set_double_attribute("hoffset",
            xoffset);
    element_sptr->set_double_attribute("voffset",
            yoffset);
    element_sptr->set_double_attribute("circular_aperture_radius",
            radius);
    element_sptr->set_double_attribute("circular_aperture_radius",
            radius);
    Lattice_element_slice_sptr slice_sptr(
            new Lattice_element_slice(element_sptr));
    Circular_aperture_operation circular_aperture_operation(slice_sptr);

    const double eps = 1.0e-8;

    // these should all be kept
    const double keeper_list[][2] = {
        {0.0, 0.0}, {radius-eps, 0.0}, {0.0, radius-eps}};

    // these should all be cut
    const double loser_list[][2] = {
        {radius+eps, 0.0}, {0.0, radius+eps},
        {radius+eps, radius+eps}};

    const int nkeepers = sizeof(keeper_list)/(2*sizeof(double));
    const int nlosers = sizeof(loser_list)/(2*sizeof(double));

    // construct particle list with all reflections
    MArray2d keep_particles(boost::extents[4*nkeepers][7]);
    MArray2d loser_particles(boost::extents[4*nlosers][7]);

    // fill the particle arrays with reflections of the particle lists
    int j = 0;
    for (int i=0; i<nkeepers; ++i) {
        keep_particles[j][0] = keeper_list[i][0]+xoffset;
        keep_particles[j][2] = keeper_list[i][1]+yoffset;
        ++j;
        keep_particles[j][0] = -keeper_list[i][0]+xoffset;
        keep_particles[j][2] = keeper_list[i][1]+yoffset;
        ++j;
        keep_particles[j][0] = keeper_list[i][0]+xoffset;
        keep_particles[j][2] = -keeper_list[i][1]+yoffset;
        ++j;
        keep_particles[j][0] = -keeper_list[i][0]+xoffset;
        keep_particles[j][2] = -keeper_list[i][1]+yoffset;
        ++j;
    }

    j = 0;
    for (int i=0; i<nlosers; ++i) {
        loser_particles[j][0] = loser_list[i][0];
        loser_particles[j][2] = loser_list[i][1];
        ++j;
        loser_particles[j][0] = -loser_list[i][0];
        loser_particles[j][2] = loser_list[i][1];
        ++j;
        loser_particles[j][0] = loser_list[i][0];
        loser_particles[j][2] = -loser_list[i][1];
        ++j;
        loser_particles[j][0] = -loser_list[i][0];
        loser_particles[j][2] = -loser_list[i][1];
        ++j;
    }

    // all the keepers better be kept and all the losers better be not
    for (int i=0; i<4*nkeepers; ++i) {
        std::cout << "egs: should be kept:" << keep_particles[i][0] << ", " << keep_particles[i][2] << std::endl;
        BOOST_CHECK(!circular_aperture_operation(keep_particles, i));
    }
    for (int i=0; i<4*nlosers; ++i) {
        std::cout << "egs: should be lost:" << loser_particles[i][0] << ", " << loser_particles[i][2] << std::endl;
        BOOST_CHECK(circular_aperture_operation(loser_particles, i));
    }
}


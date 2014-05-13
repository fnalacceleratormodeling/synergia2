#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/simulation/aperture_operation.h"
#include "lattice_fixture.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

const double tolerance = 1.0e-12;

const double width = 0.04;
const double height = 0.02;
const double ear_offset = .002;

BOOST_FIXTURE_TEST_CASE(construct, Lattice_fixture)
{
    Lattice_element_sptr element_sptr(lattice_sptr->get_elements().front());
    element_sptr->set_double_attribute("rectangular_aperture_width",
            width);
    element_sptr->set_double_attribute("rectangular_aperture_height",
            height);
    element_sptr->set_double_attribute("rectangular_aperture_ear_offset",
                                       ear_offset);
    Lattice_element_slice_sptr slice_sptr(
            new Lattice_element_slice(element_sptr));
    Rectangular_with_ears_aperture_operation rectangular_with_ears_aperture_operation(slice_sptr);
}

BOOST_FIXTURE_TEST_CASE(construct_bad, Lattice_fixture)
{
    Lattice_element_sptr element_sptr(lattice_sptr->get_elements().front());
    Lattice_element_slice_sptr slice_sptr(
            new Lattice_element_slice(element_sptr));
    bool caught = false;
    try {
        Rectangular_with_ears_aperture_operation rectangular_with_ears_aperture_operation(
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
        Rectangular_with_ears_aperture_operation rectangular_with_ears_aperture_operation(
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
        Rectangular_with_ears_aperture_operation rectangular_with_ears_aperture_operation(
                slice_sptr);
    }
    catch (std::runtime_error) {
        caught = true;
    }
    BOOST_CHECK(caught);

    caught = false;
    element_sptr->set_double_attribute("rectangular_aperture_ear_offset",
            ear_offset);
    try {
        Rectangular_with_ears_aperture_operation rectangular_with_ears_aperture_operation(
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
    element_sptr->set_double_attribute("rectangular_aperture_ear_offset",
            ear_offset);
    Lattice_element_slice_sptr slice_sptr(
            new Lattice_element_slice(element_sptr));
    Rectangular_with_ears_aperture_operation rectangular_with_ears_aperture_operation(slice_sptr);
}

BOOST_FIXTURE_TEST_CASE(operatorequals, Lattice_fixture)
{
    Lattice_elements::iterator it(lattice_sptr->get_elements().begin());
    Lattice_element_sptr element1_sptr(*it);
    element1_sptr->set_double_attribute(
            "rectangular_aperture_width", width);
    element1_sptr->set_double_attribute("rectangular_aperture_height",
            height);
    element1_sptr->set_double_attribute("rectangular_aperture_ear_offset",
            ear_offset);
    Lattice_element_slice_sptr slice1_sptr(
            new Lattice_element_slice(element1_sptr));
    Rectangular_with_ears_aperture_operation
            rectangular_with_ears_aperture_operation1(slice1_sptr);

    ++it;
    Lattice_element_sptr element2_sptr(*it);
    element2_sptr->set_double_attribute(
            "rectangular_aperture_width", width);
    element2_sptr->set_double_attribute("rectangular_aperture_height",
            height);
    element2_sptr->set_double_attribute("rectangular_aperture_ear_offset",
            ear_offset);
    Lattice_element_slice_sptr slice2_sptr(
            new Lattice_element_slice(element2_sptr));
    Rectangular_with_ears_aperture_operation
            rectangular_with_ears_aperture_operation2(slice2_sptr);
    BOOST_CHECK(rectangular_with_ears_aperture_operation1 == rectangular_with_ears_aperture_operation2);

    ++it;
    Lattice_element_sptr element3_sptr(*it);
    element3_sptr->set_double_attribute(
            "rectangular_aperture_width", 2 * width);
    element3_sptr->set_double_attribute("rectangular_aperture_height",
            height);
    element3_sptr->set_double_attribute("rectangular_aperture_ear_offset",
            ear_offset);
    Lattice_element_slice_sptr slice3_sptr(
            new Lattice_element_slice(element3_sptr));
    Rectangular_with_ears_aperture_operation
            rectangular_with_ears_aperture_operation3(slice3_sptr);
    BOOST_CHECK(!(rectangular_with_ears_aperture_operation1 == rectangular_with_ears_aperture_operation3));

    ++it;
    Lattice_element_sptr element4_sptr(*it);
    element4_sptr->set_double_attribute(
            "rectangular_aperture_width", width);
    element4_sptr->set_double_attribute("rectangular_aperture_height",
            2 * height);
    element4_sptr->set_double_attribute("rectangular_aperture_ear_offset",
            2 * ear_offset);
    Lattice_element_slice_sptr slice4_sptr(
            new Lattice_element_slice(element4_sptr));
    Rectangular_with_ears_aperture_operation
            rectangular_with_ears_aperture_operation4(slice4_sptr);
    BOOST_CHECK(!(rectangular_with_ears_aperture_operation1 == rectangular_with_ears_aperture_operation4));
}

BOOST_FIXTURE_TEST_CASE(cut, Lattice_fixture)
{
    Lattice_element_sptr element_sptr(lattice_sptr->get_elements().front());
    element_sptr->set_double_attribute("rectangular_aperture_width",
            width);
    element_sptr->set_double_attribute("rectangular_aperture_height",
            height);
    element_sptr->set_double_attribute("rectangular_aperture_ear_offset",
            ear_offset);
    Lattice_element_slice_sptr slice_sptr(
            new Lattice_element_slice(element_sptr));
    Rectangular_with_ears_aperture_operation rectangular_with_ears_aperture_operation(slice_sptr);

    double radius = height/2.0-ear_offset;
    const double eps = 1.0e-8;

    // these should all be kept
    const double keeper_list[][2] = {
        {0.0, 0.0}, {width/2.0-eps, 0.0}, {0.0, height/2.0-eps}, {width/2.0-eps, height/2.0-eps},
        {width/2.0+radius-eps, 0.0}, {width/2.0+radius-eps, ear_offset-eps},
        {width/2.0+radius*0.707, ear_offset+radius*0.707}, {width/2.0+radius*0.5, ear_offset+radius*0.866},
        {width/2.0+radius*0.866, ear_offset+radius*0.5}};

    // these should all be cut
    const double loser_list[][2] = {
        {width/2.0+radius+eps, 0.0}, {0.0, height/2.0+eps}, {width/2.0+radius+eps, ear_offset},
        {width/2.0+0.71*radius, ear_offset+0.71*radius}, {width/2.0+0.51*radius, ear_offset+0.867*radius},
        {width/2.0+0.867*radius, ear_offset+0.51*radius}};

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
        BOOST_CHECK(!rectangular_with_ears_aperture_operation(keep_particles, i));
    }
    for (int i=0; i<4*nlosers; ++i) {
        std::cout << "egs: should be lost:" << loser_particles[i][0] << ", " << loser_particles[i][2] << std::endl;
        BOOST_CHECK(rectangular_with_ears_aperture_operation(loser_particles, i));
    }
}

#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/simulation/aperture_operation.h"
#include "lattice_fixture.h"
#include "synergia/utils/boost_test_mpi_fixture.h"

#include <sstream>

BOOST_GLOBAL_FIXTURE(MPI_fixture);

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
            new Lattice_element_slice(element_sptr));
    Polygon_aperture_operation polygon_aperture_operation(slice_sptr);
}

BOOST_FIXTURE_TEST_CASE(construct_bad, Lattice_fixture)
{
    Lattice_element_sptr element_sptr(lattice_sptr->get_elements().front());
    Lattice_element_slice_sptr slice_sptr(
            new Lattice_element_slice(element_sptr));
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
            new Lattice_element_slice(element_sptr));
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
    const int verbosity = 5;
    Logger logger(0);
    polygon_aperture_operation.apply(b.bunch, verbosity, logger);
}

BOOST_FIXTURE_TEST_CASE(operatorequals, Lattice_fixture)
{
    Lattice_elements::iterator it(lattice_sptr->get_elements().begin());
    Lattice_element_sptr element1_sptr(*it);
    Lattice_element_slice_sptr slice1_sptr(
            new Lattice_element_slice(element1_sptr));
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
            new Lattice_element_slice(element2_sptr));
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
            new Lattice_element_slice(element3_sptr));
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

    ++it;
    Lattice_element_sptr element4_sptr(*it);
    Lattice_element_slice_sptr slice4_sptr(
            new Lattice_element_slice(element4_sptr));
    element4_sptr->set_double_attribute("the_number_of_vertices", 4);
    element4_sptr->set_double_attribute("pax1", pax1);
    element4_sptr->set_double_attribute("pay1", pay1);
    element4_sptr->set_double_attribute("pax2", pax2);
    element4_sptr->set_double_attribute("pay2", pay2);
    element4_sptr->set_double_attribute("pax3", pax3);
    element4_sptr->set_double_attribute("pay3", pay3);
    element4_sptr->set_double_attribute("pax4", pax4);
    element4_sptr->set_double_attribute("pay4", pay4);
    element4_sptr->set_double_attribute("min_radius2", 0.5);
    Polygon_aperture_operation polygon_aperture_operation4(slice4_sptr);
    BOOST_CHECK(!(polygon_aperture_operation1 == polygon_aperture_operation4));
}


// test aperture of a rotated square.
const double topy = 0.05;  // The maximum y (and x) of the rotated square
// points must be ordered in counter-clockwise direction

const double star_pnts[][2] = {
    {0.0, topy}, {-topy, 0.0}, {0.0, -topy}, {topy, 0.0}, {0.0, topy}};

const double xoffset = -0.051; // for offset aperture tests
const double yoffset = 0.023;

BOOST_FIXTURE_TEST_CASE(cut, Lattice_fixture)
{
    Lattice_element_sptr element_sptr(lattice_sptr->get_elements().front());
    int npnts = sizeof(star_pnts)/(2*sizeof(double));
    // paxn attributes indexed from 1
    for (int i=0; i<npnts; ++i) {
        std::stringstream paxss, payss;
        paxss << "pax";
        payss << "pay";
        paxss << i+1;
        payss << i+1;
        element_sptr->set_double_attribute(paxss.str(), star_pnts[i][0]);
        element_sptr->set_double_attribute(payss.str(), star_pnts[i][1]);
    }
    element_sptr->set_double_attribute("the_number_of_vertices", npnts);

    Lattice_element_slice_sptr slice_sptr(
            new Lattice_element_slice(element_sptr));
    Polygon_aperture_operation polygon_aperture_operation(slice_sptr);

    const double eps = 1.0e-8;

    // these should all be kept
    // The euqation for the line of the upper right edge of the square is
    //  y = topy - x.  points under that edge should be retained by the aperture cut.

    const int nkeepers = 400;
    const int nlosers = 400;
    MArray2d keeper_list(boost::extents[nkeepers][2]);
    MArray2d loser_list(boost::extents[nlosers][2]);

    // x goes from 0 to topy - eps
    double xstep = (topy - eps)/nkeepers;
    double x = 0.0;
    for (int i=0; i<nkeepers; ++i) {
        keeper_list[i][0] = x;
        keeper_list[i][1] = topy - x - eps;
        x += xstep;
    }

    // these should all be cut
    // particles outside the line y = topy - eps should be cut
    x = 0.0;
    for (int i=0; i<nlosers; ++i) {
        loser_list[i][0] = x;
        loser_list[i][1] = topy - x + eps;
        x += xstep;
    }

    // construct particle list with all rotations
    MArray2d keep_particles(boost::extents[4*nkeepers][7]);
    MArray2d loser_particles(boost::extents[4*nlosers][7]);

    // fill the particle arrays with rotateions of the particle lists
    for (int i=0; i<nkeepers; ++i) {
        keep_particles[i][0] = keeper_list[i][0];
        keep_particles[i][2] = keeper_list[i][1];
        // rotate 90 degrees clockwise
        keep_particles[i+nkeepers][0] = keeper_list[i][1];
        keep_particles[i+nkeepers][2] = -keeper_list[i][0];
        // rotate 180 degrees clockwise
        keep_particles[i+2*nkeepers][0] = -keeper_list[i][0];
        keep_particles[i+2*nkeepers][2] = -keeper_list[i][1];
        // totate 270 degrees clockwise
        keep_particles[i+3*nkeepers][0] = -keeper_list[i][1];
        keep_particles[i+3*nkeepers][2] = -keeper_list[i][0];
     }

    for (int i=0; i<nlosers; ++i) {
        loser_particles[i][0] = loser_list[i][0];
        loser_particles[i][2] = loser_list[i][1];
        //rotate 90 degrees clockwise
        loser_particles[i+nlosers][0] = loser_list[i][1];
        loser_particles[i+nlosers][2] = -loser_list[i][0];
        // rotate 180 degrees
        loser_particles[i+2*nlosers][0] = -loser_list[i][0];
        loser_particles[i+2*nlosers][2] = -loser_list[i][1];
        // rotate 270 degrees clockwise
        loser_particles[i+3*nlosers][0] = -loser_list[i][1];
        loser_particles[i+3*nlosers][2] = loser_list[i][0];
    }

    // all the keepers better be kept and all the losers better be not
    for (int i=0; i<4*nkeepers; ++i) {
        std::cout << "egs: should be kept:" << std::setprecision(12) << keep_particles[i][0] << ", " << keep_particles[i][2] << std::endl;
        BOOST_CHECK(!polygon_aperture_operation(keep_particles, i));
    }
    for (int i=0; i<4*nlosers; ++i) {
        std::cout << "egs: should be lost:" << std::setprecision(12) << loser_particles[i][0] << ", " << loser_particles[i][2] << std::endl;
        BOOST_CHECK(polygon_aperture_operation(loser_particles, i));
    }
}

BOOST_FIXTURE_TEST_CASE(cut_offset, Lattice_fixture)
{
    Lattice_element_sptr element_sptr(lattice_sptr->get_elements().front());
    int npnts = sizeof(star_pnts)/(2*sizeof(double));
    // paxn attributes indexed from 1
    for (int i=0; i<npnts; ++i) {
        std::stringstream paxss, payss;
        paxss << "pax";
        payss << "pay";
        paxss << i+1;
        payss << i+1;
        element_sptr->set_double_attribute(paxss.str(), star_pnts[i][0]);
        element_sptr->set_double_attribute(payss.str(), star_pnts[i][1]);
    }
    element_sptr->set_double_attribute("the_number_of_vertices", npnts);
    element_sptr->set_double_attribute("hoffset",
            xoffset);
    element_sptr->set_double_attribute("voffset",
            yoffset);

    Lattice_element_slice_sptr slice_sptr(
            new Lattice_element_slice(element_sptr));
    Polygon_aperture_operation polygon_aperture_operation(slice_sptr);

    const double eps = 1.0e-8;

    // these should all be kept
    // The euqation for the line of the upper right edge of the square is
    //  y = topy - x.  points under that edge should be retained by the aperture cut.

    const int nkeepers = 400;
    const int nlosers = 400;
    MArray2d keeper_list(boost::extents[nkeepers][2]);
    MArray2d loser_list(boost::extents[nlosers][2]);

    // x goes from 0 to topy - eps
    double xstep = (topy - eps)/nkeepers;
    double x = 0.0;
    for (int i=0; i<nkeepers; ++i) {
        keeper_list[i][0] = x;
        keeper_list[i][1] = topy - x - eps;
        x += xstep;
    }

    // these should all be cut
    // particles outside the line y = topy - eps should be cut
    x = 0.0;
    for (int i=0; i<nlosers; ++i) {
        loser_list[i][0] = x;
        loser_list[i][1] = topy - x + eps;
        x += xstep;
    }

    // construct particle list with all rotations
    MArray2d keep_particles(boost::extents[4*nkeepers][7]);
    MArray2d loser_particles(boost::extents[4*nlosers][7]);

    // fill the particle arrays with rotateions of the particle lists
    for (int i=0; i<nkeepers; ++i) {
        keep_particles[i][0] = keeper_list[i][0]+xoffset;
        keep_particles[i][2] = keeper_list[i][1]+yoffset;
        // rotate 90 degrees clockwise
        keep_particles[i+nkeepers][0] = keeper_list[i][1]+xoffset;
        keep_particles[i+nkeepers][2] = -keeper_list[i][0]+yoffset;
        // rotate 180 degrees clockwise
        keep_particles[i+2*nkeepers][0] = -keeper_list[i][0]+xoffset;
        keep_particles[i+2*nkeepers][2] = -keeper_list[i][1]+yoffset;
        // totate 270 degrees clockwise
        keep_particles[i+3*nkeepers][0] = -keeper_list[i][1]+xoffset;
        keep_particles[i+3*nkeepers][2] = -keeper_list[i][0]+yoffset;
     }

    for (int i=0; i<nlosers; ++i) {
        loser_particles[i][0] = loser_list[i][0]+xoffset;
        loser_particles[i][2] = loser_list[i][1]+yoffset;
        //rotate 90 degrees clockwise
        loser_particles[i+nlosers][0] = loser_list[i][1]+xoffset;
        loser_particles[i+nlosers][2] = -loser_list[i][0]+yoffset;
        // rotate 180 degrees
        loser_particles[i+2*nlosers][0] = -loser_list[i][0]+xoffset;
        loser_particles[i+2*nlosers][2] = -loser_list[i][1]+yoffset;
        // rotate 270 degrees clockwise
        loser_particles[i+3*nlosers][0] = -loser_list[i][1]+xoffset;
        loser_particles[i+3*nlosers][2] = loser_list[i][0]+yoffset;
    }

    // all the keepers better be kept and all the losers better be not
    for (int i=0; i<4*nkeepers; ++i) {
        std::cout << "egs: should be kept:" << std::setprecision(12) << keep_particles[i][0] << ", " << keep_particles[i][2] << std::endl;
        BOOST_CHECK(!polygon_aperture_operation(keep_particles, i));
    }
    for (int i=0; i<4*nlosers; ++i) {
        std::cout << "egs: should be lost:" << std::setprecision(12) << loser_particles[i][0] << ", " << loser_particles[i][2] << std::endl;
        BOOST_CHECK(polygon_aperture_operation(loser_particles, i));
    }
}


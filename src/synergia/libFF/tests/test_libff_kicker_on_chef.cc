#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/utils/floating_point.h"
#include "synergia/utils/serialization.h"

#include "synergia/lattice/lattice.h"
#include "synergia/lattice/lattice_element.h"
#include "synergia/lattice/chef_lattice.h"
#include "synergia/foundation/four_momentum.h"
#include "synergia/foundation/reference_particle.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/simulation/lattice_simulator.h"
#include "synergia/simulation/bunch_simulator.h"
#include "synergia/simulation/independent_stepper.h"
#include "synergia/simulation/independent_stepper_elements.h"
#include "synergia/simulation/propagator.h"

#include "synergia/lattice/madx_reader.h"

#include "synergia/bunch/bunch.h"

#include "beamline/beamline.h"
#include "beamline/RefRegVisitor.h"

#include "synergia/utils/boost_test_mpi_fixture.h"


#define ELEMENT_FIXTURE(ele) \
    struct ele ## _fixture : public propagator_fixture \
    { ele ## _fixture() : propagator_fixture("seq_" #ele) { } };


#define ELEMENT_FIXTURE_STEPS(ele, steps) \
    struct ele ## _fixture : public propagator_fixture \
    { ele ## _fixture() : propagator_fixture("seq_" #ele, steps) { } };



BOOST_GLOBAL_FIXTURE(MPI_fixture); // needed to initialize MPI

const double tolerance = 1.0e-8;
const double simple_tolerance = 1.0e-12;


void element_check(MArray2d_ref pff, MArray2d_ref pcf, double tolerance)
{
    BOOST_CHECK_CLOSE(pff[0][0], pcf[0][0], tolerance);
    BOOST_CHECK_CLOSE(pff[0][1], pcf[0][1], tolerance);
    BOOST_CHECK_CLOSE(pff[0][2], pcf[0][2], tolerance);
    BOOST_CHECK_CLOSE(pff[0][3], pcf[0][3], tolerance);
    BOOST_CHECK_CLOSE(pff[0][4], pcf[0][4], tolerance);
    BOOST_CHECK_CLOSE(pff[0][5], pcf[0][5], tolerance);
}


// fixture
struct propagator_fixture
{
    propagator_fixture(std::string const & sname, int steps = 1)
    : commxx(new Commxx())
    , seq_name(sname)
    {
        BOOST_TEST_MESSAGE("setup propagator fixture");

        const int map_order = 1;
        const int num_steps = steps;

        // chef propagator
        lattice_chef = reader.get_lattice_sptr(seq_name, "fodo.madx");
        lattice_chef->set_all_string_attribute("extractor_type", "chef_propagate");

        stepper_chef.reset(new Independent_stepper(lattice_chef, map_order, num_steps));
        propagator_chef = new Propagator(stepper_chef);

        Reference_particle refpart(lattice_chef->get_reference_particle());
        Four_momentum four_momentum(refpart.get_four_momentum());
        double momentum = four_momentum.get_momentum();
        four_momentum.set_momentum(momentum*0.25);
        refpart.set_four_momentum(four_momentum);

        bunch_chef.reset(new Bunch(refpart, 1, 1.0e9, commxx));
        bunch_simulator_chef = new Bunch_simulator(bunch_chef);

        // libff propgator
        lattice_ff = reader.get_lattice_sptr(seq_name, "fodo.madx");
        lattice_ff->set_all_string_attribute("extractor_type", "libff");

        stepper_ff.reset(new Independent_stepper(lattice_ff, map_order, num_steps));
        propagator_ff = new Propagator(stepper_ff);

        bunch_ff.reset(new Bunch(refpart, 1, 1.0e9, commxx));
        bunch_simulator_ff = new Bunch_simulator(bunch_ff);
    }

    virtual ~propagator_fixture()
    {
        BOOST_TEST_MESSAGE("teardown propagator fixture");

        delete propagator_chef;
        delete propagator_ff;

        delete bunch_simulator_chef;
        delete bunch_simulator_ff;
    }

    void propagate_chef()
    { propagator_chef->propagate(*bunch_simulator_chef, 1); }

    void propagate_ff()
    { propagator_ff->propagate(*bunch_simulator_ff, 1); }

    MArray2d_ref p_chef() const
    { return bunch_chef->get_local_particles(); }

    MArray2d_ref p_ff() const
    { return bunch_ff->get_local_particles(); }

    Commxx_sptr commxx;
    std::string seq_name;

    MadX_reader reader;

    Lattice_sptr lattice_chef;
    Lattice_sptr lattice_ff;

    Independent_stepper_sptr stepper_chef;
    Independent_stepper_sptr stepper_ff;

    Propagator * propagator_chef;
    Propagator * propagator_ff;

    Bunch_sptr bunch_chef;
    Bunch_sptr bunch_ff;

    Bunch_simulator * bunch_simulator_chef;
    Bunch_simulator * bunch_simulator_ff;
};

ELEMENT_FIXTURE(hkicker);
ELEMENT_FIXTURE(vkicker);
ELEMENT_FIXTURE(kicker);

ELEMENT_FIXTURE_STEPS(kicker2, 2);

ELEMENT_FIXTURE(long_hkicker);
ELEMENT_FIXTURE(long_vkicker);
ELEMENT_FIXTURE(long_kicker);

ELEMENT_FIXTURE(long_hkicker_simple);
ELEMENT_FIXTURE(long_vkicker_simple);
ELEMENT_FIXTURE(long_kicker_simple);

BOOST_FIXTURE_TEST_CASE( test_hkicker, hkicker_fixture )
{
    MArray2d_ref pcf = p_chef();
    MArray2d_ref pff = p_ff();

    pcf[0][0] = 0.1;
    pcf[0][1] = 0.1;
    pcf[0][2] = 0.1;
    pcf[0][3] = 0.1;
    pcf[0][4] = 0.1;
    pcf[0][5] = 0.1;

    pff[0][0] = 0.1;
    pff[0][1] = 0.1;
    pff[0][2] = 0.1;
    pff[0][3] = 0.1;
    pff[0][4] = 0.1;
    pff[0][5] = 0.1;

    std::cout << std::setprecision(16);
    std::cout << "\nhkicker\n";

    propagate_chef();
    propagate_ff();

    for(int i=0; i<6; ++i)
    {
        std::cout << pcf[0][i] << " <--> " << pff[0][i] << "\n";
    }

    element_check(pff, pcf, tolerance);
    BOOST_CHECK(true);
}


BOOST_FIXTURE_TEST_CASE( test_vkicker, vkicker_fixture )
{
    MArray2d_ref pcf = p_chef();
    MArray2d_ref pff = p_ff();

    pcf[0][0] = 0.1;
    pcf[0][1] = 0.1;
    pcf[0][2] = 0.1;
    pcf[0][3] = 0.1;
    pcf[0][4] = 0.1;
    pcf[0][5] = 0.1;

    pff[0][0] = 0.1;
    pff[0][1] = 0.1;
    pff[0][2] = 0.1;
    pff[0][3] = 0.1;
    pff[0][4] = 0.1;
    pff[0][5] = 0.1;

    std::cout << std::setprecision(16);
    std::cout << "\nvkicker\n";

    propagate_chef();
    propagate_ff();

    for(int i=0; i<6; ++i)
    {
        std::cout << pcf[0][i] << " <--> " << pff[0][i] << "\n";
    }

    element_check(pff, pcf, tolerance);
    BOOST_CHECK(true);
}

BOOST_FIXTURE_TEST_CASE( test_kicker, kicker_fixture )
{
    MArray2d_ref pcf = p_chef();
    MArray2d_ref pff = p_ff();

    pcf[0][0] = 0.1;
    pcf[0][1] = 0.1;
    pcf[0][2] = 0.1;
    pcf[0][3] = 0.1;
    pcf[0][4] = 0.1;
    pcf[0][5] = 0.1;

    pff[0][0] = 0.1;
    pff[0][1] = 0.1;
    pff[0][2] = 0.1;
    pff[0][3] = 0.1;
    pff[0][4] = 0.1;
    pff[0][5] = 0.1;

    std::cout << std::setprecision(16);
    std::cout << "\nkicker\n";

    propagate_chef();
    propagate_ff();

    for(int i=0; i<6; ++i)
    {
        std::cout << pcf[0][i] << " <--> " << pff[0][i] << "\n";
    }

    element_check(pff, pcf, tolerance);
    BOOST_CHECK(true);
}

// the tolerance for the long kicker tests is so loose because CHEF uses
// a drift-single kick-drift algorithm but libff is using a 6th order yoshida algorithm.
const double long_kicker_tolerance = 1.0;

BOOST_FIXTURE_TEST_CASE( test_long_hkicker, long_hkicker_fixture )
{
    MArray2d_ref pcf = p_chef();
    MArray2d_ref pff = p_ff();

    pcf[0][0] = 0.1;
    pcf[0][1] = 0.1;
    pcf[0][2] = 0.1;
    pcf[0][3] = 0.1;
    pcf[0][4] = 0.1;
    pcf[0][5] = 0.1;

    pff[0][0] = 0.1;
    pff[0][1] = 0.1;
    pff[0][2] = 0.1;
    pff[0][3] = 0.1;
    pff[0][4] = 0.1;
    pff[0][5] = 0.1;

    std::cout << std::setprecision(16);
    std::cout << "\nlong_hkicker\n";

    propagate_chef();
    propagate_ff();

    for(int i=0; i<6; ++i)
    {
        std::cout << pcf[0][i] << " <--> " << pff[0][i] << "\n";
    }

    element_check(pff, pcf, long_kicker_tolerance);
    BOOST_CHECK(true);
}

BOOST_FIXTURE_TEST_CASE( test_long_hkicker_simple, long_hkicker_simple_fixture )
{
    MArray2d_ref pcf = p_chef();
    MArray2d_ref pff = p_ff();

    pcf[0][0] = 0.1;
    pcf[0][1] = 0.1;
    pcf[0][2] = 0.1;
    pcf[0][3] = 0.1;
    pcf[0][4] = 0.1;
    pcf[0][5] = 0.1;

    pff[0][0] = 0.1;
    pff[0][1] = 0.1;
    pff[0][2] = 0.1;
    pff[0][3] = 0.1;
    pff[0][4] = 0.1;
    pff[0][5] = 0.1;

    std::cout << std::setprecision(16);
    std::cout << "\nlong_hkicker_simple\n";

    propagate_chef();
    propagate_ff();

    for(int i=0; i<6; ++i)
    {
        std::cout << pcf[0][i] << " <--> " << pff[0][i] << "\n";
    }

    element_check(pff, pcf, simple_tolerance);
    BOOST_CHECK(true);
}


BOOST_FIXTURE_TEST_CASE( test_long_vkicker, long_vkicker_fixture )
{
    MArray2d_ref pcf = p_chef();
    MArray2d_ref pff = p_ff();

    pcf[0][0] = 0.1;
    pcf[0][1] = 0.1;
    pcf[0][2] = 0.1;
    pcf[0][3] = 0.1;
    pcf[0][4] = 0.1;
    pcf[0][5] = 0.1;

    pff[0][0] = 0.1;
    pff[0][1] = 0.1;
    pff[0][2] = 0.1;
    pff[0][3] = 0.1;
    pff[0][4] = 0.1;
    pff[0][5] = 0.1;

    std::cout << std::setprecision(16);
    std::cout << "\nlong_vkicker\n";

    propagate_chef();
    propagate_ff();

    for(int i=0; i<6; ++i)
    {
        std::cout << pcf[0][i] << " <--> " << pff[0][i] << "\n";
    }

    element_check(pff, pcf, long_kicker_tolerance);
    BOOST_CHECK(true);
}

BOOST_FIXTURE_TEST_CASE( test_long_vkicker_simple, long_vkicker_simple_fixture )
{
    MArray2d_ref pcf = p_chef();
    MArray2d_ref pff = p_ff();

    pcf[0][0] = 0.1;
    pcf[0][1] = 0.1;
    pcf[0][2] = 0.1;
    pcf[0][3] = 0.1;
    pcf[0][4] = 0.1;
    pcf[0][5] = 0.1;

    pff[0][0] = 0.1;
    pff[0][1] = 0.1;
    pff[0][2] = 0.1;
    pff[0][3] = 0.1;
    pff[0][4] = 0.1;
    pff[0][5] = 0.1;

    std::cout << std::setprecision(16);
    std::cout << "\nlong_vkicker_simple\n";

    propagate_chef();
    propagate_ff();

    for(int i=0; i<6; ++i)
    {
        std::cout << pcf[0][i] << " <--> " << pff[0][i] << "\n";
    }

    element_check(pff, pcf, simple_tolerance);
    BOOST_CHECK(true);
}

BOOST_FIXTURE_TEST_CASE( test_long_kicker, long_kicker_fixture )
{
    MArray2d_ref pcf = p_chef();
    MArray2d_ref pff = p_ff();

    pcf[0][0] = 0.1;
    pcf[0][1] = 0.1;
    pcf[0][2] = 0.1;
    pcf[0][3] = 0.1;
    pcf[0][4] = 0.1;
    pcf[0][5] = 0.1;

    pff[0][0] = 0.1;
    pff[0][1] = 0.1;
    pff[0][2] = 0.1;
    pff[0][3] = 0.1;
    pff[0][4] = 0.1;
    pff[0][5] = 0.1;

    std::cout << std::setprecision(16);
    std::cout << "\nlong_kicker\n";

    propagate_chef();
    propagate_ff();

    for(int i=0; i<6; ++i)
    {
        std::cout << pcf[0][i] << " <--> " << pff[0][i] << "\n";
    }

    element_check(pff, pcf, long_kicker_tolerance);
    BOOST_CHECK(true);
}

BOOST_FIXTURE_TEST_CASE( test_long_kicker_simple, long_kicker_simple_fixture )
{
    MArray2d_ref pcf = p_chef();
    MArray2d_ref pff = p_ff();

    pcf[0][0] = 0.1;
    pcf[0][1] = 0.1;
    pcf[0][2] = 0.1;
    pcf[0][3] = 0.1;
    pcf[0][4] = 0.1;
    pcf[0][5] = 0.1;

    pff[0][0] = 0.1;
    pff[0][1] = 0.1;
    pff[0][2] = 0.1;
    pff[0][3] = 0.1;
    pff[0][4] = 0.1;
    pff[0][5] = 0.1;

    std::cout << std::setprecision(16);
    std::cout << "\nlong_kicker_simple\n";

    propagate_chef();
    propagate_ff();

    for(int i=0; i<6; ++i)
    {
        std::cout << pcf[0][i] << " <--> " << pff[0][i] << "\n";
    }

    element_check(pff, pcf, simple_tolerance);
    BOOST_CHECK(true);
}


BOOST_FIXTURE_TEST_CASE( test_kicker_steps, kicker2_fixture )
{
    MArray2d_ref pcf = p_chef();
    MArray2d_ref pff = p_ff();

    pcf[0][0] = 0.1;
    pcf[0][1] = 0.1;
    pcf[0][2] = 0.1;
    pcf[0][3] = 0.1;
    pcf[0][4] = 0.1;
    pcf[0][5] = 0.1;

    pff[0][0] = 0.1;
    pff[0][1] = 0.1;
    pff[0][2] = 0.1;
    pff[0][3] = 0.1;
    pff[0][4] = 0.1;
    pff[0][5] = 0.1;

    std::cout << std::setprecision(16);
    std::cout << "\nkicker with multiple steps\n";

    propagate_chef();
    propagate_ff();

    for(int i=0; i<6; ++i)
    {
        std::cout << pcf[0][i] << " <--> " << pff[0][i] << "\n";
    }

    element_check(pff, pcf, tolerance);
    BOOST_CHECK(true);
}



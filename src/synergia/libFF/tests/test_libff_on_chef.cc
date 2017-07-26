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

#include "synergia/libFF/ff_element_map.h"

#define ELEMENT_FIXTURE(ele) \
    struct ele ## _fixture : public propagator_fixture \
    { ele ## _fixture() : propagator_fixture("seq_" #ele) { } };


#define ELEMENT_FIXTURE_STEPS(ele, steps) \
    struct ele ## _fixture : public propagator_fixture \
    { ele ## _fixture() : propagator_fixture("seq_" #ele, steps) { } };



BOOST_GLOBAL_FIXTURE(MPI_fixture); // needed to initialize MPI

const double tolerance = 1.0e-8;

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
        four_momentum.set_momentum(momentum*0.95);
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

ELEMENT_FIXTURE(drift);
ELEMENT_FIXTURE(rbend);
ELEMENT_FIXTURE(sbend);
ELEMENT_FIXTURE(cfsbend);
ELEMENT_FIXTURE(quadrupole);
ELEMENT_FIXTURE(quadrupole2);
#if 0
ELEMENT_FIXTURE(sextupole);
ELEMENT_FIXTURE(sextupole2);
ELEMENT_FIXTURE(octupole);
ELEMENT_FIXTURE(octupole2);
#endif
ELEMENT_FIXTURE(rfc);
ELEMENT_FIXTURE(mp1);
ELEMENT_FIXTURE(mp2);
ELEMENT_FIXTURE(mp3);
ELEMENT_FIXTURE(mp4);
ELEMENT_FIXTURE(constfoc);

ELEMENT_FIXTURE(hkicker);
ELEMENT_FIXTURE(vkicker);
ELEMENT_FIXTURE(kicker);

ELEMENT_FIXTURE_STEPS(kicker2, 2);

BOOST_FIXTURE_TEST_CASE( test_drift, drift_fixture )
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
    std::cout << "\ndrift\n";

    propagate_chef();
    propagate_ff();

    for(int i=0; i<6; ++i)
    {
        //std::cout << pcf[0][i] << " <--> " << pff[0][i] << "\n";
    }

    element_check(pff, pcf, tolerance);
}

BOOST_FIXTURE_TEST_CASE( test_rbend, rbend_fixture )
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

    //FF_rbend::set_yoshida_steps(1);

    std::cout << std::setprecision(16);
    std::cout << "\nrbend\n";

    propagate_chef();
    propagate_ff();

    for(int i=0; i<6; ++i)
    {
        //std::cout << pcf[0][i] << " <--> " << pff[0][i] << "\n";
    }

    element_check(pff, pcf, tolerance);
    BOOST_CHECK(true);
}

BOOST_FIXTURE_TEST_CASE( test_sbend, sbend_fixture )
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

    //FF_rbend::set_yoshida_steps(1);

    std::cout << std::setprecision(16);
    std::cout << "\nsbend\n";

    propagate_chef();
    propagate_ff();

    for(int i=0; i<6; ++i)
    {
        //std::cout << pcf[0][i] << " <--> " << pff[0][i] << "\n";
    }

    element_check(pff, pcf, tolerance);
    BOOST_CHECK(true);
}

BOOST_FIXTURE_TEST_CASE( test_cfsbend, cfsbend_fixture )
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

    //FF_rbend::set_yoshida_steps(1);

    std::cout << std::setprecision(16);
    std::cout << "\ncf_sbend\n";

    propagate_chef();
    propagate_ff();

    for(int i=0; i<6; ++i)
    {
        //std::cout << pcf[0][i] << " <--> " << pff[0][i] << "\n";
    }

    element_check(pff, pcf, 3e-5);
    BOOST_CHECK(true);
}

BOOST_FIXTURE_TEST_CASE( test_quadrupole, quadrupole_fixture )
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
    std::cout << "\nquadrupole\n";

    propagate_chef();
    propagate_ff();

    for(int i=0; i<6; ++i)
    {
        //std::cout << pcf[0][i] << " <--> " << pff[0][i] << "\n";
    }

    element_check(pff, pcf, tolerance);
    BOOST_CHECK(true);
}

BOOST_FIXTURE_TEST_CASE( test_quadrupole_with_tilt, quadrupole2_fixture )
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
    std::cout << "\nquadrupole with tilt\n";

    propagate_chef();
    propagate_ff();

    for(int i=0; i<6; ++i)
    {
        //std::cout << pcf[0][i] << " <--> " << pff[0][i] << "\n";
    }

    element_check(pff, pcf, tolerance);
    BOOST_CHECK(true);
}

// tests comparing sextupoles and octupoles against chef have been
// removed because they make no sense.  CHEF only uses one thin kick
// in the middle of the element compared to a yoshida expansion for libff.
// The yoshida expansion agrees with the mad-x by the way.
#if 0
BOOST_FIXTURE_TEST_CASE( test_sextupole, sextupole_fixture )
{
    the_big_giant_global_ff_element_map.get_element_type("sextupole")->set_yoshida_steps(6);
    std::cout << "lattice_chef: " << lattice_chef->as_string() << std::endl;
    std::cout << "lattice_ff: " << lattice_ff->as_string() << std::endl;
    BmlPtr chef_beamline(stepper_chef->get_lattice_simulator().get_chef_lattice().get_beamline_sptr());
    std::cout << "chef_beamline: " << chef_beamline_as_string(chef_beamline) << std::endl;

    MArray2d_ref pcf = p_chef();
    MArray2d_ref pff = p_ff();

    pcf[0][0] = 0.01;
    pcf[0][1] = 0.01;
    pcf[0][2] = 0.01;
    pcf[0][3] = 0.01;
    pcf[0][4] = 0.01;
    pcf[0][5] = 0.01;

    pff[0][0] = 0.01;
    pff[0][1] = 0.01;
    pff[0][2] = 0.01;
    pff[0][3] = 0.01;
    pff[0][4] = 0.01;
    pff[0][5] = 0.01;

    std::cout << std::setprecision(16);
    std::cout << "\nsextupole\n";

    propagate_chef();
    propagate_ff();

    for(int i=0; i<6; ++i)
    {
        //std::cout << pcf[0][i] << " <--> " << pff[0][i] << "\n";
    }

    // tolerance is looser beccause libFF does yoshida but chef doesn't
    //element_check(pff, pcf, 0.001);
    element_check(pff, pcf, tolerance);
    BOOST_CHECK(true);
}

BOOST_FIXTURE_TEST_CASE( test_sextupole_with_tilt, sextupole2_fixture )
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
    std::cout << "\nsextupole with tilt\n";

    propagate_chef();
    propagate_ff();

    for(int i=0; i<6; ++i)
    {
        //std::cout << pcf[0][i] << " <--> " << pff[0][i] << "\n";
    }

    // tolerance is looser because libff does yoshida but chef doesn't
    element_check(pff, pcf, 0.001);
    BOOST_CHECK(true);
}

BOOST_FIXTURE_TEST_CASE( test_octupole, octupole_fixture )
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
    std::cout << "\noctupole\n";

    propagate_chef();
    propagate_ff();

    for(int i=0; i<6; ++i)
    {
        //std::cout << pcf[0][i] << " <--> " << pff[0][i] << "\n";
    }

    // tolerance is looser because libff does yoshida but chef doesn't
    element_check(pff, pcf, 0.001);
    BOOST_CHECK(true);
}

BOOST_FIXTURE_TEST_CASE( test_octupole_with_tilt, octupole2_fixture )
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
    std::cout << "\noctupole with tilt\n";

    propagate_chef();
    propagate_ff();

    for(int i=0; i<6; ++i)
    {
        //std::cout << pcf[0][i] << " <--> " << pff[0][i] << "\n";
    }

    // tolerance is looser because libff does yoshida but chef doesn't
    element_check(pff, pcf, 0.001);
    BOOST_CHECK(true);
}
#endif

BOOST_FIXTURE_TEST_CASE( test_rfc, rfc_fixture )
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
    std::cout << "\nrfcavity\n";

    propagate_chef();
    propagate_ff();

    for(int i=0; i<6; ++i)
    {
        //std::cout << pcf[0][i] << " <--> " << pff[0][i] << "\n";
    }

    element_check(pff, pcf, tolerance);
    BOOST_CHECK(true);
}

BOOST_FIXTURE_TEST_CASE( test_mp1, mp1_fixture )
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
    std::cout << "\nthin multipole 1 -- mad8 format\n";

    propagate_chef();
    propagate_ff();

    for(int i=0; i<6; ++i)
    {
        //std::cout << pcf[0][i] << " <--> " << pff[0][i] << "\n";
    }

    element_check(pff, pcf, tolerance);
    BOOST_CHECK(true);
}

BOOST_FIXTURE_TEST_CASE( test_mp2, mp2_fixture )
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
    std::cout << "\nthin multipole 2 -- mad8 format with tilts\n";

    propagate_chef();
    propagate_ff();

    for(int i=0; i<6; ++i)
    {
        //std::cout << pcf[0][i] << " <--> " << pff[0][i] << "\n";
    }

    element_check(pff, pcf, tolerance);
    BOOST_CHECK(true);
}

BOOST_FIXTURE_TEST_CASE( test_mp3, mp3_fixture )
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
    std::cout << "\nthin multipole 3 -- madx format\n";

    propagate_chef();
    propagate_ff();

    for(int i=0; i<6; ++i)
    {
        //std::cout << pcf[0][i] << " <--> " << pff[0][i] << "\n";
    }

    element_check(pff, pcf, tolerance);
    BOOST_CHECK(true);
}

BOOST_FIXTURE_TEST_CASE( test_mp4, mp4_fixture )
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
    std::cout << "\nthin multipole 4 -- madx format with tilts\n";

    propagate_chef();
    propagate_ff();

    for(int i=0; i<6; ++i)
    {
        //std::cout << pcf[0][i] << " <--> " << pff[0][i] << "\n";
    }

    element_check(pff, pcf, tolerance);
    BOOST_CHECK(true);
}

BOOST_FIXTURE_TEST_CASE( test_constfoc, constfoc_fixture )
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
    std::cout << "\nconstfoc\n";

    propagate_chef();
    propagate_ff();

    for(int i=0; i<6; ++i)
    {
        //std::cout << pcf[0][i] << " <--> " << pff[0][i] << "\n";
    }

    element_check(pff, pcf, tolerance);
    BOOST_CHECK(true);
}


#if 0

BOOST_FIXTURE_TEST_CASE( test_rbend, rbend_fixture )
{

    double x  = (sqrt(1-4.0*sin(0.17/2.0)*sin(0.17/2.0)) - 1.0) / (2.0 * sin(0.17/2.0));
    double px = -2.0*sin(0.17/2.0);

    MArray2d_ref pff = p_ff();

    std::cout << "steps\tex\tt\n";

    for (int i=0; i<1024; ++i)
    {
        //int steps = 1 << i;
        int steps = i;
        FF_rbend::set_yoshida_steps(steps);

        for (int j=0; j<1000; ++j)
        {
            pff[j][0] = 0.0;
            pff[j][1] = 0.0;
            pff[j][2] = 0.0;
            pff[j][3] = 0.0;
            pff[j][4] = 0.0;
            pff[j][5] = 0.0;
        }

        double start = MPI_Wtime();
        propagate_ff();
        double end = MPI_Wtime();

        double dx  = abs(pff[0][0] - x);
        double dpx = abs(pff[0][1] - px);

        double ex  = -log(abs(dx / x));
        double epx = -log(abs(dpx / px));

        std::cout << std::setprecision(16);

        std::cout << i << "\t" << ex << "\t" << end-start << "\n";

#if 0
        std::cout << "steps = " << steps << "\n";

        std::cout << "delta x  = " << dx  << ", E(x)  = " << ex  << "\n";
        std::cout << "delta px = " << dpx << ", E(px) = " << epx << "\n";
#endif
    }
}
#endif

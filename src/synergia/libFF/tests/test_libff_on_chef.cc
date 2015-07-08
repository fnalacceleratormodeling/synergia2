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



BOOST_GLOBAL_FIXTURE(MPI_fixture) // needed to initialize MPI

const double tolerance = 1.0e-12;




// fixture
struct propagator_fixture
{
    propagator_fixture(std::string const & sname)
    : commxx(new Commxx())
    , seq_name(sname)
    {
        BOOST_TEST_MESSAGE("setup propagator fixture");

        // chef propagator
        lattice_chef = reader.get_lattice_sptr(seq_name, "fodo.madx");
        lattice_chef->set_all_string_attribute("extractor_type", "chef_propagate");

        stepper_chef.reset(new Independent_stepper_elements(lattice_chef, 1, 1));
        propagator_chef = new Propagator(stepper_chef);

        bunch_chef.reset(new Bunch(lattice_chef->get_reference_particle(), 1, 1.0e9, commxx));
        bunch_simulator_chef = new Bunch_simulator(bunch_chef);
    
        // libff propgator
        lattice_ff = reader.get_lattice_sptr(seq_name, "fodo.madx");
        lattice_ff->set_all_string_attribute("extractor_type", "libff");

        stepper_ff.reset(new Independent_stepper_elements(lattice_ff, 1, 1));
        propagator_ff = new Propagator(stepper_ff);

        bunch_ff.reset(new Bunch(lattice_ff->get_reference_particle(), 1, 1.0e9, commxx));
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

    Independent_stepper_elements_sptr stepper_chef;
    Independent_stepper_elements_sptr stepper_ff;

    Propagator * propagator_chef;
    Propagator * propagator_ff;

    Bunch_sptr bunch_chef;
    Bunch_sptr bunch_ff;

    Bunch_simulator * bunch_simulator_chef;
    Bunch_simulator * bunch_simulator_ff;
};

ELEMENT_FIXTURE(drift);
ELEMENT_FIXTURE(quadrupole);
ELEMENT_FIXTURE(sextupole);

BOOST_FIXTURE_TEST_CASE( test_drift, drift_fixture )
{
    MArray2d_ref pcf = p_chef();
    MArray2d_ref pff = p_ff();

    pcf[0][0] = 0.1;
    pcf[0][1] = 0.1;
    pcf[0][2] = 0.0;
    pcf[0][3] = 0.0;
    pcf[0][4] = 0.0;
    pcf[0][5] = 0.0;

    pff[0][0] = 0.1;
    pff[0][1] = 0.1;
    pff[0][2] = 0.0;
    pff[0][3] = 0.0;
    pff[0][4] = 0.0;
    pff[0][5] = 0.0;

    propagate_chef();
    propagate_ff();

    std::cout << std::setprecision(10);

    for(int i=0; i<6; ++i)
    {
        std::cout << pcf[0][i] << " <--> " << pff[0][i] << "\n";
    }

    BOOST_CHECK(true);
}

BOOST_FIXTURE_TEST_CASE( test_quadrupole, quadrupole_fixture )
{
    MArray2d_ref pcf = p_chef();
    MArray2d_ref pff = p_ff();

    pcf[0][0] = 0.0;
    pcf[0][1] = 0.1;
    pcf[0][2] = 0.0;
    pcf[0][3] = 0.0;
    pcf[0][4] = 0.0;
    pcf[0][5] = 0.0;

    pff[0][0] = 0.0;
    pff[0][1] = 0.1;
    pff[0][2] = 0.0;
    pff[0][3] = 0.0;
    pff[0][4] = 0.0;
    pff[0][5] = 0.0;

    propagate_chef();
    propagate_ff();

    std::cout << std::setprecision(10);

    for(int i=0; i<6; ++i)
    {
        std::cout << pcf[0][i] << " <--> " << pff[0][i] << "\n";
    }

    BOOST_CHECK(true);
}

BOOST_FIXTURE_TEST_CASE( test_sextupole, sextupole_fixture )
{
    MArray2d_ref pcf = p_chef();
    MArray2d_ref pff = p_ff();

    pcf[0][0] = 0.1;
    pcf[0][1] = 0.1;
    pcf[0][2] = 0.0;
    pcf[0][3] = 0.0;
    pcf[0][4] = 0.0;
    pcf[0][5] = 0.0;

    pff[0][0] = 0.1;
    pff[0][1] = 0.1;
    pff[0][2] = 0.0;
    pff[0][3] = 0.0;
    pff[0][4] = 0.0;
    pff[0][5] = 0.0;

    propagate_chef();
    propagate_ff();

    std::cout << std::setprecision(10);

    for(int i=0; i<6; ++i)
    {
        std::cout << pcf[0][i] << " <--> " << pff[0][i] << "\n";
    }

    BOOST_CHECK(true);
}



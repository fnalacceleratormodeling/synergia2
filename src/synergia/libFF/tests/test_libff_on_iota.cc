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
#include "synergia/libFF/ff_rbend.h"

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
struct iota_fixture
{
    iota_fixture(int steps = 1)
    : commxx(new Commxx())
    , seq_name("iota")
    {
        BOOST_TEST_MESSAGE("setup propagator fixture");

        const int map_order = 1;
        const int num_steps = steps;

        // chef propagator
        lattice_chef = reader.get_lattice_sptr(seq_name, "iota_nll.madx");
        lattice_chef->set_all_string_attribute("extractor_type", "chef_propagate");

        stepper_chef.reset(new Independent_stepper_elements(lattice_chef, map_order, num_steps));
        propagator_chef = new Propagator(stepper_chef);

        Reference_particle refpart(lattice_chef->get_reference_particle());
        Four_momentum four_momentum(refpart.get_four_momentum());
        double momentum = four_momentum.get_momentum();
        refpart.set_four_momentum(four_momentum);

        bunch_chef.reset(new Bunch(refpart, 1, 1.0e9, commxx));
        bunch_simulator_chef = new Bunch_simulator(bunch_chef);

        // libff propgator
        lattice_ff = reader.get_lattice_sptr(seq_name, "iota_nll.madx");
        lattice_ff->set_all_string_attribute("extractor_type", "libff");

        stepper_ff.reset(new Independent_stepper_elements(lattice_ff, map_order, num_steps));

        propagator_ff = new Propagator(stepper_ff);

        bunch_ff.reset(new Bunch(refpart, 1, 1.0e9, commxx));
        bunch_simulator_ff = new Bunch_simulator(bunch_ff);
    }

    virtual ~iota_fixture()
    {
        BOOST_TEST_MESSAGE("teardown propagator fixture");

        delete propagator_chef;
        delete propagator_ff;

        delete bunch_simulator_chef;
        delete bunch_simulator_ff;
    }

    void propagate_chef()
    { propagator_chef->propagate(*bunch_simulator_chef, 1, 1, 2); }

    void propagate_ff()
    { propagator_ff->propagate(*bunch_simulator_ff, 1, 1, 2); }

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

BOOST_FIXTURE_TEST_CASE( test_iota, iota_fixture )
{
    MArray2d_ref pcf = p_chef();
    MArray2d_ref pff = p_ff();

    pcf[0][0] = -0.00171090887340559;
    pcf[0][1] = -0.00154341073914047;
    pcf[0][2] = -0.00407240368042679;
    pcf[0][3] = -0.0027155614225723;
    pcf[0][4] = -0.0532406406114278;
    pcf[0][5] = 0.0;

    pff[0][0] = -0.00171090887340559;
    pff[0][1] = -0.00154341073914047;
    pff[0][2] = -0.00407240368042679;
    pff[0][3] = -0.0027155614225723;
    pff[0][4] = -0.0532406406114278;
    pff[0][5] = 0.0;

    std::cout << std::setprecision(16);
    std::cout << "\niota\n";

    propagate_chef();
    propagate_ff();

    for(int i=0; i<6; ++i)
    {
        std::cout << pcf[0][i] << " <--> " << pff[0][i] << "\n";
    }

    element_check(pff, pcf, tolerance);
    BOOST_CHECK(true);
}

#define BOOST_TEST_MAIN
#include <boost/lexical_cast.hpp>
#include <boost/test/unit_test.hpp>
#include "synergia/simulation/propagator.h"
#include "synergia/simulation/lattice_simulator.h"
#include "synergia/simulation/independent_stepper.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/foundation/math_constants.h"
#include "synergia/lattice/madx_adaptor_map.h"
#include "synergia/lattice/chef_utils.h"
#include "synergia/bunch/bunch.h"
#include "synergia/bunch/diagnostics_basic.h"
#include "synergia/bunch/diagnostics_full2.h"
#include "synergia/bunch/diagnostics_particles.h"
#include "synergia/bunch/diagnostics_bulk_track.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
#include "synergia/simulation/diagnostics_actions.h"
#include "synergia/simulation/lattice_elements_actions.h"
#include "synergia/utils/floating_point.h"

#include <beamline/beamline.h>

#define DBG 1

BOOST_GLOBAL_FIXTURE(MPI_fixture);

const double tolerance = 1.0e-12;

// test whether particles with the same actual momentum but different
// bunch reference momentum get to the same spot in a bending magnet.
BOOST_AUTO_TEST_CASE(mass_spectrometer)
{
    const double angle = mconstants::pi/2;
    const double magnet_length = 0.5;
    const double kenergy = 0.1;
    const double volt = 10.0;   // MV = 0.01 GeV

    Lattice_element analyzer("sbend", "analyzer");
    analyzer.set_double_attribute("l", magnet_length);
    analyzer.set_double_attribute("angle", angle);
    analyzer.set_string_attribute("extractor_type", "chef_propagate");

    Lattice_sptr lattice1_sptr(new Lattice("spectrometer",
                                           Element_adaptor_map_sptr(new MadX_adaptor_map)));
    lattice1_sptr->append(analyzer);

    const double pfrac = 0.95;

    Four_momentum pmom(pconstants::mp);
    pmom.set_kinetic_energy(kenergy);
    Reference_particle refpart1(pconstants::proton_charge, pmom);
    double momentum = pmom.get_momentum();
    pmom.set_momentum(pfrac * momentum);
    Reference_particle refpart2(pconstants::proton_charge, pmom);

    std::cout << "refpart1 momentum: " << refpart1.get_momentum() << std::endl;
    std::cout << "refpart2 momentum: " << refpart2.get_momentum() << std::endl;

    lattice1_sptr->set_reference_particle(refpart1);

    const int macro_particles = 1;
    const double real_particles = 1.0e9;
    Commxx_sptr comm_sptr(new Commxx);
    Bunch_sptr bunch1_sptr( new Bunch(refpart1, macro_particles, real_particles, comm_sptr));
    Bunch_sptr bunch2_sptr( new Bunch(refpart2, macro_particles, real_particles, comm_sptr));

    MArray2d_ref local_particles1(bunch1_sptr->get_local_particles());
    MArray2d_ref local_particles2(bunch2_sptr->get_local_particles());

    // bunch1 will have a particle with nominal reference momentum but
    // dpop = pfrac-1

    // bunch2 will have a particle with reference momentum 95% of bunch1's
    // momentum but dpop=0

    local_particles1[0][Bunch::dpop] = pfrac-1.0;

    const int map_order = 1;
    const int num_steps = 1;

    Bunch_simulator bunch_simulator1(bunch1_sptr);
    Bunch_simulator bunch_simulator2(bunch2_sptr);

    Independent_stepper_sptr stepper1_sptr(new Independent_stepper(lattice1_sptr, map_order, num_steps));

    Propagator propagator1(stepper1_sptr);

    propagator1.propagate(bunch_simulator1, 1, 1, 1);
    propagator1.propagate(bunch_simulator2, 1, 1, 1);

    std::cout << "bunch1 x: " << local_particles1[0][Bunch::x] << std::endl;
    std::cout << "bunch2 x: " << local_particles2[0][Bunch::x] << std::endl;

    BOOST_CHECK(floating_point_equal(local_particles1[0][Bunch::x],
                                     local_particles2[0][Bunch::x], tolerance));

}

BOOST_AUTO_TEST_CASE(accelerate_particles)
{
    const double drift_length = 9.0;
    const double cav_length = 1.0;
    const double kenergy = 0.1;
    const double gvolts = 0.2;

    int p = cout.precision();
    std::cout.precision(5);

    Lattice_element d("drift", "d");
    d.set_double_attribute("l", drift_length);

    Lattice_element cav("rfcavity", "cav");
    cav.set_double_attribute("l", cav_length);
    cav.set_double_attribute("harmon", 1.0);
    cav.set_double_attribute("lag", 0.25);
    cav.set_double_attribute("volt", gvolts*1000.0); //rfcavity takes MV

    Lattice_sptr lattice1_sptr(new Lattice("channel",
                                           Element_adaptor_map_sptr(new MadX_adaptor_map)));

    lattice1_sptr->append(cav);
    lattice1_sptr->append(d);

    for (Lattice_elements::const_iterator it = lattice1_sptr->get_elements().begin();
         it != lattice1_sptr->get_elements().end(); ++it) {
        (*it)->set_string_attribute("extractor_type", "chef_propagate");
    }

    lattice1_sptr->print();

    Four_momentum pmom(pconstants::mp);
    pmom.set_kinetic_energy(kenergy);
    Reference_particle refpart1(pconstants::proton_charge, pmom);
    double momentum = pmom.get_momentum();
    lattice1_sptr->set_reference_particle(refpart1);

    const int macro_particles = 3;
    const double real_particles = 1.0e9;
    Commxx_sptr comm_sptr(new Commxx);
    Bunch_sptr bunch1_sptr( new Bunch(refpart1, macro_particles, real_particles, comm_sptr));

    double initial_energy = bunch1_sptr->get_reference_particle().get_total_energy();

    MArray2d_ref local_particles1(bunch1_sptr->get_local_particles());

    local_particles1[1][Bunch::xp] = 0.05;
    local_particles1[2][Bunch::xp] = 0.1;

    double px1_i = local_particles1[1][Bunch::xp]*momentum;
    double px2_i = local_particles1[2][Bunch::xp]*momentum;

    const int map_order = 1;
    const int num_steps = 1;

    Bunch_simulator bunch_simulator1(bunch1_sptr);

    Independent_stepper_sptr stepper1_sptr(new Independent_stepper(lattice1_sptr, map_order, num_steps));


    stepper1_sptr->force_update_operations_no_collective();

    Propagator propagator1(stepper1_sptr);

    std::cout << "before propagate particles[0][:]:";
    for (int i=0; i<6; ++i) {
        std::cout << " " << local_particles1[0][i];
    }
    std::cout << std::endl;

    propagator1.propagate(bunch_simulator1, 1, 1, 5);

    std::cout << "after propagate particles[0][:]:";
    for (int i=0; i<6; ++i) {
        std::cout << " " << local_particles1[0][i];
    }
    std::cout << std::endl;

    double final_energy = bunch1_sptr->get_reference_particle().get_total_energy();
    double final_momentum = bunch1_sptr->get_reference_particle().get_momentum();

    std::cout << "accelerate_particles: initial energy: " << initial_energy << std::endl;
    std::cout << "accelerate_particles: initial momentum: " << momentum << std::endl;
    std::cout << "accelerate_particles: final energy: " << final_energy << std::endl;
    std::cout << "accelerate_particles: final momentum: " << final_momentum << std::endl;

    BOOST_CHECK(floating_point_equal(initial_energy + gvolts,
                                     final_energy, tolerance));

#if 0
    double px1_f = final_momentum * local_particles1[1][Bunch::xp];
    double px2_f = final_momentum * local_particles1[2][Bunch::xp];
    std::cout << "particles[1][1]: " << local_particles1[1][Bunch::xp] << ", px1_f: " << px1_f << std::endl;
    std::cout << "particles[2][1]: " << local_particles1[2][Bunch::xp] << ", px2_f: " << px2_f << std::endl;
#endif
    std::cout.precision(p);

}


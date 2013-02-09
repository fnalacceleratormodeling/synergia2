#define BOOST_TEST_MAIN
#include <boost/lexical_cast.hpp>
#include <boost/test/unit_test.hpp>
#include "synergia/simulation/propagator.h"
#include "synergia/simulation/lattice_simulator.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/lattice/chef_utils.h"
#include "synergia/bunch/bunch.h"
#include "synergia/bunch/diagnostics_basic.h"
#include "synergia/bunch/diagnostics_full2.h"
#include "synergia/bunch/diagnostics_particles.h"
#include "propagator_fixture.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
#include "synergia/simulation/diagnostics_actions.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

const double tolerance = 1.0e-12;

BOOST_FIXTURE_TEST_CASE(construct, Propagator_fixture)
{
}

BOOST_FIXTURE_TEST_CASE(propagate, Propagator_fixture)
{
    Bunch_sptr bunch_sptr(new Bunch(l.b.bunch));
    Bunch_simulator bunch_simulator(bunch_sptr);

    Diagnostics_sptr first_step_full2_diag_sptr(
            new Diagnostics_full2("first_full2_per_step.h5"));
    Diagnostics_sptr second_step_full2_diag_sptr(
            new Diagnostics_full2("second_full2_per_step.h5"));
    bunch_simulator.add_per_step(first_step_full2_diag_sptr);
    bunch_simulator.add_per_step(second_step_full2_diag_sptr);

    Diagnostics_sptr first_turn_particles_diag_sptr(
            new Diagnostics_particles("first_particles_per_turn.h5"));
    Diagnostics_sptr second_turn_particles_diag_sptr(
            new Diagnostics_particles("second_particles_per_turn.h5"));
    bunch_simulator.add_per_turn(first_turn_particles_diag_sptr);
    bunch_simulator.add_per_turn(second_turn_particles_diag_sptr);

    int num_turns = 4;
    propagator.propagate(bunch_simulator, num_turns);
}

BOOST_FIXTURE_TEST_CASE(propagate_max_turns, Propagator_fixture)
{
    Bunch_sptr bunch_sptr(new Bunch(l.b.bunch));
    Bunch_simulator bunch_simulator(bunch_sptr);

    Diagnostics_sptr first_step_full2_diag_sptr(
            new Diagnostics_full2("first_full2_per_step.h5"));
    Diagnostics_sptr second_step_full2_diag_sptr(
            new Diagnostics_full2("second_full2_per_step.h5"));
    bunch_simulator.add_per_step(first_step_full2_diag_sptr);
    bunch_simulator.add_per_step(second_step_full2_diag_sptr);

    Diagnostics_sptr first_turn_particles_diag_sptr(
            new Diagnostics_particles("first_particles_per_turn.h5"));
    Diagnostics_sptr second_turn_particles_diag_sptr(
            new Diagnostics_particles("second_particles_per_turn.h5"));
    bunch_simulator.add_per_turn(first_turn_particles_diag_sptr);
    bunch_simulator.add_per_turn(second_turn_particles_diag_sptr);

    int num_turns = 4;
    int max_turns = 1;
    propagator.propagate(bunch_simulator, num_turns, max_turns);

    const char second_checkpoint[] = "second_checkpoint";
    propagator.set_checkpoint_dir(second_checkpoint);
    propagator.resume(Propagator::default_checkpoint_dir, false, 0, false, 0);

    propagator.resume(second_checkpoint, true, 2, false, 0);
}

BOOST_FIXTURE_TEST_CASE(serialize_xml, Propagator_fixture)
{
    xml_save(propagator, "propagator.xml");

    Bunch_sptr bunch_sptr(new Bunch(l.b.bunch));
    Bunch_simulator bunch_simulator(bunch_sptr);

    bunch_simulator.add_per_turn(
            Diagnostics_sptr(
                    new Diagnostics_particles("test_propagate_per_turn.h5")));
    bunch_simulator.add_per_step(
            Diagnostics_sptr(
                    new Diagnostics_basic("test_propagate_per_step.h5")));

    int num_turns = 4;
    propagator.propagate(bunch_simulator, num_turns);
    xml_save(propagator, "propagator2.xml");

    Propagator loaded;
    xml_load(propagator, "propagator.xml");
}


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

BOOST_FIXTURE_TEST_CASE(get_stepper_sptr, Propagator_fixture)
{
    BOOST_CHECK_EQUAL(propagator.get_stepper_sptr(), stepper_sptr);
}

BOOST_FIXTURE_TEST_CASE(set_get_checkpoint_period, Propagator_fixture)
{
    BOOST_CHECK_EQUAL(propagator.get_checkpoint_period(),
            Propagator::default_checkpoint_period);
    const int new_period = 77;
    propagator.set_checkpoint_period(new_period);
    BOOST_CHECK_EQUAL(propagator.get_checkpoint_period(), new_period);
}

BOOST_FIXTURE_TEST_CASE(set_get_checkpoint_dir, Propagator_fixture)
{
    BOOST_CHECK_EQUAL(propagator.get_checkpoint_dir(),
            Propagator::default_checkpoint_dir);
    std::string new_dir("garply");
    propagator.set_checkpoint_dir(new_dir);
    BOOST_CHECK_EQUAL(propagator.get_checkpoint_dir(), new_dir);
}

BOOST_FIXTURE_TEST_CASE(set_get_checkpoint_with_xml, Propagator_fixture)
{
    BOOST_CHECK_EQUAL(propagator.get_checkpoint_with_xml(),
            false);
    bool new_with_xml = true;
    propagator.set_checkpoint_with_xml(new_with_xml);
    BOOST_CHECK_EQUAL(propagator.get_checkpoint_with_xml(), new_with_xml);
}

BOOST_FIXTURE_TEST_CASE(set_get_final_checkpoint, Propagator_fixture)
{
    BOOST_CHECK_EQUAL(propagator.get_final_checkpoint(),
            false);
    bool new_final_checkpoint = true;
    propagator.set_final_checkpoint(new_final_checkpoint);
    BOOST_CHECK_EQUAL(propagator.get_final_checkpoint(), new_final_checkpoint);
}

BOOST_FIXTURE_TEST_CASE(set_get_concurrent_io, Propagator_fixture)
{
    BOOST_CHECK_EQUAL(propagator.get_concurrent_io(),
            Propagator::default_concurrent_io);
    const int new_period = 77;
    propagator.set_concurrent_io(new_period);
    BOOST_CHECK_EQUAL(propagator.get_concurrent_io(), new_period);
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

BOOST_FIXTURE_TEST_CASE(propagate_train, Propagator_fixture)
{
    const double bunch_spacing = 1.7;
    Bunch_train_sptr bunch_train_sptr(new Bunch_train(bs.bunches, bunch_spacing));
    Bunch_train_simulator bunch_train_simulator(bunch_train_sptr);

    Diagnostics_sptr step_full2_diag0_sptr(
            new Diagnostics_full2("full2_per_step0.h5"));
    bunch_train_simulator.add_per_step(0, step_full2_diag0_sptr);

    Diagnostics_sptr turn_basic_diag0_sptr(
            new Diagnostics_basic("basic_per_turn0.h5"));
    bunch_train_simulator.add_per_turn(0, turn_basic_diag0_sptr);

    Diagnostics_sptr step_full2_diag1_sptr(
            new Diagnostics_full2("full2_per_step1.h5"));
    bunch_train_simulator.add_per_step(1, step_full2_diag1_sptr);

    Diagnostics_sptr turn_basic_diag1_sptr(
            new Diagnostics_basic("basic_per_turn1.h5"));
    bunch_train_simulator.add_per_turn(1, turn_basic_diag1_sptr);

    Diagnostics_sptr step_full2_diag2_sptr(
            new Diagnostics_full2("full2_per_step2.h5"));
    bunch_train_simulator.add_per_step(2, step_full2_diag2_sptr);

    Diagnostics_sptr turn_basic_diag2_sptr(
            new Diagnostics_basic("basic_per_turn2.h5"));
    bunch_train_simulator.add_per_turn(2, turn_basic_diag2_sptr);


    int num_turns = 4;
    propagator.propagate(bunch_train_simulator, num_turns);
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

BOOST_FIXTURE_TEST_CASE(serialize_train_xml, Propagator_fixture)
{
    xml_save(propagator, "propagator_train.xml");

    const double bunch_spacing = 1.7;
    Bunch_train_sptr bunch_train_sptr(new Bunch_train(bs.bunches, bunch_spacing));
    Bunch_train_simulator bunch_train_simulator(bunch_train_sptr);

    bunch_train_simulator.add_per_turn(0,
            Diagnostics_sptr(
                    new Diagnostics_full2("test_propagate_per_turn0.h5")));
    bunch_train_simulator.add_per_step(0,
            Diagnostics_sptr(
                    new Diagnostics_basic("test_propagate_per_step0.h5")));

    bunch_train_simulator.add_per_turn(1,
            Diagnostics_sptr(
                    new Diagnostics_full2("test_propagate_per_turn1.h5")));
    bunch_train_simulator.add_per_step(1,
            Diagnostics_sptr(
                    new Diagnostics_basic("test_propagate_per_step1.h5")));

    bunch_train_simulator.add_per_turn(2,
            Diagnostics_sptr(
                    new Diagnostics_full2("test_propagate_per_turn2.h5")));
    bunch_train_simulator.add_per_step(2,
            Diagnostics_sptr(
                    new Diagnostics_basic("test_propagate_per_step2.h5")));

    int num_turns = 4;
    propagator.propagate(bunch_train_simulator, num_turns);
    xml_save(propagator, "propagator_train2.xml");

    Propagator loaded;
    xml_load(propagator, "propagator_train.xml");
}


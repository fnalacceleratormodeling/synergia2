#define BOOST_TEST_MAIN
#include <boost/lexical_cast.hpp>
#include <boost/test/unit_test.hpp>
#include "synergia/simulation/propagator.h"
#include "synergia/simulation/lattice_simulator.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/lattice/chef_utils.h"
#include "synergia/bunch/bunch.h"
#include "lattice_fixture.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
#include "synergia/simulation/diagnostics_actions.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

const double tolerance = 1.0e-12;

struct Object_to_sptr_hack
{
    void
    operator()(void const *) const
    {
    }
};

BOOST_FIXTURE_TEST_CASE(construct, Lattice_fixture)
{
    Dummy_collective_operator_sptr space_charge(
            new Dummy_collective_operator("space_charge"));
    Lattice_simulator lattice_simulator(lattice_sptr, 2);

    Split_operator_stepper_sptr stepper_sptr(
            new Split_operator_stepper(lattice_simulator, space_charge, 7));

    Propagator propagator(stepper_sptr);
}

BOOST_FIXTURE_TEST_CASE(propagate, Lattice_fixture)
{
    Dummy_collective_operator_sptr space_charge(
            new Dummy_collective_operator("space_charge"));

    Lattice_simulator lattice_simulator(lattice_sptr, 2);

    Split_operator_stepper_sptr stepper_sptr(
            new Split_operator_stepper(lattice_simulator, space_charge, 4));
    Propagator propagator(stepper_sptr);

    Bunch_sptr bunch_sptr(new Bunch(b.bunch));
    Bunch_simulator bunch_simulator(bunch_sptr);

    Diagnostics_sptr first_step_full2_diag_sptr(
            new Diagnostics_full2("first_full2_per_step.h5"));
    Diagnostics_sptr second_step_full2_diag_sptr(
            new Diagnostics_full2("second_full2_per_step.h5"));
    bunch_simulator.add_per_step(
            first_step_full2_diag_sptr);
    bunch_simulator.add_per_step(
            second_step_full2_diag_sptr);

    Diagnostics_sptr
            first_turn_particles_diag_sptr(
                    new Diagnostics_particles(
                            "first_particles_per_turn.h5"));
    Diagnostics_sptr second_turn_particles_diag_sptr(
            new Diagnostics_particles(
                    "second_particles_per_turn.h5"));
    bunch_simulator.add_per_turn(
            first_turn_particles_diag_sptr);
    bunch_simulator.add_per_turn(
            second_turn_particles_diag_sptr);

    int num_turns = 4;
    propagator.propagate(bunch_simulator, num_turns);

}

BOOST_FIXTURE_TEST_CASE(propagate_max_turns, Lattice_fixture)
{
    Dummy_collective_operator_sptr space_charge(
            new Dummy_collective_operator("space_charge"));

    Lattice_simulator lattice_simulator(lattice_sptr, 2);

    Split_operator_stepper_sptr stepper_sptr(
            new Split_operator_stepper(lattice_simulator, space_charge, 4));
    Propagator propagator(stepper_sptr);

    Bunch_sptr bunch_sptr(new Bunch(b.bunch));
    Bunch_simulator bunch_simulator(bunch_sptr);

    Diagnostics_sptr first_step_full2_diag_sptr(
            new Diagnostics_full2("first_full2_per_step.h5"));
    Diagnostics_sptr second_step_full2_diag_sptr(
            new Diagnostics_full2("second_full2_per_step.h5"));
    bunch_simulator.add_per_step(
            first_step_full2_diag_sptr);
    bunch_simulator.add_per_step(
            second_step_full2_diag_sptr);

    Diagnostics_sptr
            first_turn_particles_diag_sptr(
                    new Diagnostics_particles(
                            "first_particles_per_turn.h5"));
    Diagnostics_sptr second_turn_particles_diag_sptr(
            new Diagnostics_particles(
                    "second_particles_per_turn.h5"));
    bunch_simulator.add_per_turn(
            first_turn_particles_diag_sptr);
    bunch_simulator.add_per_turn(
            second_turn_particles_diag_sptr);

    int num_turns = 4;
    int max_turns = 1;
    propagator.propagate(bunch_simulator, num_turns, max_turns);

    const char second_checkpoint[] = "second_checkpoint";
    propagator.set_checkpoint_dir(second_checkpoint);
    propagator.resume(Propagator::default_checkpoint_dir, false, 0, false, 0);

    propagator.resume(second_checkpoint, true, 2, false, 0);
}

//BOOST_FIXTURE_TEST_CASE(propagate1, Lattice_fixture)
//{
//
//    Dummy_collective_operator_sptr space_charge(new Dummy_collective_operator(
//            "space_charge"));
//    Lattice_simulator lattice_simulator(lattice_sptr, 2);
//
//    Split_operator_stepper_sptr stepper_sptr(new Split_operator_stepper(
//            lattice_simulator, space_charge, 4));
//    Propagator propagator(stepper_sptr);
//
//    Bunch_sptr bunch_sptr(&b.bunch,Object_to_sptr_hack());
//
//    Diagnostics_full2_sptr per_step_diagnostics1(new Diagnostics_full2(bunch_sptr,
//            "test_propagate_per_step1.h5"));
//    Diagnostics_full2_sptr per_step_diagnostics2(new Diagnostics_full2(bunch_sptr,
//            "test_propagate_per_step2.h5"));
//    Multi_diagnostics multi_diagnostics_step;
//    multi_diagnostics_step.append(per_step_diagnostics1);
//    multi_diagnostics_step.append(per_step_diagnostics2);
//
//    Diagnostics_particles_sptr per_turn_diagnostics1(new Diagnostics_particles(bunch_sptr,
//            "test_propagate_per_turn1.h5"));
//    Diagnostics_particles_sptr per_turn_diagnostics2(new Diagnostics_particles(bunch_sptr,
//            "test_propagate_per_turn2.h5"));
//    Multi_diagnostics multi_diagnostics_turn;
//    multi_diagnostics_turn.append(per_turn_diagnostics1);
//    multi_diagnostics_turn.append(per_turn_diagnostics2);
//
//    int num_turns = 4;
//    propagator.propagate(b.bunch, num_turns, multi_diagnostics_step,
//            multi_diagnostics_turn);
//
//}

//BOOST_FIXTURE_TEST_CASE(propagate2, Lattice_fixture)
//{
//    Dummy_collective_operator_sptr space_charge(new Dummy_collective_operator(
//            "space_charge"));
//    Lattice_simulator lattice_simulator(lattice_sptr, 2);
//
//    Split_operator_stepper_sptr stepper_sptr(new Split_operator_stepper(
//            lattice_simulator, space_charge, 4));
//    Propagator propagator(stepper_sptr);
//
//    Bunch_sptr bunch_sptr(new Bunch(b.bunch));
//
//    Diagnostics_particles per_turn_diagnostics(bunch_sptr,
//            "test_propagate_per_turn.h5");
//    Diagnostics_basic per_step_diagnostics(bunch_sptr,
//            "test_propagate_per_step.h5");
//    int num_turns = 4;
//    propagator.propagate(*bunch_sptr, num_turns, per_step_diagnostics,
//            per_turn_diagnostics);
//}
BOOST_FIXTURE_TEST_CASE(serialize_xml, Lattice_fixture)
{
    Dummy_collective_operator_sptr space_charge(
            new Dummy_collective_operator("space_charge"));
    Lattice_simulator lattice_simulator(lattice_sptr, 2);

    Split_operator_stepper_sptr stepper_sptr(
            new Split_operator_stepper(lattice_simulator, space_charge, 7));

    Propagator propagator(stepper_sptr);
    xml_save(propagator, "propagator.xml");

    Bunch_sptr bunch_sptr(new Bunch(b.bunch));
    Bunch_simulator bunch_simulator(bunch_sptr);

    bunch_simulator.add_per_turn(
            Diagnostics_sptr(
                    new Diagnostics_particles(
                            "test_propagate_per_turn.h5")));
    bunch_simulator.add_per_step(
            Diagnostics_sptr(
                    new Diagnostics_basic(
                            "test_propagate_per_step.h5")));

    int num_turns = 4;
    propagator.propagate(bunch_simulator, num_turns);
    xml_save(propagator, "propagator2.xml");

    Propagator loaded;
    xml_load(propagator, "propagator.xml");
}


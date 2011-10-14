#define BOOST_TEST_MAIN
#include <boost/lexical_cast.hpp>
#include <boost/test/unit_test.hpp>
#include "synergia/simulation/propagator.h"
#include "synergia/simulation/lattice_simulator.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/lattice/chef_utils.h"
#include "synergia/bunch/bunch.h"
#include "synergia/bunch/bunch_with_diagnostics.h"
#include "lattice_fixture.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
#include "synergia/simulation/standard_diagnostics_actions.h"


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
    Dummy_collective_operator_sptr space_charge(new Dummy_collective_operator(
            "space_charge"));
    Lattice_simulator lattice_simulator(lattice_sptr, 2);

    Split_operator_stepper_sptr stepper_sptr(new Split_operator_stepper(
            lattice_simulator, space_charge, 7));

    Propagator propagator(stepper_sptr);
}

BOOST_FIXTURE_TEST_CASE(propagate, Lattice_fixture)
 {
  Dummy_collective_operator_sptr space_charge(new Dummy_collective_operator(
            "space_charge"));
    
    Lattice_simulator lattice_simulator(lattice_sptr, 2);

    Split_operator_stepper_sptr stepper_sptr(new Split_operator_stepper(
            lattice_simulator, space_charge, 4));
    Propagator propagator(stepper_sptr);

    Bunch_sptr bunch_sptr(new Bunch(b.bunch));
    Standard_diagnostics_actions_sptr diagnostics_actions_sptr(new Standard_diagnostics_actions);
    Bunch_with_diagnostics bunch_with_diagnostics(bunch_sptr, diagnostics_actions_sptr);
  
    Diagnostics_sptr first_step_full2_diag_sptr(new  Diagnostics_full2(bunch_sptr,"first_full2_per_step.h5"));     
    Diagnostics_sptr second_step_full2_diag_sptr(new  Diagnostics_full2(bunch_sptr,"second_full2_per_step.h5"));
    bunch_with_diagnostics.add_per_step_diagnostics(first_step_full2_diag_sptr);
    bunch_with_diagnostics.add_per_step_diagnostics(second_step_full2_diag_sptr);
      
    Diagnostics_sptr first_turn_particles_diag_sptr(new  Diagnostics_particles(bunch_sptr,"first_particles_per_turn.h5"));
    Diagnostics_sptr second_turn_particles_diag_sptr(new  Diagnostics_particles(bunch_sptr,"second_particles_per_turn.h5"));
    bunch_with_diagnostics.add_per_turn_diagnostics(first_turn_particles_diag_sptr);
    bunch_with_diagnostics.add_per_turn_diagnostics(second_turn_particles_diag_sptr);
       
    int num_turns = 4;
    propagator.propagate(bunch_with_diagnostics, num_turns);
     
 }


BOOST_FIXTURE_TEST_CASE(propagate1, Lattice_fixture)
{

    Dummy_collective_operator_sptr space_charge(new Dummy_collective_operator(
            "space_charge"));
    Lattice_simulator lattice_simulator(lattice_sptr, 2);

    Split_operator_stepper_sptr stepper_sptr(new Split_operator_stepper(
            lattice_simulator, space_charge, 4));
    Propagator propagator(stepper_sptr);

    Bunch_sptr bunch_sptr(&b.bunch,Object_to_sptr_hack());
    
    Diagnostics_full2_sptr per_step_diagnostics1(new Diagnostics_full2(bunch_sptr,
            "test_propagate_per_step1.h5"));
    Diagnostics_full2_sptr per_step_diagnostics2(new Diagnostics_full2(bunch_sptr,
            "test_propagate_per_step2.h5"));
    Multi_diagnostics multi_diagnostics_step;
    multi_diagnostics_step.append(per_step_diagnostics1);
    multi_diagnostics_step.append(per_step_diagnostics2);

    Diagnostics_particles_sptr per_turn_diagnostics1(new Diagnostics_particles(bunch_sptr,
            "test_propagate_per_turn1.h5"));
    Diagnostics_particles_sptr per_turn_diagnostics2(new Diagnostics_particles(bunch_sptr,
            "test_propagate_per_turn2.h5"));
    Multi_diagnostics multi_diagnostics_turn;
    multi_diagnostics_turn.append(per_turn_diagnostics1);
    multi_diagnostics_turn.append(per_turn_diagnostics2);

    int num_turns = 4;
    propagator.propagate(b.bunch, num_turns, multi_diagnostics_step,
            multi_diagnostics_turn);

}

BOOST_FIXTURE_TEST_CASE(propagate2, Lattice_fixture)
{
    Dummy_collective_operator_sptr space_charge(new Dummy_collective_operator(
            "space_charge"));
    Lattice_simulator lattice_simulator(lattice_sptr, 2);

    Split_operator_stepper_sptr stepper_sptr(new Split_operator_stepper(
            lattice_simulator, space_charge, 4));
    Propagator propagator(stepper_sptr);

    Bunch_sptr bunch_sptr(new Bunch(b.bunch));
    
    Diagnostics_particles per_turn_diagnostics(bunch_sptr,
            "test_propagate_per_turn.h5");
    Diagnostics_basic per_step_diagnostics(bunch_sptr,
            "test_propagate_per_step.h5");
    int num_turns = 4;
    propagator.propagate(*bunch_sptr, num_turns, per_step_diagnostics,
            per_turn_diagnostics);
}




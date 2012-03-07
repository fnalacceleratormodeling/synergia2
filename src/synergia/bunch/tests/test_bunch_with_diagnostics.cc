#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <cstdlib>
#include "synergia/foundation/physical_constants.h"
#include "synergia/bunch/bunch.h"
#include "synergia/bunch/diagnostics.h"
#include "synergia/bunch/bunch_with_diagnostics.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
#include "synergia/simulation/tests/bunch_fixture.h"
#include "synergia/simulation/diagnostics_actions.h"
#include "synergia/simulation/diagnostics_actions.cc"
//#include "synergia/simulation/propagate_actions.cc"

BOOST_GLOBAL_FIXTURE(MPI_fixture)

// const double tolerance = 1.0e-14;
//
// const double mass = 100.0;
// const double total_energy = 125.0;
// const int total_num = 100;
// const double real_num = 2.0e12;


BOOST_FIXTURE_TEST_CASE(construct,Bunch_fixture)
{
     Bunch_sptr bunch_sptr(new Bunch(bunch));
     Diagnostics_actions_sptr diagnostics_actions_sptr(new Diagnostics_actions);

      Diagnostics_sptr step_full2_diag_sptr(new  Diagnostics_full2(bunch_sptr,"full2_per_step.h5"));
      diagnostics_actions_sptr->add_per_step(step_full2_diag_sptr);
      Diagnostics_sptr turn_particles_diag_sptr(new  Diagnostics_particles(bunch_sptr,"particles_per_turn.h5"));
      diagnostics_actions_sptr->add_per_turn(turn_particles_diag_sptr);

      Bunch_with_diagnostics bunch_with_diagnostics(bunch_sptr, diagnostics_actions_sptr);
      bunch_with_diagnostics.check_bunch_pointer_in_diagnostics();

      Bunch_sptr bsptr=bunch_with_diagnostics.get_bunch_sptr();
      Diagnostics_actions_sptr  staptr=bunch_with_diagnostics.get_diagnostics_actions_sptr();

}

BOOST_FIXTURE_TEST_CASE(construct1,Bunch_fixture)
{
     Bunch_sptr bunch_sptr(new Bunch(bunch));
     Diagnostics_actions_sptr diagnostics_actions_sptr(new Diagnostics_actions);

      Diagnostics_sptr step_full2_diag_sptr(new  Diagnostics_full2(bunch_sptr,"full2_per_step.h5"));
      Diagnostics_sptr turn_particles_diag_sptr(new  Diagnostics_particles(bunch_sptr,"particles_per_turn.h5"));


      Bunch_with_diagnostics bunch_with_diagnostics(bunch_sptr, diagnostics_actions_sptr);

      bunch_with_diagnostics.add_per_step_diagnostics(step_full2_diag_sptr);
      bunch_with_diagnostics.add_per_turn_diagnostics(turn_particles_diag_sptr);
      bunch_with_diagnostics.check_bunch_pointer_in_diagnostics();

      Bunch_sptr bsptr=bunch_with_diagnostics.get_bunch_sptr();
      Diagnostics_actions_sptr  staptr=bunch_with_diagnostics.get_diagnostics_actions_sptr();

}
// BOOST_FIXTURE_TEST_CASE(construct2,Bunch_fixture)
// {    // it is supposed to fail, i.e. to throw an error
//      Bunch_sptr bunch_sptr(new Bunch(bunch));
//      Bunch_sptr bunch1_sptr(new Bunch(bunch));
//      Diagnostics_actions_sptr diagnostics_actions_sptr(new Diagnostics_actions);
//
//       Diagnostics_sptr step_full2_diag_sptr(new  Diagnostics_full2(bunch_sptr,"full2_per_step.h5"));
//       Diagnostics_sptr turn_particles_diag_sptr(new  Diagnostics_particles(bunch_sptr,"particles_per_turn.h5"));
//
//
//       Bunch_with_diagnostics bunch_with_diagnostics(bunch1_sptr, diagnostics_actions_sptr);
//
//       bunch_with_diagnostics.add_per_step_diagnostics(step_full2_diag_sptr);
//       bunch_with_diagnostics.add_per_turn_diagnostics(turn_particles_diag_sptr);
//
//
// }
